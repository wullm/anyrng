/*******************************************************************************
 * This file is part of AnyRNG.
 * Copyright (c) 2021 Willem Elbers (whe@willemelbers.com)
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/

#include "../include/random.h"
#include <stdio.h>
#include <math.h>

void generate_header(struct sampler *rng, char *fname);

double custom_pdf(double x, void *params) {
  /* Unpack the parameters */
  double *pars = (double *)params;
  double T = pars[0];   // temperature
  double mu = pars[1];  // chemical potential

  /* Calculate the unnormalized Fermi-Dirac density function */
  return (x <= 0.0) ? 0.0 : x * x / (exp((x - mu) / T) + 1.0);
}

int main(int argc, char *argv[]) {
    /* Initialize the inverse transform sampler */
    struct sampler rng;
    double pars[2] = {1.0, 0.0};
    init_sampler(&rng, custom_pdf, 1e-5, 25.0, &pars);

    if (argc < 2) {
        printf("Usage: ./anyrng [header filename]\n");
        return 0;
    }

    /* Dump the tables and an inline rng method to a header file */
    char *fname = argv[1];
    generate_header(&rng, fname);
    printf("Custom header exported to %s.\n", fname);

    /* For testing purposes, seed a xorshift random number generator */
    rng_state seed = rand_uint64_init(12345);

    /* Test the sampler */
    double u = sampleUniform(&seed); //uniform random variate
    double x = draw_sampler(&rng, u); //custom variate

    printf("u = %f\n", u);
    printf("x = F^-1(u) = %f\n", x);

    /* Clean up the sampler struct */
    clean_sampler(&rng);

    return 0;
}

/**
 * @brief Dump the transform tables and an inline rng method to a header file
 *
 * @param rng The #sampler for the custom distribution
 * @param fname File name for the header
 */
void generate_header(struct sampler *rng, char *fname) {
    /* Open the header file */
    FILE *f = fopen(fname, "w");

    /* File description */
    fprintf(f, "/**\n"
               "*  @file %s\n"
               "*  @brief Custom header file generated with AnyRNG by Willem Elbers. Allows \n"
               "*  generating pseudo-random numbers from a predefined distribution using fast \n"
               "*  numerical inversion (Hormann & Leydold, 2003). For more details, refer to\n"
               "*  https://github.com/wullm/AnyRNG.\n"
               "*/\n\n", fname);

    /* Define the structs */
    fprintf(f, "struct spline {\n"
               "  float a0, a1, a2, a3;\n"
               "};\n\n");
    fprintf(f, "struct custom_sampler {\n"
               "  float *endpoints;\n"
               "  struct spline *splines;\n"
               "  int intervalNum;\n"
               "  float *index_table;\n"
               "};\n\n");

    /* Dump the tables */
    fprintf(f, "static float endpoints[%d] = {\n  ", rng->intervalNum);
    for (int i=0; i<rng->intervalNum; i++) {
        fprintf(f, "%e%s%s", rng->intervals[i].Fl, (i < rng->intervalNum-1) ? ", " : "", (i % 5) == 4 ? "\n  " : "");
    }
    fprintf(f, "};\n");
    fprintf(f, "static struct spline splines[%d] = {\n", rng->intervalNum);
    for (int i=0; i<rng->intervalNum; i++) {
        fprintf(f, "  {%e, %e, %e, %e}%s", rng->intervals[i].a0, rng->intervals[i].a1, rng->intervals[i].a2, rng->intervals[i].a3, (i < rng->intervalNum-1) ? ",\n" : "");
    }
    fprintf(f, "};\n");
    fprintf(f, "static float index_table[%d] = {\n  ", SEARCH_TABLE_LENGTH);
    for (int i=0; i<SEARCH_TABLE_LENGTH; i++) {
        fprintf(f, "%e%s%s", rng->index[i], (i < SEARCH_TABLE_LENGTH-1) ? ", " : "", (i % 5) == 4 ? "\n  " : "");
    }
    fprintf(f, "};\n");
    fprintf(f, "static struct custom_sampler transform_rng = {endpoints, splines, %d, index_table};\n\n", rng->intervalNum);

    /* Write the transform method */
    fprintf(f, "/**\n"
               "* @brief Transform a uniform random number into a custom variate X = F^-1(u)\n"
               "*\n"
               "* @param u Random number to be transformed\n"
               "*/\n"
               "static inline double transform_variate(double u) {\n"
               "  /* Use the search table to find a nearby interval */\n"
               "  int tablength = %d;\n"
               "  int intu = (int)(u * tablength);\n"
               "  int idx = transform_rng.index_table[intu < tablength ? intu : tablength - 1];\n"
               "  int i;\n\n"
               "  /* Find the exact interval, i.e. the largest interval such that u > F(p) */\n"
               "  for (i = idx; i < transform_rng.intervalNum-1; i++) {\n"
               "    if (transform_rng.endpoints[i+1] >= u) break;\n"
               "  }\n\n"
               "  float Fl = transform_rng.endpoints[i];\n"
               "  float Fr = transform_rng.endpoints[i+1];\n"
               "  struct spline *iv = &transform_rng.splines[i];\n\n"
               "  /* Evaluate F^-1(u) using the Hermite approximation of F in this interval */\n"
               "  double u_tilde = (u - Fl) / (Fr - Fl);\n"
               "  double H = iv->a0 + iv->a1 * u_tilde + iv->a2 * u_tilde * u_tilde +\n"
               "             iv->a3 * u_tilde * u_tilde * u_tilde;\n\n"
               "  return H;\n"
               "}\n", SEARCH_TABLE_LENGTH);

    /* Close the file */
    fclose(f);
}
