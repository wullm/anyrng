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

#include <stdlib.h>
#include <math.h>

#include "../include/random.h"

/* Generate standard normal variable with Box-Mueller */
double sampleNorm(rng_state *state) {
    /* Generate random integers */
    const uint64_t A = rand_uint64(state);
    const uint64_t B = rand_uint64(state);
    const double RMax = (double) UINT64_MAX + 1;

    /* Map the random integers to the open (!) unit interval */
    const double u = ((double) A + 0.5) / RMax;
    const double v = ((double) B + 0.5) / RMax;

    /* Map to two Gaussians (the second is not used - inefficient) */
    const double z0 = sqrt(-2 * log(u)) * cos(2 * M_PI * v);
    //double z1 = sqrt(-2 * log(u)) * sin(2 * M_PI * v);

    return z0;
}

/**
 * @brief Numerical evaluation of the cumulative distribution function
 *
 * @param xl Left endpoint of the integration
 * @param xr Right endpoint of the integration
 * @param f Probability density function reference
 * @param params Parameters for the distribution function
 */
double numerical_cdf(double xl, double xr, pdf f, void *params) {
  /* Midpoint rule integration */
  double out = 0.0;
  double delta = (xr - xl) / NUMERICAL_CDF_SAMPLES;
  for (int i = 0; i < NUMERICAL_CDF_SAMPLES; i++) {
    double x = xl + (i + 0.5) * delta;
    out += delta * f(x, params);
  }
  return out;
}

/**
 * @brief Initialize the numerical inversion sampler.
 *
 * @param s The #sampler to initialize
 * @param pdf Function reference of the probability density function
 * @param xl Left endpoint of the domain
 * @param xr Right endpoint of the domain
 * @param params Parameters to be passed to the pdf
 *
 * We will compute Hermite polynomial approximations of the cdf F(X) in
 * discrete intervals, which are then used to quickly evaluate the inverse
 * transform X = F^-1(u) of a uniform random variate u.
 */
void init_sampler(struct sampler *s, pdf f, double xl, double xr,
                  void *params) {
  /* Store the parameters and endpoints */
  s->xl = xl;
  s->xr = xr;
  s->f = f;
  s->params = params;

  /* Normalization of the pdf */
  s->norm = 1.0 / numerical_cdf(xl, xr, f, params);

  /* Create the intervals, starting with just one */
  s->intervalNum = 1;
  s->intervals = malloc(s->intervalNum * sizeof(struct interval));

  /* Initially, the first interval covers the entire domain */
  s->intervals[0].id = 0;
  s->intervals[0].l = xl;
  s->intervals[0].r = xr;
  s->intervals[0].Fl = 0.0;
  s->intervals[0].Fr = 1.0;

  /* The current interval under consideration */
  int current_interval_id = 0;

  /* Split intervals until they are small enough */
  char done = 0;
  while (!done) {
    /* The interval under consideration */
    struct interval *iv = &s->intervals[current_interval_id];

    /* Check if the interval is too big (covers more than 5%) */
    if (iv->Fr - iv->Fl > 0.05) {
      /* Split the interval in half */
      split_interval(s, current_interval_id);

    } else if (iv->r >= xr) {
      /* Stop if we are at the end */
      done = 1;
    } else {
      /* Move on to the next interval */
      current_interval_id = iv->nid;
    }
  }

  /* Return to the first interval */
  current_interval_id = 0;

  /* Now calculate Hermite polynomials in intervals and split them up if
   * they are not monotonic or if the error is too big. */
  done = 0;
  while (!done) {
    /* The interval under consideration */
    struct interval *iv = &s->intervals[current_interval_id];

    /* Evaluate the normalized pdf at the endpoints */
    double fl = s->norm * f(iv->l, params);
    double fr = s->norm * f(iv->r, params);

    /* Calculate the cubic Hermite approximation */
    iv->a0 = iv->l;
    iv->a1 = (iv->Fr - iv->Fl) / fl;
    iv->a2 = 3 * (iv->r - iv->l) - (iv->Fr - iv->Fl) * (2. / fl + 1. / fr);
    iv->a3 = 2 * (iv->l - iv->r) + (iv->Fr - iv->Fl) * (1. / fl + 1. / fr);

    /* Evaluate the error at the midpoint */
    double u = 0.5 * (iv->Fr + iv->Fl);
    double H = iv->a0 + iv->a1 * 0.5 + iv->a2 * 0.25 + iv->a3 * 0.125;
    iv->error = fabs(s->norm * numerical_cdf(xl, H, f, params) - u);

    /* Monotonicity check */
    double delta = (iv->Fr - iv->Fl) / (iv->r - iv->l);
    char monotonic = (delta <= 3 * fl) && (delta <= 3 * fr);

    /* If the error is too big or if the polynomial is not monotonic */
    if (iv->error > 1e-6 || !monotonic) {
      /* Split the interval in half */
      split_interval(s, current_interval_id);

    } else if (iv->r >= xr) {
      /* Stop if we are at the end */
      done = 1;
    } else {
      /* Move on to the next interval */
      current_interval_id = iv->nid;
    }
  }

  /* Sort the intervals */
  qsort(s->intervals, s->intervalNum, sizeof(struct interval), compareByLeft);

  /* Allocate memory for the search table */
  s->index = (double *)malloc(SEARCH_TABLE_LENGTH * sizeof(double));

  /* Generate the index search table */
  for (int i = 0; i < SEARCH_TABLE_LENGTH; i++) {
    double u = (double)i / SEARCH_TABLE_LENGTH;

    /* Find the largest interval such that u > F(p) */
    double maxJ = 0;
    int int_i = 0;
    for (int j = 0; j < s->intervalNum; j++) {
      if (s->intervals[j].Fr < u && s->intervals[j].r > maxJ) {
        maxJ = s->intervals[j].r;
        int_i = j;
      }
    }
    s->index[i] = int_i;
  }
}

/**
 * @brief Split an interval in half and update the links
 *
 * @param s The #sampler containing the interval
 * @param current_interval_id Id of the interval to be split
 */
void split_interval(struct sampler *s, int current_interval_id) {
  /* The current interval that will be halved */
  struct interval *iv = &s->intervals[current_interval_id];

  /* Split the interval in half */
  double m = iv->l + 0.5 * (iv->r - iv->l);
  double Fm = s->norm * numerical_cdf(s->xl, m, s->f, s->params);

  /* ID of the new interval */
  int id = s->intervalNum;
  s->intervalNum++;

  /* Allocate memory for the new interval (inefficient, but not critical) */
  s->intervals =
      realloc(s->intervals, s->intervalNum * sizeof(struct interval));

  if (s->intervals == NULL) {
    // error("Error reallocating memory for the intervals.");
  }

  /* Update the pointer to the current interval */
  iv = &s->intervals[current_interval_id];

  /* Insert the right half as a new interval */
  s->intervals[id].id = id;
  s->intervals[id].l = m;
  s->intervals[id].r = iv->r;
  s->intervals[id].Fl = Fm;
  s->intervals[id].Fr = iv->Fr;
  s->intervals[id].nid = iv->nid;  // link to the old interval's right-neighbour

  /* Update the old interval to cover just the left half */
  iv->r = m;
  iv->Fr = Fm;
  iv->nid = id;  // link the left-half to the right-half
}

/**
 * @brief Clean up the inversion sampler
 *
 * @param s The #sampler to be cleaned
 */
void clean_sampler(struct sampler *s) {
  free(s->intervals);
  free(s->index);
}

/**
 * @brief Transform a uniform random number into a custom variate X = F^-1(u)
 *
 * @param s The #sampler for the distribution
 * @param u Random number to be transformed
 */
double draw_sampler(struct sampler *s, double u) {
  /* Use the search table to find a nearby interval */
  int If = floor(u * SEARCH_TABLE_LENGTH);
  int idx = s->index[If < SEARCH_TABLE_LENGTH ? If : SEARCH_TABLE_LENGTH - 1];

  /* Find the exact interval, i.e. the largest interval such that u > F(p) */
  double maxJ = 0;
  int int_i = idx;
  for (int i = idx; i < s->intervalNum; i++) {
    if (s->intervals[i].Fr < u && s->intervals[i].r > maxJ) {
      maxJ = s->intervals[i].r;
      int_i = i;
    } else {
      break;
    }
  }

  /* The correct interval */
  struct interval *iv = &s->intervals[int_i];

  /* Evaluate F^-1(u) using the Hermite approximation of F in this interval */
  double u_tilde = (u - iv->Fl) / (iv->Fr - iv->Fl);
  double H = iv->a0 + iv->a1 * u_tilde + iv->a2 * u_tilde * u_tilde +
             iv->a3 * u_tilde * u_tilde * u_tilde;

  return H;
}


// double draw_sampler_short(struct sampler_short *s, double u) {
  // /* Use the search table to find a nearby interval */
  // int If = floor(u * SEARCH_TABLE_LENGTH);
  // int idx = s->index[If < SEARCH_TABLE_LENGTH ? If : SEARCH_TABLE_LENGTH - 1];
  //
  // /* Find the exact interval, i.e. the largest interval such that u > F(p) */
  // double maxJ = 0;
  // int int_i = idx;
  // for (int i = idx; i < s->intervalNum; i++) {
  //   if (s->intervals[i].Fr < u && s->intervals[i].r > maxJ) {
  //     maxJ = s->intervals[i].r;
  //     int_i = i;
  //   } else {
  //     break;
  //   }
  // }
  //
  // /* The correct interval */
  // struct interval_short *iv = &s->intervals[int_i];
  //
  // /* Evaluate F^-1(u) using the Hermite approximation of F in this interval */
  // double u_tilde = (u - iv->Fl) / (iv->Fr - iv->Fl);
  // double H = iv->a0 + iv->a1 * u_tilde + iv->a2 * u_tilde * u_tilde +
  //            iv->a3 * u_tilde * u_tilde * u_tilde;
  //
  // return H;
// }
