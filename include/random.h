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

#ifndef RANDOM_H
#define RANDOM_H

/* We use the xoshiro256** pseudo-random number generator */
#include "../include/random_xorshift.h"

#define SEARCH_TABLE_LENGTH 100
#define NUMERICAL_CDF_SAMPLES 1000

/* We allow for arbitrary probability density functions */
typedef double (*pdf)(double x, void *params);

/* PDF for a Fermi-Dirac distribution */
double fd_pdf(double x, void *params);

/* PDF for a Bose-Einstein distribution */
double be_pdf(double x, void *params);

/* A numerical inversion sampler that can be used for arbitrary distributions */
struct sampler {
  /*! The normalization of the pdf */
  double norm;

  /*! The left endpoint of the domain */
  double xl;

  /*! The right endpoint of the domain */
  double xr;

  /*! Pointer to the probability density function */
  pdf f;

  /*! Array of optional parameters passed to the pdf */
  void *params;

  /*! The intervals used in the interpolation */
  struct interval *intervals;

  /*! The number of intervals */
  int intervalNum;

  /*! The indexed search table */
  double *index;
};

/* Intervals used by the numerical inversion sampler */
struct interval {
  int id;
  double l, r;            // endpoints left and right
  double Fl, Fr;          // cdf evaluations at endpoints
  double a0, a1, a2, a3;  // cubic Hermite coefficients
  double error;           // error at midpoint
  int nid;                // the next interval
};

/* Compare intervals by the value of the CDF at the left endpoint */
static inline int compareByLeft(const void *a, const void *b) {
  struct interval *ia = (struct interval *)a;
  struct interval *ib = (struct interval *)b;
  return ia->Fl >= ib->Fl;
}

/* Methods that allow one to sample from arbitrary distribution */
void init_sampler(struct sampler *s, pdf f, double xl, double xr, void *params);
void split_interval(struct sampler *s, int current_interval_id);
void clean_sampler(struct sampler *s);
double draw_sampler(struct sampler *s, double u);



#endif
