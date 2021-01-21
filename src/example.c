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

/* This custom header needs to be generated first with the main program! */
#include "../header.h"

/* A xorshift random number generator */
#include "../include/random_xorshift.h"

/* Standard headers */
#include <stdio.h>
#include <sys/time.h>

int main() {
    /* Seed a xorshift random number generator */
    rng_state seed = rand_uint64_init(10124);

    /* Start the timer */
    struct timeval time_stop, time_start;
    gettimeofday(&time_start, NULL);

    /* Generate a million random numbers */
    double tot = 0;
    long long num = 1000000;
    for (int i=0; i<num; i++) {
        /* Generate a uniform random number */
        double u = sampleUniform(&seed);
        /* Transform to custom distribution */
        double x = transform_variate(u);
        tot += x;
    }

    printf("Mean: %e\n", tot/num);

    /* End the timer */
    gettimeofday(&time_stop, NULL);
    long unsigned microsec = (time_stop.tv_sec - time_start.tv_sec) * 1000000
                           + time_stop.tv_usec - time_start.tv_usec;
    printf("\nTime elapsed: %.5f s\n", microsec/1e6);

    return 0;
}
