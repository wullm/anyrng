AnyRNG: A random number generator for arbitrary distributions  {#mainpage}
=============================================================

Author: Willem Elbers

Allows generating pseudo-random numbers from arbitrary distributions using
fast numerical inversion (Hormann & Leydold, 2003).

Getting started:
----------------

Edit src/anyrng.c to change the probability density function (pdf) and
distribution parameters. An example pdf for a Fermi-Dirac distribution is:

double custom_pdf(double x, void *params) {
  /* Unpack the parameters */
  double *pars = (double *)params;
  double T = pars[0];   // temperature
  double mu = pars[1];  // chemical potential

  /* Calculate the unnormalized Fermi-Dirac density function */
  return (x <= 0.0) ? 0.0 : x * x / (exp((x - mu) / T) + 1.0);
}

The program automatically normalizes the pdf. Once the function has been
changed, the program can be compiled and run with

make
./anyrng

This will produce a stand alone header file that can easily be included
elsewhere. An example program is provided, which can be compiled and run
- after generating the header - with

make example
./example
