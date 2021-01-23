AnyRNG: A random number generator for arbitrary distributions
=============================================================

Author: Willem Elbers

Allows generating pseudo-random numbers from arbitrary distributions using
fast numerical inversion (Hormann & Leydold, 2003).

Getting started:
----------------

Edit src/anyrng.c to change the probability density function (pdf),
distribution parameters, and the range of the quantile function F^-1. An
example pdf for a Fermi-Dirac distribution is:

```
double custom_pdf(double x, void *params) {
  /* Unpack the parameters */
  double *pars = (double *)params;
  double T = pars[0];   // temperature
  double mu = pars[1];  // chemical potential

  /* Calculate the unnormalized Fermi-Dirac density function */
  return (x <= 0.0) ? 0.0 : x * x / (exp((x - mu) / T) + 1.0);
}
```

The program ensures that the pdf is normalized.  

Once the function has been changed, the program can be compiled and run with

```console
make
./anyrng filename [tolerance]
```

This will produce a stand alone header file that can easily be included
elsewhere. An example program is provided, which can be compiled and run,
after generating the header, with

```console
make example
./example
```
