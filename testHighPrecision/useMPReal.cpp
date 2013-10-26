// test speed of high precision lib
// MPFR C++
//   http://www.holoborodko.com/pavel/mpfr/
// compile command:
//   g++ useMPReal.cpp -lmpfr

#include <iostream>
#include "mpreal.h"

void mpfr_free_cache (void);

#include <sys/time.h>
static struct timeval _tv;
void tic()
{
  gettimeofday(&_tv, NULL);
}
double toc(const char *st=NULL)
{
  struct timeval tv1;
  gettimeofday(&tv1, NULL);
  double t = (double)tv1.tv_sec + 1e-6 * tv1.tv_usec
             - ((double)_tv.tv_sec + 1e-6 * _tv.tv_usec);
  if (st) {
    printf(st, t);
  }
  return t;
}

int main(int argc, char* argv[])
{
//  using mpfr::mpreal;
//  using std::cout;
//  using std::endl;
//
//  // Required precision of computations in decimal digits
//  // Play with it to check different precisions
//  const int digits = 50;
//
//  // Setup default precision for all subsequent computations
//  // MPFR accepts precision in bits - so we do the conversion
//  mpreal::set_default_prec(mpfr::digits2bits(digits));
//
//  // Compute all the vital characteristics of mpreal (in current precision)
//  // Analogous to lamch from LAPACK
//  const mpreal one         =    1.0;
//  const mpreal zero        =    0.0;
//  const mpreal eps         =    std::numeric_limits<mpreal>::epsilon();
//  const int    base        =    std::numeric_limits<mpreal>::radix;
//  const mpreal prec        =    eps * base;
//  const int bindigits      =    std::numeric_limits<mpreal>::digits(); // eqv. to mpfr::mpreal::get_default_prec();
//  const mpreal rnd         =    std::numeric_limits<mpreal>::round_error();
//  const mpreal maxval      =    std::numeric_limits<mpreal>::max();
//  const mpreal minval      =    std::numeric_limits<mpreal>::min();
//  const mpreal small       =    one / maxval;
//  const mpreal sfmin       =    (small > minval) ? small * (one + eps) : minval;
//  const mpreal round       =    std::numeric_limits<mpreal>::round_style();
//  const int    min_exp     =    std::numeric_limits<mpreal>::min_exponent;
//  const mpreal underflow   =    std::numeric_limits<mpreal>::min();
//  const int    max_exp     =    std::numeric_limits<mpreal>::max_exponent;
//  const mpreal overflow    =    std::numeric_limits<mpreal>::max();
//
//  // Additionally compute pi with required accuracy - just for fun :)
//  const mpreal pi          =    mpfr::const_pi();
//
//  cout.precision(digits);    // Show all the digits
//  cout << "pi         =    "<<    pi          << endl;
//  cout << "eps        =    "<<    eps         << endl;
//  cout << "base       =    "<<    base        << endl;
//  cout << "prec       =    "<<    prec        << endl;
//  cout << "b.digits   =    "<<    bindigits   << endl;
//  cout << "rnd        =    "<<    rnd         << endl;
//  cout << "maxval     =    "<<    maxval      << endl;
//  cout << "minval     =    "<<    minval      << endl;
//  cout << "small      =    "<<    small       << endl;
//  cout << "sfmin      =    "<<    sfmin       << endl;
//  cout << "1/sfmin    =    "<<    1 / sfmin   << endl;
//  cout << "round      =    "<<    round       << endl;
//  cout << "max_exp    =    "<<    max_exp     << endl;
//  cout << "min_exp    =    "<<    min_exp     << endl;
//  cout << "underflow  =    "<<    underflow   << endl;
//  cout << "overflow   =    "<<    overflow    << endl;
//
//  cout << pi.toDouble() <<endl;

  using mpfr::mpreal;
  using std::cout;
  using std::endl;

  int digits = 100000;
  mpreal::set_default_prec(mpfr::digits2bits(digits));
  cout << "set precision: " << digits << " digits"  << endl;
  tic();
  mpreal pi = mpfr::const_pi();
  mpreal e  = exp(mpreal(1));
  toc("pi and e time: %9.6f\n");

  mpreal u,y,z;

  cout << "Initializing random number generator." << endl;
  tic();
  mpfr::random((unsigned)time(NULL));
  toc("%9.6f sec.\n");

  cout << "get a uniformly distributed random number." << endl;
  tic();
  u = mpfr::random();
  toc("%9.6f sec.");

  tic();
//  z = pow(e,1/pi);   // slow!  139.521635sec
  z = exp(-log(e)/pi);
  toc("z time: %9.6f\n");
  cout << "z = " << z.toDouble() <<endl;

  tic();
//  z = pow(e,1/pi);   // slow!  139.521635sec
  y = exp(log(sqrt(e))/sqrt(pi));
  toc("y time: %9.6f\n");
  cout << "y = " << y.toDouble() <<endl;

  tic();
  u = z-y;
  toc("%9.6f sec: u = z-y;\n");

  tic();
  u = z*y;
  toc("%9.6f sec: u = z*y;\n");

  tic();
  u = 1/y;
  toc("%9.6f sec: u = 1/y;\n");

  tic();
  u = z/y;
  toc("%9.6f sec: u = z/y;\n");

  tic();
  u = sqrt(z);
  toc("%9.6f sec: u = sqrt(z);\n");

  tic();
  u = rec_sqrt(z);
  toc("%9.6f sec: u = 1/sqrt(z);\n");

  tic();
  u = agm(z,y);
  toc("%9.6f sec: u = agm(z,y);\n");

  tic();
  u = log(z);
  toc("%9.6f sec: u = log(z);\n");

  tic();
  u = log2(z);
  toc("%9.6f sec: u = log2(z);\n");

  tic();
  u = log10(z);
  toc("%9.6f sec: u = log10(z);\n");

  tic();
  u = exp(z);
  toc("%9.6f sec: u = exp(z);\n");

  tic();
  u = exp2(z);
  toc("%9.6f sec: u = exp2(z);\n");

  tic();
  u = exp10(z);
  toc("%9.6f sec: u = exp10(z);\n");

  tic();
  u = sin(z);
  toc("%9.6f sec: u = sin(z);\n");

  tic();
  u = cos(z);
  toc("%9.6f sec: u = cos(z);\n");

  tic();
  u = tan(z);
  toc("%9.6f sec: u = tan(z);\n");

  tic();
  u = atan(z);
  toc("%9.6f sec: u = atan(z);\n");

  tic();
  u = asin(z);
  toc("%9.6f sec: u = asin(z);\n");

  tic();
  u = acos(z);
  toc("%9.6f sec: u = acos(z);\n");

  tic();
  u = atan2(z, y);
  toc("%9.6f sec: u = atan2(z, y);\n");

  tic();
  u = sinh(z);
  toc("%9.6f sec: u = sinh(z);\n");

  tic();
  u = cosh(z);
  toc("%9.6f sec: u = cosh(z);\n");

  tic();
  u = tanh(z);
  toc("%9.6f sec: u = tanh(z);\n");

  tic();
  u = atanh(z);
  toc("%9.6f sec: u = atanh(z);\n");

  tic();
  u = asinh(z);
  toc("%9.6f sec: u = asinh(z);\n");

  tic();
  u = acosh(y);
  toc("%9.6f sec: u = acosh(y);\n");

  // now test slow functions
  digits = 10000;
  mpreal::set_default_prec(mpfr::digits2bits(digits));
  cout << "set precision: " << digits << " digits" << endl;
  z.setPrecision(mpfr::digits2bits(digits));
  y.setPrecision(mpfr::digits2bits(digits));

  // slow
  tic();
  u = besselj0(z);
  toc("%9.6f sec: u = besselj0(z);\n");

  tic();
  u = besselj1(z);
  toc("%9.6f sec: u = besselj1(z);\n");

  tic();
  u = besseljn(1, z);
  toc("%9.6f sec: u = besseljn(1, z);\n");

  u = bessely0(z);
  toc("%9.6f sec: u = bessely0(z);\n");

  tic();
  u = bessely1(z);
  toc("%9.6f sec: u = bessely1(z);\n");

  tic();
  u = besselyn(1, z);
  toc("%9.6f sec: u = besselyn(1, z);\n");

  tic();
  u = ai(z);
  toc("%9.6f sec: u = Airy_Ai(z);\n");

  tic();
  u = erf(z);
  toc("%9.6f sec: u = erf(z);\n");

  tic();
  u = eint(z);
  toc("%9.6f sec: u = eint(z);\n");

  // now test very slow functions
  digits = 1000;
  mpreal::set_default_prec(mpfr::digits2bits(digits));
  cout << "set precision: " << digits << " digits" << endl;
  z.setPrecision(mpfr::digits2bits(digits));
  y.setPrecision(mpfr::digits2bits(digits));

  // very slow
  tic();
  u = gamma(z);
  toc("%9.6f sec: u = gamma(z);\n");

  tic();
  u = lngamma(z);
  toc("%9.6f sec: u = lngamma(z);\n");

  tic();
  u = digamma(z);
  toc("%9.6f sec: u = digamma(z);\n");

  tic();
  u = zeta(z);
  toc("%9.6f sec: u = zeta(z);\n");

  mpfr_free_cache();
  return 0;
}
