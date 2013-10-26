// matlab like timer

#include <stdio.h>
#include <sys/time.h>

static struct timeval _tv0;
static struct timeval _tv1;

void tic()
{ // rigorous speaking, we should check the return value of gettimeofday()
  gettimeofday(&_tv0, NULL);
}

double toc0()
{
  gettimeofday(&_tv1, NULL);
  return _tv1.tv_sec + 1e-6 * _tv1.tv_usec
      - (_tv0.tv_sec + 1e-6 * _tv0.tv_usec);
}

double toc(const char *st = "")
{
  double t = toc0();
  printf("%s: t = %.6f s.\n", st, t);
  return t;
}

double tocs(const char *st = "\006")
{
  double t = toc0();
  if (st) {
    if (st[0] == '\006' && st[1] == '\0' )
      printf("Elapsed time is %.6f seconds.\n", t);
    else
      printf(st, t);
  }
  return t;
}
