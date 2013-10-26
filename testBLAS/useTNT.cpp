// use TNT(Template Numerical Toolkit) and JAMA/C++ to preform
// basic matrix operation
// http://math.nist.gov/tnt/jama_doxygen/class_JAMA__LU.html

#include <sys/time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <iomanip>
#include <math.h>

#include <tnt/tnt.h>
#include <jama/jama_lu.h>

using std::cout;
using std::endl;
using std::cerr;

using namespace TNT;

double GetTimeNow()
{
  struct timeval tv;
  if (gettimeofday(&tv,NULL)) {
    cerr<<"Fail to get time! (gettimeofday)"<<endl;
    exit(1);
  }
  return tv.tv_sec + tv.tv_usec * 1e-6;
}

double MatAbsAve(const Array2D<double> A, const Array2D<double> B, size_t m, size_t n)
{
  double s=0;
  if (B.dim1()>0) {
    for (size_t i=0; i<m; i++) {
      for (size_t j=0; j<n; j++) {
        s += fabs(B[i][j] - A[i][j]);
      }
    }
  } else {
    for (size_t i=0; i<m; i++) {
      for (size_t j=0; j<n; j++) {
        s += fabs(A[i][j]);
      }
    }
  }
  return s/m/n;
}

double MatAbsMax(const Array2D<double> A, const Array2D<double> B, size_t m, size_t n)
{
  double v, s=0;
  if (B.dim1()>0) {
    for (size_t i=0; i<m; i++) {
      for (size_t j=0; j<n; j++) {
        v = fabs(B[i][j] - A[i][j]);
        if (v>s) s=v;
      }
    }
  } else {
    for (size_t i=0; i<m; i++) {
      for (size_t j=0; j<n; j++) {
        v = fabs(A[i][j]);
        if (v>s) s=v;
      }
    }
  }
  return s;
}

int main(int argc, char *argv[])
{
  size_t m = 4000;  // row
  size_t n = 4000;  // column
  double t0, t1;

  if (argc>1) {
    m=n=atoi(argv[1]);
  }

  t0 = GetTimeNow();
  cout<<"Allocating memory:"<<endl;
  Array2D<double> A(m,n);
  Array2D<double> B(m,n);
  Array2D<double> C(m,n);
  Array2D<double> D(m,n);
  Array1D<double> v(n);
  Array1D<double> y(n);
  t1 = GetTimeNow();
  printf("  time: %5.3f sec.\n", t1-t0);

  size_t A_size = A.dim1()*A.dim2()*sizeof(double);
  size_t B_size = m*n*sizeof(double);
  size_t C_size = m*n*sizeof(double);
  size_t D_size = m*n*sizeof(double);

  cout<<"Memory cost:"<<endl;
  cout<<"  A("<<m<<"*"<<n<<"): "<<A_size/1024.0/1024.0<<"MB"<<endl;
  cout<<"  B("<<m<<"*"<<n<<"): "<<B_size/1024.0/1024.0<<"MB"<<endl;
  cout<<"  C("<<m<<"*"<<n<<"): "<<C_size/1024.0/1024.0<<"MB"<<endl;
  cout<<"  D("<<m<<"*"<<n<<"): "<<C_size/1024.0/1024.0<<"MB"<<endl;
  cout<<" total: "<<(A_size+B_size+C_size+D_size)/1024.0/1024.0<<"MB"<<endl;
  cout<<endl;


  cout<<"Problem size: "<<m<<" * "<<n<<endl;
  cout<<"generating: A[i][j]"<<endl;
  t0 = GetTimeNow();
  for (size_t i = 0; i < m; ++i) {
    for (size_t j = 0; j < n; ++j) {
      A[i][j] = (rand()-((RAND_MAX)>>1)) & ((1<<9)-1);
    }
  }

  t1 = GetTimeNow();
  printf("  time: %5.3f sec.\n", t1-t0);
  t0 = GetTimeNow();

  cout<<"copying: B = A"<<endl;
  B = A.copy();

  t1 = GetTimeNow();
  printf("  time: %5.3f sec.\n", t1-t0);
  t0 = GetTimeNow();

  cout<<"calculating: C = A * B"<<endl;
  C = matmult(A,B);

  t1 = GetTimeNow();
  printf("  time: %5.3f sec.\n", t1-t0);

  t0 = GetTimeNow();
  {
  size_t n_rhs = 1;
  cout<<"solving: D: A * D = C  "<<"(NRHS = "<<n_rhs<<")"<<endl;

  JAMA::LU<double> lu(A);

  t1 = GetTimeNow();
  printf("  time: %5.3f sec. (DONE LU decomp)\n", t1-t0);

  for (int i=0; i<C.dim1(); i++)
    v[i] = C[i][0];
  y = lu.solve(v);
  for (int i=0; i<C.dim1(); i++)
    D[i][0] = y[i];

  t1 = GetTimeNow();
  printf("  time: %5.3f sec.\n", t1-t0);

  cout<<"cal: error (B-D)_Abs"<<endl;
  cout<<"  mean abs error = "<<MatAbsAve(B,D,m,n_rhs)<<endl;
  cout<<"  max  abs error = "<<MatAbsMax(B,D,m,n_rhs)<<endl;
  Array2D<double> ZeroEmpty(0,0);
  cout<<"  abs norm  = "<<MatAbsAve(B,ZeroEmpty,m,n_rhs)<<endl;
  }

  t0 = GetTimeNow();
  {
  size_t n_rhs = m;
  cout<<"solving: D: A * D = C  "<<"(NRHS = "<<n_rhs<<")"<<endl;

  JAMA::LU<double> lu(A);

  t1 = GetTimeNow();
  printf("  time: %5.3f sec. (DONE LU decomp)\n", t1-t0);

  D = lu.solve(C);

  t1 = GetTimeNow();
  printf("  time: %5.3f sec.\n", t1-t0);

  cout<<"cal: error (B-D)_Abs"<<endl;
  cout<<"  mean abs error = "<<MatAbsAve(B,D,m,n_rhs)<<endl;
  cout<<"  max  abs error = "<<MatAbsMax(B,D,m,n_rhs)<<endl;
  Array2D<double> ZeroEmpty(0,0);
  cout<<"  abs norm  = "<<MatAbsAve(B,ZeroEmpty,m,n_rhs)<<endl;
  }

  return 0;
}
