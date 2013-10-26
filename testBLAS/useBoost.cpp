// g++ useBoost.cpp

#include <sys/time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <iomanip>
#include <math.h>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/triangular.hpp>  // for solve()
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>

using std::cout;
using std::endl;
using std::cerr;

using namespace boost::numeric;

double GetTimeNow()
{
  struct timeval tv;
  if (gettimeofday(&tv,NULL)) {
    cerr<<"Fail to get time! (gettimeofday)"<<endl;
    exit(1);
  }
  return tv.tv_sec + tv.tv_usec * 1e-6;
}

double MatAbsAve(const ublas::matrix<double> A, const ublas::matrix<double> B, size_t m, size_t n)
{
  double s=0;
  if (B.size1()>0) {
    for (size_t i=0; i<m; i++) {
      for (size_t j=0; j<n; j++) {
        s += fabs(B(i,j) - A(i,j));
      }
    }
  } else {
    for (size_t i=0; i<m; i++) {
      for (size_t j=0; j<n; j++) {
        s += fabs(A(i,j));
      }
    }
  }
  return s/m/n;
}

double MatAbsMax(const ublas::matrix<double> A, const ublas::matrix<double> B, size_t m, size_t n)
{
  double v, s=0;
  if (B.size1()>0) {
    for (size_t i=0; i<m; i++) {
      for (size_t j=0; j<n; j++) {
        v = fabs(B(i,j) - A(i,j));
        if (v>s) s=v;
      }
    }
  } else {
    for (size_t i=0; i<m; i++) {
      for (size_t j=0; j<n; j++) {
        v = fabs(A(i,j));
        if (v>s) s=v;
      }
    }
  }
  return s;
}

int main(int argc, char *argv[])
{
  size_t m = 1000;  // row
  size_t n = 1000;  // column
  double t0, t1;

  if (argc>1) {
    m=n=atoi(argv[1]);
  }

  t0 = GetTimeNow();
  cout<<"Allocating memory:"<<endl;
  ublas::matrix<double> A(m,n), B(m,n), C(m,n), D(m,n);
  ublas::permutation_matrix<double> P(n);
  ublas::vector<double> v(n), y(n);
  t1 = GetTimeNow();
  printf("  time: %5.3f sec.\n", t1-t0);

  size_t A_size = A.size1()*A.size2()*sizeof(double);
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
      A(i,j) = (rand()-((RAND_MAX)>>1)) & ((1<<9)-1);
    }
  }

  t1 = GetTimeNow();
  printf("  time: %5.3f sec.\n", t1-t0);
  t0 = GetTimeNow();

  cout<<"copying: B = A"<<endl;
  B = A;

  t1 = GetTimeNow();
  printf("  time: %5.3f sec.\n", t1-t0);
  t0 = GetTimeNow();

  cout<<"calculating: C = A * B"<<endl;
  C = prod(A,B);

  t1 = GetTimeNow();
  printf("  time: %5.3f sec.\n", t1-t0);

//  ublas::column(D, 0) = solve(m, v, ublas::lower_tag());
  // bnu::project(m, bnu::range(0,m.size1()), bnu::range(2,3)) = bnu::project(m, bnu::range(0,m.size1()), bnu::range(3,4));

  t0 = GetTimeNow();
  {
  size_t n_rhs = 1;
  cout<<"solving: D: A * D = C  "<<"(NRHS = "<<n_rhs<<")"<<endl;

  lu_factorize(A,P);
  t1 = GetTimeNow();
  printf("  time: %5.3f sec. (DONE LU decomp)\n", t1-t0);
  v = ublas::column(C, 0);
  lu_substitute(A,P,v);
  ublas::column(D, 0) = v;

  t1 = GetTimeNow();
  printf("  time: %5.3f sec.\n", t1-t0);

  cout<<"cal: error (B-D)_Abs"<<endl;
  cout<<"  mean abs error = "<<MatAbsAve(B,D,m,n_rhs)<<endl;
  cout<<"  max  abs error = "<<MatAbsMax(B,D,m,n_rhs)<<endl;
  ublas::matrix<double> ZeroEmpty(0,0);
  cout<<"  abs norm  = "<<MatAbsAve(B,ZeroEmpty,m,n_rhs)<<endl;
  }

  A = B;
  t0 = GetTimeNow();
  {
  size_t n_rhs = m;
  cout<<"solving: D: A * D = C  "<<"(NRHS = "<<n_rhs<<")"<<endl;

  lu_factorize(A,P);
  t1 = GetTimeNow();
  printf("  time: %5.3f sec. (DONE LU decomp)\n", t1-t0);
  for (size_t id_rhs = 0; id_rhs<n_rhs; id_rhs++) {
    v = ublas::column(C, id_rhs);
    lu_substitute(A,P,v);
    ublas::column(D, id_rhs) = v;
  }

  t1 = GetTimeNow();
  printf("  time: %5.3f sec.\n", t1-t0);

  cout<<"cal: error (B-D)_Abs"<<endl;
  cout<<"  mean abs error = "<<MatAbsAve(B,D,m,n_rhs)<<endl;
  cout<<"  max  abs error = "<<MatAbsMax(B,D,m,n_rhs)<<endl;
  ublas::matrix<double> ZeroEmpty(0,0);
  cout<<"  abs norm  = "<<MatAbsAve(B,ZeroEmpty,m,n_rhs)<<endl;
  }

  return 0;
}
