// test ATLAS
// http://math-atlas.sourceforge.net/
// http://www.netlib.org/lapack/
// ALGLIB, http://www.alglib.net/matrixops/lu.php
//
// For ATLAS 3.8.4, link with libcblas.so and libatlas.so
// For ATLAS 3.10.0,link with libtatlas.so

//compile in Debian 7 (wheezy)
//  g++ -I/usr/include/atlas -llapack_atlas -lcblas -latlas useATLAS.cpp
//  g++ -I/usr/include/atlas -llapack_atlas -lblas useATLAS.cpp

#include <sys/time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <iomanip>
#include <math.h>

//#include <atlas/clapack.h>

// In ATLAS 3.10.0, we need this extern "C"
extern "C"
{
#include <cblas.h>
#include <clapack.h>
}

using std::cout;
using std::endl;
using std::cerr;

typedef double **TyMat;
typedef double *TyVec;

void vector_show_row(const TyVec v, size_t size)
{
  cout<<"size: "<< size<<endl;
  for (size_t i = 0; i < size; ++i) {
    cout<<v[i]<<" ";
  }
  cout<<endl<<endl;
}

void matrix_show_square(const TyMat matA, size_t m, size_t n)
{
  cout<<"size: "<<m<<" * "<<n<<endl;
  for (size_t i = 0; i < m; ++i) {
    for (size_t j = 0; j < n; ++j) {
      cout<<matA[i][j]<<" ";
    }
    cout<<endl;
  }
  cout<<endl;
}

void CreateMat(TyMat *pmat, TyVec *pmatp, size_t m, size_t n)
{
  *pmat = (double**)malloc(m*sizeof(double*));
//  *pmatp = (double*)malloc(m*n*sizeof(double));
  *pmatp = (double*)calloc(m*n,sizeof(double));
  for (size_t i=0; i<m; i++) (*pmat)[i] = &(*pmatp)[i*n];
}

double GetTimeNow()
{
  struct timeval tv;
  if (gettimeofday(&tv,NULL)) {
    cerr<<"Fail to get time! (gettimeofday)"<<endl;
    exit(1);
  }
  return tv.tv_sec + tv.tv_usec * 1e-6;
}

double VecAbsAve(const TyVec pA, const TyVec pB, size_t m)
{
  double s=0;
  if (pB) {
    for (size_t i=0; i<m; i++) {
      s += fabs(pB[i] - pA[i]);
    }
  } else {
    for (size_t i=0; i<m; i++) {
      s += fabs(pA[i]);
    }
  }
  return s/m;
}

double MatAbsAve(const TyMat A, const TyMat B, size_t m, size_t n)
{
  double s=0;
  if (B) {
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

double MatAbsMax(const TyMat A, const TyMat B, size_t m, size_t n)
{
  double v, s=0;
  if (B) {
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

void SimpleTranspose(TyMat A, size_t m, size_t n)
{
  double t;
  for (size_t i=0; i<m; i++)
    for (size_t j=i+1; j<n; j++) {
      t = A[i][j];
      A[i][j] = A[j][i];
      A[j][i] = t;
    }
}


extern "C" void ATL_buildinfo(void);

int main(int argc, char *argv[])
{
  size_t m = 4000;  // row
  size_t n = 4000;  // column
  double t0, t1;

  cout<<"Version info for ATLAS:"<<endl;
  // if link fails, you are using ATLAS version older than 3.3.6.
  ATL_buildinfo();
  cout<<endl;

  if (argc>1) {
    m=n=atoi(argv[1]);
  }

  t0 = GetTimeNow();
  cout<<"Allocating memory:"<<endl;
  TyMat matA=NULL;
  TyVec p_matA=NULL;
  CreateMat(&matA, &p_matA, m, n);
  TyMat matB=NULL;
  TyVec p_matB=NULL;
  CreateMat(&matB, &p_matB, m, n);
  TyMat matC=NULL;
  TyVec p_matC=NULL;
  CreateMat(&matC, &p_matC, m, n);
  TyMat matD=NULL;
  TyVec p_matD=NULL;
  CreateMat(&matD, &p_matD, m, n);
  int *ipiv=(int*)malloc(n*sizeof(int));
  t1 = GetTimeNow();
  printf("  time: %5.3f sec.\n", t1-t0);

  size_t A_size = m*n*sizeof(double);
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
//      matA[i][j] = i*i+j*j-i*j+i;  // low rank
      matA[i][j] = (rand()-((RAND_MAX)>>1)) & ((1<<9)-1);
    }
  }

  t1 = GetTimeNow();
  printf("  time: %5.3f sec.\n", t1-t0);
  t0 = GetTimeNow();

  cout<<"copying: B = A"<<endl;
  memcpy(p_matB, p_matA, m*n*sizeof(double));

  t1 = GetTimeNow();
  printf("  time: %5.3f sec.\n", t1-t0);
  t0 = GetTimeNow();

  cout<<"calculating: C = A * B"<<endl;
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, n,
              1.0, p_matA, m, p_matB, m, 0.0, p_matC, m);

  t1 = GetTimeNow();
  printf("  time: %5.3f sec.\n", t1-t0);

  t0 = GetTimeNow();
  {
  size_t n_rhs = 1;
  cout<<"solving: D: A * D = C  "<<"(NRHS = "<<n_rhs<<")"<<endl;
  memcpy(p_matD, p_matC, m*n*sizeof(double));  // make a copy of C, used to do solve
//  clapack_dgesv(CblasRowMajor, n, m, A, )

//  cblas_dtrsm(CblasRowMajor, CblasLeft, CblasLower, CblasNoTrans, CblasNonUnit,
//              m, n)
//  trsv

// using LAPACK
// we have to do a transpose first (you can do this while coping C to D)
  SimpleTranspose(matD, m, n);
//int clapack_♦gesv ( const enum CBLAS ORDER Order, const int N, const int NRHS,
//                     TYPE *A, const int lda, int *ipiv, TYPE *B, const int ldb )
//  RowMajor or CblasColMajor,  order of matrix A, number of right hand sides,
//  A, leading dimension of A, array that records pivoting(just leave empty),
//  result vector b on entry, x on exit,  leading dimension of b,  return INFO
  clapack_dgesv(CblasRowMajor, m, n_rhs,
                p_matA, n, ipiv, p_matD, m);  // use posv for SPD
  // Transpose back...
  SimpleTranspose(matD, m, n);

  t1 = GetTimeNow();
  printf("  time: %5.3f sec.\n", t1-t0);

  cout<<"cal: error (B-D)_Abs"<<endl;
  cout<<"  mean abs error = "<<MatAbsAve(matB,matD,m,n_rhs)<<endl;
  cout<<"  max  abs error = "<<MatAbsMax(matB,matD,m,n_rhs)<<endl;
  cout<<"  abs norm  = "<<MatAbsAve(matB,NULL,m,n_rhs)<<endl;
  }

  memcpy(p_matA, p_matB, m*n*sizeof(double));
  t0 = GetTimeNow();
  {
  size_t n_rhs = m;
  cout<<"solving: D: A * D = C  "<<"(NRHS = "<<n_rhs<<")"<<endl;
  memcpy(p_matD, p_matC, m*n*sizeof(double));  // make a copy of C, used to do solve
//  clapack_dgesv(CblasRowMajor, n, m, A, )

//  cblas_dtrsm(CblasRowMajor, CblasLeft, CblasLower, CblasNoTrans, CblasNonUnit,
//              m, n)
//  trsv

// using LAPACK
// we have to do a transpose first (you can do this while coping C to D)
  SimpleTranspose(matD, m, n);
//int clapack_♦gesv ( const enum CBLAS ORDER Order, const int N, const int NRHS,
//                     TYPE *A, const int lda, int *ipiv, TYPE *B, const int ldb )
//  RowMajor or CblasColMajor,  order of matrix A, number of right hand sides,
//  A, leading dimension of A, array that records pivoting(just leave empty),
//  result vector b on entry, x on exit,  leading dimension of b,  return INFO
  clapack_dgesv(CblasRowMajor, m, n_rhs,
                p_matA, n, ipiv, p_matD, m);  // use posv for SPD
  // Transpose back...
  SimpleTranspose(matD, m, n);

  t1 = GetTimeNow();
  printf("  time: %5.3f sec.\n", t1-t0);

  cout<<"cal: error (B-D)_Abs"<<endl;
  cout<<"  mean abs error = "<<MatAbsAve(matB,matD,m,n_rhs)<<endl;
  cout<<"  max  abs error = "<<MatAbsMax(matB,matD,m,n_rhs)<<endl;
  cout<<"  abs norm  = "<<MatAbsAve(matB,NULL,m,n_rhs)<<endl;
  }

//  matrix_show_square(matA, m, n);
//  matrix_show_square(matB, m, n);
//  matrix_show_square(matC, m, n);
//  matrix_show_square(matD, m, n);
//  cout<<matC[m-1][n-1]<<endl;

  free(ipiv);
  free(p_matD);
  free(matD);
  free(p_matC);
  free(matC);
  free(p_matB);
  free(matB);
  free(p_matA);
  free(matA);

  return 0;
}
