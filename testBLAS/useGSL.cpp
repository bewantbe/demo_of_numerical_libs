// test BLAS, and linear solver in GSL (with standard BLAS)
// compile with (Debian 7 Wheezy)
// g++ useGSL.cpp -lgsl -lblas
// g++ useGSL.cpp -lgsl -lgslcblas

#include <sys/time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <iomanip>
#include <math.h>

#define HAVE_INLINE             // use inline version gsl_vector_get etc.
#define GSL_RANGE_CHECK_OFF
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
//#include <gsl/gsl_cblas.h>
#include <gsl/gsl_linalg.h>

using std::cout;
using std::endl;
using std::cerr;

void gsl_vector_show_row(const gsl_vector *v)
{
  cout<<"size: "<<v->size<<endl;
  for (size_t i = 0; i < v->size; ++i) {
    cout<<gsl_vector_get(v, i)<<" ";
  }
  cout<<endl<<endl;
}

void gsl_matrix_show_square(const gsl_matrix *m)
{
  cout<<"size: "<<m->size1<<" "<<m->size2<<endl;
  for (size_t i = 0; i < m->size1; ++i) {
    for (size_t j = 0; j < m->size2; ++j) {
      cout<<gsl_matrix_get(m, i, j)<<" ";
    }
    cout<<endl;
  }
  cout<<endl;
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

double MatAbsAve(const gsl_matrix *A, const gsl_matrix *B, size_t m, size_t n)
{
  double s=0;
  if (B) {
    for (size_t i=0; i<m; i++) {
      for (size_t j=0; j<n; j++) {
        s += fabs(gsl_matrix_get(B,i,j) - gsl_matrix_get(A,i,j));
      }
    }
  } else {
    for (size_t i=0; i<m; i++) {
      for (size_t j=0; j<n; j++) {
        s += fabs(gsl_matrix_get(A,i,j));
      }
    }
  }
  return s/m/n;
}

double MatAbsMax(const gsl_matrix *A, const gsl_matrix *B, size_t m, size_t n)
{
  double v, s=0;
  if (B) {
    for (size_t i=0; i<m; i++) {
      for (size_t j=0; j<n; j++) {
        v = fabs(gsl_matrix_get(B,i,j) - gsl_matrix_get(A,i,j));
        if (v>s) s=v;
      }
    }
  } else {
    for (size_t i=0; i<m; i++) {
      for (size_t j=0; j<n; j++) {
        v = fabs(gsl_matrix_get(A,i,j));
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

  //cout<<"Version info from GSL:"<<endl;

  if (argc>1) {
    m=n=atoi(argv[1]);
  }

  t0 = GetTimeNow();
  cout<<"Allocating memory:"<<endl;
  gsl_matrix *A = gsl_matrix_alloc(m, n);
  gsl_matrix *B = gsl_matrix_alloc(m, n);
  gsl_matrix *C = gsl_matrix_alloc(m, n);
  gsl_matrix *D = gsl_matrix_alloc(m, n);
  gsl_permutation *p = gsl_permutation_alloc(m);
  gsl_vector *v = gsl_vector_alloc(n);
  gsl_vector *y = gsl_vector_alloc(n);
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
  for (size_t i = 0; i < A->size1; ++i) {
    for (size_t j = 0; j < A->size2; ++j) {
      gsl_matrix_set(A,i,j, (rand()-((RAND_MAX)>>1)) & ((1<<9)-1));
    }
  }

  t1 = GetTimeNow();
  printf("  time: %5.3f sec.\n", t1-t0);
  t0 = GetTimeNow();

  cout<<"copying: B = A"<<endl;
  gsl_matrix_memcpy(B, A);

  t1 = GetTimeNow();
  printf("  time: %5.3f sec.\n", t1-t0);
  t0 = GetTimeNow();

  cout<<"calculating: C = A * B"<<endl;
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, A, B, 0.0, C);  // alias is not allowed

  t1 = GetTimeNow();
  printf("  time: %5.3f sec.\n", t1-t0);

  t0 = GetTimeNow();
  {
  size_t n_rhs = 1;
  cout<<"solving: D: A * D = C  "<<"(NRHS = "<<n_rhs<<")"<<endl;
//  gsl_matrix_memcpy(D, C);  // make a copy of C, used to do solve

  int s;
  gsl_linalg_LU_decomp(A, p, &s);
  t1 = GetTimeNow();
  printf("  time: %5.3f sec. (DONE LU decomp)\n", t1-t0);
  gsl_matrix_get_col(v, C, 0);      // copy one column of D to v
  gsl_linalg_LU_solve(A, p, v, y);  // solve with triangle A
  gsl_matrix_set_col(D, 0, y);      // copy solution y to first column of D

  t1 = GetTimeNow();
  printf("  time: %5.3f sec.\n", t1-t0);

  cout<<"cal: error (B-D)_Abs"<<endl;
  cout<<"  mean abs error = "<<MatAbsAve(B,D,m,n_rhs)<<endl;
  cout<<"  max  abs error = "<<MatAbsMax(B,D,m,n_rhs)<<endl;
  cout<<"  abs norm  = "<<MatAbsAve(B,NULL,m,n_rhs)<<endl;
  }

  gsl_matrix_memcpy(A, B);
  t0 = GetTimeNow();
  {
  size_t n_rhs = m;
  cout<<"solving: D: A * D = C  "<<"(NRHS = "<<n_rhs<<")"<<endl;
//  gsl_matrix_memcpy(D, C);  // make a copy of C, used to do solve

  int s;
  gsl_linalg_LU_decomp(A, p, &s);
  t1 = GetTimeNow();
  printf("  time: %5.3f sec. (DONE LU decomp)\n", t1-t0);
  gsl_vector_view v_view, y_view;
  for (size_t id_rhs=0; id_rhs<n_rhs; id_rhs++) {
//    gsl_matrix_get_col(v, C, id_rhs);  // copy one column of D to v
//    gsl_linalg_LU_solve(A, p, v, y);   // solve with triangle A
//    gsl_matrix_set_col(D, id_rhs, y);  // copy solution y to first column of D
    v_view = gsl_matrix_column(C, id_rhs);  // get column reference of C in v
    y_view = gsl_matrix_column(D, id_rhs);  // get column reference of D in y
    gsl_linalg_LU_solve(A, p, &v_view.vector, &y_view.vector);   // solve with triangle A
  }

  t1 = GetTimeNow();
  printf("  time: %5.3f sec.\n", t1-t0);

  cout<<"cal: error (B-D)_Abs"<<endl;
  cout<<"  mean abs error = "<<MatAbsAve(B,D,m,n_rhs)<<endl;
  cout<<"  max  abs error = "<<MatAbsMax(B,D,m,n_rhs)<<endl;
  cout<<"  abs norm  = "<<MatAbsAve(B,NULL,m,n_rhs)<<endl;
  }

  gsl_matrix_free(A);
  gsl_matrix_free(B);
  gsl_matrix_free(C);
  gsl_matrix_free(D);
  gsl_vector_free(v);
  gsl_vector_free(y);
  gsl_permutation_free(p);

  return 0;
}
