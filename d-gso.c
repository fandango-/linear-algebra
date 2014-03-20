#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "flint/flint.h"
#include "flint/ulong_extras.h"
#include "flint/double_extras.h"
#include "test_helpers.c"

typedef struct
{
	double * entries;
	slong r;
	slong c;
	double ** rows;
} d_mat_struct;

typedef d_mat_struct d_mat_t[1];

#define d_mat_entry(mat,i,j) (*((mat)->rows[i] + (j)))

void _d_vec_add(double * r1, double * r2, double * r3, ulong n)
{
   ulong i;
   for (i = 0; i < n; i++)
	  r1[i] = r2[i] + r3[i];
}

void _d_vec_sub(double * r1, double * r2, double * r3, ulong n)
{
   ulong i;
   for (i = 0; i < n; i++)
	  r1[i] = r2[i] - r3[i];
}

double _d_vec_scalar_product(double * vec1, double * vec2, ulong n)
{
  double sum;

  sum = vec1[0] * vec2[0];
  long i;
  for (i = 1; i < n; i++)
     sum += vec1[i] * vec2[i];

  return sum;
} 

double _d_vec_norm(double * vec, ulong n)
{
  double sum;

  sum = vec[0] * vec[0];
  long i;
  for (i = 1 ; i < n ; i++)
     sum += vec[i] * vec[i];

  return sum;
}

void _d_vec_set(double * vec1, const double * vec2, slong len2)
{
	if (vec1 != vec2) {
		slong i;
		for (i = 0; i < len2; i++)
			vec1[i] = vec2[i];
	}
} 

void d_mat_randtest(d_mat_t mat, flint_rand_t state)
{
    slong r, c, i, j;

    r = mat->r;
    c = mat->c;

    for (i = 0; i < r; i++)
        for (j = 0; j < c; j++)
            d_mat_entry(mat, i, j) = d_randtest(state);
}

void d_mat_init(d_mat_t mat, slong rows, slong cols)
{
    if ((rows) && (cols))
    {
        slong i;
        mat->entries = flint_calloc(rows * cols, sizeof(double));
        mat->rows = flint_malloc(rows * sizeof(double *));
        
        for (i = 0; i < rows; i++)
            mat->rows[i] = mat->entries + i * cols;
    }
    else
        mat->entries = NULL;

    mat->r = rows;
    mat->c = cols;
}

void d_mat_clear(d_mat_t mat)
{
    if (mat->entries)
    {
        flint_free(mat->entries);
        flint_free(mat->rows);
    }
}

void d_mat_print(d_mat_t B) 
{
   long i, j; 

   flint_printf("[");
   for (i = 0; i < B->r; i++) 
   {
      flint_printf("[");
      for (j = 0; j < B->c; j++) 
      { 
         flint_printf("%E", d_mat_entry(B, i, j)); 
         if (j < B->c-1) flint_printf(" "); 
      }
      flint_printf("]\n");
   }  
   flint_printf("]\n"); 
}

void d_mat_swap(d_mat_t mat1, d_mat_t mat2)
{
    if (mat1 != mat2)
    {
        d_mat_struct tmp;

        tmp = *mat1;
        *mat1 = *mat2;
        *mat2 = tmp;
    }
}

void d_mat_set(d_mat_t mat1, const d_mat_t mat2)
{
    if (mat1 != mat2)
    {
        slong i;

        if (mat2->r && mat2->c)
            for (i = 0; i < mat2->r; i++)
                _d_vec_set(mat1->rows[i], mat2->rows[i], mat2->c);
    }
}

void d_mat_swap_rows(d_mat_t mat, slong r, slong s)
{
	if (mat->entries) {
    if (r != s)
    {
        double * u;
        
        u = mat->rows[s];
        mat->rows[s] = mat->rows[r];
        mat->rows[r] = u; 
    }
}
}

void d_mat_gso(d_mat_t B, const d_mat_t A)
{
	slong i, j, k, flag;
	double t, s;
	
	if (B->r != A->r || B->c != A->c) {
		flint_printf("Exception (d_mat_gso). Incompatible dimensions.\n");
        abort();
	}
	
	if (B == A) {		
		d_mat_t t;
		d_mat_init(t, A->r, A->c);
		d_mat_gso(t, A);
		d_mat_swap(B, t);
		d_mat_clear(t);
		return;
	}
	
	if (!A->r)
	{
		return;
	}
		
	for (k = 0; k < A->c; k++) {
		for (j = 0; j < A->r; j++) {
			d_mat_entry(B, j, k) = d_mat_entry(A, j, k);
		}
		flag = 1;
		while (flag) {
			t = 0;
			for (i = 0; i < k; i++) {
				s = 0;
				for (j = 0; j < A->r; j++) {
					s += d_mat_entry(B, j, i) * d_mat_entry(B, j, k);
				}
				t += s * s;
				for (j = 0; j < A->r; j++) {
					d_mat_entry(B, j, k) -= s * d_mat_entry(B, j, i);
				}
			}
			s = 0;
			for (j = 0; j < A->r; j++) {
				s += d_mat_entry(B, j, k) * d_mat_entry(B, j, k);
			}
			t += s;
			flag = 0;
			if (s < t) {
				if (s * D_EPS == 0) s = 0;
				else flag = 1;
			}
		}
		s = sqrt(s);
		if(s != 0) s = 1/s;
		for(j = 0; j < A->r; j++) {
			d_mat_entry(B, j, k) *= s;
		}
	}				
}

int main(void) 
{
	int i;
    FLINT_TEST_INIT(state);
    

    flint_printf("gso....");
    fflush(stdout);

	for (i = 0; i < 100 * flint_test_multiplier(); i++) {
        double dot;
        int j, k, l;
        d_mat_t A;

        slong m, n;

        m = n_randint(state, 10);
        n = n_randint(state, 10);
		
        d_mat_init(A, m, n);        
        
        d_mat_randtest(A, state);
                       
        d_mat_gso(A, A);
                
        for (j = 0; j < n; j++) {
			for (k = j + 1; k < n; k++) {
				
				dot = 0;
				for (l = 0; l < m; l++) {
				dot += d_mat_entry(A, l, j) * d_mat_entry(A, l, k);
				}
																
				if (fabs(dot) > D_EPS)
				{
					flint_printf("FAIL:\n");
					flint_printf("A:\n");
					d_mat_print(A);
					flint_printf("%g\n", dot);
					flint_printf("%d %d\n", j, k);
					abort();
				}	
			}
		}
				
        d_mat_clear(A);
	}
	
	FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}

