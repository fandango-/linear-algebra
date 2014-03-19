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

void d_mat_randtest(d_mat_t mat, flint_rand_t state)
{
    slong r, c, i, j;

    r = mat->r;
    c = mat->c;

    for (i = 0; i < r; i++)
        for (j = 0; j < c; j++)
            * (mat->rows[i] + j) = d_randtest(state);
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
	if (B->entries) {
   long i, j; 

   printf("[");
   for (i = 0; i < B->r; i++) 
   {
      printf("[");
      for (j = 0; j < B->c; j++) 
      { 
         printf("%E", *(B->rows[i] + j)); 
         if (j < B->c-1) printf(" "); 
      }
      printf("]\n");
   }  
   printf("]\n"); 
}
}

void d_swap(double * a, double * b)
{
	if (a != b) {
		double tmp;
		
		tmp = *a;
		*a = *b;
		*b = tmp;
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

/* void d_mat_gso(double ** B, const double ** A, ulong r, ulong c)
{
	slong i, j, k, flag;
	double num, den, mu, dot;
	double t, s;
	
	if((double **) B == (double **) A) {		
		double ** t = d_mat_init(r, c);
		d_mat_gso(t, A, r, c);
		for(i = 0; i < r; i++) {
			for(j = 0; j < c; j++) {
				d_swap(&B[i][j], &t[i][j]);
			}
		}
		d_mat_clear(t);
		return;
	}
	
	if(!r)
	{
		return;
	}
		
	for(k = 0; k < c; k++) {
		for(j = 0; j < r; j++) {
			B[j][k] = A[j][k];
		}
		flag = 1;
		while (flag) {
			t = 0;
			for(i = 0; i < k; i++) {
				s = 0;
				for(j = 0; j < r; j++) {
					s += B[j][i] * B[j][k];
				}
				t += s * s;
				for(j = 0; j < r; j++) {
					B[j][k] -= s * B[j][i];
				}
			}
			s = 0;
			for(j = 0; j < r; j++) {
				s += B[j][k] * B[j][k];
			}
			t += s;
			flag = 0;
			if(s < t) {
				if(s * D_EPS == 0) s = 0;
				else flag = 1;
			}
		}
		s = sqrt(s);
		if(s != 0) s = 1/s;
		for(j = 0; j < r; j++) {
			B[j][k] *= s;
		}
	}
				
	for(i = 0; i < c; i++) {
		for(j = 0; j < r; j++) {
			B[j][i] = A[j][i];
		}
		
		flag = 1;

		while (flag) {
		for(j = 0; j < i; j++) {
			num = B[0][i] * B[0][j];
					 
			for(k = 1; k < r; k++) {
				num += B[k][i] * B[k][j];
			}
			
			den = B[0][j] * B[0][j];
										
			for(k = 1; k < r; k++) {
				den += B[k][j] * B[k][j];
			}
			
			if(fabs(den) > 0)
			{
				mu = num / den;
			
				for(k = 0; k < r; k++) {
					B[k][i] -= mu * B[k][j];
				}
			}
		}
		for (j = 0; j < i; j++) {
			dot = 0;
			for(k = 0; k < r; k++) {
				dot += B[k][i] * B[k][j];
			}
			if(fabs(dot) > 0.5e-12) {
				break;
			}
		}
		if (j == i) {
			flag = 0;
		}
	}
			
		for(j = 0; j < i; j++) {
			num = B[0][i] * B[0][j];
					 
			for(k = 1; k < r; k++) {
				num += B[k][i] * B[k][j];
			}
			
			den = B[0][j] * B[0][j];
										
			for(k = 1; k < r; k++) {
				den += B[k][j] * B[k][j];
			}
			
			if(fabs(den) > 0)
			{
				mu = num / den;
			
				for(k = 0; k < r; k++) {
					B[k][i] -= mu * B[k][j];
				}
			}
		}
		
		for(j = 0; j < i; j++) {
			num = B[0][i] * B[0][j];
					 
			for(k = 1; k < r; k++) {
				num += B[k][i] * B[k][j];
			}
			
			den = B[0][j] * B[0][j];
										
			for(k = 1; k < r; k++) {
				den += B[k][j] * B[k][j];
			}
			
			if(fabs(den) > 0)
			{
				mu = num / den;
			
				for(k = 0; k < r; k++) {
					B[k][i] -= mu * B[k][j];
				}
			}
		}
	}
} */

int main(void) 
{
	// int i;
    FLINT_TEST_INIT(state);
    

    /* flint_printf("gso....");
    fflush(stdout);

	for(i = 0; i < 100 * flint_test_multiplier(); i++) {
        double dot;
        int j, k, l;
        d_mat_t A;

        slong m, n;

        m = n_randint(state, 10);
        n = n_randint(state, 10);
		
        d_mat_init(A, m, n);        
        
        d_mat_randtest(A, state);
                       
        // d_mat_gso(A, (void *) A, m, n);
        
        for(j = 0; j < n; j++) {
			for(k = j + 1; k < n; k++) {
				
				dot = 0;
				for(l = 0; l < m; l++) {
				dot += A[l][j] * A[l][k];
				}
																
				if (fabs(dot) > D_EPS)
				{
					flint_printf("FAIL:\n");
					flint_printf("A:\n");
					d_mat_print(A, m, n);
					printf("%g\n", dot);
					printf("%d %d\n", j, k);
					abort();
				}	
			}
		}
				
        d_mat_clear(A);
	} */
	int i;
	
	for(i = 0; i < 100 * flint_test_multiplier(); i++) {
		int m, n;
		d_mat_t A;
	
		m = n_randint(state, 10);
		n = n_randint(state, 10);
		
		d_mat_init(A, m, n);        
        
		d_mat_randtest(A, state);
				
		d_mat_swap_rows(A, 0, m-1);
		
		d_mat_clear(A);
	}
	
	
	
    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}

