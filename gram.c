#include "flint/fmpz_mat.h"

void fmpz_mat_gram(fmpz_mat_t B, const fmpz_mat_t A)
/*
 *  Sets B to the Gram matrix of the m-dimensional lattice L in 
	n-dimensional Euclidean space R^n spanned by the rows of
	the m × n matrix A 
 *  Requires B to be a m x m matrix, else an exception raised
*/
{
	slong i, j, k;
	
	if(B->r != A->r || B->c != A->r) {
		flint_printf("Exception (fmpz_mat_gram). Incompatible dimensions.\n");
		abort();
	}
	
	if(B == A) {
		fmpz_mat_t t;
		fmpz_mat_init(t, B->r, B->c);
		fmpz_mat_gram(t, A);
		fmpz_mat_swap(B, t);
		fmpz_mat_clear(t);
		return;
	}
	
	if(A->c == 0) {
		fmpz_mat_zero(B);
		return;
	}
	
	for(i = 0; i < B->r; i++) {
		for(j = 0; j < B->c; j++) {
			fmpz_mul(fmpz_mat_entry(B, i, j),
					 fmpz_mat_entry(A, i, 0),
					 fmpz_mat_entry(A, j, 0));
					 
			for (k = 1; k < A->c; k++) {
                fmpz_addmul(fmpz_mat_entry(B, i, j),
                            fmpz_mat_entry(A, i, k),
                            fmpz_mat_entry(A, j, k));
            }
		}
	}
}

int main(void)
{
	fmpz_mat_t A, B, C, D;
    slong i;
    FLINT_TEST_INIT(state);

    flint_printf("gram....");
    fflush(stdout);

    

    for (i = 0; i < 100 * 10; i++)
    {
        slong m, n;

        m = n_randint(state, 50);
        n = n_randint(state, 50);

        fmpz_mat_init(A, m, n);
        fmpz_mat_init(B, n, m);
        fmpz_mat_init(C, m, m);
        fmpz_mat_init(D, m, m);

        fmpz_mat_randtest(A, state, n_randint(state, 200) + 1);
        
        fmpz_mat_transpose(B, A);
        fmpz_mat_mul(C, A, B);
        fmpz_mat_gram(D, A);

        if (!fmpz_mat_equal(C, D))
        {
            flint_printf("FAIL: results not equal\n");
            abort();
        }

        fmpz_mat_clear(A);
        fmpz_mat_clear(B);
        fmpz_mat_clear(C);
        fmpz_mat_clear(D);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");    
    return 0;
}
