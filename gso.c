#include "flint/fmpq_mat.h"
#include "test_helpers.c"

void fmpq_mat_gso(fmpq_mat_t B, const fmpq_mat_t A)
/* Input: A basis a1, ...., an of R^n as the columns of an m x n matrix A
 * Output: An orthogonal basis of R^n as the columns of m x n matrix B
*/
{
	slong i, j, k;
	fmpq_t num, den, mu;
	fmpq_init(num);
	fmpq_init(den);
	fmpq_init(mu);
	
	if(B->r != A->r || B->c != A->c) {
		flint_printf("Exception (fmpq_mat_gso). Incompatible dimensions.\n");
        abort();
	}
	
	if(B == A) {
		fmpq_mat_t t;
		fmpq_mat_init(t, B->r, B->c);
		fmpq_mat_gso(t, A);
		for(i = 0; i < B->r; i++) {
			for(j = 0; j < B->c; j++) {
				fmpq_swap(fmpq_mat_entry(B, i, j), 
                          fmpq_mat_entry(t, i, j));
			}
		}
		fmpq_mat_clear(t);
	}
				
	for(i = 0; i < A->c; i++) {
		for(j = 0; j < A->r; j++) {
			fmpq_set(fmpq_mat_entry(B, j, i),
					 fmpq_mat_entry(A, j, i));
		}

		for(j = 0; j < i; j++) {
			fmpq_zero(num);
					 
			for(k = 0; k < A->r; k++) {
				fmpq_addmul(num,
							fmpq_mat_entry(A, k, i),
							fmpq_mat_entry(B, k, j));
			}
			
			fmpq_zero(den);
										
			for(k = 0; k < A->r; k++) {
				fmpq_addmul(den,
							fmpq_mat_entry(B, k, j),
							fmpq_mat_entry(B, k, j));
			}
			
			if(!fmpq_is_zero(den))
			{
				fmpq_div(mu, num, den);
			
				for(k = 0; k < A->r; k++) {
					fmpq_submul(fmpq_mat_entry(B, k, i),
								mu,
								fmpq_mat_entry(B, k, j));
				}
			}
		}
	}
	
	fmpq_clear(num);
	fmpq_clear(den);
	fmpq_clear(mu);
}

int main(void)
{
	int i;
    FLINT_TEST_INIT(state);
    

    flint_printf("gso....");
    fflush(stdout);

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpq_mat_t A;

        slong m, n, bits;

        m = n_randint(state, 10);
        n = n_randint(state, 10);

        bits = 1 + n_randint(state, 100);

        fmpq_mat_init(A, m, n);

        fmpq_mat_randtest(A, state, bits);
        
        fmpq_mat_gso(A, A);
        
        fmpq_t dot;
        fmpq_init(dot);
        int j, k, l;
        for(j = 0; j < n; j++) {
			for(k = j + 1; k < n; k++) {
				
				fmpq_zero(dot);
									  
				for(l = 0; l < m; l++)
				{
					fmpq_addmul(dot,
								fmpq_mat_entry(A, l, j),
								fmpq_mat_entry(A, l, k));
				}
				
				if (!fmpq_is_zero(dot))
				{
					flint_printf("FAIL:\n");
					flint_printf("A:\n");
					fmpq_mat_print(A);
					abort();
				}	
			}
		}
		
        fmpq_mat_clear(A);
        fmpq_clear(dot);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
