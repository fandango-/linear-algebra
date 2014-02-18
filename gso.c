#include "flint/fmpq_mat.h"

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
	
	if(A->r == 0) {
		fmpq_mat_zero(B);
		return;
	}
		
	for(i = 0; i < A->c; i++) {
		for(j = 0; j < A->r; j++) {
			fmpq_set(fmpq_mat_entry(B, j, i),
					 fmpq_mat_entry(A, j, i));
		}

		for(j = 0; j < i; j++) {
			fmpq_mul(num,
					 fmpq_mat_entry(A, 0, i),
					 fmpq_mat_entry(B, 0, j));
					 
			for(k = 1; k < A->r; k++) {
				fmpq_addmul(num,
							fmpq_mat_entry(A, k, i),
							fmpq_mat_entry(B, k, j));
			}
			
			fmpq_mul(den,
					 fmpq_mat_entry(B, 0, j),
					 fmpq_mat_entry(B, 0, j));
										
			for(k = 1; k < A->r; k++) {
				fmpq_addmul(den,
							fmpq_mat_entry(B, k, j),
							fmpq_mat_entry(B, k, j));
			}
			
			fmpq_div(mu, num, den);
			
			for(k = 0; k < A->r; k++) {
				fmpq_submul(fmpq_mat_entry(B, k, i),
							mu,
							fmpq_mat_entry(B, k, j));
			}
		}
	}
	
	fmpq_clear(num);
	fmpq_clear(den);
	fmpq_clear(mu);
}
