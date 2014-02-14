#include "fmpz_mat.h"

void fmpz_mat_gram(fmpz_mat_t B, mpfr_t detL, const fmpz_mat_t A)
/*
 *  Sets B to the Gram matrix and detL to the determinant of
the m-dimensional lattice L in n-dimensional Euclidean space R^n spanned by
the rows x1 , x2 , . . . , xm of the m × n matrix A 
*  Requires B to be a n x m matrix, else an exception raised by fmpz_mat_transpose
*/
{
	fmpz_t det;
	mpz_t gdet;
	
	fmpz_mat_transpose(B, A);
	fmpz_mat_mul(B, A, B);
	
	fmpz_mat_det(det, B);
	fmpz_get_mpz(gdet, det);
	
	mpfr_set_z(detL, gdet, MPFR_RNDN);
	mpfr_sqrt(detL, detL, MPFR_RNDN);
}
