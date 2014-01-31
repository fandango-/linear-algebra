/*
 * rref.c
 * 
 * Copyright 2014 aman <aman@aman-Aspire-4750>
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA 02110-1301, USA.
 * 
 * 
 */


#include <stdio.h>
#include "flint/fmpz.h"

typedef struct frac_struct {
	fmpz_t num;
	fmpz_t den;
} Fraction;

int rows;
int cols;

Fraction frac_divide(Fraction res, Fraction a, Fraction b)
{
	fmpz_mul(res.num, a.num, b.den);
	fmpz_mul(res.den, a.den, b.num);
	return res;
}

Fraction frac_subtract(Fraction res, Fraction a, Fraction b)
{
	fmpz_mul(res.num, a.num, b.den);
	fmpz_submul(res.num, b.num, a.den);
	fmpz_mul(res.den, a.den, b.den);
	return res;
}

Fraction frac_multiply(Fraction res, Fraction a, Fraction b)
{
	fmpz_mul(res.num, a.num, b.num);
	fmpz_mul(res.den, a.den, b.den);
	return res;
}

void inverse(Fraction m[rows][2*cols])
{
	int r = 0;
	int i, j, k, l;
	Fraction temp, det;
	fmpz_init_set_ui(det.num, 1);
	fmpz_init_set_ui(det.den, 1);
	fmpz_t g;
	for(j = 0; j < 2 * cols; j++) {
		l = -1;
		i = r;
		while(l == -1 && i < rows) {
			if(!fmpz_is_zero(m[i][j].num)) {
				l = i;
			}
			i++;
		}
		if(l != -1) {
			if(l != r) {
				for(k = 0; k < 2 * cols; k++) {
					temp = m[r][k];
					m[r][k] = m[l][k];
					m[l][k] = temp;
				}
				fmpz_mul_si(det.num, det.num, -1);
			}
			temp = m[r][j];
			det = frac_multiply(det, det, temp);
			fmpz_gcd(g, det.num, det.den);
			fmpz_divexact(det.num, det.num, g);
			fmpz_divexact(det.den, det.den, g);
			for(k = 0; k < 2 * cols; k++) {
				m[r][k] = frac_divide(m[r][k], m[r][k], temp);
				fmpz_gcd(g, m[r][k].num, m[r][k].den);
				fmpz_divexact(m[r][k].num, m[r][k].num, g);
				fmpz_divexact(m[r][k].den, m[r][k].den, g);
			}
			for(i = 0; i < rows; i++) {
				temp = m[i][j];
				if(i != r) {
					for(k = 0; k < 2 * cols; k++) {
						m[i][k] = frac_subtract(m[i][k], m[i][k], frac_multiply(temp, temp, m[r][k]));
						fmpz_gcd(g, m[i][k].num, m[i][k].den);
						fmpz_divexact(m[i][k].num, m[i][k].num, g);
						fmpz_divexact(m[i][k].den, m[i][k].den, g);
					}
				}
			}
			r++;
		}
	}
	for(i = 0; i < rows; i++) {
		for(j = cols; j < 2 * cols; j++) {
			if(fmpz_equal_si(m[i][j].den, -1)) {
				fmpz_neg(m[i][j].num, m[i][j].num);
				fmpz_print(m[i][j].num);
				printf("\t");
			} else if(fmpz_equal_ui(m[i][j].den, 1)) {
				fmpz_print(m[i][j].num);
				printf("\t");
			} else if(fmpz_cmp_si(m[i][j].den, 0) < 0) {
				fmpz_neg(m[i][j].num, m[i][j].num);
				fmpz_neg(m[i][j].den, m[i][j].den);
				fmpz_print(m[i][j].num);
				printf("/");
				fmpz_print(m[i][j].den);
				printf("\t");
			} else {
				fmpz_print(m[i][j].num);
				printf("/");
				fmpz_print(m[i][j].den);
				printf("\t");
			}
		}
		printf("\n");
	}
	if(fmpz_equal_si(det.den, -1)) {
		fmpz_neg(det.num, det.num);
	}
	printf("Its determinant is:\n");
	fmpz_print(det.num);
	printf("\n");
	return;
}

int rref(Fraction m[rows][cols])
{
	int r = 0;
	int i, j, k, l;
	Fraction temp;
	fmpz_t g;
	for(j = 0; j < cols; j++) {
		l = -1;
		i = r;
		while(l == -1 && i < rows) {
			if(!fmpz_is_zero(m[i][j].num)) {
				l = i;
			}
			i++;
		}
		if(l != -1) {
			if(l != r) {
				for(k = 0; k < cols; k++) {
					temp = m[r][k];
					m[r][k] = m[l][k];
					m[l][k] = temp;
				}
			}
			temp = m[r][j];
			for(k = 0; k < cols; k++) {
				m[r][k] = frac_divide(m[r][k], m[r][k], temp);
				fmpz_gcd(g, m[r][k].num, m[r][k].den);
				fmpz_divexact(m[r][k].num, m[r][k].num, g);
				fmpz_divexact(m[r][k].den, m[r][k].den, g);
			}
			for(i = 0; i < rows; i++) {
				temp = m[i][j];
				if(i != r) {
					for(k = 0; k < cols; k++) {
						m[i][k] = frac_subtract(m[i][k], m[i][k], frac_multiply(temp, temp, m[r][k]));
						fmpz_gcd(g, m[i][k].num, m[i][k].den);
						fmpz_divexact(m[i][k].num, m[i][k].num, g);
						fmpz_divexact(m[i][k].den, m[i][k].den, g);
					}
				}
			}
			r++;
		}
	}
	int flag = 1;
	for(i = 0; i < rows; i++) {
		for(j = 0; j < cols; j++) {
			if(fmpz_equal_si(m[i][j].den, -1)) {
				fmpz_neg(m[i][j].num, m[i][j].num);
				fmpz_print(m[i][j].num);
				printf("\t");
			} else if(fmpz_equal_ui(m[i][j].den, 1)) {
				fmpz_print(m[i][j].num);
				printf("\t");
			} else if(fmpz_cmp_si(m[i][j].den, 0) < 0) {
				fmpz_neg(m[i][j].num, m[i][j].num);
				fmpz_neg(m[i][j].den, m[i][j].den);
				fmpz_print(m[i][j].num);
				printf("/");
				fmpz_print(m[i][j].den);
				printf("\t");
			} else {
				fmpz_print(m[i][j].num);
				printf("/");
				fmpz_print(m[i][j].den);
				printf("\t");
			}
			if(rows == cols) {
				if(i == j) {
					if(!fmpz_equal_ui(m[i][j].num, 1)) {
						flag = 0;
					}
				}
			}
		}
		printf("\n");
	}
	return flag;
}

int main(int argc, char **argv)
{
	printf("Enter number of rows:\n");
	scanf("%d", &rows);
	printf("Enter number of columns:\n");
	scanf("%d", &cols);
	Fraction m[rows][cols];
	Fraction mi[rows][2*cols];
    int i, j;
    printf("Enter the elements:\n");
    for(i = 0; i < rows; i++) {
		for(j = 0; j < cols; j++) {
			fmpz_read(m[i][j].num);
			fmpz_init_set_ui(m[i][j].den, 1);
			if(rows == cols) {
				mi[i][j] = m[i][j];
			}
		}
	}
	printf("The reduced row echelon form of the given matrix is:\n");
    int f = rref(m);
    if(rows == cols) {
		printf("The given matrix is square.\n");
		if(f == 0) {
			printf("However, it is singular and hence not invertible.\n");
		} else {
			for(i = 0; i < rows; i++) {
				for(j = cols; j < 2 * cols; j++) {
					if(j == i + cols) {
						fmpz_init_set_ui(mi[i][j].num, 1);
					} else {
						fmpz_init(mi[i][j].num);
					}
					fmpz_init_set_ui(mi[i][j].den, 1);
				}
			}
			printf("It is non-singular and its inverse is:\n");
			inverse(mi);
		}
	}
	return 0;
}

