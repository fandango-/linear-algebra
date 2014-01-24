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

typedef struct frac_struct {
	int num;
	int den;
} Fraction;

int rows;
int cols;

Fraction frac_divide(Fraction a, Fraction b)
{
	return (Fraction) { a.num * b.den, a.den * b.num };
}

Fraction frac_subtract(Fraction a, Fraction b)
{
	return (Fraction) { a.num * b.den - b.num * a.den, a.den * b.den };
}

Fraction frac_multiply(Fraction a, Fraction b)
{
	return (Fraction) { a.num * b.num, a.den * b.den };
}

int gcd(int a, int b)
{
	return ( b == 0 ) ? a : gcd(b, a % b);
}

void inverse(Fraction m[rows][2*cols])
{
	int r = 0;
	int i, j, k, l;
	Fraction temp;
	int g;
	for(j = 0; j < 2 * cols; j++) {
		l = -1;
		i = r;
		while(l == -1 && i < rows) {
			if(m[i][j].num != 0) {
				l = i;
			}
			i++;
		}
		if(l != -1) {
			for(k = 0; k < 2 * cols; k++) {
				temp = m[r][k];
				m[r][k] = m[l][k];
				m[l][k] = temp;
			}
			temp = m[r][j];
			for(k = 0; k < 2 * cols; k++) {
				m[r][k] = frac_divide(m[r][k], temp);
				g = gcd(m[r][k].num, m[r][k].den);
				m[r][k].num /= g;
				m[r][k].den /= g;
			}
			for(i = 0; i < rows; i++) {
				temp = m[i][j];
				if(i != r) {
					for(k = 0; k < 2 * cols; k++) {
						m[i][k] = frac_subtract(m[i][k], frac_multiply(temp, m[r][k]));
						g = gcd(m[i][k].num, m[i][k].den);
						m[i][k].num /= g;
						m[i][k].den /= g;
					}
				}
			}
			r++;
		}
	}
	for(i = 0; i < rows; i++) {
		for(j = cols; j < 2 * cols; j++) {
			g = gcd(m[i][j].num, m[i][j].den);
			m[i][j].num /= g;
			m[i][j].den /= g;
			if(m[i][j].den == 1) {
				printf("%d\t", m[i][j].num);
			} else {
				printf("%d/%d\t", m[i][j].num, m[i][j].den);
			}
		}
		printf("\n");
	}
	return;
}

int rref(Fraction m[rows][cols])
{
	int r = 0;
	int i, j, k, l;
	Fraction temp;
	int g;
	for(j = 0; j < cols; j++) {
		l = -1;
		i = r;
		while(l == -1 && i < rows) {
			if(m[i][j].num != 0) {
				l = i;
			}
			i++;
		}
		if(l != -1) {
			for(k = 0; k < cols; k++) {
				temp = m[r][k];
				m[r][k] = m[l][k];
				m[l][k] = temp;
			}
			temp = m[r][j];
			for(k = 0; k < cols; k++) {
				m[r][k] = frac_divide(m[r][k], temp);
				g = gcd(m[r][k].num, m[r][k].den);
				m[r][k].num /= g;
				m[r][k].den /= g;
			}
			for(i = 0; i < rows; i++) {
				temp = m[i][j];
				if(i != r) {
					for(k = 0; k < cols; k++) {
						m[i][k] = frac_subtract(m[i][k], frac_multiply(temp, m[r][k]));
						g = gcd(m[i][k].num, m[i][k].den);
						m[i][k].num /= g;
						m[i][k].den /= g;
					}
				}
			}
			r++;
		}
	}
	int flag = 1;
	for(i = 0; i < rows; i++) {
		for(j = 0; j < cols; j++) {
			g = gcd(m[i][j].num, m[i][j].den);
			m[i][j].num /= g;
			m[i][j].den /= g;
			if(m[i][j].den == 1) {
				printf("%d\t", m[i][j].num);
			} else {
				printf("%d/%d\t", m[i][j].num, m[i][j].den);
			}
			if(rows == cols) {
				if(i == j) {
					if(m[i][j].num != 1) {
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
			scanf("%d", &m[i][j].num);
			m[i][j].den = 1;
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
						mi[i][j].num = 1;
					} else {
						mi[i][j].num = 0;
					}
					mi[i][j].den = 1;
				}
			}
			printf("It is non-singular and its inverse is:\n");
			inverse(mi);
		}
	}
	return 0;
}

