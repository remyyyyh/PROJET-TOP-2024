#include "stencil/solve.h"

#include <assert.h>
#include <math.h>

#define min(a, b)                                                                                  \
({                                                                                             \
	__typeof__(a) _a = (a);                                                                    \
	__typeof__(b) _b = (b);                                                                    \
	_a < _b ? _a : _b;                                                                         \
})

usz BI = 512;
usz BJ = 64;
usz BK = 64;

void solve_jacobi(mesh_t* A, mesh_t const* B, mesh_t* C) {
    assert(A->dim_x == B->dim_x && B->dim_x == C->dim_x);
    assert(A->dim_y == B->dim_y && B->dim_y == C->dim_y);
    assert(A->dim_z == B->dim_z && B->dim_z == C->dim_z);

    usz const dim_x = A->dim_x;
    usz const dim_y = A->dim_y;
    usz const dim_z = A->dim_z;

	usz min_k = 100000;
	usz min_j = 100000;
	usz min_i = 100000;

    for (usz kk = STENCIL_ORDER; kk < dim_z - STENCIL_ORDER; kk += BK) {
        for (usz jj = STENCIL_ORDER; jj < dim_y - STENCIL_ORDER; jj += BJ) {
            for (usz ii = STENCIL_ORDER; ii < dim_x - STENCIL_ORDER; ii += BI) {
				min_i = min(ii + BI, dim_x - STENCIL_ORDER);
				min_j = min(ii + BJ, dim_y - STENCIL_ORDER);
				min_k = min(ii + BK, dim_z - STENCIL_ORDER);


				for (usz k = kk; k < min_k; ++k) {
					for (usz j = jj; j < min_j; ++j) {
						for (usz i = ii; i < min_i; ++i) {
							C->cells[i][j][k].value = A->cells[i][j][k].value * B->cells[i][j][k].value;

							for (usz o = 1; o <= STENCIL_ORDER; ++o) {
								C->cells[i][j][k].value += A->cells[i + o][j][k].value *
														B->cells[i + o][j][k].value / pow(17.0, (f64)o);
								C->cells[i][j][k].value += A->cells[i - o][j][k].value *
														B->cells[i - o][j][k].value / pow(17.0, (f64)o);
								C->cells[i][j][k].value += A->cells[i][j + o][k].value *
														B->cells[i][j + o][k].value / pow(17.0, (f64)o);
								C->cells[i][j][k].value += A->cells[i][j - o][k].value *
														B->cells[i][j - o][k].value / pow(17.0, (f64)o);
								C->cells[i][j][k].value += A->cells[i][j][k + o].value *
														B->cells[i][j][k + o].value / pow(17.0, (f64)o);
								C->cells[i][j][k].value += A->cells[i][j][k - o].value *
														B->cells[i][j][k - o].value / pow(17.0, (f64)o);
							}
						}
					}
				}
			}
		}
	}

    mesh_copy_core(A, C);
}
