#include "stencil/solve.h"

#include <assert.h>
#include <immintrin.h>
#include <math.h>

#define min(a, b)               \
	({                          \
		__typeof__(a) _a = (a); \
		__typeof__(b) _b = (b); \
		_a < _b ? _a : _b;      \
	})

usz BI = 8;
usz BJ = 8;
usz BK = 4096;


void solve_jacobi(mesh_t *A, mesh_t const *B, mesh_t *C)
{
	assert(A->dim_x == B->dim_x && B->dim_x == C->dim_x);
	assert(A->dim_y == B->dim_y && B->dim_y == C->dim_y);
	assert(A->dim_z == B->dim_z && B->dim_z == C->dim_z);

	usz const dim_x = A->dim_x;
	usz const dim_y = A->dim_y;
	usz const dim_z = A->dim_z;

	usz min_k = 100000;
	usz min_j = 100000;
	usz min_i = 100000;

	f64(*restrict A_span_value)[A->dim_y][A->dim_z] = (f64(*)[A->dim_y][A->dim_z])A->value;
	f64(*restrict B_span_value)[A->dim_y][A->dim_z] = (f64(*)[A->dim_y][A->dim_z])B->value;
	f64(*restrict C_span_value)[A->dim_y][A->dim_z] = (f64(*)[A->dim_y][A->dim_z])C->value;

	f64 val;
	f64 pow17[STENCIL_ORDER];

	for (usz o = 0; o < STENCIL_ORDER; ++o)
		pow17[o] = 1 / pow(17.0, (f64)(o + 1));

	for (usz ii = STENCIL_ORDER; ii < dim_x - STENCIL_ORDER; ii += BI)
	{
		for (usz jj = STENCIL_ORDER; jj < dim_y - STENCIL_ORDER; jj += BJ)
		{

			for (usz kk = STENCIL_ORDER; kk < dim_z - STENCIL_ORDER; kk += BK)
			{
				min_i = min(ii + BI, dim_x - STENCIL_ORDER);
				min_j = min(jj + BJ, dim_y - STENCIL_ORDER);
				min_k = min(kk + BK, dim_z - STENCIL_ORDER);

				for (usz i = ii; i < min_i; ++i)
				{
					for (usz j = jj; j < min_j; ++j)
					{
						for (usz k = kk; k < min_k; ++k)
						{

							C_span_value[i][j][k] = A_span_value[i][j][k] * B_span_value[i][j][k];

							for (usz o = 1; o <= STENCIL_ORDER; ++o)
							{
								val = A_span_value[i + o][j][k] * B_span_value[i + o][j][k];
								val += A_span_value[i - o][j][k] * B_span_value[i - o][j][k];
								val += A_span_value[i][j + o][k] * B_span_value[i][j + o][k];
								val += A_span_value[i][j - o][k] * B_span_value[i][j - o][k];
								val += A_span_value[i][j][k + o] * B_span_value[i][j][k + o];
								val += A_span_value[i][j][k - o] * B_span_value[i][j][k - o];

								C_span_value[i][j][k] += val * pow17[o - 1];
							}
						}
					}
				}
			}
		}
	}

	mesh_copy_core(A, C);
}
