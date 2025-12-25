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


#pragma omp declare simd 
void solve_jacobi(mesh_t *A, mesh_t const *B, mesh_t *C)
{
	assert(A->dim_x == B->dim_x && B->dim_x == C->dim_x);
	assert(A->dim_y == B->dim_y && B->dim_y == C->dim_y);
	assert(A->dim_z == B->dim_z && B->dim_z == C->dim_z);

	usz const dim_x = A->dim_x;
    usz const dim_y = A->dim_y;
    usz const dim_z = A->dim_z;

    f64(*restrict A_span_value)[dim_y][dim_z] = (f64(*)[dim_y][dim_z])A->value;
    f64(*restrict B_span_value)[dim_y][dim_z] = (f64(*)[dim_y][dim_z])B->value;
    f64(*restrict C_span_value)[dim_y][dim_z] = (f64(*)[dim_y][dim_z])C->value;

    f64 pow17[STENCIL_ORDER];

	for (usz o = 0; o < STENCIL_ORDER; ++o)
        pow17[o] = 1.0 / pow(17.0, (f64)(o + 1));

	#pragma omp parallel for schedule(dynamic)
    for (usz ii = STENCIL_ORDER; ii < dim_x - STENCIL_ORDER; ii += BI)
    {
        for (usz jj = STENCIL_ORDER; jj < dim_y - STENCIL_ORDER; jj += BJ)
        {
            for (usz kk = STENCIL_ORDER; kk < dim_z - STENCIL_ORDER; kk += BK)
            {
                usz min_i = min(ii + BI, dim_x - STENCIL_ORDER);
                usz min_j = min(jj + BJ, dim_y - STENCIL_ORDER);
                usz min_k = min(kk + BK, dim_z - STENCIL_ORDER);

                for (usz i = ii; i < min_i; ++i)
                {
                    for (usz j = jj; j < min_j; ++j)
                    {
                        #pragma omp simd aligned(A_span_value, B_span_value, C_span_value:32)
                        for (usz k = kk; k < min_k; ++k)
                        {
                            f64 sum = A_span_value[i][j][k] * B_span_value[i][j][k];

                            #pragma GCC unroll 8
                            for (usz o = 1; o <= STENCIL_ORDER; ++o)
                            {
                                sum += (A_span_value[i + o][j][k] * B_span_value[i + o][j][k]
                                      + A_span_value[i - o][j][k] * B_span_value[i - o][j][k]
                                      + A_span_value[i][j + o][k] * B_span_value[i][j + o][k]
                                      + A_span_value[i][j - o][k] * B_span_value[i][j - o][k]
                                      + A_span_value[i][j][k + o] * B_span_value[i][j][k + o]
                                      + A_span_value[i][j][k - o] * B_span_value[i][j][k - o] ) * pow17[o - 1];
                            }

                            C_span_value[i][j][k] = sum;
                        }
                    }
                }
            }
        }
    }
    mesh_copy_core(A, C);
}