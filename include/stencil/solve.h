#pragma once

#include "mesh.h"

void solve_jacobi(mesh_t* A, mesh_t const* B, mesh_t* C);
void solve_jacobi_blocked(mesh_t* A, mesh_t const* B, mesh_t* C);
