#pragma once

#include "mesh.h"
#include "comm_handler.h"

void init_meshes(mesh_t* A, mesh_t* B, mesh_t* C, comm_handler_t const* comm_handler);
