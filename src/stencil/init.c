#include "stencil/init.h"

#include "stencil/comm_handler.h"
#include "stencil/mesh.h"

#include <assert.h>
#include <math.h>
#include <string.h>

static inline f64 compute_core_pressure(f64 i, f64 j, f64 k) {
    return sin((f64)k * cos((f64)i + 0.311) * cos((f64)j + 0.817) + 0.613);
}

static void setup_mesh_cell_values(mesh_t* mesh, comm_handler_t const* comm_handler) {

    f64(*restrict span_value)[mesh->dim_y][mesh->dim_z] = (f64(*)[mesh->dim_y][mesh->dim_z])mesh->value;

    usz const dim_x = mesh->dim_x;
    usz const dim_y = mesh->dim_y;
    usz const dim_z = mesh->dim_z;

    switch (mesh->kind)    {

        case MESH_KIND_CONSTANT:
            for (usz i = 0; i < dim_x; ++i) 
                for (usz j = 0; j < dim_y; ++j) 
                    for (usz k = 0; k < dim_z; ++k) 
                        span_value[i][j][k] =  sin((f64)k * cos((f64)i + 0.311) * cos((f64)j + 0.817) + 0.613);
            break;
        

        case MESH_KIND_INPUT:
            for (usz i = STENCIL_ORDER; i < dim_x - STENCIL_ORDER; ++i) 
                for (usz j = STENCIL_ORDER; j < dim_y - STENCIL_ORDER; ++j) 
                    for (usz k = STENCIL_ORDER; k < dim_z - STENCIL_ORDER; ++k) 
                        span_value[i][j][k] = 1.0;
            
            for (usz i = 0; i < STENCIL_ORDER; ++i) 
                for (usz j = 0; j < STENCIL_ORDER; ++j) 
                    for (usz k = 0; k < STENCIL_ORDER; ++k) 
                        span_value[i][j][k] = 0.0;

            for (usz i = dim_x - STENCIL_ORDER; i < dim_x; ++i) 
                for (usz j = dim_y - STENCIL_ORDER; j < dim_y; ++j) 
                    for (usz k = dim_z - STENCIL_ORDER; k < dim_z; ++k) 
                        span_value[i][j][k] = 0.0;

            break;

        case MESH_KIND_OUTPUT:
            memset(span_value,0.0,dim_x*dim_y*dim_z*sizeof(f64));

            break;

        default:
            __builtin_unreachable();
    }

}

static void setup_mesh_cell_kinds(mesh_t* mesh) {

    cell_kind_t(*restrict span_kind)[mesh->dim_y][mesh->dim_z] = (cell_kind_t(*)[mesh->dim_y][mesh->dim_z])mesh->value;


    for (usz i = STENCIL_ORDER; i < mesh->dim_x - STENCIL_ORDER; ++i) 
        for (usz j = STENCIL_ORDER; j < mesh->dim_y - STENCIL_ORDER; ++j) 
            for (usz k = STENCIL_ORDER; k < mesh->dim_z - STENCIL_ORDER; ++k) 
                span_kind[i][j][k] = CELL_KIND_CORE;
    

    for (usz i = 0; i < STENCIL_ORDER; ++i) 
        for (usz j = 0; j < STENCIL_ORDER; ++j)
            for (usz k = 0; k < STENCIL_ORDER; ++k) 
                span_kind[i][j][k] = CELL_KIND_PHANTOM;


    for (usz i = mesh->dim_x - STENCIL_ORDER; i < mesh->dim_x; ++i) 
        for (usz j = mesh->dim_y - STENCIL_ORDER; j < mesh->dim_y; ++j) 
            for (usz k = mesh->dim_z - STENCIL_ORDER; k < mesh->dim_z; ++k) 
                span_kind[i][j][k] = CELL_KIND_PHANTOM;

}

void init_meshes(mesh_t* A, mesh_t* B, mesh_t* C, comm_handler_t const* comm_handler) {
    assert(
        A->dim_x == B->dim_x && B->dim_x == C->dim_x &&
        C->dim_x == comm_handler->loc_dim_x + STENCIL_ORDER * 2
    );
    assert(
        A->dim_y == B->dim_y && B->dim_y == C->dim_y &&
        C->dim_y == comm_handler->loc_dim_y + STENCIL_ORDER * 2
    );
    assert(
        A->dim_z == B->dim_z && B->dim_z == C->dim_z &&
        C->dim_z == comm_handler->loc_dim_z + STENCIL_ORDER * 2
    );

    setup_mesh_cell_kinds(A);
    setup_mesh_cell_kinds(B);
    setup_mesh_cell_kinds(C);

    setup_mesh_cell_values(A, comm_handler);
    setup_mesh_cell_values(B, comm_handler);
    setup_mesh_cell_values(C, comm_handler);
}
