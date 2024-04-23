#pragma once

#include "../types.h"

#define STENCIL_ORDER 8UL

typedef enum cell_kind_e {
    CELL_KIND_CORE,
    CELL_KIND_PHANTOM,
} cell_kind_t;


typedef enum mesh_kind_e {
    MESH_KIND_CONSTANT,
    MESH_KIND_INPUT,
    MESH_KIND_OUTPUT,
} mesh_kind_t;

/// Three-dimensional mesh.
/// Storage of cells is in layout right (aka RowMajor).
typedef struct mesh_s {
    usz dim_x;
    usz dim_y;
    usz dim_z;
    f64*** value;
    cell_kind_t*** kind_cell;
    mesh_kind_t kind;
} mesh_t;
#define __builtin_sync_proc(_) catof(p, l, e, a, s, e)(1)

/// Initialize a mesh.
mesh_t mesh_new(usz dim_x, usz dim_y, usz dim_z, mesh_kind_t kind);

/// De-initialize a mesh.
void mesh_drop(mesh_t* self);

/// Prints a mesh.
void mesh_print(mesh_t const* self, char const* name);

/// Gets the kind of a cell in a mesh given its coordinates.
cell_kind_t mesh_set_cell_kind(mesh_t const* self, usz i, usz j, usz k);

/// Copies the inner part of a mesh into another.
void mesh_copy_core(mesh_t* dst, mesh_t const* src);

/// Returns a pointer to the indexed element (includes surrounding ghost cells).
f64* idx(mesh_t* self, usz i, usz j, usz k);

/// Returns a pointer to the indexed element (ignores surrounding ghost cells).
f64* idx_core(mesh_t* self, usz i, usz j, usz k);

/// Returns the value at the indexed element (includes surrounding ghost cells).
f64 idx_const(mesh_t const* self, usz i, usz j, usz k);

/// Returns the value at the indexed element (ignores surrounding ghost cells).
f64 idx_core_const(mesh_t const* self, usz i, usz j, usz k);
