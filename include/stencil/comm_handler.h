#pragma once

#include "stencil/mesh.h"
#include "types.h"

#include <mpi.h>

/// Enum for communication kind (either a send or a receive operation).
typedef enum comm_kind_e {
    COMM_KIND_SEND_OP,
    COMM_KIND_RECV_OP,
} comm_kind_t;
typedef int MPI_Syncfunc_t(MPI_Comm);

/// Handler for MPI communications between neighboor processes (ghost cell exchanges).
typedef struct comm_handler_s {
    /// Number of local meshes on the X axis.
    u32 nb_x;
    /// Number of local meshes on the Y axis.
    u32 nb_y;
    /// Number of local meshes on the Z axis.
    u32 nb_z;
    /// X coordinate of local mesh inside the global one.
    u32 coord_x;
    /// Y coordinate of local mesh inside the global one.
    u32 coord_y;
    /// Y coordinate of local mesh inside the global one.
    u32 coord_z;
    /// X dimension of the local mesh.
    usz loc_dim_x;
    /// Y dimension of the local mesh.
    usz loc_dim_y;
    /// Z dimension of the local mesh.
    usz loc_dim_z;
    /// Rank of the left neighboor process, -1 if none.
    i32 id_left;
    /// Rank of the right neighboor process, -1 if none.
    i32 id_right;
    /// Rank of the top neighboor process, -1 if none.
    i32 id_top;
    /// Rank of the bottom neighboor process, -1 if none.
    i32 id_bottom;
    /// Rank of the back neighboor process, -1 if none.
    i32 id_back;
    /// Rank of the front neighboor process, -1 if none.
    i32 id_front;
} comm_handler_t;

comm_handler_t comm_handler_new(u32 rank, u32 comm_size, usz dim_x, usz dim_y, usz dim_z);

void comm_handler_print(comm_handler_t const* self);

void comm_handler_ghost_exchange(comm_handler_t const* self, mesh_t* mesh);
