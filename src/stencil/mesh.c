#include "stencil/mesh.h"

#include "logging.h"

#include <assert.h>
#include <stdlib.h>

mesh_t mesh_new(usz dim_x, usz dim_y, usz dim_z, mesh_kind_t kind)
{
    usz const ghost_size = 2 * STENCIL_ORDER;

    cell_kind_t *kind_cell = malloc(sizeof(cell_kind_t) * (dim_x + ghost_size) * (dim_y + ghost_size) * (dim_z + ghost_size));
    f64 *value = malloc(sizeof(f64) * (dim_x + ghost_size) * (dim_y + ghost_size) * (dim_z + ghost_size));

    return (mesh_t){
        .dim_x = dim_x + ghost_size,
        .dim_y = dim_y + ghost_size,
        .dim_z = dim_z + ghost_size,
        .value = value,
        .kind_cell = kind_cell,
        .kind = kind,
    };
}

void mesh_drop(mesh_t *self)
{
    if (NULL != self->value)
    {
        free(self->value);
    }

    if (NULL != self->kind_cell)
    {
        free(self->kind_cell);
    }
}

static char const *mesh_kind_as_str(mesh_t const *self)
{
    static char const *MESH_KINDS_STR[] = {
        "CONSTANT",
        "INPUT",
        "OUTPUT",
    };
    return MESH_KINDS_STR[(usz)self->kind];
}

void mesh_print(mesh_t const *self, char const *name)
{
    fprintf(
        stderr,
        "****************************************\n"
        "MESH `%s`\n\tKIND: %s\n\tDIMS: %zux%zux%zu\n\tVALUES:\n",
        name,
        mesh_kind_as_str(self),
        self->dim_x,
        self->dim_y,
        self->dim_z);

	f64(*restrict span_value)[self->dim_y][self->dim_z] = (f64(*)[self->dim_y][self->dim_z])self->value;
	cell_kind_t(*restrict span_kind)[self->dim_y][self->dim_z] = (cell_kind_t(*)[self->dim_y][self->dim_z])self->kind_cell;

    for (usz i = 0; i < self->dim_x; ++i)
    {
        for (usz j = 0; j < self->dim_y; ++j)
        {
            for (usz k = 0; k < self->dim_z; ++k)
            {
                printf(
                    "%s%6.3lf%s ",
                    CELL_KIND_CORE == span_kind[i][j][k] ? "\x1b[1m" : "",
                    span_value[i][j][k],
                    "\x1b[0m");
            }
            puts("");
        }
        puts("");
    }
}

cell_kind_t mesh_set_cell_kind(mesh_t const *self, usz i, usz j, usz k)
{
    if ((i >= STENCIL_ORDER && i < self->dim_x - STENCIL_ORDER) &&
        (j >= STENCIL_ORDER && j < self->dim_y - STENCIL_ORDER) &&
        (k >= STENCIL_ORDER && k < self->dim_z - STENCIL_ORDER))
    {
        return CELL_KIND_CORE;
    }
    else
    {
        return CELL_KIND_PHANTOM;
    }
}

void mesh_copy_core(mesh_t *dst, mesh_t const *src)
{
    assert(dst->dim_x == src->dim_x);
    assert(dst->dim_y == src->dim_y);
    assert(dst->dim_z == src->dim_z);

    cell_kind_t(*restrict dst_kind)[dst->dim_y][dst->dim_z] = (cell_kind_t(*)[dst->dim_y][dst->dim_z])dst->kind_cell;
    cell_kind_t(*restrict src_kind)[dst->dim_y][dst->dim_z] = (cell_kind_t(*)[dst->dim_y][dst->dim_z])src->kind_cell;

	f64(*restrict dst_value)[dst->dim_y][dst->dim_z] = (f64(*)[dst->dim_y][dst->dim_z])dst->value;
	f64(*restrict src_value)[dst->dim_y][dst->dim_z] = (f64(*)[dst->dim_y][dst->dim_z])src->value;

    for (usz k = STENCIL_ORDER; k < dst->dim_z - STENCIL_ORDER; ++k)
    {
        for (usz j = STENCIL_ORDER; j < dst->dim_y - STENCIL_ORDER; ++j)
        {
            for (usz i = STENCIL_ORDER; i < dst->dim_x - STENCIL_ORDER; ++i)
            {
                assert(dst_kind[i][j][k] == CELL_KIND_CORE);
                assert(src_kind[i][j][k] == CELL_KIND_CORE);
                dst_value[i][j][k] = src_value[i][j][k];
            }
        }
    }
}
