#pragma once

#include "../types.h"

/// Problem configuration.
typedef struct config_s {
    usz dim_x;
    usz dim_y;
    usz dim_z;
    usz niter;
} config_t;

/// Parse configuration from a file.
config_t config_parse_from_file(char const file_name[static 1]);

/// Retrieve size of x-axis from configuration.
usz config_dim_x(config_t self);

/// Retrieve size of y-axis from configuration.
usz config_dim_y(config_t self);

/// Retrieve size of z-axis from configuration.
usz config_dim_z(config_t self);

/// Retrieve number of iterations from configuration.
usz config_niter(config_t self);

/// Prints a configuration.
void config_print(config_t const* self);
