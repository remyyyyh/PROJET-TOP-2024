#define _GNU_SOURCE

#include "stencil/config.h"

#include "logging.h"

#include <stdio.h>
#include <string.h>

static inline config_t config_default() {
    return (config_t){
        .dim_x = 100,
        .dim_y = 100,
        .dim_z = 100,
        .niter = 5,
    };
}

config_t config_parse_from_file(char const file_name[static 1]) {
    FILE* cfp = fopen(file_name, "rb");
    if (NULL == cfp) {
        warn("failed to open configuration file %s, using default", file_name);
        return config_default();
    }

    config_t self = config_default();
    usz MAX_LINE_LEN = 64;
    char* line_buf = malloc(MAX_LINE_LEN);
    usz line_num = 0;
    while (0 == feof(cfp)) {
        line_num += 1;
        getline(&line_buf, &MAX_LINE_LEN, cfp);
        if (NULL == line_buf) {
            warn("failed to read line %zu in file %s, using default", line_num, file_name);
            return config_default();
        }

        if ('#' == line_buf[0]) {
            continue;
        }

        usz const MAX_TOKEN_LEN = 5;
        char key[MAX_TOKEN_LEN + 5];
        usz val;
        sscanf(line_buf, "%5s=%zu\n", key, &val);

        if (strcmp("dim_x", key) == 0) {
            self.dim_x = val;
        } else if (strcmp("dim_y", key) == 0) {
            self.dim_y = val;
        } else if (strcmp("dim_z", key) == 0) {
            self.dim_z = val;
        } else if (strcmp("niter", key) == 0) {
            self.niter = val;
        } else {
            warn("unknown key `%s` at line %zu", key, line_num);
            return config_default();
        }
    }

    free(line_buf);
    fclose(cfp);
    return self;
}

inline usz config_dim_x(config_t self) {
    return self.dim_x;
}

inline usz config_dim_y(config_t self) {
    return self.dim_y;
}

inline usz config_dim_z(config_t self) {
    return self.dim_z;
}

inline usz config_niter(config_t self) {
    return self.niter;
}

void config_print(config_t const* self) {
    fprintf(
        stderr,
        "****************************************\n"
        "         STENCIL CONFIGURATION\n"
        "X-axis dimension ................... %zu\n"
        "Y-axis dimension ................... %zu\n"
        "Z-axis dimension ................... %zu\n"
        "Number of iterations ............... %zu\n",
        self->dim_x,
        self->dim_y,
        self->dim_z,
        self->niter
    );
}
