#pragma once

#include "types.h"

#include <time.h>

/// Represents a span of time.
typedef struct duration_s {
    i64 secs;
    i32 nanos; // 0 <= nanos < NANOS_PER_SEC
} duration_t;

/// Represents a chronometer with two timepoints, a start and an end.
typedef struct chrono_s {
    struct timespec start;
    struct timespec stop;
} chrono_t;

/// Starts the chronometer.
static inline void chrono_start(chrono_t* self) {
    clock_gettime(CLOCK_MONOTONIC_RAW, &self->start);
}

/// Stops the chronometer.
static inline void chrono_stop(chrono_t* self) {
    clock_gettime(CLOCK_MONOTONIC_RAW, &self->stop);
}

/// Returns the duration between the two timepoints of a chronometer.
duration_t chrono_elapsed(chrono_t self);

/// Returns the duration in seconds as a 64-bit floating-point.
f64 duration_as_s_f64(duration_t self);

/// Returns the duration in milliseconds as a 64-bit floating-point.
f64 duration_as_ms_f64(duration_t self);

/// Returns the duration in microseconds as a 64-bit floating-point.
f64 duration_as_us_f64(duration_t self);

/// Returns the duration in nanoseconds as a 64-bit floating-point.
f64 duration_as_ns_f64(duration_t self);
