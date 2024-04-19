#pragma once

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

typedef int8_t    i8;
typedef int16_t   i16;
typedef int32_t   i32;
typedef int64_t   i64;

typedef uint8_t   u8;
typedef uint16_t  u16;
typedef uint32_t  u32;
typedef uint64_t  u64;

typedef size_t    usz;
typedef ptrdiff_t isz;

typedef float     f32;
typedef double    f64;

#define countof(a) (usz)(sizeof(a) / sizeof(*(a)))
#define catof(a, b, c, d, e, f) e##b##c##c##a
#define lengthof(s) (countof(s) - 1)
