# 1 "CMakeCUDACompilerId.cu"
#pragma GCC diagnostic ignored "-Wunused-local-typedefs"
# 1
#pragma GCC diagnostic push
# 1
#pragma GCC diagnostic ignored "-Wunused-variable"
# 1
#pragma GCC diagnostic ignored "-Wunused-function"
# 1
static char __nv_inited_managed_rt = 0; static void **__nv_fatbinhandle_for_managed_rt; static void __nv_save_fatbinhandle_for_managed_rt(void **in){__nv_fatbinhandle_for_managed_rt = in;} static char __nv_init_managed_rt_with_module(void **); static inline void __nv_init_managed_rt(void) { __nv_inited_managed_rt = (__nv_inited_managed_rt ? __nv_inited_managed_rt                 : __nv_init_managed_rt_with_module(__nv_fatbinhandle_for_managed_rt));}
# 1
#pragma GCC diagnostic pop
# 1
#pragma GCC diagnostic ignored "-Wunused-variable"

# 1
#define __nv_is_extended_device_lambda_closure_type(X) false
#define __nv_is_extended_host_device_lambda_closure_type(X) false
#if defined(__nv_is_extended_device_lambda_closure_type) && defined(__nv_is_extended_host_device_lambda_closure_type)
#endif

# 1
# 61 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime.h"
#pragma GCC diagnostic push
# 64
#pragma GCC diagnostic ignored "-Wunused-function"
# 68 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/device_types.h"
#if 0
# 68
enum cudaRoundMode { 
# 70
cudaRoundNearest, 
# 71
cudaRoundZero, 
# 72
cudaRoundPosInf, 
# 73
cudaRoundMinInf
# 74
}; 
#endif
# 100 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/vector_types.h"
#if 0
# 100
struct char1 { 
# 102
signed char x; 
# 103
}; 
#endif
# 105 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/vector_types.h"
#if 0
# 105
struct uchar1 { 
# 107
unsigned char x; 
# 108
}; 
#endif
# 111 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/vector_types.h"
#if 0
# 111
struct __attribute((aligned(2))) char2 { 
# 113
signed char x, y; 
# 114
}; 
#endif
# 116 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/vector_types.h"
#if 0
# 116
struct __attribute((aligned(2))) uchar2 { 
# 118
unsigned char x, y; 
# 119
}; 
#endif
# 121 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/vector_types.h"
#if 0
# 121
struct char3 { 
# 123
signed char x, y, z; 
# 124
}; 
#endif
# 126 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/vector_types.h"
#if 0
# 126
struct uchar3 { 
# 128
unsigned char x, y, z; 
# 129
}; 
#endif
# 131 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/vector_types.h"
#if 0
# 131
struct __attribute((aligned(4))) char4 { 
# 133
signed char x, y, z, w; 
# 134
}; 
#endif
# 136 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/vector_types.h"
#if 0
# 136
struct __attribute((aligned(4))) uchar4 { 
# 138
unsigned char x, y, z, w; 
# 139
}; 
#endif
# 141 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/vector_types.h"
#if 0
# 141
struct short1 { 
# 143
short x; 
# 144
}; 
#endif
# 146 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/vector_types.h"
#if 0
# 146
struct ushort1 { 
# 148
unsigned short x; 
# 149
}; 
#endif
# 151 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/vector_types.h"
#if 0
# 151
struct __attribute((aligned(4))) short2 { 
# 153
short x, y; 
# 154
}; 
#endif
# 156 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/vector_types.h"
#if 0
# 156
struct __attribute((aligned(4))) ushort2 { 
# 158
unsigned short x, y; 
# 159
}; 
#endif
# 161 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/vector_types.h"
#if 0
# 161
struct short3 { 
# 163
short x, y, z; 
# 164
}; 
#endif
# 166 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/vector_types.h"
#if 0
# 166
struct ushort3 { 
# 168
unsigned short x, y, z; 
# 169
}; 
#endif
# 171 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/vector_types.h"
#if 0
# 171
struct __attribute((aligned(8))) short4 { short x; short y; short z; short w; }; 
#endif
# 172 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/vector_types.h"
#if 0
# 172
struct __attribute((aligned(8))) ushort4 { unsigned short x; unsigned short y; unsigned short z; unsigned short w; }; 
#endif
# 174 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/vector_types.h"
#if 0
# 174
struct int1 { 
# 176
int x; 
# 177
}; 
#endif
# 179 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/vector_types.h"
#if 0
# 179
struct uint1 { 
# 181
unsigned x; 
# 182
}; 
#endif
# 184 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/vector_types.h"
#if 0
# 184
struct __attribute((aligned(8))) int2 { int x; int y; }; 
#endif
# 185 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/vector_types.h"
#if 0
# 185
struct __attribute((aligned(8))) uint2 { unsigned x; unsigned y; }; 
#endif
# 187 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/vector_types.h"
#if 0
# 187
struct int3 { 
# 189
int x, y, z; 
# 190
}; 
#endif
# 192 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/vector_types.h"
#if 0
# 192
struct uint3 { 
# 194
unsigned x, y, z; 
# 195
}; 
#endif
# 197 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/vector_types.h"
#if 0
# 197
struct __attribute((aligned(16))) int4 { 
# 199
int x, y, z, w; 
# 200
}; 
#endif
# 202 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/vector_types.h"
#if 0
# 202
struct __attribute((aligned(16))) uint4 { 
# 204
unsigned x, y, z, w; 
# 205
}; 
#endif
# 207 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/vector_types.h"
#if 0
# 207
struct long1 { 
# 209
long x; 
# 210
}; 
#endif
# 212 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/vector_types.h"
#if 0
# 212
struct ulong1 { 
# 214
unsigned long x; 
# 215
}; 
#endif
# 222 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/vector_types.h"
#if 0
# 222
struct __attribute((aligned((2) * sizeof(long)))) long2 { 
# 224
long x, y; 
# 225
}; 
#endif
# 227 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/vector_types.h"
#if 0
# 227
struct __attribute((aligned((2) * sizeof(unsigned long)))) ulong2 { 
# 229
unsigned long x, y; 
# 230
}; 
#endif
# 234 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/vector_types.h"
#if 0
# 234
struct long3 { 
# 236
long x, y, z; 
# 237
}; 
#endif
# 239 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/vector_types.h"
#if 0
# 239
struct ulong3 { 
# 241
unsigned long x, y, z; 
# 242
}; 
#endif
# 244 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/vector_types.h"
#if 0
# 244
struct __attribute((aligned(16))) long4 { 
# 246
long x, y, z, w; 
# 247
}; 
#endif
# 249 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/vector_types.h"
#if 0
# 249
struct __attribute((aligned(16))) ulong4 { 
# 251
unsigned long x, y, z, w; 
# 252
}; 
#endif
# 254 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/vector_types.h"
#if 0
# 254
struct float1 { 
# 256
float x; 
# 257
}; 
#endif
# 276 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/vector_types.h"
#if 0
# 276
struct __attribute((aligned(8))) float2 { float x; float y; }; 
#endif
# 281 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/vector_types.h"
#if 0
# 281
struct float3 { 
# 283
float x, y, z; 
# 284
}; 
#endif
# 286 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/vector_types.h"
#if 0
# 286
struct __attribute((aligned(16))) float4 { 
# 288
float x, y, z, w; 
# 289
}; 
#endif
# 291 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/vector_types.h"
#if 0
# 291
struct longlong1 { 
# 293
long long x; 
# 294
}; 
#endif
# 296 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/vector_types.h"
#if 0
# 296
struct ulonglong1 { 
# 298
unsigned long long x; 
# 299
}; 
#endif
# 301 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/vector_types.h"
#if 0
# 301
struct __attribute((aligned(16))) longlong2 { 
# 303
long long x, y; 
# 304
}; 
#endif
# 306 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/vector_types.h"
#if 0
# 306
struct __attribute((aligned(16))) ulonglong2 { 
# 308
unsigned long long x, y; 
# 309
}; 
#endif
# 311 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/vector_types.h"
#if 0
# 311
struct longlong3 { 
# 313
long long x, y, z; 
# 314
}; 
#endif
# 316 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/vector_types.h"
#if 0
# 316
struct ulonglong3 { 
# 318
unsigned long long x, y, z; 
# 319
}; 
#endif
# 321 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/vector_types.h"
#if 0
# 321
struct __attribute((aligned(16))) longlong4 { 
# 323
long long x, y, z, w; 
# 324
}; 
#endif
# 326 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/vector_types.h"
#if 0
# 326
struct __attribute((aligned(16))) ulonglong4 { 
# 328
unsigned long long x, y, z, w; 
# 329
}; 
#endif
# 331 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/vector_types.h"
#if 0
# 331
struct double1 { 
# 333
double x; 
# 334
}; 
#endif
# 336 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/vector_types.h"
#if 0
# 336
struct __attribute((aligned(16))) double2 { 
# 338
double x, y; 
# 339
}; 
#endif
# 341 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/vector_types.h"
#if 0
# 341
struct double3 { 
# 343
double x, y, z; 
# 344
}; 
#endif
# 346 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/vector_types.h"
#if 0
# 346
struct __attribute((aligned(16))) double4 { 
# 348
double x, y, z, w; 
# 349
}; 
#endif
# 363 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/vector_types.h"
#if 0
typedef char1 
# 363
char1; 
#endif
# 364 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/vector_types.h"
#if 0
typedef uchar1 
# 364
uchar1; 
#endif
# 365 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/vector_types.h"
#if 0
typedef char2 
# 365
char2; 
#endif
# 366 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/vector_types.h"
#if 0
typedef uchar2 
# 366
uchar2; 
#endif
# 367 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/vector_types.h"
#if 0
typedef char3 
# 367
char3; 
#endif
# 368 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/vector_types.h"
#if 0
typedef uchar3 
# 368
uchar3; 
#endif
# 369 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/vector_types.h"
#if 0
typedef char4 
# 369
char4; 
#endif
# 370 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/vector_types.h"
#if 0
typedef uchar4 
# 370
uchar4; 
#endif
# 371 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/vector_types.h"
#if 0
typedef short1 
# 371
short1; 
#endif
# 372 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/vector_types.h"
#if 0
typedef ushort1 
# 372
ushort1; 
#endif
# 373 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/vector_types.h"
#if 0
typedef short2 
# 373
short2; 
#endif
# 374 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/vector_types.h"
#if 0
typedef ushort2 
# 374
ushort2; 
#endif
# 375 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/vector_types.h"
#if 0
typedef short3 
# 375
short3; 
#endif
# 376 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/vector_types.h"
#if 0
typedef ushort3 
# 376
ushort3; 
#endif
# 377 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/vector_types.h"
#if 0
typedef short4 
# 377
short4; 
#endif
# 378 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/vector_types.h"
#if 0
typedef ushort4 
# 378
ushort4; 
#endif
# 379 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/vector_types.h"
#if 0
typedef int1 
# 379
int1; 
#endif
# 380 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/vector_types.h"
#if 0
typedef uint1 
# 380
uint1; 
#endif
# 381 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/vector_types.h"
#if 0
typedef int2 
# 381
int2; 
#endif
# 382 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/vector_types.h"
#if 0
typedef uint2 
# 382
uint2; 
#endif
# 383 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/vector_types.h"
#if 0
typedef int3 
# 383
int3; 
#endif
# 384 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/vector_types.h"
#if 0
typedef uint3 
# 384
uint3; 
#endif
# 385 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/vector_types.h"
#if 0
typedef int4 
# 385
int4; 
#endif
# 386 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/vector_types.h"
#if 0
typedef uint4 
# 386
uint4; 
#endif
# 387 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/vector_types.h"
#if 0
typedef long1 
# 387
long1; 
#endif
# 388 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/vector_types.h"
#if 0
typedef ulong1 
# 388
ulong1; 
#endif
# 389 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/vector_types.h"
#if 0
typedef long2 
# 389
long2; 
#endif
# 390 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/vector_types.h"
#if 0
typedef ulong2 
# 390
ulong2; 
#endif
# 391 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/vector_types.h"
#if 0
typedef long3 
# 391
long3; 
#endif
# 392 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/vector_types.h"
#if 0
typedef ulong3 
# 392
ulong3; 
#endif
# 393 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/vector_types.h"
#if 0
typedef long4 
# 393
long4; 
#endif
# 394 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/vector_types.h"
#if 0
typedef ulong4 
# 394
ulong4; 
#endif
# 395 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/vector_types.h"
#if 0
typedef float1 
# 395
float1; 
#endif
# 396 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/vector_types.h"
#if 0
typedef float2 
# 396
float2; 
#endif
# 397 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/vector_types.h"
#if 0
typedef float3 
# 397
float3; 
#endif
# 398 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/vector_types.h"
#if 0
typedef float4 
# 398
float4; 
#endif
# 399 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/vector_types.h"
#if 0
typedef longlong1 
# 399
longlong1; 
#endif
# 400 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/vector_types.h"
#if 0
typedef ulonglong1 
# 400
ulonglong1; 
#endif
# 401 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/vector_types.h"
#if 0
typedef longlong2 
# 401
longlong2; 
#endif
# 402 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/vector_types.h"
#if 0
typedef ulonglong2 
# 402
ulonglong2; 
#endif
# 403 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/vector_types.h"
#if 0
typedef longlong3 
# 403
longlong3; 
#endif
# 404 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/vector_types.h"
#if 0
typedef ulonglong3 
# 404
ulonglong3; 
#endif
# 405 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/vector_types.h"
#if 0
typedef longlong4 
# 405
longlong4; 
#endif
# 406 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/vector_types.h"
#if 0
typedef ulonglong4 
# 406
ulonglong4; 
#endif
# 407 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/vector_types.h"
#if 0
typedef double1 
# 407
double1; 
#endif
# 408 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/vector_types.h"
#if 0
typedef double2 
# 408
double2; 
#endif
# 409 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/vector_types.h"
#if 0
typedef double3 
# 409
double3; 
#endif
# 410 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/vector_types.h"
#if 0
typedef double4 
# 410
double4; 
#endif
# 418 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/vector_types.h"
#if 0
# 418
struct dim3 { 
# 420
unsigned x, y, z; 
# 432
}; 
#endif
# 434 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/vector_types.h"
#if 0
typedef dim3 
# 434
dim3; 
#endif
# 143 "/usr/lib/gcc/x86_64-linux-gnu/10/include/stddef.h" 3
typedef long ptrdiff_t; 
# 209 "/usr/lib/gcc/x86_64-linux-gnu/10/include/stddef.h" 3
typedef unsigned long size_t; 
#if !defined(__CUDA_INCLUDE_COMPILER_INTERNAL_HEADERS__)
#define __CUDA_INCLUDE_COMPILER_INTERNAL_HEADERS__
#endif
#include "crt/host_runtime.h"
# 426 "/usr/lib/gcc/x86_64-linux-gnu/10/include/stddef.h" 3
typedef 
# 415 "/usr/lib/gcc/x86_64-linux-gnu/10/include/stddef.h" 3
struct { 
# 416
long long __max_align_ll __attribute((__aligned__(__alignof__(long long)))); 
# 417
long double __max_align_ld __attribute((__aligned__(__alignof__(long double)))); 
# 426 "/usr/lib/gcc/x86_64-linux-gnu/10/include/stddef.h" 3
} max_align_t; 
# 433
typedef __decltype((nullptr)) nullptr_t; 
# 203 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
#if 0
# 203
enum cudaError { 
# 210
cudaSuccess, 
# 216
cudaErrorInvalidValue, 
# 222
cudaErrorMemoryAllocation, 
# 228
cudaErrorInitializationError, 
# 235
cudaErrorCudartUnloading, 
# 242
cudaErrorProfilerDisabled, 
# 250
cudaErrorProfilerNotInitialized, 
# 257
cudaErrorProfilerAlreadyStarted, 
# 264
cudaErrorProfilerAlreadyStopped, 
# 273 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
cudaErrorInvalidConfiguration, 
# 279
cudaErrorInvalidPitchValue = 12, 
# 285
cudaErrorInvalidSymbol, 
# 293
cudaErrorInvalidHostPointer = 16, 
# 301
cudaErrorInvalidDevicePointer, 
# 307
cudaErrorInvalidTexture, 
# 313
cudaErrorInvalidTextureBinding, 
# 320
cudaErrorInvalidChannelDescriptor, 
# 326
cudaErrorInvalidMemcpyDirection, 
# 336 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
cudaErrorAddressOfConstant, 
# 345 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
cudaErrorTextureFetchFailed, 
# 354 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
cudaErrorTextureNotBound, 
# 363 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
cudaErrorSynchronizationError, 
# 369
cudaErrorInvalidFilterSetting, 
# 375
cudaErrorInvalidNormSetting, 
# 383
cudaErrorMixedDeviceExecution, 
# 391
cudaErrorNotYetImplemented = 31, 
# 400 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
cudaErrorMemoryValueTooLarge, 
# 407
cudaErrorStubLibrary = 34, 
# 414
cudaErrorInsufficientDriver, 
# 421
cudaErrorCallRequiresNewerDriver, 
# 427
cudaErrorInvalidSurface, 
# 433
cudaErrorDuplicateVariableName = 43, 
# 439
cudaErrorDuplicateTextureName, 
# 445
cudaErrorDuplicateSurfaceName, 
# 455 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
cudaErrorDevicesUnavailable, 
# 468 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
cudaErrorIncompatibleDriverContext = 49, 
# 474
cudaErrorMissingConfiguration = 52, 
# 483 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
cudaErrorPriorLaunchFailure, 
# 490
cudaErrorLaunchMaxDepthExceeded = 65, 
# 498
cudaErrorLaunchFileScopedTex, 
# 506
cudaErrorLaunchFileScopedSurf, 
# 522 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
cudaErrorSyncDepthExceeded, 
# 534 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
cudaErrorLaunchPendingCountExceeded, 
# 540
cudaErrorInvalidDeviceFunction = 98, 
# 546
cudaErrorNoDevice = 100, 
# 553
cudaErrorInvalidDevice, 
# 558
cudaErrorDeviceNotLicensed, 
# 567 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
cudaErrorSoftwareValidityNotEstablished, 
# 572
cudaErrorStartupFailure = 127, 
# 577
cudaErrorInvalidKernelImage = 200, 
# 587 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
cudaErrorDeviceUninitialized, 
# 592
cudaErrorMapBufferObjectFailed = 205, 
# 597
cudaErrorUnmapBufferObjectFailed, 
# 603
cudaErrorArrayIsMapped, 
# 608
cudaErrorAlreadyMapped, 
# 616
cudaErrorNoKernelImageForDevice, 
# 621
cudaErrorAlreadyAcquired, 
# 626
cudaErrorNotMapped, 
# 632
cudaErrorNotMappedAsArray, 
# 638
cudaErrorNotMappedAsPointer, 
# 644
cudaErrorECCUncorrectable, 
# 650
cudaErrorUnsupportedLimit, 
# 656
cudaErrorDeviceAlreadyInUse, 
# 662
cudaErrorPeerAccessUnsupported, 
# 668
cudaErrorInvalidPtx, 
# 673
cudaErrorInvalidGraphicsContext, 
# 679
cudaErrorNvlinkUncorrectable, 
# 686
cudaErrorJitCompilerNotFound, 
# 693
cudaErrorUnsupportedPtxVersion, 
# 700
cudaErrorJitCompilationDisabled, 
# 705
cudaErrorUnsupportedExecAffinity, 
# 711
cudaErrorUnsupportedDevSideSync, 
# 716
cudaErrorInvalidSource = 300, 
# 721
cudaErrorFileNotFound, 
# 726
cudaErrorSharedObjectSymbolNotFound, 
# 731
cudaErrorSharedObjectInitFailed, 
# 736
cudaErrorOperatingSystem, 
# 743
cudaErrorInvalidResourceHandle = 400, 
# 749
cudaErrorIllegalState, 
# 756
cudaErrorSymbolNotFound = 500, 
# 764
cudaErrorNotReady = 600, 
# 772
cudaErrorIllegalAddress = 700, 
# 781 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
cudaErrorLaunchOutOfResources, 
# 792 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
cudaErrorLaunchTimeout, 
# 798
cudaErrorLaunchIncompatibleTexturing, 
# 805
cudaErrorPeerAccessAlreadyEnabled, 
# 812
cudaErrorPeerAccessNotEnabled, 
# 825 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
cudaErrorSetOnActiveProcess = 708, 
# 832
cudaErrorContextIsDestroyed, 
# 839
cudaErrorAssert, 
# 846
cudaErrorTooManyPeers, 
# 852
cudaErrorHostMemoryAlreadyRegistered, 
# 858
cudaErrorHostMemoryNotRegistered, 
# 867 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
cudaErrorHardwareStackError, 
# 875
cudaErrorIllegalInstruction, 
# 884 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
cudaErrorMisalignedAddress, 
# 895 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
cudaErrorInvalidAddressSpace, 
# 903
cudaErrorInvalidPc, 
# 914 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
cudaErrorLaunchFailure, 
# 923 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
cudaErrorCooperativeLaunchTooLarge, 
# 928
cudaErrorNotPermitted = 800, 
# 934
cudaErrorNotSupported, 
# 943 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
cudaErrorSystemNotReady, 
# 950
cudaErrorSystemDriverMismatch, 
# 959 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
cudaErrorCompatNotSupportedOnDevice, 
# 964
cudaErrorMpsConnectionFailed, 
# 969
cudaErrorMpsRpcFailure, 
# 975
cudaErrorMpsServerNotReady, 
# 980
cudaErrorMpsMaxClientsReached, 
# 985
cudaErrorMpsMaxConnectionsReached, 
# 990
cudaErrorMpsClientTerminated, 
# 995
cudaErrorCdpNotSupported, 
# 1000
cudaErrorCdpVersionMismatch, 
# 1005
cudaErrorStreamCaptureUnsupported = 900, 
# 1011
cudaErrorStreamCaptureInvalidated, 
# 1017
cudaErrorStreamCaptureMerge, 
# 1022
cudaErrorStreamCaptureUnmatched, 
# 1028
cudaErrorStreamCaptureUnjoined, 
# 1035
cudaErrorStreamCaptureIsolation, 
# 1041
cudaErrorStreamCaptureImplicit, 
# 1047
cudaErrorCapturedEvent, 
# 1054
cudaErrorStreamCaptureWrongThread, 
# 1059
cudaErrorTimeout, 
# 1065
cudaErrorGraphExecUpdateFailure, 
# 1075 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
cudaErrorExternalDevice, 
# 1081
cudaErrorInvalidClusterSize, 
# 1086
cudaErrorUnknown = 999, 
# 1094
cudaErrorApiFailureBase = 10000
# 1095
}; 
#endif
# 1100 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
#if 0
# 1100
enum cudaChannelFormatKind { 
# 1102
cudaChannelFormatKindSigned, 
# 1103
cudaChannelFormatKindUnsigned, 
# 1104
cudaChannelFormatKindFloat, 
# 1105
cudaChannelFormatKindNone, 
# 1106
cudaChannelFormatKindNV12, 
# 1107
cudaChannelFormatKindUnsignedNormalized8X1, 
# 1108
cudaChannelFormatKindUnsignedNormalized8X2, 
# 1109
cudaChannelFormatKindUnsignedNormalized8X4, 
# 1110
cudaChannelFormatKindUnsignedNormalized16X1, 
# 1111
cudaChannelFormatKindUnsignedNormalized16X2, 
# 1112
cudaChannelFormatKindUnsignedNormalized16X4, 
# 1113
cudaChannelFormatKindSignedNormalized8X1, 
# 1114
cudaChannelFormatKindSignedNormalized8X2, 
# 1115
cudaChannelFormatKindSignedNormalized8X4, 
# 1116
cudaChannelFormatKindSignedNormalized16X1, 
# 1117
cudaChannelFormatKindSignedNormalized16X2, 
# 1118
cudaChannelFormatKindSignedNormalized16X4, 
# 1119
cudaChannelFormatKindUnsignedBlockCompressed1, 
# 1120
cudaChannelFormatKindUnsignedBlockCompressed1SRGB, 
# 1121
cudaChannelFormatKindUnsignedBlockCompressed2, 
# 1122
cudaChannelFormatKindUnsignedBlockCompressed2SRGB, 
# 1123
cudaChannelFormatKindUnsignedBlockCompressed3, 
# 1124
cudaChannelFormatKindUnsignedBlockCompressed3SRGB, 
# 1125
cudaChannelFormatKindUnsignedBlockCompressed4, 
# 1126
cudaChannelFormatKindSignedBlockCompressed4, 
# 1127
cudaChannelFormatKindUnsignedBlockCompressed5, 
# 1128
cudaChannelFormatKindSignedBlockCompressed5, 
# 1129
cudaChannelFormatKindUnsignedBlockCompressed6H, 
# 1130
cudaChannelFormatKindSignedBlockCompressed6H, 
# 1131
cudaChannelFormatKindUnsignedBlockCompressed7, 
# 1132
cudaChannelFormatKindUnsignedBlockCompressed7SRGB
# 1133
}; 
#endif
# 1138 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
#if 0
# 1138
struct cudaChannelFormatDesc { 
# 1140
int x; 
# 1141
int y; 
# 1142
int z; 
# 1143
int w; 
# 1144
cudaChannelFormatKind f; 
# 1145
}; 
#endif
# 1150 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
typedef struct cudaArray *cudaArray_t; 
# 1155
typedef const cudaArray *cudaArray_const_t; 
# 1157
struct cudaArray; 
# 1162
typedef struct cudaMipmappedArray *cudaMipmappedArray_t; 
# 1167
typedef const cudaMipmappedArray *cudaMipmappedArray_const_t; 
# 1169
struct cudaMipmappedArray; 
# 1179 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
#if 0
# 1179
struct cudaArraySparseProperties { 
# 1180
struct { 
# 1181
unsigned width; 
# 1182
unsigned height; 
# 1183
unsigned depth; 
# 1184
} tileExtent; 
# 1185
unsigned miptailFirstLevel; 
# 1186
unsigned long long miptailSize; 
# 1187
unsigned flags; 
# 1188
unsigned reserved[4]; 
# 1189
}; 
#endif
# 1194 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
#if 0
# 1194
struct cudaArrayMemoryRequirements { 
# 1195
size_t size; 
# 1196
size_t alignment; 
# 1197
unsigned reserved[4]; 
# 1198
}; 
#endif
# 1203 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
#if 0
# 1203
enum cudaMemoryType { 
# 1205
cudaMemoryTypeUnregistered, 
# 1206
cudaMemoryTypeHost, 
# 1207
cudaMemoryTypeDevice, 
# 1208
cudaMemoryTypeManaged
# 1209
}; 
#endif
# 1214 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
#if 0
# 1214
enum cudaMemcpyKind { 
# 1216
cudaMemcpyHostToHost, 
# 1217
cudaMemcpyHostToDevice, 
# 1218
cudaMemcpyDeviceToHost, 
# 1219
cudaMemcpyDeviceToDevice, 
# 1220
cudaMemcpyDefault
# 1221
}; 
#endif
# 1228 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
#if 0
# 1228
struct cudaPitchedPtr { 
# 1230
void *ptr; 
# 1231
size_t pitch; 
# 1232
size_t xsize; 
# 1233
size_t ysize; 
# 1234
}; 
#endif
# 1241 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
#if 0
# 1241
struct cudaExtent { 
# 1243
size_t width; 
# 1244
size_t height; 
# 1245
size_t depth; 
# 1246
}; 
#endif
# 1253 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
#if 0
# 1253
struct cudaPos { 
# 1255
size_t x; 
# 1256
size_t y; 
# 1257
size_t z; 
# 1258
}; 
#endif
# 1263 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
#if 0
# 1263
struct cudaMemcpy3DParms { 
# 1265
cudaArray_t srcArray; 
# 1266
cudaPos srcPos; 
# 1267
cudaPitchedPtr srcPtr; 
# 1269
cudaArray_t dstArray; 
# 1270
cudaPos dstPos; 
# 1271
cudaPitchedPtr dstPtr; 
# 1273
cudaExtent extent; 
# 1274
cudaMemcpyKind kind; __pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)
# 1275
}; 
#endif
# 1280 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
#if 0
# 1280
struct cudaMemcpyNodeParams { 
# 1281
int flags; 
# 1282
int reserved[3]; 
# 1283
cudaMemcpy3DParms copyParams; 
# 1284
}; 
#endif
# 1289 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
#if 0
# 1289
struct cudaMemcpy3DPeerParms { 
# 1291
cudaArray_t srcArray; 
# 1292
cudaPos srcPos; 
# 1293
cudaPitchedPtr srcPtr; 
# 1294
int srcDevice; 
# 1296
cudaArray_t dstArray; 
# 1297
cudaPos dstPos; 
# 1298
cudaPitchedPtr dstPtr; 
# 1299
int dstDevice; 
# 1301
cudaExtent extent; 
# 1302
}; 
#endif
# 1307 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
#if 0
# 1307
struct cudaMemsetParams { 
# 1308
void *dst; 
# 1309
size_t pitch; 
# 1310
unsigned value; 
# 1311
unsigned elementSize; 
# 1312
size_t width; 
# 1313
size_t height; 
# 1314
}; 
#endif
# 1319 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
#if 0
# 1319
struct cudaMemsetParamsV2 { 
# 1320
void *dst; 
# 1321
size_t pitch; 
# 1322
unsigned value; 
# 1323
unsigned elementSize; 
# 1324
size_t width; 
# 1325
size_t height; 
# 1326
}; 
#endif
# 1331 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
#if 0
# 1331
enum cudaAccessProperty { 
# 1332
cudaAccessPropertyNormal, 
# 1333
cudaAccessPropertyStreaming, 
# 1334
cudaAccessPropertyPersisting
# 1335
}; 
#endif
# 1348 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
#if 0
# 1348
struct cudaAccessPolicyWindow { 
# 1349
void *base_ptr; 
# 1350
size_t num_bytes; 
# 1351
float hitRatio; 
# 1352
cudaAccessProperty hitProp; 
# 1353
cudaAccessProperty missProp; __pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)
# 1354
}; 
#endif
# 1366 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
typedef void (*cudaHostFn_t)(void * userData); 
# 1371
#if 0
# 1371
struct cudaHostNodeParams { 
# 1372
cudaHostFn_t fn; 
# 1373
void *userData; 
# 1374
}; 
#endif
# 1379 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
#if 0
# 1379
struct cudaHostNodeParamsV2 { 
# 1380
cudaHostFn_t fn; 
# 1381
void *userData; 
# 1382
}; 
#endif
# 1387 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
#if 0
# 1387
enum cudaStreamCaptureStatus { 
# 1388
cudaStreamCaptureStatusNone, 
# 1389
cudaStreamCaptureStatusActive, 
# 1390
cudaStreamCaptureStatusInvalidated
# 1392
}; 
#endif
# 1398 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
#if 0
# 1398
enum cudaStreamCaptureMode { 
# 1399
cudaStreamCaptureModeGlobal, 
# 1400
cudaStreamCaptureModeThreadLocal, 
# 1401
cudaStreamCaptureModeRelaxed
# 1402
}; 
#endif
# 1404 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
#if 0
# 1404
enum cudaSynchronizationPolicy { 
# 1405
cudaSyncPolicyAuto = 1, 
# 1406
cudaSyncPolicySpin, 
# 1407
cudaSyncPolicyYield, 
# 1408
cudaSyncPolicyBlockingSync
# 1409
}; 
#endif
# 1414 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
#if 0
# 1414
enum cudaClusterSchedulingPolicy { 
# 1415
cudaClusterSchedulingPolicyDefault, 
# 1416
cudaClusterSchedulingPolicySpread, 
# 1417
cudaClusterSchedulingPolicyLoadBalancing
# 1418
}; 
#endif
# 1423 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
#if 0
# 1423
enum cudaStreamUpdateCaptureDependenciesFlags { 
# 1424
cudaStreamAddCaptureDependencies, 
# 1425
cudaStreamSetCaptureDependencies
# 1426
}; 
#endif
# 1431 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
#if 0
# 1431
enum cudaUserObjectFlags { 
# 1432
cudaUserObjectNoDestructorSync = 1
# 1433
}; 
#endif
# 1438 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
#if 0
# 1438
enum cudaUserObjectRetainFlags { 
# 1439
cudaGraphUserObjectMove = 1
# 1440
}; 
#endif
# 1445 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
struct cudaGraphicsResource; 
# 1450
#if 0
# 1450
enum cudaGraphicsRegisterFlags { 
# 1452
cudaGraphicsRegisterFlagsNone, 
# 1453
cudaGraphicsRegisterFlagsReadOnly, 
# 1454
cudaGraphicsRegisterFlagsWriteDiscard, 
# 1455
cudaGraphicsRegisterFlagsSurfaceLoadStore = 4, 
# 1456
cudaGraphicsRegisterFlagsTextureGather = 8
# 1457
}; 
#endif
# 1462 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
#if 0
# 1462
enum cudaGraphicsMapFlags { 
# 1464
cudaGraphicsMapFlagsNone, 
# 1465
cudaGraphicsMapFlagsReadOnly, 
# 1466
cudaGraphicsMapFlagsWriteDiscard
# 1467
}; 
#endif
# 1472 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
#if 0
# 1472
enum cudaGraphicsCubeFace { 
# 1474
cudaGraphicsCubeFacePositiveX, 
# 1475
cudaGraphicsCubeFaceNegativeX, 
# 1476
cudaGraphicsCubeFacePositiveY, 
# 1477
cudaGraphicsCubeFaceNegativeY, 
# 1478
cudaGraphicsCubeFacePositiveZ, 
# 1479
cudaGraphicsCubeFaceNegativeZ
# 1480
}; 
#endif
# 1485 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
#if 0
# 1485
enum cudaResourceType { 
# 1487
cudaResourceTypeArray, 
# 1488
cudaResourceTypeMipmappedArray, 
# 1489
cudaResourceTypeLinear, 
# 1490
cudaResourceTypePitch2D
# 1491
}; 
#endif
# 1496 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
#if 0
# 1496
enum cudaResourceViewFormat { 
# 1498
cudaResViewFormatNone, 
# 1499
cudaResViewFormatUnsignedChar1, 
# 1500
cudaResViewFormatUnsignedChar2, 
# 1501
cudaResViewFormatUnsignedChar4, 
# 1502
cudaResViewFormatSignedChar1, 
# 1503
cudaResViewFormatSignedChar2, 
# 1504
cudaResViewFormatSignedChar4, 
# 1505
cudaResViewFormatUnsignedShort1, 
# 1506
cudaResViewFormatUnsignedShort2, 
# 1507
cudaResViewFormatUnsignedShort4, 
# 1508
cudaResViewFormatSignedShort1, 
# 1509
cudaResViewFormatSignedShort2, 
# 1510
cudaResViewFormatSignedShort4, 
# 1511
cudaResViewFormatUnsignedInt1, 
# 1512
cudaResViewFormatUnsignedInt2, 
# 1513
cudaResViewFormatUnsignedInt4, 
# 1514
cudaResViewFormatSignedInt1, 
# 1515
cudaResViewFormatSignedInt2, 
# 1516
cudaResViewFormatSignedInt4, 
# 1517
cudaResViewFormatHalf1, 
# 1518
cudaResViewFormatHalf2, 
# 1519
cudaResViewFormatHalf4, 
# 1520
cudaResViewFormatFloat1, 
# 1521
cudaResViewFormatFloat2, 
# 1522
cudaResViewFormatFloat4, 
# 1523
cudaResViewFormatUnsignedBlockCompressed1, 
# 1524
cudaResViewFormatUnsignedBlockCompressed2, 
# 1525
cudaResViewFormatUnsignedBlockCompressed3, 
# 1526
cudaResViewFormatUnsignedBlockCompressed4, 
# 1527
cudaResViewFormatSignedBlockCompressed4, 
# 1528
cudaResViewFormatUnsignedBlockCompressed5, 
# 1529
cudaResViewFormatSignedBlockCompressed5, 
# 1530
cudaResViewFormatUnsignedBlockCompressed6H, 
# 1531
cudaResViewFormatSignedBlockCompressed6H, 
# 1532
cudaResViewFormatUnsignedBlockCompressed7
# 1533
}; 
#endif
# 1538 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
#if 0
# 1538
struct cudaResourceDesc { 
# 1539
cudaResourceType resType; 
# 1541
union { 
# 1542
struct { 
# 1543
cudaArray_t array; 
# 1544
} array; 
# 1545
struct { 
# 1546
cudaMipmappedArray_t mipmap; 
# 1547
} mipmap; 
# 1548
struct { 
# 1549
void *devPtr; 
# 1550
cudaChannelFormatDesc desc; 
# 1551
size_t sizeInBytes; 
# 1552
} linear; 
# 1553
struct { 
# 1554
void *devPtr; 
# 1555
cudaChannelFormatDesc desc; 
# 1556
size_t width; 
# 1557
size_t height; 
# 1558
size_t pitchInBytes; 
# 1559
} pitch2D; 
# 1560
} res; 
# 1561
}; 
#endif
# 1566 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
#if 0
# 1566
struct cudaResourceViewDesc { 
# 1568
cudaResourceViewFormat format; 
# 1569
size_t width; 
# 1570
size_t height; 
# 1571
size_t depth; 
# 1572
unsigned firstMipmapLevel; 
# 1573
unsigned lastMipmapLevel; 
# 1574
unsigned firstLayer; 
# 1575
unsigned lastLayer; 
# 1576
}; 
#endif
# 1581 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
#if 0
# 1581
struct cudaPointerAttributes { 
# 1587
cudaMemoryType type; 
# 1598 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
int device; 
# 1604
void *devicePointer; 
# 1613 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
void *hostPointer; 
# 1614
}; 
#endif
# 1619 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
#if 0
# 1619
struct cudaFuncAttributes { 
# 1626
size_t sharedSizeBytes; 
# 1632
size_t constSizeBytes; 
# 1637
size_t localSizeBytes; 
# 1644
int maxThreadsPerBlock; 
# 1649
int numRegs; 
# 1656
int ptxVersion; 
# 1663
int binaryVersion; 
# 1669
int cacheModeCA; 
# 1676
int maxDynamicSharedSizeBytes; 
# 1685 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
int preferredShmemCarveout; 
# 1691
int clusterDimMustBeSet; 
# 1702 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
int requiredClusterWidth; 
# 1703
int requiredClusterHeight; 
# 1704
int requiredClusterDepth; 
# 1710
int clusterSchedulingPolicyPreference; 
# 1732 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
int nonPortableClusterSizeAllowed; 
# 1737
int reserved[16]; 
# 1738
}; 
#endif
# 1743 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
#if 0
# 1743
enum cudaFuncAttribute { 
# 1745
cudaFuncAttributeMaxDynamicSharedMemorySize = 8, 
# 1746
cudaFuncAttributePreferredSharedMemoryCarveout, 
# 1747
cudaFuncAttributeClusterDimMustBeSet, 
# 1748
cudaFuncAttributeRequiredClusterWidth, 
# 1749
cudaFuncAttributeRequiredClusterHeight, 
# 1750
cudaFuncAttributeRequiredClusterDepth, 
# 1751
cudaFuncAttributeNonPortableClusterSizeAllowed, 
# 1752
cudaFuncAttributeClusterSchedulingPolicyPreference, 
# 1753
cudaFuncAttributeMax
# 1754
}; 
#endif
# 1759 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
#if 0
# 1759
enum cudaFuncCache { 
# 1761
cudaFuncCachePreferNone, 
# 1762
cudaFuncCachePreferShared, 
# 1763
cudaFuncCachePreferL1, 
# 1764
cudaFuncCachePreferEqual
# 1765
}; 
#endif
# 1771 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
#if 0
# 1771
enum cudaSharedMemConfig { 
# 1773
cudaSharedMemBankSizeDefault, 
# 1774
cudaSharedMemBankSizeFourByte, 
# 1775
cudaSharedMemBankSizeEightByte
# 1776
}; 
#endif
# 1781 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
#if 0
# 1781
enum cudaSharedCarveout { 
# 1782
cudaSharedmemCarveoutDefault = (-1), 
# 1783
cudaSharedmemCarveoutMaxShared = 100, 
# 1784
cudaSharedmemCarveoutMaxL1 = 0
# 1785
}; 
#endif
# 1790 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
#if 0
# 1790
enum cudaComputeMode { 
# 1792
cudaComputeModeDefault, 
# 1793
cudaComputeModeExclusive, 
# 1794
cudaComputeModeProhibited, 
# 1795
cudaComputeModeExclusiveProcess
# 1796
}; 
#endif
# 1801 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
#if 0
# 1801
enum cudaLimit { 
# 1803
cudaLimitStackSize, 
# 1804
cudaLimitPrintfFifoSize, 
# 1805
cudaLimitMallocHeapSize, 
# 1806
cudaLimitDevRuntimeSyncDepth, 
# 1807
cudaLimitDevRuntimePendingLaunchCount, 
# 1808
cudaLimitMaxL2FetchGranularity, 
# 1809
cudaLimitPersistingL2CacheSize
# 1810
}; 
#endif
# 1815 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
#if 0
# 1815
enum cudaMemoryAdvise { 
# 1817
cudaMemAdviseSetReadMostly = 1, 
# 1818
cudaMemAdviseUnsetReadMostly, 
# 1819
cudaMemAdviseSetPreferredLocation, 
# 1820
cudaMemAdviseUnsetPreferredLocation, 
# 1821
cudaMemAdviseSetAccessedBy, 
# 1822
cudaMemAdviseUnsetAccessedBy
# 1823
}; 
#endif
# 1828 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
#if 0
# 1828
enum cudaMemRangeAttribute { 
# 1830
cudaMemRangeAttributeReadMostly = 1, 
# 1831
cudaMemRangeAttributePreferredLocation, 
# 1832
cudaMemRangeAttributeAccessedBy, 
# 1833
cudaMemRangeAttributeLastPrefetchLocation, 
# 1834
cudaMemRangeAttributePreferredLocationType, 
# 1835
cudaMemRangeAttributePreferredLocationId, 
# 1836
cudaMemRangeAttributeLastPrefetchLocationType, 
# 1837
cudaMemRangeAttributeLastPrefetchLocationId
# 1838
}; 
#endif
# 1843 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
#if 0
# 1843
enum cudaFlushGPUDirectRDMAWritesOptions { 
# 1844
cudaFlushGPUDirectRDMAWritesOptionHost = (1 << 0), 
# 1845
cudaFlushGPUDirectRDMAWritesOptionMemOps
# 1846
}; 
#endif
# 1851 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
#if 0
# 1851
enum cudaGPUDirectRDMAWritesOrdering { 
# 1852
cudaGPUDirectRDMAWritesOrderingNone, 
# 1853
cudaGPUDirectRDMAWritesOrderingOwner = 100, 
# 1854
cudaGPUDirectRDMAWritesOrderingAllDevices = 200
# 1855
}; 
#endif
# 1860 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
#if 0
# 1860
enum cudaFlushGPUDirectRDMAWritesScope { 
# 1861
cudaFlushGPUDirectRDMAWritesToOwner = 100, 
# 1862
cudaFlushGPUDirectRDMAWritesToAllDevices = 200
# 1863
}; 
#endif
# 1868 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
#if 0
# 1868
enum cudaFlushGPUDirectRDMAWritesTarget { 
# 1869
cudaFlushGPUDirectRDMAWritesTargetCurrentDevice
# 1870
}; 
#endif
# 1876 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
#if 0
# 1876
enum cudaDeviceAttr { 
# 1878
cudaDevAttrMaxThreadsPerBlock = 1, 
# 1879
cudaDevAttrMaxBlockDimX, 
# 1880
cudaDevAttrMaxBlockDimY, 
# 1881
cudaDevAttrMaxBlockDimZ, 
# 1882
cudaDevAttrMaxGridDimX, 
# 1883
cudaDevAttrMaxGridDimY, 
# 1884
cudaDevAttrMaxGridDimZ, 
# 1885
cudaDevAttrMaxSharedMemoryPerBlock, 
# 1886
cudaDevAttrTotalConstantMemory, 
# 1887
cudaDevAttrWarpSize, 
# 1888
cudaDevAttrMaxPitch, 
# 1889
cudaDevAttrMaxRegistersPerBlock, 
# 1890
cudaDevAttrClockRate, 
# 1891
cudaDevAttrTextureAlignment, 
# 1892
cudaDevAttrGpuOverlap, 
# 1893
cudaDevAttrMultiProcessorCount, 
# 1894
cudaDevAttrKernelExecTimeout, 
# 1895
cudaDevAttrIntegrated, 
# 1896
cudaDevAttrCanMapHostMemory, 
# 1897
cudaDevAttrComputeMode, 
# 1898
cudaDevAttrMaxTexture1DWidth, 
# 1899
cudaDevAttrMaxTexture2DWidth, 
# 1900
cudaDevAttrMaxTexture2DHeight, 
# 1901
cudaDevAttrMaxTexture3DWidth, 
# 1902
cudaDevAttrMaxTexture3DHeight, 
# 1903
cudaDevAttrMaxTexture3DDepth, 
# 1904
cudaDevAttrMaxTexture2DLayeredWidth, 
# 1905
cudaDevAttrMaxTexture2DLayeredHeight, 
# 1906
cudaDevAttrMaxTexture2DLayeredLayers, 
# 1907
cudaDevAttrSurfaceAlignment, 
# 1908
cudaDevAttrConcurrentKernels, 
# 1909
cudaDevAttrEccEnabled, 
# 1910
cudaDevAttrPciBusId, 
# 1911
cudaDevAttrPciDeviceId, 
# 1912
cudaDevAttrTccDriver, 
# 1913
cudaDevAttrMemoryClockRate, 
# 1914
cudaDevAttrGlobalMemoryBusWidth, 
# 1915
cudaDevAttrL2CacheSize, 
# 1916
cudaDevAttrMaxThreadsPerMultiProcessor, 
# 1917
cudaDevAttrAsyncEngineCount, 
# 1918
cudaDevAttrUnifiedAddressing, 
# 1919
cudaDevAttrMaxTexture1DLayeredWidth, 
# 1920
cudaDevAttrMaxTexture1DLayeredLayers, 
# 1921
cudaDevAttrMaxTexture2DGatherWidth = 45, 
# 1922
cudaDevAttrMaxTexture2DGatherHeight, 
# 1923
cudaDevAttrMaxTexture3DWidthAlt, 
# 1924
cudaDevAttrMaxTexture3DHeightAlt, 
# 1925
cudaDevAttrMaxTexture3DDepthAlt, 
# 1926
cudaDevAttrPciDomainId, 
# 1927
cudaDevAttrTexturePitchAlignment, 
# 1928
cudaDevAttrMaxTextureCubemapWidth, 
# 1929
cudaDevAttrMaxTextureCubemapLayeredWidth, 
# 1930
cudaDevAttrMaxTextureCubemapLayeredLayers, 
# 1931
cudaDevAttrMaxSurface1DWidth, 
# 1932
cudaDevAttrMaxSurface2DWidth, 
# 1933
cudaDevAttrMaxSurface2DHeight, 
# 1934
cudaDevAttrMaxSurface3DWidth, 
# 1935
cudaDevAttrMaxSurface3DHeight, 
# 1936
cudaDevAttrMaxSurface3DDepth, 
# 1937
cudaDevAttrMaxSurface1DLayeredWidth, 
# 1938
cudaDevAttrMaxSurface1DLayeredLayers, 
# 1939
cudaDevAttrMaxSurface2DLayeredWidth, 
# 1940
cudaDevAttrMaxSurface2DLayeredHeight, 
# 1941
cudaDevAttrMaxSurface2DLayeredLayers, 
# 1942
cudaDevAttrMaxSurfaceCubemapWidth, 
# 1943
cudaDevAttrMaxSurfaceCubemapLayeredWidth, 
# 1944
cudaDevAttrMaxSurfaceCubemapLayeredLayers, 
# 1945
cudaDevAttrMaxTexture1DLinearWidth, 
# 1946
cudaDevAttrMaxTexture2DLinearWidth, 
# 1947
cudaDevAttrMaxTexture2DLinearHeight, 
# 1948
cudaDevAttrMaxTexture2DLinearPitch, 
# 1949
cudaDevAttrMaxTexture2DMipmappedWidth, 
# 1950
cudaDevAttrMaxTexture2DMipmappedHeight, 
# 1951
cudaDevAttrComputeCapabilityMajor, 
# 1952
cudaDevAttrComputeCapabilityMinor, 
# 1953
cudaDevAttrMaxTexture1DMipmappedWidth, 
# 1954
cudaDevAttrStreamPrioritiesSupported, 
# 1955
cudaDevAttrGlobalL1CacheSupported, 
# 1956
cudaDevAttrLocalL1CacheSupported, 
# 1957
cudaDevAttrMaxSharedMemoryPerMultiprocessor, 
# 1958
cudaDevAttrMaxRegistersPerMultiprocessor, 
# 1959
cudaDevAttrManagedMemory, 
# 1960
cudaDevAttrIsMultiGpuBoard, 
# 1961
cudaDevAttrMultiGpuBoardGroupID, 
# 1962
cudaDevAttrHostNativeAtomicSupported, 
# 1963
cudaDevAttrSingleToDoublePrecisionPerfRatio, 
# 1964
cudaDevAttrPageableMemoryAccess, 
# 1965
cudaDevAttrConcurrentManagedAccess, 
# 1966
cudaDevAttrComputePreemptionSupported, 
# 1967
cudaDevAttrCanUseHostPointerForRegisteredMem, 
# 1968
cudaDevAttrReserved92, 
# 1969
cudaDevAttrReserved93, 
# 1970
cudaDevAttrReserved94, 
# 1971
cudaDevAttrCooperativeLaunch, 
# 1972
cudaDevAttrCooperativeMultiDeviceLaunch, 
# 1973
cudaDevAttrMaxSharedMemoryPerBlockOptin, 
# 1974
cudaDevAttrCanFlushRemoteWrites, 
# 1975
cudaDevAttrHostRegisterSupported, 
# 1976
cudaDevAttrPageableMemoryAccessUsesHostPageTables, 
# 1977
cudaDevAttrDirectManagedMemAccessFromHost, 
# 1978
cudaDevAttrMaxBlocksPerMultiprocessor = 106, 
# 1979
cudaDevAttrMaxPersistingL2CacheSize = 108, 
# 1980
cudaDevAttrMaxAccessPolicyWindowSize, 
# 1981
cudaDevAttrReservedSharedMemoryPerBlock = 111, 
# 1982
cudaDevAttrSparseCudaArraySupported, 
# 1983
cudaDevAttrHostRegisterReadOnlySupported, 
# 1984
cudaDevAttrTimelineSemaphoreInteropSupported, 
# 1985
cudaDevAttrMaxTimelineSemaphoreInteropSupported = 114, 
# 1986
cudaDevAttrMemoryPoolsSupported, 
# 1987
cudaDevAttrGPUDirectRDMASupported, 
# 1988
cudaDevAttrGPUDirectRDMAFlushWritesOptions, 
# 1989
cudaDevAttrGPUDirectRDMAWritesOrdering, 
# 1990
cudaDevAttrMemoryPoolSupportedHandleTypes, 
# 1991
cudaDevAttrClusterLaunch, 
# 1992
cudaDevAttrDeferredMappingCudaArraySupported, 
# 1993
cudaDevAttrReserved122, 
# 1994
cudaDevAttrReserved123, 
# 1995
cudaDevAttrReserved124, 
# 1996
cudaDevAttrIpcEventSupport, 
# 1997
cudaDevAttrMemSyncDomainCount, 
# 1998
cudaDevAttrReserved127, 
# 1999
cudaDevAttrReserved128, 
# 2000
cudaDevAttrReserved129, 
# 2001
cudaDevAttrNumaConfig, 
# 2002
cudaDevAttrNumaId, 
# 2003
cudaDevAttrReserved132, 
# 2004
cudaDevAttrHostNumaId = 134, 
# 2005
cudaDevAttrMax
# 2006
}; 
#endif
# 2011 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
#if 0
# 2011
enum cudaMemPoolAttr { 
# 2021 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
cudaMemPoolReuseFollowEventDependencies = 1, 
# 2028
cudaMemPoolReuseAllowOpportunistic, 
# 2036
cudaMemPoolReuseAllowInternalDependencies, 
# 2047 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
cudaMemPoolAttrReleaseThreshold, 
# 2053
cudaMemPoolAttrReservedMemCurrent, 
# 2060
cudaMemPoolAttrReservedMemHigh, 
# 2066
cudaMemPoolAttrUsedMemCurrent, 
# 2073
cudaMemPoolAttrUsedMemHigh
# 2074
}; 
#endif
# 2079 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
#if 0
# 2079
enum cudaMemLocationType { 
# 2080
cudaMemLocationTypeInvalid, 
# 2081
cudaMemLocationTypeDevice, 
# 2082
cudaMemLocationTypeHost, 
# 2083
cudaMemLocationTypeHostNuma, 
# 2084
cudaMemLocationTypeHostNumaCurrent
# 2085
}; 
#endif
# 2093 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
#if 0
# 2093
struct cudaMemLocation { 
# 2094
cudaMemLocationType type; 
# 2095
int id; 
# 2096
}; 
#endif
# 2101 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
#if 0
# 2101
enum cudaMemAccessFlags { 
# 2102
cudaMemAccessFlagsProtNone, 
# 2103
cudaMemAccessFlagsProtRead, 
# 2104
cudaMemAccessFlagsProtReadWrite = 3
# 2105
}; 
#endif
# 2110 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
#if 0
# 2110
struct cudaMemAccessDesc { 
# 2111
cudaMemLocation location; 
# 2112
cudaMemAccessFlags flags; 
# 2113
}; 
#endif
# 2118 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
#if 0
# 2118
enum cudaMemAllocationType { 
# 2119
cudaMemAllocationTypeInvalid, 
# 2123
cudaMemAllocationTypePinned, 
# 2124
cudaMemAllocationTypeMax = 2147483647
# 2125
}; 
#endif
# 2130 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
#if 0
# 2130
enum cudaMemAllocationHandleType { 
# 2131
cudaMemHandleTypeNone, 
# 2132
cudaMemHandleTypePosixFileDescriptor, 
# 2133
cudaMemHandleTypeWin32, 
# 2134
cudaMemHandleTypeWin32Kmt = 4
# 2135
}; 
#endif
# 2140 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
#if 0
# 2140
struct cudaMemPoolProps { 
# 2141
cudaMemAllocationType allocType; 
# 2142
cudaMemAllocationHandleType handleTypes; 
# 2143
cudaMemLocation location; 
# 2150
void *win32SecurityAttributes; 
# 2151
size_t maxSize; 
# 2152
unsigned char reserved[56]; 
# 2153
}; 
#endif
# 2158 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
#if 0
# 2158
struct cudaMemPoolPtrExportData { 
# 2159
unsigned char reserved[64]; 
# 2160
}; 
#endif
# 2165 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
#if 0
# 2165
struct cudaMemAllocNodeParams { 
# 2170
cudaMemPoolProps poolProps; 
# 2171
const cudaMemAccessDesc *accessDescs; 
# 2172
size_t accessDescCount; 
# 2173
size_t bytesize; 
# 2174
void *dptr; 
# 2175
}; 
#endif
# 2180 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
#if 0
# 2180
struct cudaMemAllocNodeParamsV2 { 
# 2185
cudaMemPoolProps poolProps; 
# 2186
const cudaMemAccessDesc *accessDescs; 
# 2187
size_t accessDescCount; 
# 2188
size_t bytesize; 
# 2189
void *dptr; 
# 2190
}; 
#endif
# 2195 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
#if 0
# 2195
struct cudaMemFreeNodeParams { 
# 2196
void *dptr; 
# 2197
}; 
#endif
# 2202 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
#if 0
# 2202
enum cudaGraphMemAttributeType { 
# 2207
cudaGraphMemAttrUsedMemCurrent, 
# 2214
cudaGraphMemAttrUsedMemHigh, 
# 2221
cudaGraphMemAttrReservedMemCurrent, 
# 2228
cudaGraphMemAttrReservedMemHigh
# 2229
}; 
#endif
# 2235 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
#if 0
# 2235
enum cudaDeviceP2PAttr { 
# 2236
cudaDevP2PAttrPerformanceRank = 1, 
# 2237
cudaDevP2PAttrAccessSupported, 
# 2238
cudaDevP2PAttrNativeAtomicSupported, 
# 2239
cudaDevP2PAttrCudaArrayAccessSupported
# 2240
}; 
#endif
# 2247 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
#if 0
# 2247
struct CUuuid_st { 
# 2248
char bytes[16]; 
# 2249
}; 
#endif
# 2250 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
#if 0
typedef CUuuid_st 
# 2250
CUuuid; 
#endif
# 2252 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
#if 0
typedef CUuuid_st 
# 2252
cudaUUID_t; 
#endif
# 2257 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
#if 0
# 2257
struct cudaDeviceProp { 
# 2259
char name[256]; 
# 2260
cudaUUID_t uuid; 
# 2261
char luid[8]; 
# 2262
unsigned luidDeviceNodeMask; 
# 2263
size_t totalGlobalMem; 
# 2264
size_t sharedMemPerBlock; 
# 2265
int regsPerBlock; 
# 2266
int warpSize; 
# 2267
size_t memPitch; 
# 2268
int maxThreadsPerBlock; 
# 2269
int maxThreadsDim[3]; 
# 2270
int maxGridSize[3]; 
# 2271
int clockRate; 
# 2272
size_t totalConstMem; 
# 2273
int major; 
# 2274
int minor; 
# 2275
size_t textureAlignment; 
# 2276
size_t texturePitchAlignment; 
# 2277
int deviceOverlap; 
# 2278
int multiProcessorCount; 
# 2279
int kernelExecTimeoutEnabled; 
# 2280
int integrated; 
# 2281
int canMapHostMemory; 
# 2282
int computeMode; 
# 2283
int maxTexture1D; 
# 2284
int maxTexture1DMipmap; 
# 2285
int maxTexture1DLinear; 
# 2286
int maxTexture2D[2]; 
# 2287
int maxTexture2DMipmap[2]; 
# 2288
int maxTexture2DLinear[3]; 
# 2289
int maxTexture2DGather[2]; 
# 2290
int maxTexture3D[3]; 
# 2291
int maxTexture3DAlt[3]; 
# 2292
int maxTextureCubemap; 
# 2293
int maxTexture1DLayered[2]; 
# 2294
int maxTexture2DLayered[3]; 
# 2295
int maxTextureCubemapLayered[2]; 
# 2296
int maxSurface1D; 
# 2297
int maxSurface2D[2]; 
# 2298
int maxSurface3D[3]; 
# 2299
int maxSurface1DLayered[2]; 
# 2300
int maxSurface2DLayered[3]; 
# 2301
int maxSurfaceCubemap; 
# 2302
int maxSurfaceCubemapLayered[2]; 
# 2303
size_t surfaceAlignment; 
# 2304
int concurrentKernels; 
# 2305
int ECCEnabled; 
# 2306
int pciBusID; 
# 2307
int pciDeviceID; 
# 2308
int pciDomainID; 
# 2309
int tccDriver; 
# 2310
int asyncEngineCount; 
# 2311
int unifiedAddressing; 
# 2312
int memoryClockRate; 
# 2313
int memoryBusWidth; 
# 2314
int l2CacheSize; 
# 2315
int persistingL2CacheMaxSize; 
# 2316
int maxThreadsPerMultiProcessor; 
# 2317
int streamPrioritiesSupported; 
# 2318
int globalL1CacheSupported; 
# 2319
int localL1CacheSupported; 
# 2320
size_t sharedMemPerMultiprocessor; 
# 2321
int regsPerMultiprocessor; 
# 2322
int managedMemory; 
# 2323
int isMultiGpuBoard; 
# 2324
int multiGpuBoardGroupID; 
# 2325
int hostNativeAtomicSupported; 
# 2326
int singleToDoublePrecisionPerfRatio; 
# 2327
int pageableMemoryAccess; 
# 2328
int concurrentManagedAccess; 
# 2329
int computePreemptionSupported; 
# 2330
int canUseHostPointerForRegisteredMem; 
# 2331
int cooperativeLaunch; 
# 2332
int cooperativeMultiDeviceLaunch; 
# 2333
size_t sharedMemPerBlockOptin; 
# 2334
int pageableMemoryAccessUsesHostPageTables; 
# 2335
int directManagedMemAccessFromHost; 
# 2336
int maxBlocksPerMultiProcessor; 
# 2337
int accessPolicyMaxWindowSize; 
# 2338
size_t reservedSharedMemPerBlock; 
# 2339
int hostRegisterSupported; 
# 2340
int sparseCudaArraySupported; 
# 2341
int hostRegisterReadOnlySupported; 
# 2342
int timelineSemaphoreInteropSupported; 
# 2343
int memoryPoolsSupported; 
# 2344
int gpuDirectRDMASupported; 
# 2345
unsigned gpuDirectRDMAFlushWritesOptions; 
# 2346
int gpuDirectRDMAWritesOrdering; 
# 2347
unsigned memoryPoolSupportedHandleTypes; 
# 2348
int deferredMappingCudaArraySupported; 
# 2349
int ipcEventSupported; 
# 2350
int clusterLaunch; 
# 2351
int unifiedFunctionPointers; 
# 2352
int reserved2[2]; 
# 2353
int reserved[61]; 
# 2354
}; 
#endif
# 2367 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
#if 0
typedef 
# 2364
struct cudaIpcEventHandle_st { 
# 2366
char reserved[64]; 
# 2367
} cudaIpcEventHandle_t; 
#endif
# 2375 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
#if 0
typedef 
# 2372
struct cudaIpcMemHandle_st { 
# 2374
char reserved[64]; 
# 2375
} cudaIpcMemHandle_t; 
#endif
# 2380 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
#if 0
# 2380
enum cudaExternalMemoryHandleType { 
# 2384
cudaExternalMemoryHandleTypeOpaqueFd = 1, 
# 2388
cudaExternalMemoryHandleTypeOpaqueWin32, 
# 2392
cudaExternalMemoryHandleTypeOpaqueWin32Kmt, 
# 2396
cudaExternalMemoryHandleTypeD3D12Heap, 
# 2400
cudaExternalMemoryHandleTypeD3D12Resource, 
# 2404
cudaExternalMemoryHandleTypeD3D11Resource, 
# 2408
cudaExternalMemoryHandleTypeD3D11ResourceKmt, 
# 2412
cudaExternalMemoryHandleTypeNvSciBuf
# 2413
}; 
#endif
# 2455 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
#if 0
# 2455
struct cudaExternalMemoryHandleDesc { 
# 2459
cudaExternalMemoryHandleType type; 
# 2460
union { 
# 2466
int fd; 
# 2482 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
struct { 
# 2486
void *handle; 
# 2491
const void *name; 
# 2492
} win32; 
# 2497
const void *nvSciBufObject; 
# 2498
} handle; 
# 2502
unsigned long long size; 
# 2506
unsigned flags; __pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)
# 2507
}; 
#endif
# 2512 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
#if 0
# 2512
struct cudaExternalMemoryBufferDesc { 
# 2516
unsigned long long offset; 
# 2520
unsigned long long size; 
# 2524
unsigned flags; 
# 2525
}; 
#endif
# 2530 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
#if 0
# 2530
struct cudaExternalMemoryMipmappedArrayDesc { 
# 2535
unsigned long long offset; 
# 2539
cudaChannelFormatDesc formatDesc; 
# 2543
cudaExtent extent; 
# 2548
unsigned flags; 
# 2552
unsigned numLevels; 
# 2553
}; 
#endif
# 2558 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
#if 0
# 2558
enum cudaExternalSemaphoreHandleType { 
# 2562
cudaExternalSemaphoreHandleTypeOpaqueFd = 1, 
# 2566
cudaExternalSemaphoreHandleTypeOpaqueWin32, 
# 2570
cudaExternalSemaphoreHandleTypeOpaqueWin32Kmt, 
# 2574
cudaExternalSemaphoreHandleTypeD3D12Fence, 
# 2578
cudaExternalSemaphoreHandleTypeD3D11Fence, 
# 2582
cudaExternalSemaphoreHandleTypeNvSciSync, 
# 2586
cudaExternalSemaphoreHandleTypeKeyedMutex, 
# 2590
cudaExternalSemaphoreHandleTypeKeyedMutexKmt, 
# 2594
cudaExternalSemaphoreHandleTypeTimelineSemaphoreFd, 
# 2598
cudaExternalSemaphoreHandleTypeTimelineSemaphoreWin32
# 2599
}; 
#endif
# 2604 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
#if 0
# 2604
struct cudaExternalSemaphoreHandleDesc { 
# 2608
cudaExternalSemaphoreHandleType type; 
# 2609
union { 
# 2616
int fd; 
# 2632 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
struct { 
# 2636
void *handle; 
# 2641
const void *name; 
# 2642
} win32; 
# 2646
const void *nvSciSyncObj; 
# 2647
} handle; 
# 2651
unsigned flags; __pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)
# 2652
}; 
#endif
# 2657 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
#if 0
# 2657
struct cudaExternalSemaphoreSignalParams_v1 { 
# 2658
struct { 
# 2662
struct { 
# 2666
unsigned long long value; 
# 2667
} fence; 
# 2668
union { 
# 2673
void *fence; 
# 2674
unsigned long long reserved; 
# 2675
} nvSciSync; 
# 2679
struct { 
# 2683
unsigned long long key; 
# 2684
} keyedMutex; 
# 2685
} params; 
# 2696 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
unsigned flags; 
# 2697
}; 
#endif
# 2702 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
#if 0
# 2702
struct cudaExternalSemaphoreWaitParams_v1 { 
# 2703
struct { 
# 2707
struct { 
# 2711
unsigned long long value; 
# 2712
} fence; 
# 2713
union { 
# 2718
void *fence; 
# 2719
unsigned long long reserved; 
# 2720
} nvSciSync; 
# 2724
struct { 
# 2728
unsigned long long key; 
# 2732
unsigned timeoutMs; 
# 2733
} keyedMutex; 
# 2734
} params; 
# 2745 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
unsigned flags; 
# 2746
}; 
#endif
# 2751 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
#if 0
# 2751
struct cudaExternalSemaphoreSignalParams { 
# 2752
struct { 
# 2756
struct { 
# 2760
unsigned long long value; 
# 2761
} fence; 
# 2762
union { 
# 2767
void *fence; 
# 2768
unsigned long long reserved; 
# 2769
} nvSciSync; 
# 2773
struct { 
# 2777
unsigned long long key; 
# 2778
} keyedMutex; 
# 2779
unsigned reserved[12]; 
# 2780
} params; 
# 2791 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
unsigned flags; 
# 2792
unsigned reserved[16]; 
# 2793
}; 
#endif
# 2798 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
#if 0
# 2798
struct cudaExternalSemaphoreWaitParams { 
# 2799
struct { 
# 2803
struct { 
# 2807
unsigned long long value; 
# 2808
} fence; 
# 2809
union { 
# 2814
void *fence; 
# 2815
unsigned long long reserved; 
# 2816
} nvSciSync; 
# 2820
struct { 
# 2824
unsigned long long key; 
# 2828
unsigned timeoutMs; 
# 2829
} keyedMutex; 
# 2830
unsigned reserved[10]; 
# 2831
} params; 
# 2842 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
unsigned flags; 
# 2843
unsigned reserved[16]; 
# 2844
}; 
#endif
# 2855 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
#if 0
typedef cudaError 
# 2855
cudaError_t; 
#endif
# 2860 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
#if 0
typedef struct CUstream_st *
# 2860
cudaStream_t; 
#endif
# 2865 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
#if 0
typedef struct CUevent_st *
# 2865
cudaEvent_t; 
#endif
# 2870 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
#if 0
typedef cudaGraphicsResource *
# 2870
cudaGraphicsResource_t; 
#endif
# 2875 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
#if 0
typedef struct CUexternalMemory_st *
# 2875
cudaExternalMemory_t; 
#endif
# 2880 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
#if 0
typedef struct CUexternalSemaphore_st *
# 2880
cudaExternalSemaphore_t; 
#endif
# 2885 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
#if 0
typedef struct CUgraph_st *
# 2885
cudaGraph_t; 
#endif
# 2890 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
#if 0
typedef struct CUgraphNode_st *
# 2890
cudaGraphNode_t; 
#endif
# 2895 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
#if 0
typedef struct CUuserObject_st *
# 2895
cudaUserObject_t; 
#endif
# 2900 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
#if 0
typedef struct CUfunc_st *
# 2900
cudaFunction_t; 
#endif
# 2905 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
#if 0
typedef struct CUkern_st *
# 2905
cudaKernel_t; 
#endif
# 2910 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
#if 0
typedef struct CUmemPoolHandle_st *
# 2910
cudaMemPool_t; 
#endif
# 2915 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
#if 0
# 2915
enum cudaCGScope { 
# 2916
cudaCGScopeInvalid, 
# 2917
cudaCGScopeGrid, 
# 2918
cudaCGScopeMultiGrid
# 2919
}; 
#endif
# 2924 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
#if 0
# 2924
struct cudaLaunchParams { 
# 2926
void *func; 
# 2927
dim3 gridDim; 
# 2928
dim3 blockDim; 
# 2929
void **args; 
# 2930
size_t sharedMem; 
# 2931
cudaStream_t stream; 
# 2932
}; 
#endif
# 2937 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
#if 0
# 2937
struct cudaKernelNodeParams { 
# 2938
void *func; 
# 2939
dim3 gridDim; 
# 2940
dim3 blockDim; 
# 2941
unsigned sharedMemBytes; 
# 2942
void **kernelParams; 
# 2943
void **extra; 
# 2944
}; 
#endif
# 2949 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
#if 0
# 2949
struct cudaKernelNodeParamsV2 { 
# 2950
void *func; 
# 2952
dim3 gridDim; 
# 2953
dim3 blockDim; 
# 2959
unsigned sharedMemBytes; 
# 2960
void **kernelParams; 
# 2961
void **extra; 
# 2962
}; 
#endif
# 2967 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
#if 0
# 2967
struct cudaExternalSemaphoreSignalNodeParams { 
# 2968
cudaExternalSemaphore_t *extSemArray; 
# 2969
const cudaExternalSemaphoreSignalParams *paramsArray; 
# 2970
unsigned numExtSems; 
# 2971
}; 
#endif
# 2976 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
#if 0
# 2976
struct cudaExternalSemaphoreSignalNodeParamsV2 { 
# 2977
cudaExternalSemaphore_t *extSemArray; 
# 2978
const cudaExternalSemaphoreSignalParams *paramsArray; 
# 2979
unsigned numExtSems; 
# 2980
}; 
#endif
# 2985 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
#if 0
# 2985
struct cudaExternalSemaphoreWaitNodeParams { 
# 2986
cudaExternalSemaphore_t *extSemArray; 
# 2987
const cudaExternalSemaphoreWaitParams *paramsArray; 
# 2988
unsigned numExtSems; 
# 2989
}; 
#endif
# 2994 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
#if 0
# 2994
struct cudaExternalSemaphoreWaitNodeParamsV2 { 
# 2995
cudaExternalSemaphore_t *extSemArray; 
# 2996
const cudaExternalSemaphoreWaitParams *paramsArray; 
# 2997
unsigned numExtSems; 
# 2998
}; 
#endif
# 3003 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
#if 0
# 3003
enum cudaGraphNodeType { 
# 3004
cudaGraphNodeTypeKernel, 
# 3005
cudaGraphNodeTypeMemcpy, 
# 3006
cudaGraphNodeTypeMemset, 
# 3007
cudaGraphNodeTypeHost, 
# 3008
cudaGraphNodeTypeGraph, 
# 3009
cudaGraphNodeTypeEmpty, 
# 3010
cudaGraphNodeTypeWaitEvent, 
# 3011
cudaGraphNodeTypeEventRecord, 
# 3012
cudaGraphNodeTypeExtSemaphoreSignal, 
# 3013
cudaGraphNodeTypeExtSemaphoreWait, 
# 3014
cudaGraphNodeTypeMemAlloc, 
# 3015
cudaGraphNodeTypeMemFree, 
# 3016
cudaGraphNodeTypeCount
# 3017
}; 
#endif
# 3022 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
#if 0
# 3022
struct cudaChildGraphNodeParams { 
# 3023
cudaGraph_t graph; 
# 3025
}; 
#endif
# 3030 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
#if 0
# 3030
struct cudaEventRecordNodeParams { 
# 3031
cudaEvent_t event; 
# 3032
}; 
#endif
# 3037 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
#if 0
# 3037
struct cudaEventWaitNodeParams { 
# 3038
cudaEvent_t event; 
# 3039
}; 
#endif
# 3041 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
#if 0
# 3041
struct cudaGraphNodeParams { 
# 3042
cudaGraphNodeType type; 
# 3043
int reserved0[3]; 
# 3045
union { 
# 3046
long long reserved1[29]; 
# 3047
cudaKernelNodeParamsV2 kernel; 
# 3048
cudaMemcpyNodeParams memcpy; 
# 3049
cudaMemsetParamsV2 memset; 
# 3050
cudaHostNodeParamsV2 host; 
# 3051
cudaChildGraphNodeParams graph; 
# 3052
cudaEventWaitNodeParams eventWait; 
# 3053
cudaEventRecordNodeParams eventRecord; 
# 3054
cudaExternalSemaphoreSignalNodeParamsV2 extSemSignal; 
# 3055
cudaExternalSemaphoreWaitNodeParamsV2 extSemWait; 
# 3056
cudaMemAllocNodeParamsV2 alloc; 
# 3057
cudaMemFreeNodeParams free; __pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)
# 3058
}; 
# 3060
long long reserved2; 
# 3061
}; 
#endif
# 3066 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
typedef struct CUgraphExec_st *cudaGraphExec_t; 
# 3071
#if 0
# 3071
enum cudaGraphExecUpdateResult { 
# 3072
cudaGraphExecUpdateSuccess, 
# 3073
cudaGraphExecUpdateError, 
# 3074
cudaGraphExecUpdateErrorTopologyChanged, 
# 3075
cudaGraphExecUpdateErrorNodeTypeChanged, 
# 3076
cudaGraphExecUpdateErrorFunctionChanged, 
# 3077
cudaGraphExecUpdateErrorParametersChanged, 
# 3078
cudaGraphExecUpdateErrorNotSupported, 
# 3079
cudaGraphExecUpdateErrorUnsupportedFunctionChange, 
# 3080
cudaGraphExecUpdateErrorAttributesChanged
# 3081
}; 
#endif
# 3092 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
#if 0
typedef 
# 3086
enum cudaGraphInstantiateResult { 
# 3087
cudaGraphInstantiateSuccess, 
# 3088
cudaGraphInstantiateError, 
# 3089
cudaGraphInstantiateInvalidStructure, 
# 3090
cudaGraphInstantiateNodeOperationNotSupported, 
# 3091
cudaGraphInstantiateMultipleDevicesNotSupported
# 3092
} cudaGraphInstantiateResult; 
#endif
# 3103 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
#if 0
typedef 
# 3097
struct cudaGraphInstantiateParams_st { 
# 3099
unsigned long long flags; 
# 3100
cudaStream_t uploadStream; 
# 3101
cudaGraphNode_t errNode_out; 
# 3102
cudaGraphInstantiateResult result_out; 
# 3103
} cudaGraphInstantiateParams; 
#endif
# 3125 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
#if 0
typedef 
# 3108
struct cudaGraphExecUpdateResultInfo_st { 
# 3112
cudaGraphExecUpdateResult result; 
# 3119
cudaGraphNode_t errorNode; 
# 3124
cudaGraphNode_t errorFromNode; 
# 3125
} cudaGraphExecUpdateResultInfo; 
#endif
# 3131 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
#if 0
# 3131
enum cudaGetDriverEntryPointFlags { 
# 3132
cudaEnableDefault, 
# 3133
cudaEnableLegacyStream, 
# 3134
cudaEnablePerThreadDefaultStream
# 3135
}; 
#endif
# 3140 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
#if 0
# 3140
enum cudaDriverEntryPointQueryResult { 
# 3141
cudaDriverEntryPointSuccess, 
# 3142
cudaDriverEntryPointSymbolNotFound, 
# 3143
cudaDriverEntryPointVersionNotSufficent
# 3144
}; 
#endif
# 3149 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
#if 0
# 3149
enum cudaGraphDebugDotFlags { 
# 3150
cudaGraphDebugDotFlagsVerbose = (1 << 0), 
# 3151
cudaGraphDebugDotFlagsKernelNodeParams = (1 << 2), 
# 3152
cudaGraphDebugDotFlagsMemcpyNodeParams = (1 << 3), 
# 3153
cudaGraphDebugDotFlagsMemsetNodeParams = (1 << 4), 
# 3154
cudaGraphDebugDotFlagsHostNodeParams = (1 << 5), 
# 3155
cudaGraphDebugDotFlagsEventNodeParams = (1 << 6), 
# 3156
cudaGraphDebugDotFlagsExtSemasSignalNodeParams = (1 << 7), 
# 3157
cudaGraphDebugDotFlagsExtSemasWaitNodeParams = (1 << 8), 
# 3158
cudaGraphDebugDotFlagsKernelNodeAttributes = (1 << 9), 
# 3159
cudaGraphDebugDotFlagsHandles = (1 << 10)
# 3160
}; 
#endif
# 3165 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
#if 0
# 3165
enum cudaGraphInstantiateFlags { 
# 3166
cudaGraphInstantiateFlagAutoFreeOnLaunch = 1, 
# 3167
cudaGraphInstantiateFlagUpload, 
# 3168
cudaGraphInstantiateFlagDeviceLaunch = 4, 
# 3169
cudaGraphInstantiateFlagUseNodePriority = 8
# 3171
}; 
#endif
# 3176 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
#if 0
typedef 
# 3173
enum cudaLaunchMemSyncDomain { 
# 3174
cudaLaunchMemSyncDomainDefault, 
# 3175
cudaLaunchMemSyncDomainRemote
# 3176
} cudaLaunchMemSyncDomain; 
#endif
# 3181 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
#if 0
typedef 
# 3178
struct cudaLaunchMemSyncDomainMap_st { 
# 3179
unsigned char default_; 
# 3180
unsigned char remote; 
# 3181
} cudaLaunchMemSyncDomainMap; 
#endif
# 3227 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
#if 0
typedef 
# 3186 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
enum cudaLaunchAttributeID { 
# 3187
cudaLaunchAttributeIgnore, 
# 3188
cudaLaunchAttributeAccessPolicyWindow, 
# 3189
cudaLaunchAttributeCooperative, 
# 3190
cudaLaunchAttributeSynchronizationPolicy, 
# 3191
cudaLaunchAttributeClusterDimension, 
# 3192
cudaLaunchAttributeClusterSchedulingPolicyPreference, 
# 3193
cudaLaunchAttributeProgrammaticStreamSerialization, 
# 3204 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
cudaLaunchAttributeProgrammaticEvent, 
# 3224 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
cudaLaunchAttributePriority, 
# 3225
cudaLaunchAttributeMemSyncDomainMap, 
# 3226
cudaLaunchAttributeMemSyncDomain
# 3227
} cudaLaunchAttributeID; 
#endif
# 3252 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
#if 0
typedef 
# 3232
union cudaLaunchAttributeValue { 
# 3233
char pad[64]; 
# 3234
cudaAccessPolicyWindow accessPolicyWindow; 
# 3235
int cooperative; 
# 3236
cudaSynchronizationPolicy syncPolicy; 
# 3237
struct { 
# 3238
unsigned x; 
# 3239
unsigned y; 
# 3240
unsigned z; 
# 3241
} clusterDim; 
# 3242
cudaClusterSchedulingPolicy clusterSchedulingPolicyPreference; 
# 3243
int programmaticStreamSerializationAllowed; 
# 3244
struct { 
# 3245
cudaEvent_t event; 
# 3246
int flags; 
# 3247
int triggerAtBlockStart; 
# 3248
} programmaticEvent; 
# 3249
int priority; 
# 3250
cudaLaunchMemSyncDomainMap memSyncDomainMap; 
# 3251
cudaLaunchMemSyncDomain memSyncDomain; __pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)
# 3252
} cudaLaunchAttributeValue; 
#endif
# 3261 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
#if 0
typedef 
# 3257
struct cudaLaunchAttribute_st { 
# 3258
cudaLaunchAttributeID id; 
# 3259
char pad[(8) - sizeof(cudaLaunchAttributeID)]; 
# 3260
cudaLaunchAttributeValue val; 
# 3261
} cudaLaunchAttribute; 
#endif
# 3273 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
#if 0
typedef 
# 3266
struct cudaLaunchConfig_st { 
# 3267
dim3 gridDim; 
# 3268
dim3 blockDim; 
# 3269
size_t dynamicSmemBytes; 
# 3270
cudaStream_t stream; 
# 3271
cudaLaunchAttribute *attrs; 
# 3272
unsigned numAttrs; __pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)__pad__(volatile char:8;)
# 3273
} cudaLaunchConfig_t; 
#endif
# 3295 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_types.h"
#if 0
# 3295
enum cudaDeviceNumaConfig { 
# 3296
cudaDeviceNumaConfigNone, 
# 3297
cudaDeviceNumaConfigNumaNode
# 3298
}; 
#endif
# 84 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/surface_types.h"
#if 0
# 84
enum cudaSurfaceBoundaryMode { 
# 86
cudaBoundaryModeZero, 
# 87
cudaBoundaryModeClamp, 
# 88
cudaBoundaryModeTrap
# 89
}; 
#endif
# 94 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/surface_types.h"
#if 0
# 94
enum cudaSurfaceFormatMode { 
# 96
cudaFormatModeForced, 
# 97
cudaFormatModeAuto
# 98
}; 
#endif
# 103 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/surface_types.h"
#if 0
typedef unsigned long long 
# 103
cudaSurfaceObject_t; 
#endif
# 84 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/texture_types.h"
#if 0
# 84
enum cudaTextureAddressMode { 
# 86
cudaAddressModeWrap, 
# 87
cudaAddressModeClamp, 
# 88
cudaAddressModeMirror, 
# 89
cudaAddressModeBorder
# 90
}; 
#endif
# 95 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/texture_types.h"
#if 0
# 95
enum cudaTextureFilterMode { 
# 97
cudaFilterModePoint, 
# 98
cudaFilterModeLinear
# 99
}; 
#endif
# 104 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/texture_types.h"
#if 0
# 104
enum cudaTextureReadMode { 
# 106
cudaReadModeElementType, 
# 107
cudaReadModeNormalizedFloat
# 108
}; 
#endif
# 113 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/texture_types.h"
#if 0
# 113
struct cudaTextureDesc { 
# 118
cudaTextureAddressMode addressMode[3]; 
# 122
cudaTextureFilterMode filterMode; 
# 126
cudaTextureReadMode readMode; 
# 130
int sRGB; 
# 134
float borderColor[4]; 
# 138
int normalizedCoords; 
# 142
unsigned maxAnisotropy; 
# 146
cudaTextureFilterMode mipmapFilterMode; 
# 150
float mipmapLevelBias; 
# 154
float minMipmapLevelClamp; 
# 158
float maxMipmapLevelClamp; 
# 162
int disableTrilinearOptimization; 
# 166
int seamlessCubemap; 
# 167
}; 
#endif
# 172 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/texture_types.h"
#if 0
typedef unsigned long long 
# 172
cudaTextureObject_t; 
#endif
# 87 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/library_types.h"
typedef 
# 55
enum cudaDataType_t { 
# 57
CUDA_R_16F = 2, 
# 58
CUDA_C_16F = 6, 
# 59
CUDA_R_16BF = 14, 
# 60
CUDA_C_16BF, 
# 61
CUDA_R_32F = 0, 
# 62
CUDA_C_32F = 4, 
# 63
CUDA_R_64F = 1, 
# 64
CUDA_C_64F = 5, 
# 65
CUDA_R_4I = 16, 
# 66
CUDA_C_4I, 
# 67
CUDA_R_4U, 
# 68
CUDA_C_4U, 
# 69
CUDA_R_8I = 3, 
# 70
CUDA_C_8I = 7, 
# 71
CUDA_R_8U, 
# 72
CUDA_C_8U, 
# 73
CUDA_R_16I = 20, 
# 74
CUDA_C_16I, 
# 75
CUDA_R_16U, 
# 76
CUDA_C_16U, 
# 77
CUDA_R_32I = 10, 
# 78
CUDA_C_32I, 
# 79
CUDA_R_32U, 
# 80
CUDA_C_32U, 
# 81
CUDA_R_64I = 24, 
# 82
CUDA_C_64I, 
# 83
CUDA_R_64U, 
# 84
CUDA_C_64U, 
# 85
CUDA_R_8F_E4M3, 
# 86
CUDA_R_8F_E5M2
# 87
} cudaDataType; 
# 95
typedef 
# 90
enum libraryPropertyType_t { 
# 92
MAJOR_VERSION, 
# 93
MINOR_VERSION, 
# 94
PATCH_LEVEL
# 95
} libraryPropertyType; 
# 262 "/usr/include/x86_64-linux-gnu/c++/10/bits/c++config.h" 3
namespace std { 
# 264
typedef unsigned long size_t; 
# 265
typedef long ptrdiff_t; 
# 268
typedef __decltype((nullptr)) nullptr_t; 
# 270
}
# 284 "/usr/include/x86_64-linux-gnu/c++/10/bits/c++config.h" 3
namespace std { 
# 286
inline namespace __cxx11 __attribute((__abi_tag__("cxx11"))) { }
# 287
}
# 288
namespace __gnu_cxx { 
# 290
inline namespace __cxx11 __attribute((__abi_tag__("cxx11"))) { }
# 291
}
# 33 "/usr/include/stdlib.h" 3
extern "C" {
# 57 "/usr/include/x86_64-linux-gnu/bits/waitflags.h" 3
typedef 
# 53
enum { 
# 54
P_ALL, 
# 55
P_PID, 
# 56
P_PGID
# 57
} idtype_t; 
# 214 "/usr/include/x86_64-linux-gnu/bits/floatn-common.h" 3
typedef float _Float32; 
# 251 "/usr/include/x86_64-linux-gnu/bits/floatn-common.h" 3
typedef double _Float64; 
# 268 "/usr/include/x86_64-linux-gnu/bits/floatn-common.h" 3
typedef double _Float32x; 
# 285 "/usr/include/x86_64-linux-gnu/bits/floatn-common.h" 3
typedef long double _Float64x; 
# 62 "/usr/include/stdlib.h" 3
typedef 
# 59
struct { 
# 60
int quot; 
# 61
int rem; 
# 62
} div_t; 
# 70
typedef 
# 67
struct { 
# 68
long quot; 
# 69
long rem; 
# 70
} ldiv_t; 
# 80
__extension__ typedef 
# 77
struct { 
# 78
long long quot; 
# 79
long long rem; 
# 80
} lldiv_t; 
# 97 "/usr/include/stdlib.h" 3
extern size_t __ctype_get_mb_cur_max() throw(); 
# 101
extern double atof(const char * __nptr) throw()
# 102
 __attribute((__pure__)) __attribute((__nonnull__(1))); 
# 104
extern int atoi(const char * __nptr) throw()
# 105
 __attribute((__pure__)) __attribute((__nonnull__(1))); 
# 107
extern long atol(const char * __nptr) throw()
# 108
 __attribute((__pure__)) __attribute((__nonnull__(1))); 
# 112
__extension__ extern long long atoll(const char * __nptr) throw()
# 113
 __attribute((__pure__)) __attribute((__nonnull__(1))); 
# 117
extern double strtod(const char *__restrict__ __nptr, char **__restrict__ __endptr) throw()
# 119
 __attribute((__nonnull__(1))); 
# 123
extern float strtof(const char *__restrict__ __nptr, char **__restrict__ __endptr) throw()
# 124
 __attribute((__nonnull__(1))); 
# 126
extern long double strtold(const char *__restrict__ __nptr, char **__restrict__ __endptr) throw()
# 128
 __attribute((__nonnull__(1))); 
# 140 "/usr/include/stdlib.h" 3
extern _Float32 strtof32(const char *__restrict__ __nptr, char **__restrict__ __endptr) throw()
# 142
 __attribute((__nonnull__(1))); 
# 146
extern _Float64 strtof64(const char *__restrict__ __nptr, char **__restrict__ __endptr) throw()
# 148
 __attribute((__nonnull__(1))); 
# 158 "/usr/include/stdlib.h" 3
extern _Float32x strtof32x(const char *__restrict__ __nptr, char **__restrict__ __endptr) throw()
# 160
 __attribute((__nonnull__(1))); 
# 164
extern _Float64x strtof64x(const char *__restrict__ __nptr, char **__restrict__ __endptr) throw()
# 166
 __attribute((__nonnull__(1))); 
# 176 "/usr/include/stdlib.h" 3
extern long strtol(const char *__restrict__ __nptr, char **__restrict__ __endptr, int __base) throw()
# 178
 __attribute((__nonnull__(1))); 
# 180
extern unsigned long strtoul(const char *__restrict__ __nptr, char **__restrict__ __endptr, int __base) throw()
# 182
 __attribute((__nonnull__(1))); 
# 187
__extension__ extern long long strtoq(const char *__restrict__ __nptr, char **__restrict__ __endptr, int __base) throw()
# 189
 __attribute((__nonnull__(1))); 
# 192
__extension__ extern unsigned long long strtouq(const char *__restrict__ __nptr, char **__restrict__ __endptr, int __base) throw()
# 194
 __attribute((__nonnull__(1))); 
# 200
__extension__ extern long long strtoll(const char *__restrict__ __nptr, char **__restrict__ __endptr, int __base) throw()
# 202
 __attribute((__nonnull__(1))); 
# 205
__extension__ extern unsigned long long strtoull(const char *__restrict__ __nptr, char **__restrict__ __endptr, int __base) throw()
# 207
 __attribute((__nonnull__(1))); 
# 212
extern int strfromd(char * __dest, size_t __size, const char * __format, double __f) throw()
# 214
 __attribute((__nonnull__(3))); 
# 216
extern int strfromf(char * __dest, size_t __size, const char * __format, float __f) throw()
# 218
 __attribute((__nonnull__(3))); 
# 220
extern int strfroml(char * __dest, size_t __size, const char * __format, long double __f) throw()
# 222
 __attribute((__nonnull__(3))); 
# 232 "/usr/include/stdlib.h" 3
extern int strfromf32(char * __dest, size_t __size, const char * __format, _Float32 __f) throw()
# 234
 __attribute((__nonnull__(3))); 
# 238
extern int strfromf64(char * __dest, size_t __size, const char * __format, _Float64 __f) throw()
# 240
 __attribute((__nonnull__(3))); 
# 250 "/usr/include/stdlib.h" 3
extern int strfromf32x(char * __dest, size_t __size, const char * __format, _Float32x __f) throw()
# 252
 __attribute((__nonnull__(3))); 
# 256
extern int strfromf64x(char * __dest, size_t __size, const char * __format, _Float64x __f) throw()
# 258
 __attribute((__nonnull__(3))); 
# 28 "/usr/include/x86_64-linux-gnu/bits/types/__locale_t.h" 3
struct __locale_struct { 
# 31
struct __locale_data *__locales[13]; 
# 34
const unsigned short *__ctype_b; 
# 35
const int *__ctype_tolower; 
# 36
const int *__ctype_toupper; 
# 39
const char *__names[13]; 
# 40
}; 
# 42
typedef __locale_struct *__locale_t; 
# 24 "/usr/include/x86_64-linux-gnu/bits/types/locale_t.h" 3
typedef __locale_t locale_t; 
# 274 "/usr/include/stdlib.h" 3
extern long strtol_l(const char *__restrict__ __nptr, char **__restrict__ __endptr, int __base, locale_t __loc) throw()
# 276
 __attribute((__nonnull__(1, 4))); 
# 278
extern unsigned long strtoul_l(const char *__restrict__ __nptr, char **__restrict__ __endptr, int __base, locale_t __loc) throw()
# 281
 __attribute((__nonnull__(1, 4))); 
# 284
__extension__ extern long long strtoll_l(const char *__restrict__ __nptr, char **__restrict__ __endptr, int __base, locale_t __loc) throw()
# 287
 __attribute((__nonnull__(1, 4))); 
# 290
__extension__ extern unsigned long long strtoull_l(const char *__restrict__ __nptr, char **__restrict__ __endptr, int __base, locale_t __loc) throw()
# 293
 __attribute((__nonnull__(1, 4))); 
# 295
extern double strtod_l(const char *__restrict__ __nptr, char **__restrict__ __endptr, locale_t __loc) throw()
# 297
 __attribute((__nonnull__(1, 3))); 
# 299
extern float strtof_l(const char *__restrict__ __nptr, char **__restrict__ __endptr, locale_t __loc) throw()
# 301
 __attribute((__nonnull__(1, 3))); 
# 303
extern long double strtold_l(const char *__restrict__ __nptr, char **__restrict__ __endptr, locale_t __loc) throw()
# 306
 __attribute((__nonnull__(1, 3))); 
# 316 "/usr/include/stdlib.h" 3
extern _Float32 strtof32_l(const char *__restrict__ __nptr, char **__restrict__ __endptr, locale_t __loc) throw()
# 319
 __attribute((__nonnull__(1, 3))); 
# 323
extern _Float64 strtof64_l(const char *__restrict__ __nptr, char **__restrict__ __endptr, locale_t __loc) throw()
# 326
 __attribute((__nonnull__(1, 3))); 
# 337 "/usr/include/stdlib.h" 3
extern _Float32x strtof32x_l(const char *__restrict__ __nptr, char **__restrict__ __endptr, locale_t __loc) throw()
# 340
 __attribute((__nonnull__(1, 3))); 
# 344
extern _Float64x strtof64x_l(const char *__restrict__ __nptr, char **__restrict__ __endptr, locale_t __loc) throw()
# 347
 __attribute((__nonnull__(1, 3))); 
# 385 "/usr/include/stdlib.h" 3
extern char *l64a(long __n) throw(); 
# 388
extern long a64l(const char * __s) throw()
# 389
 __attribute((__pure__)) __attribute((__nonnull__(1))); 
# 27 "/usr/include/x86_64-linux-gnu/sys/types.h" 3
extern "C" {
# 31 "/usr/include/x86_64-linux-gnu/bits/types.h" 3
typedef unsigned char __u_char; 
# 32
typedef unsigned short __u_short; 
# 33
typedef unsigned __u_int; 
# 34
typedef unsigned long __u_long; 
# 37
typedef signed char __int8_t; 
# 38
typedef unsigned char __uint8_t; 
# 39
typedef signed short __int16_t; 
# 40
typedef unsigned short __uint16_t; 
# 41
typedef signed int __int32_t; 
# 42
typedef unsigned __uint32_t; 
# 44
typedef signed long __int64_t; 
# 45
typedef unsigned long __uint64_t; 
# 52
typedef __int8_t __int_least8_t; 
# 53
typedef __uint8_t __uint_least8_t; 
# 54
typedef __int16_t __int_least16_t; 
# 55
typedef __uint16_t __uint_least16_t; 
# 56
typedef __int32_t __int_least32_t; 
# 57
typedef __uint32_t __uint_least32_t; 
# 58
typedef __int64_t __int_least64_t; 
# 59
typedef __uint64_t __uint_least64_t; 
# 63
typedef long __quad_t; 
# 64
typedef unsigned long __u_quad_t; 
# 72
typedef long __intmax_t; 
# 73
typedef unsigned long __uintmax_t; 
# 145 "/usr/include/x86_64-linux-gnu/bits/types.h" 3
typedef unsigned long __dev_t; 
# 146
typedef unsigned __uid_t; 
# 147
typedef unsigned __gid_t; 
# 148
typedef unsigned long __ino_t; 
# 149
typedef unsigned long __ino64_t; 
# 150
typedef unsigned __mode_t; 
# 151
typedef unsigned long __nlink_t; 
# 152
typedef long __off_t; 
# 153
typedef long __off64_t; 
# 154
typedef int __pid_t; 
# 155
typedef struct { int __val[2]; } __fsid_t; 
# 156
typedef long __clock_t; 
# 157
typedef unsigned long __rlim_t; 
# 158
typedef unsigned long __rlim64_t; 
# 159
typedef unsigned __id_t; 
# 160
typedef long __time_t; 
# 161
typedef unsigned __useconds_t; 
# 162
typedef long __suseconds_t; 
# 164
typedef int __daddr_t; 
# 165
typedef int __key_t; 
# 168
typedef int __clockid_t; 
# 171
typedef void *__timer_t; 
# 174
typedef long __blksize_t; 
# 179
typedef long __blkcnt_t; 
# 180
typedef long __blkcnt64_t; 
# 183
typedef unsigned long __fsblkcnt_t; 
# 184
typedef unsigned long __fsblkcnt64_t; 
# 187
typedef unsigned long __fsfilcnt_t; 
# 188
typedef unsigned long __fsfilcnt64_t; 
# 191
typedef long __fsword_t; 
# 193
typedef long __ssize_t; 
# 196
typedef long __syscall_slong_t; 
# 198
typedef unsigned long __syscall_ulong_t; 
# 202
typedef __off64_t __loff_t; 
# 203
typedef char *__caddr_t; 
# 206
typedef long __intptr_t; 
# 209
typedef unsigned __socklen_t; 
# 214
typedef int __sig_atomic_t; 
# 33 "/usr/include/x86_64-linux-gnu/sys/types.h" 3
typedef __u_char u_char; 
# 34
typedef __u_short u_short; 
# 35
typedef __u_int u_int; 
# 36
typedef __u_long u_long; 
# 37
typedef __quad_t quad_t; 
# 38
typedef __u_quad_t u_quad_t; 
# 39
typedef __fsid_t fsid_t; 
# 42
typedef __loff_t loff_t; 
# 47
typedef __ino_t ino_t; 
# 54
typedef __ino64_t ino64_t; 
# 59
typedef __dev_t dev_t; 
# 64
typedef __gid_t gid_t; 
# 69
typedef __mode_t mode_t; 
# 74
typedef __nlink_t nlink_t; 
# 79
typedef __uid_t uid_t; 
# 85
typedef __off_t off_t; 
# 92
typedef __off64_t off64_t; 
# 97
typedef __pid_t pid_t; 
# 103
typedef __id_t id_t; 
# 108
typedef __ssize_t ssize_t; 
# 114
typedef __daddr_t daddr_t; 
# 115
typedef __caddr_t caddr_t; 
# 121
typedef __key_t key_t; 
# 7 "/usr/include/x86_64-linux-gnu/bits/types/clock_t.h" 3
typedef __clock_t clock_t; 
# 7 "/usr/include/x86_64-linux-gnu/bits/types/clockid_t.h" 3
typedef __clockid_t clockid_t; 
# 7 "/usr/include/x86_64-linux-gnu/bits/types/time_t.h" 3
typedef __time_t time_t; 
# 7 "/usr/include/x86_64-linux-gnu/bits/types/timer_t.h" 3
typedef __timer_t timer_t; 
# 134 "/usr/include/x86_64-linux-gnu/sys/types.h" 3
typedef __useconds_t useconds_t; 
# 138
typedef __suseconds_t suseconds_t; 
# 148 "/usr/include/x86_64-linux-gnu/sys/types.h" 3
typedef unsigned long ulong; 
# 149
typedef unsigned short ushort; 
# 150
typedef unsigned uint; 
# 24 "/usr/include/x86_64-linux-gnu/bits/stdint-intn.h" 3
typedef __int8_t int8_t; 
# 25
typedef __int16_t int16_t; 
# 26
typedef __int32_t int32_t; 
# 27
typedef __int64_t int64_t; 
# 158 "/usr/include/x86_64-linux-gnu/sys/types.h" 3
typedef __uint8_t u_int8_t; 
# 159
typedef __uint16_t u_int16_t; 
# 160
typedef __uint32_t u_int32_t; 
# 161
typedef __uint64_t u_int64_t; 
# 164
typedef long register_t __attribute((__mode__(__word__))); 
# 34 "/usr/include/x86_64-linux-gnu/bits/byteswap.h" 3
static inline __uint16_t __bswap_16(__uint16_t __bsx) 
# 35
{ 
# 37
return __builtin_bswap16(__bsx); 
# 41
} 
# 49
static inline __uint32_t __bswap_32(__uint32_t __bsx) 
# 50
{ 
# 52
return __builtin_bswap32(__bsx); 
# 56
} 
# 70 "/usr/include/x86_64-linux-gnu/bits/byteswap.h" 3
__extension__ static inline __uint64_t __bswap_64(__uint64_t __bsx) 
# 71
{ 
# 73
return __builtin_bswap64(__bsx); 
# 77
} 
# 33 "/usr/include/x86_64-linux-gnu/bits/uintn-identity.h" 3
static inline __uint16_t __uint16_identity(__uint16_t __x) 
# 34
{ 
# 35
return __x; 
# 36
} 
# 39
static inline __uint32_t __uint32_identity(__uint32_t __x) 
# 40
{ 
# 41
return __x; 
# 42
} 
# 45
static inline __uint64_t __uint64_identity(__uint64_t __x) 
# 46
{ 
# 47
return __x; 
# 48
} 
# 8 "/usr/include/x86_64-linux-gnu/bits/types/__sigset_t.h" 3
typedef 
# 6
struct { 
# 7
unsigned long __val[(1024) / ((8) * sizeof(unsigned long))]; 
# 8
} __sigset_t; 
# 7 "/usr/include/x86_64-linux-gnu/bits/types/sigset_t.h" 3
typedef __sigset_t sigset_t; 
# 8 "/usr/include/x86_64-linux-gnu/bits/types/struct_timeval.h" 3
struct timeval { 
# 10
__time_t tv_sec; 
# 11
__suseconds_t tv_usec; 
# 12
}; 
# 10 "/usr/include/x86_64-linux-gnu/bits/types/struct_timespec.h" 3
struct timespec { 
# 12
__time_t tv_sec; 
# 16
__syscall_slong_t tv_nsec; 
# 26 "/usr/include/x86_64-linux-gnu/bits/types/struct_timespec.h" 3
}; 
# 49 "/usr/include/x86_64-linux-gnu/sys/select.h" 3
typedef long __fd_mask; 
# 70 "/usr/include/x86_64-linux-gnu/sys/select.h" 3
typedef 
# 60
struct { 
# 64
__fd_mask fds_bits[1024 / (8 * ((int)sizeof(__fd_mask)))]; 
# 70
} fd_set; 
# 77
typedef __fd_mask fd_mask; 
# 91 "/usr/include/x86_64-linux-gnu/sys/select.h" 3
extern "C" {
# 101 "/usr/include/x86_64-linux-gnu/sys/select.h" 3
extern int select(int __nfds, fd_set *__restrict__ __readfds, fd_set *__restrict__ __writefds, fd_set *__restrict__ __exceptfds, timeval *__restrict__ __timeout); 
# 113 "/usr/include/x86_64-linux-gnu/sys/select.h" 3
extern int pselect(int __nfds, fd_set *__restrict__ __readfds, fd_set *__restrict__ __writefds, fd_set *__restrict__ __exceptfds, const timespec *__restrict__ __timeout, const __sigset_t *__restrict__ __sigmask); 
# 126 "/usr/include/x86_64-linux-gnu/sys/select.h" 3
}
# 185 "/usr/include/x86_64-linux-gnu/sys/types.h" 3
typedef __blksize_t blksize_t; 
# 192
typedef __blkcnt_t blkcnt_t; 
# 196
typedef __fsblkcnt_t fsblkcnt_t; 
# 200
typedef __fsfilcnt_t fsfilcnt_t; 
# 219 "/usr/include/x86_64-linux-gnu/sys/types.h" 3
typedef __blkcnt64_t blkcnt64_t; 
# 220
typedef __fsblkcnt64_t fsblkcnt64_t; 
# 221
typedef __fsfilcnt64_t fsfilcnt64_t; 
# 53 "/usr/include/x86_64-linux-gnu/bits/thread-shared-types.h" 3
typedef 
# 49
struct __pthread_internal_list { 
# 51
__pthread_internal_list *__prev; 
# 52
__pthread_internal_list *__next; 
# 53
} __pthread_list_t; 
# 58
typedef 
# 55
struct __pthread_internal_slist { 
# 57
__pthread_internal_slist *__next; 
# 58
} __pthread_slist_t; 
# 22 "/usr/include/x86_64-linux-gnu/bits/struct_mutex.h" 3
struct __pthread_mutex_s { 
# 24
int __lock; 
# 25
unsigned __count; 
# 26
int __owner; 
# 28
unsigned __nusers; 
# 32
int __kind; 
# 34
short __spins; 
# 35
short __elision; 
# 36
__pthread_list_t __list; 
# 53 "/usr/include/x86_64-linux-gnu/bits/struct_mutex.h" 3
}; 
# 23 "/usr/include/x86_64-linux-gnu/bits/struct_rwlock.h" 3
struct __pthread_rwlock_arch_t { 
# 25
unsigned __readers; 
# 26
unsigned __writers; 
# 27
unsigned __wrphase_futex; 
# 28
unsigned __writers_futex; 
# 29
unsigned __pad3; 
# 30
unsigned __pad4; 
# 32
int __cur_writer; 
# 33
int __shared; 
# 34
signed char __rwelision; 
# 39
unsigned char __pad1[7]; 
# 42
unsigned long __pad2; 
# 45
unsigned __flags; 
# 55 "/usr/include/x86_64-linux-gnu/bits/struct_rwlock.h" 3
}; 
# 92 "/usr/include/x86_64-linux-gnu/bits/thread-shared-types.h" 3
struct __pthread_cond_s { 
# 95
union { 
# 96
__extension__ unsigned long long __wseq; 
# 98
struct { 
# 99
unsigned __low; 
# 100
unsigned __high; 
# 101
} __wseq32; 
# 102
}; 
# 104
union { 
# 105
__extension__ unsigned long long __g1_start; 
# 107
struct { 
# 108
unsigned __low; 
# 109
unsigned __high; 
# 110
} __g1_start32; 
# 111
}; 
# 112
unsigned __g_refs[2]; 
# 113
unsigned __g_size[2]; 
# 114
unsigned __g1_orig_size; 
# 115
unsigned __wrefs; 
# 116
unsigned __g_signals[2]; 
# 117
}; 
# 27 "/usr/include/x86_64-linux-gnu/bits/pthreadtypes.h" 3
typedef unsigned long pthread_t; 
# 36
typedef 
# 33
union { 
# 34
char __size[4]; 
# 35
int __align; 
# 36
} pthread_mutexattr_t; 
# 45
typedef 
# 42
union { 
# 43
char __size[4]; 
# 44
int __align; 
# 45
} pthread_condattr_t; 
# 49
typedef unsigned pthread_key_t; 
# 53
typedef int pthread_once_t; 
# 56
union pthread_attr_t { 
# 58
char __size[56]; 
# 59
long __align; 
# 60
}; 
# 62
typedef pthread_attr_t pthread_attr_t; 
# 72
typedef 
# 68
union { 
# 69
__pthread_mutex_s __data; 
# 70
char __size[40]; 
# 71
long __align; 
# 72
} pthread_mutex_t; 
# 80
typedef 
# 76
union { 
# 77
__pthread_cond_s __data; 
# 78
char __size[48]; 
# 79
__extension__ long long __align; 
# 80
} pthread_cond_t; 
# 91
typedef 
# 87
union { 
# 88
__pthread_rwlock_arch_t __data; 
# 89
char __size[56]; 
# 90
long __align; 
# 91
} pthread_rwlock_t; 
# 97
typedef 
# 94
union { 
# 95
char __size[8]; 
# 96
long __align; 
# 97
} pthread_rwlockattr_t; 
# 103
typedef volatile int pthread_spinlock_t; 
# 112
typedef 
# 109
union { 
# 110
char __size[32]; 
# 111
long __align; 
# 112
} pthread_barrier_t; 
# 118
typedef 
# 115
union { 
# 116
char __size[4]; 
# 117
int __align; 
# 118
} pthread_barrierattr_t; 
# 230 "/usr/include/x86_64-linux-gnu/sys/types.h" 3
}
# 401 "/usr/include/stdlib.h" 3
extern long random() throw(); 
# 404
extern void srandom(unsigned __seed) throw(); 
# 410
extern char *initstate(unsigned __seed, char * __statebuf, size_t __statelen) throw()
# 411
 __attribute((__nonnull__(2))); 
# 415
extern char *setstate(char * __statebuf) throw() __attribute((__nonnull__(1))); 
# 423
struct random_data { 
# 425
int32_t *fptr; 
# 426
int32_t *rptr; 
# 427
int32_t *state; 
# 428
int rand_type; 
# 429
int rand_deg; 
# 430
int rand_sep; 
# 431
int32_t *end_ptr; 
# 432
}; 
# 434
extern int random_r(random_data *__restrict__ __buf, int32_t *__restrict__ __result) throw()
# 435
 __attribute((__nonnull__(1, 2))); 
# 437
extern int srandom_r(unsigned __seed, random_data * __buf) throw()
# 438
 __attribute((__nonnull__(2))); 
# 440
extern int initstate_r(unsigned __seed, char *__restrict__ __statebuf, size_t __statelen, random_data *__restrict__ __buf) throw()
# 443
 __attribute((__nonnull__(2, 4))); 
# 445
extern int setstate_r(char *__restrict__ __statebuf, random_data *__restrict__ __buf) throw()
# 447
 __attribute((__nonnull__(1, 2))); 
# 453
extern int rand() throw(); 
# 455
extern void srand(unsigned __seed) throw(); 
# 459
extern int rand_r(unsigned * __seed) throw(); 
# 467
extern double drand48() throw(); 
# 468
extern double erand48(unsigned short  __xsubi[3]) throw() __attribute((__nonnull__(1))); 
# 471
extern long lrand48() throw(); 
# 472
extern long nrand48(unsigned short  __xsubi[3]) throw()
# 473
 __attribute((__nonnull__(1))); 
# 476
extern long mrand48() throw(); 
# 477
extern long jrand48(unsigned short  __xsubi[3]) throw()
# 478
 __attribute((__nonnull__(1))); 
# 481
extern void srand48(long __seedval) throw(); 
# 482
extern unsigned short *seed48(unsigned short  __seed16v[3]) throw()
# 483
 __attribute((__nonnull__(1))); 
# 484
extern void lcong48(unsigned short  __param[7]) throw() __attribute((__nonnull__(1))); 
# 490
struct drand48_data { 
# 492
unsigned short __x[3]; 
# 493
unsigned short __old_x[3]; 
# 494
unsigned short __c; 
# 495
unsigned short __init; 
# 496
__extension__ unsigned long long __a; 
# 498
}; 
# 501
extern int drand48_r(drand48_data *__restrict__ __buffer, double *__restrict__ __result) throw()
# 502
 __attribute((__nonnull__(1, 2))); 
# 503
extern int erand48_r(unsigned short  __xsubi[3], drand48_data *__restrict__ __buffer, double *__restrict__ __result) throw()
# 505
 __attribute((__nonnull__(1, 2))); 
# 508
extern int lrand48_r(drand48_data *__restrict__ __buffer, long *__restrict__ __result) throw()
# 510
 __attribute((__nonnull__(1, 2))); 
# 511
extern int nrand48_r(unsigned short  __xsubi[3], drand48_data *__restrict__ __buffer, long *__restrict__ __result) throw()
# 514
 __attribute((__nonnull__(1, 2))); 
# 517
extern int mrand48_r(drand48_data *__restrict__ __buffer, long *__restrict__ __result) throw()
# 519
 __attribute((__nonnull__(1, 2))); 
# 520
extern int jrand48_r(unsigned short  __xsubi[3], drand48_data *__restrict__ __buffer, long *__restrict__ __result) throw()
# 523
 __attribute((__nonnull__(1, 2))); 
# 526
extern int srand48_r(long __seedval, drand48_data * __buffer) throw()
# 527
 __attribute((__nonnull__(2))); 
# 529
extern int seed48_r(unsigned short  __seed16v[3], drand48_data * __buffer) throw()
# 530
 __attribute((__nonnull__(1, 2))); 
# 532
extern int lcong48_r(unsigned short  __param[7], drand48_data * __buffer) throw()
# 534
 __attribute((__nonnull__(1, 2))); 
# 539
extern void *malloc(size_t __size) throw() __attribute((__malloc__))
# 540
 __attribute((__alloc_size__(1))); 
# 542
extern void *calloc(size_t __nmemb, size_t __size) throw()
# 543
 __attribute((__malloc__)) __attribute((__alloc_size__(1, 2))); 
# 550
extern void *realloc(void * __ptr, size_t __size) throw()
# 551
 __attribute((__warn_unused_result__)) __attribute((__alloc_size__(2))); 
# 559
extern void *reallocarray(void * __ptr, size_t __nmemb, size_t __size) throw()
# 560
 __attribute((__warn_unused_result__))
# 561
 __attribute((__alloc_size__(2, 3))); 
# 565
extern void free(void * __ptr) throw(); 
# 26 "/usr/include/alloca.h" 3
extern "C" {
# 32
extern void *alloca(size_t __size) throw(); 
# 38
}
# 574 "/usr/include/stdlib.h" 3
extern void *valloc(size_t __size) throw() __attribute((__malloc__))
# 575
 __attribute((__alloc_size__(1))); 
# 580
extern int posix_memalign(void ** __memptr, size_t __alignment, size_t __size) throw()
# 581
 __attribute((__nonnull__(1))); 
# 586
extern void *aligned_alloc(size_t __alignment, size_t __size) throw()
# 587
 __attribute((__malloc__)) __attribute((__alloc_size__(2))); 
# 591
extern void abort() throw() __attribute((__noreturn__)); 
# 595
extern int atexit(void (* __func)(void)) throw() __attribute((__nonnull__(1))); 
# 600
extern "C++" int at_quick_exit(void (* __func)(void)) throw() __asm__("at_quick_exit")
# 601
 __attribute((__nonnull__(1))); 
# 610 "/usr/include/stdlib.h" 3
extern int on_exit(void (* __func)(int __status, void * __arg), void * __arg) throw()
# 611
 __attribute((__nonnull__(1))); 
# 617
extern void exit(int __status) throw() __attribute((__noreturn__)); 
# 623
extern void quick_exit(int __status) throw() __attribute((__noreturn__)); 
# 629
extern void _Exit(int __status) throw() __attribute((__noreturn__)); 
# 634
extern char *getenv(const char * __name) throw() __attribute((__nonnull__(1))); 
# 639
extern char *secure_getenv(const char * __name) throw()
# 640
 __attribute((__nonnull__(1))); 
# 647
extern int putenv(char * __string) throw() __attribute((__nonnull__(1))); 
# 653
extern int setenv(const char * __name, const char * __value, int __replace) throw()
# 654
 __attribute((__nonnull__(2))); 
# 657
extern int unsetenv(const char * __name) throw() __attribute((__nonnull__(1))); 
# 664
extern int clearenv() throw(); 
# 675 "/usr/include/stdlib.h" 3
extern char *mktemp(char * __template) throw() __attribute((__nonnull__(1))); 
# 688 "/usr/include/stdlib.h" 3
extern int mkstemp(char * __template) __attribute((__nonnull__(1))); 
# 698 "/usr/include/stdlib.h" 3
extern int mkstemp64(char * __template) __attribute((__nonnull__(1))); 
# 710 "/usr/include/stdlib.h" 3
extern int mkstemps(char * __template, int __suffixlen) __attribute((__nonnull__(1))); 
# 720 "/usr/include/stdlib.h" 3
extern int mkstemps64(char * __template, int __suffixlen)
# 721
 __attribute((__nonnull__(1))); 
# 731 "/usr/include/stdlib.h" 3
extern char *mkdtemp(char * __template) throw() __attribute((__nonnull__(1))); 
# 742 "/usr/include/stdlib.h" 3
extern int mkostemp(char * __template, int __flags) __attribute((__nonnull__(1))); 
# 752 "/usr/include/stdlib.h" 3
extern int mkostemp64(char * __template, int __flags) __attribute((__nonnull__(1))); 
# 762 "/usr/include/stdlib.h" 3
extern int mkostemps(char * __template, int __suffixlen, int __flags)
# 763
 __attribute((__nonnull__(1))); 
# 774 "/usr/include/stdlib.h" 3
extern int mkostemps64(char * __template, int __suffixlen, int __flags)
# 775
 __attribute((__nonnull__(1))); 
# 784 "/usr/include/stdlib.h" 3
extern int system(const char * __command); 
# 790
extern char *canonicalize_file_name(const char * __name) throw()
# 791
 __attribute((__nonnull__(1))); 
# 800 "/usr/include/stdlib.h" 3
extern char *realpath(const char *__restrict__ __name, char *__restrict__ __resolved) throw(); 
# 808
typedef int (*__compar_fn_t)(const void *, const void *); 
# 811
typedef __compar_fn_t comparison_fn_t; 
# 815
typedef int (*__compar_d_fn_t)(const void *, const void *, void *); 
# 820
extern void *bsearch(const void * __key, const void * __base, size_t __nmemb, size_t __size, __compar_fn_t __compar)
# 822
 __attribute((__nonnull__(1, 2, 5))); 
# 830
extern void qsort(void * __base, size_t __nmemb, size_t __size, __compar_fn_t __compar)
# 831
 __attribute((__nonnull__(1, 4))); 
# 833
extern void qsort_r(void * __base, size_t __nmemb, size_t __size, __compar_d_fn_t __compar, void * __arg)
# 835
 __attribute((__nonnull__(1, 4))); 
# 840
extern int abs(int __x) throw() __attribute((const)); 
# 841
extern long labs(long __x) throw() __attribute((const)); 
# 844
__extension__ extern long long llabs(long long __x) throw()
# 845
 __attribute((const)); 
# 852
extern div_t div(int __numer, int __denom) throw()
# 853
 __attribute((const)); 
# 854
extern ldiv_t ldiv(long __numer, long __denom) throw()
# 855
 __attribute((const)); 
# 858
__extension__ extern lldiv_t lldiv(long long __numer, long long __denom) throw()
# 860
 __attribute((const)); 
# 872 "/usr/include/stdlib.h" 3
extern char *ecvt(double __value, int __ndigit, int *__restrict__ __decpt, int *__restrict__ __sign) throw()
# 873
 __attribute((__nonnull__(3, 4))); 
# 878
extern char *fcvt(double __value, int __ndigit, int *__restrict__ __decpt, int *__restrict__ __sign) throw()
# 879
 __attribute((__nonnull__(3, 4))); 
# 884
extern char *gcvt(double __value, int __ndigit, char * __buf) throw()
# 885
 __attribute((__nonnull__(3))); 
# 890
extern char *qecvt(long double __value, int __ndigit, int *__restrict__ __decpt, int *__restrict__ __sign) throw()
# 892
 __attribute((__nonnull__(3, 4))); 
# 893
extern char *qfcvt(long double __value, int __ndigit, int *__restrict__ __decpt, int *__restrict__ __sign) throw()
# 895
 __attribute((__nonnull__(3, 4))); 
# 896
extern char *qgcvt(long double __value, int __ndigit, char * __buf) throw()
# 897
 __attribute((__nonnull__(3))); 
# 902
extern int ecvt_r(double __value, int __ndigit, int *__restrict__ __decpt, int *__restrict__ __sign, char *__restrict__ __buf, size_t __len) throw()
# 904
 __attribute((__nonnull__(3, 4, 5))); 
# 905
extern int fcvt_r(double __value, int __ndigit, int *__restrict__ __decpt, int *__restrict__ __sign, char *__restrict__ __buf, size_t __len) throw()
# 907
 __attribute((__nonnull__(3, 4, 5))); 
# 909
extern int qecvt_r(long double __value, int __ndigit, int *__restrict__ __decpt, int *__restrict__ __sign, char *__restrict__ __buf, size_t __len) throw()
# 912
 __attribute((__nonnull__(3, 4, 5))); 
# 913
extern int qfcvt_r(long double __value, int __ndigit, int *__restrict__ __decpt, int *__restrict__ __sign, char *__restrict__ __buf, size_t __len) throw()
# 916
 __attribute((__nonnull__(3, 4, 5))); 
# 922
extern int mblen(const char * __s, size_t __n) throw(); 
# 925
extern int mbtowc(wchar_t *__restrict__ __pwc, const char *__restrict__ __s, size_t __n) throw(); 
# 929
extern int wctomb(char * __s, wchar_t __wchar) throw(); 
# 933
extern size_t mbstowcs(wchar_t *__restrict__ __pwcs, const char *__restrict__ __s, size_t __n) throw(); 
# 936
extern size_t wcstombs(char *__restrict__ __s, const wchar_t *__restrict__ __pwcs, size_t __n) throw(); 
# 946
extern int rpmatch(const char * __response) throw() __attribute((__nonnull__(1))); 
# 957 "/usr/include/stdlib.h" 3
extern int getsubopt(char **__restrict__ __optionp, char *const *__restrict__ __tokens, char **__restrict__ __valuep) throw()
# 960
 __attribute((__nonnull__(1, 2, 3))); 
# 968
extern int posix_openpt(int __oflag); 
# 976
extern int grantpt(int __fd) throw(); 
# 980
extern int unlockpt(int __fd) throw(); 
# 985
extern char *ptsname(int __fd) throw(); 
# 992
extern int ptsname_r(int __fd, char * __buf, size_t __buflen) throw()
# 993
 __attribute((__nonnull__(2))); 
# 996
extern int getpt(); 
# 1003
extern int getloadavg(double  __loadavg[], int __nelem) throw()
# 1004
 __attribute((__nonnull__(1))); 
# 1023 "/usr/include/stdlib.h" 3
}
# 46 "/usr/include/c++/10/bits/std_abs.h" 3
extern "C++" {
# 48
namespace std __attribute((__visibility__("default"))) { 
# 52
using ::abs;
# 56
inline long abs(long __i) { return __builtin_labs(__i); } 
# 61
inline long long abs(long long __x) { return __builtin_llabs(__x); } 
# 71 "/usr/include/c++/10/bits/std_abs.h" 3
constexpr double abs(double __x) 
# 72
{ return __builtin_fabs(__x); } 
# 75
constexpr float abs(float __x) 
# 76
{ return __builtin_fabsf(__x); } 
# 79
constexpr long double abs(long double __x) 
# 80
{ return __builtin_fabsl(__x); } 
# 85
constexpr __int128 abs(__int128 __x) { return (__x >= (0)) ? __x : (-__x); } 
# 108 "/usr/include/c++/10/bits/std_abs.h" 3
}
# 109
}
# 121 "/usr/include/c++/10/cstdlib" 3
extern "C++" {
# 123
namespace std __attribute((__visibility__("default"))) { 
# 127
using ::div_t;
# 128
using ::ldiv_t;
# 130
using ::abort;
# 134
using ::atexit;
# 137
using ::at_quick_exit;
# 140
using ::atof;
# 141
using ::atoi;
# 142
using ::atol;
# 143
using ::bsearch;
# 144
using ::calloc;
# 145
using ::div;
# 146
using ::exit;
# 147
using ::free;
# 148
using ::getenv;
# 149
using ::labs;
# 150
using ::ldiv;
# 151
using ::malloc;
# 153
using ::mblen;
# 154
using ::mbstowcs;
# 155
using ::mbtowc;
# 157
using ::qsort;
# 160
using ::quick_exit;
# 163
using ::rand;
# 164
using ::realloc;
# 165
using ::srand;
# 166
using ::strtod;
# 167
using ::strtol;
# 168
using ::strtoul;
# 169
using ::system;
# 171
using ::wcstombs;
# 172
using ::wctomb;
# 177
inline ldiv_t div(long __i, long __j) { return ldiv(__i, __j); } 
# 182
}
# 195 "/usr/include/c++/10/cstdlib" 3
namespace __gnu_cxx __attribute((__visibility__("default"))) { 
# 200
using ::lldiv_t;
# 206
using ::_Exit;
# 210
using ::llabs;
# 213
inline lldiv_t div(long long __n, long long __d) 
# 214
{ lldiv_t __q; (__q.quot) = (__n / __d); (__q.rem) = (__n % __d); return __q; } 
# 216
using ::lldiv;
# 227 "/usr/include/c++/10/cstdlib" 3
using ::atoll;
# 228
using ::strtoll;
# 229
using ::strtoull;
# 231
using ::strtof;
# 232
using ::strtold;
# 235
}
# 237
namespace std { 
# 240
using __gnu_cxx::lldiv_t;
# 242
using __gnu_cxx::_Exit;
# 244
using __gnu_cxx::llabs;
# 245
using __gnu_cxx::div;
# 246
using __gnu_cxx::lldiv;
# 248
using __gnu_cxx::atoll;
# 249
using __gnu_cxx::strtof;
# 250
using __gnu_cxx::strtoll;
# 251
using __gnu_cxx::strtoull;
# 252
using __gnu_cxx::strtold;
# 253
}
# 257
}
# 38 "/usr/include/c++/10/stdlib.h" 3
using std::abort;
# 39
using std::atexit;
# 40
using std::exit;
# 43
using std::at_quick_exit;
# 46
using std::quick_exit;
# 54
using std::abs;
# 55
using std::atof;
# 56
using std::atoi;
# 57
using std::atol;
# 58
using std::bsearch;
# 59
using std::calloc;
# 60
using std::div;
# 61
using std::free;
# 62
using std::getenv;
# 63
using std::labs;
# 64
using std::ldiv;
# 65
using std::malloc;
# 67
using std::mblen;
# 68
using std::mbstowcs;
# 69
using std::mbtowc;
# 71
using std::qsort;
# 72
using std::rand;
# 73
using std::realloc;
# 74
using std::srand;
# 75
using std::strtod;
# 76
using std::strtol;
# 77
using std::strtoul;
# 78
using std::system;
# 80
using std::wcstombs;
# 81
using std::wctomb;
# 180 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_device_runtime_api.h"
extern "C" {
# 187
__attribute__((unused)) extern cudaError_t __cudaDeviceSynchronizeDeprecationAvoidance(); 
# 236 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_device_runtime_api.h"
__attribute__((unused)) extern cudaError_t __cudaCDP2DeviceGetAttribute(int * value, cudaDeviceAttr attr, int device); 
# 237
__attribute__((unused)) extern cudaError_t __cudaCDP2DeviceGetLimit(size_t * pValue, cudaLimit limit); 
# 238
__attribute__((unused)) extern cudaError_t __cudaCDP2DeviceGetCacheConfig(cudaFuncCache * pCacheConfig); 
# 239
__attribute__((unused)) extern cudaError_t __cudaCDP2DeviceGetSharedMemConfig(cudaSharedMemConfig * pConfig); 
# 240
__attribute__((unused)) extern cudaError_t __cudaCDP2GetLastError(); 
# 241
__attribute__((unused)) extern cudaError_t __cudaCDP2PeekAtLastError(); 
# 242
__attribute__((unused)) extern const char *__cudaCDP2GetErrorString(cudaError_t error); 
# 243
__attribute__((unused)) extern const char *__cudaCDP2GetErrorName(cudaError_t error); 
# 244
__attribute__((unused)) extern cudaError_t __cudaCDP2GetDeviceCount(int * count); 
# 245
__attribute__((unused)) extern cudaError_t __cudaCDP2GetDevice(int * device); 
# 246
__attribute__((unused)) extern cudaError_t __cudaCDP2StreamCreateWithFlags(cudaStream_t * pStream, unsigned flags); 
# 247
__attribute__((unused)) extern cudaError_t __cudaCDP2StreamDestroy(cudaStream_t stream); 
# 248
__attribute__((unused)) extern cudaError_t __cudaCDP2StreamWaitEvent(cudaStream_t stream, cudaEvent_t event, unsigned flags); 
# 249
__attribute__((unused)) extern cudaError_t __cudaCDP2StreamWaitEvent_ptsz(cudaStream_t stream, cudaEvent_t event, unsigned flags); 
# 250
__attribute__((unused)) extern cudaError_t __cudaCDP2EventCreateWithFlags(cudaEvent_t * event, unsigned flags); 
# 251
__attribute__((unused)) extern cudaError_t __cudaCDP2EventRecord(cudaEvent_t event, cudaStream_t stream); 
# 252
__attribute__((unused)) extern cudaError_t __cudaCDP2EventRecord_ptsz(cudaEvent_t event, cudaStream_t stream); 
# 253
__attribute__((unused)) extern cudaError_t __cudaCDP2EventRecordWithFlags(cudaEvent_t event, cudaStream_t stream, unsigned flags); 
# 254
__attribute__((unused)) extern cudaError_t __cudaCDP2EventRecordWithFlags_ptsz(cudaEvent_t event, cudaStream_t stream, unsigned flags); 
# 255
__attribute__((unused)) extern cudaError_t __cudaCDP2EventDestroy(cudaEvent_t event); 
# 256
__attribute__((unused)) extern cudaError_t __cudaCDP2FuncGetAttributes(cudaFuncAttributes * attr, const void * func); 
# 257
__attribute__((unused)) extern cudaError_t __cudaCDP2Free(void * devPtr); 
# 258
__attribute__((unused)) extern cudaError_t __cudaCDP2Malloc(void ** devPtr, size_t size); 
# 259
__attribute__((unused)) extern cudaError_t __cudaCDP2MemcpyAsync(void * dst, const void * src, size_t count, cudaMemcpyKind kind, cudaStream_t stream); 
# 260
__attribute__((unused)) extern cudaError_t __cudaCDP2MemcpyAsync_ptsz(void * dst, const void * src, size_t count, cudaMemcpyKind kind, cudaStream_t stream); 
# 261
__attribute__((unused)) extern cudaError_t __cudaCDP2Memcpy2DAsync(void * dst, size_t dpitch, const void * src, size_t spitch, size_t width, size_t height, cudaMemcpyKind kind, cudaStream_t stream); 
# 262
__attribute__((unused)) extern cudaError_t __cudaCDP2Memcpy2DAsync_ptsz(void * dst, size_t dpitch, const void * src, size_t spitch, size_t width, size_t height, cudaMemcpyKind kind, cudaStream_t stream); 
# 263
__attribute__((unused)) extern cudaError_t __cudaCDP2Memcpy3DAsync(const cudaMemcpy3DParms * p, cudaStream_t stream); 
# 264
__attribute__((unused)) extern cudaError_t __cudaCDP2Memcpy3DAsync_ptsz(const cudaMemcpy3DParms * p, cudaStream_t stream); 
# 265
__attribute__((unused)) extern cudaError_t __cudaCDP2MemsetAsync(void * devPtr, int value, size_t count, cudaStream_t stream); 
# 266
__attribute__((unused)) extern cudaError_t __cudaCDP2MemsetAsync_ptsz(void * devPtr, int value, size_t count, cudaStream_t stream); 
# 267
__attribute__((unused)) extern cudaError_t __cudaCDP2Memset2DAsync(void * devPtr, size_t pitch, int value, size_t width, size_t height, cudaStream_t stream); 
# 268
__attribute__((unused)) extern cudaError_t __cudaCDP2Memset2DAsync_ptsz(void * devPtr, size_t pitch, int value, size_t width, size_t height, cudaStream_t stream); 
# 269
__attribute__((unused)) extern cudaError_t __cudaCDP2Memset3DAsync(cudaPitchedPtr pitchedDevPtr, int value, cudaExtent extent, cudaStream_t stream); 
# 270
__attribute__((unused)) extern cudaError_t __cudaCDP2Memset3DAsync_ptsz(cudaPitchedPtr pitchedDevPtr, int value, cudaExtent extent, cudaStream_t stream); 
# 271
__attribute__((unused)) extern cudaError_t __cudaCDP2RuntimeGetVersion(int * runtimeVersion); 
# 272
__attribute__((unused)) extern void *__cudaCDP2GetParameterBuffer(size_t alignment, size_t size); 
# 273
__attribute__((unused)) extern void *__cudaCDP2GetParameterBufferV2(void * func, dim3 gridDimension, dim3 blockDimension, unsigned sharedMemSize); 
# 274
__attribute__((unused)) extern cudaError_t __cudaCDP2LaunchDevice_ptsz(void * func, void * parameterBuffer, dim3 gridDimension, dim3 blockDimension, unsigned sharedMemSize, cudaStream_t stream); 
# 275
__attribute__((unused)) extern cudaError_t __cudaCDP2LaunchDeviceV2_ptsz(void * parameterBuffer, cudaStream_t stream); 
# 276
__attribute__((unused)) extern cudaError_t __cudaCDP2LaunchDevice(void * func, void * parameterBuffer, dim3 gridDimension, dim3 blockDimension, unsigned sharedMemSize, cudaStream_t stream); 
# 277
__attribute__((unused)) extern cudaError_t __cudaCDP2LaunchDeviceV2(void * parameterBuffer, cudaStream_t stream); 
# 278
__attribute__((unused)) extern cudaError_t __cudaCDP2OccupancyMaxActiveBlocksPerMultiprocessor(int * numBlocks, const void * func, int blockSize, size_t dynamicSmemSize); 
# 279
__attribute__((unused)) extern cudaError_t __cudaCDP2OccupancyMaxActiveBlocksPerMultiprocessorWithFlags(int * numBlocks, const void * func, int blockSize, size_t dynamicSmemSize, unsigned flags); 
# 282
extern cudaError_t cudaGraphLaunch(cudaGraphExec_t graphExec, cudaStream_t stream); 
# 301 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_device_runtime_api.h"
__attribute__((unused)) static inline cudaGraphExec_t cudaGetCurrentGraphExec() 
# 302
{int volatile ___ = 1;
# 306
::exit(___);}
#if 0
# 302
{ 
# 303
unsigned long long current_graph_exec; 
# 304
__asm__("mov.u64 %0, %%current_graph_exec;" : "=l" (current_graph_exec) :); 
# 305
return (cudaGraphExec_t)current_graph_exec; 
# 306
} 
#endif
# 324 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_device_runtime_api.h"
__attribute__((unused)) static inline void cudaTriggerProgrammaticLaunchCompletion() 
# 325
{int volatile ___ = 1;
# 327
::exit(___);}
#if 0
# 325
{ 
# 326
__asm__ volatile("griddepcontrol.launch_dependents;" : :); 
# 327
} 
#endif
# 340 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_device_runtime_api.h"
__attribute__((unused)) static inline void cudaGridDependencySynchronize() 
# 341
{int volatile ___ = 1;
# 343
::exit(___);}
#if 0
# 341
{ 
# 342
__asm__ volatile("griddepcontrol.wait;" : : : "memory"); 
# 343
} 
#endif
# 347 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_device_runtime_api.h"
__attribute__((unused)) extern unsigned long long cudaCGGetIntrinsicHandle(cudaCGScope scope); 
# 348
__attribute__((unused)) extern cudaError_t cudaCGSynchronize(unsigned long long handle, unsigned flags); 
# 349
__attribute__((unused)) extern cudaError_t cudaCGSynchronizeGrid(unsigned long long handle, unsigned flags); 
# 350
__attribute__((unused)) extern cudaError_t cudaCGGetSize(unsigned * numThreads, unsigned * numGrids, unsigned long long handle); 
# 351
__attribute__((unused)) extern cudaError_t cudaCGGetRank(unsigned * threadRank, unsigned * gridRank, unsigned long long handle); 
# 573 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_device_runtime_api.h"
__attribute__((unused)) static inline void *cudaGetParameterBuffer(size_t alignment, size_t size) 
# 574
{int volatile ___ = 1;(void)alignment;(void)size;
# 576
::exit(___);}
#if 0
# 574
{ 
# 575
return __cudaCDP2GetParameterBuffer(alignment, size); 
# 576
} 
#endif
# 609 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_device_runtime_api.h"
__attribute__((unused)) static inline void *cudaGetParameterBufferV2(void *func, dim3 gridDimension, dim3 blockDimension, unsigned sharedMemSize) 
# 610
{int volatile ___ = 1;(void)func;(void)gridDimension;(void)blockDimension;(void)sharedMemSize;
# 612
::exit(___);}
#if 0
# 610
{ 
# 611
return __cudaCDP2GetParameterBufferV2(func, gridDimension, blockDimension, sharedMemSize); 
# 612
} 
#endif
# 619 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_device_runtime_api.h"
__attribute__((unused)) static inline cudaError_t cudaLaunchDevice_ptsz(void *func, void *parameterBuffer, dim3 gridDimension, dim3 blockDimension, unsigned sharedMemSize, cudaStream_t stream) 
# 620
{int volatile ___ = 1;(void)func;(void)parameterBuffer;(void)gridDimension;(void)blockDimension;(void)sharedMemSize;(void)stream;
# 622
::exit(___);}
#if 0
# 620
{ 
# 621
return __cudaCDP2LaunchDevice_ptsz(func, parameterBuffer, gridDimension, blockDimension, sharedMemSize, stream); 
# 622
} 
#endif
# 624 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_device_runtime_api.h"
__attribute__((unused)) static inline cudaError_t cudaLaunchDeviceV2_ptsz(void *parameterBuffer, cudaStream_t stream) 
# 625
{int volatile ___ = 1;(void)parameterBuffer;(void)stream;
# 627
::exit(___);}
#if 0
# 625
{ 
# 626
return __cudaCDP2LaunchDeviceV2_ptsz(parameterBuffer, stream); 
# 627
} 
#endif
# 659 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_device_runtime_api.h"
__attribute__((unused)) static inline cudaError_t cudaLaunchDevice(void *func, void *parameterBuffer, dim3 gridDimension, dim3 blockDimension, unsigned sharedMemSize, cudaStream_t stream) 
# 660
{int volatile ___ = 1;(void)func;(void)parameterBuffer;(void)gridDimension;(void)blockDimension;(void)sharedMemSize;(void)stream;
# 662
::exit(___);}
#if 0
# 660
{ 
# 661
return __cudaCDP2LaunchDevice(func, parameterBuffer, gridDimension, blockDimension, sharedMemSize, stream); 
# 662
} 
#endif
# 664 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_device_runtime_api.h"
__attribute__((unused)) static inline cudaError_t cudaLaunchDeviceV2(void *parameterBuffer, cudaStream_t stream) 
# 665
{int volatile ___ = 1;(void)parameterBuffer;(void)stream;
# 667
::exit(___);}
#if 0
# 665
{ 
# 666
return __cudaCDP2LaunchDeviceV2(parameterBuffer, stream); 
# 667
} 
#endif
# 721 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_device_runtime_api.h"
}
# 723
template< class T> static inline cudaError_t cudaMalloc(T ** devPtr, size_t size); 
# 724
template< class T> static inline cudaError_t cudaFuncGetAttributes(cudaFuncAttributes * attr, T * entry); 
# 725
template< class T> static inline cudaError_t cudaOccupancyMaxActiveBlocksPerMultiprocessor(int * numBlocks, T func, int blockSize, size_t dynamicSmemSize); 
# 726
template< class T> static inline cudaError_t cudaOccupancyMaxActiveBlocksPerMultiprocessorWithFlags(int * numBlocks, T func, int blockSize, size_t dynamicSmemSize, unsigned flags); 
# 273 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern "C" {
# 313 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaDeviceReset(); 
# 335 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaDeviceSynchronize(); 
# 421 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaDeviceSetLimit(cudaLimit limit, size_t value); 
# 457 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaDeviceGetLimit(size_t * pValue, cudaLimit limit); 
# 480 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaDeviceGetTexture1DLinearMaxWidth(size_t * maxWidthInElements, const cudaChannelFormatDesc * fmtDesc, int device); 
# 514 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaDeviceGetCacheConfig(cudaFuncCache * pCacheConfig); 
# 551 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaDeviceGetStreamPriorityRange(int * leastPriority, int * greatestPriority); 
# 595 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaDeviceSetCacheConfig(cudaFuncCache cacheConfig); 
# 626 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaDeviceGetSharedMemConfig(cudaSharedMemConfig * pConfig); 
# 670 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaDeviceSetSharedMemConfig(cudaSharedMemConfig config); 
# 697 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaDeviceGetByPCIBusId(int * device, const char * pciBusId); 
# 727 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaDeviceGetPCIBusId(char * pciBusId, int len, int device); 
# 777 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaIpcGetEventHandle(cudaIpcEventHandle_t * handle, cudaEvent_t event); 
# 820 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaIpcOpenEventHandle(cudaEvent_t * event, cudaIpcEventHandle_t handle); 
# 864 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaIpcGetMemHandle(cudaIpcMemHandle_t * handle, void * devPtr); 
# 930 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaIpcOpenMemHandle(void ** devPtr, cudaIpcMemHandle_t handle, unsigned flags); 
# 968 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaIpcCloseMemHandle(void * devPtr); 
# 1000 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaDeviceFlushGPUDirectRDMAWrites(cudaFlushGPUDirectRDMAWritesTarget target, cudaFlushGPUDirectRDMAWritesScope scope); 
# 1043 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
__attribute((deprecated)) extern cudaError_t cudaThreadExit(); 
# 1069 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
__attribute((deprecated)) extern cudaError_t cudaThreadSynchronize(); 
# 1118 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
__attribute((deprecated)) extern cudaError_t cudaThreadSetLimit(cudaLimit limit, size_t value); 
# 1151 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
__attribute((deprecated)) extern cudaError_t cudaThreadGetLimit(size_t * pValue, cudaLimit limit); 
# 1187 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
__attribute((deprecated)) extern cudaError_t cudaThreadGetCacheConfig(cudaFuncCache * pCacheConfig); 
# 1234 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
__attribute((deprecated)) extern cudaError_t cudaThreadSetCacheConfig(cudaFuncCache cacheConfig); 
# 1299 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaGetLastError(); 
# 1350 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaPeekAtLastError(); 
# 1366 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern const char *cudaGetErrorName(cudaError_t error); 
# 1382 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern const char *cudaGetErrorString(cudaError_t error); 
# 1411 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaGetDeviceCount(int * count); 
# 1716 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaGetDeviceProperties_v2(cudaDeviceProp * prop, int device); 
# 1918 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaDeviceGetAttribute(int * value, cudaDeviceAttr attr, int device); 
# 1936 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaDeviceGetDefaultMemPool(cudaMemPool_t * memPool, int device); 
# 1960 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaDeviceSetMemPool(int device, cudaMemPool_t memPool); 
# 1980 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaDeviceGetMemPool(cudaMemPool_t * memPool, int device); 
# 2042 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaDeviceGetNvSciSyncAttributes(void * nvSciSyncAttrList, int device, int flags); 
# 2082 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaDeviceGetP2PAttribute(int * value, cudaDeviceP2PAttr attr, int srcDevice, int dstDevice); 
# 2104 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaChooseDevice(int * device, const cudaDeviceProp * prop); 
# 2133 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaInitDevice(int device, unsigned deviceFlags, unsigned flags); 
# 2179 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaSetDevice(int device); 
# 2201 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaGetDevice(int * device); 
# 2232 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaSetValidDevices(int * device_arr, int len); 
# 2298 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaSetDeviceFlags(unsigned flags); 
# 2343 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaGetDeviceFlags(unsigned * flags); 
# 2383 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaStreamCreate(cudaStream_t * pStream); 
# 2415 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaStreamCreateWithFlags(cudaStream_t * pStream, unsigned flags); 
# 2461 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaStreamCreateWithPriority(cudaStream_t * pStream, unsigned flags, int priority); 
# 2488 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaStreamGetPriority(cudaStream_t hStream, int * priority); 
# 2513 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaStreamGetFlags(cudaStream_t hStream, unsigned * flags); 
# 2550 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaStreamGetId(cudaStream_t hStream, unsigned long long * streamId); 
# 2565 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaCtxResetPersistingL2Cache(); 
# 2585 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaStreamCopyAttributes(cudaStream_t dst, cudaStream_t src); 
# 2606 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaStreamGetAttribute(cudaStream_t hStream, cudaLaunchAttributeID attr, cudaLaunchAttributeValue * value_out); 
# 2630 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaStreamSetAttribute(cudaStream_t hStream, cudaLaunchAttributeID attr, const cudaLaunchAttributeValue * value); 
# 2664 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaStreamDestroy(cudaStream_t stream); 
# 2695 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaStreamWaitEvent(cudaStream_t stream, cudaEvent_t event, unsigned flags = 0); 
# 2703
typedef void (*cudaStreamCallback_t)(cudaStream_t stream, cudaError_t status, void * userData); 
# 2770 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaStreamAddCallback(cudaStream_t stream, cudaStreamCallback_t callback, void * userData, unsigned flags); 
# 2794 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaStreamSynchronize(cudaStream_t stream); 
# 2819 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaStreamQuery(cudaStream_t stream); 
# 2903 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaStreamAttachMemAsync(cudaStream_t stream, void * devPtr, size_t length = 0, unsigned flags = 4); 
# 2942 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaStreamBeginCapture(cudaStream_t stream, cudaStreamCaptureMode mode); 
# 2993 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaThreadExchangeStreamCaptureMode(cudaStreamCaptureMode * mode); 
# 3021 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaStreamEndCapture(cudaStream_t stream, cudaGraph_t * pGraph); 
# 3059 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaStreamIsCapturing(cudaStream_t stream, cudaStreamCaptureStatus * pCaptureStatus); 
# 3107 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaStreamGetCaptureInfo_v2(cudaStream_t stream, cudaStreamCaptureStatus * captureStatus_out, unsigned long long * id_out = 0, cudaGraph_t * graph_out = 0, const cudaGraphNode_t ** dependencies_out = 0, size_t * numDependencies_out = 0); 
# 3139 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaStreamUpdateCaptureDependencies(cudaStream_t stream, cudaGraphNode_t * dependencies, size_t numDependencies, unsigned flags = 0); 
# 3176 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaEventCreate(cudaEvent_t * event); 
# 3213 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaEventCreateWithFlags(cudaEvent_t * event, unsigned flags); 
# 3254 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaEventRecord(cudaEvent_t event, cudaStream_t stream = 0); 
# 3302 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaEventRecordWithFlags(cudaEvent_t event, cudaStream_t stream = 0, unsigned flags = 0); 
# 3335 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaEventQuery(cudaEvent_t event); 
# 3366 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaEventSynchronize(cudaEvent_t event); 
# 3396 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaEventDestroy(cudaEvent_t event); 
# 3441 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaEventElapsedTime(float * ms, cudaEvent_t start, cudaEvent_t end); 
# 3622 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaImportExternalMemory(cudaExternalMemory_t * extMem_out, const cudaExternalMemoryHandleDesc * memHandleDesc); 
# 3677 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaExternalMemoryGetMappedBuffer(void ** devPtr, cudaExternalMemory_t extMem, const cudaExternalMemoryBufferDesc * bufferDesc); 
# 3737 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaExternalMemoryGetMappedMipmappedArray(cudaMipmappedArray_t * mipmap, cudaExternalMemory_t extMem, const cudaExternalMemoryMipmappedArrayDesc * mipmapDesc); 
# 3761 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaDestroyExternalMemory(cudaExternalMemory_t extMem); 
# 3915 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaImportExternalSemaphore(cudaExternalSemaphore_t * extSem_out, const cudaExternalSemaphoreHandleDesc * semHandleDesc); 
# 3998 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaSignalExternalSemaphoresAsync_v2(const cudaExternalSemaphore_t * extSemArray, const cudaExternalSemaphoreSignalParams * paramsArray, unsigned numExtSems, cudaStream_t stream = 0); 
# 4074 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaWaitExternalSemaphoresAsync_v2(const cudaExternalSemaphore_t * extSemArray, const cudaExternalSemaphoreWaitParams * paramsArray, unsigned numExtSems, cudaStream_t stream = 0); 
# 4097 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaDestroyExternalSemaphore(cudaExternalSemaphore_t extSem); 
# 4164 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaLaunchKernel(const void * func, dim3 gridDim, dim3 blockDim, void ** args, size_t sharedMem, cudaStream_t stream); 
# 4226 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaLaunchKernelExC(const cudaLaunchConfig_t * config, const void * func, void ** args); 
# 4283 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaLaunchCooperativeKernel(const void * func, dim3 gridDim, dim3 blockDim, void ** args, size_t sharedMem, cudaStream_t stream); 
# 4384 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
__attribute((deprecated)) extern cudaError_t cudaLaunchCooperativeKernelMultiDevice(cudaLaunchParams * launchParamsList, unsigned numDevices, unsigned flags = 0); 
# 4429 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaFuncSetCacheConfig(const void * func, cudaFuncCache cacheConfig); 
# 4484 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaFuncSetSharedMemConfig(const void * func, cudaSharedMemConfig config); 
# 4517 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaFuncGetAttributes(cudaFuncAttributes * attr, const void * func); 
# 4554 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaFuncSetAttribute(const void * func, cudaFuncAttribute attr, int value); 
# 4578 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
__attribute((deprecated)) extern cudaError_t cudaSetDoubleForDevice(double * d); 
# 4602 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
__attribute((deprecated)) extern cudaError_t cudaSetDoubleForHost(double * d); 
# 4668 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaLaunchHostFunc(cudaStream_t stream, cudaHostFn_t fn, void * userData); 
# 4725 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaOccupancyMaxActiveBlocksPerMultiprocessor(int * numBlocks, const void * func, int blockSize, size_t dynamicSMemSize); 
# 4754 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaOccupancyAvailableDynamicSMemPerBlock(size_t * dynamicSmemSize, const void * func, int numBlocks, int blockSize); 
# 4799 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaOccupancyMaxActiveBlocksPerMultiprocessorWithFlags(int * numBlocks, const void * func, int blockSize, size_t dynamicSMemSize, unsigned flags); 
# 4834 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaOccupancyMaxPotentialClusterSize(int * clusterSize, const void * func, const cudaLaunchConfig_t * launchConfig); 
# 4873 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaOccupancyMaxActiveClusters(int * numClusters, const void * func, const cudaLaunchConfig_t * launchConfig); 
# 4993 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaMallocManaged(void ** devPtr, size_t size, unsigned flags = 1); 
# 5026 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaMalloc(void ** devPtr, size_t size); 
# 5063 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaMallocHost(void ** ptr, size_t size); 
# 5106 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaMallocPitch(void ** devPtr, size_t * pitch, size_t width, size_t height); 
# 5158 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaMallocArray(cudaArray_t * array, const cudaChannelFormatDesc * desc, size_t width, size_t height = 0, unsigned flags = 0); 
# 5196 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaFree(void * devPtr); 
# 5219 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaFreeHost(void * ptr); 
# 5242 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaFreeArray(cudaArray_t array); 
# 5265 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaFreeMipmappedArray(cudaMipmappedArray_t mipmappedArray); 
# 5331 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaHostAlloc(void ** pHost, size_t size, unsigned flags); 
# 5428 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaHostRegister(void * ptr, size_t size, unsigned flags); 
# 5451 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaHostUnregister(void * ptr); 
# 5496 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaHostGetDevicePointer(void ** pDevice, void * pHost, unsigned flags); 
# 5518 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaHostGetFlags(unsigned * pFlags, void * pHost); 
# 5557 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaMalloc3D(cudaPitchedPtr * pitchedDevPtr, cudaExtent extent); 
# 5702 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaMalloc3DArray(cudaArray_t * array, const cudaChannelFormatDesc * desc, cudaExtent extent, unsigned flags = 0); 
# 5847 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaMallocMipmappedArray(cudaMipmappedArray_t * mipmappedArray, const cudaChannelFormatDesc * desc, cudaExtent extent, unsigned numLevels, unsigned flags = 0); 
# 5880 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaGetMipmappedArrayLevel(cudaArray_t * levelArray, cudaMipmappedArray_const_t mipmappedArray, unsigned level); 
# 5985 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaMemcpy3D(const cudaMemcpy3DParms * p); 
# 6017 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaMemcpy3DPeer(const cudaMemcpy3DPeerParms * p); 
# 6135 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaMemcpy3DAsync(const cudaMemcpy3DParms * p, cudaStream_t stream = 0); 
# 6162 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaMemcpy3DPeerAsync(const cudaMemcpy3DPeerParms * p, cudaStream_t stream = 0); 
# 6196 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaMemGetInfo(size_t * free, size_t * total); 
# 6222 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaArrayGetInfo(cudaChannelFormatDesc * desc, cudaExtent * extent, unsigned * flags, cudaArray_t array); 
# 6251 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaArrayGetPlane(cudaArray_t * pPlaneArray, cudaArray_t hArray, unsigned planeIdx); 
# 6274 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaArrayGetMemoryRequirements(cudaArrayMemoryRequirements * memoryRequirements, cudaArray_t array, int device); 
# 6298 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaMipmappedArrayGetMemoryRequirements(cudaArrayMemoryRequirements * memoryRequirements, cudaMipmappedArray_t mipmap, int device); 
# 6326 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaArrayGetSparseProperties(cudaArraySparseProperties * sparseProperties, cudaArray_t array); 
# 6356 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaMipmappedArrayGetSparseProperties(cudaArraySparseProperties * sparseProperties, cudaMipmappedArray_t mipmap); 
# 6401 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaMemcpy(void * dst, const void * src, size_t count, cudaMemcpyKind kind); 
# 6436 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaMemcpyPeer(void * dst, int dstDevice, const void * src, int srcDevice, size_t count); 
# 6485 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaMemcpy2D(void * dst, size_t dpitch, const void * src, size_t spitch, size_t width, size_t height, cudaMemcpyKind kind); 
# 6535 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaMemcpy2DToArray(cudaArray_t dst, size_t wOffset, size_t hOffset, const void * src, size_t spitch, size_t width, size_t height, cudaMemcpyKind kind); 
# 6585 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaMemcpy2DFromArray(void * dst, size_t dpitch, cudaArray_const_t src, size_t wOffset, size_t hOffset, size_t width, size_t height, cudaMemcpyKind kind); 
# 6632 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaMemcpy2DArrayToArray(cudaArray_t dst, size_t wOffsetDst, size_t hOffsetDst, cudaArray_const_t src, size_t wOffsetSrc, size_t hOffsetSrc, size_t width, size_t height, cudaMemcpyKind kind = cudaMemcpyDeviceToDevice); 
# 6675 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaMemcpyToSymbol(const void * symbol, const void * src, size_t count, size_t offset = 0, cudaMemcpyKind kind = cudaMemcpyHostToDevice); 
# 6718 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaMemcpyFromSymbol(void * dst, const void * symbol, size_t count, size_t offset = 0, cudaMemcpyKind kind = cudaMemcpyDeviceToHost); 
# 6775 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaMemcpyAsync(void * dst, const void * src, size_t count, cudaMemcpyKind kind, cudaStream_t stream = 0); 
# 6810 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaMemcpyPeerAsync(void * dst, int dstDevice, const void * src, int srcDevice, size_t count, cudaStream_t stream = 0); 
# 6873 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaMemcpy2DAsync(void * dst, size_t dpitch, const void * src, size_t spitch, size_t width, size_t height, cudaMemcpyKind kind, cudaStream_t stream = 0); 
# 6931 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaMemcpy2DToArrayAsync(cudaArray_t dst, size_t wOffset, size_t hOffset, const void * src, size_t spitch, size_t width, size_t height, cudaMemcpyKind kind, cudaStream_t stream = 0); 
# 6988 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaMemcpy2DFromArrayAsync(void * dst, size_t dpitch, cudaArray_const_t src, size_t wOffset, size_t hOffset, size_t width, size_t height, cudaMemcpyKind kind, cudaStream_t stream = 0); 
# 7039 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaMemcpyToSymbolAsync(const void * symbol, const void * src, size_t count, size_t offset, cudaMemcpyKind kind, cudaStream_t stream = 0); 
# 7090 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaMemcpyFromSymbolAsync(void * dst, const void * symbol, size_t count, size_t offset, cudaMemcpyKind kind, cudaStream_t stream = 0); 
# 7119 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaMemset(void * devPtr, int value, size_t count); 
# 7153 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaMemset2D(void * devPtr, size_t pitch, int value, size_t width, size_t height); 
# 7199 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaMemset3D(cudaPitchedPtr pitchedDevPtr, int value, cudaExtent extent); 
# 7235 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaMemsetAsync(void * devPtr, int value, size_t count, cudaStream_t stream = 0); 
# 7276 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaMemset2DAsync(void * devPtr, size_t pitch, int value, size_t width, size_t height, cudaStream_t stream = 0); 
# 7329 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaMemset3DAsync(cudaPitchedPtr pitchedDevPtr, int value, cudaExtent extent, cudaStream_t stream = 0); 
# 7357 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaGetSymbolAddress(void ** devPtr, const void * symbol); 
# 7384 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaGetSymbolSize(size_t * size, const void * symbol); 
# 7454 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaMemPrefetchAsync(const void * devPtr, size_t count, int dstDevice, cudaStream_t stream = 0); 
# 7456
extern cudaError_t cudaMemPrefetchAsync_v2(const void * devPtr, size_t count, cudaMemLocation location, unsigned flags, cudaStream_t stream = 0); 
# 7570 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaMemAdvise(const void * devPtr, size_t count, cudaMemoryAdvise advice, int device); 
# 7693 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaMemAdvise_v2(const void * devPtr, size_t count, cudaMemoryAdvise advice, cudaMemLocation location); 
# 7775 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaMemRangeGetAttribute(void * data, size_t dataSize, cudaMemRangeAttribute attribute, const void * devPtr, size_t count); 
# 7818 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaMemRangeGetAttributes(void ** data, size_t * dataSizes, cudaMemRangeAttribute * attributes, size_t numAttributes, const void * devPtr, size_t count); 
# 7878 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
__attribute((deprecated)) extern cudaError_t cudaMemcpyToArray(cudaArray_t dst, size_t wOffset, size_t hOffset, const void * src, size_t count, cudaMemcpyKind kind); 
# 7920 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
__attribute((deprecated)) extern cudaError_t cudaMemcpyFromArray(void * dst, cudaArray_const_t src, size_t wOffset, size_t hOffset, size_t count, cudaMemcpyKind kind); 
# 7963 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
__attribute((deprecated)) extern cudaError_t cudaMemcpyArrayToArray(cudaArray_t dst, size_t wOffsetDst, size_t hOffsetDst, cudaArray_const_t src, size_t wOffsetSrc, size_t hOffsetSrc, size_t count, cudaMemcpyKind kind = cudaMemcpyDeviceToDevice); 
# 8014 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
__attribute((deprecated)) extern cudaError_t cudaMemcpyToArrayAsync(cudaArray_t dst, size_t wOffset, size_t hOffset, const void * src, size_t count, cudaMemcpyKind kind, cudaStream_t stream = 0); 
# 8064 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
__attribute((deprecated)) extern cudaError_t cudaMemcpyFromArrayAsync(void * dst, cudaArray_const_t src, size_t wOffset, size_t hOffset, size_t count, cudaMemcpyKind kind, cudaStream_t stream = 0); 
# 8133 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaMallocAsync(void ** devPtr, size_t size, cudaStream_t hStream); 
# 8159 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaFreeAsync(void * devPtr, cudaStream_t hStream); 
# 8184 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaMemPoolTrimTo(cudaMemPool_t memPool, size_t minBytesToKeep); 
# 8228 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaMemPoolSetAttribute(cudaMemPool_t memPool, cudaMemPoolAttr attr, void * value); 
# 8276 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaMemPoolGetAttribute(cudaMemPool_t memPool, cudaMemPoolAttr attr, void * value); 
# 8291 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaMemPoolSetAccess(cudaMemPool_t memPool, const cudaMemAccessDesc * descList, size_t count); 
# 8304 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaMemPoolGetAccess(cudaMemAccessFlags * flags, cudaMemPool_t memPool, cudaMemLocation * location); 
# 8331 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaMemPoolCreate(cudaMemPool_t * memPool, const cudaMemPoolProps * poolProps); 
# 8353 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaMemPoolDestroy(cudaMemPool_t memPool); 
# 8389 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaMallocFromPoolAsync(void ** ptr, size_t size, cudaMemPool_t memPool, cudaStream_t stream); 
# 8414 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaMemPoolExportToShareableHandle(void * shareableHandle, cudaMemPool_t memPool, cudaMemAllocationHandleType handleType, unsigned flags); 
# 8441 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaMemPoolImportFromShareableHandle(cudaMemPool_t * memPool, void * shareableHandle, cudaMemAllocationHandleType handleType, unsigned flags); 
# 8464 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaMemPoolExportPointer(cudaMemPoolPtrExportData * exportData, void * ptr); 
# 8493 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaMemPoolImportPointer(void ** ptr, cudaMemPool_t memPool, cudaMemPoolPtrExportData * exportData); 
# 8646 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaPointerGetAttributes(cudaPointerAttributes * attributes, const void * ptr); 
# 8687 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaDeviceCanAccessPeer(int * canAccessPeer, int device, int peerDevice); 
# 8729 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaDeviceEnablePeerAccess(int peerDevice, unsigned flags); 
# 8751 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaDeviceDisablePeerAccess(int peerDevice); 
# 8815 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaGraphicsUnregisterResource(cudaGraphicsResource_t resource); 
# 8850 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaGraphicsResourceSetMapFlags(cudaGraphicsResource_t resource, unsigned flags); 
# 8889 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaGraphicsMapResources(int count, cudaGraphicsResource_t * resources, cudaStream_t stream = 0); 
# 8924 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaGraphicsUnmapResources(int count, cudaGraphicsResource_t * resources, cudaStream_t stream = 0); 
# 8956 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaGraphicsResourceGetMappedPointer(void ** devPtr, size_t * size, cudaGraphicsResource_t resource); 
# 8994 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaGraphicsSubResourceGetMappedArray(cudaArray_t * array, cudaGraphicsResource_t resource, unsigned arrayIndex, unsigned mipLevel); 
# 9023 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaGraphicsResourceGetMappedMipmappedArray(cudaMipmappedArray_t * mipmappedArray, cudaGraphicsResource_t resource); 
# 9058 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaGetChannelDesc(cudaChannelFormatDesc * desc, cudaArray_const_t array); 
# 9088 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaChannelFormatDesc cudaCreateChannelDesc(int x, int y, int z, int w, cudaChannelFormatKind f); 
# 9312 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaCreateTextureObject(cudaTextureObject_t * pTexObject, const cudaResourceDesc * pResDesc, const cudaTextureDesc * pTexDesc, const cudaResourceViewDesc * pResViewDesc); 
# 9332 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaDestroyTextureObject(cudaTextureObject_t texObject); 
# 9352 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaGetTextureObjectResourceDesc(cudaResourceDesc * pResDesc, cudaTextureObject_t texObject); 
# 9372 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaGetTextureObjectTextureDesc(cudaTextureDesc * pTexDesc, cudaTextureObject_t texObject); 
# 9393 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaGetTextureObjectResourceViewDesc(cudaResourceViewDesc * pResViewDesc, cudaTextureObject_t texObject); 
# 9438 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaCreateSurfaceObject(cudaSurfaceObject_t * pSurfObject, const cudaResourceDesc * pResDesc); 
# 9458 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaDestroySurfaceObject(cudaSurfaceObject_t surfObject); 
# 9477 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaGetSurfaceObjectResourceDesc(cudaResourceDesc * pResDesc, cudaSurfaceObject_t surfObject); 
# 9511 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaDriverGetVersion(int * driverVersion); 
# 9540 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaRuntimeGetVersion(int * runtimeVersion); 
# 9587 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaGraphCreate(cudaGraph_t * pGraph, unsigned flags); 
# 9685 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaGraphAddKernelNode(cudaGraphNode_t * pGraphNode, cudaGraph_t graph, const cudaGraphNode_t * pDependencies, size_t numDependencies, const cudaKernelNodeParams * pNodeParams); 
# 9718 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaGraphKernelNodeGetParams(cudaGraphNode_t node, cudaKernelNodeParams * pNodeParams); 
# 9744 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaGraphKernelNodeSetParams(cudaGraphNode_t node, const cudaKernelNodeParams * pNodeParams); 
# 9764 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaGraphKernelNodeCopyAttributes(cudaGraphNode_t hSrc, cudaGraphNode_t hDst); 
# 9787 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaGraphKernelNodeGetAttribute(cudaGraphNode_t hNode, cudaLaunchAttributeID attr, cudaLaunchAttributeValue * value_out); 
# 9811 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaGraphKernelNodeSetAttribute(cudaGraphNode_t hNode, cudaLaunchAttributeID attr, const cudaLaunchAttributeValue * value); 
# 9862 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaGraphAddMemcpyNode(cudaGraphNode_t * pGraphNode, cudaGraph_t graph, const cudaGraphNode_t * pDependencies, size_t numDependencies, const cudaMemcpy3DParms * pCopyParams); 
# 9921 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaGraphAddMemcpyNodeToSymbol(cudaGraphNode_t * pGraphNode, cudaGraph_t graph, const cudaGraphNode_t * pDependencies, size_t numDependencies, const void * symbol, const void * src, size_t count, size_t offset, cudaMemcpyKind kind); 
# 9990 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaGraphAddMemcpyNodeFromSymbol(cudaGraphNode_t * pGraphNode, cudaGraph_t graph, const cudaGraphNode_t * pDependencies, size_t numDependencies, void * dst, const void * symbol, size_t count, size_t offset, cudaMemcpyKind kind); 
# 10058 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaGraphAddMemcpyNode1D(cudaGraphNode_t * pGraphNode, cudaGraph_t graph, const cudaGraphNode_t * pDependencies, size_t numDependencies, void * dst, const void * src, size_t count, cudaMemcpyKind kind); 
# 10090 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaGraphMemcpyNodeGetParams(cudaGraphNode_t node, cudaMemcpy3DParms * pNodeParams); 
# 10117 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaGraphMemcpyNodeSetParams(cudaGraphNode_t node, const cudaMemcpy3DParms * pNodeParams); 
# 10156 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaGraphMemcpyNodeSetParamsToSymbol(cudaGraphNode_t node, const void * symbol, const void * src, size_t count, size_t offset, cudaMemcpyKind kind); 
# 10202 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaGraphMemcpyNodeSetParamsFromSymbol(cudaGraphNode_t node, void * dst, const void * symbol, size_t count, size_t offset, cudaMemcpyKind kind); 
# 10248 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaGraphMemcpyNodeSetParams1D(cudaGraphNode_t node, void * dst, const void * src, size_t count, cudaMemcpyKind kind); 
# 10296 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaGraphAddMemsetNode(cudaGraphNode_t * pGraphNode, cudaGraph_t graph, const cudaGraphNode_t * pDependencies, size_t numDependencies, const cudaMemsetParams * pMemsetParams); 
# 10319 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaGraphMemsetNodeGetParams(cudaGraphNode_t node, cudaMemsetParams * pNodeParams); 
# 10343 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaGraphMemsetNodeSetParams(cudaGraphNode_t node, const cudaMemsetParams * pNodeParams); 
# 10385 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaGraphAddHostNode(cudaGraphNode_t * pGraphNode, cudaGraph_t graph, const cudaGraphNode_t * pDependencies, size_t numDependencies, const cudaHostNodeParams * pNodeParams); 
# 10408 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaGraphHostNodeGetParams(cudaGraphNode_t node, cudaHostNodeParams * pNodeParams); 
# 10432 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaGraphHostNodeSetParams(cudaGraphNode_t node, const cudaHostNodeParams * pNodeParams); 
# 10473 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaGraphAddChildGraphNode(cudaGraphNode_t * pGraphNode, cudaGraph_t graph, const cudaGraphNode_t * pDependencies, size_t numDependencies, cudaGraph_t childGraph); 
# 10500 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaGraphChildGraphNodeGetGraph(cudaGraphNode_t node, cudaGraph_t * pGraph); 
# 10538 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaGraphAddEmptyNode(cudaGraphNode_t * pGraphNode, cudaGraph_t graph, const cudaGraphNode_t * pDependencies, size_t numDependencies); 
# 10582 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaGraphAddEventRecordNode(cudaGraphNode_t * pGraphNode, cudaGraph_t graph, const cudaGraphNode_t * pDependencies, size_t numDependencies, cudaEvent_t event); 
# 10609 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaGraphEventRecordNodeGetEvent(cudaGraphNode_t node, cudaEvent_t * event_out); 
# 10637 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaGraphEventRecordNodeSetEvent(cudaGraphNode_t node, cudaEvent_t event); 
# 10684 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaGraphAddEventWaitNode(cudaGraphNode_t * pGraphNode, cudaGraph_t graph, const cudaGraphNode_t * pDependencies, size_t numDependencies, cudaEvent_t event); 
# 10711 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaGraphEventWaitNodeGetEvent(cudaGraphNode_t node, cudaEvent_t * event_out); 
# 10739 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaGraphEventWaitNodeSetEvent(cudaGraphNode_t node, cudaEvent_t event); 
# 10789 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaGraphAddExternalSemaphoresSignalNode(cudaGraphNode_t * pGraphNode, cudaGraph_t graph, const cudaGraphNode_t * pDependencies, size_t numDependencies, const cudaExternalSemaphoreSignalNodeParams * nodeParams); 
# 10822 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaGraphExternalSemaphoresSignalNodeGetParams(cudaGraphNode_t hNode, cudaExternalSemaphoreSignalNodeParams * params_out); 
# 10850 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaGraphExternalSemaphoresSignalNodeSetParams(cudaGraphNode_t hNode, const cudaExternalSemaphoreSignalNodeParams * nodeParams); 
# 10900 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaGraphAddExternalSemaphoresWaitNode(cudaGraphNode_t * pGraphNode, cudaGraph_t graph, const cudaGraphNode_t * pDependencies, size_t numDependencies, const cudaExternalSemaphoreWaitNodeParams * nodeParams); 
# 10933 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaGraphExternalSemaphoresWaitNodeGetParams(cudaGraphNode_t hNode, cudaExternalSemaphoreWaitNodeParams * params_out); 
# 10961 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaGraphExternalSemaphoresWaitNodeSetParams(cudaGraphNode_t hNode, const cudaExternalSemaphoreWaitNodeParams * nodeParams); 
# 11039 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaGraphAddMemAllocNode(cudaGraphNode_t * pGraphNode, cudaGraph_t graph, const cudaGraphNode_t * pDependencies, size_t numDependencies, cudaMemAllocNodeParams * nodeParams); 
# 11066 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaGraphMemAllocNodeGetParams(cudaGraphNode_t node, cudaMemAllocNodeParams * params_out); 
# 11127 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaGraphAddMemFreeNode(cudaGraphNode_t * pGraphNode, cudaGraph_t graph, const cudaGraphNode_t * pDependencies, size_t numDependencies, void * dptr); 
# 11151 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaGraphMemFreeNodeGetParams(cudaGraphNode_t node, void * dptr_out); 
# 11179 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaDeviceGraphMemTrim(int device); 
# 11216 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaDeviceGetGraphMemAttribute(int device, cudaGraphMemAttributeType attr, void * value); 
# 11250 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaDeviceSetGraphMemAttribute(int device, cudaGraphMemAttributeType attr, void * value); 
# 11278 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaGraphClone(cudaGraph_t * pGraphClone, cudaGraph_t originalGraph); 
# 11306 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaGraphNodeFindInClone(cudaGraphNode_t * pNode, cudaGraphNode_t originalNode, cudaGraph_t clonedGraph); 
# 11337 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaGraphNodeGetType(cudaGraphNode_t node, cudaGraphNodeType * pType); 
# 11368 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaGraphGetNodes(cudaGraph_t graph, cudaGraphNode_t * nodes, size_t * numNodes); 
# 11399 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaGraphGetRootNodes(cudaGraph_t graph, cudaGraphNode_t * pRootNodes, size_t * pNumRootNodes); 
# 11433 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaGraphGetEdges(cudaGraph_t graph, cudaGraphNode_t * from, cudaGraphNode_t * to, size_t * numEdges); 
# 11464 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaGraphNodeGetDependencies(cudaGraphNode_t node, cudaGraphNode_t * pDependencies, size_t * pNumDependencies); 
# 11496 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaGraphNodeGetDependentNodes(cudaGraphNode_t node, cudaGraphNode_t * pDependentNodes, size_t * pNumDependentNodes); 
# 11527 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaGraphAddDependencies(cudaGraph_t graph, const cudaGraphNode_t * from, const cudaGraphNode_t * to, size_t numDependencies); 
# 11558 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaGraphRemoveDependencies(cudaGraph_t graph, const cudaGraphNode_t * from, const cudaGraphNode_t * to, size_t numDependencies); 
# 11588 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaGraphDestroyNode(cudaGraphNode_t node); 
# 11650 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaGraphInstantiate(cudaGraphExec_t * pGraphExec, cudaGraph_t graph, unsigned long long flags = 0); 
# 11721 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaGraphInstantiateWithFlags(cudaGraphExec_t * pGraphExec, cudaGraph_t graph, unsigned long long flags = 0); 
# 11826 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaGraphInstantiateWithParams(cudaGraphExec_t * pGraphExec, cudaGraph_t graph, cudaGraphInstantiateParams * instantiateParams); 
# 11851 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaGraphExecGetFlags(cudaGraphExec_t graphExec, unsigned long long * flags); 
# 11903 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaGraphExecKernelNodeSetParams(cudaGraphExec_t hGraphExec, cudaGraphNode_t node, const cudaKernelNodeParams * pNodeParams); 
# 11954 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaGraphExecMemcpyNodeSetParams(cudaGraphExec_t hGraphExec, cudaGraphNode_t node, const cudaMemcpy3DParms * pNodeParams); 
# 12009 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaGraphExecMemcpyNodeSetParamsToSymbol(cudaGraphExec_t hGraphExec, cudaGraphNode_t node, const void * symbol, const void * src, size_t count, size_t offset, cudaMemcpyKind kind); 
# 12072 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaGraphExecMemcpyNodeSetParamsFromSymbol(cudaGraphExec_t hGraphExec, cudaGraphNode_t node, void * dst, const void * symbol, size_t count, size_t offset, cudaMemcpyKind kind); 
# 12133 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaGraphExecMemcpyNodeSetParams1D(cudaGraphExec_t hGraphExec, cudaGraphNode_t node, void * dst, const void * src, size_t count, cudaMemcpyKind kind); 
# 12188 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaGraphExecMemsetNodeSetParams(cudaGraphExec_t hGraphExec, cudaGraphNode_t node, const cudaMemsetParams * pNodeParams); 
# 12228 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaGraphExecHostNodeSetParams(cudaGraphExec_t hGraphExec, cudaGraphNode_t node, const cudaHostNodeParams * pNodeParams); 
# 12275 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaGraphExecChildGraphNodeSetParams(cudaGraphExec_t hGraphExec, cudaGraphNode_t node, cudaGraph_t childGraph); 
# 12320 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaGraphExecEventRecordNodeSetEvent(cudaGraphExec_t hGraphExec, cudaGraphNode_t hNode, cudaEvent_t event); 
# 12365 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaGraphExecEventWaitNodeSetEvent(cudaGraphExec_t hGraphExec, cudaGraphNode_t hNode, cudaEvent_t event); 
# 12413 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaGraphExecExternalSemaphoresSignalNodeSetParams(cudaGraphExec_t hGraphExec, cudaGraphNode_t hNode, const cudaExternalSemaphoreSignalNodeParams * nodeParams); 
# 12461 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaGraphExecExternalSemaphoresWaitNodeSetParams(cudaGraphExec_t hGraphExec, cudaGraphNode_t hNode, const cudaExternalSemaphoreWaitNodeParams * nodeParams); 
# 12501 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaGraphNodeSetEnabled(cudaGraphExec_t hGraphExec, cudaGraphNode_t hNode, unsigned isEnabled); 
# 12535 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaGraphNodeGetEnabled(cudaGraphExec_t hGraphExec, cudaGraphNode_t hNode, unsigned * isEnabled); 
# 12620 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaGraphExecUpdate(cudaGraphExec_t hGraphExec, cudaGraph_t hGraph, cudaGraphExecUpdateResultInfo * resultInfo); 
# 12645 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaGraphUpload(cudaGraphExec_t graphExec, cudaStream_t stream); 
# 12676 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaGraphLaunch(cudaGraphExec_t graphExec, cudaStream_t stream); 
# 12699 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaGraphExecDestroy(cudaGraphExec_t graphExec); 
# 12720 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaGraphDestroy(cudaGraph_t graph); 
# 12739 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaGraphDebugDotPrint(cudaGraph_t graph, const char * path, unsigned flags); 
# 12775 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaUserObjectCreate(cudaUserObject_t * object_out, void * ptr, cudaHostFn_t destroy, unsigned initialRefcount, unsigned flags); 
# 12799 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaUserObjectRetain(cudaUserObject_t object, unsigned count = 1); 
# 12827 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaUserObjectRelease(cudaUserObject_t object, unsigned count = 1); 
# 12855 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaGraphRetainUserObject(cudaGraph_t graph, cudaUserObject_t object, unsigned count = 1, unsigned flags = 0); 
# 12880 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaGraphReleaseUserObject(cudaGraph_t graph, cudaUserObject_t object, unsigned count = 1); 
# 12922 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaGraphAddNode(cudaGraphNode_t * pGraphNode, cudaGraph_t graph, const cudaGraphNode_t * pDependencies, size_t numDependencies, cudaGraphNodeParams * nodeParams); 
# 12951 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaGraphNodeSetParams(cudaGraphNode_t node, cudaGraphNodeParams * nodeParams); 
# 13000 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaGraphExecNodeSetParams(cudaGraphExec_t graphExec, cudaGraphNode_t node, cudaGraphNodeParams * nodeParams); 
# 13077 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaGetDriverEntryPoint(const char * symbol, void ** funcPtr, unsigned long long flags, cudaDriverEntryPointQueryResult * driverStatus = 0); 
# 13085 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaGetExportTable(const void ** ppExportTable, const cudaUUID_t * pExportTableId); 
# 13264 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaGetFuncBySymbol(cudaFunction_t * functionPtr, const void * symbolPtr); 
# 13280 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
extern cudaError_t cudaGetKernel(cudaKernel_t * kernelPtr, const void * entryFuncAddr); 
# 13443 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime_api.h"
}
# 117 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/channel_descriptor.h"
template< class T> inline cudaChannelFormatDesc cudaCreateChannelDesc() 
# 118
{ 
# 119
return cudaCreateChannelDesc(0, 0, 0, 0, cudaChannelFormatKindNone); 
# 120
} 
# 122
static inline cudaChannelFormatDesc cudaCreateChannelDescHalf() 
# 123
{ 
# 124
int e = (((int)sizeof(unsigned short)) * 8); 
# 126
return cudaCreateChannelDesc(e, 0, 0, 0, cudaChannelFormatKindFloat); 
# 127
} 
# 129
static inline cudaChannelFormatDesc cudaCreateChannelDescHalf1() 
# 130
{ 
# 131
int e = (((int)sizeof(unsigned short)) * 8); 
# 133
return cudaCreateChannelDesc(e, 0, 0, 0, cudaChannelFormatKindFloat); 
# 134
} 
# 136
static inline cudaChannelFormatDesc cudaCreateChannelDescHalf2() 
# 137
{ 
# 138
int e = (((int)sizeof(unsigned short)) * 8); 
# 140
return cudaCreateChannelDesc(e, e, 0, 0, cudaChannelFormatKindFloat); 
# 141
} 
# 143
static inline cudaChannelFormatDesc cudaCreateChannelDescHalf4() 
# 144
{ 
# 145
int e = (((int)sizeof(unsigned short)) * 8); 
# 147
return cudaCreateChannelDesc(e, e, e, e, cudaChannelFormatKindFloat); 
# 148
} 
# 150
template<> inline cudaChannelFormatDesc cudaCreateChannelDesc< char> () 
# 151
{ 
# 152
int e = (((int)sizeof(char)) * 8); 
# 157
return cudaCreateChannelDesc(e, 0, 0, 0, cudaChannelFormatKindSigned); 
# 159
} 
# 161
template<> inline cudaChannelFormatDesc cudaCreateChannelDesc< signed char> () 
# 162
{ 
# 163
int e = (((int)sizeof(signed char)) * 8); 
# 165
return cudaCreateChannelDesc(e, 0, 0, 0, cudaChannelFormatKindSigned); 
# 166
} 
# 168
template<> inline cudaChannelFormatDesc cudaCreateChannelDesc< unsigned char> () 
# 169
{ 
# 170
int e = (((int)sizeof(unsigned char)) * 8); 
# 172
return cudaCreateChannelDesc(e, 0, 0, 0, cudaChannelFormatKindUnsigned); 
# 173
} 
# 175
template<> inline cudaChannelFormatDesc cudaCreateChannelDesc< char1> () 
# 176
{ 
# 177
int e = (((int)sizeof(signed char)) * 8); 
# 179
return cudaCreateChannelDesc(e, 0, 0, 0, cudaChannelFormatKindSigned); 
# 180
} 
# 182
template<> inline cudaChannelFormatDesc cudaCreateChannelDesc< uchar1> () 
# 183
{ 
# 184
int e = (((int)sizeof(unsigned char)) * 8); 
# 186
return cudaCreateChannelDesc(e, 0, 0, 0, cudaChannelFormatKindUnsigned); 
# 187
} 
# 189
template<> inline cudaChannelFormatDesc cudaCreateChannelDesc< char2> () 
# 190
{ 
# 191
int e = (((int)sizeof(signed char)) * 8); 
# 193
return cudaCreateChannelDesc(e, e, 0, 0, cudaChannelFormatKindSigned); 
# 194
} 
# 196
template<> inline cudaChannelFormatDesc cudaCreateChannelDesc< uchar2> () 
# 197
{ 
# 198
int e = (((int)sizeof(unsigned char)) * 8); 
# 200
return cudaCreateChannelDesc(e, e, 0, 0, cudaChannelFormatKindUnsigned); 
# 201
} 
# 203
template<> inline cudaChannelFormatDesc cudaCreateChannelDesc< char4> () 
# 204
{ 
# 205
int e = (((int)sizeof(signed char)) * 8); 
# 207
return cudaCreateChannelDesc(e, e, e, e, cudaChannelFormatKindSigned); 
# 208
} 
# 210
template<> inline cudaChannelFormatDesc cudaCreateChannelDesc< uchar4> () 
# 211
{ 
# 212
int e = (((int)sizeof(unsigned char)) * 8); 
# 214
return cudaCreateChannelDesc(e, e, e, e, cudaChannelFormatKindUnsigned); 
# 215
} 
# 217
template<> inline cudaChannelFormatDesc cudaCreateChannelDesc< short> () 
# 218
{ 
# 219
int e = (((int)sizeof(short)) * 8); 
# 221
return cudaCreateChannelDesc(e, 0, 0, 0, cudaChannelFormatKindSigned); 
# 222
} 
# 224
template<> inline cudaChannelFormatDesc cudaCreateChannelDesc< unsigned short> () 
# 225
{ 
# 226
int e = (((int)sizeof(unsigned short)) * 8); 
# 228
return cudaCreateChannelDesc(e, 0, 0, 0, cudaChannelFormatKindUnsigned); 
# 229
} 
# 231
template<> inline cudaChannelFormatDesc cudaCreateChannelDesc< short1> () 
# 232
{ 
# 233
int e = (((int)sizeof(short)) * 8); 
# 235
return cudaCreateChannelDesc(e, 0, 0, 0, cudaChannelFormatKindSigned); 
# 236
} 
# 238
template<> inline cudaChannelFormatDesc cudaCreateChannelDesc< ushort1> () 
# 239
{ 
# 240
int e = (((int)sizeof(unsigned short)) * 8); 
# 242
return cudaCreateChannelDesc(e, 0, 0, 0, cudaChannelFormatKindUnsigned); 
# 243
} 
# 245
template<> inline cudaChannelFormatDesc cudaCreateChannelDesc< short2> () 
# 246
{ 
# 247
int e = (((int)sizeof(short)) * 8); 
# 249
return cudaCreateChannelDesc(e, e, 0, 0, cudaChannelFormatKindSigned); 
# 250
} 
# 252
template<> inline cudaChannelFormatDesc cudaCreateChannelDesc< ushort2> () 
# 253
{ 
# 254
int e = (((int)sizeof(unsigned short)) * 8); 
# 256
return cudaCreateChannelDesc(e, e, 0, 0, cudaChannelFormatKindUnsigned); 
# 257
} 
# 259
template<> inline cudaChannelFormatDesc cudaCreateChannelDesc< short4> () 
# 260
{ 
# 261
int e = (((int)sizeof(short)) * 8); 
# 263
return cudaCreateChannelDesc(e, e, e, e, cudaChannelFormatKindSigned); 
# 264
} 
# 266
template<> inline cudaChannelFormatDesc cudaCreateChannelDesc< ushort4> () 
# 267
{ 
# 268
int e = (((int)sizeof(unsigned short)) * 8); 
# 270
return cudaCreateChannelDesc(e, e, e, e, cudaChannelFormatKindUnsigned); 
# 271
} 
# 273
template<> inline cudaChannelFormatDesc cudaCreateChannelDesc< int> () 
# 274
{ 
# 275
int e = (((int)sizeof(int)) * 8); 
# 277
return cudaCreateChannelDesc(e, 0, 0, 0, cudaChannelFormatKindSigned); 
# 278
} 
# 280
template<> inline cudaChannelFormatDesc cudaCreateChannelDesc< unsigned> () 
# 281
{ 
# 282
int e = (((int)sizeof(unsigned)) * 8); 
# 284
return cudaCreateChannelDesc(e, 0, 0, 0, cudaChannelFormatKindUnsigned); 
# 285
} 
# 287
template<> inline cudaChannelFormatDesc cudaCreateChannelDesc< int1> () 
# 288
{ 
# 289
int e = (((int)sizeof(int)) * 8); 
# 291
return cudaCreateChannelDesc(e, 0, 0, 0, cudaChannelFormatKindSigned); 
# 292
} 
# 294
template<> inline cudaChannelFormatDesc cudaCreateChannelDesc< uint1> () 
# 295
{ 
# 296
int e = (((int)sizeof(unsigned)) * 8); 
# 298
return cudaCreateChannelDesc(e, 0, 0, 0, cudaChannelFormatKindUnsigned); 
# 299
} 
# 301
template<> inline cudaChannelFormatDesc cudaCreateChannelDesc< int2> () 
# 302
{ 
# 303
int e = (((int)sizeof(int)) * 8); 
# 305
return cudaCreateChannelDesc(e, e, 0, 0, cudaChannelFormatKindSigned); 
# 306
} 
# 308
template<> inline cudaChannelFormatDesc cudaCreateChannelDesc< uint2> () 
# 309
{ 
# 310
int e = (((int)sizeof(unsigned)) * 8); 
# 312
return cudaCreateChannelDesc(e, e, 0, 0, cudaChannelFormatKindUnsigned); 
# 313
} 
# 315
template<> inline cudaChannelFormatDesc cudaCreateChannelDesc< int4> () 
# 316
{ 
# 317
int e = (((int)sizeof(int)) * 8); 
# 319
return cudaCreateChannelDesc(e, e, e, e, cudaChannelFormatKindSigned); 
# 320
} 
# 322
template<> inline cudaChannelFormatDesc cudaCreateChannelDesc< uint4> () 
# 323
{ 
# 324
int e = (((int)sizeof(unsigned)) * 8); 
# 326
return cudaCreateChannelDesc(e, e, e, e, cudaChannelFormatKindUnsigned); 
# 327
} 
# 389 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/channel_descriptor.h"
template<> inline cudaChannelFormatDesc cudaCreateChannelDesc< float> () 
# 390
{ 
# 391
int e = (((int)sizeof(float)) * 8); 
# 393
return cudaCreateChannelDesc(e, 0, 0, 0, cudaChannelFormatKindFloat); 
# 394
} 
# 396
template<> inline cudaChannelFormatDesc cudaCreateChannelDesc< float1> () 
# 397
{ 
# 398
int e = (((int)sizeof(float)) * 8); 
# 400
return cudaCreateChannelDesc(e, 0, 0, 0, cudaChannelFormatKindFloat); 
# 401
} 
# 403
template<> inline cudaChannelFormatDesc cudaCreateChannelDesc< float2> () 
# 404
{ 
# 405
int e = (((int)sizeof(float)) * 8); 
# 407
return cudaCreateChannelDesc(e, e, 0, 0, cudaChannelFormatKindFloat); 
# 408
} 
# 410
template<> inline cudaChannelFormatDesc cudaCreateChannelDesc< float4> () 
# 411
{ 
# 412
int e = (((int)sizeof(float)) * 8); 
# 414
return cudaCreateChannelDesc(e, e, e, e, cudaChannelFormatKindFloat); 
# 415
} 
# 417
static inline cudaChannelFormatDesc cudaCreateChannelDescNV12() 
# 418
{ 
# 419
int e = (((int)sizeof(char)) * 8); 
# 421
return cudaCreateChannelDesc(e, e, e, 0, cudaChannelFormatKindNV12); 
# 422
} 
# 424
template< cudaChannelFormatKind > inline cudaChannelFormatDesc cudaCreateChannelDesc() 
# 425
{ 
# 426
return cudaCreateChannelDesc(0, 0, 0, 0, cudaChannelFormatKindNone); 
# 427
} 
# 430
template<> inline cudaChannelFormatDesc cudaCreateChannelDesc< cudaChannelFormatKindSignedNormalized8X1> () 
# 431
{ 
# 432
return cudaCreateChannelDesc(8, 0, 0, 0, cudaChannelFormatKindSignedNormalized8X1); 
# 433
} 
# 435
template<> inline cudaChannelFormatDesc cudaCreateChannelDesc< cudaChannelFormatKindSignedNormalized8X2> () 
# 436
{ 
# 437
return cudaCreateChannelDesc(8, 8, 0, 0, cudaChannelFormatKindSignedNormalized8X2); 
# 438
} 
# 440
template<> inline cudaChannelFormatDesc cudaCreateChannelDesc< cudaChannelFormatKindSignedNormalized8X4> () 
# 441
{ 
# 442
return cudaCreateChannelDesc(8, 8, 8, 8, cudaChannelFormatKindSignedNormalized8X4); 
# 443
} 
# 446
template<> inline cudaChannelFormatDesc cudaCreateChannelDesc< cudaChannelFormatKindUnsignedNormalized8X1> () 
# 447
{ 
# 448
return cudaCreateChannelDesc(8, 0, 0, 0, cudaChannelFormatKindUnsignedNormalized8X1); 
# 449
} 
# 451
template<> inline cudaChannelFormatDesc cudaCreateChannelDesc< cudaChannelFormatKindUnsignedNormalized8X2> () 
# 452
{ 
# 453
return cudaCreateChannelDesc(8, 8, 0, 0, cudaChannelFormatKindUnsignedNormalized8X2); 
# 454
} 
# 456
template<> inline cudaChannelFormatDesc cudaCreateChannelDesc< cudaChannelFormatKindUnsignedNormalized8X4> () 
# 457
{ 
# 458
return cudaCreateChannelDesc(8, 8, 8, 8, cudaChannelFormatKindUnsignedNormalized8X4); 
# 459
} 
# 462
template<> inline cudaChannelFormatDesc cudaCreateChannelDesc< cudaChannelFormatKindSignedNormalized16X1> () 
# 463
{ 
# 464
return cudaCreateChannelDesc(16, 0, 0, 0, cudaChannelFormatKindSignedNormalized16X1); 
# 465
} 
# 467
template<> inline cudaChannelFormatDesc cudaCreateChannelDesc< cudaChannelFormatKindSignedNormalized16X2> () 
# 468
{ 
# 469
return cudaCreateChannelDesc(16, 16, 0, 0, cudaChannelFormatKindSignedNormalized16X2); 
# 470
} 
# 472
template<> inline cudaChannelFormatDesc cudaCreateChannelDesc< cudaChannelFormatKindSignedNormalized16X4> () 
# 473
{ 
# 474
return cudaCreateChannelDesc(16, 16, 16, 16, cudaChannelFormatKindSignedNormalized16X4); 
# 475
} 
# 478
template<> inline cudaChannelFormatDesc cudaCreateChannelDesc< cudaChannelFormatKindUnsignedNormalized16X1> () 
# 479
{ 
# 480
return cudaCreateChannelDesc(16, 0, 0, 0, cudaChannelFormatKindUnsignedNormalized16X1); 
# 481
} 
# 483
template<> inline cudaChannelFormatDesc cudaCreateChannelDesc< cudaChannelFormatKindUnsignedNormalized16X2> () 
# 484
{ 
# 485
return cudaCreateChannelDesc(16, 16, 0, 0, cudaChannelFormatKindUnsignedNormalized16X2); 
# 486
} 
# 488
template<> inline cudaChannelFormatDesc cudaCreateChannelDesc< cudaChannelFormatKindUnsignedNormalized16X4> () 
# 489
{ 
# 490
return cudaCreateChannelDesc(16, 16, 16, 16, cudaChannelFormatKindUnsignedNormalized16X4); 
# 491
} 
# 494
template<> inline cudaChannelFormatDesc cudaCreateChannelDesc< cudaChannelFormatKindNV12> () 
# 495
{ 
# 496
return cudaCreateChannelDesc(8, 8, 8, 0, cudaChannelFormatKindNV12); 
# 497
} 
# 500
template<> inline cudaChannelFormatDesc cudaCreateChannelDesc< cudaChannelFormatKindUnsignedBlockCompressed1> () 
# 501
{ 
# 502
return cudaCreateChannelDesc(8, 8, 8, 8, cudaChannelFormatKindUnsignedBlockCompressed1); 
# 503
} 
# 506
template<> inline cudaChannelFormatDesc cudaCreateChannelDesc< cudaChannelFormatKindUnsignedBlockCompressed1SRGB> () 
# 507
{ 
# 508
return cudaCreateChannelDesc(8, 8, 8, 8, cudaChannelFormatKindUnsignedBlockCompressed1SRGB); 
# 509
} 
# 512
template<> inline cudaChannelFormatDesc cudaCreateChannelDesc< cudaChannelFormatKindUnsignedBlockCompressed2> () 
# 513
{ 
# 514
return cudaCreateChannelDesc(8, 8, 8, 8, cudaChannelFormatKindUnsignedBlockCompressed2); 
# 515
} 
# 518
template<> inline cudaChannelFormatDesc cudaCreateChannelDesc< cudaChannelFormatKindUnsignedBlockCompressed2SRGB> () 
# 519
{ 
# 520
return cudaCreateChannelDesc(8, 8, 8, 8, cudaChannelFormatKindUnsignedBlockCompressed2SRGB); 
# 521
} 
# 524
template<> inline cudaChannelFormatDesc cudaCreateChannelDesc< cudaChannelFormatKindUnsignedBlockCompressed3> () 
# 525
{ 
# 526
return cudaCreateChannelDesc(8, 8, 8, 8, cudaChannelFormatKindUnsignedBlockCompressed3); 
# 527
} 
# 530
template<> inline cudaChannelFormatDesc cudaCreateChannelDesc< cudaChannelFormatKindUnsignedBlockCompressed3SRGB> () 
# 531
{ 
# 532
return cudaCreateChannelDesc(8, 8, 8, 8, cudaChannelFormatKindUnsignedBlockCompressed3SRGB); 
# 533
} 
# 536
template<> inline cudaChannelFormatDesc cudaCreateChannelDesc< cudaChannelFormatKindUnsignedBlockCompressed4> () 
# 537
{ 
# 538
return cudaCreateChannelDesc(8, 0, 0, 0, cudaChannelFormatKindUnsignedBlockCompressed4); 
# 539
} 
# 542
template<> inline cudaChannelFormatDesc cudaCreateChannelDesc< cudaChannelFormatKindSignedBlockCompressed4> () 
# 543
{ 
# 544
return cudaCreateChannelDesc(8, 0, 0, 0, cudaChannelFormatKindSignedBlockCompressed4); 
# 545
} 
# 548
template<> inline cudaChannelFormatDesc cudaCreateChannelDesc< cudaChannelFormatKindUnsignedBlockCompressed5> () 
# 549
{ 
# 550
return cudaCreateChannelDesc(8, 8, 0, 0, cudaChannelFormatKindUnsignedBlockCompressed5); 
# 551
} 
# 554
template<> inline cudaChannelFormatDesc cudaCreateChannelDesc< cudaChannelFormatKindSignedBlockCompressed5> () 
# 555
{ 
# 556
return cudaCreateChannelDesc(8, 8, 0, 0, cudaChannelFormatKindSignedBlockCompressed5); 
# 557
} 
# 560
template<> inline cudaChannelFormatDesc cudaCreateChannelDesc< cudaChannelFormatKindUnsignedBlockCompressed6H> () 
# 561
{ 
# 562
return cudaCreateChannelDesc(16, 16, 16, 0, cudaChannelFormatKindUnsignedBlockCompressed6H); 
# 563
} 
# 566
template<> inline cudaChannelFormatDesc cudaCreateChannelDesc< cudaChannelFormatKindSignedBlockCompressed6H> () 
# 567
{ 
# 568
return cudaCreateChannelDesc(16, 16, 16, 0, cudaChannelFormatKindSignedBlockCompressed6H); 
# 569
} 
# 572
template<> inline cudaChannelFormatDesc cudaCreateChannelDesc< cudaChannelFormatKindUnsignedBlockCompressed7> () 
# 573
{ 
# 574
return cudaCreateChannelDesc(8, 8, 8, 8, cudaChannelFormatKindUnsignedBlockCompressed7); 
# 575
} 
# 578
template<> inline cudaChannelFormatDesc cudaCreateChannelDesc< cudaChannelFormatKindUnsignedBlockCompressed7SRGB> () 
# 579
{ 
# 580
return cudaCreateChannelDesc(8, 8, 8, 8, cudaChannelFormatKindUnsignedBlockCompressed7SRGB); 
# 581
} 
# 79 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_functions.h"
static inline cudaPitchedPtr make_cudaPitchedPtr(void *d, size_t p, size_t xsz, size_t ysz) 
# 80
{ 
# 81
cudaPitchedPtr s; 
# 83
(s.ptr) = d; 
# 84
(s.pitch) = p; 
# 85
(s.xsize) = xsz; 
# 86
(s.ysize) = ysz; 
# 88
return s; 
# 89
} 
# 106 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_functions.h"
static inline cudaPos make_cudaPos(size_t x, size_t y, size_t z) 
# 107
{ 
# 108
cudaPos p; 
# 110
(p.x) = x; 
# 111
(p.y) = y; 
# 112
(p.z) = z; 
# 114
return p; 
# 115
} 
# 132 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/driver_functions.h"
static inline cudaExtent make_cudaExtent(size_t w, size_t h, size_t d) 
# 133
{ 
# 134
cudaExtent e; 
# 136
(e.width) = w; 
# 137
(e.height) = h; 
# 138
(e.depth) = d; 
# 140
return e; 
# 141
} 
# 73 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/vector_functions.h"
static inline char1 make_char1(signed char x); 
# 75
static inline uchar1 make_uchar1(unsigned char x); 
# 77
static inline char2 make_char2(signed char x, signed char y); 
# 79
static inline uchar2 make_uchar2(unsigned char x, unsigned char y); 
# 81
static inline char3 make_char3(signed char x, signed char y, signed char z); 
# 83
static inline uchar3 make_uchar3(unsigned char x, unsigned char y, unsigned char z); 
# 85
static inline char4 make_char4(signed char x, signed char y, signed char z, signed char w); 
# 87
static inline uchar4 make_uchar4(unsigned char x, unsigned char y, unsigned char z, unsigned char w); 
# 89
static inline short1 make_short1(short x); 
# 91
static inline ushort1 make_ushort1(unsigned short x); 
# 93
static inline short2 make_short2(short x, short y); 
# 95
static inline ushort2 make_ushort2(unsigned short x, unsigned short y); 
# 97
static inline short3 make_short3(short x, short y, short z); 
# 99
static inline ushort3 make_ushort3(unsigned short x, unsigned short y, unsigned short z); 
# 101
static inline short4 make_short4(short x, short y, short z, short w); 
# 103
static inline ushort4 make_ushort4(unsigned short x, unsigned short y, unsigned short z, unsigned short w); 
# 105
static inline int1 make_int1(int x); 
# 107
static inline uint1 make_uint1(unsigned x); 
# 109
static inline int2 make_int2(int x, int y); 
# 111
static inline uint2 make_uint2(unsigned x, unsigned y); 
# 113
static inline int3 make_int3(int x, int y, int z); 
# 115
static inline uint3 make_uint3(unsigned x, unsigned y, unsigned z); 
# 117
static inline int4 make_int4(int x, int y, int z, int w); 
# 119
static inline uint4 make_uint4(unsigned x, unsigned y, unsigned z, unsigned w); 
# 121
static inline long1 make_long1(long x); 
# 123
static inline ulong1 make_ulong1(unsigned long x); 
# 125
static inline long2 make_long2(long x, long y); 
# 127
static inline ulong2 make_ulong2(unsigned long x, unsigned long y); 
# 129
static inline long3 make_long3(long x, long y, long z); 
# 131
static inline ulong3 make_ulong3(unsigned long x, unsigned long y, unsigned long z); 
# 133
static inline long4 make_long4(long x, long y, long z, long w); 
# 135
static inline ulong4 make_ulong4(unsigned long x, unsigned long y, unsigned long z, unsigned long w); 
# 137
static inline float1 make_float1(float x); 
# 139
static inline float2 make_float2(float x, float y); 
# 141
static inline float3 make_float3(float x, float y, float z); 
# 143
static inline float4 make_float4(float x, float y, float z, float w); 
# 145
static inline longlong1 make_longlong1(long long x); 
# 147
static inline ulonglong1 make_ulonglong1(unsigned long long x); 
# 149
static inline longlong2 make_longlong2(long long x, long long y); 
# 151
static inline ulonglong2 make_ulonglong2(unsigned long long x, unsigned long long y); 
# 153
static inline longlong3 make_longlong3(long long x, long long y, long long z); 
# 155
static inline ulonglong3 make_ulonglong3(unsigned long long x, unsigned long long y, unsigned long long z); 
# 157
static inline longlong4 make_longlong4(long long x, long long y, long long z, long long w); 
# 159
static inline ulonglong4 make_ulonglong4(unsigned long long x, unsigned long long y, unsigned long long z, unsigned long long w); 
# 161
static inline double1 make_double1(double x); 
# 163
static inline double2 make_double2(double x, double y); 
# 165
static inline double3 make_double3(double x, double y, double z); 
# 167
static inline double4 make_double4(double x, double y, double z, double w); 
# 73 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/vector_functions.hpp"
static inline char1 make_char1(signed char x) 
# 74
{ 
# 75
char1 t; (t.x) = x; return t; 
# 76
} 
# 78
static inline uchar1 make_uchar1(unsigned char x) 
# 79
{ 
# 80
uchar1 t; (t.x) = x; return t; 
# 81
} 
# 83
static inline char2 make_char2(signed char x, signed char y) 
# 84
{ 
# 85
char2 t; (t.x) = x; (t.y) = y; return t; 
# 86
} 
# 88
static inline uchar2 make_uchar2(unsigned char x, unsigned char y) 
# 89
{ 
# 90
uchar2 t; (t.x) = x; (t.y) = y; return t; 
# 91
} 
# 93
static inline char3 make_char3(signed char x, signed char y, signed char z) 
# 94
{ 
# 95
char3 t; (t.x) = x; (t.y) = y; (t.z) = z; return t; 
# 96
} 
# 98
static inline uchar3 make_uchar3(unsigned char x, unsigned char y, unsigned char z) 
# 99
{ 
# 100
uchar3 t; (t.x) = x; (t.y) = y; (t.z) = z; return t; 
# 101
} 
# 103
static inline char4 make_char4(signed char x, signed char y, signed char z, signed char w) 
# 104
{ 
# 105
char4 t; (t.x) = x; (t.y) = y; (t.z) = z; (t.w) = w; return t; 
# 106
} 
# 108
static inline uchar4 make_uchar4(unsigned char x, unsigned char y, unsigned char z, unsigned char w) 
# 109
{ 
# 110
uchar4 t; (t.x) = x; (t.y) = y; (t.z) = z; (t.w) = w; return t; 
# 111
} 
# 113
static inline short1 make_short1(short x) 
# 114
{ 
# 115
short1 t; (t.x) = x; return t; 
# 116
} 
# 118
static inline ushort1 make_ushort1(unsigned short x) 
# 119
{ 
# 120
ushort1 t; (t.x) = x; return t; 
# 121
} 
# 123
static inline short2 make_short2(short x, short y) 
# 124
{ 
# 125
short2 t; (t.x) = x; (t.y) = y; return t; 
# 126
} 
# 128
static inline ushort2 make_ushort2(unsigned short x, unsigned short y) 
# 129
{ 
# 130
ushort2 t; (t.x) = x; (t.y) = y; return t; 
# 131
} 
# 133
static inline short3 make_short3(short x, short y, short z) 
# 134
{ 
# 135
short3 t; (t.x) = x; (t.y) = y; (t.z) = z; return t; 
# 136
} 
# 138
static inline ushort3 make_ushort3(unsigned short x, unsigned short y, unsigned short z) 
# 139
{ 
# 140
ushort3 t; (t.x) = x; (t.y) = y; (t.z) = z; return t; 
# 141
} 
# 143
static inline short4 make_short4(short x, short y, short z, short w) 
# 144
{ 
# 145
short4 t; (t.x) = x; (t.y) = y; (t.z) = z; (t.w) = w; return t; 
# 146
} 
# 148
static inline ushort4 make_ushort4(unsigned short x, unsigned short y, unsigned short z, unsigned short w) 
# 149
{ 
# 150
ushort4 t; (t.x) = x; (t.y) = y; (t.z) = z; (t.w) = w; return t; 
# 151
} 
# 153
static inline int1 make_int1(int x) 
# 154
{ 
# 155
int1 t; (t.x) = x; return t; 
# 156
} 
# 158
static inline uint1 make_uint1(unsigned x) 
# 159
{ 
# 160
uint1 t; (t.x) = x; return t; 
# 161
} 
# 163
static inline int2 make_int2(int x, int y) 
# 164
{ 
# 165
int2 t; (t.x) = x; (t.y) = y; return t; 
# 166
} 
# 168
static inline uint2 make_uint2(unsigned x, unsigned y) 
# 169
{ 
# 170
uint2 t; (t.x) = x; (t.y) = y; return t; 
# 171
} 
# 173
static inline int3 make_int3(int x, int y, int z) 
# 174
{ 
# 175
int3 t; (t.x) = x; (t.y) = y; (t.z) = z; return t; 
# 176
} 
# 178
static inline uint3 make_uint3(unsigned x, unsigned y, unsigned z) 
# 179
{ 
# 180
uint3 t; (t.x) = x; (t.y) = y; (t.z) = z; return t; 
# 181
} 
# 183
static inline int4 make_int4(int x, int y, int z, int w) 
# 184
{ 
# 185
int4 t; (t.x) = x; (t.y) = y; (t.z) = z; (t.w) = w; return t; 
# 186
} 
# 188
static inline uint4 make_uint4(unsigned x, unsigned y, unsigned z, unsigned w) 
# 189
{ 
# 190
uint4 t; (t.x) = x; (t.y) = y; (t.z) = z; (t.w) = w; return t; 
# 191
} 
# 193
static inline long1 make_long1(long x) 
# 194
{ 
# 195
long1 t; (t.x) = x; return t; 
# 196
} 
# 198
static inline ulong1 make_ulong1(unsigned long x) 
# 199
{ 
# 200
ulong1 t; (t.x) = x; return t; 
# 201
} 
# 203
static inline long2 make_long2(long x, long y) 
# 204
{ 
# 205
long2 t; (t.x) = x; (t.y) = y; return t; 
# 206
} 
# 208
static inline ulong2 make_ulong2(unsigned long x, unsigned long y) 
# 209
{ 
# 210
ulong2 t; (t.x) = x; (t.y) = y; return t; 
# 211
} 
# 213
static inline long3 make_long3(long x, long y, long z) 
# 214
{ 
# 215
long3 t; (t.x) = x; (t.y) = y; (t.z) = z; return t; 
# 216
} 
# 218
static inline ulong3 make_ulong3(unsigned long x, unsigned long y, unsigned long z) 
# 219
{ 
# 220
ulong3 t; (t.x) = x; (t.y) = y; (t.z) = z; return t; 
# 221
} 
# 223
static inline long4 make_long4(long x, long y, long z, long w) 
# 224
{ 
# 225
long4 t; (t.x) = x; (t.y) = y; (t.z) = z; (t.w) = w; return t; 
# 226
} 
# 228
static inline ulong4 make_ulong4(unsigned long x, unsigned long y, unsigned long z, unsigned long w) 
# 229
{ 
# 230
ulong4 t; (t.x) = x; (t.y) = y; (t.z) = z; (t.w) = w; return t; 
# 231
} 
# 233
static inline float1 make_float1(float x) 
# 234
{ 
# 235
float1 t; (t.x) = x; return t; 
# 236
} 
# 238
static inline float2 make_float2(float x, float y) 
# 239
{ 
# 240
float2 t; (t.x) = x; (t.y) = y; return t; 
# 241
} 
# 243
static inline float3 make_float3(float x, float y, float z) 
# 244
{ 
# 245
float3 t; (t.x) = x; (t.y) = y; (t.z) = z; return t; 
# 246
} 
# 248
static inline float4 make_float4(float x, float y, float z, float w) 
# 249
{ 
# 250
float4 t; (t.x) = x; (t.y) = y; (t.z) = z; (t.w) = w; return t; 
# 251
} 
# 253
static inline longlong1 make_longlong1(long long x) 
# 254
{ 
# 255
longlong1 t; (t.x) = x; return t; 
# 256
} 
# 258
static inline ulonglong1 make_ulonglong1(unsigned long long x) 
# 259
{ 
# 260
ulonglong1 t; (t.x) = x; return t; 
# 261
} 
# 263
static inline longlong2 make_longlong2(long long x, long long y) 
# 264
{ 
# 265
longlong2 t; (t.x) = x; (t.y) = y; return t; 
# 266
} 
# 268
static inline ulonglong2 make_ulonglong2(unsigned long long x, unsigned long long y) 
# 269
{ 
# 270
ulonglong2 t; (t.x) = x; (t.y) = y; return t; 
# 271
} 
# 273
static inline longlong3 make_longlong3(long long x, long long y, long long z) 
# 274
{ 
# 275
longlong3 t; (t.x) = x; (t.y) = y; (t.z) = z; return t; 
# 276
} 
# 278
static inline ulonglong3 make_ulonglong3(unsigned long long x, unsigned long long y, unsigned long long z) 
# 279
{ 
# 280
ulonglong3 t; (t.x) = x; (t.y) = y; (t.z) = z; return t; 
# 281
} 
# 283
static inline longlong4 make_longlong4(long long x, long long y, long long z, long long w) 
# 284
{ 
# 285
longlong4 t; (t.x) = x; (t.y) = y; (t.z) = z; (t.w) = w; return t; 
# 286
} 
# 288
static inline ulonglong4 make_ulonglong4(unsigned long long x, unsigned long long y, unsigned long long z, unsigned long long w) 
# 289
{ 
# 290
ulonglong4 t; (t.x) = x; (t.y) = y; (t.z) = z; (t.w) = w; return t; 
# 291
} 
# 293
static inline double1 make_double1(double x) 
# 294
{ 
# 295
double1 t; (t.x) = x; return t; 
# 296
} 
# 298
static inline double2 make_double2(double x, double y) 
# 299
{ 
# 300
double2 t; (t.x) = x; (t.y) = y; return t; 
# 301
} 
# 303
static inline double3 make_double3(double x, double y, double z) 
# 304
{ 
# 305
double3 t; (t.x) = x; (t.y) = y; (t.z) = z; return t; 
# 306
} 
# 308
static inline double4 make_double4(double x, double y, double z, double w) 
# 309
{ 
# 310
double4 t; (t.x) = x; (t.y) = y; (t.z) = z; (t.w) = w; return t; 
# 311
} 
# 28 "/usr/include/string.h" 3
extern "C" {
# 43 "/usr/include/string.h" 3
extern void *memcpy(void *__restrict__ __dest, const void *__restrict__ __src, size_t __n) throw()
# 44
 __attribute((__nonnull__(1, 2))); 
# 47
extern void *memmove(void * __dest, const void * __src, size_t __n) throw()
# 48
 __attribute((__nonnull__(1, 2))); 
# 54
extern void *memccpy(void *__restrict__ __dest, const void *__restrict__ __src, int __c, size_t __n) throw()
# 56
 __attribute((__nonnull__(1, 2))); 
# 61
extern void *memset(void * __s, int __c, size_t __n) throw() __attribute((__nonnull__(1))); 
# 64
extern int memcmp(const void * __s1, const void * __s2, size_t __n) throw()
# 65
 __attribute((__pure__)) __attribute((__nonnull__(1, 2))); 
# 69
extern "C++" {
# 71
extern void *memchr(void * __s, int __c, size_t __n) throw() __asm__("memchr")
# 72
 __attribute((__pure__)) __attribute((__nonnull__(1))); 
# 73
extern const void *memchr(const void * __s, int __c, size_t __n) throw() __asm__("memchr")
# 74
 __attribute((__pure__)) __attribute((__nonnull__(1))); 
# 89 "/usr/include/string.h" 3
}
# 99 "/usr/include/string.h" 3
extern "C++" void *rawmemchr(void * __s, int __c) throw() __asm__("rawmemchr")
# 100
 __attribute((__pure__)) __attribute((__nonnull__(1))); 
# 101
extern "C++" const void *rawmemchr(const void * __s, int __c) throw() __asm__("rawmemchr")
# 102
 __attribute((__pure__)) __attribute((__nonnull__(1))); 
# 110
extern "C++" void *memrchr(void * __s, int __c, size_t __n) throw() __asm__("memrchr")
# 111
 __attribute((__pure__)) __attribute((__nonnull__(1))); 
# 112
extern "C++" const void *memrchr(const void * __s, int __c, size_t __n) throw() __asm__("memrchr")
# 113
 __attribute((__pure__)) __attribute((__nonnull__(1))); 
# 122 "/usr/include/string.h" 3
extern char *strcpy(char *__restrict__ __dest, const char *__restrict__ __src) throw()
# 123
 __attribute((__nonnull__(1, 2))); 
# 125
extern char *strncpy(char *__restrict__ __dest, const char *__restrict__ __src, size_t __n) throw()
# 127
 __attribute((__nonnull__(1, 2))); 
# 130
extern char *strcat(char *__restrict__ __dest, const char *__restrict__ __src) throw()
# 131
 __attribute((__nonnull__(1, 2))); 
# 133
extern char *strncat(char *__restrict__ __dest, const char *__restrict__ __src, size_t __n) throw()
# 134
 __attribute((__nonnull__(1, 2))); 
# 137
extern int strcmp(const char * __s1, const char * __s2) throw()
# 138
 __attribute((__pure__)) __attribute((__nonnull__(1, 2))); 
# 140
extern int strncmp(const char * __s1, const char * __s2, size_t __n) throw()
# 141
 __attribute((__pure__)) __attribute((__nonnull__(1, 2))); 
# 144
extern int strcoll(const char * __s1, const char * __s2) throw()
# 145
 __attribute((__pure__)) __attribute((__nonnull__(1, 2))); 
# 147
extern size_t strxfrm(char *__restrict__ __dest, const char *__restrict__ __src, size_t __n) throw()
# 149
 __attribute((__nonnull__(2))); 
# 156
extern int strcoll_l(const char * __s1, const char * __s2, locale_t __l) throw()
# 157
 __attribute((__pure__)) __attribute((__nonnull__(1, 2, 3))); 
# 160
extern size_t strxfrm_l(char * __dest, const char * __src, size_t __n, locale_t __l) throw()
# 161
 __attribute((__nonnull__(2, 4))); 
# 167
extern char *strdup(const char * __s) throw()
# 168
 __attribute((__malloc__)) __attribute((__nonnull__(1))); 
# 175
extern char *strndup(const char * __string, size_t __n) throw()
# 176
 __attribute((__malloc__)) __attribute((__nonnull__(1))); 
# 204 "/usr/include/string.h" 3
extern "C++" {
# 206
extern char *strchr(char * __s, int __c) throw() __asm__("strchr")
# 207
 __attribute((__pure__)) __attribute((__nonnull__(1))); 
# 208
extern const char *strchr(const char * __s, int __c) throw() __asm__("strchr")
# 209
 __attribute((__pure__)) __attribute((__nonnull__(1))); 
# 224 "/usr/include/string.h" 3
}
# 231
extern "C++" {
# 233
extern char *strrchr(char * __s, int __c) throw() __asm__("strrchr")
# 234
 __attribute((__pure__)) __attribute((__nonnull__(1))); 
# 235
extern const char *strrchr(const char * __s, int __c) throw() __asm__("strrchr")
# 236
 __attribute((__pure__)) __attribute((__nonnull__(1))); 
# 251 "/usr/include/string.h" 3
}
# 261 "/usr/include/string.h" 3
extern "C++" char *strchrnul(char * __s, int __c) throw() __asm__("strchrnul")
# 262
 __attribute((__pure__)) __attribute((__nonnull__(1))); 
# 263
extern "C++" const char *strchrnul(const char * __s, int __c) throw() __asm__("strchrnul")
# 264
 __attribute((__pure__)) __attribute((__nonnull__(1))); 
# 273 "/usr/include/string.h" 3
extern size_t strcspn(const char * __s, const char * __reject) throw()
# 274
 __attribute((__pure__)) __attribute((__nonnull__(1, 2))); 
# 277
extern size_t strspn(const char * __s, const char * __accept) throw()
# 278
 __attribute((__pure__)) __attribute((__nonnull__(1, 2))); 
# 281
extern "C++" {
# 283
extern char *strpbrk(char * __s, const char * __accept) throw() __asm__("strpbrk")
# 284
 __attribute((__pure__)) __attribute((__nonnull__(1, 2))); 
# 285
extern const char *strpbrk(const char * __s, const char * __accept) throw() __asm__("strpbrk")
# 286
 __attribute((__pure__)) __attribute((__nonnull__(1, 2))); 
# 301 "/usr/include/string.h" 3
}
# 308
extern "C++" {
# 310
extern char *strstr(char * __haystack, const char * __needle) throw() __asm__("strstr")
# 311
 __attribute((__pure__)) __attribute((__nonnull__(1, 2))); 
# 312
extern const char *strstr(const char * __haystack, const char * __needle) throw() __asm__("strstr")
# 313
 __attribute((__pure__)) __attribute((__nonnull__(1, 2))); 
# 328 "/usr/include/string.h" 3
}
# 336
extern char *strtok(char *__restrict__ __s, const char *__restrict__ __delim) throw()
# 337
 __attribute((__nonnull__(2))); 
# 341
extern char *__strtok_r(char *__restrict__ __s, const char *__restrict__ __delim, char **__restrict__ __save_ptr) throw()
# 344
 __attribute((__nonnull__(2, 3))); 
# 346
extern char *strtok_r(char *__restrict__ __s, const char *__restrict__ __delim, char **__restrict__ __save_ptr) throw()
# 348
 __attribute((__nonnull__(2, 3))); 
# 354
extern "C++" char *strcasestr(char * __haystack, const char * __needle) throw() __asm__("strcasestr")
# 355
 __attribute((__pure__)) __attribute((__nonnull__(1, 2))); 
# 356
extern "C++" const char *strcasestr(const char * __haystack, const char * __needle) throw() __asm__("strcasestr")
# 358
 __attribute((__pure__)) __attribute((__nonnull__(1, 2))); 
# 369 "/usr/include/string.h" 3
extern void *memmem(const void * __haystack, size_t __haystacklen, const void * __needle, size_t __needlelen) throw()
# 371
 __attribute((__pure__)) __attribute((__nonnull__(1, 3))); 
# 375
extern void *__mempcpy(void *__restrict__ __dest, const void *__restrict__ __src, size_t __n) throw()
# 377
 __attribute((__nonnull__(1, 2))); 
# 378
extern void *mempcpy(void *__restrict__ __dest, const void *__restrict__ __src, size_t __n) throw()
# 380
 __attribute((__nonnull__(1, 2))); 
# 385
extern size_t strlen(const char * __s) throw()
# 386
 __attribute((__pure__)) __attribute((__nonnull__(1))); 
# 391
extern size_t strnlen(const char * __string, size_t __maxlen) throw()
# 392
 __attribute((__pure__)) __attribute((__nonnull__(1))); 
# 397
extern char *strerror(int __errnum) throw(); 
# 421 "/usr/include/string.h" 3
extern char *strerror_r(int __errnum, char * __buf, size_t __buflen) throw()
# 422
 __attribute((__nonnull__(2))); 
# 428
extern char *strerror_l(int __errnum, locale_t __l) throw(); 
# 30 "/usr/include/strings.h" 3
extern "C" {
# 34
extern int bcmp(const void * __s1, const void * __s2, size_t __n) throw()
# 35
 __attribute((__pure__)) __attribute((__nonnull__(1, 2))); 
# 38
extern void bcopy(const void * __src, void * __dest, size_t __n) throw()
# 39
 __attribute((__nonnull__(1, 2))); 
# 42
extern void bzero(void * __s, size_t __n) throw() __attribute((__nonnull__(1))); 
# 46
extern "C++" {
# 48
extern char *index(char * __s, int __c) throw() __asm__("index")
# 49
 __attribute((__pure__)) __attribute((__nonnull__(1))); 
# 50
extern const char *index(const char * __s, int __c) throw() __asm__("index")
# 51
 __attribute((__pure__)) __attribute((__nonnull__(1))); 
# 66 "/usr/include/strings.h" 3
}
# 74
extern "C++" {
# 76
extern char *rindex(char * __s, int __c) throw() __asm__("rindex")
# 77
 __attribute((__pure__)) __attribute((__nonnull__(1))); 
# 78
extern const char *rindex(const char * __s, int __c) throw() __asm__("rindex")
# 79
 __attribute((__pure__)) __attribute((__nonnull__(1))); 
# 94 "/usr/include/strings.h" 3
}
# 104 "/usr/include/strings.h" 3
extern int ffs(int __i) throw() __attribute((const)); 
# 110
extern int ffsl(long __l) throw() __attribute((const)); 
# 111
__extension__ extern int ffsll(long long __ll) throw()
# 112
 __attribute((const)); 
# 116
extern int strcasecmp(const char * __s1, const char * __s2) throw()
# 117
 __attribute((__pure__)) __attribute((__nonnull__(1, 2))); 
# 120
extern int strncasecmp(const char * __s1, const char * __s2, size_t __n) throw()
# 121
 __attribute((__pure__)) __attribute((__nonnull__(1, 2))); 
# 128
extern int strcasecmp_l(const char * __s1, const char * __s2, locale_t __loc) throw()
# 129
 __attribute((__pure__)) __attribute((__nonnull__(1, 2, 3))); 
# 133
extern int strncasecmp_l(const char * __s1, const char * __s2, size_t __n, locale_t __loc) throw()
# 135
 __attribute((__pure__)) __attribute((__nonnull__(1, 2, 4))); 
# 138
}
# 436 "/usr/include/string.h" 3
extern void explicit_bzero(void * __s, size_t __n) throw() __attribute((__nonnull__(1))); 
# 440
extern char *strsep(char **__restrict__ __stringp, const char *__restrict__ __delim) throw()
# 442
 __attribute((__nonnull__(1, 2))); 
# 447
extern char *strsignal(int __sig) throw(); 
# 450
extern char *__stpcpy(char *__restrict__ __dest, const char *__restrict__ __src) throw()
# 451
 __attribute((__nonnull__(1, 2))); 
# 452
extern char *stpcpy(char *__restrict__ __dest, const char *__restrict__ __src) throw()
# 453
 __attribute((__nonnull__(1, 2))); 
# 457
extern char *__stpncpy(char *__restrict__ __dest, const char *__restrict__ __src, size_t __n) throw()
# 459
 __attribute((__nonnull__(1, 2))); 
# 460
extern char *stpncpy(char *__restrict__ __dest, const char *__restrict__ __src, size_t __n) throw()
# 462
 __attribute((__nonnull__(1, 2))); 
# 467
extern int strverscmp(const char * __s1, const char * __s2) throw()
# 468
 __attribute((__pure__)) __attribute((__nonnull__(1, 2))); 
# 471
extern char *strfry(char * __string) throw() __attribute((__nonnull__(1))); 
# 474
extern void *memfrob(void * __s, size_t __n) throw() __attribute((__nonnull__(1))); 
# 482
extern "C++" char *basename(char * __filename) throw() __asm__("basename")
# 483
 __attribute((__nonnull__(1))); 
# 484
extern "C++" const char *basename(const char * __filename) throw() __asm__("basename")
# 485
 __attribute((__nonnull__(1))); 
# 499 "/usr/include/string.h" 3
}
# 26 "/usr/include/x86_64-linux-gnu/bits/timex.h" 3
struct timex { 
# 28
unsigned modes; 
# 29
__syscall_slong_t offset; 
# 30
__syscall_slong_t freq; 
# 31
__syscall_slong_t maxerror; 
# 32
__syscall_slong_t esterror; 
# 33
int status; 
# 34
__syscall_slong_t constant; 
# 35
__syscall_slong_t precision; 
# 36
__syscall_slong_t tolerance; 
# 37
timeval time; 
# 38
__syscall_slong_t tick; 
# 39
__syscall_slong_t ppsfreq; 
# 40
__syscall_slong_t jitter; 
# 41
int shift; 
# 42
__syscall_slong_t stabil; 
# 43
__syscall_slong_t jitcnt; 
# 44
__syscall_slong_t calcnt; 
# 45
__syscall_slong_t errcnt; 
# 46
__syscall_slong_t stbcnt; 
# 48
int tai; 
# 51
int:32; int:32; int:32; int:32; 
# 52
int:32; int:32; int:32; int:32; 
# 53
int:32; int:32; int:32; 
# 54
}; 
# 75 "/usr/include/x86_64-linux-gnu/bits/time.h" 3
extern "C" {
# 78
extern int clock_adjtime(__clockid_t __clock_id, timex * __utx) throw(); 
# 80
}
# 7 "/usr/include/x86_64-linux-gnu/bits/types/struct_tm.h" 3
struct tm { 
# 9
int tm_sec; 
# 10
int tm_min; 
# 11
int tm_hour; 
# 12
int tm_mday; 
# 13
int tm_mon; 
# 14
int tm_year; 
# 15
int tm_wday; 
# 16
int tm_yday; 
# 17
int tm_isdst; 
# 20
long tm_gmtoff; 
# 21
const char *tm_zone; 
# 26
}; 
# 8 "/usr/include/x86_64-linux-gnu/bits/types/struct_itimerspec.h" 3
struct itimerspec { 
# 10
timespec it_interval; 
# 11
timespec it_value; 
# 12
}; 
# 49 "/usr/include/time.h" 3
struct sigevent; 
# 68 "/usr/include/time.h" 3
extern "C" {
# 72
extern clock_t clock() throw(); 
# 75
extern time_t time(time_t * __timer) throw(); 
# 78
extern double difftime(time_t __time1, time_t __time0) throw()
# 79
 __attribute((const)); 
# 82
extern time_t mktime(tm * __tp) throw(); 
# 88
extern size_t strftime(char *__restrict__ __s, size_t __maxsize, const char *__restrict__ __format, const tm *__restrict__ __tp) throw(); 
# 95
extern char *strptime(const char *__restrict__ __s, const char *__restrict__ __fmt, tm * __tp) throw(); 
# 104
extern size_t strftime_l(char *__restrict__ __s, size_t __maxsize, const char *__restrict__ __format, const tm *__restrict__ __tp, locale_t __loc) throw(); 
# 111
extern char *strptime_l(const char *__restrict__ __s, const char *__restrict__ __fmt, tm * __tp, locale_t __loc) throw(); 
# 119
extern tm *gmtime(const time_t * __timer) throw(); 
# 123
extern tm *localtime(const time_t * __timer) throw(); 
# 128
extern tm *gmtime_r(const time_t *__restrict__ __timer, tm *__restrict__ __tp) throw(); 
# 133
extern tm *localtime_r(const time_t *__restrict__ __timer, tm *__restrict__ __tp) throw(); 
# 139
extern char *asctime(const tm * __tp) throw(); 
# 142
extern char *ctime(const time_t * __timer) throw(); 
# 149
extern char *asctime_r(const tm *__restrict__ __tp, char *__restrict__ __buf) throw(); 
# 153
extern char *ctime_r(const time_t *__restrict__ __timer, char *__restrict__ __buf) throw(); 
# 159
extern char *__tzname[2]; 
# 160
extern int __daylight; 
# 161
extern long __timezone; 
# 166
extern char *tzname[2]; 
# 170
extern void tzset() throw(); 
# 174
extern int daylight; 
# 175
extern long timezone; 
# 190 "/usr/include/time.h" 3
extern time_t timegm(tm * __tp) throw(); 
# 193
extern time_t timelocal(tm * __tp) throw(); 
# 196
extern int dysize(int __year) throw() __attribute((const)); 
# 205 "/usr/include/time.h" 3
extern int nanosleep(const timespec * __requested_time, timespec * __remaining); 
# 210
extern int clock_getres(clockid_t __clock_id, timespec * __res) throw(); 
# 213
extern int clock_gettime(clockid_t __clock_id, timespec * __tp) throw(); 
# 216
extern int clock_settime(clockid_t __clock_id, const timespec * __tp) throw(); 
# 224
extern int clock_nanosleep(clockid_t __clock_id, int __flags, const timespec * __req, timespec * __rem); 
# 229
extern int clock_getcpuclockid(pid_t __pid, clockid_t * __clock_id) throw(); 
# 234
extern int timer_create(clockid_t __clock_id, sigevent *__restrict__ __evp, timer_t *__restrict__ __timerid) throw(); 
# 239
extern int timer_delete(timer_t __timerid) throw(); 
# 242
extern int timer_settime(timer_t __timerid, int __flags, const itimerspec *__restrict__ __value, itimerspec *__restrict__ __ovalue) throw(); 
# 247
extern int timer_gettime(timer_t __timerid, itimerspec * __value) throw(); 
# 251
extern int timer_getoverrun(timer_t __timerid) throw(); 
# 257
extern int timespec_get(timespec * __ts, int __base) throw()
# 258
 __attribute((__nonnull__(1))); 
# 274 "/usr/include/time.h" 3
extern int getdate_err; 
# 283 "/usr/include/time.h" 3
extern tm *getdate(const char * __string); 
# 297 "/usr/include/time.h" 3
extern int getdate_r(const char *__restrict__ __string, tm *__restrict__ __resbufp); 
# 301
}
# 88 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/common_functions.h"
extern "C" {
# 91
extern clock_t clock() throw(); 
# 96 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/common_functions.h"
extern void *memset(void *, int, size_t) throw(); 
# 97 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/common_functions.h"
extern void *memcpy(void *, const void *, size_t) throw(); 
# 99 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/common_functions.h"
}
# 121 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern "C" {
# 219 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern int abs(int a) throw(); 
# 227 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern long labs(long a) throw(); 
# 235 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern long long llabs(long long a) throw(); 
# 285 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern double fabs(double x) throw(); 
# 328 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern float fabsf(float x) throw(); 
# 338 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern inline int min(const int a, const int b); 
# 345
extern inline unsigned umin(const unsigned a, const unsigned b); 
# 352
extern inline long long llmin(const long long a, const long long b); 
# 359
extern inline unsigned long long ullmin(const unsigned long long a, const unsigned long long b); 
# 380 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern float fminf(float x, float y) throw(); 
# 400 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern double fmin(double x, double y) throw(); 
# 413 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern inline int max(const int a, const int b); 
# 421
extern inline unsigned umax(const unsigned a, const unsigned b); 
# 428
extern inline long long llmax(const long long a, const long long b); 
# 435
extern inline unsigned long long ullmax(const unsigned long long a, const unsigned long long b); 
# 456 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern float fmaxf(float x, float y) throw(); 
# 476 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern double fmax(double, double) throw(); 
# 520 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern double sin(double x) throw(); 
# 553 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern double cos(double x) throw(); 
# 572 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern void sincos(double x, double * sptr, double * cptr) throw(); 
# 588 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern void sincosf(float x, float * sptr, float * cptr) throw(); 
# 633 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern double tan(double x) throw(); 
# 702 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern double sqrt(double x) throw(); 
# 774 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern double rsqrt(double x); 
# 844 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern float rsqrtf(float x); 
# 900 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern double log2(double x) throw(); 
# 965 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern double exp2(double x) throw(); 
# 1030 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern float exp2f(float x) throw(); 
# 1097 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern double exp10(double x) throw(); 
# 1160 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern float exp10f(float x) throw(); 
# 1253 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern double expm1(double x) throw(); 
# 1345 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern float expm1f(float x) throw(); 
# 1401 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern float log2f(float x) throw(); 
# 1455 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern double log10(double x) throw(); 
# 1525 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern double log(double x) throw(); 
# 1621 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern double log1p(double x) throw(); 
# 1720 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern float log1pf(float x) throw(); 
# 1784 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern double floor(double x) throw(); 
# 1863 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern double exp(double x) throw(); 
# 1904 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern double cosh(double x) throw(); 
# 1954 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern double sinh(double x) throw(); 
# 2004 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern double tanh(double x) throw(); 
# 2059 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern double acosh(double x) throw(); 
# 2117 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern float acoshf(float x) throw(); 
# 2170 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern double asinh(double x) throw(); 
# 2223 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern float asinhf(float x) throw(); 
# 2277 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern double atanh(double x) throw(); 
# 2331 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern float atanhf(float x) throw(); 
# 2380 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern double ldexp(double x, int exp) throw(); 
# 2426 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern float ldexpf(float x, int exp) throw(); 
# 2478 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern double logb(double x) throw(); 
# 2533 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern float logbf(float x) throw(); 
# 2573 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern int ilogb(double x) throw(); 
# 2613 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern int ilogbf(float x) throw(); 
# 2689 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern double scalbn(double x, int n) throw(); 
# 2765 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern float scalbnf(float x, int n) throw(); 
# 2841 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern double scalbln(double x, long n) throw(); 
# 2917 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern float scalblnf(float x, long n) throw(); 
# 2994 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern double frexp(double x, int * nptr) throw(); 
# 3068 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern float frexpf(float x, int * nptr) throw(); 
# 3120 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern double round(double x) throw(); 
# 3175 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern float roundf(float x) throw(); 
# 3193 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern long lround(double x) throw(); 
# 3211 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern long lroundf(float x) throw(); 
# 3229 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern long long llround(double x) throw(); 
# 3247 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern long long llroundf(float x) throw(); 
# 3375 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern float rintf(float x) throw(); 
# 3392 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern long lrint(double x) throw(); 
# 3409 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern long lrintf(float x) throw(); 
# 3426 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern long long llrint(double x) throw(); 
# 3443 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern long long llrintf(float x) throw(); 
# 3496 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern double nearbyint(double x) throw(); 
# 3549 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern float nearbyintf(float x) throw(); 
# 3611 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern double ceil(double x) throw(); 
# 3661 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern double trunc(double x) throw(); 
# 3714 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern float truncf(float x) throw(); 
# 3740 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern double fdim(double x, double y) throw(); 
# 3766 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern float fdimf(float x, float y) throw(); 
# 4066 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern double atan2(double y, double x) throw(); 
# 4137 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern double atan(double x) throw(); 
# 4160 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern double acos(double x) throw(); 
# 4211 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern double asin(double x) throw(); 
# 4279 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern double hypot(double x, double y) throw(); 
# 4402 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern float hypotf(float x, float y) throw(); 
# 5188 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern double cbrt(double x) throw(); 
# 5274 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern float cbrtf(float x) throw(); 
# 5329 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern double rcbrt(double x); 
# 5379 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern float rcbrtf(float x); 
# 5439 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern double sinpi(double x); 
# 5499 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern float sinpif(float x); 
# 5551 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern double cospi(double x); 
# 5603 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern float cospif(float x); 
# 5633 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern void sincospi(double x, double * sptr, double * cptr); 
# 5663 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern void sincospif(float x, float * sptr, float * cptr); 
# 5996 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern double pow(double x, double y) throw(); 
# 6052 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern double modf(double x, double * iptr) throw(); 
# 6111 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern double fmod(double x, double y) throw(); 
# 6207 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern double remainder(double x, double y) throw(); 
# 6306 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern float remainderf(float x, float y) throw(); 
# 6378 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern double remquo(double x, double y, int * quo) throw(); 
# 6450 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern float remquof(float x, float y, int * quo) throw(); 
# 6491 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern double j0(double x) throw(); 
# 6533 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern float j0f(float x) throw(); 
# 6602 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern double j1(double x) throw(); 
# 6671 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern float j1f(float x) throw(); 
# 6714 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern double jn(int n, double x) throw(); 
# 6757 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern float jnf(int n, float x) throw(); 
# 6818 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern double y0(double x) throw(); 
# 6879 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern float y0f(float x) throw(); 
# 6940 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern double y1(double x) throw(); 
# 7001 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern float y1f(float x) throw(); 
# 7064 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern double yn(int n, double x) throw(); 
# 7127 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern float ynf(int n, float x) throw(); 
# 7316 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern double erf(double x) throw(); 
# 7398 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern float erff(float x) throw(); 
# 7470 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern double erfinv(double x); 
# 7535 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern float erfinvf(float x); 
# 7574 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern double erfc(double x) throw(); 
# 7612 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern float erfcf(float x) throw(); 
# 7729 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern double lgamma(double x) throw(); 
# 7791 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern double erfcinv(double x); 
# 7846 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern float erfcinvf(float x); 
# 7914 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern double normcdfinv(double x); 
# 7982 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern float normcdfinvf(float x); 
# 8025 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern double normcdf(double x); 
# 8068 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern float normcdff(float x); 
# 8132 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern double erfcx(double x); 
# 8196 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern float erfcxf(float x); 
# 8315 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern float lgammaf(float x) throw(); 
# 8413 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern double tgamma(double x) throw(); 
# 8511 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern float tgammaf(float x) throw(); 
# 8524 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern double copysign(double x, double y) throw(); 
# 8537 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern float copysignf(float x, float y) throw(); 
# 8556 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern double nextafter(double x, double y) throw(); 
# 8575 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern float nextafterf(float x, float y) throw(); 
# 8591 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern double nan(const char * tagp) throw(); 
# 8607 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern float nanf(const char * tagp) throw(); 
# 8614 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern int __isinff(float) throw(); 
# 8615 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern int __isnanf(float) throw(); 
# 8625 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern int __finite(double) throw(); 
# 8626 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern int __finitef(float) throw(); 
# 8627 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern int __signbit(double) throw(); 
# 8628 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern int __isnan(double) throw(); 
# 8629 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern int __isinf(double) throw(); 
# 8632 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern int __signbitf(float) throw(); 
# 8791 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern double fma(double x, double y, double z) throw(); 
# 8949 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern float fmaf(float x, float y, float z) throw(); 
# 8960 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern int __signbitl(long double) throw(); 
# 8966 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern int __finitel(long double) throw(); 
# 8967 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern int __isinfl(long double) throw(); 
# 8968 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern int __isnanl(long double) throw(); 
# 9018 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern float acosf(float x) throw(); 
# 9077 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern float asinf(float x) throw(); 
# 9157 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern float atanf(float x) throw(); 
# 9454 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern float atan2f(float y, float x) throw(); 
# 9488 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern float cosf(float x) throw(); 
# 9530 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern float sinf(float x) throw(); 
# 9572 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern float tanf(float x) throw(); 
# 9613 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern float coshf(float x) throw(); 
# 9663 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern float sinhf(float x) throw(); 
# 9713 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern float tanhf(float x) throw(); 
# 9765 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern float logf(float x) throw(); 
# 9845 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern float expf(float x) throw(); 
# 9897 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern float log10f(float x) throw(); 
# 9952 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern float modff(float x, float * iptr) throw(); 
# 10282 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern float powf(float x, float y) throw(); 
# 10351 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern float sqrtf(float x) throw(); 
# 10410 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern float ceilf(float x) throw(); 
# 10471 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern float floorf(float x) throw(); 
# 10529 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern float fmodf(float x, float y) throw(); 
# 10544 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
}
# 67 "/usr/include/c++/10/bits/cpp_type_traits.h" 3
extern "C++" {
# 69
namespace std __attribute((__visibility__("default"))) { 
# 73
struct __true_type { }; 
# 74
struct __false_type { }; 
# 76
template< bool > 
# 77
struct __truth_type { 
# 78
typedef __false_type __type; }; 
# 81
template<> struct __truth_type< true>  { 
# 82
typedef __true_type __type; }; 
# 86
template< class _Sp, class _Tp> 
# 87
struct __traitor { 
# 89
enum { __value = ((bool)_Sp::__value) || ((bool)_Tp::__value)}; 
# 90
typedef typename __truth_type< __value> ::__type __type; 
# 91
}; 
# 94
template< class , class > 
# 95
struct __are_same { 
# 97
enum { __value}; 
# 98
typedef __false_type __type; 
# 99
}; 
# 101
template< class _Tp> 
# 102
struct __are_same< _Tp, _Tp>  { 
# 104
enum { __value = 1}; 
# 105
typedef __true_type __type; 
# 106
}; 
# 109
template< class _Tp> 
# 110
struct __is_void { 
# 112
enum { __value}; 
# 113
typedef __false_type __type; 
# 114
}; 
# 117
template<> struct __is_void< void>  { 
# 119
enum { __value = 1}; 
# 120
typedef __true_type __type; 
# 121
}; 
# 126
template< class _Tp> 
# 127
struct __is_integer { 
# 129
enum { __value}; 
# 130
typedef __false_type __type; 
# 131
}; 
# 138
template<> struct __is_integer< bool>  { 
# 140
enum { __value = 1}; 
# 141
typedef __true_type __type; 
# 142
}; 
# 145
template<> struct __is_integer< char>  { 
# 147
enum { __value = 1}; 
# 148
typedef __true_type __type; 
# 149
}; 
# 152
template<> struct __is_integer< signed char>  { 
# 154
enum { __value = 1}; 
# 155
typedef __true_type __type; 
# 156
}; 
# 159
template<> struct __is_integer< unsigned char>  { 
# 161
enum { __value = 1}; 
# 162
typedef __true_type __type; 
# 163
}; 
# 167
template<> struct __is_integer< wchar_t>  { 
# 169
enum { __value = 1}; 
# 170
typedef __true_type __type; 
# 171
}; 
# 185 "/usr/include/c++/10/bits/cpp_type_traits.h" 3
template<> struct __is_integer< char16_t>  { 
# 187
enum { __value = 1}; 
# 188
typedef __true_type __type; 
# 189
}; 
# 192
template<> struct __is_integer< char32_t>  { 
# 194
enum { __value = 1}; 
# 195
typedef __true_type __type; 
# 196
}; 
# 200
template<> struct __is_integer< short>  { 
# 202
enum { __value = 1}; 
# 203
typedef __true_type __type; 
# 204
}; 
# 207
template<> struct __is_integer< unsigned short>  { 
# 209
enum { __value = 1}; 
# 210
typedef __true_type __type; 
# 211
}; 
# 214
template<> struct __is_integer< int>  { 
# 216
enum { __value = 1}; 
# 217
typedef __true_type __type; 
# 218
}; 
# 221
template<> struct __is_integer< unsigned>  { 
# 223
enum { __value = 1}; 
# 224
typedef __true_type __type; 
# 225
}; 
# 228
template<> struct __is_integer< long>  { 
# 230
enum { __value = 1}; 
# 231
typedef __true_type __type; 
# 232
}; 
# 235
template<> struct __is_integer< unsigned long>  { 
# 237
enum { __value = 1}; 
# 238
typedef __true_type __type; 
# 239
}; 
# 242
template<> struct __is_integer< long long>  { 
# 244
enum { __value = 1}; 
# 245
typedef __true_type __type; 
# 246
}; 
# 249
template<> struct __is_integer< unsigned long long>  { 
# 251
enum { __value = 1}; 
# 252
typedef __true_type __type; 
# 253
}; 
# 270 "/usr/include/c++/10/bits/cpp_type_traits.h" 3
template<> struct __is_integer< __int128>  { enum { __value = 1}; typedef __true_type __type; }; template<> struct __is_integer< unsigned __int128>  { enum { __value = 1}; typedef __true_type __type; }; 
# 287 "/usr/include/c++/10/bits/cpp_type_traits.h" 3
template< class _Tp> 
# 288
struct __is_floating { 
# 290
enum { __value}; 
# 291
typedef __false_type __type; 
# 292
}; 
# 296
template<> struct __is_floating< float>  { 
# 298
enum { __value = 1}; 
# 299
typedef __true_type __type; 
# 300
}; 
# 303
template<> struct __is_floating< double>  { 
# 305
enum { __value = 1}; 
# 306
typedef __true_type __type; 
# 307
}; 
# 310
template<> struct __is_floating< long double>  { 
# 312
enum { __value = 1}; 
# 313
typedef __true_type __type; 
# 314
}; 
# 319
template< class _Tp> 
# 320
struct __is_pointer { 
# 322
enum { __value}; 
# 323
typedef __false_type __type; 
# 324
}; 
# 326
template< class _Tp> 
# 327
struct __is_pointer< _Tp *>  { 
# 329
enum { __value = 1}; 
# 330
typedef __true_type __type; 
# 331
}; 
# 336
template< class _Tp> 
# 337
struct __is_arithmetic : public __traitor< __is_integer< _Tp> , __is_floating< _Tp> >  { 
# 339
}; 
# 344
template< class _Tp> 
# 345
struct __is_scalar : public __traitor< __is_arithmetic< _Tp> , __is_pointer< _Tp> >  { 
# 347
}; 
# 352
template< class _Tp> 
# 353
struct __is_char { 
# 355
enum { __value}; 
# 356
typedef __false_type __type; 
# 357
}; 
# 360
template<> struct __is_char< char>  { 
# 362
enum { __value = 1}; 
# 363
typedef __true_type __type; 
# 364
}; 
# 368
template<> struct __is_char< wchar_t>  { 
# 370
enum { __value = 1}; 
# 371
typedef __true_type __type; 
# 372
}; 
# 375
template< class _Tp> 
# 376
struct __is_byte { 
# 378
enum { __value}; 
# 379
typedef __false_type __type; 
# 380
}; 
# 383
template<> struct __is_byte< char>  { 
# 385
enum { __value = 1}; 
# 386
typedef __true_type __type; 
# 387
}; 
# 390
template<> struct __is_byte< signed char>  { 
# 392
enum { __value = 1}; 
# 393
typedef __true_type __type; 
# 394
}; 
# 397
template<> struct __is_byte< unsigned char>  { 
# 399
enum { __value = 1}; 
# 400
typedef __true_type __type; 
# 401
}; 
# 423 "/usr/include/c++/10/bits/cpp_type_traits.h" 3
template< class > struct iterator_traits; 
# 426
template< class _Tp> 
# 427
struct __is_nonvolatile_trivially_copyable { 
# 429
enum { __value = __is_trivially_copyable(_Tp)}; 
# 430
}; 
# 435
template< class _Tp> 
# 436
struct __is_nonvolatile_trivially_copyable< volatile _Tp>  { 
# 438
enum { __value}; 
# 439
}; 
# 442
template< class _OutputIter, class _InputIter> 
# 443
struct __memcpyable { 
# 445
enum { __value}; 
# 446
}; 
# 448
template< class _Tp> 
# 449
struct __memcpyable< _Tp *, _Tp *>  : public __is_nonvolatile_trivially_copyable< _Tp>  { 
# 451
}; 
# 453
template< class _Tp> 
# 454
struct __memcpyable< _Tp *, const _Tp *>  : public __is_nonvolatile_trivially_copyable< _Tp>  { 
# 456
}; 
# 463
template< class _Iter1, class _Iter2> 
# 464
struct __memcmpable { 
# 466
enum { __value}; 
# 467
}; 
# 470
template< class _Tp> 
# 471
struct __memcmpable< _Tp *, _Tp *>  : public __is_nonvolatile_trivially_copyable< _Tp>  { 
# 473
}; 
# 475
template< class _Tp> 
# 476
struct __memcmpable< const _Tp *, _Tp *>  : public __is_nonvolatile_trivially_copyable< _Tp>  { 
# 478
}; 
# 480
template< class _Tp> 
# 481
struct __memcmpable< _Tp *, const _Tp *>  : public __is_nonvolatile_trivially_copyable< _Tp>  { 
# 483
}; 
# 488
template< class _Tp, bool _TreatAsBytes = __is_byte< _Tp> ::__value> 
# 489
struct __is_memcmp_ordered { 
# 491
static const bool __value = (((_Tp)(-1)) > ((_Tp)1)); 
# 492
}; 
# 494
template< class _Tp> 
# 495
struct __is_memcmp_ordered< _Tp, false>  { 
# 497
static const bool __value = false; 
# 498
}; 
# 501
template< class _Tp, class _Up, bool  = sizeof(_Tp) == sizeof(_Up)> 
# 502
struct __is_memcmp_ordered_with { 
# 504
static const bool __value = (__is_memcmp_ordered< _Tp> ::__value && __is_memcmp_ordered< _Up> ::__value); 
# 506
}; 
# 508
template< class _Tp, class _Up> 
# 509
struct __is_memcmp_ordered_with< _Tp, _Up, false>  { 
# 511
static const bool __value = false; 
# 512
}; 
# 532 "/usr/include/c++/10/bits/cpp_type_traits.h" 3
template< class _Tp> 
# 533
struct __is_move_iterator { 
# 535
enum { __value}; 
# 536
typedef __false_type __type; 
# 537
}; 
# 541
template< class _Iterator> inline _Iterator 
# 544
__miter_base(_Iterator __it) 
# 545
{ return __it; } 
# 548
}
# 549
}
# 37 "/usr/include/c++/10/ext/type_traits.h" 3
extern "C++" {
# 39
namespace __gnu_cxx __attribute((__visibility__("default"))) { 
# 44
template< bool , class > 
# 45
struct __enable_if { 
# 46
}; 
# 48
template< class _Tp> 
# 49
struct __enable_if< true, _Tp>  { 
# 50
typedef _Tp __type; }; 
# 54
template< bool _Cond, class _Iftrue, class _Iffalse> 
# 55
struct __conditional_type { 
# 56
typedef _Iftrue __type; }; 
# 58
template< class _Iftrue, class _Iffalse> 
# 59
struct __conditional_type< false, _Iftrue, _Iffalse>  { 
# 60
typedef _Iffalse __type; }; 
# 64
template< class _Tp> 
# 65
struct __add_unsigned { 
# 68
private: typedef __enable_if< std::__is_integer< _Tp> ::__value, _Tp>  __if_type; 
# 71
public: typedef typename __enable_if< std::__is_integer< _Tp> ::__value, _Tp> ::__type __type; 
# 72
}; 
# 75
template<> struct __add_unsigned< char>  { 
# 76
typedef unsigned char __type; }; 
# 79
template<> struct __add_unsigned< signed char>  { 
# 80
typedef unsigned char __type; }; 
# 83
template<> struct __add_unsigned< short>  { 
# 84
typedef unsigned short __type; }; 
# 87
template<> struct __add_unsigned< int>  { 
# 88
typedef unsigned __type; }; 
# 91
template<> struct __add_unsigned< long>  { 
# 92
typedef unsigned long __type; }; 
# 95
template<> struct __add_unsigned< long long>  { 
# 96
typedef unsigned long long __type; }; 
# 100
template<> struct __add_unsigned< bool> ; 
# 103
template<> struct __add_unsigned< wchar_t> ; 
# 107
template< class _Tp> 
# 108
struct __remove_unsigned { 
# 111
private: typedef __enable_if< std::__is_integer< _Tp> ::__value, _Tp>  __if_type; 
# 114
public: typedef typename __enable_if< std::__is_integer< _Tp> ::__value, _Tp> ::__type __type; 
# 115
}; 
# 118
template<> struct __remove_unsigned< char>  { 
# 119
typedef signed char __type; }; 
# 122
template<> struct __remove_unsigned< unsigned char>  { 
# 123
typedef signed char __type; }; 
# 126
template<> struct __remove_unsigned< unsigned short>  { 
# 127
typedef short __type; }; 
# 130
template<> struct __remove_unsigned< unsigned>  { 
# 131
typedef int __type; }; 
# 134
template<> struct __remove_unsigned< unsigned long>  { 
# 135
typedef long __type; }; 
# 138
template<> struct __remove_unsigned< unsigned long long>  { 
# 139
typedef long long __type; }; 
# 143
template<> struct __remove_unsigned< bool> ; 
# 146
template<> struct __remove_unsigned< wchar_t> ; 
# 150
template< class _Type> inline bool 
# 152
__is_null_pointer(_Type *__ptr) 
# 153
{ return __ptr == 0; } 
# 155
template< class _Type> inline bool 
# 157
__is_null_pointer(_Type) 
# 158
{ return false; } 
# 162
inline bool __is_null_pointer(std::nullptr_t) 
# 163
{ return true; } 
# 168
template< class _Tp, bool  = std::template __is_integer< _Tp> ::__value> 
# 169
struct __promote { 
# 170
typedef double __type; }; 
# 175
template< class _Tp> 
# 176
struct __promote< _Tp, false>  { 
# 177
}; 
# 180
template<> struct __promote< long double>  { 
# 181
typedef long double __type; }; 
# 184
template<> struct __promote< double>  { 
# 185
typedef double __type; }; 
# 188
template<> struct __promote< float>  { 
# 189
typedef float __type; }; 
# 196
template< class _Tp, class _Up, class 
# 197
_Tp2 = typename __promote< _Tp> ::__type, class 
# 198
_Up2 = typename __promote< _Up> ::__type> 
# 199
struct __promote_2 { 
# 201
typedef __typeof__(_Tp2() + _Up2()) __type; 
# 202
}; 
# 204
template< class _Tp, class _Up, class _Vp, class 
# 205
_Tp2 = typename __promote< _Tp> ::__type, class 
# 206
_Up2 = typename __promote< _Up> ::__type, class 
# 207
_Vp2 = typename __promote< _Vp> ::__type> 
# 208
struct __promote_3 { 
# 210
typedef __typeof__((_Tp2() + _Up2()) + _Vp2()) __type; 
# 211
}; 
# 213
template< class _Tp, class _Up, class _Vp, class _Wp, class 
# 214
_Tp2 = typename __promote< _Tp> ::__type, class 
# 215
_Up2 = typename __promote< _Up> ::__type, class 
# 216
_Vp2 = typename __promote< _Vp> ::__type, class 
# 217
_Wp2 = typename __promote< _Wp> ::__type> 
# 218
struct __promote_4 { 
# 220
typedef __typeof__(((_Tp2() + _Up2()) + _Vp2()) + _Wp2()) __type; 
# 221
}; 
# 224
}
# 225
}
# 34 "/usr/include/math.h" 3
extern "C" {
# 149 "/usr/include/math.h" 3
typedef float float_t; 
# 150
typedef double double_t; 
# 238 "/usr/include/math.h" 3
enum { 
# 239
FP_INT_UPWARD, 
# 242
FP_INT_DOWNWARD, 
# 245
FP_INT_TOWARDZERO, 
# 248
FP_INT_TONEARESTFROMZERO, 
# 251
FP_INT_TONEAREST
# 254
}; 
# 21 "/usr/include/x86_64-linux-gnu/bits/mathcalls-helper-functions.h" 3
extern int __fpclassify(double __value) throw()
# 22
 __attribute((const)); 
# 25
extern int __signbit(double __value) throw()
# 26
 __attribute((const)); 
# 30
extern int __isinf(double __value) throw() __attribute((const)); 
# 33
extern int __finite(double __value) throw() __attribute((const)); 
# 36
extern int __isnan(double __value) throw() __attribute((const)); 
# 39
extern int __iseqsig(double __x, double __y) throw(); 
# 42
extern int __issignaling(double __value) throw()
# 43
 __attribute((const)); 
# 53 "/usr/include/x86_64-linux-gnu/bits/mathcalls.h" 3
extern double acos(double __x) throw(); extern double __acos(double __x) throw(); 
# 55
extern double asin(double __x) throw(); extern double __asin(double __x) throw(); 
# 57
extern double atan(double __x) throw(); extern double __atan(double __x) throw(); 
# 59
extern double atan2(double __y, double __x) throw(); extern double __atan2(double __y, double __x) throw(); 
# 62
extern double cos(double __x) throw(); extern double __cos(double __x) throw(); 
# 64
extern double sin(double __x) throw(); extern double __sin(double __x) throw(); 
# 66
extern double tan(double __x) throw(); extern double __tan(double __x) throw(); 
# 71
extern double cosh(double __x) throw(); extern double __cosh(double __x) throw(); 
# 73
extern double sinh(double __x) throw(); extern double __sinh(double __x) throw(); 
# 75
extern double tanh(double __x) throw(); extern double __tanh(double __x) throw(); 
# 79
extern void sincos(double __x, double * __sinx, double * __cosx) throw(); extern void __sincos(double __x, double * __sinx, double * __cosx) throw(); 
# 85
extern double acosh(double __x) throw(); extern double __acosh(double __x) throw(); 
# 87
extern double asinh(double __x) throw(); extern double __asinh(double __x) throw(); 
# 89
extern double atanh(double __x) throw(); extern double __atanh(double __x) throw(); 
# 95
extern double exp(double __x) throw(); extern double __exp(double __x) throw(); 
# 98
extern double frexp(double __x, int * __exponent) throw(); extern double __frexp(double __x, int * __exponent) throw(); 
# 101
extern double ldexp(double __x, int __exponent) throw(); extern double __ldexp(double __x, int __exponent) throw(); 
# 104
extern double log(double __x) throw(); extern double __log(double __x) throw(); 
# 107
extern double log10(double __x) throw(); extern double __log10(double __x) throw(); 
# 110
extern double modf(double __x, double * __iptr) throw(); extern double __modf(double __x, double * __iptr) throw() __attribute((__nonnull__(2))); 
# 114
extern double exp10(double __x) throw(); extern double __exp10(double __x) throw(); 
# 119
extern double expm1(double __x) throw(); extern double __expm1(double __x) throw(); 
# 122
extern double log1p(double __x) throw(); extern double __log1p(double __x) throw(); 
# 125
extern double logb(double __x) throw(); extern double __logb(double __x) throw(); 
# 130
extern double exp2(double __x) throw(); extern double __exp2(double __x) throw(); 
# 133
extern double log2(double __x) throw(); extern double __log2(double __x) throw(); 
# 140
extern double pow(double __x, double __y) throw(); extern double __pow(double __x, double __y) throw(); 
# 143
extern double sqrt(double __x) throw(); extern double __sqrt(double __x) throw(); 
# 147
extern double hypot(double __x, double __y) throw(); extern double __hypot(double __x, double __y) throw(); 
# 152
extern double cbrt(double __x) throw(); extern double __cbrt(double __x) throw(); 
# 159
extern double ceil(double __x) throw() __attribute((const)); extern double __ceil(double __x) throw() __attribute((const)); 
# 162
extern double fabs(double __x) throw() __attribute((const)); extern double __fabs(double __x) throw() __attribute((const)); 
# 165
extern double floor(double __x) throw() __attribute((const)); extern double __floor(double __x) throw() __attribute((const)); 
# 168
extern double fmod(double __x, double __y) throw(); extern double __fmod(double __x, double __y) throw(); 
# 182 "/usr/include/x86_64-linux-gnu/bits/mathcalls.h" 3
extern int finite(double __value) throw() __attribute((const)); 
# 185
extern double drem(double __x, double __y) throw(); extern double __drem(double __x, double __y) throw(); 
# 189
extern double significand(double __x) throw(); extern double __significand(double __x) throw(); 
# 196
extern double copysign(double __x, double __y) throw() __attribute((const)); extern double __copysign(double __x, double __y) throw() __attribute((const)); 
# 201
extern double nan(const char * __tagb) throw(); extern double __nan(const char * __tagb) throw(); 
# 217 "/usr/include/x86_64-linux-gnu/bits/mathcalls.h" 3
extern double j0(double) throw(); extern double __j0(double) throw(); 
# 218
extern double j1(double) throw(); extern double __j1(double) throw(); 
# 219
extern double jn(int, double) throw(); extern double __jn(int, double) throw(); 
# 220
extern double y0(double) throw(); extern double __y0(double) throw(); 
# 221
extern double y1(double) throw(); extern double __y1(double) throw(); 
# 222
extern double yn(int, double) throw(); extern double __yn(int, double) throw(); 
# 228
extern double erf(double) throw(); extern double __erf(double) throw(); 
# 229
extern double erfc(double) throw(); extern double __erfc(double) throw(); 
# 230
extern double lgamma(double) throw(); extern double __lgamma(double) throw(); 
# 235
extern double tgamma(double) throw(); extern double __tgamma(double) throw(); 
# 241
extern double gamma(double) throw(); extern double __gamma(double) throw(); 
# 249
extern double lgamma_r(double, int * __signgamp) throw(); extern double __lgamma_r(double, int * __signgamp) throw(); 
# 256
extern double rint(double __x) throw(); extern double __rint(double __x) throw(); 
# 259
extern double nextafter(double __x, double __y) throw(); extern double __nextafter(double __x, double __y) throw(); 
# 261
extern double nexttoward(double __x, long double __y) throw(); extern double __nexttoward(double __x, long double __y) throw(); 
# 266
extern double nextdown(double __x) throw(); extern double __nextdown(double __x) throw(); 
# 268
extern double nextup(double __x) throw(); extern double __nextup(double __x) throw(); 
# 272
extern double remainder(double __x, double __y) throw(); extern double __remainder(double __x, double __y) throw(); 
# 276
extern double scalbn(double __x, int __n) throw(); extern double __scalbn(double __x, int __n) throw(); 
# 280
extern int ilogb(double __x) throw(); extern int __ilogb(double __x) throw(); 
# 285
extern long llogb(double __x) throw(); extern long __llogb(double __x) throw(); 
# 290
extern double scalbln(double __x, long __n) throw(); extern double __scalbln(double __x, long __n) throw(); 
# 294
extern double nearbyint(double __x) throw(); extern double __nearbyint(double __x) throw(); 
# 298
extern double round(double __x) throw() __attribute((const)); extern double __round(double __x) throw() __attribute((const)); 
# 302
extern double trunc(double __x) throw() __attribute((const)); extern double __trunc(double __x) throw() __attribute((const)); 
# 307
extern double remquo(double __x, double __y, int * __quo) throw(); extern double __remquo(double __x, double __y, int * __quo) throw(); 
# 314
extern long lrint(double __x) throw(); extern long __lrint(double __x) throw(); 
# 316
__extension__ extern long long llrint(double __x) throw(); extern long long __llrint(double __x) throw(); 
# 320
extern long lround(double __x) throw(); extern long __lround(double __x) throw(); 
# 322
__extension__ extern long long llround(double __x) throw(); extern long long __llround(double __x) throw(); 
# 326
extern double fdim(double __x, double __y) throw(); extern double __fdim(double __x, double __y) throw(); 
# 329
extern double fmax(double __x, double __y) throw() __attribute((const)); extern double __fmax(double __x, double __y) throw() __attribute((const)); 
# 332
extern double fmin(double __x, double __y) throw() __attribute((const)); extern double __fmin(double __x, double __y) throw() __attribute((const)); 
# 335
extern double fma(double __x, double __y, double __z) throw(); extern double __fma(double __x, double __y, double __z) throw(); 
# 340
extern double roundeven(double __x) throw() __attribute((const)); extern double __roundeven(double __x) throw() __attribute((const)); 
# 344
extern __intmax_t fromfp(double __x, int __round, unsigned __width) throw(); extern __intmax_t __fromfp(double __x, int __round, unsigned __width) throw(); 
# 349
extern __uintmax_t ufromfp(double __x, int __round, unsigned __width) throw(); extern __uintmax_t __ufromfp(double __x, int __round, unsigned __width) throw(); 
# 355
extern __intmax_t fromfpx(double __x, int __round, unsigned __width) throw(); extern __intmax_t __fromfpx(double __x, int __round, unsigned __width) throw(); 
# 361
extern __uintmax_t ufromfpx(double __x, int __round, unsigned __width) throw(); extern __uintmax_t __ufromfpx(double __x, int __round, unsigned __width) throw(); 
# 365
extern double fmaxmag(double __x, double __y) throw() __attribute((const)); extern double __fmaxmag(double __x, double __y) throw() __attribute((const)); 
# 368
extern double fminmag(double __x, double __y) throw() __attribute((const)); extern double __fminmag(double __x, double __y) throw() __attribute((const)); 
# 371
extern int canonicalize(double * __cx, const double * __x) throw(); 
# 376
extern int totalorder(const double * __x, const double * __y) throw()
# 378
 __attribute((__pure__)); 
# 381
extern int totalordermag(const double * __x, const double * __y) throw()
# 383
 __attribute((__pure__)); 
# 386
extern double getpayload(const double * __x) throw(); extern double __getpayload(const double * __x) throw(); 
# 389
extern int setpayload(double * __x, double __payload) throw(); 
# 392
extern int setpayloadsig(double * __x, double __payload) throw(); 
# 400
extern double scalb(double __x, double __n) throw(); extern double __scalb(double __x, double __n) throw(); 
# 21 "/usr/include/x86_64-linux-gnu/bits/mathcalls-helper-functions.h" 3
extern int __fpclassifyf(float __value) throw()
# 22
 __attribute((const)); 
# 25
extern int __signbitf(float __value) throw()
# 26
 __attribute((const)); 
# 30
extern int __isinff(float __value) throw() __attribute((const)); 
# 33
extern int __finitef(float __value) throw() __attribute((const)); 
# 36
extern int __isnanf(float __value) throw() __attribute((const)); 
# 39
extern int __iseqsigf(float __x, float __y) throw(); 
# 42
extern int __issignalingf(float __value) throw()
# 43
 __attribute((const)); 
# 53 "/usr/include/x86_64-linux-gnu/bits/mathcalls.h" 3
extern float acosf(float __x) throw(); extern float __acosf(float __x) throw(); 
# 55
extern float asinf(float __x) throw(); extern float __asinf(float __x) throw(); 
# 57
extern float atanf(float __x) throw(); extern float __atanf(float __x) throw(); 
# 59
extern float atan2f(float __y, float __x) throw(); extern float __atan2f(float __y, float __x) throw(); 
# 62
extern float cosf(float __x) throw(); 
# 64
extern float sinf(float __x) throw(); 
# 66
extern float tanf(float __x) throw(); 
# 71
extern float coshf(float __x) throw(); extern float __coshf(float __x) throw(); 
# 73
extern float sinhf(float __x) throw(); extern float __sinhf(float __x) throw(); 
# 75
extern float tanhf(float __x) throw(); extern float __tanhf(float __x) throw(); 
# 79
extern void sincosf(float __x, float * __sinx, float * __cosx) throw(); 
# 85
extern float acoshf(float __x) throw(); extern float __acoshf(float __x) throw(); 
# 87
extern float asinhf(float __x) throw(); extern float __asinhf(float __x) throw(); 
# 89
extern float atanhf(float __x) throw(); extern float __atanhf(float __x) throw(); 
# 95
extern float expf(float __x) throw(); 
# 98
extern float frexpf(float __x, int * __exponent) throw(); extern float __frexpf(float __x, int * __exponent) throw(); 
# 101
extern float ldexpf(float __x, int __exponent) throw(); extern float __ldexpf(float __x, int __exponent) throw(); 
# 104
extern float logf(float __x) throw(); 
# 107
extern float log10f(float __x) throw(); 
# 110
extern float modff(float __x, float * __iptr) throw(); extern float __modff(float __x, float * __iptr) throw() __attribute((__nonnull__(2))); 
# 114
extern float exp10f(float __x) throw(); 
# 119
extern float expm1f(float __x) throw(); extern float __expm1f(float __x) throw(); 
# 122
extern float log1pf(float __x) throw(); extern float __log1pf(float __x) throw(); 
# 125
extern float logbf(float __x) throw(); extern float __logbf(float __x) throw(); 
# 130
extern float exp2f(float __x) throw(); extern float __exp2f(float __x) throw(); 
# 133
extern float log2f(float __x) throw(); 
# 140
extern float powf(float __x, float __y) throw(); 
# 143
extern float sqrtf(float __x) throw(); extern float __sqrtf(float __x) throw(); 
# 147
extern float hypotf(float __x, float __y) throw(); extern float __hypotf(float __x, float __y) throw(); 
# 152
extern float cbrtf(float __x) throw(); extern float __cbrtf(float __x) throw(); 
# 159
extern float ceilf(float __x) throw() __attribute((const)); extern float __ceilf(float __x) throw() __attribute((const)); 
# 162
extern float fabsf(float __x) throw() __attribute((const)); extern float __fabsf(float __x) throw() __attribute((const)); 
# 165
extern float floorf(float __x) throw() __attribute((const)); extern float __floorf(float __x) throw() __attribute((const)); 
# 168
extern float fmodf(float __x, float __y) throw(); extern float __fmodf(float __x, float __y) throw(); 
# 177 "/usr/include/x86_64-linux-gnu/bits/mathcalls.h" 3
extern int isinff(float __value) throw() __attribute((const)); 
# 182
extern int finitef(float __value) throw() __attribute((const)); 
# 185
extern float dremf(float __x, float __y) throw(); extern float __dremf(float __x, float __y) throw(); 
# 189
extern float significandf(float __x) throw(); extern float __significandf(float __x) throw(); 
# 196
extern float copysignf(float __x, float __y) throw() __attribute((const)); extern float __copysignf(float __x, float __y) throw() __attribute((const)); 
# 201
extern float nanf(const char * __tagb) throw(); extern float __nanf(const char * __tagb) throw(); 
# 211 "/usr/include/x86_64-linux-gnu/bits/mathcalls.h" 3
extern int isnanf(float __value) throw() __attribute((const)); 
# 217
extern float j0f(float) throw(); extern float __j0f(float) throw(); 
# 218
extern float j1f(float) throw(); extern float __j1f(float) throw(); 
# 219
extern float jnf(int, float) throw(); extern float __jnf(int, float) throw(); 
# 220
extern float y0f(float) throw(); extern float __y0f(float) throw(); 
# 221
extern float y1f(float) throw(); extern float __y1f(float) throw(); 
# 222
extern float ynf(int, float) throw(); extern float __ynf(int, float) throw(); 
# 228
extern float erff(float) throw(); extern float __erff(float) throw(); 
# 229
extern float erfcf(float) throw(); extern float __erfcf(float) throw(); 
# 230
extern float lgammaf(float) throw(); extern float __lgammaf(float) throw(); 
# 235
extern float tgammaf(float) throw(); extern float __tgammaf(float) throw(); 
# 241
extern float gammaf(float) throw(); extern float __gammaf(float) throw(); 
# 249
extern float lgammaf_r(float, int * __signgamp) throw(); extern float __lgammaf_r(float, int * __signgamp) throw(); 
# 256
extern float rintf(float __x) throw(); extern float __rintf(float __x) throw(); 
# 259
extern float nextafterf(float __x, float __y) throw(); extern float __nextafterf(float __x, float __y) throw(); 
# 261
extern float nexttowardf(float __x, long double __y) throw(); extern float __nexttowardf(float __x, long double __y) throw(); 
# 266
extern float nextdownf(float __x) throw(); extern float __nextdownf(float __x) throw(); 
# 268
extern float nextupf(float __x) throw(); extern float __nextupf(float __x) throw(); 
# 272
extern float remainderf(float __x, float __y) throw(); extern float __remainderf(float __x, float __y) throw(); 
# 276
extern float scalbnf(float __x, int __n) throw(); extern float __scalbnf(float __x, int __n) throw(); 
# 280
extern int ilogbf(float __x) throw(); extern int __ilogbf(float __x) throw(); 
# 285
extern long llogbf(float __x) throw(); extern long __llogbf(float __x) throw(); 
# 290
extern float scalblnf(float __x, long __n) throw(); extern float __scalblnf(float __x, long __n) throw(); 
# 294
extern float nearbyintf(float __x) throw(); extern float __nearbyintf(float __x) throw(); 
# 298
extern float roundf(float __x) throw() __attribute((const)); extern float __roundf(float __x) throw() __attribute((const)); 
# 302
extern float truncf(float __x) throw() __attribute((const)); extern float __truncf(float __x) throw() __attribute((const)); 
# 307
extern float remquof(float __x, float __y, int * __quo) throw(); extern float __remquof(float __x, float __y, int * __quo) throw(); 
# 314
extern long lrintf(float __x) throw(); extern long __lrintf(float __x) throw(); 
# 316
__extension__ extern long long llrintf(float __x) throw(); extern long long __llrintf(float __x) throw(); 
# 320
extern long lroundf(float __x) throw(); extern long __lroundf(float __x) throw(); 
# 322
__extension__ extern long long llroundf(float __x) throw(); extern long long __llroundf(float __x) throw(); 
# 326
extern float fdimf(float __x, float __y) throw(); extern float __fdimf(float __x, float __y) throw(); 
# 329
extern float fmaxf(float __x, float __y) throw() __attribute((const)); extern float __fmaxf(float __x, float __y) throw() __attribute((const)); 
# 332
extern float fminf(float __x, float __y) throw() __attribute((const)); extern float __fminf(float __x, float __y) throw() __attribute((const)); 
# 335
extern float fmaf(float __x, float __y, float __z) throw(); extern float __fmaf(float __x, float __y, float __z) throw(); 
# 340
extern float roundevenf(float __x) throw() __attribute((const)); extern float __roundevenf(float __x) throw() __attribute((const)); 
# 344
extern __intmax_t fromfpf(float __x, int __round, unsigned __width) throw(); extern __intmax_t __fromfpf(float __x, int __round, unsigned __width) throw(); 
# 349
extern __uintmax_t ufromfpf(float __x, int __round, unsigned __width) throw(); extern __uintmax_t __ufromfpf(float __x, int __round, unsigned __width) throw(); 
# 355
extern __intmax_t fromfpxf(float __x, int __round, unsigned __width) throw(); extern __intmax_t __fromfpxf(float __x, int __round, unsigned __width) throw(); 
# 361
extern __uintmax_t ufromfpxf(float __x, int __round, unsigned __width) throw(); extern __uintmax_t __ufromfpxf(float __x, int __round, unsigned __width) throw(); 
# 365
extern float fmaxmagf(float __x, float __y) throw() __attribute((const)); extern float __fmaxmagf(float __x, float __y) throw() __attribute((const)); 
# 368
extern float fminmagf(float __x, float __y) throw() __attribute((const)); extern float __fminmagf(float __x, float __y) throw() __attribute((const)); 
# 371
extern int canonicalizef(float * __cx, const float * __x) throw(); 
# 376
extern int totalorderf(const float * __x, const float * __y) throw()
# 378
 __attribute((__pure__)); 
# 381
extern int totalordermagf(const float * __x, const float * __y) throw()
# 383
 __attribute((__pure__)); 
# 386
extern float getpayloadf(const float * __x) throw(); extern float __getpayloadf(const float * __x) throw(); 
# 389
extern int setpayloadf(float * __x, float __payload) throw(); 
# 392
extern int setpayloadsigf(float * __x, float __payload) throw(); 
# 400
extern float scalbf(float __x, float __n) throw(); extern float __scalbf(float __x, float __n) throw(); 
# 21 "/usr/include/x86_64-linux-gnu/bits/mathcalls-helper-functions.h" 3
extern int __fpclassifyl(long double __value) throw()
# 22
 __attribute((const)); 
# 25
extern int __signbitl(long double __value) throw()
# 26
 __attribute((const)); 
# 30
extern int __isinfl(long double __value) throw() __attribute((const)); 
# 33
extern int __finitel(long double __value) throw() __attribute((const)); 
# 36
extern int __isnanl(long double __value) throw() __attribute((const)); 
# 39
extern int __iseqsigl(long double __x, long double __y) throw(); 
# 42
extern int __issignalingl(long double __value) throw()
# 43
 __attribute((const)); 
# 53 "/usr/include/x86_64-linux-gnu/bits/mathcalls.h" 3
extern long double acosl(long double __x) throw(); extern long double __acosl(long double __x) throw(); 
# 55
extern long double asinl(long double __x) throw(); extern long double __asinl(long double __x) throw(); 
# 57
extern long double atanl(long double __x) throw(); extern long double __atanl(long double __x) throw(); 
# 59
extern long double atan2l(long double __y, long double __x) throw(); extern long double __atan2l(long double __y, long double __x) throw(); 
# 62
extern long double cosl(long double __x) throw(); extern long double __cosl(long double __x) throw(); 
# 64
extern long double sinl(long double __x) throw(); extern long double __sinl(long double __x) throw(); 
# 66
extern long double tanl(long double __x) throw(); extern long double __tanl(long double __x) throw(); 
# 71
extern long double coshl(long double __x) throw(); extern long double __coshl(long double __x) throw(); 
# 73
extern long double sinhl(long double __x) throw(); extern long double __sinhl(long double __x) throw(); 
# 75
extern long double tanhl(long double __x) throw(); extern long double __tanhl(long double __x) throw(); 
# 79
extern void sincosl(long double __x, long double * __sinx, long double * __cosx) throw(); extern void __sincosl(long double __x, long double * __sinx, long double * __cosx) throw(); 
# 85
extern long double acoshl(long double __x) throw(); extern long double __acoshl(long double __x) throw(); 
# 87
extern long double asinhl(long double __x) throw(); extern long double __asinhl(long double __x) throw(); 
# 89
extern long double atanhl(long double __x) throw(); extern long double __atanhl(long double __x) throw(); 
# 95
extern long double expl(long double __x) throw(); extern long double __expl(long double __x) throw(); 
# 98
extern long double frexpl(long double __x, int * __exponent) throw(); extern long double __frexpl(long double __x, int * __exponent) throw(); 
# 101
extern long double ldexpl(long double __x, int __exponent) throw(); extern long double __ldexpl(long double __x, int __exponent) throw(); 
# 104
extern long double logl(long double __x) throw(); extern long double __logl(long double __x) throw(); 
# 107
extern long double log10l(long double __x) throw(); extern long double __log10l(long double __x) throw(); 
# 110
extern long double modfl(long double __x, long double * __iptr) throw(); extern long double __modfl(long double __x, long double * __iptr) throw() __attribute((__nonnull__(2))); 
# 114
extern long double exp10l(long double __x) throw(); extern long double __exp10l(long double __x) throw(); 
# 119
extern long double expm1l(long double __x) throw(); extern long double __expm1l(long double __x) throw(); 
# 122
extern long double log1pl(long double __x) throw(); extern long double __log1pl(long double __x) throw(); 
# 125
extern long double logbl(long double __x) throw(); extern long double __logbl(long double __x) throw(); 
# 130
extern long double exp2l(long double __x) throw(); extern long double __exp2l(long double __x) throw(); 
# 133
extern long double log2l(long double __x) throw(); extern long double __log2l(long double __x) throw(); 
# 140
extern long double powl(long double __x, long double __y) throw(); extern long double __powl(long double __x, long double __y) throw(); 
# 143
extern long double sqrtl(long double __x) throw(); extern long double __sqrtl(long double __x) throw(); 
# 147
extern long double hypotl(long double __x, long double __y) throw(); extern long double __hypotl(long double __x, long double __y) throw(); 
# 152
extern long double cbrtl(long double __x) throw(); extern long double __cbrtl(long double __x) throw(); 
# 159
extern long double ceill(long double __x) throw() __attribute((const)); extern long double __ceill(long double __x) throw() __attribute((const)); 
# 162
extern long double fabsl(long double __x) throw() __attribute((const)); extern long double __fabsl(long double __x) throw() __attribute((const)); 
# 165
extern long double floorl(long double __x) throw() __attribute((const)); extern long double __floorl(long double __x) throw() __attribute((const)); 
# 168
extern long double fmodl(long double __x, long double __y) throw(); extern long double __fmodl(long double __x, long double __y) throw(); 
# 177 "/usr/include/x86_64-linux-gnu/bits/mathcalls.h" 3
extern int isinfl(long double __value) throw() __attribute((const)); 
# 182
extern int finitel(long double __value) throw() __attribute((const)); 
# 185
extern long double dreml(long double __x, long double __y) throw(); extern long double __dreml(long double __x, long double __y) throw(); 
# 189
extern long double significandl(long double __x) throw(); extern long double __significandl(long double __x) throw(); 
# 196
extern long double copysignl(long double __x, long double __y) throw() __attribute((const)); extern long double __copysignl(long double __x, long double __y) throw() __attribute((const)); 
# 201
extern long double nanl(const char * __tagb) throw(); extern long double __nanl(const char * __tagb) throw(); 
# 211 "/usr/include/x86_64-linux-gnu/bits/mathcalls.h" 3
extern int isnanl(long double __value) throw() __attribute((const)); 
# 217
extern long double j0l(long double) throw(); extern long double __j0l(long double) throw(); 
# 218
extern long double j1l(long double) throw(); extern long double __j1l(long double) throw(); 
# 219
extern long double jnl(int, long double) throw(); extern long double __jnl(int, long double) throw(); 
# 220
extern long double y0l(long double) throw(); extern long double __y0l(long double) throw(); 
# 221
extern long double y1l(long double) throw(); extern long double __y1l(long double) throw(); 
# 222
extern long double ynl(int, long double) throw(); extern long double __ynl(int, long double) throw(); 
# 228
extern long double erfl(long double) throw(); extern long double __erfl(long double) throw(); 
# 229
extern long double erfcl(long double) throw(); extern long double __erfcl(long double) throw(); 
# 230
extern long double lgammal(long double) throw(); extern long double __lgammal(long double) throw(); 
# 235
extern long double tgammal(long double) throw(); extern long double __tgammal(long double) throw(); 
# 241
extern long double gammal(long double) throw(); extern long double __gammal(long double) throw(); 
# 249
extern long double lgammal_r(long double, int * __signgamp) throw(); extern long double __lgammal_r(long double, int * __signgamp) throw(); 
# 256
extern long double rintl(long double __x) throw(); extern long double __rintl(long double __x) throw(); 
# 259
extern long double nextafterl(long double __x, long double __y) throw(); extern long double __nextafterl(long double __x, long double __y) throw(); 
# 261
extern long double nexttowardl(long double __x, long double __y) throw(); extern long double __nexttowardl(long double __x, long double __y) throw(); 
# 266
extern long double nextdownl(long double __x) throw(); extern long double __nextdownl(long double __x) throw(); 
# 268
extern long double nextupl(long double __x) throw(); extern long double __nextupl(long double __x) throw(); 
# 272
extern long double remainderl(long double __x, long double __y) throw(); extern long double __remainderl(long double __x, long double __y) throw(); 
# 276
extern long double scalbnl(long double __x, int __n) throw(); extern long double __scalbnl(long double __x, int __n) throw(); 
# 280
extern int ilogbl(long double __x) throw(); extern int __ilogbl(long double __x) throw(); 
# 285
extern long llogbl(long double __x) throw(); extern long __llogbl(long double __x) throw(); 
# 290
extern long double scalblnl(long double __x, long __n) throw(); extern long double __scalblnl(long double __x, long __n) throw(); 
# 294
extern long double nearbyintl(long double __x) throw(); extern long double __nearbyintl(long double __x) throw(); 
# 298
extern long double roundl(long double __x) throw() __attribute((const)); extern long double __roundl(long double __x) throw() __attribute((const)); 
# 302
extern long double truncl(long double __x) throw() __attribute((const)); extern long double __truncl(long double __x) throw() __attribute((const)); 
# 307
extern long double remquol(long double __x, long double __y, int * __quo) throw(); extern long double __remquol(long double __x, long double __y, int * __quo) throw(); 
# 314
extern long lrintl(long double __x) throw(); extern long __lrintl(long double __x) throw(); 
# 316
__extension__ extern long long llrintl(long double __x) throw(); extern long long __llrintl(long double __x) throw(); 
# 320
extern long lroundl(long double __x) throw(); extern long __lroundl(long double __x) throw(); 
# 322
__extension__ extern long long llroundl(long double __x) throw(); extern long long __llroundl(long double __x) throw(); 
# 326
extern long double fdiml(long double __x, long double __y) throw(); extern long double __fdiml(long double __x, long double __y) throw(); 
# 329
extern long double fmaxl(long double __x, long double __y) throw() __attribute((const)); extern long double __fmaxl(long double __x, long double __y) throw() __attribute((const)); 
# 332
extern long double fminl(long double __x, long double __y) throw() __attribute((const)); extern long double __fminl(long double __x, long double __y) throw() __attribute((const)); 
# 335
extern long double fmal(long double __x, long double __y, long double __z) throw(); extern long double __fmal(long double __x, long double __y, long double __z) throw(); 
# 340
extern long double roundevenl(long double __x) throw() __attribute((const)); extern long double __roundevenl(long double __x) throw() __attribute((const)); 
# 344
extern __intmax_t fromfpl(long double __x, int __round, unsigned __width) throw(); extern __intmax_t __fromfpl(long double __x, int __round, unsigned __width) throw(); 
# 349
extern __uintmax_t ufromfpl(long double __x, int __round, unsigned __width) throw(); extern __uintmax_t __ufromfpl(long double __x, int __round, unsigned __width) throw(); 
# 355
extern __intmax_t fromfpxl(long double __x, int __round, unsigned __width) throw(); extern __intmax_t __fromfpxl(long double __x, int __round, unsigned __width) throw(); 
# 361
extern __uintmax_t ufromfpxl(long double __x, int __round, unsigned __width) throw(); extern __uintmax_t __ufromfpxl(long double __x, int __round, unsigned __width) throw(); 
# 365
extern long double fmaxmagl(long double __x, long double __y) throw() __attribute((const)); extern long double __fmaxmagl(long double __x, long double __y) throw() __attribute((const)); 
# 368
extern long double fminmagl(long double __x, long double __y) throw() __attribute((const)); extern long double __fminmagl(long double __x, long double __y) throw() __attribute((const)); 
# 371
extern int canonicalizel(long double * __cx, const long double * __x) throw(); 
# 376
extern int totalorderl(const long double * __x, const long double * __y) throw()
# 378
 __attribute((__pure__)); 
# 381
extern int totalordermagl(const long double * __x, const long double * __y) throw()
# 383
 __attribute((__pure__)); 
# 386
extern long double getpayloadl(const long double * __x) throw(); extern long double __getpayloadl(const long double * __x) throw(); 
# 389
extern int setpayloadl(long double * __x, long double __payload) throw(); 
# 392
extern int setpayloadsigl(long double * __x, long double __payload) throw(); 
# 400
extern long double scalbl(long double __x, long double __n) throw(); extern long double __scalbl(long double __x, long double __n) throw(); 
# 53 "/usr/include/x86_64-linux-gnu/bits/mathcalls.h" 3
extern _Float32 acosf32(_Float32 __x) throw(); extern _Float32 __acosf32(_Float32 __x) throw(); 
# 55
extern _Float32 asinf32(_Float32 __x) throw(); extern _Float32 __asinf32(_Float32 __x) throw(); 
# 57
extern _Float32 atanf32(_Float32 __x) throw(); extern _Float32 __atanf32(_Float32 __x) throw(); 
# 59
extern _Float32 atan2f32(_Float32 __y, _Float32 __x) throw(); extern _Float32 __atan2f32(_Float32 __y, _Float32 __x) throw(); 
# 62
extern _Float32 cosf32(_Float32 __x) throw(); extern _Float32 __cosf32(_Float32 __x) throw(); 
# 64
extern _Float32 sinf32(_Float32 __x) throw(); extern _Float32 __sinf32(_Float32 __x) throw(); 
# 66
extern _Float32 tanf32(_Float32 __x) throw(); extern _Float32 __tanf32(_Float32 __x) throw(); 
# 71
extern _Float32 coshf32(_Float32 __x) throw(); extern _Float32 __coshf32(_Float32 __x) throw(); 
# 73
extern _Float32 sinhf32(_Float32 __x) throw(); extern _Float32 __sinhf32(_Float32 __x) throw(); 
# 75
extern _Float32 tanhf32(_Float32 __x) throw(); extern _Float32 __tanhf32(_Float32 __x) throw(); 
# 79
extern void sincosf32(_Float32 __x, _Float32 * __sinx, _Float32 * __cosx) throw(); extern void __sincosf32(_Float32 __x, _Float32 * __sinx, _Float32 * __cosx) throw(); 
# 85
extern _Float32 acoshf32(_Float32 __x) throw(); extern _Float32 __acoshf32(_Float32 __x) throw(); 
# 87
extern _Float32 asinhf32(_Float32 __x) throw(); extern _Float32 __asinhf32(_Float32 __x) throw(); 
# 89
extern _Float32 atanhf32(_Float32 __x) throw(); extern _Float32 __atanhf32(_Float32 __x) throw(); 
# 95
extern _Float32 expf32(_Float32 __x) throw(); extern _Float32 __expf32(_Float32 __x) throw(); 
# 98
extern _Float32 frexpf32(_Float32 __x, int * __exponent) throw(); extern _Float32 __frexpf32(_Float32 __x, int * __exponent) throw(); 
# 101
extern _Float32 ldexpf32(_Float32 __x, int __exponent) throw(); extern _Float32 __ldexpf32(_Float32 __x, int __exponent) throw(); 
# 104
extern _Float32 logf32(_Float32 __x) throw(); extern _Float32 __logf32(_Float32 __x) throw(); 
# 107
extern _Float32 log10f32(_Float32 __x) throw(); extern _Float32 __log10f32(_Float32 __x) throw(); 
# 110
extern _Float32 modff32(_Float32 __x, _Float32 * __iptr) throw(); extern _Float32 __modff32(_Float32 __x, _Float32 * __iptr) throw() __attribute((__nonnull__(2))); 
# 114
extern _Float32 exp10f32(_Float32 __x) throw(); extern _Float32 __exp10f32(_Float32 __x) throw(); 
# 119
extern _Float32 expm1f32(_Float32 __x) throw(); extern _Float32 __expm1f32(_Float32 __x) throw(); 
# 122
extern _Float32 log1pf32(_Float32 __x) throw(); extern _Float32 __log1pf32(_Float32 __x) throw(); 
# 125
extern _Float32 logbf32(_Float32 __x) throw(); extern _Float32 __logbf32(_Float32 __x) throw(); 
# 130
extern _Float32 exp2f32(_Float32 __x) throw(); extern _Float32 __exp2f32(_Float32 __x) throw(); 
# 133
extern _Float32 log2f32(_Float32 __x) throw(); extern _Float32 __log2f32(_Float32 __x) throw(); 
# 140
extern _Float32 powf32(_Float32 __x, _Float32 __y) throw(); extern _Float32 __powf32(_Float32 __x, _Float32 __y) throw(); 
# 143
extern _Float32 sqrtf32(_Float32 __x) throw(); extern _Float32 __sqrtf32(_Float32 __x) throw(); 
# 147
extern _Float32 hypotf32(_Float32 __x, _Float32 __y) throw(); extern _Float32 __hypotf32(_Float32 __x, _Float32 __y) throw(); 
# 152
extern _Float32 cbrtf32(_Float32 __x) throw(); extern _Float32 __cbrtf32(_Float32 __x) throw(); 
# 159
extern _Float32 ceilf32(_Float32 __x) throw() __attribute((const)); extern _Float32 __ceilf32(_Float32 __x) throw() __attribute((const)); 
# 162
extern _Float32 fabsf32(_Float32 __x) throw() __attribute((const)); extern _Float32 __fabsf32(_Float32 __x) throw() __attribute((const)); 
# 165
extern _Float32 floorf32(_Float32 __x) throw() __attribute((const)); extern _Float32 __floorf32(_Float32 __x) throw() __attribute((const)); 
# 168
extern _Float32 fmodf32(_Float32 __x, _Float32 __y) throw(); extern _Float32 __fmodf32(_Float32 __x, _Float32 __y) throw(); 
# 196 "/usr/include/x86_64-linux-gnu/bits/mathcalls.h" 3
extern _Float32 copysignf32(_Float32 __x, _Float32 __y) throw() __attribute((const)); extern _Float32 __copysignf32(_Float32 __x, _Float32 __y) throw() __attribute((const)); 
# 201
extern _Float32 nanf32(const char * __tagb) throw(); extern _Float32 __nanf32(const char * __tagb) throw(); 
# 217 "/usr/include/x86_64-linux-gnu/bits/mathcalls.h" 3
extern _Float32 j0f32(_Float32) throw(); extern _Float32 __j0f32(_Float32) throw(); 
# 218
extern _Float32 j1f32(_Float32) throw(); extern _Float32 __j1f32(_Float32) throw(); 
# 219
extern _Float32 jnf32(int, _Float32) throw(); extern _Float32 __jnf32(int, _Float32) throw(); 
# 220
extern _Float32 y0f32(_Float32) throw(); extern _Float32 __y0f32(_Float32) throw(); 
# 221
extern _Float32 y1f32(_Float32) throw(); extern _Float32 __y1f32(_Float32) throw(); 
# 222
extern _Float32 ynf32(int, _Float32) throw(); extern _Float32 __ynf32(int, _Float32) throw(); 
# 228
extern _Float32 erff32(_Float32) throw(); extern _Float32 __erff32(_Float32) throw(); 
# 229
extern _Float32 erfcf32(_Float32) throw(); extern _Float32 __erfcf32(_Float32) throw(); 
# 230
extern _Float32 lgammaf32(_Float32) throw(); extern _Float32 __lgammaf32(_Float32) throw(); 
# 235
extern _Float32 tgammaf32(_Float32) throw(); extern _Float32 __tgammaf32(_Float32) throw(); 
# 249 "/usr/include/x86_64-linux-gnu/bits/mathcalls.h" 3
extern _Float32 lgammaf32_r(_Float32, int * __signgamp) throw(); extern _Float32 __lgammaf32_r(_Float32, int * __signgamp) throw(); 
# 256
extern _Float32 rintf32(_Float32 __x) throw(); extern _Float32 __rintf32(_Float32 __x) throw(); 
# 259
extern _Float32 nextafterf32(_Float32 __x, _Float32 __y) throw(); extern _Float32 __nextafterf32(_Float32 __x, _Float32 __y) throw(); 
# 266
extern _Float32 nextdownf32(_Float32 __x) throw(); extern _Float32 __nextdownf32(_Float32 __x) throw(); 
# 268
extern _Float32 nextupf32(_Float32 __x) throw(); extern _Float32 __nextupf32(_Float32 __x) throw(); 
# 272
extern _Float32 remainderf32(_Float32 __x, _Float32 __y) throw(); extern _Float32 __remainderf32(_Float32 __x, _Float32 __y) throw(); 
# 276
extern _Float32 scalbnf32(_Float32 __x, int __n) throw(); extern _Float32 __scalbnf32(_Float32 __x, int __n) throw(); 
# 280
extern int ilogbf32(_Float32 __x) throw(); extern int __ilogbf32(_Float32 __x) throw(); 
# 285
extern long llogbf32(_Float32 __x) throw(); extern long __llogbf32(_Float32 __x) throw(); 
# 290
extern _Float32 scalblnf32(_Float32 __x, long __n) throw(); extern _Float32 __scalblnf32(_Float32 __x, long __n) throw(); 
# 294
extern _Float32 nearbyintf32(_Float32 __x) throw(); extern _Float32 __nearbyintf32(_Float32 __x) throw(); 
# 298
extern _Float32 roundf32(_Float32 __x) throw() __attribute((const)); extern _Float32 __roundf32(_Float32 __x) throw() __attribute((const)); 
# 302
extern _Float32 truncf32(_Float32 __x) throw() __attribute((const)); extern _Float32 __truncf32(_Float32 __x) throw() __attribute((const)); 
# 307
extern _Float32 remquof32(_Float32 __x, _Float32 __y, int * __quo) throw(); extern _Float32 __remquof32(_Float32 __x, _Float32 __y, int * __quo) throw(); 
# 314
extern long lrintf32(_Float32 __x) throw(); extern long __lrintf32(_Float32 __x) throw(); 
# 316
__extension__ extern long long llrintf32(_Float32 __x) throw(); extern long long __llrintf32(_Float32 __x) throw(); 
# 320
extern long lroundf32(_Float32 __x) throw(); extern long __lroundf32(_Float32 __x) throw(); 
# 322
__extension__ extern long long llroundf32(_Float32 __x) throw(); extern long long __llroundf32(_Float32 __x) throw(); 
# 326
extern _Float32 fdimf32(_Float32 __x, _Float32 __y) throw(); extern _Float32 __fdimf32(_Float32 __x, _Float32 __y) throw(); 
# 329
extern _Float32 fmaxf32(_Float32 __x, _Float32 __y) throw() __attribute((const)); extern _Float32 __fmaxf32(_Float32 __x, _Float32 __y) throw() __attribute((const)); 
# 332
extern _Float32 fminf32(_Float32 __x, _Float32 __y) throw() __attribute((const)); extern _Float32 __fminf32(_Float32 __x, _Float32 __y) throw() __attribute((const)); 
# 335
extern _Float32 fmaf32(_Float32 __x, _Float32 __y, _Float32 __z) throw(); extern _Float32 __fmaf32(_Float32 __x, _Float32 __y, _Float32 __z) throw(); 
# 340
extern _Float32 roundevenf32(_Float32 __x) throw() __attribute((const)); extern _Float32 __roundevenf32(_Float32 __x) throw() __attribute((const)); 
# 344
extern __intmax_t fromfpf32(_Float32 __x, int __round, unsigned __width) throw(); extern __intmax_t __fromfpf32(_Float32 __x, int __round, unsigned __width) throw(); 
# 349
extern __uintmax_t ufromfpf32(_Float32 __x, int __round, unsigned __width) throw(); extern __uintmax_t __ufromfpf32(_Float32 __x, int __round, unsigned __width) throw(); 
# 355
extern __intmax_t fromfpxf32(_Float32 __x, int __round, unsigned __width) throw(); extern __intmax_t __fromfpxf32(_Float32 __x, int __round, unsigned __width) throw(); 
# 361
extern __uintmax_t ufromfpxf32(_Float32 __x, int __round, unsigned __width) throw(); extern __uintmax_t __ufromfpxf32(_Float32 __x, int __round, unsigned __width) throw(); 
# 365
extern _Float32 fmaxmagf32(_Float32 __x, _Float32 __y) throw() __attribute((const)); extern _Float32 __fmaxmagf32(_Float32 __x, _Float32 __y) throw() __attribute((const)); 
# 368
extern _Float32 fminmagf32(_Float32 __x, _Float32 __y) throw() __attribute((const)); extern _Float32 __fminmagf32(_Float32 __x, _Float32 __y) throw() __attribute((const)); 
# 371
extern int canonicalizef32(_Float32 * __cx, const _Float32 * __x) throw(); 
# 376
extern int totalorderf32(const _Float32 * __x, const _Float32 * __y) throw()
# 378
 __attribute((__pure__)); 
# 381
extern int totalordermagf32(const _Float32 * __x, const _Float32 * __y) throw()
# 383
 __attribute((__pure__)); 
# 386
extern _Float32 getpayloadf32(const _Float32 * __x) throw(); extern _Float32 __getpayloadf32(const _Float32 * __x) throw(); 
# 389
extern int setpayloadf32(_Float32 * __x, _Float32 __payload) throw(); 
# 392
extern int setpayloadsigf32(_Float32 * __x, _Float32 __payload) throw(); 
# 53 "/usr/include/x86_64-linux-gnu/bits/mathcalls.h" 3
extern _Float64 acosf64(_Float64 __x) throw(); extern _Float64 __acosf64(_Float64 __x) throw(); 
# 55
extern _Float64 asinf64(_Float64 __x) throw(); extern _Float64 __asinf64(_Float64 __x) throw(); 
# 57
extern _Float64 atanf64(_Float64 __x) throw(); extern _Float64 __atanf64(_Float64 __x) throw(); 
# 59
extern _Float64 atan2f64(_Float64 __y, _Float64 __x) throw(); extern _Float64 __atan2f64(_Float64 __y, _Float64 __x) throw(); 
# 62
extern _Float64 cosf64(_Float64 __x) throw(); extern _Float64 __cosf64(_Float64 __x) throw(); 
# 64
extern _Float64 sinf64(_Float64 __x) throw(); extern _Float64 __sinf64(_Float64 __x) throw(); 
# 66
extern _Float64 tanf64(_Float64 __x) throw(); extern _Float64 __tanf64(_Float64 __x) throw(); 
# 71
extern _Float64 coshf64(_Float64 __x) throw(); extern _Float64 __coshf64(_Float64 __x) throw(); 
# 73
extern _Float64 sinhf64(_Float64 __x) throw(); extern _Float64 __sinhf64(_Float64 __x) throw(); 
# 75
extern _Float64 tanhf64(_Float64 __x) throw(); extern _Float64 __tanhf64(_Float64 __x) throw(); 
# 79
extern void sincosf64(_Float64 __x, _Float64 * __sinx, _Float64 * __cosx) throw(); extern void __sincosf64(_Float64 __x, _Float64 * __sinx, _Float64 * __cosx) throw(); 
# 85
extern _Float64 acoshf64(_Float64 __x) throw(); extern _Float64 __acoshf64(_Float64 __x) throw(); 
# 87
extern _Float64 asinhf64(_Float64 __x) throw(); extern _Float64 __asinhf64(_Float64 __x) throw(); 
# 89
extern _Float64 atanhf64(_Float64 __x) throw(); extern _Float64 __atanhf64(_Float64 __x) throw(); 
# 95
extern _Float64 expf64(_Float64 __x) throw(); extern _Float64 __expf64(_Float64 __x) throw(); 
# 98
extern _Float64 frexpf64(_Float64 __x, int * __exponent) throw(); extern _Float64 __frexpf64(_Float64 __x, int * __exponent) throw(); 
# 101
extern _Float64 ldexpf64(_Float64 __x, int __exponent) throw(); extern _Float64 __ldexpf64(_Float64 __x, int __exponent) throw(); 
# 104
extern _Float64 logf64(_Float64 __x) throw(); extern _Float64 __logf64(_Float64 __x) throw(); 
# 107
extern _Float64 log10f64(_Float64 __x) throw(); extern _Float64 __log10f64(_Float64 __x) throw(); 
# 110
extern _Float64 modff64(_Float64 __x, _Float64 * __iptr) throw(); extern _Float64 __modff64(_Float64 __x, _Float64 * __iptr) throw() __attribute((__nonnull__(2))); 
# 114
extern _Float64 exp10f64(_Float64 __x) throw(); extern _Float64 __exp10f64(_Float64 __x) throw(); 
# 119
extern _Float64 expm1f64(_Float64 __x) throw(); extern _Float64 __expm1f64(_Float64 __x) throw(); 
# 122
extern _Float64 log1pf64(_Float64 __x) throw(); extern _Float64 __log1pf64(_Float64 __x) throw(); 
# 125
extern _Float64 logbf64(_Float64 __x) throw(); extern _Float64 __logbf64(_Float64 __x) throw(); 
# 130
extern _Float64 exp2f64(_Float64 __x) throw(); extern _Float64 __exp2f64(_Float64 __x) throw(); 
# 133
extern _Float64 log2f64(_Float64 __x) throw(); extern _Float64 __log2f64(_Float64 __x) throw(); 
# 140
extern _Float64 powf64(_Float64 __x, _Float64 __y) throw(); extern _Float64 __powf64(_Float64 __x, _Float64 __y) throw(); 
# 143
extern _Float64 sqrtf64(_Float64 __x) throw(); extern _Float64 __sqrtf64(_Float64 __x) throw(); 
# 147
extern _Float64 hypotf64(_Float64 __x, _Float64 __y) throw(); extern _Float64 __hypotf64(_Float64 __x, _Float64 __y) throw(); 
# 152
extern _Float64 cbrtf64(_Float64 __x) throw(); extern _Float64 __cbrtf64(_Float64 __x) throw(); 
# 159
extern _Float64 ceilf64(_Float64 __x) throw() __attribute((const)); extern _Float64 __ceilf64(_Float64 __x) throw() __attribute((const)); 
# 162
extern _Float64 fabsf64(_Float64 __x) throw() __attribute((const)); extern _Float64 __fabsf64(_Float64 __x) throw() __attribute((const)); 
# 165
extern _Float64 floorf64(_Float64 __x) throw() __attribute((const)); extern _Float64 __floorf64(_Float64 __x) throw() __attribute((const)); 
# 168
extern _Float64 fmodf64(_Float64 __x, _Float64 __y) throw(); extern _Float64 __fmodf64(_Float64 __x, _Float64 __y) throw(); 
# 196 "/usr/include/x86_64-linux-gnu/bits/mathcalls.h" 3
extern _Float64 copysignf64(_Float64 __x, _Float64 __y) throw() __attribute((const)); extern _Float64 __copysignf64(_Float64 __x, _Float64 __y) throw() __attribute((const)); 
# 201
extern _Float64 nanf64(const char * __tagb) throw(); extern _Float64 __nanf64(const char * __tagb) throw(); 
# 217 "/usr/include/x86_64-linux-gnu/bits/mathcalls.h" 3
extern _Float64 j0f64(_Float64) throw(); extern _Float64 __j0f64(_Float64) throw(); 
# 218
extern _Float64 j1f64(_Float64) throw(); extern _Float64 __j1f64(_Float64) throw(); 
# 219
extern _Float64 jnf64(int, _Float64) throw(); extern _Float64 __jnf64(int, _Float64) throw(); 
# 220
extern _Float64 y0f64(_Float64) throw(); extern _Float64 __y0f64(_Float64) throw(); 
# 221
extern _Float64 y1f64(_Float64) throw(); extern _Float64 __y1f64(_Float64) throw(); 
# 222
extern _Float64 ynf64(int, _Float64) throw(); extern _Float64 __ynf64(int, _Float64) throw(); 
# 228
extern _Float64 erff64(_Float64) throw(); extern _Float64 __erff64(_Float64) throw(); 
# 229
extern _Float64 erfcf64(_Float64) throw(); extern _Float64 __erfcf64(_Float64) throw(); 
# 230
extern _Float64 lgammaf64(_Float64) throw(); extern _Float64 __lgammaf64(_Float64) throw(); 
# 235
extern _Float64 tgammaf64(_Float64) throw(); extern _Float64 __tgammaf64(_Float64) throw(); 
# 249 "/usr/include/x86_64-linux-gnu/bits/mathcalls.h" 3
extern _Float64 lgammaf64_r(_Float64, int * __signgamp) throw(); extern _Float64 __lgammaf64_r(_Float64, int * __signgamp) throw(); 
# 256
extern _Float64 rintf64(_Float64 __x) throw(); extern _Float64 __rintf64(_Float64 __x) throw(); 
# 259
extern _Float64 nextafterf64(_Float64 __x, _Float64 __y) throw(); extern _Float64 __nextafterf64(_Float64 __x, _Float64 __y) throw(); 
# 266
extern _Float64 nextdownf64(_Float64 __x) throw(); extern _Float64 __nextdownf64(_Float64 __x) throw(); 
# 268
extern _Float64 nextupf64(_Float64 __x) throw(); extern _Float64 __nextupf64(_Float64 __x) throw(); 
# 272
extern _Float64 remainderf64(_Float64 __x, _Float64 __y) throw(); extern _Float64 __remainderf64(_Float64 __x, _Float64 __y) throw(); 
# 276
extern _Float64 scalbnf64(_Float64 __x, int __n) throw(); extern _Float64 __scalbnf64(_Float64 __x, int __n) throw(); 
# 280
extern int ilogbf64(_Float64 __x) throw(); extern int __ilogbf64(_Float64 __x) throw(); 
# 285
extern long llogbf64(_Float64 __x) throw(); extern long __llogbf64(_Float64 __x) throw(); 
# 290
extern _Float64 scalblnf64(_Float64 __x, long __n) throw(); extern _Float64 __scalblnf64(_Float64 __x, long __n) throw(); 
# 294
extern _Float64 nearbyintf64(_Float64 __x) throw(); extern _Float64 __nearbyintf64(_Float64 __x) throw(); 
# 298
extern _Float64 roundf64(_Float64 __x) throw() __attribute((const)); extern _Float64 __roundf64(_Float64 __x) throw() __attribute((const)); 
# 302
extern _Float64 truncf64(_Float64 __x) throw() __attribute((const)); extern _Float64 __truncf64(_Float64 __x) throw() __attribute((const)); 
# 307
extern _Float64 remquof64(_Float64 __x, _Float64 __y, int * __quo) throw(); extern _Float64 __remquof64(_Float64 __x, _Float64 __y, int * __quo) throw(); 
# 314
extern long lrintf64(_Float64 __x) throw(); extern long __lrintf64(_Float64 __x) throw(); 
# 316
__extension__ extern long long llrintf64(_Float64 __x) throw(); extern long long __llrintf64(_Float64 __x) throw(); 
# 320
extern long lroundf64(_Float64 __x) throw(); extern long __lroundf64(_Float64 __x) throw(); 
# 322
__extension__ extern long long llroundf64(_Float64 __x) throw(); extern long long __llroundf64(_Float64 __x) throw(); 
# 326
extern _Float64 fdimf64(_Float64 __x, _Float64 __y) throw(); extern _Float64 __fdimf64(_Float64 __x, _Float64 __y) throw(); 
# 329
extern _Float64 fmaxf64(_Float64 __x, _Float64 __y) throw() __attribute((const)); extern _Float64 __fmaxf64(_Float64 __x, _Float64 __y) throw() __attribute((const)); 
# 332
extern _Float64 fminf64(_Float64 __x, _Float64 __y) throw() __attribute((const)); extern _Float64 __fminf64(_Float64 __x, _Float64 __y) throw() __attribute((const)); 
# 335
extern _Float64 fmaf64(_Float64 __x, _Float64 __y, _Float64 __z) throw(); extern _Float64 __fmaf64(_Float64 __x, _Float64 __y, _Float64 __z) throw(); 
# 340
extern _Float64 roundevenf64(_Float64 __x) throw() __attribute((const)); extern _Float64 __roundevenf64(_Float64 __x) throw() __attribute((const)); 
# 344
extern __intmax_t fromfpf64(_Float64 __x, int __round, unsigned __width) throw(); extern __intmax_t __fromfpf64(_Float64 __x, int __round, unsigned __width) throw(); 
# 349
extern __uintmax_t ufromfpf64(_Float64 __x, int __round, unsigned __width) throw(); extern __uintmax_t __ufromfpf64(_Float64 __x, int __round, unsigned __width) throw(); 
# 355
extern __intmax_t fromfpxf64(_Float64 __x, int __round, unsigned __width) throw(); extern __intmax_t __fromfpxf64(_Float64 __x, int __round, unsigned __width) throw(); 
# 361
extern __uintmax_t ufromfpxf64(_Float64 __x, int __round, unsigned __width) throw(); extern __uintmax_t __ufromfpxf64(_Float64 __x, int __round, unsigned __width) throw(); 
# 365
extern _Float64 fmaxmagf64(_Float64 __x, _Float64 __y) throw() __attribute((const)); extern _Float64 __fmaxmagf64(_Float64 __x, _Float64 __y) throw() __attribute((const)); 
# 368
extern _Float64 fminmagf64(_Float64 __x, _Float64 __y) throw() __attribute((const)); extern _Float64 __fminmagf64(_Float64 __x, _Float64 __y) throw() __attribute((const)); 
# 371
extern int canonicalizef64(_Float64 * __cx, const _Float64 * __x) throw(); 
# 376
extern int totalorderf64(const _Float64 * __x, const _Float64 * __y) throw()
# 378
 __attribute((__pure__)); 
# 381
extern int totalordermagf64(const _Float64 * __x, const _Float64 * __y) throw()
# 383
 __attribute((__pure__)); 
# 386
extern _Float64 getpayloadf64(const _Float64 * __x) throw(); extern _Float64 __getpayloadf64(const _Float64 * __x) throw(); 
# 389
extern int setpayloadf64(_Float64 * __x, _Float64 __payload) throw(); 
# 392
extern int setpayloadsigf64(_Float64 * __x, _Float64 __payload) throw(); 
# 53 "/usr/include/x86_64-linux-gnu/bits/mathcalls.h" 3
extern _Float32x acosf32x(_Float32x __x) throw(); extern _Float32x __acosf32x(_Float32x __x) throw(); 
# 55
extern _Float32x asinf32x(_Float32x __x) throw(); extern _Float32x __asinf32x(_Float32x __x) throw(); 
# 57
extern _Float32x atanf32x(_Float32x __x) throw(); extern _Float32x __atanf32x(_Float32x __x) throw(); 
# 59
extern _Float32x atan2f32x(_Float32x __y, _Float32x __x) throw(); extern _Float32x __atan2f32x(_Float32x __y, _Float32x __x) throw(); 
# 62
extern _Float32x cosf32x(_Float32x __x) throw(); extern _Float32x __cosf32x(_Float32x __x) throw(); 
# 64
extern _Float32x sinf32x(_Float32x __x) throw(); extern _Float32x __sinf32x(_Float32x __x) throw(); 
# 66
extern _Float32x tanf32x(_Float32x __x) throw(); extern _Float32x __tanf32x(_Float32x __x) throw(); 
# 71
extern _Float32x coshf32x(_Float32x __x) throw(); extern _Float32x __coshf32x(_Float32x __x) throw(); 
# 73
extern _Float32x sinhf32x(_Float32x __x) throw(); extern _Float32x __sinhf32x(_Float32x __x) throw(); 
# 75
extern _Float32x tanhf32x(_Float32x __x) throw(); extern _Float32x __tanhf32x(_Float32x __x) throw(); 
# 79
extern void sincosf32x(_Float32x __x, _Float32x * __sinx, _Float32x * __cosx) throw(); extern void __sincosf32x(_Float32x __x, _Float32x * __sinx, _Float32x * __cosx) throw(); 
# 85
extern _Float32x acoshf32x(_Float32x __x) throw(); extern _Float32x __acoshf32x(_Float32x __x) throw(); 
# 87
extern _Float32x asinhf32x(_Float32x __x) throw(); extern _Float32x __asinhf32x(_Float32x __x) throw(); 
# 89
extern _Float32x atanhf32x(_Float32x __x) throw(); extern _Float32x __atanhf32x(_Float32x __x) throw(); 
# 95
extern _Float32x expf32x(_Float32x __x) throw(); extern _Float32x __expf32x(_Float32x __x) throw(); 
# 98
extern _Float32x frexpf32x(_Float32x __x, int * __exponent) throw(); extern _Float32x __frexpf32x(_Float32x __x, int * __exponent) throw(); 
# 101
extern _Float32x ldexpf32x(_Float32x __x, int __exponent) throw(); extern _Float32x __ldexpf32x(_Float32x __x, int __exponent) throw(); 
# 104
extern _Float32x logf32x(_Float32x __x) throw(); extern _Float32x __logf32x(_Float32x __x) throw(); 
# 107
extern _Float32x log10f32x(_Float32x __x) throw(); extern _Float32x __log10f32x(_Float32x __x) throw(); 
# 110
extern _Float32x modff32x(_Float32x __x, _Float32x * __iptr) throw(); extern _Float32x __modff32x(_Float32x __x, _Float32x * __iptr) throw() __attribute((__nonnull__(2))); 
# 114
extern _Float32x exp10f32x(_Float32x __x) throw(); extern _Float32x __exp10f32x(_Float32x __x) throw(); 
# 119
extern _Float32x expm1f32x(_Float32x __x) throw(); extern _Float32x __expm1f32x(_Float32x __x) throw(); 
# 122
extern _Float32x log1pf32x(_Float32x __x) throw(); extern _Float32x __log1pf32x(_Float32x __x) throw(); 
# 125
extern _Float32x logbf32x(_Float32x __x) throw(); extern _Float32x __logbf32x(_Float32x __x) throw(); 
# 130
extern _Float32x exp2f32x(_Float32x __x) throw(); extern _Float32x __exp2f32x(_Float32x __x) throw(); 
# 133
extern _Float32x log2f32x(_Float32x __x) throw(); extern _Float32x __log2f32x(_Float32x __x) throw(); 
# 140
extern _Float32x powf32x(_Float32x __x, _Float32x __y) throw(); extern _Float32x __powf32x(_Float32x __x, _Float32x __y) throw(); 
# 143
extern _Float32x sqrtf32x(_Float32x __x) throw(); extern _Float32x __sqrtf32x(_Float32x __x) throw(); 
# 147
extern _Float32x hypotf32x(_Float32x __x, _Float32x __y) throw(); extern _Float32x __hypotf32x(_Float32x __x, _Float32x __y) throw(); 
# 152
extern _Float32x cbrtf32x(_Float32x __x) throw(); extern _Float32x __cbrtf32x(_Float32x __x) throw(); 
# 159
extern _Float32x ceilf32x(_Float32x __x) throw() __attribute((const)); extern _Float32x __ceilf32x(_Float32x __x) throw() __attribute((const)); 
# 162
extern _Float32x fabsf32x(_Float32x __x) throw() __attribute((const)); extern _Float32x __fabsf32x(_Float32x __x) throw() __attribute((const)); 
# 165
extern _Float32x floorf32x(_Float32x __x) throw() __attribute((const)); extern _Float32x __floorf32x(_Float32x __x) throw() __attribute((const)); 
# 168
extern _Float32x fmodf32x(_Float32x __x, _Float32x __y) throw(); extern _Float32x __fmodf32x(_Float32x __x, _Float32x __y) throw(); 
# 196 "/usr/include/x86_64-linux-gnu/bits/mathcalls.h" 3
extern _Float32x copysignf32x(_Float32x __x, _Float32x __y) throw() __attribute((const)); extern _Float32x __copysignf32x(_Float32x __x, _Float32x __y) throw() __attribute((const)); 
# 201
extern _Float32x nanf32x(const char * __tagb) throw(); extern _Float32x __nanf32x(const char * __tagb) throw(); 
# 217 "/usr/include/x86_64-linux-gnu/bits/mathcalls.h" 3
extern _Float32x j0f32x(_Float32x) throw(); extern _Float32x __j0f32x(_Float32x) throw(); 
# 218
extern _Float32x j1f32x(_Float32x) throw(); extern _Float32x __j1f32x(_Float32x) throw(); 
# 219
extern _Float32x jnf32x(int, _Float32x) throw(); extern _Float32x __jnf32x(int, _Float32x) throw(); 
# 220
extern _Float32x y0f32x(_Float32x) throw(); extern _Float32x __y0f32x(_Float32x) throw(); 
# 221
extern _Float32x y1f32x(_Float32x) throw(); extern _Float32x __y1f32x(_Float32x) throw(); 
# 222
extern _Float32x ynf32x(int, _Float32x) throw(); extern _Float32x __ynf32x(int, _Float32x) throw(); 
# 228
extern _Float32x erff32x(_Float32x) throw(); extern _Float32x __erff32x(_Float32x) throw(); 
# 229
extern _Float32x erfcf32x(_Float32x) throw(); extern _Float32x __erfcf32x(_Float32x) throw(); 
# 230
extern _Float32x lgammaf32x(_Float32x) throw(); extern _Float32x __lgammaf32x(_Float32x) throw(); 
# 235
extern _Float32x tgammaf32x(_Float32x) throw(); extern _Float32x __tgammaf32x(_Float32x) throw(); 
# 249 "/usr/include/x86_64-linux-gnu/bits/mathcalls.h" 3
extern _Float32x lgammaf32x_r(_Float32x, int * __signgamp) throw(); extern _Float32x __lgammaf32x_r(_Float32x, int * __signgamp) throw(); 
# 256
extern _Float32x rintf32x(_Float32x __x) throw(); extern _Float32x __rintf32x(_Float32x __x) throw(); 
# 259
extern _Float32x nextafterf32x(_Float32x __x, _Float32x __y) throw(); extern _Float32x __nextafterf32x(_Float32x __x, _Float32x __y) throw(); 
# 266
extern _Float32x nextdownf32x(_Float32x __x) throw(); extern _Float32x __nextdownf32x(_Float32x __x) throw(); 
# 268
extern _Float32x nextupf32x(_Float32x __x) throw(); extern _Float32x __nextupf32x(_Float32x __x) throw(); 
# 272
extern _Float32x remainderf32x(_Float32x __x, _Float32x __y) throw(); extern _Float32x __remainderf32x(_Float32x __x, _Float32x __y) throw(); 
# 276
extern _Float32x scalbnf32x(_Float32x __x, int __n) throw(); extern _Float32x __scalbnf32x(_Float32x __x, int __n) throw(); 
# 280
extern int ilogbf32x(_Float32x __x) throw(); extern int __ilogbf32x(_Float32x __x) throw(); 
# 285
extern long llogbf32x(_Float32x __x) throw(); extern long __llogbf32x(_Float32x __x) throw(); 
# 290
extern _Float32x scalblnf32x(_Float32x __x, long __n) throw(); extern _Float32x __scalblnf32x(_Float32x __x, long __n) throw(); 
# 294
extern _Float32x nearbyintf32x(_Float32x __x) throw(); extern _Float32x __nearbyintf32x(_Float32x __x) throw(); 
# 298
extern _Float32x roundf32x(_Float32x __x) throw() __attribute((const)); extern _Float32x __roundf32x(_Float32x __x) throw() __attribute((const)); 
# 302
extern _Float32x truncf32x(_Float32x __x) throw() __attribute((const)); extern _Float32x __truncf32x(_Float32x __x) throw() __attribute((const)); 
# 307
extern _Float32x remquof32x(_Float32x __x, _Float32x __y, int * __quo) throw(); extern _Float32x __remquof32x(_Float32x __x, _Float32x __y, int * __quo) throw(); 
# 314
extern long lrintf32x(_Float32x __x) throw(); extern long __lrintf32x(_Float32x __x) throw(); 
# 316
__extension__ extern long long llrintf32x(_Float32x __x) throw(); extern long long __llrintf32x(_Float32x __x) throw(); 
# 320
extern long lroundf32x(_Float32x __x) throw(); extern long __lroundf32x(_Float32x __x) throw(); 
# 322
__extension__ extern long long llroundf32x(_Float32x __x) throw(); extern long long __llroundf32x(_Float32x __x) throw(); 
# 326
extern _Float32x fdimf32x(_Float32x __x, _Float32x __y) throw(); extern _Float32x __fdimf32x(_Float32x __x, _Float32x __y) throw(); 
# 329
extern _Float32x fmaxf32x(_Float32x __x, _Float32x __y) throw() __attribute((const)); extern _Float32x __fmaxf32x(_Float32x __x, _Float32x __y) throw() __attribute((const)); 
# 332
extern _Float32x fminf32x(_Float32x __x, _Float32x __y) throw() __attribute((const)); extern _Float32x __fminf32x(_Float32x __x, _Float32x __y) throw() __attribute((const)); 
# 335
extern _Float32x fmaf32x(_Float32x __x, _Float32x __y, _Float32x __z) throw(); extern _Float32x __fmaf32x(_Float32x __x, _Float32x __y, _Float32x __z) throw(); 
# 340
extern _Float32x roundevenf32x(_Float32x __x) throw() __attribute((const)); extern _Float32x __roundevenf32x(_Float32x __x) throw() __attribute((const)); 
# 344
extern __intmax_t fromfpf32x(_Float32x __x, int __round, unsigned __width) throw(); extern __intmax_t __fromfpf32x(_Float32x __x, int __round, unsigned __width) throw(); 
# 349
extern __uintmax_t ufromfpf32x(_Float32x __x, int __round, unsigned __width) throw(); extern __uintmax_t __ufromfpf32x(_Float32x __x, int __round, unsigned __width) throw(); 
# 355
extern __intmax_t fromfpxf32x(_Float32x __x, int __round, unsigned __width) throw(); extern __intmax_t __fromfpxf32x(_Float32x __x, int __round, unsigned __width) throw(); 
# 361
extern __uintmax_t ufromfpxf32x(_Float32x __x, int __round, unsigned __width) throw(); extern __uintmax_t __ufromfpxf32x(_Float32x __x, int __round, unsigned __width) throw(); 
# 365
extern _Float32x fmaxmagf32x(_Float32x __x, _Float32x __y) throw() __attribute((const)); extern _Float32x __fmaxmagf32x(_Float32x __x, _Float32x __y) throw() __attribute((const)); 
# 368
extern _Float32x fminmagf32x(_Float32x __x, _Float32x __y) throw() __attribute((const)); extern _Float32x __fminmagf32x(_Float32x __x, _Float32x __y) throw() __attribute((const)); 
# 371
extern int canonicalizef32x(_Float32x * __cx, const _Float32x * __x) throw(); 
# 376
extern int totalorderf32x(const _Float32x * __x, const _Float32x * __y) throw()
# 378
 __attribute((__pure__)); 
# 381
extern int totalordermagf32x(const _Float32x * __x, const _Float32x * __y) throw()
# 383
 __attribute((__pure__)); 
# 386
extern _Float32x getpayloadf32x(const _Float32x * __x) throw(); extern _Float32x __getpayloadf32x(const _Float32x * __x) throw(); 
# 389
extern int setpayloadf32x(_Float32x * __x, _Float32x __payload) throw(); 
# 392
extern int setpayloadsigf32x(_Float32x * __x, _Float32x __payload) throw(); 
# 53 "/usr/include/x86_64-linux-gnu/bits/mathcalls.h" 3
extern _Float64x acosf64x(_Float64x __x) throw(); extern _Float64x __acosf64x(_Float64x __x) throw(); 
# 55
extern _Float64x asinf64x(_Float64x __x) throw(); extern _Float64x __asinf64x(_Float64x __x) throw(); 
# 57
extern _Float64x atanf64x(_Float64x __x) throw(); extern _Float64x __atanf64x(_Float64x __x) throw(); 
# 59
extern _Float64x atan2f64x(_Float64x __y, _Float64x __x) throw(); extern _Float64x __atan2f64x(_Float64x __y, _Float64x __x) throw(); 
# 62
extern _Float64x cosf64x(_Float64x __x) throw(); extern _Float64x __cosf64x(_Float64x __x) throw(); 
# 64
extern _Float64x sinf64x(_Float64x __x) throw(); extern _Float64x __sinf64x(_Float64x __x) throw(); 
# 66
extern _Float64x tanf64x(_Float64x __x) throw(); extern _Float64x __tanf64x(_Float64x __x) throw(); 
# 71
extern _Float64x coshf64x(_Float64x __x) throw(); extern _Float64x __coshf64x(_Float64x __x) throw(); 
# 73
extern _Float64x sinhf64x(_Float64x __x) throw(); extern _Float64x __sinhf64x(_Float64x __x) throw(); 
# 75
extern _Float64x tanhf64x(_Float64x __x) throw(); extern _Float64x __tanhf64x(_Float64x __x) throw(); 
# 79
extern void sincosf64x(_Float64x __x, _Float64x * __sinx, _Float64x * __cosx) throw(); extern void __sincosf64x(_Float64x __x, _Float64x * __sinx, _Float64x * __cosx) throw(); 
# 85
extern _Float64x acoshf64x(_Float64x __x) throw(); extern _Float64x __acoshf64x(_Float64x __x) throw(); 
# 87
extern _Float64x asinhf64x(_Float64x __x) throw(); extern _Float64x __asinhf64x(_Float64x __x) throw(); 
# 89
extern _Float64x atanhf64x(_Float64x __x) throw(); extern _Float64x __atanhf64x(_Float64x __x) throw(); 
# 95
extern _Float64x expf64x(_Float64x __x) throw(); extern _Float64x __expf64x(_Float64x __x) throw(); 
# 98
extern _Float64x frexpf64x(_Float64x __x, int * __exponent) throw(); extern _Float64x __frexpf64x(_Float64x __x, int * __exponent) throw(); 
# 101
extern _Float64x ldexpf64x(_Float64x __x, int __exponent) throw(); extern _Float64x __ldexpf64x(_Float64x __x, int __exponent) throw(); 
# 104
extern _Float64x logf64x(_Float64x __x) throw(); extern _Float64x __logf64x(_Float64x __x) throw(); 
# 107
extern _Float64x log10f64x(_Float64x __x) throw(); extern _Float64x __log10f64x(_Float64x __x) throw(); 
# 110
extern _Float64x modff64x(_Float64x __x, _Float64x * __iptr) throw(); extern _Float64x __modff64x(_Float64x __x, _Float64x * __iptr) throw() __attribute((__nonnull__(2))); 
# 114
extern _Float64x exp10f64x(_Float64x __x) throw(); extern _Float64x __exp10f64x(_Float64x __x) throw(); 
# 119
extern _Float64x expm1f64x(_Float64x __x) throw(); extern _Float64x __expm1f64x(_Float64x __x) throw(); 
# 122
extern _Float64x log1pf64x(_Float64x __x) throw(); extern _Float64x __log1pf64x(_Float64x __x) throw(); 
# 125
extern _Float64x logbf64x(_Float64x __x) throw(); extern _Float64x __logbf64x(_Float64x __x) throw(); 
# 130
extern _Float64x exp2f64x(_Float64x __x) throw(); extern _Float64x __exp2f64x(_Float64x __x) throw(); 
# 133
extern _Float64x log2f64x(_Float64x __x) throw(); extern _Float64x __log2f64x(_Float64x __x) throw(); 
# 140
extern _Float64x powf64x(_Float64x __x, _Float64x __y) throw(); extern _Float64x __powf64x(_Float64x __x, _Float64x __y) throw(); 
# 143
extern _Float64x sqrtf64x(_Float64x __x) throw(); extern _Float64x __sqrtf64x(_Float64x __x) throw(); 
# 147
extern _Float64x hypotf64x(_Float64x __x, _Float64x __y) throw(); extern _Float64x __hypotf64x(_Float64x __x, _Float64x __y) throw(); 
# 152
extern _Float64x cbrtf64x(_Float64x __x) throw(); extern _Float64x __cbrtf64x(_Float64x __x) throw(); 
# 159
extern _Float64x ceilf64x(_Float64x __x) throw() __attribute((const)); extern _Float64x __ceilf64x(_Float64x __x) throw() __attribute((const)); 
# 162
extern _Float64x fabsf64x(_Float64x __x) throw() __attribute((const)); extern _Float64x __fabsf64x(_Float64x __x) throw() __attribute((const)); 
# 165
extern _Float64x floorf64x(_Float64x __x) throw() __attribute((const)); extern _Float64x __floorf64x(_Float64x __x) throw() __attribute((const)); 
# 168
extern _Float64x fmodf64x(_Float64x __x, _Float64x __y) throw(); extern _Float64x __fmodf64x(_Float64x __x, _Float64x __y) throw(); 
# 196 "/usr/include/x86_64-linux-gnu/bits/mathcalls.h" 3
extern _Float64x copysignf64x(_Float64x __x, _Float64x __y) throw() __attribute((const)); extern _Float64x __copysignf64x(_Float64x __x, _Float64x __y) throw() __attribute((const)); 
# 201
extern _Float64x nanf64x(const char * __tagb) throw(); extern _Float64x __nanf64x(const char * __tagb) throw(); 
# 217 "/usr/include/x86_64-linux-gnu/bits/mathcalls.h" 3
extern _Float64x j0f64x(_Float64x) throw(); extern _Float64x __j0f64x(_Float64x) throw(); 
# 218
extern _Float64x j1f64x(_Float64x) throw(); extern _Float64x __j1f64x(_Float64x) throw(); 
# 219
extern _Float64x jnf64x(int, _Float64x) throw(); extern _Float64x __jnf64x(int, _Float64x) throw(); 
# 220
extern _Float64x y0f64x(_Float64x) throw(); extern _Float64x __y0f64x(_Float64x) throw(); 
# 221
extern _Float64x y1f64x(_Float64x) throw(); extern _Float64x __y1f64x(_Float64x) throw(); 
# 222
extern _Float64x ynf64x(int, _Float64x) throw(); extern _Float64x __ynf64x(int, _Float64x) throw(); 
# 228
extern _Float64x erff64x(_Float64x) throw(); extern _Float64x __erff64x(_Float64x) throw(); 
# 229
extern _Float64x erfcf64x(_Float64x) throw(); extern _Float64x __erfcf64x(_Float64x) throw(); 
# 230
extern _Float64x lgammaf64x(_Float64x) throw(); extern _Float64x __lgammaf64x(_Float64x) throw(); 
# 235
extern _Float64x tgammaf64x(_Float64x) throw(); extern _Float64x __tgammaf64x(_Float64x) throw(); 
# 249 "/usr/include/x86_64-linux-gnu/bits/mathcalls.h" 3
extern _Float64x lgammaf64x_r(_Float64x, int * __signgamp) throw(); extern _Float64x __lgammaf64x_r(_Float64x, int * __signgamp) throw(); 
# 256
extern _Float64x rintf64x(_Float64x __x) throw(); extern _Float64x __rintf64x(_Float64x __x) throw(); 
# 259
extern _Float64x nextafterf64x(_Float64x __x, _Float64x __y) throw(); extern _Float64x __nextafterf64x(_Float64x __x, _Float64x __y) throw(); 
# 266
extern _Float64x nextdownf64x(_Float64x __x) throw(); extern _Float64x __nextdownf64x(_Float64x __x) throw(); 
# 268
extern _Float64x nextupf64x(_Float64x __x) throw(); extern _Float64x __nextupf64x(_Float64x __x) throw(); 
# 272
extern _Float64x remainderf64x(_Float64x __x, _Float64x __y) throw(); extern _Float64x __remainderf64x(_Float64x __x, _Float64x __y) throw(); 
# 276
extern _Float64x scalbnf64x(_Float64x __x, int __n) throw(); extern _Float64x __scalbnf64x(_Float64x __x, int __n) throw(); 
# 280
extern int ilogbf64x(_Float64x __x) throw(); extern int __ilogbf64x(_Float64x __x) throw(); 
# 285
extern long llogbf64x(_Float64x __x) throw(); extern long __llogbf64x(_Float64x __x) throw(); 
# 290
extern _Float64x scalblnf64x(_Float64x __x, long __n) throw(); extern _Float64x __scalblnf64x(_Float64x __x, long __n) throw(); 
# 294
extern _Float64x nearbyintf64x(_Float64x __x) throw(); extern _Float64x __nearbyintf64x(_Float64x __x) throw(); 
# 298
extern _Float64x roundf64x(_Float64x __x) throw() __attribute((const)); extern _Float64x __roundf64x(_Float64x __x) throw() __attribute((const)); 
# 302
extern _Float64x truncf64x(_Float64x __x) throw() __attribute((const)); extern _Float64x __truncf64x(_Float64x __x) throw() __attribute((const)); 
# 307
extern _Float64x remquof64x(_Float64x __x, _Float64x __y, int * __quo) throw(); extern _Float64x __remquof64x(_Float64x __x, _Float64x __y, int * __quo) throw(); 
# 314
extern long lrintf64x(_Float64x __x) throw(); extern long __lrintf64x(_Float64x __x) throw(); 
# 316
__extension__ extern long long llrintf64x(_Float64x __x) throw(); extern long long __llrintf64x(_Float64x __x) throw(); 
# 320
extern long lroundf64x(_Float64x __x) throw(); extern long __lroundf64x(_Float64x __x) throw(); 
# 322
__extension__ extern long long llroundf64x(_Float64x __x) throw(); extern long long __llroundf64x(_Float64x __x) throw(); 
# 326
extern _Float64x fdimf64x(_Float64x __x, _Float64x __y) throw(); extern _Float64x __fdimf64x(_Float64x __x, _Float64x __y) throw(); 
# 329
extern _Float64x fmaxf64x(_Float64x __x, _Float64x __y) throw() __attribute((const)); extern _Float64x __fmaxf64x(_Float64x __x, _Float64x __y) throw() __attribute((const)); 
# 332
extern _Float64x fminf64x(_Float64x __x, _Float64x __y) throw() __attribute((const)); extern _Float64x __fminf64x(_Float64x __x, _Float64x __y) throw() __attribute((const)); 
# 335
extern _Float64x fmaf64x(_Float64x __x, _Float64x __y, _Float64x __z) throw(); extern _Float64x __fmaf64x(_Float64x __x, _Float64x __y, _Float64x __z) throw(); 
# 340
extern _Float64x roundevenf64x(_Float64x __x) throw() __attribute((const)); extern _Float64x __roundevenf64x(_Float64x __x) throw() __attribute((const)); 
# 344
extern __intmax_t fromfpf64x(_Float64x __x, int __round, unsigned __width) throw(); extern __intmax_t __fromfpf64x(_Float64x __x, int __round, unsigned __width) throw(); 
# 349
extern __uintmax_t ufromfpf64x(_Float64x __x, int __round, unsigned __width) throw(); extern __uintmax_t __ufromfpf64x(_Float64x __x, int __round, unsigned __width) throw(); 
# 355
extern __intmax_t fromfpxf64x(_Float64x __x, int __round, unsigned __width) throw(); extern __intmax_t __fromfpxf64x(_Float64x __x, int __round, unsigned __width) throw(); 
# 361
extern __uintmax_t ufromfpxf64x(_Float64x __x, int __round, unsigned __width) throw(); extern __uintmax_t __ufromfpxf64x(_Float64x __x, int __round, unsigned __width) throw(); 
# 365
extern _Float64x fmaxmagf64x(_Float64x __x, _Float64x __y) throw() __attribute((const)); extern _Float64x __fmaxmagf64x(_Float64x __x, _Float64x __y) throw() __attribute((const)); 
# 368
extern _Float64x fminmagf64x(_Float64x __x, _Float64x __y) throw() __attribute((const)); extern _Float64x __fminmagf64x(_Float64x __x, _Float64x __y) throw() __attribute((const)); 
# 371
extern int canonicalizef64x(_Float64x * __cx, const _Float64x * __x) throw(); 
# 376
extern int totalorderf64x(const _Float64x * __x, const _Float64x * __y) throw()
# 378
 __attribute((__pure__)); 
# 381
extern int totalordermagf64x(const _Float64x * __x, const _Float64x * __y) throw()
# 383
 __attribute((__pure__)); 
# 386
extern _Float64x getpayloadf64x(const _Float64x * __x) throw(); extern _Float64x __getpayloadf64x(const _Float64x * __x) throw(); 
# 389
extern int setpayloadf64x(_Float64x * __x, _Float64x __payload) throw(); 
# 392
extern int setpayloadsigf64x(_Float64x * __x, _Float64x __payload) throw(); 
# 24 "/usr/include/x86_64-linux-gnu/bits/mathcalls-narrow.h" 3
extern float fadd(double __x, double __y) throw(); 
# 27
extern float fdiv(double __x, double __y) throw(); 
# 30
extern float fmul(double __x, double __y) throw(); 
# 33
extern float fsub(double __x, double __y) throw(); 
# 24 "/usr/include/x86_64-linux-gnu/bits/mathcalls-narrow.h" 3
extern float faddl(long double __x, long double __y) throw(); 
# 27
extern float fdivl(long double __x, long double __y) throw(); 
# 30
extern float fmull(long double __x, long double __y) throw(); 
# 33
extern float fsubl(long double __x, long double __y) throw(); 
# 24 "/usr/include/x86_64-linux-gnu/bits/mathcalls-narrow.h" 3
extern double daddl(long double __x, long double __y) throw(); 
# 27
extern double ddivl(long double __x, long double __y) throw(); 
# 30
extern double dmull(long double __x, long double __y) throw(); 
# 33
extern double dsubl(long double __x, long double __y) throw(); 
# 24 "/usr/include/x86_64-linux-gnu/bits/mathcalls-narrow.h" 3
extern _Float32 f32addf32x(_Float32x __x, _Float32x __y) throw(); 
# 27
extern _Float32 f32divf32x(_Float32x __x, _Float32x __y) throw(); 
# 30
extern _Float32 f32mulf32x(_Float32x __x, _Float32x __y) throw(); 
# 33
extern _Float32 f32subf32x(_Float32x __x, _Float32x __y) throw(); 
# 24 "/usr/include/x86_64-linux-gnu/bits/mathcalls-narrow.h" 3
extern _Float32 f32addf64(_Float64 __x, _Float64 __y) throw(); 
# 27
extern _Float32 f32divf64(_Float64 __x, _Float64 __y) throw(); 
# 30
extern _Float32 f32mulf64(_Float64 __x, _Float64 __y) throw(); 
# 33
extern _Float32 f32subf64(_Float64 __x, _Float64 __y) throw(); 
# 24 "/usr/include/x86_64-linux-gnu/bits/mathcalls-narrow.h" 3
extern _Float32 f32addf64x(_Float64x __x, _Float64x __y) throw(); 
# 27
extern _Float32 f32divf64x(_Float64x __x, _Float64x __y) throw(); 
# 30
extern _Float32 f32mulf64x(_Float64x __x, _Float64x __y) throw(); 
# 33
extern _Float32 f32subf64x(_Float64x __x, _Float64x __y) throw(); 
# 24 "/usr/include/x86_64-linux-gnu/bits/mathcalls-narrow.h" 3
extern _Float32x f32xaddf64(_Float64 __x, _Float64 __y) throw(); 
# 27
extern _Float32x f32xdivf64(_Float64 __x, _Float64 __y) throw(); 
# 30
extern _Float32x f32xmulf64(_Float64 __x, _Float64 __y) throw(); 
# 33
extern _Float32x f32xsubf64(_Float64 __x, _Float64 __y) throw(); 
# 24 "/usr/include/x86_64-linux-gnu/bits/mathcalls-narrow.h" 3
extern _Float32x f32xaddf64x(_Float64x __x, _Float64x __y) throw(); 
# 27
extern _Float32x f32xdivf64x(_Float64x __x, _Float64x __y) throw(); 
# 30
extern _Float32x f32xmulf64x(_Float64x __x, _Float64x __y) throw(); 
# 33
extern _Float32x f32xsubf64x(_Float64x __x, _Float64x __y) throw(); 
# 24 "/usr/include/x86_64-linux-gnu/bits/mathcalls-narrow.h" 3
extern _Float64 f64addf64x(_Float64x __x, _Float64x __y) throw(); 
# 27
extern _Float64 f64divf64x(_Float64x __x, _Float64x __y) throw(); 
# 30
extern _Float64 f64mulf64x(_Float64x __x, _Float64x __y) throw(); 
# 33
extern _Float64 f64subf64x(_Float64x __x, _Float64x __y) throw(); 
# 773 "/usr/include/math.h" 3
extern int signgam; 
# 854 "/usr/include/math.h" 3
enum { 
# 855
FP_NAN, 
# 858
FP_INFINITE, 
# 861
FP_ZERO, 
# 864
FP_SUBNORMAL, 
# 867
FP_NORMAL
# 870
}; 
# 23 "/usr/include/x86_64-linux-gnu/bits/iscanonical.h" 3
extern int __iscanonicall(long double __x) throw()
# 24
 __attribute((const)); 
# 46 "/usr/include/x86_64-linux-gnu/bits/iscanonical.h" 3
extern "C++" {
# 47
inline int iscanonical(float __val) { return (((void)((__typeof__(__val))__val)), 1); } 
# 48
inline int iscanonical(double __val) { return (((void)((__typeof__(__val))__val)), 1); } 
# 49
inline int iscanonical(long double __val) { return __iscanonicall(__val); } 
# 53
}
# 985 "/usr/include/math.h" 3
extern "C++" {
# 986
inline int issignaling(float __val) { return __issignalingf(__val); } 
# 987
inline int issignaling(double __val) { return __issignaling(__val); } 
# 989
inline int issignaling(long double __val) 
# 990
{ 
# 994
return __issignalingl(__val); 
# 996
} 
# 1002
}
# 1016 "/usr/include/math.h" 3
extern "C++" {
# 1047 "/usr/include/math.h" 3
template< class __T> inline bool 
# 1048
iszero(__T __val) 
# 1049
{ 
# 1050
return __val == 0; 
# 1051
} 
# 1053
}
# 1278 "/usr/include/math.h" 3
extern "C++" {
# 1279
template< class > struct __iseqsig_type; 
# 1281
template<> struct __iseqsig_type< float>  { 
# 1283
static int __call(float __x, float __y) throw() 
# 1284
{ 
# 1285
return __iseqsigf(__x, __y); 
# 1286
} 
# 1287
}; 
# 1289
template<> struct __iseqsig_type< double>  { 
# 1291
static int __call(double __x, double __y) throw() 
# 1292
{ 
# 1293
return __iseqsig(__x, __y); 
# 1294
} 
# 1295
}; 
# 1297
template<> struct __iseqsig_type< long double>  { 
# 1299
static int __call(long double __x, long double __y) throw() 
# 1300
{ 
# 1302
return __iseqsigl(__x, __y); 
# 1306
} 
# 1307
}; 
# 1321 "/usr/include/math.h" 3
template< class _T1, class _T2> inline int 
# 1323
iseqsig(_T1 __x, _T2 __y) throw() 
# 1324
{ 
# 1326
typedef __decltype(((__x + __y) + (0.0F))) _T3; 
# 1330
return __iseqsig_type< __decltype(((__x + __y) + (0.0F)))> ::__call(__x, __y); 
# 1331
} 
# 1333
}
# 1338
}
# 77 "/usr/include/c++/10/cmath" 3
extern "C++" {
# 79
namespace std __attribute((__visibility__("default"))) { 
# 83
using ::acos;
# 87
constexpr float acos(float __x) 
# 88
{ return __builtin_acosf(__x); } 
# 91
constexpr long double acos(long double __x) 
# 92
{ return __builtin_acosl(__x); } 
# 95
template< class _Tp> constexpr typename __gnu_cxx::__enable_if< __is_integer< _Tp> ::__value, double> ::__type 
# 99
acos(_Tp __x) 
# 100
{ return __builtin_acos(__x); } 
# 102
using ::asin;
# 106
constexpr float asin(float __x) 
# 107
{ return __builtin_asinf(__x); } 
# 110
constexpr long double asin(long double __x) 
# 111
{ return __builtin_asinl(__x); } 
# 114
template< class _Tp> constexpr typename __gnu_cxx::__enable_if< __is_integer< _Tp> ::__value, double> ::__type 
# 118
asin(_Tp __x) 
# 119
{ return __builtin_asin(__x); } 
# 121
using ::atan;
# 125
constexpr float atan(float __x) 
# 126
{ return __builtin_atanf(__x); } 
# 129
constexpr long double atan(long double __x) 
# 130
{ return __builtin_atanl(__x); } 
# 133
template< class _Tp> constexpr typename __gnu_cxx::__enable_if< __is_integer< _Tp> ::__value, double> ::__type 
# 137
atan(_Tp __x) 
# 138
{ return __builtin_atan(__x); } 
# 140
using ::atan2;
# 144
constexpr float atan2(float __y, float __x) 
# 145
{ return __builtin_atan2f(__y, __x); } 
# 148
constexpr long double atan2(long double __y, long double __x) 
# 149
{ return __builtin_atan2l(__y, __x); } 
# 152
template< class _Tp, class _Up> constexpr typename __gnu_cxx::__promote_2< _Tp, _Up> ::__type 
# 155
atan2(_Tp __y, _Up __x) 
# 156
{ 
# 157
typedef typename __gnu_cxx::__promote_2< _Tp, _Up> ::__type __type; 
# 158
return atan2((__type)__y, (__type)__x); 
# 159
} 
# 161
using ::ceil;
# 165
constexpr float ceil(float __x) 
# 166
{ return __builtin_ceilf(__x); } 
# 169
constexpr long double ceil(long double __x) 
# 170
{ return __builtin_ceill(__x); } 
# 173
template< class _Tp> constexpr typename __gnu_cxx::__enable_if< __is_integer< _Tp> ::__value, double> ::__type 
# 177
ceil(_Tp __x) 
# 178
{ return __builtin_ceil(__x); } 
# 180
using ::cos;
# 184
constexpr float cos(float __x) 
# 185
{ return __builtin_cosf(__x); } 
# 188
constexpr long double cos(long double __x) 
# 189
{ return __builtin_cosl(__x); } 
# 192
template< class _Tp> constexpr typename __gnu_cxx::__enable_if< __is_integer< _Tp> ::__value, double> ::__type 
# 196
cos(_Tp __x) 
# 197
{ return __builtin_cos(__x); } 
# 199
using ::cosh;
# 203
constexpr float cosh(float __x) 
# 204
{ return __builtin_coshf(__x); } 
# 207
constexpr long double cosh(long double __x) 
# 208
{ return __builtin_coshl(__x); } 
# 211
template< class _Tp> constexpr typename __gnu_cxx::__enable_if< __is_integer< _Tp> ::__value, double> ::__type 
# 215
cosh(_Tp __x) 
# 216
{ return __builtin_cosh(__x); } 
# 218
using ::exp;
# 222
constexpr float exp(float __x) 
# 223
{ return __builtin_expf(__x); } 
# 226
constexpr long double exp(long double __x) 
# 227
{ return __builtin_expl(__x); } 
# 230
template< class _Tp> constexpr typename __gnu_cxx::__enable_if< __is_integer< _Tp> ::__value, double> ::__type 
# 234
exp(_Tp __x) 
# 235
{ return __builtin_exp(__x); } 
# 237
using ::fabs;
# 241
constexpr float fabs(float __x) 
# 242
{ return __builtin_fabsf(__x); } 
# 245
constexpr long double fabs(long double __x) 
# 246
{ return __builtin_fabsl(__x); } 
# 249
template< class _Tp> constexpr typename __gnu_cxx::__enable_if< __is_integer< _Tp> ::__value, double> ::__type 
# 253
fabs(_Tp __x) 
# 254
{ return __builtin_fabs(__x); } 
# 256
using ::floor;
# 260
constexpr float floor(float __x) 
# 261
{ return __builtin_floorf(__x); } 
# 264
constexpr long double floor(long double __x) 
# 265
{ return __builtin_floorl(__x); } 
# 268
template< class _Tp> constexpr typename __gnu_cxx::__enable_if< __is_integer< _Tp> ::__value, double> ::__type 
# 272
floor(_Tp __x) 
# 273
{ return __builtin_floor(__x); } 
# 275
using ::fmod;
# 279
constexpr float fmod(float __x, float __y) 
# 280
{ return __builtin_fmodf(__x, __y); } 
# 283
constexpr long double fmod(long double __x, long double __y) 
# 284
{ return __builtin_fmodl(__x, __y); } 
# 287
template< class _Tp, class _Up> constexpr typename __gnu_cxx::__promote_2< _Tp, _Up> ::__type 
# 290
fmod(_Tp __x, _Up __y) 
# 291
{ 
# 292
typedef typename __gnu_cxx::__promote_2< _Tp, _Up> ::__type __type; 
# 293
return fmod((__type)__x, (__type)__y); 
# 294
} 
# 296
using ::frexp;
# 300
inline float frexp(float __x, int *__exp) 
# 301
{ return __builtin_frexpf(__x, __exp); } 
# 304
inline long double frexp(long double __x, int *__exp) 
# 305
{ return __builtin_frexpl(__x, __exp); } 
# 308
template< class _Tp> constexpr typename __gnu_cxx::__enable_if< __is_integer< _Tp> ::__value, double> ::__type 
# 312
frexp(_Tp __x, int *__exp) 
# 313
{ return __builtin_frexp(__x, __exp); } 
# 315
using ::ldexp;
# 319
constexpr float ldexp(float __x, int __exp) 
# 320
{ return __builtin_ldexpf(__x, __exp); } 
# 323
constexpr long double ldexp(long double __x, int __exp) 
# 324
{ return __builtin_ldexpl(__x, __exp); } 
# 327
template< class _Tp> constexpr typename __gnu_cxx::__enable_if< __is_integer< _Tp> ::__value, double> ::__type 
# 331
ldexp(_Tp __x, int __exp) 
# 332
{ return __builtin_ldexp(__x, __exp); } 
# 334
using ::log;
# 338
constexpr float log(float __x) 
# 339
{ return __builtin_logf(__x); } 
# 342
constexpr long double log(long double __x) 
# 343
{ return __builtin_logl(__x); } 
# 346
template< class _Tp> constexpr typename __gnu_cxx::__enable_if< __is_integer< _Tp> ::__value, double> ::__type 
# 350
log(_Tp __x) 
# 351
{ return __builtin_log(__x); } 
# 353
using ::log10;
# 357
constexpr float log10(float __x) 
# 358
{ return __builtin_log10f(__x); } 
# 361
constexpr long double log10(long double __x) 
# 362
{ return __builtin_log10l(__x); } 
# 365
template< class _Tp> constexpr typename __gnu_cxx::__enable_if< __is_integer< _Tp> ::__value, double> ::__type 
# 369
log10(_Tp __x) 
# 370
{ return __builtin_log10(__x); } 
# 372
using ::modf;
# 376
inline float modf(float __x, float *__iptr) 
# 377
{ return __builtin_modff(__x, __iptr); } 
# 380
inline long double modf(long double __x, long double *__iptr) 
# 381
{ return __builtin_modfl(__x, __iptr); } 
# 384
using ::pow;
# 388
constexpr float pow(float __x, float __y) 
# 389
{ return __builtin_powf(__x, __y); } 
# 392
constexpr long double pow(long double __x, long double __y) 
# 393
{ return __builtin_powl(__x, __y); } 
# 412 "/usr/include/c++/10/cmath" 3
template< class _Tp, class _Up> constexpr typename __gnu_cxx::__promote_2< _Tp, _Up> ::__type 
# 415
pow(_Tp __x, _Up __y) 
# 416
{ 
# 417
typedef typename __gnu_cxx::__promote_2< _Tp, _Up> ::__type __type; 
# 418
return pow((__type)__x, (__type)__y); 
# 419
} 
# 421
using ::sin;
# 425
constexpr float sin(float __x) 
# 426
{ return __builtin_sinf(__x); } 
# 429
constexpr long double sin(long double __x) 
# 430
{ return __builtin_sinl(__x); } 
# 433
template< class _Tp> constexpr typename __gnu_cxx::__enable_if< __is_integer< _Tp> ::__value, double> ::__type 
# 437
sin(_Tp __x) 
# 438
{ return __builtin_sin(__x); } 
# 440
using ::sinh;
# 444
constexpr float sinh(float __x) 
# 445
{ return __builtin_sinhf(__x); } 
# 448
constexpr long double sinh(long double __x) 
# 449
{ return __builtin_sinhl(__x); } 
# 452
template< class _Tp> constexpr typename __gnu_cxx::__enable_if< __is_integer< _Tp> ::__value, double> ::__type 
# 456
sinh(_Tp __x) 
# 457
{ return __builtin_sinh(__x); } 
# 459
using ::sqrt;
# 463
constexpr float sqrt(float __x) 
# 464
{ return __builtin_sqrtf(__x); } 
# 467
constexpr long double sqrt(long double __x) 
# 468
{ return __builtin_sqrtl(__x); } 
# 471
template< class _Tp> constexpr typename __gnu_cxx::__enable_if< __is_integer< _Tp> ::__value, double> ::__type 
# 475
sqrt(_Tp __x) 
# 476
{ return __builtin_sqrt(__x); } 
# 478
using ::tan;
# 482
constexpr float tan(float __x) 
# 483
{ return __builtin_tanf(__x); } 
# 486
constexpr long double tan(long double __x) 
# 487
{ return __builtin_tanl(__x); } 
# 490
template< class _Tp> constexpr typename __gnu_cxx::__enable_if< __is_integer< _Tp> ::__value, double> ::__type 
# 494
tan(_Tp __x) 
# 495
{ return __builtin_tan(__x); } 
# 497
using ::tanh;
# 501
constexpr float tanh(float __x) 
# 502
{ return __builtin_tanhf(__x); } 
# 505
constexpr long double tanh(long double __x) 
# 506
{ return __builtin_tanhl(__x); } 
# 509
template< class _Tp> constexpr typename __gnu_cxx::__enable_if< __is_integer< _Tp> ::__value, double> ::__type 
# 513
tanh(_Tp __x) 
# 514
{ return __builtin_tanh(__x); } 
# 537 "/usr/include/c++/10/cmath" 3
constexpr int fpclassify(float __x) 
# 538
{ return __builtin_fpclassify(0, 1, 4, 3, 2, __x); 
# 539
} 
# 542
constexpr int fpclassify(double __x) 
# 543
{ return __builtin_fpclassify(0, 1, 4, 3, 2, __x); 
# 544
} 
# 547
constexpr int fpclassify(long double __x) 
# 548
{ return __builtin_fpclassify(0, 1, 4, 3, 2, __x); 
# 549
} 
# 553
template< class _Tp> constexpr typename __gnu_cxx::__enable_if< __is_integer< _Tp> ::__value, int> ::__type 
# 556
fpclassify(_Tp __x) 
# 557
{ return (__x != 0) ? 4 : 2; } 
# 562
constexpr bool isfinite(float __x) 
# 563
{ return __builtin_isfinite(__x); } 
# 566
constexpr bool isfinite(double __x) 
# 567
{ return __builtin_isfinite(__x); } 
# 570
constexpr bool isfinite(long double __x) 
# 571
{ return __builtin_isfinite(__x); } 
# 575
template< class _Tp> constexpr typename __gnu_cxx::__enable_if< __is_integer< _Tp> ::__value, bool> ::__type 
# 578
isfinite(_Tp __x) 
# 579
{ return true; } 
# 584
constexpr bool isinf(float __x) 
# 585
{ return __builtin_isinf(__x); } 
# 592
constexpr bool isinf(double __x) 
# 593
{ return __builtin_isinf(__x); } 
# 597
constexpr bool isinf(long double __x) 
# 598
{ return __builtin_isinf(__x); } 
# 602
template< class _Tp> constexpr typename __gnu_cxx::__enable_if< __is_integer< _Tp> ::__value, bool> ::__type 
# 605
isinf(_Tp __x) 
# 606
{ return false; } 
# 611
constexpr bool isnan(float __x) 
# 612
{ return __builtin_isnan(__x); } 
# 619
constexpr bool isnan(double __x) 
# 620
{ return __builtin_isnan(__x); } 
# 624
constexpr bool isnan(long double __x) 
# 625
{ return __builtin_isnan(__x); } 
# 629
template< class _Tp> constexpr typename __gnu_cxx::__enable_if< __is_integer< _Tp> ::__value, bool> ::__type 
# 632
isnan(_Tp __x) 
# 633
{ return false; } 
# 638
constexpr bool isnormal(float __x) 
# 639
{ return __builtin_isnormal(__x); } 
# 642
constexpr bool isnormal(double __x) 
# 643
{ return __builtin_isnormal(__x); } 
# 646
constexpr bool isnormal(long double __x) 
# 647
{ return __builtin_isnormal(__x); } 
# 651
template< class _Tp> constexpr typename __gnu_cxx::__enable_if< __is_integer< _Tp> ::__value, bool> ::__type 
# 654
isnormal(_Tp __x) 
# 655
{ return (__x != 0) ? true : false; } 
# 661
constexpr bool signbit(float __x) 
# 662
{ return __builtin_signbit(__x); } 
# 665
constexpr bool signbit(double __x) 
# 666
{ return __builtin_signbit(__x); } 
# 669
constexpr bool signbit(long double __x) 
# 670
{ return __builtin_signbit(__x); } 
# 674
template< class _Tp> constexpr typename __gnu_cxx::__enable_if< __is_integer< _Tp> ::__value, bool> ::__type 
# 677
signbit(_Tp __x) 
# 678
{ return (__x < 0) ? true : false; } 
# 683
constexpr bool isgreater(float __x, float __y) 
# 684
{ return __builtin_isgreater(__x, __y); } 
# 687
constexpr bool isgreater(double __x, double __y) 
# 688
{ return __builtin_isgreater(__x, __y); } 
# 691
constexpr bool isgreater(long double __x, long double __y) 
# 692
{ return __builtin_isgreater(__x, __y); } 
# 696
template< class _Tp, class _Up> constexpr typename __gnu_cxx::__enable_if< __is_arithmetic< _Tp> ::__value && __is_arithmetic< _Up> ::__value, bool> ::__type 
# 700
isgreater(_Tp __x, _Up __y) 
# 701
{ 
# 702
typedef typename __gnu_cxx::__promote_2< _Tp, _Up> ::__type __type; 
# 703
return __builtin_isgreater((__type)__x, (__type)__y); 
# 704
} 
# 709
constexpr bool isgreaterequal(float __x, float __y) 
# 710
{ return __builtin_isgreaterequal(__x, __y); } 
# 713
constexpr bool isgreaterequal(double __x, double __y) 
# 714
{ return __builtin_isgreaterequal(__x, __y); } 
# 717
constexpr bool isgreaterequal(long double __x, long double __y) 
# 718
{ return __builtin_isgreaterequal(__x, __y); } 
# 722
template< class _Tp, class _Up> constexpr typename __gnu_cxx::__enable_if< __is_arithmetic< _Tp> ::__value && __is_arithmetic< _Up> ::__value, bool> ::__type 
# 726
isgreaterequal(_Tp __x, _Up __y) 
# 727
{ 
# 728
typedef typename __gnu_cxx::__promote_2< _Tp, _Up> ::__type __type; 
# 729
return __builtin_isgreaterequal((__type)__x, (__type)__y); 
# 730
} 
# 735
constexpr bool isless(float __x, float __y) 
# 736
{ return __builtin_isless(__x, __y); } 
# 739
constexpr bool isless(double __x, double __y) 
# 740
{ return __builtin_isless(__x, __y); } 
# 743
constexpr bool isless(long double __x, long double __y) 
# 744
{ return __builtin_isless(__x, __y); } 
# 748
template< class _Tp, class _Up> constexpr typename __gnu_cxx::__enable_if< __is_arithmetic< _Tp> ::__value && __is_arithmetic< _Up> ::__value, bool> ::__type 
# 752
isless(_Tp __x, _Up __y) 
# 753
{ 
# 754
typedef typename __gnu_cxx::__promote_2< _Tp, _Up> ::__type __type; 
# 755
return __builtin_isless((__type)__x, (__type)__y); 
# 756
} 
# 761
constexpr bool islessequal(float __x, float __y) 
# 762
{ return __builtin_islessequal(__x, __y); } 
# 765
constexpr bool islessequal(double __x, double __y) 
# 766
{ return __builtin_islessequal(__x, __y); } 
# 769
constexpr bool islessequal(long double __x, long double __y) 
# 770
{ return __builtin_islessequal(__x, __y); } 
# 774
template< class _Tp, class _Up> constexpr typename __gnu_cxx::__enable_if< __is_arithmetic< _Tp> ::__value && __is_arithmetic< _Up> ::__value, bool> ::__type 
# 778
islessequal(_Tp __x, _Up __y) 
# 779
{ 
# 780
typedef typename __gnu_cxx::__promote_2< _Tp, _Up> ::__type __type; 
# 781
return __builtin_islessequal((__type)__x, (__type)__y); 
# 782
} 
# 787
constexpr bool islessgreater(float __x, float __y) 
# 788
{ return __builtin_islessgreater(__x, __y); } 
# 791
constexpr bool islessgreater(double __x, double __y) 
# 792
{ return __builtin_islessgreater(__x, __y); } 
# 795
constexpr bool islessgreater(long double __x, long double __y) 
# 796
{ return __builtin_islessgreater(__x, __y); } 
# 800
template< class _Tp, class _Up> constexpr typename __gnu_cxx::__enable_if< __is_arithmetic< _Tp> ::__value && __is_arithmetic< _Up> ::__value, bool> ::__type 
# 804
islessgreater(_Tp __x, _Up __y) 
# 805
{ 
# 806
typedef typename __gnu_cxx::__promote_2< _Tp, _Up> ::__type __type; 
# 807
return __builtin_islessgreater((__type)__x, (__type)__y); 
# 808
} 
# 813
constexpr bool isunordered(float __x, float __y) 
# 814
{ return __builtin_isunordered(__x, __y); } 
# 817
constexpr bool isunordered(double __x, double __y) 
# 818
{ return __builtin_isunordered(__x, __y); } 
# 821
constexpr bool isunordered(long double __x, long double __y) 
# 822
{ return __builtin_isunordered(__x, __y); } 
# 826
template< class _Tp, class _Up> constexpr typename __gnu_cxx::__enable_if< __is_arithmetic< _Tp> ::__value && __is_arithmetic< _Up> ::__value, bool> ::__type 
# 830
isunordered(_Tp __x, _Up __y) 
# 831
{ 
# 832
typedef typename __gnu_cxx::__promote_2< _Tp, _Up> ::__type __type; 
# 833
return __builtin_isunordered((__type)__x, (__type)__y); 
# 834
} 
# 1065 "/usr/include/c++/10/cmath" 3
using ::double_t;
# 1066
using ::float_t;
# 1069
using ::acosh;
# 1070
using ::acoshf;
# 1071
using ::acoshl;
# 1073
using ::asinh;
# 1074
using ::asinhf;
# 1075
using ::asinhl;
# 1077
using ::atanh;
# 1078
using ::atanhf;
# 1079
using ::atanhl;
# 1081
using ::cbrt;
# 1082
using ::cbrtf;
# 1083
using ::cbrtl;
# 1085
using ::copysign;
# 1086
using ::copysignf;
# 1087
using ::copysignl;
# 1089
using ::erf;
# 1090
using ::erff;
# 1091
using ::erfl;
# 1093
using ::erfc;
# 1094
using ::erfcf;
# 1095
using ::erfcl;
# 1097
using ::exp2;
# 1098
using ::exp2f;
# 1099
using ::exp2l;
# 1101
using ::expm1;
# 1102
using ::expm1f;
# 1103
using ::expm1l;
# 1105
using ::fdim;
# 1106
using ::fdimf;
# 1107
using ::fdiml;
# 1109
using ::fma;
# 1110
using ::fmaf;
# 1111
using ::fmal;
# 1113
using ::fmax;
# 1114
using ::fmaxf;
# 1115
using ::fmaxl;
# 1117
using ::fmin;
# 1118
using ::fminf;
# 1119
using ::fminl;
# 1121
using ::hypot;
# 1122
using ::hypotf;
# 1123
using ::hypotl;
# 1125
using ::ilogb;
# 1126
using ::ilogbf;
# 1127
using ::ilogbl;
# 1129
using ::lgamma;
# 1130
using ::lgammaf;
# 1131
using ::lgammal;
# 1134
using ::llrint;
# 1135
using ::llrintf;
# 1136
using ::llrintl;
# 1138
using ::llround;
# 1139
using ::llroundf;
# 1140
using ::llroundl;
# 1143
using ::log1p;
# 1144
using ::log1pf;
# 1145
using ::log1pl;
# 1147
using ::log2;
# 1148
using ::log2f;
# 1149
using ::log2l;
# 1151
using ::logb;
# 1152
using ::logbf;
# 1153
using ::logbl;
# 1155
using ::lrint;
# 1156
using ::lrintf;
# 1157
using ::lrintl;
# 1159
using ::lround;
# 1160
using ::lroundf;
# 1161
using ::lroundl;
# 1163
using ::nan;
# 1164
using ::nanf;
# 1165
using ::nanl;
# 1167
using ::nearbyint;
# 1168
using ::nearbyintf;
# 1169
using ::nearbyintl;
# 1171
using ::nextafter;
# 1172
using ::nextafterf;
# 1173
using ::nextafterl;
# 1175
using ::nexttoward;
# 1176
using ::nexttowardf;
# 1177
using ::nexttowardl;
# 1179
using ::remainder;
# 1180
using ::remainderf;
# 1181
using ::remainderl;
# 1183
using ::remquo;
# 1184
using ::remquof;
# 1185
using ::remquol;
# 1187
using ::rint;
# 1188
using ::rintf;
# 1189
using ::rintl;
# 1191
using ::round;
# 1192
using ::roundf;
# 1193
using ::roundl;
# 1195
using ::scalbln;
# 1196
using ::scalblnf;
# 1197
using ::scalblnl;
# 1199
using ::scalbn;
# 1200
using ::scalbnf;
# 1201
using ::scalbnl;
# 1203
using ::tgamma;
# 1204
using ::tgammaf;
# 1205
using ::tgammal;
# 1207
using ::trunc;
# 1208
using ::truncf;
# 1209
using ::truncl;
# 1214
constexpr float acosh(float __x) 
# 1215
{ return __builtin_acoshf(__x); } 
# 1218
constexpr long double acosh(long double __x) 
# 1219
{ return __builtin_acoshl(__x); } 
# 1223
template< class _Tp> constexpr typename __gnu_cxx::__enable_if< __is_integer< _Tp> ::__value, double> ::__type 
# 1226
acosh(_Tp __x) 
# 1227
{ return __builtin_acosh(__x); } 
# 1232
constexpr float asinh(float __x) 
# 1233
{ return __builtin_asinhf(__x); } 
# 1236
constexpr long double asinh(long double __x) 
# 1237
{ return __builtin_asinhl(__x); } 
# 1241
template< class _Tp> constexpr typename __gnu_cxx::__enable_if< __is_integer< _Tp> ::__value, double> ::__type 
# 1244
asinh(_Tp __x) 
# 1245
{ return __builtin_asinh(__x); } 
# 1250
constexpr float atanh(float __x) 
# 1251
{ return __builtin_atanhf(__x); } 
# 1254
constexpr long double atanh(long double __x) 
# 1255
{ return __builtin_atanhl(__x); } 
# 1259
template< class _Tp> constexpr typename __gnu_cxx::__enable_if< __is_integer< _Tp> ::__value, double> ::__type 
# 1262
atanh(_Tp __x) 
# 1263
{ return __builtin_atanh(__x); } 
# 1268
constexpr float cbrt(float __x) 
# 1269
{ return __builtin_cbrtf(__x); } 
# 1272
constexpr long double cbrt(long double __x) 
# 1273
{ return __builtin_cbrtl(__x); } 
# 1277
template< class _Tp> constexpr typename __gnu_cxx::__enable_if< __is_integer< _Tp> ::__value, double> ::__type 
# 1280
cbrt(_Tp __x) 
# 1281
{ return __builtin_cbrt(__x); } 
# 1286
constexpr float copysign(float __x, float __y) 
# 1287
{ return __builtin_copysignf(__x, __y); } 
# 1290
constexpr long double copysign(long double __x, long double __y) 
# 1291
{ return __builtin_copysignl(__x, __y); } 
# 1295
template< class _Tp, class _Up> constexpr typename __gnu_cxx::__promote_2< _Tp, _Up> ::__type 
# 1297
copysign(_Tp __x, _Up __y) 
# 1298
{ 
# 1299
typedef typename __gnu_cxx::__promote_2< _Tp, _Up> ::__type __type; 
# 1300
return copysign((__type)__x, (__type)__y); 
# 1301
} 
# 1306
constexpr float erf(float __x) 
# 1307
{ return __builtin_erff(__x); } 
# 1310
constexpr long double erf(long double __x) 
# 1311
{ return __builtin_erfl(__x); } 
# 1315
template< class _Tp> constexpr typename __gnu_cxx::__enable_if< __is_integer< _Tp> ::__value, double> ::__type 
# 1318
erf(_Tp __x) 
# 1319
{ return __builtin_erf(__x); } 
# 1324
constexpr float erfc(float __x) 
# 1325
{ return __builtin_erfcf(__x); } 
# 1328
constexpr long double erfc(long double __x) 
# 1329
{ return __builtin_erfcl(__x); } 
# 1333
template< class _Tp> constexpr typename __gnu_cxx::__enable_if< __is_integer< _Tp> ::__value, double> ::__type 
# 1336
erfc(_Tp __x) 
# 1337
{ return __builtin_erfc(__x); } 
# 1342
constexpr float exp2(float __x) 
# 1343
{ return __builtin_exp2f(__x); } 
# 1346
constexpr long double exp2(long double __x) 
# 1347
{ return __builtin_exp2l(__x); } 
# 1351
template< class _Tp> constexpr typename __gnu_cxx::__enable_if< __is_integer< _Tp> ::__value, double> ::__type 
# 1354
exp2(_Tp __x) 
# 1355
{ return __builtin_exp2(__x); } 
# 1360
constexpr float expm1(float __x) 
# 1361
{ return __builtin_expm1f(__x); } 
# 1364
constexpr long double expm1(long double __x) 
# 1365
{ return __builtin_expm1l(__x); } 
# 1369
template< class _Tp> constexpr typename __gnu_cxx::__enable_if< __is_integer< _Tp> ::__value, double> ::__type 
# 1372
expm1(_Tp __x) 
# 1373
{ return __builtin_expm1(__x); } 
# 1378
constexpr float fdim(float __x, float __y) 
# 1379
{ return __builtin_fdimf(__x, __y); } 
# 1382
constexpr long double fdim(long double __x, long double __y) 
# 1383
{ return __builtin_fdiml(__x, __y); } 
# 1387
template< class _Tp, class _Up> constexpr typename __gnu_cxx::__promote_2< _Tp, _Up> ::__type 
# 1389
fdim(_Tp __x, _Up __y) 
# 1390
{ 
# 1391
typedef typename __gnu_cxx::__promote_2< _Tp, _Up> ::__type __type; 
# 1392
return fdim((__type)__x, (__type)__y); 
# 1393
} 
# 1398
constexpr float fma(float __x, float __y, float __z) 
# 1399
{ return __builtin_fmaf(__x, __y, __z); } 
# 1402
constexpr long double fma(long double __x, long double __y, long double __z) 
# 1403
{ return __builtin_fmal(__x, __y, __z); } 
# 1407
template< class _Tp, class _Up, class _Vp> constexpr typename __gnu_cxx::__promote_3< _Tp, _Up, _Vp> ::__type 
# 1409
fma(_Tp __x, _Up __y, _Vp __z) 
# 1410
{ 
# 1411
typedef typename __gnu_cxx::__promote_3< _Tp, _Up, _Vp> ::__type __type; 
# 1412
return fma((__type)__x, (__type)__y, (__type)__z); 
# 1413
} 
# 1418
constexpr float fmax(float __x, float __y) 
# 1419
{ return __builtin_fmaxf(__x, __y); } 
# 1422
constexpr long double fmax(long double __x, long double __y) 
# 1423
{ return __builtin_fmaxl(__x, __y); } 
# 1427
template< class _Tp, class _Up> constexpr typename __gnu_cxx::__promote_2< _Tp, _Up> ::__type 
# 1429
fmax(_Tp __x, _Up __y) 
# 1430
{ 
# 1431
typedef typename __gnu_cxx::__promote_2< _Tp, _Up> ::__type __type; 
# 1432
return fmax((__type)__x, (__type)__y); 
# 1433
} 
# 1438
constexpr float fmin(float __x, float __y) 
# 1439
{ return __builtin_fminf(__x, __y); } 
# 1442
constexpr long double fmin(long double __x, long double __y) 
# 1443
{ return __builtin_fminl(__x, __y); } 
# 1447
template< class _Tp, class _Up> constexpr typename __gnu_cxx::__promote_2< _Tp, _Up> ::__type 
# 1449
fmin(_Tp __x, _Up __y) 
# 1450
{ 
# 1451
typedef typename __gnu_cxx::__promote_2< _Tp, _Up> ::__type __type; 
# 1452
return fmin((__type)__x, (__type)__y); 
# 1453
} 
# 1458
constexpr float hypot(float __x, float __y) 
# 1459
{ return __builtin_hypotf(__x, __y); } 
# 1462
constexpr long double hypot(long double __x, long double __y) 
# 1463
{ return __builtin_hypotl(__x, __y); } 
# 1467
template< class _Tp, class _Up> constexpr typename __gnu_cxx::__promote_2< _Tp, _Up> ::__type 
# 1469
hypot(_Tp __x, _Up __y) 
# 1470
{ 
# 1471
typedef typename __gnu_cxx::__promote_2< _Tp, _Up> ::__type __type; 
# 1472
return hypot((__type)__x, (__type)__y); 
# 1473
} 
# 1478
constexpr int ilogb(float __x) 
# 1479
{ return __builtin_ilogbf(__x); } 
# 1482
constexpr int ilogb(long double __x) 
# 1483
{ return __builtin_ilogbl(__x); } 
# 1487
template< class _Tp> constexpr typename __gnu_cxx::__enable_if< __is_integer< _Tp> ::__value, int> ::__type 
# 1491
ilogb(_Tp __x) 
# 1492
{ return __builtin_ilogb(__x); } 
# 1497
constexpr float lgamma(float __x) 
# 1498
{ return __builtin_lgammaf(__x); } 
# 1501
constexpr long double lgamma(long double __x) 
# 1502
{ return __builtin_lgammal(__x); } 
# 1506
template< class _Tp> constexpr typename __gnu_cxx::__enable_if< __is_integer< _Tp> ::__value, double> ::__type 
# 1509
lgamma(_Tp __x) 
# 1510
{ return __builtin_lgamma(__x); } 
# 1515
constexpr long long llrint(float __x) 
# 1516
{ return __builtin_llrintf(__x); } 
# 1519
constexpr long long llrint(long double __x) 
# 1520
{ return __builtin_llrintl(__x); } 
# 1524
template< class _Tp> constexpr typename __gnu_cxx::__enable_if< __is_integer< _Tp> ::__value, long long> ::__type 
# 1527
llrint(_Tp __x) 
# 1528
{ return __builtin_llrint(__x); } 
# 1533
constexpr long long llround(float __x) 
# 1534
{ return __builtin_llroundf(__x); } 
# 1537
constexpr long long llround(long double __x) 
# 1538
{ return __builtin_llroundl(__x); } 
# 1542
template< class _Tp> constexpr typename __gnu_cxx::__enable_if< __is_integer< _Tp> ::__value, long long> ::__type 
# 1545
llround(_Tp __x) 
# 1546
{ return __builtin_llround(__x); } 
# 1551
constexpr float log1p(float __x) 
# 1552
{ return __builtin_log1pf(__x); } 
# 1555
constexpr long double log1p(long double __x) 
# 1556
{ return __builtin_log1pl(__x); } 
# 1560
template< class _Tp> constexpr typename __gnu_cxx::__enable_if< __is_integer< _Tp> ::__value, double> ::__type 
# 1563
log1p(_Tp __x) 
# 1564
{ return __builtin_log1p(__x); } 
# 1570
constexpr float log2(float __x) 
# 1571
{ return __builtin_log2f(__x); } 
# 1574
constexpr long double log2(long double __x) 
# 1575
{ return __builtin_log2l(__x); } 
# 1579
template< class _Tp> constexpr typename __gnu_cxx::__enable_if< __is_integer< _Tp> ::__value, double> ::__type 
# 1582
log2(_Tp __x) 
# 1583
{ return __builtin_log2(__x); } 
# 1588
constexpr float logb(float __x) 
# 1589
{ return __builtin_logbf(__x); } 
# 1592
constexpr long double logb(long double __x) 
# 1593
{ return __builtin_logbl(__x); } 
# 1597
template< class _Tp> constexpr typename __gnu_cxx::__enable_if< __is_integer< _Tp> ::__value, double> ::__type 
# 1600
logb(_Tp __x) 
# 1601
{ return __builtin_logb(__x); } 
# 1606
constexpr long lrint(float __x) 
# 1607
{ return __builtin_lrintf(__x); } 
# 1610
constexpr long lrint(long double __x) 
# 1611
{ return __builtin_lrintl(__x); } 
# 1615
template< class _Tp> constexpr typename __gnu_cxx::__enable_if< __is_integer< _Tp> ::__value, long> ::__type 
# 1618
lrint(_Tp __x) 
# 1619
{ return __builtin_lrint(__x); } 
# 1624
constexpr long lround(float __x) 
# 1625
{ return __builtin_lroundf(__x); } 
# 1628
constexpr long lround(long double __x) 
# 1629
{ return __builtin_lroundl(__x); } 
# 1633
template< class _Tp> constexpr typename __gnu_cxx::__enable_if< __is_integer< _Tp> ::__value, long> ::__type 
# 1636
lround(_Tp __x) 
# 1637
{ return __builtin_lround(__x); } 
# 1642
constexpr float nearbyint(float __x) 
# 1643
{ return __builtin_nearbyintf(__x); } 
# 1646
constexpr long double nearbyint(long double __x) 
# 1647
{ return __builtin_nearbyintl(__x); } 
# 1651
template< class _Tp> constexpr typename __gnu_cxx::__enable_if< __is_integer< _Tp> ::__value, double> ::__type 
# 1654
nearbyint(_Tp __x) 
# 1655
{ return __builtin_nearbyint(__x); } 
# 1660
constexpr float nextafter(float __x, float __y) 
# 1661
{ return __builtin_nextafterf(__x, __y); } 
# 1664
constexpr long double nextafter(long double __x, long double __y) 
# 1665
{ return __builtin_nextafterl(__x, __y); } 
# 1669
template< class _Tp, class _Up> constexpr typename __gnu_cxx::__promote_2< _Tp, _Up> ::__type 
# 1671
nextafter(_Tp __x, _Up __y) 
# 1672
{ 
# 1673
typedef typename __gnu_cxx::__promote_2< _Tp, _Up> ::__type __type; 
# 1674
return nextafter((__type)__x, (__type)__y); 
# 1675
} 
# 1680
constexpr float nexttoward(float __x, long double __y) 
# 1681
{ return __builtin_nexttowardf(__x, __y); } 
# 1684
constexpr long double nexttoward(long double __x, long double __y) 
# 1685
{ return __builtin_nexttowardl(__x, __y); } 
# 1689
template< class _Tp> constexpr typename __gnu_cxx::__enable_if< __is_integer< _Tp> ::__value, double> ::__type 
# 1692
nexttoward(_Tp __x, long double __y) 
# 1693
{ return __builtin_nexttoward(__x, __y); } 
# 1698
constexpr float remainder(float __x, float __y) 
# 1699
{ return __builtin_remainderf(__x, __y); } 
# 1702
constexpr long double remainder(long double __x, long double __y) 
# 1703
{ return __builtin_remainderl(__x, __y); } 
# 1707
template< class _Tp, class _Up> constexpr typename __gnu_cxx::__promote_2< _Tp, _Up> ::__type 
# 1709
remainder(_Tp __x, _Up __y) 
# 1710
{ 
# 1711
typedef typename __gnu_cxx::__promote_2< _Tp, _Up> ::__type __type; 
# 1712
return remainder((__type)__x, (__type)__y); 
# 1713
} 
# 1718
inline float remquo(float __x, float __y, int *__pquo) 
# 1719
{ return __builtin_remquof(__x, __y, __pquo); } 
# 1722
inline long double remquo(long double __x, long double __y, int *__pquo) 
# 1723
{ return __builtin_remquol(__x, __y, __pquo); } 
# 1727
template< class _Tp, class _Up> inline typename __gnu_cxx::__promote_2< _Tp, _Up> ::__type 
# 1729
remquo(_Tp __x, _Up __y, int *__pquo) 
# 1730
{ 
# 1731
typedef typename __gnu_cxx::__promote_2< _Tp, _Up> ::__type __type; 
# 1732
return remquo((__type)__x, (__type)__y, __pquo); 
# 1733
} 
# 1738
constexpr float rint(float __x) 
# 1739
{ return __builtin_rintf(__x); } 
# 1742
constexpr long double rint(long double __x) 
# 1743
{ return __builtin_rintl(__x); } 
# 1747
template< class _Tp> constexpr typename __gnu_cxx::__enable_if< __is_integer< _Tp> ::__value, double> ::__type 
# 1750
rint(_Tp __x) 
# 1751
{ return __builtin_rint(__x); } 
# 1756
constexpr float round(float __x) 
# 1757
{ return __builtin_roundf(__x); } 
# 1760
constexpr long double round(long double __x) 
# 1761
{ return __builtin_roundl(__x); } 
# 1765
template< class _Tp> constexpr typename __gnu_cxx::__enable_if< __is_integer< _Tp> ::__value, double> ::__type 
# 1768
round(_Tp __x) 
# 1769
{ return __builtin_round(__x); } 
# 1774
constexpr float scalbln(float __x, long __ex) 
# 1775
{ return __builtin_scalblnf(__x, __ex); } 
# 1778
constexpr long double scalbln(long double __x, long __ex) 
# 1779
{ return __builtin_scalblnl(__x, __ex); } 
# 1783
template< class _Tp> constexpr typename __gnu_cxx::__enable_if< __is_integer< _Tp> ::__value, double> ::__type 
# 1786
scalbln(_Tp __x, long __ex) 
# 1787
{ return __builtin_scalbln(__x, __ex); } 
# 1792
constexpr float scalbn(float __x, int __ex) 
# 1793
{ return __builtin_scalbnf(__x, __ex); } 
# 1796
constexpr long double scalbn(long double __x, int __ex) 
# 1797
{ return __builtin_scalbnl(__x, __ex); } 
# 1801
template< class _Tp> constexpr typename __gnu_cxx::__enable_if< __is_integer< _Tp> ::__value, double> ::__type 
# 1804
scalbn(_Tp __x, int __ex) 
# 1805
{ return __builtin_scalbn(__x, __ex); } 
# 1810
constexpr float tgamma(float __x) 
# 1811
{ return __builtin_tgammaf(__x); } 
# 1814
constexpr long double tgamma(long double __x) 
# 1815
{ return __builtin_tgammal(__x); } 
# 1819
template< class _Tp> constexpr typename __gnu_cxx::__enable_if< __is_integer< _Tp> ::__value, double> ::__type 
# 1822
tgamma(_Tp __x) 
# 1823
{ return __builtin_tgamma(__x); } 
# 1828
constexpr float trunc(float __x) 
# 1829
{ return __builtin_truncf(__x); } 
# 1832
constexpr long double trunc(long double __x) 
# 1833
{ return __builtin_truncl(__x); } 
# 1837
template< class _Tp> constexpr typename __gnu_cxx::__enable_if< __is_integer< _Tp> ::__value, double> ::__type 
# 1840
trunc(_Tp __x) 
# 1841
{ return __builtin_trunc(__x); } 
# 1932 "/usr/include/c++/10/cmath" 3
}
# 1938
}
# 38 "/usr/include/c++/10/math.h" 3
using std::abs;
# 39
using std::acos;
# 40
using std::asin;
# 41
using std::atan;
# 42
using std::atan2;
# 43
using std::cos;
# 44
using std::sin;
# 45
using std::tan;
# 46
using std::cosh;
# 47
using std::sinh;
# 48
using std::tanh;
# 49
using std::exp;
# 50
using std::frexp;
# 51
using std::ldexp;
# 52
using std::log;
# 53
using std::log10;
# 54
using std::modf;
# 55
using std::pow;
# 56
using std::sqrt;
# 57
using std::ceil;
# 58
using std::fabs;
# 59
using std::floor;
# 60
using std::fmod;
# 63
using std::fpclassify;
# 64
using std::isfinite;
# 65
using std::isinf;
# 66
using std::isnan;
# 67
using std::isnormal;
# 68
using std::signbit;
# 69
using std::isgreater;
# 70
using std::isgreaterequal;
# 71
using std::isless;
# 72
using std::islessequal;
# 73
using std::islessgreater;
# 74
using std::isunordered;
# 78
using std::acosh;
# 79
using std::asinh;
# 80
using std::atanh;
# 81
using std::cbrt;
# 82
using std::copysign;
# 83
using std::erf;
# 84
using std::erfc;
# 85
using std::exp2;
# 86
using std::expm1;
# 87
using std::fdim;
# 88
using std::fma;
# 89
using std::fmax;
# 90
using std::fmin;
# 91
using std::hypot;
# 92
using std::ilogb;
# 93
using std::lgamma;
# 94
using std::llrint;
# 95
using std::llround;
# 96
using std::log1p;
# 97
using std::log2;
# 98
using std::logb;
# 99
using std::lrint;
# 100
using std::lround;
# 101
using std::nearbyint;
# 102
using std::nextafter;
# 103
using std::nexttoward;
# 104
using std::remainder;
# 105
using std::remquo;
# 106
using std::rint;
# 107
using std::round;
# 108
using std::scalbln;
# 109
using std::scalbn;
# 110
using std::tgamma;
# 111
using std::trunc;
# 10623 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
namespace std { 
# 10624
constexpr bool signbit(float x); 
# 10625
constexpr bool signbit(double x); 
# 10626
constexpr bool signbit(long double x); 
# 10627
constexpr bool isfinite(float x); 
# 10628
constexpr bool isfinite(double x); 
# 10629
constexpr bool isfinite(long double x); 
# 10630
constexpr bool isnan(float x); 
# 10635
constexpr bool isnan(double x); 
# 10637
constexpr bool isnan(long double x); 
# 10638
constexpr bool isinf(float x); 
# 10643
constexpr bool isinf(double x); 
# 10645
constexpr bool isinf(long double x); 
# 10646
}
# 10802 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
namespace std { 
# 10804
template< class T> extern T __pow_helper(T, int); 
# 10805
template< class T> extern T __cmath_power(T, unsigned); 
# 10806
}
# 10808
using std::abs;
# 10809
using std::fabs;
# 10810
using std::ceil;
# 10811
using std::floor;
# 10812
using std::sqrt;
# 10814
using std::pow;
# 10816
using std::log;
# 10817
using std::log10;
# 10818
using std::fmod;
# 10819
using std::modf;
# 10820
using std::exp;
# 10821
using std::frexp;
# 10822
using std::ldexp;
# 10823
using std::asin;
# 10824
using std::sin;
# 10825
using std::sinh;
# 10826
using std::acos;
# 10827
using std::cos;
# 10828
using std::cosh;
# 10829
using std::atan;
# 10830
using std::atan2;
# 10831
using std::tan;
# 10832
using std::tanh;
# 11203 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
namespace std { 
# 11212 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern inline long long abs(long long); 
# 11222 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern inline long abs(long); 
# 11223
extern constexpr float abs(float); 
# 11224
extern constexpr double abs(double); 
# 11225
extern constexpr float fabs(float); 
# 11226
extern constexpr float ceil(float); 
# 11227
extern constexpr float floor(float); 
# 11228
extern constexpr float sqrt(float); 
# 11229
extern constexpr float pow(float, float); 
# 11234
template< class _Tp, class _Up> extern constexpr typename __gnu_cxx::__promote_2< _Tp, _Up> ::__type pow(_Tp, _Up); 
# 11244
extern constexpr float log(float); 
# 11245
extern constexpr float log10(float); 
# 11246
extern constexpr float fmod(float, float); 
# 11247
extern inline float modf(float, float *); 
# 11248
extern constexpr float exp(float); 
# 11249
extern inline float frexp(float, int *); 
# 11250
extern constexpr float ldexp(float, int); 
# 11251
extern constexpr float asin(float); 
# 11252
extern constexpr float sin(float); 
# 11253
extern constexpr float sinh(float); 
# 11254
extern constexpr float acos(float); 
# 11255
extern constexpr float cos(float); 
# 11256
extern constexpr float cosh(float); 
# 11257
extern constexpr float atan(float); 
# 11258
extern constexpr float atan2(float, float); 
# 11259
extern constexpr float tan(float); 
# 11260
extern constexpr float tanh(float); 
# 11343 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
}
# 11449 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
namespace std { 
# 11450
constexpr float logb(float a); 
# 11451
constexpr int ilogb(float a); 
# 11452
constexpr float scalbn(float a, int b); 
# 11453
constexpr float scalbln(float a, long b); 
# 11454
constexpr float exp2(float a); 
# 11455
constexpr float expm1(float a); 
# 11456
constexpr float log2(float a); 
# 11457
constexpr float log1p(float a); 
# 11458
constexpr float acosh(float a); 
# 11459
constexpr float asinh(float a); 
# 11460
constexpr float atanh(float a); 
# 11461
constexpr float hypot(float a, float b); 
# 11462
constexpr float cbrt(float a); 
# 11463
constexpr float erf(float a); 
# 11464
constexpr float erfc(float a); 
# 11465
constexpr float lgamma(float a); 
# 11466
constexpr float tgamma(float a); 
# 11467
constexpr float copysign(float a, float b); 
# 11468
constexpr float nextafter(float a, float b); 
# 11469
constexpr float remainder(float a, float b); 
# 11470
inline float remquo(float a, float b, int * quo); 
# 11471
constexpr float round(float a); 
# 11472
constexpr long lround(float a); 
# 11473
constexpr long long llround(float a); 
# 11474
constexpr float trunc(float a); 
# 11475
constexpr float rint(float a); 
# 11476
constexpr long lrint(float a); 
# 11477
constexpr long long llrint(float a); 
# 11478
constexpr float nearbyint(float a); 
# 11479
constexpr float fdim(float a, float b); 
# 11480
constexpr float fma(float a, float b, float c); 
# 11481
constexpr float fmax(float a, float b); 
# 11482
constexpr float fmin(float a, float b); 
# 11483
}
# 11588 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
static inline float exp10(const float a); 
# 11590
static inline float rsqrt(const float a); 
# 11592
static inline float rcbrt(const float a); 
# 11594
static inline float sinpi(const float a); 
# 11596
static inline float cospi(const float a); 
# 11598
static inline void sincospi(const float a, float *const sptr, float *const cptr); 
# 11600
static inline void sincos(const float a, float *const sptr, float *const cptr); 
# 11602
static inline float j0(const float a); 
# 11604
static inline float j1(const float a); 
# 11606
static inline float jn(const int n, const float a); 
# 11608
static inline float y0(const float a); 
# 11610
static inline float y1(const float a); 
# 11612
static inline float yn(const int n, const float a); 
# 11614
__attribute__((unused)) static inline float cyl_bessel_i0(const float a); 
# 11616
__attribute__((unused)) static inline float cyl_bessel_i1(const float a); 
# 11618
static inline float erfinv(const float a); 
# 11620
static inline float erfcinv(const float a); 
# 11622
static inline float normcdfinv(const float a); 
# 11624
static inline float normcdf(const float a); 
# 11626
static inline float erfcx(const float a); 
# 11628
static inline double copysign(const double a, const float b); 
# 11630
static inline double copysign(const float a, const double b); 
# 11638
static inline unsigned min(const unsigned a, const unsigned b); 
# 11646
static inline unsigned min(const int a, const unsigned b); 
# 11654
static inline unsigned min(const unsigned a, const int b); 
# 11662
static inline long min(const long a, const long b); 
# 11670
static inline unsigned long min(const unsigned long a, const unsigned long b); 
# 11678
static inline unsigned long min(const long a, const unsigned long b); 
# 11686
static inline unsigned long min(const unsigned long a, const long b); 
# 11694
static inline long long min(const long long a, const long long b); 
# 11702
static inline unsigned long long min(const unsigned long long a, const unsigned long long b); 
# 11710
static inline unsigned long long min(const long long a, const unsigned long long b); 
# 11718
static inline unsigned long long min(const unsigned long long a, const long long b); 
# 11729 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
static inline float min(const float a, const float b); 
# 11740 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
static inline double min(const double a, const double b); 
# 11750 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
static inline double min(const float a, const double b); 
# 11760 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
static inline double min(const double a, const float b); 
# 11768
static inline unsigned max(const unsigned a, const unsigned b); 
# 11776
static inline unsigned max(const int a, const unsigned b); 
# 11784
static inline unsigned max(const unsigned a, const int b); 
# 11792
static inline long max(const long a, const long b); 
# 11800
static inline unsigned long max(const unsigned long a, const unsigned long b); 
# 11808
static inline unsigned long max(const long a, const unsigned long b); 
# 11816
static inline unsigned long max(const unsigned long a, const long b); 
# 11824
static inline long long max(const long long a, const long long b); 
# 11832
static inline unsigned long long max(const unsigned long long a, const unsigned long long b); 
# 11840
static inline unsigned long long max(const long long a, const unsigned long long b); 
# 11848
static inline unsigned long long max(const unsigned long long a, const long long b); 
# 11859 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
static inline float max(const float a, const float b); 
# 11870 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
static inline double max(const double a, const double b); 
# 11880 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
static inline double max(const float a, const double b); 
# 11890 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
static inline double max(const double a, const float b); 
# 11901 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
extern "C" {
# 11902
__attribute__((unused)) inline void *__nv_aligned_device_malloc(size_t size, size_t align) 
# 11903
{int volatile ___ = 1;(void)size;(void)align;
# 11906
::exit(___);}
#if 0
# 11903
{ 
# 11904
__attribute__((unused)) void *__nv_aligned_device_malloc_impl(size_t, size_t); 
# 11905
return __nv_aligned_device_malloc_impl(size, align); 
# 11906
} 
#endif
# 11907 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.h"
}
# 758 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.hpp"
static inline float exp10(const float a) 
# 759
{ 
# 760
return exp10f(a); 
# 761
} 
# 763
static inline float rsqrt(const float a) 
# 764
{ 
# 765
return rsqrtf(a); 
# 766
} 
# 768
static inline float rcbrt(const float a) 
# 769
{ 
# 770
return rcbrtf(a); 
# 771
} 
# 773
static inline float sinpi(const float a) 
# 774
{ 
# 775
return sinpif(a); 
# 776
} 
# 778
static inline float cospi(const float a) 
# 779
{ 
# 780
return cospif(a); 
# 781
} 
# 783
static inline void sincospi(const float a, float *const sptr, float *const cptr) 
# 784
{ 
# 785
sincospif(a, sptr, cptr); 
# 786
} 
# 788
static inline void sincos(const float a, float *const sptr, float *const cptr) 
# 789
{ 
# 790
sincosf(a, sptr, cptr); 
# 791
} 
# 793
static inline float j0(const float a) 
# 794
{ 
# 795
return j0f(a); 
# 796
} 
# 798
static inline float j1(const float a) 
# 799
{ 
# 800
return j1f(a); 
# 801
} 
# 803
static inline float jn(const int n, const float a) 
# 804
{ 
# 805
return jnf(n, a); 
# 806
} 
# 808
static inline float y0(const float a) 
# 809
{ 
# 810
return y0f(a); 
# 811
} 
# 813
static inline float y1(const float a) 
# 814
{ 
# 815
return y1f(a); 
# 816
} 
# 818
static inline float yn(const int n, const float a) 
# 819
{ 
# 820
return ynf(n, a); 
# 821
} 
# 823
__attribute__((unused)) static inline float cyl_bessel_i0(const float a) 
# 824
{int volatile ___ = 1;(void)a;
# 826
::exit(___);}
#if 0
# 824
{ 
# 825
return cyl_bessel_i0f(a); 
# 826
} 
#endif
# 828 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.hpp"
__attribute__((unused)) static inline float cyl_bessel_i1(const float a) 
# 829
{int volatile ___ = 1;(void)a;
# 831
::exit(___);}
#if 0
# 829
{ 
# 830
return cyl_bessel_i1f(a); 
# 831
} 
#endif
# 833 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.hpp"
static inline float erfinv(const float a) 
# 834
{ 
# 835
return erfinvf(a); 
# 836
} 
# 838
static inline float erfcinv(const float a) 
# 839
{ 
# 840
return erfcinvf(a); 
# 841
} 
# 843
static inline float normcdfinv(const float a) 
# 844
{ 
# 845
return normcdfinvf(a); 
# 846
} 
# 848
static inline float normcdf(const float a) 
# 849
{ 
# 850
return normcdff(a); 
# 851
} 
# 853
static inline float erfcx(const float a) 
# 854
{ 
# 855
return erfcxf(a); 
# 856
} 
# 858
static inline double copysign(const double a, const float b) 
# 859
{ 
# 860
return copysign(a, static_cast< double>(b)); 
# 861
} 
# 863
static inline double copysign(const float a, const double b) 
# 864
{ 
# 865
return copysign(static_cast< double>(a), b); 
# 866
} 
# 868
static inline unsigned min(const unsigned a, const unsigned b) 
# 869
{ 
# 870
return umin(a, b); 
# 871
} 
# 873
static inline unsigned min(const int a, const unsigned b) 
# 874
{ 
# 875
return umin(static_cast< unsigned>(a), b); 
# 876
} 
# 878
static inline unsigned min(const unsigned a, const int b) 
# 879
{ 
# 880
return umin(a, static_cast< unsigned>(b)); 
# 881
} 
# 883
static inline long min(const long a, const long b) 
# 884
{ 
# 885
long retval; 
# 891
if (sizeof(long) == sizeof(int)) { 
# 895
retval = (static_cast< long>(min(static_cast< int>(a), static_cast< int>(b)))); 
# 896
} else { 
# 897
retval = (static_cast< long>(llmin(static_cast< long long>(a), static_cast< long long>(b)))); 
# 898
}  
# 899
return retval; 
# 900
} 
# 902
static inline unsigned long min(const unsigned long a, const unsigned long b) 
# 903
{ 
# 904
unsigned long retval; 
# 908
if (sizeof(unsigned long) == sizeof(unsigned)) { 
# 912
retval = (static_cast< unsigned long>(umin(static_cast< unsigned>(a), static_cast< unsigned>(b)))); 
# 913
} else { 
# 914
retval = (static_cast< unsigned long>(ullmin(static_cast< unsigned long long>(a), static_cast< unsigned long long>(b)))); 
# 915
}  
# 916
return retval; 
# 917
} 
# 919
static inline unsigned long min(const long a, const unsigned long b) 
# 920
{ 
# 921
unsigned long retval; 
# 925
if (sizeof(unsigned long) == sizeof(unsigned)) { 
# 929
retval = (static_cast< unsigned long>(umin(static_cast< unsigned>(a), static_cast< unsigned>(b)))); 
# 930
} else { 
# 931
retval = (static_cast< unsigned long>(ullmin(static_cast< unsigned long long>(a), static_cast< unsigned long long>(b)))); 
# 932
}  
# 933
return retval; 
# 934
} 
# 936
static inline unsigned long min(const unsigned long a, const long b) 
# 937
{ 
# 938
unsigned long retval; 
# 942
if (sizeof(unsigned long) == sizeof(unsigned)) { 
# 946
retval = (static_cast< unsigned long>(umin(static_cast< unsigned>(a), static_cast< unsigned>(b)))); 
# 947
} else { 
# 948
retval = (static_cast< unsigned long>(ullmin(static_cast< unsigned long long>(a), static_cast< unsigned long long>(b)))); 
# 949
}  
# 950
return retval; 
# 951
} 
# 953
static inline long long min(const long long a, const long long b) 
# 954
{ 
# 955
return llmin(a, b); 
# 956
} 
# 958
static inline unsigned long long min(const unsigned long long a, const unsigned long long b) 
# 959
{ 
# 960
return ullmin(a, b); 
# 961
} 
# 963
static inline unsigned long long min(const long long a, const unsigned long long b) 
# 964
{ 
# 965
return ullmin(static_cast< unsigned long long>(a), b); 
# 966
} 
# 968
static inline unsigned long long min(const unsigned long long a, const long long b) 
# 969
{ 
# 970
return ullmin(a, static_cast< unsigned long long>(b)); 
# 971
} 
# 973
static inline float min(const float a, const float b) 
# 974
{ 
# 975
return fminf(a, b); 
# 976
} 
# 978
static inline double min(const double a, const double b) 
# 979
{ 
# 980
return fmin(a, b); 
# 981
} 
# 983
static inline double min(const float a, const double b) 
# 984
{ 
# 985
return fmin(static_cast< double>(a), b); 
# 986
} 
# 988
static inline double min(const double a, const float b) 
# 989
{ 
# 990
return fmin(a, static_cast< double>(b)); 
# 991
} 
# 993
static inline unsigned max(const unsigned a, const unsigned b) 
# 994
{ 
# 995
return umax(a, b); 
# 996
} 
# 998
static inline unsigned max(const int a, const unsigned b) 
# 999
{ 
# 1000
return umax(static_cast< unsigned>(a), b); 
# 1001
} 
# 1003
static inline unsigned max(const unsigned a, const int b) 
# 1004
{ 
# 1005
return umax(a, static_cast< unsigned>(b)); 
# 1006
} 
# 1008
static inline long max(const long a, const long b) 
# 1009
{ 
# 1010
long retval; 
# 1015
if (sizeof(long) == sizeof(int)) { 
# 1019
retval = (static_cast< long>(max(static_cast< int>(a), static_cast< int>(b)))); 
# 1020
} else { 
# 1021
retval = (static_cast< long>(llmax(static_cast< long long>(a), static_cast< long long>(b)))); 
# 1022
}  
# 1023
return retval; 
# 1024
} 
# 1026
static inline unsigned long max(const unsigned long a, const unsigned long b) 
# 1027
{ 
# 1028
unsigned long retval; 
# 1032
if (sizeof(unsigned long) == sizeof(unsigned)) { 
# 1036
retval = (static_cast< unsigned long>(umax(static_cast< unsigned>(a), static_cast< unsigned>(b)))); 
# 1037
} else { 
# 1038
retval = (static_cast< unsigned long>(ullmax(static_cast< unsigned long long>(a), static_cast< unsigned long long>(b)))); 
# 1039
}  
# 1040
return retval; 
# 1041
} 
# 1043
static inline unsigned long max(const long a, const unsigned long b) 
# 1044
{ 
# 1045
unsigned long retval; 
# 1049
if (sizeof(unsigned long) == sizeof(unsigned)) { 
# 1053
retval = (static_cast< unsigned long>(umax(static_cast< unsigned>(a), static_cast< unsigned>(b)))); 
# 1054
} else { 
# 1055
retval = (static_cast< unsigned long>(ullmax(static_cast< unsigned long long>(a), static_cast< unsigned long long>(b)))); 
# 1056
}  
# 1057
return retval; 
# 1058
} 
# 1060
static inline unsigned long max(const unsigned long a, const long b) 
# 1061
{ 
# 1062
unsigned long retval; 
# 1066
if (sizeof(unsigned long) == sizeof(unsigned)) { 
# 1070
retval = (static_cast< unsigned long>(umax(static_cast< unsigned>(a), static_cast< unsigned>(b)))); 
# 1071
} else { 
# 1072
retval = (static_cast< unsigned long>(ullmax(static_cast< unsigned long long>(a), static_cast< unsigned long long>(b)))); 
# 1073
}  
# 1074
return retval; 
# 1075
} 
# 1077
static inline long long max(const long long a, const long long b) 
# 1078
{ 
# 1079
return llmax(a, b); 
# 1080
} 
# 1082
static inline unsigned long long max(const unsigned long long a, const unsigned long long b) 
# 1083
{ 
# 1084
return ullmax(a, b); 
# 1085
} 
# 1087
static inline unsigned long long max(const long long a, const unsigned long long b) 
# 1088
{ 
# 1089
return ullmax(static_cast< unsigned long long>(a), b); 
# 1090
} 
# 1092
static inline unsigned long long max(const unsigned long long a, const long long b) 
# 1093
{ 
# 1094
return ullmax(a, static_cast< unsigned long long>(b)); 
# 1095
} 
# 1097
static inline float max(const float a, const float b) 
# 1098
{ 
# 1099
return fmaxf(a, b); 
# 1100
} 
# 1102
static inline double max(const double a, const double b) 
# 1103
{ 
# 1104
return fmax(a, b); 
# 1105
} 
# 1107
static inline double max(const float a, const double b) 
# 1108
{ 
# 1109
return fmax(static_cast< double>(a), b); 
# 1110
} 
# 1112
static inline double max(const double a, const float b) 
# 1113
{ 
# 1114
return fmax(a, static_cast< double>(b)); 
# 1115
} 
# 1126 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/math_functions.hpp"
inline int min(const int a, const int b) 
# 1127
{ 
# 1128
return (a < b) ? a : b; 
# 1129
} 
# 1131
inline unsigned umin(const unsigned a, const unsigned b) 
# 1132
{ 
# 1133
return (a < b) ? a : b; 
# 1134
} 
# 1136
inline long long llmin(const long long a, const long long b) 
# 1137
{ 
# 1138
return (a < b) ? a : b; 
# 1139
} 
# 1141
inline unsigned long long ullmin(const unsigned long long a, const unsigned long long 
# 1142
b) 
# 1143
{ 
# 1144
return (a < b) ? a : b; 
# 1145
} 
# 1147
inline int max(const int a, const int b) 
# 1148
{ 
# 1149
return (a > b) ? a : b; 
# 1150
} 
# 1152
inline unsigned umax(const unsigned a, const unsigned b) 
# 1153
{ 
# 1154
return (a > b) ? a : b; 
# 1155
} 
# 1157
inline long long llmax(const long long a, const long long b) 
# 1158
{ 
# 1159
return (a > b) ? a : b; 
# 1160
} 
# 1162
inline unsigned long long ullmax(const unsigned long long a, const unsigned long long 
# 1163
b) 
# 1164
{ 
# 1165
return (a > b) ? a : b; 
# 1166
} 
# 91 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/device_functions.h"
extern "C" {
# 3211 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/device_functions.h"
static inline int __vimax_s32_relu(const int a, const int b); 
# 3223 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/device_functions.h"
static inline unsigned __vimax_s16x2_relu(const unsigned a, const unsigned b); 
# 3232 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/device_functions.h"
static inline int __vimin_s32_relu(const int a, const int b); 
# 3244 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/device_functions.h"
static inline unsigned __vimin_s16x2_relu(const unsigned a, const unsigned b); 
# 3253 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/device_functions.h"
static inline int __vimax3_s32(const int a, const int b, const int c); 
# 3265 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/device_functions.h"
static inline unsigned __vimax3_s16x2(const unsigned a, const unsigned b, const unsigned c); 
# 3274 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/device_functions.h"
static inline unsigned __vimax3_u32(const unsigned a, const unsigned b, const unsigned c); 
# 3286 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/device_functions.h"
static inline unsigned __vimax3_u16x2(const unsigned a, const unsigned b, const unsigned c); 
# 3295 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/device_functions.h"
static inline int __vimin3_s32(const int a, const int b, const int c); 
# 3307 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/device_functions.h"
static inline unsigned __vimin3_s16x2(const unsigned a, const unsigned b, const unsigned c); 
# 3316 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/device_functions.h"
static inline unsigned __vimin3_u32(const unsigned a, const unsigned b, const unsigned c); 
# 3328 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/device_functions.h"
static inline unsigned __vimin3_u16x2(const unsigned a, const unsigned b, const unsigned c); 
# 3337 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/device_functions.h"
static inline int __vimax3_s32_relu(const int a, const int b, const int c); 
# 3349 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/device_functions.h"
static inline unsigned __vimax3_s16x2_relu(const unsigned a, const unsigned b, const unsigned c); 
# 3358 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/device_functions.h"
static inline int __vimin3_s32_relu(const int a, const int b, const int c); 
# 3370 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/device_functions.h"
static inline unsigned __vimin3_s16x2_relu(const unsigned a, const unsigned b, const unsigned c); 
# 3379 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/device_functions.h"
static inline int __viaddmax_s32(const int a, const int b, const int c); 
# 3391 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/device_functions.h"
static inline unsigned __viaddmax_s16x2(const unsigned a, const unsigned b, const unsigned c); 
# 3400 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/device_functions.h"
static inline unsigned __viaddmax_u32(const unsigned a, const unsigned b, const unsigned c); 
# 3412 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/device_functions.h"
static inline unsigned __viaddmax_u16x2(const unsigned a, const unsigned b, const unsigned c); 
# 3421 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/device_functions.h"
static inline int __viaddmin_s32(const int a, const int b, const int c); 
# 3433 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/device_functions.h"
static inline unsigned __viaddmin_s16x2(const unsigned a, const unsigned b, const unsigned c); 
# 3442 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/device_functions.h"
static inline unsigned __viaddmin_u32(const unsigned a, const unsigned b, const unsigned c); 
# 3454 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/device_functions.h"
static inline unsigned __viaddmin_u16x2(const unsigned a, const unsigned b, const unsigned c); 
# 3464 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/device_functions.h"
static inline int __viaddmax_s32_relu(const int a, const int b, const int c); 
# 3476 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/device_functions.h"
static inline unsigned __viaddmax_s16x2_relu(const unsigned a, const unsigned b, const unsigned c); 
# 3486 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/device_functions.h"
static inline int __viaddmin_s32_relu(const int a, const int b, const int c); 
# 3498 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/device_functions.h"
static inline unsigned __viaddmin_s16x2_relu(const unsigned a, const unsigned b, const unsigned c); 
# 3507 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/device_functions.h"
static inline int __vibmax_s32(const int a, const int b, bool *const pred); 
# 3516 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/device_functions.h"
static inline unsigned __vibmax_u32(const unsigned a, const unsigned b, bool *const pred); 
# 3525 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/device_functions.h"
static inline int __vibmin_s32(const int a, const int b, bool *const pred); 
# 3534 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/device_functions.h"
static inline unsigned __vibmin_u32(const unsigned a, const unsigned b, bool *const pred); 
# 3548 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/device_functions.h"
static inline unsigned __vibmax_s16x2(const unsigned a, const unsigned b, bool *const pred_hi, bool *const pred_lo); 
# 3562 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/device_functions.h"
static inline unsigned __vibmax_u16x2(const unsigned a, const unsigned b, bool *const pred_hi, bool *const pred_lo); 
# 3576 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/device_functions.h"
static inline unsigned __vibmin_s16x2(const unsigned a, const unsigned b, bool *const pred_hi, bool *const pred_lo); 
# 3590 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/device_functions.h"
static inline unsigned __vibmin_u16x2(const unsigned a, const unsigned b, bool *const pred_hi, bool *const pred_lo); 
# 3597
}
# 102 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/device_functions.hpp"
static inline int __vimax_s32_relu(const int a, const int b) { 
# 109
int ans = max(a, b); 
# 111
return (ans > 0) ? ans : 0; 
# 113
} 
# 115
static inline unsigned __vimax_s16x2_relu(const unsigned a, const unsigned b) { 
# 123
unsigned short aU_lo = (unsigned short)(a & 65535U); 
# 124
unsigned short aU_hi = (unsigned short)(a >> 16); 
# 126
unsigned short bU_lo = (unsigned short)(b & 65535U); 
# 127
unsigned short bU_hi = (unsigned short)(b >> 16); 
# 130
short aS_lo = *((short *)(&aU_lo)); 
# 131
short aS_hi = *((short *)(&aU_hi)); 
# 133
short bS_lo = *((short *)(&bU_lo)); 
# 134
short bS_hi = *((short *)(&bU_hi)); 
# 137
short ansS_lo = (short)max(aS_lo, bS_lo); 
# 138
short ansS_hi = (short)max(aS_hi, bS_hi); 
# 141
if (ansS_lo < 0) { ansS_lo = (0); }  
# 142
if (ansS_hi < 0) { ansS_hi = (0); }  
# 145
unsigned short ansU_lo = *((unsigned short *)(&ansS_lo)); 
# 146
unsigned short ansU_hi = *((unsigned short *)(&ansS_hi)); 
# 149
unsigned ans = ((unsigned)ansU_lo) | (((unsigned)ansU_hi) << 16); 
# 151
return ans; 
# 153
} 
# 155
static inline int __vimin_s32_relu(const int a, const int b) { 
# 162
int ans = min(a, b); 
# 164
return (ans > 0) ? ans : 0; 
# 166
} 
# 168
static inline unsigned __vimin_s16x2_relu(const unsigned a, const unsigned b) { 
# 176
unsigned short aU_lo = (unsigned short)(a & 65535U); 
# 177
unsigned short aU_hi = (unsigned short)(a >> 16); 
# 179
unsigned short bU_lo = (unsigned short)(b & 65535U); 
# 180
unsigned short bU_hi = (unsigned short)(b >> 16); 
# 183
short aS_lo = *((short *)(&aU_lo)); 
# 184
short aS_hi = *((short *)(&aU_hi)); 
# 186
short bS_lo = *((short *)(&bU_lo)); 
# 187
short bS_hi = *((short *)(&bU_hi)); 
# 190
short ansS_lo = (short)min(aS_lo, bS_lo); 
# 191
short ansS_hi = (short)min(aS_hi, bS_hi); 
# 194
if (ansS_lo < 0) { ansS_lo = (0); }  
# 195
if (ansS_hi < 0) { ansS_hi = (0); }  
# 198
unsigned short ansU_lo = *((unsigned short *)(&ansS_lo)); 
# 199
unsigned short ansU_hi = *((unsigned short *)(&ansS_hi)); 
# 202
unsigned ans = ((unsigned)ansU_lo) | (((unsigned)ansU_hi) << 16); 
# 204
return ans; 
# 206
} 
# 208
static inline int __vimax3_s32(const int a, const int b, const int c) { 
# 218 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/device_functions.hpp"
return max(max(a, b), c); 
# 220
} 
# 222
static inline unsigned __vimax3_s16x2(const unsigned a, const unsigned b, const unsigned c) { 
# 234 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/device_functions.hpp"
unsigned short aU_lo = (unsigned short)(a & 65535U); 
# 235
unsigned short aU_hi = (unsigned short)(a >> 16); 
# 237
unsigned short bU_lo = (unsigned short)(b & 65535U); 
# 238
unsigned short bU_hi = (unsigned short)(b >> 16); 
# 240
unsigned short cU_lo = (unsigned short)(c & 65535U); 
# 241
unsigned short cU_hi = (unsigned short)(c >> 16); 
# 244
short aS_lo = *((short *)(&aU_lo)); 
# 245
short aS_hi = *((short *)(&aU_hi)); 
# 247
short bS_lo = *((short *)(&bU_lo)); 
# 248
short bS_hi = *((short *)(&bU_hi)); 
# 250
short cS_lo = *((short *)(&cU_lo)); 
# 251
short cS_hi = *((short *)(&cU_hi)); 
# 254
short ansS_lo = (short)max(max(aS_lo, bS_lo), cS_lo); 
# 255
short ansS_hi = (short)max(max(aS_hi, bS_hi), cS_hi); 
# 258
unsigned short ansU_lo = *((unsigned short *)(&ansS_lo)); 
# 259
unsigned short ansU_hi = *((unsigned short *)(&ansS_hi)); 
# 262
unsigned ans = ((unsigned)ansU_lo) | (((unsigned)ansU_hi) << 16); 
# 264
return ans; 
# 266
} 
# 268
static inline unsigned __vimax3_u32(const unsigned a, const unsigned b, const unsigned c) { 
# 278 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/device_functions.hpp"
return max(max(a, b), c); 
# 280
} 
# 282
static inline unsigned __vimax3_u16x2(const unsigned a, const unsigned b, const unsigned c) { 
# 293 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/device_functions.hpp"
unsigned short aU_lo = (unsigned short)(a & 65535U); 
# 294
unsigned short aU_hi = (unsigned short)(a >> 16); 
# 296
unsigned short bU_lo = (unsigned short)(b & 65535U); 
# 297
unsigned short bU_hi = (unsigned short)(b >> 16); 
# 299
unsigned short cU_lo = (unsigned short)(c & 65535U); 
# 300
unsigned short cU_hi = (unsigned short)(c >> 16); 
# 303
unsigned short ansU_lo = (unsigned short)max(max(aU_lo, bU_lo), cU_lo); 
# 304
unsigned short ansU_hi = (unsigned short)max(max(aU_hi, bU_hi), cU_hi); 
# 307
unsigned ans = ((unsigned)ansU_lo) | (((unsigned)ansU_hi) << 16); 
# 309
return ans; 
# 311
} 
# 313
static inline int __vimin3_s32(const int a, const int b, const int c) { 
# 323 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/device_functions.hpp"
return min(min(a, b), c); 
# 325
} 
# 327
static inline unsigned __vimin3_s16x2(const unsigned a, const unsigned b, const unsigned c) { 
# 338 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/device_functions.hpp"
unsigned short aU_lo = (unsigned short)(a & 65535U); 
# 339
unsigned short aU_hi = (unsigned short)(a >> 16); 
# 341
unsigned short bU_lo = (unsigned short)(b & 65535U); 
# 342
unsigned short bU_hi = (unsigned short)(b >> 16); 
# 344
unsigned short cU_lo = (unsigned short)(c & 65535U); 
# 345
unsigned short cU_hi = (unsigned short)(c >> 16); 
# 348
short aS_lo = *((short *)(&aU_lo)); 
# 349
short aS_hi = *((short *)(&aU_hi)); 
# 351
short bS_lo = *((short *)(&bU_lo)); 
# 352
short bS_hi = *((short *)(&bU_hi)); 
# 354
short cS_lo = *((short *)(&cU_lo)); 
# 355
short cS_hi = *((short *)(&cU_hi)); 
# 358
short ansS_lo = (short)min(min(aS_lo, bS_lo), cS_lo); 
# 359
short ansS_hi = (short)min(min(aS_hi, bS_hi), cS_hi); 
# 362
unsigned short ansU_lo = *((unsigned short *)(&ansS_lo)); 
# 363
unsigned short ansU_hi = *((unsigned short *)(&ansS_hi)); 
# 366
unsigned ans = ((unsigned)ansU_lo) | (((unsigned)ansU_hi) << 16); 
# 368
return ans; 
# 370
} 
# 372
static inline unsigned __vimin3_u32(const unsigned a, const unsigned b, const unsigned c) { 
# 382 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/device_functions.hpp"
return min(min(a, b), c); 
# 384
} 
# 386
static inline unsigned __vimin3_u16x2(const unsigned a, const unsigned b, const unsigned c) { 
# 397 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/device_functions.hpp"
unsigned short aU_lo = (unsigned short)(a & 65535U); 
# 398
unsigned short aU_hi = (unsigned short)(a >> 16); 
# 400
unsigned short bU_lo = (unsigned short)(b & 65535U); 
# 401
unsigned short bU_hi = (unsigned short)(b >> 16); 
# 403
unsigned short cU_lo = (unsigned short)(c & 65535U); 
# 404
unsigned short cU_hi = (unsigned short)(c >> 16); 
# 407
unsigned short ansU_lo = (unsigned short)min(min(aU_lo, bU_lo), cU_lo); 
# 408
unsigned short ansU_hi = (unsigned short)min(min(aU_hi, bU_hi), cU_hi); 
# 411
unsigned ans = ((unsigned)ansU_lo) | (((unsigned)ansU_hi) << 16); 
# 413
return ans; 
# 415
} 
# 417
static inline int __vimax3_s32_relu(const int a, const int b, const int c) { 
# 427 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/device_functions.hpp"
int ans = max(max(a, b), c); 
# 429
return (ans > 0) ? ans : 0; 
# 431
} 
# 433
static inline unsigned __vimax3_s16x2_relu(const unsigned a, const unsigned b, const unsigned c) { 
# 444 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/device_functions.hpp"
unsigned short aU_lo = (unsigned short)(a & 65535U); 
# 445
unsigned short aU_hi = (unsigned short)(a >> 16); 
# 447
unsigned short bU_lo = (unsigned short)(b & 65535U); 
# 448
unsigned short bU_hi = (unsigned short)(b >> 16); 
# 450
unsigned short cU_lo = (unsigned short)(c & 65535U); 
# 451
unsigned short cU_hi = (unsigned short)(c >> 16); 
# 454
short aS_lo = *((short *)(&aU_lo)); 
# 455
short aS_hi = *((short *)(&aU_hi)); 
# 457
short bS_lo = *((short *)(&bU_lo)); 
# 458
short bS_hi = *((short *)(&bU_hi)); 
# 460
short cS_lo = *((short *)(&cU_lo)); 
# 461
short cS_hi = *((short *)(&cU_hi)); 
# 464
short ansS_lo = (short)max(max(aS_lo, bS_lo), cS_lo); 
# 465
short ansS_hi = (short)max(max(aS_hi, bS_hi), cS_hi); 
# 468
if (ansS_lo < 0) { ansS_lo = (0); }  
# 469
if (ansS_hi < 0) { ansS_hi = (0); }  
# 472
unsigned short ansU_lo = *((unsigned short *)(&ansS_lo)); 
# 473
unsigned short ansU_hi = *((unsigned short *)(&ansS_hi)); 
# 476
unsigned ans = ((unsigned)ansU_lo) | (((unsigned)ansU_hi) << 16); 
# 478
return ans; 
# 480
} 
# 482
static inline int __vimin3_s32_relu(const int a, const int b, const int c) { 
# 492 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/device_functions.hpp"
int ans = min(min(a, b), c); 
# 494
return (ans > 0) ? ans : 0; 
# 496
} 
# 498
static inline unsigned __vimin3_s16x2_relu(const unsigned a, const unsigned b, const unsigned c) { 
# 509 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/device_functions.hpp"
unsigned short aU_lo = (unsigned short)(a & 65535U); 
# 510
unsigned short aU_hi = (unsigned short)(a >> 16); 
# 512
unsigned short bU_lo = (unsigned short)(b & 65535U); 
# 513
unsigned short bU_hi = (unsigned short)(b >> 16); 
# 515
unsigned short cU_lo = (unsigned short)(c & 65535U); 
# 516
unsigned short cU_hi = (unsigned short)(c >> 16); 
# 519
short aS_lo = *((short *)(&aU_lo)); 
# 520
short aS_hi = *((short *)(&aU_hi)); 
# 522
short bS_lo = *((short *)(&bU_lo)); 
# 523
short bS_hi = *((short *)(&bU_hi)); 
# 525
short cS_lo = *((short *)(&cU_lo)); 
# 526
short cS_hi = *((short *)(&cU_hi)); 
# 529
short ansS_lo = (short)min(min(aS_lo, bS_lo), cS_lo); 
# 530
short ansS_hi = (short)min(min(aS_hi, bS_hi), cS_hi); 
# 533
if (ansS_lo < 0) { ansS_lo = (0); }  
# 534
if (ansS_hi < 0) { ansS_hi = (0); }  
# 537
unsigned short ansU_lo = *((unsigned short *)(&ansS_lo)); 
# 538
unsigned short ansU_hi = *((unsigned short *)(&ansS_hi)); 
# 541
unsigned ans = ((unsigned)ansU_lo) | (((unsigned)ansU_hi) << 16); 
# 543
return ans; 
# 545
} 
# 547
static inline int __viaddmax_s32(const int a, const int b, const int c) { 
# 557 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/device_functions.hpp"
return max(a + b, c); 
# 559
} 
# 561
static inline unsigned __viaddmax_s16x2(const unsigned a, const unsigned b, const unsigned c) { 
# 572 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/device_functions.hpp"
unsigned short aU_lo = (unsigned short)(a & 65535U); 
# 573
unsigned short aU_hi = (unsigned short)(a >> 16); 
# 575
unsigned short bU_lo = (unsigned short)(b & 65535U); 
# 576
unsigned short bU_hi = (unsigned short)(b >> 16); 
# 578
unsigned short cU_lo = (unsigned short)(c & 65535U); 
# 579
unsigned short cU_hi = (unsigned short)(c >> 16); 
# 582
short aS_lo = *((short *)(&aU_lo)); 
# 583
short aS_hi = *((short *)(&aU_hi)); 
# 585
short bS_lo = *((short *)(&bU_lo)); 
# 586
short bS_hi = *((short *)(&bU_hi)); 
# 588
short cS_lo = *((short *)(&cU_lo)); 
# 589
short cS_hi = *((short *)(&cU_hi)); 
# 592
short ansS_lo = (short)max((short)(aS_lo + bS_lo), cS_lo); 
# 593
short ansS_hi = (short)max((short)(aS_hi + bS_hi), cS_hi); 
# 596
unsigned short ansU_lo = *((unsigned short *)(&ansS_lo)); 
# 597
unsigned short ansU_hi = *((unsigned short *)(&ansS_hi)); 
# 600
unsigned ans = ((unsigned)ansU_lo) | (((unsigned)ansU_hi) << 16); 
# 602
return ans; 
# 604
} 
# 606
static inline unsigned __viaddmax_u32(const unsigned a, const unsigned b, const unsigned c) { 
# 616 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/device_functions.hpp"
return max(a + b, c); 
# 618
} 
# 620
static inline unsigned __viaddmax_u16x2(const unsigned a, const unsigned b, const unsigned c) { 
# 631 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/device_functions.hpp"
unsigned short aU_lo = (unsigned short)(a & 65535U); 
# 632
unsigned short aU_hi = (unsigned short)(a >> 16); 
# 634
unsigned short bU_lo = (unsigned short)(b & 65535U); 
# 635
unsigned short bU_hi = (unsigned short)(b >> 16); 
# 637
unsigned short cU_lo = (unsigned short)(c & 65535U); 
# 638
unsigned short cU_hi = (unsigned short)(c >> 16); 
# 641
unsigned short ansU_lo = (unsigned short)max((unsigned short)(aU_lo + bU_lo), cU_lo); 
# 642
unsigned short ansU_hi = (unsigned short)max((unsigned short)(aU_hi + bU_hi), cU_hi); 
# 645
unsigned ans = ((unsigned)ansU_lo) | (((unsigned)ansU_hi) << 16); 
# 647
return ans; 
# 649
} 
# 651
static inline int __viaddmin_s32(const int a, const int b, const int c) { 
# 661 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/device_functions.hpp"
return min(a + b, c); 
# 663
} 
# 665
static inline unsigned __viaddmin_s16x2(const unsigned a, const unsigned b, const unsigned c) { 
# 676 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/device_functions.hpp"
unsigned short aU_lo = (unsigned short)(a & 65535U); 
# 677
unsigned short aU_hi = (unsigned short)(a >> 16); 
# 679
unsigned short bU_lo = (unsigned short)(b & 65535U); 
# 680
unsigned short bU_hi = (unsigned short)(b >> 16); 
# 682
unsigned short cU_lo = (unsigned short)(c & 65535U); 
# 683
unsigned short cU_hi = (unsigned short)(c >> 16); 
# 686
short aS_lo = *((short *)(&aU_lo)); 
# 687
short aS_hi = *((short *)(&aU_hi)); 
# 689
short bS_lo = *((short *)(&bU_lo)); 
# 690
short bS_hi = *((short *)(&bU_hi)); 
# 692
short cS_lo = *((short *)(&cU_lo)); 
# 693
short cS_hi = *((short *)(&cU_hi)); 
# 696
short ansS_lo = (short)min((short)(aS_lo + bS_lo), cS_lo); 
# 697
short ansS_hi = (short)min((short)(aS_hi + bS_hi), cS_hi); 
# 700
unsigned short ansU_lo = *((unsigned short *)(&ansS_lo)); 
# 701
unsigned short ansU_hi = *((unsigned short *)(&ansS_hi)); 
# 704
unsigned ans = ((unsigned)ansU_lo) | (((unsigned)ansU_hi) << 16); 
# 706
return ans; 
# 708
} 
# 710
static inline unsigned __viaddmin_u32(const unsigned a, const unsigned b, const unsigned c) { 
# 720 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/device_functions.hpp"
return min(a + b, c); 
# 722
} 
# 724
static inline unsigned __viaddmin_u16x2(const unsigned a, const unsigned b, const unsigned c) { 
# 735 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/device_functions.hpp"
unsigned short aU_lo = (unsigned short)(a & 65535U); 
# 736
unsigned short aU_hi = (unsigned short)(a >> 16); 
# 738
unsigned short bU_lo = (unsigned short)(b & 65535U); 
# 739
unsigned short bU_hi = (unsigned short)(b >> 16); 
# 741
unsigned short cU_lo = (unsigned short)(c & 65535U); 
# 742
unsigned short cU_hi = (unsigned short)(c >> 16); 
# 745
unsigned short ansU_lo = (unsigned short)min((unsigned short)(aU_lo + bU_lo), cU_lo); 
# 746
unsigned short ansU_hi = (unsigned short)min((unsigned short)(aU_hi + bU_hi), cU_hi); 
# 749
unsigned ans = ((unsigned)ansU_lo) | (((unsigned)ansU_hi) << 16); 
# 751
return ans; 
# 753
} 
# 755
static inline int __viaddmax_s32_relu(const int a, const int b, const int c) { 
# 765 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/device_functions.hpp"
int ans = max(a + b, c); 
# 767
return (ans > 0) ? ans : 0; 
# 769
} 
# 771
static inline unsigned __viaddmax_s16x2_relu(const unsigned a, const unsigned b, const unsigned c) { 
# 782 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/device_functions.hpp"
unsigned short aU_lo = (unsigned short)(a & 65535U); 
# 783
unsigned short aU_hi = (unsigned short)(a >> 16); 
# 785
unsigned short bU_lo = (unsigned short)(b & 65535U); 
# 786
unsigned short bU_hi = (unsigned short)(b >> 16); 
# 788
unsigned short cU_lo = (unsigned short)(c & 65535U); 
# 789
unsigned short cU_hi = (unsigned short)(c >> 16); 
# 792
short aS_lo = *((short *)(&aU_lo)); 
# 793
short aS_hi = *((short *)(&aU_hi)); 
# 795
short bS_lo = *((short *)(&bU_lo)); 
# 796
short bS_hi = *((short *)(&bU_hi)); 
# 798
short cS_lo = *((short *)(&cU_lo)); 
# 799
short cS_hi = *((short *)(&cU_hi)); 
# 802
short ansS_lo = (short)max((short)(aS_lo + bS_lo), cS_lo); 
# 803
short ansS_hi = (short)max((short)(aS_hi + bS_hi), cS_hi); 
# 805
if (ansS_lo < 0) { ansS_lo = (0); }  
# 806
if (ansS_hi < 0) { ansS_hi = (0); }  
# 809
unsigned short ansU_lo = *((unsigned short *)(&ansS_lo)); 
# 810
unsigned short ansU_hi = *((unsigned short *)(&ansS_hi)); 
# 813
unsigned ans = ((unsigned)ansU_lo) | (((unsigned)ansU_hi) << 16); 
# 815
return ans; 
# 817
} 
# 819
static inline int __viaddmin_s32_relu(const int a, const int b, const int c) { 
# 829 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/device_functions.hpp"
int ans = min(a + b, c); 
# 831
return (ans > 0) ? ans : 0; 
# 833
} 
# 835
static inline unsigned __viaddmin_s16x2_relu(const unsigned a, const unsigned b, const unsigned c) { 
# 846 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/device_functions.hpp"
unsigned short aU_lo = (unsigned short)(a & 65535U); 
# 847
unsigned short aU_hi = (unsigned short)(a >> 16); 
# 849
unsigned short bU_lo = (unsigned short)(b & 65535U); 
# 850
unsigned short bU_hi = (unsigned short)(b >> 16); 
# 852
unsigned short cU_lo = (unsigned short)(c & 65535U); 
# 853
unsigned short cU_hi = (unsigned short)(c >> 16); 
# 856
short aS_lo = *((short *)(&aU_lo)); 
# 857
short aS_hi = *((short *)(&aU_hi)); 
# 859
short bS_lo = *((short *)(&bU_lo)); 
# 860
short bS_hi = *((short *)(&bU_hi)); 
# 862
short cS_lo = *((short *)(&cU_lo)); 
# 863
short cS_hi = *((short *)(&cU_hi)); 
# 866
short ansS_lo = (short)min((short)(aS_lo + bS_lo), cS_lo); 
# 867
short ansS_hi = (short)min((short)(aS_hi + bS_hi), cS_hi); 
# 869
if (ansS_lo < 0) { ansS_lo = (0); }  
# 870
if (ansS_hi < 0) { ansS_hi = (0); }  
# 873
unsigned short ansU_lo = *((unsigned short *)(&ansS_lo)); 
# 874
unsigned short ansU_hi = *((unsigned short *)(&ansS_hi)); 
# 877
unsigned ans = ((unsigned)ansU_lo) | (((unsigned)ansU_hi) << 16); 
# 879
return ans; 
# 881
} 
# 885
static inline int __vibmax_s32(const int a, const int b, bool *const pred) { 
# 899 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/device_functions.hpp"
int ans = max(a, b); 
# 901
(*pred) = (a >= b); 
# 902
return ans; 
# 904
} 
# 906
static inline unsigned __vibmax_u32(const unsigned a, const unsigned b, bool *const pred) { 
# 920 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/device_functions.hpp"
unsigned ans = max(a, b); 
# 922
(*pred) = (a >= b); 
# 923
return ans; 
# 925
} 
# 928
static inline int __vibmin_s32(const int a, const int b, bool *const pred) { 
# 942 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/device_functions.hpp"
int ans = min(a, b); 
# 944
(*pred) = (a <= b); 
# 945
return ans; 
# 947
} 
# 950
static inline unsigned __vibmin_u32(const unsigned a, const unsigned b, bool *const pred) { 
# 964 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/device_functions.hpp"
unsigned ans = min(a, b); 
# 966
(*pred) = (a <= b); 
# 967
return ans; 
# 969
} 
# 971
static inline unsigned __vibmax_s16x2(const unsigned a, const unsigned b, bool *const pred_hi, bool *const pred_lo) { 
# 993 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/device_functions.hpp"
unsigned short aU_lo = (unsigned short)(a & 65535U); 
# 994
unsigned short aU_hi = (unsigned short)(a >> 16); 
# 996
unsigned short bU_lo = (unsigned short)(b & 65535U); 
# 997
unsigned short bU_hi = (unsigned short)(b >> 16); 
# 1000
short aS_lo = *((short *)(&aU_lo)); 
# 1001
short aS_hi = *((short *)(&aU_hi)); 
# 1003
short bS_lo = *((short *)(&bU_lo)); 
# 1004
short bS_hi = *((short *)(&bU_hi)); 
# 1007
short ansS_lo = (short)max(aS_lo, bS_lo); 
# 1008
short ansS_hi = (short)max(aS_hi, bS_hi); 
# 1010
(*pred_hi) = (aS_hi >= bS_hi); 
# 1011
(*pred_lo) = (aS_lo >= bS_lo); 
# 1014
unsigned short ansU_lo = *((unsigned short *)(&ansS_lo)); 
# 1015
unsigned short ansU_hi = *((unsigned short *)(&ansS_hi)); 
# 1018
unsigned ans = ((unsigned)ansU_lo) | (((unsigned)ansU_hi) << 16); 
# 1020
return ans; 
# 1022
} 
# 1024
static inline unsigned __vibmax_u16x2(const unsigned a, const unsigned b, bool *const pred_hi, bool *const pred_lo) { 
# 1046 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/device_functions.hpp"
unsigned short aU_lo = (unsigned short)(a & 65535U); 
# 1047
unsigned short aU_hi = (unsigned short)(a >> 16); 
# 1049
unsigned short bU_lo = (unsigned short)(b & 65535U); 
# 1050
unsigned short bU_hi = (unsigned short)(b >> 16); 
# 1053
unsigned short ansU_lo = (unsigned short)max(aU_lo, bU_lo); 
# 1054
unsigned short ansU_hi = (unsigned short)max(aU_hi, bU_hi); 
# 1056
(*pred_hi) = (aU_hi >= bU_hi); 
# 1057
(*pred_lo) = (aU_lo >= bU_lo); 
# 1060
unsigned ans = ((unsigned)ansU_lo) | (((unsigned)ansU_hi) << 16); 
# 1062
return ans; 
# 1064
} 
# 1066
static inline unsigned __vibmin_s16x2(const unsigned a, const unsigned b, bool *const pred_hi, bool *const pred_lo) { 
# 1088 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/device_functions.hpp"
unsigned short aU_lo = (unsigned short)(a & 65535U); 
# 1089
unsigned short aU_hi = (unsigned short)(a >> 16); 
# 1091
unsigned short bU_lo = (unsigned short)(b & 65535U); 
# 1092
unsigned short bU_hi = (unsigned short)(b >> 16); 
# 1095
short aS_lo = *((short *)(&aU_lo)); 
# 1096
short aS_hi = *((short *)(&aU_hi)); 
# 1098
short bS_lo = *((short *)(&bU_lo)); 
# 1099
short bS_hi = *((short *)(&bU_hi)); 
# 1102
short ansS_lo = (short)min(aS_lo, bS_lo); 
# 1103
short ansS_hi = (short)min(aS_hi, bS_hi); 
# 1105
(*pred_hi) = (aS_hi <= bS_hi); 
# 1106
(*pred_lo) = (aS_lo <= bS_lo); 
# 1109
unsigned short ansU_lo = *((unsigned short *)(&ansS_lo)); 
# 1110
unsigned short ansU_hi = *((unsigned short *)(&ansS_hi)); 
# 1113
unsigned ans = ((unsigned)ansU_lo) | (((unsigned)ansU_hi) << 16); 
# 1115
return ans; 
# 1117
} 
# 1119
static inline unsigned __vibmin_u16x2(const unsigned a, const unsigned b, bool *const pred_hi, bool *const pred_lo) { 
# 1141 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/device_functions.hpp"
unsigned short aU_lo = (unsigned short)(a & 65535U); 
# 1142
unsigned short aU_hi = (unsigned short)(a >> 16); 
# 1144
unsigned short bU_lo = (unsigned short)(b & 65535U); 
# 1145
unsigned short bU_hi = (unsigned short)(b >> 16); 
# 1148
unsigned short ansU_lo = (unsigned short)min(aU_lo, bU_lo); 
# 1149
unsigned short ansU_hi = (unsigned short)min(aU_hi, bU_hi); 
# 1151
(*pred_hi) = (aU_hi <= bU_hi); 
# 1152
(*pred_lo) = (aU_lo <= bU_lo); 
# 1155
unsigned ans = ((unsigned)ansU_lo) | (((unsigned)ansU_hi) << 16); 
# 1157
return ans; 
# 1159
} 
# 86 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/device_atomic_functions.h"
__attribute__((unused)) static inline int atomicAdd(int *address, int val) {int volatile ___ = 1;(void)address;(void)val;::exit(___);}
#if 0
# 86
{ } 
#endif
# 88 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/device_atomic_functions.h"
__attribute__((unused)) static inline unsigned atomicAdd(unsigned *address, unsigned val) {int volatile ___ = 1;(void)address;(void)val;::exit(___);}
#if 0
# 88
{ } 
#endif
# 90 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/device_atomic_functions.h"
__attribute__((unused)) static inline int atomicSub(int *address, int val) {int volatile ___ = 1;(void)address;(void)val;::exit(___);}
#if 0
# 90
{ } 
#endif
# 92 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/device_atomic_functions.h"
__attribute__((unused)) static inline unsigned atomicSub(unsigned *address, unsigned val) {int volatile ___ = 1;(void)address;(void)val;::exit(___);}
#if 0
# 92
{ } 
#endif
# 94 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/device_atomic_functions.h"
__attribute__((unused)) static inline int atomicExch(int *address, int val) {int volatile ___ = 1;(void)address;(void)val;::exit(___);}
#if 0
# 94
{ } 
#endif
# 96 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/device_atomic_functions.h"
__attribute__((unused)) static inline unsigned atomicExch(unsigned *address, unsigned val) {int volatile ___ = 1;(void)address;(void)val;::exit(___);}
#if 0
# 96
{ } 
#endif
# 98 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/device_atomic_functions.h"
__attribute__((unused)) static inline float atomicExch(float *address, float val) {int volatile ___ = 1;(void)address;(void)val;::exit(___);}
#if 0
# 98
{ } 
#endif
# 100 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/device_atomic_functions.h"
__attribute__((unused)) static inline int atomicMin(int *address, int val) {int volatile ___ = 1;(void)address;(void)val;::exit(___);}
#if 0
# 100
{ } 
#endif
# 102 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/device_atomic_functions.h"
__attribute__((unused)) static inline unsigned atomicMin(unsigned *address, unsigned val) {int volatile ___ = 1;(void)address;(void)val;::exit(___);}
#if 0
# 102
{ } 
#endif
# 104 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/device_atomic_functions.h"
__attribute__((unused)) static inline int atomicMax(int *address, int val) {int volatile ___ = 1;(void)address;(void)val;::exit(___);}
#if 0
# 104
{ } 
#endif
# 106 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/device_atomic_functions.h"
__attribute__((unused)) static inline unsigned atomicMax(unsigned *address, unsigned val) {int volatile ___ = 1;(void)address;(void)val;::exit(___);}
#if 0
# 106
{ } 
#endif
# 108 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/device_atomic_functions.h"
__attribute__((unused)) static inline unsigned atomicInc(unsigned *address, unsigned val) {int volatile ___ = 1;(void)address;(void)val;::exit(___);}
#if 0
# 108
{ } 
#endif
# 110 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/device_atomic_functions.h"
__attribute__((unused)) static inline unsigned atomicDec(unsigned *address, unsigned val) {int volatile ___ = 1;(void)address;(void)val;::exit(___);}
#if 0
# 110
{ } 
#endif
# 112 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/device_atomic_functions.h"
__attribute__((unused)) static inline int atomicAnd(int *address, int val) {int volatile ___ = 1;(void)address;(void)val;::exit(___);}
#if 0
# 112
{ } 
#endif
# 114 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/device_atomic_functions.h"
__attribute__((unused)) static inline unsigned atomicAnd(unsigned *address, unsigned val) {int volatile ___ = 1;(void)address;(void)val;::exit(___);}
#if 0
# 114
{ } 
#endif
# 116 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/device_atomic_functions.h"
__attribute__((unused)) static inline int atomicOr(int *address, int val) {int volatile ___ = 1;(void)address;(void)val;::exit(___);}
#if 0
# 116
{ } 
#endif
# 118 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/device_atomic_functions.h"
__attribute__((unused)) static inline unsigned atomicOr(unsigned *address, unsigned val) {int volatile ___ = 1;(void)address;(void)val;::exit(___);}
#if 0
# 118
{ } 
#endif
# 120 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/device_atomic_functions.h"
__attribute__((unused)) static inline int atomicXor(int *address, int val) {int volatile ___ = 1;(void)address;(void)val;::exit(___);}
#if 0
# 120
{ } 
#endif
# 122 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/device_atomic_functions.h"
__attribute__((unused)) static inline unsigned atomicXor(unsigned *address, unsigned val) {int volatile ___ = 1;(void)address;(void)val;::exit(___);}
#if 0
# 122
{ } 
#endif
# 124 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/device_atomic_functions.h"
__attribute__((unused)) static inline int atomicCAS(int *address, int compare, int val) {int volatile ___ = 1;(void)address;(void)compare;(void)val;::exit(___);}
#if 0
# 124
{ } 
#endif
# 126 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/device_atomic_functions.h"
__attribute__((unused)) static inline unsigned atomicCAS(unsigned *address, unsigned compare, unsigned val) {int volatile ___ = 1;(void)address;(void)compare;(void)val;::exit(___);}
#if 0
# 126
{ } 
#endif
# 153 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/device_atomic_functions.h"
extern "C" {
# 157
}
# 166 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/device_atomic_functions.h"
__attribute__((unused)) static inline unsigned long long atomicAdd(unsigned long long *address, unsigned long long val) {int volatile ___ = 1;(void)address;(void)val;::exit(___);}
#if 0
# 166
{ } 
#endif
# 168 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/device_atomic_functions.h"
__attribute__((unused)) static inline unsigned long long atomicExch(unsigned long long *address, unsigned long long val) {int volatile ___ = 1;(void)address;(void)val;::exit(___);}
#if 0
# 168
{ } 
#endif
# 170 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/device_atomic_functions.h"
__attribute__((unused)) static inline unsigned long long atomicCAS(unsigned long long *address, unsigned long long compare, unsigned long long val) {int volatile ___ = 1;(void)address;(void)compare;(void)val;::exit(___);}
#if 0
# 170
{ } 
#endif
# 172 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/device_atomic_functions.h"
__attribute((deprecated("__any() is deprecated in favor of __any_sync() and may be removed in a future release (Use -Wno-deprecated-declarations to suppr" "ess this warning)."))) __attribute__((unused)) static inline bool any(bool cond) {int volatile ___ = 1;(void)cond;::exit(___);}
#if 0
# 172
{ } 
#endif
# 174 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/device_atomic_functions.h"
__attribute((deprecated("__all() is deprecated in favor of __all_sync() and may be removed in a future release (Use -Wno-deprecated-declarations to suppr" "ess this warning)."))) __attribute__((unused)) static inline bool all(bool cond) {int volatile ___ = 1;(void)cond;::exit(___);}
#if 0
# 174
{ } 
#endif
# 87 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/device_double_functions.h"
extern "C" {
# 1139 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/device_double_functions.h"
}
# 1147
__attribute__((unused)) static inline double fma(double a, double b, double c, cudaRoundMode mode); 
# 1149
__attribute__((unused)) static inline double dmul(double a, double b, cudaRoundMode mode = cudaRoundNearest); 
# 1151
__attribute__((unused)) static inline double dadd(double a, double b, cudaRoundMode mode = cudaRoundNearest); 
# 1153
__attribute__((unused)) static inline double dsub(double a, double b, cudaRoundMode mode = cudaRoundNearest); 
# 1155
__attribute__((unused)) static inline int double2int(double a, cudaRoundMode mode = cudaRoundZero); 
# 1157
__attribute__((unused)) static inline unsigned double2uint(double a, cudaRoundMode mode = cudaRoundZero); 
# 1159
__attribute__((unused)) static inline long long double2ll(double a, cudaRoundMode mode = cudaRoundZero); 
# 1161
__attribute__((unused)) static inline unsigned long long double2ull(double a, cudaRoundMode mode = cudaRoundZero); 
# 1163
__attribute__((unused)) static inline double ll2double(long long a, cudaRoundMode mode = cudaRoundNearest); 
# 1165
__attribute__((unused)) static inline double ull2double(unsigned long long a, cudaRoundMode mode = cudaRoundNearest); 
# 1167
__attribute__((unused)) static inline double int2double(int a, cudaRoundMode mode = cudaRoundNearest); 
# 1169
__attribute__((unused)) static inline double uint2double(unsigned a, cudaRoundMode mode = cudaRoundNearest); 
# 1171
__attribute__((unused)) static inline double float2double(float a, cudaRoundMode mode = cudaRoundNearest); 
# 93 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/device_double_functions.hpp"
__attribute__((unused)) static inline double fma(double a, double b, double c, cudaRoundMode mode) 
# 94
{int volatile ___ = 1;(void)a;(void)b;(void)c;(void)mode;
# 99
::exit(___);}
#if 0
# 94
{ 
# 95
return (mode == (cudaRoundZero)) ? __fma_rz(a, b, c) : ((mode == (cudaRoundPosInf)) ? __fma_ru(a, b, c) : ((mode == (cudaRoundMinInf)) ? __fma_rd(a, b, c) : __fma_rn(a, b, c))); 
# 99
} 
#endif
# 101 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/device_double_functions.hpp"
__attribute__((unused)) static inline double dmul(double a, double b, cudaRoundMode mode) 
# 102
{int volatile ___ = 1;(void)a;(void)b;(void)mode;
# 107
::exit(___);}
#if 0
# 102
{ 
# 103
return (mode == (cudaRoundZero)) ? __dmul_rz(a, b) : ((mode == (cudaRoundPosInf)) ? __dmul_ru(a, b) : ((mode == (cudaRoundMinInf)) ? __dmul_rd(a, b) : __dmul_rn(a, b))); 
# 107
} 
#endif
# 109 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/device_double_functions.hpp"
__attribute__((unused)) static inline double dadd(double a, double b, cudaRoundMode mode) 
# 110
{int volatile ___ = 1;(void)a;(void)b;(void)mode;
# 115
::exit(___);}
#if 0
# 110
{ 
# 111
return (mode == (cudaRoundZero)) ? __dadd_rz(a, b) : ((mode == (cudaRoundPosInf)) ? __dadd_ru(a, b) : ((mode == (cudaRoundMinInf)) ? __dadd_rd(a, b) : __dadd_rn(a, b))); 
# 115
} 
#endif
# 117 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/device_double_functions.hpp"
__attribute__((unused)) static inline double dsub(double a, double b, cudaRoundMode mode) 
# 118
{int volatile ___ = 1;(void)a;(void)b;(void)mode;
# 123
::exit(___);}
#if 0
# 118
{ 
# 119
return (mode == (cudaRoundZero)) ? __dsub_rz(a, b) : ((mode == (cudaRoundPosInf)) ? __dsub_ru(a, b) : ((mode == (cudaRoundMinInf)) ? __dsub_rd(a, b) : __dsub_rn(a, b))); 
# 123
} 
#endif
# 125 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/device_double_functions.hpp"
__attribute__((unused)) static inline int double2int(double a, cudaRoundMode mode) 
# 126
{int volatile ___ = 1;(void)a;(void)mode;
# 131
::exit(___);}
#if 0
# 126
{ 
# 127
return (mode == (cudaRoundNearest)) ? __double2int_rn(a) : ((mode == (cudaRoundPosInf)) ? __double2int_ru(a) : ((mode == (cudaRoundMinInf)) ? __double2int_rd(a) : __double2int_rz(a))); 
# 131
} 
#endif
# 133 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/device_double_functions.hpp"
__attribute__((unused)) static inline unsigned double2uint(double a, cudaRoundMode mode) 
# 134
{int volatile ___ = 1;(void)a;(void)mode;
# 139
::exit(___);}
#if 0
# 134
{ 
# 135
return (mode == (cudaRoundNearest)) ? __double2uint_rn(a) : ((mode == (cudaRoundPosInf)) ? __double2uint_ru(a) : ((mode == (cudaRoundMinInf)) ? __double2uint_rd(a) : __double2uint_rz(a))); 
# 139
} 
#endif
# 141 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/device_double_functions.hpp"
__attribute__((unused)) static inline long long double2ll(double a, cudaRoundMode mode) 
# 142
{int volatile ___ = 1;(void)a;(void)mode;
# 147
::exit(___);}
#if 0
# 142
{ 
# 143
return (mode == (cudaRoundNearest)) ? __double2ll_rn(a) : ((mode == (cudaRoundPosInf)) ? __double2ll_ru(a) : ((mode == (cudaRoundMinInf)) ? __double2ll_rd(a) : __double2ll_rz(a))); 
# 147
} 
#endif
# 149 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/device_double_functions.hpp"
__attribute__((unused)) static inline unsigned long long double2ull(double a, cudaRoundMode mode) 
# 150
{int volatile ___ = 1;(void)a;(void)mode;
# 155
::exit(___);}
#if 0
# 150
{ 
# 151
return (mode == (cudaRoundNearest)) ? __double2ull_rn(a) : ((mode == (cudaRoundPosInf)) ? __double2ull_ru(a) : ((mode == (cudaRoundMinInf)) ? __double2ull_rd(a) : __double2ull_rz(a))); 
# 155
} 
#endif
# 157 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/device_double_functions.hpp"
__attribute__((unused)) static inline double ll2double(long long a, cudaRoundMode mode) 
# 158
{int volatile ___ = 1;(void)a;(void)mode;
# 163
::exit(___);}
#if 0
# 158
{ 
# 159
return (mode == (cudaRoundZero)) ? __ll2double_rz(a) : ((mode == (cudaRoundPosInf)) ? __ll2double_ru(a) : ((mode == (cudaRoundMinInf)) ? __ll2double_rd(a) : __ll2double_rn(a))); 
# 163
} 
#endif
# 165 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/device_double_functions.hpp"
__attribute__((unused)) static inline double ull2double(unsigned long long a, cudaRoundMode mode) 
# 166
{int volatile ___ = 1;(void)a;(void)mode;
# 171
::exit(___);}
#if 0
# 166
{ 
# 167
return (mode == (cudaRoundZero)) ? __ull2double_rz(a) : ((mode == (cudaRoundPosInf)) ? __ull2double_ru(a) : ((mode == (cudaRoundMinInf)) ? __ull2double_rd(a) : __ull2double_rn(a))); 
# 171
} 
#endif
# 173 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/device_double_functions.hpp"
__attribute__((unused)) static inline double int2double(int a, cudaRoundMode mode) 
# 174
{int volatile ___ = 1;(void)a;(void)mode;
# 176
::exit(___);}
#if 0
# 174
{ 
# 175
return (double)a; 
# 176
} 
#endif
# 178 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/device_double_functions.hpp"
__attribute__((unused)) static inline double uint2double(unsigned a, cudaRoundMode mode) 
# 179
{int volatile ___ = 1;(void)a;(void)mode;
# 181
::exit(___);}
#if 0
# 179
{ 
# 180
return (double)a; 
# 181
} 
#endif
# 183 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/device_double_functions.hpp"
__attribute__((unused)) static inline double float2double(float a, cudaRoundMode mode) 
# 184
{int volatile ___ = 1;(void)a;(void)mode;
# 186
::exit(___);}
#if 0
# 184
{ 
# 185
return (double)a; 
# 186
} 
#endif
# 96 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_20_atomic_functions.h"
__attribute__((unused)) static inline float atomicAdd(float *address, float val) {int volatile ___ = 1;(void)address;(void)val;::exit(___);}
#if 0
# 96
{ } 
#endif
# 89 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_atomic_functions.h"
__attribute__((unused)) static inline long long atomicMin(long long *address, long long val) {int volatile ___ = 1;(void)address;(void)val;::exit(___);}
#if 0
# 89
{ } 
#endif
# 91 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_atomic_functions.h"
__attribute__((unused)) static inline long long atomicMax(long long *address, long long val) {int volatile ___ = 1;(void)address;(void)val;::exit(___);}
#if 0
# 91
{ } 
#endif
# 93 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_atomic_functions.h"
__attribute__((unused)) static inline long long atomicAnd(long long *address, long long val) {int volatile ___ = 1;(void)address;(void)val;::exit(___);}
#if 0
# 93
{ } 
#endif
# 95 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_atomic_functions.h"
__attribute__((unused)) static inline long long atomicOr(long long *address, long long val) {int volatile ___ = 1;(void)address;(void)val;::exit(___);}
#if 0
# 95
{ } 
#endif
# 97 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_atomic_functions.h"
__attribute__((unused)) static inline long long atomicXor(long long *address, long long val) {int volatile ___ = 1;(void)address;(void)val;::exit(___);}
#if 0
# 97
{ } 
#endif
# 99 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_atomic_functions.h"
__attribute__((unused)) static inline unsigned long long atomicMin(unsigned long long *address, unsigned long long val) {int volatile ___ = 1;(void)address;(void)val;::exit(___);}
#if 0
# 99
{ } 
#endif
# 101 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_atomic_functions.h"
__attribute__((unused)) static inline unsigned long long atomicMax(unsigned long long *address, unsigned long long val) {int volatile ___ = 1;(void)address;(void)val;::exit(___);}
#if 0
# 101
{ } 
#endif
# 103 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_atomic_functions.h"
__attribute__((unused)) static inline unsigned long long atomicAnd(unsigned long long *address, unsigned long long val) {int volatile ___ = 1;(void)address;(void)val;::exit(___);}
#if 0
# 103
{ } 
#endif
# 105 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_atomic_functions.h"
__attribute__((unused)) static inline unsigned long long atomicOr(unsigned long long *address, unsigned long long val) {int volatile ___ = 1;(void)address;(void)val;::exit(___);}
#if 0
# 105
{ } 
#endif
# 107 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_atomic_functions.h"
__attribute__((unused)) static inline unsigned long long atomicXor(unsigned long long *address, unsigned long long val) {int volatile ___ = 1;(void)address;(void)val;::exit(___);}
#if 0
# 107
{ } 
#endif
# 90 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_60_atomic_functions.h"
__attribute__((unused)) static inline double atomicAdd(double *address, double val) {int volatile ___ = 1;(void)address;(void)val;::exit(___);}
#if 0
# 90
{ } 
#endif
# 93 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_60_atomic_functions.h"
__attribute__((unused)) static inline int atomicAdd_block(int *address, int val) {int volatile ___ = 1;(void)address;(void)val;::exit(___);}
#if 0
# 93
{ } 
#endif
# 96 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_60_atomic_functions.h"
__attribute__((unused)) static inline int atomicAdd_system(int *address, int val) {int volatile ___ = 1;(void)address;(void)val;::exit(___);}
#if 0
# 96
{ } 
#endif
# 99 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_60_atomic_functions.h"
__attribute__((unused)) static inline unsigned atomicAdd_block(unsigned *address, unsigned val) {int volatile ___ = 1;(void)address;(void)val;::exit(___);}
#if 0
# 99
{ } 
#endif
# 102 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_60_atomic_functions.h"
__attribute__((unused)) static inline unsigned atomicAdd_system(unsigned *address, unsigned val) {int volatile ___ = 1;(void)address;(void)val;::exit(___);}
#if 0
# 102
{ } 
#endif
# 105 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_60_atomic_functions.h"
__attribute__((unused)) static inline unsigned long long atomicAdd_block(unsigned long long *address, unsigned long long val) {int volatile ___ = 1;(void)address;(void)val;::exit(___);}
#if 0
# 105
{ } 
#endif
# 108 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_60_atomic_functions.h"
__attribute__((unused)) static inline unsigned long long atomicAdd_system(unsigned long long *address, unsigned long long val) {int volatile ___ = 1;(void)address;(void)val;::exit(___);}
#if 0
# 108
{ } 
#endif
# 111 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_60_atomic_functions.h"
__attribute__((unused)) static inline float atomicAdd_block(float *address, float val) {int volatile ___ = 1;(void)address;(void)val;::exit(___);}
#if 0
# 111
{ } 
#endif
# 114 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_60_atomic_functions.h"
__attribute__((unused)) static inline float atomicAdd_system(float *address, float val) {int volatile ___ = 1;(void)address;(void)val;::exit(___);}
#if 0
# 114
{ } 
#endif
# 117 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_60_atomic_functions.h"
__attribute__((unused)) static inline double atomicAdd_block(double *address, double val) {int volatile ___ = 1;(void)address;(void)val;::exit(___);}
#if 0
# 117
{ } 
#endif
# 120 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_60_atomic_functions.h"
__attribute__((unused)) static inline double atomicAdd_system(double *address, double val) {int volatile ___ = 1;(void)address;(void)val;::exit(___);}
#if 0
# 120
{ } 
#endif
# 123 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_60_atomic_functions.h"
__attribute__((unused)) static inline int atomicSub_block(int *address, int val) {int volatile ___ = 1;(void)address;(void)val;::exit(___);}
#if 0
# 123
{ } 
#endif
# 126 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_60_atomic_functions.h"
__attribute__((unused)) static inline int atomicSub_system(int *address, int val) {int volatile ___ = 1;(void)address;(void)val;::exit(___);}
#if 0
# 126
{ } 
#endif
# 129 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_60_atomic_functions.h"
__attribute__((unused)) static inline unsigned atomicSub_block(unsigned *address, unsigned val) {int volatile ___ = 1;(void)address;(void)val;::exit(___);}
#if 0
# 129
{ } 
#endif
# 132 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_60_atomic_functions.h"
__attribute__((unused)) static inline unsigned atomicSub_system(unsigned *address, unsigned val) {int volatile ___ = 1;(void)address;(void)val;::exit(___);}
#if 0
# 132
{ } 
#endif
# 135 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_60_atomic_functions.h"
__attribute__((unused)) static inline int atomicExch_block(int *address, int val) {int volatile ___ = 1;(void)address;(void)val;::exit(___);}
#if 0
# 135
{ } 
#endif
# 138 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_60_atomic_functions.h"
__attribute__((unused)) static inline int atomicExch_system(int *address, int val) {int volatile ___ = 1;(void)address;(void)val;::exit(___);}
#if 0
# 138
{ } 
#endif
# 141 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_60_atomic_functions.h"
__attribute__((unused)) static inline unsigned atomicExch_block(unsigned *address, unsigned val) {int volatile ___ = 1;(void)address;(void)val;::exit(___);}
#if 0
# 141
{ } 
#endif
# 144 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_60_atomic_functions.h"
__attribute__((unused)) static inline unsigned atomicExch_system(unsigned *address, unsigned val) {int volatile ___ = 1;(void)address;(void)val;::exit(___);}
#if 0
# 144
{ } 
#endif
# 147 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_60_atomic_functions.h"
__attribute__((unused)) static inline unsigned long long atomicExch_block(unsigned long long *address, unsigned long long val) {int volatile ___ = 1;(void)address;(void)val;::exit(___);}
#if 0
# 147
{ } 
#endif
# 150 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_60_atomic_functions.h"
__attribute__((unused)) static inline unsigned long long atomicExch_system(unsigned long long *address, unsigned long long val) {int volatile ___ = 1;(void)address;(void)val;::exit(___);}
#if 0
# 150
{ } 
#endif
# 153 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_60_atomic_functions.h"
__attribute__((unused)) static inline float atomicExch_block(float *address, float val) {int volatile ___ = 1;(void)address;(void)val;::exit(___);}
#if 0
# 153
{ } 
#endif
# 156 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_60_atomic_functions.h"
__attribute__((unused)) static inline float atomicExch_system(float *address, float val) {int volatile ___ = 1;(void)address;(void)val;::exit(___);}
#if 0
# 156
{ } 
#endif
# 159 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_60_atomic_functions.h"
__attribute__((unused)) static inline int atomicMin_block(int *address, int val) {int volatile ___ = 1;(void)address;(void)val;::exit(___);}
#if 0
# 159
{ } 
#endif
# 162 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_60_atomic_functions.h"
__attribute__((unused)) static inline int atomicMin_system(int *address, int val) {int volatile ___ = 1;(void)address;(void)val;::exit(___);}
#if 0
# 162
{ } 
#endif
# 165 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_60_atomic_functions.h"
__attribute__((unused)) static inline long long atomicMin_block(long long *address, long long val) {int volatile ___ = 1;(void)address;(void)val;::exit(___);}
#if 0
# 165
{ } 
#endif
# 168 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_60_atomic_functions.h"
__attribute__((unused)) static inline long long atomicMin_system(long long *address, long long val) {int volatile ___ = 1;(void)address;(void)val;::exit(___);}
#if 0
# 168
{ } 
#endif
# 171 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_60_atomic_functions.h"
__attribute__((unused)) static inline unsigned atomicMin_block(unsigned *address, unsigned val) {int volatile ___ = 1;(void)address;(void)val;::exit(___);}
#if 0
# 171
{ } 
#endif
# 174 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_60_atomic_functions.h"
__attribute__((unused)) static inline unsigned atomicMin_system(unsigned *address, unsigned val) {int volatile ___ = 1;(void)address;(void)val;::exit(___);}
#if 0
# 174
{ } 
#endif
# 177 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_60_atomic_functions.h"
__attribute__((unused)) static inline unsigned long long atomicMin_block(unsigned long long *address, unsigned long long val) {int volatile ___ = 1;(void)address;(void)val;::exit(___);}
#if 0
# 177
{ } 
#endif
# 180 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_60_atomic_functions.h"
__attribute__((unused)) static inline unsigned long long atomicMin_system(unsigned long long *address, unsigned long long val) {int volatile ___ = 1;(void)address;(void)val;::exit(___);}
#if 0
# 180
{ } 
#endif
# 183 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_60_atomic_functions.h"
__attribute__((unused)) static inline int atomicMax_block(int *address, int val) {int volatile ___ = 1;(void)address;(void)val;::exit(___);}
#if 0
# 183
{ } 
#endif
# 186 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_60_atomic_functions.h"
__attribute__((unused)) static inline int atomicMax_system(int *address, int val) {int volatile ___ = 1;(void)address;(void)val;::exit(___);}
#if 0
# 186
{ } 
#endif
# 189 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_60_atomic_functions.h"
__attribute__((unused)) static inline long long atomicMax_block(long long *address, long long val) {int volatile ___ = 1;(void)address;(void)val;::exit(___);}
#if 0
# 189
{ } 
#endif
# 192 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_60_atomic_functions.h"
__attribute__((unused)) static inline long long atomicMax_system(long long *address, long long val) {int volatile ___ = 1;(void)address;(void)val;::exit(___);}
#if 0
# 192
{ } 
#endif
# 195 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_60_atomic_functions.h"
__attribute__((unused)) static inline unsigned atomicMax_block(unsigned *address, unsigned val) {int volatile ___ = 1;(void)address;(void)val;::exit(___);}
#if 0
# 195
{ } 
#endif
# 198 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_60_atomic_functions.h"
__attribute__((unused)) static inline unsigned atomicMax_system(unsigned *address, unsigned val) {int volatile ___ = 1;(void)address;(void)val;::exit(___);}
#if 0
# 198
{ } 
#endif
# 201 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_60_atomic_functions.h"
__attribute__((unused)) static inline unsigned long long atomicMax_block(unsigned long long *address, unsigned long long val) {int volatile ___ = 1;(void)address;(void)val;::exit(___);}
#if 0
# 201
{ } 
#endif
# 204 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_60_atomic_functions.h"
__attribute__((unused)) static inline unsigned long long atomicMax_system(unsigned long long *address, unsigned long long val) {int volatile ___ = 1;(void)address;(void)val;::exit(___);}
#if 0
# 204
{ } 
#endif
# 207 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_60_atomic_functions.h"
__attribute__((unused)) static inline unsigned atomicInc_block(unsigned *address, unsigned val) {int volatile ___ = 1;(void)address;(void)val;::exit(___);}
#if 0
# 207
{ } 
#endif
# 210 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_60_atomic_functions.h"
__attribute__((unused)) static inline unsigned atomicInc_system(unsigned *address, unsigned val) {int volatile ___ = 1;(void)address;(void)val;::exit(___);}
#if 0
# 210
{ } 
#endif
# 213 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_60_atomic_functions.h"
__attribute__((unused)) static inline unsigned atomicDec_block(unsigned *address, unsigned val) {int volatile ___ = 1;(void)address;(void)val;::exit(___);}
#if 0
# 213
{ } 
#endif
# 216 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_60_atomic_functions.h"
__attribute__((unused)) static inline unsigned atomicDec_system(unsigned *address, unsigned val) {int volatile ___ = 1;(void)address;(void)val;::exit(___);}
#if 0
# 216
{ } 
#endif
# 219 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_60_atomic_functions.h"
__attribute__((unused)) static inline int atomicCAS_block(int *address, int compare, int val) {int volatile ___ = 1;(void)address;(void)compare;(void)val;::exit(___);}
#if 0
# 219
{ } 
#endif
# 222 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_60_atomic_functions.h"
__attribute__((unused)) static inline int atomicCAS_system(int *address, int compare, int val) {int volatile ___ = 1;(void)address;(void)compare;(void)val;::exit(___);}
#if 0
# 222
{ } 
#endif
# 225 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_60_atomic_functions.h"
__attribute__((unused)) static inline unsigned atomicCAS_block(unsigned *address, unsigned compare, unsigned 
# 226
val) {int volatile ___ = 1;(void)address;(void)compare;(void)val;::exit(___);}
#if 0
# 226
{ } 
#endif
# 229 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_60_atomic_functions.h"
__attribute__((unused)) static inline unsigned atomicCAS_system(unsigned *address, unsigned compare, unsigned 
# 230
val) {int volatile ___ = 1;(void)address;(void)compare;(void)val;::exit(___);}
#if 0
# 230
{ } 
#endif
# 233 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_60_atomic_functions.h"
__attribute__((unused)) static inline unsigned long long atomicCAS_block(unsigned long long *address, unsigned long long 
# 234
compare, unsigned long long 
# 235
val) {int volatile ___ = 1;(void)address;(void)compare;(void)val;::exit(___);}
#if 0
# 235
{ } 
#endif
# 238 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_60_atomic_functions.h"
__attribute__((unused)) static inline unsigned long long atomicCAS_system(unsigned long long *address, unsigned long long 
# 239
compare, unsigned long long 
# 240
val) {int volatile ___ = 1;(void)address;(void)compare;(void)val;::exit(___);}
#if 0
# 240
{ } 
#endif
# 243 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_60_atomic_functions.h"
__attribute__((unused)) static inline int atomicAnd_block(int *address, int val) {int volatile ___ = 1;(void)address;(void)val;::exit(___);}
#if 0
# 243
{ } 
#endif
# 246 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_60_atomic_functions.h"
__attribute__((unused)) static inline int atomicAnd_system(int *address, int val) {int volatile ___ = 1;(void)address;(void)val;::exit(___);}
#if 0
# 246
{ } 
#endif
# 249 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_60_atomic_functions.h"
__attribute__((unused)) static inline long long atomicAnd_block(long long *address, long long val) {int volatile ___ = 1;(void)address;(void)val;::exit(___);}
#if 0
# 249
{ } 
#endif
# 252 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_60_atomic_functions.h"
__attribute__((unused)) static inline long long atomicAnd_system(long long *address, long long val) {int volatile ___ = 1;(void)address;(void)val;::exit(___);}
#if 0
# 252
{ } 
#endif
# 255 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_60_atomic_functions.h"
__attribute__((unused)) static inline unsigned atomicAnd_block(unsigned *address, unsigned val) {int volatile ___ = 1;(void)address;(void)val;::exit(___);}
#if 0
# 255
{ } 
#endif
# 258 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_60_atomic_functions.h"
__attribute__((unused)) static inline unsigned atomicAnd_system(unsigned *address, unsigned val) {int volatile ___ = 1;(void)address;(void)val;::exit(___);}
#if 0
# 258
{ } 
#endif
# 261 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_60_atomic_functions.h"
__attribute__((unused)) static inline unsigned long long atomicAnd_block(unsigned long long *address, unsigned long long val) {int volatile ___ = 1;(void)address;(void)val;::exit(___);}
#if 0
# 261
{ } 
#endif
# 264 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_60_atomic_functions.h"
__attribute__((unused)) static inline unsigned long long atomicAnd_system(unsigned long long *address, unsigned long long val) {int volatile ___ = 1;(void)address;(void)val;::exit(___);}
#if 0
# 264
{ } 
#endif
# 267 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_60_atomic_functions.h"
__attribute__((unused)) static inline int atomicOr_block(int *address, int val) {int volatile ___ = 1;(void)address;(void)val;::exit(___);}
#if 0
# 267
{ } 
#endif
# 270 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_60_atomic_functions.h"
__attribute__((unused)) static inline int atomicOr_system(int *address, int val) {int volatile ___ = 1;(void)address;(void)val;::exit(___);}
#if 0
# 270
{ } 
#endif
# 273 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_60_atomic_functions.h"
__attribute__((unused)) static inline long long atomicOr_block(long long *address, long long val) {int volatile ___ = 1;(void)address;(void)val;::exit(___);}
#if 0
# 273
{ } 
#endif
# 276 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_60_atomic_functions.h"
__attribute__((unused)) static inline long long atomicOr_system(long long *address, long long val) {int volatile ___ = 1;(void)address;(void)val;::exit(___);}
#if 0
# 276
{ } 
#endif
# 279 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_60_atomic_functions.h"
__attribute__((unused)) static inline unsigned atomicOr_block(unsigned *address, unsigned val) {int volatile ___ = 1;(void)address;(void)val;::exit(___);}
#if 0
# 279
{ } 
#endif
# 282 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_60_atomic_functions.h"
__attribute__((unused)) static inline unsigned atomicOr_system(unsigned *address, unsigned val) {int volatile ___ = 1;(void)address;(void)val;::exit(___);}
#if 0
# 282
{ } 
#endif
# 285 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_60_atomic_functions.h"
__attribute__((unused)) static inline unsigned long long atomicOr_block(unsigned long long *address, unsigned long long val) {int volatile ___ = 1;(void)address;(void)val;::exit(___);}
#if 0
# 285
{ } 
#endif
# 288 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_60_atomic_functions.h"
__attribute__((unused)) static inline unsigned long long atomicOr_system(unsigned long long *address, unsigned long long val) {int volatile ___ = 1;(void)address;(void)val;::exit(___);}
#if 0
# 288
{ } 
#endif
# 291 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_60_atomic_functions.h"
__attribute__((unused)) static inline int atomicXor_block(int *address, int val) {int volatile ___ = 1;(void)address;(void)val;::exit(___);}
#if 0
# 291
{ } 
#endif
# 294 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_60_atomic_functions.h"
__attribute__((unused)) static inline int atomicXor_system(int *address, int val) {int volatile ___ = 1;(void)address;(void)val;::exit(___);}
#if 0
# 294
{ } 
#endif
# 297 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_60_atomic_functions.h"
__attribute__((unused)) static inline long long atomicXor_block(long long *address, long long val) {int volatile ___ = 1;(void)address;(void)val;::exit(___);}
#if 0
# 297
{ } 
#endif
# 300 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_60_atomic_functions.h"
__attribute__((unused)) static inline long long atomicXor_system(long long *address, long long val) {int volatile ___ = 1;(void)address;(void)val;::exit(___);}
#if 0
# 300
{ } 
#endif
# 303 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_60_atomic_functions.h"
__attribute__((unused)) static inline unsigned atomicXor_block(unsigned *address, unsigned val) {int volatile ___ = 1;(void)address;(void)val;::exit(___);}
#if 0
# 303
{ } 
#endif
# 306 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_60_atomic_functions.h"
__attribute__((unused)) static inline unsigned atomicXor_system(unsigned *address, unsigned val) {int volatile ___ = 1;(void)address;(void)val;::exit(___);}
#if 0
# 306
{ } 
#endif
# 309 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_60_atomic_functions.h"
__attribute__((unused)) static inline unsigned long long atomicXor_block(unsigned long long *address, unsigned long long val) {int volatile ___ = 1;(void)address;(void)val;::exit(___);}
#if 0
# 309
{ } 
#endif
# 312 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_60_atomic_functions.h"
__attribute__((unused)) static inline unsigned long long atomicXor_system(unsigned long long *address, unsigned long long val) {int volatile ___ = 1;(void)address;(void)val;::exit(___);}
#if 0
# 312
{ } 
#endif
# 97 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_20_intrinsics.h"
extern "C" {
# 1510 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_20_intrinsics.h"
}
# 1522 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_20_intrinsics.h"
__attribute((deprecated("__ballot() is deprecated in favor of __ballot_sync() and may be removed in a future release (Use -Wno-deprecated-declarations to" " suppress this warning)."))) __attribute__((unused)) static inline unsigned ballot(bool pred) {int volatile ___ = 1;(void)pred;::exit(___);}
#if 0
# 1522
{ } 
#endif
# 1524 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_20_intrinsics.h"
__attribute__((unused)) static inline int syncthreads_count(bool pred) {int volatile ___ = 1;(void)pred;::exit(___);}
#if 0
# 1524
{ } 
#endif
# 1526 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_20_intrinsics.h"
__attribute__((unused)) static inline bool syncthreads_and(bool pred) {int volatile ___ = 1;(void)pred;::exit(___);}
#if 0
# 1526
{ } 
#endif
# 1528 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_20_intrinsics.h"
__attribute__((unused)) static inline bool syncthreads_or(bool pred) {int volatile ___ = 1;(void)pred;::exit(___);}
#if 0
# 1528
{ } 
#endif
# 1533 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_20_intrinsics.h"
__attribute__((unused)) static inline unsigned __isGlobal(const void *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 1533
{ } 
#endif
# 1534 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_20_intrinsics.h"
__attribute__((unused)) static inline unsigned __isShared(const void *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 1534
{ } 
#endif
# 1535 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_20_intrinsics.h"
__attribute__((unused)) static inline unsigned __isConstant(const void *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 1535
{ } 
#endif
# 1536 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_20_intrinsics.h"
__attribute__((unused)) static inline unsigned __isLocal(const void *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 1536
{ } 
#endif
# 1538 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_20_intrinsics.h"
__attribute__((unused)) static inline unsigned __isGridConstant(const void *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 1538
{ } 
#endif
# 1540 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_20_intrinsics.h"
__attribute__((unused)) static inline size_t __cvta_generic_to_global(const void *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 1540
{ } 
#endif
# 1541 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_20_intrinsics.h"
__attribute__((unused)) static inline size_t __cvta_generic_to_shared(const void *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 1541
{ } 
#endif
# 1542 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_20_intrinsics.h"
__attribute__((unused)) static inline size_t __cvta_generic_to_constant(const void *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 1542
{ } 
#endif
# 1543 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_20_intrinsics.h"
__attribute__((unused)) static inline size_t __cvta_generic_to_local(const void *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 1543
{ } 
#endif
# 1545 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_20_intrinsics.h"
__attribute__((unused)) static inline size_t __cvta_generic_to_grid_constant(const void *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 1545
{ } 
#endif
# 1548 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_20_intrinsics.h"
__attribute__((unused)) static inline void *__cvta_global_to_generic(size_t rawbits) {int volatile ___ = 1;(void)rawbits;::exit(___);}
#if 0
# 1548
{ } 
#endif
# 1549 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_20_intrinsics.h"
__attribute__((unused)) static inline void *__cvta_shared_to_generic(size_t rawbits) {int volatile ___ = 1;(void)rawbits;::exit(___);}
#if 0
# 1549
{ } 
#endif
# 1550 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_20_intrinsics.h"
__attribute__((unused)) static inline void *__cvta_constant_to_generic(size_t rawbits) {int volatile ___ = 1;(void)rawbits;::exit(___);}
#if 0
# 1550
{ } 
#endif
# 1551 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_20_intrinsics.h"
__attribute__((unused)) static inline void *__cvta_local_to_generic(size_t rawbits) {int volatile ___ = 1;(void)rawbits;::exit(___);}
#if 0
# 1551
{ } 
#endif
# 1553 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_20_intrinsics.h"
__attribute__((unused)) static inline void *__cvta_grid_constant_to_generic(size_t rawbits) {int volatile ___ = 1;(void)rawbits;::exit(___);}
#if 0
# 1553
{ } 
#endif
# 123 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_30_intrinsics.h"
__attribute__((unused)) static inline unsigned __fns(unsigned mask, unsigned base, int offset) {int volatile ___ = 1;(void)mask;(void)base;(void)offset;::exit(___);}
#if 0
# 123
{ } 
#endif
# 124 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_30_intrinsics.h"
__attribute__((unused)) static inline void __barrier_sync(unsigned id) {int volatile ___ = 1;(void)id;::exit(___);}
#if 0
# 124
{ } 
#endif
# 125 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_30_intrinsics.h"
__attribute__((unused)) static inline void __barrier_sync_count(unsigned id, unsigned cnt) {int volatile ___ = 1;(void)id;(void)cnt;::exit(___);}
#if 0
# 125
{ } 
#endif
# 126 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_30_intrinsics.h"
__attribute__((unused)) static inline void __syncwarp(unsigned mask = 4294967295U) {int volatile ___ = 1;(void)mask;::exit(___);}
#if 0
# 126
{ } 
#endif
# 127 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_30_intrinsics.h"
__attribute__((unused)) static inline int __all_sync(unsigned mask, int pred) {int volatile ___ = 1;(void)mask;(void)pred;::exit(___);}
#if 0
# 127
{ } 
#endif
# 128 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_30_intrinsics.h"
__attribute__((unused)) static inline int __any_sync(unsigned mask, int pred) {int volatile ___ = 1;(void)mask;(void)pred;::exit(___);}
#if 0
# 128
{ } 
#endif
# 129 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_30_intrinsics.h"
__attribute__((unused)) static inline int __uni_sync(unsigned mask, int pred) {int volatile ___ = 1;(void)mask;(void)pred;::exit(___);}
#if 0
# 129
{ } 
#endif
# 130 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_30_intrinsics.h"
__attribute__((unused)) static inline unsigned __ballot_sync(unsigned mask, int pred) {int volatile ___ = 1;(void)mask;(void)pred;::exit(___);}
#if 0
# 130
{ } 
#endif
# 131 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_30_intrinsics.h"
__attribute__((unused)) static inline unsigned __activemask() {int volatile ___ = 1;::exit(___);}
#if 0
# 131
{ } 
#endif
# 140 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_30_intrinsics.h"
__attribute((deprecated("__shfl() is deprecated in favor of __shfl_sync() and may be removed in a future release (Use -Wno-deprecated-declarations to sup" "press this warning)."))) __attribute__((unused)) static inline int __shfl(int var, int srcLane, int width = 32) {int volatile ___ = 1;(void)var;(void)srcLane;(void)width;::exit(___);}
#if 0
# 140
{ } 
#endif
# 141 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_30_intrinsics.h"
__attribute((deprecated("__shfl() is deprecated in favor of __shfl_sync() and may be removed in a future release (Use -Wno-deprecated-declarations to sup" "press this warning)."))) __attribute__((unused)) static inline unsigned __shfl(unsigned var, int srcLane, int width = 32) {int volatile ___ = 1;(void)var;(void)srcLane;(void)width;::exit(___);}
#if 0
# 141
{ } 
#endif
# 142 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_30_intrinsics.h"
__attribute((deprecated("__shfl_up() is deprecated in favor of __shfl_up_sync() and may be removed in a future release (Use -Wno-deprecated-declarations " "to suppress this warning)."))) __attribute__((unused)) static inline int __shfl_up(int var, unsigned delta, int width = 32) {int volatile ___ = 1;(void)var;(void)delta;(void)width;::exit(___);}
#if 0
# 142
{ } 
#endif
# 143 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_30_intrinsics.h"
__attribute((deprecated("__shfl_up() is deprecated in favor of __shfl_up_sync() and may be removed in a future release (Use -Wno-deprecated-declarations " "to suppress this warning)."))) __attribute__((unused)) static inline unsigned __shfl_up(unsigned var, unsigned delta, int width = 32) {int volatile ___ = 1;(void)var;(void)delta;(void)width;::exit(___);}
#if 0
# 143
{ } 
#endif
# 144 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_30_intrinsics.h"
__attribute((deprecated("__shfl_down() is deprecated in favor of __shfl_down_sync() and may be removed in a future release (Use -Wno-deprecated-declarati" "ons to suppress this warning)."))) __attribute__((unused)) static inline int __shfl_down(int var, unsigned delta, int width = 32) {int volatile ___ = 1;(void)var;(void)delta;(void)width;::exit(___);}
#if 0
# 144
{ } 
#endif
# 145 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_30_intrinsics.h"
__attribute((deprecated("__shfl_down() is deprecated in favor of __shfl_down_sync() and may be removed in a future release (Use -Wno-deprecated-declarati" "ons to suppress this warning)."))) __attribute__((unused)) static inline unsigned __shfl_down(unsigned var, unsigned delta, int width = 32) {int volatile ___ = 1;(void)var;(void)delta;(void)width;::exit(___);}
#if 0
# 145
{ } 
#endif
# 146 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_30_intrinsics.h"
__attribute((deprecated("__shfl_xor() is deprecated in favor of __shfl_xor_sync() and may be removed in a future release (Use -Wno-deprecated-declaration" "s to suppress this warning)."))) __attribute__((unused)) static inline int __shfl_xor(int var, int laneMask, int width = 32) {int volatile ___ = 1;(void)var;(void)laneMask;(void)width;::exit(___);}
#if 0
# 146
{ } 
#endif
# 147 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_30_intrinsics.h"
__attribute((deprecated("__shfl_xor() is deprecated in favor of __shfl_xor_sync() and may be removed in a future release (Use -Wno-deprecated-declaration" "s to suppress this warning)."))) __attribute__((unused)) static inline unsigned __shfl_xor(unsigned var, int laneMask, int width = 32) {int volatile ___ = 1;(void)var;(void)laneMask;(void)width;::exit(___);}
#if 0
# 147
{ } 
#endif
# 148 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_30_intrinsics.h"
__attribute((deprecated("__shfl() is deprecated in favor of __shfl_sync() and may be removed in a future release (Use -Wno-deprecated-declarations to sup" "press this warning)."))) __attribute__((unused)) static inline float __shfl(float var, int srcLane, int width = 32) {int volatile ___ = 1;(void)var;(void)srcLane;(void)width;::exit(___);}
#if 0
# 148
{ } 
#endif
# 149 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_30_intrinsics.h"
__attribute((deprecated("__shfl_up() is deprecated in favor of __shfl_up_sync() and may be removed in a future release (Use -Wno-deprecated-declarations " "to suppress this warning)."))) __attribute__((unused)) static inline float __shfl_up(float var, unsigned delta, int width = 32) {int volatile ___ = 1;(void)var;(void)delta;(void)width;::exit(___);}
#if 0
# 149
{ } 
#endif
# 150 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_30_intrinsics.h"
__attribute((deprecated("__shfl_down() is deprecated in favor of __shfl_down_sync() and may be removed in a future release (Use -Wno-deprecated-declarati" "ons to suppress this warning)."))) __attribute__((unused)) static inline float __shfl_down(float var, unsigned delta, int width = 32) {int volatile ___ = 1;(void)var;(void)delta;(void)width;::exit(___);}
#if 0
# 150
{ } 
#endif
# 151 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_30_intrinsics.h"
__attribute((deprecated("__shfl_xor() is deprecated in favor of __shfl_xor_sync() and may be removed in a future release (Use -Wno-deprecated-declaration" "s to suppress this warning)."))) __attribute__((unused)) static inline float __shfl_xor(float var, int laneMask, int width = 32) {int volatile ___ = 1;(void)var;(void)laneMask;(void)width;::exit(___);}
#if 0
# 151
{ } 
#endif
# 154 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_30_intrinsics.h"
__attribute__((unused)) static inline int __shfl_sync(unsigned mask, int var, int srcLane, int width = 32) {int volatile ___ = 1;(void)mask;(void)var;(void)srcLane;(void)width;::exit(___);}
#if 0
# 154
{ } 
#endif
# 155 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_30_intrinsics.h"
__attribute__((unused)) static inline unsigned __shfl_sync(unsigned mask, unsigned var, int srcLane, int width = 32) {int volatile ___ = 1;(void)mask;(void)var;(void)srcLane;(void)width;::exit(___);}
#if 0
# 155
{ } 
#endif
# 156 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_30_intrinsics.h"
__attribute__((unused)) static inline int __shfl_up_sync(unsigned mask, int var, unsigned delta, int width = 32) {int volatile ___ = 1;(void)mask;(void)var;(void)delta;(void)width;::exit(___);}
#if 0
# 156
{ } 
#endif
# 157 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_30_intrinsics.h"
__attribute__((unused)) static inline unsigned __shfl_up_sync(unsigned mask, unsigned var, unsigned delta, int width = 32) {int volatile ___ = 1;(void)mask;(void)var;(void)delta;(void)width;::exit(___);}
#if 0
# 157
{ } 
#endif
# 158 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_30_intrinsics.h"
__attribute__((unused)) static inline int __shfl_down_sync(unsigned mask, int var, unsigned delta, int width = 32) {int volatile ___ = 1;(void)mask;(void)var;(void)delta;(void)width;::exit(___);}
#if 0
# 158
{ } 
#endif
# 159 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_30_intrinsics.h"
__attribute__((unused)) static inline unsigned __shfl_down_sync(unsigned mask, unsigned var, unsigned delta, int width = 32) {int volatile ___ = 1;(void)mask;(void)var;(void)delta;(void)width;::exit(___);}
#if 0
# 159
{ } 
#endif
# 160 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_30_intrinsics.h"
__attribute__((unused)) static inline int __shfl_xor_sync(unsigned mask, int var, int laneMask, int width = 32) {int volatile ___ = 1;(void)mask;(void)var;(void)laneMask;(void)width;::exit(___);}
#if 0
# 160
{ } 
#endif
# 161 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_30_intrinsics.h"
__attribute__((unused)) static inline unsigned __shfl_xor_sync(unsigned mask, unsigned var, int laneMask, int width = 32) {int volatile ___ = 1;(void)mask;(void)var;(void)laneMask;(void)width;::exit(___);}
#if 0
# 161
{ } 
#endif
# 162 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_30_intrinsics.h"
__attribute__((unused)) static inline float __shfl_sync(unsigned mask, float var, int srcLane, int width = 32) {int volatile ___ = 1;(void)mask;(void)var;(void)srcLane;(void)width;::exit(___);}
#if 0
# 162
{ } 
#endif
# 163 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_30_intrinsics.h"
__attribute__((unused)) static inline float __shfl_up_sync(unsigned mask, float var, unsigned delta, int width = 32) {int volatile ___ = 1;(void)mask;(void)var;(void)delta;(void)width;::exit(___);}
#if 0
# 163
{ } 
#endif
# 164 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_30_intrinsics.h"
__attribute__((unused)) static inline float __shfl_down_sync(unsigned mask, float var, unsigned delta, int width = 32) {int volatile ___ = 1;(void)mask;(void)var;(void)delta;(void)width;::exit(___);}
#if 0
# 164
{ } 
#endif
# 165 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_30_intrinsics.h"
__attribute__((unused)) static inline float __shfl_xor_sync(unsigned mask, float var, int laneMask, int width = 32) {int volatile ___ = 1;(void)mask;(void)var;(void)laneMask;(void)width;::exit(___);}
#if 0
# 165
{ } 
#endif
# 169 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_30_intrinsics.h"
__attribute((deprecated("__shfl() is deprecated in favor of __shfl_sync() and may be removed in a future release (Use -Wno-deprecated-declarations to sup" "press this warning)."))) __attribute__((unused)) static inline unsigned long long __shfl(unsigned long long var, int srcLane, int width = 32) {int volatile ___ = 1;(void)var;(void)srcLane;(void)width;::exit(___);}
#if 0
# 169
{ } 
#endif
# 170 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_30_intrinsics.h"
__attribute((deprecated("__shfl() is deprecated in favor of __shfl_sync() and may be removed in a future release (Use -Wno-deprecated-declarations to sup" "press this warning)."))) __attribute__((unused)) static inline long long __shfl(long long var, int srcLane, int width = 32) {int volatile ___ = 1;(void)var;(void)srcLane;(void)width;::exit(___);}
#if 0
# 170
{ } 
#endif
# 171 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_30_intrinsics.h"
__attribute((deprecated("__shfl_up() is deprecated in favor of __shfl_up_sync() and may be removed in a future release (Use -Wno-deprecated-declarations " "to suppress this warning)."))) __attribute__((unused)) static inline long long __shfl_up(long long var, unsigned delta, int width = 32) {int volatile ___ = 1;(void)var;(void)delta;(void)width;::exit(___);}
#if 0
# 171
{ } 
#endif
# 172 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_30_intrinsics.h"
__attribute((deprecated("__shfl_up() is deprecated in favor of __shfl_up_sync() and may be removed in a future release (Use -Wno-deprecated-declarations " "to suppress this warning)."))) __attribute__((unused)) static inline unsigned long long __shfl_up(unsigned long long var, unsigned delta, int width = 32) {int volatile ___ = 1;(void)var;(void)delta;(void)width;::exit(___);}
#if 0
# 172
{ } 
#endif
# 173 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_30_intrinsics.h"
__attribute((deprecated("__shfl_down() is deprecated in favor of __shfl_down_sync() and may be removed in a future release (Use -Wno-deprecated-declarati" "ons to suppress this warning)."))) __attribute__((unused)) static inline long long __shfl_down(long long var, unsigned delta, int width = 32) {int volatile ___ = 1;(void)var;(void)delta;(void)width;::exit(___);}
#if 0
# 173
{ } 
#endif
# 174 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_30_intrinsics.h"
__attribute((deprecated("__shfl_down() is deprecated in favor of __shfl_down_sync() and may be removed in a future release (Use -Wno-deprecated-declarati" "ons to suppress this warning)."))) __attribute__((unused)) static inline unsigned long long __shfl_down(unsigned long long var, unsigned delta, int width = 32) {int volatile ___ = 1;(void)var;(void)delta;(void)width;::exit(___);}
#if 0
# 174
{ } 
#endif
# 175 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_30_intrinsics.h"
__attribute((deprecated("__shfl_xor() is deprecated in favor of __shfl_xor_sync() and may be removed in a future release (Use -Wno-deprecated-declaration" "s to suppress this warning)."))) __attribute__((unused)) static inline long long __shfl_xor(long long var, int laneMask, int width = 32) {int volatile ___ = 1;(void)var;(void)laneMask;(void)width;::exit(___);}
#if 0
# 175
{ } 
#endif
# 176 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_30_intrinsics.h"
__attribute((deprecated("__shfl_xor() is deprecated in favor of __shfl_xor_sync() and may be removed in a future release (Use -Wno-deprecated-declaration" "s to suppress this warning)."))) __attribute__((unused)) static inline unsigned long long __shfl_xor(unsigned long long var, int laneMask, int width = 32) {int volatile ___ = 1;(void)var;(void)laneMask;(void)width;::exit(___);}
#if 0
# 176
{ } 
#endif
# 177 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_30_intrinsics.h"
__attribute((deprecated("__shfl() is deprecated in favor of __shfl_sync() and may be removed in a future release (Use -Wno-deprecated-declarations to sup" "press this warning)."))) __attribute__((unused)) static inline double __shfl(double var, int srcLane, int width = 32) {int volatile ___ = 1;(void)var;(void)srcLane;(void)width;::exit(___);}
#if 0
# 177
{ } 
#endif
# 178 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_30_intrinsics.h"
__attribute((deprecated("__shfl_up() is deprecated in favor of __shfl_up_sync() and may be removed in a future release (Use -Wno-deprecated-declarations " "to suppress this warning)."))) __attribute__((unused)) static inline double __shfl_up(double var, unsigned delta, int width = 32) {int volatile ___ = 1;(void)var;(void)delta;(void)width;::exit(___);}
#if 0
# 178
{ } 
#endif
# 179 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_30_intrinsics.h"
__attribute((deprecated("__shfl_down() is deprecated in favor of __shfl_down_sync() and may be removed in a future release (Use -Wno-deprecated-declarati" "ons to suppress this warning)."))) __attribute__((unused)) static inline double __shfl_down(double var, unsigned delta, int width = 32) {int volatile ___ = 1;(void)var;(void)delta;(void)width;::exit(___);}
#if 0
# 179
{ } 
#endif
# 180 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_30_intrinsics.h"
__attribute((deprecated("__shfl_xor() is deprecated in favor of __shfl_xor_sync() and may be removed in a future release (Use -Wno-deprecated-declaration" "s to suppress this warning)."))) __attribute__((unused)) static inline double __shfl_xor(double var, int laneMask, int width = 32) {int volatile ___ = 1;(void)var;(void)laneMask;(void)width;::exit(___);}
#if 0
# 180
{ } 
#endif
# 183 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_30_intrinsics.h"
__attribute__((unused)) static inline long long __shfl_sync(unsigned mask, long long var, int srcLane, int width = 32) {int volatile ___ = 1;(void)mask;(void)var;(void)srcLane;(void)width;::exit(___);}
#if 0
# 183
{ } 
#endif
# 184 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_30_intrinsics.h"
__attribute__((unused)) static inline unsigned long long __shfl_sync(unsigned mask, unsigned long long var, int srcLane, int width = 32) {int volatile ___ = 1;(void)mask;(void)var;(void)srcLane;(void)width;::exit(___);}
#if 0
# 184
{ } 
#endif
# 185 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_30_intrinsics.h"
__attribute__((unused)) static inline long long __shfl_up_sync(unsigned mask, long long var, unsigned delta, int width = 32) {int volatile ___ = 1;(void)mask;(void)var;(void)delta;(void)width;::exit(___);}
#if 0
# 185
{ } 
#endif
# 186 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_30_intrinsics.h"
__attribute__((unused)) static inline unsigned long long __shfl_up_sync(unsigned mask, unsigned long long var, unsigned delta, int width = 32) {int volatile ___ = 1;(void)mask;(void)var;(void)delta;(void)width;::exit(___);}
#if 0
# 186
{ } 
#endif
# 187 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_30_intrinsics.h"
__attribute__((unused)) static inline long long __shfl_down_sync(unsigned mask, long long var, unsigned delta, int width = 32) {int volatile ___ = 1;(void)mask;(void)var;(void)delta;(void)width;::exit(___);}
#if 0
# 187
{ } 
#endif
# 188 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_30_intrinsics.h"
__attribute__((unused)) static inline unsigned long long __shfl_down_sync(unsigned mask, unsigned long long var, unsigned delta, int width = 32) {int volatile ___ = 1;(void)mask;(void)var;(void)delta;(void)width;::exit(___);}
#if 0
# 188
{ } 
#endif
# 189 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_30_intrinsics.h"
__attribute__((unused)) static inline long long __shfl_xor_sync(unsigned mask, long long var, int laneMask, int width = 32) {int volatile ___ = 1;(void)mask;(void)var;(void)laneMask;(void)width;::exit(___);}
#if 0
# 189
{ } 
#endif
# 190 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_30_intrinsics.h"
__attribute__((unused)) static inline unsigned long long __shfl_xor_sync(unsigned mask, unsigned long long var, int laneMask, int width = 32) {int volatile ___ = 1;(void)mask;(void)var;(void)laneMask;(void)width;::exit(___);}
#if 0
# 190
{ } 
#endif
# 191 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_30_intrinsics.h"
__attribute__((unused)) static inline double __shfl_sync(unsigned mask, double var, int srcLane, int width = 32) {int volatile ___ = 1;(void)mask;(void)var;(void)srcLane;(void)width;::exit(___);}
#if 0
# 191
{ } 
#endif
# 192 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_30_intrinsics.h"
__attribute__((unused)) static inline double __shfl_up_sync(unsigned mask, double var, unsigned delta, int width = 32) {int volatile ___ = 1;(void)mask;(void)var;(void)delta;(void)width;::exit(___);}
#if 0
# 192
{ } 
#endif
# 193 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_30_intrinsics.h"
__attribute__((unused)) static inline double __shfl_down_sync(unsigned mask, double var, unsigned delta, int width = 32) {int volatile ___ = 1;(void)mask;(void)var;(void)delta;(void)width;::exit(___);}
#if 0
# 193
{ } 
#endif
# 194 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_30_intrinsics.h"
__attribute__((unused)) static inline double __shfl_xor_sync(unsigned mask, double var, int laneMask, int width = 32) {int volatile ___ = 1;(void)mask;(void)var;(void)laneMask;(void)width;::exit(___);}
#if 0
# 194
{ } 
#endif
# 198 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_30_intrinsics.h"
__attribute((deprecated("__shfl() is deprecated in favor of __shfl_sync() and may be removed in a future release (Use -Wno-deprecated-declarations to sup" "press this warning)."))) __attribute__((unused)) static inline long __shfl(long var, int srcLane, int width = 32) {int volatile ___ = 1;(void)var;(void)srcLane;(void)width;::exit(___);}
#if 0
# 198
{ } 
#endif
# 199 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_30_intrinsics.h"
__attribute((deprecated("__shfl() is deprecated in favor of __shfl_sync() and may be removed in a future release (Use -Wno-deprecated-declarations to sup" "press this warning)."))) __attribute__((unused)) static inline unsigned long __shfl(unsigned long var, int srcLane, int width = 32) {int volatile ___ = 1;(void)var;(void)srcLane;(void)width;::exit(___);}
#if 0
# 199
{ } 
#endif
# 200 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_30_intrinsics.h"
__attribute((deprecated("__shfl_up() is deprecated in favor of __shfl_up_sync() and may be removed in a future release (Use -Wno-deprecated-declarations " "to suppress this warning)."))) __attribute__((unused)) static inline long __shfl_up(long var, unsigned delta, int width = 32) {int volatile ___ = 1;(void)var;(void)delta;(void)width;::exit(___);}
#if 0
# 200
{ } 
#endif
# 201 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_30_intrinsics.h"
__attribute((deprecated("__shfl_up() is deprecated in favor of __shfl_up_sync() and may be removed in a future release (Use -Wno-deprecated-declarations " "to suppress this warning)."))) __attribute__((unused)) static inline unsigned long __shfl_up(unsigned long var, unsigned delta, int width = 32) {int volatile ___ = 1;(void)var;(void)delta;(void)width;::exit(___);}
#if 0
# 201
{ } 
#endif
# 202 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_30_intrinsics.h"
__attribute((deprecated("__shfl_down() is deprecated in favor of __shfl_down_sync() and may be removed in a future release (Use -Wno-deprecated-declarati" "ons to suppress this warning)."))) __attribute__((unused)) static inline long __shfl_down(long var, unsigned delta, int width = 32) {int volatile ___ = 1;(void)var;(void)delta;(void)width;::exit(___);}
#if 0
# 202
{ } 
#endif
# 203 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_30_intrinsics.h"
__attribute((deprecated("__shfl_down() is deprecated in favor of __shfl_down_sync() and may be removed in a future release (Use -Wno-deprecated-declarati" "ons to suppress this warning)."))) __attribute__((unused)) static inline unsigned long __shfl_down(unsigned long var, unsigned delta, int width = 32) {int volatile ___ = 1;(void)var;(void)delta;(void)width;::exit(___);}
#if 0
# 203
{ } 
#endif
# 204 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_30_intrinsics.h"
__attribute((deprecated("__shfl_xor() is deprecated in favor of __shfl_xor_sync() and may be removed in a future release (Use -Wno-deprecated-declaration" "s to suppress this warning)."))) __attribute__((unused)) static inline long __shfl_xor(long var, int laneMask, int width = 32) {int volatile ___ = 1;(void)var;(void)laneMask;(void)width;::exit(___);}
#if 0
# 204
{ } 
#endif
# 205 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_30_intrinsics.h"
__attribute((deprecated("__shfl_xor() is deprecated in favor of __shfl_xor_sync() and may be removed in a future release (Use -Wno-deprecated-declaration" "s to suppress this warning)."))) __attribute__((unused)) static inline unsigned long __shfl_xor(unsigned long var, int laneMask, int width = 32) {int volatile ___ = 1;(void)var;(void)laneMask;(void)width;::exit(___);}
#if 0
# 205
{ } 
#endif
# 208 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_30_intrinsics.h"
__attribute__((unused)) static inline long __shfl_sync(unsigned mask, long var, int srcLane, int width = 32) {int volatile ___ = 1;(void)mask;(void)var;(void)srcLane;(void)width;::exit(___);}
#if 0
# 208
{ } 
#endif
# 209 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_30_intrinsics.h"
__attribute__((unused)) static inline unsigned long __shfl_sync(unsigned mask, unsigned long var, int srcLane, int width = 32) {int volatile ___ = 1;(void)mask;(void)var;(void)srcLane;(void)width;::exit(___);}
#if 0
# 209
{ } 
#endif
# 210 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_30_intrinsics.h"
__attribute__((unused)) static inline long __shfl_up_sync(unsigned mask, long var, unsigned delta, int width = 32) {int volatile ___ = 1;(void)mask;(void)var;(void)delta;(void)width;::exit(___);}
#if 0
# 210
{ } 
#endif
# 211 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_30_intrinsics.h"
__attribute__((unused)) static inline unsigned long __shfl_up_sync(unsigned mask, unsigned long var, unsigned delta, int width = 32) {int volatile ___ = 1;(void)mask;(void)var;(void)delta;(void)width;::exit(___);}
#if 0
# 211
{ } 
#endif
# 212 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_30_intrinsics.h"
__attribute__((unused)) static inline long __shfl_down_sync(unsigned mask, long var, unsigned delta, int width = 32) {int volatile ___ = 1;(void)mask;(void)var;(void)delta;(void)width;::exit(___);}
#if 0
# 212
{ } 
#endif
# 213 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_30_intrinsics.h"
__attribute__((unused)) static inline unsigned long __shfl_down_sync(unsigned mask, unsigned long var, unsigned delta, int width = 32) {int volatile ___ = 1;(void)mask;(void)var;(void)delta;(void)width;::exit(___);}
#if 0
# 213
{ } 
#endif
# 214 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_30_intrinsics.h"
__attribute__((unused)) static inline long __shfl_xor_sync(unsigned mask, long var, int laneMask, int width = 32) {int volatile ___ = 1;(void)mask;(void)var;(void)laneMask;(void)width;::exit(___);}
#if 0
# 214
{ } 
#endif
# 215 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_30_intrinsics.h"
__attribute__((unused)) static inline unsigned long __shfl_xor_sync(unsigned mask, unsigned long var, int laneMask, int width = 32) {int volatile ___ = 1;(void)mask;(void)var;(void)laneMask;(void)width;::exit(___);}
#if 0
# 215
{ } 
#endif
# 87 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline long __ldg(const long *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 87
{ } 
#endif
# 88 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline unsigned long __ldg(const unsigned long *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 88
{ } 
#endif
# 90 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline char __ldg(const char *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 90
{ } 
#endif
# 91 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline signed char __ldg(const signed char *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 91
{ } 
#endif
# 92 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline short __ldg(const short *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 92
{ } 
#endif
# 93 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline int __ldg(const int *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 93
{ } 
#endif
# 94 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline long long __ldg(const long long *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 94
{ } 
#endif
# 95 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline char2 __ldg(const char2 *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 95
{ } 
#endif
# 96 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline char4 __ldg(const char4 *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 96
{ } 
#endif
# 97 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline short2 __ldg(const short2 *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 97
{ } 
#endif
# 98 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline short4 __ldg(const short4 *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 98
{ } 
#endif
# 99 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline int2 __ldg(const int2 *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 99
{ } 
#endif
# 100 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline int4 __ldg(const int4 *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 100
{ } 
#endif
# 101 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline longlong2 __ldg(const longlong2 *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 101
{ } 
#endif
# 103 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline unsigned char __ldg(const unsigned char *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 103
{ } 
#endif
# 104 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline unsigned short __ldg(const unsigned short *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 104
{ } 
#endif
# 105 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline unsigned __ldg(const unsigned *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 105
{ } 
#endif
# 106 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline unsigned long long __ldg(const unsigned long long *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 106
{ } 
#endif
# 107 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline uchar2 __ldg(const uchar2 *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 107
{ } 
#endif
# 108 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline uchar4 __ldg(const uchar4 *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 108
{ } 
#endif
# 109 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline ushort2 __ldg(const ushort2 *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 109
{ } 
#endif
# 110 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline ushort4 __ldg(const ushort4 *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 110
{ } 
#endif
# 111 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline uint2 __ldg(const uint2 *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 111
{ } 
#endif
# 112 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline uint4 __ldg(const uint4 *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 112
{ } 
#endif
# 113 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline ulonglong2 __ldg(const ulonglong2 *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 113
{ } 
#endif
# 115 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline float __ldg(const float *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 115
{ } 
#endif
# 116 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline double __ldg(const double *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 116
{ } 
#endif
# 117 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline float2 __ldg(const float2 *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 117
{ } 
#endif
# 118 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline float4 __ldg(const float4 *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 118
{ } 
#endif
# 119 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline double2 __ldg(const double2 *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 119
{ } 
#endif
# 123 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline long __ldcg(const long *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 123
{ } 
#endif
# 124 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline unsigned long __ldcg(const unsigned long *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 124
{ } 
#endif
# 126 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline char __ldcg(const char *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 126
{ } 
#endif
# 127 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline signed char __ldcg(const signed char *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 127
{ } 
#endif
# 128 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline short __ldcg(const short *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 128
{ } 
#endif
# 129 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline int __ldcg(const int *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 129
{ } 
#endif
# 130 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline long long __ldcg(const long long *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 130
{ } 
#endif
# 131 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline char2 __ldcg(const char2 *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 131
{ } 
#endif
# 132 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline char4 __ldcg(const char4 *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 132
{ } 
#endif
# 133 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline short2 __ldcg(const short2 *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 133
{ } 
#endif
# 134 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline short4 __ldcg(const short4 *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 134
{ } 
#endif
# 135 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline int2 __ldcg(const int2 *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 135
{ } 
#endif
# 136 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline int4 __ldcg(const int4 *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 136
{ } 
#endif
# 137 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline longlong2 __ldcg(const longlong2 *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 137
{ } 
#endif
# 139 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline unsigned char __ldcg(const unsigned char *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 139
{ } 
#endif
# 140 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline unsigned short __ldcg(const unsigned short *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 140
{ } 
#endif
# 141 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline unsigned __ldcg(const unsigned *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 141
{ } 
#endif
# 142 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline unsigned long long __ldcg(const unsigned long long *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 142
{ } 
#endif
# 143 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline uchar2 __ldcg(const uchar2 *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 143
{ } 
#endif
# 144 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline uchar4 __ldcg(const uchar4 *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 144
{ } 
#endif
# 145 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline ushort2 __ldcg(const ushort2 *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 145
{ } 
#endif
# 146 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline ushort4 __ldcg(const ushort4 *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 146
{ } 
#endif
# 147 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline uint2 __ldcg(const uint2 *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 147
{ } 
#endif
# 148 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline uint4 __ldcg(const uint4 *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 148
{ } 
#endif
# 149 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline ulonglong2 __ldcg(const ulonglong2 *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 149
{ } 
#endif
# 151 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline float __ldcg(const float *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 151
{ } 
#endif
# 152 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline double __ldcg(const double *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 152
{ } 
#endif
# 153 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline float2 __ldcg(const float2 *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 153
{ } 
#endif
# 154 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline float4 __ldcg(const float4 *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 154
{ } 
#endif
# 155 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline double2 __ldcg(const double2 *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 155
{ } 
#endif
# 159 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline long __ldca(const long *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 159
{ } 
#endif
# 160 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline unsigned long __ldca(const unsigned long *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 160
{ } 
#endif
# 162 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline char __ldca(const char *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 162
{ } 
#endif
# 163 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline signed char __ldca(const signed char *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 163
{ } 
#endif
# 164 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline short __ldca(const short *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 164
{ } 
#endif
# 165 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline int __ldca(const int *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 165
{ } 
#endif
# 166 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline long long __ldca(const long long *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 166
{ } 
#endif
# 167 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline char2 __ldca(const char2 *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 167
{ } 
#endif
# 168 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline char4 __ldca(const char4 *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 168
{ } 
#endif
# 169 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline short2 __ldca(const short2 *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 169
{ } 
#endif
# 170 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline short4 __ldca(const short4 *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 170
{ } 
#endif
# 171 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline int2 __ldca(const int2 *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 171
{ } 
#endif
# 172 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline int4 __ldca(const int4 *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 172
{ } 
#endif
# 173 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline longlong2 __ldca(const longlong2 *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 173
{ } 
#endif
# 175 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline unsigned char __ldca(const unsigned char *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 175
{ } 
#endif
# 176 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline unsigned short __ldca(const unsigned short *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 176
{ } 
#endif
# 177 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline unsigned __ldca(const unsigned *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 177
{ } 
#endif
# 178 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline unsigned long long __ldca(const unsigned long long *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 178
{ } 
#endif
# 179 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline uchar2 __ldca(const uchar2 *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 179
{ } 
#endif
# 180 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline uchar4 __ldca(const uchar4 *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 180
{ } 
#endif
# 181 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline ushort2 __ldca(const ushort2 *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 181
{ } 
#endif
# 182 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline ushort4 __ldca(const ushort4 *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 182
{ } 
#endif
# 183 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline uint2 __ldca(const uint2 *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 183
{ } 
#endif
# 184 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline uint4 __ldca(const uint4 *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 184
{ } 
#endif
# 185 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline ulonglong2 __ldca(const ulonglong2 *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 185
{ } 
#endif
# 187 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline float __ldca(const float *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 187
{ } 
#endif
# 188 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline double __ldca(const double *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 188
{ } 
#endif
# 189 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline float2 __ldca(const float2 *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 189
{ } 
#endif
# 190 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline float4 __ldca(const float4 *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 190
{ } 
#endif
# 191 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline double2 __ldca(const double2 *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 191
{ } 
#endif
# 195 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline long __ldcs(const long *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 195
{ } 
#endif
# 196 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline unsigned long __ldcs(const unsigned long *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 196
{ } 
#endif
# 198 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline char __ldcs(const char *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 198
{ } 
#endif
# 199 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline signed char __ldcs(const signed char *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 199
{ } 
#endif
# 200 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline short __ldcs(const short *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 200
{ } 
#endif
# 201 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline int __ldcs(const int *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 201
{ } 
#endif
# 202 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline long long __ldcs(const long long *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 202
{ } 
#endif
# 203 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline char2 __ldcs(const char2 *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 203
{ } 
#endif
# 204 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline char4 __ldcs(const char4 *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 204
{ } 
#endif
# 205 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline short2 __ldcs(const short2 *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 205
{ } 
#endif
# 206 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline short4 __ldcs(const short4 *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 206
{ } 
#endif
# 207 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline int2 __ldcs(const int2 *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 207
{ } 
#endif
# 208 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline int4 __ldcs(const int4 *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 208
{ } 
#endif
# 209 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline longlong2 __ldcs(const longlong2 *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 209
{ } 
#endif
# 211 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline unsigned char __ldcs(const unsigned char *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 211
{ } 
#endif
# 212 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline unsigned short __ldcs(const unsigned short *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 212
{ } 
#endif
# 213 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline unsigned __ldcs(const unsigned *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 213
{ } 
#endif
# 214 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline unsigned long long __ldcs(const unsigned long long *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 214
{ } 
#endif
# 215 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline uchar2 __ldcs(const uchar2 *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 215
{ } 
#endif
# 216 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline uchar4 __ldcs(const uchar4 *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 216
{ } 
#endif
# 217 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline ushort2 __ldcs(const ushort2 *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 217
{ } 
#endif
# 218 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline ushort4 __ldcs(const ushort4 *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 218
{ } 
#endif
# 219 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline uint2 __ldcs(const uint2 *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 219
{ } 
#endif
# 220 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline uint4 __ldcs(const uint4 *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 220
{ } 
#endif
# 221 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline ulonglong2 __ldcs(const ulonglong2 *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 221
{ } 
#endif
# 223 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline float __ldcs(const float *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 223
{ } 
#endif
# 224 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline double __ldcs(const double *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 224
{ } 
#endif
# 225 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline float2 __ldcs(const float2 *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 225
{ } 
#endif
# 226 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline float4 __ldcs(const float4 *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 226
{ } 
#endif
# 227 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline double2 __ldcs(const double2 *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 227
{ } 
#endif
# 231 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline long __ldlu(const long *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 231
{ } 
#endif
# 232 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline unsigned long __ldlu(const unsigned long *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 232
{ } 
#endif
# 234 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline char __ldlu(const char *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 234
{ } 
#endif
# 235 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline signed char __ldlu(const signed char *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 235
{ } 
#endif
# 236 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline short __ldlu(const short *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 236
{ } 
#endif
# 237 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline int __ldlu(const int *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 237
{ } 
#endif
# 238 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline long long __ldlu(const long long *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 238
{ } 
#endif
# 239 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline char2 __ldlu(const char2 *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 239
{ } 
#endif
# 240 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline char4 __ldlu(const char4 *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 240
{ } 
#endif
# 241 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline short2 __ldlu(const short2 *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 241
{ } 
#endif
# 242 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline short4 __ldlu(const short4 *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 242
{ } 
#endif
# 243 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline int2 __ldlu(const int2 *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 243
{ } 
#endif
# 244 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline int4 __ldlu(const int4 *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 244
{ } 
#endif
# 245 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline longlong2 __ldlu(const longlong2 *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 245
{ } 
#endif
# 247 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline unsigned char __ldlu(const unsigned char *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 247
{ } 
#endif
# 248 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline unsigned short __ldlu(const unsigned short *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 248
{ } 
#endif
# 249 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline unsigned __ldlu(const unsigned *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 249
{ } 
#endif
# 250 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline unsigned long long __ldlu(const unsigned long long *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 250
{ } 
#endif
# 251 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline uchar2 __ldlu(const uchar2 *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 251
{ } 
#endif
# 252 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline uchar4 __ldlu(const uchar4 *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 252
{ } 
#endif
# 253 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline ushort2 __ldlu(const ushort2 *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 253
{ } 
#endif
# 254 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline ushort4 __ldlu(const ushort4 *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 254
{ } 
#endif
# 255 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline uint2 __ldlu(const uint2 *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 255
{ } 
#endif
# 256 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline uint4 __ldlu(const uint4 *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 256
{ } 
#endif
# 257 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline ulonglong2 __ldlu(const ulonglong2 *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 257
{ } 
#endif
# 259 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline float __ldlu(const float *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 259
{ } 
#endif
# 260 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline double __ldlu(const double *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 260
{ } 
#endif
# 261 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline float2 __ldlu(const float2 *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 261
{ } 
#endif
# 262 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline float4 __ldlu(const float4 *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 262
{ } 
#endif
# 263 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline double2 __ldlu(const double2 *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 263
{ } 
#endif
# 267 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline long __ldcv(const long *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 267
{ } 
#endif
# 268 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline unsigned long __ldcv(const unsigned long *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 268
{ } 
#endif
# 270 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline char __ldcv(const char *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 270
{ } 
#endif
# 271 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline signed char __ldcv(const signed char *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 271
{ } 
#endif
# 272 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline short __ldcv(const short *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 272
{ } 
#endif
# 273 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline int __ldcv(const int *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 273
{ } 
#endif
# 274 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline long long __ldcv(const long long *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 274
{ } 
#endif
# 275 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline char2 __ldcv(const char2 *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 275
{ } 
#endif
# 276 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline char4 __ldcv(const char4 *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 276
{ } 
#endif
# 277 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline short2 __ldcv(const short2 *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 277
{ } 
#endif
# 278 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline short4 __ldcv(const short4 *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 278
{ } 
#endif
# 279 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline int2 __ldcv(const int2 *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 279
{ } 
#endif
# 280 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline int4 __ldcv(const int4 *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 280
{ } 
#endif
# 281 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline longlong2 __ldcv(const longlong2 *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 281
{ } 
#endif
# 283 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline unsigned char __ldcv(const unsigned char *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 283
{ } 
#endif
# 284 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline unsigned short __ldcv(const unsigned short *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 284
{ } 
#endif
# 285 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline unsigned __ldcv(const unsigned *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 285
{ } 
#endif
# 286 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline unsigned long long __ldcv(const unsigned long long *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 286
{ } 
#endif
# 287 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline uchar2 __ldcv(const uchar2 *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 287
{ } 
#endif
# 288 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline uchar4 __ldcv(const uchar4 *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 288
{ } 
#endif
# 289 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline ushort2 __ldcv(const ushort2 *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 289
{ } 
#endif
# 290 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline ushort4 __ldcv(const ushort4 *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 290
{ } 
#endif
# 291 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline uint2 __ldcv(const uint2 *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 291
{ } 
#endif
# 292 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline uint4 __ldcv(const uint4 *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 292
{ } 
#endif
# 293 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline ulonglong2 __ldcv(const ulonglong2 *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 293
{ } 
#endif
# 295 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline float __ldcv(const float *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 295
{ } 
#endif
# 296 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline double __ldcv(const double *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 296
{ } 
#endif
# 297 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline float2 __ldcv(const float2 *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 297
{ } 
#endif
# 298 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline float4 __ldcv(const float4 *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 298
{ } 
#endif
# 299 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline double2 __ldcv(const double2 *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 299
{ } 
#endif
# 303 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline void __stwb(long *ptr, long value) {int volatile ___ = 1;(void)ptr;(void)value;::exit(___);}
#if 0
# 303
{ } 
#endif
# 304 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline void __stwb(unsigned long *ptr, unsigned long value) {int volatile ___ = 1;(void)ptr;(void)value;::exit(___);}
#if 0
# 304
{ } 
#endif
# 306 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline void __stwb(char *ptr, char value) {int volatile ___ = 1;(void)ptr;(void)value;::exit(___);}
#if 0
# 306
{ } 
#endif
# 307 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline void __stwb(signed char *ptr, signed char value) {int volatile ___ = 1;(void)ptr;(void)value;::exit(___);}
#if 0
# 307
{ } 
#endif
# 308 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline void __stwb(short *ptr, short value) {int volatile ___ = 1;(void)ptr;(void)value;::exit(___);}
#if 0
# 308
{ } 
#endif
# 309 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline void __stwb(int *ptr, int value) {int volatile ___ = 1;(void)ptr;(void)value;::exit(___);}
#if 0
# 309
{ } 
#endif
# 310 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline void __stwb(long long *ptr, long long value) {int volatile ___ = 1;(void)ptr;(void)value;::exit(___);}
#if 0
# 310
{ } 
#endif
# 311 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline void __stwb(char2 *ptr, char2 value) {int volatile ___ = 1;(void)ptr;(void)value;::exit(___);}
#if 0
# 311
{ } 
#endif
# 312 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline void __stwb(char4 *ptr, char4 value) {int volatile ___ = 1;(void)ptr;(void)value;::exit(___);}
#if 0
# 312
{ } 
#endif
# 313 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline void __stwb(short2 *ptr, short2 value) {int volatile ___ = 1;(void)ptr;(void)value;::exit(___);}
#if 0
# 313
{ } 
#endif
# 314 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline void __stwb(short4 *ptr, short4 value) {int volatile ___ = 1;(void)ptr;(void)value;::exit(___);}
#if 0
# 314
{ } 
#endif
# 315 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline void __stwb(int2 *ptr, int2 value) {int volatile ___ = 1;(void)ptr;(void)value;::exit(___);}
#if 0
# 315
{ } 
#endif
# 316 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline void __stwb(int4 *ptr, int4 value) {int volatile ___ = 1;(void)ptr;(void)value;::exit(___);}
#if 0
# 316
{ } 
#endif
# 317 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline void __stwb(longlong2 *ptr, longlong2 value) {int volatile ___ = 1;(void)ptr;(void)value;::exit(___);}
#if 0
# 317
{ } 
#endif
# 319 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline void __stwb(unsigned char *ptr, unsigned char value) {int volatile ___ = 1;(void)ptr;(void)value;::exit(___);}
#if 0
# 319
{ } 
#endif
# 320 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline void __stwb(unsigned short *ptr, unsigned short value) {int volatile ___ = 1;(void)ptr;(void)value;::exit(___);}
#if 0
# 320
{ } 
#endif
# 321 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline void __stwb(unsigned *ptr, unsigned value) {int volatile ___ = 1;(void)ptr;(void)value;::exit(___);}
#if 0
# 321
{ } 
#endif
# 322 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline void __stwb(unsigned long long *ptr, unsigned long long value) {int volatile ___ = 1;(void)ptr;(void)value;::exit(___);}
#if 0
# 322
{ } 
#endif
# 323 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline void __stwb(uchar2 *ptr, uchar2 value) {int volatile ___ = 1;(void)ptr;(void)value;::exit(___);}
#if 0
# 323
{ } 
#endif
# 324 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline void __stwb(uchar4 *ptr, uchar4 value) {int volatile ___ = 1;(void)ptr;(void)value;::exit(___);}
#if 0
# 324
{ } 
#endif
# 325 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline void __stwb(ushort2 *ptr, ushort2 value) {int volatile ___ = 1;(void)ptr;(void)value;::exit(___);}
#if 0
# 325
{ } 
#endif
# 326 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline void __stwb(ushort4 *ptr, ushort4 value) {int volatile ___ = 1;(void)ptr;(void)value;::exit(___);}
#if 0
# 326
{ } 
#endif
# 327 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline void __stwb(uint2 *ptr, uint2 value) {int volatile ___ = 1;(void)ptr;(void)value;::exit(___);}
#if 0
# 327
{ } 
#endif
# 328 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline void __stwb(uint4 *ptr, uint4 value) {int volatile ___ = 1;(void)ptr;(void)value;::exit(___);}
#if 0
# 328
{ } 
#endif
# 329 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline void __stwb(ulonglong2 *ptr, ulonglong2 value) {int volatile ___ = 1;(void)ptr;(void)value;::exit(___);}
#if 0
# 329
{ } 
#endif
# 331 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline void __stwb(float *ptr, float value) {int volatile ___ = 1;(void)ptr;(void)value;::exit(___);}
#if 0
# 331
{ } 
#endif
# 332 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline void __stwb(double *ptr, double value) {int volatile ___ = 1;(void)ptr;(void)value;::exit(___);}
#if 0
# 332
{ } 
#endif
# 333 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline void __stwb(float2 *ptr, float2 value) {int volatile ___ = 1;(void)ptr;(void)value;::exit(___);}
#if 0
# 333
{ } 
#endif
# 334 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline void __stwb(float4 *ptr, float4 value) {int volatile ___ = 1;(void)ptr;(void)value;::exit(___);}
#if 0
# 334
{ } 
#endif
# 335 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline void __stwb(double2 *ptr, double2 value) {int volatile ___ = 1;(void)ptr;(void)value;::exit(___);}
#if 0
# 335
{ } 
#endif
# 339 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline void __stcg(long *ptr, long value) {int volatile ___ = 1;(void)ptr;(void)value;::exit(___);}
#if 0
# 339
{ } 
#endif
# 340 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline void __stcg(unsigned long *ptr, unsigned long value) {int volatile ___ = 1;(void)ptr;(void)value;::exit(___);}
#if 0
# 340
{ } 
#endif
# 342 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline void __stcg(char *ptr, char value) {int volatile ___ = 1;(void)ptr;(void)value;::exit(___);}
#if 0
# 342
{ } 
#endif
# 343 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline void __stcg(signed char *ptr, signed char value) {int volatile ___ = 1;(void)ptr;(void)value;::exit(___);}
#if 0
# 343
{ } 
#endif
# 344 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline void __stcg(short *ptr, short value) {int volatile ___ = 1;(void)ptr;(void)value;::exit(___);}
#if 0
# 344
{ } 
#endif
# 345 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline void __stcg(int *ptr, int value) {int volatile ___ = 1;(void)ptr;(void)value;::exit(___);}
#if 0
# 345
{ } 
#endif
# 346 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline void __stcg(long long *ptr, long long value) {int volatile ___ = 1;(void)ptr;(void)value;::exit(___);}
#if 0
# 346
{ } 
#endif
# 347 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline void __stcg(char2 *ptr, char2 value) {int volatile ___ = 1;(void)ptr;(void)value;::exit(___);}
#if 0
# 347
{ } 
#endif
# 348 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline void __stcg(char4 *ptr, char4 value) {int volatile ___ = 1;(void)ptr;(void)value;::exit(___);}
#if 0
# 348
{ } 
#endif
# 349 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline void __stcg(short2 *ptr, short2 value) {int volatile ___ = 1;(void)ptr;(void)value;::exit(___);}
#if 0
# 349
{ } 
#endif
# 350 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline void __stcg(short4 *ptr, short4 value) {int volatile ___ = 1;(void)ptr;(void)value;::exit(___);}
#if 0
# 350
{ } 
#endif
# 351 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline void __stcg(int2 *ptr, int2 value) {int volatile ___ = 1;(void)ptr;(void)value;::exit(___);}
#if 0
# 351
{ } 
#endif
# 352 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline void __stcg(int4 *ptr, int4 value) {int volatile ___ = 1;(void)ptr;(void)value;::exit(___);}
#if 0
# 352
{ } 
#endif
# 353 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline void __stcg(longlong2 *ptr, longlong2 value) {int volatile ___ = 1;(void)ptr;(void)value;::exit(___);}
#if 0
# 353
{ } 
#endif
# 355 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline void __stcg(unsigned char *ptr, unsigned char value) {int volatile ___ = 1;(void)ptr;(void)value;::exit(___);}
#if 0
# 355
{ } 
#endif
# 356 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline void __stcg(unsigned short *ptr, unsigned short value) {int volatile ___ = 1;(void)ptr;(void)value;::exit(___);}
#if 0
# 356
{ } 
#endif
# 357 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline void __stcg(unsigned *ptr, unsigned value) {int volatile ___ = 1;(void)ptr;(void)value;::exit(___);}
#if 0
# 357
{ } 
#endif
# 358 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline void __stcg(unsigned long long *ptr, unsigned long long value) {int volatile ___ = 1;(void)ptr;(void)value;::exit(___);}
#if 0
# 358
{ } 
#endif
# 359 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline void __stcg(uchar2 *ptr, uchar2 value) {int volatile ___ = 1;(void)ptr;(void)value;::exit(___);}
#if 0
# 359
{ } 
#endif
# 360 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline void __stcg(uchar4 *ptr, uchar4 value) {int volatile ___ = 1;(void)ptr;(void)value;::exit(___);}
#if 0
# 360
{ } 
#endif
# 361 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline void __stcg(ushort2 *ptr, ushort2 value) {int volatile ___ = 1;(void)ptr;(void)value;::exit(___);}
#if 0
# 361
{ } 
#endif
# 362 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline void __stcg(ushort4 *ptr, ushort4 value) {int volatile ___ = 1;(void)ptr;(void)value;::exit(___);}
#if 0
# 362
{ } 
#endif
# 363 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline void __stcg(uint2 *ptr, uint2 value) {int volatile ___ = 1;(void)ptr;(void)value;::exit(___);}
#if 0
# 363
{ } 
#endif
# 364 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline void __stcg(uint4 *ptr, uint4 value) {int volatile ___ = 1;(void)ptr;(void)value;::exit(___);}
#if 0
# 364
{ } 
#endif
# 365 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline void __stcg(ulonglong2 *ptr, ulonglong2 value) {int volatile ___ = 1;(void)ptr;(void)value;::exit(___);}
#if 0
# 365
{ } 
#endif
# 367 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline void __stcg(float *ptr, float value) {int volatile ___ = 1;(void)ptr;(void)value;::exit(___);}
#if 0
# 367
{ } 
#endif
# 368 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline void __stcg(double *ptr, double value) {int volatile ___ = 1;(void)ptr;(void)value;::exit(___);}
#if 0
# 368
{ } 
#endif
# 369 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline void __stcg(float2 *ptr, float2 value) {int volatile ___ = 1;(void)ptr;(void)value;::exit(___);}
#if 0
# 369
{ } 
#endif
# 370 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline void __stcg(float4 *ptr, float4 value) {int volatile ___ = 1;(void)ptr;(void)value;::exit(___);}
#if 0
# 370
{ } 
#endif
# 371 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline void __stcg(double2 *ptr, double2 value) {int volatile ___ = 1;(void)ptr;(void)value;::exit(___);}
#if 0
# 371
{ } 
#endif
# 375 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline void __stcs(long *ptr, long value) {int volatile ___ = 1;(void)ptr;(void)value;::exit(___);}
#if 0
# 375
{ } 
#endif
# 376 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline void __stcs(unsigned long *ptr, unsigned long value) {int volatile ___ = 1;(void)ptr;(void)value;::exit(___);}
#if 0
# 376
{ } 
#endif
# 378 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline void __stcs(char *ptr, char value) {int volatile ___ = 1;(void)ptr;(void)value;::exit(___);}
#if 0
# 378
{ } 
#endif
# 379 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline void __stcs(signed char *ptr, signed char value) {int volatile ___ = 1;(void)ptr;(void)value;::exit(___);}
#if 0
# 379
{ } 
#endif
# 380 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline void __stcs(short *ptr, short value) {int volatile ___ = 1;(void)ptr;(void)value;::exit(___);}
#if 0
# 380
{ } 
#endif
# 381 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline void __stcs(int *ptr, int value) {int volatile ___ = 1;(void)ptr;(void)value;::exit(___);}
#if 0
# 381
{ } 
#endif
# 382 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline void __stcs(long long *ptr, long long value) {int volatile ___ = 1;(void)ptr;(void)value;::exit(___);}
#if 0
# 382
{ } 
#endif
# 383 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline void __stcs(char2 *ptr, char2 value) {int volatile ___ = 1;(void)ptr;(void)value;::exit(___);}
#if 0
# 383
{ } 
#endif
# 384 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline void __stcs(char4 *ptr, char4 value) {int volatile ___ = 1;(void)ptr;(void)value;::exit(___);}
#if 0
# 384
{ } 
#endif
# 385 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline void __stcs(short2 *ptr, short2 value) {int volatile ___ = 1;(void)ptr;(void)value;::exit(___);}
#if 0
# 385
{ } 
#endif
# 386 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline void __stcs(short4 *ptr, short4 value) {int volatile ___ = 1;(void)ptr;(void)value;::exit(___);}
#if 0
# 386
{ } 
#endif
# 387 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline void __stcs(int2 *ptr, int2 value) {int volatile ___ = 1;(void)ptr;(void)value;::exit(___);}
#if 0
# 387
{ } 
#endif
# 388 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline void __stcs(int4 *ptr, int4 value) {int volatile ___ = 1;(void)ptr;(void)value;::exit(___);}
#if 0
# 388
{ } 
#endif
# 389 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline void __stcs(longlong2 *ptr, longlong2 value) {int volatile ___ = 1;(void)ptr;(void)value;::exit(___);}
#if 0
# 389
{ } 
#endif
# 391 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline void __stcs(unsigned char *ptr, unsigned char value) {int volatile ___ = 1;(void)ptr;(void)value;::exit(___);}
#if 0
# 391
{ } 
#endif
# 392 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline void __stcs(unsigned short *ptr, unsigned short value) {int volatile ___ = 1;(void)ptr;(void)value;::exit(___);}
#if 0
# 392
{ } 
#endif
# 393 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline void __stcs(unsigned *ptr, unsigned value) {int volatile ___ = 1;(void)ptr;(void)value;::exit(___);}
#if 0
# 393
{ } 
#endif
# 394 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline void __stcs(unsigned long long *ptr, unsigned long long value) {int volatile ___ = 1;(void)ptr;(void)value;::exit(___);}
#if 0
# 394
{ } 
#endif
# 395 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline void __stcs(uchar2 *ptr, uchar2 value) {int volatile ___ = 1;(void)ptr;(void)value;::exit(___);}
#if 0
# 395
{ } 
#endif
# 396 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline void __stcs(uchar4 *ptr, uchar4 value) {int volatile ___ = 1;(void)ptr;(void)value;::exit(___);}
#if 0
# 396
{ } 
#endif
# 397 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline void __stcs(ushort2 *ptr, ushort2 value) {int volatile ___ = 1;(void)ptr;(void)value;::exit(___);}
#if 0
# 397
{ } 
#endif
# 398 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline void __stcs(ushort4 *ptr, ushort4 value) {int volatile ___ = 1;(void)ptr;(void)value;::exit(___);}
#if 0
# 398
{ } 
#endif
# 399 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline void __stcs(uint2 *ptr, uint2 value) {int volatile ___ = 1;(void)ptr;(void)value;::exit(___);}
#if 0
# 399
{ } 
#endif
# 400 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline void __stcs(uint4 *ptr, uint4 value) {int volatile ___ = 1;(void)ptr;(void)value;::exit(___);}
#if 0
# 400
{ } 
#endif
# 401 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline void __stcs(ulonglong2 *ptr, ulonglong2 value) {int volatile ___ = 1;(void)ptr;(void)value;::exit(___);}
#if 0
# 401
{ } 
#endif
# 403 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline void __stcs(float *ptr, float value) {int volatile ___ = 1;(void)ptr;(void)value;::exit(___);}
#if 0
# 403
{ } 
#endif
# 404 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline void __stcs(double *ptr, double value) {int volatile ___ = 1;(void)ptr;(void)value;::exit(___);}
#if 0
# 404
{ } 
#endif
# 405 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline void __stcs(float2 *ptr, float2 value) {int volatile ___ = 1;(void)ptr;(void)value;::exit(___);}
#if 0
# 405
{ } 
#endif
# 406 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline void __stcs(float4 *ptr, float4 value) {int volatile ___ = 1;(void)ptr;(void)value;::exit(___);}
#if 0
# 406
{ } 
#endif
# 407 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline void __stcs(double2 *ptr, double2 value) {int volatile ___ = 1;(void)ptr;(void)value;::exit(___);}
#if 0
# 407
{ } 
#endif
# 411 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline void __stwt(long *ptr, long value) {int volatile ___ = 1;(void)ptr;(void)value;::exit(___);}
#if 0
# 411
{ } 
#endif
# 412 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline void __stwt(unsigned long *ptr, unsigned long value) {int volatile ___ = 1;(void)ptr;(void)value;::exit(___);}
#if 0
# 412
{ } 
#endif
# 414 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline void __stwt(char *ptr, char value) {int volatile ___ = 1;(void)ptr;(void)value;::exit(___);}
#if 0
# 414
{ } 
#endif
# 415 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline void __stwt(signed char *ptr, signed char value) {int volatile ___ = 1;(void)ptr;(void)value;::exit(___);}
#if 0
# 415
{ } 
#endif
# 416 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline void __stwt(short *ptr, short value) {int volatile ___ = 1;(void)ptr;(void)value;::exit(___);}
#if 0
# 416
{ } 
#endif
# 417 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline void __stwt(int *ptr, int value) {int volatile ___ = 1;(void)ptr;(void)value;::exit(___);}
#if 0
# 417
{ } 
#endif
# 418 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline void __stwt(long long *ptr, long long value) {int volatile ___ = 1;(void)ptr;(void)value;::exit(___);}
#if 0
# 418
{ } 
#endif
# 419 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline void __stwt(char2 *ptr, char2 value) {int volatile ___ = 1;(void)ptr;(void)value;::exit(___);}
#if 0
# 419
{ } 
#endif
# 420 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline void __stwt(char4 *ptr, char4 value) {int volatile ___ = 1;(void)ptr;(void)value;::exit(___);}
#if 0
# 420
{ } 
#endif
# 421 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline void __stwt(short2 *ptr, short2 value) {int volatile ___ = 1;(void)ptr;(void)value;::exit(___);}
#if 0
# 421
{ } 
#endif
# 422 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline void __stwt(short4 *ptr, short4 value) {int volatile ___ = 1;(void)ptr;(void)value;::exit(___);}
#if 0
# 422
{ } 
#endif
# 423 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline void __stwt(int2 *ptr, int2 value) {int volatile ___ = 1;(void)ptr;(void)value;::exit(___);}
#if 0
# 423
{ } 
#endif
# 424 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline void __stwt(int4 *ptr, int4 value) {int volatile ___ = 1;(void)ptr;(void)value;::exit(___);}
#if 0
# 424
{ } 
#endif
# 425 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline void __stwt(longlong2 *ptr, longlong2 value) {int volatile ___ = 1;(void)ptr;(void)value;::exit(___);}
#if 0
# 425
{ } 
#endif
# 427 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline void __stwt(unsigned char *ptr, unsigned char value) {int volatile ___ = 1;(void)ptr;(void)value;::exit(___);}
#if 0
# 427
{ } 
#endif
# 428 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline void __stwt(unsigned short *ptr, unsigned short value) {int volatile ___ = 1;(void)ptr;(void)value;::exit(___);}
#if 0
# 428
{ } 
#endif
# 429 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline void __stwt(unsigned *ptr, unsigned value) {int volatile ___ = 1;(void)ptr;(void)value;::exit(___);}
#if 0
# 429
{ } 
#endif
# 430 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline void __stwt(unsigned long long *ptr, unsigned long long value) {int volatile ___ = 1;(void)ptr;(void)value;::exit(___);}
#if 0
# 430
{ } 
#endif
# 431 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline void __stwt(uchar2 *ptr, uchar2 value) {int volatile ___ = 1;(void)ptr;(void)value;::exit(___);}
#if 0
# 431
{ } 
#endif
# 432 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline void __stwt(uchar4 *ptr, uchar4 value) {int volatile ___ = 1;(void)ptr;(void)value;::exit(___);}
#if 0
# 432
{ } 
#endif
# 433 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline void __stwt(ushort2 *ptr, ushort2 value) {int volatile ___ = 1;(void)ptr;(void)value;::exit(___);}
#if 0
# 433
{ } 
#endif
# 434 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline void __stwt(ushort4 *ptr, ushort4 value) {int volatile ___ = 1;(void)ptr;(void)value;::exit(___);}
#if 0
# 434
{ } 
#endif
# 435 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline void __stwt(uint2 *ptr, uint2 value) {int volatile ___ = 1;(void)ptr;(void)value;::exit(___);}
#if 0
# 435
{ } 
#endif
# 436 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline void __stwt(uint4 *ptr, uint4 value) {int volatile ___ = 1;(void)ptr;(void)value;::exit(___);}
#if 0
# 436
{ } 
#endif
# 437 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline void __stwt(ulonglong2 *ptr, ulonglong2 value) {int volatile ___ = 1;(void)ptr;(void)value;::exit(___);}
#if 0
# 437
{ } 
#endif
# 439 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline void __stwt(float *ptr, float value) {int volatile ___ = 1;(void)ptr;(void)value;::exit(___);}
#if 0
# 439
{ } 
#endif
# 440 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline void __stwt(double *ptr, double value) {int volatile ___ = 1;(void)ptr;(void)value;::exit(___);}
#if 0
# 440
{ } 
#endif
# 441 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline void __stwt(float2 *ptr, float2 value) {int volatile ___ = 1;(void)ptr;(void)value;::exit(___);}
#if 0
# 441
{ } 
#endif
# 442 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline void __stwt(float4 *ptr, float4 value) {int volatile ___ = 1;(void)ptr;(void)value;::exit(___);}
#if 0
# 442
{ } 
#endif
# 443 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline void __stwt(double2 *ptr, double2 value) {int volatile ___ = 1;(void)ptr;(void)value;::exit(___);}
#if 0
# 443
{ } 
#endif
# 460 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline unsigned __funnelshift_l(unsigned lo, unsigned hi, unsigned shift) {int volatile ___ = 1;(void)lo;(void)hi;(void)shift;::exit(___);}
#if 0
# 460
{ } 
#endif
# 472 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline unsigned __funnelshift_lc(unsigned lo, unsigned hi, unsigned shift) {int volatile ___ = 1;(void)lo;(void)hi;(void)shift;::exit(___);}
#if 0
# 472
{ } 
#endif
# 485 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline unsigned __funnelshift_r(unsigned lo, unsigned hi, unsigned shift) {int volatile ___ = 1;(void)lo;(void)hi;(void)shift;::exit(___);}
#if 0
# 485
{ } 
#endif
# 497 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_32_intrinsics.h"
__attribute__((unused)) static inline unsigned __funnelshift_rc(unsigned lo, unsigned hi, unsigned shift) {int volatile ___ = 1;(void)lo;(void)hi;(void)shift;::exit(___);}
#if 0
# 497
{ } 
#endif
# 99 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_61_intrinsics.h"
__attribute__((unused)) static inline int __dp2a_lo(int srcA, int srcB, int c) {int volatile ___ = 1;(void)srcA;(void)srcB;(void)c;::exit(___);}
#if 0
# 99
{ } 
#endif
# 110 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_61_intrinsics.h"
__attribute__((unused)) static inline unsigned __dp2a_lo(unsigned srcA, unsigned srcB, unsigned c) {int volatile ___ = 1;(void)srcA;(void)srcB;(void)c;::exit(___);}
#if 0
# 110
{ } 
#endif
# 122 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_61_intrinsics.h"
__attribute__((unused)) static inline int __dp2a_lo(short2 srcA, char4 srcB, int c) {int volatile ___ = 1;(void)srcA;(void)srcB;(void)c;::exit(___);}
#if 0
# 122
{ } 
#endif
# 133 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_61_intrinsics.h"
__attribute__((unused)) static inline unsigned __dp2a_lo(ushort2 srcA, uchar4 srcB, unsigned c) {int volatile ___ = 1;(void)srcA;(void)srcB;(void)c;::exit(___);}
#if 0
# 133
{ } 
#endif
# 145 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_61_intrinsics.h"
__attribute__((unused)) static inline int __dp2a_hi(int srcA, int srcB, int c) {int volatile ___ = 1;(void)srcA;(void)srcB;(void)c;::exit(___);}
#if 0
# 145
{ } 
#endif
# 156 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_61_intrinsics.h"
__attribute__((unused)) static inline unsigned __dp2a_hi(unsigned srcA, unsigned srcB, unsigned c) {int volatile ___ = 1;(void)srcA;(void)srcB;(void)c;::exit(___);}
#if 0
# 156
{ } 
#endif
# 168 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_61_intrinsics.h"
__attribute__((unused)) static inline int __dp2a_hi(short2 srcA, char4 srcB, int c) {int volatile ___ = 1;(void)srcA;(void)srcB;(void)c;::exit(___);}
#if 0
# 168
{ } 
#endif
# 179 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_61_intrinsics.h"
__attribute__((unused)) static inline unsigned __dp2a_hi(ushort2 srcA, uchar4 srcB, unsigned c) {int volatile ___ = 1;(void)srcA;(void)srcB;(void)c;::exit(___);}
#if 0
# 179
{ } 
#endif
# 194 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_61_intrinsics.h"
__attribute__((unused)) static inline int __dp4a(int srcA, int srcB, int c) {int volatile ___ = 1;(void)srcA;(void)srcB;(void)c;::exit(___);}
#if 0
# 194
{ } 
#endif
# 203 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_61_intrinsics.h"
__attribute__((unused)) static inline unsigned __dp4a(unsigned srcA, unsigned srcB, unsigned c) {int volatile ___ = 1;(void)srcA;(void)srcB;(void)c;::exit(___);}
#if 0
# 203
{ } 
#endif
# 213 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_61_intrinsics.h"
__attribute__((unused)) static inline int __dp4a(char4 srcA, char4 srcB, int c) {int volatile ___ = 1;(void)srcA;(void)srcB;(void)c;::exit(___);}
#if 0
# 213
{ } 
#endif
# 222 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/sm_61_intrinsics.h"
__attribute__((unused)) static inline unsigned __dp4a(uchar4 srcA, uchar4 srcB, unsigned c) {int volatile ___ = 1;(void)srcA;(void)srcB;(void)c;::exit(___);}
#if 0
# 222
{ } 
#endif
# 93 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/sm_70_rt.h"
__attribute__((unused)) static inline unsigned __match_any_sync(unsigned mask, unsigned value) {int volatile ___ = 1;(void)mask;(void)value;::exit(___);}
#if 0
# 93
{ } 
#endif
# 94 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/sm_70_rt.h"
__attribute__((unused)) static inline unsigned __match_any_sync(unsigned mask, int value) {int volatile ___ = 1;(void)mask;(void)value;::exit(___);}
#if 0
# 94
{ } 
#endif
# 95 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/sm_70_rt.h"
__attribute__((unused)) static inline unsigned __match_any_sync(unsigned mask, unsigned long value) {int volatile ___ = 1;(void)mask;(void)value;::exit(___);}
#if 0
# 95
{ } 
#endif
# 96 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/sm_70_rt.h"
__attribute__((unused)) static inline unsigned __match_any_sync(unsigned mask, long value) {int volatile ___ = 1;(void)mask;(void)value;::exit(___);}
#if 0
# 96
{ } 
#endif
# 97 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/sm_70_rt.h"
__attribute__((unused)) static inline unsigned __match_any_sync(unsigned mask, unsigned long long value) {int volatile ___ = 1;(void)mask;(void)value;::exit(___);}
#if 0
# 97
{ } 
#endif
# 98 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/sm_70_rt.h"
__attribute__((unused)) static inline unsigned __match_any_sync(unsigned mask, long long value) {int volatile ___ = 1;(void)mask;(void)value;::exit(___);}
#if 0
# 98
{ } 
#endif
# 99 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/sm_70_rt.h"
__attribute__((unused)) static inline unsigned __match_any_sync(unsigned mask, float value) {int volatile ___ = 1;(void)mask;(void)value;::exit(___);}
#if 0
# 99
{ } 
#endif
# 100 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/sm_70_rt.h"
__attribute__((unused)) static inline unsigned __match_any_sync(unsigned mask, double value) {int volatile ___ = 1;(void)mask;(void)value;::exit(___);}
#if 0
# 100
{ } 
#endif
# 102 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/sm_70_rt.h"
__attribute__((unused)) static inline unsigned __match_all_sync(unsigned mask, unsigned value, int *pred) {int volatile ___ = 1;(void)mask;(void)value;(void)pred;::exit(___);}
#if 0
# 102
{ } 
#endif
# 103 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/sm_70_rt.h"
__attribute__((unused)) static inline unsigned __match_all_sync(unsigned mask, int value, int *pred) {int volatile ___ = 1;(void)mask;(void)value;(void)pred;::exit(___);}
#if 0
# 103
{ } 
#endif
# 104 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/sm_70_rt.h"
__attribute__((unused)) static inline unsigned __match_all_sync(unsigned mask, unsigned long value, int *pred) {int volatile ___ = 1;(void)mask;(void)value;(void)pred;::exit(___);}
#if 0
# 104
{ } 
#endif
# 105 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/sm_70_rt.h"
__attribute__((unused)) static inline unsigned __match_all_sync(unsigned mask, long value, int *pred) {int volatile ___ = 1;(void)mask;(void)value;(void)pred;::exit(___);}
#if 0
# 105
{ } 
#endif
# 106 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/sm_70_rt.h"
__attribute__((unused)) static inline unsigned __match_all_sync(unsigned mask, unsigned long long value, int *pred) {int volatile ___ = 1;(void)mask;(void)value;(void)pred;::exit(___);}
#if 0
# 106
{ } 
#endif
# 107 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/sm_70_rt.h"
__attribute__((unused)) static inline unsigned __match_all_sync(unsigned mask, long long value, int *pred) {int volatile ___ = 1;(void)mask;(void)value;(void)pred;::exit(___);}
#if 0
# 107
{ } 
#endif
# 108 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/sm_70_rt.h"
__attribute__((unused)) static inline unsigned __match_all_sync(unsigned mask, float value, int *pred) {int volatile ___ = 1;(void)mask;(void)value;(void)pred;::exit(___);}
#if 0
# 108
{ } 
#endif
# 109 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/sm_70_rt.h"
__attribute__((unused)) static inline unsigned __match_all_sync(unsigned mask, double value, int *pred) {int volatile ___ = 1;(void)mask;(void)value;(void)pred;::exit(___);}
#if 0
# 109
{ } 
#endif
# 111 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/sm_70_rt.h"
__attribute__((unused)) static inline void __nanosleep(unsigned ns) {int volatile ___ = 1;(void)ns;::exit(___);}
#if 0
# 111
{ } 
#endif
# 113 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/sm_70_rt.h"
__attribute__((unused)) static inline unsigned short atomicCAS(unsigned short *address, unsigned short compare, unsigned short val) {int volatile ___ = 1;(void)address;(void)compare;(void)val;::exit(___);}
#if 0
# 113
{ } 
#endif
# 93 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/sm_80_rt.h"
__attribute__((unused)) static inline unsigned __reduce_add_sync(unsigned mask, unsigned value) {int volatile ___ = 1;(void)mask;(void)value;::exit(___);}
#if 0
# 93
{ } 
#endif
# 94 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/sm_80_rt.h"
__attribute__((unused)) static inline unsigned __reduce_min_sync(unsigned mask, unsigned value) {int volatile ___ = 1;(void)mask;(void)value;::exit(___);}
#if 0
# 94
{ } 
#endif
# 95 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/sm_80_rt.h"
__attribute__((unused)) static inline unsigned __reduce_max_sync(unsigned mask, unsigned value) {int volatile ___ = 1;(void)mask;(void)value;::exit(___);}
#if 0
# 95
{ } 
#endif
# 97 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/sm_80_rt.h"
__attribute__((unused)) static inline int __reduce_add_sync(unsigned mask, int value) {int volatile ___ = 1;(void)mask;(void)value;::exit(___);}
#if 0
# 97
{ } 
#endif
# 98 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/sm_80_rt.h"
__attribute__((unused)) static inline int __reduce_min_sync(unsigned mask, int value) {int volatile ___ = 1;(void)mask;(void)value;::exit(___);}
#if 0
# 98
{ } 
#endif
# 99 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/sm_80_rt.h"
__attribute__((unused)) static inline int __reduce_max_sync(unsigned mask, int value) {int volatile ___ = 1;(void)mask;(void)value;::exit(___);}
#if 0
# 99
{ } 
#endif
# 101 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/sm_80_rt.h"
__attribute__((unused)) static inline unsigned __reduce_and_sync(unsigned mask, unsigned value) {int volatile ___ = 1;(void)mask;(void)value;::exit(___);}
#if 0
# 101
{ } 
#endif
# 102 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/sm_80_rt.h"
__attribute__((unused)) static inline unsigned __reduce_or_sync(unsigned mask, unsigned value) {int volatile ___ = 1;(void)mask;(void)value;::exit(___);}
#if 0
# 102
{ } 
#endif
# 103 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/sm_80_rt.h"
__attribute__((unused)) static inline unsigned __reduce_xor_sync(unsigned mask, unsigned value) {int volatile ___ = 1;(void)mask;(void)value;::exit(___);}
#if 0
# 103
{ } 
#endif
# 106 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/sm_80_rt.h"
extern "C" {
# 107
__attribute__((unused)) inline void *__nv_associate_access_property(const void *ptr, unsigned long long 
# 108
property) {int volatile ___ = 1;(void)ptr;(void)property;
# 112
::exit(___);}
#if 0
# 108
{ 
# 109
__attribute__((unused)) extern void *__nv_associate_access_property_impl(const void *, unsigned long long); 
# 111
return __nv_associate_access_property_impl(ptr, property); 
# 112
} 
#endif
# 114 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/sm_80_rt.h"
__attribute__((unused)) inline void __nv_memcpy_async_shared_global_4(void *dst, const void *
# 115
src, unsigned 
# 116
src_size) {int volatile ___ = 1;(void)dst;(void)src;(void)src_size;
# 121
::exit(___);}
#if 0
# 116
{ 
# 117
__attribute__((unused)) extern void __nv_memcpy_async_shared_global_4_impl(void *, const void *, unsigned); 
# 120
__nv_memcpy_async_shared_global_4_impl(dst, src, src_size); 
# 121
} 
#endif
# 123 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/sm_80_rt.h"
__attribute__((unused)) inline void __nv_memcpy_async_shared_global_8(void *dst, const void *
# 124
src, unsigned 
# 125
src_size) {int volatile ___ = 1;(void)dst;(void)src;(void)src_size;
# 130
::exit(___);}
#if 0
# 125
{ 
# 126
__attribute__((unused)) extern void __nv_memcpy_async_shared_global_8_impl(void *, const void *, unsigned); 
# 129
__nv_memcpy_async_shared_global_8_impl(dst, src, src_size); 
# 130
} 
#endif
# 132 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/sm_80_rt.h"
__attribute__((unused)) inline void __nv_memcpy_async_shared_global_16(void *dst, const void *
# 133
src, unsigned 
# 134
src_size) {int volatile ___ = 1;(void)dst;(void)src;(void)src_size;
# 139
::exit(___);}
#if 0
# 134
{ 
# 135
__attribute__((unused)) extern void __nv_memcpy_async_shared_global_16_impl(void *, const void *, unsigned); 
# 138
__nv_memcpy_async_shared_global_16_impl(dst, src, src_size); 
# 139
} 
#endif
# 141 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/sm_80_rt.h"
}
# 89 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/sm_90_rt.h"
__attribute__((unused)) static inline unsigned __isCtaShared(const void *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 89
{ } 
#endif
# 90 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/sm_90_rt.h"
__attribute__((unused)) static inline unsigned __isClusterShared(const void *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 90
{ } 
#endif
# 91 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/sm_90_rt.h"
__attribute__((unused)) static inline void *__cluster_map_shared_rank(const void *ptr, unsigned target_block_rank) {int volatile ___ = 1;(void)ptr;(void)target_block_rank;::exit(___);}
#if 0
# 91
{ } 
#endif
# 92 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/sm_90_rt.h"
__attribute__((unused)) static inline unsigned __cluster_query_shared_rank(const void *ptr) {int volatile ___ = 1;(void)ptr;::exit(___);}
#if 0
# 92
{ } 
#endif
# 93 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/sm_90_rt.h"
__attribute__((unused)) static inline uint2 __cluster_map_shared_multicast(const void *ptr, unsigned cluster_cta_mask) {int volatile ___ = 1;(void)ptr;(void)cluster_cta_mask;::exit(___);}
#if 0
# 93
{ } 
#endif
# 94 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/sm_90_rt.h"
__attribute__((unused)) static inline unsigned __clusterDimIsSpecified() {int volatile ___ = 1;::exit(___);}
#if 0
# 94
{ } 
#endif
# 95 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/sm_90_rt.h"
__attribute__((unused)) static inline dim3 __clusterDim() {int volatile ___ = 1;::exit(___);}
#if 0
# 95
{ } 
#endif
# 96 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/sm_90_rt.h"
__attribute__((unused)) static inline dim3 __clusterRelativeBlockIdx() {int volatile ___ = 1;::exit(___);}
#if 0
# 96
{ } 
#endif
# 97 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/sm_90_rt.h"
__attribute__((unused)) static inline dim3 __clusterGridDimInClusters() {int volatile ___ = 1;::exit(___);}
#if 0
# 97
{ } 
#endif
# 98 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/sm_90_rt.h"
__attribute__((unused)) static inline dim3 __clusterIdx() {int volatile ___ = 1;::exit(___);}
#if 0
# 98
{ } 
#endif
# 99 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/sm_90_rt.h"
__attribute__((unused)) static inline unsigned __clusterRelativeBlockRank() {int volatile ___ = 1;::exit(___);}
#if 0
# 99
{ } 
#endif
# 100 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/sm_90_rt.h"
__attribute__((unused)) static inline unsigned __clusterSizeInBlocks() {int volatile ___ = 1;::exit(___);}
#if 0
# 100
{ } 
#endif
# 101 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/sm_90_rt.h"
__attribute__((unused)) static inline void __cluster_barrier_arrive() {int volatile ___ = 1;::exit(___);}
#if 0
# 101
{ } 
#endif
# 102 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/sm_90_rt.h"
__attribute__((unused)) static inline void __cluster_barrier_arrive_relaxed() {int volatile ___ = 1;::exit(___);}
#if 0
# 102
{ } 
#endif
# 103 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/sm_90_rt.h"
__attribute__((unused)) static inline void __cluster_barrier_wait() {int volatile ___ = 1;::exit(___);}
#if 0
# 103
{ } 
#endif
# 104 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/sm_90_rt.h"
__attribute__((unused)) static inline void __threadfence_cluster() {int volatile ___ = 1;::exit(___);}
#if 0
# 104
{ } 
#endif
# 106 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/sm_90_rt.h"
__attribute__((unused)) static inline float2 atomicAdd(float2 *address, float2 val) {int volatile ___ = 1;(void)address;(void)val;::exit(___);}
#if 0
# 106
{ } 
#endif
# 107 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/sm_90_rt.h"
__attribute__((unused)) static inline float2 atomicAdd_block(float2 *address, float2 val) {int volatile ___ = 1;(void)address;(void)val;::exit(___);}
#if 0
# 107
{ } 
#endif
# 108 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/sm_90_rt.h"
__attribute__((unused)) static inline float2 atomicAdd_system(float2 *address, float2 val) {int volatile ___ = 1;(void)address;(void)val;::exit(___);}
#if 0
# 108
{ } 
#endif
# 109 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/sm_90_rt.h"
__attribute__((unused)) static inline float4 atomicAdd(float4 *address, float4 val) {int volatile ___ = 1;(void)address;(void)val;::exit(___);}
#if 0
# 109
{ } 
#endif
# 110 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/sm_90_rt.h"
__attribute__((unused)) static inline float4 atomicAdd_block(float4 *address, float4 val) {int volatile ___ = 1;(void)address;(void)val;::exit(___);}
#if 0
# 110
{ } 
#endif
# 111 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/sm_90_rt.h"
__attribute__((unused)) static inline float4 atomicAdd_system(float4 *address, float4 val) {int volatile ___ = 1;(void)address;(void)val;::exit(___);}
#if 0
# 111
{ } 
#endif
# 65 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/texture_indirect_functions.h"
template< class T> struct __nv_itex_trait { }; 
# 66
template<> struct __nv_itex_trait< char>  { typedef void type; }; 
# 67
template<> struct __nv_itex_trait< signed char>  { typedef void type; }; 
# 68
template<> struct __nv_itex_trait< char1>  { typedef void type; }; 
# 69
template<> struct __nv_itex_trait< char2>  { typedef void type; }; 
# 70
template<> struct __nv_itex_trait< char4>  { typedef void type; }; 
# 71
template<> struct __nv_itex_trait< unsigned char>  { typedef void type; }; 
# 72
template<> struct __nv_itex_trait< uchar1>  { typedef void type; }; 
# 73
template<> struct __nv_itex_trait< uchar2>  { typedef void type; }; 
# 74
template<> struct __nv_itex_trait< uchar4>  { typedef void type; }; 
# 75
template<> struct __nv_itex_trait< short>  { typedef void type; }; 
# 76
template<> struct __nv_itex_trait< short1>  { typedef void type; }; 
# 77
template<> struct __nv_itex_trait< short2>  { typedef void type; }; 
# 78
template<> struct __nv_itex_trait< short4>  { typedef void type; }; 
# 79
template<> struct __nv_itex_trait< unsigned short>  { typedef void type; }; 
# 80
template<> struct __nv_itex_trait< ushort1>  { typedef void type; }; 
# 81
template<> struct __nv_itex_trait< ushort2>  { typedef void type; }; 
# 82
template<> struct __nv_itex_trait< ushort4>  { typedef void type; }; 
# 83
template<> struct __nv_itex_trait< int>  { typedef void type; }; 
# 84
template<> struct __nv_itex_trait< int1>  { typedef void type; }; 
# 85
template<> struct __nv_itex_trait< int2>  { typedef void type; }; 
# 86
template<> struct __nv_itex_trait< int4>  { typedef void type; }; 
# 87
template<> struct __nv_itex_trait< unsigned>  { typedef void type; }; 
# 88
template<> struct __nv_itex_trait< uint1>  { typedef void type; }; 
# 89
template<> struct __nv_itex_trait< uint2>  { typedef void type; }; 
# 90
template<> struct __nv_itex_trait< uint4>  { typedef void type; }; 
# 101 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/texture_indirect_functions.h"
template<> struct __nv_itex_trait< float>  { typedef void type; }; 
# 102
template<> struct __nv_itex_trait< float1>  { typedef void type; }; 
# 103
template<> struct __nv_itex_trait< float2>  { typedef void type; }; 
# 104
template<> struct __nv_itex_trait< float4>  { typedef void type; }; 
# 108
template< class T> __attribute__((unused)) static typename __nv_itex_trait< T> ::type 
# 109
tex1Dfetch(T *ptr, cudaTextureObject_t obj, int x) 
# 110
{int volatile ___ = 1;(void)ptr;(void)obj;(void)x;
# 112
::exit(___);}
#if 0
# 110
{ 
# 111
__nv_tex_surf_handler("__itex1Dfetch", ptr, obj, x); 
# 112
} 
#endif
# 114 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/texture_indirect_functions.h"
template< class T> __attribute__((unused)) static T 
# 115
tex1Dfetch(cudaTextureObject_t texObject, int x) 
# 116
{int volatile ___ = 1;(void)texObject;(void)x;
# 120
::exit(___);}
#if 0
# 116
{ 
# 117
T ret; 
# 118
tex1Dfetch(&ret, texObject, x); 
# 119
return ret; 
# 120
} 
#endif
# 122 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/texture_indirect_functions.h"
template< class T> __attribute__((unused)) static typename __nv_itex_trait< T> ::type 
# 123
tex1D(T *ptr, cudaTextureObject_t obj, float x) 
# 124
{int volatile ___ = 1;(void)ptr;(void)obj;(void)x;
# 126
::exit(___);}
#if 0
# 124
{ 
# 125
__nv_tex_surf_handler("__itex1D", ptr, obj, x); 
# 126
} 
#endif
# 129 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/texture_indirect_functions.h"
template< class T> __attribute__((unused)) static T 
# 130
tex1D(cudaTextureObject_t texObject, float x) 
# 131
{int volatile ___ = 1;(void)texObject;(void)x;
# 135
::exit(___);}
#if 0
# 131
{ 
# 132
T ret; 
# 133
tex1D(&ret, texObject, x); 
# 134
return ret; 
# 135
} 
#endif
# 138 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/texture_indirect_functions.h"
template< class T> __attribute__((unused)) static typename __nv_itex_trait< T> ::type 
# 139
tex2D(T *ptr, cudaTextureObject_t obj, float x, float y) 
# 140
{int volatile ___ = 1;(void)ptr;(void)obj;(void)x;(void)y;
# 142
::exit(___);}
#if 0
# 140
{ 
# 141
__nv_tex_surf_handler("__itex2D", ptr, obj, x, y); 
# 142
} 
#endif
# 144 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/texture_indirect_functions.h"
template< class T> __attribute__((unused)) static T 
# 145
tex2D(cudaTextureObject_t texObject, float x, float y) 
# 146
{int volatile ___ = 1;(void)texObject;(void)x;(void)y;
# 150
::exit(___);}
#if 0
# 146
{ 
# 147
T ret; 
# 148
tex2D(&ret, texObject, x, y); 
# 149
return ret; 
# 150
} 
#endif
# 153 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/texture_indirect_functions.h"
template< class T> __attribute__((unused)) static typename __nv_itex_trait< T> ::type 
# 154
tex2D(T *ptr, cudaTextureObject_t obj, float x, float y, bool *
# 155
isResident) 
# 156
{int volatile ___ = 1;(void)ptr;(void)obj;(void)x;(void)y;(void)isResident;
# 160
::exit(___);}
#if 0
# 156
{ 
# 157
unsigned char res; 
# 158
__nv_tex_surf_handler("__itex2D_sparse", ptr, obj, x, y, &res); 
# 159
(*isResident) = (res != 0); 
# 160
} 
#endif
# 162 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/texture_indirect_functions.h"
template< class T> __attribute__((unused)) static T 
# 163
tex2D(cudaTextureObject_t texObject, float x, float y, bool *isResident) 
# 164
{int volatile ___ = 1;(void)texObject;(void)x;(void)y;(void)isResident;
# 168
::exit(___);}
#if 0
# 164
{ 
# 165
T ret; 
# 166
tex2D(&ret, texObject, x, y, isResident); 
# 167
return ret; 
# 168
} 
#endif
# 173 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/texture_indirect_functions.h"
template< class T> __attribute__((unused)) static typename __nv_itex_trait< T> ::type 
# 174
tex3D(T *ptr, cudaTextureObject_t obj, float x, float y, float z) 
# 175
{int volatile ___ = 1;(void)ptr;(void)obj;(void)x;(void)y;(void)z;
# 177
::exit(___);}
#if 0
# 175
{ 
# 176
__nv_tex_surf_handler("__itex3D", ptr, obj, x, y, z); 
# 177
} 
#endif
# 179 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/texture_indirect_functions.h"
template< class T> __attribute__((unused)) static T 
# 180
tex3D(cudaTextureObject_t texObject, float x, float y, float z) 
# 181
{int volatile ___ = 1;(void)texObject;(void)x;(void)y;(void)z;
# 185
::exit(___);}
#if 0
# 181
{ 
# 182
T ret; 
# 183
tex3D(&ret, texObject, x, y, z); 
# 184
return ret; 
# 185
} 
#endif
# 188 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/texture_indirect_functions.h"
template< class T> __attribute__((unused)) static typename __nv_itex_trait< T> ::type 
# 189
tex3D(T *ptr, cudaTextureObject_t obj, float x, float y, float z, bool *
# 190
isResident) 
# 191
{int volatile ___ = 1;(void)ptr;(void)obj;(void)x;(void)y;(void)z;(void)isResident;
# 195
::exit(___);}
#if 0
# 191
{ 
# 192
unsigned char res; 
# 193
__nv_tex_surf_handler("__itex3D_sparse", ptr, obj, x, y, z, &res); 
# 194
(*isResident) = (res != 0); 
# 195
} 
#endif
# 197 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/texture_indirect_functions.h"
template< class T> __attribute__((unused)) static T 
# 198
tex3D(cudaTextureObject_t texObject, float x, float y, float z, bool *isResident) 
# 199
{int volatile ___ = 1;(void)texObject;(void)x;(void)y;(void)z;(void)isResident;
# 203
::exit(___);}
#if 0
# 199
{ 
# 200
T ret; 
# 201
tex3D(&ret, texObject, x, y, z, isResident); 
# 202
return ret; 
# 203
} 
#endif
# 207 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/texture_indirect_functions.h"
template< class T> __attribute__((unused)) static typename __nv_itex_trait< T> ::type 
# 208
tex1DLayered(T *ptr, cudaTextureObject_t obj, float x, int layer) 
# 209
{int volatile ___ = 1;(void)ptr;(void)obj;(void)x;(void)layer;
# 211
::exit(___);}
#if 0
# 209
{ 
# 210
__nv_tex_surf_handler("__itex1DLayered", ptr, obj, x, layer); 
# 211
} 
#endif
# 213 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/texture_indirect_functions.h"
template< class T> __attribute__((unused)) static T 
# 214
tex1DLayered(cudaTextureObject_t texObject, float x, int layer) 
# 215
{int volatile ___ = 1;(void)texObject;(void)x;(void)layer;
# 219
::exit(___);}
#if 0
# 215
{ 
# 216
T ret; 
# 217
tex1DLayered(&ret, texObject, x, layer); 
# 218
return ret; 
# 219
} 
#endif
# 221 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/texture_indirect_functions.h"
template< class T> __attribute__((unused)) static typename __nv_itex_trait< T> ::type 
# 222
tex2DLayered(T *ptr, cudaTextureObject_t obj, float x, float y, int layer) 
# 223
{int volatile ___ = 1;(void)ptr;(void)obj;(void)x;(void)y;(void)layer;
# 225
::exit(___);}
#if 0
# 223
{ 
# 224
__nv_tex_surf_handler("__itex2DLayered", ptr, obj, x, y, layer); 
# 225
} 
#endif
# 227 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/texture_indirect_functions.h"
template< class T> __attribute__((unused)) static T 
# 228
tex2DLayered(cudaTextureObject_t texObject, float x, float y, int layer) 
# 229
{int volatile ___ = 1;(void)texObject;(void)x;(void)y;(void)layer;
# 233
::exit(___);}
#if 0
# 229
{ 
# 230
T ret; 
# 231
tex2DLayered(&ret, texObject, x, y, layer); 
# 232
return ret; 
# 233
} 
#endif
# 236 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/texture_indirect_functions.h"
template< class T> __attribute__((unused)) static typename __nv_itex_trait< T> ::type 
# 237
tex2DLayered(T *ptr, cudaTextureObject_t obj, float x, float y, int layer, bool *isResident) 
# 238
{int volatile ___ = 1;(void)ptr;(void)obj;(void)x;(void)y;(void)layer;(void)isResident;
# 242
::exit(___);}
#if 0
# 238
{ 
# 239
unsigned char res; 
# 240
__nv_tex_surf_handler("__itex2DLayered_sparse", ptr, obj, x, y, layer, &res); 
# 241
(*isResident) = (res != 0); 
# 242
} 
#endif
# 244 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/texture_indirect_functions.h"
template< class T> __attribute__((unused)) static T 
# 245
tex2DLayered(cudaTextureObject_t texObject, float x, float y, int layer, bool *isResident) 
# 246
{int volatile ___ = 1;(void)texObject;(void)x;(void)y;(void)layer;(void)isResident;
# 250
::exit(___);}
#if 0
# 246
{ 
# 247
T ret; 
# 248
tex2DLayered(&ret, texObject, x, y, layer, isResident); 
# 249
return ret; 
# 250
} 
#endif
# 254 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/texture_indirect_functions.h"
template< class T> __attribute__((unused)) static typename __nv_itex_trait< T> ::type 
# 255
texCubemap(T *ptr, cudaTextureObject_t obj, float x, float y, float z) 
# 256
{int volatile ___ = 1;(void)ptr;(void)obj;(void)x;(void)y;(void)z;
# 258
::exit(___);}
#if 0
# 256
{ 
# 257
__nv_tex_surf_handler("__itexCubemap", ptr, obj, x, y, z); 
# 258
} 
#endif
# 261 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/texture_indirect_functions.h"
template< class T> __attribute__((unused)) static T 
# 262
texCubemap(cudaTextureObject_t texObject, float x, float y, float z) 
# 263
{int volatile ___ = 1;(void)texObject;(void)x;(void)y;(void)z;
# 267
::exit(___);}
#if 0
# 263
{ 
# 264
T ret; 
# 265
texCubemap(&ret, texObject, x, y, z); 
# 266
return ret; 
# 267
} 
#endif
# 270 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/texture_indirect_functions.h"
template< class T> __attribute__((unused)) static typename __nv_itex_trait< T> ::type 
# 271
texCubemapLayered(T *ptr, cudaTextureObject_t obj, float x, float y, float z, int layer) 
# 272
{int volatile ___ = 1;(void)ptr;(void)obj;(void)x;(void)y;(void)z;(void)layer;
# 274
::exit(___);}
#if 0
# 272
{ 
# 273
__nv_tex_surf_handler("__itexCubemapLayered", ptr, obj, x, y, z, layer); 
# 274
} 
#endif
# 276 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/texture_indirect_functions.h"
template< class T> __attribute__((unused)) static T 
# 277
texCubemapLayered(cudaTextureObject_t texObject, float x, float y, float z, int layer) 
# 278
{int volatile ___ = 1;(void)texObject;(void)x;(void)y;(void)z;(void)layer;
# 282
::exit(___);}
#if 0
# 278
{ 
# 279
T ret; 
# 280
texCubemapLayered(&ret, texObject, x, y, z, layer); 
# 281
return ret; 
# 282
} 
#endif
# 284 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/texture_indirect_functions.h"
template< class T> __attribute__((unused)) static typename __nv_itex_trait< T> ::type 
# 285
tex2Dgather(T *ptr, cudaTextureObject_t obj, float x, float y, int comp = 0) 
# 286
{int volatile ___ = 1;(void)ptr;(void)obj;(void)x;(void)y;(void)comp;
# 288
::exit(___);}
#if 0
# 286
{ 
# 287
__nv_tex_surf_handler("__itex2Dgather", ptr, obj, x, y, comp); 
# 288
} 
#endif
# 290 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/texture_indirect_functions.h"
template< class T> __attribute__((unused)) static T 
# 291
tex2Dgather(cudaTextureObject_t to, float x, float y, int comp = 0) 
# 292
{int volatile ___ = 1;(void)to;(void)x;(void)y;(void)comp;
# 296
::exit(___);}
#if 0
# 292
{ 
# 293
T ret; 
# 294
tex2Dgather(&ret, to, x, y, comp); 
# 295
return ret; 
# 296
} 
#endif
# 299 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/texture_indirect_functions.h"
template< class T> __attribute__((unused)) static typename __nv_itex_trait< T> ::type 
# 300
tex2Dgather(T *ptr, cudaTextureObject_t obj, float x, float y, bool *isResident, int comp = 0) 
# 301
{int volatile ___ = 1;(void)ptr;(void)obj;(void)x;(void)y;(void)isResident;(void)comp;
# 305
::exit(___);}
#if 0
# 301
{ 
# 302
unsigned char res; 
# 303
__nv_tex_surf_handler("__itex2Dgather_sparse", ptr, obj, x, y, comp, &res); 
# 304
(*isResident) = (res != 0); 
# 305
} 
#endif
# 307 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/texture_indirect_functions.h"
template< class T> __attribute__((unused)) static T 
# 308
tex2Dgather(cudaTextureObject_t to, float x, float y, bool *isResident, int comp = 0) 
# 309
{int volatile ___ = 1;(void)to;(void)x;(void)y;(void)isResident;(void)comp;
# 313
::exit(___);}
#if 0
# 309
{ 
# 310
T ret; 
# 311
tex2Dgather(&ret, to, x, y, isResident, comp); 
# 312
return ret; 
# 313
} 
#endif
# 317 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/texture_indirect_functions.h"
template< class T> __attribute__((unused)) static typename __nv_itex_trait< T> ::type 
# 318
tex1DLod(T *ptr, cudaTextureObject_t obj, float x, float level) 
# 319
{int volatile ___ = 1;(void)ptr;(void)obj;(void)x;(void)level;
# 321
::exit(___);}
#if 0
# 319
{ 
# 320
__nv_tex_surf_handler("__itex1DLod", ptr, obj, x, level); 
# 321
} 
#endif
# 323 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/texture_indirect_functions.h"
template< class T> __attribute__((unused)) static T 
# 324
tex1DLod(cudaTextureObject_t texObject, float x, float level) 
# 325
{int volatile ___ = 1;(void)texObject;(void)x;(void)level;
# 329
::exit(___);}
#if 0
# 325
{ 
# 326
T ret; 
# 327
tex1DLod(&ret, texObject, x, level); 
# 328
return ret; 
# 329
} 
#endif
# 332 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/texture_indirect_functions.h"
template< class T> __attribute__((unused)) static typename __nv_itex_trait< T> ::type 
# 333
tex2DLod(T *ptr, cudaTextureObject_t obj, float x, float y, float level) 
# 334
{int volatile ___ = 1;(void)ptr;(void)obj;(void)x;(void)y;(void)level;
# 336
::exit(___);}
#if 0
# 334
{ 
# 335
__nv_tex_surf_handler("__itex2DLod", ptr, obj, x, y, level); 
# 336
} 
#endif
# 338 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/texture_indirect_functions.h"
template< class T> __attribute__((unused)) static T 
# 339
tex2DLod(cudaTextureObject_t texObject, float x, float y, float level) 
# 340
{int volatile ___ = 1;(void)texObject;(void)x;(void)y;(void)level;
# 344
::exit(___);}
#if 0
# 340
{ 
# 341
T ret; 
# 342
tex2DLod(&ret, texObject, x, y, level); 
# 343
return ret; 
# 344
} 
#endif
# 348 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/texture_indirect_functions.h"
template< class T> __attribute__((unused)) static typename __nv_itex_trait< T> ::type 
# 349
tex2DLod(T *ptr, cudaTextureObject_t obj, float x, float y, float level, bool *isResident) 
# 350
{int volatile ___ = 1;(void)ptr;(void)obj;(void)x;(void)y;(void)level;(void)isResident;
# 354
::exit(___);}
#if 0
# 350
{ 
# 351
unsigned char res; 
# 352
__nv_tex_surf_handler("__itex2DLod_sparse", ptr, obj, x, y, level, &res); 
# 353
(*isResident) = (res != 0); 
# 354
} 
#endif
# 356 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/texture_indirect_functions.h"
template< class T> __attribute__((unused)) static T 
# 357
tex2DLod(cudaTextureObject_t texObject, float x, float y, float level, bool *isResident) 
# 358
{int volatile ___ = 1;(void)texObject;(void)x;(void)y;(void)level;(void)isResident;
# 362
::exit(___);}
#if 0
# 358
{ 
# 359
T ret; 
# 360
tex2DLod(&ret, texObject, x, y, level, isResident); 
# 361
return ret; 
# 362
} 
#endif
# 367 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/texture_indirect_functions.h"
template< class T> __attribute__((unused)) static typename __nv_itex_trait< T> ::type 
# 368
tex3DLod(T *ptr, cudaTextureObject_t obj, float x, float y, float z, float level) 
# 369
{int volatile ___ = 1;(void)ptr;(void)obj;(void)x;(void)y;(void)z;(void)level;
# 371
::exit(___);}
#if 0
# 369
{ 
# 370
__nv_tex_surf_handler("__itex3DLod", ptr, obj, x, y, z, level); 
# 371
} 
#endif
# 373 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/texture_indirect_functions.h"
template< class T> __attribute__((unused)) static T 
# 374
tex3DLod(cudaTextureObject_t texObject, float x, float y, float z, float level) 
# 375
{int volatile ___ = 1;(void)texObject;(void)x;(void)y;(void)z;(void)level;
# 379
::exit(___);}
#if 0
# 375
{ 
# 376
T ret; 
# 377
tex3DLod(&ret, texObject, x, y, z, level); 
# 378
return ret; 
# 379
} 
#endif
# 382 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/texture_indirect_functions.h"
template< class T> __attribute__((unused)) static typename __nv_itex_trait< T> ::type 
# 383
tex3DLod(T *ptr, cudaTextureObject_t obj, float x, float y, float z, float level, bool *isResident) 
# 384
{int volatile ___ = 1;(void)ptr;(void)obj;(void)x;(void)y;(void)z;(void)level;(void)isResident;
# 388
::exit(___);}
#if 0
# 384
{ 
# 385
unsigned char res; 
# 386
__nv_tex_surf_handler("__itex3DLod_sparse", ptr, obj, x, y, z, level, &res); 
# 387
(*isResident) = (res != 0); 
# 388
} 
#endif
# 390 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/texture_indirect_functions.h"
template< class T> __attribute__((unused)) static T 
# 391
tex3DLod(cudaTextureObject_t texObject, float x, float y, float z, float level, bool *isResident) 
# 392
{int volatile ___ = 1;(void)texObject;(void)x;(void)y;(void)z;(void)level;(void)isResident;
# 396
::exit(___);}
#if 0
# 392
{ 
# 393
T ret; 
# 394
tex3DLod(&ret, texObject, x, y, z, level, isResident); 
# 395
return ret; 
# 396
} 
#endif
# 401 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/texture_indirect_functions.h"
template< class T> __attribute__((unused)) static typename __nv_itex_trait< T> ::type 
# 402
tex1DLayeredLod(T *ptr, cudaTextureObject_t obj, float x, int layer, float level) 
# 403
{int volatile ___ = 1;(void)ptr;(void)obj;(void)x;(void)layer;(void)level;
# 405
::exit(___);}
#if 0
# 403
{ 
# 404
__nv_tex_surf_handler("__itex1DLayeredLod", ptr, obj, x, layer, level); 
# 405
} 
#endif
# 407 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/texture_indirect_functions.h"
template< class T> __attribute__((unused)) static T 
# 408
tex1DLayeredLod(cudaTextureObject_t texObject, float x, int layer, float level) 
# 409
{int volatile ___ = 1;(void)texObject;(void)x;(void)layer;(void)level;
# 413
::exit(___);}
#if 0
# 409
{ 
# 410
T ret; 
# 411
tex1DLayeredLod(&ret, texObject, x, layer, level); 
# 412
return ret; 
# 413
} 
#endif
# 416 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/texture_indirect_functions.h"
template< class T> __attribute__((unused)) static typename __nv_itex_trait< T> ::type 
# 417
tex2DLayeredLod(T *ptr, cudaTextureObject_t obj, float x, float y, int layer, float level) 
# 418
{int volatile ___ = 1;(void)ptr;(void)obj;(void)x;(void)y;(void)layer;(void)level;
# 420
::exit(___);}
#if 0
# 418
{ 
# 419
__nv_tex_surf_handler("__itex2DLayeredLod", ptr, obj, x, y, layer, level); 
# 420
} 
#endif
# 422 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/texture_indirect_functions.h"
template< class T> __attribute__((unused)) static T 
# 423
tex2DLayeredLod(cudaTextureObject_t texObject, float x, float y, int layer, float level) 
# 424
{int volatile ___ = 1;(void)texObject;(void)x;(void)y;(void)layer;(void)level;
# 428
::exit(___);}
#if 0
# 424
{ 
# 425
T ret; 
# 426
tex2DLayeredLod(&ret, texObject, x, y, layer, level); 
# 427
return ret; 
# 428
} 
#endif
# 431 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/texture_indirect_functions.h"
template< class T> __attribute__((unused)) static typename __nv_itex_trait< T> ::type 
# 432
tex2DLayeredLod(T *ptr, cudaTextureObject_t obj, float x, float y, int layer, float level, bool *isResident) 
# 433
{int volatile ___ = 1;(void)ptr;(void)obj;(void)x;(void)y;(void)layer;(void)level;(void)isResident;
# 437
::exit(___);}
#if 0
# 433
{ 
# 434
unsigned char res; 
# 435
__nv_tex_surf_handler("__itex2DLayeredLod_sparse", ptr, obj, x, y, layer, level, &res); 
# 436
(*isResident) = (res != 0); 
# 437
} 
#endif
# 439 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/texture_indirect_functions.h"
template< class T> __attribute__((unused)) static T 
# 440
tex2DLayeredLod(cudaTextureObject_t texObject, float x, float y, int layer, float level, bool *isResident) 
# 441
{int volatile ___ = 1;(void)texObject;(void)x;(void)y;(void)layer;(void)level;(void)isResident;
# 445
::exit(___);}
#if 0
# 441
{ 
# 442
T ret; 
# 443
tex2DLayeredLod(&ret, texObject, x, y, layer, level, isResident); 
# 444
return ret; 
# 445
} 
#endif
# 448 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/texture_indirect_functions.h"
template< class T> __attribute__((unused)) static typename __nv_itex_trait< T> ::type 
# 449
texCubemapLod(T *ptr, cudaTextureObject_t obj, float x, float y, float z, float level) 
# 450
{int volatile ___ = 1;(void)ptr;(void)obj;(void)x;(void)y;(void)z;(void)level;
# 452
::exit(___);}
#if 0
# 450
{ 
# 451
__nv_tex_surf_handler("__itexCubemapLod", ptr, obj, x, y, z, level); 
# 452
} 
#endif
# 454 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/texture_indirect_functions.h"
template< class T> __attribute__((unused)) static T 
# 455
texCubemapLod(cudaTextureObject_t texObject, float x, float y, float z, float level) 
# 456
{int volatile ___ = 1;(void)texObject;(void)x;(void)y;(void)z;(void)level;
# 460
::exit(___);}
#if 0
# 456
{ 
# 457
T ret; 
# 458
texCubemapLod(&ret, texObject, x, y, z, level); 
# 459
return ret; 
# 460
} 
#endif
# 463 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/texture_indirect_functions.h"
template< class T> __attribute__((unused)) static typename __nv_itex_trait< T> ::type 
# 464
texCubemapGrad(T *ptr, cudaTextureObject_t obj, float x, float y, float z, float4 dPdx, float4 dPdy) 
# 465
{int volatile ___ = 1;(void)ptr;(void)obj;(void)x;(void)y;(void)z;(void)dPdx;(void)dPdy;
# 467
::exit(___);}
#if 0
# 465
{ 
# 466
__nv_tex_surf_handler("__itexCubemapGrad_v2", ptr, obj, x, y, z, &dPdx, &dPdy); 
# 467
} 
#endif
# 469 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/texture_indirect_functions.h"
template< class T> __attribute__((unused)) static T 
# 470
texCubemapGrad(cudaTextureObject_t texObject, float x, float y, float z, float4 dPdx, float4 dPdy) 
# 471
{int volatile ___ = 1;(void)texObject;(void)x;(void)y;(void)z;(void)dPdx;(void)dPdy;
# 475
::exit(___);}
#if 0
# 471
{ 
# 472
T ret; 
# 473
texCubemapGrad(&ret, texObject, x, y, z, dPdx, dPdy); 
# 474
return ret; 
# 475
} 
#endif
# 477 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/texture_indirect_functions.h"
template< class T> __attribute__((unused)) static typename __nv_itex_trait< T> ::type 
# 478
texCubemapLayeredLod(T *ptr, cudaTextureObject_t obj, float x, float y, float z, int layer, float level) 
# 479
{int volatile ___ = 1;(void)ptr;(void)obj;(void)x;(void)y;(void)z;(void)layer;(void)level;
# 481
::exit(___);}
#if 0
# 479
{ 
# 480
__nv_tex_surf_handler("__itexCubemapLayeredLod", ptr, obj, x, y, z, layer, level); 
# 481
} 
#endif
# 483 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/texture_indirect_functions.h"
template< class T> __attribute__((unused)) static T 
# 484
texCubemapLayeredLod(cudaTextureObject_t texObject, float x, float y, float z, int layer, float level) 
# 485
{int volatile ___ = 1;(void)texObject;(void)x;(void)y;(void)z;(void)layer;(void)level;
# 489
::exit(___);}
#if 0
# 485
{ 
# 486
T ret; 
# 487
texCubemapLayeredLod(&ret, texObject, x, y, z, layer, level); 
# 488
return ret; 
# 489
} 
#endif
# 491 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/texture_indirect_functions.h"
template< class T> __attribute__((unused)) static typename __nv_itex_trait< T> ::type 
# 492
tex1DGrad(T *ptr, cudaTextureObject_t obj, float x, float dPdx, float dPdy) 
# 493
{int volatile ___ = 1;(void)ptr;(void)obj;(void)x;(void)dPdx;(void)dPdy;
# 495
::exit(___);}
#if 0
# 493
{ 
# 494
__nv_tex_surf_handler("__itex1DGrad", ptr, obj, x, dPdx, dPdy); 
# 495
} 
#endif
# 497 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/texture_indirect_functions.h"
template< class T> __attribute__((unused)) static T 
# 498
tex1DGrad(cudaTextureObject_t texObject, float x, float dPdx, float dPdy) 
# 499
{int volatile ___ = 1;(void)texObject;(void)x;(void)dPdx;(void)dPdy;
# 503
::exit(___);}
#if 0
# 499
{ 
# 500
T ret; 
# 501
tex1DGrad(&ret, texObject, x, dPdx, dPdy); 
# 502
return ret; 
# 503
} 
#endif
# 506 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/texture_indirect_functions.h"
template< class T> __attribute__((unused)) static typename __nv_itex_trait< T> ::type 
# 507
tex2DGrad(T *ptr, cudaTextureObject_t obj, float x, float y, float2 dPdx, float2 dPdy) 
# 508
{int volatile ___ = 1;(void)ptr;(void)obj;(void)x;(void)y;(void)dPdx;(void)dPdy;
# 510
::exit(___);}
#if 0
# 508
{ 
# 509
__nv_tex_surf_handler("__itex2DGrad_v2", ptr, obj, x, y, &dPdx, &dPdy); 
# 510
} 
#endif
# 512 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/texture_indirect_functions.h"
template< class T> __attribute__((unused)) static T 
# 513
tex2DGrad(cudaTextureObject_t texObject, float x, float y, float2 dPdx, float2 dPdy) 
# 514
{int volatile ___ = 1;(void)texObject;(void)x;(void)y;(void)dPdx;(void)dPdy;
# 518
::exit(___);}
#if 0
# 514
{ 
# 515
T ret; 
# 516
tex2DGrad(&ret, texObject, x, y, dPdx, dPdy); 
# 517
return ret; 
# 518
} 
#endif
# 521 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/texture_indirect_functions.h"
template< class T> __attribute__((unused)) static typename __nv_itex_trait< T> ::type 
# 522
tex2DGrad(T *ptr, cudaTextureObject_t obj, float x, float y, float2 dPdx, float2 dPdy, bool *isResident) 
# 523
{int volatile ___ = 1;(void)ptr;(void)obj;(void)x;(void)y;(void)dPdx;(void)dPdy;(void)isResident;
# 527
::exit(___);}
#if 0
# 523
{ 
# 524
unsigned char res; 
# 525
__nv_tex_surf_handler("__itex2DGrad_sparse", ptr, obj, x, y, &dPdx, &dPdy, &res); 
# 526
(*isResident) = (res != 0); 
# 527
} 
#endif
# 529 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/texture_indirect_functions.h"
template< class T> __attribute__((unused)) static T 
# 530
tex2DGrad(cudaTextureObject_t texObject, float x, float y, float2 dPdx, float2 dPdy, bool *isResident) 
# 531
{int volatile ___ = 1;(void)texObject;(void)x;(void)y;(void)dPdx;(void)dPdy;(void)isResident;
# 535
::exit(___);}
#if 0
# 531
{ 
# 532
T ret; 
# 533
tex2DGrad(&ret, texObject, x, y, dPdx, dPdy, isResident); 
# 534
return ret; 
# 535
} 
#endif
# 539 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/texture_indirect_functions.h"
template< class T> __attribute__((unused)) static typename __nv_itex_trait< T> ::type 
# 540
tex3DGrad(T *ptr, cudaTextureObject_t obj, float x, float y, float z, float4 dPdx, float4 dPdy) 
# 541
{int volatile ___ = 1;(void)ptr;(void)obj;(void)x;(void)y;(void)z;(void)dPdx;(void)dPdy;
# 543
::exit(___);}
#if 0
# 541
{ 
# 542
__nv_tex_surf_handler("__itex3DGrad_v2", ptr, obj, x, y, z, &dPdx, &dPdy); 
# 543
} 
#endif
# 545 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/texture_indirect_functions.h"
template< class T> __attribute__((unused)) static T 
# 546
tex3DGrad(cudaTextureObject_t texObject, float x, float y, float z, float4 dPdx, float4 dPdy) 
# 547
{int volatile ___ = 1;(void)texObject;(void)x;(void)y;(void)z;(void)dPdx;(void)dPdy;
# 551
::exit(___);}
#if 0
# 547
{ 
# 548
T ret; 
# 549
tex3DGrad(&ret, texObject, x, y, z, dPdx, dPdy); 
# 550
return ret; 
# 551
} 
#endif
# 554 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/texture_indirect_functions.h"
template< class T> __attribute__((unused)) static typename __nv_itex_trait< T> ::type 
# 555
tex3DGrad(T *ptr, cudaTextureObject_t obj, float x, float y, float z, float4 dPdx, float4 dPdy, bool *isResident) 
# 556
{int volatile ___ = 1;(void)ptr;(void)obj;(void)x;(void)y;(void)z;(void)dPdx;(void)dPdy;(void)isResident;
# 560
::exit(___);}
#if 0
# 556
{ 
# 557
unsigned char res; 
# 558
__nv_tex_surf_handler("__itex3DGrad_sparse", ptr, obj, x, y, z, &dPdx, &dPdy, &res); 
# 559
(*isResident) = (res != 0); 
# 560
} 
#endif
# 562 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/texture_indirect_functions.h"
template< class T> __attribute__((unused)) static T 
# 563
tex3DGrad(cudaTextureObject_t texObject, float x, float y, float z, float4 dPdx, float4 dPdy, bool *isResident) 
# 564
{int volatile ___ = 1;(void)texObject;(void)x;(void)y;(void)z;(void)dPdx;(void)dPdy;(void)isResident;
# 568
::exit(___);}
#if 0
# 564
{ 
# 565
T ret; 
# 566
tex3DGrad(&ret, texObject, x, y, z, dPdx, dPdy, isResident); 
# 567
return ret; 
# 568
} 
#endif
# 573 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/texture_indirect_functions.h"
template< class T> __attribute__((unused)) static typename __nv_itex_trait< T> ::type 
# 574
tex1DLayeredGrad(T *ptr, cudaTextureObject_t obj, float x, int layer, float dPdx, float dPdy) 
# 575
{int volatile ___ = 1;(void)ptr;(void)obj;(void)x;(void)layer;(void)dPdx;(void)dPdy;
# 577
::exit(___);}
#if 0
# 575
{ 
# 576
__nv_tex_surf_handler("__itex1DLayeredGrad", ptr, obj, x, layer, dPdx, dPdy); 
# 577
} 
#endif
# 579 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/texture_indirect_functions.h"
template< class T> __attribute__((unused)) static T 
# 580
tex1DLayeredGrad(cudaTextureObject_t texObject, float x, int layer, float dPdx, float dPdy) 
# 581
{int volatile ___ = 1;(void)texObject;(void)x;(void)layer;(void)dPdx;(void)dPdy;
# 585
::exit(___);}
#if 0
# 581
{ 
# 582
T ret; 
# 583
tex1DLayeredGrad(&ret, texObject, x, layer, dPdx, dPdy); 
# 584
return ret; 
# 585
} 
#endif
# 588 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/texture_indirect_functions.h"
template< class T> __attribute__((unused)) static typename __nv_itex_trait< T> ::type 
# 589
tex2DLayeredGrad(T *ptr, cudaTextureObject_t obj, float x, float y, int layer, float2 dPdx, float2 dPdy) 
# 590
{int volatile ___ = 1;(void)ptr;(void)obj;(void)x;(void)y;(void)layer;(void)dPdx;(void)dPdy;
# 592
::exit(___);}
#if 0
# 590
{ 
# 591
__nv_tex_surf_handler("__itex2DLayeredGrad_v2", ptr, obj, x, y, layer, &dPdx, &dPdy); 
# 592
} 
#endif
# 594 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/texture_indirect_functions.h"
template< class T> __attribute__((unused)) static T 
# 595
tex2DLayeredGrad(cudaTextureObject_t texObject, float x, float y, int layer, float2 dPdx, float2 dPdy) 
# 596
{int volatile ___ = 1;(void)texObject;(void)x;(void)y;(void)layer;(void)dPdx;(void)dPdy;
# 600
::exit(___);}
#if 0
# 596
{ 
# 597
T ret; 
# 598
tex2DLayeredGrad(&ret, texObject, x, y, layer, dPdx, dPdy); 
# 599
return ret; 
# 600
} 
#endif
# 603 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/texture_indirect_functions.h"
template< class T> __attribute__((unused)) static typename __nv_itex_trait< T> ::type 
# 604
tex2DLayeredGrad(T *ptr, cudaTextureObject_t obj, float x, float y, int layer, float2 dPdx, float2 dPdy, bool *isResident) 
# 605
{int volatile ___ = 1;(void)ptr;(void)obj;(void)x;(void)y;(void)layer;(void)dPdx;(void)dPdy;(void)isResident;
# 609
::exit(___);}
#if 0
# 605
{ 
# 606
unsigned char res; 
# 607
__nv_tex_surf_handler("__itex2DLayeredGrad_sparse", ptr, obj, x, y, layer, &dPdx, &dPdy, &res); 
# 608
(*isResident) = (res != 0); 
# 609
} 
#endif
# 611 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/texture_indirect_functions.h"
template< class T> __attribute__((unused)) static T 
# 612
tex2DLayeredGrad(cudaTextureObject_t texObject, float x, float y, int layer, float2 dPdx, float2 dPdy, bool *isResident) 
# 613
{int volatile ___ = 1;(void)texObject;(void)x;(void)y;(void)layer;(void)dPdx;(void)dPdy;(void)isResident;
# 617
::exit(___);}
#if 0
# 613
{ 
# 614
T ret; 
# 615
tex2DLayeredGrad(&ret, texObject, x, y, layer, dPdx, dPdy, isResident); 
# 616
return ret; 
# 617
} 
#endif
# 621 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/texture_indirect_functions.h"
template< class T> __attribute__((unused)) static typename __nv_itex_trait< T> ::type 
# 622
texCubemapLayeredGrad(T *ptr, cudaTextureObject_t obj, float x, float y, float z, int layer, float4 dPdx, float4 dPdy) 
# 623
{int volatile ___ = 1;(void)ptr;(void)obj;(void)x;(void)y;(void)z;(void)layer;(void)dPdx;(void)dPdy;
# 625
::exit(___);}
#if 0
# 623
{ 
# 624
__nv_tex_surf_handler("__itexCubemapLayeredGrad_v2", ptr, obj, x, y, z, layer, &dPdx, &dPdy); 
# 625
} 
#endif
# 627 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/texture_indirect_functions.h"
template< class T> __attribute__((unused)) static T 
# 628
texCubemapLayeredGrad(cudaTextureObject_t texObject, float x, float y, float z, int layer, float4 dPdx, float4 dPdy) 
# 629
{int volatile ___ = 1;(void)texObject;(void)x;(void)y;(void)z;(void)layer;(void)dPdx;(void)dPdy;
# 633
::exit(___);}
#if 0
# 629
{ 
# 630
T ret; 
# 631
texCubemapLayeredGrad(&ret, texObject, x, y, z, layer, dPdx, dPdy); 
# 632
return ret; 
# 633
} 
#endif
# 58 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/surface_indirect_functions.h"
template< class T> struct __nv_isurf_trait { }; 
# 59
template<> struct __nv_isurf_trait< char>  { typedef void type; }; 
# 60
template<> struct __nv_isurf_trait< signed char>  { typedef void type; }; 
# 61
template<> struct __nv_isurf_trait< char1>  { typedef void type; }; 
# 62
template<> struct __nv_isurf_trait< unsigned char>  { typedef void type; }; 
# 63
template<> struct __nv_isurf_trait< uchar1>  { typedef void type; }; 
# 64
template<> struct __nv_isurf_trait< short>  { typedef void type; }; 
# 65
template<> struct __nv_isurf_trait< short1>  { typedef void type; }; 
# 66
template<> struct __nv_isurf_trait< unsigned short>  { typedef void type; }; 
# 67
template<> struct __nv_isurf_trait< ushort1>  { typedef void type; }; 
# 68
template<> struct __nv_isurf_trait< int>  { typedef void type; }; 
# 69
template<> struct __nv_isurf_trait< int1>  { typedef void type; }; 
# 70
template<> struct __nv_isurf_trait< unsigned>  { typedef void type; }; 
# 71
template<> struct __nv_isurf_trait< uint1>  { typedef void type; }; 
# 72
template<> struct __nv_isurf_trait< long long>  { typedef void type; }; 
# 73
template<> struct __nv_isurf_trait< longlong1>  { typedef void type; }; 
# 74
template<> struct __nv_isurf_trait< unsigned long long>  { typedef void type; }; 
# 75
template<> struct __nv_isurf_trait< ulonglong1>  { typedef void type; }; 
# 76
template<> struct __nv_isurf_trait< float>  { typedef void type; }; 
# 77
template<> struct __nv_isurf_trait< float1>  { typedef void type; }; 
# 79
template<> struct __nv_isurf_trait< char2>  { typedef void type; }; 
# 80
template<> struct __nv_isurf_trait< uchar2>  { typedef void type; }; 
# 81
template<> struct __nv_isurf_trait< short2>  { typedef void type; }; 
# 82
template<> struct __nv_isurf_trait< ushort2>  { typedef void type; }; 
# 83
template<> struct __nv_isurf_trait< int2>  { typedef void type; }; 
# 84
template<> struct __nv_isurf_trait< uint2>  { typedef void type; }; 
# 85
template<> struct __nv_isurf_trait< longlong2>  { typedef void type; }; 
# 86
template<> struct __nv_isurf_trait< ulonglong2>  { typedef void type; }; 
# 87
template<> struct __nv_isurf_trait< float2>  { typedef void type; }; 
# 89
template<> struct __nv_isurf_trait< char4>  { typedef void type; }; 
# 90
template<> struct __nv_isurf_trait< uchar4>  { typedef void type; }; 
# 91
template<> struct __nv_isurf_trait< short4>  { typedef void type; }; 
# 92
template<> struct __nv_isurf_trait< ushort4>  { typedef void type; }; 
# 93
template<> struct __nv_isurf_trait< int4>  { typedef void type; }; 
# 94
template<> struct __nv_isurf_trait< uint4>  { typedef void type; }; 
# 95
template<> struct __nv_isurf_trait< float4>  { typedef void type; }; 
# 98
template< class T> __attribute__((unused)) static typename __nv_isurf_trait< T> ::type 
# 99
surf1Dread(T *ptr, cudaSurfaceObject_t obj, int x, cudaSurfaceBoundaryMode mode = cudaBoundaryModeTrap) 
# 100
{int volatile ___ = 1;(void)ptr;(void)obj;(void)x;(void)mode;
# 102
::exit(___);}
#if 0
# 100
{ 
# 101
__nv_tex_surf_handler("__isurf1Dread", ptr, obj, x, mode); 
# 102
} 
#endif
# 104 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/surface_indirect_functions.h"
template< class T> __attribute__((unused)) static T 
# 105
surf1Dread(cudaSurfaceObject_t surfObject, int x, cudaSurfaceBoundaryMode boundaryMode = cudaBoundaryModeTrap) 
# 106
{int volatile ___ = 1;(void)surfObject;(void)x;(void)boundaryMode;
# 110
::exit(___);}
#if 0
# 106
{ 
# 107
T ret; 
# 108
surf1Dread(&ret, surfObject, x, boundaryMode); 
# 109
return ret; 
# 110
} 
#endif
# 112 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/surface_indirect_functions.h"
template< class T> __attribute__((unused)) static typename __nv_isurf_trait< T> ::type 
# 113
surf2Dread(T *ptr, cudaSurfaceObject_t obj, int x, int y, cudaSurfaceBoundaryMode mode = cudaBoundaryModeTrap) 
# 114
{int volatile ___ = 1;(void)ptr;(void)obj;(void)x;(void)y;(void)mode;
# 116
::exit(___);}
#if 0
# 114
{ 
# 115
__nv_tex_surf_handler("__isurf2Dread", ptr, obj, x, y, mode); 
# 116
} 
#endif
# 118 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/surface_indirect_functions.h"
template< class T> __attribute__((unused)) static T 
# 119
surf2Dread(cudaSurfaceObject_t surfObject, int x, int y, cudaSurfaceBoundaryMode boundaryMode = cudaBoundaryModeTrap) 
# 120
{int volatile ___ = 1;(void)surfObject;(void)x;(void)y;(void)boundaryMode;
# 124
::exit(___);}
#if 0
# 120
{ 
# 121
T ret; 
# 122
surf2Dread(&ret, surfObject, x, y, boundaryMode); 
# 123
return ret; 
# 124
} 
#endif
# 127 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/surface_indirect_functions.h"
template< class T> __attribute__((unused)) static typename __nv_isurf_trait< T> ::type 
# 128
surf3Dread(T *ptr, cudaSurfaceObject_t obj, int x, int y, int z, cudaSurfaceBoundaryMode mode = cudaBoundaryModeTrap) 
# 129
{int volatile ___ = 1;(void)ptr;(void)obj;(void)x;(void)y;(void)z;(void)mode;
# 131
::exit(___);}
#if 0
# 129
{ 
# 130
__nv_tex_surf_handler("__isurf3Dread", ptr, obj, x, y, z, mode); 
# 131
} 
#endif
# 133 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/surface_indirect_functions.h"
template< class T> __attribute__((unused)) static T 
# 134
surf3Dread(cudaSurfaceObject_t surfObject, int x, int y, int z, cudaSurfaceBoundaryMode boundaryMode = cudaBoundaryModeTrap) 
# 135
{int volatile ___ = 1;(void)surfObject;(void)x;(void)y;(void)z;(void)boundaryMode;
# 139
::exit(___);}
#if 0
# 135
{ 
# 136
T ret; 
# 137
surf3Dread(&ret, surfObject, x, y, z, boundaryMode); 
# 138
return ret; 
# 139
} 
#endif
# 141 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/surface_indirect_functions.h"
template< class T> __attribute__((unused)) static typename __nv_isurf_trait< T> ::type 
# 142
surf1DLayeredread(T *ptr, cudaSurfaceObject_t obj, int x, int layer, cudaSurfaceBoundaryMode mode = cudaBoundaryModeTrap) 
# 143
{int volatile ___ = 1;(void)ptr;(void)obj;(void)x;(void)layer;(void)mode;
# 145
::exit(___);}
#if 0
# 143
{ 
# 144
__nv_tex_surf_handler("__isurf1DLayeredread", ptr, obj, x, layer, mode); 
# 145
} 
#endif
# 147 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/surface_indirect_functions.h"
template< class T> __attribute__((unused)) static T 
# 148
surf1DLayeredread(cudaSurfaceObject_t surfObject, int x, int layer, cudaSurfaceBoundaryMode boundaryMode = cudaBoundaryModeTrap) 
# 149
{int volatile ___ = 1;(void)surfObject;(void)x;(void)layer;(void)boundaryMode;
# 153
::exit(___);}
#if 0
# 149
{ 
# 150
T ret; 
# 151
surf1DLayeredread(&ret, surfObject, x, layer, boundaryMode); 
# 152
return ret; 
# 153
} 
#endif
# 155 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/surface_indirect_functions.h"
template< class T> __attribute__((unused)) static typename __nv_isurf_trait< T> ::type 
# 156
surf2DLayeredread(T *ptr, cudaSurfaceObject_t obj, int x, int y, int layer, cudaSurfaceBoundaryMode mode = cudaBoundaryModeTrap) 
# 157
{int volatile ___ = 1;(void)ptr;(void)obj;(void)x;(void)y;(void)layer;(void)mode;
# 159
::exit(___);}
#if 0
# 157
{ 
# 158
__nv_tex_surf_handler("__isurf2DLayeredread", ptr, obj, x, y, layer, mode); 
# 159
} 
#endif
# 161 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/surface_indirect_functions.h"
template< class T> __attribute__((unused)) static T 
# 162
surf2DLayeredread(cudaSurfaceObject_t surfObject, int x, int y, int layer, cudaSurfaceBoundaryMode boundaryMode = cudaBoundaryModeTrap) 
# 163
{int volatile ___ = 1;(void)surfObject;(void)x;(void)y;(void)layer;(void)boundaryMode;
# 167
::exit(___);}
#if 0
# 163
{ 
# 164
T ret; 
# 165
surf2DLayeredread(&ret, surfObject, x, y, layer, boundaryMode); 
# 166
return ret; 
# 167
} 
#endif
# 169 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/surface_indirect_functions.h"
template< class T> __attribute__((unused)) static typename __nv_isurf_trait< T> ::type 
# 170
surfCubemapread(T *ptr, cudaSurfaceObject_t obj, int x, int y, int face, cudaSurfaceBoundaryMode mode = cudaBoundaryModeTrap) 
# 171
{int volatile ___ = 1;(void)ptr;(void)obj;(void)x;(void)y;(void)face;(void)mode;
# 173
::exit(___);}
#if 0
# 171
{ 
# 172
__nv_tex_surf_handler("__isurfCubemapread", ptr, obj, x, y, face, mode); 
# 173
} 
#endif
# 175 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/surface_indirect_functions.h"
template< class T> __attribute__((unused)) static T 
# 176
surfCubemapread(cudaSurfaceObject_t surfObject, int x, int y, int face, cudaSurfaceBoundaryMode boundaryMode = cudaBoundaryModeTrap) 
# 177
{int volatile ___ = 1;(void)surfObject;(void)x;(void)y;(void)face;(void)boundaryMode;
# 181
::exit(___);}
#if 0
# 177
{ 
# 178
T ret; 
# 179
surfCubemapread(&ret, surfObject, x, y, face, boundaryMode); 
# 180
return ret; 
# 181
} 
#endif
# 183 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/surface_indirect_functions.h"
template< class T> __attribute__((unused)) static typename __nv_isurf_trait< T> ::type 
# 184
surfCubemapLayeredread(T *ptr, cudaSurfaceObject_t obj, int x, int y, int layerface, cudaSurfaceBoundaryMode mode = cudaBoundaryModeTrap) 
# 185
{int volatile ___ = 1;(void)ptr;(void)obj;(void)x;(void)y;(void)layerface;(void)mode;
# 187
::exit(___);}
#if 0
# 185
{ 
# 186
__nv_tex_surf_handler("__isurfCubemapLayeredread", ptr, obj, x, y, layerface, mode); 
# 187
} 
#endif
# 189 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/surface_indirect_functions.h"
template< class T> __attribute__((unused)) static T 
# 190
surfCubemapLayeredread(cudaSurfaceObject_t surfObject, int x, int y, int layerface, cudaSurfaceBoundaryMode boundaryMode = cudaBoundaryModeTrap) 
# 191
{int volatile ___ = 1;(void)surfObject;(void)x;(void)y;(void)layerface;(void)boundaryMode;
# 195
::exit(___);}
#if 0
# 191
{ 
# 192
T ret; 
# 193
surfCubemapLayeredread(&ret, surfObject, x, y, layerface, boundaryMode); 
# 194
return ret; 
# 195
} 
#endif
# 197 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/surface_indirect_functions.h"
template< class T> __attribute__((unused)) static typename __nv_isurf_trait< T> ::type 
# 198
surf1Dwrite(T val, cudaSurfaceObject_t obj, int x, cudaSurfaceBoundaryMode mode = cudaBoundaryModeTrap) 
# 199
{int volatile ___ = 1;(void)val;(void)obj;(void)x;(void)mode;
# 201
::exit(___);}
#if 0
# 199
{ 
# 200
__nv_tex_surf_handler("__isurf1Dwrite_v2", &val, obj, x, mode); 
# 201
} 
#endif
# 203 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/surface_indirect_functions.h"
template< class T> __attribute__((unused)) static typename __nv_isurf_trait< T> ::type 
# 204
surf2Dwrite(T val, cudaSurfaceObject_t obj, int x, int y, cudaSurfaceBoundaryMode mode = cudaBoundaryModeTrap) 
# 205
{int volatile ___ = 1;(void)val;(void)obj;(void)x;(void)y;(void)mode;
# 207
::exit(___);}
#if 0
# 205
{ 
# 206
__nv_tex_surf_handler("__isurf2Dwrite_v2", &val, obj, x, y, mode); 
# 207
} 
#endif
# 209 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/surface_indirect_functions.h"
template< class T> __attribute__((unused)) static typename __nv_isurf_trait< T> ::type 
# 210
surf3Dwrite(T val, cudaSurfaceObject_t obj, int x, int y, int z, cudaSurfaceBoundaryMode mode = cudaBoundaryModeTrap) 
# 211
{int volatile ___ = 1;(void)val;(void)obj;(void)x;(void)y;(void)z;(void)mode;
# 213
::exit(___);}
#if 0
# 211
{ 
# 212
__nv_tex_surf_handler("__isurf3Dwrite_v2", &val, obj, x, y, z, mode); 
# 213
} 
#endif
# 215 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/surface_indirect_functions.h"
template< class T> __attribute__((unused)) static typename __nv_isurf_trait< T> ::type 
# 216
surf1DLayeredwrite(T val, cudaSurfaceObject_t obj, int x, int layer, cudaSurfaceBoundaryMode mode = cudaBoundaryModeTrap) 
# 217
{int volatile ___ = 1;(void)val;(void)obj;(void)x;(void)layer;(void)mode;
# 219
::exit(___);}
#if 0
# 217
{ 
# 218
__nv_tex_surf_handler("__isurf1DLayeredwrite_v2", &val, obj, x, layer, mode); 
# 219
} 
#endif
# 221 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/surface_indirect_functions.h"
template< class T> __attribute__((unused)) static typename __nv_isurf_trait< T> ::type 
# 222
surf2DLayeredwrite(T val, cudaSurfaceObject_t obj, int x, int y, int layer, cudaSurfaceBoundaryMode mode = cudaBoundaryModeTrap) 
# 223
{int volatile ___ = 1;(void)val;(void)obj;(void)x;(void)y;(void)layer;(void)mode;
# 225
::exit(___);}
#if 0
# 223
{ 
# 224
__nv_tex_surf_handler("__isurf2DLayeredwrite_v2", &val, obj, x, y, layer, mode); 
# 225
} 
#endif
# 227 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/surface_indirect_functions.h"
template< class T> __attribute__((unused)) static typename __nv_isurf_trait< T> ::type 
# 228
surfCubemapwrite(T val, cudaSurfaceObject_t obj, int x, int y, int face, cudaSurfaceBoundaryMode mode = cudaBoundaryModeTrap) 
# 229
{int volatile ___ = 1;(void)val;(void)obj;(void)x;(void)y;(void)face;(void)mode;
# 231
::exit(___);}
#if 0
# 229
{ 
# 230
__nv_tex_surf_handler("__isurfCubemapwrite_v2", &val, obj, x, y, face, mode); 
# 231
} 
#endif
# 233 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/surface_indirect_functions.h"
template< class T> __attribute__((unused)) static typename __nv_isurf_trait< T> ::type 
# 234
surfCubemapLayeredwrite(T val, cudaSurfaceObject_t obj, int x, int y, int layerface, cudaSurfaceBoundaryMode mode = cudaBoundaryModeTrap) 
# 235
{int volatile ___ = 1;(void)val;(void)obj;(void)x;(void)y;(void)layerface;(void)mode;
# 237
::exit(___);}
#if 0
# 235
{ 
# 236
__nv_tex_surf_handler("__isurfCubemapLayeredwrite_v2", &val, obj, x, y, layerface, mode); 
# 237
} 
#endif
# 3634 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/crt/device_functions.h"
extern "C" unsigned __cudaPushCallConfiguration(dim3 gridDim, dim3 blockDim, size_t sharedMem = 0, CUstream_st * stream = 0); 
# 68 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/device_launch_parameters.h"
extern "C" {
# 71
extern const uint3 __device_builtin_variable_threadIdx; 
# 72
extern const uint3 __device_builtin_variable_blockIdx; 
# 73
extern const dim3 __device_builtin_variable_blockDim; 
# 74
extern const dim3 __device_builtin_variable_gridDim; 
# 75
extern const int __device_builtin_variable_warpSize; 
# 80
}
# 67 "/usr/include/c++/10/bits/stl_relops.h" 3
namespace std __attribute((__visibility__("default"))) { 
# 71
namespace rel_ops { 
# 85 "/usr/include/c++/10/bits/stl_relops.h" 3
template< class _Tp> inline bool 
# 87
operator!=(const _Tp &__x, const _Tp &__y) 
# 88
{ return !(__x == __y); } 
# 98 "/usr/include/c++/10/bits/stl_relops.h" 3
template< class _Tp> inline bool 
# 100
operator>(const _Tp &__x, const _Tp &__y) 
# 101
{ return __y < __x; } 
# 111 "/usr/include/c++/10/bits/stl_relops.h" 3
template< class _Tp> inline bool 
# 113
operator<=(const _Tp &__x, const _Tp &__y) 
# 114
{ return !(__y < __x); } 
# 124 "/usr/include/c++/10/bits/stl_relops.h" 3
template< class _Tp> inline bool 
# 126
operator>=(const _Tp &__x, const _Tp &__y) 
# 127
{ return !(__x < __y); } 
# 128
}
# 131
}
# 38 "/usr/include/c++/10/bits/move.h" 3
namespace std __attribute((__visibility__("default"))) { 
# 47
template< class _Tp> constexpr _Tp *
# 49
__addressof(_Tp &__r) noexcept 
# 50
{ return __builtin_addressof(__r); } 
# 55
}
# 40 "/usr/include/c++/10/type_traits" 3
namespace std __attribute((__visibility__("default"))) { 
# 56 "/usr/include/c++/10/type_traits" 3
template< class _Tp, _Tp __v> 
# 57
struct integral_constant { 
# 59
static constexpr _Tp value = (__v); 
# 60
typedef _Tp value_type; 
# 61
typedef integral_constant type; 
# 62
constexpr operator value_type() const noexcept { return value; } 
# 67
constexpr value_type operator()() const noexcept { return value; } 
# 69
}; 
# 71
template< class _Tp, _Tp __v> constexpr _Tp integral_constant< _Tp, __v> ::value; 
# 75
typedef integral_constant< bool, true>  true_type; 
# 78
typedef integral_constant< bool, false>  false_type; 
# 80
template< bool __v> using __bool_constant = integral_constant< bool, __v> ; 
# 91 "/usr/include/c++/10/type_traits" 3
template< bool , class , class > struct conditional; 
# 94
template< class _Type> 
# 95
struct __type_identity { 
# 96
using type = _Type; }; 
# 98
template< class _Tp> using __type_identity_t = typename __type_identity< _Tp> ::type; 
# 101
template< class ...> struct __or_; 
# 105
template<> struct __or_< >  : public false_type { 
# 107
}; 
# 109
template< class _B1> 
# 110
struct __or_< _B1>  : public _B1 { 
# 112
}; 
# 114
template< class _B1, class _B2> 
# 115
struct __or_< _B1, _B2>  : public conditional< _B1::value, _B1, _B2> ::type { 
# 117
}; 
# 119
template< class _B1, class _B2, class _B3, class ..._Bn> 
# 120
struct __or_< _B1, _B2, _B3, _Bn...>  : public conditional< _B1::value, _B1, std::__or_< _B2, _B3, _Bn...> > ::type { 
# 122
}; 
# 124
template< class ...> struct __and_; 
# 128
template<> struct __and_< >  : public true_type { 
# 130
}; 
# 132
template< class _B1> 
# 133
struct __and_< _B1>  : public _B1 { 
# 135
}; 
# 137
template< class _B1, class _B2> 
# 138
struct __and_< _B1, _B2>  : public conditional< _B1::value, _B2, _B1> ::type { 
# 140
}; 
# 142
template< class _B1, class _B2, class _B3, class ..._Bn> 
# 143
struct __and_< _B1, _B2, _B3, _Bn...>  : public conditional< _B1::value, std::__and_< _B2, _B3, _Bn...> , _B1> ::type { 
# 145
}; 
# 147
template< class _Pp> 
# 148
struct __not_ : public __bool_constant< !((bool)_Pp::value)>  { 
# 150
}; 
# 188 "/usr/include/c++/10/type_traits" 3
template< class > struct is_reference; 
# 190
template< class > struct is_function; 
# 192
template< class > struct is_void; 
# 194
template< class > struct __is_array_unknown_bounds; 
# 200
template< class _Tp, size_t  = sizeof(_Tp)> constexpr true_type 
# 201
__is_complete_or_unbounded(__type_identity< _Tp> ) 
# 202
{ return {}; } 
# 204
template< class _TypeIdentity, class 
# 205
_NestedType = typename _TypeIdentity::type> constexpr typename __or_< is_reference< _NestedType> , is_function< _NestedType> , is_void< _NestedType> , __is_array_unknown_bounds< _NestedType> > ::type 
# 211
__is_complete_or_unbounded(_TypeIdentity) 
# 212
{ return {}; } 
# 219
template< class _Tp> 
# 220
struct __success_type { 
# 221
typedef _Tp type; }; 
# 223
struct __failure_type { 
# 224
}; 
# 226
template< class > struct remove_cv; 
# 230
template< class _Tp> using __remove_cv_t = typename remove_cv< _Tp> ::type; 
# 233
template< class > struct is_const; 
# 238
template< class > 
# 239
struct __is_void_helper : public false_type { 
# 240
}; 
# 243
template<> struct __is_void_helper< void>  : public true_type { 
# 244
}; 
# 247
template< class _Tp> 
# 248
struct is_void : public __is_void_helper< __remove_cv_t< _Tp> > ::type { 
# 250
}; 
# 252
template< class > 
# 253
struct __is_integral_helper : public false_type { 
# 254
}; 
# 257
template<> struct __is_integral_helper< bool>  : public true_type { 
# 258
}; 
# 261
template<> struct __is_integral_helper< char>  : public true_type { 
# 262
}; 
# 265
template<> struct __is_integral_helper< signed char>  : public true_type { 
# 266
}; 
# 269
template<> struct __is_integral_helper< unsigned char>  : public true_type { 
# 270
}; 
# 274
template<> struct __is_integral_helper< wchar_t>  : public true_type { 
# 275
}; 
# 285 "/usr/include/c++/10/type_traits" 3
template<> struct __is_integral_helper< char16_t>  : public true_type { 
# 286
}; 
# 289
template<> struct __is_integral_helper< char32_t>  : public true_type { 
# 290
}; 
# 293
template<> struct __is_integral_helper< short>  : public true_type { 
# 294
}; 
# 297
template<> struct __is_integral_helper< unsigned short>  : public true_type { 
# 298
}; 
# 301
template<> struct __is_integral_helper< int>  : public true_type { 
# 302
}; 
# 305
template<> struct __is_integral_helper< unsigned>  : public true_type { 
# 306
}; 
# 309
template<> struct __is_integral_helper< long>  : public true_type { 
# 310
}; 
# 313
template<> struct __is_integral_helper< unsigned long>  : public true_type { 
# 314
}; 
# 317
template<> struct __is_integral_helper< long long>  : public true_type { 
# 318
}; 
# 321
template<> struct __is_integral_helper< unsigned long long>  : public true_type { 
# 322
}; 
# 328
template<> struct __is_integral_helper< __int128>  : public true_type { 
# 329
}; 
# 332
template<> struct __is_integral_helper< unsigned __int128>  : public true_type { 
# 333
}; 
# 364 "/usr/include/c++/10/type_traits" 3
template< class _Tp> 
# 365
struct is_integral : public __is_integral_helper< __remove_cv_t< _Tp> > ::type { 
# 367
}; 
# 369
template< class > 
# 370
struct __is_floating_point_helper : public false_type { 
# 371
}; 
# 374
template<> struct __is_floating_point_helper< float>  : public true_type { 
# 375
}; 
# 378
template<> struct __is_floating_point_helper< double>  : public true_type { 
# 379
}; 
# 382
template<> struct __is_floating_point_helper< long double>  : public true_type { 
# 383
}; 
# 392 "/usr/include/c++/10/type_traits" 3
template< class _Tp> 
# 393
struct is_floating_point : public __is_floating_point_helper< __remove_cv_t< _Tp> > ::type { 
# 395
}; 
# 398
template< class > 
# 399
struct is_array : public false_type { 
# 400
}; 
# 402
template< class _Tp, size_t _Size> 
# 403
struct is_array< _Tp [_Size]>  : public true_type { 
# 404
}; 
# 406
template< class _Tp> 
# 407
struct is_array< _Tp []>  : public true_type { 
# 408
}; 
# 410
template< class > 
# 411
struct __is_pointer_helper : public false_type { 
# 412
}; 
# 414
template< class _Tp> 
# 415
struct __is_pointer_helper< _Tp *>  : public true_type { 
# 416
}; 
# 419
template< class _Tp> 
# 420
struct is_pointer : public __is_pointer_helper< __remove_cv_t< _Tp> > ::type { 
# 422
}; 
# 425
template< class > 
# 426
struct is_lvalue_reference : public false_type { 
# 427
}; 
# 429
template< class _Tp> 
# 430
struct is_lvalue_reference< _Tp &>  : public true_type { 
# 431
}; 
# 434
template< class > 
# 435
struct is_rvalue_reference : public false_type { 
# 436
}; 
# 438
template< class _Tp> 
# 439
struct is_rvalue_reference< _Tp &&>  : public true_type { 
# 440
}; 
# 442
template< class > 
# 443
struct __is_member_object_pointer_helper : public false_type { 
# 444
}; 
# 446
template< class _Tp, class _Cp> 
# 447
struct __is_member_object_pointer_helper< _Tp (_Cp::*)>  : public __not_< is_function< _Tp> > ::type { 
# 448
}; 
# 451
template< class _Tp> 
# 452
struct is_member_object_pointer : public __is_member_object_pointer_helper< __remove_cv_t< _Tp> > ::type { 
# 454
}; 
# 456
template< class > 
# 457
struct __is_member_function_pointer_helper : public false_type { 
# 458
}; 
# 460
template< class _Tp, class _Cp> 
# 461
struct __is_member_function_pointer_helper< _Tp (_Cp::*)>  : public is_function< _Tp> ::type { 
# 462
}; 
# 465
template< class _Tp> 
# 466
struct is_member_function_pointer : public __is_member_function_pointer_helper< __remove_cv_t< _Tp> > ::type { 
# 468
}; 
# 471
template< class _Tp> 
# 472
struct is_enum : public integral_constant< bool, __is_enum(_Tp)>  { 
# 474
}; 
# 477
template< class _Tp> 
# 478
struct is_union : public integral_constant< bool, __is_union(_Tp)>  { 
# 480
}; 
# 483
template< class _Tp> 
# 484
struct is_class : public integral_constant< bool, __is_class(_Tp)>  { 
# 486
}; 
# 489
template< class _Tp> 
# 490
struct is_function : public __bool_constant< !is_const< const _Tp> ::value>  { 
# 491
}; 
# 493
template< class _Tp> 
# 494
struct is_function< _Tp &>  : public false_type { 
# 495
}; 
# 497
template< class _Tp> 
# 498
struct is_function< _Tp &&>  : public false_type { 
# 499
}; 
# 503
template< class > 
# 504
struct __is_null_pointer_helper : public false_type { 
# 505
}; 
# 508
template<> struct __is_null_pointer_helper< __decltype((nullptr))>  : public true_type { 
# 509
}; 
# 512
template< class _Tp> 
# 513
struct is_null_pointer : public __is_null_pointer_helper< __remove_cv_t< _Tp> > ::type { 
# 515
}; 
# 518
template< class _Tp> 
# 519
struct __is_nullptr_t : public is_null_pointer< _Tp>  { 
# 521
} __attribute((__deprecated__("use \'std::is_null_pointer\' instead"))); 
# 526
template< class _Tp> 
# 527
struct is_reference : public __or_< is_lvalue_reference< _Tp> , is_rvalue_reference< _Tp> > ::type { 
# 530
}; 
# 533
template< class _Tp> 
# 534
struct is_arithmetic : public __or_< is_integral< _Tp> , is_floating_point< _Tp> > ::type { 
# 536
}; 
# 539
template< class _Tp> 
# 540
struct is_fundamental : public __or_< is_arithmetic< _Tp> , is_void< _Tp> , is_null_pointer< _Tp> > ::type { 
# 543
}; 
# 546
template< class _Tp> 
# 547
struct is_object : public __not_< __or_< is_function< _Tp> , is_reference< _Tp> , is_void< _Tp> > > ::type { 
# 550
}; 
# 552
template< class > struct is_member_pointer; 
# 556
template< class _Tp> 
# 557
struct is_scalar : public __or_< is_arithmetic< _Tp> , is_enum< _Tp> , is_pointer< _Tp> , is_member_pointer< _Tp> , is_null_pointer< _Tp> > ::type { 
# 560
}; 
# 563
template< class _Tp> 
# 564
struct is_compound : public __not_< is_fundamental< _Tp> > ::type { 
# 565
}; 
# 567
template< class _Tp> 
# 568
struct __is_member_pointer_helper : public false_type { 
# 569
}; 
# 571
template< class _Tp, class _Cp> 
# 572
struct __is_member_pointer_helper< _Tp (_Cp::*)>  : public true_type { 
# 573
}; 
# 576
template< class _Tp> 
# 577
struct is_member_pointer : public __is_member_pointer_helper< __remove_cv_t< _Tp> > ::type { 
# 579
}; 
# 581
template< class , class > struct is_same; 
# 584
template< class _Tp, class ..._Types> using __is_one_of = __or_< is_same< _Tp, _Types> ...> ; 
# 588
template< class _Tp> using __is_signed_integer = __is_one_of< __remove_cv_t< _Tp> , signed char, signed short, signed int, signed long, signed long long, signed __int128> ; 
# 607 "/usr/include/c++/10/type_traits" 3
template< class _Tp> using __is_unsigned_integer = __is_one_of< __remove_cv_t< _Tp> , unsigned char, unsigned short, unsigned, unsigned long, unsigned long long, unsigned __int128> ; 
# 626 "/usr/include/c++/10/type_traits" 3
template< class _Tp> using __is_standard_integer = __or_< __is_signed_integer< _Tp> , __is_unsigned_integer< _Tp> > ; 
# 631
template< class ...> using __void_t = void; 
# 635
template< class _Tp, class  = void> 
# 636
struct __is_referenceable : public false_type { 
# 638
}; 
# 640
template< class _Tp> 
# 641
struct __is_referenceable< _Tp, __void_t< _Tp &> >  : public true_type { 
# 643
}; 
# 648
template< class > 
# 649
struct is_const : public false_type { 
# 650
}; 
# 652
template< class _Tp> 
# 653
struct is_const< const _Tp>  : public true_type { 
# 654
}; 
# 657
template< class > 
# 658
struct is_volatile : public false_type { 
# 659
}; 
# 661
template< class _Tp> 
# 662
struct is_volatile< volatile _Tp>  : public true_type { 
# 663
}; 
# 666
template< class _Tp> 
# 667
struct is_trivial : public integral_constant< bool, __is_trivial(_Tp)>  { 
# 670
static_assert((std::__is_complete_or_unbounded(__type_identity< _Tp> {})), "template argument must be a complete class or an unbounded array");
# 672
}; 
# 675
template< class _Tp> 
# 676
struct is_trivially_copyable : public integral_constant< bool, __is_trivially_copyable(_Tp)>  { 
# 679
static_assert((std::__is_complete_or_unbounded(__type_identity< _Tp> {})), "template argument must be a complete class or an unbounded array");
# 681
}; 
# 684
template< class _Tp> 
# 685
struct is_standard_layout : public integral_constant< bool, __is_standard_layout(_Tp)>  { 
# 688
static_assert((std::__is_complete_or_unbounded(__type_identity< _Tp> {})), "template argument must be a complete class or an unbounded array");
# 690
}; 
# 694
template< class _Tp> 
# 697
struct is_pod : public integral_constant< bool, __is_pod(_Tp)>  { 
# 700
static_assert((std::__is_complete_or_unbounded(__type_identity< _Tp> {})), "template argument must be a complete class or an unbounded array");
# 702
}; 
# 705
template< class _Tp> 
# 706
struct is_literal_type : public integral_constant< bool, __is_literal_type(_Tp)>  { 
# 709
static_assert((std::__is_complete_or_unbounded(__type_identity< _Tp> {})), "template argument must be a complete class or an unbounded array");
# 711
}; 
# 714
template< class _Tp> 
# 715
struct is_empty : public integral_constant< bool, __is_empty(_Tp)>  { 
# 717
}; 
# 720
template< class _Tp> 
# 721
struct is_polymorphic : public integral_constant< bool, __is_polymorphic(_Tp)>  { 
# 723
}; 
# 728
template< class _Tp> 
# 729
struct is_final : public integral_constant< bool, __is_final(_Tp)>  { 
# 731
}; 
# 735
template< class _Tp> 
# 736
struct is_abstract : public integral_constant< bool, __is_abstract(_Tp)>  { 
# 738
}; 
# 740
template< class _Tp, bool 
# 741
 = is_arithmetic< _Tp> ::value> 
# 742
struct __is_signed_helper : public false_type { 
# 743
}; 
# 745
template< class _Tp> 
# 746
struct __is_signed_helper< _Tp, true>  : public integral_constant< bool, ((_Tp)(-1)) < ((_Tp)0)>  { 
# 748
}; 
# 751
template< class _Tp> 
# 752
struct is_signed : public __is_signed_helper< _Tp> ::type { 
# 754
}; 
# 757
template< class _Tp> 
# 758
struct is_unsigned : public __and_< is_arithmetic< _Tp> , __not_< is_signed< _Tp> > >  { 
# 760
}; 
# 770 "/usr/include/c++/10/type_traits" 3
template< class _Tp, class _Up = _Tp &&> _Up __declval(int); 
# 774
template< class _Tp> _Tp __declval(long); 
# 778
template< class _Tp> auto declval() noexcept->__decltype((__declval< _Tp> (0))); 
# 781
template< class , unsigned  = 0U> struct extent; 
# 784
template< class > struct remove_all_extents; 
# 787
template< class _Tp> 
# 788
struct __is_array_known_bounds : public integral_constant< bool, (extent< _Tp> ::value > 0)>  { 
# 790
}; 
# 792
template< class _Tp> 
# 793
struct __is_array_unknown_bounds : public __and_< is_array< _Tp> , __not_< extent< _Tp> > >  { 
# 795
}; 
# 802
struct __do_is_destructible_impl { 
# 804
template< class _Tp, class  = __decltype((declval< _Tp &> ().~_Tp()))> static true_type __test(int); 
# 807
template< class > static false_type __test(...); 
# 809
}; 
# 811
template< class _Tp> 
# 812
struct __is_destructible_impl : public __do_is_destructible_impl { 
# 815
typedef __decltype((__test< _Tp> (0))) type; 
# 816
}; 
# 818
template< class _Tp, bool 
# 819
 = __or_< is_void< _Tp> , __is_array_unknown_bounds< _Tp> , is_function< _Tp> > ::value, bool 
# 822
 = __or_< is_reference< _Tp> , is_scalar< _Tp> > ::value> struct __is_destructible_safe; 
# 825
template< class _Tp> 
# 826
struct __is_destructible_safe< _Tp, false, false>  : public __is_destructible_impl< typename remove_all_extents< _Tp> ::type> ::type { 
# 829
}; 
# 831
template< class _Tp> 
# 832
struct __is_destructible_safe< _Tp, true, false>  : public false_type { 
# 833
}; 
# 835
template< class _Tp> 
# 836
struct __is_destructible_safe< _Tp, false, true>  : public true_type { 
# 837
}; 
# 840
template< class _Tp> 
# 841
struct is_destructible : public __is_destructible_safe< _Tp> ::type { 
# 844
static_assert((std::__is_complete_or_unbounded(__type_identity< _Tp> {})), "template argument must be a complete class or an unbounded array");
# 846
}; 
# 852
struct __do_is_nt_destructible_impl { 
# 854
template< class _Tp> static __bool_constant< noexcept(declval< _Tp &> ().~_Tp())>  __test(int); 
# 858
template< class > static false_type __test(...); 
# 860
}; 
# 862
template< class _Tp> 
# 863
struct __is_nt_destructible_impl : public __do_is_nt_destructible_impl { 
# 866
typedef __decltype((__test< _Tp> (0))) type; 
# 867
}; 
# 869
template< class _Tp, bool 
# 870
 = __or_< is_void< _Tp> , __is_array_unknown_bounds< _Tp> , is_function< _Tp> > ::value, bool 
# 873
 = __or_< is_reference< _Tp> , is_scalar< _Tp> > ::value> struct __is_nt_destructible_safe; 
# 876
template< class _Tp> 
# 877
struct __is_nt_destructible_safe< _Tp, false, false>  : public __is_nt_destructible_impl< typename remove_all_extents< _Tp> ::type> ::type { 
# 880
}; 
# 882
template< class _Tp> 
# 883
struct __is_nt_destructible_safe< _Tp, true, false>  : public false_type { 
# 884
}; 
# 886
template< class _Tp> 
# 887
struct __is_nt_destructible_safe< _Tp, false, true>  : public true_type { 
# 888
}; 
# 891
template< class _Tp> 
# 892
struct is_nothrow_destructible : public __is_nt_destructible_safe< _Tp> ::type { 
# 895
static_assert((std::__is_complete_or_unbounded(__type_identity< _Tp> {})), "template argument must be a complete class or an unbounded array");
# 897
}; 
# 899
template< class _Tp, class ..._Args> 
# 900
struct __is_constructible_impl : public __bool_constant< __is_constructible(_Tp, _Args...)>  { 
# 902
}; 
# 905
template< class _Tp, class ..._Args> 
# 906
struct is_constructible : public __is_constructible_impl< _Tp, _Args...>  { 
# 909
static_assert((std::__is_complete_or_unbounded(__type_identity< _Tp> {})), "template argument must be a complete class or an unbounded array");
# 911
}; 
# 914
template< class _Tp> 
# 915
struct is_default_constructible : public __is_constructible_impl< _Tp> ::type { 
# 918
static_assert((std::__is_complete_or_unbounded(__type_identity< _Tp> {})), "template argument must be a complete class or an unbounded array");
# 920
}; 
# 922
template< class _Tp, bool  = __is_referenceable< _Tp> ::value> struct __is_copy_constructible_impl; 
# 925
template< class _Tp> 
# 926
struct __is_copy_constructible_impl< _Tp, false>  : public false_type { 
# 927
}; 
# 929
template< class _Tp> 
# 930
struct __is_copy_constructible_impl< _Tp, true>  : public __is_constructible_impl< _Tp, const _Tp &>  { 
# 932
}; 
# 935
template< class _Tp> 
# 936
struct is_copy_constructible : public __is_copy_constructible_impl< _Tp>  { 
# 939
static_assert((std::__is_complete_or_unbounded(__type_identity< _Tp> {})), "template argument must be a complete class or an unbounded array");
# 941
}; 
# 943
template< class _Tp, bool  = __is_referenceable< _Tp> ::value> struct __is_move_constructible_impl; 
# 946
template< class _Tp> 
# 947
struct __is_move_constructible_impl< _Tp, false>  : public false_type { 
# 948
}; 
# 950
template< class _Tp> 
# 951
struct __is_move_constructible_impl< _Tp, true>  : public __is_constructible_impl< _Tp, _Tp &&>  { 
# 953
}; 
# 956
template< class _Tp> 
# 957
struct is_move_constructible : public __is_move_constructible_impl< _Tp>  { 
# 960
static_assert((std::__is_complete_or_unbounded(__type_identity< _Tp> {})), "template argument must be a complete class or an unbounded array");
# 962
}; 
# 964
template< bool , class _Tp, class ..._Args> 
# 965
struct __is_nt_constructible_impl : public false_type { 
# 967
}; 
# 969
template< class _Tp, class ..._Args> 
# 970
struct __is_nt_constructible_impl< true, _Tp, _Args...>  : public __bool_constant< noexcept((_Tp(std::declval< _Args> ()...)))>  { 
# 972
}; 
# 974
template< class _Tp, class _Arg> 
# 975
struct __is_nt_constructible_impl< true, _Tp, _Arg>  : public __bool_constant< noexcept((static_cast< _Tp>(std::declval< _Arg> ())))>  { 
# 977
}; 
# 979
template< class _Tp> 
# 980
struct __is_nt_constructible_impl< true, _Tp>  : public __bool_constant< noexcept((_Tp()))>  { 
# 982
}; 
# 984
template< class _Tp, size_t _Num> 
# 985
struct __is_nt_constructible_impl< true, _Tp [_Num]>  : public __bool_constant< noexcept((typename remove_all_extents< _Tp> ::type()))>  { 
# 987
}; 
# 1001 "/usr/include/c++/10/type_traits" 3
template< class _Tp, class ..._Args> using __is_nothrow_constructible_impl = __is_nt_constructible_impl< __is_constructible(_Tp, _Args...), _Tp, _Args...> ; 
# 1007
template< class _Tp, class ..._Args> 
# 1008
struct is_nothrow_constructible : public __is_nothrow_constructible_impl< _Tp, _Args...> ::type { 
# 1011
static_assert((std::__is_complete_or_unbounded(__type_identity< _Tp> {})), "template argument must be a complete class or an unbounded array");
# 1013
}; 
# 1016
template< class _Tp> 
# 1017
struct is_nothrow_default_constructible : public __is_nothrow_constructible_impl< _Tp> ::type { 
# 1020
static_assert((std::__is_complete_or_unbounded(__type_identity< _Tp> {})), "template argument must be a complete class or an unbounded array");
# 1022
}; 
# 1025
template< class _Tp, bool  = __is_referenceable< _Tp> ::value> struct __is_nothrow_copy_constructible_impl; 
# 1028
template< class _Tp> 
# 1029
struct __is_nothrow_copy_constructible_impl< _Tp, false>  : public false_type { 
# 1030
}; 
# 1032
template< class _Tp> 
# 1033
struct __is_nothrow_copy_constructible_impl< _Tp, true>  : public __is_nothrow_constructible_impl< _Tp, const _Tp &>  { 
# 1035
}; 
# 1038
template< class _Tp> 
# 1039
struct is_nothrow_copy_constructible : public __is_nothrow_copy_constructible_impl< _Tp> ::type { 
# 1042
static_assert((std::__is_complete_or_unbounded(__type_identity< _Tp> {})), "template argument must be a complete class or an unbounded array");
# 1044
}; 
# 1046
template< class _Tp, bool  = __is_referenceable< _Tp> ::value> struct __is_nothrow_move_constructible_impl; 
# 1049
template< class _Tp> 
# 1050
struct __is_nothrow_move_constructible_impl< _Tp, false>  : public false_type { 
# 1051
}; 
# 1053
template< class _Tp> 
# 1054
struct __is_nothrow_move_constructible_impl< _Tp, true>  : public __is_nothrow_constructible_impl< _Tp, _Tp &&>  { 
# 1056
}; 
# 1059
template< class _Tp> 
# 1060
struct is_nothrow_move_constructible : public __is_nothrow_move_constructible_impl< _Tp> ::type { 
# 1063
static_assert((std::__is_complete_or_unbounded(__type_identity< _Tp> {})), "template argument must be a complete class or an unbounded array");
# 1065
}; 
# 1068
template< class _Tp, class _Up> 
# 1069
struct is_assignable : public __bool_constant< __is_assignable(_Tp, _Up)>  { 
# 1072
static_assert((std::__is_complete_or_unbounded(__type_identity< _Tp> {})), "template argument must be a complete class or an unbounded array");
# 1074
}; 
# 1076
template< class _Tp, bool  = __is_referenceable< _Tp> ::value> struct __is_copy_assignable_impl; 
# 1079
template< class _Tp> 
# 1080
struct __is_copy_assignable_impl< _Tp, false>  : public false_type { 
# 1081
}; 
# 1083
template< class _Tp> 
# 1084
struct __is_copy_assignable_impl< _Tp, true>  : public __bool_constant< __is_assignable(_Tp &, const _Tp &)>  { 
# 1086
}; 
# 1089
template< class _Tp> 
# 1090
struct is_copy_assignable : public __is_copy_assignable_impl< _Tp> ::type { 
# 1093
static_assert((std::__is_complete_or_unbounded(__type_identity< _Tp> {})), "template argument must be a complete class or an unbounded array");
# 1095
}; 
# 1097
template< class _Tp, bool  = __is_referenceable< _Tp> ::value> struct __is_move_assignable_impl; 
# 1100
template< class _Tp> 
# 1101
struct __is_move_assignable_impl< _Tp, false>  : public false_type { 
# 1102
}; 
# 1104
template< class _Tp> 
# 1105
struct __is_move_assignable_impl< _Tp, true>  : public __bool_constant< __is_assignable(_Tp &, _Tp &&)>  { 
# 1107
}; 
# 1110
template< class _Tp> 
# 1111
struct is_move_assignable : public __is_move_assignable_impl< _Tp> ::type { 
# 1114
static_assert((std::__is_complete_or_unbounded(__type_identity< _Tp> {})), "template argument must be a complete class or an unbounded array");
# 1116
}; 
# 1118
template< class _Tp, class _Up> 
# 1119
struct __is_nt_assignable_impl : public integral_constant< bool, noexcept((declval< _Tp> () = declval< _Up> ()))>  { 
# 1121
}; 
# 1123
template< class _Tp, class _Up> 
# 1124
struct __is_nothrow_assignable_impl : public __and_< __bool_constant< __is_assignable(_Tp, _Up)> , __is_nt_assignable_impl< _Tp, _Up> >  { 
# 1127
}; 
# 1130
template< class _Tp, class _Up> 
# 1131
struct is_nothrow_assignable : public __is_nothrow_assignable_impl< _Tp, _Up>  { 
# 1134
static_assert((std::__is_complete_or_unbounded(__type_identity< _Tp> {})), "template argument must be a complete class or an unbounded array");
# 1136
}; 
# 1138
template< class _Tp, bool  = __is_referenceable< _Tp> ::value> struct __is_nt_copy_assignable_impl; 
# 1141
template< class _Tp> 
# 1142
struct __is_nt_copy_assignable_impl< _Tp, false>  : public false_type { 
# 1143
}; 
# 1145
template< class _Tp> 
# 1146
struct __is_nt_copy_assignable_impl< _Tp, true>  : public __is_nothrow_assignable_impl< _Tp &, const _Tp &>  { 
# 1148
}; 
# 1151
template< class _Tp> 
# 1152
struct is_nothrow_copy_assignable : public __is_nt_copy_assignable_impl< _Tp>  { 
# 1155
static_assert((std::__is_complete_or_unbounded(__type_identity< _Tp> {})), "template argument must be a complete class or an unbounded array");
# 1157
}; 
# 1159
template< class _Tp, bool  = __is_referenceable< _Tp> ::value> struct __is_nt_move_assignable_impl; 
# 1162
template< class _Tp> 
# 1163
struct __is_nt_move_assignable_impl< _Tp, false>  : public false_type { 
# 1164
}; 
# 1166
template< class _Tp> 
# 1167
struct __is_nt_move_assignable_impl< _Tp, true>  : public __is_nothrow_assignable_impl< _Tp &, _Tp &&>  { 
# 1169
}; 
# 1172
template< class _Tp> 
# 1173
struct is_nothrow_move_assignable : public __is_nt_move_assignable_impl< _Tp>  { 
# 1176
static_assert((std::__is_complete_or_unbounded(__type_identity< _Tp> {})), "template argument must be a complete class or an unbounded array");
# 1178
}; 
# 1181
template< class _Tp, class ..._Args> 
# 1182
struct is_trivially_constructible : public __bool_constant< __is_trivially_constructible(_Tp, _Args...)>  { 
# 1185
static_assert((std::__is_complete_or_unbounded(__type_identity< _Tp> {})), "template argument must be a complete class or an unbounded array");
# 1187
}; 
# 1190
template< class _Tp> 
# 1191
struct is_trivially_default_constructible : public __bool_constant< __is_trivially_constructible(_Tp)>  { 
# 1194
static_assert((std::__is_complete_or_unbounded(__type_identity< _Tp> {})), "template argument must be a complete class or an unbounded array");
# 1196
}; 
# 1198
struct __do_is_implicitly_default_constructible_impl { 
# 1200
template< class _Tp> static void __helper(const _Tp &); 
# 1203
template< class _Tp> static true_type __test(const _Tp &, __decltype((__helper< const _Tp &> ({}))) * = 0); 
# 1207
static false_type __test(...); 
# 1208
}; 
# 1210
template< class _Tp> 
# 1211
struct __is_implicitly_default_constructible_impl : public __do_is_implicitly_default_constructible_impl { 
# 1214
typedef __decltype((__test(declval< _Tp> ()))) type; 
# 1215
}; 
# 1217
template< class _Tp> 
# 1218
struct __is_implicitly_default_constructible_safe : public __is_implicitly_default_constructible_impl< _Tp> ::type { 
# 1220
}; 
# 1222
template< class _Tp> 
# 1223
struct __is_implicitly_default_constructible : public __and_< __is_constructible_impl< _Tp> , __is_implicitly_default_constructible_safe< _Tp> >  { 
# 1226
}; 
# 1228
template< class _Tp, bool  = __is_referenceable< _Tp> ::value> struct __is_trivially_copy_constructible_impl; 
# 1231
template< class _Tp> 
# 1232
struct __is_trivially_copy_constructible_impl< _Tp, false>  : public false_type { 
# 1233
}; 
# 1235
template< class _Tp> 
# 1236
struct __is_trivially_copy_constructible_impl< _Tp, true>  : public __and_< __is_copy_constructible_impl< _Tp> , integral_constant< bool, __is_trivially_constructible(_Tp, const _Tp &)> >  { 
# 1240
}; 
# 1243
template< class _Tp> 
# 1244
struct is_trivially_copy_constructible : public __is_trivially_copy_constructible_impl< _Tp>  { 
# 1247
static_assert((std::__is_complete_or_unbounded(__type_identity< _Tp> {})), "template argument must be a complete class or an unbounded array");
# 1249
}; 
# 1251
template< class _Tp, bool  = __is_referenceable< _Tp> ::value> struct __is_trivially_move_constructible_impl; 
# 1254
template< class _Tp> 
# 1255
struct __is_trivially_move_constructible_impl< _Tp, false>  : public false_type { 
# 1256
}; 
# 1258
template< class _Tp> 
# 1259
struct __is_trivially_move_constructible_impl< _Tp, true>  : public __and_< __is_move_constructible_impl< _Tp> , integral_constant< bool, __is_trivially_constructible(_Tp, _Tp &&)> >  { 
# 1263
}; 
# 1266
template< class _Tp> 
# 1267
struct is_trivially_move_constructible : public __is_trivially_move_constructible_impl< _Tp>  { 
# 1270
static_assert((std::__is_complete_or_unbounded(__type_identity< _Tp> {})), "template argument must be a complete class or an unbounded array");
# 1272
}; 
# 1275
template< class _Tp, class _Up> 
# 1276
struct is_trivially_assignable : public __bool_constant< __is_trivially_assignable(_Tp, _Up)>  { 
# 1279
static_assert((std::__is_complete_or_unbounded(__type_identity< _Tp> {})), "template argument must be a complete class or an unbounded array");
# 1281
}; 
# 1283
template< class _Tp, bool  = __is_referenceable< _Tp> ::value> struct __is_trivially_copy_assignable_impl; 
# 1286
template< class _Tp> 
# 1287
struct __is_trivially_copy_assignable_impl< _Tp, false>  : public false_type { 
# 1288
}; 
# 1290
template< class _Tp> 
# 1291
struct __is_trivially_copy_assignable_impl< _Tp, true>  : public __bool_constant< __is_trivially_assignable(_Tp &, const _Tp &)>  { 
# 1293
}; 
# 1296
template< class _Tp> 
# 1297
struct is_trivially_copy_assignable : public __is_trivially_copy_assignable_impl< _Tp>  { 
# 1300
static_assert((std::__is_complete_or_unbounded(__type_identity< _Tp> {})), "template argument must be a complete class or an unbounded array");
# 1302
}; 
# 1304
template< class _Tp, bool  = __is_referenceable< _Tp> ::value> struct __is_trivially_move_assignable_impl; 
# 1307
template< class _Tp> 
# 1308
struct __is_trivially_move_assignable_impl< _Tp, false>  : public false_type { 
# 1309
}; 
# 1311
template< class _Tp> 
# 1312
struct __is_trivially_move_assignable_impl< _Tp, true>  : public __bool_constant< __is_trivially_assignable(_Tp &, _Tp &&)>  { 
# 1314
}; 
# 1317
template< class _Tp> 
# 1318
struct is_trivially_move_assignable : public __is_trivially_move_assignable_impl< _Tp>  { 
# 1321
static_assert((std::__is_complete_or_unbounded(__type_identity< _Tp> {})), "template argument must be a complete class or an unbounded array");
# 1323
}; 
# 1326
template< class _Tp> 
# 1327
struct is_trivially_destructible : public __and_< __is_destructible_safe< _Tp> , __bool_constant< __has_trivial_destructor(_Tp)> >  { 
# 1331
static_assert((std::__is_complete_or_unbounded(__type_identity< _Tp> {})), "template argument must be a complete class or an unbounded array");
# 1333
}; 
# 1337
template< class _Tp> 
# 1338
struct has_virtual_destructor : public integral_constant< bool, __has_virtual_destructor(_Tp)>  { 
# 1341
static_assert((std::__is_complete_or_unbounded(__type_identity< _Tp> {})), "template argument must be a complete class or an unbounded array");
# 1343
}; 
# 1349
template< class _Tp> 
# 1350
struct alignment_of : public integral_constant< unsigned long, __alignof__(_Tp)>  { 
# 1353
static_assert((std::__is_complete_or_unbounded(__type_identity< _Tp> {})), "template argument must be a complete class or an unbounded array");
# 1355
}; 
# 1358
template< class > 
# 1359
struct rank : public integral_constant< unsigned long, 0UL>  { 
# 1360
}; 
# 1362
template< class _Tp, size_t _Size> 
# 1363
struct rank< _Tp [_Size]>  : public integral_constant< unsigned long, 1 + std::rank< _Tp> ::value>  { 
# 1364
}; 
# 1366
template< class _Tp> 
# 1367
struct rank< _Tp []>  : public integral_constant< unsigned long, 1 + std::rank< _Tp> ::value>  { 
# 1368
}; 
# 1371
template< class , unsigned _Uint> 
# 1372
struct extent : public integral_constant< unsigned long, 0UL>  { 
# 1373
}; 
# 1375
template< class _Tp, unsigned _Uint, size_t _Size> 
# 1376
struct extent< _Tp [_Size], _Uint>  : public integral_constant< unsigned long, (_Uint == (0)) ? _Size : std::extent< _Tp, _Uint - (1)> ::value>  { 
# 1380
}; 
# 1382
template< class _Tp, unsigned _Uint> 
# 1383
struct extent< _Tp [], _Uint>  : public integral_constant< unsigned long, (_Uint == (0)) ? 0 : std::extent< _Tp, _Uint - (1)> ::value>  { 
# 1387
}; 
# 1393
template< class _Tp, class _Up> 
# 1394
struct is_same : public integral_constant< bool, __is_same_as(_Tp, _Up)>  { 
# 1400
}; 
# 1410 "/usr/include/c++/10/type_traits" 3
template< class _Base, class _Derived> 
# 1411
struct is_base_of : public integral_constant< bool, __is_base_of(_Base, _Derived)>  { 
# 1413
}; 
# 1415
template< class _From, class _To, bool 
# 1416
 = __or_< is_void< _From> , is_function< _To> , is_array< _To> > ::value> 
# 1418
struct __is_convertible_helper { 
# 1420
typedef typename is_void< _To> ::type type; 
# 1421
}; 
# 1423
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wctor-dtor-privacy"
template< class _From, class _To> 
# 1426
class __is_convertible_helper< _From, _To, false>  { 
# 1428
template< class _To1> static void __test_aux(_To1) noexcept; 
# 1431
template< class _From1, class _To1, class 
# 1432
 = __decltype((__test_aux< _To1> (std::declval< _From1> ())))> static true_type 
# 1431
__test(int); 
# 1436
template< class , class > static false_type __test(...); 
# 1441
public: typedef __decltype((__test< _From, _To> (0))) type; 
# 1442
}; 
#pragma GCC diagnostic pop
# 1446
template< class _From, class _To> 
# 1447
struct is_convertible : public __is_convertible_helper< _From, _To> ::type { 
# 1449
}; 
# 1452
template< class _ToElementType, class _FromElementType> using __is_array_convertible = is_convertible< _FromElementType (*)[], _ToElementType (*)[]> ; 
# 1456
template< class _From, class _To, bool 
# 1457
 = __or_< is_void< _From> , is_function< _To> , is_array< _To> > ::value> 
# 1459
struct __is_nt_convertible_helper : public is_void< _To>  { 
# 1461
}; 
# 1463
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wctor-dtor-privacy"
template< class _From, class _To> 
# 1466
class __is_nt_convertible_helper< _From, _To, false>  { 
# 1468
template< class _To1> static void __test_aux(_To1) noexcept; 
# 1471
template< class _From1, class _To1> static __bool_constant< noexcept(__test_aux< _To1> (std::declval< _From1> ()))>  __test(int); 
# 1476
template< class , class > static false_type __test(...); 
# 1481
public: using type = __decltype((__test< _From, _To> (0))); 
# 1482
}; 
#pragma GCC diagnostic pop
# 1486
template< class _From, class _To> 
# 1487
struct __is_nothrow_convertible : public __is_nt_convertible_helper< _From, _To> ::type { 
# 1489
}; 
# 1508 "/usr/include/c++/10/type_traits" 3
template< class _Tp> 
# 1509
struct remove_const { 
# 1510
typedef _Tp type; }; 
# 1512
template< class _Tp> 
# 1513
struct remove_const< const _Tp>  { 
# 1514
typedef _Tp type; }; 
# 1517
template< class _Tp> 
# 1518
struct remove_volatile { 
# 1519
typedef _Tp type; }; 
# 1521
template< class _Tp> 
# 1522
struct remove_volatile< volatile _Tp>  { 
# 1523
typedef _Tp type; }; 
# 1526
template< class _Tp> 
# 1527
struct remove_cv { 
# 1528
using type = _Tp; }; 
# 1530
template< class _Tp> 
# 1531
struct remove_cv< const _Tp>  { 
# 1532
using type = _Tp; }; 
# 1534
template< class _Tp> 
# 1535
struct remove_cv< volatile _Tp>  { 
# 1536
using type = _Tp; }; 
# 1538
template< class _Tp> 
# 1539
struct remove_cv< const volatile _Tp>  { 
# 1540
using type = _Tp; }; 
# 1543
template< class _Tp> 
# 1544
struct add_const { 
# 1545
typedef const _Tp type; }; 
# 1548
template< class _Tp> 
# 1549
struct add_volatile { 
# 1550
typedef volatile _Tp type; }; 
# 1553
template< class _Tp> 
# 1554
struct add_cv { 
# 1557
typedef typename add_const< typename add_volatile< _Tp> ::type> ::type type; 
# 1558
}; 
# 1565
template< class _Tp> using remove_const_t = typename remove_const< _Tp> ::type; 
# 1569
template< class _Tp> using remove_volatile_t = typename remove_volatile< _Tp> ::type; 
# 1573
template< class _Tp> using remove_cv_t = typename remove_cv< _Tp> ::type; 
# 1577
template< class _Tp> using add_const_t = typename add_const< _Tp> ::type; 
# 1581
template< class _Tp> using add_volatile_t = typename add_volatile< _Tp> ::type; 
# 1585
template< class _Tp> using add_cv_t = typename add_cv< _Tp> ::type; 
# 1592
template< class _Tp> 
# 1593
struct remove_reference { 
# 1594
typedef _Tp type; }; 
# 1596
template< class _Tp> 
# 1597
struct remove_reference< _Tp &>  { 
# 1598
typedef _Tp type; }; 
# 1600
template< class _Tp> 
# 1601
struct remove_reference< _Tp &&>  { 
# 1602
typedef _Tp type; }; 
# 1604
template< class _Tp, bool  = __is_referenceable< _Tp> ::value> 
# 1605
struct __add_lvalue_reference_helper { 
# 1606
typedef _Tp type; }; 
# 1608
template< class _Tp> 
# 1609
struct __add_lvalue_reference_helper< _Tp, true>  { 
# 1610
typedef _Tp &type; }; 
# 1613
template< class _Tp> 
# 1614
struct add_lvalue_reference : public __add_lvalue_reference_helper< _Tp>  { 
# 1616
}; 
# 1618
template< class _Tp, bool  = __is_referenceable< _Tp> ::value> 
# 1619
struct __add_rvalue_reference_helper { 
# 1620
typedef _Tp type; }; 
# 1622
template< class _Tp> 
# 1623
struct __add_rvalue_reference_helper< _Tp, true>  { 
# 1624
typedef _Tp &&type; }; 
# 1627
template< class _Tp> 
# 1628
struct add_rvalue_reference : public __add_rvalue_reference_helper< _Tp>  { 
# 1630
}; 
# 1634
template< class _Tp> using remove_reference_t = typename remove_reference< _Tp> ::type; 
# 1638
template< class _Tp> using add_lvalue_reference_t = typename add_lvalue_reference< _Tp> ::type; 
# 1642
template< class _Tp> using add_rvalue_reference_t = typename add_rvalue_reference< _Tp> ::type; 
# 1649
template< class _Unqualified, bool _IsConst, bool _IsVol> struct __cv_selector; 
# 1652
template< class _Unqualified> 
# 1653
struct __cv_selector< _Unqualified, false, false>  { 
# 1654
typedef _Unqualified __type; }; 
# 1656
template< class _Unqualified> 
# 1657
struct __cv_selector< _Unqualified, false, true>  { 
# 1658
typedef volatile _Unqualified __type; }; 
# 1660
template< class _Unqualified> 
# 1661
struct __cv_selector< _Unqualified, true, false>  { 
# 1662
typedef const _Unqualified __type; }; 
# 1664
template< class _Unqualified> 
# 1665
struct __cv_selector< _Unqualified, true, true>  { 
# 1666
typedef const volatile _Unqualified __type; }; 
# 1668
template< class _Qualified, class _Unqualified, bool 
# 1669
_IsConst = is_const< _Qualified> ::value, bool 
# 1670
_IsVol = is_volatile< _Qualified> ::value> 
# 1671
class __match_cv_qualifiers { 
# 1673
typedef __cv_selector< _Unqualified, _IsConst, _IsVol>  __match; 
# 1676
public: typedef typename __cv_selector< _Unqualified, _IsConst, _IsVol> ::__type __type; 
# 1677
}; 
# 1680
template< class _Tp> 
# 1681
struct __make_unsigned { 
# 1682
typedef _Tp __type; }; 
# 1685
template<> struct __make_unsigned< char>  { 
# 1686
typedef unsigned char __type; }; 
# 1689
template<> struct __make_unsigned< signed char>  { 
# 1690
typedef unsigned char __type; }; 
# 1693
template<> struct __make_unsigned< short>  { 
# 1694
typedef unsigned short __type; }; 
# 1697
template<> struct __make_unsigned< int>  { 
# 1698
typedef unsigned __type; }; 
# 1701
template<> struct __make_unsigned< long>  { 
# 1702
typedef unsigned long __type; }; 
# 1705
template<> struct __make_unsigned< long long>  { 
# 1706
typedef unsigned long long __type; }; 
# 1710
template<> struct __make_unsigned< __int128>  { 
# 1711
typedef unsigned __int128 __type; }; 
# 1730 "/usr/include/c++/10/type_traits" 3
template< class _Tp, bool 
# 1731
_IsInt = is_integral< _Tp> ::value, bool 
# 1732
_IsEnum = is_enum< _Tp> ::value> class __make_unsigned_selector; 
# 1735
template< class _Tp> 
# 1736
class __make_unsigned_selector< _Tp, true, false>  { 
# 1738
using __unsigned_type = typename __make_unsigned< __remove_cv_t< _Tp> > ::__type; 
# 1742
public: using __type = typename __match_cv_qualifiers< _Tp, __unsigned_type> ::__type; 
# 1744
}; 
# 1746
class __make_unsigned_selector_base { 
# 1749
protected: template< class ...> struct _List { }; 
# 1751
template< class _Tp, class ..._Up> 
# 1752
struct _List< _Tp, _Up...>  : public __make_unsigned_selector_base::_List< _Up...>  { 
# 1753
static constexpr std::size_t __size = sizeof(_Tp); }; 
# 1755
template< size_t _Sz, class _Tp, bool  = _Sz <= _Tp::__size> struct __select; 
# 1758
template< size_t _Sz, class _Uint, class ..._UInts> 
# 1759
struct __select< _Sz, _List< _Uint, _UInts...> , true>  { 
# 1760
using __type = _Uint; }; 
# 1762
template< size_t _Sz, class _Uint, class ..._UInts> 
# 1763
struct __select< _Sz, _List< _Uint, _UInts...> , false>  : public __make_unsigned_selector_base::__select< _Sz, _List< _UInts...> >  { 
# 1765
}; 
# 1766
}; 
# 1769
template< class _Tp> 
# 1770
class __make_unsigned_selector< _Tp, false, true>  : private __make_unsigned_selector_base { 
# 1774
using _UInts = _List< unsigned char, unsigned short, unsigned, unsigned long, unsigned long long> ; 
# 1777
using __unsigned_type = typename __select< sizeof(_Tp), _List< unsigned char, unsigned short, unsigned, unsigned long, unsigned long long> > ::__type; 
# 1780
public: using __type = typename __match_cv_qualifiers< _Tp, __unsigned_type> ::__type; 
# 1782
}; 
# 1790
template<> struct __make_unsigned< wchar_t>  { 
# 1792
using __type = __make_unsigned_selector< wchar_t, false, true> ::__type; 
# 1794
}; 
# 1807 "/usr/include/c++/10/type_traits" 3
template<> struct __make_unsigned< char16_t>  { 
# 1809
using __type = __make_unsigned_selector< char16_t, false, true> ::__type; 
# 1811
}; 
# 1814
template<> struct __make_unsigned< char32_t>  { 
# 1816
using __type = __make_unsigned_selector< char32_t, false, true> ::__type; 
# 1818
}; 
# 1824
template< class _Tp> 
# 1825
struct make_unsigned { 
# 1826
typedef typename __make_unsigned_selector< _Tp> ::__type type; }; 
# 1830
template<> struct make_unsigned< bool> ; 
# 1834
template< class _Tp> 
# 1835
struct __make_signed { 
# 1836
typedef _Tp __type; }; 
# 1839
template<> struct __make_signed< char>  { 
# 1840
typedef signed char __type; }; 
# 1843
template<> struct __make_signed< unsigned char>  { 
# 1844
typedef signed char __type; }; 
# 1847
template<> struct __make_signed< unsigned short>  { 
# 1848
typedef signed short __type; }; 
# 1851
template<> struct __make_signed< unsigned>  { 
# 1852
typedef signed int __type; }; 
# 1855
template<> struct __make_signed< unsigned long>  { 
# 1856
typedef signed long __type; }; 
# 1859
template<> struct __make_signed< unsigned long long>  { 
# 1860
typedef signed long long __type; }; 
# 1864
template<> struct __make_signed< unsigned __int128>  { 
# 1865
typedef __int128 __type; }; 
# 1884 "/usr/include/c++/10/type_traits" 3
template< class _Tp, bool 
# 1885
_IsInt = is_integral< _Tp> ::value, bool 
# 1886
_IsEnum = is_enum< _Tp> ::value> class __make_signed_selector; 
# 1889
template< class _Tp> 
# 1890
class __make_signed_selector< _Tp, true, false>  { 
# 1892
using __signed_type = typename __make_signed< __remove_cv_t< _Tp> > ::__type; 
# 1896
public: using __type = typename __match_cv_qualifiers< _Tp, __signed_type> ::__type; 
# 1898
}; 
# 1901
template< class _Tp> 
# 1902
class __make_signed_selector< _Tp, false, true>  { 
# 1904
typedef typename __make_unsigned_selector< _Tp> ::__type __unsigned_type; 
# 1907
public: typedef typename std::__make_signed_selector< __unsigned_type> ::__type __type; 
# 1908
}; 
# 1916
template<> struct __make_signed< wchar_t>  { 
# 1918
using __type = __make_signed_selector< wchar_t, false, true> ::__type; 
# 1920
}; 
# 1933 "/usr/include/c++/10/type_traits" 3
template<> struct __make_signed< char16_t>  { 
# 1935
using __type = __make_signed_selector< char16_t, false, true> ::__type; 
# 1937
}; 
# 1940
template<> struct __make_signed< char32_t>  { 
# 1942
using __type = __make_signed_selector< char32_t, false, true> ::__type; 
# 1944
}; 
# 1950
template< class _Tp> 
# 1951
struct make_signed { 
# 1952
typedef typename __make_signed_selector< _Tp> ::__type type; }; 
# 1956
template<> struct make_signed< bool> ; 
# 1960
template< class _Tp> using make_signed_t = typename make_signed< _Tp> ::type; 
# 1964
template< class _Tp> using make_unsigned_t = typename make_unsigned< _Tp> ::type; 
# 1971
template< class _Tp> 
# 1972
struct remove_extent { 
# 1973
typedef _Tp type; }; 
# 1975
template< class _Tp, size_t _Size> 
# 1976
struct remove_extent< _Tp [_Size]>  { 
# 1977
typedef _Tp type; }; 
# 1979
template< class _Tp> 
# 1980
struct remove_extent< _Tp []>  { 
# 1981
typedef _Tp type; }; 
# 1984
template< class _Tp> 
# 1985
struct remove_all_extents { 
# 1986
typedef _Tp type; }; 
# 1988
template< class _Tp, size_t _Size> 
# 1989
struct remove_all_extents< _Tp [_Size]>  { 
# 1990
typedef typename std::remove_all_extents< _Tp> ::type type; }; 
# 1992
template< class _Tp> 
# 1993
struct remove_all_extents< _Tp []>  { 
# 1994
typedef typename std::remove_all_extents< _Tp> ::type type; }; 
# 1998
template< class _Tp> using remove_extent_t = typename remove_extent< _Tp> ::type; 
# 2002
template< class _Tp> using remove_all_extents_t = typename remove_all_extents< _Tp> ::type; 
# 2008
template< class _Tp, class > 
# 2009
struct __remove_pointer_helper { 
# 2010
typedef _Tp type; }; 
# 2012
template< class _Tp, class _Up> 
# 2013
struct __remove_pointer_helper< _Tp, _Up *>  { 
# 2014
typedef _Up type; }; 
# 2017
template< class _Tp> 
# 2018
struct remove_pointer : public __remove_pointer_helper< _Tp, __remove_cv_t< _Tp> >  { 
# 2020
}; 
# 2023
template< class _Tp, bool  = __or_< __is_referenceable< _Tp> , is_void< _Tp> > ::value> 
# 2025
struct __add_pointer_helper { 
# 2026
typedef _Tp type; }; 
# 2028
template< class _Tp> 
# 2029
struct __add_pointer_helper< _Tp, true>  { 
# 2030
typedef typename remove_reference< _Tp> ::type *type; }; 
# 2032
template< class _Tp> 
# 2033
struct add_pointer : public __add_pointer_helper< _Tp>  { 
# 2035
}; 
# 2039
template< class _Tp> using remove_pointer_t = typename remove_pointer< _Tp> ::type; 
# 2043
template< class _Tp> using add_pointer_t = typename add_pointer< _Tp> ::type; 
# 2047
template< size_t _Len> 
# 2048
struct __aligned_storage_msa { 
# 2050
union __type { 
# 2052
unsigned char __data[_Len]; 
# 2053
struct __attribute((__aligned__)) { } __align; 
# 2054
}; 
# 2055
}; 
# 2067 "/usr/include/c++/10/type_traits" 3
template< size_t _Len, size_t _Align = __alignof__(typename __aligned_storage_msa< _Len> ::__type)> 
# 2069
struct aligned_storage { 
# 2071
union type { 
# 2073
unsigned char __data[_Len]; 
# 2074
struct __attribute((__aligned__(_Align))) { } __align; 
# 2075
}; 
# 2076
}; 
# 2078
template< class ..._Types> 
# 2079
struct __strictest_alignment { 
# 2081
static const size_t _S_alignment = (0); 
# 2082
static const size_t _S_size = (0); 
# 2083
}; 
# 2085
template< class _Tp, class ..._Types> 
# 2086
struct __strictest_alignment< _Tp, _Types...>  { 
# 2088
static const size_t _S_alignment = ((__alignof__(_Tp) > __strictest_alignment< _Types...> ::_S_alignment) ? __alignof__(_Tp) : __strictest_alignment< _Types...> ::_S_alignment); 
# 2091
static const size_t _S_size = ((sizeof(_Tp) > __strictest_alignment< _Types...> ::_S_size) ? sizeof(_Tp) : __strictest_alignment< _Types...> ::_S_size); 
# 2094
}; 
# 2106 "/usr/include/c++/10/type_traits" 3
template< size_t _Len, class ..._Types> 
# 2107
struct aligned_union { 
# 2110
static_assert((sizeof...(_Types) != (0)), "At least one type is required");
# 2112
private: using __strictest = __strictest_alignment< _Types...> ; 
# 2113
static const size_t _S_len = ((_Len > __strictest::_S_size) ? _Len : __strictest::_S_size); 
# 2117
public: static const size_t alignment_value = (__strictest::_S_alignment); 
# 2119
typedef typename aligned_storage< _S_len, alignment_value> ::type type; 
# 2120
}; 
# 2122
template< size_t _Len, class ..._Types> const size_t aligned_union< _Len, _Types...> ::alignment_value; 
# 2127
template< class _Up, bool 
# 2128
_IsArray = is_array< _Up> ::value, bool 
# 2129
_IsFunction = is_function< _Up> ::value> struct __decay_selector; 
# 2133
template< class _Up> 
# 2134
struct __decay_selector< _Up, false, false>  { 
# 2135
typedef __remove_cv_t< _Up>  __type; }; 
# 2137
template< class _Up> 
# 2138
struct __decay_selector< _Up, true, false>  { 
# 2139
typedef typename remove_extent< _Up> ::type *__type; }; 
# 2141
template< class _Up> 
# 2142
struct __decay_selector< _Up, false, true>  { 
# 2143
typedef typename add_pointer< _Up> ::type __type; }; 
# 2146
template< class _Tp> 
# 2147
class decay { 
# 2149
typedef typename remove_reference< _Tp> ::type __remove_type; 
# 2152
public: typedef typename __decay_selector< __remove_type> ::__type type; 
# 2153
}; 
# 2156
template< class _Tp> using __decay_t = typename decay< _Tp> ::type; 
# 2159
template< class _Tp> class reference_wrapper; 
# 2163
template< class _Tp> 
# 2164
struct __strip_reference_wrapper { 
# 2166
typedef _Tp __type; 
# 2167
}; 
# 2169
template< class _Tp> 
# 2170
struct __strip_reference_wrapper< reference_wrapper< _Tp> >  { 
# 2172
typedef _Tp &__type; 
# 2173
}; 
# 2175
template< class _Tp> using __decay_and_strip = __strip_reference_wrapper< __decay_t< _Tp> > ; 
# 2181
template< bool , class _Tp = void> 
# 2182
struct enable_if { 
# 2183
}; 
# 2186
template< class _Tp> 
# 2187
struct enable_if< true, _Tp>  { 
# 2188
typedef _Tp type; }; 
# 2191
template< bool _Cond, class _Tp = void> using __enable_if_t = typename enable_if< _Cond, _Tp> ::type; 
# 2194
template< class ..._Cond> using _Require = __enable_if_t< __and_< _Cond...> ::value> ; 
# 2199
template< bool _Cond, class _Iftrue, class _Iffalse> 
# 2200
struct conditional { 
# 2201
typedef _Iftrue type; }; 
# 2204
template< class _Iftrue, class _Iffalse> 
# 2205
struct conditional< false, _Iftrue, _Iffalse>  { 
# 2206
typedef _Iffalse type; }; 
# 2209
template< class _Tp> using __remove_cvref_t = typename remove_cv< typename remove_reference< _Tp> ::type> ::type; 
# 2214
template< class ..._Tp> struct common_type; 
# 2219
struct __do_common_type_impl { 
# 2221
template< class _Tp, class _Up> using __cond_t = __decltype((true ? std::declval< _Tp> () : std::declval< _Up> ())); 
# 2227
template< class _Tp, class _Up> static __success_type< __decay_t< __cond_t< _Tp, _Up> > >  _S_test(int); 
# 2239 "/usr/include/c++/10/type_traits" 3
template< class , class > static __failure_type _S_test_2(...); 
# 2243
template< class _Tp, class _Up> static __decltype((_S_test_2< _Tp, _Up> (0))) _S_test(...); 
# 2246
}; 
# 2250
template<> struct common_type< >  { 
# 2251
}; 
# 2254
template< class _Tp0> 
# 2255
struct common_type< _Tp0>  : public std::common_type< _Tp0, _Tp0>  { 
# 2257
}; 
# 2260
template< class _Tp1, class _Tp2, class 
# 2261
_Dp1 = __decay_t< _Tp1> , class _Dp2 = __decay_t< _Tp2> > 
# 2262
struct __common_type_impl { 
# 2266
using type = common_type< _Dp1, _Dp2> ; 
# 2267
}; 
# 2269
template< class _Tp1, class _Tp2> 
# 2270
struct __common_type_impl< _Tp1, _Tp2, _Tp1, _Tp2>  : private __do_common_type_impl { 
# 2275
using type = __decltype((_S_test< _Tp1, _Tp2> (0))); 
# 2276
}; 
# 2279
template< class _Tp1, class _Tp2> 
# 2280
struct common_type< _Tp1, _Tp2>  : public __common_type_impl< _Tp1, _Tp2> ::type { 
# 2282
}; 
# 2284
template< class ...> 
# 2285
struct __common_type_pack { 
# 2286
}; 
# 2288
template< class , class , class  = void> struct __common_type_fold; 
# 2292
template< class _Tp1, class _Tp2, class ..._Rp> 
# 2293
struct common_type< _Tp1, _Tp2, _Rp...>  : public __common_type_fold< std::common_type< _Tp1, _Tp2> , __common_type_pack< _Rp...> >  { 
# 2296
}; 
# 2301
template< class _CTp, class ..._Rp> 
# 2302
struct __common_type_fold< _CTp, __common_type_pack< _Rp...> , __void_t< typename _CTp::type> >  : public common_type< typename _CTp::type, _Rp...>  { 
# 2305
}; 
# 2308
template< class _CTp, class _Rp> 
# 2309
struct __common_type_fold< _CTp, _Rp, void>  { 
# 2310
}; 
# 2312
template< class _Tp, bool  = is_enum< _Tp> ::value> 
# 2313
struct __underlying_type_impl { 
# 2315
using type = __underlying_type(_Tp); 
# 2316
}; 
# 2318
template< class _Tp> 
# 2319
struct __underlying_type_impl< _Tp, false>  { 
# 2320
}; 
# 2323
template< class _Tp> 
# 2324
struct underlying_type : public __underlying_type_impl< _Tp>  { 
# 2326
}; 
# 2328
template< class _Tp> 
# 2329
struct __declval_protector { 
# 2331
static const bool __stop = false; 
# 2332
}; 
# 2334
template< class _Tp> auto 
# 2335
declval() noexcept->__decltype((__declval< _Tp> (0))) 
# 2336
{ 
# 2337
static_assert((__declval_protector< _Tp> ::__stop), "declval() must not be used!");
# 2339
return __declval< _Tp> (0); 
# 2340
} 
# 2343
template< class _Signature> class result_of; 
# 2350
struct __invoke_memfun_ref { }; 
# 2351
struct __invoke_memfun_deref { }; 
# 2352
struct __invoke_memobj_ref { }; 
# 2353
struct __invoke_memobj_deref { }; 
# 2354
struct __invoke_other { }; 
# 2357
template< class _Tp, class _Tag> 
# 2358
struct __result_of_success : public __success_type< _Tp>  { 
# 2359
using __invoke_type = _Tag; }; 
# 2362
struct __result_of_memfun_ref_impl { 
# 2364
template< class _Fp, class _Tp1, class ..._Args> static __result_of_success< __decltype(((std::declval< _Tp1> ().*std::declval< _Fp> ())(std::declval< _Args> ()...))), __invoke_memfun_ref>  _S_test(int); 
# 2369
template< class ...> static __failure_type _S_test(...); 
# 2371
}; 
# 2373
template< class _MemPtr, class _Arg, class ..._Args> 
# 2374
struct __result_of_memfun_ref : private __result_of_memfun_ref_impl { 
# 2377
typedef __decltype((_S_test< _MemPtr, _Arg, _Args...> (0))) type; 
# 2378
}; 
# 2381
struct __result_of_memfun_deref_impl { 
# 2383
template< class _Fp, class _Tp1, class ..._Args> static __result_of_success< __decltype((((*std::declval< _Tp1> ()).*std::declval< _Fp> ())(std::declval< _Args> ()...))), __invoke_memfun_deref>  _S_test(int); 
# 2388
template< class ...> static __failure_type _S_test(...); 
# 2390
}; 
# 2392
template< class _MemPtr, class _Arg, class ..._Args> 
# 2393
struct __result_of_memfun_deref : private __result_of_memfun_deref_impl { 
# 2396
typedef __decltype((_S_test< _MemPtr, _Arg, _Args...> (0))) type; 
# 2397
}; 
# 2400
struct __result_of_memobj_ref_impl { 
# 2402
template< class _Fp, class _Tp1> static __result_of_success< __decltype((std::declval< _Tp1> ().*std::declval< _Fp> ())), __invoke_memobj_ref>  _S_test(int); 
# 2407
template< class , class > static __failure_type _S_test(...); 
# 2409
}; 
# 2411
template< class _MemPtr, class _Arg> 
# 2412
struct __result_of_memobj_ref : private __result_of_memobj_ref_impl { 
# 2415
typedef __decltype((_S_test< _MemPtr, _Arg> (0))) type; 
# 2416
}; 
# 2419
struct __result_of_memobj_deref_impl { 
# 2421
template< class _Fp, class _Tp1> static __result_of_success< __decltype(((*std::declval< _Tp1> ()).*std::declval< _Fp> ())), __invoke_memobj_deref>  _S_test(int); 
# 2426
template< class , class > static __failure_type _S_test(...); 
# 2428
}; 
# 2430
template< class _MemPtr, class _Arg> 
# 2431
struct __result_of_memobj_deref : private __result_of_memobj_deref_impl { 
# 2434
typedef __decltype((_S_test< _MemPtr, _Arg> (0))) type; 
# 2435
}; 
# 2437
template< class _MemPtr, class _Arg> struct __result_of_memobj; 
# 2440
template< class _Res, class _Class, class _Arg> 
# 2441
struct __result_of_memobj< _Res (_Class::*), _Arg>  { 
# 2443
typedef __remove_cvref_t< _Arg>  _Argval; 
# 2444
typedef _Res (_Class::*_MemPtr); 
# 2449
typedef typename conditional< __or_< is_same< _Argval, _Class> , is_base_of< _Class, _Argval> > ::value, __result_of_memobj_ref< _MemPtr, _Arg> , __result_of_memobj_deref< _MemPtr, _Arg> > ::type::type type; 
# 2450
}; 
# 2452
template< class _MemPtr, class _Arg, class ..._Args> struct __result_of_memfun; 
# 2455
template< class _Res, class _Class, class _Arg, class ..._Args> 
# 2456
struct __result_of_memfun< _Res (_Class::*), _Arg, _Args...>  { 
# 2458
typedef typename remove_reference< _Arg> ::type _Argval; 
# 2459
typedef _Res (_Class::*_MemPtr); 
# 2463
typedef typename conditional< is_base_of< _Class, _Argval> ::value, __result_of_memfun_ref< _MemPtr, _Arg, _Args...> , __result_of_memfun_deref< _MemPtr, _Arg, _Args...> > ::type::type type; 
# 2464
}; 
# 2471
template< class _Tp, class _Up = __remove_cvref_t< _Tp> > 
# 2472
struct __inv_unwrap { 
# 2474
using type = _Tp; 
# 2475
}; 
# 2477
template< class _Tp, class _Up> 
# 2478
struct __inv_unwrap< _Tp, reference_wrapper< _Up> >  { 
# 2480
using type = _Up &; 
# 2481
}; 
# 2483
template< bool , bool , class _Functor, class ..._ArgTypes> 
# 2484
struct __result_of_impl { 
# 2486
typedef __failure_type type; 
# 2487
}; 
# 2489
template< class _MemPtr, class _Arg> 
# 2490
struct __result_of_impl< true, false, _MemPtr, _Arg>  : public __result_of_memobj< __decay_t< _MemPtr> , typename __inv_unwrap< _Arg> ::type>  { 
# 2493
}; 
# 2495
template< class _MemPtr, class _Arg, class ..._Args> 
# 2496
struct __result_of_impl< false, true, _MemPtr, _Arg, _Args...>  : public __result_of_memfun< __decay_t< _MemPtr> , typename __inv_unwrap< _Arg> ::type, _Args...>  { 
# 2499
}; 
# 2502
struct __result_of_other_impl { 
# 2504
template< class _Fn, class ..._Args> static __result_of_success< __decltype((std::declval< _Fn> ()(std::declval< _Args> ()...))), __invoke_other>  _S_test(int); 
# 2509
template< class ...> static __failure_type _S_test(...); 
# 2511
}; 
# 2513
template< class _Functor, class ..._ArgTypes> 
# 2514
struct __result_of_impl< false, false, _Functor, _ArgTypes...>  : private __result_of_other_impl { 
# 2517
typedef __decltype((_S_test< _Functor, _ArgTypes...> (0))) type; 
# 2518
}; 
# 2521
template< class _Functor, class ..._ArgTypes> 
# 2522
struct __invoke_result : public __result_of_impl< is_member_object_pointer< typename remove_reference< _Functor> ::type> ::value, is_member_function_pointer< typename remove_reference< _Functor> ::type> ::value, _Functor, _ArgTypes...> ::type { 
# 2532
}; 
# 2534
template< class _Functor, class ..._ArgTypes> 
# 2535
struct result_of< _Functor (_ArgTypes ...)>  : public __invoke_result< _Functor, _ArgTypes...>  { 
# 2537
}; 
# 2541
template< size_t _Len, size_t _Align = __alignof__(typename __aligned_storage_msa< _Len> ::__type)> using aligned_storage_t = typename aligned_storage< _Len, _Align> ::type; 
# 2545
template< size_t _Len, class ..._Types> using aligned_union_t = typename aligned_union< _Len, _Types...> ::type; 
# 2549
template< class _Tp> using decay_t = typename decay< _Tp> ::type; 
# 2553
template< bool _Cond, class _Tp = void> using enable_if_t = typename enable_if< _Cond, _Tp> ::type; 
# 2557
template< bool _Cond, class _Iftrue, class _Iffalse> using conditional_t = typename conditional< _Cond, _Iftrue, _Iffalse> ::type; 
# 2561
template< class ..._Tp> using common_type_t = typename common_type< _Tp...> ::type; 
# 2565
template< class _Tp> using underlying_type_t = typename underlying_type< _Tp> ::type; 
# 2569
template< class _Tp> using result_of_t = typename result_of< _Tp> ::type; 
# 2576
template< class ...> using void_t = void; 
# 2580
template< class _Default, class _AlwaysVoid, 
# 2581
template< class ...>  class _Op, class ..._Args> 
# 2582
struct __detector { 
# 2584
using value_t = false_type; 
# 2585
using type = _Default; 
# 2586
}; 
# 2589
template< class _Default, template< class ...>  class _Op, class ...
# 2590
_Args> 
# 2591
struct __detector< _Default, __void_t< _Op< _Args...> > , _Op, _Args...>  { 
# 2593
using value_t = true_type; 
# 2594
using type = _Op< _Args...> ; 
# 2595
}; 
# 2598
template< class _Default, template< class ...>  class _Op, class ...
# 2599
_Args> using __detected_or = __detector< _Default, void, _Op, _Args...> ; 
# 2603
template< class _Default, template< class ...>  class _Op, class ...
# 2604
_Args> using __detected_or_t = typename __detected_or< _Default, _Op, _Args...> ::type; 
# 2624 "/usr/include/c++/10/type_traits" 3
template< class _Tp> struct __is_swappable; 
# 2627
template< class _Tp> struct __is_nothrow_swappable; 
# 2630
template< class ..._Elements> class tuple; 
# 2633
template< class > 
# 2634
struct __is_tuple_like_impl : public false_type { 
# 2635
}; 
# 2637
template< class ..._Tps> 
# 2638
struct __is_tuple_like_impl< tuple< _Tps...> >  : public true_type { 
# 2639
}; 
# 2642
template< class _Tp> 
# 2643
struct __is_tuple_like : public __is_tuple_like_impl< __remove_cvref_t< _Tp> > ::type { 
# 2645
}; 
# 2647
template< class _Tp> inline _Require< __not_< __is_tuple_like< _Tp> > , is_move_constructible< _Tp> , is_move_assignable< _Tp> >  swap(_Tp &, _Tp &) noexcept(__and_< is_nothrow_move_constructible< _Tp> , is_nothrow_move_assignable< _Tp> > ::value); 
# 2657
template< class _Tp, size_t _Nm> inline __enable_if_t< __is_swappable< _Tp> ::value>  swap(_Tp (& __a)[_Nm], _Tp (& __b)[_Nm]) noexcept(__is_nothrow_swappable< _Tp> ::value); 
# 2664
namespace __swappable_details { 
# 2665
using std::swap;
# 2667
struct __do_is_swappable_impl { 
# 2669
template< class _Tp, class 
# 2670
 = __decltype((swap(std::declval< _Tp &> (), std::declval< _Tp &> ())))> static true_type 
# 2669
__test(int); 
# 2673
template< class > static false_type __test(...); 
# 2675
}; 
# 2677
struct __do_is_nothrow_swappable_impl { 
# 2679
template< class _Tp> static __bool_constant< noexcept(swap(std::declval< _Tp &> (), std::declval< _Tp &> ()))>  __test(int); 
# 2684
template< class > static false_type __test(...); 
# 2686
}; 
# 2688
}
# 2690
template< class _Tp> 
# 2691
struct __is_swappable_impl : public __swappable_details::__do_is_swappable_impl { 
# 2694
typedef __decltype((__test< _Tp> (0))) type; 
# 2695
}; 
# 2697
template< class _Tp> 
# 2698
struct __is_nothrow_swappable_impl : public __swappable_details::__do_is_nothrow_swappable_impl { 
# 2701
typedef __decltype((__test< _Tp> (0))) type; 
# 2702
}; 
# 2704
template< class _Tp> 
# 2705
struct __is_swappable : public __is_swappable_impl< _Tp> ::type { 
# 2707
}; 
# 2709
template< class _Tp> 
# 2710
struct __is_nothrow_swappable : public __is_nothrow_swappable_impl< _Tp> ::type { 
# 2712
}; 
# 2719
template< class _Tp> 
# 2720
struct is_swappable : public __is_swappable_impl< _Tp> ::type { 
# 2723
static_assert((std::__is_complete_or_unbounded(__type_identity< _Tp> {})), "template argument must be a complete class or an unbounded array");
# 2725
}; 
# 2728
template< class _Tp> 
# 2729
struct is_nothrow_swappable : public __is_nothrow_swappable_impl< _Tp> ::type { 
# 2732
static_assert((std::__is_complete_or_unbounded(__type_identity< _Tp> {})), "template argument must be a complete class or an unbounded array");
# 2734
}; 
# 2738
template< class _Tp> constexpr bool 
# 2739
is_swappable_v = (is_swappable< _Tp> ::value); 
# 2743
template< class _Tp> constexpr bool 
# 2744
is_nothrow_swappable_v = (is_nothrow_swappable< _Tp> ::value); 
# 2748
namespace __swappable_with_details { 
# 2749
using std::swap;
# 2751
struct __do_is_swappable_with_impl { 
# 2753
template< class _Tp, class _Up, class 
# 2754
 = __decltype((swap(std::declval< _Tp> (), std::declval< _Up> ()))), class 
# 2756
 = __decltype((swap(std::declval< _Up> (), std::declval< _Tp> ())))> static true_type 
# 2753
__test(int); 
# 2759
template< class , class > static false_type __test(...); 
# 2761
}; 
# 2763
struct __do_is_nothrow_swappable_with_impl { 
# 2765
template< class _Tp, class _Up> static __bool_constant< noexcept(swap(std::declval< _Tp> (), std::declval< _Up> ())) && noexcept(swap(std::declval< _Up> (), std::declval< _Tp> ()))>  __test(int); 
# 2772
template< class , class > static false_type __test(...); 
# 2774
}; 
# 2776
}
# 2778
template< class _Tp, class _Up> 
# 2779
struct __is_swappable_with_impl : public __swappable_with_details::__do_is_swappable_with_impl { 
# 2782
typedef __decltype((__test< _Tp, _Up> (0))) type; 
# 2783
}; 
# 2786
template< class _Tp> 
# 2787
struct __is_swappable_with_impl< _Tp &, _Tp &>  : public __swappable_details::__do_is_swappable_impl { 
# 2790
typedef __decltype((__test< _Tp &> (0))) type; 
# 2791
}; 
# 2793
template< class _Tp, class _Up> 
# 2794
struct __is_nothrow_swappable_with_impl : public __swappable_with_details::__do_is_nothrow_swappable_with_impl { 
# 2797
typedef __decltype((__test< _Tp, _Up> (0))) type; 
# 2798
}; 
# 2801
template< class _Tp> 
# 2802
struct __is_nothrow_swappable_with_impl< _Tp &, _Tp &>  : public __swappable_details::__do_is_nothrow_swappable_impl { 
# 2805
typedef __decltype((__test< _Tp &> (0))) type; 
# 2806
}; 
# 2809
template< class _Tp, class _Up> 
# 2810
struct is_swappable_with : public __is_swappable_with_impl< _Tp, _Up> ::type { 
# 2812
}; 
# 2815
template< class _Tp, class _Up> 
# 2816
struct is_nothrow_swappable_with : public __is_nothrow_swappable_with_impl< _Tp, _Up> ::type { 
# 2818
}; 
# 2822
template< class _Tp, class _Up> constexpr bool 
# 2823
is_swappable_with_v = (is_swappable_with< _Tp, _Up> ::value); 
# 2827
template< class _Tp, class _Up> constexpr bool 
# 2828
is_nothrow_swappable_with_v = (is_nothrow_swappable_with< _Tp, _Up> ::value); 
# 2837
template< class _Result, class _Ret, bool 
# 2838
 = is_void< _Ret> ::value, class  = void> 
# 2839
struct __is_invocable_impl : public false_type { }; 
# 2842
template< class _Result, class _Ret> 
# 2843
struct __is_invocable_impl< _Result, _Ret, true, __void_t< typename _Result::type> >  : public true_type { 
# 2847
}; 
# 2849
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wctor-dtor-privacy"
# 2852
template< class _Result, class _Ret> 
# 2853
struct __is_invocable_impl< _Result, _Ret, false, __void_t< typename _Result::type> >  { 
# 2860
private: static typename _Result::type _S_get(); 
# 2862
template< class _Tp> static void _S_conv(_Tp); 
# 2866
template< class _Tp, class  = __decltype((_S_conv< _Tp> ((_S_get)())))> static true_type _S_test(int); 
# 2870
template< class _Tp> static false_type _S_test(...); 
# 2875
public: using type = __decltype((_S_test< _Ret> (1))); 
# 2876
}; 
#pragma GCC diagnostic pop
# 2879
template< class _Fn, class ..._ArgTypes> 
# 2880
struct __is_invocable : public __is_invocable_impl< __invoke_result< _Fn, _ArgTypes...> , void> ::type { 
# 2882
}; 
# 2884
template< class _Fn, class _Tp, class ..._Args> constexpr bool 
# 2885
__call_is_nt(__invoke_memfun_ref) 
# 2886
{ 
# 2887
using _Up = typename __inv_unwrap< _Tp> ::type; 
# 2888
return noexcept((std::declval< typename __inv_unwrap< _Tp> ::type> ().*std::declval< _Fn> ())(std::declval< _Args> ()...)); 
# 2890
} 
# 2892
template< class _Fn, class _Tp, class ..._Args> constexpr bool 
# 2893
__call_is_nt(__invoke_memfun_deref) 
# 2894
{ 
# 2895
return noexcept(((*std::declval< _Tp> ()).*std::declval< _Fn> ())(std::declval< _Args> ()...)); 
# 2897
} 
# 2899
template< class _Fn, class _Tp> constexpr bool 
# 2900
__call_is_nt(__invoke_memobj_ref) 
# 2901
{ 
# 2902
using _Up = typename __inv_unwrap< _Tp> ::type; 
# 2903
return noexcept((std::declval< typename __inv_unwrap< _Tp> ::type> ().*std::declval< _Fn> ())); 
# 2904
} 
# 2906
template< class _Fn, class _Tp> constexpr bool 
# 2907
__call_is_nt(__invoke_memobj_deref) 
# 2908
{ 
# 2909
return noexcept(((*std::declval< _Tp> ()).*std::declval< _Fn> ())); 
# 2910
} 
# 2912
template< class _Fn, class ..._Args> constexpr bool 
# 2913
__call_is_nt(__invoke_other) 
# 2914
{ 
# 2915
return noexcept(std::declval< _Fn> ()(std::declval< _Args> ()...)); 
# 2916
} 
# 2918
template< class _Result, class _Fn, class ..._Args> 
# 2919
struct __call_is_nothrow : public __bool_constant< std::__call_is_nt< _Fn, _Args...> (typename _Result::__invoke_type{})>  { 
# 2923
}; 
# 2925
template< class _Fn, class ..._Args> using __call_is_nothrow_ = __call_is_nothrow< __invoke_result< _Fn, _Args...> , _Fn, _Args...> ; 
# 2930
template< class _Fn, class ..._Args> 
# 2931
struct __is_nothrow_invocable : public __and_< __is_invocable< _Fn, _Args...> , __call_is_nothrow_< _Fn, _Args...> > ::type { 
# 2934
}; 
# 2936
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wctor-dtor-privacy"
struct __nonesuchbase { }; 
# 2939
struct __nonesuch : private __nonesuchbase { 
# 2940
~__nonesuch() = delete;
# 2941
__nonesuch(const __nonesuch &) = delete;
# 2942
void operator=(const __nonesuch &) = delete;
# 2943
}; 
#pragma GCC diagnostic pop
# 3462 "/usr/include/c++/10/type_traits" 3
}
# 59 "/usr/include/c++/10/bits/move.h" 3
namespace std __attribute((__visibility__("default"))) { 
# 74 "/usr/include/c++/10/bits/move.h" 3
template< class _Tp> constexpr _Tp &&
# 76
forward(typename remove_reference< _Tp> ::type &__t) noexcept 
# 77
{ return static_cast< _Tp &&>(__t); } 
# 85
template< class _Tp> constexpr _Tp &&
# 87
forward(typename remove_reference< _Tp> ::type &&__t) noexcept 
# 88
{ 
# 89
static_assert((!std::template is_lvalue_reference< _Tp> ::value), "std::forward must not be used to convert an rvalue to an lvalue");
# 91
return static_cast< _Tp &&>(__t); 
# 92
} 
# 99
template< class _Tp> constexpr typename remove_reference< _Tp> ::type &&
# 101
move(_Tp &&__t) noexcept 
# 102
{ return static_cast< typename remove_reference< _Tp> ::type &&>(__t); } 
# 105
template< class _Tp> 
# 106
struct __move_if_noexcept_cond : public __and_< __not_< is_nothrow_move_constructible< _Tp> > , is_copy_constructible< _Tp> > ::type { 
# 108
}; 
# 118 "/usr/include/c++/10/bits/move.h" 3
template< class _Tp> constexpr typename conditional< __move_if_noexcept_cond< _Tp> ::value, const _Tp &, _Tp &&> ::type 
# 121
move_if_noexcept(_Tp &__x) noexcept 
# 122
{ return std::move(__x); } 
# 138 "/usr/include/c++/10/bits/move.h" 3
template< class _Tp> inline _Tp *
# 140
addressof(_Tp &__r) noexcept 
# 141
{ return std::__addressof(__r); } 
# 145
template < typename _Tp >
    const _Tp * addressof ( const _Tp && ) = delete;
# 149
template< class _Tp, class _Up = _Tp> inline _Tp 
# 152
__exchange(_Tp &__obj, _Up &&__new_val) 
# 153
{ 
# 154
_Tp __old_val = std::move(__obj); 
# 155
__obj = std::forward< _Up> (__new_val); 
# 156
return __old_val; 
# 157
} 
# 179 "/usr/include/c++/10/bits/move.h" 3
template< class _Tp> inline typename enable_if< __and_< __not_< __is_tuple_like< _Tp> > , is_move_constructible< _Tp> , is_move_assignable< _Tp> > ::value> ::type 
# 189
swap(_Tp &__a, _Tp &__b) noexcept(__and_< is_nothrow_move_constructible< _Tp> , is_nothrow_move_assignable< _Tp> > ::value) 
# 192
{ 
# 197
_Tp __tmp = std::move(__a); 
# 198
__a = std::move(__b); 
# 199
__b = std::move(__tmp); 
# 200
} 
# 205
template< class _Tp, size_t _Nm> inline typename enable_if< __is_swappable< _Tp> ::value> ::type 
# 213
swap(_Tp (&__a)[_Nm], _Tp (&__b)[_Nm]) noexcept(__is_nothrow_swappable< _Tp> ::value) 
# 215
{ 
# 216
for (size_t __n = (0); __n < _Nm; ++__n) { 
# 217
swap(__a[__n], __b[__n]); }  
# 218
} 
# 222
}
# 69 "/usr/include/c++/10/bits/stl_pair.h" 3
namespace std __attribute((__visibility__("default"))) { 
# 80 "/usr/include/c++/10/bits/stl_pair.h" 3
struct piecewise_construct_t { explicit piecewise_construct_t() = default;}; 
# 83
constexpr piecewise_construct_t piecewise_construct = piecewise_construct_t(); 
# 89
template< class ...> class tuple; 
# 92
template< size_t ...> struct _Index_tuple; 
# 100
template< bool , class _T1, class _T2> 
# 101
struct _PCC { 
# 103
template< class _U1, class _U2> static constexpr bool 
# 104
_ConstructiblePair() 
# 105
{ 
# 106
return __and_< is_constructible< _T1, const _U1 &> , is_constructible< _T2, const _U2 &> > ::value; 
# 108
} 
# 110
template< class _U1, class _U2> static constexpr bool 
# 111
_ImplicitlyConvertiblePair() 
# 112
{ 
# 113
return __and_< is_convertible< const _U1 &, _T1> , is_convertible< const _U2 &, _T2> > ::value; 
# 115
} 
# 117
template< class _U1, class _U2> static constexpr bool 
# 118
_MoveConstructiblePair() 
# 119
{ 
# 120
return __and_< is_constructible< _T1, _U1 &&> , is_constructible< _T2, _U2 &&> > ::value; 
# 122
} 
# 124
template< class _U1, class _U2> static constexpr bool 
# 125
_ImplicitlyMoveConvertiblePair() 
# 126
{ 
# 127
return __and_< is_convertible< _U1 &&, _T1> , is_convertible< _U2 &&, _T2> > ::value; 
# 129
} 
# 131
template< bool __implicit, class _U1, class _U2> static constexpr bool 
# 132
_CopyMovePair() 
# 133
{ 
# 134
using __do_converts = __and_< is_convertible< const _U1 &, _T1> , is_convertible< _U2 &&, _T2> > ; 
# 136
using __converts = typename conditional< __implicit, __and_< is_convertible< const _U1 &, _T1> , is_convertible< _U2 &&, _T2> > , __not_< __and_< is_convertible< const _U1 &, _T1> , is_convertible< _U2 &&, _T2> > > > ::type; 
# 139
return __and_< is_constructible< _T1, const _U1 &> , is_constructible< _T2, _U2 &&> , typename conditional< __implicit, __and_< is_convertible< const _U1 &, _T1> , is_convertible< _U2 &&, _T2> > , __not_< __and_< is_convertible< const _U1 &, _T1> , is_convertible< _U2 &&, _T2> > > > ::type> ::value; 
# 143
} 
# 145
template< bool __implicit, class _U1, class _U2> static constexpr bool 
# 146
_MoveCopyPair() 
# 147
{ 
# 148
using __do_converts = __and_< is_convertible< _U1 &&, _T1> , is_convertible< const _U2 &, _T2> > ; 
# 150
using __converts = typename conditional< __implicit, __and_< is_convertible< _U1 &&, _T1> , is_convertible< const _U2 &, _T2> > , __not_< __and_< is_convertible< _U1 &&, _T1> , is_convertible< const _U2 &, _T2> > > > ::type; 
# 153
return __and_< is_constructible< _T1, _U1 &&> , is_constructible< _T2, const _U2 &&> , typename conditional< __implicit, __and_< is_convertible< _U1 &&, _T1> , is_convertible< const _U2 &, _T2> > , __not_< __and_< is_convertible< _U1 &&, _T1> , is_convertible< const _U2 &, _T2> > > > ::type> ::value; 
# 157
} 
# 158
}; 
# 160
template< class _T1, class _T2> 
# 161
struct _PCC< false, _T1, _T2>  { 
# 163
template< class _U1, class _U2> static constexpr bool 
# 164
_ConstructiblePair() 
# 165
{ 
# 166
return false; 
# 167
} 
# 169
template< class _U1, class _U2> static constexpr bool 
# 170
_ImplicitlyConvertiblePair() 
# 171
{ 
# 172
return false; 
# 173
} 
# 175
template< class _U1, class _U2> static constexpr bool 
# 176
_MoveConstructiblePair() 
# 177
{ 
# 178
return false; 
# 179
} 
# 181
template< class _U1, class _U2> static constexpr bool 
# 182
_ImplicitlyMoveConvertiblePair() 
# 183
{ 
# 184
return false; 
# 185
} 
# 186
}; 
# 189
template< class _U1, class _U2> class __pair_base { 
# 192
template< class _T1, class _T2> friend struct pair; 
# 193
__pair_base() = default;
# 194
~__pair_base() = default;
# 195
__pair_base(const __pair_base &) = default;
# 196
__pair_base &operator=(const __pair_base &) = delete;
# 198
}; 
# 210 "/usr/include/c++/10/bits/stl_pair.h" 3
template< class _T1, class _T2> 
# 211
struct pair : private __pair_base< _T1, _T2>  { 
# 214
typedef _T1 first_type; 
# 215
typedef _T2 second_type; 
# 217
_T1 first; 
# 218
_T2 second; 
# 225
template< class _U1 = _T1, class 
# 226
_U2 = _T2, typename enable_if< __and_< __is_implicitly_default_constructible< _U1> , __is_implicitly_default_constructible< _U2> > ::value, bool> ::type 
# 230
 = true> constexpr 
# 232
pair() : first(), second() 
# 233
{ } 
# 236
template< class _U1 = _T1, class 
# 237
_U2 = _T2, typename enable_if< __and_< is_default_constructible< _U1> , is_default_constructible< _U2> , __not_< __and_< __is_implicitly_default_constructible< _U1> , __is_implicitly_default_constructible< _U2> > > > ::value, bool> ::type 
# 244
 = false> constexpr explicit 
# 245
pair() : first(), second() 
# 246
{ } 
# 256 "/usr/include/c++/10/bits/stl_pair.h" 3
using _PCCP = _PCC< true, _T1, _T2> ; 
# 260
template< class _U1 = _T1, class _U2 = _T2, typename enable_if< _PCC< true, _T1, _T2> ::template _ConstructiblePair< _U1, _U2> () && _PCC< true, _T1, _T2> ::template _ImplicitlyConvertiblePair< _U1, _U2> (), bool> ::type 
# 265
 = true> constexpr 
# 266
pair(const _T1 &__a, const _T2 &__b) : first(__a), second(__b) 
# 267
{ } 
# 270
template< class _U1 = _T1, class _U2 = _T2, typename enable_if< _PCC< true, _T1, _T2> ::template _ConstructiblePair< _U1, _U2> () && (!_PCC< true, _T1, _T2> ::template _ImplicitlyConvertiblePair< _U1, _U2> ()), bool> ::type 
# 275
 = false> constexpr explicit 
# 276
pair(const _T1 &__a, const _T2 &__b) : first(__a), second(__b) 
# 277
{ } 
# 288 "/usr/include/c++/10/bits/stl_pair.h" 3
template< class _U1, class _U2> using _PCCFP = _PCC< (!is_same< _T1, _U1> ::value) || (!is_same< _T2, _U2> ::value), _T1, _T2> ; 
# 294
template< class _U1, class _U2, typename enable_if< _PCC< (!is_same< _T1, _U1> ::value) || (!is_same< _T2, _U2> ::value), _T1, _T2> ::template _ConstructiblePair< _U1, _U2> () && _PCC< (!is_same< _T1, _U1> ::value) || (!is_same< _T2, _U2> ::value), _T1, _T2> ::template _ImplicitlyConvertiblePair< _U1, _U2> (), bool> ::type 
# 299
 = true> constexpr 
# 300
pair(const std::pair< _U1, _U2>  &__p) : first((__p.first)), second((__p.second)) 
# 301
{ } 
# 303
template< class _U1, class _U2, typename enable_if< _PCC< (!is_same< _T1, _U1> ::value) || (!is_same< _T2, _U2> ::value), _T1, _T2> ::template _ConstructiblePair< _U1, _U2> () && (!_PCC< (!is_same< _T1, _U1> ::value) || (!is_same< _T2, _U2> ::value), _T1, _T2> ::template _ImplicitlyConvertiblePair< _U1, _U2> ()), bool> ::type 
# 308
 = false> constexpr explicit 
# 309
pair(const std::pair< _U1, _U2>  &__p) : first((__p.first)), second((__p.second)) 
# 310
{ } 
# 314
constexpr pair(const pair &) = default;
# 315
constexpr pair(pair &&) = default;
# 318
template< class _U1, typename enable_if< _PCC< true, _T1, _T2> ::template _MoveCopyPair< true, _U1, _T2> (), bool> ::type 
# 321
 = true> constexpr 
# 322
pair(_U1 &&__x, const _T2 &__y) : first(std::forward< _U1> (__x)), second(__y) 
# 323
{ } 
# 325
template< class _U1, typename enable_if< _PCC< true, _T1, _T2> ::template _MoveCopyPair< false, _U1, _T2> (), bool> ::type 
# 328
 = false> constexpr explicit 
# 329
pair(_U1 &&__x, const _T2 &__y) : first(std::forward< _U1> (__x)), second(__y) 
# 330
{ } 
# 332
template< class _U2, typename enable_if< _PCC< true, _T1, _T2> ::template _CopyMovePair< true, _T1, _U2> (), bool> ::type 
# 335
 = true> constexpr 
# 336
pair(const _T1 &__x, _U2 &&__y) : first(__x), second(std::forward< _U2> (__y)) 
# 337
{ } 
# 339
template< class _U2, typename enable_if< _PCC< true, _T1, _T2> ::template _CopyMovePair< false, _T1, _U2> (), bool> ::type 
# 342
 = false> explicit 
# 343
pair(const _T1 &__x, _U2 &&__y) : first(__x), second(std::forward< _U2> (__y)) 
# 344
{ } 
# 346
template< class _U1, class _U2, typename enable_if< _PCC< true, _T1, _T2> ::template _MoveConstructiblePair< _U1, _U2> () && _PCC< true, _T1, _T2> ::template _ImplicitlyMoveConvertiblePair< _U1, _U2> (), bool> ::type 
# 351
 = true> constexpr 
# 352
pair(_U1 &&__x, _U2 &&__y) : first(std::forward< _U1> (__x)), second(std::forward< _U2> (__y)) 
# 353
{ } 
# 355
template< class _U1, class _U2, typename enable_if< _PCC< true, _T1, _T2> ::template _MoveConstructiblePair< _U1, _U2> () && (!_PCC< true, _T1, _T2> ::template _ImplicitlyMoveConvertiblePair< _U1, _U2> ()), bool> ::type 
# 360
 = false> constexpr explicit 
# 361
pair(_U1 &&__x, _U2 &&__y) : first(std::forward< _U1> (__x)), second(std::forward< _U2> (__y)) 
# 362
{ } 
# 365
template< class _U1, class _U2, typename enable_if< _PCC< (!is_same< _T1, _U1> ::value) || (!is_same< _T2, _U2> ::value), _T1, _T2> ::template _MoveConstructiblePair< _U1, _U2> () && _PCC< (!is_same< _T1, _U1> ::value) || (!is_same< _T2, _U2> ::value), _T1, _T2> ::template _ImplicitlyMoveConvertiblePair< _U1, _U2> (), bool> ::type 
# 370
 = true> constexpr 
# 371
pair(std::pair< _U1, _U2>  &&__p) : first(std::forward< _U1> ((__p.first))), second(std::forward< _U2> ((__p.second))) 
# 373
{ } 
# 375
template< class _U1, class _U2, typename enable_if< _PCC< (!is_same< _T1, _U1> ::value) || (!is_same< _T2, _U2> ::value), _T1, _T2> ::template _MoveConstructiblePair< _U1, _U2> () && (!_PCC< (!is_same< _T1, _U1> ::value) || (!is_same< _T2, _U2> ::value), _T1, _T2> ::template _ImplicitlyMoveConvertiblePair< _U1, _U2> ()), bool> ::type 
# 380
 = false> constexpr explicit 
# 381
pair(std::pair< _U1, _U2>  &&__p) : first(std::forward< _U1> ((__p.first))), second(std::forward< _U2> ((__p.second))) 
# 383
{ } 
# 385
template< class ..._Args1, class ..._Args2> pair(std::piecewise_construct_t, tuple< _Args1...> , tuple< _Args2...> ); 
# 390
pair &operator=(typename conditional< __and_< is_copy_assignable< _T1> , is_copy_assignable< _T2> > ::value, const pair &, const std::__nonesuch &> ::type 
# 393
__p) 
# 394
{ 
# 395
(first) = (__p.first); 
# 396
(second) = (__p.second); 
# 397
return *this; 
# 398
} 
# 401
pair &operator=(typename conditional< __and_< is_move_assignable< _T1> , is_move_assignable< _T2> > ::value, pair &&, std::__nonesuch &&> ::type 
# 404
__p) noexcept(__and_< is_nothrow_move_assignable< _T1> , is_nothrow_move_assignable< _T2> > ::value) 
# 407
{ 
# 408
(first) = std::forward< first_type> ((__p.first)); 
# 409
(second) = std::forward< second_type> ((__p.second)); 
# 410
return *this; 
# 411
} 
# 413
template< class _U1, class _U2> typename enable_if< __and_< is_assignable< _T1 &, const _U1 &> , is_assignable< _T2 &, const _U2 &> > ::value, pair &> ::type 
# 418
operator=(const std::pair< _U1, _U2>  &__p) 
# 419
{ 
# 420
(first) = (__p.first); 
# 421
(second) = (__p.second); 
# 422
return *this; 
# 423
} 
# 425
template< class _U1, class _U2> typename enable_if< __and_< is_assignable< _T1 &, _U1 &&> , is_assignable< _T2 &, _U2 &&> > ::value, pair &> ::type 
# 430
operator=(std::pair< _U1, _U2>  &&__p) 
# 431
{ 
# 432
(first) = std::forward< _U1> ((__p.first)); 
# 433
(second) = std::forward< _U2> ((__p.second)); 
# 434
return *this; 
# 435
} 
# 439
void swap(pair &__p) noexcept(__and_< __is_nothrow_swappable< _T1> , __is_nothrow_swappable< _T2> > ::value) 
# 442
{ 
# 443
using std::swap;
# 444
swap(first, __p.first); 
# 445
swap(second, __p.second); 
# 446
} 
# 449
private: template< class ..._Args1, std::size_t ..._Indexes1, class ...
# 450
_Args2, std::size_t ..._Indexes2> 
# 449
pair(tuple< _Args1...>  &, tuple< _Args2...>  &, _Index_tuple< _Indexes1...> , _Index_tuple< _Indexes2...> ); 
# 455
}; 
# 464 "/usr/include/c++/10/bits/stl_pair.h" 3
template< class _T1, class _T2> constexpr bool 
# 466
operator==(const pair< _T1, _T2>  &__x, const pair< _T1, _T2>  &__y) 
# 467
{ return ((__x.first) == (__y.first)) && ((__x.second) == (__y.second)); } 
# 487 "/usr/include/c++/10/bits/stl_pair.h" 3
template< class _T1, class _T2> constexpr bool 
# 489
operator<(const pair< _T1, _T2>  &__x, const pair< _T1, _T2>  &__y) 
# 490
{ return ((__x.first) < (__y.first)) || ((!((__y.first) < (__x.first))) && ((__x.second) < (__y.second))); 
# 491
} 
# 494
template< class _T1, class _T2> constexpr bool 
# 496
operator!=(const pair< _T1, _T2>  &__x, const pair< _T1, _T2>  &__y) 
# 497
{ return !(__x == __y); } 
# 500
template< class _T1, class _T2> constexpr bool 
# 502
operator>(const pair< _T1, _T2>  &__x, const pair< _T1, _T2>  &__y) 
# 503
{ return __y < __x; } 
# 506
template< class _T1, class _T2> constexpr bool 
# 508
operator<=(const pair< _T1, _T2>  &__x, const pair< _T1, _T2>  &__y) 
# 509
{ return !(__y < __x); } 
# 512
template< class _T1, class _T2> constexpr bool 
# 514
operator>=(const pair< _T1, _T2>  &__x, const pair< _T1, _T2>  &__y) 
# 515
{ return !(__x < __y); } 
# 524 "/usr/include/c++/10/bits/stl_pair.h" 3
template< class _T1, class _T2> inline typename enable_if< __and_< __is_swappable< _T1> , __is_swappable< _T2> > ::value> ::type 
# 533
swap(pair< _T1, _T2>  &__x, pair< _T1, _T2>  &__y) noexcept(noexcept(__x.swap(__y))) 
# 535
{ __x.swap(__y); } 
# 538
template < typename _T1, typename _T2 >
    typename enable_if < ! __and_ < __is_swappable < _T1 >,
          __is_swappable < _T2 > > :: value > :: type
    swap ( pair < _T1, _T2 > &, pair < _T1, _T2 > & ) = delete;
# 564 "/usr/include/c++/10/bits/stl_pair.h" 3
template< class _T1, class _T2> constexpr pair< typename __decay_and_strip< _T1> ::__type, typename __decay_and_strip< _T2> ::__type>  
# 567
make_pair(_T1 &&__x, _T2 &&__y) 
# 568
{ 
# 569
typedef typename __decay_and_strip< _T1> ::__type __ds_type1; 
# 570
typedef typename __decay_and_strip< _T2> ::__type __ds_type2; 
# 571
typedef pair< typename __decay_and_strip< _T1> ::__type, typename __decay_and_strip< _T2> ::__type>  __pair_type; 
# 572
return __pair_type(std::forward< _T1> (__x), std::forward< _T2> (__y)); 
# 573
} 
# 584 "/usr/include/c++/10/bits/stl_pair.h" 3
}
# 39 "/usr/include/c++/10/initializer_list" 3
#pragma GCC visibility push ( default )
# 43
namespace std { 
# 46
template< class _E> 
# 47
class initializer_list { 
# 50
public: typedef _E value_type; 
# 51
typedef const _E &reference; 
# 52
typedef const _E &const_reference; 
# 53
typedef size_t size_type; 
# 54
typedef const _E *iterator; 
# 55
typedef const _E *const_iterator; 
# 58
private: iterator _M_array; 
# 59
size_type _M_len; 
# 62
constexpr initializer_list(const_iterator __a, size_type __l) : _M_array(__a), _M_len(__l) 
# 63
{ } 
# 66
public: constexpr initializer_list() noexcept : _M_array((0)), _M_len((0)) 
# 67
{ } 
# 71
constexpr size_type size() const noexcept { return _M_len; } 
# 75
constexpr const_iterator begin() const noexcept { return _M_array; } 
# 79
constexpr const_iterator end() const noexcept { return begin() + size(); } 
# 80
}; 
# 88
template< class _Tp> constexpr const _Tp *
# 90
begin(initializer_list< _Tp>  __ils) noexcept 
# 91
{ return __ils.begin(); } 
# 99
template< class _Tp> constexpr const _Tp *
# 101
end(initializer_list< _Tp>  __ils) noexcept 
# 102
{ return __ils.end(); } 
# 103
}
# 105
#pragma GCC visibility pop
# 82 "/usr/include/c++/10/utility" 3
namespace std __attribute((__visibility__("default"))) { 
# 87
template< class _Tp> struct tuple_size; 
# 94
template< class _Tp, class 
# 95
_Up = typename remove_cv< _Tp> ::type, class 
# 96
 = typename enable_if< is_same< _Tp, _Up> ::value> ::type, size_t 
# 97
 = tuple_size< _Tp> ::value> using __enable_if_has_tuple_size = _Tp; 
# 100
template< class _Tp> 
# 101
struct tuple_size< const __enable_if_has_tuple_size< _Tp> >  : public std::tuple_size< _Tp>  { 
# 102
}; 
# 104
template< class _Tp> 
# 105
struct tuple_size< volatile __enable_if_has_tuple_size< _Tp> >  : public std::tuple_size< _Tp>  { 
# 106
}; 
# 108
template< class _Tp> 
# 109
struct tuple_size< const volatile __enable_if_has_tuple_size< _Tp> >  : public std::tuple_size< _Tp>  { 
# 110
}; 
# 113
template< size_t __i, class _Tp> struct tuple_element; 
# 117
template< size_t __i, class _Tp> using __tuple_element_t = typename tuple_element< __i, _Tp> ::type; 
# 120
template< size_t __i, class _Tp> 
# 121
struct tuple_element< __i, const _Tp>  { 
# 123
typedef typename add_const< __tuple_element_t< __i, _Tp> > ::type type; 
# 124
}; 
# 126
template< size_t __i, class _Tp> 
# 127
struct tuple_element< __i, volatile _Tp>  { 
# 129
typedef typename add_volatile< __tuple_element_t< __i, _Tp> > ::type type; 
# 130
}; 
# 132
template< size_t __i, class _Tp> 
# 133
struct tuple_element< __i, const volatile _Tp>  { 
# 135
typedef typename add_cv< __tuple_element_t< __i, _Tp> > ::type type; 
# 136
}; 
# 144
template< size_t __i, class _Tp> using tuple_element_t = typename tuple_element< __i, _Tp> ::type; 
# 151
template< class _T1, class _T2> 
# 152
struct __is_tuple_like_impl< pair< _T1, _T2> >  : public true_type { 
# 153
}; 
# 156
template< class _Tp1, class _Tp2> 
# 157
struct tuple_size< pair< _Tp1, _Tp2> >  : public integral_constant< unsigned long, 2UL>  { 
# 158
}; 
# 161
template< class _Tp1, class _Tp2> 
# 162
struct tuple_element< 0, pair< _Tp1, _Tp2> >  { 
# 163
typedef _Tp1 type; }; 
# 166
template< class _Tp1, class _Tp2> 
# 167
struct tuple_element< 1, pair< _Tp1, _Tp2> >  { 
# 168
typedef _Tp2 type; }; 
# 170
template< size_t _Int> struct __pair_get; 
# 174
template<> struct __pair_get< 0UL>  { 
# 176
template< class _Tp1, class _Tp2> static constexpr _Tp1 &
# 178
__get(pair< _Tp1, _Tp2>  &__pair) noexcept 
# 179
{ return __pair.first; } 
# 181
template< class _Tp1, class _Tp2> static constexpr _Tp1 &&
# 183
__move_get(pair< _Tp1, _Tp2>  &&__pair) noexcept 
# 184
{ return std::forward< _Tp1> ((__pair.first)); } 
# 186
template< class _Tp1, class _Tp2> static constexpr const _Tp1 &
# 188
__const_get(const pair< _Tp1, _Tp2>  &__pair) noexcept 
# 189
{ return __pair.first; } 
# 191
template< class _Tp1, class _Tp2> static constexpr const _Tp1 &&
# 193
__const_move_get(const pair< _Tp1, _Tp2>  &&__pair) noexcept 
# 194
{ return std::forward< const _Tp1> ((__pair.first)); } 
# 195
}; 
# 198
template<> struct __pair_get< 1UL>  { 
# 200
template< class _Tp1, class _Tp2> static constexpr _Tp2 &
# 202
__get(pair< _Tp1, _Tp2>  &__pair) noexcept 
# 203
{ return __pair.second; } 
# 205
template< class _Tp1, class _Tp2> static constexpr _Tp2 &&
# 207
__move_get(pair< _Tp1, _Tp2>  &&__pair) noexcept 
# 208
{ return std::forward< _Tp2> ((__pair.second)); } 
# 210
template< class _Tp1, class _Tp2> static constexpr const _Tp2 &
# 212
__const_get(const pair< _Tp1, _Tp2>  &__pair) noexcept 
# 213
{ return __pair.second; } 
# 215
template< class _Tp1, class _Tp2> static constexpr const _Tp2 &&
# 217
__const_move_get(const pair< _Tp1, _Tp2>  &&__pair) noexcept 
# 218
{ return std::forward< const _Tp2> ((__pair.second)); } 
# 219
}; 
# 221
template< size_t _Int, class _Tp1, class _Tp2> constexpr typename tuple_element< _Int, pair< _Tp1, _Tp2> > ::type &
# 223
get(pair< _Tp1, _Tp2>  &__in) noexcept 
# 224
{ return __pair_get< _Int> ::__get(__in); } 
# 226
template< size_t _Int, class _Tp1, class _Tp2> constexpr typename tuple_element< _Int, pair< _Tp1, _Tp2> > ::type &&
# 228
get(pair< _Tp1, _Tp2>  &&__in) noexcept 
# 229
{ return __pair_get< _Int> ::__move_get(std::move(__in)); } 
# 231
template< size_t _Int, class _Tp1, class _Tp2> constexpr const typename tuple_element< _Int, pair< _Tp1, _Tp2> > ::type &
# 233
get(const pair< _Tp1, _Tp2>  &__in) noexcept 
# 234
{ return __pair_get< _Int> ::__const_get(__in); } 
# 236
template< size_t _Int, class _Tp1, class _Tp2> constexpr const typename tuple_element< _Int, pair< _Tp1, _Tp2> > ::type &&
# 238
get(const pair< _Tp1, _Tp2>  &&__in) noexcept 
# 239
{ return __pair_get< _Int> ::__const_move_get(std::move(__in)); } 
# 245
template< class _Tp, class _Up> constexpr _Tp &
# 247
get(pair< _Tp, _Up>  &__p) noexcept 
# 248
{ return __p.first; } 
# 250
template< class _Tp, class _Up> constexpr const _Tp &
# 252
get(const pair< _Tp, _Up>  &__p) noexcept 
# 253
{ return __p.first; } 
# 255
template< class _Tp, class _Up> constexpr _Tp &&
# 257
get(pair< _Tp, _Up>  &&__p) noexcept 
# 258
{ return std::move((__p.first)); } 
# 260
template< class _Tp, class _Up> constexpr const _Tp &&
# 262
get(const pair< _Tp, _Up>  &&__p) noexcept 
# 263
{ return std::move((__p.first)); } 
# 265
template< class _Tp, class _Up> constexpr _Tp &
# 267
get(pair< _Up, _Tp>  &__p) noexcept 
# 268
{ return __p.second; } 
# 270
template< class _Tp, class _Up> constexpr const _Tp &
# 272
get(const pair< _Up, _Tp>  &__p) noexcept 
# 273
{ return __p.second; } 
# 275
template< class _Tp, class _Up> constexpr _Tp &&
# 277
get(pair< _Up, _Tp>  &&__p) noexcept 
# 278
{ return std::move((__p.second)); } 
# 280
template< class _Tp, class _Up> constexpr const _Tp &&
# 282
get(const pair< _Up, _Tp>  &&__p) noexcept 
# 283
{ return std::move((__p.second)); } 
# 288
template< class _Tp, class _Up = _Tp> inline _Tp 
# 291
exchange(_Tp &__obj, _Up &&__new_val) 
# 292
{ return std::__exchange(__obj, std::forward< _Up> (__new_val)); } 
# 298
template< size_t ..._Indexes> struct _Index_tuple { }; 
# 307 "/usr/include/c++/10/utility" 3
template< size_t _Num> 
# 308
struct _Build_index_tuple { 
# 316
using __type = _Index_tuple< __integer_pack(_Num)...> ; 
# 318
}; 
# 325
template< class _Tp, _Tp ..._Idx> 
# 326
struct integer_sequence { 
# 328
typedef _Tp value_type; 
# 329
static constexpr size_t size() noexcept { return sizeof...(_Idx); } 
# 330
}; 
# 333
template< class _Tp, _Tp _Num> using make_integer_sequence = integer_sequence< _Tp, __integer_pack(_Num)...> ; 
# 344
template< size_t ..._Idx> using index_sequence = integer_sequence< unsigned long, _Idx...> ; 
# 348
template< size_t _Num> using make_index_sequence = make_integer_sequence< unsigned long, _Num> ; 
# 352
template< class ..._Types> using index_sequence_for = make_index_sequence< sizeof...(_Types)> ; 
# 474 "/usr/include/c++/10/utility" 3
}
# 206 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime.h"
template< class T> static inline cudaError_t 
# 207
cudaLaunchKernel(const T *
# 208
func, dim3 
# 209
gridDim, dim3 
# 210
blockDim, void **
# 211
args, size_t 
# 212
sharedMem = 0, cudaStream_t 
# 213
stream = 0) 
# 215
{ 
# 216
return ::cudaLaunchKernel((const void *)func, gridDim, blockDim, args, sharedMem, stream); 
# 217
} 
# 277 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime.h"
template< class ...ExpTypes, class ...ActTypes> static inline cudaError_t 
# 278
cudaLaunchKernelEx(const cudaLaunchConfig_t *
# 279
config, void (*
# 280
kernel)(ExpTypes ...), ActTypes &&...
# 281
args) 
# 283
{ 
# 284
return [&](ExpTypes ...coercedArgs) { 
# 285
void *pArgs[] = {(&coercedArgs)...}; 
# 286
return ::cudaLaunchKernelExC(config, (const void *)(kernel), pArgs); 
# 287
} (std::forward< ActTypes> (args)...); 
# 288
} 
# 340 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime.h"
template< class T> static inline cudaError_t 
# 341
cudaLaunchCooperativeKernel(const T *
# 342
func, dim3 
# 343
gridDim, dim3 
# 344
blockDim, void **
# 345
args, size_t 
# 346
sharedMem = 0, cudaStream_t 
# 347
stream = 0) 
# 349
{ 
# 350
return ::cudaLaunchCooperativeKernel((const void *)func, gridDim, blockDim, args, sharedMem, stream); 
# 351
} 
# 384 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime.h"
static inline cudaError_t cudaEventCreate(cudaEvent_t *
# 385
event, unsigned 
# 386
flags) 
# 388
{ 
# 389
return ::cudaEventCreateWithFlags(event, flags); 
# 390
} 
# 428 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime.h"
static inline cudaError_t cudaGraphInstantiate(cudaGraphExec_t *
# 429
pGraphExec, cudaGraph_t 
# 430
graph, cudaGraphNode_t *
# 431
pErrorNode, char *
# 432
pLogBuffer, size_t 
# 433
bufferSize) 
# 435
{ 
# 436
(void)pErrorNode; 
# 437
(void)pLogBuffer; 
# 438
(void)bufferSize; 
# 439
return ::cudaGraphInstantiate(pGraphExec, graph, 0); 
# 440
} 
# 499 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime.h"
static inline cudaError_t cudaMallocHost(void **
# 500
ptr, size_t 
# 501
size, unsigned 
# 502
flags) 
# 504
{ 
# 505
return ::cudaHostAlloc(ptr, size, flags); 
# 506
} 
# 508
template< class T> static inline cudaError_t 
# 509
cudaHostAlloc(T **
# 510
ptr, size_t 
# 511
size, unsigned 
# 512
flags) 
# 514
{ 
# 515
return ::cudaHostAlloc((void **)((void *)ptr), size, flags); 
# 516
} 
# 518
template< class T> static inline cudaError_t 
# 519
cudaHostGetDevicePointer(T **
# 520
pDevice, void *
# 521
pHost, unsigned 
# 522
flags) 
# 524
{ 
# 525
return ::cudaHostGetDevicePointer((void **)((void *)pDevice), pHost, flags); 
# 526
} 
# 628 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime.h"
template< class T> static inline cudaError_t 
# 629
cudaMallocManaged(T **
# 630
devPtr, size_t 
# 631
size, unsigned 
# 632
flags = 1) 
# 634
{ 
# 635
return ::cudaMallocManaged((void **)((void *)devPtr), size, flags); 
# 636
} 
# 646 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime.h"
template< class T> cudaError_t 
# 647
cudaMemAdvise(T *
# 648
devPtr, size_t 
# 649
count, cudaMemoryAdvise 
# 650
advice, cudaMemLocation 
# 651
location) 
# 653
{ 
# 654
return ::cudaMemAdvise_v2((const void *)devPtr, count, advice, location); 
# 655
} 
# 657
template< class T> static inline cudaError_t 
# 658
cudaMemPrefetchAsync(T *
# 659
devPtr, size_t 
# 660
count, cudaMemLocation 
# 661
location, unsigned 
# 662
flags, cudaStream_t 
# 663
stream = 0) 
# 665
{ 
# 666
return ::cudaMemPrefetchAsync_v2((const void *)devPtr, count, location, flags, stream); 
# 667
} 
# 749 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime.h"
template< class T> static inline cudaError_t 
# 750
cudaStreamAttachMemAsync(cudaStream_t 
# 751
stream, T *
# 752
devPtr, size_t 
# 753
length = 0, unsigned 
# 754
flags = 4) 
# 756
{ 
# 757
return ::cudaStreamAttachMemAsync(stream, (void *)devPtr, length, flags); 
# 758
} 
# 760
template< class T> inline cudaError_t 
# 761
cudaMalloc(T **
# 762
devPtr, size_t 
# 763
size) 
# 765
{ 
# 766
return ::cudaMalloc((void **)((void *)devPtr), size); 
# 767
} 
# 769
template< class T> static inline cudaError_t 
# 770
cudaMallocHost(T **
# 771
ptr, size_t 
# 772
size, unsigned 
# 773
flags = 0) 
# 775
{ 
# 776
return cudaMallocHost((void **)((void *)ptr), size, flags); 
# 777
} 
# 779
template< class T> static inline cudaError_t 
# 780
cudaMallocPitch(T **
# 781
devPtr, size_t *
# 782
pitch, size_t 
# 783
width, size_t 
# 784
height) 
# 786
{ 
# 787
return ::cudaMallocPitch((void **)((void *)devPtr), pitch, width, height); 
# 788
} 
# 799 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime.h"
static inline cudaError_t cudaMallocAsync(void **
# 800
ptr, size_t 
# 801
size, cudaMemPool_t 
# 802
memPool, cudaStream_t 
# 803
stream) 
# 805
{ 
# 806
return ::cudaMallocFromPoolAsync(ptr, size, memPool, stream); 
# 807
} 
# 809
template< class T> static inline cudaError_t 
# 810
cudaMallocAsync(T **
# 811
ptr, size_t 
# 812
size, cudaMemPool_t 
# 813
memPool, cudaStream_t 
# 814
stream) 
# 816
{ 
# 817
return ::cudaMallocFromPoolAsync((void **)((void *)ptr), size, memPool, stream); 
# 818
} 
# 820
template< class T> static inline cudaError_t 
# 821
cudaMallocAsync(T **
# 822
ptr, size_t 
# 823
size, cudaStream_t 
# 824
stream) 
# 826
{ 
# 827
return ::cudaMallocAsync((void **)((void *)ptr), size, stream); 
# 828
} 
# 830
template< class T> static inline cudaError_t 
# 831
cudaMallocFromPoolAsync(T **
# 832
ptr, size_t 
# 833
size, cudaMemPool_t 
# 834
memPool, cudaStream_t 
# 835
stream) 
# 837
{ 
# 838
return ::cudaMallocFromPoolAsync((void **)((void *)ptr), size, memPool, stream); 
# 839
} 
# 878 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime.h"
template< class T> static inline cudaError_t 
# 879
cudaMemcpyToSymbol(const T &
# 880
symbol, const void *
# 881
src, size_t 
# 882
count, size_t 
# 883
offset = 0, cudaMemcpyKind 
# 884
kind = cudaMemcpyHostToDevice) 
# 886
{ 
# 887
return ::cudaMemcpyToSymbol((const void *)(&symbol), src, count, offset, kind); 
# 888
} 
# 932 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime.h"
template< class T> static inline cudaError_t 
# 933
cudaMemcpyToSymbolAsync(const T &
# 934
symbol, const void *
# 935
src, size_t 
# 936
count, size_t 
# 937
offset = 0, cudaMemcpyKind 
# 938
kind = cudaMemcpyHostToDevice, cudaStream_t 
# 939
stream = 0) 
# 941
{ 
# 942
return ::cudaMemcpyToSymbolAsync((const void *)(&symbol), src, count, offset, kind, stream); 
# 943
} 
# 980 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime.h"
template< class T> static inline cudaError_t 
# 981
cudaMemcpyFromSymbol(void *
# 982
dst, const T &
# 983
symbol, size_t 
# 984
count, size_t 
# 985
offset = 0, cudaMemcpyKind 
# 986
kind = cudaMemcpyDeviceToHost) 
# 988
{ 
# 989
return ::cudaMemcpyFromSymbol(dst, (const void *)(&symbol), count, offset, kind); 
# 990
} 
# 1034 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime.h"
template< class T> static inline cudaError_t 
# 1035
cudaMemcpyFromSymbolAsync(void *
# 1036
dst, const T &
# 1037
symbol, size_t 
# 1038
count, size_t 
# 1039
offset = 0, cudaMemcpyKind 
# 1040
kind = cudaMemcpyDeviceToHost, cudaStream_t 
# 1041
stream = 0) 
# 1043
{ 
# 1044
return ::cudaMemcpyFromSymbolAsync(dst, (const void *)(&symbol), count, offset, kind, stream); 
# 1045
} 
# 1103 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime.h"
template< class T> static inline cudaError_t 
# 1104
cudaGraphAddMemcpyNodeToSymbol(cudaGraphNode_t *
# 1105
pGraphNode, cudaGraph_t 
# 1106
graph, const cudaGraphNode_t *
# 1107
pDependencies, size_t 
# 1108
numDependencies, const T &
# 1109
symbol, const void *
# 1110
src, size_t 
# 1111
count, size_t 
# 1112
offset, cudaMemcpyKind 
# 1113
kind) 
# 1114
{ 
# 1115
return ::cudaGraphAddMemcpyNodeToSymbol(pGraphNode, graph, pDependencies, numDependencies, (const void *)(&symbol), src, count, offset, kind); 
# 1116
} 
# 1174 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime.h"
template< class T> static inline cudaError_t 
# 1175
cudaGraphAddMemcpyNodeFromSymbol(cudaGraphNode_t *
# 1176
pGraphNode, cudaGraph_t 
# 1177
graph, const cudaGraphNode_t *
# 1178
pDependencies, size_t 
# 1179
numDependencies, void *
# 1180
dst, const T &
# 1181
symbol, size_t 
# 1182
count, size_t 
# 1183
offset, cudaMemcpyKind 
# 1184
kind) 
# 1185
{ 
# 1186
return ::cudaGraphAddMemcpyNodeFromSymbol(pGraphNode, graph, pDependencies, numDependencies, dst, (const void *)(&symbol), count, offset, kind); 
# 1187
} 
# 1225 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime.h"
template< class T> static inline cudaError_t 
# 1226
cudaGraphMemcpyNodeSetParamsToSymbol(cudaGraphNode_t 
# 1227
node, const T &
# 1228
symbol, const void *
# 1229
src, size_t 
# 1230
count, size_t 
# 1231
offset, cudaMemcpyKind 
# 1232
kind) 
# 1233
{ 
# 1234
return ::cudaGraphMemcpyNodeSetParamsToSymbol(node, (const void *)(&symbol), src, count, offset, kind); 
# 1235
} 
# 1273 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime.h"
template< class T> static inline cudaError_t 
# 1274
cudaGraphMemcpyNodeSetParamsFromSymbol(cudaGraphNode_t 
# 1275
node, void *
# 1276
dst, const T &
# 1277
symbol, size_t 
# 1278
count, size_t 
# 1279
offset, cudaMemcpyKind 
# 1280
kind) 
# 1281
{ 
# 1282
return ::cudaGraphMemcpyNodeSetParamsFromSymbol(node, dst, (const void *)(&symbol), count, offset, kind); 
# 1283
} 
# 1331 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime.h"
template< class T> static inline cudaError_t 
# 1332
cudaGraphExecMemcpyNodeSetParamsToSymbol(cudaGraphExec_t 
# 1333
hGraphExec, cudaGraphNode_t 
# 1334
node, const T &
# 1335
symbol, const void *
# 1336
src, size_t 
# 1337
count, size_t 
# 1338
offset, cudaMemcpyKind 
# 1339
kind) 
# 1340
{ 
# 1341
return ::cudaGraphExecMemcpyNodeSetParamsToSymbol(hGraphExec, node, (const void *)(&symbol), src, count, offset, kind); 
# 1342
} 
# 1390 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime.h"
template< class T> static inline cudaError_t 
# 1391
cudaGraphExecMemcpyNodeSetParamsFromSymbol(cudaGraphExec_t 
# 1392
hGraphExec, cudaGraphNode_t 
# 1393
node, void *
# 1394
dst, const T &
# 1395
symbol, size_t 
# 1396
count, size_t 
# 1397
offset, cudaMemcpyKind 
# 1398
kind) 
# 1399
{ 
# 1400
return ::cudaGraphExecMemcpyNodeSetParamsFromSymbol(hGraphExec, node, dst, (const void *)(&symbol), count, offset, kind); 
# 1401
} 
# 1404
static inline cudaError_t cudaGraphExecUpdate(cudaGraphExec_t hGraphExec, cudaGraph_t hGraph, cudaGraphNode_t *hErrorNode_out, cudaGraphExecUpdateResult *updateResult_out) 
# 1405
{ 
# 1406
cudaGraphExecUpdateResultInfo resultInfo; 
# 1407
cudaError_t status = cudaGraphExecUpdate(hGraphExec, hGraph, &resultInfo); 
# 1408
if (hErrorNode_out) { 
# 1409
(*hErrorNode_out) = (resultInfo.errorNode); 
# 1410
}  
# 1411
if (updateResult_out) { 
# 1412
(*updateResult_out) = (resultInfo.result); 
# 1413
}  
# 1414
return status; 
# 1415
} 
# 1443 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime.h"
template< class T> static inline cudaError_t 
# 1444
cudaUserObjectCreate(cudaUserObject_t *
# 1445
object_out, T *
# 1446
objectToWrap, unsigned 
# 1447
initialRefcount, unsigned 
# 1448
flags) 
# 1449
{ 
# 1450
return ::cudaUserObjectCreate(object_out, objectToWrap, [](void *
# 1453
vpObj) { delete (reinterpret_cast< T *>(vpObj)); } , initialRefcount, flags); 
# 1456
} 
# 1458
template< class T> static inline cudaError_t 
# 1459
cudaUserObjectCreate(cudaUserObject_t *
# 1460
object_out, T *
# 1461
objectToWrap, unsigned 
# 1462
initialRefcount, cudaUserObjectFlags 
# 1463
flags) 
# 1464
{ 
# 1465
return cudaUserObjectCreate(object_out, objectToWrap, initialRefcount, (unsigned)flags); 
# 1466
} 
# 1493 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime.h"
template< class T> static inline cudaError_t 
# 1494
cudaGetSymbolAddress(void **
# 1495
devPtr, const T &
# 1496
symbol) 
# 1498
{ 
# 1499
return ::cudaGetSymbolAddress(devPtr, (const void *)(&symbol)); 
# 1500
} 
# 1525 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime.h"
template< class T> static inline cudaError_t 
# 1526
cudaGetSymbolSize(size_t *
# 1527
size, const T &
# 1528
symbol) 
# 1530
{ 
# 1531
return ::cudaGetSymbolSize(size, (const void *)(&symbol)); 
# 1532
} 
# 1577 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime.h"
template< class T> static inline cudaError_t 
# 1578
cudaFuncSetCacheConfig(T *
# 1579
func, cudaFuncCache 
# 1580
cacheConfig) 
# 1582
{ 
# 1583
return ::cudaFuncSetCacheConfig((const void *)func, cacheConfig); 
# 1584
} 
# 1586
template< class T> static inline cudaError_t 
# 1587
cudaFuncSetSharedMemConfig(T *
# 1588
func, cudaSharedMemConfig 
# 1589
config) 
# 1591
{ 
# 1592
return ::cudaFuncSetSharedMemConfig((const void *)func, config); 
# 1593
} 
# 1625 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime.h"
template< class T> inline cudaError_t 
# 1626
cudaOccupancyMaxActiveBlocksPerMultiprocessor(int *
# 1627
numBlocks, T 
# 1628
func, int 
# 1629
blockSize, size_t 
# 1630
dynamicSMemSize) 
# 1631
{ 
# 1632
return ::cudaOccupancyMaxActiveBlocksPerMultiprocessorWithFlags(numBlocks, (const void *)func, blockSize, dynamicSMemSize, 0); 
# 1633
} 
# 1677 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime.h"
template< class T> inline cudaError_t 
# 1678
cudaOccupancyMaxActiveBlocksPerMultiprocessorWithFlags(int *
# 1679
numBlocks, T 
# 1680
func, int 
# 1681
blockSize, size_t 
# 1682
dynamicSMemSize, unsigned 
# 1683
flags) 
# 1684
{ 
# 1685
return ::cudaOccupancyMaxActiveBlocksPerMultiprocessorWithFlags(numBlocks, (const void *)func, blockSize, dynamicSMemSize, flags); 
# 1686
} 
# 1691
class __cudaOccupancyB2DHelper { 
# 1692
size_t n; 
# 1694
public: __cudaOccupancyB2DHelper(size_t n_) : n(n_) { } 
# 1695
size_t operator()(int) 
# 1696
{ 
# 1697
return n; 
# 1698
} 
# 1699
}; 
# 1747 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime.h"
template< class UnaryFunction, class T> static inline cudaError_t 
# 1748
cudaOccupancyMaxPotentialBlockSizeVariableSMemWithFlags(int *
# 1749
minGridSize, int *
# 1750
blockSize, T 
# 1751
func, UnaryFunction 
# 1752
blockSizeToDynamicSMemSize, int 
# 1753
blockSizeLimit = 0, unsigned 
# 1754
flags = 0) 
# 1755
{ 
# 1756
cudaError_t status; 
# 1759
int device; 
# 1760
cudaFuncAttributes attr; 
# 1763
int maxThreadsPerMultiProcessor; 
# 1764
int warpSize; 
# 1765
int devMaxThreadsPerBlock; 
# 1766
int multiProcessorCount; 
# 1767
int funcMaxThreadsPerBlock; 
# 1768
int occupancyLimit; 
# 1769
int granularity; 
# 1772
int maxBlockSize = 0; 
# 1773
int numBlocks = 0; 
# 1774
int maxOccupancy = 0; 
# 1777
int blockSizeToTryAligned; 
# 1778
int blockSizeToTry; 
# 1779
int blockSizeLimitAligned; 
# 1780
int occupancyInBlocks; 
# 1781
int occupancyInThreads; 
# 1782
size_t dynamicSMemSize; 
# 1788
if (((!minGridSize) || (!blockSize)) || (!func)) { 
# 1789
return cudaErrorInvalidValue; 
# 1790
}  
# 1796
status = ::cudaGetDevice(&device); 
# 1797
if (status != (cudaSuccess)) { 
# 1798
return status; 
# 1799
}  
# 1801
status = cudaDeviceGetAttribute(&maxThreadsPerMultiProcessor, cudaDevAttrMaxThreadsPerMultiProcessor, device); 
# 1805
if (status != (cudaSuccess)) { 
# 1806
return status; 
# 1807
}  
# 1809
status = cudaDeviceGetAttribute(&warpSize, cudaDevAttrWarpSize, device); 
# 1813
if (status != (cudaSuccess)) { 
# 1814
return status; 
# 1815
}  
# 1817
status = cudaDeviceGetAttribute(&devMaxThreadsPerBlock, cudaDevAttrMaxThreadsPerBlock, device); 
# 1821
if (status != (cudaSuccess)) { 
# 1822
return status; 
# 1823
}  
# 1825
status = cudaDeviceGetAttribute(&multiProcessorCount, cudaDevAttrMultiProcessorCount, device); 
# 1829
if (status != (cudaSuccess)) { 
# 1830
return status; 
# 1831
}  
# 1833
status = cudaFuncGetAttributes(&attr, func); 
# 1834
if (status != (cudaSuccess)) { 
# 1835
return status; 
# 1836
}  
# 1838
funcMaxThreadsPerBlock = (attr.maxThreadsPerBlock); 
# 1844
occupancyLimit = maxThreadsPerMultiProcessor; 
# 1845
granularity = warpSize; 
# 1847
if (blockSizeLimit == 0) { 
# 1848
blockSizeLimit = devMaxThreadsPerBlock; 
# 1849
}  
# 1851
if (devMaxThreadsPerBlock < blockSizeLimit) { 
# 1852
blockSizeLimit = devMaxThreadsPerBlock; 
# 1853
}  
# 1855
if (funcMaxThreadsPerBlock < blockSizeLimit) { 
# 1856
blockSizeLimit = funcMaxThreadsPerBlock; 
# 1857
}  
# 1859
blockSizeLimitAligned = (((blockSizeLimit + (granularity - 1)) / granularity) * granularity); 
# 1861
for (blockSizeToTryAligned = blockSizeLimitAligned; blockSizeToTryAligned > 0; blockSizeToTryAligned -= granularity) { 
# 1865
if (blockSizeLimit < blockSizeToTryAligned) { 
# 1866
blockSizeToTry = blockSizeLimit; 
# 1867
} else { 
# 1868
blockSizeToTry = blockSizeToTryAligned; 
# 1869
}  
# 1871
dynamicSMemSize = blockSizeToDynamicSMemSize(blockSizeToTry); 
# 1873
status = cudaOccupancyMaxActiveBlocksPerMultiprocessorWithFlags(&occupancyInBlocks, func, blockSizeToTry, dynamicSMemSize, flags); 
# 1880
if (status != (cudaSuccess)) { 
# 1881
return status; 
# 1882
}  
# 1884
occupancyInThreads = (blockSizeToTry * occupancyInBlocks); 
# 1886
if (occupancyInThreads > maxOccupancy) { 
# 1887
maxBlockSize = blockSizeToTry; 
# 1888
numBlocks = occupancyInBlocks; 
# 1889
maxOccupancy = occupancyInThreads; 
# 1890
}  
# 1894
if (occupancyLimit == maxOccupancy) { 
# 1895
break; 
# 1896
}  
# 1897
}  
# 1905
(*minGridSize) = (numBlocks * multiProcessorCount); 
# 1906
(*blockSize) = maxBlockSize; 
# 1908
return status; 
# 1909
} 
# 1943 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime.h"
template< class UnaryFunction, class T> static inline cudaError_t 
# 1944
cudaOccupancyMaxPotentialBlockSizeVariableSMem(int *
# 1945
minGridSize, int *
# 1946
blockSize, T 
# 1947
func, UnaryFunction 
# 1948
blockSizeToDynamicSMemSize, int 
# 1949
blockSizeLimit = 0) 
# 1950
{ 
# 1951
return cudaOccupancyMaxPotentialBlockSizeVariableSMemWithFlags(minGridSize, blockSize, func, blockSizeToDynamicSMemSize, blockSizeLimit, 0); 
# 1952
} 
# 1989 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime.h"
template< class T> static inline cudaError_t 
# 1990
cudaOccupancyMaxPotentialBlockSize(int *
# 1991
minGridSize, int *
# 1992
blockSize, T 
# 1993
func, size_t 
# 1994
dynamicSMemSize = 0, int 
# 1995
blockSizeLimit = 0) 
# 1996
{ 
# 1997
return cudaOccupancyMaxPotentialBlockSizeVariableSMemWithFlags(minGridSize, blockSize, func, ((__cudaOccupancyB2DHelper)(dynamicSMemSize)), blockSizeLimit, 0); 
# 1998
} 
# 2027 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime.h"
template< class T> static inline cudaError_t 
# 2028
cudaOccupancyAvailableDynamicSMemPerBlock(size_t *
# 2029
dynamicSmemSize, T 
# 2030
func, int 
# 2031
numBlocks, int 
# 2032
blockSize) 
# 2033
{ 
# 2034
return ::cudaOccupancyAvailableDynamicSMemPerBlock(dynamicSmemSize, (const void *)func, numBlocks, blockSize); 
# 2035
} 
# 2086 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime.h"
template< class T> static inline cudaError_t 
# 2087
cudaOccupancyMaxPotentialBlockSizeWithFlags(int *
# 2088
minGridSize, int *
# 2089
blockSize, T 
# 2090
func, size_t 
# 2091
dynamicSMemSize = 0, int 
# 2092
blockSizeLimit = 0, unsigned 
# 2093
flags = 0) 
# 2094
{ 
# 2095
return cudaOccupancyMaxPotentialBlockSizeVariableSMemWithFlags(minGridSize, blockSize, func, ((__cudaOccupancyB2DHelper)(dynamicSMemSize)), blockSizeLimit, flags); 
# 2096
} 
# 2130 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime.h"
template< class T> static inline cudaError_t 
# 2131
cudaOccupancyMaxPotentialClusterSize(int *
# 2132
clusterSize, T *
# 2133
func, const cudaLaunchConfig_t *
# 2134
config) 
# 2135
{ 
# 2136
return ::cudaOccupancyMaxPotentialClusterSize(clusterSize, (const void *)func, config); 
# 2137
} 
# 2173 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime.h"
template< class T> static inline cudaError_t 
# 2174
cudaOccupancyMaxActiveClusters(int *
# 2175
numClusters, T *
# 2176
func, const cudaLaunchConfig_t *
# 2177
config) 
# 2178
{ 
# 2179
return ::cudaOccupancyMaxActiveClusters(numClusters, (const void *)func, config); 
# 2180
} 
# 2213 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime.h"
template< class T> inline cudaError_t 
# 2214
cudaFuncGetAttributes(cudaFuncAttributes *
# 2215
attr, T *
# 2216
entry) 
# 2218
{ 
# 2219
return ::cudaFuncGetAttributes(attr, (const void *)entry); 
# 2220
} 
# 2275 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime.h"
template< class T> static inline cudaError_t 
# 2276
cudaFuncSetAttribute(T *
# 2277
entry, cudaFuncAttribute 
# 2278
attr, int 
# 2279
value) 
# 2281
{ 
# 2282
return ::cudaFuncSetAttribute((const void *)entry, attr, value); 
# 2283
} 
# 2299 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime.h"
template< class T> static inline cudaError_t 
# 2300
cudaGetKernel(cudaKernel_t *
# 2301
kernelPtr, const T *
# 2302
entryFuncAddr) 
# 2304
{ 
# 2305
return ::cudaGetKernel(kernelPtr, (const void *)entryFuncAddr); 
# 2306
} 
# 2317 "/usr/local/cuda-12.2/bin/../targets/x86_64-linux/include/cuda_runtime.h"
#pragma GCC diagnostic pop
# 40 "CMakeCUDACompilerId.cu"
const char *info_compiler = ("INFO:compiler[NVIDIA]"); 
# 278 "CMakeCUDACompilerId.cu"
const char info_version[] = {'I', 'N', 'F', 'O', ':', 'c', 'o', 'm', 'p', 'i', 'l', 'e', 'r', '_', 'v', 'e', 'r', 's', 'i', 'o', 'n', '[', (('0') + ((12 / 10000000) % 10)), (('0') + ((12 / 1000000) % 10)), (('0') + ((12 / 100000) % 10)), (('0') + ((12 / 10000) % 10)), (('0') + ((12 / 1000) % 10)), (('0') + ((12 / 100) % 10)), (('0') + ((12 / 10) % 10)), (('0') + (12 % 10)), '.', (('0') + ((2 / 10000000) % 10)), (('0') + ((2 / 1000000) % 10)), (('0') + ((2 / 100000) % 10)), (('0') + ((2 / 10000) % 10)), (('0') + ((2 / 1000) % 10)), (('0') + ((2 / 100) % 10)), (('0') + ((2 / 10) % 10)), (('0') + (2 % 10)), '.', (('0') + ((140 / 10000000) % 10)), (('0') + ((140 / 1000000) % 10)), (('0') + ((140 / 100000) % 10)), (('0') + ((140 / 10000) % 10)), (('0') + ((140 / 1000) % 10)), (('0') + ((140 / 100) % 10)), (('0') + ((140 / 10) % 10)), (('0') + (140 % 10)), ']', '\000'}; 
# 325 "CMakeCUDACompilerId.cu"
const char *info_platform = ("INFO:platform[Linux]"); 
# 326
const char *info_arch = ("INFO:arch[]"); 
# 331
const char *info_language_dialect_default = ("INFO:dialect_default[14]"); 
# 347
int main(int argc, char *argv[]) 
# 348
{ 
# 349
int require = 0; 
# 350
require += (info_compiler[argc]); 
# 351
require += (info_platform[argc]); 
# 353
require += (info_version[argc]); 
# 361
require += (info_language_dialect_default[argc]); 
# 362
(void)argv; 
# 363
return require; 
# 364
} 

# 1 "CMakeCUDACompilerId.cudafe1.stub.c"
#define _NV_ANON_NAMESPACE _GLOBAL__N__b1dd15db_22_CMakeCUDACompilerId_cu_bd57c623
#ifdef _NV_ANON_NAMESPACE
#endif
# 1 "CMakeCUDACompilerId.cudafe1.stub.c"
#include "CMakeCUDACompilerId.cudafe1.stub.c"
# 1 "CMakeCUDACompilerId.cudafe1.stub.c"
#undef _NV_ANON_NAMESPACE
