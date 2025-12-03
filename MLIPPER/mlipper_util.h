#pragma once
#include <cuda_runtime.h>
#include <cstdint>
#include <string>
#include <stdexcept>
#include <cstddef>

// Common CUDA error checking macro used across modules.
#ifndef CHECK_CUDA
#define CHECK_CUDA(call)                                                                  \
    do {                                                                                  \
        cudaError_t err__ = (call);                                                       \
        if (err__ != cudaSuccess) {                                                       \
            throw std::runtime_error(std::string("CUDA error: ") + cudaGetErrorString(err__)); \
        }                                                                                 \
    } while (0)
#endif

// Fast integer log2 ceil for small unsigned values.
inline unsigned int ceil_log2_u32(unsigned int x) {
    if (x <= 1u) return 0u;
    unsigned int v = x - 1u;
    unsigned int r = 0u;
    while (v) { v >>= 1u; ++r; }
    return r;
}

// RAII wrapper for device buffers allocated with cudaMalloc.
struct DeviceBuffer {
    double* ptr{nullptr};

    DeviceBuffer() = default;
    explicit DeviceBuffer(std::size_t count) { CHECK_CUDA(cudaMalloc(&ptr, sizeof(double) * count)); }
    ~DeviceBuffer() { if (ptr) cudaFree(ptr); }

    DeviceBuffer(const DeviceBuffer&) = delete;
    DeviceBuffer& operator=(const DeviceBuffer&) = delete;
    DeviceBuffer(DeviceBuffer&& other) noexcept : ptr(other.ptr) { other.ptr = nullptr; }
    DeviceBuffer& operator=(DeviceBuffer&& other) noexcept {
        if (this != &other) {
            if (ptr) cudaFree(ptr);
            ptr = other.ptr;
            other.ptr = nullptr;
        }
        return *this;
    }

    double* get() const { return ptr; }
};

// Bounds check helper for states / rate categories.
inline void validate_states_rate(int states, int rate_cats, int max_states, int max_rate_cats) {
    if (states <= 0 || rate_cats <= 0) throw std::runtime_error("Invalid states or rate categories.");
    if (states > max_states) throw std::runtime_error("states exceeds MAX_STATES.");
    if (rate_cats > max_rate_cats) throw std::runtime_error("rate_cats exceeds MAX_RATECATS.");
}
