// f64 STREAM-triad bandwidth probe for the S12 GPU spike.
// The whole design premise (PLAN ruling #9) is that our kernels are
// BANDWIDTH-bound, not FLOP-bound — the 4070 Ti SUPER runs f64 at ~1/64 its
// FP32 rate, so what governs throughput is how fast we can stream f64 cell
// state through memory, NOT how many f64 FLOPs the SMs can retire. This
// measures the ceiling that governs the THE-RUN envelope: sustained f64
// device-memory bandwidth (triad = 2 reads + 1 write per element).
#include <cstdio>
#include <cuda_runtime.h>

__global__ void triad(double* a, const double* b, const double* c, double s, size_t n) {
    size_t i = (size_t)blockIdx.x * blockDim.x + threadIdx.x;
    size_t stride = (size_t)gridDim.x * blockDim.x;
    for (; i < n; i += stride) a[i] = b[i] + s * c[i];
}

__global__ void init(double* a, double v, size_t n) {
    size_t i = (size_t)blockIdx.x * blockDim.x + threadIdx.x;
    size_t stride = (size_t)gridDim.x * blockDim.x;
    for (; i < n; i += stride) a[i] = v;
}

#define CK(x) do { cudaError_t e = (x); if (e != cudaSuccess) { \
    printf("CUDA error %s at %s:%d\n", cudaGetErrorString(e), __FILE__, __LINE__); return 1; } } while (0)

int main() {
    int dev = 0;
    cudaDeviceProp p;
    CK(cudaGetDeviceProperties(&p, dev));
    int memClkKHz = 0, busBits = 0;   // CUDA 13: read via the attribute API
    CK(cudaDeviceGetAttribute(&memClkKHz, cudaDevAttrMemoryClockRate, dev));
    CK(cudaDeviceGetAttribute(&busBits, cudaDevAttrGlobalMemoryBusWidth, dev));
    // Advertised peak = 2 (DDR) * memClk(kHz)*1e3 * busWidth(bits)/8 / 1e9 GB/s.
    double peak = 2.0 * (double)memClkKHz * 1e3 * (busBits / 8.0) / 1e9;
    printf("Device: %s | %.0f GB VRAM | bus %d-bit | advertised peak ~%.0f GB/s\n",
           p.name, p.totalGlobalMem / 1e9, busBits, peak);

    size_t n = (size_t)1 << 26;          // 67.1M doubles per array
    size_t bytes = n * sizeof(double);   // 536 MB per array, 1.6 GB total
    double *a, *b, *c;
    CK(cudaMalloc(&a, bytes)); CK(cudaMalloc(&b, bytes)); CK(cudaMalloc(&c, bytes));

    int block = 256, grid = 4096;        // grid-stride, occupancy-saturating
    init<<<grid, block>>>(b, 1.0, n);
    init<<<grid, block>>>(c, 2.0, n);
    CK(cudaDeviceSynchronize());

    triad<<<grid, block>>>(a, b, c, 3.0, n);   // warm up
    CK(cudaDeviceSynchronize());

    cudaEvent_t t0, t1; cudaEventCreate(&t0); cudaEventCreate(&t1);
    int iters = 100;
    CK(cudaEventRecord(t0));
    for (int it = 0; it < iters; ++it) triad<<<grid, block>>>(a, b, c, 3.0, n);
    CK(cudaEventRecord(t1)); CK(cudaEventSynchronize(t1));
    float ms; cudaEventElapsedTime(&ms, t0, t1);
    double sec = (ms / 1e3) / iters;
    double moved = 3.0 * (double)bytes;              // 2 read + 1 write
    double gbs = moved / 1e9 / sec;
    printf("triad f64: %.3f ms/iter, %.1f GB/s sustained (%.0f%% of peak)\n",
           sec * 1e3, gbs, 100.0 * gbs / peak);

    // Sanity: a[0] should be b + 3*c = 1 + 6 = 7.
    double h; CK(cudaMemcpy(&h, a, sizeof(double), cudaMemcpyDeviceToHost));
    printf("check a[0]=%.1f (expect 7.0)\n", h);
    cudaFree(a); cudaFree(b); cudaFree(c);
    return 0;
}
