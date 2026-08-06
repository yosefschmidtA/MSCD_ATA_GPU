// Scratch for Phase 3

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <cuda_runtime.h>

typedef struct { float re,im; } Gcplx;

__constant__ int c_lamda[64];
static int *d_tevendim = NULL;
static int *d_tevenadd = NULL;
static float *d_tevenpar = NULL;
static float *d_talpha = NULL;
static float *d_tgamma = NULL;

static short2 *d_surviving_pairs = NULL;
static int h_pair_offset[16] = {0};
static int h_pair_count[16] = {0};

static Gcplx *d_asum = NULL;
static Gcplx *d_bsum = NULL;
static Gcplx *d_tevenelem = NULL;
static int g_tevenelem_size = 0;

extern "C" int mscdgpu_setup_summation(
    const int *tevencut, const int *tevendim, const int *tevenadd, 
    const float *tevenpar, const float *talpha, const float *tgamma,
    int ntrieven)
{
    // ... setup code ...
    return 0;
}

__global__ void k_summation(
    int m,
    int pair_count,
    const short2 *pairs,
    const Gcplx *in_sum,
    Gcplx *out_sum,
    const Gcplx *devendetec,
    const Gcplx *tevenelem,
    const int *tevendim,
    const int *tevenadd,
    const float *tevenpar,
    const float *talpha,
    const float *tgamma,
    const Gcplx *cexpix,
    int natoms,
    int radim,
    int exndata,
    int exmdata,
    int sizeint)
{
    // ... kernel code ...
}

extern "C" int mscdgpu_summation(float akin, const Gcplx *tevenelem, int ntrielem, Gcplx *asum_out, const float *patom)
{
    // ... host code ...
    return 0;
}
