/**
 **
 ** TileProcessor.cuh
 **
 ** Copyright (C) 2018 Elphel, Inc.
 **
 ** -----------------------------------------------------------------------------**
 **
 **  TileProcessor.cuh is free software: you can redistribute it and/or modify
 **  it under the terms of the GNU General Public License as published by
 **  the Free Software Foundation, either version 3 of the License, or
 **  (at your option) any later version.
 **
 **  This program is distributed in the hope that it will be useful,
 **  but WITHOUT ANY WARRANTY; without even the implied warranty of
 **  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 **  GNU General Public License for more details.
 **
 **  You should have received a copy of the GNU General Public License
 **  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 **
 **  Additional permission under GNU GPL version 3 section 7
 **
 **  If you modify this Program, or any covered work, by linking or
 **  combining it with NVIDIA Corporation's CUDA libraries from the
 **  NVIDIA CUDA Toolkit (or a modified version of those libraries),
 **  containing parts covered by the terms of NVIDIA CUDA Toolkit
 **  EULA, the licensors of this Program grant you additional
 **  permission to convey the resulting work.
 ** -----------------------------------------------------------------------------**
 */

/**
**************************************************************************
* \file TileProcessor.cuh
* \brief Top level of the Tile Processor for frequency domain

*/
// Avoiding includes in jcuda, all source files will be merged
#ifndef JCUDA
#pragma once
#include "dtt8x8.cuh"
#define THREADSX         (DTT_SIZE)
#define IMG_WIDTH       2592
#define IMG_HEIGHT      1936
#define KERNELS_HOR      164
#define KERNELS_VERT     123
#define NUM_CAMS           4
#define NUM_COLORS         3
#define KERNELS_LSTEP      4
#define THREADS_PER_TILE   8
#define TILES_PER_BLOCK    4
#define IMCLT_THREADS_PER_TILE 16
#define IMCLT_TILES_PER_BLOCK   4

#endif
//#define IMCLT14
//#define NOICLT 1
//#define TEST_IMCLT
//#define SAVE_CLT
// Not enough shared memory to have more threads per block,even just for the result clt tiles
// What to do:
// 1) make single image aberration correction: 1/4 of the result tiles
// With 4 cameras = calculate correlations (9x9), reusing kernel or just clt ones after color reducing, then output them to device memory
//Average run time =1308.124146 ms
//#define TILES_PER_BLOCK    2
//Average run time =12502.638672 - with 2 tiles/block it is longer!
///12129.268555 ms
//Average run time =4704.506348 ms (syncwarp)
//Average run time =4705.612305 ms (syncthreads)
//Average run time =1051.411255 ms
//Average run time =861.866577 ms
//Average run time =850.871277 ms had bugs
//Average run time =857.947632 ms fixed bugs


// Something broke, even w/o LPF: Average run time =1093.115112 ms
// without clt copying to device memory - Average run time =965.342407 ms - still worse
//Average run time =965.880554 ms
// combined tx and ty into a single int : Average run time =871.017944 ms
//Average run time =873.386597 ms (reduced number of registers)
//__umul24 : Average run time =879.125122 ms
// without __umul24 - back to Average run time =871.315552 ms
// Added copying clt to device memory - Average run time =942.071960 ms
// Removed rest of NOICLT : Average run time =943.456177 ms
// Added lpf: Average run time =1046.101318 ms (0.1 sec, 10%) - can be combined with the PSF kernel
//#define USE_UMUL24
////#define TILES_PER_BLOCK    4
//Average run time =5155.922852 ms
//Average run time =1166.388306 ms
//Average run time =988.750977 ms
//#define TILES_PER_BLOCK    8
//Average run time =9656.743164 ms
// Average run time =9422.057617 ms (reducing divergence)
//#define TILES_PER_BLOCK    1

//#define THREADS_PER_TILE   8
//#define IMCLT_THREADS_PER_TILE 16
//#define IMCLT_TILES_PER_BLOCK   4


#define KERNELS_STEP  (1 << KERNELS_LSTEP)
#define TILESX        (IMG_WIDTH / DTT_SIZE)
#define TILESY        (IMG_HEIGHT / DTT_SIZE)
// increase row length by 1 so vertical passes will use different ports
#define DTT_SIZE1        (DTT_SIZE + 1)
#define DTT_SIZE2        (2 * DTT_SIZE)
#define DTT_SIZE21       (DTT_SIZE2 + 1)

#define BAYER_RED   0
#define BAYER_BLUE  1
#define BAYER_GREEN 2
// assuming GR/BG as now
#define BAYER_RED_ROW 0
#define BAYER_RED_COL 1
//#define BAYER_BLUE_ROW (1 - BAYER_RED_ROW)
//#define BAYER_BLUE_COL (1 - BAYER_RED_COL)


#define DBG_TILE_X     174
#define DBG_TILE_Y     118

//#define DBG_TILE     (DBG_TILE_Y * 324 + DBG_TILE_X)
//#define DEBUG1 1
//#define DEBUG2 1
//#define DEBUG3 1
//#define DEBUG4 1
//#define DEBUG5 1
//56494
// struct tp_task
//#define TASK_SIZE      12
struct tp_task {
	int   task;
	int   txy;
//	short ty;
//	short tx;
	float xy[NUM_CAMS][2];
};
struct CltExtra{
	float data_x;   // kernel data is relative to this displacement X (0.5 pixel increments)
	float data_y;   // kernel data is relative to this displacement Y (0.5 pixel increments)
	float center_x; // actual center X (use to find derivatives)
	float center_y; // actual center X (use to find derivatives)
	float dxc_dx;   // add this to data_x per each pixel X-shift relative to the kernel center location
	float dxc_dy;   // same per each Y-shift pixel
	float dyc_dx;
	float dyc_dy;
};
/*
 Python code to generate constant coefficients:
def setup_hwindow(n=8, l=4):
    hwindow = [math.sin(math.pi*((1.0+2*i)/(4*n))) for i in range(2*n)]
    print("__constant__ float HWINDOW[] = {", end="") #
    for i in range (n):
        print("%ff"%(hwindow[i]), end ="")
        if i == (n-1):
            print("};")
        elif ((i + 1) % l) == 0:
            print(",")
            print("                                ", end ="")
        else:
            print(", ",end="")

def setup_hwindow2(n=8, l=4):
    hwindow = [0.5*math.sin(math.pi*((1.0+2*i)/(4*n))) for i in range(2*n)]
    print("__constant__ float HWINDOW2[] = {", end="") #
    for i in range (n):
        print("%ff"%(hwindow[i]), end ="")
        if i == (n-1):
            print("};")
        elif ((i + 1) % l) == 0:
            print(",")
            print("                                 ", end ="")
        else:
            print(", ",end="")

def get_fold_rindices(n=8):
    n1 = n>>1;
    rind = [0] * (2 * n) # reverse indices
    rcs =  [0] * (2 * n) # reverse signs for cosine term
    rss =  [0] * (2 * n) # reverse signs for sine term
    for x in range (n1):
        ri0 = n + n1 - x - 1
        ri1 = n + n1 + x
        ri2 =          x
        ri3 = n      - x - 1
        rind[ri0] = x
        rind[ri1] = x
        rind[ri2] = x + n1
        rind[ri3] = x + n1
        rcs[ri0] = -1
        rss[ri0] =  1
        rcs[ri1] = -1
        rss[ri1] = -1
        rcs[ri2] =  1
        rss[ri2] =  1
        rcs[ri3] = -1
        rss[ri3] =  1
    rind0 = []
    rind1 = []
    # generate start indices for the first 2 bayer rows
    for a in rind:
        rind0.append(a+rind[0]*n)
        rind1.append(a+rind[1]*n)
    #column increments for odd/even bayer rows
    inc_even = []
    inc_odd = []
    for i in range (n-1):
        inc_even.append(rind[2*i+2]-rind[2*i])
        inc_odd.append (rind[2*i+3]-rind[2*i+1])
    inc_even.reverse()
    inc_odd.reverse()
    # combine increments into int data
    inc_e = 0
    inc_o = 0
    for d in inc_even:
        inc_e = ((inc_e) << 4) | (d & 0xf)
    for d in inc_odd:
        inc_o = ((inc_o) << 4) | (d & 0xf)
    print("__constant__ int fold_indx2[2][%d] = {{"%(2*n), end="") #
    for d in rind0[:-1]:
        print('0x%2x,'%(d), end="")
    print('0x%2x},'%(rind0[-1]))
    print("                                      {", end="") #
    for d in rind1[:-1]:
        print('0x%2x,'%(d), end="")
    print('0x%2x}};'%(rind1[-1]))
    print("__constant__ int fold_inc[]=          {0x%08x, 0x%08x};"%(inc_e, inc_o))

def set_imclt_sa(stride=9):
    sa8 =[0x24,0x2c,0x34,0x3c,0x3c,0x34,0x2c,0x24,0x1c,0x14,0x0c,0x04,0x04,0x0c,0x14,0x1c]
    sa8s = [d // 8 + (d % 8) * stride for d in sa8]
    print("__constant__ int imclt_indx9[16] = {", end="") #
    for d in sa8s[:-1]:
        print('0x%02x,'%(d), end="")
    print('0x%2x};'%(sa8s[-1]))
*/


__constant__ float HWINDOW[] = {0.098017f, 0.290285f, 0.471397f, 0.634393f,
                                0.773010f, 0.881921f, 0.956940f, 0.995185f};

__constant__ float HWINDOW2[] = {0.049009f, 0.145142f, 0.235698f, 0.317197f,
                                 0.386505f, 0.440961f, 0.478470f, 0.497592f};


// Offsets in 8x8 DCT_CC/DST_SC tile for the first 2 lines of the 16x16 bayer image
__constant__ int fold_indx2[2][16] = {{0x24,0x25,0x26,0x27,0x27,0x26,0x25,0x24,0x23,0x22,0x21,0x20,0x20,0x21,0x22,0x23},
                                      {0x2c,0x2d,0x2e,0x2f,0x2f,0x2e,0x2d,0x2c,0x2b,0x2a,0x29,0x28,0x28,0x29,0x2a,0x2b}};

// increments of the offsets in 8x8 tile when going down, jumping two lines (same Bayer). Each 4 bits have to be <<3,
// addd to the current index and result should be AND-ed with 0x3f. inc_e is for even rows (0,2, ...) while inc_o - for odd ones (1,3,)
__constant__ int fold_inc[]=          {0x02feee12, 0x021eeef2};

//__constant__ int imclt_indx[16] = {0x24,0x2c,0x34,0x3c,0x3c,0x34,0x2c,0x24,0x1c,0x22,0x21,0x20,0x20,0x21,0x22,0x23};
//__constant__ int imclt_indx9[16] = {0x28,0x31,0x3a,0x43,0x43,0x3a,0x31,0x28,0x1f,0x16,0x0d,0x04,0x04,0x0d,0x16,0x1f};
__constant__ int imclt_indx9[16] = {0x28,0x29,0x2a,0x2b,0x2b,0x2a,0x29,0x28,0x27,0x26,0x25,0x24,0x24,0x25,0x26,0x27};


// Hope that if 2 outer indices are known at compile time there will be no integer multiplications
__constant__ float idct_signs[4][4][4] ={
		{ // quadrant 0, each elements corresponds to 4x4 pixel output, covering altogether 16x16
				{ 1,-1,-1,-1},
				{-1, 1, 1, 1},
				{-1, 1, 1, 1},
				{-1, 1, 1, 1}
		},{ // quadrant 1, each elements corresponds to 4x4 pixel output, covering altogether 16x16
				{ 1, 1, 1,-1},
				{-1,-1,-1, 1},
				{-1,-1,-1, 1},
				{-1,-1,-1, 1}
		},{ // quadrant 2, each elements corresponds to 4x4 pixel output, covering altogether 16x16
				{ 1,-1,-1,-1},
				{ 1,-1,-1,-1},
				{ 1,-1,-1,-1},
				{-1, 1, 1, 1}
		},{ // quadrant 3, each elements corresponds to 4x4 pixel output, covering altogether 16x16
				{ 1, 1, 1,-1},
				{ 1, 1, 1,-1},
				{ 1, 1, 1,-1},
				{-1,-1,-1, 1}
		}};
// LPF for sigma 0.9 each color (modify through cudaMemcpyToSymbol() or similar in Driver API
//#ifndef NOICLT
__constant__ float lpf_data[3][64]={
		{
				1.00000000f, 0.87041007f, 0.65943687f, 0.43487258f, 0.24970076f, 0.12518080f, 0.05616371f, 0.02728573f,
				0.87041007f, 0.75761368f, 0.57398049f, 0.37851747f, 0.21734206f, 0.10895863f, 0.04888546f, 0.02374977f,
				0.65943687f, 0.57398049f, 0.43485698f, 0.28677101f, 0.16466189f, 0.08254883f, 0.03703642f, 0.01799322f,
				0.43487258f, 0.37851747f, 0.28677101f, 0.18911416f, 0.10858801f, 0.05443770f, 0.02442406f, 0.01186582f,
				0.24970076f, 0.21734206f, 0.16466189f, 0.10858801f, 0.06235047f, 0.03125774f, 0.01402412f, 0.00681327f,
				0.12518080f, 0.10895863f, 0.08254883f, 0.05443770f, 0.03125774f, 0.01567023f, 0.00703062f, 0.00341565f,
				0.05616371f, 0.04888546f, 0.03703642f, 0.02442406f, 0.01402412f, 0.00703062f, 0.00315436f, 0.00153247f,
				0.02728573f, 0.02374977f, 0.01799322f, 0.01186582f, 0.00681327f, 0.00341565f, 0.00153247f, 0.00074451f
		},{
				1.00000000f, 0.87041007f, 0.65943687f, 0.43487258f, 0.24970076f, 0.12518080f, 0.05616371f, 0.02728573f,
				0.87041007f, 0.75761368f, 0.57398049f, 0.37851747f, 0.21734206f, 0.10895863f, 0.04888546f, 0.02374977f,
				0.65943687f, 0.57398049f, 0.43485698f, 0.28677101f, 0.16466189f, 0.08254883f, 0.03703642f, 0.01799322f,
				0.43487258f, 0.37851747f, 0.28677101f, 0.18911416f, 0.10858801f, 0.05443770f, 0.02442406f, 0.01186582f,
				0.24970076f, 0.21734206f, 0.16466189f, 0.10858801f, 0.06235047f, 0.03125774f, 0.01402412f, 0.00681327f,
				0.12518080f, 0.10895863f, 0.08254883f, 0.05443770f, 0.03125774f, 0.01567023f, 0.00703062f, 0.00341565f,
				0.05616371f, 0.04888546f, 0.03703642f, 0.02442406f, 0.01402412f, 0.00703062f, 0.00315436f, 0.00153247f,
				0.02728573f, 0.02374977f, 0.01799322f, 0.01186582f, 0.00681327f, 0.00341565f, 0.00153247f, 0.00074451f
		},{
				1.00000000f, 0.87041007f, 0.65943687f, 0.43487258f, 0.24970076f, 0.12518080f, 0.05616371f, 0.02728573f,
				0.87041007f, 0.75761368f, 0.57398049f, 0.37851747f, 0.21734206f, 0.10895863f, 0.04888546f, 0.02374977f,
				0.65943687f, 0.57398049f, 0.43485698f, 0.28677101f, 0.16466189f, 0.08254883f, 0.03703642f, 0.01799322f,
				0.43487258f, 0.37851747f, 0.28677101f, 0.18911416f, 0.10858801f, 0.05443770f, 0.02442406f, 0.01186582f,
				0.24970076f, 0.21734206f, 0.16466189f, 0.10858801f, 0.06235047f, 0.03125774f, 0.01402412f, 0.00681327f,
				0.12518080f, 0.10895863f, 0.08254883f, 0.05443770f, 0.03125774f, 0.01567023f, 0.00703062f, 0.00341565f,
				0.05616371f, 0.04888546f, 0.03703642f, 0.02442406f, 0.01402412f, 0.00703062f, 0.00315436f, 0.00153247f,
				0.02728573f, 0.02374977f, 0.01799322f, 0.01186582f, 0.00681327f, 0.00341565f, 0.00153247f, 0.00074451f
		}};
//#endif
__device__ void convertCorrectTile(
		struct CltExtra     * gpu_kernel_offsets, // [tileY][tileX][color]
		float               * gpu_kernels,        // [tileY][tileX][color]
		float               * gpu_images,
		float               * gpu_clt,
		const int             color,
		const int             lpf_mask,
		const float           centerX,
		const float           centerY,
//		const short          tx,
//		const short          ty,
		const int             txy,
		const size_t          dstride, // in floats (pixels)
		float               * clt_tile, //        [4][DTT_SIZE][DTT_SIZE1], // +1 to alternate column ports
		float               * clt_kernels, //      [4][DTT_SIZE][DTT_SIZE1], // +1 to alternate column ports
		int   int_topleft     [2],
		float residual_shift  [2],
	    float window_hor_cos  [2*DTT_SIZE],
	    float window_hor_sin  [2*DTT_SIZE],
	    float window_vert_cos [2*DTT_SIZE]);

// Fractional pixel shift (phase rotation), horizontal. In-place.
__device__ void shiftTileHor(
		float * clt_tile, //        [4][DTT_SIZE][DTT_SIZE1], // +1 to alternate column ports
		float residual_shift                         );
// Fractional pixel shift (phase rotation), vertical. In-place.
__device__ void shiftTileVert(
		float *clt_tile, //        [4][DTT_SIZE][DTT_SIZE1], // +1 to alternate column ports
		float residual_shift                         );
__device__ void convolveTiles(
		float* clt_tile, //    [4][DTT_SIZE][DTT_SIZE1], // 4 quadrants of the clt data, rows extended to optimize shared ports
		float* kernel); //      [4][DTT_SIZE][DTT_SIZE1]) // 4 quadrants of the CLT kernel (DTT3 converted)
__device__ void imclt(
		float * clt_tile,   //        [4][DTT_SIZE][DTT_SIZE1], // +1 to alternate column ports [4][8][9]
		float * mclt_tile ); //           [2* DTT_SIZE][DTT_SIZE1+ DTT_SIZE], // +1 to alternate column ports[16][17]
__device__ void imclt(
		float * clt_tile,   //        [4][DTT_SIZE][DTT_SIZE1], // +1 to alternate column ports [4][8][9]
		float * mclt_tile ); //           [2* DTT_SIZE][DTT_SIZE1+ DTT_SIZE], // +1 to alternate column ports[16][17]
__device__ void imclt_plane(
		int               color,
		float           * gpu_clt,   // [TILESY][TILESX][NUM_COLORS][DTT_SIZE*DTT_SIZE]
		float           * gpu_rbg,            // WIDTH, HEIGHT
		const size_t      dstride);            // in floats (pixels)

extern "C"
__global__ void convert_correct_tiles(
//		struct CltExtra ** gpu_kernel_offsets, // [NUM_CAMS], // changed for jcuda to avoid struct paraeters
		float           ** gpu_kernel_offsets, // [NUM_CAMS],
		float           ** gpu_kernels,        // [NUM_CAMS],
		float           ** gpu_images,         // [NUM_CAMS],
		struct tp_task  * gpu_tasks,
		float           ** gpu_clt,            // [NUM_CAMS][TILESY][TILESX][NUM_COLORS][DTT_SIZE*DTT_SIZE]
		size_t            dstride,             // in floats (pixels)
		int               num_tiles,           // number of tiles in task
		int               lpf_mask)            // apply lpf to colors : bit 0 - red, bit 1 - blue, bit2 - green

{
//	struct CltExtra* gpu_kernel_offsets = (struct CltExtra*) vgpu_kernel_offsets;
	dim3 t = threadIdx;
	int tile_in_block = threadIdx.y;
	int task_num = blockIdx.x * TILES_PER_BLOCK + tile_in_block;
	if (task_num >= num_tiles) return; // nothing to do
	struct tp_task  * gpu_task = &gpu_tasks[task_num];
	if (!gpu_task->task)       return; // NOP tile
	__shared__ struct tp_task tt [TILES_PER_BLOCK];
	// Copy task data to shared memory
	tt[tile_in_block].task =          gpu_task -> task;
//	tt[tile_in_block].tx =            gpu_task -> tx;
//	tt[tile_in_block].ty =            gpu_task -> ty;
	tt[tile_in_block].txy =           gpu_task -> txy;
	int thread0 =  threadIdx.x & 1;
	int thread12 = threadIdx.x >>1;
	if (thread12 < NUM_CAMS) {
		tt[tile_in_block].xy[thread12][thread0] = gpu_task -> xy[thread12][thread0];
	}
	if (NUM_CAMS > 4){ // unlikely
#pragma unroll
		for (int nc0 = 4; nc0 < NUM_CAMS; nc0 += 4){
			int nc = nc0 + thread12;
			if (nc < NUM_CAMS) {
				tt[tile_in_block].xy[nc][thread0] = gpu_task -> xy[nc][thread0];
			}
		}
	}
#pragma unroll
	for (int i = 0; i < (NUM_CAMS / 4); i++){
		int nc = (threadIdx.x >> 1) + (i << 2);
		if (nc < NUM_CAMS) {
			tt[tile_in_block].xy[nc][0] = gpu_task -> xy[nc][0];
			tt[tile_in_block].xy[nc][1] = gpu_task -> xy[nc][1];
		}

	}
     __syncthreads();// __syncwarp();
    __shared__ float clt_tile        [TILES_PER_BLOCK][4][DTT_SIZE][DTT_SIZE1];
    __shared__ float clt_kernels     [TILES_PER_BLOCK][4][DTT_SIZE][DTT_SIZE1]; // +1 to alternate column ports
    __shared__ int   int_topleft     [TILES_PER_BLOCK][2];
    __shared__ float residual_shift  [TILES_PER_BLOCK][2];
    __shared__ float window_hor_cos  [TILES_PER_BLOCK][2*DTT_SIZE];
    __shared__ float window_hor_sin  [TILES_PER_BLOCK][2*DTT_SIZE];
    __shared__ float window_vert_cos [TILES_PER_BLOCK][2*DTT_SIZE];

    // process each camera,l each color in series (to reduce shared memory)
    for (int ncam = 0; ncam <  NUM_CAMS; ncam++){
    	for (int color = 0; color <  NUM_COLORS; color++){
    		convertCorrectTile(
    				(struct CltExtra*)(gpu_kernel_offsets[ncam]),        // struct CltExtra* gpu_kernel_offsets,
					gpu_kernels[ncam],               // float           * gpu_kernels,
					gpu_images[ncam],                // float           * gpu_images,
					gpu_clt[ncam],                   // float           * gpu_clt,
					color,                           // const int         color,
					lpf_mask,                        // const int         lpf_mask,
					tt[tile_in_block].xy[ncam][0],   // const float       centerX,
					tt[tile_in_block].xy[ncam][1],   // const float       centerY,
//					tt[tile_in_block].tx | (tt[tile_in_block].ty <<16), //  const int txy,
					tt[tile_in_block].txy,           //  const int txy,
					dstride,                         // size_t            dstride, // in floats (pixels)
					(float * )(clt_tile [tile_in_block]),        // float clt_tile [TILES_PER_BLOCK][NUM_CAMS][NUM_COLORS][4][DTT_SIZE][DTT_SIZE])
					(float * )(clt_kernels[tile_in_block]),      // float clt_tile    [NUM_COLORS][4][DTT_SIZE][DTT_SIZE],
					int_topleft[tile_in_block],      // int   int_topleft  [NUM_COLORS][2],
					residual_shift[tile_in_block],   // float frac_topleft [NUM_COLORS][2],
					window_hor_cos[tile_in_block],   // float window_hor_cos  [NUM_COLORS][2*DTT_SIZE],
					window_hor_sin[tile_in_block],   //float window_hor_sin  [NUM_COLORS][2*DTT_SIZE],
					window_vert_cos[tile_in_block]); //float window_vert_cos [NUM_COLORS][2*DTT_SIZE]);
    		 __syncthreads();// __syncwarp();
    	}
    }
}

// Fractional pixel shift (phase rotation), horizontal. In-place. uses 8 threads (.x)
__device__ void shiftTileHor(
		float * clt_tile, //       [4][DTT_SIZE][DTT_SIZE1], // +1 to alternate column ports
		float residual_shift )
{
	int joffs = threadIdx.x; // * DTT_SIZE1;
	float * clt_tile_j0 = clt_tile +    joffs;                // ==&clt_tile[0][j][0]
	float * clt_tile_j1 = clt_tile_j0 + (DTT_SIZE1*DTT_SIZE); // ==&clt_tile[1][j][0]
	float * clt_tile_j2 = clt_tile_j1 + (DTT_SIZE1*DTT_SIZE); // ==&clt_tile[2][j][0]
	float * clt_tile_j3 = clt_tile_j2 + (DTT_SIZE1*DTT_SIZE); // ==&clt_tile[3][j][0]
	float x = residual_shift * ((threadIdx.x << 1 ) +1) * (0.5f/ DTT_SIZE);
	float ch =  cospif(x);
	float sh =  sinpif(x);
#pragma unroll
	for (int i = 0; i < DTT_SIZE; i++) {
		float clt_tile_a = *clt_tile_j0;
		float clt_tile_b = *clt_tile_j1;
		*clt_tile_j0 =  clt_tile_a * ch - clt_tile_b * sh;
		*clt_tile_j1 =  clt_tile_a * sh + clt_tile_b * ch;

		clt_tile_a = *clt_tile_j2;
		clt_tile_b = *clt_tile_j3;
		*clt_tile_j2 =  clt_tile_a * ch - clt_tile_b * sh;
		*clt_tile_j3 =  clt_tile_a * sh + clt_tile_b * ch;

		clt_tile_j0 +=DTT_SIZE1;
		clt_tile_j1 +=DTT_SIZE1;
		clt_tile_j2 +=DTT_SIZE1;
		clt_tile_j3 +=DTT_SIZE1;
	}
}


__device__ void shiftTileVert(
		float * clt_tile, //       [4][DTT_SIZE][DTT_SIZE1], // +1 to alternate column ports
		float residual_shift)
{
	int joffs = threadIdx.x * DTT_SIZE1;
	float * clt_tile_j0 = clt_tile +    joffs;                // ==&clt_tile[0][j][0]
	float * clt_tile_j1 = clt_tile_j0 + (DTT_SIZE1*DTT_SIZE); // ==&clt_tile[1][j][0]
	float * clt_tile_j2 = clt_tile_j1 + (DTT_SIZE1*DTT_SIZE); // ==&clt_tile[2][j][0]
	float * clt_tile_j3 = clt_tile_j2 + (DTT_SIZE1*DTT_SIZE); // ==&clt_tile[3][j][0]
	float x = residual_shift * ((threadIdx.x << 1 ) +1) * (0.5f/ DTT_SIZE);
	float ch =  cospif(x);
	float sh =  sinpif(x);
#pragma unroll
	for (int i = 0; i < DTT_SIZE; i++) {
		float clt_tile_a = *clt_tile_j0;
		float clt_tile_b = *clt_tile_j2;
		*clt_tile_j0 =  clt_tile_a * ch - clt_tile_b * sh;
		*clt_tile_j2 =  clt_tile_a * sh + clt_tile_b * ch;

		clt_tile_a = *clt_tile_j1;
		clt_tile_b = *clt_tile_j3;
		*clt_tile_j1 =  clt_tile_a * ch - clt_tile_b * sh;
		*clt_tile_j3 =  clt_tile_a * sh + clt_tile_b * ch;

		clt_tile_j0 ++;
		clt_tile_j1 ++;
		clt_tile_j2 ++;
		clt_tile_j3 ++;
	}
}

__device__ void convolveTiles(
		float* clt_tile, //    [4][DTT_SIZE][DTT_SIZE1], // 4 quadrants of the clt data, rows extended to optimize shared ports
		float* kernel) //      [4][DTT_SIZE][DTT_SIZE1]) // 4 quadrants of the CLT kernel (DTT3 converted)
{
	int joffs = threadIdx.x * DTT_SIZE1;
	float * kernel_j; //  =   kernel +      joffs;                // ==&kernel[0][j][0]
	float * clt_tile_j0 = clt_tile +    joffs;                // ==&clt_tile[0][j][0]
	float * clt_tile_j1 = clt_tile_j0 + (DTT_SIZE1*DTT_SIZE); // ==&clt_tile[1][j][0]
	float * clt_tile_j2 = clt_tile_j1 + (DTT_SIZE1*DTT_SIZE); // ==&clt_tile[2][j][0]
	float * clt_tile_j3 = clt_tile_j2 + (DTT_SIZE1*DTT_SIZE); // ==&clt_tile[3][j][0]
//#pragma unroll
	for (int i = 0; i < DTT_SIZE; i++){
		// k=0
		kernel_j =   kernel + joffs + i;
		float krn = *(kernel_j);
		float r0 =  *(clt_tile_j0) * krn;
		float r1 =  *(clt_tile_j1) * krn;
		float r2 =  *(clt_tile_j2) * krn;
		float r3 =  *(clt_tile_j3) * krn;
		// k = 1
		kernel_j += (DTT_SIZE1*DTT_SIZE);
		krn = *(kernel_j);
		r0 -=  *(clt_tile_j1) * krn;
		r1 +=  *(clt_tile_j0) * krn;
		r2 -=  *(clt_tile_j3) * krn;
		r3 +=  *(clt_tile_j2) * krn;
		// k=2
		kernel_j += (DTT_SIZE1*DTT_SIZE);
		krn = *(kernel_j);
		r0 -=  *(clt_tile_j2) * krn;
		r1 -=  *(clt_tile_j3) * krn;
		r2 +=  *(clt_tile_j0) * krn;
		r3 +=  *(clt_tile_j1) * krn;
		// k=3
		kernel_j += (DTT_SIZE1*DTT_SIZE);
		krn = *(kernel_j);
		r0 +=  *(clt_tile_j3) * krn;
		r1 -=  *(clt_tile_j2) * krn;
		r2 -=  *(clt_tile_j1) * krn;
		r3 +=  *(clt_tile_j0) * krn;
		*(clt_tile_j0)= r0;
		*(clt_tile_j1)= r1;
		*(clt_tile_j2)= r2;
		*(clt_tile_j3)= r3;
		clt_tile_j0 ++;
		clt_tile_j1 ++;
		clt_tile_j2 ++;
		clt_tile_j3 ++;
	}
}

__device__ void debug_print_clt1(
		float * clt_tile, //         [4][DTT_SIZE][DTT_SIZE1], // +1 to alternate column ports)
		const int color,
		int mask)
{
		if (color >= 0) printf("----------- Color = %d -----------\n",color);
		for (int dbg_quadrant = 0; dbg_quadrant < 4; dbg_quadrant++){
			printf("----------- Quadrant (c(h)-c(v), s-c, c-s, s-s) = %d -----------\n",dbg_quadrant);
			if ((mask >> dbg_quadrant) & 1) {
				for (int dbg_row = 0; dbg_row < DTT_SIZE; dbg_row++){
					for (int dbg_col = 0; dbg_col < DTT_SIZE; dbg_col++){
						printf ("%10.5f ", clt_tile[(dbg_quadrant*DTT_SIZE + dbg_row)*DTT_SIZE1 + dbg_col]);
					}
					printf("\n");
				}
			}
			printf("\n");
		}
}

__device__ void debug_print_mclt(
		float * mclt_tile, //         [4][DTT_SIZE][DTT_SIZE1], // +1 to alternate column ports)
		const int color)
{

	if (color >= 0) printf("----------- Color = %d -----------\n",color);
	for (int dbg_row = 0; dbg_row < DTT_SIZE2; dbg_row++){
		for (int dbg_col = 0; dbg_col < DTT_SIZE2; dbg_col++){
			printf ("%10.5f ", mclt_tile[dbg_row *DTT_SIZE21 + dbg_col]);
		}
		printf("\n");
	}
	printf("\n");
}

__device__ void convertCorrectTile(
		struct CltExtra     * gpu_kernel_offsets, // [tileY][tileX][color]
		float               * gpu_kernels,        // [tileY][tileX][color]
		float               * gpu_images,
		float               * gpu_clt,
		const int             color,
		const int             lpf_mask,
		const float           centerX,
		const float           centerY,
		const int             txy,
		const size_t          dstride, // in floats (pixels)
		float               * clt_tile, //        [4][DTT_SIZE][DTT_SIZE1], // +1 to alternate column ports
		float               * clt_kernels, //      [4][DTT_SIZE][DTT_SIZE1], // +1 to alternate column ports
		int   int_topleft     [2],
		float residual_shift  [2],
	    float window_hor_cos  [2*DTT_SIZE],
	    float window_hor_sin  [2*DTT_SIZE],
	    float window_vert_cos [2*DTT_SIZE])

{

	int   ktileX, ktileY;
	int   kernel_index; // common for all coors
	float kdx, kdy;
	if (threadIdx.x == 0){
		ktileX = min(KERNELS_HOR-1,  max(0, ((int) lrintf(centerX * (1.0/KERNELS_STEP)+1))));
		ktileY = min(KERNELS_VERT-1, max(0, ((int) lrintf(centerY * (1.0/KERNELS_STEP)+1))));
		kdx =    centerX - (ktileX << KERNELS_LSTEP) + (1 << (KERNELS_LSTEP -1)); // difference in pixel
		kdy =    centerY - (ktileY << KERNELS_LSTEP) + (1 << (KERNELS_LSTEP -1)); // difference in pixel
#ifdef USE_UMUL24
    	kernel_index = __umul24((ktileX + __umul24(ktileY, KERNELS_HOR)), NUM_COLORS);
#else
    	kernel_index = (ktileX + ktileY * KERNELS_HOR) * NUM_COLORS;
#endif
	}
////     __syncthreads();// __syncwarp();
    // broadcast kernel_index
    kernel_index =  __shfl_sync(
    		0xffffffff,        // unsigned mask,
			kernel_index,      // T var,
			0,                 // int srcLane,
			THREADS_PER_TILE); // int width=warpSize);
////     __syncthreads();// __syncwarp(); // is it needed?
    kdx =  __shfl_sync(
    		0xffffffff,        // unsigned mask,
			kdx,               // T var,
			0,                 // int srcLane,
			THREADS_PER_TILE); // int width=warpSize);
    kdy =  __shfl_sync(
    		0xffffffff,        // unsigned mask,
			kdy,               // T var,
			0,                 // int srcLane,
			THREADS_PER_TILE); // int width=warpSize);

     __syncthreads();// __syncwarp(); // is it needed?
    float px, py;
    // copy kernel
	int kernel_full_index = kernel_index + color;
#ifdef USE_UMUL24
    float * kernel_src = gpu_kernels + __umul24(kernel_full_index, (DTT_SIZE * DTT_SIZE * 4));
#else
    float * kernel_src = gpu_kernels + kernel_full_index* (DTT_SIZE * DTT_SIZE * 4);
#endif

    float * kernelp =  clt_kernels;

    kernel_src += threadIdx.x; // lsb;
    kernelp +=    threadIdx.x; // lsb;
#pragma unroll
    for (int j = 0; j < DTT_SIZE * 4; j++){ // all 4 components, 8 rows
    	// shared memory kernels use DTT_SIZE1 (same as image data)
    	*kernelp = *kernel_src;
    	kernelp+=DTT_SIZE1;
    	kernel_src+=THREADSX;
    }

    // Calculate offsets and prepare windows (all colors):
	struct CltExtra * clt_extra = &gpu_kernel_offsets[kernel_full_index];

	px = centerX - DTT_SIZE - (clt_extra->data_x + clt_extra->dxc_dx * kdx + clt_extra->dxc_dy * kdy) ; // fractional left corner
	int itlx = (int) floorf(px +0.5f);
	int_topleft [0] = itlx;
	float shift_hor =  itlx - px;
	residual_shift[0] = shift_hor;
	float x = shift_hor *(1.0f/16);
	float ahc = cospif(x);
	float ahs = sinpif(x);
	int i1 = DTT_SIZE;
	int i = 0;
	// embed sign for cosine and sine branches into window coefficients
#pragma unroll
	for (; i < (DTT_SIZE/2); i++ ){
		int ri = (DTT_SIZE-1) - i;
		window_hor_cos[i] =   HWINDOW[i ]*ahc + HWINDOW[ri]*ahs;
		window_hor_cos[i1] =  HWINDOW[ i]*ahs - HWINDOW[ri]*ahc;
		if (color == BAYER_GREEN){
			window_hor_sin[i] =    HWINDOW[i ]*ahc + HWINDOW[ri]*ahs; // bayer_color== 2
			window_hor_sin[i1] =   HWINDOW[ri]*ahc - HWINDOW[ i]*ahs;
		}
		i1++;
	}
	// embed sign for cosine and sine branches into window coefficients
#pragma unroll
	for (; i < DTT_SIZE; i++ ){
		int ri = (DTT_SIZE-1) - i;
		window_hor_cos[i] =   -HWINDOW[i ]*ahc - HWINDOW[ri]*ahs;
		window_hor_cos[i1] =   HWINDOW[ i]*ahs - HWINDOW[ri]*ahc;
		if (color == BAYER_GREEN){
			window_hor_sin[i] =    HWINDOW[i ]*ahc + HWINDOW[ri]*ahs;
			window_hor_sin[i1] =   HWINDOW[ i]*ahs - HWINDOW[ri]*ahc;
		}
		i1++;
	}

	py = centerY - DTT_SIZE - (clt_extra->data_y + clt_extra->dyc_dx * kdx + clt_extra->dyc_dy * kdy) ; // fractional top corner
	int itly = (int) floorf(py +0.5f);
	int_topleft[1] = itly;

	float shift_vert =  itly - py;
	residual_shift[1] = shift_vert;
	x = shift_vert *(1.0f/16);
	float avc = cospif(x);
	float avs = sinpif(x);
	i1 = DTT_SIZE; //**** Was commented out
	// embed sign for cosine branch only into window coefficients (for R,B only CC is needed, for G - CC and SC
	i = 0;
#pragma unroll
	for (; i < DTT_SIZE/2; i++ ){
		int ri = (DTT_SIZE-1) - i;
		window_vert_cos[i] =    HWINDOW[i ]*avc + HWINDOW[ri]*avs;
		window_vert_cos[i1++] = HWINDOW[ i]*avs - HWINDOW[ri]*avc;
	}
#pragma unroll
	for (; i < DTT_SIZE; i++ ){
		int ri = (DTT_SIZE-1) - i;
		window_vert_cos[i] =  -(HWINDOW[i ]*avc + HWINDOW[ri]*avs);
		window_vert_cos[i1++] = HWINDOW[ i]*avs - HWINDOW[ri]*avc;
	}


//    } //  if (color < 3) else
     __syncthreads();// __syncwarp();
#ifdef DEBUG1
    if ((threadIdx.x) == 0){
		printf("COLOR=%d\n",color);
		printf("centerX=%f,    centerY=%f\n",centerX, centerY);
		printf("ktileX=%d,     ktileY=%d\n", ktileX,  ktileY);
		printf("kdx=%f,        kdy=%f\n",    kdx, kdy);
		printf("int_topleft[%d][0]=%d,    int_topleft[%d][1]=%d\n",i,int_topleft[0],i,int_topleft[1]);
		printf("residual_shift[%d][0]=%f, residual_shift[%d][1]=%f\n",i,residual_shift[0],i,residual_shift[1]);
    }
     __syncthreads();// __syncwarp();
#endif



    // prepare, fold and write data to DTT buffers
    int dstride2 = dstride << 1; // in floats (pixels)
    int color0 = color & 1;
    int color1 = (color >>1) & 1;

    for (int gpass = 0; gpass < (color1 + 1); gpass++) { // Only once for R, B, twice - for G
    	int col_tl = int_topleft[0]; //  + (threadIdx.x << 1);
    	int row_tl = int_topleft[1];
    	// for red, blue and green, pass 0
    	int local_col = ((col_tl & 1) ^ (BAYER_RED_COL ^ color0 ^ color1 ^ gpass)) + (threadIdx.x << 1); // green red row: invert column from red
//    	int local_row = ((row_tl & 1) ^ BAYER_RED_ROW ^ gpass);                            // use red row
//    	int local_row = ((row_tl & 1) ^ BAYER_RED_ROW ^ color0 ^ color1 ^ gpass);                            // use red row
    	int local_row = ((row_tl & 1) ^ BAYER_RED_ROW ^ color0          ^ gpass);                            // use red row
    	float hwind_cos = window_hor_cos[local_col];
    	float hwind_sin = window_hor_sin[local_col]; // **** only used for green

    	int dtt_offset =     fold_indx2[local_row][local_col];
    	int dtt_offset_inc = fold_inc[local_row];
#ifdef USE_UMUL24
    	float *dct_buf = clt_tile + __umul24(gpass << 1 , DTT_SIZE * DTT_SIZE1);
    	float *dst_buf = clt_tile + __umul24((gpass << 1) + 1 , DTT_SIZE * DTT_SIZE1);   // **** only used for green
#else
    	float *dct_buf = clt_tile + ((gpass << 1) * (DTT_SIZE * DTT_SIZE1));
    	float *dst_buf = clt_tile + (((gpass << 1) + 1) * (DTT_SIZE * DTT_SIZE1));   // **** only used for green
#endif

    	if ((col_tl >= 0) && ((col_tl < (IMG_WIDTH - DTT_SIZE * 2))) && (row_tl >= 0) && ((row_tl < (IMG_HEIGHT - DTT_SIZE * 2)))) {
#ifdef USE_UMUL24
    		float *image_p = gpu_images + __umul24(row_tl + local_row, dstride)+ col_tl + local_col;
#else
    		float *image_p = gpu_images + dstride * (row_tl + local_row)+ col_tl + local_col;
#endif
#pragma unroll
    		for (int i = 0; i < 8; i++) {
//    			float d = (*image_p) * window_vert_cos[local_row]; //warp illegal address (0,2,1)
    			float d = (*image_p);
    			d *= window_vert_cos[local_row]; //warp illegal address (0,2,1)
    			int dtt_offset1 = dtt_offset + (dtt_offset >> 3); // converting for 9-long rows (DTT_SIZE1)
    			dct_buf[dtt_offset1] = d * hwind_cos;
    			dst_buf[dtt_offset1] = d * hwind_sin; // **** only used for green
    			dtt_offset = ( dtt_offset + ((dtt_offset_inc & 0xf) << 3)) & 0x3f;
    			dtt_offset_inc >>= 4;
    			local_row += 2;
    			image_p += dstride2;
    		}
    	} else { // handling border tiles (slower)
    		int eff_col = (min(IMG_HEIGHT/2 -1, max(0, col_tl >> 1)) << 1) + (col_tl & 1);
    		int row_lsb =  row_tl & 1;
    		int row_pair = row_tl >> 1;
#ifdef USE_UMUL24
    		float *image_p = gpu_images + __umul24(local_row, dstride) + (eff_col + local_col);
#else
    		float *image_p = gpu_images + dstride * local_row+ (eff_col + local_col);
#endif
#pragma unroll
    		for (int i = 0; i < 8; i++) {
    			int eff_row = (min(IMG_WIDTH/2 - 1, max(0, row_pair + i)) << 1) + row_lsb;
#ifdef USE_UMUL24
    			float d =  image_p[__umul24(eff_row,dstride)] * window_vert_cos[local_row];
#else
    			float d =  image_p[dstride * eff_row] * window_vert_cos[local_row];
#endif

    			int dtt_offset1 = dtt_offset + (dtt_offset >> 3); // converting for 9-long rows (DTT_SIZE1)
    			dct_buf[dtt_offset1] = d * hwind_cos;
    			dst_buf[dtt_offset1] = d * hwind_sin; // **** only used for green

    			dtt_offset = ( dtt_offset + ((dtt_offset_inc & 0xf) << 3)) & 0x3f;
    			dtt_offset_inc >>= 4;
    			local_row += 2;
    		}
    	}
    }
     __syncthreads();// __syncwarp();
#ifdef DEBUG2
    if ((threadIdx.x == 0) && (color == BAYER_GREEN)){
        printf("\nFOLDED DTT Tiles Green before reduction\n");
    	debug_print_clt1(clt_tile, color, 0xf); // all quadrants for green only
    }
     __syncthreads();// __syncwarp();
#endif
/*
    if (color == BAYER_GREEN) {
    	// reduce 4 green DTT buffers into 2 (so free future rotated green that were borrowed)
//    	float *dtt_buf =  ((float *) clt_tile[0]) + threadIdx.x;
//    	float *dtt_buf1 = ((float *) clt_tile[2]) + threadIdx.x;
    	float *dtt_buf =  clt_tile + threadIdx.x;
    	float *dtt_buf1 = dtt_buf+ (2 * DTT_SIZE1 * DTT_SIZE); // ((float *) clt_tile[2]) + threadIdx.x;

    	(*dtt_buf) += (*dtt_buf1);
    	dtt_buf +=    (4 * DTT_SIZE1);
    	dtt_buf1 +=   (4 * DTT_SIZE1);
    	(*dtt_buf) += (*dtt_buf1);

    	dtt_buf =         clt_tile + (DTT_SIZE1 * DTT_SIZE) + threadIdx.x; // ((float *) clt_tile[1]) + threadIdx.x;
    	dtt_buf1 =        dtt_buf +  (2 * DTT_SIZE1 * DTT_SIZE);           // ((float *) clt_tile[3]) + threadIdx.x;
    	(*dtt_buf) += (*dtt_buf1);
    	dtt_buf += (4 * DTT_SIZE1);
    	dtt_buf1 += (4 * DTT_SIZE1);
    	(*dtt_buf) += (*dtt_buf1);
    	 __syncthreads();// __syncwarp();
    }
*/
     if (color == BAYER_GREEN) {
    	 // reduce 4 green DTT buffers into 2 (so free future rotated green that were borrowed)
    	 float *dtt_buf =  clt_tile + threadIdx.x;
    	 float *dtt_buf1 = dtt_buf+ (2 * DTT_SIZE1 * DTT_SIZE); // ((float *) clt_tile[2]) + threadIdx.x;
#pragma unroll
    	 for (int i = 0; i < 2*DTT_SIZE; i++) {
    		 (*dtt_buf) += (*dtt_buf1);
    		 dtt_buf +=    DTT_SIZE1;
    		 dtt_buf1 +=   DTT_SIZE1;
    	 }
    	 __syncthreads();// __syncwarp();
     }

#ifdef DEBUG2
    if ((threadIdx.x) == 0){
        printf("\nFOLDED DTT Tiles,color=%d\n", color);
    	debug_print_clt1(clt_tile, color, (color== BAYER_GREEN)?3:1); // only 1 quadrant for R,B and 2 - for G
    }
     __syncthreads();// __syncwarp();
#endif

    dctiv_nodiverg( // all colors
#ifdef USE_UMUL24
    		clt_tile + __umul24(threadIdx.x,DTT_SIZE1), // [0][threadIdx.x], // pointer to start of row
#else
			clt_tile + (DTT_SIZE1 * threadIdx.x), // [0][threadIdx.x], // pointer to start of row
#endif
			1); //int inc);
    if (color == BAYER_GREEN){
        dstiv_nodiverg( // all colors
#ifdef USE_UMUL24
        		clt_tile + __umul24(threadIdx.x + DTT_SIZE, DTT_SIZE1), // clt_tile[1][threadIdx.x], // pointer to start of row
#else
				clt_tile + DTT_SIZE1 * (threadIdx.x + DTT_SIZE), // clt_tile[1][threadIdx.x], // pointer to start of row
#endif
    			1); //int inc);

    }
  	 __syncthreads();// __syncwarp();

#ifdef DEBUG2
    if ((threadIdx.x) == 0){
        printf("\nDTT Tiles after horizontal pass, color=%d\n",color);
    	debug_print_clt1(clt_tile, color, (color== BAYER_GREEN)?3:1); // only 1 quadrant for R,B and 2 - for G
    }
     __syncthreads();// __syncwarp();
#endif
    dctiv_nodiverg( // all colors
    		clt_tile + threadIdx.x, //  &clt_tile[0][0][threadIdx.x], // pointer to start of column
			DTT_SIZE1); // int inc,
    if (color == BAYER_GREEN){
//        dstiv_nodiverg( // all colors
          dctiv_nodiverg( // all colors
        		clt_tile + threadIdx.x + (DTT_SIZE1 * DTT_SIZE), // &clt_tile[1][0][threadIdx.x], // pointer to start of column
    			DTT_SIZE1); // int inc,
    }
  	 __syncthreads();// __syncwarp();

#ifdef DEBUG2
    if ((threadIdx.x) == 0){
        printf("\nDTT Tiles after vertical pass (both passes), color = %d\n",color);
    	debug_print_clt1(clt_tile, color, (color== BAYER_GREEN)?3:1); // only 1 quadrant for R,B and 2 - for G
    }
     __syncthreads();// __syncwarp();
#endif

// Replicate DTT, so non-bayer can still use same in-place rotation code
    float *src, *dst;
    int   negate; // , dst_inc;

    // Replicate horizontally (for R and B only):
    if (color != BAYER_GREEN) {
    	negate = 1-(((int_topleft[0] & 1) ^ (BAYER_RED_COL ^ color)) << 1); // +1/-1
    	src = clt_tile + threadIdx.x; // &clt_tile[0][0][threadIdx.x    ];
    	dst = clt_tile + (DTT_SIZE1 * DTT_SIZE) + (threadIdx.x ^ 7); // &clt_tile[1][0][threadIdx.x ^ 7];
#pragma unroll
    	for (int i = 0; i < DTT_SIZE; i++){
    		*dst = negate*(*src);
    		src += DTT_SIZE1;
    		dst += DTT_SIZE1;
    	}
    	 __syncthreads();// __syncwarp();
#ifdef DEBUG2
    	if ((threadIdx.x) == 0){
    		printf("\nDTT Tiles after first replicating, color = %d\n",color);
    		debug_print_clt1(clt_tile, color, 0x3);
    	}
    	 __syncthreads();// __syncwarp();
#endif

    }
    // replicate all colors down diagonal
	negate = 1-(((int_topleft[0] & 1) ^ (int_topleft[1] & 1) ^ (BAYER_RED_COL ^ BAYER_RED_ROW ^ (color >> 1))) << 1); // +1/-1 // 1 -
	// CC -> SS
	src = clt_tile + threadIdx.x; // &clt_tile[0][0][threadIdx.x    ];
	dst = clt_tile + (DTT_SIZE1 * (DTT_SIZE * 3 + 7)) +  (threadIdx.x ^ 7); // &clt_tile[3][7][threadIdx.x ^ 7];
#pragma unroll
    for (int i = 0; i < DTT_SIZE; i++){
    	*dst = negate*(*src);
    	src += DTT_SIZE1;
    	dst -= DTT_SIZE1;
    }
    //SC -> CS
	src = clt_tile + (DTT_SIZE1 * DTT_SIZE) + threadIdx.x; // &clt_tile[1][0][threadIdx.x    ];
	dst = clt_tile + (DTT_SIZE1 * (DTT_SIZE * 2 + 7)) +  (threadIdx.x ^ 7); // &clt_tile[2][7][threadIdx.x    ];
#pragma unroll
    for (int i = 0; i < DTT_SIZE; i++){
    	*dst = negate*(*src);
    	src += DTT_SIZE1;
    	dst -= DTT_SIZE1;
    }
#ifdef DEBUG2
    	if ((threadIdx.x) == 0){
    		printf("\nDTT Tiles after all replicating, color = %d\n",color);
    		debug_print_clt1(clt_tile, color, 0xf);
    	}
    	 __syncthreads();// __syncwarp();
#endif



#ifdef DEBUG2
    if ((threadIdx.x) == 0){
        printf("\nKernel tiles to convolve, color = %d\n",color);
    	debug_print_clt1(clt_kernels, color, 0xf); // all colors, all quadrants
    }
     __syncthreads();// __syncwarp();
#endif


    // convolve first, then rotate to match Java and make it easier to verify
    convolveTiles(
    		clt_tile,      // float clt_tile   [4][DTT_SIZE][DTT_SIZE1], // 4 quadrants of the clt data, rows extended to optimize shared ports
			clt_kernels); // float kernel     [4][DTT_SIZE][DTT_SIZE1]); // 4 quadrants of the CLT kernel (DTT3 converted)
     __syncthreads();// __syncwarp();
#ifdef DEBUG2
    if ((threadIdx.x) == 0){
        printf("\nDTT Tiles after convolution, color = %d\n",color);
    	debug_print_clt1(clt_tile, color, 0xf); // all colors, all quadrants
    }
     __syncthreads();// __syncwarp();
#endif

    // rotate phases: first horizontal, then vertical
    shiftTileHor(
    		clt_tile,           // float clt_tile       [4][DTT_SIZE][DTT_SIZE1], // +1 to alternate column ports
			residual_shift[0]); // float residual_shift);
     __syncthreads();// __syncwarp();
#ifdef DEBUG2
    if ((threadIdx.x) == 0){
        printf("\nDTT Tiles after horizontal shift, color = %d\n",color);
    	debug_print_clt1(clt_tile, color,  0xf); // only 1 quadrant for R,B and 2 - for G
    }
     __syncthreads();// __syncwarp();
#endif


    shiftTileVert(
    		clt_tile,           // float clt_tile       [4][DTT_SIZE][DTT_SIZE1], // +1 to alternate column ports
			residual_shift[1]); // float residual_shift);
     __syncthreads();// __syncwarp();
#ifdef DEBUG1
    if ((threadIdx.x) == 0){
        printf("\nDTT Tiles after vertical shift, color = %d\n",color);
    	debug_print_clt1(clt_tile, color,  0xf); // only 1 quadrant for R,B and 2 - for G
        printf("\nDTT All done\n");
    }
     __syncthreads();// __syncwarp();
#endif


#ifdef DBG_TILE
#ifdef DEBUG3
    if ((threadIdx.x) == 0){
        printf("\nDTT Tiles after vertical shift, color = %d\n",color);
    	debug_print_clt1(clt_tile, color,  0xf); // only 1 quadrant for R,B and 2 - for G
        printf("\nDTT All done\n");
    }
     __syncthreads();// __syncwarp();
#endif
#endif



     // optionally apply LF
     if ((lpf_mask >> color) & 1){
    	 float * clt = clt_tile + threadIdx.x;
#pragma unroll
    	 for (int q = 0; q < 4; q++) {
    		 float *lpf = lpf_data[color] + threadIdx.x;
#pragma unroll
    		 for (int i = 0; i <8; i++){
    			 (*clt) *= (*lpf);
    			 clt   += DTT_SIZE1;
    			 lpf   += DTT_SIZE;
    		 }
    	 }
         __syncthreads();// __syncwarp();
#ifdef DBG_TILE
#ifdef DEBUG3
         if ((threadIdx.x) == 0){
        	 printf("\nDTT Tiles after LPF, color = %d\n",color);
        	 debug_print_clt1(clt_tile, color,  0xf); // only 1 quadrant for R,B and 2 - for G
        	 printf("\nDTT All done\n");
         }
     __syncthreads();// __syncwarp();
#endif
#endif
     }
     //	const int         tx = txy & 0xffff; // slow again
     //	const int         ty = txy >> 16;

     int offset_src = threadIdx.x;
//     int offset_dst = ((ty * TILESX + tx)*NUM_COLORS + color)* ( 4 * DTT_SIZE * DTT_SIZE) + threadIdx.x; // gpu_kernels + kernel_full_index* (DTT_SIZE * DTT_SIZE * 4);
//     int offset_dst = ((ty * TILESX + tx)*NUM_COLORS + color)* ( 4 * DTT_SIZE * DTT_SIZE) + threadIdx.x; // gpu_kernels + kernel_full_index* (DTT_SIZE * DTT_SIZE * 4);
#ifdef USE_UMUL24
     int offset_dst = __umul24(    __umul24( __umul24(txy >> 16, TILESX) + (txy & 0xfff)  , NUM_COLORS) + color ,  4 * DTT_SIZE * DTT_SIZE) + threadIdx.x;
#else
     int offset_dst = (((txy >> 16) * TILESX + (txy & 0xfff))*NUM_COLORS + color)* ( 4 * DTT_SIZE * DTT_SIZE) + threadIdx.x;
#endif

     float * clt_src = clt_tile + offset_src; // threadIdx.x;
     float * clt_dst = gpu_clt +  offset_dst; // ((ty * TILESX + tx)*NUM_COLORS + color)* ( 4 * DTT_SIZE * DTT_SIZE1) + threadIdx.x; // gpu_kernels + kernel_full_index* (DTT_SIZE * DTT_SIZE * 4);
//#ifndef NOICLT

#ifdef DBG_TILE
#ifdef DEBUG3
    if ((threadIdx.x) == 0){
        printf("clt_src = 0x%lx\n",clt_src);
        printf("clt_dst = 0x%lx\n",clt_dst);
    }
#endif
#endif



#pragma unroll
    for (int j = 0; j < DTT_SIZE * 4; j++){ // all 4 components, 8 rows
    	// shared memory tiles use DTT_SIZE1
    	*clt_dst =  *clt_src;
    	clt_src   += DTT_SIZE1;
    	clt_dst   += DTT_SIZE;
    }
    __syncthreads();// __syncwarp();
    // just for testing perform imclt, save result to clt_kernels

//#endif
}

#ifndef NOICLT1


extern "C"
__global__ void test_imclt(
		float           * gpu_clt,            // [TILESY][TILESX][NUM_COLORS][DTT_SIZE*DTT_SIZE]
		int             ncam) // just for debug print
// Initially - no output, will add later

{
//	dim3 t = threadIdx;
	int tile_in_block = threadIdx.y;
	int tile_num = blockIdx.x * IMCLT_TILES_PER_BLOCK + tile_in_block;
	if (tile_num >= 1) return; // just testing with a single tile
	int thr3 =    threadIdx.x >> 3;
	int column =  threadIdx.x; // modify to use 2*8 threads, if needed.
//	int thr012 =  threadIdx.x & 7;

	// Read clt tile to
    __shared__ float clt_tiles        [IMCLT_TILES_PER_BLOCK][4][DTT_SIZE][DTT_SIZE1];
    __shared__ float mclt_tiles       [IMCLT_TILES_PER_BLOCK][DTT_SIZE2][DTT_SIZE21];

    // Read clt tile from device memory
    for (int color = 0; color < NUM_COLORS; color++) {
    	float * clt_tile = ((float *) clt_tiles) +  tile_in_block * (4 * DTT_SIZE * DTT_SIZE1); // top left quadrant0
    	float * gpu_tile = ((float *) gpu_clt) +  ((DBG_TILE_Y * TILESX + DBG_TILE_X) * NUM_COLORS + color) * (4 * DTT_SIZE * DTT_SIZE); // top left quadrant0

#ifdef DEBUG3
    	if ((threadIdx.x) == 0){
    		printf("\n\n\n================== gpu_tile = 0x%lx, clt_tile = 0x%lx, COLOR=%d, ncam = %d  ======================\n",gpu_tile,clt_tile,color,ncam);
    	}
#endif
    	clt_tile += column + thr3; // first 2 rows
    	gpu_tile += column;  // first 2 rows
#pragma unroll
    	for (int i = 0; i < DTT_SIZE2; i++){
    		*clt_tile= *gpu_tile;
    		clt_tile += (2 * DTT_SIZE1);
    		gpu_tile += (2 * DTT_SIZE);
    	}
    	// reset mclt tile to zero


    	float * mclt_tile = ((float*) mclt_tiles) +  tile_in_block * (DTT_SIZE2 * DTT_SIZE21) + column;
#pragma unroll
    	for (int i = 0; i < DTT_SIZE2; i++){
    		*mclt_tile= 0.0f;
    		mclt_tile += DTT_SIZE21;
    	}


    	__syncthreads();// __syncwarp();
    	imclt(
    			((float*) clt_tiles) +  tile_in_block * (4 * DTT_SIZE * DTT_SIZE1), // float * clt_tile,   //        [4][DTT_SIZE][DTT_SIZE1], // +1 to alternate column ports [4][8][9]
				((float*) mclt_tiles) +  tile_in_block * (DTT_SIZE2 * DTT_SIZE21)); // float * mclt_tile )
    	__syncthreads();// __syncwarp();
    }
}

extern "C"
__global__ void imclt_rbg(
		float           * gpu_clt,            // [TILESY][TILESX][NUM_COLORS][DTT_SIZE*DTT_SIZE]
		float           * gpu_rbg,            // WIDTH, 3 * HEIGHT
		int               color,
		int               v_offset,
		int               h_offset,
		const size_t      dstride)            // in floats (pixels)
{
	float *color_plane = gpu_rbg + dstride * (IMG_HEIGHT + DTT_SIZE) * color;
	int pass =           (v_offset << 1) + h_offset;     		// 0..3 to correctly acummulate 16x16 tiles stride 8
	int tile_in_block = threadIdx.y;
	int tile_num = blockIdx.x * IMCLT_TILES_PER_BLOCK + tile_in_block;
//	if (tile_num >= (TILESY * TILESX)) {
//		return; // just testing with a single tile
//	}
//	int tilesy_half = (TILESY + (v_offset ^ 1)) >> 1;
	int tilesx_half = (TILESX + (h_offset ^ 1)) >> 1;
	int tileY_half =  tile_num / tilesx_half;
	int tileX_half =  tile_num - tileY_half * tilesx_half;
	int tileY = (tileY_half << 1) + v_offset;
	int tileX = (tileX_half << 1) + h_offset;
	if (tileY >= TILESY) {
		return; // just testing with a single tile
	}
#ifdef DEBUG4
	if (threadIdx.x == 0) {
		if (tileY == DBG_TILE_Y) {
			printf("tileX == %d, tileY = %d\n",tileX, tileY);
		}
		if (tileX == DBG_TILE_X) {
			printf("tileX == %d, tileY = %d\n",tileX, tileY);
		}
		if ((tileX == DBG_TILE_X) && (tileY == DBG_TILE_Y)) {
			printf("tileX == %d, tileY = %d\n",tileX, tileY);
		}
	}
#endif

	int thr3 =    threadIdx.x >> 3;
	int column =  threadIdx.x; // modify to use 2 * 8 threads, if needed.

    __shared__ float clt_tiles        [IMCLT_TILES_PER_BLOCK][4][DTT_SIZE][DTT_SIZE1];
    __shared__ float mclt_tiles       [IMCLT_TILES_PER_BLOCK][DTT_SIZE2][DTT_SIZE21];

    // copy clt (frequency domain data)
    float * clt_tile = ((float *) clt_tiles) +  tile_in_block * (4 * DTT_SIZE * DTT_SIZE1); // top left quadrant0
//    float * gpu_tile = ((float *) gpu_clt) +  ((DBG_TILE_Y * TILESX + DBG_TILE_X) * NUM_COLORS + color) * (4 * DTT_SIZE * DTT_SIZE); // top left quadrant0
    float * gpu_tile = ((float *) gpu_clt) +  ((tileY * TILESX + tileX) * NUM_COLORS + color) * (4 * DTT_SIZE * DTT_SIZE); // top left quadrant0

    clt_tile += column + thr3; // first 2 rows
    gpu_tile += column;  // first 2 rows
#pragma unroll
    for (int i = 0; i < DTT_SIZE2; i++){
    	*clt_tile= *gpu_tile;
    	clt_tile += (2 * DTT_SIZE1);
    	gpu_tile += (2 * DTT_SIZE);
    }

	float * mclt_top = ((float*) mclt_tiles) +  tile_in_block * (DTT_SIZE2 * DTT_SIZE21) + column;
	float * rbg_top = color_plane + (tileY * DTT_SIZE)* dstride + (tileX * DTT_SIZE) + column;
	float * mclt_tile = mclt_top;

	if (pass == 0){ // just set mclt tile to all 0
#pragma unroll
		for (int i = 0; i < DTT_SIZE2; i++){
			*mclt_tile= 0.0f;
			mclt_tile += DTT_SIZE21;
		}
	} else {
		float * rbg_p = rbg_top;
#pragma unroll
			for (int i = 0; i < DTT_SIZE2; i++){
				*mclt_tile= *rbg_p;
				mclt_tile += DTT_SIZE21;
				rbg_p +=     dstride; // DTT_SIZE2;
			}
	}
	__syncthreads();// __syncwarp();
	imclt(
			((float*) clt_tiles) +  tile_in_block * (4 * DTT_SIZE * DTT_SIZE1), // float * clt_tile,   //        [4][DTT_SIZE][DTT_SIZE1], // +1 to alternate column ports [4][8][9]
			((float*) mclt_tiles) +  tile_in_block * (DTT_SIZE2 * DTT_SIZE21)); // float * mclt_tile )
	__syncthreads();// __syncwarp();

#ifdef DEBUG5
    if (((threadIdx.x) == 0) &&(tileX == DBG_TILE_X)  && (tileY == DBG_TILE_Y)){
//        printf("\nMCLT Tiles after IMCLT\n");
		printf("tileX == %d, tileY = %d\n",tileX, tileY);
    	debug_print_mclt(mclt_tile, -1); // only 1 quadrant for R,B and 2 - for G
    }
    __syncthreads();// __syncwarp();
#endif



//	save result (back)
	float * rbg_p = rbg_top;
	mclt_tile =     mclt_top;
	if ((tileX == 0)  && (tileY == 0)){
#pragma unroll
		for (int i = 0; i < DTT_SIZE2; i++){
			*rbg_p = 100.0f; // just testing
			mclt_tile += DTT_SIZE21;
			rbg_p +=     dstride; // DTT_SIZE2; // FIXME
		}
	} else if ((tileX == DBG_TILE_X)  && (tileY == DBG_TILE_Y)){
#pragma unroll
		for (int i = 0; i < DTT_SIZE2; i++){
			*rbg_p = (*mclt_tile) * 2.0; // just testing
			mclt_tile += DTT_SIZE21;
			rbg_p +=     dstride; // DTT_SIZE2; // FIXME
		}
	} else {
#pragma unroll
		for (int i = 0; i < DTT_SIZE2; i++){
			*rbg_p = *mclt_tile;
			mclt_tile += DTT_SIZE21;
			rbg_p +=     dstride; // DTT_SIZE2; // FIXME
		}
	}
}
/*

//	int margins = (tileX == 0) | ((tileY == 0) << 1) | ((tileX == (TILESX - 1)) << 2)| ((tileY == (TILESY - 1)) << 3); // bits 0 - left, 1 - top, 2 - right, 3 - bottom
//	int thr012 =  threadIdx.x & 7;
	// shift up/left by 4 pixels if no margins are used
//	float * rbg_tl = color_plane + (tileY * DTT_SIZE - (DTT_SIZE/2))* dstride + (tileX * DTT_SIZE - (DTT_SIZE/2));

		} else { // marginal tile
			int i = 0;
			int bottom = DTT_SIZE2;
			if (margins & 4){
				bottom -= DTT_SIZE2 - DTT_SIZE /2;
			}
			if (margins & 2) {
#pragma unroll
				for (i=0; i < (DTT_SIZE /2); i++){
					*mclt_tile= 0.0f;
					mclt_tile += DTT_SIZE21;
					rbg_p +=     DTT_SIZE2;
				}
			}
			if (margins & 1){
#pragma unroll
				for (; i < bottom; i++){
					if (column < (DTT_SIZE /2)) *mclt_tile= 0.0f;
					else                        *mclt_tile= *rbg_p;
					mclt_tile += DTT_SIZE21;
					rbg_p +=     DTT_SIZE2;
				}

			} else if (margins & 4){
#pragma unroll
				for (; i < bottom; i++){
					if (column >= (DTT_SIZE + DTT_SIZE /2)) *mclt_tile= 0.0f;
					else                                    *mclt_tile= *rbg_p;
					mclt_tile += DTT_SIZE21;
					rbg_p +=     DTT_SIZE2;
				}
			} else {
#pragma unroll
				for (; i < bottom; i++){
					*mclt_tile= *rbg_p;
					mclt_tile += DTT_SIZE21;
					rbg_p +=     DTT_SIZE2;
				}
			}
			if (margins & 8) {
#pragma unroll
				for (int i = 0; i < (DTT_SIZE /2); i++){
					*mclt_tile= 0.0f;
					mclt_tile += DTT_SIZE21;
					rbg_p +=     DTT_SIZE2;
				}
			}
		}
__device__ void imclt_plane(
		int               color,
		float           * gpu_clt,   // [TILESY][TILESX][NUM_COLORS][DTT_SIZE*DTT_SIZE]
		float           * gpu_rbg,            // WIDTH, HEIGHT
		const size_t      dstride)            // in floats (pixels)
{
	for (int v_offset = 0; v_offset < 2; v_offset++){
		for (int h_offset = 0; h_offset < 2; v_offset++){

		}

	}

}
	for (int color = 0; color < NUM_COLORS; color++){
		float *color_plane = gpu_rbg + dstride * IMG_HEIGHT * color;
		imclt_plane(
				color,        // int               color,
		        gpu_clt,      // float           * gpu_clt,   // [TILESY][TILESX][NUM_COLORS][DTT_SIZE*DTT_SIZE]
				color_plane,  // float           * gpu_rbg,            // WIDTH, HEIGHT
				dstride);     // const size_t      dstride)
	}

*/

//
// Uses 16 threads, gets 4*8*8 clt tiles, performs idtt-iv (swapping 1 and 2 quadrants) and then unfolds with window,
// adding to the output 16x16 tile (to use Read-modify-write with 4 passes over the frame. Shuld be zeroed before the
// first pass
//__constant__ int imclt_indx9[16] = {0x28,0x31,0x3a,0x43,0x43,0x3a,0x31,0x28,0x1f,0x16,0x0d,0x04,0x04,0x0d,0x16,0x1f};
__device__ void imclt(
		float * clt_tile,   //        [4][DTT_SIZE][DTT_SIZE1], // +1 to alternate column ports [4][8][9]
		float * mclt_tile ) //           [2* DTT_SIZE][DTT_SIZE1+ DTT_SIZE], // +1 to alternate column ports[16][17]
{
	int thr3 =    threadIdx.x >> 3;
	int column =  threadIdx.x; // modify to use 2*8 threads, if needed.
	int thr012 =  threadIdx.x & 7;
	int column4 = threadIdx.x >> 2;
//	int wcolumn =column ^ (7 * thr3); //0..7,7,..0
//	int wcolumn = ((thr3 << 3) -1) ^ thr3; //0..7,7,..0
	int wcolumn = ((thr3 << 3) - thr3) ^ thr012; //0..7,7,..0
	float * clt_tile1 = clt_tile +  (DTT_SIZE1 * DTT_SIZE);
	float * clt_tile2 = clt_tile1 + (DTT_SIZE1 * DTT_SIZE);
	float * clt_tile3 = clt_tile2 + (DTT_SIZE1 * DTT_SIZE);
#ifdef DEBUG3
    if ((threadIdx.x) == 0){
        printf("\nDTT Tiles before IDTT\n");
    	debug_print_clt1(clt_tile, -1,  0xf); // only 1 quadrant for R,B and 2 - for G
    }
     __syncthreads();// __syncwarp();
#endif

	// perform horizontal dct-iv on quadrants 0 and 1
    dctiv_nodiverg(
//    		clt_tile +  DTT_SIZE1 * (thr012 +   DTT_SIZE * thr3), // pointer to start of row for quadrants 0 and 1
    		clt_tile +  DTT_SIZE1 * (thr012 + 2*DTT_SIZE * thr3), // pointer to start of row for quadrants 0 and 2
			1);
	// perform horizontal dst-iv on quadrants 2 and 3
    dstiv_nodiverg( // all colors
//    		clt_tile2 + DTT_SIZE1 * (thr012 +   DTT_SIZE * thr3), // pointer to start of row for quadrants 2 and 3
    		clt_tile1 + DTT_SIZE1 * (thr012 + 2*DTT_SIZE * thr3), // pointer to start of row for quadrants 1 and 3
			1);
    __syncthreads();// __syncwarp();
	// perform vertical   dct-iv on quadrants 0 and 2
    dctiv_nodiverg(
//    		clt_tile +  thr012 + (DTT_SIZE1 * 2*DTT_SIZE) * thr3, // pointer to start of row for quadrants 0 and 2
    		clt_tile +  thr012 + (DTT_SIZE1 *   DTT_SIZE) * thr3, // pointer to start of row for quadrants 0 and 1
			DTT_SIZE1);
	// perform vertical   dst-iv on quadrants 1 and 3
    dstiv_nodiverg(
//    		clt_tile1 + thr012 + (DTT_SIZE1 * 2*DTT_SIZE) * thr3, // pointer to start of row for quadrants 1 and 3
    		clt_tile2 + thr012 + (DTT_SIZE1 *   DTT_SIZE) * thr3, // pointer to start of row for quadrants 2 and 3
			DTT_SIZE1);
    __syncthreads();// __syncwarp();

#ifdef DEBUG3
    if ((threadIdx.x) == 0){
        printf("\nDTT Tiles after IDTT\n");
    	debug_print_clt1(clt_tile, -1,  0xf); // only 1 quadrant for R,B and 2 - for G
    }
     __syncthreads();// __syncwarp();
#endif


    float hw = HWINDOW2[wcolumn];
    int clt_offset = imclt_indx9[column]; // index in each of the 4 iclt quadrants, accounting for stride=9
    float * rslt = mclt_tile + column;
#pragma unroll
    for (int i = 0; i < 4; i++){
    	float val = *rslt;
    	float w = HWINDOW2[i] * hw;
    	float d0 = idct_signs[0][0][column4] * (*(clt_tile +  clt_offset));
    	float d1 = idct_signs[1][0][column4] * (*(clt_tile1 + clt_offset));
    	float d2 = idct_signs[2][0][column4] * (*(clt_tile2 + clt_offset));
    	float d3 = idct_signs[3][0][column4] * (*(clt_tile3 + clt_offset));
    	d0+=d1;
    	d2+=d3;
    	d0+= d2;
    	if (i < 3){
    		clt_offset +=  DTT_SIZE1;
    	}
//    	*rslt = __fmaf_rd(w,d0,val); // w*d0 + val
    	val = __fmaf_rd(w,d0,val); // w*d0 + val
    	*rslt = val;
    	rslt += DTT_SIZE21;
    }
#pragma unroll
    for (int i = 4; i < 8; i++){
    	float val = *rslt;
    	float w = HWINDOW2[i] * hw;
    	float d0 = idct_signs[0][1][column4] * (*(clt_tile +  clt_offset));
    	float d1 = idct_signs[1][1][column4] * (*(clt_tile1 + clt_offset));
    	float d2 = idct_signs[2][1][column4] * (*(clt_tile2 + clt_offset));
    	float d3 = idct_signs[3][1][column4] * (*(clt_tile3 + clt_offset));
    	d0+=d1;
    	d2+=d3;
    	d0+= d2;
//    	if (i < 7){
   		clt_offset -=  DTT_SIZE1;
//    	}
    	*rslt = __fmaf_rd(w,d0,val); // w*d0 + val
    	rslt += DTT_SIZE21;
    }
#pragma unroll
    for (int i = 7; i >= 4; i--){
    	float val = *rslt;
    	float w = HWINDOW2[i] * hw;
    	float d0 = idct_signs[0][2][column4] * (*(clt_tile +  clt_offset));
    	float d1 = idct_signs[1][2][column4] * (*(clt_tile1 + clt_offset));
    	float d2 = idct_signs[2][2][column4] * (*(clt_tile2 + clt_offset));
    	float d3 = idct_signs[3][2][column4] * (*(clt_tile3 + clt_offset));
    	d0+=d1;
    	d2+=d3;
    	d0+= d2;
    	if (i > 4){
    		clt_offset -=  DTT_SIZE1;
    	}
    	*rslt = __fmaf_rd(w,d0,val); // w*d0 + val
    	rslt += DTT_SIZE21;
    }
#pragma unroll
    for (int i = 3; i >= 0; i--){
    	float val = *rslt;
    	float w = HWINDOW2[i] * hw;
    	float d0 = idct_signs[0][3][column4] * (*(clt_tile +  clt_offset));
    	float d1 = idct_signs[1][3][column4] * (*(clt_tile1 + clt_offset));
    	float d2 = idct_signs[2][3][column4] * (*(clt_tile2 + clt_offset));
    	float d3 = idct_signs[3][3][column4] * (*(clt_tile3 + clt_offset));
    	d0+=d1;
    	d2+=d3;
    	d0+= d2;
    	if (i > 0){
    		clt_offset +=  DTT_SIZE1;
    	}
    	*rslt = __fmaf_rd(w,d0,val); // w*d0 + val
    	rslt += DTT_SIZE21;
    }
#ifdef DEBUG3
    __syncthreads();// __syncwarp();
    if ((threadIdx.x) == 0){
        printf("\nMCLT Tiles after IMCLT\n");
    	debug_print_mclt(mclt_tile, -1); // only 1 quadrant for R,B and 2 - for G
    }
    __syncthreads();// __syncwarp();
#endif
}
#endif


