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

#pragma once
#include "dtt8x8.cuh"
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
#define TILES_PER_BLOCK    4
//Average run time =5155.922852 ms
//Average run time =1166.388306 ms
//Average run time =988.750977 ms
//#define TILES_PER_BLOCK    8
//Average run time =9656.743164 ms
// Average run time =9422.057617 ms (reducing divergence)
//#define TILES_PER_BLOCK    1
#define THREADS_PER_TILE   8
#define IMG_WIDTH       2592
#define IMG_HEIGHT      1936
#define NUM_CAMS           4
#define NUM_COLORS         3
#define KERNELS_LSTEP      4
#define KERNELS_HOR      164
#define KERNELS_VERT     123
#define IMAGE_TILE_SIDE  18
#define KERNELS_STEP  (1 << KERNELS_LSTEP)
#define TILESX        (IMG_WIDTH / DTT_SIZE)
#define TILESY        (IMG_HEIGHT / DTT_SIZE)
// increase row length by 1 so vertical passes will use different ports
#define THREADSX         (DTT_SIZE)
#define DTT_SIZE1        (DTT_SIZE + 1)

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
//#undef  DEBUG2
//56494
// struct tp_task
//#define TASK_SIZE      12
struct tp_task {
	long task;
	short ty;
	short tx;
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
        if i == (2*n-1):
            print("};")
        elif ((i + 1) % l) == 0:
            print(",")
            print("                                ", end ="")
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

*/


__constant__ float HWINDOW[] = {0.098017f, 0.290285f, 0.471397f, 0.634393f,
                                0.773010f, 0.881921f, 0.956940f, 0.995185f};

// Offsets in 8x8 DCT_CC/DST_SC tile for the first 2 lines of the 16x16 bayer image
__constant__ int fold_indx2[2][16] = {{0x24,0x25,0x26,0x27,0x27,0x26,0x25,0x24,0x23,0x22,0x21,0x20,0x20,0x21,0x22,0x23},
                                      {0x2c,0x2d,0x2e,0x2f,0x2f,0x2e,0x2d,0x2c,0x2b,0x2a,0x29,0x28,0x28,0x29,0x2a,0x2b}};

// increments of the offsets in 8x8 tile when going down, jumping two lines (same Bayer). Each 4 bits have to be <<3,
// addd to the current index and result should be AND-ed with 0x3f. inc_e is for even rows (0,2, ...) while inc_o - for odd ones (1,3,)
__constant__ int fold_inc[]=          {0x02feee12, 0x021eeef2};

// index table for convolutions
__constant__ int zi[4][4] =  {{ 0, -1, -2,  3},
							  { 1,  0, -3, -2},
							  { 2, -3,  0, -1},
							  { 3,  2,  1,  0}};

__constant__ int za[4][4] =  {{ 0,  1,  2,  3},
							  { 1,  0,  3,  2},
							  { 2,  3,  0,  1},
							  { 3,  2,  1,  0}};
__constant__ int zs[4][4] =  {{ 0, -1, -1,  1},
							  { 1,  0, -1, -1},
							  { 1, -1,  0, -1},
							  { 1,  1,  1,  0}};



__device__ void convertCorrectTile(
		struct CltExtra     * gpu_kernel_offsets, // [tileY][tileX][color]
		float               * gpu_kernels,        // [tileY][tileX][color]
		float               * gpu_images,
		float               * gpu_clt,
		const int             color,
		const float           centerX,
		const float           centerY,
		const size_t          dstride, // in floats (pixels)
		// clt_tile[0] - before rotation, [0][0] - R:DCT/DCT, [0][1] - B:DCT/DCT, [0][2] - G:DCT/DCT, [0][3] - G:DST/DCT,
		// clt_tile[1], clt_tile[2], and clt_tile[3] - after rotation, 4 quadrants each
	    // changed, above is wrong now
		//		float clt_tile        [4][DTT_SIZE][DTT_SIZE1], // +1 to alternate column ports
		float               * clt_tile, //        [4][DTT_SIZE][DTT_SIZE1], // +1 to alternate column ports
		//		float clt_kernels     [4][DTT_SIZE][DTT_SIZE1], // +1 to alternate column ports
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
__device__ void shiftTileHor1(
		float clt_tile       [4][DTT_SIZE][DTT_SIZE1], // +1 to alternate column ports
		float residual_shift                         );

// Fractional pixel shift (phase rotation), vertical. In-place.
__device__ void shiftTileVert(
		float *clt_tile, //        [4][DTT_SIZE][DTT_SIZE1], // +1 to alternate column ports
		float residual_shift                         );
__device__ void shiftTileVert1(
		float clt_tile       [4][DTT_SIZE][DTT_SIZE1], // +1 to alternate column ports
		float residual_shift                         );

__device__ void convolveTiles(
		float* clt_tile, //    [4][DTT_SIZE][DTT_SIZE1], // 4 quadrants of the clt data, rows extended to optimize shared ports
		float* kernel); //      [4][DTT_SIZE][DTT_SIZE1]) // 4 quadrants of the CLT kernel (DTT3 converted)

__device__ void convolveTiles0(
		float clt_tile   [4][DTT_SIZE][DTT_SIZE1], // 4 quadrants of the clt data, rows extended to optimize shared ports
		float kernel     [4][DTT_SIZE][DTT_SIZE1]); // 4 quadrants of the CLT kernel (DTT3 converted)


extern "C"
__global__ void tileProcessor(
		struct CltExtra ** gpu_kernel_offsets, // [NUM_CAMS],
		float           ** gpu_kernels,        // [NUM_CAMS],
		float           ** gpu_images,         // [NUM_CAMS],
		struct tp_task  * gpu_tasks,
		float           ** gpu_clt,            // [NUM_CAMS][TILESY][TILESX][NUM_COLORS][DTT_SIZE*DTT_SIZE]
		size_t            dstride, // // in floats (pixels)
		int               num_tiles) // number of tiles in task

{
	dim3 t = threadIdx;
	int tile_in_block = threadIdx.y;
	int task_num = blockIdx.x * TILES_PER_BLOCK + tile_in_block;
	if (task_num >= num_tiles) return; // nothing to do
	struct tp_task  * gpu_task = &gpu_tasks[task_num];
	if (!gpu_task->task)       return; // NOP tile
	__shared__ struct tp_task tt [TILES_PER_BLOCK];
	// Copy task data to shared memory
	tt[tile_in_block].task =          gpu_task -> task;
	tt[tile_in_block].tx =            gpu_task -> tx;
	tt[tile_in_block].ty =            gpu_task -> ty;
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


    // set memory for CLT result (per tile, per camera, per color, per clt, per row, per column
	// clt_tile[][0] - before rotation, [][0][0] - R:DCT/DCT, [][0][1] - B:DCT/DCT, [][0][2] - G:DCT/DCT, [][0][3] - G:DST/DCT,
	// clt_tile[][1], clt_tile[][2], and clt_tile[][3] - after rotation, 4 quadrants each
    // changed, above is wrong now
///    __shared__ float clt_tile        [TILES_PER_BLOCK][NUM_CAMS][NUM_COLORS][4][DTT_SIZE][DTT_SIZE1];
    __shared__ float clt_tile        [TILES_PER_BLOCK][4][DTT_SIZE][DTT_SIZE1];
    // sharing shared memory for cameras as they are corrected one after another
    // TODO: evaluate total shared memory usage, maybe this sharing is not needed

    __shared__ float clt_kernels     [TILES_PER_BLOCK][4][DTT_SIZE][DTT_SIZE1]; // +1 to alternate column ports
    __shared__ int   int_topleft     [TILES_PER_BLOCK][2];
    __shared__ float residual_shift  [TILES_PER_BLOCK][2];

    __shared__ float window_hor_cos  [TILES_PER_BLOCK][2*DTT_SIZE];
    __shared__ float window_hor_sin  [TILES_PER_BLOCK][2*DTT_SIZE];
    __shared__ float window_vert_cos [TILES_PER_BLOCK][2*DTT_SIZE];

//IMAGE_TILE_SIDE
    // process each camera in series
    for (int ncam = 0; ncam <  NUM_CAMS; ncam++){
    	for (int color = 0; color <  NUM_COLORS; color++){
    		convertCorrectTile(
    				gpu_kernel_offsets[ncam],        // float           * gpu_kernel_offsets,
					gpu_kernels[ncam],               // float           * gpu_kernels,
					gpu_images[ncam],                // float           * gpu_images,
					gpu_clt[ncam],                   // float           * gpu_clt,
					color,                           // const int         color,
					tt[tile_in_block].xy[ncam][0],   // const float       centerX,
					tt[tile_in_block].xy[ncam][1],   // const float       centerY,
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



__device__ void shiftTileHor2(
		float *fclt_tile, //   [4][DTT_SIZE][DTT_SIZE1], // 4 quadrants of the clt data, rows extended to optimize shared ports
		float residual_shift )
{
	float (*clt_tile) [4][DTT_SIZE][DTT_SIZE1] =  (float(*)[4][DTT_SIZE][DTT_SIZE1]) fclt_tile;
	int j =   threadIdx.x;
	float x = residual_shift * ((j << 1 ) +1) * (0.5f/ DTT_SIZE);
	float ch =  cospif(x);
	float sh =  sinpif(x);
#pragma unroll
	for (int i = 0; i < DTT_SIZE; i++) {
		float t =           (*clt_tile)[0][i][j] * ch - (*clt_tile)[1][i][j] * sh;
		(*clt_tile)[1][i][j] = (*clt_tile)[0][i][j] * sh + (*clt_tile)[1][i][j] * ch;
		(*clt_tile)[0][i][j] = t;

		t =                 (*clt_tile)[2][i][j] * ch - (*clt_tile)[3][i][j] * sh;
		(*clt_tile)[3][i][j] = (*clt_tile)[2][i][j] * sh + (*clt_tile)[3][i][j] * ch;
		(*clt_tile)[2][i][j] = t;
	}
}

__device__ void shiftTileHor1(
		float clt_tile       [4][DTT_SIZE][DTT_SIZE1], // +1 to alternate column ports
		float residual_shift )
{
	int j =   threadIdx.x;
	float x = residual_shift * ((j << 1 ) +1) * (0.5f/ DTT_SIZE);
	float ch =  cospif(x);
	float sh =  sinpif(x);
#pragma unroll
	for (int i = 0; i < DTT_SIZE; i++) {
		float t =           clt_tile[0][i][j] * ch - clt_tile[1][i][j] * sh;
		clt_tile[1][i][j] = clt_tile[0][i][j] * sh + clt_tile[1][i][j] * ch;
		clt_tile[0][i][j] = t;

		t =                 clt_tile[2][i][j] * ch - clt_tile[3][i][j] * sh;
		clt_tile[3][i][j] = clt_tile[2][i][j] * sh + clt_tile[3][i][j] * ch;
		clt_tile[2][i][j] = t;
	}
}







// Fractional pixel shift (phase rotation), vertical. In-place.
__device__ void shiftTileVert0(
		float clt_tile       [4][DTT_SIZE][DTT_SIZE1], // +1 to alternate column ports
		float residual_shift)
{
	int j =   threadIdx.x;
	float x = residual_shift * ((j << 1 ) +1) * (0.5f/ DTT_SIZE);
	float ch =  cospif(x);
	float sh =  sinpif(x);
#pragma unroll
	for (int i = 0; i < DTT_SIZE; i++) {
		float t =           clt_tile[0][j][i] * ch - clt_tile[1][j][i] * sh;
		clt_tile[1][j][i] = clt_tile[0][j][i] * sh + clt_tile[1][j][i] * ch;
		clt_tile[0][j][i] = t;

		t =                 clt_tile[2][j][i] * ch - clt_tile[3][j][i] * sh;
		clt_tile[3][j][i] = clt_tile[2][j][i] * sh + clt_tile[3][j][i] * ch;
		clt_tile[2][j][i] = t;
	}
}

__device__ void shiftTileVert1(
		float clt_tile       [4][DTT_SIZE][DTT_SIZE1], // +1 to alternate column ports
		float residual_shift)
{
	int j =   threadIdx.x;
	float x = residual_shift * ((j << 1 ) +1) * (0.5f/ DTT_SIZE);
	float ch =  cospif(x);
	float sh =  sinpif(x);
#pragma unroll
	for (int i = 0; i < DTT_SIZE; i++) {
		float t =           clt_tile[0][j][i] * ch - clt_tile[2][j][i] * sh;
		clt_tile[2][j][i] = clt_tile[0][j][i] * sh + clt_tile[2][j][i] * ch;
		clt_tile[0][j][i] = t;

		t =                 clt_tile[1][j][i] * ch - clt_tile[3][j][i] * sh;
		clt_tile[3][j][i] = clt_tile[1][j][i] * sh + clt_tile[3][j][i] * ch;
		clt_tile[1][j][i] = t;
	}
}

__device__ void shiftTileVert2(
		float *fclt_tile, //   [4][DTT_SIZE][DTT_SIZE1], // 4 quadrants of the clt data, rows extended to optimize shared ports
		float residual_shift)
{
	float (*clt_tile) [4][DTT_SIZE][DTT_SIZE1] =  (float(*)[4][DTT_SIZE][DTT_SIZE1]) fclt_tile;
	int j =   threadIdx.x;
	float x = residual_shift * ((j << 1 ) +1) * (0.5f/ DTT_SIZE);
	float ch =  cospif(x);
	float sh =  sinpif(x);
#pragma unroll
	for (int i = 0; i < DTT_SIZE; i++) {
		float t =           (*clt_tile)[0][j][i] * ch - (*clt_tile)[2][j][i] * sh;
		(*clt_tile)[2][j][i] = (*clt_tile)[0][j][i] * sh + (*clt_tile)[2][j][i] * ch;
		(*clt_tile)[0][j][i] = t;

		t =                 (*clt_tile)[1][j][i] * ch - (*clt_tile)[3][j][i] * sh;
		(*clt_tile)[3][j][i] = (*clt_tile)[1][j][i] * sh + (*clt_tile)[3][j][i] * ch;
		(*clt_tile)[1][j][i] = t;
	}
}


// Fractional pixel shift (phase rotation), vertical. In-place.


__device__ void convolveTiles1(
		float clt_tile   [4][DTT_SIZE][DTT_SIZE1], // 4 quadrants of the clt data, rows extended to optimize shared ports
		float kernel     [4][DTT_SIZE][DTT_SIZE1]) // 4 quadrants of the CLT kernel (DTT3 converted)
{
	int j = threadIdx.x;
	for (int i = 0; i < DTT_SIZE; i++){
		float r0 = 0;
		float r1 = 0;
		float r2 = 0;
		float r3 = 0;
		for (int k = 0; k < 4; k++){
			r0 += zs[0][k]*clt_tile[za[0][k]][j][i] * kernel[k][j][i];
			r1 += zs[1][k]*clt_tile[za[1][k]][j][i] * kernel[k][j][i];
			r2 += zs[2][k]*clt_tile[za[2][k]][j][i] * kernel[k][j][i];
			r3 += zs[3][k]*clt_tile[za[3][k]][j][i] * kernel[k][j][i];
		}
		clt_tile[0][j][i]= r0;
		clt_tile[1][j][i]= r1;
		clt_tile[2][j][i]= r2;
		clt_tile[3][j][i]= r3;
	}
}

__device__ void convolveTiles0(
		float clt_tile   [4][DTT_SIZE][DTT_SIZE1], // 4 quadrants of the clt data, rows extended to optimize shared ports
		float kernel     [4][DTT_SIZE][DTT_SIZE1]) // 4 quadrants of the CLT kernel (DTT3 converted)
{
	int j = threadIdx.x;
	for (int i = 0; i < DTT_SIZE; i++){
		float r0 = 0;
		float r1 = 0;
		float r2 = 0;
		float r3 = 0;
		for (int k = 0; k < 4; k++){
			if (zi[0][k] < 0) r0 -= clt_tile[-zi[0][k]][j][i] * kernel[k][j][i];
			else			  r0 += clt_tile[ zi[0][k]][j][i] * kernel[k][j][i];

			if (zi[1][k] < 0) r1 -= clt_tile[-zi[1][k]][j][i] * kernel[k][j][i];
			else			  r1 += clt_tile[ zi[1][k]][j][i] * kernel[k][j][i];

			if (zi[2][k] < 0) r2 -= clt_tile[-zi[2][k]][j][i] * kernel[k][j][i];
			else			  r2 += clt_tile[ zi[2][k]][j][i] * kernel[k][j][i];

			if (zi[3][k] < 0) r3 -= clt_tile[-zi[3][k]][j][i] * kernel[k][j][i];
			else			  r3 += clt_tile[ zi[3][k]][j][i] * kernel[k][j][i];
		}
		clt_tile[0][j][i]= r0;
		clt_tile[1][j][i]= r1;
		clt_tile[2][j][i]= r2;
		clt_tile[3][j][i]= r3;
	}
}

__device__ void convolveTiles2(
		float *fclt_tile, //   [4][DTT_SIZE][DTT_SIZE1], // 4 quadrants of the clt data, rows extended to optimize shared ports
		float *fkernel) //      [4][DTT_SIZE][DTT_SIZE1]) // 4 quadrants of the CLT kernel (DTT3 converted)
{
	float (*clt_tile) [4][DTT_SIZE][DTT_SIZE1] =  (float(*)[4][DTT_SIZE][DTT_SIZE1]) fclt_tile;
	float (*kernel)    [4][DTT_SIZE][DTT_SIZE1] = (float(*)[4][DTT_SIZE][DTT_SIZE1]) fkernel;
	int j = threadIdx.x;
	for (int i = 0; i < DTT_SIZE; i++){
		float r0 = 0;
		float r1 = 0;
		float r2 = 0;
		float r3 = 0;
		for (int k = 0; k < 4; k++){
			if (zi[0][k] < 0) r0 -= (*clt_tile)[-zi[0][k]][j][i] * (*kernel)[k][j][i];
			else			  r0 += (*clt_tile)[ zi[0][k]][j][i] * (*kernel)[k][j][i];

			if (zi[1][k] < 0) r1 -= (*clt_tile)[-zi[1][k]][j][i] * (*kernel)[k][j][i];
			else			  r1 += (*clt_tile)[ zi[1][k]][j][i] * (*kernel)[k][j][i];

			if (zi[2][k] < 0) r2 -= (*clt_tile)[-zi[2][k]][j][i] * (*kernel)[k][j][i];
			else			  r2 += (*clt_tile)[ zi[2][k]][j][i] * (*kernel)[k][j][i];

			if (zi[3][k] < 0) r3 -= (*clt_tile)[-zi[3][k]][j][i] * (*kernel)[k][j][i];
			else			  r3 += (*clt_tile)[ zi[3][k]][j][i] * (*kernel)[k][j][i];
		}
		(*clt_tile)[0][j][i]= r0;
		(*clt_tile)[1][j][i]= r1;
		(*clt_tile)[2][j][i]= r2;
		(*clt_tile)[3][j][i]= r3;
	}
}




__device__ void debug_print_clt(
		float clt_tile        [4][DTT_SIZE][DTT_SIZE1], // +1 to alternate column ports)
		const int color,
		int mask)
{
		printf("----------- Color = %d -----------\n",color);
		for (int dbg_quadrant = 0; dbg_quadrant < 4; dbg_quadrant++){
			printf("----------- Quadrant (c(h)-c(v), s-c, c-s, s-s) = %d -----------\n",dbg_quadrant);
			if ((mask >> dbg_quadrant) & 1) {
				for (int dbg_row = 0; dbg_row < DTT_SIZE; dbg_row++){
					for (int dbg_col = 0; dbg_col < DTT_SIZE; dbg_col++){
						printf ("%10.5f ", clt_tile[dbg_quadrant][dbg_row][dbg_col]);
					}
					printf("\n");
				}
			}
			printf("\n");
		}
}
__device__ void debug_print_clt1(
		float * clt_tile, //         [4][DTT_SIZE][DTT_SIZE1], // +1 to alternate column ports)
		const int color,
		int mask)
{
		printf("----------- Color = %d -----------\n",color);
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



// Uses 32 threads
__device__ void convertCorrectTile(
		struct CltExtra     * gpu_kernel_offsets, // [tileY][tileX][color]
		float               * gpu_kernels,        // [tileY][tileX][color]
		float               * gpu_images,
		float               * gpu_clt,
		const int             color,
		const float           centerX,
		const float           centerY,
		const size_t          dstride, // in floats (pixels)
		// clt_tile[0] - before rotation, [0][0] - R:DCT/DCT, [0][1] - B:DCT/DCT, [0][2] - G:DCT/DCT, [0][3] - G:DST/DCT,
		// clt_tile[1], clt_tile[2], and clt_tile[3] - after rotation, 4 quadrants each
	    // changed, above is wrong now
//		float clt_tile        [4][DTT_SIZE][DTT_SIZE1], // +1 to alternate column ports
		float               * clt_tile, //        [4][DTT_SIZE][DTT_SIZE1], // +1 to alternate column ports
//		float clt_kernels     [4][DTT_SIZE][DTT_SIZE1], // +1 to alternate column ports
		float               * clt_kernels, //      [4][DTT_SIZE][DTT_SIZE1], // +1 to alternate column ports
		int   int_topleft     [2],
		float residual_shift  [2],
	    float window_hor_cos  [2*DTT_SIZE],
	    float window_hor_sin  [2*DTT_SIZE],
	    float window_vert_cos [2*DTT_SIZE])

{

///    __shared__ float window_hor_cos  [NUM_COLORS][2*DTT_SIZE];
///    __shared__ float window_hor_sin  [NUM_COLORS][2*DTT_SIZE];
///    __shared__ float window_vert_cos [NUM_COLORS][2*DTT_SIZE];
	// get correct kernel tile, then use 2 threads per kernel and image
	int   ktileX, ktileY;
	int   kernel_index; // common for all coors
	float kdx, kdy;
	if (threadIdx.x == 0){
		ktileX = min(KERNELS_HOR-1,  max(0, ((int) lrintf(centerX * (1.0/KERNELS_STEP)+1))));
		ktileY = min(KERNELS_VERT-1, max(0, ((int) lrintf(centerY * (1.0/KERNELS_STEP)+1))));
		kdx =    centerX - (ktileX << KERNELS_LSTEP) + (1 << (KERNELS_LSTEP -1)); // difference in pixel
		kdy =    centerY - (ktileY << KERNELS_LSTEP) + (1 << (KERNELS_LSTEP -1)); // difference in pixel
    	kernel_index = (ktileX + ktileY * KERNELS_HOR) * NUM_COLORS;
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
    float * kernel_src = gpu_kernels + kernel_full_index* (DTT_SIZE * DTT_SIZE * 4);
    float * kernelp =  clt_kernels;

    kernel_src += threadIdx.x; // lsb;
    kernelp +=    threadIdx.x; // lsb;
#pragma unroll
    for (int j = 0; j < DTT_SIZE * 4; j++){ // all 4 components, 8 rows
    	// shared memory kernels use DTT_SIZE1 (same as image data)
    	*kernelp = *kernel_src;
    	kernelp+=DTT_SIZE1;
    	kernel_src+=THREADSX;
/*
    	*kernelp = *kernel_src;
    	kernelp+=DTT_SIZE1;
    	kernel_src+=THREADSX;
    	*kernelp = *kernel_src;
    	kernelp+=DTT_SIZE1;
    	kernel_src+=THREADSX;
    	*kernelp = *kernel_src;
    	kernelp+=DTT_SIZE1;
    	kernel_src+=THREADSX;
    	*/
    }

    // Calculate offsets and prepare windows (all colors):
//	int kernel_full_index = kernel_index + color;
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



    // threads 0..23 loaded 3 color kernels, threads 24-27 - prepared hor and vert windows for R and B, threads 28..31 - for G
    // prepare, fold and write data to DTT buffers
    int dstride2 = dstride << 1; // in floats (pixels)
    int color0 = color & 1;
    int color1 = (color >>1) & 1;

    for (int gpass = 0; gpass < (color0 + 1); gpass++) { // Only once for R, B, twice - for G
    	int col_tl = int_topleft[0]; //  + (threadIdx.x << 1);
    	int row_tl = int_topleft[1];
    	//	int local_col = ((col_tl & 1) ^ BAYER_RED_COL ^ color0) + (threadIdx.x << 1);
    	//	int local_row = ((row_tl & 1) ^ BAYER_RED_ROW ^ color0);

    	// for red, blue and green, pass 0
    	int local_col = ((col_tl & 1) ^ (BAYER_RED_COL ^ color0 ^ color1 ^ gpass)) + (threadIdx.x << 1); // green red row: invert column from red
    	int local_row = ((row_tl & 1) ^ BAYER_RED_ROW ^ gpass);                            // use red row
    	float hwind_cos = window_hor_cos[local_col];
    	float hwind_sin = window_hor_sin[local_col]; // **** only used for green

    	int dtt_offset =     fold_indx2[local_row][local_col];
    	int dtt_offset_inc = fold_inc[local_row];
//    	float *dct_buf = (float *) clt_tile[ gpass << 1];
//    	float *dst_buf = (float *) clt_tile[(gpass << 1)+1];            // **** only used for green
    	float *dct_buf = clt_tile + ((gpass << 1) * (DTT_SIZE * DTT_SIZE1));
    	float *dst_buf = clt_tile + (((gpass << 1) + 1) * (DTT_SIZE * DTT_SIZE1));   // **** only used for green

    	if ((col_tl >= 0) && ((col_tl < (IMG_WIDTH - DTT_SIZE * 2))) && (row_tl >= 0) && ((row_tl < (IMG_HEIGHT - DTT_SIZE * 2)))) {
    		float *image_p = gpu_images + dstride * (row_tl + local_row)+ col_tl + local_col;
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
    		float *image_p = gpu_images + dstride * local_row+ (eff_col + local_col);
#pragma unroll
    		for (int i = 0; i < 8; i++) {
    			int eff_row = (min(IMG_WIDTH/2 - 1, max(0, row_pair + i)) << 1) + row_lsb;
    			float d =  image_p[dstride * eff_row] * window_vert_cos[local_row];

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
    	(*dtt_buf) += (*dtt_buf1); dtt_buf += (4 * DTT_SIZE1); dtt_buf1 += (4 * DTT_SIZE1);
    	(*dtt_buf) += (*dtt_buf1);
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
//    		clt_tile[0][threadIdx.x], // pointer to start of row
    		clt_tile + (DTT_SIZE1 * threadIdx.x), // [0][threadIdx.x], // pointer to start of row
			1); //int inc);
    if (color == BAYER_GREEN){
        dstiv_nodiverg( // all colors
        		clt_tile + DTT_SIZE1 * (threadIdx.x + DTT_SIZE), // clt_tile[1][threadIdx.x], // pointer to start of row
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
        dstiv_nodiverg( // all colors
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




}


