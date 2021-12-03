/**
 **
 ** dtt8x8.h
 **
 ** Copyright (C) 2018 Elphel, Inc.
 **
 ** -----------------------------------------------------------------------------**
 **
 **  dtt8x8.cuh is free software: you can redistribute it and/or modify
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
* \file dtt8x8.h
* \brief DCT-II, DST-II, DCT-IV and DST-IV for Complex Lapped Transform of 16x16 (stride 8)
*        in GPU
* This file contains building blocks for the 16x16 stride 8 COmplex Lapped Transform (CLT)
* implementation. DTT-IV are used for forward and inverse 2D CLT, DTT-II - to convert correlation
* results from the frequency to pixel domain. DTT-III (inverse of DTT-II) is not implemented
* here it is used to convert convolution kernels and LPF to the frequency domain - done in
* software.
*
* This file is cpompatible with both runtime and driver API, runtime is used for development
* with Nvidia Nsight, driver API when calling these kernels from Java
*/
#ifndef JCUDA
#define DTT_SIZE_LOG2                 3
#endif

#pragma once
#define DTT_SIZE                     (1 << DTT_SIZE_LOG2)
#define DTT_SIZE1        (DTT_SIZE + 1)
#define DTT_SIZE2        (2 * DTT_SIZE)
#define DTT_SIZE21       (DTT_SIZE2 + 1)
#define DTT_SIZE4        (4 * DTT_SIZE)
#define DTT_SIZE2M1      (DTT_SIZE2 - 1)
#define BAYER_RED   0
#define BAYER_BLUE  1
#define BAYER_GREEN 2
// assuming GR/BG as now
#define BAYER_RED_ROW 0
#define BAYER_RED_COL 1

#define DTTTEST_BLOCK_WIDTH          32
#define DTTTEST_BLOCK_HEIGHT         16
#define DTTTEST_BLK_STRIDE     (DTTTEST_BLOCK_WIDTH+1)

//extern __constant__ float idct_signs[4][4][4];
//extern __constant__ int imclt_indx9[16];
//extern __constant__ float HWINDOW2[];


// kernels (not used so far)
#if 0
extern "C" __global__ void GPU_DTT24_DRV(float *dst, float *src, int src_stride, int dtt_mode);
#endif// #if 0

//=========================== 2D functions ===============
extern __device__ void corrUnfoldTile(
		int corr_radius,
		float* qdata0, //    [4][DTT_SIZE][DTT_SIZE1], // 4 quadrants of the clt data, rows extended to optimize shared ports
		float* rslt);  //   [DTT_SIZE2M1][DTT_SIZE2M1]) // 15x15

extern __device__ void dttii_2d(
		float * clt_corr); // shared memory, [4][DTT_SIZE1][DTT_SIZE]

extern __device__ void dttiv_color_2d(
		float * clt_tile,
		int color);
extern __device__ void dttiv_mono_2d(
		float * clt_tile);
extern __device__ void imclt(
		float * clt_tile,   //        [4][DTT_SIZE][DTT_SIZE1], // +1 to alternate column ports [4][8][9]
		float * mclt_tile );

extern __device__ void imclt8threads(
		int     do_acc,     // 1 - add to previous value, 0 - overwrite
		float * clt_tile,   //        [4][DTT_SIZE][DTT_SIZE1], // +1 to alternate column ports [4][8][9]
		float * mclt_tile,  //           [2* DTT_SIZE][DTT_SIZE1+ DTT_SIZE], // +1 to alternate column ports[16][17]
		int     debug);
