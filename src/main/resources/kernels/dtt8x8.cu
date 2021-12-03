/**
 **
 ** dtt8x8.cuh
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
* \file dtt8x8.cuh
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
#include "dtt8x8.h"
#endif

//#define CUDART_INF_F            __int_as_float(0x7f800000)
/*
 Python code to generate constant coefficients:
def dct_constants():
    COSPI_1_8_SQRT2 = math.cos(math.pi/8)*math.sqrt(2.0)
    COSPI_3_8_SQRT2 = math.cos(3*math.pi/8)*math.sqrt(2.0)
    SQRT_2 = math.sqrt(2.0)
    SQRT1_2 = 1/math.sqrt(2.0)
    SQRT1_8 = 1/math.sqrt(8.0)
    CN = [[math.cos((2*k+1)*(math.pi/(8*(2 << t))))  for k in range (2 << t)] for t in range (2)]
    SN = [[math.sin((2*k+1)*(math.pi/(8*(2 << t))))  for k in range (2 << t)] for t in range (2)]
    print("__constant__ float COSPI_1_8_SQRT2 = %ff;"%(COSPI_1_8_SQRT2))
    print("__constant__ float COSPI_3_8_SQRT2 = %ff;"%(COSPI_3_8_SQRT2))
    print("__constant__ float SQRT_2 = %ff;"%         (SQRT_2))
    print("__constant__ float SQRT1_2 = %ff;"%        (SQRT1_2))
    print("__constant__ float SQRT1_8 = %ff;"%        (SQRT1_8))
    print("__constant__ float COSN1[] = {%ff,%ff};"%         (CN[0][0],CN[0][1]))
    print("__constant__ float COSN2[] = {%ff,%ff,%ff,%ff};"% (CN[1][0],CN[1][1],CN[1][2],CN[1][3]))
    print("__constant__ float SINN1[] = {%ff,%ff};"%         (SN[0][0],SN[0][1]))
    print("__constant__ float SINN2[] = {%ff,%ff,%ff,%ff};"% (SN[1][0],SN[1][1],SN[1][2],SN[1][3]))
*/
__constant__ float COSPI_1_8_SQRT2 = 1.306563f;
__constant__ float COSPI_3_8_SQRT2 = 0.541196f;
__constant__ float SQRT_2 = 1.414214f;
__constant__ float SQRT1_2 = 0.707107f;
__constant__ float SQRT1_8 = 0.353553f;
__constant__ float COSN1[] = {0.980785f,0.831470f};
__constant__ float COSN2[] = {0.995185f,0.956940f,0.881921f,0.773010f};
__constant__ float SINN1[] = {0.195090f,0.555570f};
__constant__ float SINN2[] = {0.098017f,0.290285f,0.471397f,0.634393f};
__constant__ int imclt_indx9[16] = {0x28,0x29,0x2a,0x2b,0x2b,0x2a,0x29,0x28,0x27,0x26,0x25,0x24,0x24,0x25,0x26,0x27};
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
__constant__ float HWINDOW2[] =  {0.049009f, 0.145142f, 0.235698f, 0.317197f,
                                  0.386505f, 0.440961f, 0.478470f, 0.497592f};


inline __device__ void dttii_shared_mem_nonortho(float * x0,  int inc, int dst_not_dct); // does not scale by y[0] (y[7]) by 1/sqrt[0]
inline __device__ void dttii_shared_mem(float * x0,  int inc, int dst_not_dct);   // used in GPU_DTT24_DRV
inline __device__ void dttiv_shared_mem(float * x0,  int inc, int dst_not_dct);   // used in GPU_DTT24_DRV
inline __device__ void dttiv_nodiverg  (float * x,   int inc, int dst_not_dct);   // not used
inline __device__ void dctiv_nodiverg  (float * x0,  int inc);                    // used in TP
inline __device__ void dstiv_nodiverg  (float * x0,  int inc);                    // used in TP

inline __device__ void dct_ii8         ( float x[8], float y[8]); // x,y point to 8-element arrays each // not used
inline __device__ void dct_iv8         ( float x[8], float y[8]); // x,y point to 8-element arrays each // not used
inline __device__ void dst_iv8         ( float x[8], float y[8]); // x,y point to 8-element arrays each // not used
inline __device__ void _dctii_nrecurs8 ( float x[8], float y[8]); // x,y point to 8-element arrays each // not used
inline __device__ void _dctiv_nrecurs8 ( float x[8], float y[8]); // x,y point to 8-element arrays each // not used


/**
**************************************************************************
*  Converts 2D image (in the GPU memory) using 8x8 DTT 8x8 tiles.
*  Mostly for testing and profiling individual conversions
*
* \param dst                        [OUT] - Coefficients as 8x8 tiles
* \param src                         [IN] - Source image of floats
* \param src_stride                  [IN] - Source image stride
* \param mode                        [IN] - DTT mode:
*     0 - horizontal DCT-IV followed by vertical DCT-IV
*     1 - horizontal DST-IV followed by vertical DCT-IV
*     2 - horizontal DCT-IV followed by vertical DST-IV
*     3 - horizontal DST-IV followed by vertical DST-IV
*     4 - horizontal DCT-II followed by vertical DCT-II
*     5 - horizontal DST-II followed by vertical DCT-II
*     6 - horizontal DCT-II followed by vertical DST-II
*     7 - horizontal DST-II followed by vertical DST-II
*
* \return None
*/
#ifdef BBBB
extern "C"
__global__ void GPU_DTT24_DRV(float *dst, float *src, int src_stride, int dtt_mode)
{
	int dtt_mode0 = dtt_mode & 1;
	int dtt_mode1 = (dtt_mode >>1) & 1;

    __shared__ float block[DTTTEST_BLOCK_HEIGHT * DTTTEST_BLK_STRIDE];

    int OffsThreadInRow = threadIdx.y * DTT_SIZE + threadIdx.x;
    int OffsThreadInCol = threadIdx.z * DTT_SIZE;
    src += ((blockIdx.y * DTTTEST_BLOCK_HEIGHT + OffsThreadInCol) * src_stride) + blockIdx.x * DTTTEST_BLOCK_WIDTH + OffsThreadInRow;
    dst += ((blockIdx.y * DTTTEST_BLOCK_HEIGHT + OffsThreadInCol) * src_stride) + blockIdx.x * DTTTEST_BLOCK_WIDTH + OffsThreadInRow;
    float *bl_ptr = block + OffsThreadInCol * DTTTEST_BLK_STRIDE + OffsThreadInRow;

#pragma unroll

    for (unsigned int i = 0; i < DTT_SIZE; i++)
        bl_ptr[i * DTTTEST_BLK_STRIDE] = src[i * src_stride];

    __syncthreads();
    // horizontal pass
    if (dtt_mode > 3) {
    	dttii_shared_mem                   (block + (OffsThreadInCol + threadIdx.x) * DTTTEST_BLK_STRIDE + OffsThreadInRow - threadIdx.x, 1, dtt_mode0);
    } else {
    	dttiv_shared_mem                   (block + (OffsThreadInCol + threadIdx.x) * DTTTEST_BLK_STRIDE + OffsThreadInRow - threadIdx.x, 1, dtt_mode0);
    }

    __syncthreads();
    // vertical pass
    if (dtt_mode > 3) {
    	dttii_shared_mem                    (bl_ptr, DTTTEST_BLK_STRIDE, dtt_mode1);
    } else {
    	dttiv_shared_mem                    (bl_ptr, DTTTEST_BLK_STRIDE, dtt_mode1);
    }
    __syncthreads();
    for (unsigned int i = 0; i < DTT_SIZE; i++)
        dst[i * src_stride] = bl_ptr[i * DTTTEST_BLK_STRIDE];
}
#endif //#ifdef BBBB



inline __device__ void _dctiv_nrecurs8( float x[8], float y[8]) // x,y point to 8-element arrays each
{
	float u00=            ( COSN2[0] * x[0] + SINN2[0] * x[7]);
	float u10=            (-SINN2[3] * x[3] + COSN2[3] * x[4]);

	float u01=            ( COSN2[1] * x[1] + SINN2[1] * x[6]);
	float u11=           -(-SINN2[2] * x[2] + COSN2[2] * x[5]);

	float u02=            ( COSN2[2] * x[2] + SINN2[2] * x[5]);
	float u12=            (-SINN2[1] * x[1] + COSN2[1] * x[6]);

	float u03=            ( COSN2[3] * x[3] + SINN2[3] * x[4]);
	float u13=           -(-SINN2[0] * x[0] + COSN2[0] * x[7]);

//	_dctii_nrecurs4(u00, u01, u02, u03, &v00, &v01, &v02, &v03);

	float ua00= u00 + u03;
	float ua10= u00 - u03;

	float ua01= u01 + u02;
	float ua11= u01 - u02;

	float v00= ua00 + ua01;
	float v02= ua00 - ua01;

	float v01= COSPI_1_8_SQRT2 * ua10 + COSPI_3_8_SQRT2 * ua11;
	float v03= COSPI_3_8_SQRT2 * ua10 - COSPI_1_8_SQRT2 * ua11;

//	_dctii_nrecurs4(u10, u11, u12, u13, &v10, &v11, &v12, &v13);

	float ub00= u10 + u13;
	float ub10= u10 - u13;

	float ub01= u11 + u12;
	float ub11= u11 - u12;

	float vb00= ub00 + ub01;
	float vb01= ub00 - ub01;

	float vb10= COSPI_1_8_SQRT2*ub10 + COSPI_3_8_SQRT2*ub11;
	float vb11= COSPI_3_8_SQRT2*ub10 - COSPI_1_8_SQRT2*ub11;


	y[0] =  SQRT_2 * v00;    // w0[0];
	y[1] =  v01 -  vb11;    // w1[0];
	// j == 1
	y[2] =  v01 +  vb11;    // w0[1];
	y[3] =  v02 +  vb01;    // w1[1];
	// j == 2
	y[4] =  v02 -  vb01;    // w0[2];
	y[5] =  v03 -  vb10;    // w1[2]; - same as y[3]
	// j == 3
	y[6] =  v03 +  vb10;    // w0[3];
	y[7] =  SQRT_2 * vb00;    // w1[3];
}

__device__ void _dttiv(float x0, float x1,float x2, float x3,float x4, float x5,float x6, float x7,
		float *y0, float *y1, float *y2, float *y3, float *y4, float *y5, float *y6, float *y7, int dst_not_dct)
{
	float u00, u01, u02, u03, u10, u11, u12, u13;
	if (dst_not_dct) { // DSTIV
		u00=  ( COSN2[0] * x7 + SINN2[0] * x0);
		u10=  (-SINN2[3] * x4 + COSN2[3] * x3);

		u01=  ( COSN2[1] * x6 + SINN2[1] * x1);
		u11= -(-SINN2[2] * x5 + COSN2[2] * x2);

		u02=  ( COSN2[2] * x5 + SINN2[2] * x2);
		u12=  (-SINN2[1] * x6 + COSN2[1] * x1);

		u03=  ( COSN2[3] * x4 + SINN2[3] * x3);
		u13= -(-SINN2[0] * x7 + COSN2[0] * x0);
	} else { // DCTIV
		u00=  ( COSN2[0] * x0 + SINN2[0] * x7);
		u10=  (-SINN2[3] * x3 + COSN2[3] * x4);

		u01=  ( COSN2[1] * x1 + SINN2[1] * x6);
		u11= -(-SINN2[2] * x2 + COSN2[2] * x5);

		u02=  ( COSN2[2] * x2 + SINN2[2] * x5);
		u12=  (-SINN2[1] * x1 + COSN2[1] * x6);

		u03=  ( COSN2[3] * x3 + SINN2[3] * x4);
		u13= -(-SINN2[0] * x0 + COSN2[0] * x7);
	}

//	_dctii_nrecurs4(u00, u01, u02, u03, &v00, &v01, &v02, &v03);

	float ua00= u00 + u03;
	float ua10= u00 - u03;

	float ua01= u01 + u02;
	float ua11= u01 - u02;

	float v00= ua00 + ua01;
	float v02= ua00 - ua01;

	float v01= COSPI_1_8_SQRT2 * ua10 + COSPI_3_8_SQRT2 * ua11;
	float v03= COSPI_3_8_SQRT2 * ua10 - COSPI_1_8_SQRT2 * ua11;

//	_dctii_nrecurs4(u10, u11, u12, u13, &v10, &v11, &v12, &v13);

	float ub00= u10 + u13;
	float ub10= u10 - u13;

	float ub01= u11 + u12;
	float ub11= u11 - u12;

	float vb00= ub00 + ub01;
	float vb01= ub00 - ub01;

	float vb10= COSPI_1_8_SQRT2*ub10 + COSPI_3_8_SQRT2*ub11;
	float vb11= COSPI_3_8_SQRT2*ub10 - COSPI_1_8_SQRT2*ub11;

	*y0 =  v00 * 0.5f;              // w0[0];
	// j == 1
	*y2 =  (v01 +  vb11) * SQRT1_8; // w0[1];
	// j == 2
	*y4 =  (v02 -  vb01) * SQRT1_8; // w0[2];
	// j == 3
	*y6 =  (v03 +  vb10) * SQRT1_8; // w0[3];
	if (dst_not_dct) { // DSTIV
		*y1 =  (vb11 - v01)  * SQRT1_8; // w1[0];
		*y3 = -(v02 +  vb01) * SQRT1_8; // w1[1];
		*y5 =  (vb10 - v03)  * SQRT1_8; // w1[2]; - same as y[3]
		*y7 = -vb00 * 0.5f;             // w1[3];
	} else {
		*y1 =  (v01 -  vb11) * SQRT1_8; // w1[0];
		*y3 =  (v02 +  vb01) * SQRT1_8; // w1[1];
		*y5 =  (v03 -  vb10) * SQRT1_8; // w1[2]; - same as y[3]
		*y7 =  vb00 * 0.5f;             // w1[3];
	}
}

inline __device__ void dttii_shared_mem(float * x0,  int inc, int dst_not_dct)
{
	float *x1 = x0 + inc;
	float *x2 = x1 + inc;
	float *x3 = x2 + inc;
	float *x4 = x3 + inc;
	float *x5 = x4 + inc;
	float *x6 = x5 + inc;
	float *x7 = x6 + inc;
	float u00, u01, u02, u03, u10, u11, u12, u13;
	if (dst_not_dct) { // DSTII
		// invert odd input samples
		u00= ( (*x0) - (*x7));
		u10= ( (*x0) + (*x7));

		u01= (-(*x1) + (*x6));
		u11= (-(*x1) - (*x6));

		u02= ( (*x2) - (*x5));
		u12= ( (*x2) + (*x5));

		u03= (-(*x3) + (*x4));
		u13= (-(*x3) - (*x4));
	} else { // DCTII
		u00= ( (*x0) + (*x7));
		u10= ( (*x0) - (*x7));

		u01= ( (*x1) + (*x6));
		u11= ( (*x1) - (*x6));

		u02= ( (*x2) + (*x5));
		u12= ( (*x2) - (*x5));

		u03= ( (*x3) + (*x4));
		u13= ( (*x3) - (*x4));
	}
	//	_dctii_nrecurs4(u00,u01, u02, u03, &v00, &v01, &v02, &v03);

		float w00= u00 + u03;
		float w10= u00 - u03;

		float w01= (u01 + u02);
		float w11= (u01 - u02);

		float v01= COSPI_1_8_SQRT2 * w10 + COSPI_3_8_SQRT2 * w11;
		float v03= COSPI_3_8_SQRT2 * w10 - COSPI_1_8_SQRT2 * w11;
	//	_dctiv_nrecurs4(u10, u11, u12, u13, &v10, &v11, &v12, &v13);
		float w20=            ( COSN1[0] * u10 + SINN1[0] * u13);
		float w30=            (-SINN1[1] * u11 + COSN1[1] * u12);

		float w21=            ( COSN1[1] * u11 + SINN1[1] * u12);
		float w31=           -(-SINN1[0] * u10 + COSN1[0] * u13);
		float v11 = w20 - w21 - w30 + w31;
		float v12 = w20 - w21 + w30 - w31;

	if (dst_not_dct) { // DSTII
		// Invert output sequence
		*x0 =   (w30 + w31)*  0.5f;    // v13 * SQRT1_8; z10 * 0.5f
		*x1 =   v03 *         SQRT1_8;

		*x2 =   v12 *         SQRT1_8;
		*x3 =   (w00 - w01) * SQRT1_8; // v02 * SQRT1_8

		*x4 =   v11 *         SQRT1_8;
		*x5 =   v01 *         SQRT1_8;

		*x6 =   (w20 + w21) * 0.5f;    // v10 * SQRT1_8; z00 * 0.5f;
		*x7 =   (w00 + w01) * SQRT1_8; // v00 * SQRT1_8
	} else {
		*x0 =   (w00 + w01) * SQRT1_8; // v00 * SQRT1_8
		*x1 =   (w20 + w21) * 0.5f;    // v10 * SQRT1_8; z00 * 0.5f;

		*x2 =   v01 *         SQRT1_8;
		*x3 =   v11 *         SQRT1_8;

		*x4 =   (w00 - w01) * SQRT1_8; // v02 * SQRT1_8
		*x5 =   v12 *         SQRT1_8;

		*x6 =   v03 *         SQRT1_8;
		*x7 =   (w30 + w31)*  0.5f;    // v13 * SQRT1_8; z10 * 0.5f
	}
}

inline __device__ void dttii_shared_mem_nonortho(float * x0,  int inc, int dst_not_dct)
{
	float *x1 = x0 + inc;
	float *x2 = x1 + inc;
	float *x3 = x2 + inc;
	float *x4 = x3 + inc;
	float *x5 = x4 + inc;
	float *x6 = x5 + inc;
	float *x7 = x6 + inc;
	float u00, u01, u02, u03, u10, u11, u12, u13;
	if (dst_not_dct) { // DSTII
		// invert odd input samples
		u00= ( (*x0) - (*x7));
		u10= ( (*x0) + (*x7));

		u01= (-(*x1) + (*x6));
		u11= (-(*x1) - (*x6));

		u02= ( (*x2) - (*x5));
		u12= ( (*x2) + (*x5));

		u03= (-(*x3) + (*x4));
		u13= (-(*x3) - (*x4));
	} else { // DCTII
		u00= ( (*x0) + (*x7));
		u10= ( (*x0) - (*x7));

		u01= ( (*x1) + (*x6));
		u11= ( (*x1) - (*x6));

		u02= ( (*x2) + (*x5));
		u12= ( (*x2) - (*x5));

		u03= ( (*x3) + (*x4));
		u13= ( (*x3) - (*x4));
	}
	//	_dctii_nrecurs4(u00,u01, u02, u03, &v00, &v01, &v02, &v03);

		float w00= u00 + u03;
		float w10= u00 - u03;

		float w01= (u01 + u02);
		float w11= (u01 - u02);

		float v01= COSPI_1_8_SQRT2 * w10 + COSPI_3_8_SQRT2 * w11;
		float v03= COSPI_3_8_SQRT2 * w10 - COSPI_1_8_SQRT2 * w11;
	//	_dctiv_nrecurs4(u10, u11, u12, u13, &v10, &v11, &v12, &v13);
		float w20=            ( COSN1[0] * u10 + SINN1[0] * u13);
		float w30=            (-SINN1[1] * u11 + COSN1[1] * u12);

		float w21=            ( COSN1[1] * u11 + SINN1[1] * u12);
		float w31=           -(-SINN1[0] * u10 + COSN1[0] * u13);
		float v11 = w20 - w21 - w30 + w31;
		float v12 = w20 - w21 + w30 - w31;

	if (dst_not_dct) { // DSTII
		// Invert output sequence
		*x0 =   (w30 + w31)*  0.5f;    // v13 * SQRT1_8; z10 * 0.5f
		*x1 =   v03 *         SQRT1_8;

		*x2 =   v12 *         SQRT1_8;
		*x3 =   (w00 - w01) * SQRT1_8; // v02 * SQRT1_8

		*x4 =   v11 *         SQRT1_8;
		*x5 =   v01 *         SQRT1_8;

		*x6 =   (w20 + w21) * 0.5f;    // v10 * SQRT1_8; z00 * 0.5f;
		*x7 =   (w00 + w01) * 0.5f;    // SQRT1_8; // v00 * SQRT1_8 //*** no 1/sqrt(2)!
	} else {
		*x0 =   (w00 + w01) * 0.5f;    // SQRT1_8; // v00 * SQRT1_8 //*** no 1/sqrt(2)!
		*x1 =   (w20 + w21) * 0.5f;    // v10 * SQRT1_8; z00 * 0.5f;

		*x2 =   v01 *         SQRT1_8;
		*x3 =   v11 *         SQRT1_8;

		*x4 =   (w00 - w01) * SQRT1_8; // v02 * SQRT1_8
		*x5 =   v12 *         SQRT1_8;

		*x6 =   v03 *         SQRT1_8;
		*x7 =   (w30 + w31)*  0.5f;    // v13 * SQRT1_8; z10 * 0.5f
	}
}

inline __device__ void dttiv_shared_mem(float * x0,  int inc, int dst_not_dct)
{
	float *x1 = x0 + inc;
	float *x2 = x1 + inc;
	float *x3 = x2 + inc;
	float *x4 = x3 + inc;
	float *x5 = x4 + inc;
	float *x6 = x5 + inc;
	float *x7 = x6 + inc;
	float u00, u01, u02, u03, u10, u11, u12, u13;
	if (dst_not_dct) { // DSTIV
		u00=  ( COSN2[0] * (*x7) + SINN2[0] * (*x0));
		u10=  (-SINN2[3] * (*x4) + COSN2[3] * (*x3));

		u01=  ( COSN2[1] * (*x6) + SINN2[1] * (*x1));
		u11= -(-SINN2[2] * (*x5) + COSN2[2] * (*x2));

		u02=  ( COSN2[2] * (*x5) + SINN2[2] * (*x2));
		u12=  (-SINN2[1] * (*x6) + COSN2[1] * (*x1));

		u03=  ( COSN2[3] * (*x4) + SINN2[3] * (*x3));
		u13= -(-SINN2[0] * (*x7) + COSN2[0] * (*x0));
	} else { // DCTIV
		u00=  ( COSN2[0] * (*x0) + SINN2[0] * (*x7));
		u10=  (-SINN2[3] * (*x3) + COSN2[3] * (*x4));

		u01=  ( COSN2[1] * (*x1) + SINN2[1] * (*x6));
		u11= -(-SINN2[2] * (*x2) + COSN2[2] * (*x5));

		u02=  ( COSN2[2] * (*x2) + SINN2[2] * (*x5));
		u12=  (-SINN2[1] * (*x1) + COSN2[1] * (*x6));

		u03=  ( COSN2[3] * (*x3) + SINN2[3] * (*x4));
		u13= -(-SINN2[0] * (*x0) + COSN2[0] * (*x7));
	}

//	_dctii_nrecurs4(u00, u01, u02, u03, &v00, &v01, &v02, &v03);

	float ua00= u00 + u03;
	float ua10= u00 - u03;

	float ua01= u01 + u02;
	float ua11= u01 - u02;

	float v00= ua00 + ua01;
	float v02= ua00 - ua01;

	float v01= COSPI_1_8_SQRT2 * ua10 + COSPI_3_8_SQRT2 * ua11;
	float v03= COSPI_3_8_SQRT2 * ua10 - COSPI_1_8_SQRT2 * ua11;

//	_dctii_nrecurs4(u10, u11, u12, u13, &v10, &v11, &v12, &v13);

	float ub00= u10 + u13;
	float ub10= u10 - u13;

	float ub01= u11 + u12;
	float ub11= u11 - u12;

	float vb00= ub00 + ub01;
	float vb01= ub00 - ub01;

	float vb10= COSPI_1_8_SQRT2*ub10 + COSPI_3_8_SQRT2*ub11;
	float vb11= COSPI_3_8_SQRT2*ub10 - COSPI_1_8_SQRT2*ub11;


	*x0 =  v00 * 0.5f;              // w0[0];
	*x2 =  (v01 +  vb11) * SQRT1_8; // w0[1];
	*x4 =  (v02 -  vb01) * SQRT1_8; // w0[2];
	*x6 =  (v03 +  vb10) * SQRT1_8; // w0[3];
	if (dst_not_dct) { // DSTIV
		*x1 =  (vb11 - v01)  * SQRT1_8; // w1[0];
		*x3 = -(v02 +  vb01) * SQRT1_8; // w1[1];
		*x5 =  (vb10 - v03)  * SQRT1_8; // w1[2]; - same as y[3]
		*x7 = -vb00 * 0.5f;             // w1[3];
	} else {
		*x1 =  (v01 -  vb11) * SQRT1_8; // w1[0];
		*x3 =  (v02 +  vb01) * SQRT1_8; // w1[1];
		*x5 =  (v03 -  vb10) * SQRT1_8; // w1[2]; - same as y[3]
		*x7 =  vb00 * 0.5f;             // w1[3];
	}
}

inline __device__ void dttiv_nodiverg(float * x,  int inc, int dst_not_dct)
{
	float sgn = 1 - 2* dst_not_dct;
	float *y0 = x;
	float *y1 = y0 + inc;
	float *y2 = y1 + inc;
	float *y3 = y2 + inc;
	float *y4 = y3 + inc;
	float *y5 = y4 + inc;
	float *y6 = y5 + inc;
	float *y7 = y6 + inc;

	float *x0 =  x + dst_not_dct * 7 * inc;
	// negate inc, replace
	inc *= sgn;
	float *x1 = x0 + inc;
	float *x2 = x1 + inc;
	float *x3 = x2 + inc;
	float *x4 = x3 + inc;
	float *x5 = x4 + inc;
	float *x6 = x5 + inc;
	float *x7 = x6 + inc;
	float u00, u01, u02, u03, u10, u11, u12, u13;
	u00=  ( COSN2[0] * (*x0) + SINN2[0] * (*x7));
	u10=  (-SINN2[3] * (*x3) + COSN2[3] * (*x4));

	u01=  ( COSN2[1] * (*x1) + SINN2[1] * (*x6));
	u11= -(-SINN2[2] * (*x2) + COSN2[2] * (*x5));

	u02=  ( COSN2[2] * (*x2) + SINN2[2] * (*x5));
	u12=  (-SINN2[1] * (*x1) + COSN2[1] * (*x6));

	u03=  ( COSN2[3] * (*x3) + SINN2[3] * (*x4));
	u13= -(-SINN2[0] * (*x0) + COSN2[0] * (*x7));

//	_dctii_nrecurs4(u00, u01, u02, u03, &v00, &v01, &v02, &v03);

	float ua00= u00 + u03;
	float ua10= u00 - u03;

	float ua01= u01 + u02;
	float ua11= u01 - u02;

	float v00= ua00 + ua01;
	float v02= ua00 - ua01;

	float v01= COSPI_1_8_SQRT2 * ua10 + COSPI_3_8_SQRT2 * ua11;
	float v03= COSPI_3_8_SQRT2 * ua10 - COSPI_1_8_SQRT2 * ua11;

//	_dctii_nrecurs4(u10, u11, u12, u13, &v10, &v11, &v12, &v13);

	float ub00= u10 + u13;
	float ub10= u10 - u13;

	float ub01= u11 + u12;
	float ub11= u11 - u12;

	float vb00= ub00 + ub01;
	float vb01= ub00 - ub01;

	float vb10= COSPI_1_8_SQRT2*ub10 + COSPI_3_8_SQRT2*ub11;
	float vb11= COSPI_3_8_SQRT2*ub10 - COSPI_1_8_SQRT2*ub11;


	*y0 =  v00 * 0.5f;              // w0[0];
	*y2 =  (v01 +  vb11) * SQRT1_8; // w0[1];
	*y4 =  (v02 -  vb01) * SQRT1_8; // w0[2];
	*y6 =  (v03 +  vb10) * SQRT1_8; // w0[3];
	*y1 =  sgn * (v01 -  vb11) * SQRT1_8; // w1[0];
	*y3 =  sgn * (v02 +  vb01) * SQRT1_8; // w1[1];
	*y5 =  sgn * (v03 -  vb10) * SQRT1_8; // w1[2]; - same as y[3]
	*y7 =  sgn * vb00 * 0.5f;             // w1[3];
}

inline __device__ void dctiv_nodiverg(float * x0,  int inc)
{
	float *x1 = x0 + inc;
	float *x2 = x1 + inc;
	float *x3 = x2 + inc;
	float *x4 = x3 + inc;
	float *x5 = x4 + inc;
	float *x6 = x5 + inc;
	float *x7 = x6 + inc;
	float u00, u01, u02, u03, u10, u11, u12, u13;
	u00=  ( COSN2[0] * (*x0) + SINN2[0] * (*x7));
	u10=  (-SINN2[3] * (*x3) + COSN2[3] * (*x4));

	u01=  ( COSN2[1] * (*x1) + SINN2[1] * (*x6));
	u11= -(-SINN2[2] * (*x2) + COSN2[2] * (*x5));

	u02=  ( COSN2[2] * (*x2) + SINN2[2] * (*x5));
	u12=  (-SINN2[1] * (*x1) + COSN2[1] * (*x6));

	u03=  ( COSN2[3] * (*x3) + SINN2[3] * (*x4));
	u13= -(-SINN2[0] * (*x0) + COSN2[0] * (*x7));

//	_dctii_nrecurs4(u00, u01, u02, u03, &v00, &v01, &v02, &v03);

	float ua00= u00 + u03;
	float ua10= u00 - u03;

	float ua01= u01 + u02;
	float ua11= u01 - u02;

	float v00= ua00 + ua01;
	float v02= ua00 - ua01;

	float v01= COSPI_1_8_SQRT2 * ua10 + COSPI_3_8_SQRT2 * ua11;
	float v03= COSPI_3_8_SQRT2 * ua10 - COSPI_1_8_SQRT2 * ua11;

//	_dctii_nrecurs4(u10, u11, u12, u13, &v10, &v11, &v12, &v13);

	float ub00= u10 + u13;
	float ub10= u10 - u13;

	float ub01= u11 + u12;
	float ub11= u11 - u12;

	float vb00= ub00 + ub01;
	float vb01= ub00 - ub01;

	float vb10= COSPI_1_8_SQRT2*ub10 + COSPI_3_8_SQRT2*ub11;
	float vb11= COSPI_3_8_SQRT2*ub10 - COSPI_1_8_SQRT2*ub11;


	*x0 =  v00 * 0.5f;              // w0[0];
	*x2 =  (v01 +  vb11) * SQRT1_8; // w0[1];
	*x4 =  (v02 -  vb01) * SQRT1_8; // w0[2];
	*x6 =  (v03 +  vb10) * SQRT1_8; // w0[3];
	*x1 =  (v01 -  vb11) * SQRT1_8; // w1[0];
	*x3 =  (v02 +  vb01) * SQRT1_8; // w1[1];
	*x5 =  (v03 -  vb10) * SQRT1_8; // w1[2]; - same as y[3]
	*x7 =   vb00 * 0.5f;             // w1[3];
}

inline __device__ void dstiv_nodiverg(float * x,  int inc)
{
	float *x0 =  x +  7 * inc;
	// negate inc, replace
	inc = -inc;
	float *x1 = x0 + inc;
	float *x2 = x1 + inc;
	float *x3 = x2 + inc;
	float *x4 = x3 + inc;
	float *x5 = x4 + inc;
	float *x6 = x5 + inc;
	float *x7 = x6 + inc;
	float u00, u01, u02, u03, u10, u11, u12, u13;
	u00=  ( COSN2[0] * (*x0) + SINN2[0] * (*x7));
	u10=  (-SINN2[3] * (*x3) + COSN2[3] * (*x4));

	u01=  ( COSN2[1] * (*x1) + SINN2[1] * (*x6));
	u11= -(-SINN2[2] * (*x2) + COSN2[2] * (*x5));

	u02=  ( COSN2[2] * (*x2) + SINN2[2] * (*x5));
	u12=  (-SINN2[1] * (*x1) + COSN2[1] * (*x6));

	u03=  ( COSN2[3] * (*x3) + SINN2[3] * (*x4));
	u13= -(-SINN2[0] * (*x0) + COSN2[0] * (*x7));

//	_dctii_nrecurs4(u00, u01, u02, u03, &v00, &v01, &v02, &v03);

	float ua00= u00 + u03;
	float ua10= u00 - u03;

	float ua01= u01 + u02;
	float ua11= u01 - u02;

	float v00= ua00 + ua01;
	float v02= ua00 - ua01;

	float v01= COSPI_1_8_SQRT2 * ua10 + COSPI_3_8_SQRT2 * ua11;
	float v03= COSPI_3_8_SQRT2 * ua10 - COSPI_1_8_SQRT2 * ua11;

//	_dctii_nrecurs4(u10, u11, u12, u13, &v10, &v11, &v12, &v13);

	float ub00= u10 + u13;
	float ub10= u10 - u13;

	float ub01= u11 + u12;
	float ub11= u11 - u12;

	float vb00= ub00 + ub01;
	float vb01= ub00 - ub01;

	float vb10= COSPI_1_8_SQRT2*ub10 + COSPI_3_8_SQRT2*ub11;
	float vb11= COSPI_3_8_SQRT2*ub10 - COSPI_1_8_SQRT2*ub11;


	*x7 =  v00 * 0.5f;              // w0[0];
	*x5 =  (v01 +  vb11) * SQRT1_8; // w0[1];
	*x3 =  (v02 -  vb01) * SQRT1_8; // w0[2];
	*x1 =  (v03 +  vb10) * SQRT1_8; // w0[3];

	*x6 =  (vb11 - v01)  * SQRT1_8; // w1[0];
	*x4 = -(v02 +  vb01) * SQRT1_8; // w1[1];
	*x2 =  (vb10 - v03)  * SQRT1_8; // w1[2]; - same as y[3]
	*x0 = -vb00 * 0.5f;             // w1[3];
}



inline  __device__ void _dctii_nrecurs8( float x[8], float y[8]) // x,y point to 8-element arrays each
{
	float u00= (x[0] + x[7]);
	float u10= (x[0] - x[7]);

	float u01= (x[1] + x[6]);
	float u11= (x[1] - x[6]);

	float u02= (x[2] + x[5]);
	float u12= (x[2] - x[5]);

	float u03= (x[3] + x[4]);
	float u13= (x[3] - x[4]);

//	_dctii_nrecurs4(u00, u01, u02, u03, &v00, &v01, &v02, &v03);

	float w00= u00 + u03;
	float w10= u00 - u03;

	float w01= (u01 + u02);
	float w11= (u01 - u02);

	float v00= w00 + w01;
	float v02= w00 - w01;
	float v01= COSPI_1_8_SQRT2 * w10 + COSPI_3_8_SQRT2 * w11;
	float v03= COSPI_3_8_SQRT2 * w10 - COSPI_1_8_SQRT2 * w11;

//	_dctiv_nrecurs4(u10, u11, u12, u13, &v10, &v11, &v12, &v13);
	float w20=            ( COSN1[0] * u10 + SINN1[0] * u13);
	float w30=            (-SINN1[1] * u11 + COSN1[1] * u12);

	float w21=            ( COSN1[1] * u11 + SINN1[1] * u12);
	float w31=           -(-SINN1[0] * u10 + COSN1[0] * u13);

//	_dctii_nrecurs2(u00, u01, &v00, &v01);
	float z00= w20 + w21;
	float z01= w20 - w21;

//	_dctii_nrecurs2(u10, u11, &v10, &v11);
	float z10= w30 + w31;
	float z11= w30 - w31;

	float v10 = SQRT_2 * z00;
	float v11 = z01 - z11;

	float v12 = z01 + z11;
	float v13 = SQRT_2 * z10;

	y[0] =   v00;
	y[1] =   v10;

	y[2] =   v01;
	y[3] =   v11;

	y[4] =   v02;
	y[5] =   v12;

	y[6] =   v03;
	y[7] =   v13;
}

inline  __device__ void dct_ii8( float x[8], float y[8]) // x,y point to 8-element arrays each
{
	_dctii_nrecurs8(x, y);
#pragma unroll
	for (int i = 0; i < 8 ; i++) {
		y[i] *= SQRT1_8;
	}
}


__device__ void dct_iv8( float x[8], float y[8]) // x,y point to 8-element arrays each
{
	_dctiv_nrecurs8(x, y);
#pragma unroll
	for (int i = 0; i < 8 ; i++) {
		y[i] *= SQRT1_8;
	}

}

inline __device__ void dst_iv8( float x[8], float y[8]) // x,y point to 8-element arrays each
{
	float xr[8];
#pragma unroll
	for (int i=0; i < 8;i++){
		xr[i] = x[7 - i];
	}
	_dctiv_nrecurs8(xr, y);
#pragma unroll
	for (int i=0; i < 8;i+=2){
		y[i]   *=  SQRT1_8;
		y[i+1] *= -SQRT1_8;
	}
}


//=========================== 2D functions ===============
__device__ void corrUnfoldTile(
		int corr_radius,
		float* qdata0, //    [4][DTT_SIZE][DTT_SIZE1], // 4 quadrants of the clt data, rows extended to optimize shared ports
		float* rslt)  //   [DTT_SIZE2M1][DTT_SIZE2M1]) // 15x15
{
	int size2r1 = 2 * corr_radius + 1; // 15
	int crp1 = corr_radius + 1;        //8
///	const int rslt_base_index = DTT_SIZE2M1 * (DTT_SIZE) - DTT_SIZE; // offset of the center
	int rslt_base_index = size2r1 * crp1 - crp1; // offset of the center

	float * qdata1 = qdata0 + (DTT_SIZE * DTT_SIZE1);
	float * qdata2 = qdata1 + (DTT_SIZE * DTT_SIZE1);
	float * qdata3 = qdata2 + (DTT_SIZE * DTT_SIZE1);
	int i = threadIdx.x;
	if (i > corr_radius) {
		return; // not needed, only use inner
	}
//	printf("\corrUnfoldTile() corr_radius=%d, i=%d\n",corr_radius,i);
	float corr_pixscale = 0.25f;
	int i_transform_size = i * DTT_SIZE1; // used to address source rows which are 9 long
	int im1_transform_size = i_transform_size - DTT_SIZE1; // negative for i = 0, use only after divergence
///	int rslt_row_offs = i * DTT_SIZE2M1;
	int rslt_row_offs = i * size2r1;
	int rslt_base_index_p = rslt_base_index + rslt_row_offs; // i * DTT_SIZE2M1;
	int rslt_base_index_m = rslt_base_index - rslt_row_offs; // i * DTT_SIZE2M1;
	rslt[rslt_base_index_p] = corr_pixscale * qdata0[i_transform_size]; // incomplete, will only be used for thread i=0
	rslt[rslt_base_index_m] = rslt[rslt_base_index_p];                  // nop for i=0 incomplete, will only be used for thread i=0
///	for (int j = 1; j < DTT_SIZE; j++) {
	for (int j = 1; j <= corr_radius; j++) {
		int rslt_base_index_pp = rslt_base_index_p + j;
		int rslt_base_index_pm = rslt_base_index_p - j;
		rslt[rslt_base_index_pp] = corr_pixscale * (
				 qdata0[i_transform_size + j] +
				 qdata1[i_transform_size + j -1]); // incomplete, will only be used for thread i=0
		rslt[rslt_base_index_pm] = corr_pixscale * (
				 qdata0[i_transform_size + j] +
				-qdata1[i_transform_size + j -1]); // incomplete, will only be used for thread i=0
	}
	if (i == 0) {
		return;
	}
///	im1_transform_size = i_transform_size - DTT_SIZE1; // already is calculated
	float d = corr_pixscale * qdata2[im1_transform_size];
	rslt[rslt_base_index_p] += d;
	rslt[rslt_base_index_m] -= d;
	for (int j = 1; j <= corr_radius; j++) {
		int rslt_base_index_pp = rslt_base_index_p + j;
		int rslt_base_index_pm = rslt_base_index_p - j;
		int rslt_base_index_mp = rslt_base_index_m + j;
		int rslt_base_index_mm = rslt_base_index_m - j;
		float d2 = corr_pixscale * qdata2[im1_transform_size + j];
		float d3 = corr_pixscale * qdata3[im1_transform_size + j -1];
		//rslt[rslt_base_index_mp], rslt[rslt_base_index_mp] are partially calculated in the cycle common with i=0
		rslt[rslt_base_index_mp] = rslt[rslt_base_index_pp] - d2 - d3;
		rslt[rslt_base_index_mm] = rslt[rslt_base_index_pm] - d2 + d3;
		rslt[rslt_base_index_pp] += d2 + d3;
		rslt[rslt_base_index_pm] += d2 - d3;
	}
}

__device__ void dttii_2d(
		float * clt_corr) // shared memory, [4][DTT_SIZE1][DTT_SIZE]
{
    // change to 16-32 threads?? in next iteration
    // vert pass (hor pass in Java, before transpose. Here transposed, no transform needed)
    for (int q = 0; q < 4; q++){
    	int is_sin = (q >> 1) & 1;
    	dttii_shared_mem_nonortho(clt_corr + q * (DTT_SIZE1 * DTT_SIZE) + threadIdx.x , DTT_SIZE1, is_sin); // vertical pass, thread is column
    }
    __syncthreads();

    // hor pass, corresponding to vert pass in Java
    for (int q = 0; q < 4; q++){
    	int is_sin = q & 1;
    	dttii_shared_mem_nonortho(clt_corr + (q * DTT_SIZE + threadIdx.x) * DTT_SIZE1 ,  1, is_sin); // horizontal pass, tread is row
    }
    __syncthreads();

}

__device__ void dttiv_color_2d(
		float * clt_tile,
		int color)
{
    dctiv_nodiverg( // all colors
			clt_tile + (DTT_SIZE1 * threadIdx.x), // [0][threadIdx.x], // pointer to start of row
			1); //int inc);
//	__syncthreads();// worsened
    if (color == BAYER_GREEN){
        dstiv_nodiverg( // all colors
				clt_tile + DTT_SIZE1 * threadIdx.x + DTT_SIZE1 * DTT_SIZE, // clt_tile[1][threadIdx.x], // pointer to start of row
    			1); //int inc);

    }
  	 __syncthreads();// __syncwarp();

#ifdef DEBUG222
    if ((threadIdx.x) == 0){
        printf("\nDTT Tiles after horizontal pass, color=%d\n",color);
    	debug_print_clt1(clt_tile, color, (color== BAYER_GREEN)?3:1); // only 1 quadrant for R,B and 2 - for G
    }
     __syncthreads();// __syncwarp();
#endif
    dctiv_nodiverg( // all colors
    		clt_tile + threadIdx.x, //  &clt_tile[0][0][threadIdx.x], // pointer to start of column
			DTT_SIZE1); // int inc,
//	__syncthreads();// worsened
    if (color == BAYER_GREEN){
          dctiv_nodiverg( // all colors
        		clt_tile + threadIdx.x + (DTT_SIZE1 * DTT_SIZE), // &clt_tile[1][0][threadIdx.x], // pointer to start of column
    			DTT_SIZE1); // int inc,
    }
  	 __syncthreads();// __syncwarp();
}

__device__ void dttiv_mono_2d(
		float * clt_tile)
{
	// Copy 0-> 1

    dctiv_nodiverg(
			clt_tile + (DTT_SIZE1 * threadIdx.x) + (0 * DTT_SIZE1 * DTT_SIZE),
			1); //int inc);
    dstiv_nodiverg(
    		clt_tile + (DTT_SIZE1 * threadIdx.x) + (1 * DTT_SIZE1 * DTT_SIZE),
			1); //int inc);
    dctiv_nodiverg(
			clt_tile + (DTT_SIZE1 * threadIdx.x) + (2 * DTT_SIZE1 * DTT_SIZE),
			1); //int inc);
    dstiv_nodiverg(
    		clt_tile + (DTT_SIZE1 * threadIdx.x) + (3 * DTT_SIZE1 * DTT_SIZE),
			1); //int inc);
	__syncthreads();// __syncwarp();

#ifdef DEBUG222
    if ((threadIdx.x) == 0){
        printf("\nDTT Tiles after horizontal pass, color=%d\n",color);
    	debug_print_clt1(clt_tile, color, (color== BAYER_GREEN)?3:1); // only 1 quadrant for R,B and 2 - for G
    }
     __syncthreads();// __syncwarp();
#endif

 	dctiv_nodiverg( // CC
 			clt_tile + threadIdx.x,
			DTT_SIZE1); // int inc,
 	dctiv_nodiverg( // SC
 			clt_tile + threadIdx.x + 1 * (DTT_SIZE1 * DTT_SIZE),
			DTT_SIZE1); // int inc,
 	dstiv_nodiverg( // CS
 			clt_tile + threadIdx.x + 2 * (DTT_SIZE1 * DTT_SIZE), // &clt_tile[1][0][threadIdx.x], // pointer to start of column
			DTT_SIZE1); // int inc,
 	dstiv_nodiverg( // SS
 			clt_tile + threadIdx.x + 3 * (DTT_SIZE1 * DTT_SIZE), // &clt_tile[1][0][threadIdx.x], // pointer to start of column
			DTT_SIZE1); // int inc,
  	 __syncthreads();// __syncwarp();
}



//
// Uses 16 threads, gets 4*8*8 clt tiles, performs idtt-iv (swapping 1 and 2 quadrants) and then unfolds with window,
// adding to the output 16x16 tile (to use Read-modify-write with 4 passes over the frame. Should be zeroed before the
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
    		clt_tile +  DTT_SIZE1 * (thr012 + 2*DTT_SIZE * thr3), // pointer to start of row for quadrants 0 and 2
			1);
	// perform horizontal dst-iv on quadrants 2 and 3
    dstiv_nodiverg( // all colors
    		clt_tile1 + DTT_SIZE1 * (thr012 + 2*DTT_SIZE * thr3), // pointer to start of row for quadrants 1 and 3
			1);
    __syncthreads();// __syncwarp();
	// perform vertical   dct-iv on quadrants 0 and 2
    dctiv_nodiverg(
    		clt_tile +  thr012 + (DTT_SIZE1 *   DTT_SIZE) * thr3, // pointer to start of row for quadrants 0 and 1
			DTT_SIZE1);
	// perform vertical   dst-iv on quadrants 1 and 3
    dstiv_nodiverg(
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


// Uses 8 threads, gets 4*8*8 clt tiles, performs idtt-iv (swapping 1 and 2 quadrants) and then unfolds to the 16x16
// adding to the output 16x16 tile (to use Read-modify-write with 4 passes over the frame. Should be zeroed before the
// first pass
//__constant__ int imclt_indx9[16] = {0x28,0x31,0x3a,0x43,0x43,0x3a,0x31,0x28,0x1f,0x16,0x0d,0x04,0x04,0x0d,0x16,0x1f};

__device__ void imclt8threads(
		int     do_acc,     // 1 - add to previous value, 0 - overwrite
		float * clt_tile,   //        [4][DTT_SIZE][DTT_SIZE1], // +1 to alternate column ports [4][8][9]
		float * mclt_tile,  //           [2* DTT_SIZE][DTT_SIZE1+ DTT_SIZE], // +1 to alternate column ports[16][17]
		int debug)
{
//	int thr3 =    threadIdx.x >> 3;
//	int column =  threadIdx.x; // modify to use 2*8 threads, if needed.
//	int thr012 =  threadIdx.x & 7;
//	int column4 = threadIdx.x >> 2;
//	int wcolumn = ((thr3 << 3) - thr3) ^ thr012; //0..7,7,..0
	float * clt_tile1 = clt_tile +  (DTT_SIZE1 * DTT_SIZE);
	float * clt_tile2 = clt_tile1 + (DTT_SIZE1 * DTT_SIZE);
	float * clt_tile3 = clt_tile2 + (DTT_SIZE1 * DTT_SIZE);
#ifdef DEBUG7
    if (debug && (threadIdx.x == 0) && (threadIdx.y == 0)){
        printf("\nDTT Tiles before IDTT\n");
        debug_print_clt_scaled(clt_tile, -1,  0xf, 0.25); // only 1 quadrant for R,B and 2 - for G
    }
     __syncthreads();// __syncwarp();
#endif

	// perform horizontal dct-iv on quadrants 0 and 1
    dctiv_nodiverg( // quadrant 0
    		clt_tile +  threadIdx.x,                              // pointer to start of row for quadrant 0
			DTT_SIZE1);
    dctiv_nodiverg( // quadrant 1
    		clt_tile +  threadIdx.x + (1 * DTT_SIZE * DTT_SIZE1), // pointer to start of row for quadrant 1
			DTT_SIZE1);
	// perform horizontal dst-iv on quadrants 2 and 3
    dstiv_nodiverg( // quadrant 2
    		clt_tile +  threadIdx.x + (2 * DTT_SIZE * DTT_SIZE1), // pointer to start of row for quadrant 2
			DTT_SIZE1);
    dstiv_nodiverg( // quadrant 3
    		clt_tile +  threadIdx.x + (3 * DTT_SIZE * DTT_SIZE1), // pointer to start of row for quadrant 3
			DTT_SIZE1);
    __syncthreads();// __syncwarp();
	// perform vertical   dct-iv on quadrants 0 and 2
    dctiv_nodiverg( // quadrant 0
    		clt_tile +  DTT_SIZE1 * threadIdx.x,                              // pointer to start of row for quadrant 0
			1);
    dctiv_nodiverg( // quadrant 2
    		clt_tile +  DTT_SIZE1 * threadIdx.x + (2 * DTT_SIZE * DTT_SIZE1), // pointer to start of row for quadrant 2
			1);
    // perform vertical   dst-iv on quadrants 1 and 3
    dstiv_nodiverg( // quadrant 1
    		clt_tile +  DTT_SIZE1 * threadIdx.x + (1 * DTT_SIZE * DTT_SIZE1), // pointer to start of row for quadrant 1
			1);
    dstiv_nodiverg( // quadrant 3
    		clt_tile +  DTT_SIZE1 * threadIdx.x + (3 * DTT_SIZE * DTT_SIZE1), // pointer to start of row for quadrant 3
			1);
    __syncthreads();// __syncwarp();

#ifdef DEBUG7
    if (debug && (threadIdx.x == 0) && (threadIdx.y == 0)){
    	printf("\nDTT Tiles after IDTT\n");
    	debug_print_clt_scaled(clt_tile, -1,  0xf, 0.25); // only 1 quadrant for R,B and 2 - for G
    }
    __syncthreads();// __syncwarp();
#endif
    // re-using 16-thread code (thr3 was bit 3 of threadIdx.x).
    for (int thr3 = 0; thr3 < 2; thr3++){
    	int thr3m = (thr3 << 3);
    	int column =  threadIdx.x + thr3m; // modify to use 2*8 threads, if needed.
    	int thr012 =  threadIdx.x & 7; // == threadIdx.x
    	int column4 = column >> 2; // (threadIdx.x >> 2) | (thr3 << 1) ; // different !
    	int wcolumn = (thr3m - thr3) ^ thr012; //0..7,7,..0

    	float hw = HWINDOW2[wcolumn];
    	int clt_offset = imclt_indx9[column]; // index in each of the 4 iclt quadrants, accounting for stride=9
    	float * rslt = mclt_tile + column;
#ifdef DEBUG7
        if (debug && (threadIdx.x == 0) && (threadIdx.y == 0)){
    	printf("\nUnrolling: thr3=%d, thr3m=%d, column=%d, thr012=%d, column4=%d, wcolumn=%d, hw=%f, clt_offset=%d\n",
    			thr3, thr3m, column, thr012, column4, wcolumn, hw, clt_offset);
    	debug_print_clt1(clt_tile, -1,  0xf); // only 1 quadrant for R,B and 2 - for G
    }
    __syncthreads();// __syncwarp();
#endif

#pragma unroll
    	for (int i = 0; i < 4; i++){
    		float val = *rslt;
    		// facc
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
    		// val =__fmaf_rd(w,d0,val); // w*d0 + val
    		// *rslt = val;
    		*rslt = do_acc? __fmaf_rd(w,d0,val) : w * d0; // w*d0 + val do_acc - common for all thereads
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
//    		*rslt = __fmaf_rd(w,d0,val); // w*d0 + val
    		*rslt = do_acc? __fmaf_rd(w,d0,val) : w * d0; // w*d0 + val do_acc - common for all thereads

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
    		//*rslt = __fmaf_rd(w,d0,val); // w*d0 + val
    		*rslt = do_acc? __fmaf_rd(w,d0,val) : w * d0; // w*d0 + val do_acc - common for all thereads
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
    		//*rslt = __fmaf_rd(w,d0,val); // w*d0 + val
    		*rslt = do_acc? __fmaf_rd(w,d0,val) : w * d0; // w*d0 + val do_acc - common for all thereads
    		rslt += DTT_SIZE21;
    	}
    }
#ifdef DEBUG7
    __syncthreads();// __syncwarp();
	for (int ccam = 0; ccam < NUM_CAMS; ccam++) {
		if (debug  && (threadIdx.x == 0) && (threadIdx.y == ccam)){
			printf("\nMCLT Tiles after IMCLT, cam=%d\n", threadIdx.y);
			debug_print_mclt(
					mclt_tile, //         [4][DTT_SIZE][DTT_SIZE1], // +1 to alternate column ports)
					-1);
		}
		__syncthreads();// __syncwarp();
	}
    __syncthreads();// __syncwarp();
#endif
}




//#endif

