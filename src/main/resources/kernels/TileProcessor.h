/**
 **
 ** TileProcessor.h
 **
 ** Copyright (C) 2020 Elphel, Inc.
 **
 ** -----------------------------------------------------------------------------**
 **
 **  TileProcessor.h is free software: you can redistribute it and/or modify
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
 * \file TileProcessor.h
 * \brief header file for  the Tile Processor for frequency domain

 */
#pragma once
#ifndef NUM_CAMS
#include "tp_defines.h"
#endif


extern "C" __global__ void convert_direct( // called with a single block, single thread
		//		struct CltExtra ** gpu_kernel_offsets, // [NUM_CAMS], // changed for jcuda to avoid struct parameters
		int                num_cams,           // actual number of cameras
		int                num_colors,         // actual number of colors: 3 for RGB, 1 for LWIR/mono
		float           ** gpu_kernel_offsets, // [NUM_CAMS],
		float           ** gpu_kernels,        // [NUM_CAMS],
		float           ** gpu_images,         // [NUM_CAMS],
		float            * gpu_ftasks,         // flattened tasks, 29 floats for quad EO, 101 floats for LWIR16
//		struct tp_task   * gpu_tasks,
		float           ** gpu_clt,            // [NUM_CAMS][TILES-Y][TILES-X][NUM_COLORS][DTT_SIZE*DTT_SIZE]
		size_t             dstride,            // in floats (pixels)
		int                num_tiles,          // number of tiles in task
		int                lpf_mask,           // apply lpf to colors : bit 0 - red, bit 1 - blue, bit2 - green. Now - always 0 !
		int                woi_width,
		int                woi_height,
		int                kernels_hor,
		int                kernels_vert,
		int *              gpu_active_tiles,      // pointer to the calculated number of non-zero tiles
		int *              pnum_active_tiles,  //  indices to gpu_tasks
		int                tilesx);

extern "C" __global__ void correlate2D(
		int               num_cams,
//		int *             sel_pairs,
		int               sel_pairs0,
		int               sel_pairs1,
		int               sel_pairs2,
		int               sel_pairs3,
		float          ** gpu_clt,            // [NUM_CAMS] ->[TILES-Y][TILES-X][NUM_COLORS][DTT_SIZE*DTT_SIZE]
		int               colors,             // number of colors (3/1)
		float             scale0,             // scale for R
		float             scale1,             // scale for B
		float             scale2,             // scale for G
		float             fat_zero,           // here - absolute
		float            * gpu_ftasks,         // flattened tasks, 29 floats for quad EO, 101 floats for LWIR16
//		struct tp_task  * gpu_tasks,          // array of per-tile tasks (now bits 4..9 - correlation pairs)
		int               num_tiles,          // number of tiles in task
		int               tilesx,             // number of tile rows
		int             * gpu_corr_indices,   // packed tile+pair
		int             * pnum_corr_tiles,    // pointer to a number of correlation tiles to process
		const size_t      corr_stride,        // in floats
		int               corr_radius,        // radius of the output correlation (7 for 15x15)
		float           * gpu_corrs);          // correlation output data

extern "C" __global__ void corr2D_normalize(
		int               num_corr_tiles,     // number of correlation tiles to process
		const size_t      corr_stride_td,     // in floats
		float           * gpu_corrs_td,       // correlation tiles in transform domain
		const size_t      corr_stride,        // in floats
		float           * gpu_corrs,          // correlation output data (either pixel domain or transform domain
		float             fat_zero,           // here - absolute
		int               corr_radius);        // radius of the output correlation (7 for 15x15)

extern "C" __global__ void corr2D_combine(
		int               num_tiles,          // number of tiles to process (each with num_pairs)
		int               num_pairs,          // num pairs per tile (should be the same)
		int               init_output,        // !=0 - reset output tiles to zero before accumulating
		int               pairs_mask,         // selected pairs (0x3 - horizontal, 0xc - vertical, 0xf - quad, 0x30 - cross)
		int             * gpu_corr_indices,   // packed tile+pair
		int             * gpu_combo_indices,  // output if noty null: packed tile+pairs_mask (will point to the first used pair
		const size_t      corr_stride,        // (in floats) stride for the input TD correlations
		float           * gpu_corrs,          // input correlation tiles
		const size_t      corr_stride_combo,  // (in floats) stride for the output TD correlations (same as input)
		float           * gpu_corrs_combo);   // combined correlation output (one per tile)

extern "C" __global__ void textures_nonoverlap(
		int               num_cams,           // number of cameras
		float            * gpu_ftasks,         // flattened tasks, 29 floats for quad EO, 101 floats
		//		struct tp_task  * gpu_tasks,
		int               num_tiles,          // number of tiles in task list
//		int               num_tilesx,         // number of tiles in a row
// declare arrays in device code?
		int             * gpu_texture_indices,// packed tile + bits (now only (1 << 7)
		int             * pnum_texture_tiles,  // returns total number of elements in gpu_texture_indices array
		float          ** gpu_clt,            // [NUM_CAMS] ->[TILES-Y][TILES-X][NUM_COLORS][DTT_SIZE*DTT_SIZE]
		// TODO: use geometry_correction rXY !
		struct gc       * gpu_geometry_correction,
		int               colors,             // number of colors (3/1)
		int               is_lwir,            // do not perform shot correction
		float             params[5],
		float             weights[3],         // scale for R,B,G
		int               dust_remove,        // Do not reduce average weight when only one image differs much from the average
// combining both non-overlap and overlap (each calculated if pointer is not null )
		size_t            texture_stride,     // in floats (now 256*4 = 1024)  // may be 0 if not needed
		float           * gpu_texture_tiles,  // (number of colors +1 + ?)*16*16 rgba texture tiles    // may be 0 if not needed
		int               linescan_order,     // 0 low-res tiles have tghe same order, as gpu_texture_indices, 1 - in linescan order
		float           * gpu_diff_rgb_combo, //); // diff[NUM_CAMS], R[NUM_CAMS], B[NUM_CAMS],G[NUM_CAMS] // may be 0 if not needed
		int               num_tilesx);

extern "C"
__global__ void imclt_rbg_all(
		int                num_cams,
		float           ** gpu_clt,            // [NUM_CAMS][TILES-Y][TILES-X][NUM_COLORS][DTT_SIZE*DTT_SIZE]
		float           ** gpu_corr_images,    // [NUM_CAMS][WIDTH, 3 * HEIGHT]
		int                apply_lpf,
		int                colors,
		int                woi_twidth,
		int                woi_theight,
		const size_t       dstride);            // in floats (pixels)

extern "C" __global__ void erase8x8(
		float           * gpu_top_left,
		const size_t      dstride);

extern "C" __global__ void imclt_rbg(
		float           * gpu_clt,            // [TILES-Y][TILES-X][NUM_COLORS][DTT_SIZE*DTT_SIZE]
		float           * gpu_rbg,            // WIDTH, 3 * HEIGHT
		int               apply_lpf,
		int               mono,               // defines lpf filter
		int               color,              // defines location of clt data
		int               v_offset,
		int               h_offset,
		int               woi_twidth,
		int               woi_theight,
		const size_t      dstride);            // in floats (pixels)

extern "C" __global__ void generate_RBGA(
		int                num_cams,           // number of cameras used
		// Parameters to generate texture tasks
		float            * gpu_ftasks,         // flattened tasks, 29 floats for quad EO, 101 floats for LWIR16
//		struct tp_task   * gpu_tasks,
		int                num_tiles,          // number of tiles in task list
		// declare arrays in device code?
		int              * gpu_texture_indices,// packed tile + bits (now only (1 << 7)
		int              * num_texture_tiles,  // number of texture tiles to process  (8 separate elements for accumulation)
		int              * woi,                // x,y,width,height of the woi
		int                width,  // <= TILES-X, use for faster processing of LWIR images (should be actual + 1)
		int                height, // <= TILES-Y, use for faster processing of LWIR images
		// Parameters for the texture generation
		float          ** gpu_clt,            // [NUM_CAMS] ->[TILES-Y][TILES-X][NUM_COLORS][DTT_SIZE*DTT_SIZE]
		// TODO: use geometry_correction rXY !
		struct gc       * gpu_geometry_correction,
		int               colors,             // number of colors (3/1)
		int               is_lwir,            // do not perform shot correction
		float             params[5],          // mitigating CUDA_ERROR_INVALID_PTX
		float             weights[3],         // scale for R,B,G
		int               dust_remove,        // Do not reduce average weight when only one image differs much from the average
		int               keep_weights,       // return channel weights after A in RGBA (was removed)
		const size_t      texture_rbga_stride,     // in floats
		float           * gpu_texture_tiles);  // (number of colors +1 + ?)*16*16 rgba texture tiles

