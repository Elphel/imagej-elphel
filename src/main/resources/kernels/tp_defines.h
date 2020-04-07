/**
 **
 ** tp_defines.h
 **
 ** Copyright (C) 2020 Elphel, Inc.
 **
 ** -----------------------------------------------------------------------------**
 **
 **  tp_defines.h is free software: you can redistribute it and/or modify
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
* \file tp_defines.h
* \brief Defines for running in C++ environment, replaced when called from Java

*/
// Avoiding includes in jcuda, all source files will be merged
#pragma once
#ifndef JCUDA
#define THREADSX         (DTT_SIZE)
#define NUM_CAMS                  4
#define NUM_PAIRS                 6
#define NUM_COLORS                3
#define IMG_WIDTH              2592
#define IMG_HEIGHT             1936
#define KERNELS_HOR             164
#define KERNELS_VERT            123
#define KERNELS_LSTEP             4
#define THREADS_PER_TILE          8
#define TILES_PER_BLOCK           4
#define CORR_THREADS_PER_TILE     8
#define CORR_TILES_PER_BLOCK      4
#define TEXTURE_THREADS_PER_TILE  8
#define TEXTURE_TILES_PER_BLOCK   1
#define IMCLT_THREADS_PER_TILE   16
#define IMCLT_TILES_PER_BLOCK     4
#define CORR_NTILE_SHIFT          8 // higher bits - number of a pair, other bits tile number
#define CORR_PAIRS_MASK        0x3f// lower bits used to address correlation pair for the selected tile
#define CORR_TEXTURE_BIT          7 // bit 7 used to request texture for the tile
#define TASK_CORR_BITS            4
#define TASK_TEXTURE_N_BIT        0 // Texture with North neighbor
#define TASK_TEXTURE_E_BIT        1 // Texture with East  neighbor
#define TASK_TEXTURE_S_BIT        2 // Texture with South neighbor
#define TASK_TEXTURE_W_BIT        3 // Texture with West  neighbor
#define TASK_TEXTURE_BIT          3 // bit to request texture calculation int task field of struct tp_task
#define LIST_TEXTURE_BIT          7 // bit to request texture calculation
#define CORR_OUT_RAD              4
#define FAT_ZERO_WEIGHT           0.0001 // add to port weights to avoid nan

#define THREADS_DYNAMIC_BITS      5 // treads in block for CDP creation of the texture list

//#undef HAS_PRINTF
#define HAS_PRINTF
//7
//#define DEBUG1 1
//#define DEBUG2 1
//#define DEBUG3 1
//#define DEBUG4 1
//#define DEBUG5 1
//#define DEBUG6 1
/*
#define DEBUG7 1
#define DEBUG8 1
#define DEBUG9 1
*/
#define DEBUG10 1
#define DEBUG11 1
#define DEBUG12 1
//#define USE_textures_gen
#define DEBUG_OOB1 1
#endif //#ifndef JCUDA

