/**
 **
 ** geometry_correction.h
 **
 ** Copyright (C) 2020 Elphel, Inc.
 **
 ** -----------------------------------------------------------------------------**
 **
 **  geometry_correction.h is free software: you can redistribute it and/or modify
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
* \file geometry_correction.h
* \brief header file for geometry correction - per-tile/per camera calculation of the tile offset

*/
#pragma once
#ifndef NUM_CAMS
#include "tp_defines.h"
#endif

#define NVRTC_BUG 1
#ifndef M_PI
#define M_PI  3.14159265358979323846 /* pi */
#endif
#ifndef offsetof
#define offsetof(st, m) \
    ((size_t)&(((st *)0)->m))
//#define offsetof(TYPE, MEMBER) __builtin_offsetof (TYPE, MEMBER)
#endif


#define SCENE_UNITS_SCALE  0.001 // meters from mm
#define MIN_DISPARITY      0.01  // minimal disparity to try to convert to world coordinates
struct tp_task {
	int   task;
	union {
		int      txy;
		unsigned short sxy[2];
	};
	float xy[NUM_CAMS][2];
	float target_disparity;
	float disp_dist[NUM_CAMS][4]; // calculated with getPortsCoordinates()
};

struct corr_vector{
	float tilt    [NUM_CAMS-1]; // 0..2
	float azimuth [NUM_CAMS-1]; // 3..5
	float roll    [NUM_CAMS];   // 6..9
	float zoom    [NUM_CAMS-1]; // 10..12
	// for ERS correction:
	float imu_rot [3]; // d_tilt/dt (rad/s), d_az/dt, d_roll/dt 13..15
	float imu_move[3]; // dx/dt, dy/dt, dz/dt 16..19
};
#ifdef NVRTC_BUG
struct trot_deriv{
	float rots    [NUM_CAMS][3][3];
	float d_daz   [NUM_CAMS][3][3];
	float d_tilt  [NUM_CAMS][3][3];
	float d_roll  [NUM_CAMS][3][3];
	float d_zoom  [NUM_CAMS][3][3];
};
#else
union trot_deriv{
	struct {
		float rots    [NUM_CAMS][3][3];
		float d_daz   [NUM_CAMS][3][3];
		float d_tilt  [NUM_CAMS][3][3];
		float d_roll  [NUM_CAMS][3][3];
		float d_zoom  [NUM_CAMS][3][3];
	};
	float matrices [5][NUM_CAMS][3][3];
};
#endif

struct gc {
	float pixelCorrectionWidth; //  =2592;   // virtual camera center is at (pixelCorrectionWidth/2, pixelCorrectionHeight/2)
	float pixelCorrectionHeight; // =1936;
	float line_time;        // duration of one scan line readout (for ERS)
	float focalLength;      // =FOCAL_LENGTH;
	float pixelSize;        // =  PIXEL_SIZE; //um
	float distortionRadius; // =  DISTORTION_RADIUS; // mm - half width of the sensor
#ifndef	NVRTC_BUG
	union {
		struct {
#endif
			float distortionC;      // r^2
			float distortionB;      // r^3
			float distortionA;      // r^4 (normalized to focal length or to sensor half width?)
			float distortionA5;     //r^5 (normalized to focal length or to sensor half width?)
			float distortionA6;     //r^6 (normalized to focal length or to sensor half width?)
			float distortionA7;     //r^7 (normalized to focal length or to sensor half width?)
			float distortionA8;     //r^8 (normalized to focal length or to sensor half width?)
#ifndef	NVRTC_BUG
//		};
//		float rad_coeff [7];
//	};
#endif
	// parameters, common for all sensors
	float    elevation;     // degrees, up - positive;
	float    heading;       // degrees, CW (from top) - positive

	float forward    [NUM_CAMS];
	float right      [NUM_CAMS];
	float height     [NUM_CAMS];
	float roll       [NUM_CAMS];    // degrees, CW (to target) - positive
	float pXY0       [NUM_CAMS][2];
	float common_right;    // mm right, camera center
	float common_forward;  // mm forward (to target), camera center
	float common_height;   // mm up, camera center
	float common_roll;     // degrees CW (to target) camera as a whole
//	float [][] XYZ_he;     // all cameras coordinates transformed to eliminate heading and elevation (rolls preserved)
//	float [][] XYZ_her = null; // XYZ of the lenses in a corrected CCS (adjusted for to elevation, heading,  common_roll)
	float rXY        [NUM_CAMS][2]; // XY pairs of the in a normal plane, relative to disparityRadius
//	float [][] rXY_ideal = {{-0.5, -0.5}, {0.5,-0.5}, {-0.5, 0.5}, {0.5,0.5}};
// only used for the multi-quad systems
	float cameraRadius; // =0; // average distance from the "mass center" of the sensors to the sensors
	float disparityRadius; // =150.0; // distance between cameras to normalize disparity units to. sqrt(2)*disparityRadius for quad
};
#define RAD_COEFF_LEN 7
extern "C" __global__ void get_tiles_offsets(
		struct tp_task     * gpu_tasks,
		int                  num_tiles,          // number of tiles in task
		struct gc          * gpu_geometry_correction,
		struct corr_vector * gpu_correction_vector,
		float *              gpu_rByRDist, // length should match RBYRDIST_LEN
		trot_deriv   * gpu_rot_deriv);

#if 0
// uses 3 threadIdx.x, 3 - threadIdx.y, 4 - threadIdx.z
extern "C" __global__ void calc_rot_matrices(
		struct corr_vector * gpu_correction_vector);
#endif
// uses NUM_CAMS blocks, (3,3,3) threads
extern "C" __global__ void calc_rot_deriv(
		struct corr_vector * gpu_correction_vector,
		trot_deriv   * gpu_rot_deriv);


