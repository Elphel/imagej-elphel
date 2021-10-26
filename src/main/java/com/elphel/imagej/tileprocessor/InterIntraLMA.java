package com.elphel.imagej.tileprocessor;

import java.util.Arrays;

/**
 **
 ** InterIntraLMA - Comparing contrast for different interscene and intrascene
 ** connfifurations
 **
 ** Copyright (C) 2021 Elphel, Inc.
 **
 ** -----------------------------------------------------------------------------**
 **
 **  InterIntraLMA.java is free software: you can redistribute it and/or modify
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
 ** -----------------------------------------------------------------------------**
 **
 */


public class InterIntraLMA {
	/**
	 * Estimate noise threshold for each combination of inter/intra and number of
	 * sensors.
	 * @param noise_file synthetic noise level per file (either rnd or fpn)
	 * @param senosor_mode_file sensor mode (0 - 16, 1 - 8, 2 - 4, 3 - 2 sensors) per file
	 * @param inter_file true - interscene , false - intrascene only
	 * @param good_file_tile true if tile of this file is considered "good"
	 * @param min_modes minimal number of modes the tile should be defined for (undefine if less than)
	 * @return per tile per mode average between highest noise of the good tile and lowest noise
	 * of the bad tile. Double.NaN if noise threshold can not be determined. null if the tile is
	 * undefined for all modes   
	 */
	public static double [][] getNoiseThreshold(
			double []    noise_file, //  = new double [noise_files.length];
			int []       sensor_mode_file,
			boolean []   inter_file,
			boolean [][] good_file_tile,
			int          min_modes)
	{
		int num_sensor_modes = 0;
		int num_tiles = good_file_tile[0].length;
		for (int i = 0; i < sensor_mode_file.length; i++) {
			if (sensor_mode_file[i] > num_sensor_modes) {
				num_sensor_modes = sensor_mode_file[i];
			}
		}
		num_sensor_modes ++;
		int num_modes = 2 * num_sensor_modes;
		double [][] rslt = new double [num_tiles][]; // number of tiles  
		double [][][] noise_interval = new double[num_modes][num_tiles][2]; // [modes][tiles]
		for (int i = 0; i < num_modes; i++) {
			for (int j= 0; j < num_tiles; j++) {
				noise_interval[i][j][0] =Double.NaN;
				noise_interval[i][j][1] =Double.NaN;
			}
		}
		for (int nf = 0; nf < noise_file.length; nf++) {
			double noise = noise_file[nf];
			int mode = sensor_mode_file[nf] + (inter_file[nf] ? 0: num_sensor_modes); 
			for (int ntile = 0; ntile < num_tiles; ntile++) {
				if (good_file_tile[nf][ntile]) { // good tile
					if (!(noise <= noise_interval[mode][ntile][0])){ // including Double.isNaN(noise_interval[mode][ntile][0]
						noise_interval[mode][ntile][0] = noise;
					}
				} else { // bad tile
					if (!(noise >= noise_interval[mode][ntile][1])){ // including Double.isNaN(noise_interval[mode][ntile][1]
						noise_interval[mode][ntile][1] = noise;
					}
				}
			}
		}
		for (int ntile = 0; ntile < num_tiles; ntile++){ 
			int num_defined = 0;
			for (int mode = 0; mode < num_modes; mode++) {
				if (!Double.isNaN(noise_interval[mode][ntile][0]) && !Double.isNaN(noise_interval[mode][ntile][0])) {
					num_defined++;
				}
			}
			if (num_defined >= min_modes) {
				rslt[ntile] = new double [num_modes];
				for (int mode = 0; mode < num_modes; mode++) {
					if (!Double.isNaN(noise_interval[mode][ntile][0]) && !Double.isNaN(noise_interval[mode][ntile][0])) {
						rslt[ntile][mode] = 0.5 * (noise_interval[mode][ntile][0] + noise_interval[mode][ntile][1]);
					} else {
						rslt[ntile][mode] = Double.NaN;
					}
				}
			}
		}
		return rslt;
	}

}
