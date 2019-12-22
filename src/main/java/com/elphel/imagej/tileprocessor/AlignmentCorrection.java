package com.elphel.imagej.tileprocessor;
/**
 **
 ** AlignmentCorrection - try to apply minor adjustments to the misaligned camera
 **
 ** Copyright (C) 2017 Elphel, Inc.
 **
 ** -----------------------------------------------------------------------------**
 **
 **  AlignmentCorrection.java is free software: you can redistribute it and/or modify
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
import java.util.ArrayList;
import java.util.Arrays;

import com.elphel.imagej.cameras.CLTParameters;
import com.elphel.imagej.common.DoubleGaussianBlur;
import com.elphel.imagej.common.PolynomialApproximation;
import com.elphel.imagej.common.ShowDoubleFloatArrays;

//import GeometryCorrection.CorrVector;
import Jama.Matrix;
import ij.ImagePlus;
import ij.ImageStack;

public class AlignmentCorrection {
	static final int NUM_SLICES =      10; // disp, strength, dx0, dy0, dx1, dy1, dx2, dy2, dx3, dy3)
	static final int NUM_ALL_SLICES =  14; // disp, strength, dx0, dy0, str0, dx1, dy1, str1, dx2, dy2, str2, dx3, dy3 str3,)
	static final int NUM_SENSORS = 4;
	static final int [] INDICES_14_10 = {0,1,2,3,5,6,8,9,11,12};
	static final int [] INDICES_10_DISP = {2,4,7,9}; // which indices need to add disparity to restore full offsets
	static final int [] INDICES_14_WEIGHTS = {4,7,10,13}; // now pair weights can be zeros where common weight is not (and respective values are NaNs)
	static final int    INDEX_14_WEIGHT = 1;
	static final int    INDEX_10_WEIGHT = 1;
	static final int    INDEX_10_DISPARITY = 0;
	static final int    INDEX_14_DISPARITY = 0;
	static final double DISP_SCALE = 2.0;

	QuadCLT qc;

	public class Sample{
		public int series;
		public int tile;
		public double weight;
		Sample(int series, int tile, double weight) {
			this.series = series;
			this.tile =  tile;
			this.weight = weight;

		}
	}

	public class Mismatch{
		public boolean   use_disparity; // adjust dx0+dx1+dy0+dy1 == 0
		public double [] pXY; // tile center x,y
		public double    disparity_task;
		public double    disparity_meas;
		public double    strength;
		public double [][] offsets; // measured {{dx0,dy0}, {dx1,dy1},{dx2,dy2},{dx3,dy3}};

		public Mismatch()
		{
			pXY = new double[2];
			disparity_task = 0.0;
			disparity_meas = 0.0;
			strength =  0.0;
			offsets = new double[NUM_SENSORS][2];
		}
		public Mismatch (
				boolean   use_disparity, // adjust dx0+dx1+dy0+dy1 == 0
				double [] pXY, // tile center x,y
				double    disparity_task,
				double    disparity_meas,
				double    strength,
				double [][] offsets )
		{
			this.use_disparity = use_disparity;
			this.pXY = pXY;
			//FIXME:  2 next fields are not used
			this.disparity_task = disparity_task; // currently wrong, series is 0/1, not measurement number
			this.disparity_meas = disparity_meas;
			this.strength =  (Double.isNaN(strength) || (offsets==null) || Double.isNaN(pXY[0]) || Double.isNaN(pXY[1]))?0.0:strength;
			this.offsets = offsets;
			if (strength != 0.0) {
				boolean bug = (offsets == null);
				if (!bug) {
					for (int i = 0; i < offsets.length; i++) {
						bug |= (offsets[i] == null) || Double.isNaN(offsets[i][0])  || Double.isNaN(offsets[i][1]);
						if (bug) break;
					}
				}
				if (bug) {
					System.out.println("***** BUG: NaN offsets, but strength != 0: x= "+pXY[0]+", y = "+pXY[1]+" ******");
					System.out.println("***** BUG: NaN offsets, but strength != 0: x= "+pXY[0]+", y = "+pXY[1]+" ******");
					this.strength = 0;
				}
			}
		}

		public boolean usesDisparity()
		{
			return use_disparity;
		}
		public double [] getPXY(){
			return pXY;
		}

		public double getDisparityTask() {
			return disparity_task;
		}
		public double getDisparityMeas() {
			return disparity_meas;
		}
		public double getStrength() {
			return strength;
		}

		public double [] getOffsets(){
			double [] offsets1d = new double [offsets.length * offsets[0].length];
			for (int i = 0; i < offsets.length; i++){
				for (int j = 0; j < offsets[i].length; j++){
					offsets1d[i * offsets[i].length + j] = offsets[i][j];
				}
			}
			return offsets1d;
		}

		/**
		 * Calculate error vector of 8 components, the last one will be optionally masked by weight - used only when measured
		 * disparity should be 0 (as at infinity)
		 * @return 8 components
		 */
		public double [] getY()
		{
			double [] y = {
					offsets[0][1], //dy0;
					offsets[1][1], //dy1;
					offsets[2][0], //dx2;
					offsets[3][0], //dx3;
					0.5*(offsets[1][0] - offsets[0][0]), //(dx1 - dx0)/2; 0.5/magic is already applied
					0.5*(offsets[3][1] - offsets[2][1]), //(dy3 - dy2)/2;
					0.25* (offsets[0][0] + offsets[1][0] - offsets[2][1] - offsets[3][1]), //(dx0 + dx1 -dy2 - dy3)/4;
					use_disparity ? (DISP_SCALE* 0.125* (offsets[0][0] + offsets[1][0] + offsets[2][1] + offsets[3][1])) : 0.0 // only when disparity is known to be 0
			};
			return y;
		}
/*
		public double [] getY()
		{
			double [] y = {
					offsets[0][1], //dy0;
					offsets[1][1], //dy1;
					offsets[2][0], //dx2;
					offsets[3][0], //dx3;
					(offsets[1][0] - offsets[0][0]), //(dx1 - dx0)/2; 0.5/magic is already applied
					(offsets[3][1] - offsets[2][1]), //(dy3 - dy2)/2;
					0.5* (offsets[0][0] + offsets[1][0] - offsets[2][1] - offsets[3][1]), //(dx0 + dx1 -dy2 - dy3)/4;
					use_disparity ? (0.25* (offsets[0][0] + offsets[1][0] + offsets[2][1] + offsets[3][1])) : 0.0 // only when disparity is known to be 0
			};
			return y;
		}

 */
		public void copyToY(
				double [] y,
				int n_sample)
		{
			double [] sub_y = getY();
			System.arraycopy(sub_y, 0, y, n_sample * (2 * NUM_SENSORS), (2 * NUM_SENSORS));
		}

		public void copyToW(
				double [] w,
				int n_sample)
		{
			for (int i = 0; i < 2 * NUM_SENSORS; i++){
//				w[n_sample * (2 * NUM_SENSORS) + i] = (!Double.isNaN(strength) && (use_disparity || (i < 7))) ? strength : 0.0;
				// TODO: find out about (i < 7) - it zeroes out y[3]?
				w[n_sample * (2 * NUM_SENSORS) + i] =(use_disparity || (i < 7)) ? strength : 0.0; // Last element corresponds to disparity,
				// set weight = 0 for samples that do not have it.
			}
		}

// measurement vectors and correction vectors - seem to mismatch for disparity (sym0)
/*
		public double [][] get_dMismatch_dXY()
		{
			double [][] dMismatch_dXY = { // extra 0.5 is because differences dxi, dyi are already *= 0.5/magic
					// FIXME: 0.5/magic is removed, but what is measured is (x1-x0)/2m ...
					//x0       y0      x1       y1      x2      y2      x3      y3
					{ 0.0 ,   -0.5,    0.0 ,   0.5 ,   0.0 ,   0.0 ,   0.0 ,   0.0   }, // mv0 = dy0 = y1 - y0
					{ 0.0 ,    0.0 ,   0.0 ,   0.0 ,   0.0 ,  -0.5 ,   0.0 ,   0.5   }, // mv1 = dy1 = y3 - y2
					{-0.5 ,    0.0 ,   0.0 ,   0.0 ,   0.5 ,   0.0 ,   0.0 ,   0.0   }, // mv2 = dx2 = x2 - x0
					{ 0.0 ,    0.0 ,  -0.5 ,   0.0 ,   0.0 ,   0.0 ,   0.5 ,   0.0   }, // mv3 = dx3 = x3 - x1
					{-0.25,    0.0 ,   0.25,   0.0 ,  -0.25,   0.0 ,   0.25,   0.0   }, // mv4 = (dx1 - dx0)/2 = (x3 - x2 - x0 + x1) / 2
					{ 0.0 ,   -0.25,   0.0 ,  -0.25,   0.0 ,   0.25,   0.0 ,   0.25  }, // mv5 = (dy3 - dy2)/2 = (y3 - y1 - y0 + y2) / 2
					{-0.125,   0.125,  0.125,  0.125, -0.125, -0.125,  0.125, -0.125 }, // mv6 = (dx0 + dx1 -dy2 - dy3)/4 = (x1 - x0 + x3 - x2 - y2 + y0 - y3 + y1)/4
					{-0.0625, -0.0625, 0.0625,-0.0625,-0.0625, 0.0625, 0.0625, 0.0625}};// mv7 = (dx0 + dx1 +dy2 + dy3)/8=  (x1 - x0 + x3 - x2 + y2 - y0 + y3 - y1)/8
			return dMismatch_dXY;
		}
 */

		public double [][] get_dMismatch_dXY()
		{
			double [][] dMismatch_dXY = { // extra 0.5 is because differences dxi, dyi are already *= 0.5/magic
					// FIXME: 0.5/magic is removed, but what is measured is (x1-x0)/2m ...
					//x0       y0      x1       y1      x2      y2      x3      y3
					{ 0.0 ,   -0.5,    0.0 ,   0.5 ,   0.0 ,   0.0 ,   0.0 ,   0.0   }, // mv0 = dy0 = y1 - y0
					{ 0.0 ,    0.0 ,   0.0 ,   0.0 ,   0.0 ,  -0.5 ,   0.0 ,   0.5   }, // mv1 = dy1 = y3 - y2
					{-0.5 ,    0.0 ,   0.0 ,   0.0 ,   0.5 ,   0.0 ,   0.0 ,   0.0   }, // mv2 = dx2 = x2 - x0
					{ 0.0 ,    0.0 ,  -0.5 ,   0.0 ,   0.0 ,   0.0 ,   0.5 ,   0.0   }, // mv3 = dx3 = x3 - x1
					{ 0.25,    0.0 ,  -0.25,   0.0 ,  -0.25,   0.0 ,   0.25,   0.0   }, // mv4 = (dx1 - dx0)/2 = (x3 - x2 + x0 - x1) / 2
					{ 0.0 ,    0.25,   0.0 ,  -0.25,   0.0 ,  -0.25,   0.0 ,   0.25  }, // mv5 = (dy3 - dy2)/2 = (y3 - y1 + y0 - y2) / 2
					{-0.125,   0.125,  0.125,  0.125, -0.125, -0.125,  0.125, -0.125 }, // mv6 = (dx0 + dx1 -dy2 - dy3)/4 = (x1 - x0 + x3 - x2 - y2 + y0 - y3 + y1)/4
					{-0.0625*DISP_SCALE, -0.0625*DISP_SCALE, 0.0625*DISP_SCALE,-0.0625*DISP_SCALE,-0.0625*DISP_SCALE, 0.0625*DISP_SCALE, 0.0625*DISP_SCALE, 0.0625*DISP_SCALE}};// mv7 = (dx0 + dx1 +dy2 + dy3)/8=  (x1 - x0 + x3 - x2 + y2 - y0 + y3 - y1)/8
			return dMismatch_dXY;
		}



		/**
		 * Convert transposed jacobian from {d_dx0,d_dy0, ...,d_dy3} to d_mvi (measurement vectors),
		 * where sum of measurement vectors squared is minimized. Same matrix multiplications
		 * is applied to each group of 8 columns. last column in each group is only non-zero if
		 * disparity is known to be 0;
		 * @param jt transposed Jacobian of 13/10/9 rows and 8*n columns
		 * @return converted transposed Jacobian of the same dimensions
		 */
		double [][] convertJt_mv(
				double [][] jt)
		{
			double [][] dMismatch_dXY = get_dMismatch_dXY();

			double [][] jt_conv = new double [jt.length][jt[0].length/dMismatch_dXY[0].length*dMismatch_dXY.length]; // now dMismatch_dXY is square
			// multiplying by transposed dMismatch_dXY
			for (int g = 0; g < jt[0].length/dMismatch_dXY[0].length; g++) {
				int indx_in =  dMismatch_dXY[0].length * g;
				int indx_out = dMismatch_dXY.length * g;
				for (int np = 0; np < jt.length; np++){
					for (int oindx = 0; oindx < dMismatch_dXY.length; oindx++){
						for (int k = 0; k <  dMismatch_dXY[0].length; k++) {
							jt_conv[np][indx_out + oindx] += jt[np][indx_in + k] * dMismatch_dXY[oindx][k];
						}
					}
				}
			}
			return jt_conv;
		}

	}

	//System.arraycopy(dpixels, (tileY*width+tileX)*dct_size + i*width, tile_in, i*n2, n2);
	AlignmentCorrection (QuadCLT qc){
		this.qc = qc;
	}

	public double [][][] infinityCorrection(
			final boolean    use_poly,
			final double     min_strength_in,
			final double     max_diff,
			final int        max_iterations,
			final double     max_coeff_diff,
			final double     far_pull, //  = 0.2; // 1; //  0.5;
			final double     strength_pow,
			final int        smplSide, //        = 2;      // Sample size (side of a square)
			final int        smplNum, //         = 3;      // Number after removing worst (should be >1)
			final double     smplRms, //         = 0.1;    // Maximal RMS of the remaining tiles in a sample
			// histogram parameters
			final int        hist_smpl_side, // 8 x8 masked, 16x16 sampled
			final double     hist_disp_min,
			final double     hist_disp_step,
			final int        hist_num_bins,
			final double     hist_sigma,
			final double     hist_max_diff,
			final int        hist_min_samples,
			final boolean    hist_norm_center, // if there are more tiles that fit than min_samples, replace with

			CLTParameters           clt_parameters,
			double [][]      disp_strength_in,
			int              tilesX,
			double           magic_coeff, // still not understood coefficient that reduces reported disparity value.  Seems to be around 8.5
			int              debugLevel)
	{
		double [][] disp_strength;
		double min_strength;
		if (smplSide > 1){
			disp_strength = filterDisparityStrength (
					disp_strength_in,
					min_strength_in,    // final double     strength_floor,
					strength_pow,
					smplSide, //        = 2;      // Sample size (side of a square)
					smplNum, //         = 3;      // Number after removing worst (should be >1)
					smplRms, //         = 0.1;    // Maximal RMS of the remaining tiles in a sample
					tilesX);
			min_strength = 0; // all > 0
			if (debugLevel > 0){
				String [] titles= {"disp","strength","dx0","dy0","dx1","dy1","dx2","dy2","dx3","dy3"};
				double [][] dbg_img = disp_strength.clone();
				for (int n = 0; n < disp_strength.length; n++){
					dbg_img[n] = disp_strength[n].clone();
				}
				for (int n = 0; n < dbg_img.length; n++) if (n != 1){
					for (int i = 0; i < dbg_img[n].length; i++) {
						if (dbg_img[1][i] == 0.0){
							dbg_img[n][i] = Double.NaN;
						}
					}
				}
				(new ShowDoubleFloatArrays()).showArrays(dbg_img, tilesX, disp_strength[0].length/tilesX, true, "filtered_ds", titles);

			}

		} else {
			disp_strength = disp_strength_in;
			min_strength = min_strength_in;
			if (debugLevel > 0){
				String [] titles= {"disp","strength","dx0","dy0","dx1","dy1","dx2","dy2","dx3","dy3"};

				double [][] dbg_img = disp_strength.clone();
				for (int n = 0; n < disp_strength.length; n++){
					dbg_img[n] = disp_strength[n].clone();
				}
				for (int n = 0; n < dbg_img.length; n++) if (n != 1){
					for (int i = 0; i < dbg_img[n].length; i++) {
						if (dbg_img[1][i] == 0.0){
							dbg_img[n][i] = Double.NaN;
						}
					}
				}
				(new ShowDoubleFloatArrays()).showArrays(dbg_img, tilesX, disp_strength[0].length/tilesX, true, "inp_ds" , titles);

			}

		}
		if (hist_smpl_side > 0) { // 0 to bypass histogram filtering
			disp_strength = filterHistogramFar (
					disp_strength,    // final double[][] disp_strength_in,
					hist_smpl_side,   // final int        smpl_side, // 8 x8 masked, 16x16 sampled
					hist_disp_min,    // final double     disp_min,
					hist_disp_step,   // final double     disp_step,
					hist_num_bins,    // final int        num_bins,
					hist_sigma,       // final double     sigma,
					hist_max_diff,    // final double     max_diff,
					hist_min_samples, // final int        min_samples,
					hist_norm_center, // final boolean    norm_center, // if there are more tiles that fit than minsamples, replace with
					tilesX);          // final int        tilesX)
			if (debugLevel > 0){
				String [] titles= {"disp","strength","dx0","dy0","dx1","dy1","dx2","dy2","dx3","dy3"};
				double [][] dbg_img = disp_strength.clone();
				for (int n = 0; n < disp_strength.length; n++){
					dbg_img[n] = disp_strength[n].clone();
				}
//				for (int n = 0; n < dbg_img.length; n+=2){
//					for (int i = 0; i < dbg_img[n].length; i++) {
//						if (dbg_img[n+1][i] == 0.0){
//							dbg_img[n][i] = Double.NaN;
//						}
//					}
//				}
				for (int n = 0; n < dbg_img.length; n++) if (n != 1){
					for (int i = 0; i < dbg_img[n].length; i++) {
						if (dbg_img[1][i] == 0.0){
							dbg_img[n][i] = Double.NaN;
						}
					}
				}

				(new ShowDoubleFloatArrays()).showArrays(dbg_img, tilesX, disp_strength[0].length/tilesX, true, "hist_filt_ds", titles);

			}
		}

		ArrayList<Sample> samples_list = selectInfinityTiles(
				clt_parameters.fcorr_radius,       // 	final double fcorr_radius,
				clt_parameters.fcorr_inf_vert,// final boolean use_vertical,
				min_strength,
				max_diff,
				max_iterations,
				max_coeff_diff,
				far_pull, //  = 0.2; // 1; //  0.5;
				clt_parameters,
				disp_strength,
				tilesX,
				magic_coeff, // still not understood coefficient that reduces reported disparity value.  Seems to be around 8.5
				debugLevel);

		double [][][] mismatch_corr_coefficients = null;
		ArrayList<Mismatch> mismatch_list = use_poly? null : (new ArrayList<Mismatch>());

		mismatch_corr_coefficients = infinityMismatchCorrection(
				clt_parameters.disp_scan_start, // final double  disp_scan_start,
				clt_parameters.disp_scan_step,  // final double  disp_scan_step,
				use_poly,                       // final boolean use_poly,
				clt_parameters.fcorr_inf_quad,  // final boolean use_quadratic,
				clt_parameters.fcorr_inf_vert,  // final boolean use_vertical,
				// now disparity is already restored
				false, //clt_parameters.ly_inf_en,       // final boolean use_disparity, // for infinity - if true, restores differences in the direction of disparity that was subtracted during measurement)
				// For ly_inf_en need to make sure that programmed disparity was 0.0, so
				clt_parameters.ly_inf_disp,     //final boolean allow_dispatity,
				clt_parameters,
				disp_strength,
				samples_list,
				tilesX,
				magic_coeff, // still not understood coefficient that reduces reported disparity value.  Seems to be around 8.5
				mismatch_list,                  // ArrayList<Mismatch> mismatch_list,
				debugLevel);
		if (debugLevel > -1){
			System.out.println("infinityCorrection(): coefficient increments from infinityMismatchCorrection");
			if (mismatch_corr_coefficients == null) { // non-null only for poly !
				System.out.println("imismatch_corr_coefficients == null (non-null for polynomial correction only)");
				return mismatch_corr_coefficients;
			}
			show_fine_corr(
					mismatch_corr_coefficients, // double [][][] corr,
					"");// String prefix)
		}
		return mismatch_corr_coefficients;
	}

	/**
	 * Discard correction data outside of the center image area
	 * @param fcorr_radius fraction of the image to use (1.0 - 100%)
	 * @param tilesX width in tiles
	 * @param tilesY height in tiles
	 * @return boolean array in linescan order
	 */
	public boolean[] getCenterMask(
			double     fcorr_radius,
			int tilesX,
			int tilesY)
	{
		boolean [] mask = new boolean [tilesX * tilesY];

		int y0 = (int)  (0.5 * tilesY*(1.0 - fcorr_radius));
		int y1 = (int)  (0.5 * tilesY*(1.0 + fcorr_radius));
		int x0 = (int)  (0.5 * tilesX*(1.0 - fcorr_radius));
		int x1 = (int)  (0.5 * tilesX*(1.0 + fcorr_radius));
		if (y0 < 0)      y0 = 0;
		if (y1 > tilesY) y1 = tilesY;
		if (x0 < 0)      x0 = 0;
		if (x1 > tilesX) x1 = tilesX;
		for (int ty = y0; ty < y1; ty++){
			for (int tx = x0; tx < x1; tx++){
				mask[tx + tilesX * ty] = true;
			}
		}
		return mask;
	}


	/**
	 * Select infinity tiles from a single or multiple image sets
	 * Next parameters are made separate to be able to modify them between different runs keeping clt_parameters
	 * @param fcorr_radius do not use peripheral tiles
	 * @param min_strength minimal correlation strength to use tile
	 * @param max_diff maximal disparity difference between tiles and previous approximation to use tiles
	 * @param max_iterations maximal number of iterations to find disparity surface
	 * @param max_coeff_diff coefficients maximal change to continue iterations
	 * @param far_pull weight coefficient for the tile outside of the range, but farther
	 * @param clt_parameters CLT parameters
	 * @param disp_strength array of a single or multiple disparity/strength pairs (0,2, .. - disparity,
	 *  1,3,.. - corresponding strengths
	 * @param tilesX number of tiles in each data line
	 * @param magic_coeff still not understood coefficient that reduces reported disparity value.  Seems to be around 0.85
	 * @param debugLevel debug level
	 * @return per sub-camera, per direction (x,y) 6 quadratic polynomial coefficients, same format as fine_geometry_correction()
	 */
	public ArrayList<Sample> selectInfinityTiles(
			final double fcorr_radius,
			final boolean use_vertical,
			final double min_strength,
			final double max_diff, // also temporarily maximal difference from 0.0
			final int max_iterations,
			final double max_coeff_diff,
			final double far_pull, //  = 0.2; // 1; //  0.5;
			CLTParameters           clt_parameters,
			double [][] disp_strength,
			int         tilesX,
			double      magic_coeff, // still not understood coefficient that reduces reported disparity value.  Seems to be around 0.85
			int debugLevel)
	{
		final int numTiles =            disp_strength[0].length;
		final int tilesY =              numTiles/tilesX;
		final boolean []  center_mask = getCenterMask(fcorr_radius, tilesX, tilesY);
		double [] disparity_poly = new double[6];
		PolynomialApproximation pa = new PolynomialApproximation();
		double thresholdLin = 1.0E-20;  // threshold ratio of matrix determinant to norm for linear approximation (det too low - fail)
		double thresholdQuad = 1.0E-30; // threshold ratio of matrix determinant to norm for quadratic approximation (det too low - fail)
		double [] disp_surface = new double[numTiles];
		ArrayList<Sample> samples_list = new ArrayList<Sample>();
		// start with average disparity:
		double sdw= 0.0, sw = 0.0;
		int num_samples = 0;
		for (int num_set = 0; num_set < disp_strength.length/NUM_SLICES; num_set++){
			int disp_index = NUM_SLICES * num_set;
			int str_index = NUM_SLICES * num_set + 1;
			for (int nTile = 0; nTile < numTiles; nTile++) if (center_mask[nTile]){
				if (disp_strength[str_index][nTile] > min_strength) {
					//clt_parameters.fcorr_inf_diff
					if (Math.abs(disp_strength[disp_index][nTile]) <= max_diff) {
						double weight= disp_strength[str_index][nTile];
						sdw += weight * disp_strength[disp_index][nTile];
						sw +=  weight;
						num_samples ++;
					} else {
						disp_strength[str_index][nTile] = 0.0;
					}
				}
			}
		}
		if (sw != 0) {
			sdw /= sw;
		}
		disparity_poly[5] = sdw;
		if (debugLevel > 0) {
			System.out.println("selectInfinityTiles(): Number of input samples exceeding "+min_strength+" strength is "+num_samples+", average disparity is "+sdw);
		}



		for (int pass = 0; pass < max_iterations; pass++){
			for (int nTile = 0; nTile < numTiles; nTile++){
				int tileX = nTile % tilesX;
				int tileY = nTile / tilesX;
				double x =  (2.0 * tileX)/tilesX - 1.0; // -1.0 to +1.0;
				double y =  (2.0 * tileY)/tilesY - 1.0; // -1.0 to +1.0
				disp_surface[nTile] =
						disparity_poly[0] * x * x +
						disparity_poly[1] * y * y +
						disparity_poly[2] * x * y +
						disparity_poly[3] * x +
						disparity_poly[4] * y +
						disparity_poly[5];
			}

			samples_list.clear();
			for (int num_set = 0; num_set < disp_strength.length/NUM_SLICES; num_set++){
				int disp_index = NUM_SLICES * num_set;
				int str_index = NUM_SLICES * num_set + 1;
				for (int nTile = 0; nTile < numTiles; nTile++) if (center_mask[nTile]){
					if ((disp_strength[str_index][nTile] > min_strength) &&
							//							  (Math.abs(disp_strength[disp_index][nTile] - disp_surface[nTile]) < max_diff)){
							((disp_strength[disp_index][nTile] - disp_surface[nTile]) < max_diff)){
						double weight= disp_strength[str_index][nTile];
						// next are for far tiles, that differ by more than max_diff
						if (Math.abs(disp_strength[disp_index][nTile] - disp_surface[nTile]) > max_diff){
							weight *= far_pull;
						}
						samples_list.add(new Sample(num_set, nTile, weight));
					}
				}
			}
			double [][][] mdata;
			if (use_vertical) {
				mdata = new double[samples_list.size()][3][];
			} else {
				mdata = new double [1][samples_list.size()][3];
			}

			int indx = 0;
			for (Sample s: samples_list){
				int tileX = s.tile % tilesX;
				int tileY = s.tile / tilesX;
				if (use_vertical) {
					mdata[indx][0] = new double [2];
					mdata[indx][0][0] =  (2.0 * tileX)/tilesX - 1.0; // -1.0 to +1.0;
					mdata[indx][0][1] =  (2.0 * tileY)/tilesY - 1.0; // -1.0 to +1.0
					mdata[indx][1] = new double [1];
					mdata[indx][1][0] =  (disp_strength[2 * s.series + 0][s.tile]/magic_coeff - disp_surface[s.tile]); // disparity residual
					mdata[indx][2][0] =  s.weight; // disp_strength[2 * s.series + 1][s.tile]; // strength
					//					  if (Math.abs( mdata[indx][1][0]) > max_diff) { // far tiles
					//						  mdata[indx][2][0] *= far_pull;
					//					  }

				} else {
					mdata[0][indx][0] = (2.0 * tileX)/tilesX - 1.0; // -1.0 to +1.0;
					mdata[0][indx][1] = (disp_strength[2 * s.series + 0][s.tile]/magic_coeff - disp_surface[s.tile]); // disparity residual
					mdata[0][indx][2] = s.weight; // disp_strength[2 * s.series + 1][s.tile]; // strength
					//						  if (Math.abs( mdata[0][indx][1]) > max_diff) { // far tiles
					//							  mdata[0][indx][2] *= far_pull;
					//						  }
				}
				indx ++;
			}
			if ((debugLevel > 2) && (pass < 21)){
				String [] titles = {"disparity","approx","diff", "strength"};
				double [][] dbg_img = new double [titles.length][numTiles];
				for (int nTile = 0; nTile < numTiles; nTile++){
					int tileX = nTile % tilesX;
					int tileY = nTile / tilesX;
					double x =  (2.0 * tileX)/tilesX - 1.0; // -1.0 to +1.0;
					double y =  (2.0 * tileY)/tilesY - 1.0; // -1.0 to +1.0
					dbg_img[1][nTile] =
							disparity_poly[0] * x * x +
							disparity_poly[1] * y * y +
							disparity_poly[2] * x * y +
							disparity_poly[3] * x +
							disparity_poly[4] * y +
							disparity_poly[5];
				}

				for (Sample s: samples_list){
					dbg_img[0][s.tile] += disp_strength[2 * s.series][s.tile] * s.weight; // disp_strength[2 * p.x + 1][p.y];
					dbg_img[3][s.tile] += s.weight; // disp_strength[2 * p.x + 1][p.y];
				}
				for (int nTile = 0; nTile < numTiles; nTile++) {
					if (dbg_img[3][nTile] > 0.0) {
						dbg_img[0][nTile] /=dbg_img[3][nTile];
					} else {
						dbg_img[0][nTile] = Double.NaN;
					}
				}

				for (int nTile = 0; nTile < numTiles; nTile++) {
					if (dbg_img[3][nTile] > 0.0) {
						dbg_img[2][nTile] = dbg_img[0][nTile] - dbg_img[1][nTile];
					} else {
						dbg_img[2][nTile] = Double.NaN;
					}
				}

				(new ShowDoubleFloatArrays()).showArrays(dbg_img, tilesX, tilesY, true, "infinity_"+pass, titles);
			}





			double [][] coeffs = new double[1][6];
			if (use_vertical){
				double[][] approx2d = pa.quadraticApproximation(
						mdata,
						!clt_parameters.fcorr_inf_quad, // boolean forceLinear,  // use linear approximation
						thresholdLin,  // threshold ratio of matrix determinant to norm for linear approximation (det too low - fail)
						thresholdQuad,  // threshold ratio of matrix determinant to norm for quadratic approximation (det too low - fail)
						debugLevel); {

						}
						if (approx2d[0].length == 6) {
							coeffs = approx2d;
						} else {
							for (int i = 0; i < 3; i++){
								coeffs[0][3+i] = approx2d[0][i];
							}
						}
			} else {
				double [] approx1d = pa.polynomialApproximation1d(mdata[0], clt_parameters.fcorr_inf_quad ? 2 : 1);
				//			  disparity_approximation = new double[6];
				coeffs[0][5] = approx1d[0];
				coeffs[0][3] = approx1d[1];
				if (approx1d.length > 2){
					coeffs[0][0] = approx1d[2];
				}
			}
			if (debugLevel > -1){
				System.out.println(String.format(
						"infinityCorrection() disparity pass=%03d A=%8.5f B=%8.5f C=%8.5f D=%8.5f E=%8.5f F=%8.5f",
						pass, disparity_poly[0], disparity_poly[1], disparity_poly[2],
						disparity_poly[3], disparity_poly[4], disparity_poly[5]) );
			}
			boolean all_done = true;
			for (int i = 0; i < disparity_poly.length; i++){
				if (Math.abs(coeffs[0][i]) > max_coeff_diff){
					all_done = false;
					break;
				}
			}

			if (all_done) break;
			for (int i = 0; i < disparity_poly.length; i++){
				disparity_poly[i] += coeffs[0][i];
			}
		}
		return samples_list;
		// use last generated samples_list;
	}

	/**
	 * Calculate quadratic polynomials for each subcamera X/Y correction to match disparity = 0 at infinity
	 * Next parameters are made separate to be able to modify them between different runs keeping clt_parameters
	 * @param clt_parameters CLT parameters
	 * @param disp_strength array of a single or multiple disparity/strength pairs (0,2, .. - disparity,
	 *  1,3,.. - corresponding strengths
	 * @param samples_list sample list generated by selectInfinityTiles method, each element references measurement series,
	 * tile index and (possibly modified) weight of each tile
	 * @param tilesX number of tiles in each data line
	 * @param magic_coeff still not understood coefficient that reduces reported disparity value.  Seems to be around 0.85
	 * @param debugLevel debug level
	 * @return per sub-camera, per direction (x,y) 6 quadratic polynomial coefficients, same format as fine_geometry_correction()
	 */
	public double [][][] infinityCorrection_old(
			//			final double min_strength0,
			//			final double max_diff0,
			//			final int max_iterations0,
			//			final double max_coeff_diff0,
			//			final double far_pull0, //  = 0.2; // 1; //  0.5;
			final boolean use_vertical,
			CLTParameters           clt_parameters,
			double [][] disp_strength,
			ArrayList<Sample> samples_list,
			int         tilesX,
			double      magic_coeff, // still not understood coefficient that reduces reported disparity value.  Seems to be around 8.5
			int debugLevel)
	{
		final int numTiles =            disp_strength[0].length;
		final int tilesY =              numTiles/tilesX;
		double [] disparity_poly = new double[6];
		PolynomialApproximation pa = new PolynomialApproximation();
		double thresholdLin = 1.0E-20;  // threshold ratio of matrix determinant to norm for linear approximation (det too low - fail)
		double thresholdQuad = 1.0E-30; // threshold ratio of matrix determinant to norm for quadratic approximation (det too low - fail)

		// use last generated samples_list;

		double [][][] mdata;
		if (use_vertical) {
			mdata = new double[samples_list.size()][3][];
		} else {
			mdata = new double [8][samples_list.size()][3];
		}
		int indx = 0;
		final Matrix [] corr_rots = qc.geometryCorrection.getCorrVector().getRotMatrices(); // get array of per-sensor rotation matrices
		for (Sample s: samples_list){
			int tileX = s.tile % tilesX;
			int tileY = s.tile / tilesX;
			double centerX = tileX * qc.tp.getTileSize() + qc.tp.getTileSize()/2;// - shiftX;
			double centerY = tileY * qc.tp.getTileSize() + qc.tp.getTileSize()/2;//- shiftY;

			double [][] centersXY_disp = qc.geometryCorrection.getPortsCoordinatesAndDerivatives(
					qc.geometryCorrection, //			GeometryCorrection gc_main,
					false,          // boolean use_rig_offsets,
					corr_rots, // Matrix []   rots,
					null,      //  Matrix [][] deriv_rots,
					null,      // double [][] pXYderiv, // if not null, should be double[8][]
					null,      // double [][] disp_dist used to correct 3D correlations
					centerX,
					centerY,
					disp_strength[2 * s.series + 0][s.tile]/magic_coeff); // disparity
			double [][] centersXY_inf = qc.geometryCorrection.getPortsCoordinatesAndDerivatives(
					qc.geometryCorrection, //			GeometryCorrection gc_main,
					false,          // boolean use_rig_offsets,
					corr_rots, // Matrix []   rots,
					null,      //  Matrix [][] deriv_rots,
					null,      // double [][] pXYderiv, // if not null, should be double[8][]
					null,      // double [][] disp_dist used to correct 3D correlations
					centerX,
					centerY,
					0.0); // disparity

			for (int i = 0; i < centersXY_disp.length;i++){
				centersXY_disp[i][0] -= centersXY_inf[i][0];
				centersXY_disp[i][1] -= centersXY_inf[i][1];
			}
			if (use_vertical) {
				mdata[indx][0] = new double [2];
				mdata[indx][0][0] =  (2.0 * tileX)/tilesX - 1.0; // -1.0 to +1.0;
				mdata[indx][0][1] =  (2.0 * tileY)/tilesY - 1.0; // -1.0 to +1.0
				mdata[indx][1] = new double [8];
				for (int n = 0; n < 8; n++){
					mdata[indx][1][n] =   -centersXY_disp[n / 2][n % 2];
				}
				mdata[indx][2] = new double [1];
				mdata[indx][2][0] = s.weight; //  disp_strength[2 * p.x + 1][p.y]; // strength

			} else {
				for (int n = 0; n < 8; n++){
					mdata[n][indx][0] = (2.0 * tileX)/tilesX - 1.0; // -1.0 to +1.0;
					mdata[n][indx][1] = -centersXY_disp[n / 2][n % 2];
					mdata[n][indx][2] = s.weight; // disp_strength[2 * p.x + 1][p.y]; // strength
				}
			}
			indx ++;
		}
		double [][] coeffs = new double[8][6];
		if (use_vertical){
			double [][] approx2d = pa.quadraticApproximation(
					mdata,
					!clt_parameters.fcorr_inf_quad, // boolean forceLinear,  // use linear approximation
					thresholdLin,  // threshold ratio of matrix determinant to norm for linear approximation (det too low - fail)
					thresholdQuad,  // threshold ratio of matrix determinant to norm for quadratic approximation (det too low - fail)
					debugLevel);
			for (int n = 0; n < 8; n++){
				if (approx2d[n].length == 6) {
					coeffs[n] = approx2d[n];
				} else {
					for (int i = 0; i < 3; i++){
						coeffs[n][3+i] = approx2d[n][i];
					}
				}
			}
		} else {
			for (int n = 0; n < 8; n++){
				double [] approx1d = pa.polynomialApproximation1d(mdata[n], clt_parameters.fcorr_inf_quad ? 2 : 1);
				//			  disparity_approximation = new double[6];
				coeffs[n][5] = approx1d[0];
				coeffs[n][3] = approx1d[1];
				if (approx1d.length > 2){
					coeffs[n][0] = approx1d[2];
				}
			}
		}

		double [][][] inf_corr = new double [4][2][];
		for (int n = 0; n < 8; n++){
			inf_corr[n / 2][n % 2] = coeffs[n];
		}

		if (debugLevel > 0) {
			String [] titles = {"disparity","approx","diff", "strength"};
			double [][] dbg_img = new double [titles.length][numTiles];
			for (int nTile = 0; nTile < numTiles; nTile++){
				int tileX = nTile % tilesX;
				int tileY = nTile / tilesX;
				double x =  (2.0 * tileX)/tilesX - 1.0; // -1.0 to +1.0;
				double y =  (2.0 * tileY)/tilesY - 1.0; // -1.0 to +1.0
				dbg_img[1][nTile] =
						disparity_poly[0] * x * x +
						disparity_poly[1] * y * y +
						disparity_poly[2] * x * y +
						disparity_poly[3] * x +
						disparity_poly[4] * y +
						disparity_poly[5];
			}

			for (Sample s: samples_list){
				dbg_img[0][s.tile] += disp_strength[2 * s.series][s.tile] * s.weight; // disp_strength[2 * p.x + 1][p.y];
				dbg_img[3][s.tile] += s.weight; // disp_strength[2 * s.series + 1][s.tile];
			}

			for (int nTile = 0; nTile < numTiles; nTile++) {
				if (dbg_img[3][nTile] > 0.0) {
					dbg_img[0][nTile] /=dbg_img[3][nTile];
				} else {
					dbg_img[0][nTile] = Double.NaN;
				}
			}

			for (int nTile = 0; nTile < numTiles; nTile++) {
				if (dbg_img[3][nTile] > 0.0) {
					dbg_img[2][nTile] = dbg_img[0][nTile] - dbg_img[1][nTile];
				} else {
					dbg_img[2][nTile] = Double.NaN;
				}
			}

			(new ShowDoubleFloatArrays()).showArrays(dbg_img, tilesX, tilesY, true, "AC_infinityCorrection", titles);

		}
		return inf_corr;
	}

	/**
	 * Correct channel mismatch (preserving disparity) using the same tiles as those for correcting disparity
	 * at infinity
	 * Next parameters are made separate to be able to modify them between different runs keeping clt_parameters
	 * @param use_poly calculate quadratic/linear correction (can not work with high disparity).
	 *                 False - adjust extrinsics (tilt, azimuth, roll)
	 * @param clt_parameters CLT parameters
	 * @param disp_strength array of a single or multiple disparity/strength pairs (0,2, .. - disparity,
	 *  1,3,.. - corresponding strengths
	 * @param samples_list sample list generated by selectInfinityTiles method, each element references measurement series,
	 * tile index and (possibly modified) weight of each tile
	 * @param tilesX number of tiles in each data line
	 * @param magic_coeff understood - interaction of the CM maximum and correlation window
	 * @param mismatch_list data to calculate extrinsic corrections or null
	 * @param debugLevel debug level
	 * @return per sub-camera, per direction (x,y) 6 quadratic polynomial coefficients, same format as fine_geometry_correction()
	 */
	public double [][][] infinityMismatchCorrection(
			final double  disp_scan_start,
			final double  disp_scan_step,
			final boolean use_poly,
			final boolean use_quadratic,
			final boolean use_vertical,
			final boolean use_disparity, // for infinity // now disabled?
			final boolean allow_dispatity,
			CLTParameters           clt_parameters,
			double [][] disp_strength,
			ArrayList<Sample> samples_list,
			int         tilesX,
			double      magic_coeff, // still not understood coefficient that reduces reported disparity value.  Seems to be around 8.5
			ArrayList<Mismatch> mismatch_list,
			int debugLevel)
	{
		int dbgTileX = 100;
		int dbgTileY = 100;
		// Mismatch data has disparity values already subtracted, so to correct disparity at infinity, disparity values should be restored
		final int num_tiles =            disp_strength[0].length;
		final int tilesY =              num_tiles/tilesX;
		PolynomialApproximation pa = new PolynomialApproximation();
		double thresholdLin = 1.0E-20;  // threshold ratio of matrix determinant to norm for linear approximation (det too low - fail)
		double thresholdQuad = 1.0E-30; // threshold ratio of matrix determinant to norm for quadratic approximation (det too low - fail)
		//		final int [] indices_mismatch = {1,4,6,9}; //   "dy0", "dy1", "dx2", "dx3"
		//		final int [] indices_mismatch = {2+1, 2+4, 2+6, 2+9}; //   "dy0", "dy1", "dx2", "dx3"
		//		final int [] indices_mismatch = {2+1, 2+3, 2+4, 2+6}; //   "dy0", "dy1", "dx2", "dx3"
		final int [][] indices_mismatch = {
				{2+0, 2+2, 2+4, 2+6},  // "dx0", "dx1", "dx2", "dx3"
				{2+1, 2+3, 2+5, 2+7}}; // "dy0", "dy1", "dy2", "dy3"

		// use last generated samples_list;

		double [][][] mdata;
		if (use_vertical) {
			mdata = new double[samples_list.size()][3][];
		} else {
			mdata = new double [8][samples_list.size()][3];
		}
		int indx = 0;
		double [][] A_arr = {
				{2.0, 1.0, 1.0},
				{1.0, 3.0, 2.0},
				{1.0, 2.0, 3.0}};
		Matrix A = new Matrix(A_arr);
		Matrix AINV = A.inverse();
		double scale = 1.0; // Why was it here? 0.5/magic_coeff;

		double [][] dbg_xy = null;
		if (clt_parameters.show_extrinsic && (debugLevel > -2)) { // TODO: Add clt_parameters
//			dbg_xy = new double [9][num_tiles];
			dbg_xy = new double [18][num_tiles];
		}
		for (Sample s: samples_list){
			int tileX = s.tile % tilesX;
			int tileY = s.tile / tilesX;
			if ((debugLevel > 0) && (tileX == dbgTileX)  && (tileY == dbgTileY)) {
				System.out.println("infinityMismatchCorrection(): tileX = "+tileX+", tileY = "+tileY);
			}
			double [] xy = new double[8]; // same as coefficients: x0,y0,x1,y1,x2,y2,x3,y3
			// Calculate x0,x1,x2,x3 and y0,y1,y2,y3 assuming x0+x1+x2+x3 = 0,y0+y1+y2+y3 = 0 and minimizing squares of errors
			// as each each 4: "dx0", "dx1", "dx2", "dx3" and "dy0", "dy1", "dy2", "dy3" are over-defined
			double [][] dxy = new double[4][2];
			for (int dir = 0; dir < 2; dir++) { // 0 - X, 1 - Y
//				double [] dxy = new double[4];
				for (int i = 0; i < 4; i++){
					dxy[i][dir] = scale * disp_strength[indices_mismatch[dir][i] + (s.series * NUM_SLICES)][s.tile];
					if (Double.isNaN(dxy[i][dir])) {
						int ii = indices_mismatch[dir][i] + (s.series * NUM_SLICES);
						System.out.println("**** BUG: infinityMismatchCorrection() tileX="+tileX+", tileY="+tileY+" dxy["+i+"]["+dir+"]= "+dxy[i][dir]+
								", index="+ii+", s.tile = "+s.tile);
						System.out.println("**** BUG: infinityMismatchCorrection() tileX="+tileX+", tileY="+tileY+" dxy["+i+"]["+dir+"]= "+dxy[i][dir]+
								", index="+ii+", s.tile = "+s.tile);
					}
				}

				/*
use best fit as they are overdefined
y1 - y0 = 2 * dy0
y3 - y1 = 2 * dy3
y3 - y2 = 2 * dy1
y2 - y0 = 2 * dy2

x2 - x0 = 2 * dx2
x3 - x2 = 2 * dx1
x3 - x1 = 2 * dx3
x1 - x0 = 2 * dx0

x1 - x0 = 2 * dx0
x3 - x1 = 2 * dx3
x3 - x2 = 2 * dx1
x2 - x0 = 2 * dx2

y0 + y1 + y2 + y3 = 0

x0 + x1 + x2 + x3 = 0

e1 = (y1 - y0  - 2*dy0)^2 +
     (y2 - y0  - 2*dy2)^2 +
     (y3 - y1  - 2*dy3)^2 +
     (y3 - y2  - 2*dy1)^2 =

     (y1 - y0  - 2*dy0)^2 +
     (y2 - y0  - 2*dy2)^2 +
     (-(y0+y1+y2) - y1  - 2*dy3)^2 +
     (-(y0+y1+y2) - y2  - 2*dy1)^2 =

     (y1 - y0  - 2*dy0)^2 +
     (y2 - y0  - 2*dy2)^2 +
     (y0 + 2 * y1 + y2  + 2*dy3)^2 +
     (y0 + y1 + 2 * y2  + 2*dy1)^2

de1/dy0 = -2 * (y1 - y0  - 2*dy0)
          -2 * (y2 - y0  - 2*dy2)
          +2 * (y0 + 2 * y1 + y2  + 2*dy3)
          +2 * (y0 + y1 + 2 * y2  + 2*dy1) =

          2 * (-(y1 - y0  - 2*dy0)-(y2 - y0  - 2*dy2)+ (y0 + 2 * y1 + y2  + 2*dy3) + (y0 + y1 + 2 * y2  + 2*dy1)) =

          2 * ((2*dy0 + y0 - y1) + (2*dy2 + y0 - y2) + (y0 + 2 * y1 + y2  + 2*dy3) + (y0 + y1 + 2 * y2  + 2*dy1)) =
          2 * (2 * (dy0+dy1+dy2+dy3) + y0 * 4 + 2 * y1 + 2 * y2)

2 * y0 + 1 * y1 + 1 * y2 + (dy0+dy1+dy2+dy3) = 0;

de1/dy1 = +2 * (y1 - y0  - 2*dy0)
          +2 * (2*y0 + 4 * y1 + 2*y2  + 4*dy3)
          +2 * (y0 + y1 + 2 * y2  + 2*dy1) =
          2* ( 2*(-dy0 +2*dy3 +dy1) + y0 * 2 + y1 * 6 + y2 * 4
1 * y0 + 3 * y1 + 2 * y2 + (-dy0 +2*dy3 +dy1)  = 0

de1/dy2 = +2 * (y2 - y0  - 2*dy2)
          +2 * (y0 + 2 * y1 + y2  + 2*dy3)
          +2 * (2*y0 + 2*y1 + 4 * y2  + 4*dy1) =

          2*(2*(-dy2 + dy3 +2*dy1) + y0 * 2 + y1 * 4 + y2 * 6
1 * y0 + 2 * y1 + 3 * y2 + (-dy2 + dy3 +2*dy1) = 0;





2 * y0 + 1 * y1 + 1 * y2 + (dy0+dy1+dy2+dy3) = 0;
1 * y0 + 3 * y1 + 2 * y2 + (-dy0 +2*dy3 +dy1)  = 0
1 * y0 + 2 * y1 + 3 * y2 + (-dy2 +dy3 + 2*dy1) = 0;
    | 2  1  1 |
A = | 1  3  2 |
    | 1  2  3 |

    |-dy0   -dy1 -dy2   -dy3 |
B = |+dy0   -dy1      -2*dy3 |
    |+dy2 -2*dy1        -dy3 |
				 */

				/*
				 *			    |-dy0   -dy1 -dy2   -dy3 |
				 *			B = |+dy0   -dy1      -2*dy3 |
				 *			    |+dy2 -2*dy1        -dy3 |
				 */

				double [] B_arr = {
					   -dxy[0][dir]      -dxy[1][dir] -dxy[2][dir] -dxy[3][dir],
						dxy[0][dir]      -dxy[1][dir]     -2 * dxy[3][dir],
						dxy[2][dir]  -2 * dxy[1][dir]         -dxy[3][dir]};
				Matrix B = new Matrix(B_arr, 3); // 3 rows
				Matrix X = AINV.times(B);
				for (int i = 0; i < 3; i++) {
					xy[2 * i + dir] =  X.get(i, 0);
					xy[2 * 3 + dir] -= X.get(i, 0);
				}
			}
			if (use_disparity) {
				double d = disp_strength[s.series * NUM_SLICES + 0][s.tile];
				xy[0] -= d;	xy[1] -= d;
				xy[2] += d;	xy[3] -= d;
				xy[4] -= d;	xy[5] += d;
				xy[6] += d;	xy[7] += d;
			}
			if (dbg_xy != null){
				for (int i = 0; i < xy.length; i++){
					dbg_xy[i][s.tile] += xy[i] * s.weight;
				}
				dbg_xy[9][s.tile] += s.weight*(-xy[0]-xy[1]+xy[2]-xy[3]-xy[4]+xy[5]+xy[6]+xy[7])/8;
				dbg_xy[8][s.tile] += s.weight;
				for (int i = 0; i < dxy.length;i++) {
					dbg_xy[2* i +10][s.tile] += dxy[i][0] * s.weight;
					dbg_xy[2* i +11][s.tile] += dxy[i][1] * s.weight;
				}
			}


			if (use_poly) {
				if (use_vertical) {
					// 2d optimization for 4 functions of x, y
					mdata[indx][0] = new double [2];
					mdata[indx][0][0] =  (2.0 * tileX)/tilesX - 1.0; // -1.0 to +1.0;
					mdata[indx][0][1] =  (2.0 * tileY)/tilesY - 1.0; // -1.0 to +1.0
					mdata[indx][1] = xy; // new double [8]; // 4

					mdata[indx][2] = new double [1];
					mdata[indx][2][0] = s.weight; //  disp_strength[2 * p.x + 1][p.y]; // strength

				} else {
					// 4 individual 1d optimization for functions of x only
					for (int n = 0; n < xy.length; n++){
						mdata[n][indx][0] = (2.0 * tileX)/tilesX - 1.0; // -1.0 to +1.0;
						mdata[n][indx][1] = xy[n];
						mdata[n][indx][2] = s.weight; // disp_strength[2 * p.x + 1][p.y]; // strength
					}
				}
			} else if (mismatch_list != null){
				//				double centerX = tileX * transform_size + transform_size/2 - shiftX;
				//				double centerY = tileY * transform_size + transform_size/2 - shiftY;
				double [] centerXY = {
						tileX * qc.tp.getTileSize() + qc.tp.getTileSize()/2,// - shiftX;
						tileY * qc.tp.getTileSize() + qc.tp.getTileSize()/2};//- shiftY;
				double disparity_task = 	disp_scan_start + disp_scan_step * s.series; // Not needed even with known disparity
				double disparity_meas = 	disp_strength[s.series * NUM_SLICES + 0][s.tile];
				double strength =    	disp_strength[s.series * NUM_SLICES + 1][s.tile];
				if (Double.isNaN(disparity_meas)) {
					System.out.println("infinityMismatchCorrection(): disparity_meas=NaN: s.tile= "+s.tile);
				}
				mismatch_list.add(new Mismatch(
						allow_dispatity && (s.series == 0), // true,  //false, // public boolean   use_disparity; // adjust dx0+dx1+dy0+dy1 == 0
						centerXY,
						disparity_task,  // not used
						disparity_meas,  // not used
						strength,
						dxy)); // xy));
			}
			indx ++;
		}
		if (dbg_xy != null){
			for (int nTile = 0; nTile < num_tiles; nTile++){
				if (dbg_xy[8][nTile] > 0.0){
					for (int i = 0; i < dbg_xy.length; i++) if (i != 8){
						dbg_xy[i][nTile] /= dbg_xy[8][nTile];
					}
				} else {
					for (int i = 0; i < dbg_xy.length; i++) if (i != 8){
						dbg_xy[i][nTile] = Double.NaN;
						dbg_xy[9][nTile] = Double.NaN;
					}
				}
			}
			String [] titles = {"x0", "y0", "x1", "y1", "x2", "y2", "x3","y3","weight","~disp","dx0","dy0","dx1","dy1","dx2","dy2","dx3","dy3"};
			(new ShowDoubleFloatArrays()).showArrays(
					dbg_xy,
					tilesX,
					tilesY,
					true,
					"xy_mismatch",
					titles);
		}

		if (!use_poly){
			return null; // for extrinsics - use mismatch_list
		}

		double [][] coeffs = new double[8][6];
		if (use_vertical){
			double [][] approx2d = pa.quadraticApproximation(
					mdata,
					!use_quadratic, // boolean forceLinear,  // use linear approximation
					thresholdLin,  // threshold ratio of matrix determinant to norm for linear approximation (det too low - fail)
					thresholdQuad,  // threshold ratio of matrix determinant to norm for quadratic approximation (det too low - fail)
					debugLevel);
			for (int n = 0; n < approx2d.length; n++){
				if (approx2d[n].length == 6) {
					coeffs[n] = approx2d[n];
				} else {
					for (int i = 0; i < 3; i++){
						coeffs[n][3+i] = approx2d[n][i];
					}
				}
			}
		} else {
			for (int n = 0; n < mdata.length; n++){
				double [] approx1d = pa.polynomialApproximation1d(mdata[n], use_quadratic ? 2 : 1);
				coeffs[n][5] = approx1d[0];
				coeffs[n][3] = approx1d[1];
				if (approx1d.length > 2){
					coeffs[n][0] = approx1d[2];
				}
			}
		}

		// convert to 8 sets of coefficient for px0, py0, px1, py1, ... py3.
		double [][][] coeff_full = new double [4][2][6];
		for (int j = 0; j < coeffs[0].length; j++){
			coeff_full[0][0][j] = coeffs[0][j];
			coeff_full[0][1][j] = coeffs[1][j];
			coeff_full[1][0][j] = coeffs[2][j];
			coeff_full[1][1][j] = coeffs[3][j];
			coeff_full[2][0][j] = coeffs[4][j];
			coeff_full[2][1][j] = coeffs[5][j];
			coeff_full[3][0][j] = coeffs[6][j];
			coeff_full[3][1][j] = coeffs[7][j];
		}
		return coeff_full;
	}

	public double[][] filterLazyEyePairs (
			final double[][] samples_in,
			final int        smpl_side, // 8 x8 masked, 16x16 sampled
			final double     rms_max,
			final double     frac_keep,
			final int        min_samples,
			final boolean    norm_center, // if there are more tiles that fit than minsamples, replace with a single equal weight
			final int        tilesX)
	{
		if (samples_in.length > NUM_SLICES){
			final double [][] samples = new double [samples_in.length][];
			for (int nfirst = 0; nfirst < samples_in.length; nfirst += NUM_SLICES){
				double [][] ds = new double[NUM_SLICES][]; // {disp_strength_in[i],disp_strength_in[i+1]};
				for (int slice = 0; slice < NUM_SLICES; slice++){
					ds[slice] = samples[nfirst+slice];
				}
				double [][] ds_rslt = filterLazyEyePairs (
						ds,
						smpl_side, // 8 x8 masked, 16x16 sampled
						rms_max,
						frac_keep,
						min_samples,
						norm_center,
						tilesX);
				for (int slice = 0; slice < NUM_SLICES; slice++){
					samples[nfirst+slice] = ds_rslt[slice];
				}
			}
			return samples;
		}
		final int num_pairs = 4;
		final int num_tiles = samples_in[0].length;
		int tilesY = num_tiles/tilesX;
		final double [][] samples = new double [NUM_SLICES][num_tiles];
		//		final int low_lim = -smpl_side/2;
		//		final int high_lim = 3 * smpl_side/2;
		final int step = norm_center ? smpl_side : 1;
		final int start = norm_center ? smpl_side/2 : 0;
		for (int tY = start; tY < tilesY; tY += step){
			for (int tX = start; tX < tilesX ; tX += step){
				ArrayList<Integer> sample_list = new ArrayList<Integer>();
				for (int sY = -smpl_side; sY < smpl_side; sY++){
					int y = tY + sY;
					if ((y >= 0) && (y <tilesY)) {
						for (int sX = -smpl_side; sX < smpl_side; sX++){
							int x = tX + sX;
							if ((x >= 0) && (x < tilesX)) {
								int nTile = y * tilesX + x;
								double w = samples_in[1][nTile];
								if (w > 0.0){
									sample_list.add(nTile);
								}
							}
						}
					}
				}
				int samples_left = (int) (frac_keep * sample_list.size());
				if (samples_left >=  min_samples){
					double [] s1 = new double [2 * num_pairs];
					double [] s2 = new double [2 * num_pairs];
					double sw = 0;
					for (Integer nTile: sample_list){
						double w = samples_in[1][nTile];
						for (int i = 0; i < s1.length; i++){
							double d = samples_in[2 + i][nTile];
							s1[i] += w * d;
							s2[i] += w *d *d;
						}
						sw += w;
					}
					double rms = 0.0;
					for (int i = 0; i < s1.length; i++){
						s1[i] /= sw; // average
						s2[i] /= sw;
						rms += s2[i] - s1[i] * s1[i];
					}
					rms = Math.sqrt(rms/s1.length);
					//					if (rms_max > 0.0){
					//						System.out.println("filterLazyEyePairs(): tX="+tX+", tY="+ tY+" rms = "+rms);
					//					}
					if (rms <= rms_max){ // otherwise discard all the tile
						while (sample_list.size() > samples_left){
							int worst_index = -1;
							double worst_var  = 0;
							for (int i = 0; i < sample_list.size(); i++){
								int nTile = sample_list.get(i);
								double sv = 0;
								for (int n = 0; n < s1.length; n++){
									double d = samples_in[2 + n][nTile] - s1[n];
									sv += d * d;
								}
								if (sv >  worst_var) {
									worst_var = sv;
									worst_index = i;
								}
							}
							// remove worst
							int nTile = sample_list.remove(worst_index);
							double w = samples_in[1][nTile];
							for (int n = 0; n < s1.length; n++){
								s1[n] *= sw; // now it is weighted sum
								double d = samples_in[2 + n][nTile];
								s1[n] -= d * w;
								s1[n] /= (sw - w); // average again
							}
							sw -= w;

						}
						if (norm_center) {
							double sd = 0.0;
							//							double sw1 = 0.0;
							for (Integer nTile:sample_list) {
								double w = samples_in[1][nTile];
								double d = samples_in[0][nTile];
								sd += w*d;
								//								sw1 += w;
							}
							sd /= sw;
							int nTile =tY * tilesX + tX;
							samples[0][nTile] = sd;
							samples[1][nTile] = 1.0; // sw; // equal weight
							for (int j = 0; j < 2 * num_pairs; j++) {
								samples[2 + j][nTile] = s1[j];
							}
						} else { // if (norm_center)
							int nTile =tY * tilesX + tX;
							if (sample_list.contains(nTile)){
								for (int j = 0; j < samples.length; j++) {
									samples[j][nTile] = samples_in[j][nTile];
								}
							}

						} // if (norm_center) else
					}
				}
			}
		}
		return samples;

	}


	public double[][] filterHistogramFar (
			final double[][] disp_strength_in,
			final int        smpl_side, // 8 x8 masked, 16x16 sampled
			final double     disp_min,
			final double     disp_step,
			final int        num_bins,
			final double     sigma,
			final double     max_diff,
			final int        min_samples,
			final boolean    norm_center, // if there are more tiles that fit than minsamples, replace with a single equal weight
			final int        tilesX)
	{
		if (disp_strength_in.length > NUM_SLICES){
			final double [][] disp_strength = new double [disp_strength_in.length][];
			for (int nfirst = 0; nfirst < disp_strength_in.length; nfirst += NUM_SLICES){
				double [][] ds = new double[NUM_SLICES][]; // {disp_strength_in[i],disp_strength_in[i+1]};
				for (int slice = 0; slice < NUM_SLICES; slice++){
					ds[slice] = disp_strength[nfirst+slice];
				}
				double [][] ds_rslt = filterHistogramFar (
						ds,
						smpl_side, // 8 x8 masked, 16x16 sampled
						disp_min,
						disp_step,
						num_bins,
						sigma,
						max_diff,
						min_samples,
						norm_center,
						tilesX);
				for (int slice = 0; slice < NUM_SLICES; slice++){
					disp_strength[nfirst+slice] = ds_rslt[slice];
				}
			}
			return disp_strength;
		}

		final int num_tiles = disp_strength_in[0].length;
		int tilesY = num_tiles/tilesX;
		final double [][] disp_strength = new double [NUM_SLICES][num_tiles];
		//		final int low_lim = -smpl_side/2;
		//		final int high_lim = 3 * smpl_side/2;
		final DoubleGaussianBlur gb=new DoubleGaussianBlur();
		final int step = norm_center ? smpl_side : 1;
		final int start = norm_center ? smpl_side/2 : 0;
		for (int tY = start; tY < tilesY; tY += step){
			for (int tX = start; tX < tilesX ; tX += step){
				double [] hist = new double [num_bins];
				int num_samples = 0;
				for (int sY = -smpl_side; sY < smpl_side; sY++){
					int y = tY + sY;
					if ((y >= 0) && (y <tilesY)) {
						for (int sX = -smpl_side; sX < smpl_side; sX++){
							int x = tX + sX;
							if ((x >= 0) && (x < tilesX)) {
								int nTile = y * tilesX + x;
								double w = disp_strength_in[1][nTile];
								if (w > 0.0){
									int bin = (int) ((disp_strength_in[0][nTile] - disp_min)/disp_step);
									if ((bin >= 0) && (bin < num_bins)){
										hist[bin] += w;
										num_samples ++;
									}
								}
							}
						}
					}
				}
				if (num_samples >= min_samples) {
					if (sigma > 0.0){
						gb.blur1Direction(hist, num_bins, 1, sigma/disp_step, 0.01,true);
					}
				}
				// find first increase
				int incr = -1;
				for (int i = 0; i < (num_bins-2); i++) if (hist[i+1] > hist[i]){
					incr = i;
					break;
				}

				if (incr < 0) continue;

				// find first decrease after increase

				int decr = -1;
				for (int i = incr+1; i < (num_bins-2); i++) if (hist[i+1] < hist[i]){
					decr = i;
					break;
				}
				if (decr < 0) continue;

				// now hist[decr] is a local max, lowest in the range

				double a = 0.5*(hist[decr + 1] -2*hist[decr] + hist[decr - 1]);
				double b = 0.5*(hist[decr + 1] - hist[decr - 1]);
				double dx =  - b/(2*a);
				if (dx > 1.0) dx = 1.0;
				else if (dx < -1.0) dx = -1.0;
				double disp_inf = disp_min + disp_step * (decr + dx);
				if (norm_center) {
					num_samples = 0;
					double sdw = 0.0, sw = 0.0;
					double [] sm = new double [NUM_SLICES - 2];

					for (int sY = -smpl_side/2; sY < smpl_side/2; sY++){
						int y = tY + sY;
						if ((y >= 0) && (y <tilesY)) {
							for (int sX = -smpl_side/2; sX < smpl_side/2; sX++){
								int x = tX + sX;
								if ((x >= 0) && (x <tilesX)) {
									int nTile = y * tilesX + x;
									double w = disp_strength_in[1][nTile];
									if (w > 0.0){
										if (Math.abs(disp_strength_in[0][nTile] - disp_inf) < max_diff) {
											sw += disp_strength_in[1][nTile];
											sdw += disp_strength_in[0][nTile] * disp_strength_in[1][nTile];
											for (int j = 0; j < sm.length; j++){
												sm[j] += disp_strength_in[2 + j][nTile] * disp_strength_in[1][nTile];
											}
											num_samples++;
										}
									}
								}
							}
						}
					}
					if ((num_samples > min_samples) && (sw > 0.0)){
						int nTile = tY * tilesX + tX;
						disp_strength[0][nTile] = sdw/sw; // weighted average disparity
						disp_strength[1][nTile] = 1.0; // equal weight
						for (int j = 0; j < sm.length; j++){
							disp_strength[2 + j][nTile] = sm[j]/sw;
						}
					}
				} else { // if (norm_center)
					// just center tile - fits or not
					int nTile = tY * tilesX + tX;
					double w = disp_strength_in[1][nTile];
					if (w > 0.0){
						if (Math.abs(disp_strength_in[0][nTile] - disp_inf) < max_diff) {
							for (int j = 0; j < disp_strength.length; j++) {
								disp_strength[j][nTile] = disp_strength_in[j][nTile];
							}
						}
					}
				} //if (norm_center) else

			}
		}
		return disp_strength;
	}

	public double[][] filterDisparityStrength (
			final double[][] disp_strength_in,
			final double     strength_floor,
			final double     strength_pow,
			final int        smplSide, //        = 2;      // Sample size (side of a square)
			final int        smplNum, //         = 3;      // Number after removing worst (should be >1)
			final double     smplRms, //         = 0.1;    // Maximal RMS of the remaining tiles in a sample
			final int        tilesX)
	{
		final int dbg_tile = -34145; // 37005;
		if (disp_strength_in.length > NUM_SLICES){
			final double [][] disp_strength = new double [disp_strength_in.length][];
			for (int nfirst = 0; nfirst < disp_strength_in.length; nfirst += NUM_SLICES){
				double [][] ds = new double[NUM_SLICES][]; // {disp_strength_in[i],disp_strength_in[i+1]};
				for (int slice = 0; slice < NUM_SLICES; slice++){
					ds[slice] = disp_strength_in[nfirst+slice];
				}
				double [][] ds_rslt = filterDisparityStrength (
						ds,
						strength_floor,
						strength_pow,
						smplSide, //        = 2;      // Sample size (side of a square)
						smplNum, //         = 3;      // Number after removing worst (should be >1)
						smplRms, //         = 0.1;    // Maximal RMS of the remaining tiles in a sample
						tilesX);
				for (int slice = 0; slice < NUM_SLICES; slice++){
					disp_strength[nfirst+slice] = ds_rslt[slice];
				}
			}
			return disp_strength;
		}

		final int num_tiles = disp_strength_in[0].length;
		int tilesY = num_tiles/tilesX;
		final double [][] disp_strength = new double [NUM_SLICES][num_tiles]; // disp, strength, mismatches
		final int index_shift  = (smplSide /2) * (tilesX + 1) ; // diagonal shift
		final int smplLen = smplSide*smplSide;

		final double [] weight = new double [num_tiles]; // modified weight
		final double [] disp = disp_strength_in[0];
		final double [][] mismatch = new double [NUM_SLICES-2][];
		for (int i = 0; i < mismatch.length; i++) {
			mismatch[i] = disp_strength_in[i + 2];
		}
		final double smlVar = smplRms * smplRms; // maximal variance (weighted average of the squared difference from the mean)
		for (int nTile = 0; nTile < num_tiles; nTile++){
			if (nTile == dbg_tile) {
				System.out.println("filterDisparityStrength().1: nTile = dbg_tile = "+dbg_tile);
			}
			double w = disp_strength_in[1][nTile] - strength_floor;
			if (w > 0){
				if (strength_pow != 1.0) w = Math.pow(w, strength_pow);
				weight[nTile] = w;
			}
		}
		for (int tY = 0; tY < (tilesY - smplSide); tY++){
			for (int tX = 0; tX < (tilesX - smplSide); tX++){
				if ((tY*tilesX + tX + index_shift) == dbg_tile) {
					System.out.println("filterDisparityStrength().2: nTile = dbg_tile = "+dbg_tile);
				}

				int num_in_sample = 0;
				boolean [] smpl_sel = new boolean [smplLen];
				double [] smpl_d =  new double [smplLen];
				double [] smpl_w =  new double [smplLen];
				double [][] smpl_mismatch = new double [NUM_SLICES-2][smplLen];
				for (int sy = 0; sy < smplSide; sy++){
					int y = tY + sy; //  - smpl_center;
					for (int sx = 0; sx < smplSide; sx++){
						int x = tX + sx; // - smpl_center;
						int indx = y * tilesX + x; // absolute center
						if (weight[indx] > 0.0){
							int indxs = sy * smplSide + sx;
							smpl_sel[indxs] = true;
							smpl_d[indxs] = disp[indx];
							smpl_w[indxs] = weight[indx];
							for (int i = 0; i < smpl_mismatch.length; i++){
								smpl_mismatch[i][indxs] = mismatch[i][indx];
							}
							num_in_sample ++;
						}
					}
				}
				if (num_in_sample >= smplNum){ // try, remove worst
					// calculate
					double sd=0.0, sd2 = 0.0, sw = 0.0;
					double [] sm = new double [mismatch.length];
					for (int i = 0; i < smplLen; i++) if (smpl_sel[i]) {
						double dw = smpl_d[i] * smpl_w[i];
						sd += dw;
						sd2 += dw * smpl_d[i];
						sw +=       smpl_w[i];
						for (int j = 0; j < sm.length; j++){
							sm[j] += smpl_mismatch[j][i] * smpl_w[i];
						}
					}
					// remove worst, update sd2, sd and sw
					while ((num_in_sample > smplNum) && (sw > 0)){ // try, remove worst
						double d_mean = sd/sw;
						//						double [] mismatch_mean = new double[sm.length];
						//						for (int j = 0; j < sm.length; j++){
						//							mismatch_mean[j] = sm[j]/sw;
						//						}
						int iworst = -1;
						double dworst2 = 0.0;
						for (int i = 0; i < smplLen; i++) if (smpl_sel[i]) {
							double d2 = (smpl_d[i] - d_mean);
							d2 *=d2;
							if (d2 > dworst2) {
								iworst = i;
								dworst2 = d2;
							}
						}
						if (iworst < 0){
							System.out.println("**** this is a BUG in filterDisparityStrength() ****"); // happened when smpl_d = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
							break;
						}
						// remove worst sample
						smpl_sel[iworst] = false;
						double dw = smpl_d[iworst] * smpl_w[iworst];
						sd -= dw;
						sd2 -= dw * smpl_d[iworst];
						sw -=       smpl_w[iworst];

						for (int j = 0; j < sm.length; j++){
							sm[j] -= smpl_mismatch[j][iworst] * smpl_w[iworst];
						}
						num_in_sample --;
					}
					// calculate variance of the remaining set
					if (sw > 0.0) {
						sd /= sw;
						sd2 /= sw;
						double var = sd2 - sd * sd;
						if (var < smlVar) { // good, save in the result array
							int indx = tY*tilesX + tX + index_shift;
							disp_strength[0][indx] = sd;
							disp_strength[1][indx] = sw;
							for (int j = 0; j < sm.length; j++){
								disp_strength[2 + j][indx] = sm[j]/sw ;
							}
						}
					} else {
						num_in_sample = 0;
						System.out.println("**** this is a BUG in getDisparityStrength(), shoud not happen ? ****");
					}
				}
			}
		}
		return disp_strength;
	}







	public void show_fine_corr(
			double [][][] corr,
			String prefix)
	{
		String sadd = (prefix.length() > 0)?(prefix+" "):"";

		for (int n = 0; n< corr.length; n++){
			for (int i = 0; i< corr[n].length; i++){
				System.out.print(sadd+"port"+n+": "+QuadCLT.fine_corr_dir_names[i]+": ");
				for (int j = 0; j< corr[n][i].length; j++){
					System.out.print(QuadCLT.fine_corr_coeff_names[j]+"="+corr[n][i][j]+" ");
				}
				System.out.println();
			}
		}
	}

	public double [][] getDoubleFromImage(
			ImagePlus imp_src,
			int debugLevel)
	{
		//		  double min_max_ratio = 1.3;
		ImageStack clt_mismatches_stack= imp_src.getStack();
		final int tilesX = clt_mismatches_stack.getWidth(); // tp.getTilesX();
		final int tilesY = clt_mismatches_stack.getHeight(); // tp.getTilesY();
		final int nTiles =tilesX * tilesY;
		final double [][] data = new double [clt_mismatches_stack.getSize()][nTiles];
		for (int n = 0; n < data.length; n++) {
			float [] fpixels = (float[]) clt_mismatches_stack.getPixels(n +1);
			for (int i = 0; i < nTiles; i++){
				data[n][i] = fpixels[i];
			}
		}
		return data;
	}
	public double [][] getFineCorrFromDoubleArray(
			double [][] data,
			int         tilesX,
			int debugLevel)
	{
		//		  double min_max_ratio = 1.3;

		final int tilesY = data[0].length/tilesX; // tp.getTilesY();
		final int nTiles =tilesX * tilesY;
		final int num_scans = data.length/NUM_ALL_SLICES;

		final double [][] scans = new double [num_scans * NUM_ALL_SLICES][nTiles];

		for (int ns = 0; ns < num_scans; ns++){
			//		    	double [][]min_max_strength = new double[2][nTiles];
			for (int pair = 0; pair < 4; pair++){
				scans[ns * NUM_ALL_SLICES + 0] = data[12 * num_scans + ns];
				scans[ns * NUM_ALL_SLICES + 1] = data[13 * num_scans + ns];
				scans[ns * NUM_ALL_SLICES + pair * 3 + 2] = data[(2 * pair + 0) * num_scans + ns];
				scans[ns * NUM_ALL_SLICES + pair * 3 + 3] = data[(2 * pair + 1) * num_scans + ns];
				scans[ns * NUM_ALL_SLICES + pair * 3 + 4] = data[(8 + pair   ) *  num_scans + ns];
			}
		}
		if (debugLevel > 0) {
			String [] prefixes = {"disparity", "strength", "dx0", "dy0", "str0", "dx1", "dy1", "str1", "dx2", "dy2", "str2", "dx3", "dy3", "str3"};
			String [] titles = new String [num_scans * NUM_ALL_SLICES];
			for (int ns = 0; ns < num_scans; ns++){
				for (int i = 0; i < NUM_ALL_SLICES; i++){
					titles[ns * NUM_ALL_SLICES + i] = prefixes[i]+"_"+ns;
				}
			}
			(new ShowDoubleFloatArrays()).showArrays(scans, tilesX, tilesY, true, "scans" , titles);
		}
		return scans;
	}

	public double [][] getFineCorrFromImage_old(
			ImagePlus imp_src,
			int debugLevel)
	{
		//		  double min_max_ratio = 1.3;
		ImageStack clt_mismatches_stack= imp_src.getStack();
		final int tilesX = clt_mismatches_stack.getWidth(); // tp.getTilesX();
		final int tilesY = clt_mismatches_stack.getHeight(); // tp.getTilesY();
		final int nTiles =tilesX * tilesY;
		final int num_scans =  clt_mismatches_stack.getSize()/NUM_ALL_SLICES;

		final double [][] scans = new double [num_scans * NUM_ALL_SLICES][nTiles];

		for (int ns = 0; ns < num_scans; ns++){
			//		    	double [][]min_max_strength = new double[2][nTiles];
			for (int pair = 0; pair < 4; pair++){
				float [][] fset = new float [3][];
				fset[0] = (float[]) clt_mismatches_stack.getPixels(12 * num_scans + ns +1);
				fset[1] = (float[]) clt_mismatches_stack.getPixels(13 * num_scans + ns +1); //
				for (int i = 0; i < nTiles; i++){
					scans[ns * NUM_ALL_SLICES + 0][i] = fset[0][i]; // disparity
					scans[ns * NUM_ALL_SLICES + 1][i] = fset[1][i]; // strength
				}
				fset[0] = (float[]) clt_mismatches_stack.getPixels((2 * pair + 0) * num_scans + ns +1);
				fset[1] = (float[]) clt_mismatches_stack.getPixels((2 * pair + 1) * num_scans + ns +1); //
				fset[2] = (float[]) clt_mismatches_stack.getPixels((8 + pair   ) *  num_scans + ns +1);
				for (int i = 0; i < nTiles; i++){
					scans[ns * NUM_ALL_SLICES + pair * 3 + 2][i] = fset[0][i]; // dx_i
					scans[ns * NUM_ALL_SLICES + pair * 3 + 3][i] = fset[1][i]; // dy_i
					scans[ns * NUM_ALL_SLICES + pair * 3 + 4][i] = fset[2][i]; // str_i
				}
			}
		}
		if (debugLevel > -1) {
			String [] prefixes = {"disparity", "strength", "dx0", "dy0", "str0", "dx1", "dy1", "str1", "dx2", "dy2", "str2", "dx3", "dy3", "str3"};
			String [] titles = new String [num_scans * NUM_ALL_SLICES];
			for (int ns = 0; ns < num_scans; ns++){
				for (int i = 0; i < NUM_ALL_SLICES; i++){
					titles[ns * NUM_ALL_SLICES + i] = prefixes[i]+"_"+ns;
				}
			}
			(new ShowDoubleFloatArrays()).showArrays(scans, tilesX, tilesY, true, "scans" , titles);
		}
		return scans;
	}

	public double [][][] lazyEyeCorrection(
			final boolean    use_poly, // Use polynomial correction, false - correct tilt/azimuth/roll of each sensor
			final boolean    restore_disp_inf, // Restore subtracted disparity for scan #0 (infinity)
			final double     fcorr_radius,
			final double     min_strength_in,
			final double     max_diff,
			//				final double comp_strength_var,
			final int        max_iterations,
			final double     max_coeff_diff,
			final double     far_pull, //  = 0.2; // 1; //  0.5;
			final double     strength_pow,
			final double     lazyEyeCompDiff, // clt_parameters.fcorr_disp_diff
			final int        lazyEyeSmplSide, //        = 2;      // Sample size (side of a square)
			final int        lazyEyeSmplNum, //         = 3;      // Number after removing worst (should be >1)
			final double     lazyEyeSmplRms, //         = 0.1;    // Maximal RMS of the remaining tiles in a sample
			final double     lazyEyeDispVariation, // maximal full disparity difference between the tile and 8 neighborxs
			final double     lazyEyeDispRelVariation,
			final double     ly_norm_disp, //  =    5.0;     // Reduce weight of higher disparity tiles
			final int        smplSide, //        = 2;      // Sample size (side of a square)
			final int        smplNum, //         = 3;      // Number after removing worst (should be >1)
			final double     smplRms, //         = 0.1;    // Maximal RMS of the remaining tiles in a sample
			// histogram parameters
			final int        hist_smpl_side, // 8 x8 masked, 16x16 sampled
			final double     hist_disp_min,
			final double     hist_disp_step,
			final int        hist_num_bins,
			final double     hist_sigma,
			final double     hist_max_diff,
			final int        hist_min_samples,
			final boolean    hist_norm_center, // if there are more tiles that fit than min_samples, replace with
			final double     inf_fraction,    // fraction of the weight for the infinity tiles

			final int        min_per_quadrant, // minimal tiles per quadrant (not counting the worst) tp proceed
			final int        min_inf,          // minimal number of tiles at infinity to proceed
			final int        min_inf_to_scale, // minimal number of tiles at infinity to apply weight scaling

			final boolean    right_left, // equalize weights of right/left FoV (use with horizon in both halves and gross infinity correction)
			CLTParameters    clt_parameters,
			double [][]      scans_14,
			double [][]      target_disparity, // null or programmed disparity (1 per each 14 entries of scans_14)
			int              tilesX,
			double           magic_coeff, // still not understood coefficient that reduces reported disparity value.  Seems to be around 8.5
			int              debugLevel){

//		final double lazyEyeDispRelVariation = 0.02;

		final int dbg_nTile = -34145; // 37005; // -59038;
		final int num_scans = scans_14.length/NUM_ALL_SLICES;
		final int num_tiles = scans_14[0].length;
		final int tilesY = num_tiles/tilesX;
		final boolean []  center_mask = getCenterMask(fcorr_radius, tilesX, tilesY);
		final double [][] scans = new double [num_scans * NUM_SLICES][];
		final double [][] comp_strength_rms = new double [num_scans][num_tiles];
		for (int ns = 0; ns < num_scans; ns++){
			final double [] min_weights = new double [num_tiles];
			for (int nTile = 0; nTile < num_tiles; nTile++){
				if (nTile == dbg_nTile) {
					System.out.println("lazyEyeCorrection(), nTile="+nTile);
				}
				double w = scans_14[ns * NUM_ALL_SLICES + INDEX_14_WEIGHT][nTile];
				for (int i = 0; i < INDICES_14_WEIGHTS.length; i++) {
					w = Math.min(w, scans_14[ns * NUM_ALL_SLICES + INDICES_14_WEIGHTS[i]][nTile]);
				}
				min_weights[nTile] = w;
			}
			for (int i = 0; i < INDICES_14_10.length; i++){
				if (i == INDEX_10_WEIGHT) {
					scans[ns * NUM_SLICES + i] = min_weights;
				} else {
					scans[ns * NUM_SLICES + i] = scans_14[ns * NUM_ALL_SLICES + INDICES_14_10[i]];
				}
			}
		}
		if (debugLevel > -2) { // -2) { //  100) {
			if (debugLevel > -3) { // -1) { // -2) { //  100) {
				(new ShowDoubleFloatArrays()).showArrays(scans, tilesX, tilesY, true, "scans_pre-disp");
			}
		}

		// Add disparity to dx0, dx1, dy2, dy3 pairs
		if (  restore_disp_inf) { // ==clt_parameters.inf_restore_disp
			for (int nTile = 0; nTile < num_tiles; nTile++) if (scans[INDEX_10_WEIGHT][nTile] > 0){
				for (int i = 0; i < INDICES_10_DISP.length; i++) {
					scans[INDICES_10_DISP[i]][nTile] += scans[INDEX_10_DISPARITY][nTile];
				}
			}
		}

		if (debugLevel > -2) { // -2) { //  100) {
			if (debugLevel > -3) { // -1) { // -2) { //  100) {
				(new ShowDoubleFloatArrays()).showArrays(scans, tilesX, tilesY, true, "scans_post-disp");
			}
		}
//	static final int [] INDICES_10_DISP = {2,4,7,9}; // which indices need to add disparity to restore full offsets
		//INDEX_10_WEIGHT

		for (int ns = 0; ns < num_scans; ns++){
			for (int nTile = 0; nTile < num_tiles; nTile++){
				double s1=0.0, s2=0.0;
				for (int pair = 0; pair <4; pair++){
					double s = scans_14[ns * NUM_ALL_SLICES + 4 + 3 * pair][nTile];
					s1 += s;
					s2 +=s * s;
				}
				s1 /= 4;
				s2 /= 4;
				comp_strength_rms[ns][nTile] = Math.sqrt(s2 - s1*s1);
			}
		}
		/*
		 * None of comp_strength_rms methods works to detect potential outliers for horizontal/vertical features
		 */

		if (debugLevel > 0) {
			String [] titles = new String [num_scans];
			for (int ns = 0; ns < num_scans; ns++){
				titles[ns] = "scan_" + ns;
			}
			(new ShowDoubleFloatArrays()).showArrays(comp_strength_rms, tilesX, tilesY, true, "comp_strength_rms" , titles);
		}

// FIXME: Seems that disparity should be combined with dxy for BG scan before that
		double[][] filtered_scans = filterDisparityStrength (
				scans, // final double[][] disp_strength_in, // [1][37006] >0, [21][37006] = NaN
				min_strength_in, // final double     strength_floor,
				strength_pow, // final double     strength_pow,
				lazyEyeSmplSide, // final int        smplSide, //        = 2;      // Sample size (side of a square)
				lazyEyeSmplNum, // final int        smplNum, //         = 3;      // Number after removing worst (should be >1)
				lazyEyeSmplRms, // final double     smplRms, //         = 0.1;    // Maximal RMS of the remaining tiles in a sample
				tilesX);// final int        tilesX);

		if (debugLevel > -2) { // -2) { //  100) {
			(new ShowDoubleFloatArrays()).showArrays(filtered_scans, tilesX, tilesY, true, "filtered_scans");
			if (debugLevel > -3) { // -1) { // -2) { //  100) {
				(new ShowDoubleFloatArrays()).showArrays(scans, tilesX, tilesY, true, "scans");
				(new ShowDoubleFloatArrays()).showArrays(target_disparity, tilesX, tilesY, true, "target_disparity");
			}
		}

		if (debugLevel > -2) {
			System.out.println("lazyEyeCorrection().a 1: removing tiles with residual disparity absolute value > "+lazyEyeCompDiff);
//			System.out.println("lazyEyeCorrection().a 1: removing tiles with residual disparity absolute value > "+(2*lazyEyeCompDiff));
		}

		double [][] combo_mismatch = new double [NUM_SLICES][num_tiles];
		double [] combo_comp_rms = new double [num_tiles];
		for (int ns = 0; ns < num_scans; ns++){
			for (int nTile = 0; nTile < num_tiles; nTile++) {
				if (nTile == dbg_nTile){
					System.out.println("lazyEyeCorrection().1: nTile="+nTile); // filtered_scans[2][37005] = NaN
				}
				double w = filtered_scans[ns * NUM_SLICES + 1][nTile];
				if (w > 0.0){
					double disp = filtered_scans[ns * NUM_SLICES + 0][nTile];
//					if (Math.abs(disp) <= 2.0*lazyEyeCompDiff) {
					if (Math.abs(disp) <= lazyEyeCompDiff) {
						for (int i = 2; i < NUM_SLICES; i++) if (i != 1){
							combo_mismatch[i][nTile] += filtered_scans[ns * NUM_SLICES + i][nTile] * w;
						}
						//FIXME: ???? target_disparity is not 0 for bg
						// combo_mismatch combines both infinity and regular for the same tile, mixing "disparity" and "target disparity" with weights and magic_scale
						// Seems to be wrong, as target_disparity is only estimated disparity, not measured. Or is it measured for non-infinity?
						// At least bg scan is measured with disparity =0, even as target_disparity is not 0
						// combo data is later used as a non-infinity to correct all but disparity

						if (target_disparity != null){
							combo_mismatch[0][nTile] += (
									filtered_scans[ns * NUM_SLICES + 0][nTile]/clt_parameters.corr_magic_scale +
//									target_disparity[ns][nTile])* w;
									// disabling for bg - already combined
							((ns==0)?0.0:target_disparity[ns][nTile])* w);
						} else {
							combo_mismatch[0][nTile] += (
									filtered_scans[ns * NUM_SLICES + 0][nTile]/clt_parameters.corr_magic_scale +
									(clt_parameters.disp_scan_start + clt_parameters.disp_scan_step * ns ) )* w;
						}
						combo_mismatch[1][nTile] += w;
						combo_comp_rms[nTile] += w * comp_strength_rms[ns][nTile];
					}
				}
			}
		}

		for (int nTile = 0; nTile < num_tiles; nTile++) {
			if (nTile == dbg_nTile){
				System.out.println("lazyEyeCorrection().2: nTile="+nTile);
			}

			double w = combo_mismatch[1][nTile];
			if (w > 0.0){
				for (int i = 0; i < NUM_SLICES; i++) if (i != 1){
					combo_mismatch[i][nTile] /= w;
				}
				combo_comp_rms[nTile] /= w;
			} else {
				for (int i = 0; i < NUM_SLICES; i++) if (i != 1){
					combo_mismatch[i][nTile] = Double.NaN;
				}
				combo_comp_rms[nTile] = Double.NaN;
			}
		}

		// reduce influence of high disparity,  using combined disparity
//		double norm_ly_disparity = 100.0; // disabling
		for (int nTile = 0; nTile < num_tiles; nTile++) {
			if ((combo_mismatch[0][nTile] > 0) && (combo_mismatch[0][nTile] > ly_norm_disp)) {
				combo_mismatch[1][nTile] *= ly_norm_disp/combo_mismatch[0][nTile];
			}
		}

		if (debugLevel > 0) { // -2) { //  100) {
			(new ShowDoubleFloatArrays()).showArrays(combo_comp_rms, tilesX, tilesY, "combo_comp_rms");
		}

		// instance of class to operate navigation over tiles
		// compare tile disparity (combo) with those of neighbors, discard if too different
		final TileNeibs tnImage = new TileNeibs(tilesX, tilesY); // num_tiles/tilesX);
		for (int nTile = 0; nTile < num_tiles; nTile++) if (combo_mismatch[1][nTile] > 0.0){
			if (nTile == dbg_nTile){
				System.out.println("lazyEyeCorrection().3: nTile="+nTile);
			}
			double d = combo_mismatch[0][nTile];
			double lev = lazyEyeDispVariation + lazyEyeDispRelVariation * d;
			for (int dir = 0; dir <8; dir++){
				int nTile1 = tnImage.getNeibIndex(nTile, dir);
				if ((nTile1 >= 0) && (combo_mismatch[1][nTile1] > 0.0)){
					if (Math.abs(combo_mismatch[0][nTile1] - d) > lev) { // azyEyeDispVariation){
						combo_mismatch[1][nTile] = 0.0;
						combo_mismatch[0][nTile] = Double.NaN;
						for (int i = 2; i < NUM_SLICES; i++) if (i != 1){
							combo_mismatch[i][nTile] = Double.NaN;
						}
						break;
					}
				}
			}
		}


		// here combo_mismatch[2][37005] = Double.NaN,combo_mismatch[1][37005] != 0.0, combo_mismatch[0][37005] = 0.0
		if (debugLevel > 0) { // 0) {
			String [] prefixes = {"disparity", "strength", "dx0", "dy0", "dx1", "dy1", "dx2", "dy2", "dx3", "dy3"};
			(new ShowDoubleFloatArrays()).showArrays(combo_mismatch, tilesX, combo_mismatch[0].length/tilesX, true, "combo_mismatch" , prefixes);
		}
		if (clt_parameters.lyf_filter) {
			combo_mismatch =  filterLazyEyePairs (
					combo_mismatch, // final double[][] samples_in,
					clt_parameters.lyf_smpl_side ,   // 8,              // final int        smpl_side, // 8 x8 masked, 16x16 sampled
					clt_parameters.lyf_rms_max ,     // 0.25,           // final double     rms_max, TODO: find reasonable one not critical?
					clt_parameters.lyf_frac_keep ,   // 0.5,            // final double     frac_keep,
					clt_parameters.lyf_min_samples , // 5,              // final int        min_samples,
					clt_parameters.lyf_norm_center , // true,           // final boolean    norm_center, // if there are more tiles that fit than minsamples, replace with a single equal weight
					tilesX);        // final int        tilesX);
		}
		if (debugLevel > 0) {
			String [] prefixes = {"disparity", "strength", "dx0", "dy0", "dx1", "dy1", "dx2", "dy2", "dx3", "dy3"};
			(new ShowDoubleFloatArrays()).showArrays(combo_mismatch, tilesX, combo_mismatch[0].length/tilesX, true, "filtered_mismatch" , prefixes);
		}



		// extract infinity data to be processed as infinity

		double [][] inf_scan = new double [NUM_SLICES][];
		for (int i = 0; i < NUM_SLICES; i++){
			inf_scan[i] = scans[i]; // just copy first half (infinity data - now it has deltas x/y with restored disparity (no magic_scale)
		}

		// Optionally filter infinity scan data
		if (smplSide > 1){
			inf_scan = filterDisparityStrength (
					inf_scan,
					min_strength_in,    // final double     strength_floor,
					strength_pow,
					smplSide, //        = 2;      // Sample size (side of a square)
					smplNum, //         = 3;      // Number after removing worst (should be >1)
					smplRms, //         = 0.1;    // Maximal RMS of the remaining tiles in a sample
					tilesX);

			if (debugLevel > 0){
				double [][] dbg_img = inf_scan.clone();
				for (int n = 0; n < inf_scan.length; n++){
					dbg_img[n] = inf_scan[n].clone();
				}
				for (int n = 0; n < dbg_img.length; n+=2){
					for (int i = 0; i < dbg_img[n].length; i++) {
						if (dbg_img[n+1][i] == 0.0){
							dbg_img[n][i] = Double.NaN;
						}
					}
				}

				(new ShowDoubleFloatArrays()).showArrays(dbg_img, tilesX, tilesY, true, "filtered_infinity_ds"); // , titles);

			}
		}

		// Optionally filter infinity scan data with a histogram filter

		if (hist_smpl_side > 0) { // 0 to bypass histogram filtering
			inf_scan = filterHistogramFar (
					inf_scan,         // final double[][] disp_strength_in,
					hist_smpl_side,   // final int        smpl_side, // 8 x8 masked, 16x16 sampled
					hist_disp_min,    // final double     disp_min,
					hist_disp_step,   // final double     disp_step,
					hist_num_bins,    // final int        num_bins,
					hist_sigma,       // final double     sigma,
					hist_max_diff,    // final double     max_diff,
					hist_min_samples, // final int        min_samples,
					hist_norm_center, // final boolean    norm_center, // if there are more tiles that fit than minsamples, replace with
					tilesX);          // final int        tilesX)
			if (debugLevel > 1){
				double [][] dbg_img = inf_scan.clone();
				for (int n = 0; n < inf_scan.length; n++){
					dbg_img[n] = inf_scan[n].clone();
				}
				for (int n = 0; n < dbg_img.length; n+=2){
					for (int i = 0; i < dbg_img[n].length; i++) {
						if (dbg_img[n+1][i] == 0.0){
							dbg_img[n][i] = Double.NaN;
						}
					}
				}
				(new ShowDoubleFloatArrays()).showArrays(dbg_img, tilesX, tilesY, true, "hist_filt_ds"); // , titles);
			}
		}
		// combine infinity and lazy eye scan data into a single array
		double [][] inf_and_ly = new double [2 * NUM_SLICES][];
		for (int i = 0; i < NUM_SLICES; i++){
			inf_and_ly[i] = inf_scan[i];
			inf_and_ly[i + NUM_SLICES] = combo_mismatch[i];
		}

// make all zero strength tiles to have NaN values to use histrograms in ImageJ
		for (int ns = 0; ns < 2; ns++) {
			for (int nt = 0; nt < inf_and_ly[INDEX_10_WEIGHT + ns * NUM_SLICES].length; nt++ ) {
				if (inf_and_ly[INDEX_10_WEIGHT +  ns * NUM_SLICES][nt] == 0.0) {
					for (int i = 0; i < NUM_SLICES; i++) if (i != INDEX_10_WEIGHT){
						inf_and_ly[i +  ns * NUM_SLICES][nt] = Double.NaN;
					}
				}

			}
		}

//	static final int    INDEX_10_WEIGHT = 1;
		if (debugLevel > 0) {
			System.out.println("test123");
		}
		if ((debugLevel > -1) && (hist_smpl_side > 0)) { // 0) {
			String [] prefixes = {"disparity", "strength", "dx0", "dy0", "dx1", "dy1", "dx2", "dy2", "dx3", "dy3"};
			String [] titles = new String [2 * NUM_SLICES];
			for (int i = 0; i < NUM_SLICES; i++){
				titles[i] = prefixes[i]+"-inf";
				titles[i + NUM_SLICES] = prefixes[i]+"-ly";
			}
			(new ShowDoubleFloatArrays()).showArrays(inf_and_ly, tilesX, tilesY, true, "inf_and_ly",titles);

			int step = hist_smpl_side; // should be the same for both filters
			int tilesX1 = tilesX/step;
			int tilesY1 = tilesY/step;
			int num_tiles1 = tilesX1 * tilesY1;
			double [][] dbg_img = new double [inf_and_ly.length][num_tiles1];
			for (int tY = 0; tY < tilesY1; tY++) {
				for (int tX = 0; tX < tilesX1; tX++) {
					//		    			if (tY == 14) {
					//		    				System.out.println("lazyEyeCorrection(): tX="+tX+", tY="+tY);
					//		    			}
					int nTile1 = tX + tY*tilesX1;
					for (int sY = 0; sY < step; sY ++) {
						for (int sX = 0; sX < step; sX ++) {
							int nTile = (sX + step * tX) + (sY + step * tY) * tilesX;
							for (int ns = 0; ns <2; ns++){
								double w = inf_and_ly[ns*NUM_SLICES + 1][nTile];
								if (w > 0.0){
									for (int i = 0; i < NUM_SLICES; i++) if (i != 1) {
										dbg_img[ns * NUM_SLICES + i][nTile1] += w * inf_and_ly[ns * NUM_SLICES + i][nTile];
									}
									dbg_img[ns * NUM_SLICES + 1][nTile1] += w;
								}
							}
						}
					}

					for (int ns = 0; ns < 2; ns++){
						double w = dbg_img[ns * NUM_SLICES + 1][nTile1];
						if (w > 0.0){
							for (int i = 0; i < NUM_SLICES; i++) if (i != 1) {
								dbg_img[ns * NUM_SLICES + i][nTile1] /= w;
							}
						} else {
							for (int i = 0; i < NUM_SLICES; i++) if (i != 1) {
								dbg_img[ns * NUM_SLICES + i][nTile1] = Double.NaN;
							}
						}
					}
				}
			}
			(new ShowDoubleFloatArrays()).showArrays(dbg_img, tilesX1, tilesY1, true, "inf_and_ly8",titles);
		}
		if (debugLevel > 0) {
			System.out.println("test1234");
		}
		// create list for infinity data
//		/clt_parameters.ly_inf_en,
		ArrayList<Sample> inf_samples_list;
//		if (clt_parameters.ly_inf_en) {
		if (true) { //clt_parameters.ly_inf_en) {
			inf_samples_list = selectInfinityTiles(
					clt_parameters.fcorr_radius,       // 	final double fcorr_radius,
					clt_parameters.fcorr_inf_vert,// final boolean use_vertical,
					0.0, // any > 0.0
					max_diff, // max_diff, //clt_parameters.fcorr_inf_diff
					max_iterations, // max_iterations, // clt_parameters.inf_iters
					max_coeff_diff, // max_coeff_diff, // clt_parameters.inf_final_diff
					far_pull, // far_pull, // clt_parameters.inf_far_pull,  = 0.2; // 1; //  0.5;
					clt_parameters,
					inf_scan,
					tilesX,
					magic_coeff, // magic_coeff, // still not understood coefficient that reduces reported disparity value.  Seems to be around 8.5
					debugLevel);

			if (debugLevel > -1) {
				double inf_weight = 0.0;
				for (Sample s: inf_samples_list) {
					inf_weight += s.weight;
				}
				System.out.println("lazyEyeCorrection(): number of infinity samples="+inf_samples_list.size()+", total weight = "+inf_weight);

			}
			// adjust weight to balance infinity data and lazy eye one. As some tiles were discarded by selectInfinityTiles() list and not the original
			// array has to be used to find the total weight of the infinity tile. Other ones will be used with no extra filtering
			double [] total_weights = new double[2];
			int num_inf = 0;
			int [] per_quad = new int[4];
			for (Sample s: inf_samples_list) {
				total_weights[0] += s.weight;
				if (s.weight > 0.0) {
					num_inf ++;
					int hx = (s.tile % tilesX) / (tilesX/2);
					int hy = (s.tile / tilesX) / (tilesY/2);
					per_quad[hx + 2*hy]++;
				}
			}

			for (int nTile = 0; nTile < num_tiles; nTile++) if (center_mask[nTile]){ // calculate total weight of non-infinity
				total_weights[1]+= inf_and_ly[1 * NUM_SLICES + 1][nTile];
				if (inf_and_ly[1 * NUM_SLICES + 1][nTile] > 0.0) {
					int hx = (nTile % tilesX) / (tilesX/2);
					int hy = (nTile / tilesX) / (tilesY/2);
					per_quad[hx + 2*hy]++;
				}
			}

			int [] pq = per_quad.clone();
			Arrays.sort(per_quad);
			if (debugLevel > -20) {
				System.out.println(String.format("Tiles per quadrants :[%d, %d, %d, %d], tiles at infinity %d", pq[0],pq[1],pq[2],pq[3],num_inf));
			}
			if (per_quad[1] < min_per_quadrant) {
				if (debugLevel > -20) {
					System.out.println(String.format("Too few tiles in quadrants :[%d, %d, %d, %d], minimum for the second worst is %d", pq[0],pq[1],pq[2],pq[3],min_per_quadrant));
				}
				return null;
			}

			if (num_inf < min_inf) {
				if (debugLevel > -20) {
					System.out.println(String.format("Too few tiles at infinity: %d minimum is %d", num_inf, min_inf));
				}
				return null;
			}

			double inf_fraction_limited =  (inf_fraction >= 0.0) ?((inf_fraction > 1.0) ? 1.0 : inf_fraction):0.0;

			double [] weights = {
					inf_fraction_limited *         (total_weights[0] + total_weights[1]) / total_weights[0],
					(1.0 - inf_fraction_limited) * (total_weights[0] + total_weights[1]) / total_weights[1],
			};

			if (num_inf < min_inf_to_scale) {
				if (debugLevel>-1) {
					System.out.println("Too few infinity tiles to boost ("+num_inf+" < "+min_inf_to_scale+", keeping original weights");
				}
			} else if (weights[0] > weights[1]) {
				if (debugLevel>-1) {
					System.out.println("Boosting weights of far tiles (weights[0]="+weights[0]+", weights[1]="+weights[1]);
				}

				for (int ns = 0; ns <2; ns++) {
					for (int nTile = 0; nTile < num_tiles; nTile++) {
						inf_and_ly[ns * NUM_SLICES + 1][nTile] *= weights[ns];
					}
				}
				for (Sample s: inf_samples_list) {
					s.weight *= weights[0];
				}

			} else {
				if (debugLevel>-1) {
					System.out.println("There are already more far tiles than requested (weights[0]="+weights[0]+", weights[1]="+weights[1]+", so keeping original weights");
				}
			}
		}
		///-----


		// Supplement list with the lazy eye scans data - use all tiles
		for (int nTile = 0; nTile < num_tiles; nTile++) if (center_mask[nTile]) {
			double w = inf_and_ly[1 * NUM_SLICES + 1][nTile];
			if (w > 0.0) {
				inf_samples_list.add(new Sample(1,nTile,w));
			}
		}

		if (debugLevel > 0) {
			double inf_weight = 0.0;
			for (Sample s: inf_samples_list) {
				inf_weight += s.weight;
			}
			System.out.println("lazyEyeCorrection(): number of all samples="+inf_samples_list.size()+", total weight = "+inf_weight);

		}
		if (right_left) {
			System.out.println("Balancing right/left part of FoV weights, width = "+tilesX+" tiles");
			double [] right_left_weights = {0.0, 0.0};
			for (Sample s: inf_samples_list) {
				int rl = (s.tile % tilesX) / (tilesX/2);
				right_left_weights[rl] += s.weight;
			}
			System.out.println("Weights: left:"+ right_left_weights[0]+", right:"+ right_left_weights[1]);
			for (Sample s: inf_samples_list) {
				int rl = (s.tile % tilesX) / (tilesX/2);
				s.weight *= right_left_weights[1 - rl];
			}
			// just verifying
			right_left_weights[0] = 0.0;
			right_left_weights[1] = 0.0;
			for (Sample s: inf_samples_list) {
				int rl = (s.tile % tilesX) / (tilesX/2);
				right_left_weights[rl] += s.weight;
			}
			System.out.println("Weights after balancing: left:"+ right_left_weights[0]+", right:"+ right_left_weights[1]);
		}



		if (debugLevel > 0) {
			String [] prefixes = {"disparity", "strength", "dx0", "dy0", "dx1", "dy1", "dx2", "dy2", "dx3", "dy3"};
			String [] titles = new String [num_scans * NUM_SLICES];
			for (int ns = 0; ns < num_scans; ns++){
				for (int i = 0; i < NUM_SLICES; i++){
					titles[ns * NUM_SLICES + i] = prefixes[i]+"_"+ns;
				}
			}
			(new ShowDoubleFloatArrays()).showArrays(filtered_scans, tilesX, tilesY, true, "filtered_scans_a" , titles);
		}


		if (debugLevel > -1) {
			String [] prefixes = {"disparity", "strength", "dx0", "dy0", "dx1", "dy1", "dx2", "dy2", "dx3", "dy3"};
			(new ShowDoubleFloatArrays()).showArrays(combo_mismatch, tilesX, combo_mismatch[0].length/tilesX, true, "combo_mismatch" , prefixes);
		}

		if (debugLevel > -1) { // 0) {
			String [] prefixes = {"disparity", "strength", "dx0", "dy0", "dx1", "dy1", "dx2", "dy2", "dx3", "dy3"};
			String [] titles = new String [2 * NUM_SLICES];
			for (int i = 0; i < NUM_SLICES; i++){
				titles[i] = prefixes[i]+"-inf";
				titles[i + NUM_SLICES] = prefixes[i]+"-ly";
			}
			(new ShowDoubleFloatArrays()).showArrays(inf_and_ly, tilesX, tilesY, true, "inf_and_ly_last",titles);
		}
		ArrayList<Mismatch> mismatch_list = use_poly? null : (new ArrayList<Mismatch>());
		// inf_and_ly here has filtered disparity and offsets, should be process clt_parameters.ly_inf_disp before filters
		// for rig with known disparity - use series = 0 - it will allow disparity adjustment
		double [][][] mismatch_corr_coefficients = infinityMismatchCorrection(
				clt_parameters.disp_scan_start, // final double  disp_scan_start,
				clt_parameters.disp_scan_step,  // final double  disp_scan_step,
				use_poly,                       // final boolean use_poly,
				clt_parameters.fcorr_quadratic, // final boolean use_quadratic,
				true, // clt_parameters.fcorr_inf_vert,  // final boolean use_vertical,
				// too late to restore disparity - should be dome earlier
				false,                          // final boolean use_disparity, // for infinity
				true, // clt_parameters.ly_inf_disp,     //final boolean allow_dispatity,
				clt_parameters,                 // EyesisCorrectionParameters.CLTParameters           clt_parameters,
				inf_and_ly,                     // double [][] disp_strength,
				inf_samples_list,               // ArrayList<Sample> samples_list,
				tilesX,                         // int         tilesX,
				magic_coeff,                    // double      , // still not understood coefficient that reduces reported disparity value.  Seems to be around 8.5
				mismatch_list,                  // ArrayList<Mismatch> mismatch_list,
				debugLevel);                    // int debugLevel)
		if (debugLevel > -2) {
			System.out.println("===== lazyEyeCorrection(): correction coefficients =====");
			if (mismatch_corr_coefficients != null) {
				show_fine_corr(
						mismatch_corr_coefficients,
						"mismatch_corr_coefficients");
			} else {
				System.out.println("Are null - non-null are for poly correction only");
			}
		}
		if (!use_poly && (mismatch_list != null)){
			double [] old_new_rms = new double[1];
			boolean apply_extrinsic = true;
			int solveCorr_debug =  ((clt_parameters.lym_iter == 1) && (clt_parameters.ly_par_sel != 0))? 2 : debugLevel;
			GeometryCorrection.CorrVector corr_vector = solveCorr (
					clt_parameters.ly_inf_en,      // boolean use_disparity,     // if true will ignore disparity data even if available (was false)
//					clt_parameters.ly_combo_en,    // boolean use_other_extr,    // adjust other extrinsic parameters that do not influence disparity, common roll and zoom
					clt_parameters.ly_aztilt_en,// boolean use_aztilts,       // Adjust azimuths and tilts excluding disparity
					clt_parameters.ly_diff_roll_en,// boolean use_diff_rolls,    // Adjust differential rolls (3 of 4 angles)
					clt_parameters.ly_inf_force,   // boolean force_convergence, // if true try to adjust convergence (disparity, symmetrical parameter 0) even with no disparity
					clt_parameters.ly_com_roll,    // boolean    common_roll,    // Enable common roll (valid for high disparity range only)
					clt_parameters.ly_focalLength, // boolean    corr_focalLength,     // Correct scales (focal length temperature? variations)
					clt_parameters.ly_par_sel,     //int     manual_par_sel,    // Manually select the parameter mask bit 0 - sym0, bit1 - sym1, ... (0 - use boolean flags, != 0 - ignore boolean flags)
					mismatch_list,                          // ArrayList<Mismatch> mismatch_list,
					qc.geometryCorrection,                  // GeometryCorrection geometryCorrection,
///					null, // GeometryCorrection geometryCorrection_main, // if is aux camera using main cameras' coordinates. Disparity is still in aux camera pixels
					qc.geometryCorrection.getCorrVector(),  // GeometryCorrection.CorrVector corr_vector,
					old_new_rms,                            // double [] old_new_rms, // should be double[2]
//					2); // debugLevel); // 2); // 1); // int debugLevel)
					solveCorr_debug); // debugLevel); // 2); // 1); // int debugLevel)
//TODO: ** Put 2 here to debug derivative images (diff_dmv_dsym - does not match yet, probably different "addition" of angles)

			if (debugLevel > -1){
				System.out.println("Old extrinsic corrections:");
				System.out.println(qc.geometryCorrection.getCorrVector().toString());
				System.out.println("Delta extrinsic corrections:");
				System.out.println(corr_vector.toString());
			}
			if (apply_extrinsic){
				boolean ok = qc.geometryCorrection.getCorrVector().incrementVector(corr_vector, clt_parameters.ly_corr_scale);
				if (!ok) {
					System.out.println("Failed to solve correction, corr_vector:"+corr_vector.toString());
					return null;
				}
				if (debugLevel > -1){
					System.out.println("New extrinsic corrections:");
					System.out.println(qc.geometryCorrection.getCorrVector().toString());
				}
			}
			mismatch_corr_coefficients = new double [1][2][];
			mismatch_corr_coefficients[0][0] = corr_vector.toSymArray(null);
//			mismatch_corr_coefficients[0][1] = new double[1];
//			old_new_rms[1] = getRMS(getYminusFx(mismatch_list), getWeights(mismatch_list));
			mismatch_corr_coefficients[0][1] = old_new_rms;
		} else {
			if (debugLevel > -2){
				System.out.println("Extrinsic parameters (tilt, azimuth, roll) of subcameras is disabled, use_poly="+
						use_poly+" (should be false for extrinsics)");
				System.out.println(qc.geometryCorrection.getCorrVector().toString());
			}
			return mismatch_corr_coefficients;
		}
		return mismatch_corr_coefficients;
	}





/**
 *
 * @param use_poly
 * @param restore_disp_inf
 * @param fcorr_radius
 * @param min_strength_in
 * @param strength_pow
 * @param lazyEyeCompDiff
 * @param lazyEyeSmplSide
 * @param lazyEyeSmplNum
 * @param lazyEyeSmplRms
 * @param lazyEyeDispVariation
 * @param lazyEyeDispRelVariation
 * @param ly_norm_disp
 * @param smplSide
 * @param smplNum
 * @param smplRms
 * @param hist_smpl_side
 * @param hist_disp_min
 * @param hist_disp_step
 * @param hist_num_bins
 * @param hist_sigma
 * @param hist_max_diff
 * @param hist_min_samples
 * @param hist_norm_center
 * @param inf_fraction
 * @param min_per_quadrant
 * @param min_inf
 * @param min_inf_to_scale
 * @param inf_max_disparity
 * @param clt_parameters
 * @param scans_14
 * @param gt_disparity_strength
 * @param filter_ds
 * @param filter_lyf
 * @param tilesX
 * @param magic_coeff
 * @param debugLevel
 * @return will return null if can not adjust
 */







	public double [][][] lazyEyeCorrectionFromGT(
//			final GeometryCorrection geometryCorrection_main, // if not null - this is an AUX camera of a rig
			final boolean    use_poly, // Use polynomial correction, false - correct tilt/azimuth/roll of each sensor
			final boolean    restore_disp_inf, // Restore subtracted disparity for scan #0 (infinity) always true
			final double     fcorr_radius,
			final double     min_strength_in,

			final double     strength_pow,
			final double     lazyEyeCompDiff, // clt_parameters.fcorr_disp_diff
			final int        lazyEyeSmplSide, //        = 2;      // Sample size (side of a square)
			final int        lazyEyeSmplNum, //         = 3;      // Number after removing worst (should be >1)
			final double     lazyEyeSmplRms, //         = 0.1;    // Maximal RMS of the remaining tiles in a sample
			final double     lazyEyeDispVariation, // maximal full disparity difference between the tile and 8 neighborxs
			final double     lazyEyeDispRelVariation,
			final double     ly_norm_disp, //  =    5.0;     // Reduce weight of higher disparity tiles
			final int        smplSide, //        = 2;      // Sample size (side of a square)
			final int        smplNum, //         = 3;      // Number after removing worst (should be >1)
			final double     smplRms, //         = 0.1;    // Maximal RMS of the remaining tiles in a sample
			// histogram parameters
			final int        hist_smpl_side, // 8 x8 masked, 16x16 sampled
			final double     hist_disp_min,
			final double     hist_disp_step,
			final int        hist_num_bins,
			final double     hist_sigma,
			final double     hist_max_diff,
			final int        hist_min_samples,
			final boolean    hist_norm_center, // if there are more tiles that fit than min_samples, replace with
			final double     inf_fraction,     // fraction of the weight for the infinity tiles
			final int        min_per_quadrant, // minimal tiles per quadrant (not counting the worst) tp proceed
			final int        min_inf,          // minimal number of tiles at infinity to proceed
			final int        min_inf_to_scale, // minimal number of tiles at infinity to apply weight scaling
			final double     inf_max_disparity, // use all smaller disparities as inf_fraction
			CLTParameters    clt_parameters,
			double [][]      scans_14, // here - always 14 - infinity and non-infinity
			double [][][]    gt_disparity_strength, // 1 pair for each 14 entries of scans_14 (normally - just 1 scan
			final boolean    filter_ds, //
			final boolean    filter_lyf, // ~clt_parameters.lyf_filter, but may be different, now off for a single cameras
			int              tilesX,
			double           magic_coeff, // still not understood coefficient that reduces reported disparity value.  Seems to be around 8.5
			int              debugLevel){

//		final double lazyEyeDispRelVariation = 0.02;
		final int dbg_nTile = -34145; // 37005; // -59038;
		final int num_scans = scans_14.length/NUM_ALL_SLICES;
		final int num_tiles = scans_14[0].length;
		final int tilesY = num_tiles/tilesX;
		final boolean []  center_mask = getCenterMask(fcorr_radius, tilesX, tilesY);
		final double [][] scans = new double [num_scans * NUM_SLICES][];
//		final double [][] comp_strength_rms = new double [num_scans][num_tiles];


		for (int ns = 0; ns < num_scans; ns++){
			final double [] min_weights = new double [num_tiles];
			for (int nTile = 0; nTile < num_tiles; nTile++){
				if (nTile == dbg_nTile) {
					System.out.println("lazyEyeCorrectionFromGT(), nTile="+nTile);
				}
				double w = scans_14[ns * NUM_ALL_SLICES + INDEX_14_WEIGHT][nTile];
				for (int i = 0; i < INDICES_14_WEIGHTS.length; i++) {
					w = Math.min(w, scans_14[ns * NUM_ALL_SLICES + INDICES_14_WEIGHTS[i]][nTile]);
				}
				min_weights[nTile] = w * gt_disparity_strength[ns][1][nTile];
			}
			for (int i = 0; i < INDICES_14_10.length; i++){
				if (i == INDEX_10_WEIGHT) {
					scans[ns * NUM_SLICES + i] = min_weights;
				} else {
					scans[ns * NUM_SLICES + i] = scans_14[ns * NUM_ALL_SLICES + INDICES_14_10[i]];
				}
			}
		}
		if (debugLevel > 0) { // -2) { //  100) {
			if (debugLevel > -3) { // -1) { // -2) { //  100) {
				double [][] dbg_scans = new double[scans.length][];
				for (int i = 0; i < dbg_scans.length;i++) {
					dbg_scans[i] = scans[i].clone();
					if (i != 1) {
						for (int j = 0; j < dbg_scans[i].length; j++) if (scans[1][j]<=0.0) {
							dbg_scans[i][j] = Double.NaN;
						}
					}
				}
				(new ShowDoubleFloatArrays()).showArrays(dbg_scans, tilesX, tilesY, true, "scans_pre-disp");
			}
		}

		// Add disparity to dx0, dx1, dy2, dy3 pairs (here - always)
		if (  restore_disp_inf) { //  && false) { // ==clt_parameters.inf_restore_disp
			for (int nTile = 0; nTile < num_tiles; nTile++) if (scans[INDEX_10_WEIGHT][nTile] > 0){
				for (int i = 0; i < INDICES_10_DISP.length; i++) {
					scans[INDICES_10_DISP[i]][nTile] += scans[INDEX_10_DISPARITY][nTile];
				}
			}
		}

		if (debugLevel > 0) { // -2) { //  100) {
			if (debugLevel > -3) { // -1) { // -2) { //  100) {
				double [][] dbg_scans = new double[scans.length][];
				for (int i = 0; i < dbg_scans.length;i++) {
					dbg_scans[i] = scans[i].clone();
					if (i != 1) {
						for (int j = 0; j < dbg_scans[i].length; j++) if (scans[1][j]<=0.0) {
							dbg_scans[i][j] = Double.NaN;
						}
					}
				}
				(new ShowDoubleFloatArrays()).showArrays(dbg_scans, tilesX, tilesY, true, "scans_post-disp");
			}
		}


// FIXME: Seems that disparity should be combined with dxy for BG scan before that
		// For GT - keep it here or remove?
		double[][] filtered_scans;
		if (filter_ds) {
			filtered_scans = filterDisparityStrength (
					scans, // final double[][] disp_strength_in, // [1][37006] >0, [21][37006] = NaN
					min_strength_in, // final double     strength_floor,
					strength_pow, // final double     strength_pow,
					lazyEyeSmplSide, // final int        smplSide, //        = 2;      // Sample size (side of a square)
					lazyEyeSmplNum, // final int        smplNum, //         = 3;      // Number after removing worst (should be >1)
					lazyEyeSmplRms, // final double     smplRms, //         = 0.1;    // Maximal RMS of the remaining tiles in a sample
					tilesX);// final int        tilesX);
		} else {
			filtered_scans = scans;
		}
		if (debugLevel > 0) { // -2) { //  100) {
			{
				double [][] dbg_scans = new double[scans.length][];
				for (int i = 0; i < dbg_scans.length;i++) {
					dbg_scans[i] = filtered_scans[i].clone();
					if (i != 1) {
						for (int j = 0; j < dbg_scans[i].length; j++) if (filtered_scans[1][j]<=0.0) {
							dbg_scans[i][j] = Double.NaN;
						}
					}
				}
				(new ShowDoubleFloatArrays()).showArrays(dbg_scans, tilesX, tilesY, true, "filtered_scans");
			}
			if (debugLevel > -3) { // -1) { // -2) { //  100) {
				double [][] dbg_scans = new double[scans.length][];
				for (int i = 0; i < dbg_scans.length;i++) {
					dbg_scans[i] = scans[i].clone();
					if (i != 1) {
						for (int j = 0; j < dbg_scans[i].length; j++) if (scans[1][j]<=0.0) {
							dbg_scans[i][j] = Double.NaN;
						}
					}
				}
				(new ShowDoubleFloatArrays()).showArrays(dbg_scans, tilesX, tilesY, true, "scans-after");
				(new ShowDoubleFloatArrays()).showArrays(gt_disparity_strength[0], tilesX, tilesY, true, "gt_disparity_strength");
			}
		}

		if (debugLevel > -2) {
			System.out.println("lazyEyeCorrectionFromGT() 1: removing tiles with residual disparity absoulte value > "+lazyEyeCompDiff);
		}

		double [][] combo_mismatch = new double [NUM_SLICES][num_tiles];
		for (int ns = 0; ns < num_scans; ns++){
			for (int nTile = 0; nTile < num_tiles; nTile++) {
				if (nTile == dbg_nTile) { // || (nTile == 24971)){
					System.out.println("lazyEyeCorrectionFromGT().1: nTile="+nTile); // filtered_scans[2][37005] = NaN
				}
//				double w = filtered_scans[ns * NUM_SLICES + 1][nTile];
//				if ((w > 0.0) && (gt_disparity_strength[ns][1][nTile] > 0.0)){ // filtered strength may be non-zero where gt_disparity_strength[ns][1][nTile] is -> NaN
				// reversing - use GT strength, but skip if there is no filtered?
				double w = gt_disparity_strength[ns][1][nTile]; // GT data
				if ((w > 0.0) && (filtered_scans[ns * NUM_SLICES + 1][nTile] > 0.0)){ // filtered strength may be non-zero where gt_disparity_strength[ns][1][nTile] is -> NaN
					double disp = filtered_scans[ns * NUM_SLICES + 0][nTile];
					if (Math.abs(disp) <= lazyEyeCompDiff) {
						for (int i = 2; i < NUM_SLICES; i++) if (i != 1){
							combo_mismatch[i][nTile] += filtered_scans[ns * NUM_SLICES + i][nTile] * w;
						}
						//FIXME: ???? target_disparity is not 0 for bg
						// combo_mismatch combines both infinity and regular for the same tile, mixing "disparity" and "target disparity" with weights and magic_scale
						// Seems to be wrong, as target_disparity is only estimated disparity, not measured. Or is it measured for non-infinity?
						// At least bg scan is measured with disparity =0, even as target_disparity is not 0
						// combo data is later used as a non-infinity to correct all but disparity

						combo_mismatch[0][nTile] += gt_disparity_strength[ns][0][nTile] * w;
						combo_mismatch[1][nTile] += w;
					}
				}
			}
		}
		if (debugLevel > 0) { // 0) {
			String [] prefixes = {"disparity", "strength", "dx0", "dy0", "dx1", "dy1", "dx2", "dy2", "dx3", "dy3"};
			(new ShowDoubleFloatArrays()).showArrays(combo_mismatch, tilesX, combo_mismatch[0].length/tilesX, true, "removed_residual"+lazyEyeCompDiff, prefixes);
		}



		for (int nTile = 0; nTile < num_tiles; nTile++) {
			if (nTile == dbg_nTile){
				System.out.println("lazyEyeCorrectionFromGT().2: nTile="+nTile);
			}

			double w = combo_mismatch[1][nTile];
			if (w > 0.0){
				for (int i = 0; i < NUM_SLICES; i++) if (i != 1){
					combo_mismatch[i][nTile] /= w;
				}
			} else {
				for (int i = 0; i < NUM_SLICES; i++) if (i != 1){
					combo_mismatch[i][nTile] = Double.NaN;
				}
			}
		}

		if (debugLevel > 0) { // 0) {
			String [] prefixes = {"disparity", "strength", "dx0", "dy0", "dx1", "dy1", "dx2", "dy2", "dx3", "dy3"};
			(new ShowDoubleFloatArrays()).showArrays(combo_mismatch, tilesX, combo_mismatch[0].length/tilesX, true, "removed_residual-1"+lazyEyeCompDiff, prefixes);
		}

// reduce influence of high disparity,  using combined disparity
//		double norm_ly_disparity = 100.0; // disabling
		for (int nTile = 0; nTile < num_tiles; nTile++) {
			if ((combo_mismatch[0][nTile] > 0) && (combo_mismatch[0][nTile] > ly_norm_disp)) { // why 1-st term?
				combo_mismatch[1][nTile] *= ly_norm_disp/combo_mismatch[0][nTile];
			}
		}
		if (debugLevel > 0) { // 0) {
			String [] prefixes = {"disparity", "strength", "dx0", "dy0", "dx1", "dy1", "dx2", "dy2", "dx3", "dy3"};
			(new ShowDoubleFloatArrays()).showArrays(combo_mismatch, tilesX, combo_mismatch[0].length/tilesX, true, "removed_residual-2"+lazyEyeCompDiff, prefixes);
		}


		// instance of class to operate navigation over tiles
		// compare tile disparity (combo) with those of neighbors, discard if too different
		final TileNeibs tnImage = new TileNeibs(tilesX, tilesY); // num_tiles/tilesX);
		for (int nTile = 0; nTile < num_tiles; nTile++) if (combo_mismatch[1][nTile] > 0.0){
			if (nTile == dbg_nTile){
				System.out.println("lazyEyeCorrectionFromGT().3: nTile="+nTile);
			}
			double d = combo_mismatch[0][nTile];
			double lev = lazyEyeDispVariation + lazyEyeDispRelVariation * d;
			for (int dir = 0; dir <8; dir++){
				int nTile1 = tnImage.getNeibIndex(nTile, dir);
				if ((nTile1 >= 0) && (combo_mismatch[1][nTile1] > 0.0)){
					if (Math.abs(combo_mismatch[0][nTile1] - d) > lev) {
						combo_mismatch[1][nTile] = 0.0;
						combo_mismatch[0][nTile] = Double.NaN;
						for (int i = 2; i < NUM_SLICES; i++) if (i != 1){
							combo_mismatch[i][nTile] = Double.NaN;
						}
						break;
					}
				}
			}
		}


		// here combo_mismatch[2][37005] = Double.NaN,combo_mismatch[1][37005] != 0.0, combo_mismatch[0][37005] = 0.0
		if (debugLevel > 0) { // 0) {
			String [] prefixes = {"disparity", "strength", "dx0", "dy0", "dx1", "dy1", "dx2", "dy2", "dx3", "dy3"};
			(new ShowDoubleFloatArrays()).showArrays(combo_mismatch, tilesX, combo_mismatch[0].length/tilesX, true, "combo_mismatch" , prefixes);
		}
//		if (clt_parameters.lyf_filter) {
		if (filter_lyf) {
			combo_mismatch =  filterLazyEyePairs (
					combo_mismatch, // final double[][] samples_in,
					clt_parameters.lyf_smpl_side ,   // 8,              // final int        smpl_side, // 8 x8 masked, 16x16 sampled
					clt_parameters.lyf_rms_max ,     // 0.25,           // final double     rms_max, TODO: find reasonable one not critical?
					clt_parameters.lyf_frac_keep ,   // 0.5,            // final double     frac_keep,
					clt_parameters.lyf_min_samples , // 5,              // final int        min_samples,
					clt_parameters.lyf_norm_center , // true,           // final boolean    norm_center, // if there are more tiles that fit than minsamples, replace with a single equal weight
					tilesX);        // final int        tilesX);
		}
		if (debugLevel > 0) {
			String [] prefixes = {"disparity", "strength", "dx0", "dy0", "dx1", "dy1", "dx2", "dy2", "dx3", "dy3"};
			(new ShowDoubleFloatArrays()).showArrays(combo_mismatch, tilesX, combo_mismatch[0].length/tilesX, true, "filtered_mismatch" , prefixes);
		}

		// no need to extract and filter infinity data

// make all zero strength tiles to have NaN values to use histrograms in ImageJ
		for (int nt = 0; nt < combo_mismatch[INDEX_10_WEIGHT].length; nt++ ) {
			if (combo_mismatch[INDEX_10_WEIGHT][nt] == 0.0) {
				for (int i = 0; i < NUM_SLICES; i++) if (i != INDEX_10_WEIGHT){
					combo_mismatch[i][nt] = Double.NaN;
				}
			}
		}

//	static final int    INDEX_10_WEIGHT = 1;
		if (debugLevel > 0) {
			System.out.println("test123");
		}
		if ((debugLevel > -1) && (hist_smpl_side > 0)) { // 0) {
			String [] titles = {"disparity", "strength", "dx0", "dy0", "dx1", "dy1", "dx2", "dy2", "dx3", "dy3"};
			(new ShowDoubleFloatArrays()).showArrays(combo_mismatch, tilesX, tilesY, true, "inf_and_ly",titles);

			int step = hist_smpl_side; // should be the same for both filters
			int tilesX1 = tilesX/step;
			int tilesY1 = tilesY/step;
			int num_tiles1 = tilesX1 * tilesY1;
			double [][] dbg_img = new double [combo_mismatch.length][num_tiles1];
			for (int tY = 0; tY < tilesY1; tY++) {
				for (int tX = 0; tX < tilesX1; tX++) {
					int nTile1 = tX + tY*tilesX1;
					for (int sY = 0; sY < step; sY ++) {
						for (int sX = 0; sX < step; sX ++) {
							int nTile = (sX + step * tX) + (sY + step * tY) * tilesX;
							double w = combo_mismatch[1][nTile];
							if (w > 0.0){
								for (int i = 0; i < NUM_SLICES; i++) if (i != 1) {
									dbg_img[i][nTile1] += w * combo_mismatch[i][nTile];
								}
								dbg_img[1][nTile1] += w;
							}
						}
					}

					double w = dbg_img[1][nTile1];
					if (w > 0.0){
						for (int i = 0; i < NUM_SLICES; i++) if (i != 1) {
							dbg_img[i][nTile1] /= w;
						}
					} else {
						for (int i = 0; i < NUM_SLICES; i++) if (i != 1) {
							dbg_img[i][nTile1] = Double.NaN;
						}
					}
				}
			}
			(new ShowDoubleFloatArrays()).showArrays(dbg_img, tilesX1, tilesY1, true, "inf_and_ly8",titles);
		}
		if (debugLevel > 0) {
			System.out.println("test1234a");
		}
		// create list for infinity data
//		/clt_parameters.ly_inf_en,

// adjust weight to balance infinity data and lazy eye one, so "infinity" (or really far) tiles impact is not too small even
// if there is little of infinity in the scene. As the ground truth (rig data) is known, infinity does not need to be guessed
// from the images
//	final double     inf_fraction,    // fraction of the weight for the infinity tiles
// final double     inf_max_disparity, // use all smaller disparities as inf_fraction
		double [] total_weights = new double[2];
		int num_inf = 0;
		int [] per_quad = new int[4];
		for (int nTile = 0; nTile < combo_mismatch[INDEX_10_WEIGHT].length; nTile++ ) if (center_mask[nTile] && (combo_mismatch[INDEX_10_WEIGHT][nTile] > 0.0)){
			if (combo_mismatch[INDEX_10_DISPARITY][nTile] <= inf_max_disparity) {
				total_weights[0] += combo_mismatch[INDEX_10_WEIGHT][nTile];
				num_inf ++;
			} else {
				total_weights[1] += combo_mismatch[INDEX_10_WEIGHT][nTile];
			}
			int hx = (nTile % tilesX) / (tilesX/2);
			int hy = (nTile / tilesX) / (tilesY/2);
			per_quad[hx + 2*hy]++;
		}
		int [] pq = per_quad.clone();
		Arrays.sort(per_quad);
		if (per_quad[1] < min_per_quadrant) {
			if (debugLevel > -20) {
				System.out.println(String.format("Too few tiles in quadrants :[%d, %d, %d, %d], minimum for the second worst is %d", pq[0],pq[1],pq[2],pq[3],min_per_quadrant));
			}
			return null;
		}

		if (num_inf < min_inf) {
			if (debugLevel > -20) {
				System.out.println(String.format("Too few tiles at infinity (<%4f): %d minimum is %d", inf_max_disparity, num_inf, min_inf));
			}
			return null;
		}
		if (debugLevel > -20) {
			System.out.println(String.format("Tiles per quadrants :[%d, %d, %d, %d], tiles at infinity %d", pq[0],pq[1],pq[2],pq[3],num_inf));
		}

		double inf_fraction_limited =  (inf_fraction >= 0.0) ?((inf_fraction > 1.0) ? 1.0 : inf_fraction):0.0;

		double [] weights = {
				inf_fraction_limited *         (total_weights[0] + total_weights[1]) / total_weights[0],
				(1.0 - inf_fraction_limited) * (total_weights[0] + total_weights[1]) / total_weights[1],
		};
		if (num_inf < min_inf_to_scale) {
			if (debugLevel>-1) {
				System.out.println("Too few infinity tiles to boost ("+num_inf+" < "+min_inf_to_scale+", keeping original weights");
			}
		} else if (weights[0] > weights[1]) {
			if (debugLevel>-1) {
				System.out.println("Boosting weights of far tiles (weights[0]="+weights[0]+", weights[1]="+weights[1]);
			}
			for (int nTile = 0; nTile < num_tiles; nTile++) {
				if (combo_mismatch[INDEX_10_DISPARITY][nTile] <= inf_max_disparity) {
					combo_mismatch[1][nTile] *= weights[0];
				} else {
					combo_mismatch[1][nTile] *= weights[1];
				}
			}
		} else {
			if (debugLevel>-1) {
				System.out.println("There are already more far tiles than requested (weights[0]="+weights[0]+", weights[1]="+weights[1]+", so keeping original weights");
			}
		}

		ArrayList<Sample> samples_list = new ArrayList<Sample>();

		for (int nTile = 0; nTile < num_tiles; nTile++) if (combo_mismatch[INDEX_10_WEIGHT][nTile] > 0.0) {
			samples_list.add(new Sample(0, nTile, combo_mismatch[INDEX_10_WEIGHT][nTile])); // first should be 0 to use disparity
			if (Double.isNaN(combo_mismatch[INDEX_10_DISPARITY][nTile] )) {
				System.out.println("lazyEyeCorrectionFromGT(): Double.isNaN(combo_mismatch["+INDEX_10_DISPARITY+"]["+nTile+"])");
			}
		}

		if (debugLevel > 1) {
			double inf_weight = 0.0;
			for (Sample s: samples_list) {
				inf_weight += s.weight;
			}
			System.out.println("lazyEyeCorrectionFromGT(): number of all samples="+samples_list.size()+", total weight = "+inf_weight);

		}

		if (debugLevel > 1) {
			String [] titles = {"disparity", "strength", "dx0", "dy0", "dx1", "dy1", "dx2", "dy2", "dx3", "dy3"};
			(new ShowDoubleFloatArrays()).showArrays(filtered_scans, tilesX, tilesY, true, "filtered_scans_a" , titles);
		}


		if (debugLevel > 1) {
			String [] prefixes = {"disparity", "strength", "dx0", "dy0", "dx1", "dy1", "dx2", "dy2", "dx3", "dy3"};
			(new ShowDoubleFloatArrays()).showArrays(combo_mismatch, tilesX, combo_mismatch[0].length/tilesX, true, "combo_mismatch4" , prefixes);
		}


		ArrayList<Mismatch> mismatch_list = use_poly? null : (new ArrayList<Mismatch>());
		// inf_and_ly here has filtered disparity and offsets, should be process clt_parameters.ly_inf_disp before filters
		// for rig with known disparity - use series = 0 - it will allow disparity adjustment
		double [][][] mismatch_corr_coefficients = infinityMismatchCorrection(
				clt_parameters.disp_scan_start, // final double  disp_scan_start,
				clt_parameters.disp_scan_step,  // final double  disp_scan_step,
				use_poly,                       // final boolean use_poly,
				clt_parameters.fcorr_quadratic, // final boolean use_quadratic,
				true, // clt_parameters.fcorr_inf_vert,  // final boolean use_vertical,
				// too late to restore disparity - should be dome earlier
				false,                          // final boolean use_disparity, // for infinity
				true, // clt_parameters.ly_inf_disp,     //final boolean allow_dispatity,
				clt_parameters,                 // EyesisCorrectionParameters.CLTParameters           clt_parameters,
				combo_mismatch,                 // double [][] disp_strength,
				samples_list,                   // ArrayList<Sample> samples_list,
				tilesX,                         // int         tilesX,
				magic_coeff,                    // double      , // still not understood coefficient that reduces reported disparity value.  Seems to be around 8.5
				mismatch_list,                  // ArrayList<Mismatch> mismatch_list,
				debugLevel);                    // int debugLevel)
		if (debugLevel > -2) {
			System.out.println("===== lazyEyeCorrectionFromGT(): correction coefficients =====");
			if (mismatch_corr_coefficients != null) {
				show_fine_corr(
						mismatch_corr_coefficients,
						"mismatch_corr_coefficients");
			} else {
				System.out.println("Are null - non-null are for poly correction only");
			}
		}
		// TODO: use geometryCorrection_main (if not null)
		if (!use_poly && (mismatch_list != null)){
			double [] old_new_rms = new double[1];
			boolean apply_extrinsic = true;
			int solveCorr_debug =  ((clt_parameters.lym_iter == 1) && (clt_parameters.ly_par_sel != 0))? 2 : debugLevel;
			GeometryCorrection.CorrVector corr_vector = solveCorr (
					clt_parameters.ly_inf_en,      // boolean use_disparity,     // if true will ignore disparity data even if available (was false)
					clt_parameters.ly_aztilt_en,// boolean use_aztilts,       // Adjust azimuths and tilts excluding disparity
					clt_parameters.ly_diff_roll_en,// boolean use_diff_rolls,    // Adjust differential rolls (3 of 4 angles)
					clt_parameters.ly_inf_force,   // boolean force_convergence, // if true try to adjust convergence (disparity, symmetrical parameter 0) even with no disparity
					clt_parameters.ly_com_roll,    // boolean    common_roll,    // Enable common roll (valid for high disparity range only)
					clt_parameters.ly_focalLength, // boolean    corr_focalLength,     // Correct scales (focal length temperature? variations)
					clt_parameters.ly_par_sel,     // int     manual_par_sel,    // Manually select the parameter mask bit 0 - sym0, bit1 - sym1, ... (0 - use boolean flags, != 0 - ignore boolean flags)
					mismatch_list,                          // ArrayList<Mismatch> mismatch_list,
					qc.geometryCorrection,                  // GeometryCorrection geometryCorrection,
///					geometryCorrection_main, //  GeometryCorrection geometryCorrection_main, // if is aux camera using main cameras' coordinates. Disparity is still in aux camera pixels
					qc.geometryCorrection.getCorrVector(),  // GeometryCorrection.CorrVector corr_vector,
					old_new_rms,                            // double [] old_new_rms, // should be double[2]
//					2); // debugLevel); // 2); // 1); // int debugLevel)
					solveCorr_debug); // debugLevel); // 2); // 1); // int debugLevel)
//TODO: ** Put 2 here to debug derivative images (diff_dmv_dsym - does not match yet, probably different "addition" of angles)

			if (debugLevel > -1){
				System.out.println("Old extrinsic corrections:");
				System.out.println(qc.geometryCorrection.getCorrVector().toString());
				System.out.println("Delta extrinsic corrections:");
				System.out.println(corr_vector.toString());
			}
			if (apply_extrinsic){
				qc.geometryCorrection.getCorrVector().incrementVector(corr_vector, clt_parameters.ly_corr_scale);
				if (debugLevel > -1){
					System.out.println("New extrinsic corrections:");
					System.out.println(qc.geometryCorrection.getCorrVector().toString());
				}
			}
			mismatch_corr_coefficients = new double [1][2][];
			mismatch_corr_coefficients[0][0] = corr_vector.toSymArray(null);
			mismatch_corr_coefficients[0][1] = old_new_rms;
		} else {
			if (debugLevel > -2){
				System.out.println("Extrinsic parameters (tilt, azimuth, roll) of subcameras is disabled, use_poly="+
						use_poly+" (should be false for extrinsics)");
				System.out.println(qc.geometryCorrection.getCorrVector().toString());
			}
			return mismatch_corr_coefficients;
		}
		return mismatch_corr_coefficients;
	}



	public double [][] combineCltMismatches(
			CLTParameters clt_parameters,
			double [][][]                            clt_mismatches,
			double [][][]                            disparity_maps,
			int                                      disparity_index,
			int                                      strength_index)
	{
		int n = clt_mismatches.length;
		double [][] combo = new double [clt_parameters.disp_scan_count * AlignmentCorrection.NUM_ALL_SLICES][];
		for (int pair = 0; pair < 4; pair++){
			for (int i = 0; i < n; i++){
				combo[(2 * pair + 0) * n + i] = clt_mismatches[i][3 * pair + 0];
				combo[(2 * pair + 1) * n + i] = clt_mismatches[i][3 * pair + 1];
				combo[(2 * 4 + pair) * n + i] = clt_mismatches[i][3 * pair + 2];
			}
		}
		for (int i = 0; i < n; i++){
			combo[12 * n + i] = disparity_maps[i][disparity_index];
			combo[13 * n + i] = disparity_maps[i][strength_index];
		}
		return combo;
	}

	public void showCltMismatches(
			String                                   title,
			CLTParameters clt_parameters,
			double [][]                              combo_data,
			int tilesX,
			int tilesY)
	{
		ShowDoubleFloatArrays sdfa_instance = new ShowDoubleFloatArrays();
		String [] titles = new String [combo_data.length];
		int num_scans = combo_data.length / AlignmentCorrection.NUM_ALL_SLICES;
		for (int pair = 0; pair < 4; pair++){
			for (int i = 0; i < num_scans; i++){
				double disparity = clt_parameters.disp_scan_start + i * clt_parameters.disp_scan_step;

				titles[(2 * pair + 0) * num_scans + i] = "dx_"+pair+"_"+disparity;
				titles[(2 * pair + 1) * num_scans + i] = "dy_"+pair+"_"+disparity;
				titles[(2 * 4 + pair) * num_scans + i] = "strength_"+pair+"_"+disparity;
			}
		}
		for (int i = 0; i < num_scans; i++){
			double disparity = clt_parameters.disp_scan_start + i * clt_parameters.disp_scan_step;
			titles[ 12 * num_scans + i] = "disp_"+disparity;
			titles[ 13 * num_scans + i] = "strength_"+disparity;
		}
		sdfa_instance.showArrays(combo_data, tilesX,tilesY, true, title, titles);
	}

	public void showCltMismatch(
			String                                   title,
			CLTParameters clt_parameters,
			double [][]                              clt_mismatch,
			int tilesX,
			int tilesY)

	{
		ShowDoubleFloatArrays sdfa_instance = new ShowDoubleFloatArrays();
		String [] titles = new String [12];
		double [][] dbg_clt_mismatch= new double [12][];
		for (int pair = 0; pair < 4; pair++){
			titles[2 * pair + 0] = "dx_"+pair;
			titles[2 * pair + 1] = "dy_"+pair;
			titles[2 * 4 + pair] = "strength_"+pair;
			dbg_clt_mismatch[(2 * pair + 0)] = clt_mismatch[3 * pair + 0].clone();
			dbg_clt_mismatch[(2 * pair + 1)] = clt_mismatch[3 * pair + 1].clone();
			dbg_clt_mismatch[(2 * 4 + pair)] = clt_mismatch[3 * pair + 2];
			for (int i = 0; i < dbg_clt_mismatch[(2 * 4 + pair)].length; i++ ){
				if (dbg_clt_mismatch[(2 * 4 + pair)][i] == 0.0){
					dbg_clt_mismatch[(2 * pair + 0)][i] = Double.NaN;
					dbg_clt_mismatch[(2 * pair + 1)][i] = Double.NaN;
				}
			}
		}
		sdfa_instance.showArrays(dbg_clt_mismatch, tilesX, tilesY, true, title, titles);
	}
	/**
	 * Calculates transposed Jacobian for 8*<number of tiles> points (port x,y coordinates) for each of the symmetrical parameters (starting with zoom)
	 * @param par_mask           bitmask of selected parameters, starting with sym0 (convergence)
	 * @param mismatch_list      list of samples
	 * @param geometryCorrection instance of the camera geometry class
	 * @param corr_vector        current correction vector
	 * @param debugLevel         debug level
	 * @return                   transposed Jacobian
	 */

	double [][] getJacobianTransposed(
			boolean [] par_mask,
			ArrayList<Mismatch> mismatch_list,
			GeometryCorrection geometryCorrection,
			GeometryCorrection.CorrVector corr_vector,
			int debugLevel)
	{
		boolean dbg_images = debugLevel>1;
		int dbg_decimate = 64; // just for the debug image
		int dbg_width =  qc.tp.getTilesX()*qc.tp.getTileSize();
		int dbg_height = qc.tp.getTilesY()*qc.tp.getTileSize();
		int dbg_owidth = dbg_width/dbg_decimate;
		int dbg_oheight = dbg_height/dbg_decimate;
		int dbg_length = dbg_owidth*dbg_oheight;
		String [] dbg_titles_tar=GeometryCorrection.CORR_NAMES;
		String [] dbg_titles_sym= {"sym0","sym1","sym2","sym3","sym4","sym5","sroll0","sroll1","sroll2","sroll3"};
		String [] dbg_titles_xy=  {"x0","y0","x1","y1","x2","y2","x3","y3"};
		double [][] dbg_img_deriv = null; // compare derivatives with delta-diffs
		double [][] dbg_dxy_dsym = null;  // jacobian dxy/dsym
		if (dbg_images) {
			dbg_img_deriv = doubleNaN(dbg_titles_xy.length * dbg_titles_tar.length *2, dbg_length); // compare derivatives with delta-diffs
			dbg_dxy_dsym =  doubleNaN(dbg_titles_xy.length * dbg_titles_sym.length,    dbg_length); // jacobian dxy/dsym
		}

		int num_pars = 0;
		for (int i = 0; i < par_mask.length; i++) if (par_mask[i]) num_pars ++;

		double [][] jt = new double  [num_pars][2 * NUM_SENSORS * mismatch_list.size()];
		double [][] jt_dbg = null;
		if (debugLevel > 0){
			jt_dbg = new double  [num_pars][2 * NUM_SENSORS * mismatch_list.size()];
		}

		Matrix []   corr_rots =  corr_vector.getRotMatrices(); // get array of per-sensor rotation matrices
		Matrix [][] deriv_rots = corr_vector.getRotDeriveMatrices();
		for (int indx = 0; indx<mismatch_list.size(); indx++){ // need indx value
			Mismatch mm = mismatch_list.get(indx);
			double [] pXY = mm.getPXY();
			double [][] deriv = new double [2 * NUM_SENSORS][];
			int dbg_index =dbg_index (pXY, dbg_decimate);
			double [][] disp_dist = new double[dbg_titles_xy.length][]; // used to correct 3D correlations
			geometryCorrection.getPortsCoordinatesAndDerivatives(
					geometryCorrection, //			GeometryCorrection gc_main,
					false,       // boolean use_rig_offsets,
					corr_rots,   //  Matrix []   rots,
					deriv_rots,  //  Matrix [][] deriv_rots,
					deriv,       // 	boolean calc_deriv,
					disp_dist,   // used to correct 3D correlations
					pXY[0],      // double px,
					pXY[1],      // double py,
					mm.getDisparityMeas()); // getDisparityTask()); // double disparity)

			// convert to symmetrical coordinates
			// derivatives of each port coordinates (in pixels) for each of selected symmetric all parameters (sym0 is convergence for disparity)
			 double [][] jt_partial = corr_vector.getJtPartial(
						deriv, // double [][] port_coord_deriv,
						par_mask); // boolean [] par_mask

			// put partial (for 1 tile - 8 port coordinates) transposed jacobian into full (all tiles) transposed Jacobian
			for (int npar = 0; npar < jt.length; npar++){
				for (int n = 0; n < 2* NUM_SENSORS; n++){
//						jt[npar][2 * NUM_SENSORS * indx + n] = j_partial[n][npar]; // here Jacobian was not transposed
						jt[npar][2 * NUM_SENSORS * indx + n] = jt_partial[npar][n];
						if (Double.isNaN(jt_partial[npar][n])) {
							System.out.println("getJacobianTransposed(): npar="+npar+", indx="+indx+", n="+n);
						}
				}
			}
			if (debugLevel > 0){
				double [][] deriv_dbg = new double [2 * NUM_SENSORS][];
				double [] dbg_a_vector= null;
				geometryCorrection.getPortsCoordinatesAndDerivatives(
						false,          // boolean use_rig_offsets,
						dbg_a_vector, // double [] dbg_a_vector, // replace actual radial distortion coefficients
						1E-9, // 1E-8, //6,    // double delta, // 1e-6
						corr_vector, // CorrVector corr_vector,
						deriv_dbg, // j_partial_debug, //
						null, // disp_dist,       // used to correct 3D correlations
						pXY[0],      // double px,
						pXY[1],      // double py,
						mm.getDisparityMeas()); // Task()); // double disparity)

				// convert to symmetrical coordinates
				 double [][] jt_partial_dbg = corr_vector.getJtPartial(
							deriv_dbg, // double [][] port_coord_deriv,
							par_mask); // boolean [] par_mask

				double max_rdiff = 0;
				for (int npar = 0; npar < jt.length; npar++){
					for (int n = 0; n < 2 * NUM_SENSORS; n++){
//						jt_dbg[npar][2 * NUM_SENSORS * indx + n] = j_partial_debug[n][npar];  // here Jacobian was not transposed
						jt_dbg[npar][2 * NUM_SENSORS * indx + n] = jt_partial_dbg[npar][n];
						double avg  = 0.5*(jt_dbg[npar][2 * NUM_SENSORS * indx + n] + jt[npar][2 * NUM_SENSORS * indx + n]);
						double rdiff= Math.abs(0.5*(jt_dbg[npar][2 * NUM_SENSORS * indx + n] - jt[npar][2 * NUM_SENSORS * indx + n]));
						if (avg != 0.0){
							rdiff /= avg;
						}
						if (rdiff > max_rdiff){
							max_rdiff = rdiff;
						}
					}
				}
				if (debugLevel > 2) {
					System.out.println(String.format("px = %5.0f py = %5.0f disp_task = %7.3f disp_meas = %7.3f strength = %7.3f max rdiff = %f",
							mm.getPXY()[0], mm.getPXY()[1], mm.getDisparityTask(), mm.getDisparityMeas(), mm.getStrength(), max_rdiff));
				}
				if (dbg_images) {
					for (int i = 0; i < dbg_titles_xy.length; i++){
						for (int j = 0; j < dbg_titles_tar.length; j++){
							dbg_img_deriv[2 * (i * dbg_titles_tar.length + j) + 0][dbg_index] = deriv[i][j];
							dbg_img_deriv[2 * (i * dbg_titles_tar.length + j) + 1][dbg_index] = deriv_dbg[i][j];
						}
					}
					for (int i = 0; i < dbg_titles_xy.length; i++){
						int oj = 0;
						for (int j = 0; j < dbg_titles_sym.length; j++) if (par_mask[j]){
							dbg_dxy_dsym[i * dbg_titles_sym.length + j][dbg_index] =  jt_partial[oj][i];
							oj++;
						}
					}
				}
			}
		}
		if (dbg_images) {
			String [] dbg_img_deriv_titles = new String [dbg_titles_xy.length * dbg_titles_tar.length *2];
			for (int i = 0; i < dbg_titles_xy.length; i++){
				for (int j = 0; j < dbg_titles_tar.length; j++){
					dbg_img_deriv_titles[2 * (i * dbg_titles_tar.length + j) + 0]= dbg_titles_xy[i] + "_" +dbg_titles_tar[j];
					dbg_img_deriv_titles[2 * (i * dbg_titles_tar.length + j) + 1]= dbg_titles_xy[i] + "_" +dbg_titles_tar[j] + "delta";
				}
			}

			String [] dbg_dxy_dsym_titles = new String [dbg_titles_xy.length * dbg_titles_sym.length];
			for (int i = 0; i < dbg_titles_xy.length; i++){
				for (int j = 0; j < dbg_titles_sym.length; j++){
					dbg_dxy_dsym_titles[i * dbg_titles_sym.length + j]= dbg_titles_xy[i] + "_" +dbg_titles_sym[j];
				}
			}

			dbgImgRemoveEmpty(dbg_img_deriv);
			dbgImgRemoveEmpty(dbg_dxy_dsym);

			ShowDoubleFloatArrays sdfa_instance = new ShowDoubleFloatArrays();
			sdfa_instance.showArrays(dbg_img_deriv, dbg_owidth, dbg_oheight, true, "dbg_img_deriv", dbg_img_deriv_titles);
			sdfa_instance.showArrays(dbg_dxy_dsym,  dbg_owidth, dbg_oheight, true, "dbg_dxy_dsym",  dbg_dxy_dsym_titles);
		}
		return jt;
	}


	double [][] debug_mv_from_sym_jacobian(
			double delta,
			boolean [] par_mask,
			ArrayList<Mismatch> mismatch_list,
			GeometryCorrection geometryCorrection,
			GeometryCorrection.CorrVector corr_vector,
			int debugLevel)
	{
		int num_pars = 0;
		for (int i = 0; i < par_mask.length; i++) if (par_mask[i]) num_pars ++;

		double [][] jt_mv = new double  [num_pars][2 * NUM_SENSORS * mismatch_list.size()];
		double [] sym_par_0 = corr_vector.toSymArray(par_mask);
		for (int sym_par = 0; sym_par < num_pars; sym_par++ ) {
			double [] sym_par_p = sym_par_0.clone();
			double [] sym_par_m = sym_par_0.clone();
			sym_par_p[sym_par] += 0.5 * delta;
			sym_par_m[sym_par] -= 0.5 * delta;
			GeometryCorrection.CorrVector corr_p = geometryCorrection.getCorrVector(sym_par_p, par_mask);
			GeometryCorrection.CorrVector corr_m = geometryCorrection.getCorrVector(sym_par_m, par_mask);
			double [] mv_p = debug_mv_from_sym(
					mismatch_list,
					geometryCorrection,
					corr_p,
					debugLevel);
			double [] mv_m = debug_mv_from_sym(
					mismatch_list,
					geometryCorrection,
					corr_m,
					debugLevel);
			for (int i = 0; i < jt_mv[sym_par].length; i++){
				jt_mv[sym_par][i] = (mv_p[i]-mv_m[i])/delta;
			}

		}
		return jt_mv;
	}


	/**
	 * Debugging jacobian with two coordinate transformations - input and output. Calculating output mv vector
	 * for all coordinate points for current corr_vector (to use it with delta corr_vecotr)
	 * @param mismatch_list
	 * @param geometryCorrection
	 * @param corr_vector
	 * @param debugLevel
	 * @return
	 */
	double [] debug_mv_from_sym(
//			boolean [] par_mask,
			ArrayList<Mismatch> mismatch_list,
			GeometryCorrection geometryCorrection,
			GeometryCorrection.CorrVector corr_vector,
			int debugLevel)
	{
		double [][] dMismatch_dXY = (new Mismatch()).get_dMismatch_dXY(); // just a static array
		double [] mv = new double [2 * NUM_SENSORS * mismatch_list.size()];
		Matrix [] corr_rots = corr_vector.getRotMatrices(); // get array of per-sensor rotation matrices
		for (int indx = 0; indx<mismatch_list.size(); indx++){ // need indx value
			Mismatch mm = mismatch_list.get(indx);
			double [] pXY = mm.getPXY();
			double [][] f = geometryCorrection.getPortsCoordinatesAndDerivatives( // 4x2
					geometryCorrection, //			GeometryCorrection gc_main,
					false,          // boolean use_rig_offsets,
					corr_rots,   //  Matrix []   rots,
					null,        //  Matrix [][] deriv_rots,
					null,        //	boolean calc_deriv,
					null,        // disp_dist,       // used to correct 3D correlations
					pXY[0],      // double px,
					pXY[1],      // double py,
					mm.getDisparityMeas()); // getDisparityTask()); // double disparity)

			// convert to symmetrical coordinates
			// f is [4][2] array of port x,y coordinates - convert them to mv (linear array)

			double [] mv_partial = new double [dMismatch_dXY.length];
			for (int i = 0; i < mv_partial.length; i++){
				for (int nsens = 0; nsens < NUM_SENSORS; nsens++){
					for (int dir = 0; dir <2; dir++){
						mv_partial[i]+=f[nsens][dir]* dMismatch_dXY[i][2*nsens+dir];
					}
				}
			}
			for (int n = 0; n < 2* NUM_SENSORS; n++){
				mv[2 * NUM_SENSORS * indx + n] = mv_partial[n];
			}
		}
		return mv;
	}

	public int dbg_index(double [] pXY, int decimate)
	{
		int width =  qc.tp.getTilesX()*qc.tp.getTileSize()/decimate;
		int height = qc.tp.getTilesY()*qc.tp.getTileSize()/decimate;
		int x = (int) Math.round(pXY[0]/decimate);
		int y = (int) Math.round(pXY[1]/decimate);
		if (x < 0) x = 0;
		else if (x >= width) x = width - 1;
		if (y < 0) y = 0;
		else if (y >= height) y = height - 1;
		return x + width * y;
	}
	public double [][] doubleNaN(int layers, int leng){
		double [][] dbg_img = new double [layers][leng];
		for (int nimg = 0; nimg < dbg_img.length; nimg++) if (dbg_img[nimg] != null){
			for (int i = 0; i < dbg_img[nimg].length; i++){
				dbg_img[nimg][i] = Double.NaN;
			}
		}
		return dbg_img;
	}

	public void dbgImgRemoveEmpty(double [][] dbg_img){
		for (int nimg = 0; nimg < dbg_img.length; nimg++) if (dbg_img[nimg] != null){
			boolean has_data = false;
			for (int i = 0; i < dbg_img[nimg].length; i++){
				if (!Double.isNaN(dbg_img[nimg][i]) && (dbg_img[nimg][i] != 0.0)){
					has_data = true;
					break;
				}
			}
			if (!has_data){
				dbg_img[nimg] = null;
			}
		}
	}

	double [] getJtJTrace( // just debugging
			double [][] jt,
			double [] w)
	{
		double [] jtj_trace = new double [jt.length];
		if (w == null){
			w = new double[jt[0].length];
			for (int i = 0; i < w.length; i++){
				w[i] = 1.0/w.length;
			}
		}

		for (int i = 0; i < jt.length; i++){
			for (int k = 0; k < jt[i].length; k++){
				jtj_trace[i] += w[k] * jt[i][k]*jt[i][k];
			}
		}
		return jtj_trace;

	}

	double [][] getJTJ(
			double [][] jt,
			double [] w)
	{
		double [][] jtj = new double [jt.length][jt.length];
		for (int i = 0; i < jt.length; i++){
			for (int j = 0; j < i; j++){
				jtj[i][j] = jtj[j][i];
			}
			for (int j = i; j < jt.length; j++){
				for (int k = 0; k < jt[0].length; k++){
					jtj[i][j] += jt[i][k] * jt[j][k] * w[k];
					if (Double.isNaN(jtj[i][j])) {
						System.out.println("i="+i+", j="+j+", k="+k);
					}
				}
			}
		}
		return jtj;
	}

	double [] mulWeight(double [] y, double [] w){
		double [] yw = y.clone();
		for (int i = 0; i < yw.length; i++){
			yw[i] *= w[i];
		}
		return yw;
	}

	double getSumWeight(double [] w){
		double sum_w = 0.0;
		for (int i = 0; i < w.length; i++){
			sum_w += w[i];
		}
		return sum_w;
	}

	double getRMS(double [] y, double [] w){
		double rms = 0.0;
		double sw = 0.0;
		for (int i = 0; i < w.length; i++){
			if (Double.isNaN(y[i])) {
				// TODO: Fix source !
				System.out.println("getRMS(): y["+i+"]="+y[i]+", w[i]="+w[i]);
				w[i] = 0.0;
			}
			rms += w[i]*y[i]*y[i];
			sw += w[i];
		}
		if (sw != 0.0){
			rms = Math.sqrt(rms/sw);
		}
		return rms;
	}

	/**
	 * Create Y-F(x) vector for 8 points per tile data from the mismatch_list. Groups of 8 uses linear combinations of the measured {dx0, dy0,..dx3,dy3},
	 * and the first 7 points in each group are invariant of disparity, the 8-th is optionally disabled by the corresponding weight (for points with unknown disparity)
	 *
	 * Actually there is no explicit F(x), because after each step, the correlations are recalculated instead, and the measured are already differences from "F(x)"
	 *
	 * @param mismatch_list list of the tile samples
	 * @return array of 8*<number of tiles> values
	 */

	double [] getYminusFx(
			ArrayList<Mismatch> mismatch_list)
	{
		double [] yMinusFx = new double [2 * NUM_SENSORS * mismatch_list.size()];
		double [] w_DBG =    new double [2 * NUM_SENSORS * mismatch_list.size()];

		for (int indx = 0; indx<mismatch_list.size(); indx++){ // need indx value
			Mismatch mm = mismatch_list.get(indx);
			mm.copyToY(yMinusFx, indx);

			mm.copyToW(w_DBG, indx); // temporarily

		}
		// check for NaN:
		// temporary
		for (int i = 0; i < yMinusFx.length; i++) {
			if (Double.isNaN(yMinusFx[i]) && (w_DBG[i] != 0 )) {
				System.out.println("**** getYminusFx BUG! y["+i+"] = NaN, while w is != 0 *****");
				System.out.println("**** getYminusFx BUG! y["+i+"] = NaN, while w is != 0 *****");
			}
		}
		return yMinusFx;
	}

	double [] getWeights(
			ArrayList<Mismatch> mismatch_list)
	{
		double [] w = new double [2 * NUM_SENSORS * mismatch_list.size()];
		for (int indx = 0; indx < mismatch_list.size(); indx++){ // need indx value
			Mismatch mm = mismatch_list.get(indx);
			mm.copyToW(w, indx);
		}
		double sumw = 0.0;
		for (int i = 0; i < w.length; i++) sumw += w[i];
		if (sumw > 0.0){
			for (int i = 0; i < w.length; i++) w[i]/=sumw;
		}

		return w;
	}

	public GeometryCorrection.CorrVector  solveCorr (
			boolean use_disparity,     // adjust disparity-related extrinsics
//			boolean use_other_extr,    // adjust other extrinsic parameters that do not influence disparity, common roll and zoom
			boolean use_aztilts,       // Adjust azimuths and tilts excluding disparity
			boolean use_diff_rolls,    // Adjust differential rolls (3 of 4 angles)
			boolean force_convergence, // if true try to adjust convergence (disparity, symmetrical parameter 0) even with no disparity
			                           // data, using just radial distortions
	  		boolean common_roll,       // Enable common roll (valid for high disparity range only)
			boolean corr_focalLength,  // Correct scales (focal length temperature? variations)
	  		int     manual_par_sel,    // Manually select the parameter mask bit 0 - sym0, bit1 - sym1, ... (0 - use boolean flags, != 0 - ignore boolean flags)
			ArrayList<Mismatch> mismatch_list,
			GeometryCorrection geometryCorrection,
///			GeometryCorrection geometryCorrection_main, // if is aux camera using main cameras' coordinates. Disparity is still in aux camera pixels
			GeometryCorrection.CorrVector corr_vector,
			double [] old_new_rms, // should be double[2]
			int debugLevel)
	{
		boolean dbg_images = debugLevel > 0; // 1;
		boolean has_disparity = force_convergence; // force false;
		// See if at least some data has disparities to be adjusted
		if (use_disparity) {
			for (Mismatch mm: mismatch_list){
				if (mm.use_disparity) has_disparity = true;
				break;
			}
		}

		boolean [] par_mask = geometryCorrection.getParMask(
// temporary - just for testing
//				force_convergence, // boolean disparity_only,
//				force_convergence && has_disparity, // boolean use_disparity,
				has_disparity, // boolean use_disparity,
//				use_other_extr,    // boolean use_other_extr,
				use_aztilts,       // Adjust azimuths and tilts excluding disparity
				use_diff_rolls,    // Adjust differential rolls (3 of 4 angles)

				common_roll,// boolean common_roll,
				corr_focalLength, // boolean corr_focalLength);
		  		manual_par_sel);    // Manually select the parameter mask bit 0 - sym0, bit1 - sym1, ... (0 - use boolean flags, != 0 - ignore boolean flags)


		double [][] jta = getJacobianTransposed( // gets transposed jacobian for 8*num_tiles points (port coordinates) and symmetrical parameters
				par_mask,           // boolean [] par_mask,
				mismatch_list,      // ArrayList<Mismatch> mismatch_list,
				geometryCorrection, // GeometryCorrection geometryCorrection,
				corr_vector,        // GeometryCorrection.CorrVector corr_vector)
				debugLevel);		// int debugLevel)

//		debugLevel = 2;

		// Convert Jacobian outputs to symmetrical measurement vectors (last one is non-zero only if disparity should be adjusted)
		// so now each group of 8 points (individual coordinate pairs for each sub-camera) is replaced by 8-point "measurement vector", where only the last component
		// depends on disparity. That last component is optionally disabled by the corresponding weight for the samples that do not have absolute disparity data
		double [][] jta_mv =  (new Mismatch()).convertJt_mv (jta); //double [][] jt)

		Matrix jt = new Matrix(jta_mv);


		// Extract measured offsets differences from the list and convert them to the measured vectors for which the jta is intended
		double [] y_minus_fx_a = getYminusFx( // mv[0]..mv[7], not the measured data (dx0, dy0, ... dx3, dy3)
				mismatch_list); // ArrayList<Mismatch> mismatch_list)

		double [] weights = getWeights(
				mismatch_list); // ArrayList<Mismatch> mismatch_list)
		double [] y_minus_fx_a_weighted = mulWeight(y_minus_fx_a, weights);
		double rms0 = getRMS	(y_minus_fx_a, weights);
		if (debugLevel > -3){
			System.out.println("--- solveCorr(): initial RMS = " + rms0);
		}

		if (old_new_rms != null){
			old_new_rms[0] =rms0;
		}
		Matrix y_minus_fx_weighted = new Matrix(y_minus_fx_a_weighted, y_minus_fx_a_weighted.length);
		double [][] jtja = getJTJ(jta_mv, weights);
		Matrix jtj = new Matrix(jtja); // getJTJ(jta, weights)); // less operations than jt.times(jt.transpose());
		int dbg_decimate = 64; // just for the debug image
		int dbg_width =  qc.tp.getTilesX()*qc.tp.getTileSize();
		int dbg_height = qc.tp.getTilesY()*qc.tp.getTileSize();
		int dbg_owidth = dbg_width/dbg_decimate;
		int dbg_oheight = dbg_height/dbg_decimate;
		while (dbg_owidth < 40) {
			dbg_decimate /= 2;
			dbg_owidth = dbg_width/dbg_decimate;
			dbg_oheight = dbg_height/dbg_decimate;
		}
		int dbg_length = dbg_owidth*dbg_oheight;
		String [] dbg_titles_sym= {"sym0","sym1","sym2","sym3","sym4","sym5","sroll0","sroll1","sroll2","sroll3", "zoom0", "zoom1", "zoom2"};
		String [] dbg_titles_xy=  {"dx0","dy0","dx1","dy1","dx2","dy2","dx3","y3"};
		String [] dbg_titles_mv=  {"dy0","dy1","dx2","dx3","dx1-dx0","dy3-dy2","dh-dv","dh+dv"};
		double [][] dbg_xy = null;  // jacobian dmv/dsym
		double [][] dbg_mv = null;  // jacobian dmv/dsym
		double [][] dbg_dmv_dsym = null;  // jacobian dmv/dsym
		double [][] dbg_dmv_dsym_delta = null;  // jacobian dmv/dsym
		double [][] dbg_dmv_dsym_diff =  null;  // jacobian dmv/dsym
		if (dbg_images) {
			double [][] jta_mv_delta = debug_mv_from_sym_jacobian(
					1e-9, // 6,               // double delta,
					par_mask,           // boolean [] par_mask,
					mismatch_list,      // ArrayList<Mismatch> mismatch_list,
					geometryCorrection, // GeometryCorrection geometryCorrection,
					corr_vector,        // GeometryCorrection.CorrVector corr_vector)
					debugLevel);		// int debugLevel)

			dbg_xy =              doubleNaN(dbg_titles_xy.length,                            dbg_length); // jacobian dmv/dsym
			dbg_mv =              doubleNaN(dbg_titles_mv.length,                            dbg_length); // jacobian dmv/dsym
			dbg_dmv_dsym =        doubleNaN(dbg_titles_mv.length * dbg_titles_sym.length,    dbg_length); // jacobian dmv/dsym
			dbg_dmv_dsym_delta =  doubleNaN(dbg_titles_mv.length * dbg_titles_sym.length,    dbg_length); // jacobian dmv/dsym
			dbg_dmv_dsym_diff =   doubleNaN(dbg_titles_mv.length * dbg_titles_sym.length,    dbg_length); // jacobian dmv/dsym

			String [] dbg_dmv_dsym_titles = new String [dbg_titles_mv.length * dbg_titles_sym.length];
			for (int i = 0; i < dbg_titles_mv.length; i++){
				for (int j = 0; j < dbg_titles_sym.length; j++){
					dbg_dmv_dsym_titles[i * dbg_titles_sym.length + j]= dbg_titles_mv[i] + "_" +dbg_titles_sym[j];
				}
			}
			for (int indx = 0; indx < mismatch_list.size(); indx++){
				Mismatch mm = mismatch_list.get(indx);
				double [] pXY = mm.getPXY();
				int dbg_index =dbg_index (pXY, dbg_decimate);
				double [] xy =  mm.getOffsets();
				double [] mv =  mm.getY();
				for (int i = 0; i < xy.length; i++){
					dbg_xy[i][dbg_index] = xy[i];
				}
				for (int i = 0; i < mv.length; i++){
					dbg_mv[i][dbg_index] = mv[i];
				}
				// Now Jacobian - get from full one
				//System.out.println("solveCorr(): dbg_dmv_dsym.length="+dbg_dmv_dsym.length+ ", dbg_dmv_dsym[0].length="+dbg_dmv_dsym[0].length);
				//System.out.println("solveCorr(): jta_mv.length="+jta_mv.length+", jta_mv[0].length="+jta_mv[0].length+" dbg_index="+dbg_index);
				for (int i = 0; i < dbg_titles_mv.length; i++){
					int oj = 0;
					for (int j = 0; j < dbg_titles_sym.length; j++)  if (par_mask[j]){
						if ((i * dbg_titles_sym.length + j) >= dbg_dmv_dsym.length) {
							System.out.println("solveCorr(): dbg_dmv_dsym.length="+dbg_dmv_dsym.length+ ", dbg_dmv_dsym[0].length="+dbg_dmv_dsym[0].length);
						} else if ((dbg_titles_mv.length * indx + i) >= jta_mv [oj].length){
							System.out.println("solveCorr(): dbg_dmv_dsym.length="+dbg_dmv_dsym.length+ ", dbg_dmv_dsym[0].length="+dbg_dmv_dsym[0].length);
						}
						dbg_dmv_dsym[i * dbg_titles_sym.length + j][dbg_index] =       jta_mv [oj][dbg_titles_mv.length * indx + i];  //java.lang.ArrayIndexOutOfBoundsException: 3552
						dbg_dmv_dsym_delta[i * dbg_titles_sym.length + j][dbg_index] = jta_mv_delta [oj][dbg_titles_mv.length * indx + i];
						dbg_dmv_dsym_diff[i * dbg_titles_sym.length + j][dbg_index] =  jta_mv_delta [oj][dbg_titles_mv.length * indx + i] - jta_mv [oj][dbg_titles_mv.length * indx + i];
						oj++;
					}
				}
			}

			dbgImgRemoveEmpty(dbg_xy);
			dbgImgRemoveEmpty(dbg_mv);
			dbgImgRemoveEmpty(dbg_dmv_dsym);
			dbgImgRemoveEmpty(dbg_dmv_dsym_delta);
			dbgImgRemoveEmpty(dbg_dmv_dsym_diff);

			ShowDoubleFloatArrays sdfa_instance = new ShowDoubleFloatArrays();
			sdfa_instance.showArrays(dbg_xy,              dbg_owidth, dbg_oheight, true, "dbg_xy",        dbg_titles_xy);
			sdfa_instance.showArrays(dbg_mv,              dbg_owidth, dbg_oheight, true, "dbg_mv",        dbg_titles_mv);
			sdfa_instance.showArrays(dbg_dmv_dsym,        dbg_owidth, dbg_oheight, true, "dbg_dmv_dsym",  dbg_dmv_dsym_titles);
			sdfa_instance.showArrays(dbg_dmv_dsym_delta,  dbg_owidth, dbg_oheight, true, "delta_dmv_dsym",dbg_dmv_dsym_titles);
			sdfa_instance.showArrays(dbg_dmv_dsym_diff,   dbg_owidth, dbg_oheight, true, "diff_dmv_dsym", dbg_dmv_dsym_titles);
		}
		if (debugLevel>-1) {
			jtj.print(18, 6);
		}
		Matrix jtj_inv = jtj.inverse();
		Matrix jty = jt.times(y_minus_fx_weighted);
		Matrix mrslt = jtj_inv.times(jty);
		double []  drslt = mrslt.getColumnPackedCopy();
		// wrong sign?
		for (int i = 0; i < drslt.length; i++){
			drslt[i] *= -1.0;
		}
		//if (par_mask[0]) drslt[0] *= -1.0; //FIXME: Find actual bug, sym[0] corrects in opposite way

		GeometryCorrection.CorrVector rslt = geometryCorrection.getCorrVector(drslt, par_mask);
		if (debugLevel > -3){ // change to >0) {
			System.out.println("solveCorr() rslt (increment):");
			System.out.println(rslt.toString());
			System.out.println("--- end of increment ---");
		}

		return rslt;
	}


}
