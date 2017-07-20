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

import Jama.Matrix;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.WindowManager;

public class AlignmentCorrection {
	static int NUM_SLICES =      10; // disp, strength, dx0, dy0, dx1, dy1, dx2, dy2, dx3, dy3)
	static int NUM_ALL_SLICES =  14; // disp, strength, dx0, dy0, str0, dx1, dy1, str1, dx2, dy2, str2, dx3, dy3 str3,)
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
	
	
	
	AlignmentCorrection (QuadCLT qc){
		this.qc = qc;
	}

	public double [][][] infinityCorrection(
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

			EyesisCorrectionParameters.CLTParameters           clt_parameters,
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
				double [][] dbg_img = disp_strength.clone();
				for (int n = 0; n < disp_strength.length; n++){
					dbg_img[n] = disp_strength[n].clone();
				}
				for (int n = 0; n < dbg_img.length; n+=2){
					for (int i = 0; i < dbg_img[n].length; i++) {
						if (dbg_img[n+1][i] == 0.0){
							dbg_img[n][i] = Double.NaN;
						}
					}
				}

				(new showDoubleFloatArrays()).showArrays(dbg_img, tilesX, disp_strength[0].length/tilesX, true, "filtered_ds"); // , titles);

			}

		} else {
			disp_strength = disp_strength_in;
			min_strength = min_strength_in;
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
				double [][] dbg_img = disp_strength.clone();
				for (int n = 0; n < disp_strength.length; n++){
					dbg_img[n] = disp_strength[n].clone();
				}
				for (int n = 0; n < dbg_img.length; n+=2){
					for (int i = 0; i < dbg_img[n].length; i++) {
						if (dbg_img[n+1][i] == 0.0){
							dbg_img[n][i] = Double.NaN;
						}
					}
				}
				(new showDoubleFloatArrays()).showArrays(dbg_img, tilesX, disp_strength[0].length/tilesX, true, "hist_filt_ds"); // , titles);

			}
		}

		ArrayList<Sample> samples_list = selectInfinityTiles(
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
/*		
		double [][][] disparity_corr_coefficiants = null;
		if (clt_parameters.inf_disp_apply) {
			disparity_corr_coefficiants = infinityCorrection(
					clt_parameters.fcorr_inf_vert,// final boolean use_vertical,
					//				min_strength,
					//				max_diff,
					//				max_iterations,
					//				max_coeff_diff,
					//				far_pull, //  = 0.2; // 1; //  0.5;
					clt_parameters,
					disp_strength,
					
					samples_list,
					tilesX,
					magic_coeff, // still not understood coefficient that reduces reported disparity value.  Seems to be around 8.5  
					debugLevel);
			if (debugLevel > -1){
				System.out.println("infinityCorrection(): coefficient increments from infinityCorrection");
				show_fine_corr(
						disparity_corr_coefficiants, // double [][][] corr,
						"");// String prefix)
			}
		}
*/		
		double [][][] mismatch_corr_coefficiants = null;
//		if (clt_parameters.inf_disp_apply) {
			mismatch_corr_coefficiants = infinityMismatchCorrection(
					clt_parameters.fcorr_inf_quad, // final boolean use_quadratic,
					clt_parameters.fcorr_inf_vert, // final boolean use_vertical,
					clt_parameters.ly_inf_en,      // final boolean use_disparity, // for infinity 
					clt_parameters,
					disp_strength,
					samples_list,
					tilesX,
					magic_coeff, // still not understood coefficient that reduces reported disparity value.  Seems to be around 8.5  
					debugLevel);
			if (debugLevel > -1){
				System.out.println("infinityCorrection(): coefficient increments from infinityMismatchCorrection");
				show_fine_corr(
						mismatch_corr_coefficiants, // double [][][] corr,
						"");// String prefix)
			}
/*			
			if (disparity_corr_coefficiants == null) {
				disparity_corr_coefficiants = mismatch_corr_coefficiants;
				if (debugLevel > -1){
					System.out.println("infinityCorrection(): using only coefficient increments from infinityMismatchCorrection");
				}

			} else { 
				for (int i = 0; i < disparity_corr_coefficiants.length; i++){
					for (int j = 0; j < disparity_corr_coefficiants[i].length; j++){
						for (int k = 0; k < disparity_corr_coefficiants[i][j].length; k++){
							disparity_corr_coefficiants[i][j][k] += mismatch_corr_coefficiants[i][j][k];
						}
					}
				}
				if (debugLevel > -1){
					System.out.println("infinityCorrection(): combining coefficient increments from infinityCorrection and infinityMismatchCorrection");
				}
			}
*/			
//		}
		return mismatch_corr_coefficiants;
	}
	
	/**
	 * Select infinity tiles from a single or multiple image sets
	 * Next parameters are made separate to be able to modify them between different runs keeping clt_parameters
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
			final boolean use_vertical,
			final double min_strength,
			final double max_diff,
			final int max_iterations,
			final double max_coeff_diff,
			final double far_pull, //  = 0.2; // 1; //  0.5;
			EyesisCorrectionParameters.CLTParameters           clt_parameters,
			double [][] disp_strength,
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
		double [] disp_surface = new double[numTiles];
		ArrayList<Sample> samples_list = new ArrayList<Sample>();
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
				for (int nTile = 0; nTile < numTiles; nTile++){
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
			if ((debugLevel > 2) && (pass < 20)){
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

				(new showDoubleFloatArrays()).showArrays(dbg_img, tilesX, tilesY, true, "infinity_"+pass, titles);
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
	public double [][][] infinityCorrection(
//			final double min_strength0,
//			final double max_diff0,
//			final int max_iterations0,
//			final double max_coeff_diff0,
//			final double far_pull0, //  = 0.2; // 1; //  0.5;
			final boolean use_vertical,
			EyesisCorrectionParameters.CLTParameters           clt_parameters,
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
		for (Sample s: samples_list){
			int tileX = s.tile % tilesX;
			int tileY = s.tile / tilesX;
			double centerX = tileX * qc.tp.getTileSize() + qc.tp.getTileSize()/2;// - shiftX;
			double centerY = tileY * qc.tp.getTileSize() + qc.tp.getTileSize()/2;//- shiftY;
			double [][] centersXY_disp = qc.geometryCorrection.getPortsCoordinates(
					centerX,
					centerY,
					disp_strength[2 * s.series + 0][s.tile]/magic_coeff); // disparity
			double [][] centersXY_inf = qc.geometryCorrection.getPortsCoordinates(
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

			(new showDoubleFloatArrays()).showArrays(dbg_img, tilesX, tilesY, true, "AC_infinityCorrection", titles);

		}
		return inf_corr;
	}	  

	/**
	 * Correct channel mismatch (preserving disparity) using the same tiles as those for correcting disparity
	 * at infinity
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
	public double [][][] infinityMismatchCorrection(
			final boolean use_quadratic,
			final boolean use_vertical,
			final boolean use_disparity, // for infinity 
			EyesisCorrectionParameters.CLTParameters           clt_parameters,
			double [][] disp_strength,
			ArrayList<Sample> samples_list,
			int         tilesX,
			double      magic_coeff, // still not understood coefficient that reduces reported disparity value.  Seems to be around 8.5  
			int debugLevel)
	{
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
		double scale = 0.5/magic_coeff;
		double [][] dbg_xy = null;
		if (debugLevel > 0) {
			dbg_xy = new double [9][num_tiles];
		}
		for (Sample s: samples_list){
			int tileX = s.tile % tilesX;
			int tileY = s.tile / tilesX;
			double [] xy = new double[8]; // same as coefficients: x0,y0,x1,y1,x2,y2,x3,y3
			// Calculate x0,x1,x2,x3 and y0,y1,y2,y3 assuming x0+x1+x2+x3 = 0,y0+y1+y2+y3 = 0 and minimizing squares of errors
			// as each each 4: "dx0", "dx1", "dx2", "dx3" and "dy0", "dy1", "dy2", "dy3" are over-defined
			for (int dir = 0; dir < 2; dir++) { // 0 - X, 1 - Y
				double [] dxy = new double[4];
				for (int i = 0; i < 4; i++){
					dxy[i] = scale * disp_strength[indices_mismatch[dir][i] + (s.series * NUM_SLICES)][s.tile];
				}
				/*
				 *			    |-dy0   -dy1 -dy2   -dy3 |
				 *			B = |+dy0   -dy1      -2*dy3 |
				 *			    |+dy2 -2*dy1        -dy3 |
				 */					
			
				double [] B_arr = {
						-dxy[0]     -dxy[1] -dxy[2] -dxy[3],
						dxy[0]      -dxy[1]     -2 * dxy[3],
						dxy[2]  -2 * dxy[1]         -dxy[3]};
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
				dbg_xy[8][s.tile] += s.weight;
			}
				
			
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
			indx ++;
		}
		if (dbg_xy != null){
			for (int nTile = 0; nTile < num_tiles; nTile++){
				if (dbg_xy[8][nTile] > 0.0){
					for (int i = 0; i< 8; i++) {
						dbg_xy[i][nTile] /= dbg_xy[8][nTile];
					}
				} else {
					for (int i = 0; i< 8; i++) {
						dbg_xy[i][nTile] = Double.NaN;
					}
				}
			}
			String [] titles = {"x0", "y0", "x1", "y1", "x2", "y2", "x3","y3","weight"};
			(new showDoubleFloatArrays()).showArrays(
					dbg_xy,
					tilesX,
					tilesY,
					true,
					"xy_mismatch",
					titles);
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

	
	public double [][][] infinityMismatchCorrectionOld(
			final boolean use_quadratic,
			final boolean use_vertical,
			EyesisCorrectionParameters.CLTParameters           clt_parameters,
			double [][] disp_strength,
			ArrayList<Sample> samples_list,
			int         tilesX,
			double      magic_coeff, // still not understood coefficient that reduces reported disparity value.  Seems to be around 8.5  
			int debugLevel)
	{
		final int num_tiles =            disp_strength[0].length;
		final int tilesY =              num_tiles/tilesX;
		PolynomialApproximation pa = new PolynomialApproximation();
		double thresholdLin = 1.0E-20;  // threshold ratio of matrix determinant to norm for linear approximation (det too low - fail)
		double thresholdQuad = 1.0E-30; // threshold ratio of matrix determinant to norm for quadratic approximation (det too low - fail)
		//		final int [] indices_mismatch = {1,4,6,9}; //   "dy0", "dy1", "dx2", "dx3"
//		final int [] indices_mismatch = {2+1, 2+4, 2+6, 2+9}; //   "dy0", "dy1", "dx2", "dx3"
		final int [] indices_mismatch = {2+1, 2+3, 2+4, 2+6}; //   "dy0", "dy1", "dx2", "dx3"

		// use last generated samples_list;

		double [][][] mdata;
		if (use_vertical) {
			mdata = new double[samples_list.size()][3][];
		} else {
			mdata = new double [indices_mismatch.length][samples_list.size()][3];
		}
		int indx = 0;
		for (Sample s: samples_list){
			int tileX = s.tile % tilesX;
			int tileY = s.tile / tilesX;
			if (use_vertical) {
				// 2d optimization for 4 functions of x, y
				mdata[indx][0] = new double [2];
				mdata[indx][0][0] =  (2.0 * tileX)/tilesX - 1.0; // -1.0 to +1.0;
				mdata[indx][0][1] =  (2.0 * tileY)/tilesY - 1.0; // -1.0 to +1.0
				mdata[indx][1] = new double [indices_mismatch.length]; // 4
				for (int ip = 0; ip < indices_mismatch.length; ip++) {
					mdata[indx][1][ip] =  disp_strength[indices_mismatch[ip] + (s.series * NUM_SLICES)][s.tile];
				}

				mdata[indx][2] = new double [1];
				mdata[indx][2][0] = s.weight; //  disp_strength[2 * p.x + 1][p.y]; // strength

			} else {
				// 4 individual 1d optimization for functions of x only
				for (int n = 0; n < indices_mismatch.length; n++){
					mdata[n][indx][0] = (2.0 * tileX)/tilesX - 1.0; // -1.0 to +1.0;
					mdata[n][indx][1] = disp_strength[indices_mismatch[n] + (s.series * NUM_SLICES)][s.tile];
					mdata[n][indx][2] = s.weight; // disp_strength[2 * p.x + 1][p.y]; // strength
				}
			}
			indx ++;
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
				//			  disparity_approximation = new double[6];
				coeffs[n][5] = approx1d[0];
				coeffs[n][3] = approx1d[1];
				if (approx1d.length > 2){
					coeffs[n][0] = approx1d[2];
				}
			}
		}

		// convert to 8 sets of coefficient for px0, py0, px1, py1, ... py3.  
		double [][][] coeff_full = new double [4][2][6];
		double scale = 0.5/magic_coeff;
		for (int j = 0; j < coeffs[0].length; j++){
			coeff_full[0][0][j] = -coeffs[2][j] * scale;
			coeff_full[0][1][j] = -coeffs[0][j] * scale;
			coeff_full[1][0][j] = -coeffs[3][j] * scale;
			coeff_full[1][1][j] =  coeffs[0][j] * scale;
			coeff_full[2][0][j] =  coeffs[2][j] * scale;
			coeff_full[2][1][j] = -coeffs[1][j] * scale;
			coeff_full[3][0][j] =  coeffs[3][j] * scale;
			coeff_full[3][1][j] =  coeffs[1][j] * scale;
		}
		if (debugLevel > -1){
			double [] strength = new double [num_tiles];
			for (Sample s: samples_list){
				strength[s.tile] = s.weight;
			}
			int tilesX8 = (tilesX - 1) / 8 + 1;
			int tilesY8 = (tilesY - 1) / 8 + 1;
			int len8 = tilesX8*tilesY8;
			double [][] mismatch_8 = new double[11][len8];
			for (int i = 0; i < num_tiles; i++) {
				if (strength[i] > 0) {
					int tileX = i % tilesX;
					int tileY = i / tilesX;
					int tX8 = tileX/8;
					int tY8 = tileY/8;
					int indx8 = tY8*tilesX8+tX8;
					double tX =  (2.0 * tileX)/tilesX - 1.0; // -1.0 to +1.0;
					double tY =  (2.0 * tileY)/tilesY - 1.0; // -1.0 to +1.0
					double w = strength[i];
					for (int ip = 0; ip < 4; ip++) {
						mismatch_8[ip][indx8] += w * disp_strength[indices_mismatch[ip]][i];
					}
					mismatch_8[8][indx8] += w * tX;
					mismatch_8[9][indx8] += w * tY;
					mismatch_8[10][indx8] += w;
				}
			}
			for (int i = 0; i < len8; i++) {
				if (mismatch_8[10][i] > 0){
					for (int n = 0; n<4; n++){
						mismatch_8[n][i] /= mismatch_8[10][i]; 
					}
					mismatch_8[8][i] /= mismatch_8[10][i]; 
					mismatch_8[9][i] /= mismatch_8[10][i]; 
				} else {
					for (int n = 0; n<4; n++){
						mismatch_8[n][i] = Double.NaN; 
					}
					mismatch_8[8][i] = Double.NaN; 
					mismatch_8[9][i] = Double.NaN; 

				}
				for (int n = 0; n<4; n++){
					double tX = mismatch_8[8][i];
					double tY = mismatch_8[9][i];
					//f(x,y)=A*x^2+B*y^2+C*x*y+D*x+E*y+F
					mismatch_8[4+n][i] = (
							coeffs[n][0]*tX*tX+
							coeffs[n][1]*tY*tY+
							coeffs[n][2]*tX*tY+
							coeffs[n][3]*tX+
							coeffs[n][4]*tY+
							coeffs[n][5]);
				}
			}
			String [] titles = {"dy0","dy1","dx2","dx3","cy0","cy1","cx2","cx3","tx","ty","w"};
			showDoubleFloatArrays sdfa_instance = new showDoubleFloatArrays(); // just for debugging?
			sdfa_instance.showArrays(
					mismatch_8,
					tilesX8,
					tilesY8,
					true,
					"mismatch_8",
					titles);
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
			double w = disp_strength_in[1][nTile] - strength_floor;
			if (w > 0){
				if (strength_pow != 1.0) w = Math.pow(w, strength_pow);
				weight[nTile] = w;
			}
		}
		for (int tY = 0; tY < (tilesY - smplSide); tY++){
			for (int tX = 0; tX < (tilesX - smplSide); tX++){
				int num_in_sample = 0;
				boolean [] smpl_sel = new boolean [smplLen];
				double [] smpl_d =  new double [smplLen];
				double [] smpl_w =  new double [smplLen];
				double [][] smpl_mismatch = new double [NUM_SLICES-2][smplLen];
				for (int sy = 0; sy < smplSide; sy++){
					int y = tY + sy; //  - smpl_center;
					for (int sx = 0; sx < smplSide; sx++){
						int x = tX + sx; // - smpl_center;
						int indx = y * tilesX + x;
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
							System.out.println("**** this is a BUG in filterDisparityStrength() ****");
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
		    		float [][] fset = new float [3][];
//		    		fset[0] = (float[]) clt_mismatches_stack.getPixels(12 * num_scans + ns +1);
//		    		fset[1] = (float[]) clt_mismatches_stack.getPixels(13 * num_scans + ns +1); //
		    		scans[ns * NUM_ALL_SLICES + 0] = data[12 * num_scans + ns];
		    		scans[ns * NUM_ALL_SLICES + 1] = data[13 * num_scans + ns];
//		    		for (int i = 0; i < nTiles; i++){
//		    			scans[ns * NUM_ALL_SLICES + 0][i] = fset[0][i]; // disparity 
//		    			scans[ns * NUM_ALL_SLICES + 1][i] = fset[1][i]; // strength 
//		    		}
//		    		fset[0] = (float[]) clt_mismatches_stack.getPixels((2 * pair + 0) * num_scans + ns +1);
//		    		fset[1] = (float[]) clt_mismatches_stack.getPixels((2 * pair + 1) * num_scans + ns +1); //
//		    		fset[2] = (float[]) clt_mismatches_stack.getPixels((8 + pair   ) *  num_scans + ns +1);
	    			scans[ns * NUM_ALL_SLICES + pair * 3 + 2] = data[(2 * pair + 0) * num_scans + ns]; 
	    			scans[ns * NUM_ALL_SLICES + pair * 3 + 3] = data[(2 * pair + 1) * num_scans + ns];
	    			scans[ns * NUM_ALL_SLICES + pair * 3 + 4] = data[(8 + pair   ) *  num_scans + ns];

//		    		for (int i = 0; i < nTiles; i++){
//		    			scans[ns * NUM_ALL_SLICES + pair * 3 + 2][i] = fset[0][i]; // dx_i 
//		    			scans[ns * NUM_ALL_SLICES + pair * 3 + 3][i] = fset[1][i]; // dy_i
//		    			scans[ns * NUM_ALL_SLICES + pair * 3 + 4][i] = fset[2][i]; // str_i
//		    		}
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
				(new showDoubleFloatArrays()).showArrays(scans, tilesX, tilesY, true, "scans" , titles);
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
				(new showDoubleFloatArrays()).showArrays(scans, tilesX, tilesY, true, "scans" , titles);
		    }
		  return scans;
	  }
	  
	  
	  
		public double [][][] lazyEyeCorrection(
				final double min_strength_in,
				final double max_diff,
//				final double comp_strength_var,
				final int max_iterations,
				final double max_coeff_diff,
				final double far_pull, //  = 0.2; // 1; //  0.5;
				final double     strength_pow,
				final double     lazyEyeCompDiff, // clt_parameters.fcorr_disp_diff
				final int        lazyEyeSmplSide, //        = 2;      // Sample size (side of a square)
				final int        lazyEyeSmplNum, //         = 3;      // Number after removing worst (should be >1)
				final double     lazyEyeSmplRms, //         = 0.1;    // Maximal RMS of the remaining tiles in a sample
				final double     lazyEyeDispVariation, // maximal full disparity difference between tgh tile and 8 neighborxs
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
				EyesisCorrectionParameters.CLTParameters           clt_parameters,
				double [][]      scans_14,
				int              tilesX,
				double           magic_coeff, // still not understood coefficient that reduces reported disparity value.  Seems to be around 8.5  
				int debugLevel){
			final int num_scans = scans_14.length/NUM_ALL_SLICES;
		    final int num_tiles = scans_14[0].length;
		    final int tilesY = num_tiles/tilesX; 
		    final double [][] scans = new double [num_scans * NUM_SLICES][];
		    final int [] indices_14_10 = {0,1,2,3,5,6,8,9,11,12};
		    final double [][] comp_strength_rms = new double [num_scans][num_tiles];
		    for (int ns = 0; ns < num_scans; ns++){
		    	for (int i = 0; i < indices_14_10.length; i++){
		    		scans[ns * NUM_SLICES + i] = scans_14[ns * NUM_ALL_SLICES + indices_14_10[i]];
		    	}
		    }
//		    (new showDoubleFloatArrays()).showArrays(scans_14, tilesX, tilesY, true, "scans14");
//		    (new showDoubleFloatArrays()).showArrays(scans, tilesX, tilesY, true, "scans10");
		    
		    
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
		    for (int ns = 0; ns < num_scans; ns++){
		    	for (int nTile = 0; nTile < num_tiles; nTile++){
		    		double s1=0.0, s2=0.0;
		    		for (int pair = 0; pair <4; pair++){
		    			double s = scans_14[ns * NUM_ALL_SLICES + 4 + 3 * pair][nTile];
		    			if (pair == 0){
		    				s1 = s;
		    				s2 = s;
		    			} else {
		    				s1= Math.min(s, s1);
		    				s2= Math.max(s, s2);
		    			}
		    		}
		    		comp_strength_rms[ns][nTile] = s2 - s1;
		    	}
		    }
*/
	/* 
	 * None of comp_strength_rms methods works to detect potential outliers for horizontal/vertical features
	 */
		    
		    if (debugLevel > 0) {
		    	String [] titles = new String [num_scans]; 
			    for (int ns = 0; ns < num_scans; ns++){
		    		titles[ns] = "scan_" + ns;
			    }
				(new showDoubleFloatArrays()).showArrays(comp_strength_rms, tilesX, tilesY, true, "comp_strength_rms" , titles);
		    }

		    
			double[][] filtered_scans = filterDisparityStrength (
					scans, // final double[][] disp_strength_in,
					min_strength_in, // final double     strength_floor,
					strength_pow, // final double     strength_pow,
					lazyEyeSmplSide, // final int        smplSide, //        = 2;      // Sample size (side of a square)
					lazyEyeSmplNum, // final int        smplNum, //         = 3;      // Number after removing worst (should be >1)
					lazyEyeSmplRms, // final double     smplRms, //         = 0.1;    // Maximal RMS of the remaining tiles in a sample
					tilesX);// final int        tilesX);
		    
			if (debugLevel > -1) {
				System.out.println("lazyEyeCorrection() 1: removing tiles with residual disparity absoulte value > "+lazyEyeCompDiff);
			}

			double [][] combo_mismatch = new double [NUM_SLICES][num_tiles];
		    double [] combo_comp_rms = new double [num_tiles];
		    for (int ns = 0; ns < num_scans; ns++){
		    	for (int nTile = 0; nTile < num_tiles; nTile++) {
		    		double w = filtered_scans[ns * NUM_SLICES + 1][nTile];
		    		if (w > 0.0){
		    			double disp = filtered_scans[ns * NUM_SLICES + 0][nTile];
			    		if (Math.abs(disp) <= lazyEyeCompDiff) {
			    			for (int i = 2; i < NUM_SLICES; i++) if (i != 1){
			    				combo_mismatch[i][nTile] += filtered_scans[ns * NUM_SLICES + i][nTile] * w;
			    			}
		    				combo_mismatch[0][nTile] += (
		    						filtered_scans[ns * NUM_SLICES + 0][nTile]/clt_parameters.corr_magic_scale +
		    						clt_parameters.disp_scan_start + clt_parameters.disp_scan_step * ns)* w;
			    			combo_mismatch[1][nTile] += w;
			    			combo_comp_rms[nTile] += w * comp_strength_rms[ns][nTile];
			    		}
		    		}
		    	}
		    }
		    
	    	for (int nTile = 0; nTile < num_tiles; nTile++) {
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
	    	
	    	if (debugLevel > 100) {
	    		(new showDoubleFloatArrays()).showArrays(combo_comp_rms, tilesX, tilesY, "combo_comp_rms");
	    	}
	    	
			final TileNeibs tnImage = new TileNeibs(tilesX, tilesY); // num_tiles/tilesX);
	    	for (int nTile = 0; nTile < num_tiles; nTile++) if (combo_mismatch[1][nTile] > 0.0){
	    		double d = combo_mismatch[0][nTile];
	    		for (int dir = 0; dir <8; dir++){
	    			int nTile1 = tnImage.getNeibIndex(nTile, dir);
	    			if ((nTile1 >= 0) && (combo_mismatch[1][nTile1] > 0.0)){
	    				if (Math.abs(combo_mismatch[0][nTile1] - d) > lazyEyeDispVariation){
	    					combo_mismatch[1][nTile] = 0.0;
	    	    			for (int i = 0; i < NUM_SLICES; i++) if (i != 1){
	    	    				combo_mismatch[i][nTile] = Double.NaN;
	    	    			}
	    					break;
	    				}
	    			}
	    		}
	    	}			
	    	
		    if (debugLevel > 0) {
		    	String [] prefixes = {"disparity", "strength", "dx0", "dy0", "dx1", "dy1", "dx2", "dy2", "dx3", "dy3"};
				(new showDoubleFloatArrays()).showArrays(combo_mismatch, tilesX, combo_mismatch[0].length/tilesX, true, "combo_mismatch" , prefixes);
		    }
	    	
	    	combo_mismatch =  filterLazyEyePairs (
	    			combo_mismatch, // final double[][] samples_in,
	    			8,              // final int        smpl_side, // 8 x8 masked, 16x16 sampled
	    			0.25,           // final double     rms_max, TODO: find reasonable one mot critical?
	    			0.5,            // final double     frac_keep,
	    			5,              // final int        min_samples,
	    			true,           // final boolean    norm_center, // if there are more tiles that fit than minsamples, replace with a single equal weight
	    			tilesX);        // final int        tilesX);

		    if (debugLevel > 0) {
		    	String [] prefixes = {"disparity", "strength", "dx0", "dy0", "dx1", "dy1", "dx2", "dy2", "dx3", "dy3"};
				(new showDoubleFloatArrays()).showArrays(combo_mismatch, tilesX, combo_mismatch[0].length/tilesX, true, "filtered_mismatch" , prefixes);
		    }
	    	
	    	
	    	
	    	// extract infinity data to be processed as infinity
	    	
			double [][] inf_scan = new double [NUM_SLICES][];
			for (int i = 0; i < NUM_SLICES; i++){
				inf_scan[i] = scans[i];
			}
	    	
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

					(new showDoubleFloatArrays()).showArrays(dbg_img, tilesX, tilesY, true, "filtered_infinity_ds"); // , titles);

				}
			}
			
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
					(new showDoubleFloatArrays()).showArrays(dbg_img, tilesX, tilesY, true, "hist_filt_ds"); // , titles);
				}
			}
			// combine infinity and lasy eye scan data into a single array
			double [][] inf_and_ly = new double [2 * NUM_SLICES][];
			for (int i = 0; i < NUM_SLICES; i++){
				inf_and_ly[i] = inf_scan[i];
				inf_and_ly[i + NUM_SLICES] = combo_mismatch[i];
			}
			
			
			if (debugLevel > 0) {
		    	String [] prefixes = {"disparity", "strength", "dx0", "dy0", "dx1", "dy1", "dx2", "dy2", "dx3", "dy3"};
		    	String [] titles = new String [2 * NUM_SLICES]; 
		    	for (int i = 0; i < NUM_SLICES; i++){
		    		titles[i] = prefixes[i]+"-inf";
		    		titles[i + NUM_SLICES] = prefixes[i]+"-ly";
		    	}
		    	(new showDoubleFloatArrays()).showArrays(inf_and_ly, tilesX, tilesY, true, "inf_and_ly",titles);

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
		    	(new showDoubleFloatArrays()).showArrays(dbg_img, tilesX1, tilesY1, true, "inf_and_ly8",titles);
			}			
			
			// create list for infinity data
			ArrayList<Sample> inf_samples_list = selectInfinityTiles(
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
	    	for (Sample s: inf_samples_list) {
	    		total_weights[0] += s.weight;
	    	}
			
	    	for (int nTile = 0; nTile < num_tiles; nTile++) {
	    		total_weights[1]+= inf_and_ly[1 * NUM_SLICES + 1][nTile];
	    	}
			
			double [] weights = {
					inf_fraction *         (total_weights[0] + total_weights[1]) / total_weights[0], 
					(1.0 - inf_fraction) * (total_weights[0] + total_weights[1]) / total_weights[1], 
			};

			for (int ns = 0; ns <2; ns++) {
				for (int nTile = 0; nTile < num_tiles; nTile++) {
					inf_and_ly[ns * NUM_SLICES + 1][nTile] *= weights[ns];
				}
			}
	    	for (Sample s: inf_samples_list) {
	    		s.weight *= weights[0];
	    	}
			
			// Supplement list with the lazy eye scans data - use all tiles
			for (int nTile = 0; nTile < num_tiles; nTile++) {
				double w = inf_and_ly[1 * NUM_SLICES + 1][nTile];
				if (w > 0.0) {
					inf_samples_list.add(new Sample(1,nTile,w));
				}
			}
			
			if (debugLevel > -1) {
		    	double inf_weight = 0.0;
		    	for (Sample s: inf_samples_list) {
		    		inf_weight += s.weight;
		    	}
		    	System.out.println("lazyEyeCorrection(): number of all samples="+inf_samples_list.size()+", total weight = "+inf_weight);

		    }			
			
		    if (debugLevel > 0) {
		    	String [] prefixes = {"disparity", "strength", "dx0", "dy0", "dx1", "dy1", "dx2", "dy2", "dx3", "dy3"};
		    	String [] titles = new String [num_scans * NUM_SLICES]; 
			    for (int ns = 0; ns < num_scans; ns++){
			    	for (int i = 0; i < NUM_SLICES; i++){
			    		titles[ns * NUM_SLICES + i] = prefixes[i]+"_"+ns;
			    	}
			    }
				(new showDoubleFloatArrays()).showArrays(filtered_scans, tilesX, tilesY, true, "filtered_scans" , titles);
		    }
		    
		    
		    if (debugLevel > 0) {
		    	String [] prefixes = {"disparity", "strength", "dx0", "dy0", "dx1", "dy1", "dx2", "dy2", "dx3", "dy3"};
				(new showDoubleFloatArrays()).showArrays(combo_mismatch, tilesX, combo_mismatch[0].length/tilesX, true, "combo_mismatch" , prefixes);
		    }

		    double [][][] mismatch_corr_coefficiants = infinityMismatchCorrection(
		    		clt_parameters.fcorr_quadratic, // final boolean use_quadratic,
		    		true, // clt_parameters.fcorr_inf_vert,  // final boolean use_vertical,
					false,                          // final boolean use_disparity, // for infinity 
					clt_parameters,                 // EyesisCorrectionParameters.CLTParameters           clt_parameters,
					inf_and_ly,                     // double [][] disp_strength,
					inf_samples_list,               // ArrayList<Sample> samples_list,
					tilesX,                         // int         tilesX,
					magic_coeff,                    // double      , // still not understood coefficient that reduces reported disparity value.  Seems to be around 8.5  
					debugLevel);                    // int debugLevel)
		    if (debugLevel > -1) {
		    	System.out.println("===== lazyEyeCorrection(): correction coefficients =====");
		    	show_fine_corr(
		    			mismatch_corr_coefficiants,
		    			"mismatch_corr_coefficiants");
		    }
		    
		    return mismatch_corr_coefficiants;			
		}
	
		  public double [][] combineCltMismatches(
				  EyesisCorrectionParameters.CLTParameters clt_parameters,
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
				  EyesisCorrectionParameters.CLTParameters clt_parameters,
				  double [][]                              combo_data,
				  int tilesX,
				  int tilesY)
		  {
			  showDoubleFloatArrays sdfa_instance = new showDoubleFloatArrays();
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
				  EyesisCorrectionParameters.CLTParameters clt_parameters,
				  double [][]                              clt_mismatch,
				  int tilesX,
				  int tilesY)

		  {
			  showDoubleFloatArrays sdfa_instance = new showDoubleFloatArrays();
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

	  
	  public void process_fine_corr(
			  boolean dry_run,
			  EyesisCorrectionParameters.CLTParameters clt_parameters,
			  int debugLevel
			  ) {
		    final double disp_variation = 0.2; // 15; // 5; // 0.2 ?

	        ImagePlus imp_src = WindowManager.getCurrentImage();
	        if (imp_src==null){
	            IJ.showMessage("Error","12*n-layer file clt_mismatches is required");
	            return;
	        }
	        ImageStack clt_mismatches_stack= imp_src.getStack();
		    final int tilesX = clt_mismatches_stack.getWidth(); // tp.getTilesX();
		    final int tilesY = clt_mismatches_stack.getHeight(); // tp.getTilesY();
		    final int nTiles =tilesX * tilesY;
		    final int num_scans =  clt_mismatches_stack.getSize()/14;

		    final double [][] scans = new double [num_scans * NUM_SLICES][nTiles];
		    
		    for (int ns = 0; ns < num_scans; ns++){
		    	for (int pair = 0; pair < 4; pair++){
		    		float [][] fset = new float [2][];
		    		fset[0] = (float[]) clt_mismatches_stack.getPixels((2 * pair + 0) * num_scans + ns +1);
		    		fset[1] = (float[]) clt_mismatches_stack.getPixels((2 * pair + 1) * num_scans + ns +1); //
		    		for (int i = 0; i < nTiles; i++){
		    			scans[ns * NUM_SLICES + pair * 2 + 2][i] = fset[0][i]; // dxi 
		    			scans[ns * NUM_SLICES + pair * 2 + 3][i] = fset[1][i]; // dyi
		    		}
		    	}
	    		float [][] fset = new float [2][];
	    		fset[0] = (float[]) clt_mismatches_stack.getPixels(12 * num_scans + ns +1);
	    		fset[1] = (float[]) clt_mismatches_stack.getPixels(13 * num_scans + ns +1); //
	    		for (int i = 0; i < nTiles; i++){
	    			scans[ns * NUM_SLICES + 0][i] = fset[0][i]; // disparity 
	    			scans[ns * NUM_SLICES + 1][i] = fset[1][i]; // strength 
	    		}
		    }
		    if (debugLevel > -1) {
		    	String [] prefixes = {"disparity", "strength", "dx0", "dy0", "dx1", "dy1", "dx2", "dy2", "dx3", "dy3"};
		    	String [] titles = new String [num_scans * NUM_SLICES]; 
			    for (int ns = 0; ns < num_scans; ns++){
			    	for (int i = 0; i < NUM_SLICES; i++){
			    		titles[ns * NUM_SLICES + i] = prefixes[i]+"_"+ns;
			    	}
			    }
				(new showDoubleFloatArrays()).showArrays(scans, tilesX, scans[0].length/tilesX, true, "scans" , titles);
		    }
		    
//  		public double     fcorr_disp_diff =   1.5;   // consider only tiles with absolute residual disparity lower than
// separate here from reading image		    
		    
		    final int num_tiles = scans[0].length;
		    
		    
			double[][] filtered_scans = filterDisparityStrength (
					scans, // final double[][] disp_strength_in,
					clt_parameters.fcorr_inf_strength, // final double     strength_floor,
					clt_parameters.inf_str_pow, // final double     strength_pow,
					clt_parameters.inf_smpl_side, // final int        smplSide, //        = 2;      // Sample size (side of a square)
					clt_parameters.inf_smpl_num, // final int        smplNum, //         = 3;      // Number after removing worst (should be >1)
					clt_parameters.inf_smpl_rms, // final double     smplRms, //         = 0.1;    // Maximal RMS of the remaining tiles in a sample
					tilesX);// final int        tilesX);
		    
//  		public double     fcorr_disp_diff =   1.5;   // consider only tiles with absolute residual disparity lower than
			if (debugLevel > -1) {
				System.out.println("process_fine_corr() 2: removing tile with residual disparity absoulte value > "+ clt_parameters.fcorr_disp_diff);
			}
		    for (int ns = 0; ns < num_scans; ns++){
		    	for (int i = 0; i < num_tiles; i++){
		    		double disp = filtered_scans[ns * NUM_SLICES + 0][i];
		    		if (Math.abs(disp) > clt_parameters.fcorr_disp_diff) {
			    		filtered_scans[ns * NUM_SLICES + 1][i] = 0.0;
		    		}
		    	}
		    }
		    double [][] combo_mismatch = new double [NUM_SLICES][num_tiles];
		    for (int ns = 0; ns < num_scans; ns++){
		    	for (int nTile = 0; nTile < num_tiles; nTile++) {
		    		double w = filtered_scans[ns * NUM_SLICES + 1][nTile];
		    		if (w > 0.0){
		    			double disp = filtered_scans[ns * NUM_SLICES + 0][nTile];
			    		if (Math.abs(disp) <= clt_parameters.fcorr_disp_diff) {
			    			for (int i = 2; i < NUM_SLICES; i++) if (i != 1){
			    				combo_mismatch[i][nTile] += filtered_scans[ns * NUM_SLICES + i][nTile] * w;
			    			}
		    				combo_mismatch[0][nTile] += (
		    						filtered_scans[ns * NUM_SLICES + 0][nTile]/clt_parameters.corr_magic_scale +
		    						clt_parameters.disp_scan_start + clt_parameters.disp_scan_step * ns)* w;
			    			combo_mismatch[1][nTile] += w;
			    		}
		    		}
		    	}
		    }
		    
	    	for (int nTile = 0; nTile < num_tiles; nTile++) {
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
		    
			final TileNeibs tnImage = new TileNeibs(tilesX, tilesY); // num_tiles/tilesX);
	    	for (int nTile = 0; nTile < num_tiles; nTile++) if (combo_mismatch[1][nTile] > 0.0){
	    		double d = combo_mismatch[0][nTile];
	    		for (int dir = 0; dir <8; dir++){
	    			int nTile1 = tnImage.getNeibIndex(nTile, dir);
	    			if ((nTile1 >= 0) && (combo_mismatch[1][nTile1] > 0.0)){
	    				if (Math.abs(combo_mismatch[0][nTile1] - d) > disp_variation){
	    					combo_mismatch[1][nTile] = 0.0;
	    	    			for (int i = 0; i < NUM_SLICES; i++) if (i != 1){
	    	    				combo_mismatch[i][nTile] = Double.NaN;
	    	    			}
	    					break;
	    				}
	    			}
	    		}
	    	}			
		    
	    	// extract infinity data to be processed as infinity
	    	
			double [][] inf_scan = new double [NUM_SLICES][];
			for (int i = 0; i < NUM_SLICES; i++){
				inf_scan[i] = scans[i];
			}
	    	
			// Need to filter first!
			
			
			ArrayList<Sample> inf_samples_list = selectInfinityTiles(
					clt_parameters.fcorr_inf_vert,// final boolean use_vertical,
					0.0, // any > 0.0
					clt_parameters.fcorr_inf_diff, // max_diff, //clt_parameters.fcorr_inf_diff
					clt_parameters.inf_iters, // max_iterations, // clt_parameters.inf_iters
					clt_parameters.inf_final_diff, // max_coeff_diff, // clt_parameters.inf_final_diff
					clt_parameters.inf_far_pull, // far_pull, // clt_parameters.inf_far_pull,  = 0.2; // 1; //  0.5;
					clt_parameters,
					inf_scan,
					tilesX,
					clt_parameters.corr_magic_scale, // magic_coeff, // still not understood coefficient that reduces reported disparity value.  Seems to be around 8.5  
					debugLevel);
			
		    if (debugLevel > -1) {
		    	double inf_weight = 0.0;
		    	for (Sample s: inf_samples_list) {
		    		inf_weight += s.weight;
		    	}
		    	System.out.println("process_fine_corr(): number of infinity samples="+inf_samples_list.size()+", total weight = "+inf_weight);

		    }			
	    	
			
		    if (debugLevel > -1) {
		    	String [] prefixes = {"disparity", "strength", "dx0", "dy0", "dx1", "dy1", "dx2", "dy2", "dx3", "dy3"};
		    	String [] titles = new String [num_scans * NUM_SLICES]; 
			    for (int ns = 0; ns < num_scans; ns++){
			    	for (int i = 0; i < NUM_SLICES; i++){
			    		titles[ns * NUM_SLICES + i] = prefixes[i]+"_"+ns;
			    	}
			    }
				(new showDoubleFloatArrays()).showArrays(filtered_scans, tilesX, tilesY, true, "filtered_scans" , titles);
		    }
		    
		    
		    if (debugLevel > -1) {
		    	String [] prefixes = {"disparity", "strength", "dx0", "dy0", "dx1", "dy1", "dx2", "dy2", "dx3", "dy3"};
				(new showDoubleFloatArrays()).showArrays(combo_mismatch, tilesX, combo_mismatch[0].length/tilesX, true, "combo_mismatch" , prefixes);
		    }
		    
		    
	  }	  
}
