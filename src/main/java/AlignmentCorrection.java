import java.util.ArrayList;

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

public class AlignmentCorrection {
	QuadCLT qc;
	AlignmentCorrection (QuadCLT qc){
		this.qc = qc;
	}

	public double [][][] infinityCorrection(
			final double min_strength_in,
			final double max_diff,
			final int max_iterations,
			final double max_coeff_diff,
			final double far_pull, //  = 0.2; // 1; //  0.5;
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
			double [][] disp_strength_in,
			int         tilesX,
			double      magic_coeff, // still not understood coefficient that reduces reported disparity value.  Seems to be around 8.5  
			int debugLevel)
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
		return infinityCorrection(
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
	}


	/**
	 * Calculate quadratic polynomials for each subcamera X/Y correction to match disparity = 0 at infinity
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
	public double [][][] infinityCorrection(
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
		//		  double [] corr_weight = new double [numTiles]; // reduce weight for far tiles (by more than range farther than surface)
		class Sample{
			int series;
			int tile;
			double weight;
			Sample(int series, int tile, double weight) {
				this.series = series;
				this.tile =  tile;
				this.weight = weight;

			}
		}
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
			for (int num_set = 0; num_set < disp_strength.length/2; num_set++){
				int disp_index = 2 * num_set;
				int str_index = 2 * num_set+1;
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
			if (clt_parameters.fcorr_inf_vert) {
				mdata = new double[samples_list.size()][3][];
			} else {
				mdata = new double [1][samples_list.size()][3];
			}

			int indx = 0;
			for (Sample s: samples_list){
				int tileX = s.tile % tilesX;
				int tileY = s.tile / tilesX;
				if (clt_parameters.fcorr_inf_vert) {
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
			if ((debugLevel > -1) && (pass < 20)){
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
			if (clt_parameters.fcorr_inf_vert){
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

		// use last generated samples_list;

		double [][][] mdata;
		if (clt_parameters.fcorr_inf_vert) {
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
			if (clt_parameters.fcorr_inf_vert) {
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
		if (clt_parameters.fcorr_inf_vert){
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

			(new showDoubleFloatArrays()).showArrays(dbg_img, tilesX, tilesY, true, "infinityCorrection", titles);

		}
		return inf_corr;
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
		if (disp_strength_in.length > 2){
			final double [][] disp_strength = new double [disp_strength_in.length][];
			for (int i = 0; i < disp_strength_in.length; i+=2){
				double [][] ds = {disp_strength_in[i],disp_strength_in[i+1]};
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
				disp_strength[i] =   ds_rslt[0];
				disp_strength[i+1] = ds_rslt[1];
			}
			return disp_strength;
		}
		final int num_tiles = disp_strength_in[0].length;
		int tilesY = num_tiles/tilesX;
		final double [][] disp_strength = new double [2][num_tiles];
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
					for (int sY = 0; sY < smpl_side; sY++){
						int y = tY + sY;
						for (int sX = 0; sX < smpl_side; sX++){
							int x = tX + sX;
							int nTile = y * tilesX + x;
							double w = disp_strength_in[1][nTile];
							if (w > 0.0){
								if (Math.abs(disp_strength_in[0][nTile] - disp_inf) < max_diff) {
									sw += disp_strength_in[1][nTile];
									sdw +=disp_strength_in[0][nTile] * disp_strength_in[1][nTile];
									num_samples++;
								}
							}
						}
					}
					if ((num_samples > min_samples) && (sw > 0.0)){
						int nTile = tY * tilesX + tX;
						disp_strength[0][nTile] = sdw/sw; // weighted average disparity
						disp_strength[1][nTile] = 1.0; // equal weight
					}
				} else {
					// just center tile - fits or not
					int nTile = tY * tilesX + tX;
					double w = disp_strength_in[1][nTile];
					if (w > 0.0){
						if (Math.abs(disp_strength_in[0][nTile] - disp_inf) < max_diff) {
							disp_strength[0][nTile] = disp_strength_in[0][nTile];
							disp_strength[1][nTile] = disp_strength_in[1][nTile];
						}
					}
					
/*					
					for (int sY = 0; sY < smpl_side; sY++){
						int y = tY + sY;
						for (int sX = 0; sX < smpl_side; sX++){
							int x = tX + sX;
							int nTile = y * tilesX + x;
							double w = disp_strength_in[1][nTile];
							if (w > 0.0){
								if (Math.abs(disp_strength_in[0][nTile] - disp_inf) < max_diff) {
									disp_strength[0][nTile] = disp_strength_in[0][nTile];
									disp_strength[1][nTile] = disp_strength_in[1][nTile];
								}
							}
						}
					}
*/					
				}
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
		if (disp_strength_in.length > 2){
			final double [][] disp_strength = new double [disp_strength_in.length][];
			for (int i = 0; i < disp_strength_in.length; i+=2){
				double [][] ds = {disp_strength_in[i],disp_strength_in[i+1]};
				double [][] ds_rslt = filterDisparityStrength (
						ds,
						strength_floor,
						strength_pow,
						smplSide, //        = 2;      // Sample size (side of a square)
						smplNum, //         = 3;      // Number after removing worst (should be >1)
						smplRms, //         = 0.1;    // Maximal RMS of the remaining tiles in a sample
						tilesX);
				disp_strength[i] =   ds_rslt[0];
				disp_strength[i+1] = ds_rslt[1];
			}
			return disp_strength;
		}

		final int num_tiles = disp_strength_in[0].length;
		int tilesY = num_tiles/tilesX;
		final double [][] disp_strength = new double [2][num_tiles];
		final int index_shift  = (smplSide /2) * (tilesX + 1) ; // diagonal shift
		final int smplLen = smplSide*smplSide;

		final double [] weight = new double [num_tiles]; // modified weight
		final double [] disp = disp_strength_in[0];
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
							num_in_sample ++;
						}
					}
				}
				if (num_in_sample >= smplNum){ // try, remove worst
					// calculate 
					double sd=0.0, sd2 = 0.0, sw = 0.0;
					for (int i = 0; i < smplLen; i++) if (smpl_sel[i]) {
						double dw = smpl_d[i] * smpl_w[i];
						sd += dw;
						sd2 += dw * smpl_d[i];
						sw +=       smpl_w[i];
					}
					// remove worst, update sd2, sd and sw
					while ((num_in_sample > smplNum) && (sw > 0)){ // try, remove worst
						double d_mean = sd/sw;
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

}
