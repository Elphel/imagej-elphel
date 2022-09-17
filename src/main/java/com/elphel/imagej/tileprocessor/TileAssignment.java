package com.elphel.imagej.tileprocessor;
/**
 **
 ** TileAssignment - handle tile surfaces
 **
 ** Copyright (C) 2017 Elphel, Inc.
 **
 ** -----------------------------------------------------------------------------**
 **
 **  TileAssignment.java is free software: you can redistribute it and/or modify
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
import java.awt.Point;
import java.util.HashMap;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.concurrent.atomic.AtomicInteger;

import com.elphel.imagej.cameras.CLTParameters;
import com.elphel.imagej.common.DoubleGaussianBlur;
import com.elphel.imagej.common.ShowDoubleFloatArrays;

//import ij.IJ;


public class TileAssignment {
	private TileSurface ts;
	private CLTPass3d p3d; // last/ combined FPGA data
	private CLTParameters clt_parameters;
	private double [][][] dispStrength; // indexed as a surface (full supertiles), not as image
	private boolean [] valid_ml;
	private double [][] tile_tones;
	private double [][][] tone_diff_weight; // for each surface (not image) tile, for each of 8 directions,
	//a pair of normalized tone difference and weight (both weight and diff == 0.0 for connections to/from undefined tiles
	private double kR = 1.0; // relative weight of red to green ratio
	private double kB = 1.0; // relative weight of blue to green ratio
	private double fatZero = 0.01; // for colors, normalized to 0..1.0 range
	private int surfTilesX;
	private int surfTilesY;
	private int imgTilesX;
	private int imgTilesY;
	private TileNeibs tnImage;
	private TileNeibs tnSurface;
	private boolean [][][] valid_surf; // per measured layer, per surface tile, per surface layer

	// from clt_parameters
	private double dispNorm;
	private TACosts cost_coeff;
	private double minFgBg; //  =    0.1;  // Minimal foreground/ background separation to look for weak FG edge
	private double minFgEdge; //  =  0.2;  // Minimal foreground edge strength (stronger edges will have proportionally smaller costs)
	private double minColSep; //  =  0.05; // Minimal surface separation that requires color change
	private double minColDiff; //  = 0.01; // Minimal color variation (larger proportionally reduces cost)
	private double dispOutlier; //             = 1.0;    // Disparity difference limit (to handle outliers)
	private double strengthDiffPwr; //             = 0.25;   // Strength power when calculating disparity error
	private double strengthBestPwr; //             = 0.0;    // Strength power when calculating disparity error over best
	private double strengthDiff9Pwr; //            = 0.5;    // Strength power when calculating disparity error for group of 9
	private double taColSigma; //           = 1.5;    // Gaussian sigma to blur color difference between tiles along each direction
	private double taColFraction; //        = 0.3;    // Relative amount of the blurred color difference in the mixture


	private double shrinkWeakFgnd = 0.5; //  0.5; //reduce cost of multiple weak_fgnd links of the same tile
	private double shrinkColor = 0.5; // 0.0;    //reduce cost of surface transitions w/o color change

	static double [] NEIB_WEIGHTS = {0.5,0.25, 0.5, 0.25, 0.5,0.25, 0.5, 0.25, 1.0};

	static int INDEX_R = 0;
	static int INDEX_B = 1;
	static int INDEX_G = 2;
	static int INDEX_A = 3;

	public boolean mutual_weak_fgnd = false; // set to true when using 3x3 grid, false for 5x5;


	public class TACosts{
		public double empty;     // Cost of a tile that is not assigned
		public double nolink;    // Cost of a tile not having any neighbor in particular direction
		public double swtch;     // Cost of a tile switching to a neighbor that does not have a link
		public double color;     // Cost of a tile switching to a disconnected neighbor divided by a color mismatch
		public double diff;      // Cost of a weighted normalized tile disparity error
		public double diff_best; // Cost of a weighted normalized tile disparity error above best surface
		public double diff9;     // Cost of a weighted normalized tile disparity error for tile and 8 neighbors (DC)
  		public double weak_fgnd; // Cost of a weak foreground edge
		public double flaps;     // Cost of using supertile "flaps" (not in the center 8x8 tiles area)
  		public double ml_mismatch;    // Cost of a measurement layer not having same layer in the same location or near

		public TACosts (int mode)
		{
			this.empty =       1.0/0.35833; // 0.71715;
			this.nolink =      1.0/0.83739; // 0.95952; // 1.50705;
			this.swtch =       1.0/0.07604; // 0.05172; // 0.59474;
			this.color =       1.0/6.69936; // 3.34968; // 0.21458; // 0.19294; // 2.25763;
			this.diff =        1.0/0.11879; // 0.11039; // 1.94213;
			this.diff_best =   1.0/0.02831; // 0.00628; // 0.06731;
			this.diff9 =       1.0/0.04064; // 0.01350; // 1.09087;
			this.weak_fgnd =   1.0/2.355492; // 1.177746; // 0.02172; // 0.01726; // 0.22250;
			this.flaps =       1.0/0.1; // 0.00056; // 0.07229;
			this.ml_mismatch = 1.0;
			// ml_mismatch not yet implemented - is it needed - maybe some see-through data appears on one layer only
		}
		/*
empty=     0.00000 nolink= 1.50705 swtch=     0.59474 color= 2.25763 diff=        1.94213
diff_best= 0.06731 diff9=  1.09087 weak_fgnd= 0.22250 flaps= 0.07229 ml_mismatch= 0.00000
		 */
		public TACosts (CLTParameters clt_parameters)
		{
			this.empty =       clt_parameters.taEnEmpty?    clt_parameters.taCostEmpty    : 0.0;
			this.nolink =      clt_parameters.taEnNoLink?   clt_parameters.taCostNoLink   : 0.0;
			this.swtch =       clt_parameters.taEnSwitch?   clt_parameters.taCostSwitch   : 0.0;
			this.color =       clt_parameters.taEnColor?    clt_parameters.taCostColor    : 0.0;
			this.diff =        clt_parameters.taEnDiff?     clt_parameters.taCostDiff     : 0.0;
			this.diff_best =   clt_parameters.taEnDiffBest? clt_parameters.taCostDiffBest : 0.0;
			this.diff9 =       clt_parameters.taEnDiff9?    clt_parameters.taCostDiff9    : 0.0;
			this.weak_fgnd =   clt_parameters.taEnWeakFgnd? clt_parameters.taCostWeakFgnd : 0.0;
			this.flaps =       clt_parameters.taEnFlaps?    clt_parameters.taCostFlaps    : 0.0;
			this.ml_mismatch = clt_parameters.taEnMismatch? clt_parameters.taCostMismatch : 0.0;
			// ml_mismatch not yet implemented - is it needed - maybe some see-through data appears on one layer only
		}

		public TACosts(){}
		public TACosts(double [] vals)
		{
			set(vals);
		}

		public void set(double [] vals)
		{
			this.empty =       vals[0];
			this.nolink =      vals[1];
			this.swtch =       vals[2];
			this.color =       vals[3];
			this.diff =        vals[4];
			this.diff_best =   vals[5];
			this.diff9 =       vals[6];
			this.weak_fgnd =   vals[7];
			this.flaps =       vals[8];
			this.ml_mismatch = vals[9];
		}


		public double [] toArray(){
			double [] rslt = {
					this.empty,
					this.nolink,
					this.swtch,
					this.color,
					this.diff,
					this.diff_best,
					this.diff9,
					this.weak_fgnd,
					this.flaps,
					this.ml_mismatch
			};
			return rslt;
		}

		public String [] getTitles()
		{
			String [] titles = {
					"empty",
					"nolink",
					"swtch",
					"color",
					"diff",
					"diff_best",
					"diff9",
					"weak_fgnd",
					"flaps",
					"ml_mismatch"};
			return titles;
		}

		@Override
		public String toString()
		{
			String s = "";
			s+= String.format("empty=    %8.5f nolink=%8.5f swtch=    %8.5f color=%8.5f diff=       %8.5f\n",empty, nolink, swtch, color, diff);
			s+= String.format("diff_best=%8.5f diff9= %8.5f weak_fgnd=%8.5f flaps=%8.5f ml_mismatch=%8.5f\n",diff_best, diff9, weak_fgnd, flaps, ml_mismatch);
			return s;
		}

		public double dotProd (TACosts costs){
			double [] these = toArray();
			double [] other = costs.toArray();
			double prod = 0.0;
			for (int i = 0; i < these.length; i++){
				prod += these[i]*other[i];
			}
			return prod;
		}

		public void add (TACosts costs){
			double [] these = toArray();
			double [] other = costs.toArray();
			for (int i = 0; i < these.length; i++){
				these[i] += other[i];
			}
			set (these);
		}

		public void mul (TACosts costs){
			double [] these = toArray();
			double [] other = costs.toArray();
			for (int i = 0; i < these.length; i++){
				these[i] *= other[i];
			}
			set (these);
		}

		public void max (TACosts costs){
			double [] these = toArray();
			double [] other = costs.toArray();
			for (int i = 0; i < these.length; i++){
				these[i] = Math.max(these[i],other[i]);
			}
			set (these);
		}

		public void scale (double a){
			double [] these = toArray();
			for (int i = 0; i < these.length; i++){
				these[i] *= a;
			}
			set (these);
		}
	}


	/**
	 *
	 * @param ts
	 * @param p3d
	 * @param tile_sel
	 * @param dispStrength per measurement layer, combined disparity and strength array ([num_ml [2][])
	 * @param kR
	 * @param kB
	 * @param fatZero
	 */
	public TileAssignment (
			CLTParameters clt_parameters,
			TileSurface ts,
			CLTPass3d p3d,
			double [][][] dispStrength,
			double kR,
			double kB,
			double fatZero)
	{
		this.ts = ts;
		this.p3d = p3d;
		copyParams(clt_parameters);
		this.clt_parameters = clt_parameters; // replace with copying specific ones
		valid_ml = new boolean [dispStrength.length];
		for (int ml = 0; ml < dispStrength.length; ml++){
			valid_ml[ml] = dispStrength[ml] != null;
		}
		if (!Double.isNaN(kR))      this.kR = kR;
		if (!Double.isNaN(kB))      this.kB = kB;
		if (!Double.isNaN(fatZero)) this.fatZero = fatZero;
		this.surfTilesX = ts.getSTilesX() * ts.getSuperTileSize();
		this.surfTilesY = ts.getSTilesY() * ts.getSuperTileSize();
		this.imgTilesX = ts.getImageTilesX();
		this.imgTilesY = ts.getImageTilesY();
		this.tnImage =   new TileNeibs(imgTilesX, imgTilesY);
		this.tnSurface = new TileNeibs(surfTilesX, surfTilesY);
		setDispStrength(dispStrength);
		setTones(p3d);
	}


	public void copyParams(CLTParameters clt_parameters)
	{
		this.dispNorm =          clt_parameters.plDispNorm;
		this.minFgBg =           clt_parameters.taMinFgBg;
		this.minFgEdge =         clt_parameters.taMinFgEdge;
		this.minColSep =         clt_parameters.taMinColSep;
		this.minColDiff =        clt_parameters.taMinColDiff;

		this.dispOutlier =       clt_parameters.taOutlier;
		this.strengthDiffPwr =   clt_parameters.taDiffPwr;
		this.strengthBestPwr =   clt_parameters.taBestPwr;
		this.strengthDiff9Pwr =  clt_parameters.taDiff9Pwr;

		this.taColSigma =        clt_parameters.taColSigma;
		this.taColFraction =     clt_parameters.taColFraction;

		this.cost_coeff = new TACosts (clt_parameters);
		this.cost_coeff.mul(new TACosts (0)); // make average ~=1.0 for each used component

	}
	public void setDispStrength(
			double [][][] ds)
	{
		this.dispStrength = new double [ds.length][][];
		for (int ml = 0; ml < ds.length; ml++) if (ds[ml] != null){
			this.dispStrength[ml] = new double [ surfTilesX * surfTilesY][];
			for (int nTile = 0; nTile < ds[ml][0].length; nTile++){
				int nSurfTile = ts.getSurfaceTileIndex(nTile);
				this.dispStrength[ml][nSurfTile] = new double[2]; // ds[ml][nTile];
				this.dispStrength[ml][nSurfTile][0] = ds[ml][0][nTile];
				this.dispStrength[ml][nSurfTile][1] = ds[ml][1][nTile];
			}
		}
	}

	public int getSurfTilesX(){
		return surfTilesX;
	}
	public int getSurfTilesY(){
		return surfTilesY;
	}
//	private int surfTilesX;
//	private int surfTilesY;


	public void setTones(CLTPass3d p3d)
	{
		double [][] tt = p3d.getTileRBGA(4);

		this.tile_tones = new double [surfTilesX * surfTilesY][tt.length];
		for (int nTile = 0; nTile < tt[0].length; nTile++){
			int nSurfTile = ts.getSurfaceTileIndex(nTile);
			for (int j = 0; j < tt.length; j++){
				this.tile_tones[nSurfTile][j] = tt[j][nTile];
			}
		}

		this.tone_diff_weight = new double [surfTilesX * surfTilesY][8][2];
//		final int numTiles = imgTilesX * imgTilesY;
		final Thread[] threads = ImageDtt.newThreadArray(ts.getThreadsMax());
		final AtomicInteger ai = new AtomicInteger(0);
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				@Override
				public void run() {
					for (int nSurfTile = ai.getAndIncrement(); nSurfTile < tone_diff_weight.length; nSurfTile = ai.getAndIncrement()) {
						if (nSurfTile >= 0){
							for (int dir = 0; dir < 4; dir++){
								int nSurfTile1 = tnSurface.getNeibIndex(nSurfTile, dir);
								if (nSurfTile1 >= 0){
									tone_diff_weight[nSurfTile][dir] =
											getDiffWeight (
													tile_tones[nSurfTile],
													tile_tones[nSurfTile1]);
									tone_diff_weight[nSurfTile1][(dir + 4) % 8] = tone_diff_weight[nSurfTile][dir];
								}
							}
						}
					}
				}
			};
		}
		ImageDtt.startAndJoin(threads);
	}

	public void blurMixTones(
			final boolean use_sqrt,
			final boolean weighted) // blur weighted
	{
		this.tone_diff_weight = blurMixTones(
				this.tone_diff_weight, // final double [][][] tone_diffs,
				this.taColSigma, // final double sigma,
				this.taColFraction, // final double fraction,
				use_sqrt,
				weighted); // final boolean weighted)
	}

	public double [][][] blurMixTones(
			final double [][][] tone_diffs,
			final double sigma,
			final double fraction,
			final boolean use_sqrt,
			final boolean weighted) // blur weighted
	{
		final int num_tiles = surfTilesX * surfTilesY;
		final double [][] dir_blur =   new double [4][num_tiles];
		final double [][] dir_weight = new double [4][num_tiles];
		final Thread[] threads = ImageDtt.newThreadArray(ts.getThreadsMax());
		final AtomicInteger ai = new AtomicInteger(0);
		final double [][][] mixed_diff_weight = new double [surfTilesX * surfTilesY][8][2];

		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				@Override
				public void run() {
					for (int dir = ai.getAndIncrement(); dir < 4; dir = ai.getAndIncrement()) {
						for (int nSurfTile = 0; nSurfTile < num_tiles; nSurfTile++){
							double w = tone_diffs[nSurfTile][dir][1];
							double dw = tone_diffs[nSurfTile][dir][0] ;
							if (use_sqrt) dw = Math.sqrt(dw);
							if (weighted) dw *= w;
							dir_blur[dir][nSurfTile] = dw;
							dir_weight[dir][nSurfTile] = w;
						}
						DoubleGaussianBlur gb =new DoubleGaussianBlur();
						gb.blurDouble(dir_blur[dir],   surfTilesX, surfTilesY, sigma, sigma, 0.01);
						gb.blurDouble(dir_weight[dir], surfTilesX, surfTilesY, sigma, sigma, 0.01);
						if (weighted){
							for (int nSurfTile = 0; nSurfTile < num_tiles; nSurfTile++){
								if (dir_weight[dir][nSurfTile] != 0.0) dir_blur[dir][nSurfTile] /= dir_weight[dir][nSurfTile];
							}
						}
						int rdir = (dir + 4) % 8;
						double rfract = 1.0 - fraction;
						for (int nSurfTile = 0; nSurfTile < num_tiles; nSurfTile++){
							int nSurfTile1 = tnSurface.getNeibIndex(nSurfTile, dir);
							if (nSurfTile1 >= 0){
								double d;
								if (use_sqrt){
									d = rfract * Math.sqrt(tone_diffs[nSurfTile][dir][0]) + fraction * dir_blur[dir][nSurfTile];
									d *= d;
								} else {
									d =rfract * tone_diffs[nSurfTile][dir][0] + fraction * dir_blur[dir][nSurfTile];
								}
								mixed_diff_weight[nSurfTile][dir][0] = d;
								mixed_diff_weight[nSurfTile][dir][1] = rfract * tone_diffs[nSurfTile][dir][1] + fraction * dir_weight[dir][nSurfTile];
								mixed_diff_weight[nSurfTile1][rdir][0] = mixed_diff_weight[nSurfTile][dir][0];
								mixed_diff_weight[nSurfTile1][rdir][1] = mixed_diff_weight[nSurfTile][dir][1];
							}
						}
					}
				}
			};
		}
		ImageDtt.startAndJoin(threads);
		return mixed_diff_weight;
	}



	public void showToneDiffWeights3(String prefix){
		String [] titles = {"diffs","sqrt","centers","weights"};
		double [][] img_data = new double[titles.length][9 * surfTilesX * surfTilesY];
		TileNeibs tnSurface3 = new TileNeibs( 3 * surfTilesX, 3 * surfTilesY);
		for (int nSurfTile = 0; nSurfTile <  tone_diff_weight.length;  nSurfTile++){
			int [] txy = tnSurface.getXY(nSurfTile);
			int [] txy3 = {3 * txy[0] + 1, 3 * txy[1] + 1};
			double sdw = 0.0, sw = 0.0;
			int num_neibs = 0;
			for (int dir = 0; dir < 8; dir++){
				sw  += tone_diff_weight[nSurfTile][dir][1];
				sdw += tone_diff_weight[nSurfTile][dir][0] * tone_diff_weight[nSurfTile][dir][1];
				int nSurfTile3 = tnSurface3.getNeibIndex(tnSurface3.getIndex(txy3[0],txy3[1]),dir);
				img_data[0][nSurfTile3] = tone_diff_weight[nSurfTile][dir][0];
				img_data[1][nSurfTile3] = Math.sqrt(tone_diff_weight[nSurfTile][dir][0]);
				img_data[3][nSurfTile3] = tone_diff_weight[nSurfTile][dir][1];
				if (tone_diff_weight[nSurfTile][dir][1] > 0.0) num_neibs ++;
			}
			if (sw > 0.0){
				sdw /= sw;
			}
			if (num_neibs > 0){
				sw /= num_neibs;
			}
			int nSurfTile3 = tnSurface3.getIndex(txy3[0],txy3[1]);
			img_data[0][nSurfTile3] = sdw;
			img_data[1][nSurfTile3] = Math.sqrt(sdw);
			img_data[2][nSurfTile3] = sdw;
			img_data[3][nSurfTile3] = sw;
		}
		ShowDoubleFloatArrays.showArrays(img_data,  3 * surfTilesX, 3 * surfTilesY, true, prefix+"tone_diffs3", titles);
	}

	public void showToneDiffWeights1(String prefix){
		double [][] img_data = new double[2][surfTilesX * surfTilesY];
		for (int nSurfTile = 0; nSurfTile <  tone_diff_weight.length;  nSurfTile++){
			double sdw = 0.0, sw = 0.0;
			int num_neibs = 0;
			for (int dir = 0; dir < 8; dir++){
				sw  += tone_diff_weight[nSurfTile][dir][1];
				sdw += tone_diff_weight[nSurfTile][dir][0] * tone_diff_weight[nSurfTile][dir][1];
				if (tone_diff_weight[nSurfTile][dir][1] > 0.0) num_neibs ++;
			}
			if (sw > 0.0){
				sdw /= sw;
			}
			if (num_neibs > 0){
				sw /= num_neibs;
			}
			img_data[0][nSurfTile] = sdw;
			img_data[1][nSurfTile] = sw;
		}
		String [] titles = {"diffs","weights"};
		ShowDoubleFloatArrays.showArrays(img_data,  surfTilesX, surfTilesY, true, prefix+"tone_diffs1", titles);
	}



	private double [] getDiffWeight(
			double [] tone1,
			double [] tone2)
	{
		double [] tone_weight = {0.0, Math.max(tone1[INDEX_A] * tone2[INDEX_A], 0.0)};
		if (tone_weight[1] > 0.0){
			double [][] scaled = {
					{(Math.max(tone1[INDEX_R], 0.0) + fatZero) * kR, (Math.max(tone1[INDEX_B], 0.0) + fatZero) * kB, (Math.max(tone1[INDEX_G], 0.0) + fatZero) },
					{(Math.max(tone2[INDEX_R], 0.0) + fatZero) * kR, (Math.max(tone2[INDEX_B], 0.0) + fatZero) * kB, (Math.max(tone2[INDEX_G], 0.0) + fatZero)}};
			double rdiff2 = 0.0;
			for (int i = 0; i < 3; i++){
				double rd = (scaled[0][i] - scaled[1][i])/(scaled[0][i] + scaled[1][i]);
				rdiff2 += rd*rd;
			}
			tone_weight[0] = rdiff2/3.0;
		}
		return tone_weight;
	}

	public void intValidSurf()
	{
		TileSurface.TileData [][] tileData = ts.getTileData();
		int nSurfTiles = surfTilesX * surfTilesY;
		valid_surf = new boolean [valid_ml.length][][]; // surfTilesX * surfTilesY][];
		for (int ml = 0; ml < valid_ml.length; ml++) if (valid_ml[ml]) {
			valid_surf[ml] = new boolean [nSurfTiles][];
		}
		for (int nSurfTile = 0; nSurfTile < nSurfTiles; nSurfTile++) if (tileData[nSurfTile] != null){
			boolean [] surf_en = new boolean [tileData[nSurfTile].length];
			for (int i = 0; i < surf_en.length; i++) {
				surf_en[i] = true;
			}
			for (int ml = 0; ml < valid_ml.length; ml++) if (valid_surf[ml] != null){
				valid_surf[ml][nSurfTile] = surf_en.clone();
			}
		}
	}

	/**
	 * limit layer options for assigned cells to a single (selected) option, disable all negative ones
	 * @param tileLayers
	 */

	public void restrictSingle (
			final int [][] tileLayers)
	{
		for (int ml = 0; ml < valid_ml.length; ml++) if ((valid_surf[ml] != null) && (tileLayers[ml] != null)){
			for (int nTile = 0; nTile < tileLayers[ml].length; nTile++) if (tileLayers[ml][nTile] != 0){
				int nSurfTile = ts.getSurfaceTileIndex(nTile);
				for (int ns = 0; ns < valid_surf[ml][nSurfTile].length; ns++){
					valid_surf[ml][nSurfTile][ns] = false;
				}
				if (tileLayers[ml][nTile] > 0) valid_surf[ml][nSurfTile][tileLayers[ml][nTile] - 1] = true;
			}
		}
	}
	/**
	 * Limit layer options for assigned cells to a multiple options per tile cell
	 * @param options
	 */

	public void restrictMulti (
			final int [][][] options)
	{
		for (int ml = 0; ml < valid_ml.length; ml++) if ((valid_surf[ml] != null) && (options[ml] != null)){
			for (int nTile = 0; nTile < options[ml].length; nTile++) if ((options[ml][nTile] != null) && (options[ml][nTile].length > 0)){
				int nSurfTile = ts.getSurfaceTileIndex(nTile);
				for (int ns = 0; ns < valid_surf[ml][nSurfTile].length; ns++){
					valid_surf[ml][nSurfTile][ns] = false;
				}
				for (int i = 0; i < options[ml][nTile].length; i++) {
					valid_surf[ml][nSurfTile][options[ml][nTile][i] - 1] = true;
				}
			}
		}
	}
	/**
	 * Get signed normalized difference from the measured disparity to the selected surface,
	 * limit it to +/-dispOutlier to reduce cost of the outliers
	 * @param ml measured layer (combo, quad, hor, vert)
	 * @param nSurfTile tile index in surface array (includes full supertiles)
	 * @param ns number of the surface
	 * @return signed disparity error, scaled down for large disparities
	 */
	public double dispDiff(
			int ml,
			int nSurfTile,
			int ns)
	{
		if ((dispStrength[ml] == null) || (dispStrength[ml][nSurfTile] == null)) {
				String msg="dispDiff(): no data for ml="+ml+", nSurfTile="+nSurfTile+", ns="+ns;
				System.out.println(msg);
				return Double.NaN;
//				IJ.showMessage(msg);
//				throw new IllegalArgumentException (msg);
		}
		double dm = dispStrength[ml][nSurfTile][0];
		double ds = ts.getTileData()[nSurfTile][ns].getDisparity();
		double diff = dm - ds;
		double d_av = 0.5 * (dm + ds);
		if (d_av > dispNorm){
			diff *= dispNorm/d_av;
		}
		if      (diff >  dispOutlier) diff =  dispOutlier;
		else if (diff < -dispOutlier) diff = -dispOutlier;
		return diff;

	}

	/**
	 * Get absolute value of the normalized disparity difference from the measured value to the nearest surface
	 * @param ml measured layer (combo, quad, hor, vert)
	 * @param nSurfTile tile index in surface array (includes full supertiles)
	 * @return absolute value of the disparity error (scaled down for large disparities) to the nearest surface
	 */
	public double dispDiffBest(
			int ml,
			int nSurfTile)
	{
		double best = Double.NaN;
		if ((ts.getTileData() == null) || (ts.getTileData()[nSurfTile] == null)){
//			System.out.println("dispDiffBest("+ml+","+nSurfTile+")");
			return best;
		}

		for (int ns = 0; ns < ts.getTileData()[nSurfTile].length; ns++){
			double diff = dispDiff(ml, nSurfTile, ns);
			double adiff = Math.abs(diff);
			if (Double.isNaN(best) || (adiff < best)) {
				best = adiff;
			}
		}
		return best;
	}

	/**
	 * Get absolute value of the normalized disparity difference from the measured value to the farthest surface
	 * Used to estimate cost of un-assigned tiles
	 * @param ml measured layer (combo, quad, hor, vert)
	 * @param nSurfTile tile index in surface array (includes full supertiles)
	 * @return absolute value of the disparity error (scaled down for large disparities) to the farthest surface
	 */
	public double dispDiffWorst(
			int ml,
			int nSurfTile)
	{
		double worst = Double.NaN;
		if ((ts.getTileData() == null) || (ts.getTileData()[nSurfTile] == null)){
//			System.out.println("dispDiffWorst("+ml+","+nSurfTile+")");
			return worst;
		}
		for (int ns = 0; ns < ts.getTileData()[nSurfTile].length; ns++){
			double diff = dispDiff(ml, nSurfTile, ns);
			double adiff = Math.abs(diff);
			if (Double.isNaN(worst) || (adiff < worst)) {
				worst = adiff;
			}
		}
		return worst;
	}

	public TACosts [] getTileCosts(
			boolean   all_tiles,
			int [][]  tileLayers)
	{
		TACosts [] ta_costs = new TACosts [surfTilesX * surfTilesY];
		for (int nSurfTile = 0; nSurfTile < ta_costs.length; nSurfTile++){
			boolean has_tile = false;
			boolean has_block = false;
			for (int ml = 0; ml <  tileLayers.length; ml++) if(tileLayers[ml] != null) {
				if (tileLayers[ml][nSurfTile] < 0) {
					has_block = true;
					break;
				}
				if (tileLayers[ml][nSurfTile] > 0) {
					has_tile = true;
				}
			}
			if (!has_block && (all_tiles || has_tile)){
				ta_costs[nSurfTile] = getTileCosts(
						nSurfTile, // int                    nSurfTile,
						tileLayers, // int [][]               tileLayers,
						null); // HashMap<Point,Integer> replacements)
			}
		}		//for (int nSurfTile)

		return ta_costs;
	}

	public double [] calcTileCosts(
			int [][]  tileLayers)
	{
		double [] costs = new double [surfTilesX * surfTilesY];
		for (int nSurfTile = 0; nSurfTile < costs.length; nSurfTile++){
			boolean has_tile = false;
			boolean has_block = false;
			for (int ml = 0; ml <  tileLayers.length; ml++) if(tileLayers[ml] != null) {
				if (tileLayers[ml][nSurfTile] < 0) {
					has_block = true;
					break;
				}
				if (tileLayers[ml][nSurfTile] > 0) {
					has_tile = true;
				}
			}
			if (!has_block && has_tile){
				TACosts ta_cost = getTileCosts(
						nSurfTile, // int                    nSurfTile,
						tileLayers, // int [][]               tileLayers,
						null); // HashMap<Point,Integer> replacements)
				costs[nSurfTile] = cost_coeff.dotProd(ta_cost);

			}
		}
		return costs;
	}

	public void showTileCost(
			String prefix,
			int [][]  tileLayers)
	{
		double [] composite_costs = calcTileCosts(tileLayers);
		ShowDoubleFloatArrays.showArrays(composite_costs, surfTilesX, surfTilesY, prefix+"composite_costs");
	}

	public void showTileCosts(
			String prefix,
			int [][]  tileLayers)
	{
		String [] titles = (new TACosts()).getTitles();
		int num_stiles = surfTilesX*surfTilesY;
		double [][] cost_components = new double[titles.length][surfTilesX * surfTilesY];
		double [] NaNs = new double[titles.length];
		for (int i = 0; i < NaNs.length; i++) NaNs[i] = Double.NaN;
		for (int nSurfTile = 0; nSurfTile < num_stiles; nSurfTile++){
			boolean has_tile = false;
			boolean has_block = false;
			for (int ml = 0; ml <  tileLayers.length; ml++) if(tileLayers[ml] != null) {
				if (tileLayers[ml][nSurfTile] < 0) {
					has_block = true;
					break;
				}
				if (tileLayers[ml][nSurfTile] > 0) {
					has_tile = true;
				}
			}
			double [] costs = NaNs;
			if (!has_block) { // && has_tile){
				TACosts ta_cost = getTileCosts(
						nSurfTile, // int                    nSurfTile,
						tileLayers, // int [][]               tileLayers,
						null); // HashMap<Point,Integer> replacements)
				ta_cost.mul(cost_coeff);
				costs = ta_cost.toArray();
			}
			for (int i =0; i < cost_components.length; i++){
				cost_components[i][nSurfTile] = costs[i];
			}
		}
		ShowDoubleFloatArrays.showArrays(cost_components, surfTilesX, surfTilesY, true, prefix+"cost_components", titles);
	}


	public TACosts [] statTileCosts(
			int [][]  tileLayers)
	{
		TACosts [] ta_stats = {new TACosts(), new TACosts()}; // average, max
		int num_tiles = 0;
		int numSurfTiles =surfTilesX * surfTilesY;
		for (int nSurfTile = 0; nSurfTile < numSurfTiles; nSurfTile++){
			if (nSurfTile == 47459){
				System.out.println("statTileCosts() nSurfTile="+nSurfTile);
			}

			boolean has_tile = false;
			boolean has_block = false;
			for (int ml = 0; ml <  tileLayers.length; ml++) if(tileLayers[ml] != null) {
				if (tileLayers[ml][nSurfTile] < 0) {
					has_block = true;
					break;
				}
				if (tileLayers[ml][nSurfTile] > 0) {
					has_tile = true;
				}
			}
			if (!has_block) { // && has_tile){
				TACosts ta_cost = getTileCosts(
						nSurfTile, // int                    nSurfTile,
						tileLayers, // int [][]               tileLayers,
						null); // HashMap<Point,Integer> replacements)
				ta_stats[0].max(ta_cost);
				ta_stats[1].add(ta_cost);
				num_tiles++;

			}
		}
		if (num_tiles > 0){
			ta_stats[1].scale(1.0/num_tiles);
		}
		return ta_stats;
	}

	/**
	 * Get costs (not scaled) for the particular tile, add for all defined layers
	 * @param nSurfTile tile index in surface array (includes full supertiles)
	 * @param tileLayers current assignment
	 * @param replacements optional assignment modification (or null)
	 * @return TACosts instance with assigned values
	 */
	public TACosts getTileCosts(
			int                    nSurfTile,
			int [][]               tileLayers,
			HashMap<Point,Integer> replacements)
	{
		TACosts costs = new TACosts();
		int debugLevel = 0;
		if (nSurfTile == -51360) { //  44831) { // 47459){
			System.out.println("getTileCosts() nSurfTile="+nSurfTile);
			debugLevel = 1;
		}
		int [][] around = new int [valid_ml.length][]; // contains layer + 1
		for (int ml = 0; ml <  valid_ml.length; ml++) if(valid_ml[ml]) {
			around[ml] = new int [9];
		}
		int [] nSurfTiles = new int[9];

		for (int dir = 0; dir < 9; dir++){
			nSurfTiles[dir] = tnSurface.getNeibIndex(nSurfTile, dir);
			if (nSurfTiles[dir] < 0) {
				for (int ml = 0; ml < around.length; ml++) if (around[ml] != null){
					around[ml][dir] = TileSurface.PROHOBITED;
				}
			} else {
				for (int ml = 0; ml < around.length; ml++) if (around[ml] != null){
					around[ml][dir] = tileLayers[ml][nSurfTiles[dir]];
					if (replacements != null){
						Integer isurf = replacements.get(new Point(ml, nSurfTiles[dir]));
						if (isurf != null) {
							around[ml][dir] = isurf; // substitute
						}
					}
				}
			}
		}
		//
//		public double empty;     // Cost of a tile that is not assigned
		for (int ml = 0; ml < around.length; ml++) if (around[ml] != null){
			if (around[ml][8] == 0) costs.empty +=1.0; // Cost of a tile that is not assigned
		}

//		public double nolink;    // Cost of a tile not having any neighbor in particular direction
		int num_weak_fgnd = 0; // to reduce multiple for teh same tile
		int num_color_sep = 0; // to reduce multiple for teh same tile

		for (int ml = 0; ml < around.length; ml++) if ((around[ml] != null) && (around[ml][8] > 0)){
			if ((around[ml][8] - 1) >= ts.getTileData()[nSurfTile].length){
				System.out.println("getTileCosts() BUG: nSurfTile="+nSurfTile);
			}
			int [] neibs = ts.getTileData()[nSurfTile][around[ml][8] - 1].getNeighbors();
			for (int dir = 0; dir < 8; dir++) if (neibs[dir] >= 0){ // do not count non-existing connections
				boolean link_this =  false;
				boolean link_other = false;
				int ml1 = ml;
				int neibp1 = neibs[dir] + 1;
				if (around[ml][dir] > 0) {
					if (around[ml][dir] == neibp1) {
						link_this = true;

					} else {
						link_other = true;
					}
				} else if (around[ml][dir] == 0){
					// see if some other ml contains this surface (ignore other surfaces
					for (ml1 = 0; ml1 < around.length; ml1++) if ((ml1 != ml) && (around[ml1] != null)) {
						if (around[ml1][dir] == neibp1){
							link_this = true;
							break; // ml1 will keep other ml
						}
					}
				} else { // negative/prohibited
					continue;
				}
				if (!link_this && !link_other){
					costs.nolink +=1.0;
				} else if (link_other) {
					costs.swtch += 1.0; // cost for any switch
					// check if it is fb->bg transition and fg is weak
					double d_this =  ts.getTileData()[nSurfTiles[dir]][neibs[dir]].getDisparity();
					double d_other = ts.getTileData()[nSurfTiles[dir]][around[ml][dir] -1].getDisparity();
					double disp_diff = d_this - d_other;
					double d_av = 0.5 * (d_this + d_other);
					if (d_av > dispNorm){
						disp_diff *= dispNorm/d_av;
					}
					// disp_diff here is signed, positive if center is FG, other is BG
					if (disp_diff > minFgBg) {
						double strength = dispStrength[ml][nSurfTile][1];
						strength = Math.max(strength, minFgEdge);
						costs.weak_fgnd += minFgEdge / strength;
						num_weak_fgnd++;
					} else if (mutual_weak_fgnd && (disp_diff < -minFgBg)) {
						double [] dsmeas_other = dispStrength[ml][nSurfTiles[dir]];
						double strength = (dsmeas_other != null)? dsmeas_other[1]: 0.0; // measured strength on the other end or 0.0 if nothing there
						strength = Math.max(strength, minFgEdge);
						costs.weak_fgnd += minFgEdge / strength;
						num_weak_fgnd++;
					}
					if (Math.abs(disp_diff) > minColSep){
						double col_diff = tone_diff_weight[nSurfTile][dir][0];
						col_diff= Math.max(col_diff,minColDiff);
						costs.color += tone_diff_weight[nSurfTile][dir][1] * minColSep/col_diff;
						num_color_sep++;
					}

				} else { // v, both can not coexist
					// Anything to cost here?
				}
			}
//			if (around[ml][8] == 0) costs.empty +=1.0; // each existing measurement layer that is not assigned

		} //for (int ml = 0; ml < around.length; ml++) if ((around[ml] != null) && (around[ml][8] > 0))

		// using /4.0 to maintain same value (==1.0) for half neighbors (straight edge)
		if ((num_weak_fgnd > 0) && (shrinkWeakFgnd > 0.0)){
			costs.weak_fgnd /= Math.pow(num_weak_fgnd/4.0, shrinkWeakFgnd);
		}

		if ((num_color_sep > 0) && (shrinkColor > 0.0)){
			costs.color /= Math.pow(num_color_sep/4.0, shrinkColor);
		}

		double disp_diff_lpf = 0.0, disp_diff_weight = 0.0;  // calculate LPF of the disparity signed error over all ML and 9 cells
		double disp_diff2 = 0.0, disp_weight = 0.0;  // calculate disparity error in the center, weighted
		double disp_diff_over_best = 0.0; // , disp_weight_best = 0.0;  // calculate disparity error in the center, weighted
		for (int ml = 0; ml < around.length; ml++) if (around[ml] != null){
			double diff=0.0, weight=0.0;
			for (int dir = 0; dir < 9; dir++){
				if (around[ml][dir] > 0){ // assigned
					diff = dispDiff(
							ml, // int ml,
							nSurfTiles[dir], // int nSurfTile,
							around[ml][dir] - 1); // int ns)
				} else if ((dir == 8) && (around[ml][dir] == 0) ) { //no assigned (not prohibited) - center only, because worst does not have sign
					diff = dispDiffWorst( // worst for all surfaces
							ml, // int ml,
							nSurfTiles[dir]); // int nSurfTile,
				}
				if (dispStrength[ml][nSurfTile] != null) { // not a bug
					if (strengthDiff9Pwr > 0.0) {
						weight = dispStrength[ml][nSurfTile][1];
						if (strengthDiff9Pwr != 1.0) {
							weight = Math.pow(weight, strengthDiff9Pwr);
						}
					} else {
						weight = 1.0;
					}
					weight *= NEIB_WEIGHTS[dir];
					disp_diff_lpf += diff * weight;
					disp_diff_weight += weight;
					if (debugLevel > 0){
						System.out.println("getTileCosts() nSurfTile = "+nSurfTile+"->"+dir+" disp_diff_lpf="+disp_diff_lpf+
								" disp_diff_weight="+disp_diff_weight+" weight="+weight+ " diff="+diff+" around["+ml+"]["+dir+"]="+around[ml][dir]);
					}
				}
			}
			// now diff is for the center, weight needs to be re-calculated
			if (strengthDiffPwr > 0.0) {
				if ((dispStrength[ml] == null) || (dispStrength[ml][nSurfTile] == null)){
					System.out.println("getTileCosts() nSurfTile = "+nSurfTile+" ml = "+ml+" is it really a BUG - null pointer here?");
					weight = 1.0;
				} else {
					weight = dispStrength[ml][nSurfTile][1]; // null pointer
					if (strengthDiffPwr != 1.0) {
						weight = Math.pow(weight, strengthDiffPwr);
					}
				}
			} else {
				weight = 1.0;
			}
			disp_diff2 += diff * diff * weight;
			disp_weight += weight;

			// and once more for the disparity over best
			if (strengthBestPwr > 0.0) {
				weight = dispStrength[ml][nSurfTile][1];
				if (strengthBestPwr != 1.0) {
					weight = Math.pow(weight, strengthBestPwr);
				}
			} else {
				weight = 1.0;
			}
			disp_diff_over_best += (Math.abs(diff) -  dispDiffBest( ml, nSurfTiles[8])) * weight; // this one - not squared

		}
		if (disp_diff_weight > 0.0) {
			disp_diff_lpf /= disp_diff_weight;
			costs.diff9 += disp_diff_lpf * disp_diff_lpf;
			if (debugLevel > 0){
				System.out.println("getTileCosts() nSurfTile = "+nSurfTile+" disp_diff_lpf="+disp_diff_lpf+
						" disp_diff_weight="+disp_diff_weight+" costs.diff9="+costs.diff9);
			}
		}
		if (disp_weight > 0.0) {
			costs.diff +=      disp_diff2 / disp_weight;
			costs.diff_best += disp_diff_over_best / disp_weight;
		}

		int [] txy = tnSurface.getXY(nSurfTile);
		int nSTile = (txy[0] / ts.getSuperTileSize()) + (txy[1] / ts.getSuperTileSize()) *  ts.getSTilesX();
		for (int ml = 0; ml < around.length; ml++) if ((around[ml] != null) && (around[ml][8] > 0)){
			int parentSTile =ts.getTileData()[nSurfTile][around[ml][8] -1].getParentNsTile();
			if (parentSTile != nSTile){
				costs.flaps += 1.0;
			}
		}
		return costs;
	}

	public int [][] imgToSurf (
			int [][] tileLayersImg)
	{
		int [][] tileLayersSurf = new int[tileLayersImg.length][];
		for (int ml = 0; ml < tileLayersImg.length; ml++) if (tileLayersImg[ml] != null){
			tileLayersSurf[ml] = new int [tnSurface.getLength()];
			for (int i = 0; i < tileLayersSurf[ml].length; i++){
				tileLayersSurf[ml][i] = -1;
			}
			for (int nTile = 0; nTile < tileLayersImg[ml].length; nTile++){
				int nSurfTile = tnSurface.getIndex(tnImage.getXY(nTile));
				if (nSurfTile >= 0) {
					tileLayersSurf[ml][nSurfTile] = tileLayersImg[ml][nTile];
				}
			}
		}
		return tileLayersSurf;
	}

	public int [][] surfToImg (
			int [][] tileLayersSurf)
	{
		int [][] tileLayersImg = new int[tileLayersSurf.length][];
		for (int ml = 0; ml < tileLayersSurf.length; ml++) if (tileLayersSurf[ml] != null){
			tileLayersImg[ml] = new int [tnImage.getLength()];
			for (int i = 0; i < tileLayersImg[ml].length; i++){
				tileLayersImg[ml][i] = -1;
			}
			for (int nSurfTile = 0; nSurfTile < tileLayersSurf[ml].length; nSurfTile++){
				int nImgTile = tnImage.getIndex(tnSurface.getXY(nSurfTile));
				if (nImgTile >= 0) {
					tileLayersImg[ml][nImgTile] = tileLayersSurf[ml][nSurfTile];
				}
			}
		}
		return tileLayersImg;
	}

	public void optimizeAssignment9(
	        final boolean   noEdge,
			final int [][]  tileLayers,
            final int       debugLevel,
			final int       dbg_X,
			final int       dbg_Y)
	{
		mutual_weak_fgnd = true;
		final int step = 3;
		final int tries = 1000;
//		final int dbg_tile = dbg_X + dbg_Y * surfTilesX;
		final int dbg_tile = 47779; // 27083; // 44493;
		final int num_tiles = surfTilesX * surfTilesY;
		final int [][] tile_indices = new int [step*step][];
		for (int sty = 0; sty < step; sty ++){
			int num_y = (surfTilesY + step -1 - sty) / step;
			for (int stx = 0; stx < step; stx ++){
				int num_x = (surfTilesX + step -1 - stx) / step;
				int indx1 = sty * step + stx;
				int l = num_y * num_x;
				tile_indices[indx1] = new int [l];
				int indx2 = 0;
				for (int y = 0; y < num_y; y++){
					for (int x = 0; x < num_x; x++){
						tile_indices[indx1][indx2++] = (sty + step * y) * surfTilesX + (stx + step * x);
					}
				}
			}
		}
		final Thread[] threads = ImageDtt.newThreadArray((debugLevel > 1)? 1 : ts.getThreadsMax());
		final int numThreads =   threads.length;
		final int [] improved =  new int [numThreads];
		final AtomicInteger ai_numThread = new AtomicInteger(0);
		final AtomicInteger ai = new AtomicInteger(0);
		final AtomicInteger ai_series = new AtomicInteger(0);
		final AtomicBoolean [] dirty = new AtomicBoolean[surfTilesX * surfTilesY];
		for (int nSurfTile = 0; nSurfTile < dirty.length; nSurfTile++){
			boolean valid_tile = false; // bad may be only some ml?
			for (int ml = 0; ml < tileLayers.length; ml++) if ((tileLayers[ml] != null) && (tileLayers[ml][nSurfTile] >= 0)){
				valid_tile = true;
				break;
			}
			if (valid_tile) dirty[nSurfTile] = new AtomicBoolean(true);
		}

		final boolean [][] bad_surface = new boolean[num_tiles][];
		final TileSurface.TileData [][] tileData = ts.getTileData();
		if (noEdge) {
			for (int ithread = 0; ithread < threads.length; ithread++) {
				threads[ithread] = new Thread() {
					@Override
					public void run() {
						for (int nSurfTile = ai.getAndIncrement(); nSurfTile < num_tiles; nSurfTile = ai.getAndIncrement()) {
							if (tileData[nSurfTile] != null){
								bad_surface[nSurfTile] = new boolean [tileData[nSurfTile].length];
								for (int ns = 0; ns < tileData[nSurfTile].length; ns++) {
									int []neibs = tileData[nSurfTile][ns].getNeighbors();
									for (int i = 0; i < neibs.length; i++){
										if (neibs[i] < 0){
											bad_surface[nSurfTile][ns] = true;
										}
									}
								}
							}
						}
					}
				};
			}
			ImageDtt.startAndJoin(threads);
		}

		for (int nTry = 0 ; nTry < tries; nTry++) {
			for (int i = 0; i < improved.length; i++) improved[i] = 0;
			int this_improved = 0;
			for (int nSeries = 0; nSeries < tile_indices.length; nSeries++) {
				final int fnSeries = nSeries;
				ai.set(0);
				ai_numThread.set(0);
				for (int ithread = 0; ithread < threads.length; ithread++) {
					threads[ithread] = new Thread() {
						@Override
						public void run() {
							int numThread = ai_numThread.getAndIncrement(); // unique number of thread to write to rslt_diffs[numThread]
							for (int iTile = ai.getAndIncrement(); iTile < tile_indices[fnSeries].length; iTile = ai.getAndIncrement()) {
								int nSurfTile = tile_indices[fnSeries][iTile];
	                            int dl = ((debugLevel > 1) && (nSurfTile == dbg_tile)) ? 3: debugLevel;
								if (dl > 2){
									System.out.println("optimizeAssignment9(), nSurfTile = "+nSurfTile+" dl = "+dl);
								}
//								if (dirty[nSurfTile].get()) {
								if ((dirty[nSurfTile] != null) && dirty[nSurfTile].getAndSet(false)) {
									int num_surf = tileData[nSurfTile].length;
									int lowest_surf = 0;
									for (; lowest_surf < num_surf; lowest_surf ++){
										if ((bad_surface[nSurfTile] != null) && !bad_surface[nSurfTile][lowest_surf]) {
											break;
										}
									}
									if (lowest_surf >=  num_surf) {
										continue; // no valid surfaces at this location
									}
									double best_cost = cost_coeff.dotProd(getTileCosts(
											nSurfTile,  // int                    nSurfTile,
											tileLayers, // int [][]               tileLayers,
											null));     // HashMap<Point,Integer> replacements);
//									int [] initial_indices = tileLayers[nSurfTile].clone();
									int [] initial_indices = new int [valid_ml.length]; // 1-based
									for (int ml = 0; ml < valid_ml.length; ml ++ ){
										if (tileLayers[ml] != null) {
											initial_indices[ml] = tileLayers[ml][nSurfTile];  // 1-based
										}
									}
									int [] best_surf = null;
									int [] surfaces = null; // new int [valid_ml.length];
									int max_reset = 0; //  = valid_ml.length; // maximal ml to reset
									while (true) {
										if (surfaces == null) {
											surfaces = new int [valid_ml.length];
											for (int ml = 0; ml < valid_ml.length; ml ++){
												if  (!valid_ml[ml] || (tileLayers[ml][nSurfTile] < 0)) {
													surfaces[ml] = -1;
												}
											}
											max_reset = valid_ml.length;
										} else { // find ml to increase surface
											for (max_reset = 0; max_reset <  surfaces.length; max_reset++) if (surfaces[max_reset] >= 0){
												for (surfaces[max_reset]++; surfaces[max_reset] < num_surf; surfaces[max_reset]++) {
													if ((bad_surface[nSurfTile] == null) || !bad_surface[nSurfTile][surfaces[max_reset]]) {
														break;
													}

												}
												if (surfaces[max_reset] < num_surf) break;
											}
											if (max_reset >=  surfaces.length){
												break; // while (true) {
											}
										}
										// reset all surfaces[] with indices < max_reset to lowest
										for (int ml = 0; ml < max_reset; ml++ ) if (valid_ml[ml] && (surfaces[ml] >= 0)){
											surfaces[ml] = lowest_surf;
										}
										// now surfaces[] contain next combination of surfaces to try

										// tileLayers[nSurfTile] = surfaces;
										for (int ml = 0; ml < valid_ml.length; ml ++ ){
											if (tileLayers[ml] != null) {
												tileLayers[ml][nSurfTile] = surfaces[ml] + 1; // 1-based from 0-based
											}
										}

										double cost = cost_coeff.dotProd(getTileCosts( //
												nSurfTile,  // int                    nSurfTile,
												tileLayers, // int [][]               tileLayers,
												null));     // HashMap<Point,Integer> replacements);
										if (cost < best_cost) {
											best_cost = cost;
											best_surf = surfaces.clone();
										}

									} // while (true)
									if (best_surf != null){ // update
//										tileLayers[nSurfTile] = best_surf;
										for (int ml = 0; ml < valid_ml.length; ml ++ ){
											if (tileLayers[ml] != null) {
												tileLayers[ml][nSurfTile] = best_surf[ml] + 1;  // 1-based from 0-based
											}
										}


										for (int dir = 0; dir <8; dir++) {
											int nSurfTile1 = tnSurface.getNeibIndex(nSurfTile, dir);
											if ((nSurfTile1 >= 0) && (dirty[nSurfTile1] != null)){
												dirty[nSurfTile1].set(true);
											}
										}
										improved[numThread]++;
									} else { // restore initial data
//										tileLayers[nSurfTile] = initial_indices;
										for (int ml = 0; ml < valid_ml.length; ml ++ ){
											if (tileLayers[ml] != null) {
												tileLayers[ml][nSurfTile] = initial_indices[ml];  // 1-based from 1-based
											}
										}
									}
								}
							}
						}
					};
				}
				ImageDtt.startAndJoin(threads);
				if (debugLevel > -1){
					int num_better = 0;
					for (int i = 0; i < improved.length; i++){
						num_better += improved[i];
					}
					System.out.println("optimizeAssignment9(): pass = "+nTry+ ":"+nSeries+" improved:" + (num_better - this_improved));
					this_improved = num_better;
				}
			} // for (int nSeries = 0; nSeries < tile_indices.length; nSeries++) {
			// should be checked only after all series (now 9 passes) are finished - if anything was added - continue
			int num_improved = 0;
			for (int i = 0; i < improved.length; i++){
				num_improved += improved[i];
			}
			if (debugLevel > -1){
				System.out.println("optimizeAssignment9(): pass = "+nTry+ " improved:"+num_improved);
			}

			if (num_improved == 0) break;
		}
	}


	public void optimizeAssignment25(
	        final boolean   noEdge,
			final int [][]  tileLayers,
            final int       debugLevel,
			final int       dbg_X,
			final int       dbg_Y)
	{
		final int step = 5;
		mutual_weak_fgnd = false;
		final int tries = 1000;
//		final int dbg_tile = dbg_X + dbg_Y * surfTilesX;
		final int dbg_tile = 51360; // 51; // 44831; // 27083; // 44493;
		final int num_tiles = surfTilesX * surfTilesY;
		final int [][] tile_indices = new int [step*step][];
		for (int sty = 0; sty < step; sty ++){
			int num_y = (surfTilesY + step -1 - sty) / step;
			for (int stx = 0; stx < step; stx ++){
				int num_x = (surfTilesX + step -1 - stx) / step;
				int indx1 = sty * step + stx;
				int l = num_y * num_x;
				tile_indices[indx1] = new int [l];
				int indx2 = 0;
				for (int y = 0; y < num_y; y++){
					for (int x = 0; x < num_x; x++){
						tile_indices[indx1][indx2++] = (sty + step * y) * surfTilesX + (stx + step * x);
					}
				}
			}
		}
		final Thread[] threads = ImageDtt.newThreadArray((debugLevel > 1)? 1 : ts.getThreadsMax());
		final int numThreads =   threads.length;
		final int [] improved =  new int [numThreads];
		final AtomicInteger ai_numThread = new AtomicInteger(0);
		final AtomicInteger ai = new AtomicInteger(0);
		final AtomicInteger ai_series = new AtomicInteger(0);
		final AtomicBoolean [] dirty = new AtomicBoolean[surfTilesX * surfTilesY];
		for (int nSurfTile = 0; nSurfTile < dirty.length; nSurfTile++){
			boolean valid_tile = false; // bad may be only some ml?
			for (int ml = 0; ml < tileLayers.length; ml++) if ((tileLayers[ml] != null) && (tileLayers[ml][nSurfTile] >= 0)){
				valid_tile = true;
				break;
			}
			if (valid_tile) dirty[nSurfTile] = new AtomicBoolean(true);
		}

		final boolean [][] bad_surface = new boolean[num_tiles][];
		final TileSurface.TileData [][] tileData = ts.getTileData();
		if (noEdge) {
			for (int ithread = 0; ithread < threads.length; ithread++) {
				threads[ithread] = new Thread() {
					@Override
					public void run() {
						for (int nSurfTile = ai.getAndIncrement(); nSurfTile < num_tiles; nSurfTile = ai.getAndIncrement()) {
							if (tileData[nSurfTile] != null){
								bad_surface[nSurfTile] = new boolean [tileData[nSurfTile].length];
								for (int ns = 0; ns < tileData[nSurfTile].length; ns++) {
									int []neibs = tileData[nSurfTile][ns].getNeighbors();
									for (int i = 0; i < neibs.length; i++){
										if (neibs[i] < 0){
											bad_surface[nSurfTile][ns] = true;
										}
									}
								}
							}
						}
					}
				};
			}
			ImageDtt.startAndJoin(threads);
		}

		for (int nTry = 0 ; nTry < tries; nTry++) {
			for (int i = 0; i < improved.length; i++) improved[i] = 0;
			int this_improved = 0;
			for (int nSeries = 0; nSeries < tile_indices.length; nSeries++) {
				final int fnSeries = nSeries;
				ai.set(0);
				ai_numThread.set(0);
				for (int ithread = 0; ithread < threads.length; ithread++) {
					threads[ithread] = new Thread() {
						@Override
						public void run() {
							int numThread = ai_numThread.getAndIncrement(); // unique number of thread to write to rslt_diffs[numThread]
							for (int iTile = ai.getAndIncrement(); iTile < tile_indices[fnSeries].length; iTile = ai.getAndIncrement()) {
								int nSurfTile = tile_indices[fnSeries][iTile];
	                            int dl = ((debugLevel > 1) && (nSurfTile == dbg_tile)) ? 3: debugLevel;
								if (dl > 2){
									System.out.println("optimizeAssignment25(), nSurfTile = "+nSurfTile+" dl = "+dl);
								}
//								if (dirty[nSurfTile].get()) {
								if ((dirty[nSurfTile] != null) && dirty[nSurfTile].getAndSet(false)) {
									int num_surf = tileData[nSurfTile].length;
									int lowest_surf = 0;
									for (; lowest_surf < num_surf; lowest_surf ++){
										if ((bad_surface[nSurfTile] != null) && !bad_surface[nSurfTile][lowest_surf]) {
											break;
										}
									}
									if (lowest_surf >=  num_surf) {
										continue; // no valid surfaces at this location
									}
									double best_cost = 0.0;
									if (dl > 2) {
										System.out.println("optimizeAssignment25(), tileLayers[0]["+nSurfTile+"] = " + tileLayers[0][nSurfTile]);
//										int [] dbg_layers = new int [tileLayers.length];
//										for (int ml = 0; ml < tileLayers.length; ml++){
//											dbg_layers[ml] = (tileLayers[ml] == null) ? -0:  tileLayers[ml][nSurfTile];
//										}
//										System.out.println("optimizeAssignment25(), tileLayers["+nSurfTile+"] = " + dbg_layers);
									}
									for (int dir = 0; dir < 9; dir++)  {
										int nSurfTile1 = tnSurface.getNeibIndex(nSurfTile, dir);
										if (nSurfTile1 >= 0){
//											best_cost += cost_coeff.dotProd(getTileCosts(
//													nSurfTile1,  // int                    nSurfTile,
//													tileLayers, // int [][]               tileLayers,
//													null));     // HashMap<Point,Integer> replacements);
											TACosts ta_costs = getTileCosts(
													nSurfTile1,  // int                    nSurfTile,
													tileLayers, // int [][]               tileLayers,
													null);     // HashMap<Point,Integer> replacements);
											double ccost = cost_coeff.dotProd(ta_costs);
											best_cost += ccost;
											if (dl > 2){
												System.out.println("optimizeAssignment25(), nSurfTile = "+nSurfTile+" dir = "+dir +
														" nSurfTile1="+nSurfTile1+" ccost = "+ccost+" best_cost="+best_cost);
												System.out.println("ta_costs["+nSurfTile1+"]="+ta_costs.toString());
											}
										}
									}
//									int [] initial_indices = tileLayers[nSurfTile].clone();
									int [] initial_indices = new int [valid_ml.length]; // 1-based
									for (int ml = 0; ml < valid_ml.length; ml ++ ){
										if (tileLayers[ml] != null) {
											initial_indices[ml] = tileLayers[ml][nSurfTile];  // 1-based
										}
									}
									int [] best_surf = null;
									int [] surfaces = null; // new int [valid_ml.length];
									int max_reset = 0; //  = valid_ml.length; // maximal ml to reset
									while (true) {
										if (surfaces == null) {
											surfaces = new int [valid_ml.length];
											for (int ml = 0; ml < valid_ml.length; ml ++){
												if  (!valid_ml[ml] || (tileLayers[ml][nSurfTile] < 0)) {
													surfaces[ml] = -1;
												}
											}
											max_reset = valid_ml.length;
										} else { // find ml to increase surface
											for (max_reset = 0; max_reset <  surfaces.length; max_reset++) if (surfaces[max_reset] >= 0){
												for (surfaces[max_reset]++; surfaces[max_reset] < num_surf; surfaces[max_reset]++) {
													if ((bad_surface[nSurfTile] == null) || !bad_surface[nSurfTile][surfaces[max_reset]]) {
														break;
													}

												}
												if (surfaces[max_reset] < num_surf) break;
											}
											if (max_reset >=  surfaces.length){
												break; // while (true) {
											}
										}
										// reset all surfaces[] with indices < max_reset to lowest
										for (int ml = 0; ml < max_reset; ml++ ) if (valid_ml[ml] && (surfaces[ml] >= 0)){
											surfaces[ml] = lowest_surf;
										}
										// now surfaces[] contain next combination of surfaces to try

										// tileLayers[nSurfTile] = surfaces;
										for (int ml = 0; ml < valid_ml.length; ml ++ ){
											if (tileLayers[ml] != null) {
												tileLayers[ml][nSurfTile] = surfaces[ml] + 1; // 1-based from 0-based
											}
										}
										double cost = 0.0;
										if (dl > 2) {
											System.out.println("optimizeAssignment25(), surfaces[0]= " + surfaces[0]);
										}
										for (int dir = 0; dir < 9; dir++)  {
											int nSurfTile1 = tnSurface.getNeibIndex(nSurfTile, dir);
											if (nSurfTile1 >= 0){
												TACosts ta_costs = getTileCosts(
														nSurfTile1,  // int                    nSurfTile,
														tileLayers, // int [][]               tileLayers,
														null);     // HashMap<Point,Integer> replacements);
												double ccost = cost_coeff.dotProd(ta_costs);
												cost+=ccost;
												if (dl > 2){
													System.out.println("optimizeAssignment25(), nSurfTile = "+nSurfTile+" dir = "+dir +
															" nSurfTile1="+nSurfTile1+" ccost = "+ccost+" cost="+cost);
													System.out.println("ta_costs["+nSurfTile1+"]="+ta_costs.toString());
												}
											}
										}
										if (cost < best_cost) {
											best_cost = cost;
											best_surf = surfaces.clone();
										}

									} // while (true)
									if (best_surf != null){ // update
//										tileLayers[nSurfTile] = best_surf;
										for (int ml = 0; ml < valid_ml.length; ml ++ ){
											if (tileLayers[ml] != null) {
												tileLayers[ml][nSurfTile] = best_surf[ml] + 1;  // 1-based from 0-based
											}
										}


										for (int dir = 0; dir <8; dir++) {
											int nSurfTile1 = tnSurface.getNeibIndex(nSurfTile, dir);
											if ((nSurfTile1 >= 0) && (dirty[nSurfTile1] != null)){
												dirty[nSurfTile1].set(true);
											}
										}
										improved[numThread]++;
									} else { // restore initial data
//										tileLayers[nSurfTile] = initial_indices;
										for (int ml = 0; ml < valid_ml.length; ml ++ ){
											if (tileLayers[ml] != null) {
												tileLayers[ml][nSurfTile] = initial_indices[ml];  // 1-based from 1-based
											}
										}
									}
								}
							}
						}
					};
				}
				ImageDtt.startAndJoin(threads);
				if (debugLevel > -1){
					int num_better = 0;
					for (int i = 0; i < improved.length; i++){
						num_better += improved[i];
					}
					System.out.println("optimizeAssignment25(): pass = "+nTry+ ":"+nSeries+" improved:" + (num_better - this_improved));
					this_improved = num_better;
				}
			} // for (int nSeries = 0; nSeries < tile_indices.length; nSeries++) {
			// should be checked only after all series (now 9 passes) are finished - if anything was added - continue
			int num_improved = 0;
			for (int i = 0; i < improved.length; i++){
				num_improved += improved[i];
			}
			if (debugLevel > -1){
				System.out.println("optimizeAssignment25(): pass = "+nTry+ " improved:"+num_improved);
			}

			if (num_improved == 0) break;
		}
	}


}
