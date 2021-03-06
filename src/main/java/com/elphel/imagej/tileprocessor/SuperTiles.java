package com.elphel.imagej.tileprocessor;
/**
 **
 ** SuperTiles - Process tiles resulted from quad image sets
 **
 ** Copyright (C) 2017 Elphel, Inc.
 **
 ** -----------------------------------------------------------------------------**
 **
 **  SuperTiles.java is free software: you can redistribute it and/or modify
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
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.concurrent.atomic.AtomicInteger;

import com.elphel.imagej.common.DoubleGaussianBlur;
import com.elphel.imagej.common.ShowDoubleFloatArrays;

public class SuperTiles{
	public TileProcessor tileProcessor;
	double      step_far;
	double      step_near;
	double      step_threshold_far;    // relative to min_disparity
	double      step_threshold_near;   // relative to min_disparity
	double      bin_far;               // bin value (before rounding) matching  (min_disparity + step_threshold_far)
	double      bin_near;              // bin value (before rounding) matching  (min_disparity + step_threshold_near)
	double      min_disparity;
	double      max_disparity;
	double      stBlurSigma;
	int         numBins;
	double []   bin_centers;
	double [][] disparityHistograms = null;
	double []   stStrength =          null; // per super-tile correlation strength
	double [][][] maxMinMax = null;
	double []   bgDisparity = null;
	double []   bgStrength = null;
	int         measSel = 1; // bitmask of the selected measurements for supertiles : +1 - combo, +2 - quad +4 - hor +8 - vert

	boolean     smplMode        = true;   // Use sample mode (false - regular tile mode)

	/*
	double      strength_floor;
	double      strength_pow;
	int         smplSide        = 2;      // Sample size (side of a square)
	int         smplNum         = 3;      // Number after removing worst
	double      smplRms         = 0.1;    // Maximal RMS of the remaining tiles in a sample
	boolean     smplWnd         = false;  // final boolean    smplWnd,  // use window functions for the samples
	double      max_abs_tilt  = 2.0; // Maximal absolute tilt in pixels/tile
	double      max_rel_tilt  = 0.2; // Maximal relative tilt in pixels/tile/disparity
	double      damp_tilt  =    0.001; // Damp tilt to handle insufficient  (co-linear)data
	double      min_tilt_disp = 4.0; // Disparity switch between filtering modes - near objects use tilts, far - use max disparity
	double      transition    = 1.0; // Mode transition range (between tilted and maximal disparity)
	int         far_mode  =     1;   // Far objects filtering mode (0 - off, 1 - power of disparity)
	double      far_power =     1.0; // Raise disparity to this power before averaging for far objects
*/
	MeasuredLayersFilterParameters mlfp = new MeasuredLayersFilterParameters();     // filter parameters
	MeasuredLayers measuredLayers = null;

	CLTPass3d cltPass3d;
	TilePlanes.PlaneData [][] planes =  null;
	TilePlanes.PlaneData [][] planes_mod = null;
	double [][][] fuse_coeff = null; // coefficients to fuse planes in neighbor supertiles

	int [][]    shell_map = null;   // per supertile, per disparity plane - shell index + 1 (0 - none)
	double [][] surfaces;  // per shell, per tile (linescan order) disparity value or NaN in missing supertiles

	TileSurface tileSurface = null;


	/**
	 * currently lowest plane for each includes all tiles, so do not use it. May change in the future
	 * @param np total number of planes in a supertile
	 * @return lowest plane to iterate through
	 */

	static int LOWEST_PLANE(int np){ //
		return (np >1) ? 1 : 0;
	}

	public SuperTiles(
			CLTPass3d cltPass3d,
			double                  step_near,
			double                  step_far,
			double                  step_threshold,
			double                  min_disparity,
			double                  max_disparity,

			double                  stBlurSigma,
			boolean                 smplMode, //        = true;   // Use sample mode (false - regular tile mode)
			MeasuredLayersFilterParameters  mlfp,

//			double                  strength_floor,
//			double                  strength_pow,
//			boolean                 smplMode, //        = true;   // Use sample mode (false - regular tile mode)
//			int                     smplSide, //        = 2;      // Sample size (side of a square)
//			int                     smplNum, //         = 3;      // Number after removing worst
//			double                  smplRms, //         = 0.1;    // Maximal RMS of the remaining tiles in a sample
//			boolean                 smplWnd,  // use window functions for the samples

//  			double     max_abs_tilt,  //  2.0;   // pix per tile
//			double     max_rel_tilt,  //  0.2;   // (pix / disparity) per tile
//			double     damp_tilt,     //  0.001; // Damp tilt to handle insufficient  (co-linear)data
//			double     min_tilt_disp, //  4.0;   // Disparity switch between filtering modes - near objects use tilts, far - use max disparity
//			double     transition,    //  1.0;   // Mode transition range (between tilted and maximal disparity)
//			int        far_mode,      //  1;     // Far objects filtering mode (0 - off, 1 - power of disparity)
//			double     far_power,     //  3.0;   // Raise disparity to this power before averaging for far objects
//			boolean    null_if_none,

			int                     measSel)
	{
		this.cltPass3d =           cltPass3d;
		this.tileProcessor =       cltPass3d.getTileProcessor();
		this.step_near =           step_near;
		this.step_far =            step_far;
		this.step_threshold_far = step_threshold;
		this.min_disparity =  min_disparity;
		this.max_disparity =  max_disparity;
		this.stBlurSigma =    stBlurSigma;
		this.smplMode        = smplMode;   // Use sample mode (false - regular tile mode)
		this.mlfp = mlfp.clone();

//		this.strength_floor = strength_floor;
//		this.strength_pow =   strength_pow;
//		this.smplSide        = smplSide;   // Sample size (side of a square)
//		this.smplNum         = smplNum;    // Number after removing worst
//		this.smplRms         = smplRms;    // Maximal RMS of the remaining tiles in a sample
//		this.max_abs_tilt =  max_abs_tilt;
//		this.max_rel_tilt =  max_rel_tilt;
//		this.damp_tilt =     damp_tilt;
//		this.min_tilt_disp = min_tilt_disp;
//		this.transition =    transition;
//		this.far_mode =      far_mode;
//		this.far_power =     far_power;
//		this.smplWnd         = smplWnd;    // Use window functions for the samples


		this.measSel =        measSel;
		this.step_threshold_near = this.step_threshold_far * step_near / this.step_far ;
		this.bin_far =             this.step_threshold_far / this.step_far;
		this.bin_near =            this.step_threshold_far / this.step_far * (Math.log(this.step_near/this.step_far) + 1);
		int bin_max = disparityToBin(max_disparity);
		numBins = bin_max + 1;
		bin_centers = new double [numBins];
		int bin = 0;
		for (bin = 0; bin < numBins; bin ++){
			bin_centers[bin] = binToDisparity(bin);
		}
//		initFuseCoeff(0.5, true); // true);
		initFuseCoeff(0.5, false); // true);

// Set up MeasuredLayers
		measuredLayers = new MeasuredLayers(
				4, // combo, quad, hor, vert
				tileProcessor.getTilesX(),
				tileProcessor.getTilesY(),
				tileProcessor.getSuperTileSize());
// MeasuredLayer.setLayer(int, double[], double[], boolean[])
		measuredLayers.setLayer( // combo disparity
				0,                                // int       num_layer,
				cltPass3d.getDisparity(0),        // double [] disparity,
				cltPass3d.getStrength(),          // double [] strength,
				cltPass3d.getSelected());         // boolean [] selection) // may be null
		measuredLayers.setLayer(                  // quad disparity
				1,                                // int       num_layer,
				cltPass3d.getDisparity(1),        // double [] disparity,
				cltPass3d.getOriginalStrength(),  // double [] strength,
				cltPass3d.getSelected());         // boolean [] selection) // may be null
		measuredLayers.setLayer(                  // hor disparity
				2,                                // int       num_layer,
				cltPass3d.getDisparity(2),        // double [] disparity,
				cltPass3d.getHorStrength(),       // double [] strength,
				cltPass3d.getSelected());         // boolean [] selection) // may be null
		measuredLayers.setLayer(                  // vert disparity
				3,                                // int       num_layer,
				cltPass3d.getDisparity(3),        // double [] disparity,
				cltPass3d.getVertStrength(),      // double [] strength,
				cltPass3d.getSelected());         // boolean [] selection) // may be null

		getDisparityHistograms(
				null,       // double  []   world_plane, // tilt equi-disparity planes to match real world planes (usually horizontal (or null)
				null,       // boolean [][] tile_sel, // null  or per-measurement layer, per-tile selection. For each layer null - do not use, {} - use all
				smplMode,   // final boolean    smplMode, //        = true;   // Use sample mode (false - regular tile mode)
				mlfp,
//				smplSide,   // final int        smplSide, //        = 2;      // Sample size (side of a square)
//				smplNum,    // final int        smplNum,  //         = 3;      // Number after removing worst
//				smplRms,    // final double     smplRms,  //         = 0.1;    // Maximal RMS of the remaining tiles in a sample
//				smplWnd,    // final boolean    smplWnd,  // use window functions for the samples

//				max_abs_tilt,  // 2.0; // Maximal absolute tilt in pixels/tile
//				max_rel_tilt,  // 0.2; // Maximal relative tilt in pixels/tile/disparity
//				damp_tilt,     //    0.001; // Damp tilt to handle insufficient  (co-linear)data
//				min_tilt_disp, // 4.0; // Disparity switch between filtering modes - near objects use tilts, far - use max disparity
//				transition,    // 1.0; // Mode transition range (between tilted and maximal disparity)
//				far_mode,      //     1;   // Far objects filtering mode (0 - off, 1 - power of disparity)
//				far_power,     //    1.0; // Raise disparity to this power before averaging for far objects

				measSel);   // calculate and blur supertiles (for all, not just selected?)
		if (tileProcessor.globalDebugLevel > 0){
			System.out.println("SuperTiles(): min_disparity = "+min_disparity+", max_disparity="+max_disparity);
			System.out.println("SuperTiles(): step_far = "+step_far+", step_near="+step_near);
			System.out.println("SuperTiles(): numBins = "+numBins+", bin_far="+bin_far+", bin_near = "+bin_near);
			System.out.println("SuperTiles(): step_threshold_far = "+step_threshold_far+", step_threshold_near="+step_threshold_near);
			for (int i = 0; i<numBins; i++){
				System.out.println(i+": "+bin_centers[i]+", "+disparityToBin(bin_centers[i]));
			}
			//
		}
		if (tileProcessor.globalDebugLevel > 0){
			String [] titles = {"d0","s0","d1","s1","d2","s2","d3","s3"};
			double [][] dbg_img = new double [titles.length][];
			for (int i = 0; i < measuredLayers.getNumLayers(); i++){
				dbg_img[2 * i] =     measuredLayers.getDisparity(i);
				dbg_img[2 * i + 1] = measuredLayers.getStrength(i);
			}
			ShowDoubleFloatArrays sdfa_instance = new ShowDoubleFloatArrays();
			sdfa_instance.showArrays(dbg_img,  tileProcessor.getTilesX(), tileProcessor.getTilesY(), true, "measuredLayers",titles);
		}
	}
	public TileProcessor getTileProcessor(){
		return  tileProcessor;
	}

	public void setTileSurface( TileSurface tileSurface)
	{
		this.tileSurface = tileSurface;
	}
	public void setTileSurface(
			GeometryCorrection   geometryCorrection)
	{
		this.tileSurface = new TileSurface(
				tileProcessor.getTileSize(),      // int tileSize,
				tileProcessor.getSuperTileSize(), // int superTileSize,
				tileProcessor.getTilesX(),        // int tilesX,
				tileProcessor.getTilesY(),        // int tilesY,
				geometryCorrection,               // GeometryCorrection geometryCorrection,
				tileProcessor.threadsMax);        // int threadsMax);
	}


	public TileSurface getTileSurface()
	{
		return tileSurface;
	}

	public void initFuseCoeff(
			double scale_diag,
			boolean debug)
	{
		final int superTileSize = tileProcessor.getSuperTileSize();

		fuse_coeff = new double [8][2][superTileSize * superTileSize];
		double [] sin2orthotab = new double [2* superTileSize];
		for (int i = 0; i < 2 * superTileSize; i++){
			sin2orthotab[i] = 0.5 * (1.0 - Math.cos(Math.PI*(i + 0.5)/superTileSize));
		}
		double [] sin2diagtab = new double [4* superTileSize];
		for (int i = 0; i < 4 * superTileSize; i++){
			sin2diagtab[i] = 0.5 * (1.0 - Math.cos(Math.PI*i/(2 * superTileSize)));
		}
		//  ortho neighbors N, E, S, W
		for (int y = 0; y < superTileSize /2 ; y++){
			for (int x = 0; x < superTileSize; x++){
				int [] indx = {
						y * superTileSize +                       x,
						x * superTileSize +                       (superTileSize - y - 1),
						(superTileSize - y - 1) * superTileSize + (superTileSize - x - 1),
						(superTileSize - x - 1) * superTileSize + y,
				};
				for (int dir = 0; dir <4; dir++){
					fuse_coeff[2*dir][0][indx[dir]] = sin2orthotab[x + superTileSize/2] * sin2orthotab[y + superTileSize / 2];
					fuse_coeff[2*dir][1][indx[dir]] = sin2orthotab[x + superTileSize/2] * sin2orthotab[y + 3 * superTileSize / 2];
				}
			}
		}
		// Diagonal: NE, SE, SW, NW
		for (int y = 0; y < superTileSize ; y++){
			for (int x = y; x < superTileSize; x++){
				int [] indx = {
						y * superTileSize +                       x,
						x * superTileSize +                       (superTileSize - y - 1),
						(superTileSize - y - 1) * superTileSize + (superTileSize - x - 1),
						(superTileSize - x - 1) * superTileSize + y,
				};

				for (int dir = 0; dir <4; dir++){
					fuse_coeff[2*dir + 1][0][indx[dir]] = scale_diag *
							sin2diagtab[x + y + superTileSize + 1] *
							sin2diagtab[ x - y + 2* superTileSize];
					fuse_coeff[2*dir + 1][1][indx[dir]] = scale_diag *
							sin2diagtab[x + y + superTileSize + 1] *
							sin2diagtab[ x - y];
				}
			}
		}
		if (debug) {
			String [] titles = {
					"N_this",  "N-other",
					"NE_this", "NE-other",
					"E_this",  "E-other",
					"SE_this", "SE-other",
					"S_this",  "S-other",
					"SW_this", "SW-other",
					"W_this",  "W-other",
					"NW_this", "NW-other"};
			double [][] dbg_img = {
					fuse_coeff[0][0],fuse_coeff[0][1],
					fuse_coeff[1][0],fuse_coeff[1][1],
					fuse_coeff[2][0],fuse_coeff[2][1],
					fuse_coeff[3][0],fuse_coeff[3][1],
					fuse_coeff[4][0],fuse_coeff[4][1],
					fuse_coeff[5][0],fuse_coeff[5][1],
					fuse_coeff[6][0],fuse_coeff[6][1],
					fuse_coeff[7][0],fuse_coeff[7][1]};
			ShowDoubleFloatArrays sdfa_instance = new ShowDoubleFloatArrays(); // just for debugging?
			sdfa_instance.showArrays(dbg_img, superTileSize, superTileSize, true, "fuse_coeff",titles);
		}
	}

	public double [] getFuseThis(int dir){
		return fuse_coeff[dir][0];
	}
	public double [] getFuseOther(int dir){
		return fuse_coeff[dir][1];
	}

	public void setShellMap(int [][] shells){
		shell_map = shells;
	}

	public int [][] getShellMap(){
		return shell_map;
	}

	public void setSurfaces(double [][] surfaces){
		this.surfaces = surfaces;
	}

	public double [][] getSurfaces(){
		return surfaces;
	}

	public int disparityToBin(double disparity){ // not filtered
		double x = disparity - min_disparity;
		int bin;
		if (x < step_threshold_far){
			bin = (int) Math.round(x/step_far);
		} else if (x < step_threshold_near){
			bin = (int) Math.round(step_threshold_far / step_far * (Math.log(x / step_threshold_far) + 1));
		} else  {
			bin = (int) Math.round(bin_near + (x - step_threshold_near) / step_near);
		}
		//				if (bin < 0) bin = 0;
		//				else if (bin >= numBins) bin = numBins - 1;
		return bin;
	}

	public double binToDisparity(double bin){
		double d;
		if (bin < bin_far ) {
			d = bin * step_far;
		} else if (bin < bin_near){
			d =  step_threshold_far * Math.exp(bin * step_far / step_threshold_far - 1.0);
		} else {
			d = step_threshold_near + (bin - bin_near) * step_near;
		}
		return min_disparity + d;
	}

	public TilePlanes.PlaneData [][] getPlanes()
	{
		return planes;
	}

	public TilePlanes.PlaneData [][] getPlanesMod()
	{
		return planes_mod;
	}

	public void resetPlanesMod()
	{
		planes_mod = null;
	}

	private double [][] getLapWeights(){
		final int superTileSize = tileProcessor.superTileSize;
		final double [][] lapWeight = new double [2 * superTileSize][2 * superTileSize];
		final double [] lapWeight1d = new double [superTileSize];
		final int superTileSize2 = 2 * superTileSize;
		for (int i = 0; i < superTileSize; i++){
			lapWeight1d[i] = 0.5*(1.0 - Math.cos((i + 0.5)* Math.PI/superTileSize));
		}
		for (int i = 0; i < superTileSize; i++){
			for (int j = 0; j < superTileSize; j++){
				lapWeight[i]                     [                     j] = lapWeight1d[i]*lapWeight1d[j];
				lapWeight[superTileSize2 - 1 - i][                     j] = lapWeight[i][j];
				lapWeight[i]                     [superTileSize2 - 1 - j] = lapWeight[i][j];
				lapWeight[superTileSize2 - 1 - i][superTileSize2 - 1 - j] = lapWeight[i][j];
			}
		}
		double s = 0.0;
		for (int i = 0; i < superTileSize2; i++){
			for (int j = 0; j < superTileSize2; j++){
				s+=lapWeight[i][j];
			}
		}
		System.out.println("getLapWeights: sum = "+s);
		return lapWeight;
	}

	public void resetDisparityHistograms()
	{
		this.disparityHistograms = null;
	}

	public void setBlurSigma(double blurSigma) // resets histograms, as it needs re-calculating
	{
		this.stBlurSigma =         blurSigma;
		this.disparityHistograms = null;
	}

	public double [][] getDisparityHistograms(
			final double [][][][] disparity_strength, // pre-calculated disparity/strength [per super-tile][per-measurement layer][2][tiles] or null
			final boolean [][]    tile_sel, // null  or per-measurement layer, per-tile selection. For each layer null - do not use, {} - use all
			final boolean         smplMode, //        = true;   // Use sample mode (false - regular tile mode)

			final MeasuredLayersFilterParameters mlfp,
//			final int             smplSide, //        = 2;      // Sample size (side of a square)
//			final int             smplNum,  //         = 3;      // Number after removing worst
//			final double          smplRms,  //         = 0.1;    // Maximal RMS of the remaining tiles in a sample
//			final boolean         smplWnd, //

//			final double          max_abs_tilt,  //  2.0;   // pix per tile
//			final double          max_rel_tilt,  //  0.2;   // (pix / disparity) per tile
//			final double          damp_tilt,     //  0.001; // Damp tilt to handle insufficient  (co-linear)data
//			final double          min_tilt_disp, //  4.0;   // Disparity switch between filtering modes - near objects use tilts, far - use max disparity
//			final double          transition,    //  1.0;   // Mode transition range (between tilted and maximal disparity)
//			final int             far_mode,      //  1;     // Far objects filtering mode (0 - off, 1 - power of disparity)
//			final double          far_power,     //  3.0;   // Raise disparity to this power before averaging for far objects

			final int             measSel)  //
	{
		if (disparity_strength != null) {
			System.out.println("getDisparityHistograms() with non-null disparity_strength");
		}
		if ((this.disparityHistograms != null) &&
				(smplMode == this.smplMode) &&
				(mlfp.equals(this.mlfp)) &&
//				(smplSide == this.smplSide) &&
//				(smplNum  == this.smplNum) &&
//				(smplRms  == this.smplRms) &&

//				(max_abs_tilt  == this.max_abs_tilt) &&
//				(max_rel_tilt  == this.max_rel_tilt) &&
//				(damp_tilt  == this.damp_tilt) &&
//				(min_tilt_disp  == this.min_tilt_disp) &&
//				(transition  == this.transition) &&
//				(far_mode  == this.far_mode) &&
//				(far_power  == this.far_power) &&
				(measSel  == this.measSel)){
				return this.disparityHistograms;
		}
		this.smplMode =      smplMode;   // Use sample mode (false - regular tile mode)
		this.mlfp = mlfp.clone();
//		this.smplSide =      smplSide;   // Sample size (side of a square)
//		this.smplNum  =      smplNum;    // Number after removing worst
//		this.smplRms  =      smplRms;    // Maximal RMS of the remaining tiles in a sample

//		this.max_abs_tilt =  max_abs_tilt;
//		this.max_rel_tilt =  max_rel_tilt;
//		this.damp_tilt =     damp_tilt;
//		this.min_tilt_disp = min_tilt_disp;
//		this.transition =    transition;
//		this.far_mode =      far_mode;
//		this.far_power =     far_power;

		this.measSel  =      measSel;

		final int tilesX =        tileProcessor.getTilesX();
		final int tilesY =        tileProcessor.getTilesY();
		final int superTileSize = tileProcessor.getSuperTileSize();
		final int stilesX = (tilesX + superTileSize -1)/superTileSize;
		final int stilesY = (tilesY + superTileSize -1)/superTileSize;
		final int nStiles = stilesX * stilesY;
		final double [][] dispHist =     new double [nStiles][]; // now it will be sparse
		final double []   strengthHist = new double [nStiles];
		final Thread[] threads = ImageDtt.newThreadArray(tileProcessor.threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
//		final TilePlanes tpl = new TilePlanes(tileProcessor.getTileSize(),superTileSize);

		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				@Override
				public void run() {

					for (int nsTile = ai.getAndIncrement(); nsTile < nStiles; nsTile = ai.getAndIncrement()) {
						int stileY = nsTile / stilesX;
						int stileX = nsTile % stilesX;
						double sw = 0.0; // sum weights
						double [] hist = new double [numBins];
						//						double [][][] disp_strength = new double [measuredLayers.getNumLayers()][][] ;
						for (int nl = 0; nl < measuredLayers.getNumLayers(); nl ++) {
							if (((measSel & (1 << nl)) != 0) && ((tile_sel == null) || ((tile_sel[nl] != null)))) {
								double [][] disp_strength;
								if (disparity_strength == null) {
									if (smplMode) {
										disp_strength =  measuredLayers.getDisparityStrengthMLTilted(
												nl,             // int num_layer,
												stileX,         // int stX,
												stileY,         // int stY,
												(((tile_sel == null) || (tile_sel[nl].length == 0))? null:tile_sel[nl]), // boolean [] sel_in,
												//											null,           // boolean [] sel_in,
												mlfp,
//												strength_floor, // double strength_floor,
//												strength_pow,   // double strength_pow,
//												smplSide,       // int        smplSide, // = 2;   // Sample size (side of a square)
//												smplNum,        //int        smplNum,   // = 3;   // Number after removing worst (should be >1)
//												smplRms,        //double     smplRms,   // = 0.1; // Maximal RMS of the remaining tiles in a sample
//												smplWnd,


//												max_abs_tilt,   //  = 2.0; // Maximal absolute tilt in pixels/tile
//												max_rel_tilt,   //  = 0.2; // Maximal relative tilt in pixels/tile/disparity
//												damp_tilt,      //  =    0.001; // Damp tilt to handle insufficient  (co-linear)data
//												min_tilt_disp,  // = 4.0; // Disparity switch between filtering modes - near objects use tilts, far - use max disparity
//												transition,     //    = 1.0; // Mode transition range (between tilted and maximal disparity)
//												far_mode,       //  =     1;   // Far objects filtering mode (0 - off, 1 - power of disparity)
//												far_power,      // =     1.0; // Raise disparity to this power before averaging for far objects

												true,          // boolean null_if_none);
												-1); // int debugLevel
									} else {
										disp_strength =  measuredLayers.getDisparityStrengthML(
												nl,             // int num_layer,
												stileX,         // int stX,
												stileY,         // int stY,
												(((tile_sel == null) || (tile_sel[nl].length == 0))? null:tile_sel[nl]), // boolean [] sel_in,
												//											null,           // boolean [] sel_in,
												mlfp.strength_floor, // double strength_floor,
												mlfp.strength_pow,   // double strength_pow,
												true);          // boolean null_if_none);
									}
								} else {
									disp_strength =  (disparity_strength[nsTile] != null) ? disparity_strength[nsTile][nl] : null;
								}
								if (disp_strength != null) {
									for (int indx = 0; indx < disp_strength[1].length; indx++) {
										double w = disp_strength[1][indx];
										if ( w > 0.0){
											double d = disp_strength[0][indx];
											int bin = disparityToBin(d);
											if ((bin >= 0) && (bin < numBins)){ // maybe collect below min and above max somewhere?
												hist[bin] += w; // +1]
												sw +=w;
											}
										}
									}
								}
							}
						}
						strengthHist[nsTile] = sw / superTileSize / superTileSize; // average strength per tile in the super-tile
						if (sw > 0){
							for (int i = 0; i<numBins; i++){
								hist[i] /= sw;
							}
							dispHist[nsTile] = hist;
						} else {
							dispHist[nsTile] = null;
						}
					}
				}
			};
		}
		ImageDtt.startAndJoin(threads);
		this.disparityHistograms = dispHist;
		this.stStrength =          strengthHist;
		if (this.stBlurSigma > 0.0) {
			blurDisparityHistogram(0); // debugLevel);
		}
		return this.disparityHistograms; // dispHist;
	}





	public void blurDisparityHistogram( // in-place
			final int  debugLevel)
	{

		final double [][] dispHist = this.disparityHistograms;
		final double      sigma =    this.stBlurSigma;
		final int sTiles = dispHist.length;
		//				final int numBins = dispHist[0].length;
		final Thread[] threads = ImageDtt.newThreadArray(tileProcessor.threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				@Override
				public void run() {
					DoubleGaussianBlur gb=new DoubleGaussianBlur();
					for (int nsTile = ai.getAndIncrement(); nsTile < sTiles; nsTile = ai.getAndIncrement()) {
						if (dispHist[nsTile] != null) {
							gb.blur1Direction(dispHist[nsTile], numBins, 1, sigma, 0.01,true);
						}
					}
				}
			};
		}
		ImageDtt.startAndJoin(threads);
	}
	// returns odd-length array of max/min (x, strength) pairs

	public double [][][] getMaxMinMax(
			final double [][][][] disparity_strength, // pre-calculated disparity/strength [per super-tile][per-measurement layer][2][tiles] or null
			final boolean [][] tile_sel // null  or per-measurement layer, per-tile selection. For each layer null - do not use, {} - use all
			){
		// first find all integer maximums, and if the top is flat - use the middle. If not flat - use 2-nd degree polynomial
		if (disparityHistograms == null) getDisparityHistograms(
//				world_plane, // double  []   world_plane, // tilt equi-disparity planes to match real world planes (usually horizontal (or null)
				disparity_strength, // pre-calculated disparity/strength [per super-tile][per-measurement layer][2][tiles] or null
				tile_sel,   // boolean [][] tile_sel, // null  or per-measurement layer, per-tile selection. For each layer null - do not use, {} - use all
				this.smplMode,   // final boolean    smplMode, //        = true;   // Use sample mode (false - regular tile mode)
				this.mlfp,
//				this.smplSide,   // final int        smplSide, //        = 2;      // Sample size (side of a square)
//				this.smplNum,    // final int        smplNum,  //         = 3;      // Number after removing worst
//				this.smplRms,    // final double     smplRms,  //         = 0.1;    // Maximal RMS of the remaining tiles in a sample
//				this.smplWnd,    // final boolean    smplWnd,  // use window functions for the samples

//				this.max_abs_tilt,  // 2.0; // Maximal absolute tilt in pixels/tile
//				this.max_rel_tilt,  // 0.2; // Maximal relative tilt in pixels/tile/disparity
//				this.damp_tilt,     //    0.001; // Damp tilt to handle insufficient  (co-linear)data
//				this.min_tilt_disp, // 4.0; // Disparity switch between filtering modes - near objects use tilts, far - use max disparity
//				this.transition,    // 1.0; // Mode transition range (between tilted and maximal disparity)
//				this.far_mode,      //     1;   // Far objects filtering mode (0 - off, 1 - power of disparity)
//				this.far_power,     //    1.0; // Raise disparity to this power before averaging for far objects

				this.measSel);
		final int globalDebugLevel = tileProcessor.globalDebugLevel;
		maxMinMax = new double [disparityHistograms.length][][];
		final Thread[] threads = ImageDtt.newThreadArray(tileProcessor.threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				@Override
				public void run() {
//					DoubleGaussianBlur gb=new DoubleGaussianBlur();
					for (int nsTile = ai.getAndIncrement(); nsTile < maxMinMax.length; nsTile = ai.getAndIncrement()) {
						if (disparityHistograms[nsTile] != null) {
							double [] dh = disparityHistograms[nsTile];
							double [][] mmm = new double [numBins][2]; // will definitely be shorter
							int numMax = 0;
							int lo = 0;
							int hi = 1;
//							if ((globalDebugLevel > -1 ) && (nsTile == 795)) {
//								System.out.println("getMaxMinMax(): nsTile="+nsTile);
//							}
							while (hi < numBins) {
								// looking for next max
								while ((hi < numBins) && (dh[hi] >= dh[hi - 1])) hi++; // flat or higher - continue
								if (hi == numBins){ // last
									if (dh[hi - 1] == dh[lo]) break; // no maximums till the very end
									if (dh[hi - 1] > dh[hi-2]) {//  and is higher than previous
										mmm[numMax * 2][0] = hi - 1;
									} else { // flat top, but higher than [lo]
										int i = hi - 3;
										while (dh[i] == dh[hi - 1]) i --;
										mmm[numMax * 2][0] = 0.5 * (hi + i); // middle
									}
									mmm[numMax * 2][1] = dh[hi - 1];
									numMax++;
									break;
								} else { // d[hi] < d[hi-1] // can be 2-nd after higher first, has sharp max at hi-1, have flat max (including from start)
									if ((dh[hi - 1] == dh[lo]) || (dh[hi - 1] == dh[hi-2])){ // just a second with higher 1-st, or a flat initial top
										// flat top
										while (dh[lo] < dh[hi - 1]) lo++;
										mmm[numMax * 2][0] = 0.5 *(hi - 1 + lo);
										mmm[numMax * 2][1] = dh[hi - 1];
									} else { // normal max, use parabola for approximation and max
										double a = 0.5*(dh[hi] -2*dh[hi-1] + dh[hi-2]);
										double b = 0.5*(dh[hi] - dh[hi-2]);
										double dx =  - b/(2*a);
										// protect agains very low a,b
										if (dx > 1.0){
											dx = 1.0;
											mmm[numMax * 2][1] = dh[hi];
											System.out.println("Max: insufficient precision, dx = "+dx+" > 1.0");
										} else if (dx < -1.0){
											dx = -1.0;
											mmm[numMax * 2][1] = dh[hi - 2];
											System.out.println("Max: insufficient precision, dx = "+dx+" < -1.0");
										} else {
											mmm[numMax * 2][1] = dh[hi-1] - b*b / (4*a);
										}
										mmm[numMax * 2][0] = (hi -1) + dx;
									}
									lo = hi-1;
									numMax++;
								}
								// look for next minimum after maximum
								while ((hi < numBins) && (dh[hi] <= dh[hi - 1])) hi++; // flat or lower - continue
								if (hi == numBins){ // last
									break;// do not need to register minimum after maximum
								} else if (dh[hi - 1] == dh[hi-2]) { // flat top
									while (dh[lo] > dh[hi - 1]) lo++;
									mmm[numMax * 2 - 1][0] = 0.5 *(hi - 1 + lo);
									mmm[numMax * 2 - 1][1] = dh[hi - 1];
								} else { // normal min - use parabola
									double a = 0.5*(dh[hi] -2*dh[hi-1] + dh[hi-2]);
									double b = 0.5*(dh[hi] - dh[hi-2]);
									double dx =  - b/(2*a);
									// protect against very low a,b
									if (dx > 1.0){
										dx = 1.0;
										mmm[numMax * 2 - 1][1] = dh[hi];
										System.out.println("Min: insufficient precision, dx = "+dx+" > 1.0");
									} else if (dx < -1.0){
										dx = -1.0;
										mmm[numMax * 2 - 1][1] = dh[hi - 2];
										System.out.println("Min: insufficient precision, dx = "+dx+" < -1.0");
									} else {
										mmm[numMax * 2 -1][1] = dh[hi-1] - b*b / (4*a);
									}
									mmm[numMax * 2 - 1][0] = (hi -1) + dx;
									lo = hi -1;
								}
							}
							if (numMax > 0) {
								maxMinMax[nsTile] = new double[numMax * 2 - 1][2];
								for (int i = 0; i < maxMinMax[nsTile].length; i++){
									maxMinMax[nsTile][i][0] = binToDisparity(mmm[i][0]); // convert to actual disparity !
									maxMinMax[nsTile][i][1] = mmm[i][1] * stStrength[nsTile];
								}
								if (globalDebugLevel > 0 ) { //************
									System.out.println(nsTile+": "+dh[0]+" "+dh[1]+" "+dh[2]+" "+dh[3]+" "+dh[4]+" "+
											dh[5]+" "+dh[6]+" "+dh[7]+" "+dh[8]+" "+dh[9]+ " "+dh[10]+" "+dh[11]+" "+dh[12]+" "+dh[13]+" "+dh[14]+" "+
											dh[15]+" "+dh[16]+" "+dh[17]+" "+dh[18]+" "+dh[19]+ " "+dh[20]);
									String dbg_str = ""+nsTile+": ";
									for (int i = 0; i < maxMinMax[nsTile].length; i++){
										dbg_str += " "+maxMinMax[nsTile][i][0]+":"+maxMinMax[nsTile][i][1];
									}
									System.out.println(dbg_str);
								}
							} else {
								maxMinMax[nsTile] = null;
							}

						} else {
							maxMinMax[nsTile] = null;
						}
					}
				}
			};
		}
		ImageDtt.startAndJoin(threads);
		return maxMinMax;
	}


	public int showDisparityHistogramWidth()
	{
		final int superTileSize = tileProcessor.superTileSize;
		int sTilesX = (tileProcessor.getTilesX() + superTileSize -1)/superTileSize;
		return sTilesX * (numBins + 1) + 1;
	}
	public double [] showDisparityHistogram()
	{
		if (disparityHistograms == null){
			getDisparityHistograms(
					null, // double  []   world_plane, // tilt equi-disparity planes to match real world planes (usually horizontal (or null)
					null,   // boolean [][] tile_sel, // null  or per-measurement layer, per-tile selection. For each layer null - do not use, {} - use all
					this.smplMode,   // final boolean    smplMode, //        = true;   // Use sample mode (false - regular tile mode)
					this.mlfp,

//					this.smplSide,   // final int        smplSide, //        = 2;      // Sample size (side of a square)
//					this.smplNum,    // final int        smplNum,  //         = 3;      // Number after removing worst
//					this.smplRms,    // final double     smplRms,  //         = 0.1;    // Maximal RMS of the remaining tiles in a sample
//					this.smplWnd,    // final boolean         smplWnd, //

//					this.max_abs_tilt,  // 2.0; // Maximal absolute tilt in pixels/tile
//					this.max_rel_tilt,  // 0.2; // Maximal relative tilt in pixels/tile/disparity
//					this.damp_tilt,     //    0.001; // Damp tilt to handle insufficient  (co-linear)data
//					this.min_tilt_disp, // 4.0; // Disparity switch between filtering modes - near objects use tilts, far - use max disparity
//					this.transition,    // 1.0; // Mode transition range (between tilted and maximal disparity)
//					this.far_mode,      //     1;   // Far objects filtering mode (0 - off, 1 - power of disparity)
//					this.far_power,     //    1.0; // Raise disparity to this power before averaging for far objects

					this.measSel);   // calculate and blur with the current settings, specified at instantiation
		}
		return showDisparityHistogram(disparityHistograms);
	}


	public double [] showDisparityHistogram(
//			double  []   world_plane, // tilt equi-disparity planes to match real world planes (usually horizontal (or null)
			final double [][][][] disparity_strength, // pre-calculated disparity/strength [per super-tile][per-measurement layer][2][tiles] or null
			boolean [][] tile_sel, // null  or per-measurement layer, per-tile selection. For each layer null - do not use, {} - use all
			boolean    smplMode, //        = true;   // Use sample mode (false - regular tile mode)
			MeasuredLayersFilterParameters mlfp, // parameters
//			int        smplSide, //        = 2;      // Sample size (side of a square)
//			int        smplNum,  //         = 3;      // Number after removing worst
//			double     smplRms,  //         = 0.1;    // Maximal RMS of the remaining tiles in a sample
//			boolean    smplWnd, //

//  			double     max_abs_tilt,  //  2.0;   // pix per tile
//			double     max_rel_tilt,  //  0.2;   // (pix / disparity) per tile
//			double     damp_tilt,     //  0.001; // Damp tilt to handle insufficient  (co-linear)data
//			double     min_tilt_disp, //  4.0;   // Disparity switch between filtering modes - near objects use tilts, far - use max disparity
//			double     transition,    //  1.0;   // Mode transition range (between tilted and maximal disparity)
//			int        far_mode,      //  1;     // Far objects filtering mode (0 - off, 1 - power of disparity)
//			double     far_power,     //  3.0;   // Raise disparity to this power before averaging for far objects

			int        measSel) // bitmask of the selected measurements for supertiles : +1 - combo, +2 - quad +4 - hor +8 - vert
	{
		getDisparityHistograms( // will recalculate if does not exist or some parameters changed
//				world_plane, // double  []   world_plane, // tilt equi-disparity planes to match real world planes (usually horizo
				disparity_strength, // pre-calculated disparity/strength [per super-tile][per-measurement layer][2][tiles] or nullntal (or null)
				tile_sel,   // boolean [][] tile_sel, // null  or per-measurement layer, per-tile selection. For each layer null - do not use, {} - use all
				smplMode,   // final boolean    smplMode, //        = true;   // Use sample mode (false - regular tile mode)
				mlfp,
//				smplSide,   // final int        smplSide, //        = 2;      // Sample size (side of a square)
//				smplNum,    // final int        smplNum,  //         = 3;      // Number after removing worst
//				smplRms,    // final double     smplRms,  //         = 0.1;    // Maximal RMS of the remaining tiles in a sample
//				smplWnd,    // final boolean         smplWnd, //

//				max_abs_tilt,  // 2.0; // Maximal absolute tilt in pixels/tile
//				max_rel_tilt,  // 0.2; // Maximal relative tilt in pixels/tile/disparity
//				damp_tilt,     //    0.001; // Damp tilt to handle insufficient  (co-linear)data
//				min_tilt_disp, // 4.0; // Disparity switch between filtering modes - near objects use tilts, far - use max disparity
//				transition,    // 1.0; // Mode transition range (between tilted and maximal disparity)
//				far_mode,      //     1;   // Far objects filtering mode (0 - off, 1 - power of disparity)
//				far_power,     //    1.0; // Raise disparity to this power before averaging for far objects

				measSel); // calculate and blur with the current settings, specified at instantiation
		return showDisparityHistogram(disparityHistograms);
	}

	public double [] showDisparityHistogram(final double [][] dispHist){
		final int superTileSize = tileProcessor.superTileSize;
		int         sTileHeight0 = 0; // vertical pixels for each  histogram (excluding borders). if <= 0 make square cells
		final double [] strengthHist = stStrength;
		final int sTilesX = (tileProcessor.getTilesX() + superTileSize -1)/superTileSize;

		final int sTiles = dispHist.length;
		final int sTilesY = sTiles / sTilesX;
		//				final int numBins = dispHist[0].length; // [0] - weight
		final int sTileHeight = (sTileHeight0 > 0)? sTileHeight0 : numBins;

		final double [] maxHist = new double [sTiles];
		final int width = sTilesX * (numBins + 1) + 1;
		final int height = sTilesY * (sTileHeight + 1) +1;
		double [] rslt = new double [width*height];
		double maxW = 0.0; // use for borders between cells
		for (int i = 0; i < sTiles; i++){
			if (maxW < strengthHist[i]) maxW = strengthHist[i];
			if (dispHist[i] != null) {
				for (int j = 0; j< numBins; j++){
					if (!Double.isNaN( dispHist[i][j]) && (maxHist[i] < dispHist[i][j])) maxHist[i] = dispHist[i][j];
				}
			}
		}
		for (int nsTile = 0; nsTile < sTiles; nsTile++){
			int stileY = nsTile / sTilesX;
			int stileX = nsTile % sTilesX;
			int x0 = stileX * (numBins + 1);
			int y0 = stileY * (numBins + 1);
			int indx0 = x0 + y0*width;

			// draw rectangular frame - horisontal dotted lines
			for (int j = 0; j < numBins + 2; j+=2) {
				rslt [indx0 + j] =                            maxW;
				rslt [indx0 + (sTileHeight + 1) * width+ j] = maxW;
			}
			// vertical dotted lines
			for (int j = 0; j < sTileHeight + 2; j+=2) {
				rslt [indx0 + j * width] =                    maxW;
				rslt [indx0 + (numBins + 1)+  j * width] =    maxW;
			}
			// draw normalized histograms, using overall weight as intensity
			if (dispHist[nsTile] != null) {
				for (int bin = 0; bin <numBins; bin ++){
					int h = (int) (sTileHeight * dispHist[nsTile][bin] / maxHist[nsTile]);
					int x = bin + 1;
					for (int j = 0; j <= h;  j++) {
						int y =  sTileHeight + 1 - j;
						rslt [indx0 + y * width + x] =           strengthHist[nsTile];
					}
				}
			}
		}
		return rslt;
	}

	public double [] showMaxMinMax(){
			return showMaxMinMax( // calculate and blur with the current settings, specified at instantiation
					null,
					null);
	}
	public double [] showMaxMinMax(
//			double  []   world_plane, // tilt equi-disparity planes to match real world planes (usually horizontal (or null)
			double [][][][] disparity_strength, // pre-calculated disparity/strength [per super-tile][per-measurement layer][2][tiles] or null

			boolean [][] tile_sel // null  or per-measurement layer, per-tile selection. For each layer null - do not use, {} - use all
		){
		if (maxMinMax == null){
			getMaxMinMax( // calculate and blur with the current settings, specified at instantiation
//					world_plane, // tilt equi-disparity planes to match real world planes (usually horizontal (or null)
					disparity_strength, // pre-calculated disparity/strength [per super-tile][per-measurement layer][2][tiles] or null
					tile_sel); // null  or per-measurement layer, per-tile selection. For each layer null - do not use, {} - use all

		}
		final int superTileSize = tileProcessor.superTileSize;

		int         sTileHeight0 = 0; // vertical pixels for each  histogram (excluding borders). if <= 0 make square cells
		final double [] strengthHist = stStrength;
		final int sTilesX = (tileProcessor.getTilesX() + superTileSize -1)/superTileSize;

		//				final int sTiles = disparityHistograms.length;
		final int sTiles = maxMinMax.length;
		final int sTilesY = sTiles / sTilesX;
		//				final int numBins = dispHist[0].length; // [0] - weight
		final int sTileHeight = (sTileHeight0 > 0)? sTileHeight0 : numBins;

		final double [] maxHist = new double [sTiles];
		final int width = sTilesX * (numBins + 1) + 1;
		final int height = sTilesY * (sTileHeight + 1) +1;
		final boolean [] isMax = new boolean [numBins];
		final boolean [] isMin = new boolean [numBins];

		double [] rslt = new double [width*height];
		double maxW = 0.0; // use for borders between cells
		for (int i = 0; i < sTiles; i++){
			if (maxW < strengthHist[i]) maxW = strengthHist[i];
			if (disparityHistograms[i] != null) {
				for (int j = 0; j< numBins; j++){
					if (!Double.isNaN( disparityHistograms[i][j]) && (maxHist[i] < disparityHistograms[i][j])) maxHist[i] = disparityHistograms[i][j];
				}
			}
		}
		for (int nsTile = 0; nsTile < sTiles; nsTile++){
//			if (nsTile == 795){
//				System.out.println("showMaxMinMax(), nsTile="+nsTile);
//			}
			int stileY = nsTile / sTilesX;
			int stileX = nsTile % sTilesX;
			int x0 = stileX * (numBins + 1);
			int y0 = stileY * (numBins + 1);
			int indx0 = x0 + y0*width;

			// draw rectangular frame - horizontal dotted lines
			for (int j = 0; j < numBins + 2; j+=2) {
				rslt [indx0 + j] =                            maxW;
				rslt [indx0 + (sTileHeight + 1) * width+ j] = maxW;
			}
			// vertical dotted lines
			for (int j = 0; j < sTileHeight + 2; j+=2) {
				rslt [indx0 + j * width] =                    maxW;
				rslt [indx0 + (numBins + 1)+  j * width] =    maxW;
			}
			// draw normalized histograms, using overall weight as intensity
			if ((disparityHistograms[nsTile] != null) && (maxMinMax[nsTile] != null)) {
				for (int bin = 0; bin <numBins; bin ++){
					isMax[bin] = false;
					isMin[bin] = false;
				}
				for (int i = 0; i <maxMinMax[nsTile].length; i++){
//					int imm = (int) Math.round(maxMinMax[nsTile][i][0]);
					int imm = disparityToBin(maxMinMax[nsTile][i][0]);
					if      (imm <0 )        imm = 0;
					else if (imm >= numBins) imm = numBins - 1;
					if ((i & 1) == 0) isMax[imm] = true;
					else              isMin[imm] = true;
				}

				for (int bin = 0; bin <numBins; bin ++){
					int h = (int) (sTileHeight * disparityHistograms[nsTile][bin] / maxHist[nsTile]);
					int x = bin + 1;
					int jMin = isMin[bin] ? h : 0;
					int jMax = isMax[bin] ? h : sTileHeight;

					if (isMax[bin] || isMin[bin]) {
						for (int j = jMin; j <= jMax;  j++) {
							int y =  sTileHeight + 1 - j;
							rslt [indx0 + y * width + x] = strengthHist[nsTile];
						}
					}
				}
			}
		}
		return rslt;
	}

	// updates bgDisparity, bgStrength
	public double [][] getBgDispStrength(
			final double minBgDisparity,
			final double minBgFract)
	{
////		final double step_disparity = step_near; // TODO: implement


		if (maxMinMax == null) return null;
		final int superTileSize = tileProcessor.superTileSize;
		final int globalDebugLevel = tileProcessor.globalDebugLevel;
		final int sTiles = maxMinMax.length;
		final Thread[] threads = ImageDtt.newThreadArray(tileProcessor.threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		bgDisparity = new double[sTiles];
		bgStrength =  new double[sTiles];
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				@Override
				public void run() {
					for (int nsTile = ai.getAndIncrement(); nsTile < sTiles; nsTile = ai.getAndIncrement()) {
						if (nsTile == 49) { // 414){ // 331){
							System.out.println("getBgDispStrength(); nsTIle="+nsTile);
						}
						double [][] mmm = maxMinMax[nsTile];
						bgDisparity[nsTile] = Double.NaN;
						bgStrength[nsTile] =  0.0;
						if (mmm != null){
							int maxNum = 0;
							int numMax = (maxMinMax[nsTile].length + 1) / 2;
							double selStrenth = 0.0;
							int startIndex = 0;
							for (maxNum = 0; maxNum < numMax; maxNum++){
								if (mmm[2 * maxNum][0] >= minBgDisparity){
									if (selStrenth == 0.0) startIndex = maxNum; // will keep first non-zero maximum number
									selStrenth += mmm[2 * maxNum][1];
								}
							}
							if (selStrenth > 0.0){
								selStrenth *= minBgFract;
								double accumStrength = 0.0;
								for (maxNum = startIndex; maxNum < numMax; maxNum++){
									accumStrength += mmm[2 * maxNum][1];
									if (accumStrength >= selStrenth){
										break;
									}
								}
								if (maxNum >= numMax){
									maxNum = numMax - 1; // probably just wrong fraction (>1.0)
								}
								// if unlikely there are several maximums before minBgFract - use the strongest
								int maxIndex = startIndex;
								if (startIndex < maxNum){
									for (startIndex++; startIndex < numMax ;startIndex++){
										if (mmm[2 * startIndex][1] > mmm[2 * maxIndex][1]) maxIndex = startIndex;
									}
								}
								//maxIndex is what we need. Which strength to use - individual or accumulated? Use individual.
								///										bgDisparity[nsTile] = min_disparity + mmm[2 * maxIndex][0] * step_disparity;
								bgDisparity[nsTile] = binToDisparity (mmm[2 * maxIndex][0]);
								bgStrength[nsTile] =  mmm[2 * maxIndex][1];
							}
						}
					}
				}
			};
		}
		ImageDtt.startAndJoin(threads);
		final double [][] bgDispStrength = {bgDisparity, bgStrength};
		if (globalDebugLevel > 0) {
			ShowDoubleFloatArrays sdfa_instance = new ShowDoubleFloatArrays(); // just for debugging?
			final int stilesX = (tileProcessor.getTilesX() + superTileSize -1)/superTileSize;
			final int stilesY = (tileProcessor.getTilesY()  + superTileSize -1)/superTileSize;
			sdfa_instance.showArrays(bgDispStrength, stilesX, stilesY, true, "bgDispStrength");
		}
		return bgDispStrength;
	}

	// from per-super-tile disparity/strength interpolate per-tile disparity/strength using same sine-based window
	public double [][] getBgTileDispStrength()
	{
		if ((this.bgDisparity == null) || (this.bgStrength == null)) return null;
		final int tilesX = tileProcessor.getTilesX();
		final int tilesY = tileProcessor.getTilesY();
		final int superTileSize = tileProcessor.superTileSize;
		final int stilesX = (tilesX + superTileSize -1)/superTileSize;
		final int stilesY = (tilesY + superTileSize -1)/superTileSize;
		final int nStiles = stilesX * stilesY;
		final Thread[] threads = ImageDtt.newThreadArray(tileProcessor.threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		final int superTileSize2 = 2 * superTileSize;
		final double [][] lapWeight = getLapWeights();
		final double [] tileDisparity = new double [tilesY * tilesX]; // assuming all 0.0
		final double [] tileStrength =  new double [tilesY * tilesX]; // assuming all 0.0

		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				@Override
				public void run() {
					for (int nsTile = ai.getAndIncrement(); nsTile < nStiles; nsTile = ai.getAndIncrement()) {
						int stileY = nsTile / stilesX;
						int stileX = nsTile % stilesX;
						int tY0 = stileY * superTileSize;
						int tX0 = stileX * superTileSize;
						double [][] lapWeight_dbg = lapWeight;
						if (nsTile ==  414) { //317) { // 755){
							System.out.println("getBgTileDispStrength(x): stileY="+stileY+", stileX="+stileX+" lapWeight_dbg.length="+lapWeight_dbg.length);
						}
						for (int tY = 0; tY < superTileSize; tY++){
							int tileY = tY0 +tY;
							if (tileY < tilesY) {
								for (int tX = 0; tX < superTileSize; tX++){
									int tileX = tX0 +tX;
									if (tileX < tilesX) {
										int tIndex =  tileY * tilesX + tileX;
										// iterate for +/- 1 supterile around current, acummulate disparity/strength
										double sd = 0.0, sw = 0.0;
										for (int stDY = -1; stDY <=1; stDY ++){
											int stY =  stileY + stDY;
											int dtY =tY + superTileSize/2 -superTileSize * stDY;
											if ((stY >= 0) && (stY < stilesY) && (dtY >= 0) && (dtY < superTileSize2)){
												for (int stDX = -1; stDX <=1; stDX ++){
													int stX =  stileX + stDX;
													int dtX =tX + superTileSize/2 -superTileSize * stDX;
													if ((stX >= 0) && (stX < stilesX) && (dtX >= 0) && (dtX < superTileSize2)){
														int stIndex = stY * stilesX + stX;
														double w = bgStrength[stIndex] * lapWeight[dtY][dtX];
														if (nsTile ==  415) {
															System.out.println("tX="+tX+", tY="+tY+", tileX="+tileX+", tileY="+tileY+" stDX="+stDX+" stDY="+stDY+
																	", bgStrength["+stIndex+"] = "+bgStrength[stIndex]+
																	", lapWeight["+dtY+"]["+dtX+"]="+lapWeight[dtY][dtX]+" stY="+stY+" stX="+stX+" w="+w);
														}
														sw += w;
														if (w >0.0) sd += w * bgDisparity[stIndex]; // so NaN will be OK
													}
												}
											}
										}
										tileStrength[tIndex] = sw;
										if (sw > 0.0) {
											tileDisparity[tIndex] = sd/sw;
										} else {
											tileDisparity[tIndex] = Double.NaN;
										}
										if (nsTile ==  415) {
											System.out.println("sw= "+sw+", sd="+sd);
										}
									}
								}
							}
						}
					}
				}
			};
		}
		ImageDtt.startAndJoin(threads);
		final double [][] bgTileDispStrength = {tileDisparity, tileStrength};
		return bgTileDispStrength;
	}


	public double corrMaxEigen(
			double maxEigen,
			double dispNorm,
			TilePlanes.PlaneData pd)
	{
		return corrMaxEigen(
				maxEigen,
				dispNorm,
				pd.getZxy()[0]);
	}

	public double corrMaxEigen(
			double maxEigen,
			double dispNorm,
			double disparity)
	{
		double corrV = maxEigen;
		if ((dispNorm > 0.0) && (disparity > dispNorm)) {
			//					double dd = (dispNorm + z0)/ dispNorm; // > 1
			double dd = disparity/ dispNorm; // > 1
			corrV *= dd * dd; // > original
		}
		return corrV;

	}

	public int [][] getTransMatrix(
			boolean [][][] selections ) // for each plane should have the same non-null ml
	{
		int num_ml = measuredLayers.getNumLayers();
		int num_p = selections.length;
		int superTileSize = tileProcessor.getSuperTileSize();
		int st2 = 2 * superTileSize;
		int [][] trans_mat = new int [num_p][num_p];
		for (int ml = 0; ml < num_ml; ml++) if (selections[0][ml] != null){
			for (int y = 0; y < st2; y++){
				for (int x = 0; x < (st2 - 1); x++){
					int indx1 = y * st2 + x;
					int indx2 = y * st2 + x + 1;
					for (int np1 = 0; np1 < num_p; np1++){
						for (int np2 = 0; np2 < num_p; np2++){
							if (selections[np1][ml][indx1] && selections[np2][ml][indx2]){
								trans_mat[np1][np2]++;
							}
						}
					}
				}
			}
			for (int x = 0; x < st2; x++){
				for (int y = 0; y < (st2 - 1); y++){
					int indx1 = y * st2 + x;
					int indx2 = y * st2 + x + st2;
					for (int np1 = 0; np1 < num_p; np1++){
						for (int np2 = 0; np2 < num_p; np2++){
							if (selections[np1][ml][indx1] && selections[np2][ml][indx2]){
								trans_mat[np1][np2]++;
							}
						}
					}
				}
			}
		}
		return trans_mat;
	}

	public double [][] getTransRel(
			int [][] trans_matrix)
	{
		double [][] trans_rel = new double [trans_matrix.length][trans_matrix.length];
		for (int i = 0; i < trans_matrix.length; i++){
			for (int j = 0; j < trans_matrix[i].length; j++){
				if ((trans_matrix[i][i] + trans_matrix[j][j]) > 0){
					trans_rel[i][j] = (trans_matrix[i][j] + trans_matrix[j][i]);
					trans_rel[i][j] /= (trans_matrix[i][i] + trans_matrix[j][j]);
				}
			}
		}
		return trans_rel;
	}

	public String [] showSupertileSeparationTitles(
			double [][][]  disp_strength,
			boolean [][][] selections)
	{
		return 	showSupertileSeparationTitles(disp_strength, selections, null);
	}

	public double [][] showSupertileSeparation(
			boolean        useWorld,
			double [][][]  disp_strength,
			boolean [][][] selections)
	{
		return 	showSupertileSeparation(useWorld, disp_strength, selections, null);
	}

	public String [] showSupertileSeparationTitles(
			double [][][]  disp_strength,
			boolean [][][] selections,
			TilePlanes.PlaneData [] planes)
	{
		int num_ml = disp_strength.length;
		int num_p  = (selections == null) ? 0: selections.length;
		int num_pm = num_ml * num_p;
		int num_pd = (planes != null) ? (planes.length - LOWEST_PLANE(planes.length)) : 0;
		String [] titles = new String [num_pm + 3 * num_ml + 8 * num_pd];
		for (int np = 0; np < num_p; np++){
			for (int ml = 0; ml < num_ml; ml++){
				titles [np * num_ml + ml] = "p"+np+"_l"+ml;
			}
		}
		for (int ml = 0; ml < num_ml; ml++){
			titles [num_pm + 0 * num_ml + ml] = "disp_"+ml;
			titles [num_pm + 1 * num_ml + ml] = "strn_"+ml;
			titles [num_pm + 2 * num_ml + ml] = "sel_"+ml;
		}
		for (int npd = 0; npd < num_pd; npd++){
			titles [num_pm + 3 * num_ml + 0 * num_pd + npd] = "pd_"+npd;
			titles [num_pm + 3 * num_ml + 1 * num_pd + npd] = "pdm_"+npd;
			titles [num_pm + 3 * num_ml + 2 * num_pd + npd] = "pdd_"+npd;
			titles [num_pm + 3 * num_ml + 3 * num_pd + npd] = "pwd_"+npd;
			titles [num_pm + 3 * num_ml + 4 * num_pd + npd] = "pds_"+npd;
			titles [num_pm + 3 * num_ml + 5 * num_pd + npd] = "pws_"+npd;
			titles [num_pm + 3 * num_ml + 6 * num_pd + npd] = "-pdd_"+npd;
			titles [num_pm + 3 * num_ml + 7 * num_pd + npd] = "-pwd_"+npd;
		}

	// TilePlanes.PlaneData
		return titles;
	}

	public double [][] showSupertileSeparation(
			boolean useWorld,
			double [][][]  disp_strength,
			boolean [][][] selections,
			TilePlanes.PlaneData [] planes
			)
	{
		int superTileSize = tileProcessor.getSuperTileSize();
		int num_ml = disp_strength.length;
		int num_p  = (selections == null) ? 0: selections.length;
		int num_pm = num_ml * num_p;
		int num_pd = (planes != null) ? (planes.length - LOWEST_PLANE(planes.length)) : 0;
		final double [][] lapWeight = getLapWeights();
		double [][] data = new double [num_pm + 3 * num_ml + 8 * num_pd ][]; // 4* superTileSize*superTileSize];
		for (int np = 0; np < num_p; np++) if (selections [np] != null){
			for (int ml = 0; ml < num_ml; ml++) if ((disp_strength[ml]!=null) && (selections[np][ml] != null)){
				int nd = np * num_ml + ml;
				data[nd] = new double[4 * superTileSize * superTileSize];
				for (int i = 0; i < data[nd].length; i++){
					if (selections[np][ml][i] && (disp_strength[ml][1][i] > 0.0)){
						data[nd][i] = disp_strength[ml][0][i];
					} else {
						data[nd][i] = Double.NaN;
					}
				}
			}
		}
		for (int ml = 0; ml < num_ml; ml++) if (disp_strength[ml]!=null){
			int nd = num_pm + 2 * num_ml + ml;

			data [nd] = new double [4* superTileSize*superTileSize];
			int [] multi_sel = new int [4* superTileSize*superTileSize];
			for (int i = 0; i < data[nd].length; i++){
				data [nd][i] = Double.NaN;
				for (int np = 0; np < num_p; np++) if ((selections [np] != null) && (selections [np][ml] != null) && selections [np][ml][i]){
//					data [nd][i] = np + 1;
					multi_sel[i] |= (1 << np);
					data [nd][i] = multi_sel[i];
//					break;
				}
			}
			data [num_pm + 0 * num_ml + ml] = disp_strength[ml][0];
			data [num_pm + 1 * num_ml + ml] = disp_strength[ml][1].clone();
			// undo lapweight to show
			for (int sty = 0; sty < 2 * superTileSize; sty++){
				for (int stx = 0; stx < 2 * superTileSize; stx++){
					data [num_pm + 1 * num_ml + ml][stx + 2 * superTileSize * sty] /= lapWeight[sty][stx];
				}
			}

		}
		for (int npd = 0; npd < num_pd; npd++) if (planes[npd +LOWEST_PLANE(planes.length)] != null){
/*			double [][] ellipsoids = planes[npd +LOWEST_PLANE(planes.length)].getDoublePlaneDisparityStrength(
					useWorld,
					null, // double [] window,
					true, // boolean   use_sel,
					true, // boolean   divide_by_area,
					1.5, // double   scale_projection,
					1); // int       debugLevel) */
			double [][] ellipsoids = planes[npd +LOWEST_PLANE(planes.length)].getDoublePlaneDisparityStrength(
					useWorld,
					null, // double [] window,
					-1, //
					true, // boolean   use_sel,
					true, // boolean   divide_by_area,
					1.5, // double   scale_projection,
					0.0, // double    fraction_uni,
					1); // int       debugLevel)
			double [][] ellipsoidsW = planes[npd +LOWEST_PLANE(planes.length)].getDoublePlaneWorldDisparityStrength(
					null, // double [] window,
					-1, //
					true, // boolean   use_sel,
					true, // boolean   divide_by_area,
					1.5, // double   scale_projection,
					0.0, // double    fraction_uni,
					1); // int       debugLevel)
			data [num_pm + 3 * num_ml + 0 * num_pd + npd] = planes[npd +LOWEST_PLANE(planes.length)].getDoublePlaneDisparity(
					useWorld,
					false);

			boolean [][] msel = planes[npd +LOWEST_PLANE(planes.length)].getMeasSelection();
			int ntiles = data [num_pm + 3 * num_ml + 0 * num_pd + npd].length;
			data [num_pm + 3 * num_ml + 1 * num_pd + npd] = new double [ntiles];
			for (int i = 0; i < ntiles; i++) {
				data [num_pm + 3 * num_ml + 1 * num_pd + npd][i] = Double.NaN;
				for (int ml = 0; ml < num_ml; ml++) if ((disp_strength[ml]!=null) && (msel[ml] != null)){
					if (msel[ml][i] && (disp_strength[ml][1][i] > 0.0)){
						data [num_pm + 3 * num_ml + 1 * num_pd + npd][i] = data [num_pm + 3 * num_ml + 0 * num_pd + npd][i];
					}
				}
			}
			if (ellipsoids != null) {
				data [num_pm + 3 * num_ml + 2 * num_pd + npd] = ellipsoids[0];
				data [num_pm + 3 * num_ml + 4 * num_pd + npd] = ellipsoids[1];
				data [num_pm + 3 * num_ml + 6 * num_pd + npd] = new double[ntiles];
				for (int i = 0; i < ntiles; i++) {
					for (int ml = 0; ml < num_ml; ml++) if ((disp_strength[ml]!=null) && (msel[ml] != null)){
						if (msel[ml][i] && (disp_strength[ml][1][i] > 0.0)){
							data [num_pm + 3 * num_ml + 6 * num_pd + npd][i] = disp_strength[ml][0][i] - ellipsoids[0][i];
						} else {
							data [num_pm + 3 * num_ml + 6 * num_pd + npd][i] = Double.NaN;
						}
						break;
					}
				}
			}
			if (ellipsoidsW != null) {
				data [num_pm + 3 * num_ml + 3 * num_pd + npd] = ellipsoidsW[0];
				data [num_pm + 3 * num_ml + 5 * num_pd + npd] = ellipsoidsW[1];
				data [num_pm + 3 * num_ml + 7 * num_pd + npd] = new double[ntiles];
				for (int i = 0; i < ntiles; i++) {
					for (int ml = 0; ml < num_ml; ml++) if ((disp_strength[ml]!=null) && (msel[ml] != null)){
						if (msel[ml][i] && (disp_strength[ml][1][i] > 0.0)){
							data [num_pm + 3 * num_ml + 7 * num_pd + npd][i] = disp_strength[ml][0][i] - ellipsoidsW[0][i];
						} else {
							data [num_pm + 3 * num_ml + 7 * num_pd + npd][i] = Double.NaN;
						}
						break;
					}
				}
			}
		}
		return data;
	}

	public String [] showSupertileWorldTitles(
			double [][][]  disp_strength,
			boolean [][][] selections,
			TilePlanes.PlaneData [] planes)
	{
//		int num_p  = selections.length;
		int num_pd = (planes != null) ? (planes.length - LOWEST_PLANE(planes.length)) : 0;
		String [] titles = new String [9 * num_pd];
		for (int npd = 0; npd < num_pd; npd++){
			titles [0 * num_pd + npd] = "tile_x_"+npd;
			titles [1 * num_pd + npd] = "tile_y_"+npd;
			titles [2 * num_pd + npd] = "tile_z_"+npd;
			titles [3 * num_pd + npd] = "plane_x_"+npd;
			titles [4 * num_pd + npd] = "plane_y_"+npd;
			titles [5 * num_pd + npd] = "plane_z_"+npd;
			titles [6 * num_pd + npd] = "diff_x_"+npd;
			titles [7 * num_pd + npd] = "diff_y_"+npd;
			titles [8 * num_pd + npd] = "diff_z_"+npd;
		}

	// TilePlanes.PlaneData
		return titles;
	}
/*
	public double [][] showSupertileWorld(
			double [][][]  disp_strength,
			boolean [][][] selections,
			TilePlanes.PlaneData [] planes
			)
	{
		int superTileSize = tileProcessor.getSuperTileSize();
		int num_ml = disp_strength.length;
		int num_p  = (selections == null) ? 0: selections.length;
		int num_pm = num_ml * num_p;
		int num_pd = (planes != null) ? (planes.length - LOWEST_PLANE(planes.length)) : 0;
		int ml = 0;
		for (ml = 0; ml < num_ml; ml++) if (disp_strength[ml]!=null){
			break;
		}
		double [][] data = new double [num_pd][];
		for (int npd = 0; npd < num_pd; npd++) if (planes[npd +LOWEST_PLANE(planes.length)] != null){
			double [][] ellipsoids = planes[npd +LOWEST_PLANE(planes.length)].getDoublePlaneDisparityStrength(
					useWorld,
					null, // double [] window,
					-1, //
					true, // boolean   use_sel,
					true, // boolean   divide_by_area,
					1.5, // double   scale_projection,
					0.0, // double    fraction_uni,
					1); // int       debugLevel)
			double [][] ellipsoidsW = planes[npd +LOWEST_PLANE(planes.length)].getDoublePlaneWorldDisparityStrength(
					null, // double [] window,
					-1, //
					true, // boolean   use_sel,
					true, // boolean   divide_by_area,
					1.5, // double   scale_projection,
					0.0, // double    fraction_uni,
					1); // int       debugLevel)
			data [num_pm + 3 * num_ml + 0 * num_pd + npd] = planes[npd +LOWEST_PLANE(planes.length)].getDoublePlaneDisparity(
					useWorld,
					false);

			boolean [][] msel = planes[npd +LOWEST_PLANE(planes.length)].getMeasSelection();
			int ntiles = data [num_pm + 3 * num_ml + 0 * num_pd + npd].length;
			data [num_pm + 3 * num_ml + 1 * num_pd + npd] = new double [ntiles];
			for (int i = 0; i < ntiles; i++) {
				data [num_pm + 3 * num_ml + 1 * num_pd + npd][i] = Double.NaN;
				for (int ml = 0; ml < num_ml; ml++) if ((disp_strength[ml]!=null) && (msel[ml] != null)){
					if (msel[ml][i] && (disp_strength[ml][1][i] > 0.0)){
						data [num_pm + 3 * num_ml + 1 * num_pd + npd][i] = data [num_pm + 3 * num_ml + 0 * num_pd + npd][i];
					}
				}
			}
			if (ellipsoids != null) {
				data [num_pm + 3 * num_ml + 2 * num_pd + npd] = ellipsoids[0];
				data [num_pm + 3 * num_ml + 4 * num_pd + npd] = ellipsoids[1];
				data [num_pm + 3 * num_ml + 6 * num_pd + npd] = new double[ntiles];
				for (int i = 0; i < ntiles; i++) {
					for (int ml = 0; ml < num_ml; ml++) if ((disp_strength[ml]!=null) && (msel[ml] != null)){
						if (msel[ml][i] && (disp_strength[ml][1][i] > 0.0)){
							data [num_pm + 3 * num_ml + 6 * num_pd + npd][i] = disp_strength[ml][0][i] - ellipsoids[0][i];
						} else {
							data [num_pm + 3 * num_ml + 6 * num_pd + npd][i] = Double.NaN;
						}
						break;
					}
				}
			}
			if (ellipsoidsW != null) {
				data [num_pm + 3 * num_ml + 3 * num_pd + npd] = ellipsoidsW[0];
				data [num_pm + 3 * num_ml + 5 * num_pd + npd] = ellipsoidsW[1];
				data [num_pm + 3 * num_ml + 7 * num_pd + npd] = new double[ntiles];
				for (int i = 0; i < ntiles; i++) {
					for (int ml = 0; ml < num_ml; ml++) if ((disp_strength[ml]!=null) && (msel[ml] != null)){
						if (msel[ml][i] && (disp_strength[ml][1][i] > 0.0)){
							data [num_pm + 3 * num_ml + 7 * num_pd + npd][i] = disp_strength[ml][0][i] - ellipsoidsW[0][i];
						} else {
							data [num_pm + 3 * num_ml + 7 * num_pd + npd][i] = Double.NaN;
						}
						break;
					}
				}
			}
		}
		return data;
	}
*/





	// calculate "tilted" disparity, so planes parallel to the same world plane would have the same disparity
	// also produces non-tilted, if world_plane_norm == null
	// Protecting from behind the horizon - set strength of all tiles in the negative disparity area to 0
	public double [][][][] getPlaneDispStrengthsST(
			final double []  world_plane_norm, // real world normal vector to a suggested plane family (0,1,0) for horizontal planes, (0,0,1) const disparity, null - floating
			final int        stMeasSel, //            = 1;      // Select measurements for supertiles : +1 - combo, +2 - quad +4 - hor +8 - vert
			final boolean    plPreferDisparity, // Always start with disparity-most axis (false - lowest eigenvalue)
			final GeometryCorrection geometryCorrection,
			final boolean    correct_distortions,
			final boolean    smplMode, //        = true;   // Use sample mode (false - regular tile mode)
			final MeasuredLayersFilterParameters mlfp,
			final int        debugLevel,
			final int        dbg_X,
			final int        dbg_Y)
	{
		final boolean const_disparity = (world_plane_norm != null) && (world_plane_norm[0] == 0.0)  && (world_plane_norm[1] == 0.0);
		final boolean floating =  (world_plane_norm == null);
		final int tilesX =        tileProcessor.getTilesX();
		final int tilesY =        tileProcessor.getTilesY();
		final int superTileSize = tileProcessor.getSuperTileSize();
		final int tileSize =      tileProcessor.getTileSize();

		final int stilesX = (tilesX + superTileSize -1)/superTileSize;
		final int stilesY = (tilesY + superTileSize -1)/superTileSize;
		final int nStiles = stilesX * stilesY;
//		final Thread[] threads = ImageDtt.newThreadArray(tileProcessor.threadsMax);
		final Thread[] threads = ImageDtt.newThreadArray((debugLevel > 1)? 1 : tileProcessor.threadsMax);

		final AtomicInteger ai = new AtomicInteger(0);
		final int debug_stile = (debugLevel > -1)? (dbg_Y * stilesX + dbg_X):-1;
		// TODO: Remove when promoting PlaneData
		final TilePlanes tpl = new TilePlanes(tileSize,superTileSize, geometryCorrection);
		final double [][][][] plane_disp_strength = new double [nStiles][][][];
		// DEBUG feature:
		final double [][] zero_tilts = {{0.0,0.0}}; // set to null for float
		if (debugLevel > -1) {
			System.out.println("getPlaneDispStrengthsST(x): debugLevel == "+debugLevel);
		}
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				@Override
				public void run() {
					double [][][] plane_tilts = null; // used only for world_plane_norm != null
					// to get tile disparities needed to calculate tilts
					for (int nsTile = ai.getAndIncrement(); nsTile < nStiles; nsTile = ai.getAndIncrement()) {
                        int dl = ((debugLevel > 1) && (nsTile == debug_stile)) ? 3: debugLevel;
						if (dl > 2){
							System.out.println("getPlaneDispStrengths(): nsTile="+nsTile);
						}
						int stileY = nsTile / stilesX;
						int stileX = nsTile % stilesX;
						int [] sTiles = {stileX, stileY};
						// first make a plane from all tiles
						TilePlanes.PlaneData pd0 = tpl.new  PlaneData (
								sTiles, // int [] sTileXY,
								tileSize, // int tileSize,
								geometryCorrection, // GeometryCorrection   geometryCorrection,
								correct_distortions,
								measuredLayers,     // MeasuredLayers measuredLayers,
								plPreferDisparity);   // boolean preferDisparity)

//						if (smplMode && (world_plane_norm != null)) {
						if (smplMode && !floating) {
							plane_tilts = new double [measuredLayers.getNumLayers()][][];
							for (int ml = 0; ml < plane_tilts.length; ml++) if ((stMeasSel & ( 1 << ml)) != 0){
								plane_tilts[ml] = 	zero_tilts;
								if (!const_disparity) {
									double [][] tile_disp_strengths = measuredLayers.getDisparityStrengthML(
											ml,             // int num_layer,
											stileX,         // int stX,
											stileY,         // int stY,
											null,           // boolean [] sel_in,
											mlfp.strength_floor, // double strength_floor,
											mlfp.strength_pow,   // double strength_pow,
											true);          // boolean null_if_none);
									// if failed - keep constant disparity (plane_tilts[ml] = 	zero_tilts)

									if (tile_disp_strengths != null){
										plane_tilts[ml] = pd0.getDisparityTilts(
												world_plane_norm,     // double []     world_normal_xyz,
												tile_disp_strengths,  // double [][]   tile_disp_strengths,
												debugLevel);         // int           debugLevel);
										if(dl > 2){
											if (plane_tilts[ml] != null){
												for (int k = 0; k < 2; k++) {
													System.out.println((k > 0) ? "tilt Y, pix/tile" : "tilt X, pix/tile");
													for (int i = 0; i < 2 * superTileSize; i++){
														System.out.print(i+",");
														for (int j = 0; j < 2 * superTileSize; j++){
															int it = j + 2 * superTileSize * i;
															if (plane_tilts[ml][it] !=null){
																System.out.print(String.format("%f7.4", plane_tilts[ml][it][k]));
															}
															if (j < (2 * superTileSize - 1)){
																System.out.print(", ");
															}
														}
														System.out.println();
													}
													System.out.println();
												}
											}
										}
									}
								}
							}
						}

						plane_disp_strength[nsTile] = new double[measuredLayers.getNumLayers()][][];

						for (int ml = 0; ml < plane_disp_strength[nsTile].length; ml++) if ((stMeasSel & ( 1 << ml)) != 0){
							// TODO": apply tilt before/with getDisparityStrength()
							if (smplMode) {
								plane_disp_strength[nsTile][ml] =  measuredLayers.getDisparityStrengthMLTilted(
										ml,             // int num_layer,
										stileX,         // int stX,
										stileY,         // int stY,
										null,           // boolean [] sel_in,
										(floating? null: plane_tilts[ml]),      // double []  tiltXY, // null - free with limit on both absolute (2.0?) and relative (0.2) values
										mlfp,
										true,           // boolean null_if_none);
										dl);
							} else {
								plane_disp_strength[nsTile][ml] =  measuredLayers.getDisparityStrengthML(
										ml,             // int num_layer,
										stileX,         // int stX,
										stileY,         // int stY,
										null,           // boolean [] sel_in,
										mlfp.strength_floor, // double strength_floor,
										mlfp.strength_pow,   // double strength_pow,
										true);          // boolean null_if_none);
							}
						}
						if (world_plane_norm != null) {
							// find average disparity for the supertile (improve?)
							double sd = 0.0, sw = 0.0;
							for (int ml = 0; ml < plane_disp_strength[nsTile].length; ml++) if ((stMeasSel & ( 1 << ml)) != 0){
								if ((plane_disp_strength[nsTile] !=null) && (plane_disp_strength[nsTile][ml] !=null)  && (plane_disp_strength[nsTile][ml][1] !=null)) {
									for (int i = 0; i < plane_disp_strength[nsTile][ml][1].length; i++){
										double w = plane_disp_strength[nsTile][ml][1][i];
										double d = plane_disp_strength[nsTile][ml][0][i];
										sd += w * d;
										sw += w;
									}
								}
							}
							if (sw > 0) {
								if (dl > 0) {
									System.out.println("Plane tilted disparity for stileX = "+stileX+" stileY="+stileY+", average disparity "+(sd/sw));
								}
								if (dl>2) {
									String [] dbg_titles = showSupertileSeparationTitles( plane_disp_strength[nsTile], null);
									double [][] dbg_img = showSupertileSeparation(false, plane_disp_strength[nsTile], null); // plane_sels);
									ShowDoubleFloatArrays sdfa_instance = new ShowDoubleFloatArrays();
									sdfa_instance.showArrays(dbg_img, 2 * superTileSize, 2* superTileSize, true, "plane_separation_original_disp"+nsTile,dbg_titles);
									dbg_img = showSupertileSeparation(true, plane_disp_strength[nsTile], null); // plane_sels);
									sdfa_instance.showArrays(dbg_img, 2 * superTileSize, 2* superTileSize, true, "plane_separation_original_world"+nsTile,dbg_titles);
								}
								plane_disp_strength[nsTile] = pd0.getDisparityToPlane(
										world_plane_norm, // double []     world_normal_xyz,
										sd / sw, // average disparity // double        disp_center,
										null, // boolean [][]  tile_sel, // null - do not use, {} use all (will be modified)
										plane_disp_strength[nsTile], // double [][][] disp_str, // calculate just once if null
										dl); // 1); // int           debugLevel);
								if (dl>2) {
									String [] dbg_titles = showSupertileSeparationTitles( plane_disp_strength[nsTile], null);
									double [][] dbg_img = showSupertileSeparation(false, plane_disp_strength[nsTile], null); // plane_sels);
									ShowDoubleFloatArrays sdfa_instance = new ShowDoubleFloatArrays();
									sdfa_instance.showArrays(dbg_img, 2 * superTileSize, 2* superTileSize, true, "plane_separation_disp"+nsTile,dbg_titles);
									dbg_img = showSupertileSeparation(true, plane_disp_strength[nsTile], null); // plane_sels);
									sdfa_instance.showArrays(dbg_img, 2 * superTileSize, 2* superTileSize, true, "plane_separation_world"+nsTile,dbg_titles);
								}
							} else {
								plane_disp_strength[nsTile] = null;
							}
						}
					}
				}
			};
		}
		ImageDtt.startAndJoin(threads);
		return plane_disp_strength;
	}

	public boolean [][][][] dispClusterizeHighest(
			final double [][][][] disparity_strengths, // either normal or tilted disparity/strengths
			final boolean [][][] selected,   // tiles OK to be assigned [supertile][measurement layer] [tile index] or null (or null or per-measurement layer)
			final boolean [][][] prohibited, // already assigned tiles [supertile][measurement layer] [tile index] or null
			final boolean        search_min,
			final int            stMeasSel,  //            = 1;      // Select measurements for supertiles : +1 - combo, +2 - quad +4 - hor +8 - vert
			final double         plDispNorm,         // to increase weight of nearer planes
			final double         sigma,
			final double         disp_arange,
			final double         disp_rrange,
			final double         tolerance_above,
			final double         tolerance_below,
			final int            plMinPoints, //          =     5;  // Minimal number of points for plane detection
			final String         suffix,
			final int            debugLevel,
			final int            dbg_X,
			final int            dbg_Y)
	{
		final int tilesX =        tileProcessor.getTilesX();
		final int tilesY =        tileProcessor.getTilesY();
		final int superTileSize = tileProcessor.getSuperTileSize();
		final int stSize2 =       2 * superTileSize;
		final int stLen2 =        stSize2 * stSize2;


		final int stilesX = (tilesX + superTileSize -1)/superTileSize;
		final int stilesY = (tilesY + superTileSize -1)/superTileSize;
		final int nStiles = stilesX * stilesY;
		final Thread[] threads = ImageDtt.newThreadArray((debugLevel > 1)? 1 :tileProcessor.threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		this.planes = new TilePlanes.PlaneData[nStiles][];
		final int debug_stile = (debugLevel > -1)? (dbg_Y * stilesX + dbg_X):-1;
//		final boolean [][][] enabled = (selected == null) ? (new boolean [nStiles][][]) : selected;

		final boolean [][][][] plane_selections = new boolean [nStiles][][][]; // [supertile] [ plane number] [measurement layer] [tile index]
		class Selections2{
			boolean [][] sel;
			Selections2(boolean [][] sel) {
				this.sel = sel;
			}
			boolean [][] getSel(){
				return this.sel;
			}
		}
//		final double max_diff2 = Double.isNaN(max_diff)? Double.NaN: (max_diff*max_diff);
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				@Override
				public void run() {
					for (int nsTile = ai.getAndIncrement(); nsTile < nStiles; nsTile = ai.getAndIncrement()) {
						if (disparity_strengths[nsTile] != null){
                            int dl = ((debugLevel > 0) && (nsTile == debug_stile)) ? (debugLevel + 2) : debugLevel;
							if (dl > 2){
								System.out.println("dispClusterizeHighest(): nsTile="+nsTile);
							}
//							int stileY = nsTile / stilesX;
//							int stileX = nsTile % stilesX;
							double[][][] disp_strength = new double[measuredLayers.getNumLayers()][][];
							for (int ml = 0; ml < disp_strength.length; ml++) if ((stMeasSel & ( 1 << ml)) != 0){
								disp_strength[ml] = disparity_strengths[nsTile][ml]; // will we change it  - no, no need to clone
							}
							boolean [][] enabled = (selected == null) ? (new boolean [disp_strength.length][]) : selected[nsTile];
							for (int ml = 0; ml < disp_strength.length; ml++) if (disp_strength[ml] != null) {
								if (selected != null) {
									enabled[ml] = selected[nsTile][ml].clone();
								}else {
									enabled[ml] = new boolean [stLen2];
									for (int i = 0; i < stLen2; i++) {
										enabled[ml][i] = true;
									}
								}
								if ((prohibited != null) && (prohibited[nsTile] != null) && (prohibited[nsTile][ml] != null)) {
									for (int i = 0; i < stLen2; i++) {
										enabled[ml][i] &= !prohibited[nsTile][ml][i];
									}
								}
							}
							ArrayList<Selections2> sel_list = new ArrayList<Selections2>();
							for (int ntry = 0; ntry < 20; ntry++) {
								boolean [][] new_sel = getHighestPlaneSelection(
										enabled,         //boolean  [][]  enabled, // will not be modified
										disp_strength,   // double [][][]  disp_str,
										//FIXME: use this argument
										search_min, // boolean        search_min,
										2,     // int            grow,
										sigma,           // double         sigma,
										disp_arange,     // double         disp_arange,
										disp_rrange,     // double         disp_rrange,
										tolerance_above, // double         tolerance_above,
										tolerance_below, // double         tolerance_below,
										dl);             // int debugLevel)
								int num_tiles = 0;
								for (int ml = 0; ml < enabled.length; ml++) if (new_sel[ml] != null) {
									for (int i = 0; i < new_sel[ml].length; i++) {
										if (new_sel[ml][i]) {
											num_tiles++;
											enabled[ml][i]=false;
										}
									}
								}
								if (num_tiles < plMinPoints) {
									break;
								}
								sel_list.add(new Selections2(new_sel));
							}
							if (sel_list.isEmpty()) {
								if (dl > 2){
									System.out.println("dispClusterizeHighest(): nsTile="+nsTile+": nothing found");
								}
								continue;
							}
							plane_selections[nsTile] = new boolean [sel_list.size()][][]; //plane_sels;
							for (int np = 0; np < plane_selections[nsTile].length; np++) {
								plane_selections[nsTile][np] = sel_list.get(np).getSel();
							}
							if (dl > 2){
								System.out.println("dispClusterizeHighest(): nsTile="+nsTile);
								double [][]  dbg_img = new double [plane_selections[nsTile].length + 2][];
								String []  dbg_titles = new String [plane_selections[nsTile].length + 2];
								for (int ml = 0; ml < disp_strength.length; ml++) if (disp_strength[ml] != null){
									dbg_img[0] =                disp_strength[ml][0];
									dbg_img[dbg_img.length-1] = disp_strength[ml][1];
									break;
								}
								dbg_titles[0] =                "disparity";
								dbg_titles[dbg_img.length-1] = "strength";
								for (int np = 0; np < plane_selections[nsTile].length; np++) {
									double plane_strength = getSelectionStrength(
											true, // mod_strength, // boolean mod_strength,
											plDispNorm,          // double plDispNorm,
											plane_selections[nsTile][np], // boolean [][] sels,
											disp_strength); // double [][][] ds);
									dbg_titles[1 + np] = ""+plane_strength;
									dbg_img[1 + np] = new double [stLen2];
									for (int i = 0; i < stLen2; i++) {
										boolean has_tile = false;
										for (int ml = 0; ml < plane_selections[nsTile][np].length; ml++) if (plane_selections[nsTile][np][ml] != null){
											if (plane_selections[nsTile][np][ml][i]) {
												dbg_img[1 + np][i] +=  disp_strength[ml][0][i];
												has_tile = true;
											}
											if (!has_tile) {
												dbg_img[1 + np][i] = Double.NaN;
											}
										}
									}
								}
								(new ShowDoubleFloatArrays()).showArrays(dbg_img, stSize2, stSize2, true,
										"sel-"+nsTile+"-"+suffix,dbg_titles);
							}

						}
					}
				}
			};
		}
		ImageDtt.startAndJoin(threads);
		return plane_selections;
	}


	public boolean [][] getHighestPlaneSelection(
			boolean  [][]  enabled_ml, // will not be modified
			double [][][]  disp_str,
			boolean        search_min,
			int            grow,
			double         sigma,
			double         disp_arange,
			double         disp_rrange,
			double         tolerance_above,
			double         tolerance_below,
			int debugLevel){
		final int min_en_neib = 4;
		final boolean [][] plane_selection = new boolean[disp_str.length][];
		final double min_weight = 1.0E-20;
		final int stSize2 = 2 * tileProcessor.getSuperTileSize();
		final int stLen2 = stSize2 * stSize2;
		DoubleGaussianBlur gb=new DoubleGaussianBlur();
		TileNeibs tileNeibs = new TileNeibs(stSize2, stSize2);
		int expand = (grow + 1) & (~1);
		for (int ml = 0; ml < disp_str.length; ml++) if (disp_str[ml] != null) {
			double [] disparity =   disp_str[ml][0].clone();
			double [] strength =    disp_str[ml][1].clone();
			boolean [] enabled = enabled_ml[ml].clone();
			for (int i = 0; i < stLen2; i++) {
				if ((strength [i] > 0.0) && enabled[i]) {
					if (min_en_neib > 0) {
						int nn = 0;
						for (int dir = 0; dir < 8; dir++) {
							int i1 = tileNeibs.getNeibIndex(i, dir);
							if ((i1 >=0) && enabled_ml[ml][i1]) {
								nn++;
							}
						}
						if (nn < min_en_neib) {
							enabled[i] = false;
							strength[i] = 0.0;
							disparity[i] = 0.0;
						}
					}
					disparity[i] *= strength[i];
				} else {
					strength[i] = 0.0;
					disparity[i] = 0.0;
				}
			}
			gb.blurDouble(disparity,   stSize2, stSize2, sigma, sigma, 0.01);
			gb.blurDouble(strength,    stSize2, stSize2, sigma, sigma, 0.01);
			// TODO: handle large gaps with two sigmas?
			for (int i = 0; i < stLen2; i++) {
				if (strength [i] > min_weight) {
					disparity[i] /= strength [i];
				} else {
					strength [i] = 0.0;
					disparity[i] = 0.0;
				}
			}
			if (debugLevel > 2){
				System.out.println("getHighestPlaneSelection() debugLevel="+debugLevel);
			}
			boolean [] best_sel = null;
			double best_weight = 0.0;
			int best_height = 0;
//			int best_bottom;
			int bottom = stSize2 - 1;
			for (; bottom >= (best_height - 1); bottom --) {
				double line_max_max = 0.0, line_min_max= 0.0;
				int [] tops = new int [stSize2];
				int [][] line_se = new int [stSize2][2];
				int top;
				double disp_range = disp_arange;
				for (top = bottom; top >= 0; top--) {
					// find argmax for the top line;
					tops[top] = -1; //top *  stSize2;
					for (int ix = 0; ix < stSize2; ix++) {
						int indx = top *  stSize2 + ix;
						if (search_min) {
							if ((strength[indx] > 0.0) && ((tops[top] < 0) || (disparity[indx] < disparity[tops[top]]))){
								tops[top] = indx;
							}
						} else {
							if ((strength[indx] > 0.0) && ((tops[top] < 0) || (disparity[indx] > disparity[tops[top]]))){
								tops[top] = indx;
							}
						}
					}
					if (tops[top] > 0) {
						if (line_max_max == 0.0) {
							line_max_max = disparity[tops[top]]; // for search_min it will be line_max_min
							line_min_max = disparity[tops[top]]; // for search_min it will be line_min_min
							disp_range = disp_arange + 0.5 *(line_max_max + line_min_max) * disp_rrange;
						} else {
							if (disparity[tops[top]] > line_max_max) {
								if (disparity[tops[top]] > (line_min_max + disp_range)) {
									break; // difference more than allowed range
								}
								line_max_max = disparity[tops[top]];
								disp_range = disp_arange + 0.5 *(line_max_max + line_min_max) * disp_rrange;
							} else if (disparity[tops[top]] < line_min_max) {
								if (disparity[tops[top]] < (line_max_max - disp_range)) {
									break; // difference more than allowed range
								}
								line_min_max = disparity[tops[top]];
								disp_range = disp_arange + 0.5 *(line_max_max + line_min_max) * disp_rrange;
							}
						}
					} else { // no valid cells in a line (too far from enabled)
						break;
					}
				}
				if ((bottom - top) < (best_height - 2 - expand)) { //
					continue; // no sense even to try
				}
				disp_range = disp_arange + 0.5 *(line_max_max + line_min_max) * disp_rrange;
				double min_disp = line_max_max - disp_range; // used when looking for maximums
				double max_disp = line_min_max + disp_range; // used when looking for minimums
				// for each line set start and end, verify lines overlap
				int line = bottom;
				for (; line > top; line --) {
					// do not check first in the row
					for (line_se[line][0] = tops[line] ; line_se[line][0] >  stSize2 * line; line_se[line][0]--) {
						if (search_min) {
							if ((strength[line_se[line][0]] > 0.0) && (disparity[line_se[line][0]] > max_disp)) {
								//line_se[line][0]++; allow one extra
								break;
							}
						} else {
							if ((strength[line_se[line][0]] > 0.0) && (disparity[line_se[line][0]] < min_disp)) {
								//line_se[line][0]++; allow one extra
								break;
							}
						}
					}
					// do not check last in the row
					for (line_se[line][1] = tops[line] ; line_se[line][1] <  (stSize2 * (line + 1) -1); line_se[line][1]++) {
						if (search_min) {
							if ((strength[line_se[line][1]] > 0.0) && (disparity[line_se[line][1]] > max_disp)) {
								//line_se[line][1]--;
								break;
							}
						} else {
							if ((strength[line_se[line][1]] > 0.0) && (disparity[line_se[line][1]] < min_disp)) {
								//line_se[line][1]--;
								break;
							}
						}
					}
					// verify the lines overlap with previous
					if ((line < bottom) && (
							(line_se[line][0] >= (line_se[line + 1][1] - stSize2)) ||
							(line_se[line][1] <= (line_se[line + 1][0] - stSize2)))) {
						break;
					}
				}

				// now line - first invalid line
				// calculate actual selection and weight
				// duplicate first and last line to mitigate Gaussian blur - there may be some FG tiles cut off
				int bottom_ext = bottom;
				if (bottom_ext < (stSize2 - 1)) {
					bottom_ext++;
					line_se[bottom_ext][0] = line_se[bottom][0] + stSize2;
					line_se[bottom_ext][1] = line_se[bottom][1] + stSize2;
				}
				if ((line >= 0) && (line < (stSize2 - 1))) {
					line_se[line][0] = line_se[line + 1][0] - stSize2 ;
					line_se[line][1] = line_se[line + 1][1] - stSize2 ;
					line--;
				}

				if ((bottom_ext - line) < (best_height - expand)) {
					continue; // no sense even to try
				}

				boolean [] this_sel = new boolean [stLen2];
				double sw = 0.0;
				double min_tile_disp = search_min ? (line_min_max - tolerance_below): (min_disp - tolerance_below);
				double max_tile_disp = search_min ? (max_disp + tolerance_above) :    (line_max_max + tolerance_above);
				int first = -1, last = -1;

				for (int i = bottom_ext; i > line; i--) {
					for (int j = line_se[i][0]; j <= line_se[i][1]; j++) {
						if (  enabled[j] &&
								(disp_str[ml][1][j] > 0.0) &&
								(disp_str[ml][0][j] >= min_tile_disp) &&
								(disp_str[ml][0][j] <= max_tile_disp)) {
							sw += disp_str[ml][1][j];
							this_sel[j] = true;
							if (last < 0) last = i;
							first = i;
						}
					}
				}
				int sel_height = last - first + 1;
				if (sel_height < (best_height - expand)) {
					continue;
				}
				if (expand > 0) {
					boolean [] prohibit = new boolean [enabled.length];
					for (int i = 0; i < prohibit.length; i++) {
						prohibit[i] = !enabled[i];
					}
					boolean [] exp_sel = this_sel.clone();
					tileNeibs.growSelection(
							grow,
							exp_sel,
							prohibit);
					int i0 = -1, i1 = -1;
					for (int i = 0; i < stLen2; i++) if (exp_sel[i]){
						if (!this_sel[i]) {
							if ((disp_str[ml][1][i] > 0.0) &&
							(disp_str[ml][0][i] >= min_tile_disp) &&
							(disp_str[ml][0][i] <= max_tile_disp)) {
								this_sel[i] = true;
							}
						}
						if (this_sel[i]) {
							sw += disp_str[ml][1][i];
							if (i0 < 0) i0 = i;
							i1 = i;
						}
					}
					first = i0 / stSize2;
					last =  i1 / stSize2;
					sel_height = last - first + 1;
					if (sel_height < (best_height - expand)) {
						continue;
					}
				}
				if ((sel_height > best_height) ||
						(sw > best_weight)) {
					best_sel = this_sel;
					best_weight = sw;
					best_height = sel_height;
				}
			} // bottom
			plane_selection[ml] = best_sel;
			// add debug image here
			if ((debugLevel > 2) && (best_sel != null)) {
				double [] disp_sel = disp_str[ml][0].clone();
				for (int i = 0; i < plane_selection[ml].length; i++) {
					if (!plane_selection[ml][i]) {
						disp_sel[i] = Double.NaN;
					}
				}
				String [] dbg_titles = {"disp_blur","disp_sel","disp","str_blur","str"};
				double [][] dbg_img =
					{   disparity,
						disp_sel,
						disp_str[ml][0],
						strength,
						disp_str[ml][1]};
				(new ShowDoubleFloatArrays()).showArrays(dbg_img, stSize2, stSize2, true, "disp_blured",dbg_titles);
				System.out.println("getHighestPlaneSelection() done");
			}

		} // ml
		return plane_selection;
	}




	// use histogram data (min/max/strength) assign tiles to clusters
	// returns [supertile] [ plane number] [measurement layer] [tile index]
	public boolean [][][][] dispClusterize(
			final double [][][][] disparity_strengths, // either normal or tilted disparity/strengths
			final double [][][] hist_max_min_max, // histogram data: per tile array of odd number of disparity/strengths pairs, starting with first maximum
			final boolean [][][] selected,   // tiles OK to be assigned [supertile][measurement layer] [tile index] or null (or null or per-measurement layer)
			final boolean [][][] prohibited, // already assigned tiles [supertile][measurement layer] [tile index] or null
			final int        stMeasSel, //            = 1;      // Select measurements for supertiles : +1 - combo, +2 - quad +4 - hor +8 - vert
			// TODO: Scale max_diff, smallDiff for large disparities
			final double     max_diff,  // maximal disparity difference (to assign to a cluster (of Double.NaN)
			final int        plMinPoints, //          =     5;  // Minimal number of points for plane detection
			final double     smallDiff,  //       = 0.4;   // Consider merging initial planes if disparity difference below
			final double     highMix,    //stHighMix         = 0.4;   // Consider merging initial planes if jumps between ratio above
			final boolean    norm_max_diff,   // scale  max_diff for large (> dispNorm) average disparities
			final boolean    norm_small_diff, // scale  max_diff for large (> dispNorm) average disparities
			final double     dispNorm,
			final int        debugLevel,
			final int        dbg_X,
			final int        dbg_Y)
	{
		final int tilesX =        tileProcessor.getTilesX();
		final int tilesY =        tileProcessor.getTilesY();
		final int superTileSize = tileProcessor.getSuperTileSize();

		final int stilesX = (tilesX + superTileSize -1)/superTileSize;
		final int stilesY = (tilesY + superTileSize -1)/superTileSize;
		final int nStiles = stilesX * stilesY;
//		final Thread[] threads = ImageDtt.newThreadArray(tileProcessor.threadsMax);
		final Thread[] threads = ImageDtt.newThreadArray((debugLevel > 1)? 1 :tileProcessor.threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		this.planes = new TilePlanes.PlaneData[nStiles][];
		final int debug_stile = (debugLevel > -1)? (dbg_Y * stilesX + dbg_X):-1;

		final boolean [][][][] plane_selections = new boolean [nStiles][][][]; // [supertile] [ plane number] [measurement layer] [tile index]
		final double max_diff2 = Double.isNaN(max_diff)? Double.NaN: (max_diff*max_diff);
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				@Override
				public void run() {
					for (int nsTile = ai.getAndIncrement(); nsTile < nStiles; nsTile = ai.getAndIncrement()) {
						if (disparity_strengths[nsTile] != null){
                            int dl = ((debugLevel > 1) && (nsTile == debug_stile)) ? 3: debugLevel;
							if (dl > 1){
								System.out.println("dispClusterize(): nsTile="+nsTile);
							}
							int stileY = nsTile / stilesX;
							int stileX = nsTile % stilesX;
//							int dl =  (nsTile == debug_stile) ? 3 : 0;

							double[][][] disp_strength = new double[measuredLayers.getNumLayers()][][];

							for (int ml = 0; ml < disp_strength.length; ml++) if ((stMeasSel & ( 1 << ml)) != 0){
								disp_strength[ml] = disparity_strengths[nsTile][ml]; // will we change it  - no, no need to clone
							}

							if (hist_max_min_max[nsTile] == null){
								continue;
							}
							double [][] max_only = new double [(hist_max_min_max[nsTile].length + 1)/2][2];
							for (int i = 0; i < max_only.length; i++){
								max_only[i] = hist_max_min_max[nsTile][2 * i];
							}
							boolean [][][] plane_sels = null;
							int num_ml = disp_strength.length;
							int num_tiles = 4 * superTileSize * superTileSize;
							int [] num_sel;
							boolean [][] tile_en = new boolean [num_ml][];
							for (int ml = 0; ml < num_ml; ml++) if (
									(disp_strength[ml] != null) &&
									((stMeasSel & ( 1 << ml)) != 0) &&
									((selected == null) || ((selected[nsTile] != null) && (selected[nsTile][ml] != null)))
									) {
								tile_en[ml] = new boolean [num_tiles];
								// apply all restrictions here
								for (int i = 0; i < num_tiles; i++){
									if (    (disp_strength[ml][1][i] > 0.0) &&
											((selected == null) || (selected[nsTile][ml].length == 0) || selected[nsTile][ml][i]) &&
											((prohibited == null) || (prohibited[nsTile] == null) || (prohibited[nsTile][ml] == null) || !prohibited[nsTile][ml][i])
											) {
										tile_en[ml][i] = true;
									}
								}
							}

							for (int iter = 0; iter < 2; iter ++){
								int num_p  = max_only.length;
								plane_sels = new boolean[num_p][num_ml][];
								num_sel = new int [num_p];
								for (int np = 0; np < num_p; np++) {
									for (int ml = 0; ml < num_ml; ml++) if (
											(disp_strength[ml] != null) &&
											((stMeasSel & ( 1 << ml)) != 0) &&
											((selected == null) || ((selected[nsTile] != null) && (selected[nsTile][ml] != null)))
											) {
										plane_sels[np][ml] = new boolean[num_tiles];
									}
								}
								// compare closest to be able to use tilted planes later
								for (int ml = 0; ml < num_ml; ml++) if (tile_en[ml] != null) {
									if ( dl >2) {
										System.out.println("dispClusterize(), nsTile="+nsTile+", tile_en["+ml+"] + iter = "+iter);
										for (int i = 0; i < 2 * superTileSize; i ++){
											for (int j = 0; j < 2 * superTileSize; j ++){
												int indx = 2 * superTileSize * i + j;
												System.out.print(tile_en[ml][indx]?" +":" .");
											}
											System.out.println();
										}
									}
									for (int indx = 0; indx < num_tiles; indx++) if (tile_en[ml][indx]){
										int best_plane = -1;
										double best_d2 = Double.NaN;
										for (int np = 0; np < num_p; np++) {
											double d2 = max_only[np][0] - disp_strength[ml][0][indx];
											// add disp_norm correction here? Yes!
											if (norm_max_diff) {
												double d_avg = 0.5 * (max_only[np][0] + disp_strength[ml][0][indx]);
												if (d_avg > dispNorm){
													d2 *= dispNorm / d_avg;
												}
											}
											d2 *= d2;
											if (!(d2 >= best_d2)){
												best_d2 = d2;
												best_plane = np;
											}
										}
										if (best_plane >= 0){ // compare to max diff here too
											if  (!(best_d2 > max_diff2)) { // works if max_diff2 is Double.NaN
												plane_sels[best_plane][ml][indx] = true; // so far exclusive
											}
										}
									}
									if ( dl > 2) {
										System.out.println("dispClusterize(), nsTile="+nsTile+", plane_sels[?]["+ml+"] ,iter = "+iter);
										for (int i = 0; i < 2 * superTileSize; i ++){
											for (int j = 0; j < 2 * superTileSize; j ++){
												int indx = 2 * superTileSize * i + j;
												label1: {
													for (int np = 0; np < num_p; np++) {
														if (plane_sels[np][ml][indx]) {
															System.out.print(" " +np);
															break label1;
														}
													}
													System.out.print(" .");
												}
											}
											System.out.println();
										}
									}
								} //for (int ml = 0; ml < num_ml; ml++) if (tile_en[ml] != null) {
								// recalculate average disparities for each plane and show number of tiles in each in debug mode
								for (int np = 0; np < num_p; np++) {
									double sd = 0.0, sw = 0.0;
									for (int ml = 0; ml < num_ml; ml++) if (disp_strength[ml] != null) {
										for (int indx = 0; indx < num_tiles; indx++) if (plane_sels[np][ml][indx]){
											double w = disp_strength[ml][1][indx];
											sd += w * disp_strength[ml][0][indx];
											sw += w;
											num_sel[np]++;
										}
									}
									if (sw > 0) {
										sd /= sw;
									}
									if (dl > 2) {
										System.out.println("plane num_sel["+np+"] = "+num_sel[np]+" disp "+max_only[np][0]+"->"+sd+
												", weight "+max_only[np][1]+"->"+sw);
									}
									max_only[np][0] = sd;
									max_only[np][1] = sw;
								}
								// calculate transitions matrix (to find candidates for merge
								int [][]    trans_mat =  getTransMatrix(plane_sels);
								double [][] rel_trans =  getTransRel(trans_mat);


								if (dl > 1) {
									System.out.println("trans_mat = ");
									for (int i = 0; i < trans_mat.length; i++){
										System.out.print(i+": ");
										for (int j = 0; j < trans_mat[i].length; j++){
											System.out.print(trans_mat[i][j]+" ");
										}
										System.out.println();
									}
									System.out.println("rel_trans = ");
									for (int i = 0; i < rel_trans.length; i++){
										System.out.print(i+": ");
										for (int j = 0; j < rel_trans[i].length; j++){
											System.out.print(rel_trans[i][j]+" ");
										}
										System.out.println();
									}
								}
								if ((iter > 0 ) && (num_p > 1)){ // remove /join bad
									int windx = 0;
									int remove_indx = -1;
									for (int i = 1; i < num_p; i++)	if (num_sel[i] < num_sel[windx]) windx = i;
									if (num_sel[windx] < plMinPoints) {
										if ((debugLevel > 2) || (dl > 2)){
											System.out.println ("dispClusterize(): stileX = "+stileX+" stileY="+stileY+
													": removing plane "+windx+" with "+num_sel[windx]+" tiles ( <"+plMinPoints+")");
										}
										remove_indx = windx;
									}
									if (remove_indx < 0) {
										// find candidates for merge
										windx = -1;
										for (int i = 0; i < (num_p - 1); i++)	{
											double diff_disp = max_only[i+1][0] - max_only[i][0];
											if (norm_small_diff) {
												double d_avg = 0.5 * (max_only[i+1][0] + max_only[i][0]);
												if (d_avg > dispNorm){
													diff_disp *= dispNorm / d_avg;
												}
											}

											if ((diff_disp < smallDiff) &&  // close enough to consider merging
													(rel_trans[i][i+1] > highMix)) {
												if ((windx < 0) || (rel_trans[i][i+1] > rel_trans[windx][windx+1])) windx = i;
											}
										}
										if (windx >=0 ) {
											if ((debugLevel > 2) || (dl > 2)) {
												System.out.println ("dispClusterize(): stileX = "+stileX+" stileY="+stileY+
														": merging plane "+windx+" with " + (windx + 1)+": "+
														num_sel[windx] + " and "+num_sel[windx+1]+" tiles, "+
														" rel_trans="+rel_trans[windx][windx + 1]+ " ( > " + highMix+"),"+
														" diff="+ (max_only[windx + 1][0]- max_only[windx][0]) + " ( < " + smallDiff+" ),"+
														" disp1 = "+max_only[windx][0]+" disp2 = "+max_only[windx + 1][0]);
											}
											double sum_w = max_only[windx][1] + max_only[windx + 1][1];
											max_only[windx+1][0] = (max_only[windx][0]*max_only[windx][1] +  max_only[windx+1][0]*max_only[windx+1][1]) / sum_w;
											max_only[windx+1][1] = sum_w;
											remove_indx = windx;
										}
									}
									if (remove_indx >= 0){
										double [][] max_only_copy = max_only.clone();
										for (int i = 0; i < max_only.length; i++) max_only_copy[i] = max_only[i];
										max_only = new double [max_only.length - 1][];
										int indx = 0;
										for (int i = 0; i < max_only_copy.length; i++) if (i != remove_indx) max_only[indx++] =max_only_copy[i];
										iter = 0;
										continue; // restart from 0

									}

									// Show other candidates for merge
									if ((debugLevel > 2) || (dl > 0)){
										double max_sep = 0.2;
										if (iter  > 0) {
											for (int i = 0; i < (num_p-1); i++){
												if (rel_trans[i][i+1] > max_sep) {
													System.out.println("dispClusterize() stileX = "+stileX+" stileY="+stileY+" lowplane = "+i+
															" num_sel1 = "+num_sel[i] + " num_sel2 = "+num_sel[i+1] +
															" rel_trans="+rel_trans[i][i+1]+
															" diff="+ (max_only[i+1][0]- max_only[i][0]) +
															" disp1 = "+max_only[i][0]+" disp2 = "+max_only[i+1][0]);
												}
											}
										}
									}
								}
							} // for (int iter = 0; iter < 2; iter ++){

							if (dl > 3) {
								String [] dbg_titles = showSupertileSeparationTitles( disp_strength, plane_sels);
								double [][] dbg_img = showSupertileSeparation(false,disp_strength, plane_sels);
								ShowDoubleFloatArrays sdfa_instance = new ShowDoubleFloatArrays();
								sdfa_instance.showArrays(dbg_img, 2 * superTileSize, 2* superTileSize, true, "disp_clusterize_disp"+nsTile+"-"+debugLevel,dbg_titles);
								dbg_img = showSupertileSeparation(true, disp_strength, plane_sels);
								sdfa_instance.showArrays(dbg_img, 2 * superTileSize, 2* superTileSize, true, "disp_clusterize_world"+nsTile+"-"+debugLevel,dbg_titles);
							}
							plane_selections[nsTile] = plane_sels;
						}
					}
				}
			};
		}
		ImageDtt.startAndJoin(threads);
		return plane_selections;
	}

	private double getSelectionStrength(
			boolean mod_strength,
			double plDispNorm,
			boolean [][] sels,
			double [][][] ds)
	{
		int num_layers = ds.length;
		int superTileSize = tileProcessor.getSuperTileSize();
		if (sels == null) return 0.0;
		double sw = 0.0;
		int num_tiles = 0;
		if (mod_strength) {
			TileNeibs tileNeibs = new TileNeibs(2 * superTileSize, 2 * superTileSize);
			for (int ml = 0; ml < num_layers; ml++) if (sels[ml] != null) {
				for (int i = 0; i < sels[ml].length; i++) {
					if (sels[ml][i]) {
						int nn = 0;
						for (int dir = 0; dir < 8; dir++) {
							int i1 = tileNeibs.getNeibIndex(i, dir);
							if ((i1 >=0) && sels[ml][i1]) {
								nn++;
							}
						}
						double w = ds[ml][1][i];
						double d = ds[ml][0][i];
						if ((plDispNorm > 0.0) && (d > 0.0) && (d < plDispNorm)) {
							w *=  d / plDispNorm;
						}
						sw +=   w * nn;
					}
				}
			}
			// increasing weight of strong compact clusters - see if it improves results
			if (num_tiles > 0) {
				sw /= Math.sqrt(num_tiles);
			}
		} else { // old way

			for (int ml = 0; ml < num_layers; ml++) if (sels[ml] != null) {
				for (int i = 0; i < sels[ml].length; i++) {
					if (sels[ml][i]) {
						sw +=  ds[ml][1][i];
					}
				}
			}
		}
		return sw;
	}


	/**
	 * Use both horizontal and const disparity tiles to create tile clusters
	 * Add max_diff (maximal disparity difference while extracting initial tile selection) and max_tries (2..3) parameters
	 * Add separate method to create + remove outliers from all planes (2 different ones)?
	 * TODO later re-assign pixels according to existing plane parameters
	 * Sort plane data by center (plane or supertile) disparity
	 *
	 * @param growSelection
	 * @param stMeasSel Select measurements for supertiles : +1 - combo, +2 - quad +4 - hor +8 - vert
	 * @param plDispNorm
	 * @param plMinPoints
	 * @param plPreferDisparity
	 * @param geometryCorrection
	 * @param correct_distortions
	 * @param smplMode
	 * @param smplSide
	 * @param smplNum
	 * @param smplRms
	 * @param bin_blur_hor
	 * @param bin_blur_vert
	 * @param max_diff_hor
	 * @param max_diff_vert
	 * @param max_tries
	 * @param smallDiff
	 * @param highMix
	 * @param world_hor
	 * @param debugLevel
	 * @param dbg_X
	 * @param dbg_Y
	 * @return
	 */

	public boolean [][][][]  initialDiscriminateTiles(
			final int        growSelection,                     // grow initial selection before processing
			final int        stMeasSel,          //      = 1;      // Select measurements for supertiles : +1 - combo, +2 - quad +4 - hor +8 - vert
			final double     plDispNorm,
			final int        plMinPoints,       //          =     5;  // Minimal number of points for plane detection
			final boolean    plPreferDisparity, // Always start with disparity-most axis (false - lowest eigenvalue)
			final GeometryCorrection geometryCorrection,
			final boolean    correct_distortions,

			final boolean    smplMode,       //        = true;   // Use sample mode (false - regular tile mode)

			final MeasuredLayersFilterParameters mlfp,

			final double     bin_blur_hor,   // Blur disparity histograms for horizontal clusters by this sigma (in bins)
			final double     bin_blur_vert,  // Blur disparity histograms for constant disparity clusters by this sigma (in bins)
			// TODO: scale down max_diff_hor, max_diff_vert  for large disparities?
			final double     max_diff_hor,   // maximal disparity difference (to assign to a cluster (of Double.NaN) at first run for horizontal planes
			final double     max_diff_vert,  // maximal disparity difference (to assign to a cluster (of Double.NaN) at first run for vertical plane

			final int        max_tries,       // on last run - assign all remaining pixels to some cluster (disregard max_diff)
			final double     smallDiff,  //       = 0.4;   // Consider merging initial planes if disparity difference below
			// TODO: scale down smallDiff for large disparities?
			final double     highMix,    //stHighMix         = 0.4;   // Consider merging initial planes if jumps between ratio above
			final double []  world_hor, // horizontal plane normal (default [0.0, 1.0, 0.0])
			final boolean    show_histograms,
			final boolean [][] hor_planes, // returns plane types (hor/vert)

// Parameters for alternative initial planes that use lowest disparity for fronto planes, and farthest - for horizontal
			final boolean    mod_strength,          //  = true; // FIXME: make a parameter. when set, multiply each tile strength by the number of selected neighbors
			final boolean    clusterize_by_highest, //  = true;
			final double     clust_sigma,           //  = 0.7;
			final double     disp_arange_vert,      //  = 0.07;
			final double     disp_rrange_vert,      //  = 0.01;
			final double     disp_arange_hor,       //  =  0.035;
			final double     disp_rrange_hor,       //  =   0.005;
			final double     tolerance_above_near,  //  =  100.0; // 0.07; any?
			final double     tolerance_below_near,  //  =  -0.01;
			final double     tolerance_above_far,   //  =    0.07;
			final double     tolerance_below_far,   //  =    0.1; // 100.0; // any farther
			final int        hor_vert_overlap,      //  =       2;
			final int        used_companions,       //  =      5; // cell that has this many new used companions is considered used (borders and already use3d are considered used too)
			final int        used_true_companions,  //  = 1; // there should be at least this many new selected tiles among neighbors.,

			final int        debugLevel,
			final int        dbg_X,
			final int        dbg_Y)
	{
		// TODO: Make a configurable parameters (0 should also work)
		final int min_fresh_used = 1; // minimal number of new tiles not already used by the other direction (hor/vert)

		final int tilesX =        tileProcessor.getTilesX();
		final int tilesY =        tileProcessor.getTilesY();
		final int superTileSize = tileProcessor.getSuperTileSize();
		final int stilesX = (tilesX + superTileSize -1)/superTileSize;
		final int stilesY = (tilesY + superTileSize -1)/superTileSize;
		final int nStiles = stilesX * stilesY;
		final Thread[] threads = ImageDtt.newThreadArray((debugLevel > 1)? 1 : tileProcessor.threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		this.planes = new TilePlanes.PlaneData[nStiles][];
		final int debug_stile = (debugLevel > -1)? (dbg_Y * stilesX + dbg_X):-1;
		// TODO: Remove when promoting PlaneData
		final int num_layers = measuredLayers.getNumLayers();
		final int num_tiles = 4 * superTileSize * superTileSize;
		boolean [] grown_selection = null;
		if (growSelection >= 0) {
			grown_selection = cltPass3d.getSelected();
			if (growSelection > 0) {
				grown_selection = grown_selection.clone();
				measuredLayers.growSelection(
						growSelection, // grow,
						grown_selection, // tiles,
						null); // prohibit);
			}
		}
		measuredLayers.setLayer ( // overwrite?
				0, // int       num_layer,
				cltPass3d.getDisparity(), // double [] disparity,
				cltPass3d.getStrength(), // double [] strength,
				grown_selection); // null); // boolean [] selection) // may be null
		if (show_histograms && (debugLevel > 0)) {
			String [] titles = {"d0","s0","d1","s1","d2","s2","d3","s3","d","s","selection"};
			boolean [] dbg_sel= grown_selection; // cltPass3d.getSelected();
			double [][] dbg_img = new double [titles.length][];
			for (int i = 0; i < measuredLayers.getNumLayers(); i++){
				dbg_img[2 * i] =     measuredLayers.getDisparity(i);
				dbg_img[2 * i + 1] = measuredLayers.getStrength(i);
			}
			dbg_img[8] = cltPass3d.getDisparity();
			dbg_img[9] = cltPass3d.getStrength();
			dbg_img[10]= new double [dbg_sel.length];
			for (int i = 0; i < dbg_sel.length; i++){
				dbg_img[10][i] = ((dbg_sel == null) || dbg_sel[i])? 1.0:0.0;
			}
			ShowDoubleFloatArrays sdfa_instance = new ShowDoubleFloatArrays();
			sdfa_instance.showArrays(dbg_img,  tileProcessor.getTilesX(), tileProcessor.getTilesY(), true, "measuredLayers",titles);
		}

//		double [] world_hor = {0.0, 1.0, 0.0};
		double [] const_disp = {0.0, 0.0, 1.0}; // constant z in world coordinates, same as constant disparity

		final double [][][][] hor_disp_strength = getPlaneDispStrengthsST(
				world_hor,           // final double []  world_plane_norm, // real world normal vector to a suggested plane family (0,1,0) for horizontal planes
				stMeasSel,           //final int        stMeasSel, //            = 1;      // Select measurements for supertiles : +1 - combo, +2 - quad +4 - hor +8 - vert
				plPreferDisparity,   // final boolean    plPreferDisparity, // Always start with disparity-most axis (false - lowest eigenvalue)
				geometryCorrection,  // final GeometryCorrection geometryCorrection,
				correct_distortions, // final boolean    correct_distortions,
				smplMode,            // final boolean    smplMode, //        = true;   // Use sample mode (false - regular tile mode)
				mlfp,
				debugLevel,
				dbg_X,
				dbg_Y);

		final double [][][][] vert_disp_strength = getPlaneDispStrengthsST(
				const_disp,          // final double []  world_plane_norm, // real world normal vector to a suggested plane family (0,1,0) for horizontal planes
				stMeasSel,           //final int        stMeasSel, //            = 1;      // Select measurements for supertiles : +1 - combo, +2 - quad +4 - hor +8 - vert
				plPreferDisparity,   // final boolean    plPreferDisparity, // Always start with disparity-most axis (false - lowest eigenvalue)
				geometryCorrection,  // final GeometryCorrection geometryCorrection,
				correct_distortions, // final boolean    correct_distortions,
				smplMode,            // final boolean    smplMode, //        = true;   // Use sample mode (false - regular tile mode)
				mlfp,
				debugLevel,
				dbg_X,
				dbg_Y);


		String [] dbg_hist_titles = {"all","hor","mm_vert","mm_hor"};
		double [][] dbg_hist = new double [dbg_hist_titles.length][];
		final boolean [][][] used_hor = new boolean [nStiles][num_layers][]; // num_tiles
		final boolean [][][] used_vert = (hor_vert_overlap == 0) ? used_hor : (new boolean [nStiles][num_layers][]);
		final boolean [][][][] planes_selections = new boolean [nStiles][][][]; // num_tiles
		for (int pass = 0; pass < max_tries; pass ++) {
//			final int fpass = pass;
			// get both horizontal planes and constant disparity planes histograms,
//			resetDisparityHistograms();
			setBlurSigma(bin_blur_hor);
			final double [][][] mmm_hor = getMaxMinMax(
					hor_disp_strength,  // final double [][][][] disparity_strength, // pre-calculated disparity/strength [per super-tile][per-measurement layer][2][tiles] or null
					null); // final boolean [][] tile_sel // null  or per-measurement layer, per-tile selection. For each layer null - do not use, {} - use all
			if (show_histograms && (debugLevel > 0)) {
				dbg_hist[1] = showDisparityHistogram();
				dbg_hist[3] = showMaxMinMax();
			}
//			resetDisparityHistograms();
			setBlurSigma(bin_blur_vert);
			final double [][][] mmm_vert = getMaxMinMax(
					vert_disp_strength,  // final double [][][][] disparity_strength, // pre-calculated disparity/strength [per super-tile][per-measurement layer][2][tiles] or null
					null); // final boolean [][] tile_sel // null  or per-measurement layer, per-tile selection. For each layer null - do not use, {} - use all

			if (show_histograms && (debugLevel > 0)) {
				dbg_hist[0] = showDisparityHistogram();
				dbg_hist[2] = showMaxMinMax();
			}
			if (show_histograms && (debugLevel > 0)) {
				int hist_width0 =  showDisparityHistogramWidth();
				int hist_height0 = dbg_hist[0].length/hist_width0;
				ShowDoubleFloatArrays sdfa_instance = new ShowDoubleFloatArrays(); // just for debugging?
				sdfa_instance.showArrays(dbg_hist, hist_width0, hist_height0, true, "vert_hor_histograms_"+pass,dbg_hist_titles);
			}
			// try to independently (same selections) clusterize both ways
			if (debugLevel > 0){
				System.out.println("initialDiscriminateTiles(): before new_planes_hor, pass =" + (pass + 1) + " ( of "+max_tries+" )");
			}
//			if (clusterize_by_highest) {
			final boolean [][][][] new_planes_hor = clusterize_by_highest ?
					dispClusterizeHighest(
							hor_disp_strength,  // final double [][][][] disparity_strengths, // either normal or tilted disparity/strengths
							null,               // final boolean [][][] selected,   // tiles OK to be assigned [supertile][measurement layer] [tile index] or null (or null or per-measurement layer)
							used_hor,           // final boolean [][][] prohibited, // already assigned tiles [supertile][measurement layer] [tile index] or null
							true,               // final boolean        search_min,
							stMeasSel,          // final int            stMeasSel, //            = 1;      // Select measurements for supertiles : +1 - combo, +2 - quad +4 - hor +8 - vert
							plDispNorm,         // final double         plDispNorm,         // to increase weight of nearer planes
							clust_sigma,        // final double         sigma,
							disp_arange_hor,  // final double         disp_arange,
							disp_rrange_hor,  // final double         disp_rrange,
							tolerance_above_far,// final double         tolerance_above,
							tolerance_below_far,// final double         tolerance_below,
							plMinPoints,        // final int            plMinPoints, //          =     5;  // Minimal number of points for plane detection
							"hor",              // final String         suffix,
							debugLevel + 0,         // final int            debugLevel,
							dbg_X,              // final int            dbg_X,
							dbg_Y):             // final int            dbg_Y);

								dispClusterize(
										hor_disp_strength, // final double [][][][] disparity_strengths, // either normal or tilted disparity/strengths
										mmm_hor,           // final double [][][] hist_max_min_max, // histogram data: per tile array of odd number of disparity/strengths pairs, starting with first maximum
										null,              // final boolean [][][] selected,   // tiles OK to be assigned [supertile][measurement layer] [tile index] or null (or null or per-measurement layer)
										used_hor,          // final boolean [][][] prohibited, // already assigned tiles [supertile][measurement layer] [tile index] or null
										stMeasSel,         // final int        stMeasSel, //            = 1;      // Select measurements for supertiles : +1 - combo, +2 - quad +4 - hor +8 - vert
										((pass < (max_tries - 1)) ? max_diff_hor : Double.NaN), // final double     max_diff,  // maximal disparity difference (to assign to a cluster (of Double.NaN)
										plMinPoints,       // final int        plMinPoints, //          =     5;  // Minimal number of points for plane detection
										smallDiff,         // final double     smallDiff,  //       = 0.4;   // Consider merging initial planes if disparity difference below
										highMix,           // final double     highMix,    //stHighMix         = 0.4;   // Consider merging initial planes if jumps between ratio above
										true,              // final boolean    norm_max_diff,   // scale  max_diff for large (> dispNorm) average disparities
										true,              // final boolean    norm_small_diff, // scale  max_diff for large (> dispNorm) average disparities
										plDispNorm,        // final double     dispNorm, // TODO: make a separate variable?
										debugLevel, // 1, // debugLevel,
										dbg_X,
										dbg_Y);

			if (debugLevel > 0){
				System.out.println("initialDiscriminateTiles(): before new_planes_vert, pass =" + (pass + 1) + " ( of "+max_tries+" )");
			}
			final boolean [][][][] new_planes_vert = clusterize_by_highest ?
					dispClusterizeHighest(
							vert_disp_strength,   // final double [][][][] disparity_strengths, // either normal or tilted disparity/strengths
							null,                 // final boolean [][][] selected,   // tiles OK to be assigned [supertile][measurement layer] [tile index] or null (or null or per-measurement layer)
							used_vert,            // final boolean [][][] prohibited, // already assigned tiles [supertile][measurement layer] [tile index] or null
							false,                // final boolean        search_min,
							stMeasSel,            // final int            stMeasSel, //            = 1;      // Select measurements for supertiles : +1 - combo, +2 - quad +4 - hor +8 - vert
							plDispNorm,           // final double         plDispNorm,         // to increase weight of nearer planes
							clust_sigma,          // final double         sigma,
							disp_arange_vert,     // final double         disp_arange,
							disp_rrange_vert,     // final double         disp_rrange,
							tolerance_above_near, // final double         tolerance_above,
							tolerance_below_near, // final double         tolerance_below,
							plMinPoints,          // final int            plMinPoints, //          =     5;  // Minimal number of points for plane detection
							"vert",               // final String         suffix,
							debugLevel + 0,       // final int            debugLevel,
							dbg_X,                // final int            dbg_X,
							dbg_Y):               // final int            dbg_Y)
								dispClusterize(
										vert_disp_strength, // final double [][][][] disparity_strengths, // either normal or tilted disparity/strengths
										mmm_vert,           // final double [][][] hist_max_min_max, // histogram data: per tile array of odd number of disparity/strengths pairs, starting with first maximum
										null,              // final boolean [][][] selected,   // tiles OK to be assigned [supertile][measurement layer] [tile index] or null (or null or per-measurement layer)
										used_vert,              // final boolean [][][] prohibited, // already assigned tiles [supertile][measurement layer] [tile index] or null
										stMeasSel,         // final int        stMeasSel, //            = 1;      // Select measurements for supertiles : +1 - combo, +2 - quad +4 - hor +8 - vert
										((pass < (max_tries - 1)) ? max_diff_vert : Double.NaN), // final double     max_diff,  // maximal disparity difference (to assign to a cluster (of Double.NaN)
										plMinPoints,       // final int        plMinPoints, //          =     5;  // Minimal number of points for plane detection
										smallDiff,         // final double     smallDiff,  //       = 0.4;   // Consider merging initial planes if disparity difference below
										highMix,          // final double     highMix,    //stHighMix         = 0.4;   // Consider merging initial planes if jumps between ratio above
										true,              // final boolean    norm_max_diff,   // scale  max_diff for large (> dispNorm) average disparities
										true,              // final boolean    norm_small_diff, // scale  max_diff for large (> dispNorm) average disparities
										plDispNorm,        // final double     dispNorm, // TODO: make a separate variable?
										debugLevel, // 2, // debugLevel,
										dbg_X,
										dbg_Y);

			// compare which vert or hor provide stronger clusters, add them to selections, recalculate "used"
			ai.set(0);
			for (int ithread = 0; ithread < threads.length; ithread++) {
				threads[ithread] = new Thread() {
					@Override
					public void run() {
						TileNeibs tnStile = new TileNeibs(2 * superTileSize, 2* superTileSize);
						for (int nsTile = ai.getAndIncrement(); nsTile < nStiles; nsTile = ai.getAndIncrement()) {
                            int dl = ((debugLevel > 0) && (nsTile == debug_stile)) ? (debugLevel + 2) : debugLevel;
							if (dl > 2){
								System.out.println("initialDiscriminateTiles() selecting: nsTile="+nsTile);
							}
							double [][][][] ds = {vert_disp_strength[nsTile],hor_disp_strength[nsTile]};
							boolean [][][][] sels_all = {new_planes_vert[nsTile],new_planes_hor[nsTile]}; // make possible to iterate
							class SelStrength{
								int type;
								int indx;
								double strength;
								SelStrength(int type, int indx, double strength){
									this.type = type;
									this.indx = indx;
									this.strength = strength;
								}
							}
							ArrayList<SelStrength> selStrengthList = new ArrayList <SelStrength>();
							for (int pType = 0; pType < sels_all.length; pType++){
								if (sels_all[pType] == null) sels_all[pType] = new boolean[0][][];
								// calculate strength of each selection;
								for (int np = 0; np < sels_all[pType].length; np ++) {
									double sw = getSelectionStrength(
											mod_strength,        // boolean mod_strength,
											plDispNorm,          // double plDispNorm,
											sels_all[pType][np], // boolean [][] sels,
											ds[pType]);          // double [][][] ds);
									selStrengthList.add(new SelStrength(pType, np, sw));
								}
							}


							if (dl > 2){
								System.out.println("initialDiscriminateTiles() got list of clusters for "+nsTile);
								for (int i = 0; i < selStrengthList.size(); i++){
									System.out.println(i+": type = "+selStrengthList.get(i).type+
											" index = "+selStrengthList.get(i).indx +
											" strength = "+selStrengthList.get(i).strength);
								}
							}
							if (selStrengthList.size() > 0) {
								Collections.sort(selStrengthList, new Comparator<SelStrength>() {
									@Override
									public int compare(SelStrength lhs, SelStrength rhs) {
										// -1 - less than, 1 - greater than, 0 - equal, all inverted for descending
										return (lhs.strength > rhs.strength) ? -1 : (lhs.strength < rhs.strength ) ? 1 : 0;
									}
								});
								if (dl > 1){
									System.out.println("initialDiscriminateTiles() sorted list of clusters for "+nsTile);
									for (int i = 0; i < selStrengthList.size(); i++){
										System.out.println(i+": type = "+selStrengthList.get(i).type+
												" index = "+selStrengthList.get(i).indx +
												" strength = "+selStrengthList.get(i).strength);
									}
								}

								int best_type = selStrengthList.get(0).type;
								int num_planes = 0;
								for (; num_planes < selStrengthList.size(); num_planes++) {
									if (selStrengthList.get(num_planes).type != best_type) {
										break;
									}
								}
								int num_old_planes = (planes_selections[nsTile] == null) ? 0: planes_selections[nsTile].length;
								if (dl > 1){
									System.out.println("initialDiscriminateTiles() num_old_planes= "+num_old_planes);
								}
								boolean [][][] new_planes_selections = new boolean [num_old_planes + num_planes][][];
								boolean []     new_plane_types = new boolean [num_old_planes + num_planes];
								int np = 0;
								for (; np < num_old_planes; np++){
									new_planes_selections[np] = planes_selections[nsTile][np];
									new_plane_types[np] = hor_planes[nsTile][np];
								}
								int planes_added = 0;
								for (int nnp = 0; nnp < num_planes; nnp++) { // SelStrength ss:selStrengthList){
									SelStrength ss = selStrengthList.get(nnp);
									new_plane_types[np] =         ss.type != 0; // 0 - vert, 1 - hor;
									new_planes_selections[np] =   sels_all[ss.type][ss.indx];
									for (int ml = 0; ml < num_layers; ml++) if (sels_all[ss.type][ss.indx][ml] != null){
										if (used_hor[nsTile][ml] == null){
											used_hor[nsTile][ml] = new boolean[num_tiles];
										}
										// used_vert may be the same array as used_hor
										if (used_vert[nsTile][ml] == null){
											used_vert[nsTile][ml] = new boolean[num_tiles];
										}
										// Filling small gaps in new selections before marking them as "used"
										boolean [] used_before = (ss.type == 0) ? used_vert[nsTile][ml] : used_hor[nsTile][ml];
										boolean [] this_used = sels_all[ss.type][ss.indx][ml].clone();
										// See how many really new tiles has this_used;
										int num_really_new = 0;
										for (int i = 0; i < num_tiles; i++) if (this_used[i] && !used_vert[nsTile][ml][i] && !used_hor[nsTile][ml][i]) {
											num_really_new++;
										}
										// see if is enough
										if (num_really_new < min_fresh_used) {
											for (int i = 0; i < num_tiles; i++){
												used_before[i] |= this_used[i]; // Still mark them as used on this orientation
											}
											this_used = new boolean[this_used.length]; // unselect all

											if (dl > 1){
												System.out.println("initialDiscriminateTiles():"+nsTile+" selection does not provide enough previously unused tiles");
											}
											// disable these tiles in this orientation used

										} else {
											planes_added++;
										}

										boolean [] this_used_expanded = this_used.clone();
										tnStile.growSelection(
												2,
												this_used_expanded,
												null);
										for (int i = 0; i < num_tiles; i++) if (this_used_expanded[i] && !this_used[i]) {
											int n_new_used = 0, n_blocked = 0;
											for (int dir = 0; dir < 8; dir++) {
												int i1 = tnStile.getNeibIndex(i, dir);
												if (i1 < 0) {
													n_blocked++;
												} else if (sels_all[ss.type][ss.indx][ml][i1]) { // this_used is modified
													n_new_used++;
												} else if (used_before[i1]) {
													n_blocked++;
												}
											}
											if (((n_new_used + n_blocked) >= used_companions) && (n_new_used >= used_true_companions)) {
												this_used[i] = true;
											}
										}

										// first apply to the same direction (or both if it is the same array)
										if (ss.type == 0) {
											for (int i = 0; i < num_tiles; i++){
												used_vert[nsTile][ml][i] |= this_used[i]; // sels_all[ss.type][ss.indx][ml][i];
											}
										} else {
											for (int i = 0; i < num_tiles; i++){
												used_hor[nsTile][ml][i] |= this_used[i]; // sels_all[ss.type][ss.indx][ml][i];
											}
										}
										if (hor_vert_overlap > 0) { // shrink selection and apply to the other direction
//											boolean [] newsel = sels_all[ss.type][ss.indx][ml].clone();
											tnStile.shrinkSelection(
													hor_vert_overlap,
													this_used,
													null);
											if (ss.type != 0) {
												for (int i = 0; i < num_tiles; i++){
													used_vert[nsTile][ml][i] |= this_used[i];
												}
											} else {
												for (int i = 0; i < num_tiles; i++){
													used_hor[nsTile][ml][i] |= this_used[i];
												}
											}
										}
									}
									np++;
								}
								if (planes_added == 0) {
									if (dl > 1){
										System.out.println("initialDiscriminateTiles():"+nsTile+": No new planes added");
									}
									continue;
								}

								planes_selections[nsTile] = new_planes_selections;
								hor_planes[nsTile] =        new_plane_types;
								if (dl > 2){
									System.out.println("initialDiscriminateTiles() new_planes_selections.length= "+new_planes_selections.length);
								}
							}
						}
					}
				};
			}
			ImageDtt.startAndJoin(threads);
		}
		return planes_selections;
	}

	/**
	 * Re-assign tiles to the planes according to fitted planes. Scans each known plane parallel with specified
	 * size steps and specified number of steps, maximizing overall fit quality.
	 * Total weight of the fitted tiles? RMS?
	 * @param planes per-supertile, per plane array of plane data instances. Should have
	 *        nonexclusiveStar and nonexclusiveStarEq calculated
	 * @param stMeasSel - select measurements for supertiles : +1 - combo, +2 - quad +4 - hor +8 - vert
	 * @param plDispNorm - normalization of the measured disparity precision - closer objects will have disparity
	 *        difference proportionally reduced
	 * @param plPreferDisparity kept for historical reasons - when true, select disparity-most vector even if it has higher eigenvalue
	 * @param geometryCorrection GeometryCorrection instance to use
	 * @param correct_distortions correct geometrical distortions when converting to/from world coordinates
	 * @param smplMode sample mode
	 * @param smplSide sample side for averaging/filtering tile disparities
	 * @param smplNum number of best samples used fro averaging
	 * @param smplRms maximal sample disparity rms to consider sample valid
	 * @param smplWnd use window functions for the samples
	 * @param max_abs_tilt pix per tile
	 * @param max_rel_tilt pix / disparity) per tile
	 * @param damp_tilt Damp tilt to handle insufficient  (co-linear)data
	 * @param min_tilt_disp Disparity switch between filtering modes - near objects use tilts, far - use max disparity
	 * @param transition Mode transition range (between tilted and maximal disparity)
	 * @param far_mode Far objects filtering mode (0 - off, 1 - power of disparity)
	 * @param far_power Raise disparity to this power before averaging for far objects
	 * @param plDiscrTolerance maximal disparity difference from the plane to consider tile
	 * @param plDiscrDispRange parallel move known planes around original know value for the best overall fit
	 * @param plDiscrSteps number of steps (each direction) for each plane to search for the best fit (0 - single, 1 - 1 each side)
	 * @param plDiscrVariants total number of variants to try (protect from too many planes)
	 * @param plDiscrMode what plane to use as a hint: 0 - weighted, 1 - equalized, 2 - best, 3 - combined
	 * @param plDiscrVarFloor squared add to variance to calculate reverse flatness (used mostly for single-cell clusters)
	 * @param plDiscrSigma Gaussian sigma to compare how measured data is attracted to planes
	 * @param plDiscrBlur sigma to blur histograms while re-discriminating
	 * @param plDiscrExclusivity tile exclusivity: 1.0 - tile belongs to one plane only, 0.0 - regardless of others
	 * @param debugLevel debug level
	 * @param dbg_X supertile X for elevated debug level
	 * @param dbg_Y supertile X for elevated debug level
	 * @return per-tile, per plane, per measurement layer (type of correlation - combo, 4, hor, vert), per tile selections of
	 * filtered tiles
	 */

	public boolean [][][][]  refineDiscriminateTiles(
			final TilePlanes.PlaneData [][] planes,
			final int        stMeasSel,            //      = 1;      // Select measurements for supertiles : +1 - combo, +2 - quad +4 - hor +8 - vert
			final double     plDispNorm,
			final boolean    plPreferDisparity,    // Always start with disparity-most axis (false - lowest eigenvalue)
			final GeometryCorrection geometryCorrection,
			final boolean    correct_distortions,

			final boolean    smplMode,             //        = true;   // Use sample mode (false - regular tile mode)

			final MeasuredLayersFilterParameters mlfp,

			final double     plDiscrTolerance,     //     =   0.4;  // Maximal disparity difference from the plane to consider tile
			final double     plDiscrDispRange,     //     =   0.6;  // Parallel move known planes around original know value for the best overall fit
			final int        plDiscrSteps,         //         =   3;    // Number of steps (each direction) for each plane to search for the best fit (0 - single, 1 - 1 each side)
			final int        plDiscrMode,          //          =   3;    // What plane to use as a hint: 0 - weighted, 1 - equalized, 2 - best, 3 - combined

			final double     plDiscrVarFloor,    //       =   0.03;  // Squared add to variance to calculate reverse flatness (used mostly for single-cell clusters)
			final double     plDiscrSigma,       //          =   0.05;  // Gaussian sigma to compare how measured data is attracted to planes
			final double     plDiscrBlur,        //           =   0.1;   // Sigma to blur histograms while re-discriminating
			final double     plDiscrExclusivity, //    =   1.5;   // Tile exclusivity: 1.0 - tile belongs to one plane only, 0.0 - regardless of others
			final double     plDiscrExclus2,     //        =   0.8;   // For second pass if exclusivity > 1.0 - will assign only around strong neighbors
			final boolean    plDiscrStrict,      //         = true;   // When growing selection do not allow any offenders around (false - more these than others)
			final double     plDiscrCorrMax,     //        = 0.7;   // Attraction to different planes correlation that is too high for re-discrimination.
			final double     plDiscrCorrMerge,   //     = 0.85;  // Attraction to different planes correlation that is high enough to merge planes
			final int        plDiscrSteal,       //         =   4;     // If offender has this number of tiles (including center) the cell can not be used
			final int        plDiscrGrown,       //         =   0;     // Only use tiles within this range from original selection
			final double     plDiscrXMedian,     //        = 1.5;   // Remove outliers from the final selection that have distance more than scaled median
			final int        debugLevel,
			final int        dbg_X,
			final int        dbg_Y)
	{
		// create a list of usable planes according to the mode
		final int tilesX =        tileProcessor.getTilesX();
		final int tilesY =        tileProcessor.getTilesY();
		final int superTileSize = tileProcessor.getSuperTileSize();
		final int tileSize =      tileProcessor.getTileSize();

		final int stilesX = (tilesX + superTileSize -1)/superTileSize;
		final int stilesY = (tilesY + superTileSize -1)/superTileSize;
		final int nStiles = stilesX * stilesY;
		final TilePlanes tpl = new TilePlanes(tileSize,superTileSize, geometryCorrection);
		final Thread[] threads = ImageDtt.newThreadArray((debugLevel > 1)? 1 : tileProcessor.threadsMax);

		final AtomicInteger ai = new AtomicInteger(0);
		this.planes = new TilePlanes.PlaneData[nStiles][];
		final int debug_stile = (debugLevel > -1)? (dbg_Y * stilesX + dbg_X):-1;
		final boolean [][][][] planes_selections = new boolean [nStiles][][][]; // num_tiles
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				@Override
				public void run() {
					for (int nsTile = ai.getAndIncrement(); nsTile < nStiles; nsTile = ai.getAndIncrement()) {
//						int dl =  ((debugLevel > -1) && (nsTile == debug_stile)) ? 3 : 1;
                        int dl = ((debugLevel > 1) && (nsTile == debug_stile)) ? 3: debugLevel;

						if (dl > 1){
							System.out.println("refineDiscriminateTiles() selecting: nsTile="+nsTile);
						}
						if (planes[nsTile] != null) {
							int stileY = nsTile / stilesX;
							int stileX = nsTile % stilesX;
							int [] sTiles = {stileX, stileY};
							TilePlanes.PlaneData pd0 = tpl.new  PlaneData (
									sTiles, // int [] sTileXY,
									tileSize, // int tileSize,
									geometryCorrection, // GeometryCorrection   geometryCorrection,
									correct_distortions,
									measuredLayers,     // MeasuredLayers measuredLayers,
									plPreferDisparity);   // boolean preferDisparity)

							TilePlanes.PlaneData [] these_planes = planes[nsTile].clone();  // null pointer

							for (int ntry = 0; ntry < planes[nsTile].length; ntry ++) { // can be while(true), just for debugging
								int [] merge_planes = {-1, -1};
								planes_selections[nsTile] = pd0.reDiscriminateTiles(
										""+nsTile,         // String           prefix,
										these_planes,  //planes[nsTile],    // final PlaneData  [] planes,
										merge_planes,      // final int []     merge_planes, // indices of planes suggested to be merged
										stMeasSel,         // final int        stMeasSel, //            = 1;      // Select measurements for supertiles : +1 - combo, +2 - quad +4 - hor +8 - vert
										plDispNorm,        // final double     dispNorm,   //  Normalize disparities to the average if above
										smplMode,          // final boolean    smplMode, //        = true;   // Use sample mode (false - regular tile mode)
										mlfp,

										plDiscrTolerance,  // final double     disp_tolerance,   // maximal disparity difference from the plane to consider tile
										plDiscrVarFloor,   // final double     disp_var_floor,   // squared add to variance to calculate reverse flatness (used mostly for single-cell clusters)
										plDiscrSigma,      // final double     disp_sigma,       // G.sigma to compare how measured data is attracted to planes
										plDiscrDispRange,  // final double     disp_range,       // parallel move known planes around original know value for the best overall fit
										plDiscrSteps,      // final int        amplitude_steps,  // number of steps (each direction) for each plane to search for the best fit (0 - single, 1 - 1 each side)
										plDiscrBlur,       // final double     hist_blur,        // Sigma to blur histogram
										plDiscrExclusivity, // final double    exclusivity,      // 1.0 - tile belongs to one plane only, 0.0 - regardless of others
										plDiscrExclus2,    // final double     excluisivity2, //        =   0.8;   // For second pass if exclusivity > 1.0 - will assign only around strong neighbors
										plDiscrStrict,     // final boolean    exclusivity_strict//         = true;   // When growing selection do not allow any offenders around (false - more these than others)
										plDiscrCorrMax,    // final double     plDiscrCorrMax, //         = 0.7;   // Attraction to different planes correlation that is too high for re-discrimination.
										plDiscrCorrMerge,  // final double     attractionCorrMerge, //         = 0.85;  // Attraction to different planes correlation that is high enough to merge planes
										plDiscrSteal,      // final int        plDiscrSteal,        //         =   4;     // If offender has this number of tiles (including center) the cell can not be used
										plDiscrGrown,      // final int        plDiscrGrown,        //         =   0;     // Only use tiles within this range from original selection
										plDiscrXMedian,    // final double     outliersXMedian, //         = 1.5;   // Remove outliers from the final selection that have distance more than scaled median
										plDiscrMode,       // final int        mode, // 0 - weighted, 1 - equalized, 2 - best, 3 - combined
										dl);               // debugLevel); // final int        debugLevel)
								if ((planes_selections[nsTile] != null) || (merge_planes[0] < 0)){ // either OK, or does not know what to do (keep old)
									break;
								}
								// merging suggested plane pair
								if (debugLevel > 0) {
									System.out.println("refineDiscriminateTiles(): nsTile="+nsTile+" merging pair ["+merge_planes[0]+","+merge_planes[1]+"]");
								}
								TilePlanes.PlaneData [] new_planes = new TilePlanes.PlaneData [these_planes.length -1];
								int np1 = 0;
								for (int np = 0; np < these_planes.length; np++) {
									if (np != merge_planes[1]){
										new_planes[np1++] = these_planes[np];
									}
								}
								TilePlanes.PlaneData plane1 = these_planes[merge_planes[0]].mergePlaneToThis(
										these_planes[merge_planes[1]], // PlaneData otherPd,
										1.0,         // double    scale_other,
										1.0,         // double     starWeightPwr,    // Use this power of tile weight when calculating connection cost
										false,       // boolean   ignore_weights,
										true, // boolean   sum_weights,
										these_planes[merge_planes[0]].getPreferDisparity(), // preferDisparity,
										dl-2); // int       debugLevel)
								// combine tile selection - if next time pd0.reDiscriminateTiles() will fail, it will
								// use old selections, we need to provide them (otherwise will use selection from the first plane)
								if (plane1 == null){
									System.out.println("refineDiscriminateTiles() nsTile="+nsTile+" plane1 = null");
									break;
								}
								plane1.orMeasSelection(these_planes[merge_planes[1]].getMeasSelection());

								// separately merge corresponding nonexclusiveStar and nonexclusiveStarEq of these planes - kit is not exact,
								// but is needed just for a hint and is compatible with multithreading without recalculating other planes

								TilePlanes.PlaneData    plane1Ex =    these_planes[merge_planes[0]].getNonexclusiveStarFb();
								TilePlanes.PlaneData    plane1ExEq =  these_planes[merge_planes[0]].getNonexclusiveStarEqFb();

								TilePlanes.PlaneData    plane2Ex =    these_planes[merge_planes[1]].getNonexclusiveStarFb();
								TilePlanes.PlaneData    plane2ExEq =  these_planes[merge_planes[1]].getNonexclusiveStarEqFb();

								TilePlanes.PlaneData plane1NonExcl = plane1Ex.mergePlaneToThis(
										plane2Ex, // PlaneData otherPd,
										1.0,         // double    scale_other,
										1.0,         // double     starWeightPwr,    // Use this power of tile weight when calculating connection cost
										false,       // boolean   ignore_weights,
										true, // boolean   sum_weights,
										these_planes[merge_planes[0]].getPreferDisparity(), // preferDisparity,
										dl-2); // int       debugLevel)

								TilePlanes.PlaneData plane1NonExclEq = plane1ExEq.mergePlaneToThis(
										plane2ExEq, // PlaneData otherPd,
										1.0,         // double    scale_other,
										1.0,         // double     starWeightPwr,    // Use this power of tile weight when calculating connection cost
										true,        // boolean   ignore_weights,
										true,        // boolean   sum_weights,
										these_planes[merge_planes[0]].getPreferDisparity(), // preferDisparity,
										dl-2); // int       debugLevel)
								plane1.setNonexclusiveStar  (plane1NonExcl);
								plane1.setNonexclusiveStarEq(plane1NonExclEq);
								new_planes[merge_planes[0]] = plane1;
								if (dl > 1){
									System.out.println("First plane ("+merge_planes[0]+"):\n"+(these_planes[merge_planes[0]].toString()));
									System.out.println("Second plane ("+merge_planes[1]+"):\n"+(these_planes[merge_planes[1]].toString()));
									System.out.println("combined plane:\n"+(plane1.toString()));
								}
								these_planes = new_planes.clone();
							}
						}
					}
				}
			};
		}
		ImageDtt.startAndJoin(threads);
		return planes_selections;
	}


	public TilePlanes.PlaneData [][] createPlanesFromSelections( // never finished?
			final boolean [][][][] plane_selections, //  = new boolean [nStiles][][][]; // num_tiles
			final double  [][][][] disp_strength,
			final double     plDispNorm,
			final int        plMinPoints, //          =     5;  // Minimal number of points for plane detection
			final double     plTargetEigen, //        =   0.1;  // Remove outliers until main axis eigenvalue (possibly scaled by plDispNorm) gets below
			final double     plFractOutliers, //      =   0.3;  // Maximal fraction of outliers to remove
			final int        plMaxOutliers, //        =    20;  // Maximal number of outliers to remove
//			final double     plVertWors,    //        =    1.5  // if rotating plane vertical does not increase 'eigenvalue' more, use vertical
			final boolean    plPreferDisparity, // Always start with disparity-most axis (false - lowest eigenvalue)
			final GeometryCorrection geometryCorrection,
			final boolean    correct_distortions,
			final boolean    smplMode, //        = true;   // Use sample mode (false - regular tile mode)

			final MeasuredLayersFilterParameters mlfp,
			final double     fronto_tol, // fronto tolerance (pix) - treat almost fronto as fronto (constant disparity). <= 0 - disable this feature
			//FIXME: use following 2 parameters
			final double     fronto_rms,    // Target rms for the fronto planes - same as sqrt(plMaxEigen) for other planes
			final double     fronto_offs,   //        =   0.2;  // increasing weight of the near tiles by using difference between the reduced average as weight. <= 0 - disable
			final double     fronto_pow,    //        =   1.0;  // increase weight even more
			// now for regenerated planes - just null as it is not known if it is hor or vert
			final boolean [][] hor_planes, // plane types (hor/vert)
			final int        debugLevel,
			final int        dbg_X,
			final int        dbg_Y)
	{
		final int tilesX =        tileProcessor.getTilesX();
		final int tilesY =        tileProcessor.getTilesY();
		final int superTileSize = tileProcessor.getSuperTileSize();
		final int tileSize =      tileProcessor.getTileSize();

		final int stilesX = (tilesX + superTileSize -1)/superTileSize;
		final int stilesY = (tilesY + superTileSize -1)/superTileSize;
		final int nStiles = stilesX * stilesY;
//		final Thread[] threads = ImageDtt.newThreadArray(tileProcessor.threadsMax);
		final Thread[] threads = ImageDtt.newThreadArray((debugLevel > 1)? 1 : tileProcessor.threadsMax);

		final AtomicInteger ai = new AtomicInteger(0);
		final TilePlanes.PlaneData [][] result_planes = new TilePlanes.PlaneData[nStiles][];
		//		this.planes = new TilePlanes.PlaneData[nStiles][];
		final int debug_stile = (debugLevel > -1)? (dbg_Y * stilesX + dbg_X):-1;
		// TODO: Remove when promoting PlaneData
		final TilePlanes tpl = new TilePlanes(tileSize,superTileSize, geometryCorrection);

		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				@Override
				public void run() {
					for (int nsTile = ai.getAndIncrement(); nsTile < nStiles; nsTile = ai.getAndIncrement()) {
                        int dl = ((debugLevel > 1) && (nsTile == debug_stile)) ? 3: debugLevel;
						if (dl > 2){
							System.out.println("createPlanesFromSelections(): nsTile="+nsTile);
						}
						if (plane_selections[nsTile] != null) {

							int stileY = nsTile / stilesX;
							int stileX = nsTile % stilesX;
							int [] sTiles = {stileX, stileY};

							result_planes[nsTile] = null;
							// first make a plane from all tiles
							TilePlanes.PlaneData pd0 = tpl.new  PlaneData (
									sTiles, // int [] sTileXY,
									tileSize, // int tileSize,
									geometryCorrection, // GeometryCorrection   geometryCorrection,
									correct_distortions,
									measuredLayers,     // MeasuredLayers measuredLayers,
									plPreferDisparity);   // boolean preferDisparity)
							ArrayList<TilePlanes.PlaneData> st_planes = pd0.createTilePlanesFromSelections(
											"" + nsTile, // String        suffix,
											plane_selections[nsTile], // boolean [][][] plane_selections, //  = new boolean [nStiles][][][]; // num_tiles
											((hor_planes == null)? null:hor_planes[nsTile]), // boolean []     hor_planes,
											disp_strength[nsTile],    // double  [][][] disp_strength,
											plDispNorm,               // double       dispNorm,   //  Normalize disparities to the average if above
											plMinPoints,              // int          min_tiles,
											plTargetEigen,            // double       plTargetEigen, //        =   0.1;  // Remove outliers until main axis eigenvalue (possibly scaled by plDispNorm) gets below
											plFractOutliers,          // double       plFractOutliers, //      =   0.3;  // Maximal fraction of outliers to remove
											plMaxOutliers,            // int          plMaxOutliers, //        =    20;  // Maximal number of outliers to remove
											correct_distortions,      // boolean      correct_distortions,
											smplMode,                 // boolean      smplMode, //        = true;   // Use sample mode (false - regular tile mode)
											mlfp,
											fronto_tol,               // double       fronto_tol, // fronto tolerance (pix) - treat almost fronto as fronto (constant disparity). <= 0 - disable this feature
											fronto_rms,               // double       fronto_rms,  // Target rms for the fronto planes - same as sqrt(plMaxEigen) for other planes. May be tighter
											fronto_offs,              // double       fronto_offs,   //        =   0.2;  // increasing weight of the near tiles by using difference between the reduced average as weight. <= 0 - disable
											fronto_pow,               // double       fronto_pow,    //        =   1.0;  // increase weight even more
											dl);                      // int          debugLevel);

							if ((st_planes != null) && (!st_planes.isEmpty())){
								if (dl > 2){
									System.out.println("======= createPlanesFromSelections(): nsTile="+nsTile+" detecting bridges ==========");
								}
								boolean [][] hs = new boolean [1][];
								boolean [][][] split_sels = pd0.filterBridges(
										st_planes, // ArrayList<PlaneData> tilePlanes,
										hs, // ((hor_planes == null)? null:hor_planes[nsTile]), // boolean []     hor_planes,
										3, // int max_grow, // make configurable?
										3, // int max_grow_far,
										dl); // int debugLevel)
								if (split_sels !=null){
									if (hor_planes != null) {
										hor_planes[nsTile] = hs[0];
									}
									if (dl > 2){
										System.out.println("======= createPlanesFromSelections(): nsTile="+nsTile+" removing bridges ==========");
									}
									if (dl > 2) {
										String [] dbg_titles = showSupertileSeparationTitles( disp_strength[nsTile], plane_selections[nsTile], result_planes[nsTile]);
										double [][] dbg_img =  showSupertileSeparation(false, disp_strength[nsTile], plane_selections[nsTile], result_planes[nsTile]);
										ShowDoubleFloatArrays sdfa_instance = new ShowDoubleFloatArrays();
										sdfa_instance.showArrays(dbg_img, 2 * superTileSize, 2* superTileSize, true, "pre_bridge_disp-"+nsTile+"-"+debugLevel,dbg_titles);
										dbg_img =  showSupertileSeparation(true, disp_strength[nsTile], plane_selections[nsTile], result_planes[nsTile]);
										sdfa_instance.showArrays(dbg_img, 2 * superTileSize, 2* superTileSize, true, "pre_bridge_world-"+nsTile+"-"+debugLevel,dbg_titles);
									}

									st_planes = pd0.createTilePlanesFromSelections(
											"" + nsTile+" bridges ",  // String        suffix,
											split_sels,               // boolean [][][] plane_selections, //  = new boolean [nStiles][][][]; // num_tiles
											// FIXME: replace null (from what type it was before detecting bridges)
											((hor_planes == null)? null:hor_planes[nsTile]), // boolean []     hor_planes,
//											null,                     // boolean []     hor_planes,
											disp_strength[nsTile],    // double  [][][] disp_strength,
											plDispNorm,               // double       dispNorm,   //  Normalize disparities to the average if above
											plMinPoints,              // int          min_tiles,
											plTargetEigen,            // double       plTargetEigen, //        =   0.1;  // Remove outliers until main axis eigenvalue (possibly scaled by plDispNorm) gets below
											0.0,                      // double       plFractOutliers, //      =   0.3;  // Maximal fraction of outliers to remove
											0,                        // plMaxOutliers,            // int          plMaxOutliers, //        =    20;  // Maximal number of outliers to remove
											correct_distortions,      // boolean      correct_distortions,
											smplMode,                 // boolean      smplMode, //        = true;   // Use sample mode (false - regular tile mode)
											mlfp,
											fronto_tol,               // double       fronto_tol, // fronto tolerance (pix) - treat almost fronto as fronto (constant disparity). <= 0 - disable this feature
											fronto_rms,               // double       fronto_rms,  // Target rms for the fronto planes - same as sqrt(plMaxEigen) for other planes. May be tighter
											fronto_offs,              // double       fronto_offs,   //        =   0.2;  // increasing weight of the near tiles by using difference between the reduced average as weight. <= 0 - disable
											fronto_pow,               // double       fronto_pow,    //        =   1.0;  // increase weight even more
											dl - 1);                  // int          debugLevel);
								}

							}

							if ((st_planes != null) && (!st_planes.isEmpty())){
								if (dl > 2) {
									for (TilePlanes.PlaneData plane:st_planes){
										plane.getWorldXYZ(0);
										System.out.println(plane.toString());
									}
									// Calculate planes and print results
								}
//						this_new_plane.getWorldXYZ(0);

								if (LOWEST_PLANE(2) > 0) st_planes.add(0, st_planes.get(0)); // insert dummy at pos 0;
								result_planes[nsTile] = st_planes.toArray(new TilePlanes.PlaneData[0] );
								if (LOWEST_PLANE(2) > 0) result_planes[nsTile][0] = null; // remove dummy
								if (dl > 2){
									System.out.println("createPlanesFromSelections(): nsTile="+nsTile);
								}
								if (dl > 3) { // 2) {
									String [] dbg_titles = showSupertileSeparationTitles( disp_strength[nsTile], plane_selections[nsTile], result_planes[nsTile]);
									double [][] dbg_img =  showSupertileSeparation(false, disp_strength[nsTile], plane_selections[nsTile], result_planes[nsTile]);
									ShowDoubleFloatArrays sdfa_instance = new ShowDoubleFloatArrays();
									sdfa_instance.showArrays(dbg_img, 2 * superTileSize, 2* superTileSize, true, "create_planes_disp-"+nsTile+"-"+debugLevel,dbg_titles);
									dbg_img =  showSupertileSeparation(true, disp_strength[nsTile], plane_selections[nsTile], result_planes[nsTile]);
									sdfa_instance.showArrays(dbg_img, 2 * superTileSize, 2* superTileSize, true, "create_planes_world-"+nsTile+"-"+debugLevel,dbg_titles);
									System.out.println("createPlanesFromSelections(): disp_strength["+nsTile+"][0][0]:");
									for (int iy = 0; iy < 2 * superTileSize; iy++){
										for (int ix = 0; ix < 2 * superTileSize; ix++){
											System.out.print(disp_strength[nsTile][0][0][2 * superTileSize * iy + ix]);
											if (ix < (2 * superTileSize-1)){
												System.out.print(", ");
											} else {
												System.out.println();
											}
										}
									}
									System.out.println("createPlanesFromSelections(): disp_strength["+nsTile+"][0][1]:");
									for (int iy = 0; iy < 2 * superTileSize; iy++){
										for (int ix = 0; ix < 2 * superTileSize; ix++){
											System.out.print(disp_strength[nsTile][0][1][2 * superTileSize * iy + ix]);
											if (ix < (2 * superTileSize - 1)){
												System.out.print(", ");
											} else {
												System.out.println();
											}
										}
									}
									System.out.println();

								}
							}
						}
					}
				}
			};
		}
		ImageDtt.startAndJoin(threads); // Never finished?
		return result_planes;
	}

	public void processPlanes5(
			final int        growSelection,                     // grow initial selection before processing
			final int        stMeasSel, //            = 1;      // Select measurements for supertiles : +1 - combo, +2 - quad +4 - hor +8 - vert
			final double     plDispNorm,
			final int        plMinPoints, //          =     5;  // Minimal number of points for plane detection
			final double     plTargetEigen, //        =   0.1;  // Remove outliers until main axis eigenvalue (possibly scaled by plDispNorm) gets below
			final double     plFractOutliers, //      =   0.3;  // Maximal fraction of outliers to remove
			final int        plMaxOutliers, //        =    20;  // Maximal number of outliers to remove
			final boolean    plPreferDisparity, // Always start with disparity-most axis (false - lowest eigenvalue)
			final double     plFrontoTol,                       // fronto tolerance (pix) - treat almost fronto as fronto (constant disparity). <= 0 - disable
			final double     plFrontoRms,                       // Target rms for the fronto planes - same as sqrt(plMaxEigen) for other planes
			final double     plFrontoOffs,                      // increasing weight of the near tiles by using difference between the reduced average as weight. <= 0 - disable
			final double     PlFrontoPow,    //        =   1.0;  // increase weight even more

			final GeometryCorrection geometryCorrection,
			final boolean    correct_distortions,

			final boolean    smplMode, //        = true;   // Use sample mode (false - regular tile mode)

			final MeasuredLayersFilterParameters mlfp,

			final double     bin_blur_hor,   // Blur disparity histograms for horizontal clusters by this sigma (in bins)
			final double     bin_blur_vert,  // Blur disparity histograms for constant disparity clusters by this sigma (in bins)
			final double     max_diff_hor,   // maximal disparity difference (to assign to a cluster (of Double.NaN) at first run for horizontal planes
			final double     max_diff_vert,  // maximal disparity difference (to assign to a cluster (of Double.NaN) at first run for vertical plane
			final int        max_tries,       // on last run - assign all remaining pixels to some cluster (disregard max_diff)

			final boolean    msUseSel,        // final boolean                   use_sel,
			final boolean    msDivideByArea,  // final boolean                   divide_by_area,
			final double     msScaleProj,     // final double                    scale_projection,

			final double     smallDiff,  //       = 0.4;   // Consider merging initial planes if disparity difference below
			final double     highMix,    //stHighMix         = 0.4;   // Consider merging initial planes if jumps between ratio above
			final double []  world_hor, // horizontal plane normal (default [0.0, 1.0, 0.0])
			final boolean    show_histograms,

			// Parameters for alternative initial planes that use lowest disparity for fronto planes, and farthest - for horizontal
			final boolean    mod_strength,          //  = true; // FIXME: make a parameter. when set, multiply each tile strength by the number of selected neighbors
			final boolean    clusterize_by_highest, //  = true;
			final double     clust_sigma,           //  = 0.7;
			final double     disp_arange_vert,      //  = 0.07;
			final double     disp_rrange_vert,      //  = 0.01;
			final double     disp_arange_hor,       //  =  0.035;
			final double     disp_rrange_hor,       //  =   0.005;
			final double     tolerance_above_near,  //  =  100.0; // 0.07; any?
			final double     tolerance_below_near,  //  =  -0.01;
			final double     tolerance_above_far,   //  =    0.07;
			final double     tolerance_below_far,   //  =    0.1; // 100.0; // any farther
			final int        hor_vert_overlap,      //  =       2;
			final int        used_companions,       //  =      5; // cell that has this many new used companions is considered used (borders and already use3d are considered used too)
			final int        used_true_companions,  //  = 1; // there should be at least this many new selected tiles among neighbors.,
			final boolean    debug_initial_discriminate,
			final int        debugLevel,
			final int        dbg_X,
			final int        dbg_Y)
	{
		// use both horizontal and const disparity tiles to create tile clusters
		// Add max_diff (maximal disparity difference while extracting initial tile selection) and max_tries (2..3) parameters

		// Add separate method to create + remove outliers from all planes (2 different ones)?
		// TODO later re-assign pixels according to existing plane parameters
		// Sort plane data by center (plane or supertile) disparity
		final int superTileSize = tileProcessor.getSuperTileSize();
		final int nStiles = ((tileProcessor.getTilesX() + superTileSize -1)/superTileSize) * ((tileProcessor.getTilesY() + superTileSize -1)/superTileSize);

		boolean [][] hor_planes = new boolean [nStiles][];

		boolean [][][][]  plane_selections = initialDiscriminateTiles(
				growSelection,        // final int        growSelection,                     // grow initial selection before processing
				stMeasSel,            // final int        stMeasSel,       //      = 1;      // Select measurements for supertiles : +1 - combo, +2 - quad +4 - hor +8 - vert
				plDispNorm,           // final double     plDispNorm,
				plMinPoints,          // final int        plMinPoints, //          =     5;  // Minimal number of points for plane detection
				plPreferDisparity,    // final boolean    plPreferDisparity, // Always start with disparity-most axis (false - lowest eigenvalue)
				geometryCorrection,  // final GeometryCorrection geometryCorrection,
				correct_distortions, // final boolean    correct_distortions,

				smplMode,            // final boolean    smplMode, //        = true;   // Use sample mode (false - regular tile mode)
				mlfp,
				bin_blur_hor,   // final double     bin_blur_hor,   // Blur disparity histograms for horizontal clusters by this sigma (in bins)
				bin_blur_vert,  // final double     bin_blur_vert,  // Blur disparity histograms for constant disparity clusters by this sigma (in bins)
				max_diff_hor,   // final double     max_diff_hor,   // maximal disparity difference (to assign to a cluster (of Double.NaN) at first run for horizontal planes
				max_diff_vert,  // final double     max_diff_vert,  // maximal disparity difference (to assign to a cluster (of Double.NaN) at first run for vertical plane
				max_tries,       //final int        max_tries,       // on last run - assign all rfemaining pixels to some cluster (disregard max_diff)
				smallDiff,      // final double     smallDiff,  //       = 0.4;   // Consider merging initial planes if disparity difference below
				highMix,    //final double     highMix,    // stHighMix         = 0.4;   // Consider merging initial planes if jumps between ratio above
				world_hor, // final double []  world_hor, // horizontal plane normal (default [0.0, 1.0, 0.0])
				show_histograms, // final boolean    show_histograms,
				hor_planes, // final boolean [][] hor_planes,

				// Parameters for alternative initial planes that use lowest disparity for fronto planes, and farthest - for horizontal
				mod_strength,          // final boolean    mod_strength,          //  = true; // FIXME: make a parameter. when set, multiply each tile strength by the number of selected neighbors
				clusterize_by_highest, // final boolean    clusterize_by_highest, //  = true;
				clust_sigma,           // final double     clust_sigma,           //  = 0.7;
				disp_arange_vert,      // final double     disp_arange_vert,      //  = 0.07;
				disp_rrange_vert,      // final double     disp_rrange_vert,      //  = 0.01;
				disp_arange_hor,       // final double     disp_arange_hor,       //  =  0.035;
				disp_rrange_hor,       // final double     disp_rrange_hor,       //  =   0.005;
				tolerance_above_near,  // final double     tolerance_above_near,  //  =  100.0; // 0.07; any?
				tolerance_below_near,  // final double     tolerance_below_near,  //  =  -0.01;
				tolerance_above_far,   // final double     tolerance_above_far,   //  =    0.07;
				tolerance_below_far,   // final double     tolerance_below_far,   //  =    0.1; // 100.0; // any farther
				hor_vert_overlap,      // final int        hor_vert_overlap,      //  =       2;
				used_companions,       // final int        used_companions,       //  =      5; // cell that has this many new used companions is considered used (borders and already use3d are considered used too)
				used_true_companions,  // final int        used_true_companions,  //  = 1; // there should be at least this many new selected tiles among neighbors.,

				debugLevel+(debug_initial_discriminate? 2:0), // final int        debugLevel,
				dbg_X, // final int        dbg_X,
				dbg_Y); // final int        dbg_Y)

		// get per-tile disparity strength again (may consider using non-filtered data here)

		double [][][][] disp_strength = getPlaneDispStrengthsST( // here use actual disparity, not tilted

				null,                // final double []  world_plane_norm, // null - floating, (0,0,1) - const disparity, (0,1,0) - horizontal
				stMeasSel,           // final int        stMeasSel, //            = 1;      // Select measurements for supertiles : +1 - combo, +2 - quad +4 - hor +8 - vert

				plPreferDisparity,   // final boolean    plPreferDisparity, // Always start with disparity-most axis (false - lowest eigenvalue)
				geometryCorrection,  // final GeometryCorrection geometryCorrection,
				correct_distortions, // final boolean    correct_distortions,

				smplMode,            // final boolean    smplMode, //        = true;   // Use sample mode (false - regular tile mode)
				mlfp,
				debugLevel,          // final int        debugLevel,
				dbg_X,               // final int        dbg_X,
				dbg_Y);              // final int        dbg_Y)
		// Create plane data (ellipsoids) from the disparity/strength data and tile selections,
		// remove outliers, order each supertile result in ascending disparity (from far to near)
		// TODO: consider 1) assigning of the non-assigned tiles to clusters and 2) re-assigning all clusters one by one to the "best" plane
		TilePlanes.PlaneData [][] new_planes = createPlanesFromSelections(
				plane_selections,   // final boolean [][][][] plane_selections, //  = new boolean [nStiles][][][]; // num_tiles
				disp_strength,       // final double  [][][][] disp_strength,
				plDispNorm,          // final double     plDispNorm,
				plMinPoints,         // final int        plMinPoints, //          =     5;  // Minimal number of points for plane detection
				plTargetEigen,       // final double     plTargetEigen, //        =   0.1;  // Remove outliers until main axis eigenvalue (possibly scaled by plDispNorm) gets below
				plFractOutliers,     // final double     plFractOutliers, //      =   0.3;  // Maximal fraction of outliers to remove
				plMaxOutliers,       // final int        plMaxOutliers, //        =    20;  // Maximal number of outliers to remove
				plPreferDisparity,   // final boolean    plPreferDisparity, // Always start with disparity-most axis (false - lowest eigenvalue)
				geometryCorrection,  // final GeometryCorrection geometryCorrection,
				correct_distortions, // final boolean    correct_distortions,

				smplMode,            // final boolean    smplMode, //        = true;   // Use sample mode (false - regular tile mode)
				mlfp,
				plFrontoTol,         // final double     plFrontoTol,         // fronto tolerance (pix) - treat almost fronto as fronto (constant disparity). <= 0 - disable
				plFrontoRms,         // final double     plFrontoRms,                       // Target rms for the fronto planes - same as sqrt(plMaxEigen) for other planes
				plFrontoOffs,        // final double     plFrontoOffs,                      // increasing weight of the near tiles by using difference between the reduced average as weight. <= 0 - disable
				PlFrontoPow,      // double       fronto_pow,    //        =   1.0;  // increase weight even more
				hor_planes,          // final boolean [][] hor_planes, // returns plane types (hor/vert)
				debugLevel+(debug_initial_discriminate? 3:0), // 1, // 0, // 1, //  + 2, // 1,          // final int        debugLevel,
				dbg_X,               // final int        dbg_X,
				dbg_Y);              // final int        dbg_Y)
		this.planes = new_planes; // save as "measured" (as opposed to "smoothed" by neighbors) planes

	}

	/**
	 * Re-assign tiles to the planes according to fitted planes. Scans each known plane parallel with specified
	 * size steps and specified number of steps, maximizing overall fit quality.
	 * Total weight of the fitted tiles? RMS? (now using number of fitted planes, if the same - than rms)
	 * @param stMeasSel - select measurements for supertiles : +1 - combo, +2 - quad +4 - hor +8 - vert
	 * @param plDispNorm - normalization of the measured disparity precision - closer objects will have disparity
	 *        difference proportionally reduced
	 * @param plMinPoints Minimal number of plane tiles to remain
	 * @param plTargetEigen target smallest eigenvalue (try to remove outliers if above)
	 * @param plFractOutliers maximal fraction of all tiles to remove as outliers
	 * @param plMaxOutliers maximal absolute number of tiles to  remove as outliers
	 * @param plPreferDisparity kept for historical reasons - when true, select disparity-most vector even if it has higher eigenvalue
	 * @param geometryCorrection GeometryCorrection instance to use
	 * @param correct_distortions correct geometrical distortions when converting to/from world coordinates
	 * @param smplMode sample mode
	 * @param smplSide sample side for averaging/filtering tile disparities
	 * @param smplNum number of best samples used fro averaging
	 * @param smplRms maximal sample disparity rms to consider sample valid
	 * @param smplWnd use window functions for the samples
	 * @param max_abs_tilt pix per tile
	 * @param max_rel_tilt pix / disparity) per tile
	 * @param damp_tilt Damp tilt to handle insufficient  (co-linear)data
	 * @param min_tilt_disp Disparity switch between filtering modes - near objects use tilts, far - use max disparity
	 * @param transition Mode transition range (between tilted and maximal disparity)
	 * @param far_mode Far objects filtering mode (0 - off, 1 - power of disparity)
	 * @param far_power Raise disparity to this power before averaging for far objects
	 * @param plDiscrTolerance maximal disparity difference from the plane to consider tile
	 * @param plDiscrDispRange parallel move known planes around original know value for the best overall fit
	 * @param plDiscrSteps number of steps (each direction) for each plane to search for the best fit (0 - single, 1 - 1 each side)
	 * @param plDiscrVariants total number of variants to try (protect from too many planes)
	 * @param plDiscrMode what plane to use as a hint: 0 - weighted, 1 - equalized, 2 - best, 3 - combined
	 * @param plDiscrVarFloor squared add to variance to calculate reverse flatness (used mostly for single-cell clusters)
	 * @param plDiscrSigma Gaussian sigma to compare how measured data is attracted to planes
	 * @param plDiscrBlur sigma to blur histograms while re-discriminating
	 * @param plDiscrExclusivity tile exclusivity: 1.0 - tile belongs to one plane only, 0.0 - regardless of others
	 * @param debugLevel debug level
	 * @param dbg_X supertile X for elevated debug level
	 * @param dbg_Y supertile X for elevated debug level
	 */

	public void regeneratePlanes(
			final TilePlanes.PlaneData [][] planes,
			final int        stMeasSel, //            = 1;      // Select measurements for supertiles : +1 - combo, +2 - quad +4 - hor +8 - vert
			final double     plDispNorm,
			final int        plMinPoints, //          =     5;  // Minimal number of points for plane detection
			final double     plTargetEigen, //        =   0.1;  // Remove outliers until main axis eigenvalue (possibly scaled by plDispNorm) gets below
			final double     plFractOutliers, //      =   0.3;  // Maximal fraction of outliers to remove
			final int        plMaxOutliers, //        =    20;  // Maximal number of outliers to remove
			final boolean    plPreferDisparity,                 // Always start with disparity-most axis (false - lowest eigenvalue)
			final double     plFrontoTol,                       // fronto tolerance (pix) - treat almost fronto as fronto (constant disparity). <= 0 - disable
			final double     plFrontoRms,                       // Target rms for the fronto planes - same as sqrt(plMaxEigen) for other planes
			final double     plFrontoOffs,                      // increasing weight of the near tiles by using difference between the reduced average as weight. <= 0 - disable
			final double     PlFrontoPow,    //        =   1.0;  // increase weight even more
			final GeometryCorrection geometryCorrection,
			final boolean    correct_distortions,

			final boolean    smplMode, //        = true;   // Use sample mode (false - regular tile mode)

			final MeasuredLayersFilterParameters mlfp,

			final double     plDiscrTolerance,     //     =   0.4;  // Maximal disparity difference from the plane to consider tile
			final double     plDiscrDispRange,     //     =   0.6;  // Parallel move known planes around original know value for the best overall fit
			final int        plDiscrSteps,         //         =   3;    // Number of steps (each direction) for each plane to search for the best fit (0 - single, 1 - 1 each side)
			final int        plDiscrMode,          //          =   3;    // What plane to use as a hint: 0 - weighted, 1 - equalized, 2 - best, 3 - combined

			final double     plDiscrVarFloor, //       =   0.03;  // Squared add to variance to calculate reverse flatness (used mostly for single-cell clusters)
			final double     plDiscrSigma, //          =   0.05;  // Gaussian sigma to compare how measured data is attracted to planes
			final double     plDiscrBlur, //           =   0.1;   // Sigma to blur histograms while re-discriminating
			final double     plDiscrExclusivity, //    =   1.5;   // Tile exclusivity: 1.0 - tile belongs to one plane only, 0.0 - regardless of others
			final double     plDiscrExclus2, //        =   0.8;   // For second pass if exclusivity > 1.0 - will assign only around strong neighbors
			final boolean    plDiscrStrict,     //         = true;   // When growing selection do not allow any offenders around (false - more these than others)
			final double     plDiscrCorrMax,    //         = 0.7;   // Attraction to different planes correlation that is too high for re-discrimination.
			final double     plDiscrCorrMerge,  //         = 0.85;  // Attraction to different planes correlation that is high enough to merge planes
			final int        plDiscrSteal,      //         = 4;     // If offender has this number of tiles (including center) the cell can not be used
			final int        plDiscrGrown,      //         = 0;     // Only use tiles within this range from original selection
			final double     plDiscrXMedian,    //         = 1.5;   // Remove outliers from the final selection that have distance more than scaled median

			final int        debugLevel,
			final int        dbg_X,
			final int        dbg_Y)
	{
		// use both horizontal and const disparity tiles to create tile clusters
		// Add max_diff (maximal disparity difference while extracting initial tile selection) and max_tries (2..3) parameters

		// Add separate method to create + remove outliers from all planes (2 different ones)?
		// TODO later re-assign pixels according to existing plane parameters
		// Sort plane data by center (plane or supertile) disparity

		boolean [][][][]  plane_selections = refineDiscriminateTiles(
				planes, // final TilePlanes.PlaneData [][] planes,
				stMeasSel,            // final int        stMeasSel,       //      = 1;      // Select measurements for supertiles : +1 - combo, +2 - quad +4 - hor +8 - vert
				plDispNorm,           // final double     plDispNorm,
				plPreferDisparity,    // final boolean    plPreferDisparity, // Always start with disparity-most axis (false - lowest eigenvalue)
				geometryCorrection,   // final GeometryCorrection geometryCorrection,
				correct_distortions,  // final boolean    correct_distortions,

				smplMode,             // final boolean    smplMode, //        = true;   // Use sample mode (false - regular tile mode)
				mlfp,

				plDiscrTolerance,     //final double     plDiscrTolerance,     //     =   0.4;  // Maximal disparity difference from the plane to consider tile
				plDiscrDispRange,     // final double     plDiscrDispRange,     //     =   0.6;  // Parallel move known planes around original know value for the best overall fit
				plDiscrSteps,         // final int        plDiscrSteps,         //         =   3;    // Number of steps (each direction) for each plane to search for the best fit (0 - single, 1 - 1 each side)
				plDiscrMode,          // final int        plDiscrMode,          //          =   3;    // What plane to use as a hint: 0 - weighted, 1 - equalized, 2 - best, 3 - combined

				plDiscrVarFloor,   // final double     plDiscrVarFloor, //       =   0.03;  // Squared add to variance to calculate reverse flatness (used mostly for single-cell clusters)
				plDiscrSigma,      // final double     plDiscrSigma, //          =   0.05;  // Gaussian sigma to compare how measured data is attracted to planes
				plDiscrBlur,       // final double     plDiscrBlur, //           =   0.1;   // Sigma to blur histograms while re-discriminating
				plDiscrExclusivity,// final double     plDiscrExclusivity, //    =   0.5;   // Tile exclusivity: 1.0 - tile belongs to one plane only, 0.0 - regardless of others
				plDiscrExclus2,    // final double     plDiscrExclus2, //    =   0.5;   // For second pass if exclusivity > 1.0 - will assign only around strong neighbors
			    plDiscrStrict,     // final boolean    plDiscrStrict, //         = true;   // When growing selection do not allow any offenders around (false - more these than others)
				plDiscrCorrMax,    // final double     plDiscrCorrMax, //         = 0.7;   // Attraction to different planes correlation that is too high for re-discrimination.
				plDiscrCorrMerge,  // final double     plDiscrCorrMerge,  //     = 0.85;  // Attraction to different planes correlation that is high enough to merge planes
				plDiscrSteal,      // final int        plDiscrSteal,       //         =   4;     // If offender has this number of tiles (including center) the cell can not be used
				plDiscrGrown,      // final int        plDiscrGrown,       //         =   0;     // Only use tiles within this range from original selection

				plDiscrXMedian,    //final double    plDiscrXMedian, //         = 1.5;   // Remove outliers from the final selection that have distance more than scaled median

				debugLevel,           // final int        debugLevel,
				dbg_X,                // final int        dbg_X,
				dbg_Y);               // final int        dbg_Y)

		// get per-tile disparity strength again (may consider using non-filtered data here) ?

		double [][][][] disp_strength = getPlaneDispStrengthsST( // here use actual disparity, not tilted
				null,                // final double []  world_plane_norm, // null - floating, (0,0,1) - constant disparity,  (0,1,0) for horizontal planes
				stMeasSel,           // final int        stMeasSel, //            = 1;      // Select measurements for supertiles : +1 - combo, +2 - quad +4 - hor +8 - vert

				plPreferDisparity,   // final boolean    plPreferDisparity, // Always start with disparity-most axis (false - lowest eigenvalue)
				geometryCorrection,  // final GeometryCorrection geometryCorrection,
				correct_distortions, // final boolean    correct_distortions,

				smplMode,            // final boolean    smplMode, //        = true;   // Use sample mode (false - regular tile mode)
				mlfp,
				debugLevel,          // final int        debugLevel,
				dbg_X,               // final int        dbg_X,
				dbg_Y);              // final int        dbg_Y)
        // Create plane data (ellipsoids) from the disparity/strength data and tile selections,
		// remove outliers, order each supertile result in ascending disparity (from far to near)
		// TODO: consider 1) assigning of the non-assigned tiles to clusters and 2) re-assigning all clusters one by one to the "best" plane
		TilePlanes.PlaneData [][] new_planes = createPlanesFromSelections(
				plane_selections,   // final boolean [][][][] plane_selections, //  = new boolean [nStiles][][][]; // num_tiles
				disp_strength,       // final double  [][][][] disp_strength,
				plDispNorm,          // final double     plDispNorm,
				plMinPoints,         // final int        plMinPoints, //          =     5;  // Minimal number of points for plane detection
				plTargetEigen,       // final double     plTargetEigen, //        =   0.1;  // Remove outliers until main axis eigenvalue (possibly scaled by plDispNorm) gets below
				plFractOutliers,     // final double     plFractOutliers, //      =   0.3;  // Maximal fraction of outliers to remove
				plMaxOutliers,       // final int        plMaxOutliers, //        =    20;  // Maximal number of outliers to remove
				plPreferDisparity,   // final boolean    plPreferDisparity, // Always start with disparity-most axis (false - lowest eigenvalue)
				geometryCorrection,  // final GeometryCorrection geometryCorrection,
				correct_distortions, // final boolean    correct_distortions,

				smplMode,            // final boolean    smplMode, //        = true;   // Use sample mode (false - regular tile mode)
				mlfp,
				plFrontoTol,         // fronto tolerance (pix) - treat almost fronto as fronto (constant disparity). <= 0 - disable
				plFrontoRms,         // final double     plFrontoRms,                       // Target rms for the fronto planes - same as sqrt(plMaxEigen) for other planes
				plFrontoOffs,        // final double     plFrontoOffs,                      // increasing weight of the near tiles by using difference between the reduced average as weight. <= 0
				PlFrontoPow,         // double       fronto_pow,    //        =   1.0;  // increase weight even more
				null,                // final boolean [][] hor_planes, // plane types (hor/vert)
				debugLevel, //  + 2, // 1,          // final int        debugLevel,
				dbg_X,               // final int        dbg_X,
				dbg_Y);              // final int        dbg_Y)
		// combine old and new planes (refineDiscriminateTiles will return null for the supertile if failed to re-disciminate)

		for (int nsTile = 0; nsTile < this.planes.length; nsTile++){
			this.planes[nsTile] = (new_planes[nsTile] != null) ? new_planes[nsTile] : planes[nsTile];
		}
	}



	public int [] getShowPlanesWidthHeight()
	{
		final int tilesX =        tileProcessor.getTilesX();
		final int tilesY =        tileProcessor.getTilesY();
		final int superTileSize = tileProcessor.getSuperTileSize();
		//				final int tileSize =      tileProcessor.getTileSize();
		final int stilesX = (tilesX + superTileSize -1)/superTileSize;
		final int stilesY = (tilesY + superTileSize -1)/superTileSize;
		int [] wh = {
				stilesX * superTileSize,// * tileSize,
				stilesY * superTileSize};// * tileSize};
		return wh;
	}

	public int fillSquares()
	{
		final int tilesX =        tileProcessor.getTilesX();
		final int tilesY =        tileProcessor.getTilesY();
		final int superTileSize = tileProcessor.getSuperTileSize();
		//				final int tileSize =      tileProcessor.getTileSize();
		final int stilesX =       (tilesX + superTileSize -1)/superTileSize;
		final int stilesY =       (tilesY + superTileSize -1)/superTileSize;
//		final int nStiles =       stilesX * stilesY;
		int num_added = 0;
		for (int stY = 0; stY < (stilesY - 1); stY++ ) {
			for (int stX = 0; stX < (stilesX - 1); stX++ ) {
				int nsTile = stY * stilesX + stX;
				TilePlanes.PlaneData [][] planes4 = {
						planes[nsTile],
						planes[nsTile + 1],
						planes[nsTile + stilesX],
						planes[nsTile + stilesX + 1]};
				if ((planes4[0] != null) && (planes4[1] != null) && (planes4[2] != null) && (planes4[3] != null)){
					for (int np = 0; np < planes4[0].length; np++){
						if (planes4[0][np] == null) continue;
						int [] neibs0 = planes4[0][np].getNeibBest();
						if ((neibs0[2] < 0) || (neibs0[4] < 0)) continue;
						if (planes4[1][neibs0[2]] == null) continue;
						if (planes4[2][neibs0[4]] == null) continue;
						int [] neibs1 = planes4[1][neibs0[2]].getNeibBest();
						int [] neibs2 = planes4[2][neibs0[4]].getNeibBest();
						if ((neibs1[4] < 0) || (neibs2[2] < 0) || (neibs1[4] != neibs2[2])) continue;
						int [] neibs3 = planes4[3][neibs1[4]].getNeibBest();
						// got full square, does it already have diagonals
						if (neibs0[3] < 0) {
							neibs0[3] = neibs1[4];
							neibs3[7] = np;
							num_added++;
						}
						if (neibs1[5] < 0){
							neibs1[5] = neibs0[4];
							neibs2[1] = neibs0[2];
							num_added++;
						}
					}
				}
			}
		}
		return num_added;
	}

	public int cutCorners()
	{
		final int tilesX =        tileProcessor.getTilesX();
		final int tilesY =        tileProcessor.getTilesY();
		final int superTileSize = tileProcessor.getSuperTileSize();
		//				final int tileSize =      tileProcessor.getTileSize();
		final int stilesX =       (tilesX + superTileSize -1)/superTileSize;
		final int stilesY =       (tilesY + superTileSize -1)/superTileSize;
//		final int nStiles =       stilesX * stilesY;
		int num_added = 0;
		for (int stY = 0; stY < (stilesY - 1); stY++ ) {
			for (int stX = 0; stX < (stilesX - 1); stX++ ) {
				int nsTile = stY * stilesX + stX;
				TilePlanes.PlaneData [][] planes4 = {
						planes[nsTile],
						planes[nsTile + 1],
						planes[nsTile + stilesX],
						planes[nsTile + stilesX + 1]};
				if (planes4[0] != null) { // && (planes4[1] != null) && (planes4[2] != null)){
					for (int np = 0; np < planes4[0].length; np++){
						if (planes4[0][np] != null) {
							int [] neibs0 = planes4[0][np].getNeibBest();

							if ((neibs0[2] >= 0) && (neibs0[4] < 0)){
								if (planes4[1][neibs0[2]] != null) {
									int [] neibs1 = planes4[1][neibs0[2]].getNeibBest();
									if (neibs1[5] >= 0) {
										neibs0[4] = neibs1[5];
										int [] neibs2 = planes4[2][neibs1[5]].getNeibBest();
										neibs2[0] = np;  //?
										num_added++;
									}
								}
							}

							if ((neibs0[2] < 0) && (neibs0[4] >= 0)){
								if (planes4[2][neibs0[4]] != null) {
									int [] neibs2 = planes4[2][neibs0[4]].getNeibBest();
									if (neibs2[1] >= 0) {
										neibs0[2] = neibs2[1];
										int [] neibs1 = planes4[1][neibs2[1]].getNeibBest();
										neibs1[6] = np;
										num_added++;
									}
								}
							}
							if ((neibs0[2] >= 0) && (neibs0[3] >= 0)){
								int [] neibs1 = planes4[1][neibs0[2]].getNeibBest();
								int [] neibs3 = planes4[3][neibs0[3]].getNeibBest();
								if (neibs1[4] < 0) {
									neibs1[4] = neibs0[3];
									neibs3[0] = neibs0[2];  //?
									num_added++;
								}
							}
							if ((neibs0[3] >= 0) && (neibs0[4] >= 0)){
								int [] neibs3 = planes4[3][neibs0[3]].getNeibBest();
								int [] neibs2 = planes4[2][neibs0[4]].getNeibBest();
								if (neibs2[2] < 0) {
									neibs2[2] = neibs0[3];
									neibs3[6] = neibs0[4];
									num_added++;
								}
							}
						}
					}
				}
				if (planes4[3] != null) { // && (planes4[1] != null) && (planes4[2] != null)){
					for (int np = 0; np < planes4[3].length; np++){
						if (planes4[3][np] != null){
							int [] neibs3 = planes4[3][np].getNeibBest();

							if ((neibs3[0] >= 0) && (neibs3[7] >= 0)){
								int [] neibs0 = planes4[0][neibs3[7]].getNeibBest();
								int [] neibs1 = planes4[1][neibs3[0]].getNeibBest();
								if (neibs0[2] < 0) {
									neibs0[2] = neibs3[0];
									neibs1[6] = neibs3[7];
									num_added++;
								}
							}


							if ((neibs3[6] >= 0) && (neibs3[7] >= 0)){
								int [] neibs0 = planes4[0][neibs3[7]].getNeibBest();
								int [] neibs2 = planes4[2][neibs3[6]].getNeibBest();
								if (neibs0[4] < 0) {
									neibs0[4] = neibs3[6];
									neibs2[0] = neibs3[7];  //?
									num_added++;
								}
							}


							if ((neibs3[7] >= 0) && (neibs3[6] < 0)){
								int [] neibs0 = planes4[0][neibs3[7]].getNeibBest();

								if (neibs0[4] >= 0) {
									int [] neibs2 = planes4[2][neibs0[4]].getNeibBest();
									neibs3[6] = neibs0[4];
									neibs2[2] = np;  //?
									num_added++;
								}
							}

							if ((neibs3[6] >= 0) && (neibs3[0] < 0)){
								int [] neibs2 = planes4[2][neibs3[6]].getNeibBest();

								if (neibs2[1] >= 0) {
									int [] neibs1 = planes4[1][neibs2[1]].getNeibBest();
									neibs3[0] = neibs2[1]; //?
									neibs1[4] = np;
									num_added++;
								}
							}
						}
					}
				}
			}
		}
		return num_added;

	}

	public int [] resolveDiagonalTriangularConflicts(
			int [][][] conflicts,
			Conflicts conflict_stats, // to be updated after applying resolution
			int        starSteps, // How far to look around when calculationg connection cost
			double     orthoWeight,
			double     diagonalWeight,
			double     starPwr, // Divide cost by number of connections to this power
			double     starWeightPwr,    // Use this power of tile weight when calculating connection cost
			double     weightToDens,    // // Balance weighted density against density. 0.0 - density, 1.0 - weighted density
			double     starValPwr, //  Raise value of each tile before averaging
			double     dblTriLoss, //  When resolving double triangles allow minor degradation (0.0 - strict)
			boolean    preferDisparity,
			int        debugLevel,
			int        dbg_X,
			int        dbg_Y)
	{
		final int tilesX =        tileProcessor.getTilesX();
		final int tilesY =        tileProcessor.getTilesY();
		final int superTileSize = tileProcessor.getSuperTileSize();
		final int stilesX =       (tilesX + superTileSize -1)/superTileSize;
		final int stilesY =       (tilesY + superTileSize -1)/superTileSize;
		final int dbgTile = dbg_Y * stilesX + dbg_X;
		final TileNeibs tnSurface = new TileNeibs(stilesX, stilesY);
		int [] rslt = {0,0};
		for (int nsTile = 0; nsTile < conflicts.length; nsTile++) if (conflicts[nsTile] != null) {
			// conflicts may disappear after being fixed, recheck for null
			for (int nConfl = 0; (conflicts[nsTile] != null) && (nConfl < conflicts[nsTile].length); nConfl++){
				int dl = ((debugLevel > -1) && (nsTile == dbgTile)) ? 3 : 0;
				Conflict confl = new Conflict(nsTile, conflicts[nsTile][nConfl]);
				for (int dir_orient = 0; dir_orient < 16; dir_orient++){
					boolean right = (dir_orient & 8) != 0;
					boolean reverse = (dir_orient & 4) != 0;
					int start_dir4 = dir_orient & 3;
					if (confl.conflictExists(
							start_dir4,
							right,
							reverse)){
						boolean OK = resolveDiagonalTriangularConflict(
								confl,           // Conflict conflict,
								start_dir4,      // int start_dir4,
								right,           // boolean right,
								reverse,         // boolean reverse,
								tnSurface,       // TileNeibs tnSurface,
								conflicts,       // int [][][] conflicts,
								conflict_stats,  // Conflicts conflict_stats, // to be updated after applying resolution
								starSteps, // How far to look around when calculationg connection cost
								orthoWeight,     // double     orthoWeight,
								diagonalWeight,  // double     diagonalWeight,
								starPwr,         // double     starPwr, // Divide cost by number of connections to this power
								starWeightPwr,    // Use this power of tile weight when calculating connection cost
								weightToDens,    // // Balance weighted density against density. 0.0 - density, 1.0 - weighted density
								starValPwr,      //double     starValPwr, //  Raise value of each tile before averaging
								dblTriLoss,      //  double     dblTriLoss, //  When resolving double triangles allow minor degradation (0.0 - strict)
								preferDisparity, // boolean    preferDisparity,
								dl);             // int        debugLevel)
						if (OK) rslt[0]++; // will not count if some conflicts disappeared as a result
						else    rslt[1]++;
					}
				}
			}
		}
		return rslt;
	}



	public boolean resolveDiagonalTriangularConflict(

			Conflict conflict,
			int start_dir4,
			boolean right,
			boolean reverse,
			TileNeibs tnSurface,
			int [][][] conflicts,
			Conflicts conflict_stats, // to be updated after applying resolution
			int        starSteps, // How far to look around when calculationg connection cost
			double     orthoWeight,
			double     diagonalWeight,
			double     starPwr, // Divide cost by number of connections to this power
			double     starWeightPwr,    // Use this power of tile weight when calculating connection cost
			double     weightToDens,    // // Balance weighted density against density. 0.0 - density, 1.0 - weighted density
			double     starValPwr, //  Raise value of each tile before averaging
			double     dblTriLoss, //  When resolving double triangles allow minor degradation (0.0 - strict)
			boolean    preferDisparity,
			int        debugLevel)
	{
		if (debugLevel > 1) {
			System.out.println("resolveDiagonalTriangularConflict(): "+conflict.toString()+": start_dir4 = "+ start_dir4 + ", right = "+right+", reverse = "+reverse);
		}

		if (!conflict.conflictExists(
				start_dir4,
				right,
				reverse)) {
			System.out.println("conflict does not exist");
			return false; // should be tested before to distinguish between non-existent and failure to resolve
		}
		Conflicts iconflicts = new Conflicts(this);
		int dir01 = 2 * start_dir4;
		int dir12 = 2 * ((start_dir4 + (right? 1 : 3)) % 4);
		int dir20 = (dir12 + (right? 3 : 5)) % 8;
		int dir02 = (dir20 + 4) % 8;


		int nsTile = conflict.getSTile();
		int nl = conflict.getStartLayer(reverse);
		int nl_end = conflict.getEndLayer(reverse);
		int nsTile1 = tnSurface.getNeibIndex(nsTile, dir01);
		int nl1 = planes[nsTile][nl].getNeibBest(dir01);
		int nsTile2 = tnSurface.getNeibIndex(nsTile1, dir12);
		if (nl1 <0){
			System.out.println(conflict.toString()+" : start_dir4 = "+ start_dir4 + ", right = "+right+", reverse = "+reverse+" nl1 < 0 !");
			return false;
		}

		if (debugLevel > 1) {
			if (debugLevel > 1) {
				System.out.println("resolveDiagonalTriangularConflict(): nsTile = "+nsTile+" nsTile1 = "+nsTile1+" nsTile2="+nsTile2);
			}
		}

		int nl2 = planes[nsTile1][nl1].getNeibBest(dir12);
		// now nsTile2:nl2 is supposed to be connected to nsTile:nl_end != nl

		int [] mod_supertiles = {nsTile,nsTile2};
		int [] nsTiles = getInvolvedSupertiles( // first mod_supertiles.length entries will be mod_supertiles[]
				mod_supertiles,
				tnSurface);
		if (debugLevel > 1) {
			System.out.println("Involved supertiles:");
			for (int i = 0; i < nsTiles.length; i++){
				System.out.println(i+":"+nsTiles[i]);
			}
		}
		HashMap<Integer,Integer> replacement_tiles = new HashMap<Integer,Integer>();
		for (int i = 0; i < mod_supertiles.length; i++){
			replacement_tiles.put(mod_supertiles[i], new Integer(i));
		}
		ConnectionCosts connectionCosts = new ConnectionCosts(
				orthoWeight,
				diagonalWeight,
				starPwr,    // Divide cost by number of connections to this power
				starWeightPwr,    // Use this power of tile weight when calculating connection cost
				weightToDens,    // // Balance weighted density against density. 0.0 - density, 1.0 - weighted density
				starValPwr,      //double     starValPwr, //  Raise value of each tile before averaging
				starSteps,
				this.planes,
				tnSurface,
				preferDisparity);

		int [][][] neibs_prev = connectionCosts.initConnectionCosts(mod_supertiles, debugLevel);

		int [][][] neibs_prev_old = new int [mod_supertiles.length][][];
		double [][][]  val_weights = new double [mod_supertiles.length][][];

		// Calculate original costs and neighhbors
		updateConnectionsCost_old (
				mod_supertiles,      // int []         nsTiles,
				null,         // int [][][]     neibs_prev,
				neibs_prev_old,   // int [][][]     neibs,		// should be initialized at top dimension if neibs_prev==null
				val_weights,  // double [][][]  val_weights, // should be initialized at top dimension if neibs_prev==null
				orthoWeight,
				diagonalWeight,
				starPwr,        // double     starPwr, // Divide cost by number of connections to this power
				starWeightPwr,    // Use this power of tile weight when calculating connection cost
				tnSurface,
				preferDisparity,
				debugLevel);

		if ((neibs_prev == null) || (neibs_prev[0] == null) || (neibs_prev[1] == null)){
			System.out.println("resolveDiagonalTriangularConflict(): BUG");
			System.out.println("neibs_prev="+neibs_prev);
			System.out.println("neibs_prev[0]="+neibs_prev[0]);
			System.out.println("neibs_prev[1]="+neibs_prev[1]);
			return false;
		}
		int [][][] neibs = {neibs_prev[0].clone(),neibs_prev[1].clone()}; // both non-null
		for (int nt = 0; nt < neibs.length; nt++){
			for (int l = 0; l < neibs_prev[nt].length; l++) if (neibs_prev[nt][l]!= null){
				neibs[nt][l] = neibs_prev[nt][l].clone();
			}
		}
		int [][][] neibs_old = {neibs_prev[0].clone(),neibs_prev[1].clone()}; // both non-null
		for (int nt = 0; nt < neibs.length; nt++){
			for (int l = 0; l < neibs_prev[nt].length; l++) if (neibs_prev[nt][l]!= null){
				neibs_old[nt][l] = neibs_prev_old[nt][l].clone();
			}
		}
		int [][][] conflicts_old = new int [nsTiles.length][][];
		for (int isTile = 0; isTile < nsTiles.length; isTile++){
			conflicts_old[isTile] = iconflicts.detectTriangularTileConflicts(
					nsTiles[isTile],   // int nsTile0,
					replacement_tiles, //HashMap<Integer,Integer> replacement_tiles, //
					neibs_prev, // int [][][] replacement_neibs,
					tnSurface); // TileNeibs tnSurface)
		}
		if (debugLevel > 1) {
			System.out.println("Calculating original conflicts");
		}
		iconflicts.addConflicts(conflicts_old,
				debugLevel - 1); // debugLevel);

		if (nl2 <0){
			System.out.println(conflict.toString()+" : start_dir4 = "+ start_dir4 + ", right = "+right+", reverse = "+reverse+" nl2 < 0 !");
			return false;
		}

		// now swap diagonal connection
		neibs[0][nl][dir02] = nl2;
		neibs[1][nl2][dir20] = nl;
		neibs[0][nl_end][dir02] = neibs_prev[0][nl][dir02];
		if (neibs_prev[0][nl][dir02] >= 0) {
			neibs[1][neibs_prev[0][nl][dir02]][dir20] = nl_end;
		}

		neibs_old[0][nl][dir02] = nl2;
		neibs_old[1][nl2][dir20] = nl;
		neibs_old[0][nl_end][dir02] = neibs_prev_old[0][nl][dir02];
		if (neibs_prev_old[0][nl][dir02] >= 0) {
			neibs_old[1][neibs_prev_old[0][nl][dir02]][dir20] = nl_end;
		}



		int [][][] new_conflicts = new int [nsTiles.length][][];

		double new_costs_diff = 	connectionCosts.getConnectionsCostDiff(
				neibs,
				debugLevel);

		if (debugLevel > -1) {
			System.out.println("resolveDiagonalTriangularConflict(): resolving conflict "+conflict.toString()+
					", start dir4 = "+start_dir4+
					", orientation = "+(right?"right":"left") +
					", improvement (negative diff) = "+new_costs_diff);

		}
		double new_costs_diff_old = 	updateConnectionsCost_old (
				mod_supertiles,      // int []         nsTiles,
				neibs_prev_old,          // int [][][]     neibs_prev,
				neibs_old,               // int [][][]     neibs,		// should be initialized at top dimension if neibs_prev==null
				val_weights,         // double [][][]  val_weights, // should be initialized at top dimension if neibs_prev==null
				orthoWeight,
				diagonalWeight,
				starPwr,        // double     starPwr, // Divide cost by number of connections to this power
				starWeightPwr,    // Use this power of tile weight when calculating connection cost
				tnSurface,
				preferDisparity,
				debugLevel);

		for (int isTile = 0; isTile < nsTiles.length; isTile++){
			new_conflicts[isTile] = iconflicts.detectTriangularTileConflicts(
					nsTiles[isTile],   // int nsTile0,
					replacement_tiles, //HashMap<Integer,Integer> replacement_tiles, //
					neibs, // neibs, // int [][][] replacement_neibs,
					tnSurface); // TileNeibs tnSurface)
		}
		Conflicts new_conflicts_stats = new Conflicts(
				new_conflicts,
				this,
				debugLevel - 1); // debugLevel);
		new_conflicts_stats.subConflicts(iconflicts); // subtract old number of different types of conflicts

		if (debugLevel > -1) {
			System.out.println("resolveDiagonalTriangularConflict() OLD: resolving conflict "+conflict.toString()+
					", start dir4 = "+start_dir4+
					", orientation = "+(right?"right":"left") +
					", improvement (negative diff) = "+new_costs_diff_old);

		}
		if (debugLevel > 0) { // -1) {
			new_conflicts_stats.printConflictSummary(
					"Conflicts difference after resolution:", true, true, true);
		}
		// See if it is good
		if (debugLevel > 0) { // -1) {
			if (new_conflicts_stats.getNumConflicts() > 0) {
				if (debugLevel > -1) System.out.println("FAILURE: number of conflicts increased");
			} else if (new_costs_diff > dblTriLoss){
				if (debugLevel > -1) System.out.println("FAILURE: cost increased too much");
			} else if ((new_costs_diff >= 0.0) && (new_conflicts_stats.getNumConflicts() == 0)){
				if (debugLevel > -1) System.out.println("FAILURE: cost increased, but conflicts number is the same");
			} else if (new_costs_diff >= 0.0){
				if (debugLevel > -1) System.out.println("SUCCESS: cost did not decrease, but conflicts did");
			} else if (new_conflicts_stats.getNumConflicts() == 0){
				if (debugLevel > -1) System.out.println("SUCCESS: conflicts number is the same, but connection cost decreased");
			} else {
				if (debugLevel > -1) System.out.println("SUCCESS: both conflicts and connection cost decreased");
			}
			if (debugLevel > 1) {
				System.out.println("resolveDiagonalTriangularConflict(): nsTile = "+conflict.getSTile());
			}
		}

		if (new_conflicts_stats.getNumConflicts() > 0) {
			return false;
		}
		if (new_costs_diff > dblTriLoss){
			return false;
		}
		if ((new_costs_diff >= 0.0) && (new_conflicts_stats.getNumConflicts() == 0)){
			return false;
		}
		boolean apply = true;
		// apply changes
		if (apply) {
			// update statistics
			conflict_stats.addConflicts(new_conflicts_stats);
			// update conflict
			for (int i = 0; i < nsTiles.length; i++){
				conflicts[nsTiles[i]]=  new_conflicts[i];
			}

			// apply resolution
			for (int i = 0; i < mod_supertiles.length; i++){
				for (int l = 0; l < neibs[i].length; l ++) if (neibs[i][l] != null){
					if (planes[mod_supertiles[i]].length <= l) {
						System.out.println("BUG: (planes["+mod_supertiles[i]+"].length="+planes[mod_supertiles[i]].length+" <= "+l);
						return false;
					}
					planes[mod_supertiles[i]][l].setNeibBest(neibs[i][l]);
				}
			}
		}
		return true;
	}


	public boolean resolveMultiTriangularConflict(
			int nsTile,
			int nl1,
			int nl2,
			int dir_mask,
			TileNeibs tnSurface,
			int [][][] conflicts,
			Conflicts conflict_stats, // to be updated after applying resolution
			int        starSteps, // How far to look around when calculationg connection cost
			double     orthoWeight,
			double     diagonalWeight,
			double     starPwr, // Divide cost by number of connections to this power
			double     starWeightPwr,    // Use this power of tile weight when calculating connection cost
			double     weightToDens,    // Balance weighted density against density. 0.0 - density, 1.0 - weighted density
			double     starValPwr, //  Raise value of each tile before averaging
			double     dblTriLoss, //  When resolving double triangles allow minor degradation (0.0 - strict)
			boolean    preferDisparity,
			int        debugLevel)
	{
		Conflicts iconflicts = new Conflicts(this);
		// neibs_raw may have unused directions
		int [][][][] neibs_raw = getTriangularResolutionVariants(
				nsTile,
				nl1,
				nl2,
				dir_mask,
				tnSurface,
				debugLevel);
		if (neibs_raw.length == 0){
			System.out.println("resolveMultiTriangularConflict() - zero variants. A bug?");
			return false;
		}

		int num_tiles = 0;
		for (int i = 0; i < neibs_raw[0].length; i++){
			if (neibs_raw[0][i] != null) num_tiles++;
 		}
		int [] mod_supertiles = new int [num_tiles];
		int [] indices = new int[num_tiles];
		int indx = 0;
		indices[indx++] = 4;
		for (int i = 0; (i < neibs_raw[0].length) && (indx < indices.length); i++){
			if (neibs_raw[0][i] != null) indices[indx++] = i;
 		}
		for (int i = 0; i < indices.length; i++){
			int dir4 = indices[i];
			if (dir4 > 3) {
				mod_supertiles[i] = nsTile;
			} else {
				mod_supertiles[i] = tnSurface.getNeibIndex(nsTile, 2 * dir4);
			}
		}
		int [][][][] neibs_vars = new int [neibs_raw.length][mod_supertiles.length][][];
		for (int ivar = 0; ivar < neibs_raw.length; ivar++){
			for (int i = 0; i < indices.length; i++){
				neibs_vars[ivar][i] = neibs_raw[ivar][indices[i]];
			}
		}



		// See how this application will influence number of conflicts
		// All supertiles that may have different conflicts
		int [] nsTiles = getInvolvedSupertiles( // first mod_supertiles.length entries will be mod_supertiles[]
				mod_supertiles,
				tnSurface);
		HashMap<Integer,Integer> replacement_tiles = new HashMap<Integer,Integer>();
		for (int i = 0; i < mod_supertiles.length; i++){
			replacement_tiles.put(mod_supertiles[i], new Integer(i));
		}

		ConnectionCosts connectionCosts = new ConnectionCosts(
				orthoWeight,
				diagonalWeight,
				starPwr,    // Divide cost by number of connections to this power
				starWeightPwr,    // Use this power of tile weight when calculating connection cost
				weightToDens,    // Balance weighted density against density. 0.0 - density, 1.0 - weighted density
				starValPwr,      //double     starValPwr, //  Raise value of each tile before averaging
				starSteps,
				this.planes,
				tnSurface,
				preferDisparity);

		int [][][] neibs_prev = connectionCosts.initConnectionCosts(mod_supertiles, debugLevel);

/** */
		int [][][] neibs_prev_old = new int [mod_supertiles.length][][];
		double [][][]  val_weights_original = new double [mod_supertiles.length][][];
		updateConnectionsCost_old (
				mod_supertiles,      // int []         nsTiles,
				null,         // int [][][]     neibs_prev,
				neibs_prev_old,   // int [][][]     neibs,		// should be initialized at top dimension if neibs_prev==null
				val_weights_original,  // double [][][]  val_weights, // should be initialized at top dimension if neibs_prev==null
				orthoWeight,
				diagonalWeight,
				starPwr,        // double     starPwr, // Divide cost by number of connections to this power
				starWeightPwr,    // Use this power of tile weight when calculating connection cost
				tnSurface,
				preferDisparity,
				debugLevel);
/** */
		int [][][] conflicts_old = new int [nsTiles.length][][];
		for (int isTile = 0; isTile < nsTiles.length; isTile++){
			conflicts_old[isTile] = iconflicts.detectTriangularTileConflicts(
					nsTiles[isTile],   // int nsTile0,
					replacement_tiles, //HashMap<Integer,Integer> replacement_tiles, //
					neibs_prev, // int [][][] replacement_neibs,
					tnSurface); // TileNeibs tnSurface)
		}
		if (debugLevel > 1) {
			System.out.println("Involved supertiles:");
			for (int i = 0; i < nsTiles.length; i++){
				System.out.println(i+":"+nsTiles[i]);
			}
		}

		if (debugLevel > 1) {
			System.out.println("Calculating original conflicts");
		}

		iconflicts.addConflicts(conflicts_old,
				debugLevel - 1); // debugLevel);
		// After getting old (referfence data) iterate through all variants, for each calculate cost and number of conflicts.
		// First - just collect data - cost and full statistics of the conflicts, print them
		// Then select the best one (not incrementing number of conflicts? and reducing cost)
		int [][][][] variant_conflicts = new int [neibs_vars.length][nsTiles.length][][];

		Conflicts [] variant_conflicts_stats = new Conflicts [neibs_vars.length];


		double [] variant_costs_diff = new double [neibs_vars.length];
		for (int variant = 0; variant < neibs_vars.length; variant ++){
			if (debugLevel > 0) {
				System.out.println("resolveMultiTriangularConflict(): resolving conflict for tile "+nsTile+
						", nl1 = "+nl1+
						", nl2 = "+nl2 +
						", dir_mask = "+dir_mask+" variant = "+variant);

			}
			variant_costs_diff[variant] = 	connectionCosts.getConnectionsCostDiff(
					neibs_vars[variant],
					debugLevel);

			if (debugLevel > -1) {
				System.out.println("resolveMultiTriangularConflict(): resolving conflict for tile "+nsTile+
						", nl1 = "+nl1+
						", nl2 = "+nl2 +
						", dir_mask = "+dir_mask+" variant = "+variant+" improvement (negative diff) = "+variant_costs_diff[variant]);

			}

			double [] variant_costs_diff_old = new double [neibs_vars.length];
			double [][][]  val_weights = val_weights_original.clone();
			for (int nt = 0; nt < val_weights_original.length; nt++) if (val_weights_original[nt] != null){
				val_weights[nt] = val_weights_original[nt].clone();
				for (int nl = 0; nl < val_weights_original[nt].length; nl++) if (val_weights_original[nt][nl] != null){
					val_weights[nt][nl] = val_weights_original[nt][nl].clone();
				}
			}
			variant_costs_diff_old[variant] = 	updateConnectionsCost_old (
					mod_supertiles,      // int []         nsTiles,
					neibs_prev_old,          // int [][][]     neibs_prev,
					neibs_vars[variant], // int [][][]     neibs,		// should be initialized at top dimension if neibs_prev==null
					val_weights,         // double [][][]  val_weights, // should be initialized at top dimension if neibs_prev==null
					orthoWeight,
					diagonalWeight,
					starPwr,        // double     starPwr, // Divide cost by number of connections to this power
					starWeightPwr,    // Use this power of tile weight when calculating connection cost
					tnSurface,
					preferDisparity,
					debugLevel);



			for (int isTile = 0; isTile < nsTiles.length; isTile++){
				variant_conflicts[variant][isTile] = iconflicts.detectTriangularTileConflicts(
						nsTiles[isTile],   // int nsTile0,
						replacement_tiles, //HashMap<Integer,Integer> replacement_tiles, //
						neibs_vars[variant], // neibs, // int [][][] replacement_neibs,
						tnSurface); // TileNeibs tnSurface)
			}
			variant_conflicts_stats[variant] = new Conflicts(
					variant_conflicts[variant],
					this,
					debugLevel - 1); // debugLevel);
			variant_conflicts_stats[variant].subConflicts(iconflicts); // subtract old number of different types of conflicts

			if (debugLevel > -1) {
				System.out.println("resolveMultiTriangularConflict() OLD: resolving conflict for tile "+nsTile+
						", nl1 = "+nl1+
						", nl2 = "+nl2 +
						", dir_mask = "+dir_mask+" variant = "+variant+" improvement (negative diff) = "+variant_costs_diff_old[variant]);

			}
			if (debugLevel > -1) {
				variant_conflicts_stats[variant].printConflictSummary(
						"Conflicts difference after resolution:", true, true, false);
			}
		}
		// How to compare? 1 attempt: none of the conflicts get worse, some get better or cost goes down
		int best_variant = -1;
		int [][] num_better_worse = new int [neibs_vars.length][2];
		for (int variant = 0; variant < neibs_vars.length; variant ++){
			int num_worse = variant_conflicts_stats[variant].numBetterWorse(
					false, // boolean better,
					false, // boolean use_all,
					true, // boolean use_odo,
					false); // ); // boolean use_ood)

			int num_better = variant_conflicts_stats[variant].numBetterWorse(
					true, // boolean better,
					false, // boolean use_all,
					true, // boolean use_odo,
					false); // ); // boolean use_ood)
			num_better_worse[variant][0] = num_better;
			num_better_worse[variant][1] = num_worse;
			if ((num_worse == 0) &&
					(variant_costs_diff[variant] <= dblTriLoss) && // not too worse
					((variant_costs_diff[variant] < 0) || (num_better > 0)) && // either
					((best_variant < 0) ||
							( num_better_worse[variant][0] > num_better_worse[best_variant][0]) ||
							(( num_better_worse[variant][0] == num_better_worse[best_variant][0]) &&
									(variant_costs_diff[variant] < variant_costs_diff[best_variant])))){
				best_variant = variant;
			}
		}

		if (debugLevel > 1){
			System.out.println("resolveMultiTriangularConflict(): for tile "+nsTile);
		}

		if ((best_variant < 0) || (variant_costs_diff[best_variant] > dblTriLoss)){
			if (debugLevel > -1) {
				System.out.println("resolveMultiTriangularConflict(): FAILED find a sutable solution for tile "+nsTile+
						", nl1 = "+nl1+
						", nl2 = "+nl2 +
						", dir_mask = "+dir_mask+" of "+ neibs_vars.length+" variants");
				return false;
			}
		} else {
			if (debugLevel > -1) {
				System.out.println("resolveMultiTriangularConflict(): SUCCESS to find a sutable solution for tile "+nsTile+
						", nl1 = "+nl1+
						", nl2 = "+nl2 +
						", dir_mask = "+dir_mask+". Of "+ neibs_vars.length+" variants - use variant # " + best_variant+
						" cost difference (negative) = "+variant_costs_diff[best_variant]+" num conflict reductions = "+num_better_worse[best_variant][0]);
				variant_conflicts_stats[best_variant].printConflictSummary(
						"Conflicts number change per type: ",
						true, // use_all,
						true, //use_odo,
						true); // use_ood);
				iconflicts.printConflictSummary(
						"Conflicts before resolution: ",
						true, // use_all,
						true, //use_odo,
						true); // use_ood);
				// update statistics
				conflict_stats.addConflicts(variant_conflicts_stats[best_variant]);


				// update conflict
				for (int i = 0; i < nsTiles.length; i++){
					conflicts[nsTiles[i]]=  variant_conflicts[best_variant][i];
				}


				// apply resolution
				for (int i = 0; i < mod_supertiles.length; i++){
					for (int nl = 0; nl < neibs_vars[best_variant][i].length; nl ++) if (neibs_vars[best_variant][i][nl] != null){
						planes[mod_supertiles[i]][nl].setNeibBest(neibs_vars[best_variant][i][nl]);
					}
				}
			}
		}
		return true;
	}


	public int [] resolveMultiTriangularConflicts(
			int [][][] conflicts,
			Conflicts conflict_stats, // to be updated after applying resolution
			int        starSteps, // How far to look around when calculationg connection cost
			double     orthoWeight,
			double     diagonalWeight,
			double     starPwr, // Divide cost by number of connections to this power
			double     starWeightPwr,    // Use this power of tile weight when calculating connection cost
			double     weightToDens,    // Balance weighted density against density. 0.0 - density, 1.0 - weighted density
			double     starValPwr, //  Raise value of each tile before averaging
			double     dblTriLoss, //  When resolving double triangles allow minor degradation (0.0 - strict)
			boolean    preferDisparity,
			int        debugLevel,
			int        dbg_X,
			int        dbg_Y)
	{
		final int tilesX =        tileProcessor.getTilesX();
		final int tilesY =        tileProcessor.getTilesY();
		final int superTileSize = tileProcessor.getSuperTileSize();
		//				final int tileSize =      tileProcessor.getTileSize();
		final int stilesX =       (tilesX + superTileSize -1)/superTileSize;
		final int stilesY =       (tilesY + superTileSize -1)/superTileSize;
		final int dbgTile = dbg_Y * stilesX + dbg_X;
		final TileNeibs tnSurface = new TileNeibs(stilesX, stilesY);
		int [] rslt = {0,0};
		for (int nsTile = 0; nsTile < conflicts.length; nsTile++) if (conflicts[nsTile] != null) {
			// conflicts may disappear after being fixed, recheck for null
			for (int nConfl = 0; (conflicts[nsTile] != null) && (nConfl < conflicts[nsTile].length); nConfl++){
				int dl = ((debugLevel > 0) && (nsTile == dbgTile)) ? 3 : 0;
				boolean OK = resolveMultiTriangularConflict(
						nsTile,
						conflicts[nsTile][nConfl][0], // int nl1,
						conflicts[nsTile][nConfl][1], // int nl2,
						conflicts[nsTile][nConfl][2], // int dir_mask,
						tnSurface,
						conflicts,
						conflict_stats, // to be updated after applying resolution
//						maxEigen, // maximal eigenvalue of planes to consider
						starSteps, // How far to look around when calculationg connection cost
						orthoWeight,
						diagonalWeight,
						starPwr,        // double     starPwr, // Divide cost by number of connections to this power
						starWeightPwr,    // Use this power of tile weight when calculating connection cost
						weightToDens,    // Balance weighted density against density. 0.0 - density, 1.0 - weighted density
						starValPwr,      //double     starValPwr, //  Raise value of each tile before averaging
						dblTriLoss, //  When resolving double triangles allow minor degradation (0.0 - strict)
						preferDisparity,
						dl); // debugLevel,
				if (OK) rslt[0]++;
				else    rslt[1]++;
			}
		}
		return rslt;
	}


	public int [] resolveStarConflicts(
			int [][][] conflicts,
			Conflicts conflict_stats, // to be updated after applying resolution
			LinkPlanes lp,
			int        starSteps, // How far to look around when calculationg connection cost
			double     orthoWeight,
			double     diagonalWeight,
			double     starPwr, // Divide cost by number of connections to this power
			double     starWeightPwr,    // Use this power of tile weight when calculating connection cost
			double     weightToDens,    // Balance weighted density against density. 0.0 - density, 1.0 - weighted density
			double     starValPwr, //  Raise value of each tile before averaging
			double     dblTriLoss, //  When resolving double triangles allow minor degradation (0.0 - strict)
			boolean    newConfl, // Allow more conflicts if overall cost is reduced
			int        maxChanges,  // Maximal number of simultaneous connection changes around one tile (0 - any)
			boolean    preferDisparity,
			int        debugLevel,
			int        dbg_X,
			int        dbg_Y)
	{
		final int tilesX =        tileProcessor.getTilesX();
		final int tilesY =        tileProcessor.getTilesY();
		final int superTileSize = tileProcessor.getSuperTileSize();
		//				final int tileSize =      tileProcessor.getTileSize();
		final int stilesX =       (tilesX + superTileSize -1)/superTileSize;
		final int stilesY =       (tilesY + superTileSize -1)/superTileSize;
		final int dbgTile = dbg_Y * stilesX + dbg_X;
		final TileNeibs tnSurface = new TileNeibs(stilesX, stilesY);
		int [] rslt = {0,0};
		for (int nsTile = 0; nsTile < conflicts.length; nsTile++) if (conflicts[nsTile] != null) {
			// conflicts may disappear after being fixed, recheck for null
			for (int nConfl = 0; (conflicts[nsTile] != null) && (nConfl < conflicts[nsTile].length); nConfl++){
				int dl = ((debugLevel > 0) && (nsTile == dbgTile)) ? 3 : 0;
				boolean OK = resolveStarConflict(
						nsTile,
						conflicts[nsTile][nConfl][0], // int nl1,
						conflicts[nsTile][nConfl][1], // int nl2,
						lp, // LinkPlanes lp,
						tnSurface,
						conflicts,
						conflict_stats, // to be updated after applying resolution
						starSteps, // How far to look around when calculating connection cost
						orthoWeight,
						diagonalWeight,
						starPwr,        // double     starPwr, // Divide cost by number of connections to this power
						starWeightPwr,    // Use this power of tile weight when calculating connection cost
						weightToDens,    // Balance weighted density against density. 0.0 - density, 1.0 - weighted density
						starValPwr,     //double     starValPwr, //  Raise value of each tile before averaging
						dblTriLoss,     //  When resolving double triangles allow minor degradation (0.0 - strict)
						newConfl,       // boolean    newConfl, // Allow more conflicts if overall cost is reduced
						maxChanges,     // int        maxChanges,  // Maximal number of simultaneous connection changes around one tile (0 - any)
						preferDisparity,
						dl); // debugLevel,
				if (OK) rslt[0]++;
				else    rslt[1]++;
			}
		}
		return rslt;
	}


	public boolean resolveStarConflict0(
			int nsTile,
			int nl1,
			int nl2,
			TileNeibs tnSurface,
			int [][][] conflicts,
			Conflicts conflict_stats, // to be updated after applying resolution
			int        starSteps, // How far to look around when calculationg connection cost
			double     orthoWeight,
			double     diagonalWeight,
			double     starPwr, // Divide cost by number of connections to this power
			double     starWeightPwr,    // Use this power of tile weight when calculating connection cost
			double     weightToDens,    // Balance weighted density against density. 0.0 - density, 1.0 - weighted density
			double     starValPwr, //  Raise value of each tile before averaging
			double     dblTriLoss, //  When resolving double triangles allow minor degradation (0.0 - strict)
			boolean    newConfl, // Allow more conflicts if overall cost is reduced
			int        maxChanges,  // Maximal number of simultaneous connection changes around one tile (0 - any)
			boolean    preferDisparity,
			int        debugLevel)
	{

		if (newConfl && (dblTriLoss > 0.0)){
			dblTriLoss = 0.0; // require strict reducing of the cost if conflict increase is allowed
		}
		Conflicts iconflicts = new Conflicts(this);

		TwoLayerNeighbors twoLayerNeighbors = new TwoLayerNeighbors();

		for (int dir = -1; dir < 8; dir++){
			int nt = tnSurface.getNeibIndex(nsTile, dir);
			if ((nt >= 0) && (planes[nt] != null)) {
				int [][] neibs =  new int [planes[nt].length][];
				for (int nl = 0; nl < planes[nt].length; nl++) if ( planes[nt][nl] != null){
					neibs[nl] =  planes[nt][nl].getNeibBest();
				}
				twoLayerNeighbors.setNeighbors(neibs,dir);
			}
		}
		twoLayerNeighbors.setLayers(nl1, nl2);
		if (debugLevel > 1) {
			System.out.println("resolveStarConflict(): nsTile ="+nsTile+" nl1="+nl1+" nl2="+nl2);
		}

		int [][][][] neibs_vars_dir = twoLayerNeighbors.getNeighborVariants(
				maxChanges,
				debugLevel);

		int [] mod_supertiles = {nsTile};
		mod_supertiles = getInvolvedSupertiles( // first mod_supertiles.length entries will be mod_supertiles[]
				mod_supertiles,
				tnSurface);

		HashMap<Integer,Integer> replacement_tiles = new HashMap<Integer,Integer>();
		for (int i = 0; i < mod_supertiles.length; i++){
			replacement_tiles.put(mod_supertiles[i], new Integer(i));
		}

		// up to 9 tiles
		int [] indexToDir = new int[mod_supertiles.length];
		for (int i = 0; i < indexToDir.length; i++) indexToDir[i] = -1;
		for (int dir = -1; dir < 8; dir++){
			int nindx = (dir < 0) ? 8 : dir;
			int nt = tnSurface.getNeibIndex(nsTile, dir);
			int indx = -1;
			Integer Indx = replacement_tiles.get(nt);
			if (Indx == null){
				System.out.println("resolveStarConflict(): nsTile = "+nsTile+" nindx="+nindx+" Indx == null ");
			} else {
				indx = Indx;
			}
			if (indx >= 0) {
				indexToDir[indx] = nindx;
			}
		}
		int [][][][] neibs_vars = new int [neibs_vars_dir.length][][][];
		for (int variant = 0; variant < neibs_vars_dir.length; variant ++){
			neibs_vars[variant] = new int [indexToDir.length][][];
			for (int i = 0; i <indexToDir.length; i++){
				if (indexToDir[i] >=0){
					neibs_vars[variant][i] = neibs_vars_dir[variant][indexToDir[i]];
				} else {
					System.out.println("resolveStarConflict(): a BUG:  indexToDir["+i+"] = "+indexToDir[i]);
				}
			}
		}

		// See how this application will influence number of conflicts
		// All supertiles that may have different conflicts
		int [] nsTiles = getInvolvedSupertiles( // first mod_supertiles.length entries will be mod_supertiles[]
				mod_supertiles,
				tnSurface);

		ConnectionCosts connectionCosts = new ConnectionCosts(
				orthoWeight,
				diagonalWeight,
				starPwr,    // Divide cost by number of connections to this power
				starWeightPwr,    // Use this power of tile weight when calculating connection cost
				weightToDens,    // Balance weighted density against density. 0.0 - density, 1.0 - weighted density
				starValPwr,      //double     starValPwr, //  Raise value of each tile before averaging
				starSteps,
				this.planes,
				tnSurface,
				preferDisparity);

		int [][][] neibs_prev = connectionCosts.initConnectionCosts(mod_supertiles, debugLevel);

		int [][][] conflicts_old = new int [nsTiles.length][][];
		for (int isTile = 0; isTile < nsTiles.length; isTile++){
			conflicts_old[isTile] = iconflicts.detectTriangularTileConflicts(
					nsTiles[isTile],   // int nsTile0,
					replacement_tiles, //HashMap<Integer,Integer> replacement_tiles, //
					neibs_prev, // int [][][] replacement_neibs,
					tnSurface); // TileNeibs tnSurface)
		}
		if (debugLevel > 1) {
			System.out.println("Involved supertiles:");
			for (int i = 0; i < nsTiles.length; i++){
				System.out.println(i+":"+nsTiles[i]);
			}
		}

		if (debugLevel > 1) {
			System.out.println("Calculating original conflicts");
		}

		iconflicts.addConflicts(conflicts_old,
				debugLevel - 1); // debugLevel);
		// After getting old (referfence data) iterate through all variants, for each calculate cost and number of conflicts.
		// First - just collect data - cost and full statistics of the conflicts, print them
		// Then select the best one (not incrementing number of conflicts? and reducing cost)
		int [][][][] variant_conflicts = new int [neibs_vars.length][nsTiles.length][][];

		Conflicts [] variant_conflicts_stats = new Conflicts [neibs_vars.length];


		double [] variant_costs_diff = new double [neibs_vars.length];
		for (int variant = 0; variant < neibs_vars.length; variant ++){
			for (int isTile = 0; isTile < nsTiles.length; isTile++){
				variant_conflicts[variant][isTile] = iconflicts.detectTriangularTileConflicts(
						nsTiles[isTile],   // int nsTile0,
						replacement_tiles, //HashMap<Integer,Integer> replacement_tiles, //
						neibs_vars[variant], // neibs, // int [][][] replacement_neibs,
						tnSurface); // TileNeibs tnSurface)
			}
			variant_conflicts_stats[variant] = new Conflicts(
					variant_conflicts[variant],
					this,
					debugLevel - 1); // debugLevel);
			variant_conflicts_stats[variant].subConflicts(iconflicts); // subtract old number of different types of conflicts
			if (debugLevel > 1) {
				System.out.println("resolveStarConflict(): resolving conflict for tile "+nsTile+
						", nl1 = "+nl1+
						", nl2 = "+nl2 +
						", variant = "+variant+
						", num_conflicts diff = "+variant_conflicts_stats[variant].getNumConflicts()
						);

			}
			if (!newConfl && (variant_conflicts_stats[variant].getNumConflicts() > 0)){
				continue;
			}
			variant_costs_diff[variant] = 	connectionCosts.getConnectionsCostDiff(
					neibs_vars[variant],
					debugLevel);

			if (debugLevel > 0) {
				System.out.println("resolveStarConflict(): resolving conflict for tile "+nsTile+
						", nl1 = "+nl1+
						", nl2 = "+nl2 +
						", variant = "+variant+" improvement (negative diff) = "+variant_costs_diff[variant]+
						", num_conflicts diff = "+variant_conflicts_stats[variant].getNumConflicts());
			}


			if (debugLevel > 0) { // -1) {
				variant_conflicts_stats[variant].printConflictSummary(
						"Conflicts difference after resolution:", true, true, false);
			}
		}
		// How to compare? 1 attempt: none of the conflicts get worse, some get better or cost goes down
		int best_variant = -1;
		int best_ignore_conflicts = -1;
		int [][] num_better_worse = new int [neibs_vars.length][2];
		/*
		for (int variant = 0; variant < neibs_vars.length; variant ++){
			int num_worse = variant_conflicts_stats[variant].numBetterWorse(
					false, // boolean better,
					false, // boolean use_all,
					true, // boolean use_odo,
					false); // ); // boolean use_ood)

			int num_better = variant_conflicts_stats[variant].numBetterWorse(
					true, // boolean better,
					false, // boolean use_all,
					true, // boolean use_odo,
					false); // ); // boolean use_ood)
			num_better_worse[variant][0] = num_better;
			num_better_worse[variant][1] = num_worse;
			if ((num_worse == 0) &&
					(variant_costs_diff[variant] <= dblTriLoss) && // not too worse
					((variant_costs_diff[variant] < 0) || (num_better > 0)) && // either
					((best_variant < 0) ||
							( num_better_worse[variant][0] > num_better_worse[best_variant][0]) ||
							(( num_better_worse[variant][0] == num_better_worse[best_variant][0]) &&
									(variant_costs_diff[variant] < variant_costs_diff[best_variant])))){
				best_variant = variant;
			}
			if ((best_ignore_conflicts <0) || (variant_costs_diff[variant] < variant_costs_diff[best_ignore_conflicts])){
				best_ignore_conflicts = variant;
			}
		}
*/
		int [] num_var_conflicts = new int [neibs_vars.length];
		for (int variant = 0; variant < neibs_vars.length; variant ++){
//			int num_conflicts = variant_conflicts_stats[variant].getNumConflicts();
			num_var_conflicts[variant] = variant_conflicts_stats[variant].getNumConflicts();
			if ((num_var_conflicts[variant] <= 0) &&
					(variant_costs_diff[variant] <= dblTriLoss) && // not too worse
					((variant_costs_diff[variant] < 0) || (num_var_conflicts[variant] < 0)) && // either
					((best_variant < 0) ||
							(num_var_conflicts[variant] < num_var_conflicts[best_variant]) ||
							((num_var_conflicts[variant] == num_var_conflicts[best_variant]) &&
									(variant_costs_diff[variant] < variant_costs_diff[best_variant])))){
				best_variant = variant;
			}
			if ((best_ignore_conflicts <0) || (variant_costs_diff[variant] < variant_costs_diff[best_ignore_conflicts])){
				best_ignore_conflicts = variant;
			}
		}

		if (debugLevel > 1){
			System.out.println("resolveStarConflict(): for tile "+nsTile);
		}

		if ((best_variant < 0) && newConfl && (variant_costs_diff[best_ignore_conflicts] < 0)){ // should be cost improvement
			best_variant = best_ignore_conflicts;
			if (debugLevel > -1) {
				System.out.println("resolveMultiTriangularConflict(): conflicts increase but cost decreases "+nsTile+
						", nl1 = "+nl1+
						", nl2 = "+nl2 +" of "+ neibs_vars.length+" variants");
			}
		}
		if ((best_variant < 0) || (variant_costs_diff[best_variant] > dblTriLoss)){
			if (debugLevel > -1) {
				System.out.println("resolveMultiTriangularConflict(): FAILED find a sutable solution for tile "+nsTile+
						", nl1 = "+nl1+
						", nl2 = "+nl2 +" of "+ neibs_vars.length+" variants");
				return false;
			}
		} else {
			if (debugLevel > -1) {
				System.out.println("resolveStarConflict(): SUCCESS to find a sutable solution for tile "+nsTile+
						", nl1 = "+nl1+
						", nl2 = "+nl2 +". Of "+ neibs_vars.length+" variants - use variant # " + best_variant+
						" cost difference (negative) = "+variant_costs_diff[best_variant]+
						" num conflict reductions = "+(-num_var_conflicts[best_variant])+" (old "+num_better_worse[best_variant][0]+")");
				variant_conflicts_stats[best_variant].printConflictSummary(
						"Conflicts number change per type: ",
						true, // use_all,
						true, //use_odo,
						true); // use_ood);
				iconflicts.printConflictSummary(
						"Conflicts before resolution: ",
						true, // use_all,
						true, //use_odo,
						true); // use_ood);
				// update statistics
				conflict_stats.addConflicts(variant_conflicts_stats[best_variant]);


				// update conflict
				for (int i = 0; i < nsTiles.length; i++){
					conflicts[nsTiles[i]]=  variant_conflicts[best_variant][i];
				}


				// apply resolution
				for (int i = 0; i < mod_supertiles.length; i++){
					if (debugLevel > 1){
						System.out.println("resolveStarConflict(): nsTile = "+nsTile+ "mod_supertiles["+i+"]="+mod_supertiles[i]);
					}
					if (neibs_vars[best_variant][i] != null) {
						for (int nl = 0; nl < neibs_vars[best_variant][i].length; nl ++){
							if (debugLevel > 1){
								System.out.println("resolveStarConflict(): nl= = "+nl);
							}
							if (neibs_vars[best_variant][i][nl] != null){
								planes[mod_supertiles[i]][nl].setNeibBest(neibs_vars[best_variant][i][nl]);
							}
						}
					}
				}
			}
		}
		return true;
	}

	public boolean resolveStarConflict(
			int nsTile,
			int nl1,
			int nl2,
			LinkPlanes lp,
			TileNeibs tnSurface,
			int [][][] conflicts,
			Conflicts conflict_stats, // to be updated after applying resolution
			int        starSteps, // How far to look around when calculationg connection cost
			double     orthoWeight,
			double     diagonalWeight,
			double     starPwr, // Divide cost by number of connections to this power
			double     starWeightPwr,    // Use this power of tile weight when calculating connection cost
			double     weightToDens,    // Balance weighted density against density. 0.0 - density, 1.0 - weighted density
			double     starValPwr, //  Raise value of each tile before averaging
			double     dblTriLoss, //  When resolving double triangles allow minor degradation (0.0 - strict)
			boolean    newConfl, // Allow more conflicts if overall cost is reduced
			int        maxChanges,  // Maximal number of simultaneous connection changes around one tile (0 - any)
			boolean    preferDisparity,
			int        debugLevel)
	{

		if (newConfl && (dblTriLoss > 0.0)){
			dblTriLoss = 0.0; // require strict reducing of the cost if conflict increase is allowed
		}
		Conflicts iconflicts = new Conflicts(this);

		TwoLayerNeighbors twoLayerNeighbors = new TwoLayerNeighbors();

		for (int dir = -1; dir < 8; dir++){
			int nt = tnSurface.getNeibIndex(nsTile, dir);
			if ((nt >= 0) && (planes[nt] != null)) {
				int [][] neibs =  new int [planes[nt].length][];
				boolean [][][] merge_valid = new boolean [planes[nt].length][][];
				for (int nl = 0; nl < planes[nt].length; nl++) if ( planes[nt][nl] != null){
					neibs[nl] =  planes[nt][nl].getNeibBest();
					merge_valid[nl] = planes[nt][nl].getMergedValid();
				}
				twoLayerNeighbors.setNeighbors(neibs,dir);
				twoLayerNeighbors.setMergeValid(merge_valid,dir);
			}
		}
		twoLayerNeighbors.setLayers(nl1, nl2);
		if (debugLevel > 1) {
			System.out.println("resolveStarConflict(): nsTile ="+nsTile+" nl1="+nl1+" nl2="+nl2);
		}

		int [][][][] neibs_vars_dir = twoLayerNeighbors.getNeighborVariants(
				maxChanges,
				debugLevel);

		int [] mod_supertiles = {nsTile};
		mod_supertiles = getInvolvedSupertiles( // first mod_supertiles.length entries will be mod_supertiles[]
				mod_supertiles,
				tnSurface);

		HashMap<Integer,Integer> replacement_tiles = new HashMap<Integer,Integer>();
		for (int i = 0; i < mod_supertiles.length; i++){
			replacement_tiles.put(mod_supertiles[i], new Integer(i));
		}

		// up to 9 tiles
		int [] indexToDir = new int[mod_supertiles.length];
		for (int i = 0; i < indexToDir.length; i++) indexToDir[i] = -1;
		for (int dir = -1; dir < 8; dir++){
			int nindx = (dir < 0) ? 8 : dir;
			int nt = tnSurface.getNeibIndex(nsTile, dir);
			int indx = -1;
			Integer Indx = replacement_tiles.get(nt);
			if (Indx == null){
				System.out.println("resolveStarConflict(): nsTile = "+nsTile+" nindx="+nindx+" Indx == null ");
			} else {
				indx = Indx;
			}
			if (indx >= 0) {
				indexToDir[indx] = nindx;
			}
		}
		int [][][][] neibs_vars = new int [neibs_vars_dir.length][][][];
		for (int variant = 0; variant < neibs_vars_dir.length; variant ++){
			neibs_vars[variant] = new int [indexToDir.length][][];
			for (int i = 0; i <indexToDir.length; i++){
				if (indexToDir[i] >=0){
					neibs_vars[variant][i] = neibs_vars_dir[variant][indexToDir[i]];
				} else {
					System.out.println("resolveStarConflict(): a BUG:  indexToDir["+i+"] = "+indexToDir[i]);
				}
			}
		}

		// See how this application will influence number of conflicts
		// All supertiles that may have different conflicts
		int [] nsTiles = getInvolvedSupertiles( // first mod_supertiles.length entries will be mod_supertiles[]
				mod_supertiles,
				tnSurface);

		ConnectionCosts connectionCosts = new ConnectionCosts(
				orthoWeight,
				diagonalWeight,
				starPwr,    // Divide cost by number of connections to this power
				starWeightPwr,    // Use this power of tile weight when calculating connection cost
				weightToDens,    // Balance weighted density against density. 0.0 - density, 1.0 - weighted density
				starValPwr,      //double     starValPwr, //  Raise value of each tile before averaging
				starSteps,
				this.planes,
				tnSurface,
				preferDisparity);

		int [][][] neibs_prev = connectionCosts.initConnectionCosts(mod_supertiles, debugLevel);
// connectionCosts now contains last calculated	val/weight pairs for broader array of tile data

		int [][][] conflicts_old = new int [nsTiles.length][][];
		double conflicts_old_cost = 0.0;
		for (int isTile = 0; isTile < nsTiles.length; isTile++){
			conflicts_old[isTile] = iconflicts.detectTriangularTileConflicts(
					nsTiles[isTile],   // int nsTile0,
					replacement_tiles, // HashMap<Integer,Integer> replacement_tiles, //
					neibs_prev, // int [][][] replacement_neibs,
					tnSurface); // TileNeibs tnSurface)
// calculate cost of conflicts

			conflicts_old_cost += iconflicts.getConflictsCost(
					nsTiles[isTile],       // int           nsTile0,
					1.0,                   // double        scaleStartEnd, // include start and and layer tiles in the center in overall cost for each triangle (1.0)
					conflicts_old[isTile], // int [][]      conflicts, //
					replacement_tiles,     // HashMap<Integer,Integer> replacement_tiles, // null is OK
					neibs_prev,            // int [][][]    replacement_neibs,               // null OK if  replacement_tiles == null
					connectionCosts,       // ConnectionCosts connectionCosts,
					tnSurface); // TileNeibs tnSurface)


		}
		if (debugLevel > 1) {
			System.out.println("Sum of old conflict costs around " + nsTile + " is "+conflicts_old_cost);
			System.out.println("Involved supertiles:");
			for (int i = 0; i < nsTiles.length; i++){
				System.out.println(i+":"+nsTiles[i]);
			}
		}

		if (debugLevel > 1) {
			System.out.println("Calculating original conflicts");
		}

		iconflicts.addConflicts(conflicts_old,
				debugLevel - 1); // debugLevel);
		// After getting old (referfence data) iterate through all variants, for each calculate cost and number of conflicts.
		// First - just collect data - cost and full statistics of the conflicts, print them
		// Then select the best one (not incrementing number of conflicts? and reducing cost)
		int [][][][] variant_conflicts = new int [neibs_vars.length][nsTiles.length][][];

		Conflicts [] variant_conflicts_stats = new Conflicts [neibs_vars.length];


		double [] variant_costs_diff = new double [neibs_vars.length];
		double [] conflicts_var_cost = new double [neibs_vars.length];
		for (int variant = 0; variant < neibs_vars.length; variant ++){
			// connections costs are needed before conflicts
			variant_costs_diff[variant] = 	connectionCosts.getConnectionsCostDiff(
					neibs_vars[variant],
					debugLevel);
			conflicts_var_cost[variant] = 0.0;
			for (int isTile = 0; isTile < nsTiles.length; isTile++){
				variant_conflicts[variant][isTile] = iconflicts.detectTriangularTileConflicts(
						nsTiles[isTile],   // int nsTile0,
						replacement_tiles, //HashMap<Integer,Integer> replacement_tiles, //
						neibs_vars[variant], // neibs, // int [][][] replacement_neibs,
						tnSurface); // TileNeibs tnSurface)
				// connectionCosts now contains last calculated	val/weight pairs for broader array of tile data
				conflicts_var_cost[variant] += iconflicts.getConflictsCost(
						nsTiles[isTile],                    // int           nsTile0,
						1.0,                                // double        scaleStartEnd, // include start and and layer tiles in the center in overall cost for each triangle (1.0)
						variant_conflicts[variant][isTile], // int [][]      conflicts, //
						replacement_tiles,                  // HashMap<Integer,Integer> replacement_tiles, // null is OK
						neibs_vars[variant],                // int [][][]    replacement_neibs,               // null OK if  replacement_tiles == null
						connectionCosts,                    // ConnectionCosts connectionCosts,
						tnSurface);                         // TileNeibs tnSurface)
			}
			variant_conflicts_stats[variant] = new Conflicts(
					variant_conflicts[variant],
					this,
					debugLevel - 1); // debugLevel);
			variant_conflicts_stats[variant].subConflicts(iconflicts); // subtract old number of different types of conflicts
			if (debugLevel > 1) {
				System.out.println("resolveStarConflict(): resolving conflict for tile "+nsTile+
						", nl1 = "+nl1+
						", nl2 = "+nl2 +
						", variant = "+variant+
						"  conflit cost change = "+ (conflicts_var_cost[variant] - conflicts_old_cost)+
						", num_conflicts diff = "+variant_conflicts_stats[variant].getNumConflicts()
						);

// Initially do not process conflicts cost difference, just view it
			}
/*
			if (!newConfl && (variant_conflicts_stats[variant].getNumConflicts() > 0)){
				continue;
			}
*/

			if (debugLevel > 0) {
				System.out.println("resolveStarConflict(): resolving conflict for tile "+nsTile+
						", nl1 = "+nl1+
						", nl2 = "+nl2 +
						", variant = "+variant+" improvement (negative diff) = "+variant_costs_diff[variant]+
						"  conflit cost change = "+ (conflicts_var_cost[variant] - conflicts_old_cost)+
						", num_conflicts diff = "+variant_conflicts_stats[variant].getNumConflicts());
			}


			if (debugLevel > 0) { // -1) {
				variant_conflicts_stats[variant].printConflictSummary(
						"Conflicts difference after resolution:", true, true, false);
			}
		}
		// How to compare? 1 attempt: none of the conflicts get worse, some get better or cost goes down
		int best_variant = -1;
		int best_ignore_conflicts = -1;
//		int [][] num_better_worse = new int [neibs_vars.length][2];
		int [] num_var_conflicts = new int [neibs_vars.length];


		for (int variant = 0; variant < neibs_vars.length; variant ++){
			if (variant_costs_diff[variant] <= dblTriLoss) {
				if  (conflicts_var_cost[variant] < conflicts_old_cost) {
					if ((best_variant < 0) || (conflicts_var_cost[variant] < conflicts_var_cost[best_variant])) {
						best_variant = variant;
					}
				}
				if ((best_ignore_conflicts < 0) || (variant_costs_diff[variant] < variant_costs_diff[best_ignore_conflicts])){
					best_ignore_conflicts = variant;
				}
			}
		}

		if (debugLevel > 1){
			System.out.println("resolveStarConflict(): for tile "+nsTile);
		}
/*
		if ((best_variant < 0) && newConfl && (best_ignore_conflicts >= 0) && (variant_costs_diff[best_ignore_conflicts] < 0)){ // should be cost improvement
			best_variant = best_ignore_conflicts;
			if (debugLevel > -1) {
				System.out.println("resolveMultiTriangularConflict(): conflicts increase but cost decreases "+nsTile+
						", nl1 = "+nl1+
						", nl2 = "+nl2 +" of "+ neibs_vars.length+" variants");
			}
		}
*/
		if (newConfl){ // should be cost improvement
			best_variant = best_ignore_conflicts;
			if (debugLevel > -1) {
				System.out.println("resolveStarConflict(): ignoring conflicts if any for "+nsTile+
						", nl1 = "+nl1+
						", nl2 = "+nl2 +" of "+ neibs_vars.length+" variants");
			}
		}


		if ((best_variant < 0) || (variant_costs_diff[best_variant] > dblTriLoss) ||
//				((variant_costs_diff[best_variant] >= 0.0) && (conflicts_var_cost[best_variant] >= conflicts_old_cost) &&
				((conflicts_var_cost[best_variant] >= conflicts_old_cost) && // try always prohibit worse conflicts
						(num_var_conflicts[best_variant] >= 0))){
			if (debugLevel > -1) {
				System.out.println("resolveStarConflict(): FAILED find a sutable solution for tile "+nsTile+
						", nl1 = "+nl1+
						", nl2 = "+nl2 +" of "+ neibs_vars.length+" variants");
				if (debugLevel > 1){
					System.out.println("resolveStarConflict(): for tile "+nsTile+" FAILURE");
				}
				return false;
			}
		} else {
			if (debugLevel > -1) {
				System.out.println("resolveStarConflict(): SUCCESS to find a sutable solution for tile "+nsTile+
						", nl1 = "+nl1+
						", nl2 = "+nl2 +". Of "+ neibs_vars.length+" variants - use variant # " + best_variant+
						"  conflit cost change = "+ (conflicts_var_cost[best_variant] - conflicts_old_cost)+
						" cost difference (negative) = "+variant_costs_diff[best_variant]+
						" num conflict reductions = "+(-num_var_conflicts[best_variant]));
				variant_conflicts_stats[best_variant].printConflictSummary(
						"Conflicts number change per type: ",
						true, // use_all,
						true, //use_odo,
						true); // use_ood);
				iconflicts.printConflictSummary(
						"Conflicts before resolution: ",
						true, // use_all,
						true, //use_odo,
						true); // use_ood);
			}
			// update statistics
			conflict_stats.addConflicts(variant_conflicts_stats[best_variant]);

			// update conflict
			for (int i = 0; i < nsTiles.length; i++){
				conflicts[nsTiles[i]]=  variant_conflicts[best_variant][i];
			}


			// apply resolution
			for (int i = 0; i < mod_supertiles.length; i++){
				if (debugLevel > 1){
					System.out.println("resolveStarConflict(): nsTile = "+nsTile+ "mod_supertiles["+i+"]="+mod_supertiles[i]);
				}
				if (neibs_vars[best_variant][i] != null) {
					for (int nl = 0; nl < neibs_vars[best_variant][i].length; nl ++){
						if (debugLevel > 1){
							System.out.println("resolveStarConflict(): nl= = "+nl);
						}
						if (neibs_vars[best_variant][i][nl] != null){
							planes[mod_supertiles[i]][nl].setNeibBest(neibs_vars[best_variant][i][nl]);
						}
					}
				}
			}
			// recalculate starValueWeights for and around tiles with modified  neighbors (no outside connections changed )nsTiles
			lp.updateStarValueStrength(
					nsTiles,          // final int []         mod_supertiles,
					orthoWeight,      // final double         orthoWeight,
					diagonalWeight,   // final double         diagonalWeight,
					starPwr,          // final double         starPwr,    // Divide cost by number of connections to this power
					starWeightPwr,    // Use this power of tile weight when calculating connection cost
					weightToDens,    // Balance weighted density against density. 0.0 - density, 1.0 - weighted density
					starValPwr,       //double     starValPwr, //  Raise value of each tile before averaging
					starSteps,        // final int            steps,
					planes,           // final TilePlanes.PlaneData [][] planes,
					preferDisparity, // final boolean        preferDisparity)
					debugLevel);
		}
		if (debugLevel > 1){
			System.out.println("resolveStarConflict(): for tile "+nsTile+" OK");
		}

		return true;
	}

	public void calcStarValueStrengthOld(
			final double         orthoWeight,
			final double         diagonalWeight,
			final double         starPwr,    // Divide cost by number of connections to this power
			final double         starWeightPwr,    // Use this power of tile weight when calculating connection cost
			final double         weightToDens,    // Balance weighted density against density. 0.0 - density, 1.0 - weighted density
			final double         starValPwr, //  Raise value of each tile before averaging
			final int            steps,
			final TilePlanes.PlaneData [][] planes,
			final boolean        preferDisparity,
			final int            debugLevel)
	{
		final int tilesX =        tileProcessor.getTilesX();
		final int tilesY =        tileProcessor.getTilesY();
		final int superTileSize = tileProcessor.getSuperTileSize();
		//				final int tileSize =      tileProcessor.getTileSize();
		final int stilesX =       (tilesX + superTileSize -1)/superTileSize;
		final int stilesY =       (tilesY + superTileSize -1)/superTileSize;
		final int nStiles =       stilesX * stilesY;
		final TileNeibs tnSurface = new TileNeibs(stilesX, stilesY);
		final Thread[] threads = ImageDtt.newThreadArray(tileProcessor.threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				@Override
				public void run() {
					ConnectionCosts connectionCosts = new ConnectionCosts(
							orthoWeight,    // double         orthoWeight,
							diagonalWeight, // double         diagonalWeight,
							starPwr,        // double         starPwr,    // Divide cost by number of connections to this power
							starWeightPwr,  // double         starWeightPwr,    // Use this power of tile weight when calculating connection cost
							weightToDens,    // Balance weighted density against density. 0.0 - density, 1.0 - weighted density
							starValPwr,     //double          starValPwr, //  Raise value of each tile before averaging
							steps,          // int            steps,
							planes,         // TilePlanes.PlaneData [][] planes,
							tnSurface,      // TileNeibs tnSurface,
							preferDisparity); // boolean preferDisparity)
					int [] mod_supertiles = new int[1];
					for (int nsTile = ai.getAndIncrement(); nsTile < nStiles; nsTile = ai.getAndIncrement()) {
						if ( planes[nsTile] != null) {
							mod_supertiles[0] = nsTile;
							connectionCosts.initConnectionCosts(mod_supertiles, debugLevel);
							double [][][] val_weights = connectionCosts.getValWeights();
							for (int np = 0; np < planes[nsTile].length; np++){ // nu
								if (planes[nsTile][np] != null) {
									planes[nsTile][np].setStarValueWeight(val_weights[0][np]);
								}
							}
						}
					}
				}
			};
		}
		ImageDtt.startAndJoin(threads);
	}

	public void updateStarValueStrengthOld(
			final int []         mod_supertiles,
			final double         orthoWeight,
			final double         diagonalWeight,
			final double         starPwr,    // Divide cost by number of connections to this power
			final double         starWeightPwr,    // Use this power of tile weight when calculating connection cost
			final double         weightToDens,    // Balance weighted density against density. 0.0 - density, 1.0 - weighted density
			final double         starValPwr, //  Raise value of each tile before averaging
			final int            steps,
			final TilePlanes.PlaneData [][] planes,
			final boolean        preferDisparity,
			final int            debugLevel)
	{
		final int tilesX =        tileProcessor.getTilesX();
		final int tilesY =        tileProcessor.getTilesY();
		final int superTileSize = tileProcessor.getSuperTileSize();
		//				final int tileSize =      tileProcessor.getTileSize();
		final int stilesX =       (tilesX + superTileSize -1)/superTileSize;
		final int stilesY =       (tilesY + superTileSize -1)/superTileSize;
		final TileNeibs tnSurface = new TileNeibs(stilesX, stilesY);
		final Thread[] threads = ImageDtt.newThreadArray(tileProcessor.threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				@Override
				public void run() {
					ConnectionCosts connectionCosts = new ConnectionCosts(
							orthoWeight,    // double         orthoWeight,
							diagonalWeight, // double         diagonalWeight,
							starPwr,        // double         starPwr,    // Divide cost by number of connections to this power
							starWeightPwr,    // Use this power of tile weight when calculating connection cost
							weightToDens,    // Balance weighted density against density. 0.0 - density, 1.0 - weighted density
							starValPwr,      //double     starValPwr, //  Raise value of each tile before averaging
							steps,          // int            steps,
							planes,         // TilePlanes.PlaneData [][] planes,
							tnSurface,      // TileNeibs tnSurface,
							preferDisparity); // boolean preferDisparity)
					int [] supertiles = new int[1];
					for (int isTile = ai.getAndIncrement(); isTile < mod_supertiles.length; isTile = ai.getAndIncrement()) {
						int nsTile = mod_supertiles[isTile];
						if ((nsTile >= 0) && ( planes[nsTile] != null)) {
							supertiles[0] = nsTile;
							connectionCosts.initConnectionCosts(supertiles, debugLevel - 2);
							double [][][] val_weights = connectionCosts.getValWeights();
							for (int np = 0; np < planes[nsTile].length; np++){ // nu
								if (planes[nsTile][np] != null) {
									planes[nsTile][np].setStarValueWeight(val_weights[0][np]);
								}
							}
						}
					}
				}
			};
		}
		ImageDtt.startAndJoin(threads);
	}


	//
	/**
	 * Generate variants for changing connections while preserving number of connections between each pair of tiles
	 * each variant is free of own conflicts, but may lead to conflicts on other layers or tiles
	 * These variants should be checked on number of conflicts and overall cost to make sure conflict resolution leads
	 * to improvement - either in number of conflicts or global connection cost
	 * @param nsTile linescan tile index of the supertile
	 * @param nl1 first layer at nsTile to have triangle (2 ortho and a diagonal between them) conflict
	 * @param nl2 second layer at nsTile to have a triangle conflict, nl2 > nl1
	 * @param dir_mask start directions for the triangle conflicts 1 - N, 2 - E, 4 - S, 8 - W (starting from nl1),
	 * next 4 bits - similar fro starting from nl2
	 * @param tnSurface ileSurface.TileNeibs instance to generate indices of the tiles at certain direction from the current
	 * @param debugLevel debug level
	 * @return [variant][5][nl][8] [5] = 4 directions and center ([4]), some may be null - should be removed by the caller
	 */

	public int [][][][] getTriangularResolutionVariants(
			int nsTile,
			int nl1,
			int nl2,
			int dir_mask,
			TileNeibs tnSurface,
			int        debugLevel)
	{
		int [][] layers0 = new int[5][2]; // 4 - center, 0 - N, 1 - E, 2 - S, 3 - w. 4 has 2, each other 0,1 or 2 layers
		for (int i = 0; i < (layers0.length - 1); i++){
			layers0[i][0] = -1;
			layers0[i][1] = -1;
		}
		layers0[4][0] = nl1;
		layers0[4][1] = nl2;
		int [] neibs1 = planes[nsTile][nl1].getNeibBest();
		int [] neibs2 = planes[nsTile][nl2].getNeibBest();
		for (int dir4 = 0; dir4 < 4; dir4++){
			if ((dir_mask & (1 << dir4)) != 0) {
				layers0[dir4][0] = neibs1[2 * dir4];
				layers0[(dir4 + 1) % 4][1] = neibs2[(2 * (dir4 + 1)) % 8];
			}
			if ((dir_mask & (16 << dir4)) != 0) {
				layers0[(dir4 + 1) % 4][0] = neibs1[(2 * (dir4 + 1)) % 8];
				layers0[dir4][1] = neibs2[2 * dir4];
			}
		}
		int [][] layers = new int[5][];
		int [] nsTiles = new int [5];
		nsTiles[4] = nsTile;
		for (int i = 0; i < 5; i++){
			if ((layers0[i][0] >=0 ) || (layers0[i][1] >=0 )) {
				if ((layers0[i][0] >=0 ) && (layers0[i][1] >=0 )){
					layers[i] = layers0[i].clone();
				} else {
					layers[i] = new int[1];
					layers[i][0] = (layers0[i][0] >=0 ) ? layers0[i][0] : layers0[i][1];
				}

			} else {
				layers[i] = new int[0];
				nsTiles[i] = -1;
			}
		}
		for (int dir4 = 0; dir4 < 4; dir4++){
			nsTiles[dir4] = tnSurface.getNeibIndex(nsTile, 2 * dir4);
		}
		int [] diag_connections = new int [4];// number of variants for diagonal connections (1 or 2 where applicable)
		for (int dir4 = 0; dir4 < 4; dir4 ++){
			int dir4_next = (dir4 + 1) % 4;
			int dir_diag = (2 * dir4 + 3) % 8;
			diag_connections[dir4] = 0;
			for (int il1 = 0; il1 < layers[dir4].length; il1++){
				int [] neibs_s = planes[nsTiles[dir4]][layers[dir4][il1]].getNeibBest();
				for (int il2 = 0; il2 < layers[dir4_next].length; il2++){
					if (neibs_s[dir_diag] == layers[dir4_next][il2]) diag_connections[dir4] ++;
				}
			}
		}
		int num_vars = 1;
		for (int dir4 = 0; dir4 < 4; dir4 ++){
			if (layers[dir4].length > 0) num_vars *=2;
			// number of ways to connect diagonals, preserving same number of links
			if ((layers[dir4].length == 2) && (layers[(dir4 + 1) % 4].length == 2) && (diag_connections[dir4] == 1)) num_vars *= 2;
		}
		if (debugLevel > -1){
			System.out.println("getTriangularResolutionVariants() nsTile = "+nsTile+" variants: "+num_vars+" (may be actually less)");
		}
		int [][][][] neib_vars = new int [num_vars][][][]; // [5][][];
		int [][][] neib_var_zero = new int [nsTiles.length][][];
		for (int i = 0; i < nsTiles.length; i++){
			if (layers[i].length > 0){
				neib_var_zero[i] = new int [planes[nsTiles[i]].length][];
				for (int nl = 0; nl < planes[nsTiles[i]].length; nl++){
					if (planes[nsTiles[i]][nl] != null){
						neib_var_zero[i][nl] = planes[nsTiles[i]][nl].getNeibBest();
					}
				}
			}
		}
		 // [i][0] - which of the neighbors layers to connect to nl1, [i][1] which of the diagonal connections (if choice) to keep
		boolean [][] selection = new boolean [4][2];
		int var_num = 0;
		while(true){
			int [][][] neib_var = neib_var_zero.clone();
			for (int i = 0; i < neib_var_zero.length; i++) if (neib_var_zero[i] != null){
				neib_var[i] = neib_var_zero[i].clone();
				for (int nl = 0; nl < neib_var_zero[i].length; nl++) if (neib_var_zero[i][nl] != null){
					neib_var[i][nl] = neib_var_zero[i][nl].clone();
				}
			}
			// use current selection to build neib_vars[var_num]. Increment [var_num] if OK, not to skip invalid variants
			// start with connecting no nl1, then - to nl2 and finally - diagonals
			good_variant:
			{
				for (int dir4 = 0; dir4 < 4; dir4 ++){
					int dir4_next = (dir4 + 1) % 4;
					int dir_f = 2 * dir4;
					int dir_r = (dir_f +4 ) % 8;
					int dir12_f = (2 * dir4 + 3 ) % 8;
					int dir12_r = (dir12_f + 4 ) % 8;

					if (layers[dir4].length > 0) {
						if (layers[dir4].length > 1){
							neib_var[4][nl1][dir_f] = selection[dir4][0] ? layers[dir4][1] : layers[dir4][0];
							neib_var[4][nl2][dir_f] = selection[dir4][0] ? layers[dir4][0] : layers[dir4][1];
							neib_var[dir4][layers[dir4][0]][dir_r] = selection[dir4][0] ? nl2 : nl1;
							neib_var[dir4][layers[dir4][1]][dir_r] = selection[dir4][0] ? nl1 : nl2;

						} else { // only 1 connection
							int other_layer = (neib_var_zero[4][nl1][dir_f] == layers[dir4][0]) ? neib_var_zero[4][nl2][dir_f] : neib_var_zero[4][nl1][dir_f];
							neib_var[4][nl1][dir_f] = selection[dir4][0] ? other_layer : layers[dir4][0];
							neib_var[4][nl2][dir_f] = selection[dir4][0] ? layers[dir4][0] : other_layer;
							neib_var[dir4][layers[dir4][0]][dir_r] = selection[dir4][0] ? nl2 : nl1;
							if (other_layer >= 0) { // only if it is connected
								neib_var[dir4][other_layer][dir_r] = selection[dir4][0] ? nl1 : nl2; // null
							}
						}
						// diagonal connections - some make this variant invalid
						if (layers[dir4_next].length > 0) { // only add and verify connections if there are both radial connections
							int last_index_1 = (layers[dir4].length > 1) ? 1 :0;
							int last_index_2 = (layers[dir4_next].length > 1) ? 1 :0;
							int [][] neibs1_x = {neib_var_zero[dir4][layers[dir4][0]], // layers connected the first/single tile at dir4 from the center
									neib_var_zero[dir4][layers[dir4][last_index_1]]}; // layers connected the second/single tile at dir4 from the center
							int [][] neibs2_x = {
									neib_var_zero[dir4_next][layers[dir4_next][0]], // layers connected the first/single tile at dir4_next from the center
									neib_var_zero[dir4_next][layers[dir4_next][last_index_2]]}; // layers connected the last/single tile at dir4_next from the center
							// if there is a single variant. ( already connected)
							//if ((layers[dir4].length => 1)))
							if ((layers[dir4].length == 1) && (layers[dir4_next].length == 1)){ // 1 + 1
								if (selection[dir4][0] == selection[dir4_next][0]) continue; // no more checks - connected or not
								if (neibs1_x[0][dir12_f] != layers[dir4_next][0]) continue;     // not connected - keep whatever connections they had
								// connected, but radials are to the to different layers - invalid
								break good_variant;
							} else if ((layers[dir4].length == 2) && (layers[dir4_next].length == 2)){ // 2 + 2
								if      (diag_connections[dir4] == 0) continue;
								// want to connect layers[dir4][0] (in the direction dir12_f) to :
								int want_10 = (selection[dir4][0] ^ selection[dir4_next][0]) ? layers[dir4_next][1] : layers[dir4_next][0];
								// want to connect layers[dir4][1] (in the direction dir12_f) to :
								int want_11 = (selection[dir4][0] ^ selection[dir4_next][0]) ? layers[dir4_next][0] : layers[dir4_next][1];
								// want to connect layers[dir4_next][0] (in the direction dir12_r) to :
								int want_20 = (selection[dir4][0] ^ selection[dir4_next][0]) ? layers[dir4][1] : layers[dir4][0];
								// want to connect layers[dir4_next][1] (in the direction dir12_rf) to :
								int want_21 = (selection[dir4][0] ^ selection[dir4_next][0]) ? layers[dir4][0] : layers[dir4][1];

								if (diag_connections[dir4] == 2) {
									neib_var[dir4][layers[dir4][0]][dir12_f] = want_10;
									neib_var[dir4][layers[dir4][1]][dir12_f] = want_11;
									neib_var[dir4_next][layers[dir4_next][0]][dir12_r] = want_20;
									neib_var[dir4_next][layers[dir4_next][1]][dir12_r] = want_21;
								} else { // 1 connection, (possibly) 2 variants
									// if correct connection, skip second variant. If not - selection[dir4][1]==false - swap first pair, true - the second
									if ((want_10 == neibs1_x[0][dir12_f]) || (want_11 == neibs1_x[1][dir12_f])){
										// Already connected correctly, so just a single variant - skip the second one
										if (!selection[dir4_next][1]) continue;
										break good_variant; // prevent duplicate
									}
									if (!selection[dir4_next][1]) {
										neib_var[dir4][layers[dir4][0]][dir12_f] = neibs1_x[1][dir12_f];
										neib_var[dir4][layers[dir4][1]][dir12_f] = neibs1_x[0][dir12_f];
										// back links only if they existed
										if (neibs1_x[0][dir12_f] >= 0) neib_var[dir4_next][neibs1_x[0][dir12_f]][dir12_r] = layers[dir4][1];
										if (neibs1_x[1][dir12_f] >= 0) neib_var[dir4_next][neibs1_x[1][dir12_f]][dir12_r] = layers[dir4][0];
									} else {
										neib_var[dir4_next][layers[dir4_next][0]][dir12_r] = neibs2_x[1][dir12_r];
										neib_var[dir4_next][layers[dir4_next][1]][dir12_r] = neibs2_x[0][dir12_r];
										// back links only if they existed
										if (neibs2_x[0][dir12_r] >= 0) neib_var[dir4_next][neibs2_x[0][dir12_r]][dir12_f] = layers[dir4_next][1];
										if (neibs2_x[1][dir12_r] >= 0) neib_var[dir4_next][neibs2_x[1][dir12_r]][dir12_f] = layers[dir4_next][0];
									}
								}
							}  else {// remaining 2 + 1 and 1 + 2
								if (diag_connections[dir4] == 0) continue; // no links existed - nothing to do
								// here there was one link
								int conn_indx = (selection[dir4][0] ^ selection[dir4_next][0]) ? 1 : 0;
								int other_indx = 1 - conn_indx;
								if (layers[dir4].length == 1){
									if (neibs1_x[0][dir12_f] == layers[dir4_next][conn_indx]){
										continue; // already matched
									} else { // swap at dir4_next
										neib_var[dir4]     [layers[dir4][0]][dir12_f] = layers[dir4_next][conn_indx];
										neib_var[dir4_next][layers[dir4_next][conn_indx]][dir12_r] = layers[dir4][0];
										neib_var[dir4_next][layers[dir4_next][other_indx]][dir12_r] = neibs2_x[conn_indx][dir12_r];
										if (neibs2_x[conn_indx][dir12_r] >= 0) neib_var[dir4][neibs2_x[conn_indx][dir12_r]][dir12_f] = layers[dir4_next][other_indx];
									}
								} else {
									if (neibs2_x[0][dir12_r] == layers[dir4][conn_indx]){
										continue; // already matched
									} else { // swap at dir4
										neib_var[dir4_next][layers[dir4_next][0]]         [dir12_r] = layers[dir4]     [conn_indx];
										neib_var[dir4]     [layers[dir4]     [conn_indx]] [dir12_f] = layers[dir4_next][0];
										neib_var[dir4]     [layers[dir4]     [other_indx]][dir12_f] = neibs1_x[conn_indx][dir12_f];
										if (neibs1_x[conn_indx][dir12_f] >= 0) neib_var[dir4_next][neibs1_x[conn_indx][dir12_f]][dir12_r] = layers[dir4][other_indx];
									}
								}
							}

						}
					}
				}
				neib_vars[var_num++] = neib_var;
			}
			// calculate next selection
			next_sel:{
				for (int dir4 = 3; dir4 >= 0; dir4--){
					if (!selection[dir4][1] && (layers[dir4].length == 2) && (layers[(dir4 + 1) % 4].length == 2) && (diag_connections[dir4] == 1)){
						selection[dir4][1] = true;
						break next_sel;
//					} else if (!selection[dir4][0] && (layers[dir4].length > 1)){
					} else if (!selection[dir4][0] && (layers[dir4].length > 0)){
						selection[dir4][0] = true;
						selection[dir4][1] = false;
						break next_sel;
					} else {
						selection[dir4][0] = false;
						selection[dir4][1] = false;
					}
				}
				break; // while(true) - all done
			}
		}
// trim neib_vars
		if (debugLevel > -1){
			System.out.println("getTriangularResolutionVariants() nsTile = "+nsTile+" ACTUAL variants: "+var_num);
		}
		if (var_num == neib_vars.length) return neib_vars;
		int [][][][] neib_vars_trimmed = new int [var_num][][][];
		for (int i = 0; i < var_num; i++ ) neib_vars_trimmed[i] = neib_vars[i];

		return neib_vars_trimmed;
	}

	public void resolveConflicts(
			LinkPlanes lp,
			double     maxEigen,      // maximal eigenvalue of planes to consider
	  		boolean    conflDualTri,  // Resolve dual triangles conflict (odoodo)
	  		boolean    conflMulti,    // Resolve multiple odo triangles conflicts
	  		boolean    conflDiag,     // Resolve diagonal (ood) conflicts
	  		boolean    conflStar,     // Resolve all conflicts around a supertile
			int        starSteps, // How far to look around when calculationg connection cost
			double     orthoWeight,
			double     diagonalWeight,
			double     starPwr,         // Divide cost by number of connections to this power
			double     starWeightPwr,   // Use this power of tile weight when calculating connection cost
			double     weightToDens,    // // Balance weighted density against density. 0.0 - density, 1.0 - weighted density
			double     starValPwr, //  Raise value of each tile before averaging
			double     dblTriLoss, //  When resolving double triangles allow minor degradation (0.0 - strict)
			boolean    newConfl, // Allow more conflicts if overall cost is reduced
			int        maxChanges,  // Maximal number of simultaneous connection changes around one tile (0 - any)
			boolean    preferDisparity,
			int        debugLevel,
			int        dbg_X,
			int        dbg_Y)

	{
		lp.calcStarValueStrength(
				true, // boolean set_start_planes,
				orthoWeight,      // final double         orthoWeight,
				diagonalWeight,   // final double         diagonalWeight,
				starPwr,          // final double         starPwr,    // Divide cost by number of connections to this power
				starWeightPwr,    // final double         starWeightPwr,    // Use this power of tile weight when calculating connection cost
				weightToDens,     // Balance weighted density against density. 0.0 - density, 1.0 - weighted density
				starValPwr,      //double     starValPwr, //  Raise value of each tile before averaging
				starSteps,        // final int            steps,
				this.planes,      // final TilePlanes.PlaneData [][] planes,
				preferDisparity,  // final boolean        preferDisparity)
				debugLevel);

		Conflicts iconflicts0 = new Conflicts(this);
		int [][][] conflicts0 = iconflicts0.detectTriangularConflicts(
				1); // final int debugLevel)
		Conflicts conflicts0_stats =  new Conflicts(
				conflicts0,
				this,
				-1); // debugLevel);

		for (int pass = 0; pass < 10; pass ++) {
			int [] dual_tri_results =            {0,0};
			int [] multi_resoultion_results =    {0,0};
			int [] diagonal_resoultion_results = {0,0};
			int [] conflict_star_results = 	     {0,0};
			if (conflDualTri) {
				dual_tri_results = resolveDualTriangularConflicts(
						conflicts0, // int [][][] conflicts,
						conflicts0_stats,
						maxEigen,
						starSteps, // How far to look around when calculating connection cost
						orthoWeight, // double     orthoWeight,
						diagonalWeight, // double     diagonalWeight,
						starPwr,        // double     starPwr, // Divide cost by number of connections to this power
						starWeightPwr,    // Use this power of tile weight when calculating connection cost
						weightToDens,    // // Balance weighted density against density. 0.0 - density, 1.0 - weighted density
						starValPwr,      //double     starValPwr, //  Raise value of each tile before averaging
						dblTriLoss, // double     diagonalWeight,
						preferDisparity,
						debugLevel, // 1, //  final int debugLevel)
						dbg_X,
						dbg_Y);
				System.out.println("Pass "+(pass+1)+": dual_tri_results (success/failures) = "+dual_tri_results[0]+" / "+dual_tri_results[1]);
			}
			if (conflMulti) {
				multi_resoultion_results = resolveMultiTriangularConflicts(
						conflicts0, // int [][][] conflicts,
						conflicts0_stats,
						starSteps, // How far to look around when calculationg connection cost
						orthoWeight, // double     orthoWeight,
						diagonalWeight, // double     diagonalWeight,
						starPwr,        // double     starPwr, // Divide cost by number of connections to this power
						starWeightPwr,    // Use this power of tile weight when calculating connection cost
						weightToDens,    // // Balance weighted density against density. 0.0 - density, 1.0 - weighted density
						starValPwr,      //double     starValPwr, //  Raise value of each tile before averaging
						dblTriLoss, // double     diagonalWeight,
						preferDisparity,
						debugLevel, // 1, //  final int debugLevel)
						dbg_X,
						dbg_Y);
				System.out.println("Pass "+(pass+1)+": multi_tri_results (success/failures) = "+multi_resoultion_results[0]+" / "+multi_resoultion_results[1]);
			}
			if (conflDiag) {
				diagonal_resoultion_results = resolveDiagonalTriangularConflicts(
						conflicts0, // int [][][] conflicts,
						conflicts0_stats,
						starSteps, // How far to look around when calculationg connection cost
						orthoWeight, // double     orthoWeight,
						diagonalWeight, // double     diagonalWeight,
						starPwr,        // double     starPwr, // Divide cost by number of connections to this power
						starWeightPwr,    // Use this power of tile weight when calculating connection cost
						weightToDens,    // // Balance weighted density against density. 0.0 - density, 1.0 - weighted density
						starValPwr,      //double     starValPwr, //  Raise value of each tile before averaging
						dblTriLoss, // double     diagonalWeight,
						preferDisparity,
						((pass == 3)? 1: 0), // debugLevel, // 1, //  final int debugLevel)
						dbg_X,
						dbg_Y);
				System.out.println("Pass "+(pass+1)+": resolveDiagonalTriangularConflicts (success/failures) = "+diagonal_resoultion_results[0]+" / "+diagonal_resoultion_results[1]);
			}
			if (conflStar) {
				conflict_star_results = 	resolveStarConflicts(
						conflicts0, // int [][][] conflicts,
						conflicts0_stats,
						lp, // LinkPlanes lp,
						starSteps, // How far to look around when calculationg connection cost
						orthoWeight, // double     orthoWeight,
						diagonalWeight, // double     diagonalWeight,
						starPwr,        // double     starPwr, // Divide cost by number of connections to this power
						starWeightPwr,    // Use this power of tile weight when calculating connection cost
						weightToDens,    // // Balance weighted density against density. 0.0 - density, 1.0 - weighted density
						starValPwr,      //double     starValPwr, //  Raise value of each tile before averaging
						dblTriLoss, // double     diagonalWeight,
						newConfl, // Allow more conflicts if overall cost is reduced
						maxChanges,     // int        maxChanges,  // Maximal number of simultaneous connection changes around one tile (0 - any)
						preferDisparity,
						debugLevel, // 1, //  final int debugLevel)
						dbg_X,
						dbg_Y);
				System.out.println("Pass "+(pass+1)+": resolveStarConflicts (success/failures) = "+conflict_star_results[0]+" / "+conflict_star_results[1]);
			}

			if (    (dual_tri_results[0] == 0) &&
					(multi_resoultion_results[0] == 0) &&
					(diagonal_resoultion_results[0] == 0) &&
					(conflict_star_results[0] ==       0)) {
				System.out.println("No more improvements");
				break;
			}
		}
		Conflicts conflicts1_stats =  new Conflicts(
				conflicts0,
				this,
				1); // -1); // debugLevel);

		conflicts1_stats.printConflictSummary("Detected conflicts (all):", true,false,false);
		conflicts1_stats.printConflictSummary("Detected conflicts (ortho-diag-ortho):", false, true,false);
		conflicts1_stats.printConflictSummary("Detected conflicts(ortho-ortho-diag):", false, false, true);

		iconflicts0.resetConflicts();
		conflicts0 = iconflicts0.detectTriangularConflicts(
				1); // final int debugLevel)
		conflicts1_stats.printConflictSummary("Recounted conflicts (all):", true,false,false);
		conflicts1_stats.printConflictSummary("Recounted conflicts (ortho-diag-ortho):", false, true,false);
		conflicts1_stats.printConflictSummary("Recounted conflicts (ortho-ortho-diag):", false, false, true);


	}








	public int [] resolveDualTriangularConflicts(
			int [][][] conflicts,
			Conflicts conflict_stats, // to be updated after applying resolution
			double     maxEigen, // maximal eigenvalue of planes to consider
			int        starSteps, // How far to look around when calculationg connection cost
			double     orthoWeight,
			double     diagonalWeight,
			double     starPwr, // Divide cost by number of connections to this power
			double     starWeightPwr,    // Use this power of tile weight when calculating connection cost
			double     weightToDens,    // Balance weighted density against density. 0.0 - density, 1.0 - weighted density
			double     starValPwr, //  Raise value of each tile before averaging
			double     dblTriLoss, //  When resolving double triangles allow minor degradation (0.0 - strict)
			boolean    preferDisparity,
			int        debugLevel,
			int        dbg_X,
			int        dbg_Y)
	{
		final int tilesX =        tileProcessor.getTilesX();
		final int tilesY =        tileProcessor.getTilesY();
		final int superTileSize = tileProcessor.getSuperTileSize();
		//				final int tileSize =      tileProcessor.getTileSize();
		final int stilesX =       (tilesX + superTileSize -1)/superTileSize;
		final int stilesY =       (tilesY + superTileSize -1)/superTileSize;
		final int dbgTile = dbg_Y * stilesX + dbg_X;
		final TileNeibs tnSurface = new TileNeibs(stilesX, stilesY);
		final Conflicts iconflicts = new Conflicts(this);

		int [] rslt = {0,0};
		for (int nsTile = 0; nsTile < conflicts.length; nsTile++) if (conflicts[nsTile] != null) {
			// conflicts may disappear after being fixed, recheck for null
			for (int nConfl = 0; (conflicts[nsTile] != null) && (nConfl < conflicts[nsTile].length); nConfl++){
				Conflict confl = new Conflict(nsTile, conflicts[nsTile][nConfl]); // getDualTriOrthoDiagOrthoConflicts
//				if (numDualTri(conflicts[nsTile][nConfl][2]) > 0) {
				if (confl.getDualTriOrthoDiagOrthoConflicts() > 0) {
					for (int dir4 = 0; dir4 < 4; dir4 ++ ){
						// conflicts may disappear after being fixed, recheck for null
						if ((conflicts[nsTile] != null) && (conflicts[nsTile][nConfl] != null) && ((conflicts[nsTile][nConfl][2] & (1 << dir4)) > 0) && ((conflicts[nsTile][nConfl][2] & (16 << dir4)) > 0)) {
							if (debugLevel > 0) { //1){
								System.out.println("resolveDualTriangularConflicts(): resolving dual triangular conflict for tile "+nsTile+
										", nl1 = "+conflicts[nsTile][nConfl][0]+
										", nl2 = "+conflicts[nsTile][nConfl][1]+
										", dir4 = "+dir4);
								if (nsTile == dbgTile){
									System.out.println("resolveDualTriangularConflicts(): nsTile == dbgTile");
								}
							}
							boolean [] dual_resol = resolveDualTriangularConflict(
									nsTile,
									conflicts[nsTile][nConfl][0],
									conflicts[nsTile][nConfl][1],
									dir4,
									maxEigen, // maximal eigenvalue of planes to consider
									starWeightPwr,    // Use this power of tile weight when calculating connection cost
									tnSurface,
									preferDisparity,
									debugLevel + ((nsTile == dbgTile)? 1 : 0));
							if (dual_resol == null){
								if (debugLevel > 1) {
									System.out.println("resolveDualTriangularConflicts(): resolving dual triangular conflict for tile "+nsTile+
											", nl1 = "+conflicts[nsTile][nConfl][0]+
											", nl2 = "+conflicts[nsTile][nConfl][1]+
											", dir4 = "+dir4+"FAILED");
									rslt[1]++;
								}
							} else {
								if (debugLevel > 1) {
									System.out.println("resolveDualTriangularConflicts(): resolving dual triangular conflict for tile "+nsTile+
											", nl1 = "+conflicts[nsTile][nConfl][0]+
											", nl2 = "+conflicts[nsTile][nConfl][1]+
											", dir4 = "+dir4+" RESOLVED:");
									if (dual_resol[0])	System.out.println("swap links in dir4 = "+dir4);
									if (dual_resol[1])	System.out.println("swap links in dir4 = "+((dir4 + 1) % 4));
									if (dual_resol[2])	System.out.println("swap links in a digonal between dir4 = "+dir4+" and "+((dir4 + 1) % 4));
								}
								// calculate cost of swap
								int [] mod_supertiles = {
										nsTile,
										tnSurface.getNeibIndex(nsTile, 2 * dir4),
										tnSurface.getNeibIndex(nsTile, 2 * ((dir4 + 1) % 4))};
								int nl1 = conflicts[nsTile][nConfl][0];
								int nl2 = conflicts[nsTile][nConfl][1];

								ConnectionCosts connectionCosts = new ConnectionCosts(
										orthoWeight,
										diagonalWeight,
										starPwr,    // Divide cost by number of connections to this power
										starWeightPwr,    // Use this power of tile weight when calculating connection cost
										weightToDens,    // // Balance weighted density against density. 0.0 - density, 1.0 - weighted density
										starValPwr,      //double     starValPwr, //  Raise value of each tile before averaging
										starSteps,
										this.planes,
										tnSurface,
										preferDisparity);
								int [][][] neibs_prev = connectionCosts.initConnectionCosts(mod_supertiles, debugLevel);

								int [][][] neibs_prev_old = new int [mod_supertiles.length][][];
								double [][][]  val_weights = new double [mod_supertiles.length][][];
								updateConnectionsCost_old (
										mod_supertiles,      // int []         nsTiles,
										null,         // int [][][]     neibs_prev,
										neibs_prev_old,   // int [][][]     neibs,		// should be initialized at top dimension if neibs_prev==null
										val_weights,  // double [][][]  val_weights, // should be initialized at top dimension if neibs_prev==null
										orthoWeight,
										diagonalWeight,
										starPwr,        // double     starPwr, // Divide cost by number of connections to this power
										starWeightPwr,    // Use this power of tile weight when calculating connection cost
										tnSurface,
										preferDisparity,
										debugLevel);

								int [][][] neibs = neibs_prev.clone();
								for (int i = 0; i < neibs.length; i++) if (neibs_prev[i] != null){
									neibs[i] = neibs_prev[i].clone();
									for (int j = 0; j < neibs[i].length; j++) if (neibs_prev[i][j] != null){
										neibs[i][j] = neibs_prev[i][j].clone();
									}
								}
								//apply swaps:
//								int [] indices = {tile_index.get
								int dir1f = 2 * dir4;
								int dir1r = (2 * dir4 + 4) % 8;
								int dir2f = (2 * dir4 + 3) % 8;
								int dir2r = (2 * dir4 + 7) % 8;
								int dir3f = (2 * dir4 + 6) % 8;
								int dir3r = (2 * dir4 + 2) % 8;
								int nl1_1 = neibs_prev[0][nl1][dir1f]; // layer on the first that was connected from nl1->N
								int nl2_1 = neibs_prev[0][nl2][dir1f]; // layer on the first that was connected from nl2->N
								int nl1_2 = neibs_prev[1][nl1_1][dir2f]; // layer on second that was connected from nl1->N->SE
								int nl2_2 = neibs_prev[1][nl2_1][dir2f]; // layer on second that was connected from nl2->N->SE
								//  2->0 was already swapped

								if (dual_resol[0])	{
									if (debugLevel > 0) System.out.println("swapping links in dir4 = "+dir4);
									neibs[0][nl1][dir1f] = neibs_prev[0][nl2][dir1f];
									neibs[0][nl2][dir1f] = neibs_prev[0][nl1][dir1f];
									neibs[1][nl1_1][dir1r] = neibs_prev[1][nl2_1][dir1r];
									neibs[1][nl2_1][dir1r] = neibs_prev[1][nl1_1][dir1r];

								}
								if (dual_resol[1])	{
									if (debugLevel > 0) System.out.println("swapping links in dir4 = "+((dir4 + 1) % 4));
									neibs[0][nl1][dir3r] = neibs_prev[0][nl2][dir3r];
									neibs[0][nl2][dir3r] = neibs_prev[0][nl1][dir3r];
									neibs[2][nl1_2][dir3f] = neibs_prev[2][nl2_2][dir3f];
									neibs[2][nl2_2][dir3f] = neibs_prev[2][nl1_2][dir3f];
								}
								if (dual_resol[2]){
									if (debugLevel > 0) System.out.println("swapping links in a digonal between dir4 = "+dir4+" and "+((dir4 + 1) % 4));
									neibs[1][nl1_1][dir2f] = neibs_prev[1][nl2_1][dir2f];
									neibs[1][nl2_1][dir2f] = neibs_prev[1][nl1_1][dir2f];

									neibs[2][nl1_2][dir2r] = neibs_prev[2][nl2_2][dir2r];
									neibs[2][nl2_2][dir2r] = neibs_prev[2][nl1_2][dir2r];
								}

								double cost_diff = 	connectionCosts.getConnectionsCostDiff(
										neibs,
										debugLevel);

								if (debugLevel > 0) {
									System.out.println("resolveDualTriangularConflicts(): resolving dual triangular conflict for tile "+nsTile+
											", nl1 = "+conflicts[nsTile][nConfl][0]+
											", nl2 = "+conflicts[nsTile][nConfl][1]+
											", dir4 = "+dir4+" improvement (negative diff) = "+cost_diff);
								}

								double cost_diff_old = 	updateConnectionsCost_old (
										mod_supertiles, // int []         nsTiles,
										neibs_prev,     // int [][][]     neibs_prev,
										neibs,          // int [][][]     neibs,		// should be initialized at top dimension if neibs_prev==null
										val_weights,    // double [][][]  val_weights, // should be initialized at top dimension if neibs_prev==null
										orthoWeight,
										diagonalWeight,
										starPwr,        // double     starPwr, // Divide cost by number of connections to this power
										starWeightPwr,    // Use this power of tile weight when calculating connection cost
										tnSurface,
										preferDisparity,
										debugLevel);

								if (debugLevel > 0) {
									System.out.println("resolveDualTriangularConflicts()_OLD: resolving dual triangular conflict for tile "+nsTile+
											", nl1 = "+conflicts[nsTile][nConfl][0]+
											", nl2 = "+conflicts[nsTile][nConfl][1]+
											", dir4 = "+dir4+" improvement (negative diff) = "+cost_diff_old);
								}
								if (cost_diff < dblTriLoss) {
									if (debugLevel > 1) {
										if (cost_diff > 0 ) {
											System.out.println("resolveDualTriangularConflicts(): Score is WORSE, but degrading is less than "+
													dblTriLoss+", considered OK.");
										} else {
											System.out.println("resolveDualTriangularConflicts(): SUCCESS");
										}
									}
									// See how this application will influence number of conflicts
									// All supertiles that may have different conflicts
									int [] nsTiles = getInvolvedSupertiles( // first 3 entries will be mod_supertiles[]
											mod_supertiles,
											tnSurface);
									HashMap<Integer,Integer> replacement_tiles = new HashMap<Integer,Integer>();
									for (int i = 0; i < mod_supertiles.length; i++){
										replacement_tiles.put(mod_supertiles[i], new Integer(i));
									}
									int [][][] conflicts_old = new int [nsTiles.length][][];
									int [][][] conflicts_new = new int [nsTiles.length][][];
									for (int isTile = 0; isTile < nsTiles.length; isTile++){
										conflicts_old[isTile] = iconflicts.detectTriangularTileConflicts(
												nsTiles[isTile],   // int nsTile0,
												replacement_tiles, //HashMap<Integer,Integer> replacement_tiles, //
												neibs_prev, // int [][][] replacement_neibs,
												tnSurface); // TileNeibs tnSurface)
										conflicts_new[isTile] = iconflicts.detectTriangularTileConflicts(
												nsTiles[isTile],   // int nsTile0,
												replacement_tiles, //HashMap<Integer,Integer> replacement_tiles, //
												neibs, // int [][][] replacement_neibs,
												tnSurface); // TileNeibs tnSurface)
									}
									if (debugLevel > 1) {
										System.out.println("Involved supertiles:");
										for (int i = 0; i < nsTiles.length; i++){
											System.out.println(i+":"+nsTiles[i]);
										}
									}

									if (debugLevel > 1) {
										System.out.println("Calculating original conflicts");
									}

									iconflicts.addConflicts(conflicts_old,
											debugLevel - 1); // debugLevel);
									if (debugLevel > 1) {
										System.out.println("Calculating resolved conflicts");
									}
									Conflicts conflict_stats_new = new Conflicts(
											conflicts_new,
											this,
											debugLevel - 1); // debugLevel);
									conflict_stats_new.subConflicts(iconflicts);

									if (debugLevel > 0) {
										/*
										System.out.print("Conflicts difference after resolution:");
										printConflictSummary(conflict_stats_diff);
										*/
										conflict_stats_new.printConflictSummary(
												"Conflicts difference after resolution:", true, true, false);

									}
									if ((conflict_stats_new.getNumOrthoIncompat() >= 0) && (cost_diff > -dblTriLoss)){
										//getNumOrthoIncompat
										if (debugLevel > 1) {
											System.out.println("Number of incompatible triangles is not reduced and no significant cost reduction, resolution WILL NOT BE APPLIED");
										}
										rslt[1]++;
									} else {
										if (debugLevel > 0) {
											if (conflict_stats_new.getNumOrthoIncompat() < 0) {
												System.out.println("Number of incompatible triangles is reduced, resolution WILL BE APPLIED");
											} else {
												System.out.println("Number of incompatible triangles is not reduced but a good cost reduction, resolution WILL BE APPLIED");
											}
										}
										// update statistics
										conflict_stats.addConflicts(conflict_stats_new);

										// update conflict
										for (int i = 0; i < nsTiles.length; i++){
											conflicts[nsTiles[i]]=  conflicts_new[i];
										}
										// apply resolution
										for (int i = 0; i < mod_supertiles.length; i++){
											for (int nl = 0; nl < neibs[i].length; nl ++) if (neibs[i][nl] != null){
												planes[mod_supertiles[i]][nl].setNeibBest(neibs[i][nl]);
											}
										}
										rslt[0]++;
									}
								} else {
									if (debugLevel > 0) {
										System.out.println("resolveDualTriangularConflicts(): FAILURE");
									}
									rslt[1]++;
								}
							}
						}
					}
				}
			}
		}
		return rslt;
	}


	/**
	 *  Resolve dual triangle, where NL1 -> N -> SE -> W -> NL2 -> N -> SE -> W -> NL1 (or similar)
	 * @param nsTile supertile index
	 * @param nl1 first conflicting layer (lower of the too
	 * @param nl2 second conflicting layer
	 * @param dir4 start direction to (turn right)  0 - N, 1 - E, 2 - S, 3 - W
	 * @param maxEigen maximal sum of 2 eigenvalues of the tri-plane combinations
	 * @param tnSurface TileNeibs instance to navigate tile index and control array borders
	 * @param preferDisparity - the first eigenvalue/vector is the most disparity-like
	 *                          (false - smallest eigenvalue)
	 * @param debugLevel
	 * @return links to swap : [0] first in dir, [1] - first in dir + 1, [2] - other diagonal. Or null if failed
	 * total number of bits set is odd
	 */
	public boolean [] resolveDualTriangularConflict(
			int nsTile,
			int nl1,
			int nl2,
			int dir4,
			double maxEigen, // maximal eigenvalue of planes to consider
			double     starWeightPwr,    // Use this power of tile weight when calculating connection cost
			TileNeibs tnSurface,
			boolean preferDisparity,
			int debugLevel)
	{
		int [] neibs1 = planes[nsTile][nl1].getNeibBest();
		int [] neibs2 = planes[nsTile][nl2].getNeibBest();
		int dir1 = 2 * dir4;
		int dir2 = (2 * dir4 + 2) % 8;
		int nsTile1 = tnSurface.getNeibIndex(nsTile, dir1);
		int nsTile2 = tnSurface.getNeibIndex(nsTile, dir2);

		TilePlanes.PlaneData [][] tri_pd = {
				{planes[nsTile][nl1], planes[nsTile][nl2]},
				{planes[nsTile1][neibs1[dir1]], planes[nsTile1][neibs2[dir1]]},
				{planes[nsTile2][neibs1[dir2]], planes[nsTile2][neibs2[dir2]]}};
		TilePlanes.PlaneData [][] tri_merged = new TilePlanes.PlaneData [4][2]; // keep to the very end for debugging
		// initial connections: [0][0] -> [1][0] -> [2][1] -> [0][1] -> [1][1]-> [2][0] -> [0][0]
		// 4 variants for 3-supertile sets. get with smallest sum of the 2 main eigenvalues, as long as each
		// is below maxEigen
		if (debugLevel > 1) {
			System.out.println ("resolveDualTriangularConflict() + nsTile = "+nsTile+", nl1 = "+nl1+", nl2 = "+nl2+" dir4="+dir4);
		}
		double [][] ev_variants = new double [4][2];
		int best_pair = -1;
		double best_value = Double.NaN;
		for (int variant = 0; variant < 4; variant++){
			TilePlanes.PlaneData [][] tri_sets = {
					{tri_pd[0][0], tri_pd[1][    (variant & 1)], tri_pd[2][    ((variant >> 1) & 1)]},
					{tri_pd[0][1], tri_pd[1][1 - (variant & 1)], tri_pd[2][1 - ((variant >> 1) & 1)]}};
			for (int nt = 0; nt < tri_sets.length; nt++){
				TilePlanes.PlaneData tri_plane = tri_sets[nt][0];
				for (int indx = 1; indx <  tri_sets[nt].length; indx++){
					TilePlanes.PlaneData other_plane = tri_sets[nt][0].getPlaneToThis(  // layer here does not matter
							tri_sets[nt][indx],
							debugLevel - 2); // debugLevel);
					tri_plane = tri_plane.mergePlaneToThis(
							other_plane,     // PlaneData otherPd,
							1.0,             // double    scale_other,
							starWeightPwr,    // Use this power of tile weight when calculating connection cost
							false,           // boolean   ignore_weights,
							true,            // boolean   sum_weights,
							preferDisparity,
							debugLevel - 2); // int       debugLevel)
				}
//				ev_variants[variant][nt] = tri_plane.getValue();
				ev_variants[variant][nt] = tri_plane.getNormValue();
				tri_merged[variant][nt] = tri_plane;
			}
			if ((best_pair < 0) || ((ev_variants[variant][0] + ev_variants[variant][1]) < best_value)){
				best_value = ev_variants[variant][0] + ev_variants[variant][1];
				best_pair = variant;
			}
			if (debugLevel > 1) {
				System.out.println ("variant="+variant+" value1="+ev_variants[variant][0]+
						" value2="+ev_variants[variant][1]+" sum value = "+(ev_variants[variant][0]+ev_variants[variant][1]));
			}
		}
		if (debugLevel > 1) {
			System.out.println ("resolveDualTriangularConflict() + nsTile = "+nsTile+", nl1 = "+nl1+", nl2 = "+nl2+" dir4="+dir4+" best variant ="+ best_pair);
		}
		if (best_value < maxEigen) {
			boolean [] rslt = {
					((best_pair & 1) != 0),
					((best_pair & 2) != 0),
					false};
			rslt[2] = !rslt[0] ^ rslt[1];
			return rslt;
		}
		return null;
	}

	/**
	 * Calculate main eigenvalue of the current plane and all connected ones - used to estimate advantage of connection swap
	 * @param nsTile supertile index
	 * @param nl surface layer
	 * @param neibs array of 8 neighbors layers (N,NE,...NW), -1 - not connected
	 * @param orthoWeight multiply contribution of ortho neighbors
	 * @param diagonalWeight  multiply contribution of diagonal neighbors
	 * @param diagonalWeight  divide value by number of connections to this power (if !=0)
	 * @param tnSurface TileNeibs instance to navigate tile index and control array borders
	 * @param preferDisparity - the first eigenvalue/vector is the most disparity-like
	 *                          (false - smallest eigenvalue)
	 * @param debugLevel
	 * @return a pair of eigenvalue of the combine plane and its weight
	 */

	public double [] getStarValueWeight_old(
			int    nsTile,
			int    nl,
			int [] neibs,
			double orthoWeight,
			double diagonalWeight,
			double starPwr,    // Divide cost by number of connections to this power
			double     starWeightPwr,    // Use this power of tile weight when calculating connection cost
			TileNeibs tnSurface,
			boolean preferDisparity,
			int    debugLevel)
	{
		TilePlanes.PlaneData merged_plane = planes[nsTile][nl];
		for (int dir = 0; dir < 8; dir++){
			if (neibs[dir] >= 0){
				double other_weight = ((dir & 1) != 0) ? diagonalWeight : orthoWeight;
				TilePlanes.PlaneData other_plane = merged_plane.getPlaneToThis(  // layer here does not matter
						planes[tnSurface.getNeibIndex(nsTile, dir)][neibs[dir]],
						debugLevel - 1); // debugLevel);
				merged_plane = merged_plane.mergePlaneToThis(
						other_plane,     // PlaneData otherPd,
						other_weight,    // double    scale_other,
						starWeightPwr,   // Use this power of tile weight when calculating connection cost
						false,           // boolean   ignore_weights,
						true,            // boolean   sum_weights,
						preferDisparity,
						debugLevel - 1); // int       debugLevel)
			}
		}
//		double [] value_weight = {merged_plane.getValue(),merged_plane.getWeight()};
		double [] value_weight = {merged_plane.getNormValue(),merged_plane.getWeight()};
		if (starPwr != 0){
			value_weight[0] /= (Math.pow((planes[nsTile][nl].getNumNeibBest() + 1.0), starPwr));
		}
		return value_weight;
	}

	/**
	 * Calculate array of supertile indices that need to have connection cost recalculated when they are updated
	 * first entries of the result will be the same in input array
	 * @param mod_supertiles array of supertile indices that will be modified
	 * @param tnSurface TileNeibs instance to navigate tile index and control array borders
	 * @return array of supertile indices to watch connection cost
	 */


	public int [] getInvolvedSupertiles(
			int [] mod_supertiles,
			TileNeibs tnSurface)
	{
		HashSet<Integer> stiles_set = new HashSet<Integer>();
		for (int i = 0; i < mod_supertiles.length; i++){
			stiles_set.add(new Integer(mod_supertiles[i]));
			for (int dir = 0; dir < 8; dir++){
				Integer nsTile = tnSurface.getNeibIndex(mod_supertiles[i], dir);
				if (nsTile >= 0) stiles_set.add (nsTile);
			}
		}
		int [] stiles = new int [stiles_set.size()];
		int indx = 0;
		for (; indx < mod_supertiles.length; indx++){
			stiles_set.remove(new Integer(mod_supertiles[indx]));
			stiles[indx] = mod_supertiles[indx];
		}
		for (Integer nsTile: stiles_set){
			stiles[indx++] = nsTile;
		}
		return stiles;
	}

	/**
	 * Update cost of involved connections after swapping some of them. Total number of links should probably stay
	 *  the same.
	 * @param nsTiles array of indices of the supertiles to watch
	 * @param neibs_prev 3-d array [per supertile index][per layer][per direction] of previous connections for which
	 * the cost was already calculated. Or null - in that case established connections will be used and returned in
	 * neibs array
	 * @param neibs new connections to calculate cost for. If neibs_prev == null, this array should be initialized
	 * by the caller as new int [nsTiles.length][][]
	 * @param val_weights array of values*weights and weights per supertile index, per layer of {val*weight, weight}
	 * @param orthoWeight multiply contribution of ortho neighbors
	 * @param diagonalWeight  multiply contribution of diagonal neighbors
	 * @param tnSurface TileNeibs instance to navigate tile index and control array borders
	 * @param preferDisparity - the first eigenvalue/vector is the most disparity-like
	 *                          (false - smallest eigenvalue)
	 * @param debugLevel
	 * @return difference between new cost and old cost, negative means improvement
	 */


	public double updateConnectionsCost_old (
			int []         nsTiles,
			int [][][]     neibs_prev,
			int [][][]     neibs,		// should be initialized at top dimension if neibs_prev==null
			double [][][]  val_weights, // should be initialized at top dimension if neibs_prev==null
			double         orthoWeight,
			double         diagonalWeight,
			double         starPwr,    // Divide cost by number of connections to this power
			double     starWeightPwr,    // Use this power of tile weight when calculating connection cost
			TileNeibs tnSurface,
			boolean preferDisparity,
			int    debugLevel)
	{
		if (neibs_prev == null){
//			neibs_prev = new int [nsTiles.length][][];
			for (int isTile = 0; isTile < nsTiles.length; isTile++){
				int nsTile = nsTiles[isTile];
				if (planes[nsTile] != null){
					neibs[isTile] = new int [planes[nsTile].length][];
					val_weights[isTile] = new double [planes[nsTile].length][];
					for (int nl = 0; nl < planes[nsTile].length; nl++) if ( planes[nsTile][nl] != null){
						neibs[isTile][nl] =  planes[nsTile][nl].getNeibBest();
						//getNumNeibBest()
						val_weights[isTile][nl] = new double[2];
					}
				}
			}
		}
		// calculate old cost
		double old_value = 0.0;
		double old_weight = 0.0; // should not change during update
		for (int isTile = 0; isTile < nsTiles.length; isTile++){
			if (val_weights[isTile] != null){
				for (int nl = 0; nl < val_weights[isTile].length; nl++) if ( val_weights[isTile][nl] != null){
					old_value += val_weights[isTile][nl][0] * val_weights[isTile][nl][1];
					old_weight += val_weights[isTile][nl][1];
				}
			}
		}
		// now re-calculate val_weights where neibs are different from neibs_prev
		for (int isTile = 0; isTile < nsTiles.length; isTile++){
			int nsTile = nsTiles[isTile];
			if (planes[nsTile] != null) {
				if (neibs[isTile] != null){
					for (int nl = 0; nl < planes[nsTile].length; nl++) if (planes[nsTile][nl] != null){
						if  (neibs[isTile][nl] != null) {
							boolean neibs_changed = false;
							if ((neibs_prev == null) || (neibs_prev[isTile] == null)  || (neibs_prev[isTile][nl] == null)){
								neibs_changed = true;
							} else {
								for (int dir = 0; dir < 8; dir++) if (neibs[isTile][nl][dir] != neibs_prev[isTile][nl][dir]){
									neibs_changed = true;
									break;
								}
							}
							if (neibs_changed){
								val_weights[isTile][nl] = getStarValueWeight_old(
										nsTile,
										nl,
										neibs[isTile][nl],
										orthoWeight,
										diagonalWeight,
										starPwr, // double         starPwr,    // Divide cost by number of connections to this power
										starWeightPwr,    // Use this power of tile weight when calculating connection cost
										tnSurface,
										preferDisparity,
										-1); // debugLevel);
							}
						} else {
							val_weights[isTile][nl] = null;
						}
					}
				} else {
					val_weights[isTile] = null;
				}
			}
		}
		// calculate new cost
		double new_value = 0.0;
		double new_weight = 0.0; // should not change during update
		for (int isTile = 0; isTile < nsTiles.length; isTile++){
			if (val_weights[isTile] != null){
				for (int nl = 0; nl < val_weights[isTile].length; nl++) if ( val_weights[isTile][nl] != null){
					new_value += val_weights[isTile][nl][0] * val_weights[isTile][nl][1];
					new_weight += val_weights[isTile][nl][1];
				}
			}
		}
		// Likely weight will never change except first run, but we will still normalize by weight
		if (old_weight != 0.0) old_value /= old_weight;
		if (new_weight != 0.0) new_value /= new_weight;
		return new_value - old_value; // negative - improvement
	}



	public double [][] getShowPlanes(
			TilePlanes.PlaneData [][] planes,
			double minWeight,
			double maxEigen,
			double dispNorm,
			boolean use_NaN,
			double arrow_dark,
			double arrow_white)
	{
		final int tilesX =        tileProcessor.getTilesX();
		final int tilesY =        tileProcessor.getTilesY();
		final int superTileSize = tileProcessor.getSuperTileSize();
		//				final int tileSize =      tileProcessor.getTileSize();
		final int stilesX = (tilesX + superTileSize -1)/superTileSize;
		final int stilesY = (tilesY + superTileSize -1)/superTileSize;
		final int width =  stilesX * superTileSize; // * tileSize;
		final int height = stilesY * superTileSize; // * tileSize;
		final double [] nan_plane = new double [superTileSize*superTileSize];
		for (int i = 0; i < nan_plane.length; i++) nan_plane[i] = Double.NaN;
		final int centerIndex = (superTileSize+1) * superTileSize / 2;
		final int [] dirs = {-superTileSize,
				-superTileSize + 1,
				1,
				superTileSize +  1,
				superTileSize,
				superTileSize -  1,
				-1,
				-superTileSize - 1};
		//				final int debug_stile = 17 * stilesX + 27;
		final int debug_stile = 20 * stilesX + 27;
		//				final int debug_stile = 19 * stilesX + 27;
		// scan all get maximal number of satisfying planes
		int num_planes = 0;
		for (int i = 0; i < planes.length; i++) {
			if (i == debug_stile) {
				System.out.println("getShowPlanes():1, tile = "+i);
			}
			if ((planes[i] != null) && (planes[i].length > num_planes)){
				int num_sat = 0;
				for (int j = 0; j < planes[i].length; j++) if (planes[i][j] != null){
//					if (planes[i][j].getWeight() >= minWeight){ /// Later does not compare too minWeight
//						double eigVal = planes[i][j].getValue(); // ************** generated does not have value ????????????
						double eigVal = planes[i][j].getNormValue(); // ************** generated does not have value ????????????
						double disp = planes[i][j].getSinglePlaneDisparity(false)[centerIndex];
						/*
								if (disp > dispNorm) {
									eigVal *= dispNorm / disp;
								}
								if (eigVal < maxEigen) num_sat ++;
						 */
						if (eigVal < corrMaxEigen(
								maxEigen,
								dispNorm,
								disp)) num_sat ++;
//					}
				}
				if (num_sat > num_planes) num_planes = num_sat;
			}
		}
		double [][] data = new double[num_planes][width*height];
		for (int sty = 0; sty < stilesY; sty++){
			for (int stx = 0; stx < stilesX; stx++){
				int nsTile = sty * stilesX + stx;
				if (nsTile == debug_stile) {
					System.out.println("getShowPlanes():2, tile = "+nsTile);
				}
				int ns = 0;
				double [] plane = nan_plane;
				int indx = (sty * superTileSize) * width + (stx * superTileSize);
				boolean debug = (nsTile == debug_stile) ; // (sty == 18) && (stx == 28);
				if (debug) {
					System.out.println("getShowPlanes():"+((planes[nsTile] != null)? planes[nsTile].length : "null"));
				}
				if (planes[nsTile] != null) {
//					for (int np = 0; np < planes[nsTile].length; np++) if (planes[nsTile][np] != null){
					for (int np = 0; np < planes[nsTile].length; np++) if ((planes[nsTile][np] != null) && (planes[nsTile][np].getWeight() > 0.0)){ // disregard 0-weight planes

						if (planes[nsTile][np].getWeight() >= minWeight){
							//&& (planes[nsTile][np].getValue() < maxEigen)){
//							double eigVal = planes[nsTile][np].getValue();
							double eigVal = planes[nsTile][np].getNormValue();
							double disp = planes[nsTile][np].getZxy()[0];
							/*
									if (disp > dispNorm) {
										eigVal *= dispNorm / disp;
									}
							 */
							if (eigVal < corrMaxEigen(
									maxEigen,
									dispNorm,
									disp)) {
								plane = planes[nsTile][np].getSinglePlaneDisparity(use_NaN);
								// add connections to neighbors if available
								int [] neib_best = planes[nsTile][np].getNeibBest();
								if (neib_best != null){
									for (int dir = 0; dir < neib_best.length; dir ++) if (neib_best[dir] >= 0){
										// draw a line
										//int center_index = indx + centerIndex;
										for (int i = 0; i < (superTileSize / 2); i++){
											plane[centerIndex + i * dirs[dir]] = ((i & 1) == 0) ? arrow_white:arrow_white; //  arrow_dark : arrow_white;
										}
									}
								}
								//superTileSize
								// mark fronto tiles with NaN - center 3x3 tiles
								//superTileSize
								if (planes[nsTile][np].fronto) {
									plane[centerIndex + 0] =             Double.NaN;
									plane[centerIndex + superTileSize] = Double.NaN;
									plane[centerIndex - superTileSize] = Double.NaN;
								}
								if (planes[nsTile][np].horizontal) {
									plane[centerIndex + 0] =             Double.NaN;
									plane[centerIndex + 1] =             Double.NaN;
									plane[centerIndex - 1] =             Double.NaN;
								}
								for (int y = 0; y < superTileSize; y++) {
									if (data.length < (ns+1)){
										System.out.println("BUG in getShowPlanes()!, ns = "+ns+", data.length="+data.length);
									} else {
										System.arraycopy(plane, superTileSize * y, data[ns], indx + width * y, superTileSize); //java.lang.ArrayIndexOutOfBoundsException: 5
									}
								}
								ns ++;
							}
						}

					}
				}
				// fill same (nearest) plane or NaN if none was available
				for (; ns < num_planes; ns++){
					for (int y = 0; y < superTileSize; y++) {
						System.arraycopy(plane, superTileSize * y, data[ns], indx + width * y, superTileSize);
					}
				}

			}
		}
		return data;
	}

	public double [][] getShowShells(
			int  nlayer, // over multi-layer - do not render more than nlayer on top of each other
			TilePlanes.PlaneData [][] planes,
			int [][] shells,
			int max_shells,
			boolean fuse,
			boolean show_connections,
			boolean use_NaN,
			double arrow_dark,
			double arrow_white)
	{

		final int tilesX =        tileProcessor.getTilesX();
		final int tilesY =        tileProcessor.getTilesY();
		final int superTileSize = tileProcessor.getSuperTileSize();
		//				final int tileSize =      tileProcessor.getTileSize();
		final int stilesX = (tilesX + superTileSize -1)/superTileSize;
		final int stilesY = (tilesY + superTileSize -1)/superTileSize;
		final int width =  stilesX * superTileSize; // * tileSize;
		final int height = stilesY * superTileSize; // * tileSize;
		final double [] nan_plane = new double [superTileSize*superTileSize];
		for (int i = 0; i < nan_plane.length; i++) nan_plane[i] = Double.NaN;
		final int centerIndex = (superTileSize+1) * superTileSize / 2;
		final int [] dirs = {-superTileSize,
				-superTileSize + 1,
				1,
				superTileSize +  1,
				superTileSize,
				superTileSize -  1,
				-1,
				-superTileSize - 1};
		final int [] neib_offsets={
				-stilesX,
				-stilesX + 1,
				1,
				stilesX +  1,
				stilesX,
				stilesX -  1,
				-1,
				-stilesX - 1};
		//				final int debug_stile = 17 * stilesX + 27;
		final int debug_stile = -1; // 18 * stilesX + 30;

		int num_shells = 0;
		for (int nsTile = 0; nsTile < shells.length; nsTile++) if (shells[nsTile] != null){
			for (int np = 0; np < shells[nsTile].length; np++){
				if (shells[nsTile][np] > num_shells) {
					num_shells = shells[nsTile][np]; // starts with 1, not 0
				}
			}
		}
		if (num_shells > max_shells){
			num_shells = max_shells;
		}
		double [][] data = new double [num_shells][width*height];
		for (int nShell = 0; nShell < num_shells; nShell++) {
			for (int sty = 0; sty < stilesY; sty++){
				for (int stx = 0; stx < stilesX; stx++){
					int nsTile = sty * stilesX + stx;
					boolean debug = (nsTile == debug_stile) ; // (sty == 18) && (stx == 28);
					String [] dbg_titles = {"Center","N","NE","E","SE","S","SW","W","NW"};
					double [][] dbg_img3 = null;
					double [][] dbg_img = null;
					double [][] dbg_this = null;
					double [][] dbg_other = null;
					if (debug) {
						System.out.println("getShowShells():2, tile = "+nsTile);
					}
					double [] plane = nan_plane;
					int indx = (sty * superTileSize) * width + (stx * superTileSize);
					if (debug) {
						System.out.println("getShowShells():"+((planes[nsTile] != null)? planes[nsTile].length : "null"));
					}
					if (planes[nsTile] != null) {
						int num_rendered = 0;
						for (int np = 0; np < planes[nsTile].length; np++){
							if ((shells[nsTile][np] == (nShell + 1)) && ((nlayer == 0) || (num_rendered < nlayer)) ) {
								num_rendered ++;
								plane = planes[nsTile][np].getSinglePlaneDisparity(use_NaN);
								int [] neib_best = planes[nsTile][np].getNeibBest();
								if (fuse && (neib_best != null)){
									if (debug) {
										dbg_img3 = new double [dbg_titles.length][];
										dbg_img = new double [dbg_titles.length][];
										dbg_this = new double [dbg_titles.length][];
										dbg_other = new double [dbg_titles.length][];
										dbg_img3[0] = planes[nsTile][np].getTriplePlaneDisparity();
										dbg_img[0] = planes[nsTile][np].getSinglePlaneDisparity(false);
									}


									double [] weight_this =  new double [superTileSize*superTileSize];
									double [] weight_other = new double [superTileSize*superTileSize];
									double [] plane_combo =  new double [superTileSize*superTileSize];
									for (int dir = 0; dir < 8; dir++){
										if (neib_best[dir] >=0){
											int nsTileOther = nsTile + neib_offsets[dir];
											int npOther = neib_best[dir];
											int dirOther = (dir + 4) % 8;
											double [] plane_other =	planes[nsTileOther][npOther].getTriplePlaneDisparity(dirOther);
											double [] fuse_this =  getFuseThis(dir);
											double [] fuse_other = getFuseOther(dir);
											if (debug) {
												dbg_img3[dir + 1] =   planes[nsTileOther][npOther].getTriplePlaneDisparity();
												dbg_img[dir + 1] =    planes[nsTileOther][npOther].getTriplePlaneDisparity(dirOther);
												dbg_this[dir + 1] =   fuse_this.clone();
												dbg_other[dir + 1] =  fuse_other.clone();
											}
											for (int i = 0; i < weight_this.length; i++){
												weight_this[i] +=  fuse_this[i];
												weight_other[i] += fuse_other[i];
												plane_combo[i] +=  fuse_other[i] * plane_other[i];
											}
										}
									}
									if (debug) {
										ShowDoubleFloatArrays sdfa_instance = new ShowDoubleFloatArrays(); // just for debugging?
										sdfa_instance.showArrays(dbg_img, superTileSize, superTileSize, true, "plane_"+nsTile,dbg_titles);
										sdfa_instance.showArrays(dbg_img3, 3*superTileSize, 3*superTileSize, true, "plane3_"+nsTile,dbg_titles);
										sdfa_instance.showArrays(dbg_this, superTileSize, superTileSize, true, "weights_this_"+nsTile,dbg_titles);
										sdfa_instance.showArrays(dbg_other, superTileSize, superTileSize, true, "weights_other_"+nsTile,dbg_titles);
									}
									double [] plane_other = planes[nsTile][np].getSinglePlaneDisparity(false);
									for (int i = 0; i < weight_this.length; i++){
										plane_combo[i] += weight_this[i] * plane_other[i];
										double w = weight_this[i] + weight_other[i];
										if (w > 0){
											plane_combo[i] /= w;
										} else {
											plane_combo[i] = plane_other[i];
										}
									}
									if (use_NaN){
										for (int i = 0; i < plane_combo.length; i++){
											if (Double.isNaN(plane[i])){
												plane_combo[i] = Double.NaN;
											}
										}
									}
									plane = plane_combo;
								}

								if (show_connections) {
									if (neib_best != null){
										for (int dir = 0; dir < neib_best.length; dir ++) if (neib_best[dir] >= 0){
											// draw a line
											//int center_index = indx + centerIndex;
											for (int i = 0; i < (superTileSize / 2); i++){
												plane[centerIndex + i * dirs[dir]] = ((i & 1) == 0) ? arrow_dark : arrow_white;
											}
										}
									}
								}
								for (int y = 0; y < superTileSize; y++) {
									System.arraycopy(plane, superTileSize * y, data[nShell], indx + width * y, superTileSize);
								}
							}
						}
					}
				}
			}
		}
		return data;
	}


	public double [][][] getOverlappingShells(
			int  nlayer, // over multi-layer - do not render more than nlayer on top of each other
			TilePlanes.PlaneData [][] planes,
			int [][] shells,
			int max_shells,
			boolean fuse,
			boolean show_connections,
			boolean use_NaN,
			double arrow_dark,
			double arrow_white)
	{

		final int tilesX =        tileProcessor.getTilesX();
		final int tilesY =        tileProcessor.getTilesY();
		final int superTileSize = tileProcessor.getSuperTileSize();
		//				final int tileSize =      tileProcessor.getTileSize();
		final int stilesX = (tilesX + superTileSize -1)/superTileSize;
		final int stilesY = (tilesY + superTileSize -1)/superTileSize;
		final int width =  stilesX * superTileSize; // * tileSize;
		final int height = stilesY * superTileSize; // * tileSize;
		final double [] nan_plane = new double [superTileSize*superTileSize];
		for (int i = 0; i < nan_plane.length; i++) nan_plane[i] = Double.NaN;
		final int centerIndex = (superTileSize+1) * superTileSize / 2;
		final int [] dirs = {-superTileSize,
				-superTileSize + 1,
				1,
				superTileSize +  1,
				superTileSize,
				superTileSize -  1,
				-1,
				-superTileSize - 1};
		final int [] neib_offsets={
				-stilesX,
				-stilesX + 1,
				1,
				stilesX +  1,
				stilesX,
				stilesX -  1,
				-1,
				-stilesX - 1};
		//				final int debug_stile = 17 * stilesX + 27;
		final int debug_stile = -1; // 18 * stilesX + 30;

		int num_shells = 0;
		for (int nsTile = 0; nsTile < shells.length; nsTile++) if (shells[nsTile] != null){
			for (int np = 0; np < shells[nsTile].length; np++){
				if (shells[nsTile][np] > num_shells) {
					num_shells = shells[nsTile][np]; // starts with 1, not 0
				}
			}
		}
		if (num_shells > max_shells){
			num_shells = max_shells;
		}
		double [][] data = new double [num_shells][width*height];
		for (int nShell = 0; nShell < num_shells; nShell++) {
			for (int sty = 0; sty < stilesY; sty++){
				for (int stx = 0; stx < stilesX; stx++){
					int nsTile = sty * stilesX + stx;
					boolean debug = (nsTile == debug_stile) ; // (sty == 18) && (stx == 28);
					String [] dbg_titles = {"Center","N","NE","E","SE","S","SW","W","NW"};
					double [][] dbg_img3 = null;
					double [][] dbg_img = null;
					double [][] dbg_this = null;
					double [][] dbg_other = null;
					if (debug) {
						System.out.println("getShowShells():2, tile = "+nsTile);
					}
					double [] plane = nan_plane;
					int indx = (sty * superTileSize) * width + (stx * superTileSize);
					if (debug) {
						System.out.println("getShowShells():"+((planes[nsTile] != null)? planes[nsTile].length : "null"));
					}
					if (planes[nsTile] != null) {
						int num_rendered = 0;
						for (int np = 0; np < planes[nsTile].length; np++){
							if ((shells[nsTile][np] == (nShell + 1)) && ((nlayer == 0) || (num_rendered < nlayer)) ) {
								num_rendered ++;
								plane = planes[nsTile][np].getSinglePlaneDisparity(use_NaN);
								int [] neib_best = planes[nsTile][np].getNeibBest();
								if (fuse && (neib_best != null)){
									if (debug) {
										dbg_img3 = new double [dbg_titles.length][];
										dbg_img = new double [dbg_titles.length][];
										dbg_this = new double [dbg_titles.length][];
										dbg_other = new double [dbg_titles.length][];
										dbg_img3[0] = planes[nsTile][np].getTriplePlaneDisparity();
										dbg_img[0] = planes[nsTile][np].getSinglePlaneDisparity(false);
									}


									double [] weight_this =  new double [superTileSize*superTileSize];
									double [] weight_other = new double [superTileSize*superTileSize];
									double [] plane_combo =  new double [superTileSize*superTileSize];
									for (int dir = 0; dir < 8; dir++){
										if (neib_best[dir] >=0){
											int nsTileOther = nsTile + neib_offsets[dir];
											int npOther = neib_best[dir];
											int dirOther = (dir + 4) % 8;
											double [] plane_other =	planes[nsTileOther][npOther].getTriplePlaneDisparity(dirOther);
											double [] fuse_this =  getFuseThis(dir);
											double [] fuse_other = getFuseOther(dir);
											if (debug) {
												dbg_img3[dir + 1] =   planes[nsTileOther][npOther].getTriplePlaneDisparity();
												dbg_img[dir + 1] =    planes[nsTileOther][npOther].getTriplePlaneDisparity(dirOther);
												dbg_this[dir + 1] =   fuse_this.clone();
												dbg_other[dir + 1] =  fuse_other.clone();
											}
											for (int i = 0; i < weight_this.length; i++){
												weight_this[i] +=  fuse_this[i];
												weight_other[i] += fuse_other[i];
												plane_combo[i] +=  fuse_other[i] * plane_other[i];
											}
										}
									}
									if (debug) {
										ShowDoubleFloatArrays sdfa_instance = new ShowDoubleFloatArrays(); // just for debugging?
										sdfa_instance.showArrays(dbg_img, superTileSize, superTileSize, true, "plane_"+nsTile,dbg_titles);
										sdfa_instance.showArrays(dbg_img3, 3*superTileSize, 3*superTileSize, true, "plane3_"+nsTile,dbg_titles);
										sdfa_instance.showArrays(dbg_this, superTileSize, superTileSize, true, "weights_this_"+nsTile,dbg_titles);
										sdfa_instance.showArrays(dbg_other, superTileSize, superTileSize, true, "weights_other_"+nsTile,dbg_titles);
									}
									double [] plane_other = planes[nsTile][np].getSinglePlaneDisparity(false);
									for (int i = 0; i < weight_this.length; i++){
										plane_combo[i] += weight_this[i] * plane_other[i];
										double w = weight_this[i] + weight_other[i];
										if (w > 0){
											plane_combo[i] /= w;
										} else {
											plane_combo[i] = plane_other[i];
										}
									}
									if (use_NaN){
										for (int i = 0; i < plane_combo.length; i++){
											if (Double.isNaN(plane[i])){
												plane_combo[i] = Double.NaN;
											}
										}
									}
									plane = plane_combo;
								}

								if (show_connections) {
									if (neib_best != null){
										for (int dir = 0; dir < neib_best.length; dir ++) if (neib_best[dir] >= 0){
											// draw a line
											//int center_index = indx + centerIndex;
											for (int i = 0; i < (superTileSize / 2); i++){
												plane[centerIndex + i * dirs[dir]] = ((i & 1) == 0) ? arrow_dark : arrow_white;
											}
										}
									}
								}
								for (int y = 0; y < superTileSize; y++) {
									System.arraycopy(plane, superTileSize * y, data[nShell], indx + width * y, superTileSize);
								}
							}
						}
					}
				}
			}
		}
		return  null; // data;
	}




	/**
	 * Prepare visualization of the plane separation lines
	 * @param split_planes per supertile, per plane sets of 3 plane data instances
	 *  (first, second, merged)
	 * @param debugLevel debug level
	 * @param dbg_X supertile horizontal index to show debug information
	 * @param dbg_Y supertile vertical index to show debug information
	 * @return array to be visualized as width=superTileSize*stilesX,
	 *  height = superTileSize*stilesY data. Lines value match plane index,
	 *  length corresponds to triple-sized supertiles
	 */

	public double [] showSplitLines(
			final TilePlanes.PlaneData[][][] split_planes,
			final int                        debugLevel,
			final int                        dbg_X,
			final int                        dbg_Y)
	{
		final int [][] dirsYX = {{-1, 0},{-1,1},{0,1},{1,1},{1,0},{1,-1},{0,-1},{-1,-1}};
		final int tilesX =        tileProcessor.getTilesX();
		final int tilesY =        tileProcessor.getTilesY();
		final int superTileSize = tileProcessor.getSuperTileSize();
		final int stilesX = (tilesX + superTileSize -1)/superTileSize;
		final int stilesY = (tilesY + superTileSize -1)/superTileSize;
		final int width = superTileSize * stilesX;
		final int height = superTileSize * stilesY;
		final int debug_stile = dbg_Y * stilesX + dbg_X;
		double [] split_lines = new double[width * height];
		final int  ss3 = 3 * superTileSize;
		for (int nsTile = 0; nsTile < split_planes.length; nsTile++) {
			if (split_planes[nsTile] != null) {
				int sty = nsTile / stilesX;
				int stx = nsTile % stilesX;
				int dl = ((debugLevel > -1) && (nsTile == debug_stile)) ? 1:0;
				int min_np = (split_planes[nsTile].length > 1)? 1 : 0;
				for (int np = min_np; np < split_planes[nsTile].length; np++){
					if (split_planes[nsTile][np] != null){
						TilePlanes.PlaneData split_plane1 = split_planes[nsTile][np][0];
						TilePlanes.PlaneData split_plane2 = split_planes[nsTile][np][1];
						if ((split_plane1 != null) && ( split_plane2 != null)){
							double [] tpd =  split_plane1.getTriplePlaneDisparity();
							double [] tpd1 = split_plane2.getTriplePlaneDisparity();
							// difference of the 2 split planes disparities, show where
							// they have lowest value in any of the 2 orthogonal directions
							for (int i = 0; i < tpd.length; i++) {
								tpd[i] = Math.abs(tpd[i] - tpd1[i]);
							}

							int [] vby = new int [ss3];
							for (int x = 0; x < ss3; x ++){
								vby[x] = 0;
								for (int y = 1; y < ss3; y ++){
									if (tpd[y * ss3 + x] < tpd[vby[x] * ss3 + x]){
										vby[x] = y;
									}
								}
							}
							int [] vbx = new int [ss3];
							for (int y = 0; y < ss3; y ++){
								vbx[y] = 0;
								for (int x = 1; x < ss3; x ++){
									if (tpd[y * ss3 + x] < tpd[y * ss3 + vbx[y]]){
										vbx[y] = x;
									}
								}
							}

							for (int dx = 0; dx < ss3; dx ++){
								if ((vby[dx] > 0) && (vby[dx] < (ss3-1))) {
									int y = vby[dx] + superTileSize * (sty - 1);
									int x = dx + superTileSize * (stx - 1);
									if ((y >= 0) && (y < height) && (x >= 0) && (x < width)){
										split_lines[ y* width + x]	= np - min_np + 1;
									}
								}
							}

							for (int dy = 0; dy < ss3; dy ++){
								if ((vbx[dy] > 0) && (vbx[dy] < (ss3-1))) {
									int x = vbx[dy] + superTileSize * (stx - 1);
									int y = dy + superTileSize * (sty - 1);
									if ((y >= 0) && (y < height) && (x >= 0) && (x < width)){
										split_lines[ y* width + x]	= np - min_np + 1;
									}
								}
							}
						}
					}
				}
			}
		}
		return split_lines;
	}

	/**
	 * Apply same supertile planes merge by combining tiles and re-generating ellipsoids by diagonalizing
	 * covariance matrices. Some outliers may be removed after merge
	 * @param planes per supertile, per plane - array of supertile instances - will be modified =by merge
	 * @param merge_groups per-supertile group sets for merging. Each group set is an array of groups. Each group is an array
	 * of plane indices
	 * Parameters to generate planes (ellipsoids):
	 * @param disp_far disparity lower limit (Double.NaN - any)
	 * @param disp_near disparity upper limit (Double.NaN - any)
	 * @param dispNorm disparity normalization value (when average disparity is above, difference is proportionally reduced)
	 * @param min_weight minimal tile strength to be used
	 * @param min_tiles minimal number of tiles to generate ellipsoid
	 * Parameters for outlier removal:
	 * @param targetEigen target main eigenvalue (thickness in disparity space)
	 * @param fractOutliers maximal fraction of all tiles to be removed as outliers
	 * @param maxOutliers maximal absolute number of outliers to be removed from each plane (ellipsoid)
	 * @param debugLevel debug level
	 * @param dbg_X tile x-index for detailed debug data
	 * @param dbg_Y tile y-index for detailed debug data
	 * @return total number of plane groups merged
	 */
	public int applyMergePlanes(
			final TilePlanes.PlaneData[][]   planes,
			final int [][][]                 merge_groups,
			// parameters to generate ellipsoids
			final double                     disp_far, // minimal disparity to select (or NaN)
			final double                     disp_near, // maximal disparity to select (or NaN)
			final double                     dispNorm,   //  Normalize disparities to the average if above
			final double                     min_weight,
			final int                        min_tiles,
			// parameters to reduce outliers
			final double                     targetEigen,   //     =   0.1;  // Remove outliers until main axis eigenvalue (possibly scaled by plDispNorm) gets below
			final double                     inFractOutliers, //     =   0.3;  // Maximal fraction of outliers to remove
			final int                        maxOutliers,   //     =   20;  // Maximal number of outliers to remove
			final double                     fronto_tol,    // fronto tolerance (pix) - treat almost fronto as fronto (constant disparity). <= 0 - disable this feature
// FIXME: the following 2 parameters are not yet used
			final double                     fronto_rms,    // Target rms for the fronto planes - same as sqrt(plMaxEigen) for other planes
			final double                     fronto_offs,   //        =   0.2;  // increasing weight of the near tiles by using difference between the reduced average as weight. <= 0 - disable
			final double                     fronto_pow,    //        =   1.0;  // increase weight even more
			final int                        debugLevel,
			final int                        dbg_X,
			final int                        dbg_Y)
	{
		final int tilesX =        tileProcessor.getTilesX();
		final int superTileSize = tileProcessor.getSuperTileSize();
		final int stilesX = (tilesX + superTileSize -1)/superTileSize;
		final int debug_stile = dbg_Y * stilesX + dbg_X;
		final Thread[] threads = ImageDtt.newThreadArray((debugLevel > 1)? 1 : tileProcessor.threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);

		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				@Override
				public void run() {
					for (int nsTile = ai.getAndIncrement(); nsTile < planes.length; nsTile = ai.getAndIncrement()) {
//						int dl = ((debugLevel > -1) && (nsTile == debug_stile)) ? 4:0;
                        int dl = ((debugLevel > 1) && (nsTile == debug_stile)) ? 4: debugLevel;

						if (merge_groups[nsTile] != null){
							// first merge all to the lowest plane (they are ordered), then re-order remaining planes
							for (int ng = 0; ng < merge_groups[nsTile].length; ng++) {
								int np0 = merge_groups[nsTile][ng][0];
								 TilePlanes.PlaneData this_pd = planes[nsTile][np0];
								for (int i = 1; i < merge_groups[nsTile][ng].length; i++ ){
									int np = merge_groups[nsTile][ng][i];
									this_pd.orMeasSelection(planes[nsTile][np].getMeasSelection());
									planes[nsTile][np] = null;
								}
								this_pd.invalidateCalculated(); // is it needed?
								double [][][] disp_strength = this_pd.getPlaneFromMeas(
										null,           // boolean [][] tile_sel,null - use measuredSelection already OR-ed
										null,           //double [][][] disp_str, null - calculate
										disp_far,       // double       disp_far, // minimal disparity to select (or NaN)
										disp_near,      // double       disp_near, // maximal disparity to select (or NaN)
										dispNorm,       // double       dispNorm,   //  Normalize disparities to the average if above
										min_weight,     // double       min_weight,
										min_tiles,      // int          min_tiles,
//OK?
										smplMode,
										mlfp,
										dl - 2);            // int          debugLevel)
								if (disp_strength == null) {
									System.out.println("=== BUG in applyMergePlanes(): failed to getPlaneFromMeas() for merged planes");
									break;
								}

								this_pd.getWorldPlaneFromMeas( // re-calculate world-based planes too
										null,           // tile_sel,       // boolean [][] tile_sel, // null - do not use, {} use all (will be modified)
										disp_strength,
										disp_far,       // double       disp_far, // minimal disparity to select (or NaN)
										disp_near,      // double       disp_near, // maximal disparity to select (or NaN)
										dispNorm,       // double       dispNorm,   //  Normalize disparities to the average if above
										min_weight,     // double       min_weight,
										min_tiles,      // int          min_tiles,
//OK?
										smplMode,
										mlfp,
										dl - 2);            // int          debugLevel)



								// remove outliers //removeOutliers
								// now try to remove outliers

								double targetV = corrMaxEigen(
										targetEigen,
										dispNorm,
										planes[nsTile][np0]);

								boolean almost_fronto = this_pd.isFronto(
										fronto_tol,    // double fronto_tol,
										disp_strength,      // double [][][] disp_str,
										debugLevel);   // int debugLevel);
								// seems it is OK not to normalize fronto_rms, as it is only needed for far surfaces
								double targetEigen = almost_fronto ? (fronto_rms*fronto_rms): targetV;
//								if (this_pd.getValue() > targetV) {
								// FIXME: temporary - make a configurable parameter
								double fronto_fract_outliers = 0.8; // 2 * inFractOutliers;
								double fractOutliers = almost_fronto ? (fronto_fract_outliers) : inFractOutliers;
								int max_outliers = (int) Math.round(this_pd.getNumPoints() * fractOutliers);
								if (max_outliers > maxOutliers) max_outliers = maxOutliers;
								// in fronto mode - run at least once
								if (almost_fronto) {
									this_pd.growFrontoSelection(
											5,              // int           min_neib,
											0.1,            // double        max_over_avg,
											0.00,           // double        max_under_avg,
											disp_strength); // double [][][] disp_str)
								}
								if ((this_pd.getNormValue() > targetEigen) || almost_fronto) { // targetV) {
									boolean OK = this_pd.removeOutliers( // getPlaneFromMeas should already have run
											fronto_tol, // fronto tolerance (pix) - treat almost fronto as fronto (constant disparity). <= 0 - disable this feature
											fronto_rms,    // double        fronto_rms,  // Target rms for the fronto planes - same as sqrt(plMaxEigen) for other planes. May be tighter
											fronto_offs,   // double        fronto_offs,   //        =   0.2;  // increasing weight of the near tiles by using difference between the reduced average as weight.
											fronto_pow,    // double        fronto_pow,    //        =   1.0;  // increase weight even more
											disp_strength,
											targetV,      // double     targetEigen, // target eigenvalue for primary axis (is disparity-dependent, so is non-constant)
											max_outliers, // int        maxRemoved,  // maximal number of tiles to remove (not a constant)
											dl); // int        debugLevel)
									if (!OK) {
										System.out.println("=== BUG in applyMergePlanes(): failed to removeOutliers() for merged planes");
										break;
									}
								}
								if (this_pd.fronto) {
									//?
								}
							}
// remove all null planes (but number zero - it may be null for historical reasons)
							ArrayList<Integer> remaining_planes = new ArrayList<Integer>();
							remaining_planes.add(0);
							for (int np = 1; np < planes[nsTile].length; np++){
								if (planes[nsTile][np] != null) {
									remaining_planes.add(np);
								}
							}
							TilePlanes.PlaneData[] stack_pd = new TilePlanes.PlaneData[remaining_planes.size()];
							int indx = 0;
							for (Integer np : remaining_planes){
								stack_pd[indx++] = planes[nsTile][np];
							}
							planes[nsTile] = stack_pd;
						}
					}
				}
			};
		}
		ImageDtt.startAndJoin(threads);
		int num_merged = 0;
		for (int nsTile = 0; nsTile < merge_groups.length; nsTile++) if (merge_groups[nsTile] != null){
			for (int ng = 0; ng < merge_groups[nsTile].length; ng++){
				num_merged += merge_groups[nsTile][ng].length - 1;
			}
		}
		return num_merged;
	}




	/**
	 * Re-generate "planes" (ellipsoids) from the measured data according to provided
	 * plane pairs for some supertiles.
	 * Current version destroys connections between neighbor planes, they need to be
	 * re-calculated
	 * TODO: consider amplifying difference (angle) between planes before/during/after
	 * running this method. Like 3d "sharpening" to compensate for being washed out.
	 * @param planes per-supertile (in linescan order), per-plane data to be replaced
	 * @param brokenPd per-supertile , per-plane 3 PD instances (two first make the pair).
	 *        generated by breakPlanesToPairs method. This data will not be modified.
	 * @param max_diff maximal normalized disparity difference from the plane to consider
	 * @param other_diff maximal difference of the added tile ratio to the average
	 *  disparity difference of the exclusively selected tiles
	 * @param non_exclusive allow tiles to belong to both planes of the pair
	 * @param use_other_planes allow other tiles from the same supertile
	 * @param measSel (with use_other_planes) select measurements for supertiles:
	 *  +1 - combo, +2 - quad +4 - hor +8 - vert
	 * @param allow_parallel allow parallel shift of the specified planes before adding
	 *  more tiles (from the same pair and/or other tiles of the supertile (see
	 *  non_exclusive, use_other_planes
	 * @param splitXY generate XY separation masks for split planes, so tiles will not
	 *        be allowed to snap to the other plane in pair
	 * @param splitXYTolerance XY separation disparity tolerance (there will be some
	 *  overlapping if > 0.0
	 * @param disp_far disparity low limit
	 * @param disp_near disparity high limit (Double.NaN - no limit)
	 * @param dispNorm normalize disparity difference (keep relative) for disparities above
	 * @param min_weight minimal total weight of the plane
	 * @param min_tiles minimal number of tiles in a plane
	 * @param targetEigen target primary eigenvalue when removing outliers
	 * @param fractOutliers maximal fraction of outliers from all tiles to remove
	 * @param maxOutliers maximal number of outliers to remove
	 * @param strength_floor subtract from correlation strength (and limit by 0) to use
	 *  correlation strength as weight
	 * @param strength_pow raise correlation strength (after subtracting strength_floor) to
	 *  this power before using as weight
	 * @param smplMode use square sample mode, false - single-tile samples
	 * @param smplSide size of the square sample side
	 * @param smplNum number of averaged samples (should be <= smplSide * smplSide and > 1)
	 * @param smplRms maximal square root of variance (in disparity pixels) to accept the result
	 * @param debugLevel debug level
	 * @param dbg_X supertile horizontal index to show debug information
	 * @param dbg_Y supertile vertical index to show debug information
	 * @return number of planes replaced
	 */
	// TODO: Mark replaced planes and reduce eigenvalue requirements for them (in case they are small)
	// TODO:  planesSmooth() add power function to modify pull/neighbors relation (8 pull stronger than 1, but not 8 times stronger
	public int replaceBrokenPlanes(
			final TilePlanes.PlaneData[][]   planes,
			final TilePlanes.PlaneData[][][] brokenPd,
			final double                     max_diff, // maximal disparity difference (0 - any), will be normalized by dispNorm
			final double                     other_diff,
			final boolean                    non_exclusive,
			final boolean                    use_other_planes, // TODO:
			final int                        measSel, // Select measurements for supertiles : +1 - combo, +2 - quad +4 - hor +8 - vert
			final boolean                    allow_parallel,

			final boolean                    splitXY,
			final double                     splitXYTolerance,

			// parameters to generate ellipsoids
			final double                     disp_far, // minimal disparity to select (or NaN)
			final double                     disp_near, // maximal disparity to select (or NaN)
			final double                     dispNorm,   //  Normalize disparities to the average if above
			final double                     min_weight,
			final int                        min_tiles,
			// parameters to reduce outliers
			final double                     targetEigen,   //     =   0.1;  // Remove outliers until main axis eigenvalue (possibly scaled by plDispNorm) gets below
			final double                     inFractOutliers, //     =   0.3;  // Maximal fraction of outliers to remove
			final int                        maxOutliers,   //     =   20;  // Maximal number of outliers to remove
			final boolean                    smplMode, //        = true;   // Use sample mode (false - regular tile mode)
			final MeasuredLayersFilterParameters mlfp,
			final double                     fronto_tol, // fronto tolerance (pix) - treat almost fronto as fronto (constant disparity). <= 0 - disable this feature
//FIXME: use following 2 parameters
			final double                     fronto_rms,    // Target rms for the fronto planes - same as sqrt(plMaxEigen) for other planes
			final double                     fronto_offs,   //        =   0.2;  // increasing weight of the near tiles by using difference between the reduced average as weight. <= 0 - disable
			final double                     fronto_pow,    //        =   1.0;  // increase weight even more
			final int                        debugLevel,
			final int                        dbg_X,
			final int                        dbg_Y)
	{
		final int tilesX =        tileProcessor.getTilesX();
		final int superTileSize = tileProcessor.getSuperTileSize();
		final int stilesX = (tilesX + superTileSize -1)/superTileSize;
		final int debug_stile = dbg_Y * stilesX + dbg_X;
		final Thread[] threads = ImageDtt.newThreadArray(tileProcessor.threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		final AtomicInteger ai_numThread = new AtomicInteger(0);
		final int [] replaced = new int[threads.length];
		final int [][] dbg_dirsYX = {{-1, 0},{-1,1},{0,1},{1,1},{1,0},{1,-1},{0,-1},{-1,-1}};

		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				@Override
				public void run() {
					double [][] dbg_img = null;
					String [] dbg_titles = null;
					TilePlanes.PlaneData[][]   dbg_planes = planes;
					int numThread = ai_numThread.getAndIncrement(); // unique number of thread to write to rslt_diffs[numThread]
					for (int nsTile = ai.getAndIncrement(); nsTile < planes.length; nsTile = ai.getAndIncrement()) {
						int sty = nsTile / stilesX;
						int stx = nsTile % stilesX;
						int dl = ((debugLevel > -1) && (nsTile == debug_stile)) ? 4:0;
						int num_bplanes = 2; // brokenPd[nsTile].length;
						int num_meas_layers = measuredLayers.getNumLayers();
						int [] dbg_indices = {
								0,
								num_bplanes,
								num_bplanes + 3 * num_meas_layers,
								num_bplanes + 3 * num_meas_layers + 2 * (num_bplanes * num_meas_layers),
								num_bplanes + 3 * num_meas_layers + 2 * (num_bplanes * num_meas_layers) + num_bplanes,
						};

						if ( brokenPd[nsTile] != null) {
							if (dl > 0){
								System.out.println("replaceBrokenPlanes nsTile="+nsTile);
								dbg_titles = new String[dbg_indices[4]];
								for (int i = 0; i < num_bplanes; i++){
									dbg_titles[i + dbg_indices[0]]="src_plane_"+i;
								}
								for (int j = 0; j < num_meas_layers; j++){
									dbg_titles[3* j + dbg_indices[1] + 0] =     "disp_"+j;
									dbg_titles[3* j + dbg_indices[1] + 1] =     "mask_disp_"+j;
									dbg_titles[3* j + dbg_indices[1] + 2] =     "str_"+j;
								}
								for (int i = 0; i < num_bplanes; i++){
									for (int j = 0; j < num_meas_layers; j++){
										dbg_titles[2 * (num_meas_layers * i +j) + dbg_indices[2] + 0] =  "t"+i+"_m"+j;
										dbg_titles[2 * (num_meas_layers * i +j) + dbg_indices[2] + 1] =  "ft"+i+"_m"+j;
									}
								}
								for (int i = 0; i < num_bplanes; i++){
									dbg_titles[i + dbg_indices[3]]="plane_"+i;
								}
							}

							TilePlanes.PlaneData[][] bpd = new TilePlanes.PlaneData[brokenPd[nsTile].length][];
							for (int np = 0; np < bpd.length; np++){
								if (brokenPd[nsTile][np] != null){
									bpd[np] = new TilePlanes.PlaneData[2];
									bpd[np][0] = brokenPd[nsTile][np][0].clone();
									bpd[np][1] = brokenPd[nsTile][np][1].clone(); // no need to clone [2]
								}
							}
							if (splitXY) {
								for (int np = 0; np < bpd.length; np++){
									if (brokenPd[nsTile][np] != null){
										// if it will fail (return false) no selection masks, still possible to separate
										planes[nsTile][np].calcSelMasks(
												bpd[np][0],
												bpd[np][1],
												splitXYTolerance, // double    tolerance,
												min_tiles,
												dl); // int       debugLevel
									}
								}
							}

							// now works on the pairs in bpd, null them out in case of failure
							for (int np = 0; np < bpd.length; np++){
								if ((dl > 1) && (bpd[np] !=null)) {
									int ss2 = 2 * superTileSize;
									dbg_img = new double [dbg_titles.length][];
									for (int ni = 0; ni < num_bplanes ; ni++) {
										dbg_img[ni + dbg_indices[0]] =  bpd[np][ni].getDoublePlaneDisparity(false);
										boolean [] sm = bpd[np][ni].getSelMask();
										// checkerboard pattern over disabled part of the plane
										if (sm != null) {
											for (int i = 0; i < sm.length; i++){
												if (!sm[i] && ((((i % ss2)+ (i / ss2)) & 1) != 0)){
													dbg_img[ni + dbg_indices[0]][i] = Double.NaN;
												}
											}
										}
									}
									// add connections to neighbors
									for (int dir = 0; dir <8; dir++){
										int indx = (ss2 / 2) * (ss2 + 1);
										int dindx = dbg_dirsYX[dir][0] * ss2 + dbg_dirsYX[dir][1];
										for (int l = 0; l < (ss2 / 2 - 2); l++){
											if ( l > 2){ // keep center not modified
												for (int ni = 0; ni < num_bplanes ; ni++) {
													if (bpd[np][ni].getNeibBest(dir) >= 0){
														dbg_img[ni + dbg_indices[0]][indx] = Double.NaN;
													}
												}
											}
											indx += dindx;
										}
									}
									// show disparity (all and masked for this plane) and strength
									for (int ml = 0; ml < num_meas_layers; ml++) {
										double [][] dbg_disp_str = planes[nsTile][np].getMeasuredLayers().getDisparityStrengthML(
												ml,                     // int num_layer,
												stx,        // int stX,
												sty,        // int stY,
												null,       // boolean [] sel_in,
												mlfp,
//												strength_floor,         //  double strength_floor,
//												strength_pow,  // double strength_pow,
												true);                  // boolean null_if_none);
										if (dbg_disp_str != null){
											boolean [] dbg_msel = planes[nsTile][np].getMeasSelection(ml);
											dbg_img[3 * ml + dbg_indices[1] + 0] = dbg_disp_str[0];
											dbg_img[3 * ml + dbg_indices[1] + 2] = dbg_disp_str[1];
											dbg_img[3 * ml + dbg_indices[1] + 1] = dbg_disp_str[0].clone();
											if (dbg_msel != null) {
												for (int i = 0; i < dbg_msel.length; i++){
													if (!dbg_msel[i]) dbg_img[3 * ml + dbg_indices[1] + 1][i] = Double.NaN;
												}
											}
										}
									}

								} // (dl > 1)

								if (brokenPd[nsTile][np] != null){
									// Is it a single set to replace current plane? If yes, use all tiles, not just original selection
									int np_min = LOWEST_PLANE(planes[nsTile].length);
									boolean single_plane = (np_min == (planes[nsTile].length - 1));

									Boolean OK = planes[nsTile][np].splitPlaneTiles ( // now always OK TODO: add layer mask or make use_other_planes int
											bpd[np],          // PlaneData [] pd_set,
											single_plane,     // boolean      single_plane,
											max_diff,         // double       max_diff, // maximal disparity difference (0 - any), will be normalized by dispNorm
											other_diff,
											non_exclusive,    // boolean      non_exclusive,
											use_other_planes, // boolean      use_other_planes,
											smplMode,
											mlfp,
											measSel,          // int          measSel, // Select measurements for supertiles : +1 - combo, +2 - quad +4 - hor +8 - vert
											allow_parallel,   //boolean      allow_parallel,
											dl); // int          debugLevel)
									if (!OK) {
										brokenPd[nsTile][np] = null;
										continue;
									}
									if (dl > 1) {// show initial selection (made of disparity clone - dbg_img[3 * ml + dbg_indices[1] + 0] = dbg_disp_str[0];
										for (int ni = 0; ni < num_bplanes ; ni++) {
											boolean [][] msel = new boolean [num_meas_layers][];
											for (int ml = 0; ml < num_meas_layers; ml++) {
												msel[ml] = bpd[np][ni].getMeasSelection(ml);
												if (msel[ml] != null) {
													dbg_img[2 * (num_meas_layers * ni +ml) + dbg_indices[2] + 0] =  dbg_img[3 * ml + dbg_indices[1] + 0].clone();
													for (int i = 0; i < msel[ml].length; i++){
														if (!msel[ml][i]) dbg_img[2 * (num_meas_layers * ni +ml) + dbg_indices[2] + 0][i] = Double.NaN;
													}
												}
											}
										}
									}

									for (int npip = 0; npip < bpd[np].length; npip++) {
										double [][][] disp_strength = bpd[np][npip].getPlaneFromMeas(
												null,           // boolean [][] tile_sel, // null - do not use, {} use all (will be modified)
												null,
												disp_far,       // double       disp_far, // minimal disparity to select (or NaN)
												disp_near,      // double       disp_near, // maximal disparity to select (or NaN)
												dispNorm,       // double       dispNorm,   //  Normalize disparities to the average if above
												min_weight,     // double       min_weight,
												min_tiles,      // int          min_tiles,
// OK?
												smplMode,
												mlfp,
												dl);            // int          debugLevel)
										if (disp_strength == null) break;
										// remove outliers //removeOutliers
										// now try to remove outliers
										double targetV = corrMaxEigen(
												targetEigen,
												dispNorm,
												bpd[np][npip]);
//										if (bpd[np][npip].getValue() > targetV) {
										boolean almost_fronto = bpd[np][npip].isFronto(
												fronto_tol,    // double fronto_tol,
												disp_strength,      // double [][][] disp_str,
												debugLevel);   // int debugLevel);
										// seems it is OK not to normalize fronto_rms, as it is only needed for far surfaces
										double targetEigen = almost_fronto ? (fronto_rms*fronto_rms): targetV;
										// FIXME: temporary - make a configurable parameter
										double fronto_fract_outliers = 0.8; // 2 * inFractOutliers;
										double fractOutliers = almost_fronto ? (fronto_fract_outliers) : inFractOutliers;
										int max_outliers = (int) Math.round(bpd[np][npip].getNumPoints() * fractOutliers);
										if (max_outliers > maxOutliers) max_outliers = maxOutliers;


//										if (bpd[np][npip].getNormValue() > targetV) {
										// in fronto mode - run at least once
										if (almost_fronto) {
											bpd[np][npip].growFrontoSelection(
													4,              // int           min_neib,
													0.1,            // double        max_over_avg,
													0.00,           // double        max_under_avg,
													disp_strength); // double [][][] disp_str)
										}

										if ((bpd[np][npip].getNormValue() > targetEigen) || almost_fronto) { // targetV) {
											OK = bpd[np][npip].removeOutliers( // getPlaneFromMeas should already have run
													fronto_tol, // fronto tolerance (pix) - treat almost fronto as fronto (constant disparity). <= 0 - disable this feature
													fronto_rms,    // double        fronto_rms,  // Target rms for the fronto planes - same as sqrt(plMaxEigen) for other planes. May be tighter
													fronto_offs,   // double        fronto_offs,   //        =   0.2;  // increasing weight of the near tiles by using difference between the reduced
													fronto_pow,    // double        fronto_pow,    //        =   1.0;  // increase weight even more
													disp_strength,
													targetV,      // double     targetEigen, // target eigenvalue for primary axis (is disparity-dependent, so is non-constant)
													max_outliers, // int        maxRemoved,  // maximal number of tiles to remove (not a constant)
													dl); // int        debugLevel)
											if (!OK) break;
										}

										bpd[np][npip].getWorldPlaneFromMeas(
												null,           // boolean [][] tile_sel, // null - do not use, {} use all (will be modified)
												disp_strength,
												disp_far,       // double       disp_far, // minimal disparity to select (or NaN)
												disp_near,      // double       disp_near, // maximal disparity to select (or NaN)
												dispNorm,       // double       dispNorm,   //  Normalize disparities to the average if above
												min_weight,     // double       min_weight,
												min_tiles,      // int          min_tiles,
//												strength_floor, // double       strength_floor,
//												strength_pow,   // double       strength_pow,
// OK?
												smplMode,
												mlfp,
//												smplSide,
//												smplNum,
//												smplRms,
//												smplWnd,           // was not here: final boolean    smplWnd,  // use window functions for the samples

//												max_abs_tilt,  // 2.0; // Maximal absolute tilt in pixels/tile
//												max_rel_tilt,  // 0.2; // Maximal relative tilt in pixels/tile/disparity
//												damp_tilt,     //    0.001; // Damp tilt to handle insufficient  (co-linear)data
//												min_tilt_disp, // 4.0; // Disparity switch between filtering modes - near objects use tilts, far - use max disparity
//												transition,    // 1.0; // Mode transition range (between tilted and maximal disparity)
//												far_mode,      //     1;   // Far objects filtering mode (0 - off, 1 - power of disparity)
//												far_power,     //    1.0; // Raise disparity to this power before averaging for far objects

												dl);            // int          debugLevel)
									}
									if (!OK) {
										brokenPd[nsTile][np] = null;
										continue;
									}

									if (dl > 1) {// filtered selection (made of disparity clone - dbg_img[3 * ml + dbg_indices[1] + 0] = dbg_disp_str[0];
										for (int ni = 0; ni < num_bplanes ; ni++) {
											boolean [][] msel = new boolean [num_meas_layers][];
											for (int ml = 0; ml < num_meas_layers; ml++) {
												msel[ml] = bpd[np][ni].getMeasSelection(ml);
												if (msel[ml] != null) {
													dbg_img[2 * (num_meas_layers * ni +ml) + dbg_indices[2] +1] =  dbg_img[3 * ml + dbg_indices[1] + 0].clone();
													for (int i = 0; i < msel[ml].length; i++){
														if (!msel[ml][i]) dbg_img[2 * (num_meas_layers * ni +ml) + dbg_indices[2] + 1][i] = Double.NaN;
													}
												}
											}
										}
										// show new planes
										for (int ni = 0; ni < num_bplanes ; ni++) {
											dbg_img[ni + dbg_indices[3]] =  bpd[np][ni].getDoublePlaneDisparity(false);
										}
										// actually show the image

										if (debugLevel > -1){
											ShowDoubleFloatArrays sdfa_instance = new ShowDoubleFloatArrays();
											sdfa_instance.showArrays(dbg_img, 2 * superTileSize, 2 * superTileSize, true, "replaceBrokenPlanes"+stx+"_y"+sty+"_p"+np, dbg_titles);
										}
									}
								}
							}
							// now see if anything is left
							int npairs = 0;
							for (int np = 0; np < bpd.length; np++){
								if (bpd[np] != null) npairs++;
							}

							// Fixing that first plane in multi=-plane is not processed
							int np_start = LOWEST_PLANE(2); // planes[nsTile].length);

							if (npairs > 0){
								TilePlanes.PlaneData[] old_planes = planes[nsTile];
								int old_len = old_planes.length - LOWEST_PLANE(old_planes.length);
								planes[nsTile] = new TilePlanes.PlaneData[old_len + npairs + np_start];
								int npr = 0;
								for (int np = LOWEST_PLANE(old_planes.length); np < np_start;  np++) {
									planes[nsTile][npr++] = null; // see if there will be any conflicts, fix or replace with old_planes[0];
								}
								for (int np = 0; np < bpd.length; np++){
									if (bpd[np] != null) {
										planes[nsTile][npr++] = bpd[np][0];
										planes[nsTile][npr++] = bpd[np][1];
									} else {
										planes[nsTile][npr++] = old_planes[np];
									}
								}
							}
							replaced[numThread] += npairs;
						}
					}
				}
			};
		}
		ImageDtt.startAndJoin(threads);
		int replaced_planes = 0;
		for (int i = 0; i < replaced.length; i++) replaced_planes += replaced[i];
		return replaced_planes;
	}

	/**
	 * Create candidate planes to break a single plane in 2 by splitting consecutive connected
	 * neighbors in 2 groups that make the smallest weighted sum of the eigenvalues
	 * Groups may be optionally supplemented by the center supertile
	 * Result weights are sums of participating supertile weights, for normalization  (if needed)
	 * divide by number of neighbors plus center_pull
	 * @param center_planes [per supertile][per plane] array of plane objects to use as
	 *        the source of connections and optionally other data for the center tiles
	 * @param neib_planes  [per supertile][per plane] array of plane objects to use for
	 *        neighbor data. May be the same as center_planes or different
	 * @param center_pull - merge with center plane when calculating half-planes with this
	 *  relative weight: 0.0 - only neighbors, 1.0 - same weight of the center as each neighbor.
	 * @param min_neibs - minimal number of connected neighbors to work with (>=2)
	 * @param splitMinWeight -  minimal weight of the tile to split
	 * @param splitMinQuality - minimal quality of split (EV of the single plane over weighted average
	 *        of each half-plane
	 * @param preferDisparity - the first eigenvalue/vector is the most disparity-like
	 *                          (false - smallest eigenvalue)
	 * @param debugLevel debug level
	 * @param dbg_X supertile horizontal index to show debug information
	 * @param dbg_Y supertile vertical index to show debug information
	 * @return a pair of plane objects for each [supertile][plane][3] and a combination of both
	 *         (to compare eigenvalues)
	 */
	public TilePlanes.PlaneData[][][] breakPlanesToPairs(
			final TilePlanes.PlaneData[][] center_planes,   // measured_planes,
			final TilePlanes.PlaneData[][] neib_planes,     //mod_planes,
			final double                   center_pull,
			final int                      min_neibs,       // 2
			final double                   splitMinWeight,  //     =  2.0;  // Minimal weight of split plains to show
			final double                   splitMinQuality, //    =  1.1;  // Minimal split quality to show
			final boolean                  preferDisparity, // Always start with disparity-most axis (false - lowest eigenvalue)
			final int                      debugLevel,
			final int                      dbg_X,
			final int                      dbg_Y)
	{
		final int [][] dirsYX = {{-1, 0},{-1,1},{0,1},{1,1},{1,0},{1,-1},{0,-1},{-1,-1}};
		final int tilesX =        tileProcessor.getTilesX();
//		final int tilesY =        tileProcessor.getTilesY();
		final int superTileSize = tileProcessor.getSuperTileSize();
		final int stilesX = (tilesX + superTileSize -1)/superTileSize;
//		final int stilesY = (tilesY + superTileSize -1)/superTileSize;
		final int debug_stile = dbg_Y * stilesX + dbg_X;
		final TilePlanes.PlaneData[][][] rslt_planes = new TilePlanes.PlaneData[center_planes.length][][];
		final Thread[] threads = ImageDtt.newThreadArray(tileProcessor.threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		final String [] titles= {"center", "first", "second", "all", "f-s"};

		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				@Override
				public void run() {
					double [][] dbg_img=null;
					for (int nsTile0 = ai.getAndIncrement(); nsTile0 < center_planes.length; nsTile0 = ai.getAndIncrement()) {
						int sty0 = nsTile0 / stilesX;
						int stx0 = nsTile0 % stilesX;
						int dl = ((debugLevel > -1) && (nsTile0 == debug_stile)) ? 1:0;
						if ( center_planes[nsTile0] != null) {
							rslt_planes[nsTile0] = new TilePlanes.PlaneData[center_planes[nsTile0].length][];
							if (dl > 0){
								System.out.println("breakPlanesToPairs nsTile0="+nsTile0);
								dbg_img = new double [titles.length][];
							}
							int np0_min = (center_planes[nsTile0].length > 1) ? 1:0; // Modify if overall plane will be removed
							for (int np0 = np0_min; np0 < center_planes[nsTile0].length; np0 ++){
								TilePlanes.PlaneData center_plane = center_planes[nsTile0][np0];
								if (center_plane.getWeight() >= splitMinWeight) {
									if (dl > 0) dbg_img[ 0] = center_plane. getTriplePlaneDisparity();
									int [] neibs = center_plane.getNeibBest();
									int num_neibs = 0;
									for (int i = 0; i < neibs.length; i++) if (neibs[i] >= 0) num_neibs++;
									if (num_neibs >= min_neibs) {
										// create all pairs to test
										int num_splits = num_neibs * (num_neibs - 1) / 2;
										double [] split_quality = new double [num_splits]; // weighted sum of eigenvalues of merged
										int [][][] neibs12 = new int [2][num_splits][];
										TilePlanes.PlaneData[][] plane_triads = new TilePlanes.PlaneData[num_splits][3];
										int [] neib_index = new int [num_neibs];
										int nn = 0;
										for (int i = 0; i < neibs.length; i++) if (neibs[i] >= 0) neib_index[nn++] = i;
										nn = 0;
										for (int nStart = 1; nStart < num_neibs; nStart++) {
											for (int nEnd = nStart+1; nEnd <= num_neibs; nEnd++){
												neibs12[0][nn] = neibs.clone();
												neibs12[1][nn] = neibs.clone();
												for (int i = 0; i < neib_index.length; i++) {
													if ((i < nStart) || (i >= nEnd)){
														neibs12[0][nn][neib_index[i]] = -1;
													} else {
														neibs12[1][nn][neib_index[i]] = -1;
													}
												}
												nn++;
											}
										}

										for (int nSplit = 0; nSplit < num_splits; nSplit++) {
											// Calculate merged plane for each selection (and other), weighted average and store to split_quality
											for (int nis = 0; nis < 2; nis++) {
												plane_triads[nSplit][nis] = center_planes[nsTile0][np0].clone();
												plane_triads[nSplit][nis].scaleWeight(center_pull);
												for (int dir = 0; dir <8; dir++){
													int np = neibs12[nis][nSplit][dir];
													if (np >= 0){
														int stx = stx0 + dirsYX[dir][1];
														int sty = sty0 + dirsYX[dir][0];
														int nsTile = sty * stilesX + stx; // from where to get
														TilePlanes.PlaneData other_plane = plane_triads[nSplit][nis].getPlaneToThis(
																neib_planes[nsTile][np],
																dl-1); // debugLevel);
														if (other_plane != null){
															if (plane_triads[nSplit][nis].getWeight() > 0.0){
																plane_triads[nSplit][nis] = plane_triads[nSplit][nis].mergePlaneToThis(
																		other_plane, // PlaneData otherPd,
																		1.0,         // double    scale_other,
																		1.0,         // double     starWeightPwr,    // Use this power of tile weight when calculating connection cost
																		false,       // boolean   ignore_weights,
																		true, // boolean   sum_weights,
																		preferDisparity,
																		dl-1); // int       debugLevel)
															} else {
																plane_triads[nSplit][nis] = other_plane;
															}
														} else {
															plane_triads[nSplit][nis] = null;
															break;
														}
													}
												}
											} // for (int nis = 0; nis < 2; nis++) {
											if ((plane_triads[nSplit][0] != null) && (plane_triads[nSplit][1] != null)){
												double w1 = plane_triads[nSplit][0].getWeight();
												double w2 = plane_triads[nSplit][1].getWeight();
//												split_quality[nSplit] = (w1 * plane_triads[nSplit][0].getValue() + w2 * plane_triads[nSplit][1].getValue())/
												split_quality[nSplit] = (w1 * plane_triads[nSplit][0].getNormValue() + w2 * plane_triads[nSplit][1].getNormValue())/
														(w1 + w2);
												if (dl >0){
													plane_triads[nSplit][2] = plane_triads[nSplit][0].clone().mergePlaneToThis(
															plane_triads[nSplit][1], // PlaneData otherPd,
															1.0,                    // double    scale_other,
															1.0,                    // double     starWeightPwr,    // Use this power of tile weight when calculating connection cost
															false,                  // boolean   ignore_weights,
															true, // boolean   sum_weights,
															preferDisparity,
															dl-1);
												}
											} else {
												split_quality[nSplit] = Double.NaN;
											}
											if (dl >0){
												System.out.println("breakPlanesToPairs(): split_quality["+nSplit+"] = "+split_quality[nSplit]);
											}
										}
										// find minimum in split_quality and generate a pair of plane objects, setting neighbors for each
										// later use these plane pairs to assign tiles to each and generate a new eigenvalues/vectors
										// TODO: How to handle a pair? Any special treatment (like fusing?),
										// if the plane intersection line is inside supertile - use min/max for snapping
										int best_index = -1;
										for (int nSplit = 0; nSplit < num_splits; nSplit++) {
											if (!Double.isNaN(split_quality[nSplit]) && ((best_index < 0) || (split_quality[nSplit] < split_quality[best_index]))){
												best_index = nSplit;
											}
										}
										if (best_index >= 0) {
											plane_triads[best_index][2] = plane_triads[best_index][0].clone().mergePlaneToThis(
													plane_triads[best_index][1], // PlaneData otherPd,
													1.0,                     // double    scale_other,
													1.0, // double     starWeightPwr,    // Use this power of tile weight when calculating connection cost
													false,                   // boolean   ignore_weights,
													true,                    // boolean   sum_weights,
													preferDisparity,
													dl-1);
											//  check quality
											if ((plane_triads[best_index][2] != null) &&
//													(plane_triads[best_index][2].getValue()/split_quality[best_index] > splitMinQuality)) {
													(plane_triads[best_index][2].getNormValue()/split_quality[best_index] > splitMinQuality)) {
												for (int nis = 0; nis < 2; nis++) {
													plane_triads[best_index][nis].setNeibBest(neibs12[nis][best_index]);
												}
												rslt_planes[nsTile0][np0] = plane_triads[best_index];
												if (dl > 0) {
													int ss3 = 3 * superTileSize;
													dbg_img[ 1] = plane_triads[best_index][0].getTriplePlaneDisparity();
													dbg_img[ 2] = plane_triads[best_index][1].getTriplePlaneDisparity();
													dbg_img[ 3] = plane_triads[best_index][2].getTriplePlaneDisparity();
													dbg_img[ 4] = new double [dbg_img[ 1].length];
													for (int i = 0; i< dbg_img[ 1].length; i++){
														dbg_img[ 4][i] = dbg_img[ 1][i] - dbg_img[ 2][i];
													}
													int [] vby = new int [ss3];
													double [] dbg_split = dbg_img[ 4];
													for (int dx = 0; dx < ss3; dx ++){
														vby[dx] = 0;
														for (int dy = 1; dy < ss3; dy ++){
															int indx = dy * ss3 + dx;
															if (Math.abs(dbg_split[indx]) < Math.abs(dbg_split[vby[dx] * ss3 + dx])){
																vby[dx] = dy; // indx;
															}
														}
													}
													for (int dy = 0; dy < ss3; dy ++){
														int bx = 0;
														for (int dx = 1; dx < ss3; dx ++){
															if (Math.abs(dbg_split[dy * ss3 + dx]) < Math.abs(dbg_split[dy * ss3 + bx])){
																bx = dx;
															}
														}
														if ((bx > 0) && (bx < (ss3-1))) {
															dbg_split[dy * ss3 + bx] = Double.NaN;
															dbg_img[ 3][dy * ss3 + bx] = Double.NaN;
														}
													}

													for (int dx = 0; dx < ss3; dx ++){
														int by = vby[dx];
														if ((by > 0) && (by < (ss3-1))) {
															dbg_split[by * ss3 + dx] = Double.NaN;
															dbg_img[ 3][by * ss3 + dx] = Double.NaN;
														}
													}

													for (int dir = 0; dir <8; dir++){
														int indx = (ss3 / 2) * (ss3 + 1);
														int dindx = dirsYX[dir][0] * ss3 + dirsYX[dir][1];

														for (int l = 0; l < (ss3 / 2 - 2); l++){
															if ( l > 2){ // keep center not modified
																if (plane_triads[best_index][0].getNeibBest(dir) >= 0){
																	dbg_img[ 1][indx] = Double.NaN;
																}
																if (plane_triads[best_index][1].getNeibBest(dir) >= 0){
																	dbg_img[ 2][indx] = Double.NaN;
																}
															}
															indx += dindx;
														}
													}
													System.out.println("breakPlanesToPairs(): num_splits="+num_splits+", best_index="+best_index);
													for (int nSplit = 0; nSplit < num_splits; nSplit++) {
														if ((plane_triads[nSplit][0] != null) && (plane_triads[nSplit][1] != null) && (plane_triads[nSplit][2] != null)){
															System.out.print("breakPlanesToPairs():"+nSplit+" ");
															for (int dir = 0; dir < 8; dir++){
																System.out.print((neibs12[0][nSplit][dir] >=0)?"+":"-");
															}
															System.out.print(" ");
															for (int dir = 0; dir < 8; dir++){
																System.out.print((neibs12[1][nSplit][dir] >=0)?"+":"-");
															}
															System.out.print(" ");
															for (int dir = 0; dir < 8; dir++){
																System.out.print((plane_triads[nSplit][2].getNeibBest(dir) >=0)?"+":"-");
															}
															System.out.println(" : "+split_quality[nSplit]+
//																	" "+plane_triads[nSplit][0].getValue()+
//																	" "+plane_triads[nSplit][1].getValue()+
//																	" "+plane_triads[nSplit][2].getValue()+
																	" "+plane_triads[nSplit][0].getNormValue()+
																	" "+plane_triads[nSplit][1].getNormValue()+
																	" "+plane_triads[nSplit][2].getNormValue()+
																	" w"+plane_triads[nSplit][0].getWeight()+
																	" w"+plane_triads[nSplit][1].getWeight()+
																	" w"+plane_triads[nSplit][2].getWeight()
																	);
														}
													}

													// merged_eig_val are not set - will have to be recalculated when updating after changing tile selections
													// to make each plane
													// TODO save gain from splitting to planes

													if (debugLevel > 3){
														ShowDoubleFloatArrays sdfa_instance = new ShowDoubleFloatArrays();
														sdfa_instance.showArrays(dbg_img, 3*superTileSize, 3*superTileSize, true, "breakPlanesToPairs"+stx0+"_y"+sty0+"_p"+np0, titles);
													}
												}
											}
										}
									}
								}
							}
						}
					}
				}
			};
		}
		ImageDtt.startAndJoin(threads);
		return rslt_planes;
	}


	public boolean [][] selectPlanes(
			final double dispNorm,
			final double maxEigen, // maximal eigenvalue of planes to consider
			final double minWeight, // minimal pain weight to consider
			final boolean weak_connected, // select connected planes even if they are weak/thick
			TilePlanes.PlaneData [][] planes)
	{
		final int tilesX =        tileProcessor.getTilesX();
		final int tilesY =        tileProcessor.getTilesY();
		final int superTileSize = tileProcessor.getSuperTileSize();
		//				final int tileSize =      tileProcessor.getTileSize();
		final int stilesX = (tilesX + superTileSize -1)/superTileSize;
		final int stilesY = (tilesY + superTileSize -1)/superTileSize;
		boolean [][] selection = new boolean [planes.length][];
		for (int sty = 0; sty < stilesY; sty++){
			for (int stx = 0; stx < stilesX; stx++){
				int nsTile = sty * stilesX + stx;
				if (planes[nsTile] != null) {
					selection[nsTile] = new boolean [planes[nsTile].length];
					for (int np = 0; np < planes[nsTile].length; np++){
//						if ((planes[nsTile][np] != null) && (planes[nsTile][np].getWeight() >= minWeight)){
						if (planes[nsTile][np] != null){
							boolean has_conn = false;
							for (int dir = 0; dir < 8; dir++){
								if (planes[nsTile][np].getNeibBest(dir) >= 0){
									has_conn = true;
									break;
								}
							}
							if (weak_connected && has_conn) {
								selection[nsTile][np] = true;
							} else if  (planes[nsTile][np].getWeight() >= minWeight) {
								//&& (planes[nsTile][np].getValue() < maxEigen)){
//								double eigVal = planes[nsTile][np].getValue();
								double eigVal = planes[nsTile][np].getNormValue();
								double disp = planes[nsTile][np].getZxy()[0];
								selection[nsTile][np] = (eigVal < corrMaxEigen(
										maxEigen,
										dispNorm,
										disp));
							}
						}
					}
				}
			}
		}
		return selection;
	}

	public int [][] createSupertileShells(
			boolean[][]                     selection, // may be null
			boolean                         use_all, // use plane 0 even if there are more than 1
			boolean                         keep_orphans, // single-cell shells
			double                          orphan_strength, // minimal strength no keep orphans (single-tile planes)
			final TilePlanes.PlaneData [][] planes,
			final int         debugLevel,
			final int         dbg_X,
			final int         dbg_Y)
	{
		class TPW {
			int t;
			int p;
			double w;
			TPW (int t, int p){
				this.t= t;
				this.p= p;
				this.w = planes[t][p].getWeight();
			}
		}
		final int tilesX =        tileProcessor.getTilesX();
		final int superTileSize = tileProcessor.getSuperTileSize();
		final int stilesX = (tilesX + superTileSize -1)/superTileSize;
		final int [] dirs= {-stilesX, -stilesX+1, 1, stilesX+1, stilesX, stilesX-1, -1, -stilesX-1};
		boolean [][] field = new boolean [planes.length][];
		for (int nsTile = 0; nsTile < planes.length; nsTile++){
			if (nsTile==724){
				System.out.println("createSupertileShells()1 nsTile="+nsTile);
			}
			if (planes[nsTile] != null){
				field[nsTile] = new boolean [planes[nsTile].length];
				if (selection != null) {
					for (int i = 0; i < selection[nsTile].length; i++){
						if (!selection[nsTile][i]) field[nsTile][i] = true; // used/prohibited
					}
				}
				if (!use_all && (planes[nsTile].length > 1)){
					field[nsTile][0] = true; // used/prohibited
				}
			}
		}
		ArrayList<ArrayList<TPW>> listShells = new ArrayList<ArrayList<TPW>>();
		for (int nsTile = 0; nsTile < field.length; nsTile++){
			if (nsTile==724){
				System.out.println("createSupertileShells()2 nsTile="+nsTile);
			}

			if (field[nsTile] != null){
				for (int np = 0; np < field[nsTile].length; np++){
					if (!field[nsTile][np]){
						// see if
						ArrayList<TPW> listWave = new ArrayList<TPW>();
						int read_p = 0;
						//								int num_cluster = listShells.size()+1;
						TPW sp = new TPW(nsTile,np);
						field[sp.t][sp.p] = true; // just 1 would be OK
						listWave.add(sp);
						while (read_p < listWave.size()) { // keeping all list, not just front
							TPW tail = listWave.get(read_p++);
							int [] neibs = planes[tail.t][tail.p].getNeibBest();
							for (int dir = 0; dir < neibs.length; dir++){
								if (neibs[dir] >= 0) {
									sp = new TPW(
											tail.t + dirs[dir],
											neibs[dir]);
									if (!field[sp.t][sp.p] && ((selection == null) || (selection[sp.t][sp.p]))){
										field[sp.t][sp.p] = true;
										listWave.add(sp);
									}
								}
							}

						}
						if (listWave.size() > 1) {
							listShells.add(listWave);
						} else if (keep_orphans && (listWave.size() == 1)){ // so far always >= 1
							TPW tail = listWave.get(0);
							if (planes[tail.t][tail.p].getWeight() >= orphan_strength){
								listShells.add(listWave);
							}
						}
						if (debugLevel > 1){
							System.out.println("CreateSupertileShells(): "+(listShells.size())+":"+listWave.size());
						}
					}
				}
			}
		}
		// sort list
		Collections.sort(listShells, new Comparator<ArrayList<TPW>>() {
			@Override
			public int compare(ArrayList<TPW> lhs, ArrayList<TPW> rhs) {
				// -1 - less than, 1 - greater than, 0 - equal, all inverted for descending
				int sgn = lhs.size() > rhs.size() ? -1 : (lhs.size() < rhs.size() ) ? 1 : 0;
				if (sgn != 0) return sgn;
				double wl = 0.0, wr = 0.0;
				for (TPW p:lhs) wl += p.w;
				for (TPW p:rhs) wr += p.w;
				return wl > wr ? -1 : (wl < wr ) ? 1 : 0;
			}
		});


		int [][] shells = new int [planes.length][];
		for (int nsTile = 0; nsTile < planes.length; nsTile++){
			if (planes[nsTile] != null){
				shells[nsTile] = new int [planes[nsTile].length];
			}
		}
		for (int num_shell = 0; num_shell < listShells.size(); num_shell++){
			for (TPW p:listShells.get(num_shell)){
				shells[p.t][p.p] = num_shell+1;
			}
			if (debugLevel > 0){
				double w = 0;
				for (TPW p:listShells.get(num_shell)){
					w += p.w;
				}
				System.out.println("CreateSupertileShells(): "+(num_shell + 1)+":"+listShells.get(num_shell).size()+", total strength = "+w);
			}

		}
		return shells;
	}

	public double normDispADiff(
			double disp1,
			double disp2,
			double disp_norm
			)
	{
		double disp_avg = 0.5* (disp1 + disp2);
		double adiff = Math.abs(disp1 - disp2);
		return (disp_avg <= disp_norm) ? adiff : (adiff * disp_norm / disp_avg);
	}

	public double [] snapDisparity(
			final double []  disparity,
			final double []  strength,
			final boolean [] selection,       // can be null
			final double     dispNorm, // plDispNorm
			final double     snapDispAny,    //    =  .2,   // Maximal (scaled by plDispNorm) disparity difference to snap to plane at any strength
			final double     snapStrengthAny,    //    =  .2,   //  Maximal strength to fit any distance (if does not fit otherwise - treat as zero strength
			final double     snapNegAny,     //    =  .2,   // Maximal negative disparity difference from the best match
			final double     snapDispMax,    //    =  .5,   // Maximal (scaled by plDispNorm) disparity difference to snap to plane at low strength
			final double     snapDispWeight, //    =  .5,   // Maximal disparity diff. by weight product to snap to plane
			final int        snapZeroMode,   // Zero strength snap mode: 0: no special treatment, 1 - strongest, 2 - farthest
			final int        debugLevel)
	{
		final int tilesX =        tileProcessor.getTilesX();
//		final int tilesY =        tileProcessor.getTilesY();
		final int superTileSize = tileProcessor.getSuperTileSize();
		final int stilesX = (tilesX + superTileSize -1)/superTileSize;
//		final int stilesY = (tilesY + superTileSize -1)/superTileSize;
		final double [][] surfaces = getSurfaces();
		final int [][] shell_map =   getShellMap();
		final TilePlanes.PlaneData [][] planes = getPlanesMod(); // to get strength;
		final double [] snap_disp =  disparity.clone();

		final Thread[] threads = ImageDtt.newThreadArray(tileProcessor.threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				@Override
				public void run() {
					for (int nTile = ai.getAndIncrement(); nTile < disparity.length; nTile = ai.getAndIncrement()) {
						if ((selection == null) || selection[nTile]) {
							int tX = nTile % tilesX;
							int tY = nTile / tilesX;
							int stX = tX / superTileSize;
							int stY = tY / superTileSize;
							int nsTile = stY * stilesX + stX;
							int nTile_st = (tY * stilesX * superTileSize) + tX; // surfaces maybe (is) larger than tilesX *tilesY if
							// tilesX, tilesY are not multiples of superTileSize !
							if (nsTile == 724){
//								System.out.println("snapDisparity(): tX:tY = "+tX+":"+tY+", stX:stY="+stX+":"+stY+
//										", nTile="+nTile+", nsTile="+nsTile+", nTile_st = "+nTile_st);
							}

							// is there any surface to snap to?
							boolean no_surf = true;
							if (shell_map[nsTile] != null){
								for (int np = 0; no_surf && (np < shell_map[nsTile].length); np++){
									if (shell_map[nsTile][np] > 0) no_surf = false;
								}
							}
							if (!no_surf) {
								double [] surf_disparity = new double[shell_map[nsTile].length];
								for (int np = 0; np < surf_disparity.length; np++){
									if (shell_map[nsTile][np] > 0) {
										surf_disparity[np] = surfaces[shell_map[nsTile][np] -1][nTile_st];
									} else {
										surf_disparity[np] = Double.NaN;
									}
								}
								int best_indx = -1;
								double best_disp_diff = Double.NaN;
								double this_disp = disparity[nTile];
								double this_strength = strength[nTile];
								if (Double.isNaN(this_disp)){
									this_disp = 0.0;
									this_strength = 0.0;
								}
								for (int np = 0; np < surf_disparity.length; np++){
									if (shell_map[nsTile][np] > 0) {
										double disp_diff = normDispADiff (this_disp,surf_disparity[np], dispNorm);
										if ((best_indx < 0) || (disp_diff < best_disp_diff)) {
											best_disp_diff = disp_diff;
											best_indx =      np;
										}
									}
								}
								if ((this_strength > 0.0) || (snapZeroMode == 0)){
									if ((best_disp_diff <= snapDispAny) ||
											((best_disp_diff <= snapDispMax) && (best_disp_diff * this_strength < snapDispWeight )) ||
											((best_disp_diff <= snapNegAny) && (this_disp <  surf_disparity[best_indx]))  // farther than surface - higher tolerance
											){
										snap_disp[nTile] = surf_disparity[best_indx]; // snapped (default - keep source disparity
									} else if ((this_strength < snapStrengthAny) && (snapZeroMode != 0)){ // treat as zero strength?
										this_strength = 0.0;
									}
								}

								if ((this_strength == 0) && (snapZeroMode > 0)) {
									if (snapZeroMode == 1) {  // zero strength and non-zero snapZeroMode, consider 1 (strongest shell)
										best_indx = -1;
										double best_strength = Double.NaN;
										for (int np = 0; np < surf_disparity.length; np++){
											if (shell_map[nsTile][np] > 0) {
												double w = planes[nsTile][np].getWeight();
												if ((best_indx < 0) || (w >  best_strength)) {
													best_strength = w;
													best_indx =     np;
												}
											}
										}
										snap_disp[nTile] = surf_disparity[best_indx]; // snapped (default - keep source disparity
									} else { // zero strength and snapZeroMode > 1 (now only 2) - farthest surface (smallest disparity)
										best_indx = -1;
										double smallest_disparity = Double.NaN;
										for (int np = 0; np < surf_disparity.length; np++){
											if (shell_map[nsTile][np] > 0) {
												if ((best_indx < 0) || ( surf_disparity[np] <  smallest_disparity)) {
													smallest_disparity = surf_disparity[np];
													best_indx =     np;
												}
											}
										}
										snap_disp[nTile] = surf_disparity[best_indx]; // snapped (default - keep source disparity
									}
								}
							}
						}
					}
				}
			};
		}
		ImageDtt.startAndJoin(threads);
		return snap_disp;
	}

	/**
	 * Get number of used surfaces (some may be empty)
	 * @param snap_surf surface indices (per-tile), 0 - unmatched
	 * @param selected boolean tile selection or null (use all tiles)
	 * @return maximal number of surface used
	 */
	public int getNumSurf(
			int []     snap_surf,
			boolean [] selected) // or null
	{
		final int tilesX =        tileProcessor.getTilesX();
		final int superTileSize = tileProcessor.getSuperTileSize();
		final int stilesX = (tilesX + superTileSize -1)/superTileSize;

		int num_surf = 0;
		for (int nTile = 0; nTile < snap_surf.length; nTile ++) if (snap_surf[nTile] > 0){
			int tX = nTile % tilesX;
			int tY = nTile / tilesX;
			int stX = tX / superTileSize;
			int stY = tY / superTileSize;
			int nsTile = stY * stilesX + stX;
			int np = snap_surf[nTile] - 1;
			if ((shell_map[nsTile][np] > num_surf) && ((selected == null) || selected[nTile])) {
				num_surf = shell_map[nsTile][np];
			}
		}
		return num_surf;
	}

	/**
	 * Get number of unmatched tiles
	 * @param snap_surf surface indices (per-tile), 0 - unmatched
	 * @param selected boolean tile selection or null (use all tiles)
	 * @return number of unmatched tiles
	 */
	public int getNumNotMatched(
			int []     snap_surf,
			boolean [] selected) // or null
	{
		int num_notmatched = 0;
		for (int i = 0; i < snap_surf.length; i++){
			if ((snap_surf[i] == 0) && ((selected == null) || selected[i])) {
				num_notmatched ++;
			}
		}
		return num_notmatched;
	}

	public int getNeibMask(
			int x,
			int y,
			int width,
			int height)
	{
		final int [] corn_side_neib = { // of +1 - up (N), +2 - up-right - NE, ... +0x80 - NW
				0b00011100,    // top left corner
				0b01111100,    // top middle
				0b01110000,    // top right
				0b00011111,    // middle left
				0b11111111,    // middle
				0b11110001,    // middle right
				0b00000111,    // bottom left
				0b11000111,    // bottom middle
				0b11000001};   // bottom right
		int tileType = getNeibType(
				x,
				y,
				width,
				height);
		return corn_side_neib[tileType];
	}

	/**
	 * Returns 8-element array - which supertile neighbor enables this tile neighbor (-1 - none)
	 * @param x horizontal position of the tile in supertile
	 * @param y vertical position of the tile in supertile
	 * @param width supertile width
	 * @param height supertile height
	 * @return 8-element array
	 */

	public int[] getNeibBordDir(
			int x,
			int y,
			int width,
			int height)
	{
		final int [][] dir_dir = { // of +1 - up (N), +2 - up-right - NE, ... +0x80 - NW
			    //N  NE   E  SE   S  SW   W  NW
				{ 0,  0, -1, -1, -1,  6,  6,  7},    // top left corner
				{ 0,  0, -1, -1, -1, -1, -1,  0},    // top middle
				{ 0,  1,  2,  2, -1, -1, -1, -1},    // top right
				{-1, -1, -1, -1, -1,  6,  6,  6},    // middle left
				{-1, -1, -1, -1, -1, -1, -1, -1},    // middle
				{ 0,  2,  2,  2, -1, -1, -1, -1},    // middle right
				{-1, -1, -1,  4,  4,  5,  6,  6},    // bottom left
				{-1, -1, -1,  4,  4,  4, -1, -1},    // bottom middle
				{-1,  2,  2,  3,  4,  4, -1, -1}};   // bottom right
		int tileType = getNeibType(
				x,
				y,
				width,
				height);
		return dir_dir[tileType];
	}

	public int getNeibType(
			int x,
			int y,
			int width,
			int height)
	{
		int tileType = 4;
		int hm1 = height -1;
		int wm1 = width - 1;
		if (y == 0){
			if (x == 0){
				tileType = 0;
			} else if (x == wm1) {
				tileType = 2;
			} else {
				tileType = 1;
			}
		} else if (y == hm1) {
			if (x == 0){
				tileType = 6;
			} else if (x == wm1) {
				tileType = 8;
			} else {
				tileType = 7;
			}
		} else {
			if (x == 0){
				tileType = 3;
			} else if (x == wm1) {
				tileType = 5;
			}
		}
		return tileType;

	}

/**
 * Build tile connections to neighbors form optional selection, per-tile plane indices,
 * and a surface index.
 * @param selected tile selection (or null)
 * @param snap_surf per tile array of supertile plane index plus 1 or 0 for unmatched tiles
 * @param surf_index number of the surface to extract
 * @return per-tile array of -256 for unselected tiles or 8-bit bit mask of the neighbors (0 - N, 1 - NE, ...)
 */
	public int [] getNeighbors( // creates neighbors mask from bitmask
			final boolean [] selected, // or null
			final int []     snap_surf, // use this size - it matches image, not supertiles
			final int        surf_index)
	{
// TODO: see which are needed
		final int tilesX =        tileProcessor.getTilesX();
		final int tilesY =        tileProcessor.getTilesY();
		final int superTileSize = tileProcessor.getSuperTileSize();
		final int stilesX = (tilesX + superTileSize -1)/superTileSize;
//		final int stilesY = (tilesY + superTileSize -1)/superTileSize;
//		final double [][] surfaces = getSurfaces();
		final int [][] shell_map =   getShellMap();
		final TilePlanes.PlaneData [][] planes = getPlanesMod(); // to get strength and connections
		final int [] neibs = new int [selected.length];
		final int [] dirs8 =   {-tilesX,   -tilesX +  1, 1, tilesX + 1,  tilesX,  tilesX - 1,  -1, -tilesX -  1};
//		final int [] dirs8st = {-stilesX,  -stilesX + 1, 1, stilesX + 1, stilesX, stilesX - 1, -1, -stilesX - 1};
		final Thread[] threads = ImageDtt.newThreadArray(tileProcessor.threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		final int surf_index_plus1 = surf_index + 1;
		final int emptyTile = -256;
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				@Override
				public void run() {
					for (int nTile = ai.getAndIncrement(); nTile < snap_surf.length; nTile = ai.getAndIncrement()) {
						neibs[nTile] = emptyTile;
						if ((selected == null) || selected[nTile]) {
							int np = snap_surf[nTile] - 1;
							if (np >= 0) {
								int tX = nTile % tilesX;
								int tY = nTile / tilesX;
								int stX = tX / superTileSize;
								int stY = tY / superTileSize;
								int fX = tX % superTileSize;
								int fY = tY % superTileSize;
								int nsTile = stY * stilesX + stX;
								if (nsTile == 724){
//									System.out.println("getNeighbors(): tX:tY = "+tX+":"+tY+", stX:stY="+stX+":"+stY+
//											", nTile="+nTile+", nsTile="+nsTile);
								}
								if (shell_map[nsTile][np] == surf_index_plus1) {
									int tileMask = getNeibMask(
											tX,
											tY,
											tilesX,
											tilesY);
									int localMask = getNeibMask(
											fX,
											fY,
											superTileSize,
											superTileSize);
									int [] dir_dir = getNeibBordDir(
											fX,
											fY,
											superTileSize,
											superTileSize);
									neibs[nTile] = tileMask;
									for (int i = 0; i < 8; i++) {
										int b = 1 << i;
										if ((neibs[nTile] & b) != 0){ // that will prevent tile array bounds violation
											int nTile1 = nTile + dirs8[i];
											if (!(((selected == null) || selected[nTile1]) && (snap_surf[nTile1] == snap_surf[nTile]))){
												neibs[nTile] &= ~b; // same supertile, different index
											} else if ((localMask & b) == 0) {// going to the other supertile
												int [] neib_best = planes[nsTile][np].getNeibBest(); // array of connected planes indices/ -1s
												if ((dir_dir[i] < 0) || (neib_best[dir_dir[i]] != (snap_surf[nTile1] -1))) {
													neibs[nTile] &= ~b;
												}
											}
										}
									}
								}
							}
						}
					}
				}
			};
		}
		ImageDtt.startAndJoin(threads);
		return neibs;

	}




	/**
	 * Sort per-tile disparity map values to one of the surfaces
	 * @param disparity per-tile array of disparity values
	 * @param strength  per-tile array of correlation strength values
	 * @param selection tile selection (or null to process all tiles)
	 * @param dispNorm  disparity difference normalization - differences for average disparities
	 *                  below this value are unchanged, above - proportionally scaled down)
	 * @param snapDispAny tiles with disparity difference below this will snap to the nearest
	 *                  surface for any correlation strength
	 * @param snapStrengthAny tiles that have lower strength and can not be snapped otherwise are
	 *                  treated as zero-strength and obey snapZeroMode parameter
	 * @param snapNegAny allow snap for tiles that are below (farther) than the closest surface,
	 *                  if the disparity difference is below this value
	 * @param snapDispMax do not snap low-strength tiles if they are farther from the nearest surface
	 *                  than this.
	 * @param snapDispWeight snap tiles in the range snapDispAny and snapDispMax if the product
	 * 					of the disparity difference by the strength is below this
	 * @param snapZeroMode : 0 - no special treatment, 1 snap to strongest surface,
	 *                  2 - snap to the farthest (lowest disparity) surface
	 * @param debugLevel debug level
	 * @return integer array with per-tile elements. 0 - no surface candidate found (keep original
	 *                  disparity), > 0 - number of the plane plus 1 (use shell_map to get
	 *                  the surface number)
	 */

	public int [] snapSort(
			final double []  disparity,
			final double []  strength,
			final boolean [] selection,       // can be null
			final double     dispNorm, // plDispNorm
			final double     snapDispAny,    //    =  .2,   // Maximal (scaled by plDispNorm) disparity difference to snap to plane at any strength
			final double     snapStrengthAny,    //    =  .2,   //  Maximal strength to fit any distance (if does not fit otherwise - treat as zero strength
			final double     snapNegAny,     //    =  .2,   // Maximal negative disparity difference from the best match
			final double     snapDispMax,    //    =  .5,   // Maximal (scaled by plDispNorm) disparity difference to snap to plane at low strength
			final double     snapDispWeight, //    =  .5,   // Maximal disparity diff. by weight product to snap to plane
			final int        snapZeroMode,   // Zero strength snap mode: 0: no special treatment, 1 - strongest, 2 - farthest
			final int        debugLevel)
	{
		final int tilesX =        tileProcessor.getTilesX();
//		final int tilesY =        tileProcessor.getTilesY();
		final int superTileSize = tileProcessor.getSuperTileSize();
		final int stilesX = (tilesX + superTileSize -1)/superTileSize;
//		final int stilesY = (tilesY + superTileSize -1)/superTileSize;
		final double [][] surfaces = getSurfaces();
		final int [][] shell_map =   getShellMap();
		final TilePlanes.PlaneData [][] planes = getPlanesMod(); // to get strength;
//		final double [] snap_disp =  disparity.clone();
		final int [] snap_sort =  new int [disparity.length];

		final Thread[] threads = ImageDtt.newThreadArray(tileProcessor.threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				@Override
				public void run() {
					for (int nTile = ai.getAndIncrement(); nTile < disparity.length; nTile = ai.getAndIncrement()) {
						if ((selection == null) || selection[nTile]) {
							int tX = nTile % tilesX;
							int tY = nTile / tilesX;
							int stX = tX / superTileSize;
							int stY = tY / superTileSize;
							int nsTile = stY * stilesX + stX;
							int nTile_st = (tY * stilesX * superTileSize) + tX; // surfaces maybe (is) larger than tilesX *tilesY if
							// tilesX, tilesY are not multiples of superTileSize !
							if (nsTile == 724){
//								System.out.println("snapDisparity(): tX:tY = "+tX+":"+tY+", stX:stY="+stX+":"+stY+
//										", nTile="+nTile+", nsTile="+nsTile+", nTile_st = "+nTile_st);
							}

							// is there any surface to snap to?
							boolean no_surf = true;
							if (shell_map[nsTile] != null){
								for (int np = 0; no_surf && (np < shell_map[nsTile].length); np++){
									if (shell_map[nsTile][np] > 0) no_surf = false;
								}
							}
							if (!no_surf) {
								double [] surf_disparity = new double[shell_map[nsTile].length];
								for (int np = 0; np < surf_disparity.length; np++){
									if (shell_map[nsTile][np] > 0) {
										surf_disparity[np] = surfaces[shell_map[nsTile][np] -1][nTile_st];
									} else {
										surf_disparity[np] = Double.NaN;
									}
								}
								int best_indx = -1;
								double best_disp_diff = Double.NaN;
								double this_disp = disparity[nTile];
								double this_strength = strength[nTile];
								if (Double.isNaN(this_disp)){
									this_disp = 0.0;
									this_strength = 0.0;
								}
								for (int np = 0; np < surf_disparity.length; np++){
									if (shell_map[nsTile][np] > 0) {
										double disp_diff = normDispADiff (this_disp,surf_disparity[np], dispNorm);
										if ((best_indx < 0) || (disp_diff < best_disp_diff)) {
											best_disp_diff = disp_diff;
											best_indx =      np;
										}
									}
								}
								if ((this_strength > 0.0) || (snapZeroMode == 0)){
									if ((best_disp_diff <= snapDispAny) ||
											((best_disp_diff <= snapDispMax) && (best_disp_diff * this_strength < snapDispWeight )) ||
											((best_disp_diff <= snapNegAny) && (this_disp <  surf_disparity[best_indx]))  // farther than surface - higher tolerance
											){
										snap_sort[nTile] = best_indx+1; // snapped (default - keep source disparity
									} else if ((this_strength < snapStrengthAny) && (snapZeroMode != 0)){ // treat as zero strength?
										this_strength = 0.0;
									}
								}

								if ((this_strength == 0) && (snapZeroMode > 0)) {
									if (snapZeroMode == 1) {  // zero strength and non-zero snapZeroMode, consider 1 (strongest shell)
										best_indx = -1;
										double best_strength = Double.NaN;
										for (int np = 0; np < surf_disparity.length; np++){
											if (shell_map[nsTile][np] > 0) {
												double w = planes[nsTile][np].getWeight();
												if ((best_indx < 0) || (w >  best_strength)) {
													best_strength = w;
													best_indx =     np;
												}
											}
										}
										snap_sort[nTile] = best_indx+1; // snapped (default - keep source disparity
									} else { // zero strength and snapZeroMode > 1 (now only 2) - farthest surface (smallest disparity)
										best_indx = -1;
										double smallest_disparity = Double.NaN;
										for (int np = 0; np < surf_disparity.length; np++){
											if (shell_map[nsTile][np] > 0) {
												if ((best_indx < 0) || ( surf_disparity[np] <  smallest_disparity)) {
													smallest_disparity = surf_disparity[np];
													best_indx =     np;
												}
											}
										}
										snap_sort[nTile] = best_indx+1;
									}
								}
							}
						}
					}
				}
			};
		}
		ImageDtt.startAndJoin(threads);
		return snap_sort;
	}

	public double [][][] getDisparityStrengths(
			int        stMeasSel) //            = 1;      // Select measurements for supertiles : +1 - combo, +2 - quad +4 - hor +8 - vert)
	{
		int numMeasLayers = measuredLayers.getNumLayers();
		double [][][] ds = new double[numMeasLayers][][];
		for (int ml = 0; ml <  numMeasLayers; ml++) if ((stMeasSel & ( 1 << ml)) != 0) {
			ds[ml] = 	measuredLayers.getDisparityStrengthML (
					ml, // int num_layer,
					this.mlfp);
//					this.mlfp.strength_floor, // double strength_floor,
//					this.mlfp.strength_pow); //  double strength_pow)
		}
		return ds;
	}

	public boolean [][] getMeasurementSelections(
			int        stMeasSel) //            = 1;      // Select measurements for supertiles : +1 - combo, +2 - quad +4 - hor +8 - vert)
	{
		int numMeasLayers = measuredLayers.getNumLayers();
		boolean [][] sels = new boolean[numMeasLayers][];
		for (int ml = 0; ml <  numMeasLayers; ml++) if ((stMeasSel & ( 1 << ml)) != 0) {
			sels[ml] = 	measuredLayers.getSelection (ml);
		}
		return sels;
	}

} // end of class SuperTiles
