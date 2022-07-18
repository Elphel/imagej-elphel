/**
 **
 ** OpticalFlow - Process scene pairs
 **
 ** Copyright (C) 2020 Elphel, Inc.
 **
 ** -----------------------------------------------------------------------------**
 **
 **  OpticalFlow.java is free software: you can redistribute it and/or modify
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
package com.elphel.imagej.tileprocessor;

import java.awt.Color;
import java.awt.Font;
import java.awt.Rectangle;
import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.IOException;
import java.io.ObjectOutputStream;
import java.io.Serializable;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Properties;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.atomic.DoubleAccumulator;
import javax.xml.bind.DatatypeConverter;

import org.apache.commons.math3.complex.Quaternion;

import java.util.concurrent.ThreadLocalRandom;

import com.elphel.imagej.cameras.CLTParameters;
import com.elphel.imagej.cameras.ColorProcParameters;
import com.elphel.imagej.cameras.EyesisCorrectionParameters;
import com.elphel.imagej.common.DoubleGaussianBlur;
import com.elphel.imagej.common.ShowDoubleFloatArrays;
import com.elphel.imagej.correction.CorrectionColorProc;
import com.elphel.imagej.gpu.GpuQuad;
import com.elphel.imagej.gpu.TpTask;
import com.elphel.imagej.jp4.JP46_Reader_camera;

import ij.ImagePlus;
import ij.ImageStack;
import ij.Prefs;
import ij.gui.Line;
import ij.io.FileSaver;
import ij.process.FloatProcessor;
import ij.process.ImageConverter;
import ij.process.ImageProcessor;
import ij.plugin.filter.AVI_Writer;
import ij.plugin.filter.GaussianBlur;

public class OpticalFlow {
	public static String [] COMBO_DSN_TITLES = {"disp", "strength","disp_lma","num_valid","change",
			"disp_bg", "strength_bg","disp_lma_bg","change_bg","disp_fg","disp_bg_all"};
	public static int COMBO_DSN_INDX_DISP =        0; // cumulative disparity (from CM or POLY), FG
	public static int COMBO_DSN_INDX_STRENGTH =    1; // strength, FG
	public static int COMBO_DSN_INDX_LMA =         2; // masked copy from 0 - cumulative disparity
	public static int COMBO_DSN_INDX_VALID =       3; // initial only
	public static int COMBO_DSN_INDX_CHANGE =      4; // increment
	public static int COMBO_DSN_INDX_DISP_BG =     5; // cumulative BG disparity (from CM or POLY)
	public static int COMBO_DSN_INDX_STRENGTH_BG = 6; // background strength
	public static int COMBO_DSN_INDX_LMA_BG =      7; // masked copy from BG disparity 
	public static int COMBO_DSN_INDX_CHANGE_BG =   8; // increment, BG
	public static int COMBO_DSN_INDX_DISP_FG =     9; // cumulative disparity (from CM or POLY), FG
	public static int COMBO_DSN_INDX_DISP_BG_ALL =10; // cumulative BG disparity (Use FG where no BG is available)

	
	
	public static double [] ZERO3 = {0.0,0.0,0.0};
	public static double  LINE_ERR = 0.1;
	public int            threadsMax = 100;  // maximal number of threads to launch
	public boolean        updateStatus = true;
	public int            numSens;
    // not configured yet
    public double         scale_no_lma_disparity = 1.0; // multiply strength were disparity_lma = NaN; 
	public long                                            startTime;     // start of batch processing

	public OpticalFlow (
			int            numSens,
			double         scale_no_lma_disparity,
			int            threadsMax,  // maximal number of threads to launch
			boolean        updateStatus) {
		this.numSens =                numSens;
		this.threadsMax =             threadsMax;
		this.updateStatus =           updateStatus;
		this.scale_no_lma_disparity = scale_no_lma_disparity; // currently not used (applied directly to setDSRBG) May be removed
	}

	/**
	 * Fill gaps in scene tile values (encoded as Double.NaN) from neighbors
	 * @param nan_tiles [macrotiles][layer(disparity, strength, r,b,g)][tile-in-macrotile], has nulls at first index
	 * @param qthis scene (QuadCLT instance)
	 * @param num_passes maximal number of passes to run Laplacian
	 * @param max_change threshold change of the tile value
	 * @param debug_level debug level
	 */
	
	public void fillTilesNans(
			final double [][][] nan_tiles,
			final QuadCLT     qthis,
			final int         num_passes,
			final double      max_change,
			final int         debug_level)
	{
		final double   diagonal_weight = 0.5 * Math.sqrt(2.0); // relative to ortho
		final int margin =               0; // 1; // extra margins over 16x16 tiles to accommodate distorted destination tiles
		final TileProcessor tp =         qthis.getTileProcessor();
		final int transform_size =       tp.getTileSize();
		final int fullTileSize =         2 * (transform_size + margin);
		double wdiag = 0.25 *diagonal_weight / (diagonal_weight + 1.0);
		double wortho = 0.25 / (diagonal_weight + 1.0);
		final double [] neibw = {wortho, wdiag, wortho, wdiag, wortho, wdiag, wortho, wdiag}; 
		
		final Thread[] threads = ImageDtt.newThreadArray(threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		final int dbg_mtile = (debug_level > 1)? 203 : -1;
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				public void run() {
					for (int iMTile = ai.getAndIncrement(); iMTile < nan_tiles.length; iMTile = ai.getAndIncrement()) {
						if (iMTile == dbg_mtile) {
							System.out.println("fillTilesNans (): iMTile = "+iMTile);
						}
						if (nan_tiles[iMTile] != null) {
							tilesFillNaN(
									neibw, // final double []   neibw,
									nan_tiles[iMTile], // final double [][] slices,
									num_passes, // final int         num_passes,
									max_change, // final double      max_change,
									fullTileSize); // final int         width
						}
					}
				}
			};
		}		      
		ImageDtt.startAndJoin(threads);
		
		if (debug_level > 1) {
			// show debug image
			String title = qthis.getImageName()+"-NO-NaN";
			showMacroTiles(
					title,        // String title,
					nan_tiles, // double [][][] source_tiles,
					qthis,        // final QuadCLT qthis,
					margin);      // final int     margin); // extra margins over 16x16 tiles to accommodate distorted destination tiles
		}
		if (debug_level > 0) {
			System.out.println("fillTilesNans() DONE.");
		}
	}

	// 
	/**
	 * Helper to be called from thread of fillTilesNans()
	 * @param neibw array of 8 weights of neighbors, CW, starting with up, sum = 1.0
	 * @param slices per channel (disparity, strength, r, b, g) array of macrotile slices 
	 * @param num_passes number of times to replace value by a weighted average of 8 neighbors
	 * @param max_change break if the absolute value of the change falls below this threshold
	 * @param width width of a macrotile in tiles (height is .length/width)
	 */
	
	private void tilesFillNaN(
			final double []   neibw,
			final double [][] slices,
			final int         num_passes,
			final double      max_change,
			final int         width
			) {
		final int tiles = slices[0].length;
		final int height = tiles/width;
		double [] strength = slices[QuadCLT.DSRBG_STRENGTH]; 
		final TileNeibs tn =  new TileNeibs(width, height);
		double [] slice_in =  new double [tiles];
		double [] slice_out = new double [tiles];
		//first - process strength, then use calculated strength to fill other slices
		boolean [] fixed = new boolean [tiles]; 
		for (int i = 0; i < tiles; i++) {
			if (strength[i] > 0.0) {
				fixed[i] = true;
			} else {
				strength[i] = 0.0; // get rid of NaN; for non-strength will use average
			}
		}
		System.arraycopy(strength, 0, slice_in, 0, tiles);
		for (int pass = 0; pass < num_passes; pass ++) {
			double pass_change = 0.0;
			for (int nt = 0; nt < tiles; nt++) if (!fixed[nt]){
				double s = 0.0;
				double sw = 0.0;
				double d;
				for (int dir = 0; dir < 8; dir++) {
					int nt1 = tn.getNeibIndex(nt, dir);
					if (nt1 >=0) {
						if (fixed[nt1]) {
							d = strength[nt1];
						}else {
							d = slice_in[nt1];
						}
						s += d * neibw[dir];
						sw += neibw[dir];
					}
				}
				// sw should never be 0;
				s /= sw;
				pass_change = Math.max(pass_change, Math.abs(slice_out[nt] - s));
				slice_out[nt] = s;
			}
			if (pass_change < max_change) {
				break;
			}
			System.arraycopy(slice_out, 0, slice_in, 0, tiles);
		}
		for (int i = 0; i < fixed.length; i++) if (!fixed[i]){
			strength[i] = slice_out[i];
		}
		//non-strength							
		for (int iSlice = 0; iSlice < slices.length; iSlice++) if (iSlice != QuadCLT.DSRBG_STRENGTH){
			double [] slice =    slices[iSlice];
			System.arraycopy(slice, 0, slice_in, 0, tiles);
			double fs =0.0;
			double fsw = 0.0;
			for (int i = 0; i < fixed.length; i++) {
				if (!Double.isNaN(slice[i])  &&  (strength[i] > 0.0)) { //  - now already non-null
					fixed[i] = true;
					fs +=  slice[i] * strength[i];
					fsw += strength[i];
				}
			}
			if (fsw <= 0.0) {
				continue; // should not happen
			}
			fs /= fsw; // average value
			for (int i = 0; i < fixed.length; i++) if (! fixed[i]){
				slice_in[i] = fs;
			}								
			for (int pass = 0; pass < num_passes; pass ++) {
				double pass_change = 0.0;
				for (int nt = 0; nt < tiles; nt++) if (!fixed[nt]){
					double s = 0.0;
					double sw = 0.0;
					double d;
					for (int dir = 0; dir < 8; dir++) {
						int nt1 = tn.getNeibIndex(nt, dir);
						if (nt1 >=0) {
							if (fixed[nt1]) {
								d = slice[nt1];
							}else {
								d = slice_in[nt1];
							}
							double w = neibw[dir]; //  * strength[nt1];
							s += d * w ;
							sw += w;
						}
					}
					if (sw > 0) {
						s /= sw;
					}
					pass_change = Math.max(pass_change, Math.abs(slice_out[nt] - s));
					slice_out[nt] = s;
				}
				if (pass_change < max_change) {
					break;
				}
				System.arraycopy(slice_out, 0, slice_in, 0, tiles);
			}
			for (int i = 0; i < fixed.length; i++) if (!fixed[i]){
				slice[i] = slice_out[i];
			}
		}
	}
	
	
	
	
	/**
	 * Calculate confidence for the interscene X,Y correlation
	 * @param flowXY per-macrotile array of per-tile X,Y of the optical flow vectors. May have nulls
	 * @param width width of the macrotile array 
	 * @param best_num select this number of tghe closest matches among 8 neighbors
	 * @param ref_stdev  confidence formula: (ref_stdev ^ 2) / (neib_std_dev^2 + ref_stdev^2)  
	 * @param debug_title debug image title null - no image)
	 * @return per-tile array of triplets {x,y, confidence}. May have nulls (original and those without enough neighbors
	 */
	
	public double [][] attachVectorConfidence(
			final double [][] flowXY,
			final int         width,
			final int         best_num,
			final double      ref_stdev,
			final String      debug_title)
	{
		int height = flowXY.length/width;
		final double [][] flowXYS = new double[flowXY.length][];
		final TileNeibs tn =  new TileNeibs(width, height);
		final double ref_stdev2 =ref_stdev * ref_stdev;
		final Thread[] threads = ImageDtt.newThreadArray(threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		final int dbg_mtile = -1; // 203;
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				public void run() {
					for (int iMTile = ai.getAndIncrement(); iMTile < flowXY.length; iMTile = ai.getAndIncrement()) if (flowXY[iMTile] != null){
						if (iMTile == dbg_mtile) {
							System.out.println("attachVectorConfidence(): iMTile = "+iMTile);
						}
						double [] r2s = new double [8];
						for (int dir = 0; dir < r2s.length; dir++) {
							int indx =  tn.getNeibIndex(iMTile, dir);
							if ((indx >= 0) &&  (flowXY[indx] != null)){
								double  dx = flowXY[indx][0] - flowXY[iMTile][0];
								double  dy = flowXY[indx][1] - flowXY[iMTile][1];
								r2s[dir] = dx*dx + dy*dy;
							} else {
								r2s[dir] =Double.NaN;
							}
						}
						Arrays.sort(r2s); //  Double.NaN is considered greater than any other value and all Double.NaN values are considered equal.
						if (!Double.isNaN(r2s[best_num-1])) {
							double s1=0.0, s2 =0.0;
							for (int i = 0; i < best_num; i++) {
								s1 += r2s[i];
								s2 += r2s[i] * r2s[i];
							}
							double sd2 = (best_num * s2 - s1*s1)/(best_num * best_num);
							double confidence =  (ref_stdev * ref_stdev) / (sd2 + ref_stdev2);
							flowXYS[iMTile] = new double[] {flowXY[iMTile][0],flowXY[iMTile][1],confidence};
						}
					}
				}
			};
		}
		
		ImageDtt.startAndJoin(threads);
		if (debug_title != null) {
			showVectorXYConfidence(
					debug_title, // String      title,
					flowXYS, // double [][] flowXYS,
					width); // int         width)	
		}
		return flowXYS;
	}
	
	private int removeOutliers(
			double nsigma,
			double [][] flowXYS)
	{
		if (nsigma < 0.0) {
			return 0;
		}
		double swx = 0.0, swy = 0.0, sw = 0.0;
		for (int i = 0; i < flowXYS.length; i++) if (flowXYS[i] != null) {
			double w = flowXYS[i][2];
			swx += flowXYS[i][0]* w;
			swy += flowXYS[i][1]* w;
			sw += w;
		}
		if (sw > 0.0) {
			swx /= sw;
			swy /= sw;
			// calculate deviation regardless of weight
			int n = 0;
			double s2 = 0.0;
			for (int i = 0; i < flowXYS.length; i++) if (flowXYS[i] != null) {
				double dx = flowXYS[i][0] - swx;
				double dy = flowXYS[i][1] - swy;
				s2 += dx*dx+dy*dy;
				n++;
			}
			s2/= n;
			n=0;
			double s2_max = s2 * nsigma * nsigma;
			for (int i = 0; i < flowXYS.length; i++) if (flowXYS[i] != null) {
				double dx = flowXYS[i][0] - swx;
				double dy = flowXYS[i][1] - swy;
				if ((dx*dx+dy*dy) > s2_max) {
					flowXYS[i] = null;
					n++;
				}
			}
			return n;
			/*
			double err = 4.0;
			for (int i = 0; i < flowXYS.length; i++) if (flowXYS[i] != null){
				flowXYS[i][0] += err;
				flowXYS[i][1] += err;
			}
			*/
			
		}
		return -1;
	}
	
	
	/**
	 * Show a 2/3-slice image for the optical flow (only X,Y or X,Y,Confidence)
	 * @param title image title (null or empty - will not show)
	 * @param flowXYS
	 * @param width number of macrotiles in a row
	 */
	private void showVectorXYConfidence(
			String      title,
			double [][] flowXYS,
			int         width)
	{
		if ((title != null) && !title.equals("")) {
			int height = flowXYS.length/width;
//			String [] titles0 ={"dX","dY","Strength"};
			String [] titles ={"dX","dY","Strength","Err","dX-weighted","dY-weighted","Werr"};
//			int nslices = titles0.length;
			int nslices = 0;
			for (int i = 0; i < flowXYS.length; i++) if (flowXYS[i] != null) {
				nslices = flowXYS[i].length;
				break;
			}
			/*
			if (nslices > titles0.length) {
				nslices = titles0.length;
			}

			String [] titles =new String [nslices];
			for (int i = 0; i < nslices; i++) {
				titles[i] = titles0[i]; 
			}
			*/
			final double [][] dbg_img = new double [titles.length][width * height];
			for (int l = 0; l < dbg_img.length; l++) {
				Arrays.fill(dbg_img[l],  Double.NaN);
			}
			for (int mtile = 0; mtile < flowXYS.length; mtile++) if (flowXYS[mtile] != null){
//				for (int l = 0; l < dbg_img.length; l++) {
				for (int l = 0; l < nslices; l++) {
					dbg_img[l][mtile] = flowXYS[mtile][l];
				}
				if (nslices > 2) {
					dbg_img[3][mtile] = Math.sqrt(
							flowXYS[mtile][0] * flowXYS[mtile][0] + 
							flowXYS[mtile][1] * flowXYS[mtile][1]);
					dbg_img[4][mtile] = flowXYS[mtile][0] * flowXYS[mtile][2];
					dbg_img[5][mtile] = flowXYS[mtile][1] * flowXYS[mtile][2];
					dbg_img[6][mtile] = Math.sqrt(
							flowXYS[mtile][0] * flowXYS[mtile][0] + 
							flowXYS[mtile][1] * flowXYS[mtile][1]) * flowXYS[mtile][2];
				}
			}
			(new ShowDoubleFloatArrays()).showArrays(
					dbg_img,
					width,
					height,
					true,
					title,
					titles);
		}
	}
	private void showVectorXYConfidence(
			String        title,
			double [][][] flowXYSs,
			int           width)
	{
		if ((title != null) && !title.equals("")) {
			int height = flowXYSs[0].length/width;
			String [] titles0 ={"dX","dY","S"};
			int nslices = titles0.length;
			for (int i = 0; i < flowXYSs[0].length; i++) if (flowXYSs[0][i] != null) {
				nslices = flowXYSs[0][i].length;
				break;
			}
			if (nslices > titles0.length) {
				nslices = titles0.length;
			}

			String [] titles =new String [nslices];
			for (int i = 0; i < nslices; i++) {
				titles[i] = titles0[i]; 
			}
			
			int slices = flowXYSs.length;
			String [] dbg_titles = new String[slices * nslices];
			double [][] dbg_img = new double [dbg_titles.length][width*height];
			for (int slice = 0; slice < flowXYSs.length; slice++) {
				for (int n = 0; n < nslices; n++) {
					dbg_titles[slice + n * slices] = titles0[n]+slice;
					Arrays.fill(dbg_img[slice + slices * n], Double.NaN);
					for (int i = 0; i < flowXYSs[slice].length; i++) if (flowXYSs[slice][i] != null){
						dbg_img[slice + n * slices][i] = flowXYSs[slice][i][n];
					}
				}

			}
			(new ShowDoubleFloatArrays()).showArrays(
					dbg_img,
					width,
					height,
					true,
					title,
					dbg_titles);			
		}
	}
	
	
	
	
	/**
	 * Calculate optical flow vectors for macrotiles (vX, vY, confidence) for each non-null macrotile
	 * by multiple iterations of 2D phase correlations, finding 2D argmax (currently by windowed center of masses)
	 * and adding the corrections to the initial offsets (flowXY). Correlation takes place between reference and scene tiles,
	 * where reference macrotiles use original scene data without any interpolation, and scene macrotiles try to
	 * minimize interpolation by finding the best-fit offset in the [-0.5,0.5) range for each of X and Y directions,
	 * and then applying residual fractional shifts (flowXY_frac) as rotations in the frequency domain.
	 *   
	 * @param scene_xyz Scene X (right),Y (up), Z (negative away from camera) in the reference camera coordinates
	 *        or null to use scene instance coordinates.
	 * @param scene_atr Scene azimuth, tilt and roll (or null to use scene instance).
	 * @param scene_QuadClt Scene QuadCLT instance.
	 * @param reference_QuadClt Reference QuadCLT instance.
	 * @param reference_tiles_macro 
	 * @param reference_center_occupancy Fraction of non-null tiles in the center 8x8 area of the reference macrotiles after disparity
	 *        filtering (see tolerance_absolute,  tolerance_relative). Below this threshold - skip that macrotile.
	 * @param flowXY0 Initial offset of scene tiles (in image pixels) in x (right) and y (down) directions or null (to use all zeros)
	 * @param tolerance_absolute Filter reference macrotiles by same disparity (within a disparity range) consisting of the sum 
	 *        of absolute disparity (tolerance_absolute) and a proportional to the average disparity (tolerance_relative).  
	 * @param tolerance_relative Relative to the average disparity part of the disparity filtering.
	 * @param scene_macrotile_occupancy Skip scene macrotile if less than this fraction of all tiles in a macrotile remain
	 *        after filtering.
	 * @param num_laplassian Number of Laplassian passes while replacing undefined (NaN) tiles from neighbors for reference and scene
	 *        macrotiles.
	 * @param change_laplassian Break the loop of Laplassian passes if the maximal absolute value of the last pass changes falls below
	 *        this threshold.
	 * @param chn_weights A 4-element array of the correlation weights for strength, red, blue and green channels
	 * @param corr_sigma A low-pass sigma for the 2-d correlation (in tiles)
	 * @param fat_zero 2D correlation relative "fat zero" to damp phase correlation normalization.
	 * @param late_normalize True - normalize after combining all channels, false - normalize each channel separately.
	 *        When consolidating multiple tile late_normalize is considered true. 
	 * @param iradius Used for argmax() center of mass window (1 - 3x3, 2 - 5x5)
	 * @param dradius Radius argmax() window radius
	 * @param refine_num For argmax() based on center of masses. In each iteration center window around previously found
	 *        fractional-pixel argmax. 
	 * @param num_ignore_worsening During first num_ignore_worsening iterations, do not reduce applied correction
	 * @param max_tries Limit of the number of correlation iterations.
	 * @param magic_scale 0.85 - measured argmax has a bias caused by fading of the 2D correlation away from the center.
	 * @param min_change Stop refining offset vector when the correction falls below this value (in image pixels)
	 * @param best_num When calculating the confidence level, calculate correction vector standard deviation among
	 *        8 neighbors and use best_num of the best (closest to the current macrotile) of them. Disregard macrotile
	 *        if number of available neighbors is less than best_num. 
	 * @param ref_stdev Average/expected standard deviation of differences between the current macrotile and its best
	 *        neighbors. Confidence formula is confidence= (ref_stdev ^ 2)/(ref_stdev ^2 + stdev^2)
	 * @param debug_level
	 * @return An array of per-macrotile triplets: {X, Y, Confidence}, where X and Y are expressed in image pixels.
	 *         some macrotiles may have nulls.
	 *         Return depends on threads(?), has random variations when restarted with the same data
	 */
	public double [][] correlate2DIterate( // returns optical flow and confidence
			final ImageDttParameters  imgdtt_params,
			// for prepareSceneTiles()			
			final double []   scene_xyz,     // camera center in world coordinates
			final double []   scene_atr,     // camera orientation relative to world frame
			final QuadCLT     scene_QuadClt,
			final QuadCLT     reference_QuadClt,
			final double [][] reference_tiles_macro,
			final double      reference_center_occupancy,   // fraction of remaining  tiles in the center 8x8 area (<1.0)
//			final double [][][] reference_tiles, // prepared with prepareReferenceTiles() + fillTilesNans();
			// flowXY should be initialized to all pairs of zeros (or deliberate pixel offset pairs if initial error is too high, will be modified with each iteration
			final double [][] flowXY0, // per macro tile initial {mismatch in image pixels in X and Y directions} or null
			final double      tolerance_absolute, // absolute disparity half-range in each tile
			final double      tolerance_relative, // relative disparity half-range in each tile
			final double      scene_macrotile_occupancy,          // fraction of remaining  tiles (<1.0)
			final int         num_laplassian,
			final double      change_laplassian,
			// for correlate2DSceneToReference ()
			final double []   chn_weights, // absolute, starting from strength (strength,r,b,g)
			final double      corr_sigma,
			final double      fat_zero,
			final boolean     late_normalize,
			// for correlation2DToVectors_CM()
			final int         iradius,      // half-size of the square to process 
			final double      dradius,      // weight calculation (1/(r/dradius)^2 + 1)
			final int         refine_num,   // number of iterations to apply weights around new center
			
			final int         num_ignore_worsening, // run all tiles for few iterations before filtering
			final int         max_tries,
			// for recalculateFlowXY()
			final double      magic_scale, // 0.85 for CM
			final double      min_change,
			
			final int         best_num,
			final double      ref_stdev,
			final int         debug_level, //1
			final boolean     enable_debug_images) // true
	{
		boolean debug_mismatch = (debug_level > -10); // && enable_debug_images);
//		int debug_corr2d = 0;// 10
		
		final TileProcessor tp =         reference_QuadClt.getTileProcessor();
		final int transform_size =       tp.getTileSize();

		int dbg_mtilesX = tp.getTilesX()/transform_size;
		int dbg_mtilesY = tp.getTilesY()/transform_size;
		int dbg_mtiles = dbg_mtilesX * dbg_mtilesY; 
		int dbg_width = 3*256; // largest
		int dbg_height= dbg_mtiles * 5;
		
//		double [][] dbg_img = debug_mismatch ? (new double [max_tries][dbg_width*dbg_height]) : null; 
		double [][] dbg_img = debug_mismatch ? (new double [max_tries][]) : null; 
		
		
		double [][][] reference_tiles = prepareReferenceTiles(
				reference_QuadClt,        // final QuadCLT     qthis,
				tolerance_absolute, // final double      tolerance_absolute, // absolute disparity half-range in each tile
				tolerance_relative, // final double      tolerance_relative, // relative disparity half-range in each tile
				reference_center_occupancy,   // final double      center_occupancy,   // fraction of remaining  tiles in the center 8x8 area (<1.0)
				-1); // -1); // 2); // final int         debug_level)
		
		fillTilesNans(
				reference_tiles,          // final double [][][] nan_tiles,
				reference_QuadClt,                 // final QuadCLT     qthis,
				num_laplassian,            // final int         num_passes,
				change_laplassian,            // final double      max_change,
				-1); //-1); // 2);                    // final int         debug_level)
		if (reference_tiles_macro != null) {
			double [][] macro_centers =  getMacroPxPyDisp(
					reference_QuadClt, // final QuadCLT     reference_QuadClt,
					reference_tiles    //final double [][][] reference_tiles // prepared with prepareReferenceTiles() + fillTilesNans();
					);
			for (int i = 0; i < reference_tiles_macro.length; i++) {
				reference_tiles_macro[i] = macro_centers[i];
			}
		}
		
		final double [][] flowXY = (flowXY0 == null) ? (new double [reference_tiles.length][2]):flowXY0;
		final double [][] flowXY_frac = new double [reference_tiles.length][]; // Will contain fractional X/Y shift for CLT
		double [][] flowXY_run = flowXY; // only non-nulls for the tiles to correlate
		final double []   abs_change = new double [reference_tiles.length]; // updated 
		Arrays.fill(abs_change, Double.NaN);
		final double [] step_scale =  new double [reference_tiles.length]; // multiply increment if change exceeds previous
		Arrays.fill(step_scale, 1.0);
		final double [][] flowXY_prev =  new double [reference_tiles.length][]; // multiply increment if change exceeds previous
		
		
		for (int ntry = 0; ntry < max_tries; ntry++) {
			double [][][] scene_tiles = prepareSceneTiles(// to match to reference
					// null for {scene,reference}{xyz,atr} uses instances globals 
					scene_xyz,                // final double []   scene_xyz,     // camera center in world coordinates
					scene_atr,                // final double []   scene_atr,     // camera orientation relative to world frame
					scene_QuadClt,            // final QuadCLT     scene_QuadClt,
					reference_QuadClt,        // final QuadCLT     reference_QuadClt,
					reference_tiles,          // final double [][][] reference_tiles, // prepared with prepareReferenceTiles() + fillTilesNans();
					flowXY_run,               // final double [][] flowXY, // per macro tile {mismatch in image pixels in X and Y directions
					flowXY_frac,              // final double [][] flowXY_frac, // should be initialized as [number of macro tiles][] - returns fractional shifts [-0.5, 0.5)
					tolerance_absolute,       // final double      tolerance_absolute, // absolute disparity half-range in each tile
					tolerance_absolute,       // final double      tolerance_relative, // relative disparity half-range in each tile
					tolerance_relative,       // final double      occupancy,          // fraction of remaining  tiles (<1.0)
					num_laplassian,               // final int         num_passes,
					change_laplassian,               // final double      max_change,
					-1); //-1); // 1); // 2);                       // final int         debug_level)
			// undefine tiles in flowXY that are never used
			if (ntry == 0) {
				for (int i = 0; i <flowXY.length; i++) {
					if ((scene_tiles[i] == null) || (reference_tiles[i] == null)) {
						flowXY[i] = null;	
					}
				}
			}
			double [][] corr2dscene_ref = correlate2DSceneToReference(// to match to reference
					imgdtt_params, // final ImageDttParameters  imgdtt_params,   // Now just extra correlation parameters, later will include, most others
					scene_QuadClt,          // final QuadCLT     scene_QuadClt,
					reference_QuadClt,      // final QuadCLT     reference_QuadClt,
					scene_tiles,            // final double [][][] scene_tiles,     // prepared with prepareSceneTiles()
					reference_tiles,        // final double [][][] reference_tiles, // prepared with prepareReferenceTiles() + fillTilesNans(); - combine?
					flowXY_frac,            // final double [][] flowXY_frac, // X, YH fractional shift [-0.5,0.5) to implement with FD rotations
					chn_weights,            // final double []   chn_weights, // absolute, starting from strength (strength,r,b,g)
					corr_sigma,             // final double      corr_sigma,
					fat_zero,               // final double      fat_zero,
					late_normalize,         // final boolean     late_normalize,
					false,                  // final boolean     combine_empty_only, // only use neighbor correlations for empty tiles (false - any)
					0.0,                    // final double      combine_dradius
					0.0,                    // final double      tolerance_absolute, // absolute disparity half-range to consolidate tiles
					0.0,                    // final double      tolerance_relative, // relative disparity half-range to consolidate tiles
					-1); // 1); // final int         debug_level)
			
			double [][] vectorsXYS = correlation2DToVectors_CM(
					corr2dscene_ref,        // final double [][] corr2d_tiles, // per 2d calibration tiles (or nulls)
					transform_size,         // final int         transform_size,
					iradius,                // final int         iradius,      // half-size of the square to process 
					dradius,                // final double      dradius,      // weight calculation (1/(r/dradius)^2 + 1)
					refine_num,             // final int         refine_num,   // number of iterations to apply weights around new center
					-1);                    //final int         debug_level)
			double      this_min_change = min_change; //  (ntry < num_run_all)? 0.0: min_change;
			boolean     ignore_worsening = ntry < num_ignore_worsening; // (num_run_all + 10);
			if (debug_level > 0) { //-2) { // was >0
				System.out.println("======== NTRY "+ntry +" ========");
			}
			double [][] flowXY_dbg=null;
			double [][] flowXY_prev_dbg=null;
			double [][] vectorsXYS_dbg=null;
			double [] abs_change_dbg=null;
			double [] step_scale_dbg=null;
			if (dbg_img != null) {
				flowXY_dbg =      new double [flowXY.length][];
				flowXY_prev_dbg = new double [flowXY_prev.length][];
				vectorsXYS_dbg =  new double [vectorsXYS.length][];
				for (int ii = 0; ii < vectorsXYS.length; ii++) {
					if (flowXY[ii] != null) flowXY_dbg[ii] = flowXY[ii].clone();
					if (flowXY_prev[ii] != null) flowXY_prev_dbg[ii] = flowXY_prev[ii].clone();
					if (vectorsXYS[ii] != null) vectorsXYS_dbg[ii] = vectorsXYS[ii].clone();
				}
				abs_change_dbg = abs_change.clone();
				step_scale_dbg = step_scale.clone();
			}
			if (dbg_img != null) {//
				System.out.println("Before recalculateFlowXY() ntry = "+ntry+" debug_level="+debug_level);
			}
			flowXY_run = recalculateFlowXY(
					flowXY,                     // final double [][] flowXY, // will update
					flowXY_prev,                // final double [][] flowXY_prev, // previous flowXY (may be null for tiles)   
					vectorsXYS,                 // final double [][] corr_vectorsXY,
					abs_change,                 // final double []   abs_change, // updated
					step_scale,                 // final double []   step_scale, // multiply increment if change exceeds previous
					ignore_worsening,           // final boolean     boolean     ignore_worsening 
					magic_scale/transform_size, // final double      magic_scale, // 0.85 for CM
					this_min_change,            // final double      min_change,
					((dbg_img != null)? 2: debug_level));                         // final int         debug_level);
			if (dbg_img != null) {
				dbg_img[ntry] = new double [dbg_width*dbg_height];
				Arrays.fill(dbg_img[ntry], Double.NaN);
				for (int ii = 0; ii < scene_tiles.length; ii++) if (scene_tiles[ii] != null){
					for (int jj = 0; jj < scene_tiles[ii].length; jj++) {
						System.arraycopy(
								scene_tiles[ii][jj],
								0,
								dbg_img[ntry],
								((dbg_mtiles * 0) + ii)* dbg_width + jj * scene_tiles[ii][jj].length,
								scene_tiles[ii][jj].length); // 256
					}
				}
				for (int ii = 0; ii < corr2dscene_ref.length; ii++) if (corr2dscene_ref[ii] != null){
					System.arraycopy(
							corr2dscene_ref[ii],
							0,
							dbg_img[ntry],
							((dbg_mtiles * 1) + ii)* dbg_width,
							 corr2dscene_ref[ii].length); // 225
				}
				for (int ii = 0; ii < vectorsXYS.length; ii++) if (vectorsXYS[ii] != null){
					for (int jj = 0; jj < vectorsXYS[ii].length; jj++) {
						System.arraycopy(
								vectorsXYS[ii],
								0,
								dbg_img[ntry],
								((dbg_mtiles * 2) + ii)* dbg_width,
								vectorsXYS[ii].length); // 3
					}
				}
				if (flowXY_run != null) {
					for (int ii = 0; ii < flowXY_run.length; ii++) if (flowXY_run[ii] != null){
						System.arraycopy(
								flowXY_run[ii],
								0,
								dbg_img[ntry],
								((dbg_mtiles * 3) + ii)* dbg_width,
								flowXY_run[ii].length); // 2
					}
				}
				for (int ii = 0; ii < abs_change.length; ii++){
					if (flowXY_frac[ii] != null){
						System.arraycopy(
								flowXY_frac[ii],
								0,
								dbg_img[ntry],
								((dbg_mtiles * 4) + ii)* dbg_width,
								flowXY_frac[ii].length); // 2
					}					
					if (flowXY_prev[ii] != null){
						System.arraycopy(
								flowXY_prev[ii],
								0,
								dbg_img[ntry],
								((dbg_mtiles * 4) + ii)* dbg_width + 2, // 2 pixels right 
								flowXY_prev[ii].length); // 2
					}						
					dbg_img[ntry][((dbg_mtiles * 4) + ii)* dbg_width + 4] = abs_change[ii];
					dbg_img[ntry][((dbg_mtiles * 4) + ii)* dbg_width + 5] = step_scale[ii];

					if (flowXY_dbg[ii] != null){
						System.arraycopy(
								flowXY_dbg[ii],
								0,
								dbg_img[ntry],
								((dbg_mtiles * 4) + ii)* dbg_width + 6, // 2 pixels right 
								flowXY_dbg[ii].length); // 2
					}
					if (flowXY_prev_dbg[ii] != null){
						System.arraycopy(
								flowXY_prev_dbg[ii],
								0,
								dbg_img[ntry],
								((dbg_mtiles * 4) + ii)* dbg_width + 8, // 2 pixels right 
								flowXY_prev_dbg[ii].length); // 2
					}
					if (vectorsXYS_dbg[ii] != null){
						System.arraycopy(
								vectorsXYS_dbg[ii],
								0,
								dbg_img[ntry],
								((dbg_mtiles * 4) + ii)* dbg_width + 10, // 2 pixels right
								vectorsXYS_dbg[ii].length); // 2
					}
					dbg_img[ntry][((dbg_mtiles * 4) + ii)* dbg_width + 13] = abs_change_dbg[ii];
					dbg_img[ntry][((dbg_mtiles * 4) + ii)* dbg_width + 14] = step_scale_dbg[ii];
				}
			}
			if (flowXY_run == null) { // nothing to do left
				break; 
			}
		}
		if (dbg_img != null) {
			(new ShowDoubleFloatArrays()).showArrays(
					dbg_img,
					dbg_width,
					dbg_height,
					true,
					"lma_data_"+debug_level);
		}
		final int macroTilesX =          tp.getTilesX()/transform_size;

		String flowXYS_title =  (enable_debug_images && (debug_level > 0))?("vectorXYS_"+scene_QuadClt.getImageName()+"-ref"+reference_QuadClt.getImageName()):null;
		
		double [][] vectorXYConfidence =  attachVectorConfidence(
				flowXY,         // final double [][] flowXY,
				macroTilesX,    // final int         width,
				best_num,       // final int         best_num,
				ref_stdev,      // final double      ref_stdev,
				flowXYS_title); // final String      debug_title);    
	
		return vectorXYConfidence; // it is also in input arguments
	}
	
	/**
	 * Recalculate optical flow vectors from the new 2D correlation results 
	 * @param currentFlowXY Previous optical flow vectors (are not modified) in image pixels. May have null-s.
	 * @param corr_vectorsXY Results of the 2D correlation.
	 * @param magic_scale 0.85 - measured argmax has a bias caused by fading of the 2D correlation away from the center.
	 * @return Updated optical flow vectors in image pixels. May have null-s.
	 */
	
	double [][] recalculateFlowXY(
			final double [][] currentFlowXY,
			final double [][] corr_vectorsXY,
			final double      magic_scale) // 0.85 for CM
	{
		
		final Thread[] threads = ImageDtt.newThreadArray(threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		final double [][]   flowXY =  new double [currentFlowXY.length][];
		final int dbg_mtile = -620; // 453; // 500;
		final double rmagic_scale = 1.0/magic_scale;

		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				public void run() {
					for (int iMTile = ai.getAndIncrement(); iMTile < currentFlowXY.length; iMTile = ai.getAndIncrement())
						if ((currentFlowXY[iMTile] != null) && (corr_vectorsXY[iMTile] != null)){
							if (iMTile == dbg_mtile) {
								System.out.println("recalculateFlowXY(): iMTile = "+iMTile);
							}
							flowXY[iMTile]= new double[] {
									currentFlowXY[iMTile][0] + rmagic_scale * corr_vectorsXY[iMTile][0],
									currentFlowXY[iMTile][1] + rmagic_scale * corr_vectorsXY[iMTile][1]};
						}
				}
			};
		}		      
		ImageDtt.startAndJoin(threads);
		return flowXY;
	}
	/**
	 * Recalculate optical flow vector (in image pixels)
	 * @param flowXY current per-tile vectors (null for undefined), updated
	 * @param flowXY_prev  previous flowXY (may be null for tiles) 
	 * @param corr_vectorsXY correction vector from correlation to apply
	 * @param abs_change absolute value of last coordinate change for each tile
	 * @param step_scale  multiply increment if change exceeds previous
	 * @param ignore_worsening continue even if the change exceeds previous
	 * @param magic_scale divide correlation vector (typically 0.85/8 for CM argmax) 
	 * @param min_change minimal vector coordinate difference to repeat correlations
	 * @param debug_level if > 0; print number of tiles to correlate
	 * @return flowXY vectors only for tiles to be updated or null if no tiles left
	 */
	double [][] recalculateFlowXY(
			final double [][] flowXY, // will update
			final double [][] flowXY_prev, // previous flowXY (may be null for tiles)   
			final double [][] corr_vectorsXY,
			final double []   abs_change, // updated
			final double []   step_scale, // multiply increment if change exceeds previous
			final boolean     ignore_worsening, 
			final double      magic_scale, // 0.85 for CM
			final double      min_change,
			final int         debug_level)  
	{
		
		final Thread[] threads = ImageDtt.newThreadArray(threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		final double [][]   flowXY_task =  new double [flowXY.length][];
		final int dbg_mtile = (debug_level > 1)? 40 : -1; // 473; // 295; // 15/7 620; // 453; // 500;
		final double rmagic_scale = 1.0/magic_scale;
		final AtomicInteger aupdate = new AtomicInteger(0); //number of tiles to recalculate
		final double reduce_step = 0.5; //multiply step if calculated difference is larger thart the previous
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				public void run() {
					for (int iMTile = ai.getAndIncrement(); iMTile < flowXY.length; iMTile = ai.getAndIncrement()) {
  						if (iMTile == dbg_mtile) {
							System.out.println("recalculateFlowXY() 1: iMTile = "+iMTile);
						}
  						if (flowXY[iMTile] != null){
  							if (corr_vectorsXY[iMTile] == null) {
  								if (min_change <= 0.0) { // ==0.1
  									if (flowXY_prev[iMTile] != null) {
  										flowXY_prev[iMTile][0] = flowXY[iMTile][0];
  										flowXY_prev[iMTile][1] = flowXY[iMTile][1];
  									} else {
  										flowXY_prev[iMTile] = flowXY[iMTile].clone();
  									}
  									flowXY_task[iMTile] = flowXY[iMTile];
  									abs_change[iMTile] = Double.NaN;
  									aupdate.getAndIncrement();
  								}
  							} else { // if (corr_vectorsXY[iMTile] == null)
  								double dx = step_scale[iMTile] * rmagic_scale * corr_vectorsXY[iMTile][0];
  								double dy = step_scale[iMTile] * rmagic_scale * corr_vectorsXY[iMTile][1];
  								double new_diff = Math.sqrt(dx*dx + dy*dy);
  								
  								double last_change = abs_change[iMTile]; // may be NaN;
  								abs_change[iMTile] = new_diff;
  								
  								if ((debug_level >2) && (new_diff > last_change) && (min_change > 0.0)) {
  									System.out.println("recalculateFlowXY() 2: iMTile="+iMTile+", new_diff="+ new_diff+", last_change="+last_change);
  								}

  								
  								if ((debug_level > 1) && (iMTile == dbg_mtile))  {
//  									System.out.println(String.format("recalculateFlowXY() 2A: iMTile = %4d (%2d / %2d) flowXY = [%f/%f] step_scale = %8.6f dx = %f dy =  %f  abs= %f previous = %f, ignore_worsening = %b",
//  											iMTile, (iMTile %10), (iMTile / 10), flowXY[iMTile][0], flowXY[iMTile][1], step_scale[iMTile], dx,dy,new_diff, last_change,ignore_worsening));
  									System.out.println("recalculateFlowXY() 2A: iMTile = "+iMTile+
  											" flowXY = "+flowXY[iMTile][0]+"/"+flowXY[iMTile][0] +
  											" step_scale = "+step_scale[iMTile]+
  											"dx = "+dx+" dy =  "+dy+" abs= "+new_diff +
  											"previous = "+ last_change+" ignore_worsening ="+ignore_worsening+
  											" corr_vectorsXY[iMTile][0]="+corr_vectorsXY[iMTile][0]+" corr_vectorsXY[iMTile][1]="+corr_vectorsXY[iMTile][1]);
  								}
  								
  								if (new_diff < min_change) {
  									if (flowXY_prev[iMTile] != null) {
  										flowXY_prev[iMTile][0] = flowXY[iMTile][0];
  										flowXY_prev[iMTile][1] = flowXY[iMTile][1];
  									} else {
  										flowXY_prev[iMTile] = flowXY[iMTile].clone();
  									}
  									flowXY[iMTile][0] += dx;
  									flowXY[iMTile][1] += dy;
  	  								if ((debug_level > 1) && (iMTile == dbg_mtile))  {
  	  									System.out.println(String.format("recalculateFlowXY() 2B: iMTile = %4d (%2d / %2d) flowXY = [%f/%f]",
  	  											iMTile, (iMTile %10), (iMTile / 10), flowXY[iMTile][0], flowXY[iMTile][1]));
  	  								}
  								} else {
  									if (ignore_worsening || !(new_diff >= last_change)) { // better or ignore - continue iterations
  										//
  										if ((debug_level > 1) && (iMTile == dbg_mtile))  {
  											System.out.println(String.format("recalculateFlowXY() 3: iMTile = %4d (%2d / %2d) flowXY = [%f/%f] step_scale = %f dx = %f dy =  %f  abs= %f previous = %f CONTINUE",
  													iMTile, (iMTile %10), (iMTile / 10), flowXY[iMTile][0], flowXY[iMTile][1], step_scale[iMTile], dx,dy,new_diff, last_change));
  										}
  	  									if (flowXY_prev[iMTile] != null) {
  	  										flowXY_prev[iMTile][0] = flowXY[iMTile][0];
  	  										flowXY_prev[iMTile][1] = flowXY[iMTile][1];
  	  									} else {
  	  										flowXY_prev[iMTile] = flowXY[iMTile].clone();
  	  									}
  										flowXY[iMTile][0] += dx;
  										flowXY[iMTile][1] += dy;
  										flowXY_task[iMTile] = flowXY[iMTile]; // set to measure
  										abs_change[iMTile] = new_diff;
  										aupdate.getAndIncrement();
  	  	  								if ((debug_level > 1) && (iMTile == dbg_mtile))  {
  	  	  									System.out.println(String.format("recalculateFlowXY() 3A: iMTile = %4d (%2d / %2d) flowXY = [%f/%f]",
  	  	  											iMTile, (iMTile %10), (iMTile / 10), flowXY[iMTile][0], flowXY[iMTile][1]));
  	  	  								}
  									} else if ((new_diff >= last_change) && (min_change > 0)) { // worse - reduce step, but still apply
  	  	  								if ((debug_level > 1) && (iMTile == dbg_mtile))  {
  	  	  									System.out.println(String.format("recalculateFlowXY() 4A: iMTile = %4d (%2d / %2d) flowXY = [%f/%f]",
  	  	  											iMTile, (iMTile %10), (iMTile / 10), flowXY[iMTile][0], flowXY[iMTile][1]));
  	  	  								}
  										
   										if (debug_level > 2) {
   											System.out.println(String.format("recalculateFlowXY() 4: iMTile = %4d (%2d / %2d) flowXY = [%f/%f] step_scale = %f dx = %f dy =  %f  abs= %f previous = %f REDUCED STEP",
   													iMTile, (iMTile %40), (iMTile / 40), flowXY[iMTile][0], flowXY[iMTile][1], step_scale[iMTile], dx,dy,new_diff, last_change));
   										}
   										// do not update previous (it should be not null
  	  									if (flowXY_prev[iMTile] == null) { // should not happen
  	  										System.out.println("BUG!");
  	  										flowXY_prev[iMTile] = flowXY[iMTile].clone();
  	  									}
  	  									dx = flowXY[iMTile][0] - flowXY_prev[iMTile][0];
  	  									dy = flowXY[iMTile][1] - flowXY_prev[iMTile][1];
  	  									flowXY[iMTile][0] = flowXY_prev[iMTile][0];
  	  									flowXY[iMTile][1] = flowXY_prev[iMTile][1];
  	    								if ((debug_level > 1) && (iMTile == dbg_mtile))  {
  	    									System.out.println(String.format("recalculateFlowXY() 4B: iMTile = %4d (%2d / %2d) flowXY = [%f/%f] step_scale = %f dx = %f dy =  %f  abs= %f previous = %f",
  	    											iMTile, (iMTile %10), (iMTile / 10), flowXY[iMTile][0], flowXY[iMTile][1], step_scale[iMTile], dx,dy,new_diff, last_change));
  	    								}
  	  									
  	  									step_scale[iMTile] *= reduce_step;
   										dx *= reduce_step;
   										dy *= reduce_step; // was wrong 03/16/2022: // dx *= reduce_step;
   										
   										flowXY[iMTile][0] += dx;
   										flowXY[iMTile][1] += dy;
   		  								if ((debug_level > 1) && (iMTile == dbg_mtile))  {
   		  									System.out.println(String.format("recalculateFlowXY() 4C: iMTile = %4d (%2d / %2d) flowXY = [%f/%f] step_scale = %f dx = %f dy =  %f  abs= %f previous = %f",
   		  											iMTile, (iMTile %10), (iMTile / 10), flowXY[iMTile][0], flowXY[iMTile][1], step_scale[iMTile], dx,dy,new_diff, last_change));
   		  								}
   										flowXY_task[iMTile] = flowXY[iMTile]; // set to measure
   										abs_change[iMTile] = last_change; // restore previous step
   										aupdate.getAndIncrement();
  									}
  								}
  							}
						}
					}
				}
			};
		}		      
		ImageDtt.startAndJoin(threads);
		if (debug_level > 0) {
			System.out.println("  recalculateFlowXY(): tiles to correlate: "+aupdate.get());
		}
		if (aupdate.get() > 0) {
			return flowXY_task;
		} else {
			return null; // nothing to measure left
		}
	}
	
	/**
	 * Convert 2D correlation tiles to 2D argmax using center of masses (CM) method.
	 * @param corr2d_tiles Array of 2d correlation tiles (or nulls). Each tile is typically 225
	 *        (2 * transform_size-1) * (2 * transform_size-1). 
	 * @param transform_size CLT transform size (8)
	 * @param iradius Used for argmax() center of mass window (1 - 3x3, 2 - 5x5)
	 * @param dradius Radius argmax() window radius
	 * @param refine_num For argmax() based on center of masses. In each iteration center window around previously found
	 *        fractional-pixel argmax. 
	 * @param debug_level Debug level (now > 0 print for programmed tile, can have a breakpoint).
	 * @return A per-macrotile array of {X,Y,strength} triplets. X,Y are in tiles (not image pixels), may have nulls.
	 *         Strength is a ratio of (max - average)/stdev. There is no interpolation, so strength is influenced by
	 *         a fractional part of argmax.
	 */
	public double [][] correlation2DToVectors_CM(
			final double [][] corr2d_tiles, // per 2d correlation tiles (or nulls)
			final int         transform_size,
			final int         iradius,      // half-size of the square to process 
			final double      dradius,      // weight calculation (1/(r/dradius)^2 + 1)
			final int         refine_num,   // number of iterations to apply weights around new center
			final int         debug_level)
	{
		final Thread[] threads = ImageDtt.newThreadArray(threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		final double [][]   vectors_xys =    new double [corr2d_tiles.length][];
		final int dbg_mtile = (debug_level > 0) ? 620 : -1; // 453; // 500;

		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				public void run() {
					for (int iMTile = ai.getAndIncrement(); iMTile < corr2d_tiles.length; iMTile = ai.getAndIncrement()) if (corr2d_tiles[iMTile] != null) {
  						if (iMTile == dbg_mtile) {
							System.out.println("correlation2DToVectors_CM (): iMTile = "+iMTile);
						}
						vectors_xys[iMTile] = getCorrCenterXYS_CM(
								corr2d_tiles[iMTile], // double []   corr2d_tile,
								transform_size,       // int         transform_size,
								iradius,              // int         iradius,
								dradius,              // double      dradius,
								refine_num);          // int         refine_num); 
					}
				}
			};
		}		      
		ImageDtt.startAndJoin(threads);
		return vectors_xys;
	}
	
	
	/**
	 * Single-tile 2D correlation tile to 2D argmax pair using center of masses (CM) method.
	 * @param corr2d_tile A 2D correlation tile, typically of 225 elements: 
	 *        (2 * transform_size-1) * (2 * transform_size-1).
	 * @param transform_size CLT transform size (8)
	 * @param iradius Used for argmax() center of mass window (1 - 3x3, 2 - 5x5)
	 * @param dradius Radius argmax() window radius
	 * @param refine_num For argmax() based on center of masses. In each iteration center window around previously found
	 *        fractional-pixel argmax. 
	 * @return A {X,Y,strength} triplet. X,Y are in tiles (not image pixels), may have nulls.
	 *         Strength is a ratio of (max - average)/stdev. There is no interpolation, so strength is influenced by
	 *         a fractional part of argmax.
	 */
	private double [] getCorrCenterXYS_CM(
			double []   corr2d_tile,
			int         transform_size,
			int         iradius,
			double      dradius,
			int         refine_num) // [2 * iradius + 1][2 * iradius + 1] 
	{
		// strength - (maximum - average)/stdev?
		int corr_size = 2* transform_size - 1;
		int imax = 0;
		for (int i = 1; i < corr2d_tile.length; i++) {
			if (corr2d_tile[i] > corr2d_tile[imax]) {
				imax = i;
			}
		}
		double xMax = imax % corr_size;
		double yMax = imax / corr_size;
		double k2 = 1.0/dradius*dradius;
		for (int pass = 0; pass < refine_num; pass ++) {
			int iXMax = (int) Math.floor(xMax);
			int iYMax = (int) Math.floor(yMax);
			int iY0 = iYMax - iradius;     if (iY0 < 0) iY0 = 0;
			int iY1 = iYMax + iradius + 1; if (iY1 >= corr_size) iY1 = corr_size -1;
			int iX0 = iXMax - iradius;     if (iX0 < 0) iX0 = 0;
			int iX1 = iXMax + iradius + 1; if (iX1 >= corr_size) iX1 = corr_size -1;
			double s = 0.0, sx = 0.0, sy = 0.0;
			for (int iy =  iY0; iy <= iY1; iy++) {
				double r2y = (iy - yMax)*(iy - yMax); 
				for (int ix =  iX0; ix <= iX1; ix++) {
					double d = corr2d_tile[ix + iy * corr_size];
					double r2 = r2y + (ix - xMax)*(ix - xMax);
					double w = 1.0/(k2*r2 + 1);
					double wd = w * d;
					s += wd;
					sx += wd * ix;
					sy += wd * iy;
				}
			}
			xMax = sx/s;
			yMax = sy/s;
		}
		int iYMmax = (int) Math.round(yMax); 
		int iXMmax = (int) Math.round(xMax); 
		if (iYMmax < 0)	iYMmax = 0;
		if (iYMmax >= transform_size) iYMmax = transform_size -1;
		if (iXMmax < 0)	iXMmax = 0;
		if (iXMmax >= transform_size) iXMmax = transform_size -1;

		double dMax = corr2d_tile[iYMmax * corr_size + iXMmax]; // negative
		double s1=0.0, s2 =0.0;
		for (int i = 0; i < corr2d_tile.length; i++) {
			s1 += corr2d_tile[i];
			s2 += corr2d_tile[i] * corr2d_tile[i];
		}
		double avg = s1/corr2d_tile.length;
		double sd = Math.sqrt(corr2d_tile.length * s2 - s1*s1)/corr2d_tile.length;
		double strength = (dMax - avg)/sd;
		
		return new double [] {xMax - transform_size +1, yMax - transform_size +1, strength};
	}
	
	/**
	 * 2D correlation of scene to reference tiles.
	 *  
	 * @param scene_QuadClt scene QuadCLT instance.
	 * @param reference_QuadClt reference QuadCLT instance.
	 * @param scene_tiles Scene tiles (per macrotile, per channel (disparity, strength, r, b, g)
	 *  (2*transform_size) *  (2*transform_size), currently 256. 
	 * @param reference_tiles Reference tiles, same format as scene_tiles.
	 * @param flowXY_frac Per-tile fractional X,Y offsets in the range of [-0.5, 0.5)
	 * @param chn_weights A 4-element array of the correlation weights for strength, red, blue and green channels
	 * @param corr_sigma A low-pass sigma for the 2-d correlation (in tiles)
	 * @param fat_zero 2D correlation relative "fat zero" to damp phase correlation normalization.
	 * @param late_normalize True - normalize after combining all channels, false - normalize each channel separately.
	 *        When consolidating multiple tile late_normalize is considered true. 
	 * @param combine_empty_only If true, use neighbors consolidation for undefined tiles only.
	 * @param combine_dradius  Radius for the consolidation weights half-cosine (weight is zero outside of dradius. 
	 * @param tolerance_absolute Used for consolidation only absolute tolerance to the difference between the disparity
	 *        of the central reference macrotile and the scene macrotiles to consolidate.  
	 * @param tolerance_relative Relative disparity tolerance - add a product of tolerance_relative by the average
	 *        macrotile disparity to the tolerance_absolute for disparity filtering of the macrotiles.
	 * @param debug_level Debug level.
	 * @return Per macrotile correlation tiles (now 225-long). May have nulls for the empty tiles.
	 */
	public double [][] correlate2DSceneToReference(// to match to reference
			final ImageDttParameters  imgdtt_params,   // Now just extra correlation parameters, later will include, most others
			final QuadCLT     scene_QuadClt,
			final QuadCLT     reference_QuadClt,
			final double [][][] scene_tiles,     // prepared with prepareSceneTiles()
			final double [][][] reference_tiles, // prepared with prepareReferenceTiles() + fillTilesNans(); - combine?
			final double [][] flowXY_frac, // X, YH fractional shift [-0.5,0.5) to implement with FD rotations
			final double []   chn_weights0, // absolute, starting from strength (strength,r,b,g)
			final double      corr_sigma,
			final double      fat_zero,
			final boolean     late_normalize,
			final boolean     combine_empty_only, // only use neighbor correlations for empty corr tiles (false - any)
			// reference tile should still be defined
			final double      combine_dradius, // 1 - 3x3, 2 - 5x5
			final double      tolerance_absolute, // absolute disparity half-range to consolidate tiles
			final double      tolerance_relative, // relative disparity half-range to consolidate tiles
			final int         debug_level)
	// returns per-tile 2-d correlations (15x15)
	{
		final Thread[] threads = ImageDtt.newThreadArray(threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		final TileProcessor tp =         reference_QuadClt.getTileProcessor();
		final int tilesX =               tp.getTilesX();
		final int tilesY =               tp.getTilesY();
		final int transform_size =       tp.getTileSize();
		final int tile_length =          transform_size * transform_size;
		final int macroTilesX =          tilesX/transform_size;
		final int macroTilesY =          tilesY/transform_size;
		final int macroTiles =           macroTilesX * macroTilesY; 
		
		final double [][]   corr_tiles =    new double [macroTiles][];
		final double [][][] corr_tiles_TD = new double [macroTiles][][];
		
		final int dbg_mtile = (debug_level >0) ? 203 : -1;
		
		final int chn_offset = QuadCLT.DSRBG_STRENGTH; // start scene_tiles, reference tiles with this 2-nd index
		int dsrbg_len = reference_QuadClt.getDSRBG().length;
		final int num_channels = (chn_weights0.length < (dsrbg_len - QuadCLT.DSRBG_STRENGTH))? chn_weights0.length: (dsrbg_len - QuadCLT.DSRBG_STRENGTH);
		final double []   chn_weights = new double  [num_channels];
		for (int i = 0; i < num_channels; i++) {
			chn_weights[i] = chn_weights0[i];
		}
		final ImageDtt image_dtt = new ImageDtt(
				numSens,
				transform_size,
				imgdtt_params, // null, // FIXME: Needs  ImageDttParameters (clt_parameters.img_dtt),
				reference_QuadClt.isAux(),
				reference_QuadClt.isMonochrome(),
				reference_QuadClt.isLwir(),
				1.0);

		final double [] filter =     image_dtt.doubleGetCltLpfFd(corr_sigma);
		final int      combine_radius = (int) Math.floor(combine_dradius); // 1 - 3x3, 2 - 5x5
		final double [][] rad_weights = new double [2 * combine_radius + 1][2 * combine_radius + 1];
		for (int dY = -combine_radius; dY <= combine_radius; dY ++) {
			for (int dX = -combine_radius; dX <= combine_radius; dX ++) {
				rad_weights[dY + combine_radius][dX + combine_radius] =
						Math.cos(0.5 * Math.PI * dY / combine_dradius) *
						Math.cos(0.5 * Math.PI * dX / combine_dradius);
			}
		}
		final double [] avg_disparity_ref =   new double [macroTiles];
		final double [] avg_disparity_scene = new double [macroTiles];
		
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				public void run() {
					DttRad2 dtt = new DttRad2(transform_size);
					dtt.set_window(1);
					double [][][] clt_tiles_ref =   new double [num_channels][4][];
					double [][][] clt_tiles_scene = new double [num_channels][4][];
					
					Correlation2d corr2d = new Correlation2d(
							2, // numSens,// does it matter?
							transform_size,             // int transform_size,
							false,                      // boolean monochrome,
							false);                     //   boolean debug)
					
					for (int iMTile = ai.getAndIncrement(); iMTile < macroTiles; iMTile = ai.getAndIncrement()) {
						if (reference_tiles[iMTile] != null) { // to calculate average reference disparity
							double sw = 0.0, sdw=0.0;
							double [] disparity = reference_tiles[iMTile][QuadCLT.DSRBG_DISPARITY]; 
							double [] strength =  reference_tiles[iMTile][QuadCLT.DSRBG_STRENGTH]; 
							for (int i = 0; i < strength.length; i++) {
								if (!Double.isNaN(disparity[i]) && (strength[i] > 0.0)) {
									sw += strength[i];
									sdw += strength[i]*(disparity[i]);
								}
							}
							if (sw > 0.0) avg_disparity_ref[iMTile] = sdw/sw;
						}
						if (scene_tiles[iMTile] != null) { // to calculate average scene disparity
							double sw = 0.0, sdw=0.0;
							double [] disparity = scene_tiles[iMTile][QuadCLT.DSRBG_DISPARITY]; 
							double [] strength =  scene_tiles[iMTile][QuadCLT.DSRBG_STRENGTH]; 
							for (int i = 0; i < strength.length; i++) {
								if (!Double.isNaN(disparity[i]) && (strength[i] > 0.0)) {
									sw += strength[i];
									sdw += strength[i]*(disparity[i]);
								}
							}
							if (sw > 0.0) avg_disparity_scene[iMTile] = sdw/sw;
						}
						if ((scene_tiles[iMTile] != null) && (reference_tiles[iMTile] != null)) {
							if (iMTile == dbg_mtile) {
								System.out.println("correlate2DSceneToReference(): iMTile = "+iMTile);
							}
							// convert reference tile
							double [][][] fold_coeff_ref = dtt.get_shifted_fold_2d ( // get_shifted_fold_2d(
									transform_size,
									0.0,
									0.0,
									0); // debug level
							double [][][] fold_coeff_scene = dtt.get_shifted_fold_2d ( // get_shifted_fold_2d(
									transform_size,
									flowXY_frac[iMTile][0],
									flowXY_frac[iMTile][1],
									0); // debug level
							for (int chn = 0; chn < num_channels; chn++) {
								double [] tile_in_ref =   reference_tiles[iMTile][chn + chn_offset];
								double [] tile_in_scene = scene_tiles[iMTile][chn + chn_offset];
								// unfold and convert both reference and scene
								for (int dct_mode = 0; dct_mode < 4; dct_mode++) {
									clt_tiles_ref[chn][dct_mode] = dtt.fold_tile (tile_in_ref, transform_size, dct_mode, fold_coeff_ref);
									clt_tiles_ref[chn][dct_mode] = dtt.dttt_iv   (clt_tiles_ref[chn][dct_mode], dct_mode, transform_size);
									clt_tiles_scene[chn][dct_mode] = dtt.fold_tile (tile_in_scene, transform_size, dct_mode, fold_coeff_scene);
									clt_tiles_scene[chn][dct_mode] = dtt.dttt_iv   (clt_tiles_scene[chn][dct_mode], dct_mode, transform_size);
								}
								// Apply shift to scene only (reference is not shifted)
						        image_dtt.fract_shift(    // fractional shift in transform domain. Currently uses sin/cos - change to tables with 2? rotations
						        		clt_tiles_scene[chn], // double  [][]  clt_tile,
						        		flowXY_frac[iMTile][0],            // double        shiftX,
						        		flowXY_frac[iMTile][1],            // double        shiftY,
						                false); // debug);
							}
							if (late_normalize) {
								corr_tiles_TD[iMTile] = corr2d.correlateCompositeTD( // correlate, do not normalize, stay in TD
										clt_tiles_ref,   // double [][][] clt_data1,
										clt_tiles_scene, // double [][][] clt_data2,
										null,            // double []     lpf,
										1.0,             // double        scale_value, // scale correlation value
										chn_weights);    // double []     col_weights_in, // should have the same dimension as clt_data1 and clt_data2
							} else {
								corr_tiles[iMTile] = corr2d.correlateCompositeFD( // 
										clt_tiles_ref,   // double [][][] clt_data1,
										clt_tiles_scene, // double [][][] clt_data2,
										filter,          // double []     lpf,
										1.0,             // double        scale_value, // scale correlation value
										chn_weights,     // double []     col_weights_in, // should have the same dimension as clt_data1 and clt_data2
										fat_zero);       // double        fat_zero)

							}
						}
					}
				}
			};
		}		      
		ImageDtt.startAndJoin(threads);
		
		if (late_normalize) {
			ai.set(0);
			for (int ithread = 0; ithread < threads.length; ithread++) {
				threads[ithread] = new Thread() {
					public void run() {
						final TileNeibs tn =  new TileNeibs(macroTilesX, macroTilesY);

						double [][] corr_tile_2D = new double [4][tile_length];
						Correlation2d corr2d = new Correlation2d(
								numSens,
								transform_size,             // int transform_size,
								false,                      // boolean monochrome,
								false);                     //   boolean debug)
						// reference tile should not be null, scene = may be
						for (int iMTile = ai.getAndIncrement(); iMTile < macroTiles; iMTile = ai.getAndIncrement()) if (reference_tiles[iMTile] != null){ 
							if (true) { // !combine_empty_only  || (corr_tiles_TD[iMTile] == null)) {
								if (iMTile == dbg_mtile) {
									System.out.println("correlate2DSceneToReference() 2: iMTile = "+iMTile);
								}
								if ((combine_radius > 0) && (!combine_empty_only  || (corr_tiles_TD[iMTile] == null))) { // 
									for (int q = 0; q< 4; q++) {
										Arrays.fill(corr_tile_2D[q], 0.0);
									}
									double disp_tol = tolerance_absolute + tolerance_relative * avg_disparity_ref[iMTile];
									double sw = 0;
									int iMX = iMTile % macroTilesX;
									int iMY = iMTile / macroTilesX;
									for (int dY = -combine_radius; dY <= combine_radius; dY ++) {
										for (int dX = -combine_radius; dX <= combine_radius; dX ++) {
											int indx = tn.getIndex(iMX+dX, iMY + dY);
											if ((indx >= 0) && (corr_tiles_TD[indx] != null)) {
												if ((Math.abs(avg_disparity_scene[iMTile] - avg_disparity_ref[iMTile])) <= disp_tol) {
													double w = rad_weights[dY + combine_radius][dX + combine_radius];
													for (int q = 0; q < 4; q++) {
														for (int i = 0; i < tile_length; i++) {
															corr_tile_2D[q][i] += w * corr_tiles_TD[indx][q][i];
														}
													}
													sw+= w;
												}
											}
										}
									}
									if (sw  <= 0.0) {
										continue; // no non-null tiles around
									}
									double a = 1.0/sw;
									for (int q = 0; q < 4; q++) {
										for (int i = 0; i < tile_length; i++) {
											corr_tile_2D[q][i] *= a;
										}
									}
								} else {
									corr_tile_2D = corr_tiles_TD[iMTile]; // no need to clone, reference OK
								}
								corr2d.normalize_TD(
										corr_tile_2D,         // double [][] td,
										filter,          // double []   lpf, // or null
										fat_zero);       // double      fat_zero);

								corr_tiles[iMTile] = corr2d.convertCorrToPD(
										corr_tile_2D);     // double [][] td);
							}
						}
					}
				};
			}		      
			ImageDtt.startAndJoin(threads);
		}
		return corr_tiles;
	}
	
	/**
	 * Get width of the macrotiles array
	 * @param reference_QuadClt scene instance
	 * @return width of a macrotile array
	 */
	public int getMacroWidth(final QuadCLT     reference_QuadClt) {
		final TileProcessor tp =         reference_QuadClt.getTileProcessor();
		return   tp.getTilesX()/tp.getTileSize();
	}
	
	/**
	 * Get triplets of {pX, pY, disparity} for each reference macrotile to use with LMA fitting
	 * @param reference_QuadClt reference scene instance
	 * @param reference_tiles reference tiles prepared with prepareReferenceTiles()
	 * @return Array of [macrotile]{pX, pY, disparity}, some macrotiles may be null
	 */
	public double [][] getMacroPxPyDisp(
			final QuadCLT     reference_QuadClt,
			final double [][][] reference_tiles // prepared with prepareReferenceTiles() + fillTilesNans();
			)
	{
		final int         margin = 0; // 1; // extra margins over 16x16 tiles to accommodate distorted destination tiles
		final TileProcessor tp =         reference_QuadClt.getTileProcessor();
		final int tilesX =               tp.getTilesX();
		final int tilesY =               tp.getTilesY();
		final int transform_size =       tp.getTileSize();
		final int macroTilesX =          tilesX/transform_size;
		final int macroTilesY =          tilesY/transform_size;
		final int macroX0 =              (tilesX - macroTilesX * transform_size)/2; // distribute extra tiles symmetrically ==0 for 324
		final int macroY0 =              (tilesY - macroTilesY * transform_size)/2; // distribute extra tiles symmetrically (242 - 1 tile above and below)
		final double [][] pXpYD = new double [macroTilesX*macroTilesY][];
		final int fullTileSize =         2 * (transform_size + margin);
		final int fullTileLen =          fullTileSize * fullTileSize; 
		final Thread[] threads = ImageDtt.newThreadArray(threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				public void run() {
					for (int iMTile = ai.getAndIncrement(); iMTile < pXpYD.length; iMTile = ai.getAndIncrement()) {
						if (reference_tiles[iMTile] != null) {
							
							int mtileY = iMTile / macroTilesX; 
							int mtileX = iMTile % macroTilesX;
							double pY = transform_size * (mtileY * transform_size + macroY0 + transform_size/2);
							double pX = transform_size * (mtileX * transform_size + macroX0 + transform_size/2);
							// find average disparity
							double sw = 0.0, swd = 0.0;
							for (int iTile = 0; iTile < fullTileLen; iTile++) {
								double disparity = reference_tiles[iMTile][QuadCLT.DSRBG_DISPARITY][iTile];
								double strength =  reference_tiles[iMTile][QuadCLT.DSRBG_STRENGTH][iTile];
								if (!Double.isNaN(disparity) && (strength> 0.0)) {
									sw  += strength;
									swd += strength * disparity;
								}
							}
							if (sw > 0.0) {
								pXpYD[iMTile] = new double[] {pX, pY, swd/sw};
							}
						}
					}
				}
			};
		}		      
		ImageDtt.startAndJoin(threads);
		return pXpYD;
	}
	/**
	 * Prepare scene tiles for correlation with the reference ones. Tiles include 5 layers: disparity,
	 * strength and 3 average color components (red, blue and green). 
	 * @param scene_xyz scene X (right),Y (up), Z (negative away form camera) in the reference camera coordinates
	 *        or null to use scene instance coordinates.
	 * @param scene_atr Scene azimuth, tilt and roll (or null to use scene instance).
	 * @param scene_QuadClt Scene QuadCLT instance.
	 * @param reference_QuadClt Reference QuadCLT instance.
	 * @param reference_tiles Reference tiles prepared for correlation.
	 * @param flowXY0 Initial offset of scene tiles (in image pixels) in x (right) and y (down) directions or null (to use all zeros).
	 * @param flowXY_frac Per macrotile residual [-0.5, 0.5) offsets that will be calculated. Should be initialized to
	 *        new double [number_of_macrotiles][] by the caller.
	 * @param tolerance_absolute Filter tiles by having disparity close to the average disparity of the corresponding reference
	 *        ones. tolerance_absolute is the absolute part of the disparity difference. 
	 * @param tolerance_relative Additional component of the disparity tolerance proportional to the reference macrotile disparity.
	 * @param occupancy Skip scene macrotiles having less remaining tiles fraction of all tiles.
	 * @param num_passes Number of Laplacian passes to fill filtered by disparity tiles. 
	 * @param max_change Break the loop of Laplassian passes if the maximal absolute value of the last pass changes falls below
	 *        this threshold.
	 * @param debug_level Debug level.
	 * @return Scene macrotiles - double array [number_of_macrotiles][number_of_channels][numer_of_tiles_per_macrotile], typically
	 *         [][5][256]
	 */
	public double [][][] prepareSceneTiles(// to match to reference
			// null for {scene,reference}{xyz,atr} uses instances globals 
			final double []   scene_xyz,     // camera center in world coordinates
			final double []   scene_atr,     // camera orientation relative to world frame
			final QuadCLT     scene_QuadClt,
			final QuadCLT     reference_QuadClt,
			final double [][][] reference_tiles, // prepared with prepareReferenceTiles() + fillTilesNans();
			final double [][] flowXY, // per macro tile {mismatch in image pixels in X and Y directions
			final double [][] flowXY_frac, // should be initialized as [number of macro tiles][] - returns fractional tile shifts [-0.5, 0.5)
			final double      tolerance_absolute, // absolute disparity half-range in each tile
			final double      tolerance_relative, // relative disparity half-range in each tile
			final double      occupancy,          // fraction of remaining  tiles (<1.0)
			final int         num_passes,
			final double      max_change,
			final int         debug_level)
	{
		final double   diagonal_weight = 0.5 * Math.sqrt(2.0); // relative to ortho
		final int         margin = 0; // 1; // extra margins over 16x16 tiles to accommodate distorted destination tiles
		final TileProcessor tp =         reference_QuadClt.getTileProcessor();
		final double [][] dsrbg_scene =  scene_QuadClt.getDSRBG();
		final int tilesX =               tp.getTilesX();
		final int tilesY =               tp.getTilesY();
		final int transform_size =       tp.getTileSize();
		final int macroTilesX =          tilesX/transform_size;
		final int macroTilesY =          tilesY/transform_size;
		final int macroX0 =              (tilesX - macroTilesX * transform_size)/2; // distribute extra tiles symmetrically ==0 for 324
		final int macroY0 =              (tilesY - macroTilesY * transform_size)/2; // distribute extra tiles symmetrically (242 - 1 tile above and below)
		final double [][][] scene_tiles = new double [macroTilesX*macroTilesY][][];
		final int fullTileSize =         2 * (transform_size + margin);
		final int fullTileLen =          fullTileSize * fullTileSize; 
		final Thread[] threads = ImageDtt.newThreadArray(threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		final double []  hist_weights = new double [fullTileSize];
		for (int i = 0; i < transform_size; i++) {
			hist_weights[margin + transform_size/2+ i] = Math.sin(Math.PI * (i +0.5) / transform_size);
		}
		final ErsCorrection ers_reference = reference_QuadClt.getErsCorrection();
		final ErsCorrection ers_scene =     scene_QuadClt.getErsCorrection();
		ers_reference.setupERS(); // just in case - setUP using instance paRAMETERS
		ers_scene.setupERS();     // just in case - setUP using instance paRAMETERS
		final int hist_len = 8;
		double wdiag = 0.25 *diagonal_weight / (diagonal_weight + 1.0);
		double wortho = 0.25 / (diagonal_weight + 1.0);
		final double [] neibw = {wortho, wdiag, wortho, wdiag, wortho, wdiag, wortho, wdiag}; 
		final int dbg_mtile =  (debug_level > 0)? 54: -1; // 54 : -1; // 453; // 500;//250;
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				public void run() {
					//					double [][] pXpYD = new double [fullTileLen][];
					double [][] tXtYD = new double [fullTileLen][]; // measured in tiles, not pixels (disparity - still pixels)
					for (int iMTile = ai.getAndIncrement(); iMTile < reference_tiles.length; iMTile = ai.getAndIncrement())
						if ((reference_tiles[iMTile] != null) && (flowXY[iMTile] != null)){
							if (iMTile == dbg_mtile) {
								System.out.println("prepareSceneTiles(): iMTile = "+iMTile);
							}
							int mtileY = iMTile / macroTilesX; 
							int mtileX = iMTile % macroTilesX;
							int tY0 = mtileY * transform_size + macroY0 -transform_size/2 - margin;
							int tX0 = mtileX * transform_size + macroX0 -transform_size/2 - margin;
							//						Arrays.fill(pXpYD, null);
							Arrays.fill(tXtYD, null);
							for (int iY = 0; iY < fullTileSize; iY++) {
								int tileY = tY0 + iY;
								if ((tileY >= 0) && (tileY < tilesY)) {
									for (int iX = 0; iX < fullTileSize; iX++) {
										int tileX = tX0 + iX;
										if ((tileX >= 0) && (tileX < tilesX)) {
											//										int nTile = tileX + tileY * tilesX;
											int iTile = iX + iY * fullTileSize;

											double disparity = reference_tiles[iMTile][QuadCLT.DSRBG_DISPARITY][iTile];
											if (disparity < 0) {
												disparity = 0.0; // is it needed or should it work with negative too?
											}
											double centerX = tileX * transform_size + transform_size/2; //  - shiftX;
											double centerY = tileY * transform_size + transform_size/2; //  - shiftY;
											if (!Double.isNaN(disparity)) {
												double [] xyd = ers_reference.getImageCoordinatesERS(
														scene_QuadClt, // QuadCLT cameraQuadCLT, // camera station that got image to be to be matched 
														centerX,        // double px,                // pixel coordinate X in this camera view
														centerY,        //double py,                // pixel coordinate Y in this camera view
														disparity,      // double disparity,         // this view disparity 
														true,           // boolean distortedView,    // This camera view is distorted (diff.rect), false - rectilinear
														ZERO3,          // double [] reference_xyz,  // this view position in world coordinates (typically ZERO3)
														ZERO3,          // double [] reference_atr,  // this view orientation relative to world frame  (typically ZERO3)
														true,           // boolean distortedCamera,  // camera view is distorted (false - rectilinear)
														scene_xyz,     // double [] camera_xyz,     // camera center in world coordinates
														scene_atr,     // double [] camera_atr,     // camera orientation relative to world frame
														LINE_ERR);       // double    LINE_ERR)       // threshold error in scan lines (1.0)
												// xyd[0], xyd[1] are here offset by transform_size/2 (to the center of the tile)
												if (xyd != null) {
													tXtYD[iTile] = new double [] {
															((xyd[0] + flowXY[iMTile][0])/transform_size) - 0.5, // moving 0.5 tiles back
															((xyd[1] + flowXY[iMTile][1])/transform_size) - 0.5, // moving 0.5 tiles back
															xyd[2]};
												} else {
													tXtYD[iTile] = null;
												}
											} else {
												tXtYD[iTile] = null;
											}
										}
									}
								}
							}
							if ((debug_level>1) && (iMTile == dbg_mtile)) {
								String [] dbg_titles= {"tX","tY","disp"};
								String dbg_title= "tXtYD-MX"+mtileX+"_MY"+mtileY;
								double [][] dbg_img = new double [3][fullTileLen];
								for (int nt =0; nt < fullTileLen; nt++) {
									if(tXtYD[nt] != null) {
										for (int i=0; i < dbg_img.length; i++) dbg_img[i][nt] = tXtYD[nt][i];
									} else {
										for (int i=0; i < dbg_img.length; i++) dbg_img[i][nt] = Double.NaN;
									}
								}
								(new ShowDoubleFloatArrays()).showArrays(
										dbg_img,
										fullTileSize,
										fullTileSize,
										true,
										dbg_title,
										dbg_titles);
							}

							// Find best fractional pixel offset
							double [] hist_fx = new double [hist_len];
							double [] hist_fy = new double [hist_len];
							for (int iY = 0; iY < fullTileSize; iY++) if (hist_weights[iY] > 0.0){
								for (int iX = 0; iX < fullTileSize; iX++) if (hist_weights[iX] > 0.0){
									int iTile = iX + iY * fullTileSize;
									if (tXtYD[iTile] != null) {
										int hidx = (int) Math.floor((tXtYD[iTile][0] - Math.floor(tXtYD[iTile][0])) * hist_len);
										if (hidx < 0) hidx = 0; 
										else if (hidx >= hist_len) hidx = hist_len - 1;
										hist_fx[hidx] += hist_weights[iX];
										int hidy = (int) Math.floor((tXtYD[iTile][1] - Math.floor(tXtYD[iTile][1])) * hist_len);
										if (hidy < 0) hidy = 0; 
										else if (hidy >= hist_len) hidy = hist_len - 1;
										hist_fy[hidy] += hist_weights[iY];
									}
								}
							}
							int hist_fx_mx = 0;
							int hist_fy_mx = 0;
							for (int i = 1; i < hist_len; i++) {
								if (hist_fx[i] > hist_fx[hist_fx_mx]) hist_fx_mx = i;
								if (hist_fy[i] > hist_fy[hist_fy_mx]) hist_fy_mx = i;
							}
							double offsX = (0.5 + hist_fx_mx) / hist_len;
							double offsY = (0.5 + hist_fy_mx) / hist_len;
							double swx = 0.0, swy = 0.0, sx =  0.0, sy = 0.0;
							for (int iY = 0; iY < fullTileSize; iY++) {
								for (int iX = 0; iX < fullTileSize; iX++) {
									int iTile = iX + iY * fullTileSize;
									if (tXtYD[iTile] != null) {
										double wx = hist_weights[iX];
										double wy = hist_weights[iX];
										double dx = tXtYD[iTile][0] - offsX;
										double dy = tXtYD[iTile][1] - offsY;
										dx = dx - Math.round(dx);
										dy = dy - Math.round(dy);
										swx += wx;
										swy += wy;
										sx += wx * dx;
										sy += wy * dy;
									}
								}
							}
							offsX += sx/swx;
							offsY += sy/swy;
							// use offsX, offsY as fractional shift and for data interpolation
							if (offsX >= .5) offsX -= 1.0;
							if (offsY >= .5) offsY -= 1.0;
//							flowXY_frac[iMTile] = new double [] {offsX, offsY};
							flowXY_frac[iMTile] = new double [] {-offsX, -offsY};
							double min_tX = Double.NaN, max_tX = Double.NaN, min_tY = Double.NaN, max_tY = Double.NaN;
							for (int iY = 0; iY < fullTileSize; iY++) {
								for (int iX = 0; iX < fullTileSize; iX++) {
									int iTile = iX + iY * fullTileSize;
									if (tXtYD[iTile] != null) {
										tXtYD[iTile][0]-=offsX;
										tXtYD[iTile][1]-=offsY;
										if (Double.isNaN(min_tX)) {
											min_tX =tXtYD[iTile][0]; 
											min_tY =tXtYD[iTile][1];
											max_tX = min_tX;
											max_tY = min_tY;
										}
										if (min_tX > tXtYD[iTile][0]) min_tX = tXtYD[iTile][0];
										if (min_tY > tXtYD[iTile][1]) min_tY = tXtYD[iTile][1];
										if (max_tX < tXtYD[iTile][0]) max_tX = tXtYD[iTile][0];
										if (max_tY < tXtYD[iTile][1]) max_tY = tXtYD[iTile][1];
									}
								}
							}
							int imin_tX = (int) Math.floor(min_tX);
							int imin_tY = (int) Math.floor(min_tY);
							int imax_tX = (int) Math.ceil (max_tX);
							int imax_tY = (int) Math.ceil (max_tY);
							// See if at least some of fits into the frame
							if ((imin_tX >= tilesX) || (imin_tY >= tilesY) || (imax_tX < 0)  || (imax_tX < 0)) {
								continue; // no overlap at all
							}
							int iwidth =  imax_tX - imin_tX + 1;
							int iheight = imax_tY - imin_tY + 1;
////							if ((iwidth <= 0) || (iheight <= 0)) {
							if ((iwidth <= 1) || (iheight <= 1)) {
								System.out.println ("prepareSceneTiles(): iwidth ="+iwidth+", iheight ="+iheight+", min_tX="+min_tX+", imin_tY="+imin_tY+", max_tX="+max_tX+", imax_tY="+imax_tY);
								continue;
							}
							double [][] scene_slices = new double [dsrbg_scene.length][iwidth*iheight]; //OOM here
							for (int iY = 0; iY < iheight; iY++) {
								int tY = imin_tY + iY;
								if ((tY >= 0) && (tY < tilesY)) {
									for (int iX = 0; iX < iwidth; iX++) {
										int tX = imin_tX + iX;
										if ((tX >= 0) && (tX < tilesX)) {
											int iTile = iX + iY * iwidth;
											int tile = tX + tilesX * tY;
											if (dsrbg_scene[QuadCLT.DSRBG_STRENGTH][tile] > 0.0) {
												double d = dsrbg_scene[QuadCLT.DSRBG_DISPARITY][tile];
												scene_slices[QuadCLT.DSRBG_DISPARITY][iTile] = d;
												if (!Double.isNaN(d)) {
													scene_slices[QuadCLT.DSRBG_STRENGTH][iTile] = dsrbg_scene[QuadCLT.DSRBG_STRENGTH][tile];
												}
											}
										}
									}
								}
							}
							if ((debug_level>1) && (iMTile == dbg_mtile)) {
								String [] dbg_titles= scene_QuadClt.getDSRGGTitles(); //{"d","s","r","b","g"};
								String dbg_title= "before_disparity-MX"+mtileX+"_MY"+mtileY;
								(new ShowDoubleFloatArrays()).showArrays(
										scene_slices,
										iwidth,
										iheight,
										true,
										dbg_title,
										dbg_titles);
							}
							// filter by disparity - 1 tile around rounded
							final TileNeibs tn =  new TileNeibs(iwidth, iheight);
							for (int iY = 0; iY < fullTileSize; iY++) {
								for (int iX = 0; iX < fullTileSize; iX++) {
									int iTile = iX + iY * fullTileSize;
									if (tXtYD[iTile] != null) {
										double disp_tolerance = tolerance_absolute + tXtYD[iTile][2] * tolerance_relative;
										double disp_min = tXtYD[iTile][2] - disp_tolerance;
										double disp_max = tXtYD[iTile][2] + disp_tolerance;
										int nt = tn.getIndex(
												(int) Math.round(tXtYD[iTile][0]) - imin_tX,
												(int) Math.round(tXtYD[iTile][1]) - imin_tY);
										if (nt >= 0) { // should always be
											for (int dir = 0; dir <9; dir++) {
												int nt1 = tn.getNeibIndex(nt, dir);
												if ((nt1 >= 0) && (scene_slices[QuadCLT.DSRBG_STRENGTH][nt1] > 0.0)) {
													if (    (scene_slices[QuadCLT.DSRBG_DISPARITY][nt1] < disp_min) ||
															(scene_slices[QuadCLT.DSRBG_DISPARITY][nt1] > disp_max)) {
														scene_slices[QuadCLT.DSRBG_STRENGTH][nt1] = 0.0; // disable tile
													}
												}
											}
										}
									}
								}
							}
							if ((debug_level>1) && (iMTile == dbg_mtile)) {
								String [] dbg_titles= scene_QuadClt.getDSRGGTitles(); //{"d","s","r","b","g"};
								String dbg_title= "after_disparity-MX"+mtileX+"_MY"+mtileY;
								(new ShowDoubleFloatArrays()).showArrays(
										scene_slices,
										iwidth,
										iheight,
										true,
										dbg_title,
										dbg_titles);
							}

							// copy rest of the data (all but strength) where strength > 0
							for (int iY = 0; iY < iheight; iY++) {
								int tY = imin_tY + iY;
								if ((tY >= 0) && (tY < tilesY)) {
									for (int iX = 0; iX < iwidth; iX++) {
										int tX = imin_tX + iX;
										if ((tX >= 0) && (tX < tilesX)) {
											int iTile = iX + iY * iwidth;
											int tile = tX + tilesX * tY;
											if (scene_slices[QuadCLT.DSRBG_STRENGTH][iTile] > 0.0) {
												for (int i = 0; i < scene_slices.length; i++) if (i != QuadCLT.DSRBG_STRENGTH) {
													scene_slices[i][iTile] = dsrbg_scene[i][tile];
												}
											}
										}
									}
								}
							}
							if ((debug_level>1) && (iMTile == dbg_mtile)) {
								String [] dbg_titles=scene_QuadClt.getDSRGGTitles(); //{"d","s","r","b","g"};
								String dbg_title= "scene1-MX"+mtileX+"_MY"+mtileY;
								(new ShowDoubleFloatArrays()).showArrays(
										scene_slices,
										iwidth,
										iheight,
										true,
										dbg_title,
										dbg_titles);
							}

							// set NAN everywhere where strength is 0 (including border tiles
							int num_dead = 0;
							for (int iTile = 0; iTile < scene_slices[0].length; iTile++) {
								if (scene_slices[QuadCLT.DSRBG_STRENGTH][iTile] <=0.0) {
									for (int i = 0; i < scene_slices.length; i++) {
										scene_slices[i][iTile] = Double.NaN;
									}
									num_dead++;
								}
							}
							if ((debug_level>1) && (iMTile == dbg_mtile)) {
								System.out.println("scene2-MX"+mtileX+"_MY"+mtileY+" num_dead="+num_dead);
								String [] dbg_titles= scene_QuadClt.getDSRGGTitles(); //{"d","s","r","b","g"};
								String dbg_title= "scene2-MX"+mtileX+"_MY"+mtileY;
								(new ShowDoubleFloatArrays()).showArrays(
										scene_slices,
										iwidth,
										iheight,
										true,
										dbg_title,
										dbg_titles);
							}


							//center_occupancy -> all occupancy						
							double fract_active = 1.0*(scene_slices[0].length - num_dead)/scene_slices[0].length; // all tiles, not center
							if (fract_active < occupancy) {
								flowXY_frac[iMTile] = null;
								continue;
							}
							// Here need to fill NaNs, then

							tilesFillNaN(
									neibw,        // final double []   neibw,
									scene_slices, // final double [][] slices,
									num_passes,   // final int         num_passes,
									max_change,   // final double      max_change,
									iwidth);      // final int         width // got zero!

							if ((debug_level>1) && (iMTile == dbg_mtile)) {
								String [] dbg_titles= scene_QuadClt.getDSRGGTitles(); //{"d","s","r","b","g"};
								String dbg_title= "scene3NaN-MX"+mtileX+"_MY"+mtileY;
								(new ShowDoubleFloatArrays()).showArrays(
										scene_slices,
										iwidth,
										iheight,
										true,
										dbg_title,
										dbg_titles);
							}

							// bi-linear interpolate						
							double [][] scene_mapped = new double [dsrbg_scene.length][fullTileLen];
							boolean need_nan_filter = false;;
							for (int iY = 0; iY < fullTileSize; iY++) {
								for (int iX = 0; iX < fullTileSize; iX++) {
									int iTile = iX + iY * fullTileSize;
									if (tXtYD[iTile] != null) {
										int itX = (int) Math.floor(tXtYD[iTile][0]);
										int itY = (int) Math.floor(tXtYD[iTile][1]);
										double kX = tXtYD[iTile][0]-itX;
										double kY = tXtYD[iTile][1]-itY;
										itX -= imin_tX; // relative to scene_slices = new double [dsrbg_scene.length][iwidth*iheight];
										itY -= imin_tY; // relative to scene_slices = new double [dsrbg_scene.length][iwidth*iheight];
										int indx = itX + itY * iwidth; 
										for (int i = 0; i < scene_mapped.length; i++) {
											scene_mapped[i][iTile] = 
													(1.0 - kY) * ((1.0-kX) * scene_slices[i][indx] +          kX * scene_slices[i][indx +          1])+
													(      kY) * ((1.0-kX) * scene_slices[i][indx + iwidth] + kX * scene_slices[i][indx + iwidth + 1]);
										}									
									} else {
										for (int i = 0; i < scene_mapped.length; i++) {
											scene_mapped[i][iTile] = Double.NaN; // will need to filter?
										}
										need_nan_filter=true;
									}
								}
							}
							if ((debug_level>1) && (iMTile == dbg_mtile)) {
								String [] dbg_titles= scene_QuadClt.getDSRGGTitles(); //{"d","s","r","b","g"};
								String dbg_title= "mapped-MX"+mtileX+"_MY"+mtileY;
								(new ShowDoubleFloatArrays()).showArrays(
										scene_mapped,
										fullTileSize,
										fullTileSize,
										true,
										dbg_title,
										dbg_titles);
							}

							if (need_nan_filter) {
								tilesFillNaN(
										neibw,         // final double []   neibw,
										scene_mapped,  // final double [][] slices,
										num_passes,    // final int         num_passes,
										max_change,    // final double      max_change,
										fullTileSize); // final int         width
								if ((debug_level>1) && (iMTile == dbg_mtile)) {
									String [] dbg_titles= scene_QuadClt.getDSRGGTitles(); // {"d","s","r","b","g"};
									String dbg_title= "mappedNaN-MX"+mtileX+"_MY"+mtileY;
									(new ShowDoubleFloatArrays()).showArrays(
											scene_mapped,
											fullTileSize,
											fullTileSize,
											true,
											dbg_title,
											dbg_titles);
								}

							}
							scene_tiles[iMTile] = scene_mapped;
						}
				}
			};
		}		      
		ImageDtt.startAndJoin(threads);
		
		if (debug_level > 0) {
			
			
			// show debug image
			String title =  reference_QuadClt.getImageName() + "-" + scene_QuadClt.getImageName()+"ref-scene";
			/*
			showMacroTiles(
					title,        // String title,
					scene_tiles, // double [][][] source_tiles,
					scene_QuadClt,        // final QuadCLT qthis,
					margin);      // final int     margin); // extra margins over 16x16 tiles to accommodate distorted destination tiles
			*/
			showCompareMacroTiles(
					title,        // String title,
					new double [][][][] {reference_tiles, scene_tiles}, // double [][][][] source_tiles_sets,
					scene_QuadClt,        // final QuadCLT qthis,
					margin);      // final int     margin); // extra margins over 16x16 tiles to accommodate distorted destination tiles
			
			String [] dbg_titles= {"dX","dY"};
			String dbg_title= "flowXY_frac"; // TODO: Visualize RMS of individual tiles fitting
			double [][] dbg_img = new double [2][macroTilesX*macroTilesY];
			for (int nt =0; nt < flowXY_frac.length; nt++) {
				if(flowXY_frac[nt] != null) {
					for (int i=0; i < dbg_img.length; i++) dbg_img[i][nt] = flowXY_frac[nt][i];
				} else {
					for (int i=0; i < dbg_img.length; i++) dbg_img[i][nt] = Double.NaN;
				}
			}
			if (debug_level > 1+0) {
			(new ShowDoubleFloatArrays()).showArrays(
					dbg_img,
					macroTilesX,
					macroTilesY,
					true,
					dbg_title,
					dbg_titles);
			}
			
			
		}
//		System.out.println("fillTilesNans() DONE.");
		return scene_tiles;
		
	}
	
	/**
	 * Prepare reference tiles for correlation with the scene ones. Tiles include 5 layers: disparity,
	 * strength and 3 average color components (red, blue and green). 
	 * @param qthis Reference scene QuadCLT instance.
	 * @param tolerance_absolute Filter reference macrotiles by same disparity (within a disparity range) consisting of the sum 
	 *        of absolute disparity (tolerance_absolute) and a proportional to the average disparity (tolerance_relative).  
	 * @param tolerance_relative Relative to the average disparity part of the disparity filtering.
	 * @param center_occupancy Fraction of non-null tiles in the center 8x8 area of the reference macrotiles after disparity
	 *        filtering (see tolerance_absolute,  tolerance_relative). Below this threshold - skip that macrotile.
	 * @param debug_level Debug level.
	 * @return Reference macrotiles - double array [number_of_macrotiles][number_of_channels][numer_of_tiles_per_macrotile], typically
	 *         [][5][256]
	 */
	public double [][][] prepareReferenceTiles(
			final QuadCLT     qthis,
			final double      tolerance_absolute, // absolute disparity half-range in each tile
			final double      tolerance_relative, // relative disparity half-range in each tile
			final double      center_occupancy,   // fraction of remaining  tiles in the center 8x8 area (<1.0)
			final int         debug_level)
	{
		final int margin =               0; // 1; // extra margins over 16x16 tiles to accommodate distorted destination tiles
		final TileProcessor tp =         qthis.getTileProcessor();
		final double [][] dsrbg =        qthis.getDSRBG();
		final int tilesX =               tp.getTilesX();
		final int tilesY =               tp.getTilesY();
		final int transform_size =       tp.getTileSize();
		final int macroTilesX =          tilesX/transform_size;
		final int macroTilesY =          tilesY/transform_size;
		final int macroX0 =              (tilesX - macroTilesX * transform_size)/2; // distribute extra tiles symmetrically ==0 for 324
		final int macroY0 =              (tilesY - macroTilesY * transform_size)/2; // distribute extra tiles symmetrically (242 - 1 tile above and below)
		final double [][][] source_tiles = new double [macroTilesX*macroTilesY][][];
		final int fullTileSize =         2 * (transform_size + margin);
		final int fullTileLen =          fullTileSize * fullTileSize; 
		final int min_remain_center = (int) Math.round(center_occupancy * transform_size * transform_size);
		final Thread[] threads = ImageDtt.newThreadArray(threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		// indices of the center 8x8 of the full 20x20 tile (assuming margin = 4) 
		final Integer [] order_indices = new Integer [transform_size*transform_size];
		for (int i = 0; i <transform_size; i++) {
			int i1 = i + margin + transform_size/2;
			for (int j = 0; j <transform_size; j++) {
				int j1 = j + margin + transform_size/2;
				order_indices[i * transform_size + j] = i1 * fullTileSize + j1;
			}
		}
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				public void run() {
					final double [] disparity = new double [fullTileLen];
					final double [] strength =  new double [fullTileLen];
					Integer [] local_indices = new Integer [transform_size*transform_size];
					for (int iMTile = ai.getAndIncrement(); iMTile < source_tiles.length; iMTile = ai.getAndIncrement()) {
						int mtileY = iMTile / macroTilesX; 
						int mtileX = iMTile % macroTilesX;
						int tY0 = mtileY * transform_size + macroY0 -transform_size/2 - margin;
						int tX0 = mtileX * transform_size + macroX0 -transform_size/2 - margin;
						Arrays.fill(strength,  0.0);
						Arrays.fill(disparity, 0.0); // Double.NaN);
						System.arraycopy(order_indices,0,local_indices,0,local_indices.length);
						
						for (int iY = 0; iY < fullTileSize; iY++) {
							int tileY = tY0 + iY;
							if ((tileY >= 0) && (tileY < tilesY)) {
								for (int iX = 0; iX < fullTileSize; iX++) {
									int tileX = tX0 + iX;
									if ((tileX >= 0) && (tileX < tilesX)) {
										int nTile = tileX + tileY * tilesX;
										int iTile = iX + iY * fullTileSize;
										double d = dsrbg[QuadCLT.DSRBG_DISPARITY][nTile];
										double s = dsrbg[QuadCLT.DSRBG_STRENGTH][nTile];
										if (!Double.isNaN(d) && (s > 0.0)){
											disparity[iTile] = d;
											strength[iTile]  = s;
										}
									}
								}
							}
						}
						double sw =   0.0;
						double swd =  0.0;
						int num_remain = 0;
						for (int i = 0; i < local_indices.length; i++) {
							double d = disparity[local_indices[i]];
							double s = strength[local_indices[i]];
							sw += s;
							swd += s*d;
							num_remain++;
						}
						if (num_remain < min_remain_center) {
							continue; // already too few tiles
						}
						Arrays.sort(local_indices, new Comparator<Integer>() {
						    @Override
						    public int compare(Integer lhs, Integer rhs) {
						        // -1 - less than, 1 - greater than, 0 - equal, not inverted for ascending disparity
						        //										return lhs.disparity > rhs.disparity ? -1 : (lhs.disparity < rhs.disparity ) ? 1 : 0;
//						        return disparity[lhs] < disparity[rhs] ? -1 : (disparity[lhs] > disparity[rhs] ) ? 1 : 0;
						    	// modifying to make  NaN greater than all non-NaN and equal to each other
						        return (!(disparity[lhs] >= disparity[rhs])) ? -1 : (!(disparity[lhs] <= disparity[rhs] )) ? 1 : 0;
						    }
						});
						int indx_min = 0;
						int indx_max = num_remain - 1;

						double d_low = Double.NaN;
						double d_high = Double.NaN; 

						while (true) {
							double disp_avg = swd / sw;
							double d_tol = tolerance_absolute + Math.max(disp_avg * tolerance_relative, 0.0);
							d_low =  disp_avg - d_tol;
							d_high = disp_avg + d_tol;
							// see if both min and max are within tolerance
							double d_min = disp_avg - disparity[local_indices[indx_min]];
							double d_max = disparity[local_indices[indx_max]] - disp_avg; 
							
							if ((d_min <= d_tol) && (d_max <= d_tol)) {
								break; // already OK
							}
							num_remain --;
							if (num_remain < min_remain_center) {
								break;
							}
							int indx_gone = -1;
							if (d_min > d_max) {
								indx_gone = indx_min;
								indx_min++;
							} else {
								indx_gone = indx_max;
								indx_max--;
							}
							double d = disparity[local_indices[indx_gone]];
							double s = strength[local_indices[indx_gone]];
							sw  -= s;
							swd -= d * s;
						}
						if (num_remain < min_remain_center) {
							continue; // too few remains in this tile
						}
						// now put all tiles in fullTileSize * fullTileSize that fit in disparity tolerance and positive strength
						source_tiles[iMTile] = new double[dsrbg.length][fullTileLen];
						for (int l = 0; l < dsrbg.length; l++) {
							Arrays.fill(source_tiles[iMTile][l], Double.NaN);
						}
						for (int i = 0; i < fullTileLen; i++) {
							if (!Double.isNaN(disparity[i]) && (strength[i] > 0.0) && (disparity[i] >= d_low) && (disparity[i] <= d_high)) {
								int tileY = tY0 + i / fullTileSize;
								int tileX = tX0 + i % fullTileSize;
								if ((tileY >= 0) && (tileY < tilesY) && (tileX >= 0) && (tileX < tilesX)) {
									int tile = tileY * tilesX + tileX; 
									for (int l = 0; l < dsrbg.length; l++) {
										source_tiles[iMTile][l][i] = dsrbg[l][tile];
									}
								}
							}
						}
					}
				}
			};
		}		      
		ImageDtt.startAndJoin(threads);
		if (debug_level > 0) {
			// show debug image
			String title = qthis.getImageName()+"-F"+center_occupancy+"-A"+tolerance_absolute+"-R"+tolerance_relative;
			showMacroTiles(
					title,        // String title,
					source_tiles, // double [][][] source_tiles,
					qthis,        // final QuadCLT qthis,
					margin);      // final int     margin); // extra margins over 16x16 tiles to accommodate distorted destination tiles
			
		}

		return source_tiles;
	}
	
	/**
	 * Show macrotiles as an image stack.
	 * @param title Image title to use.
	 * @param macro_tiles macrotiles array as generated by prepareSceneTiles() or prepareReferenceTiles().
	 * @param qthis Scene instance to extract dimensions
	 * @param margin Extra margin around the tiles (not used currently, always 0)
	 */
	public void showMacroTiles(
			String title,
			double [][][] macro_tiles,
			final QuadCLT qthis,
			final int     margin) // extra margins over 16x16 tiles to accommodate distorted destination tiles
	{
		final TileProcessor tp =         qthis.getTileProcessor();
		final double [][] dsrbg =        qthis.getDSRBG();
		final int tilesX =               tp.getTilesX();
		final int tilesY =               tp.getTilesY();
		final int transform_size =       tp.getTileSize();
		final int macroTilesX =          tilesX/transform_size;
		final int macroTilesY =          tilesY/transform_size;
		final int fullTileSize =         2 * (transform_size + margin);
		
		// show debug image
		final int dbg_with =   macroTilesX * (fullTileSize +1) - 1;
		final int dbg_height = macroTilesY * (fullTileSize +1) - 1;
		final double [][] dbg_img = new double [dsrbg.length][dbg_with * dbg_height];
		for (int l = 0; l < dbg_img.length; l++) {
			Arrays.fill(dbg_img[l],  Double.NaN);
		}
		for (int mtile = 0; mtile < macro_tiles.length; mtile++) if (macro_tiles[mtile] != null){
			int mTileY = mtile / macroTilesX;
			int mTileX = mtile % macroTilesX;
			for (int iY = 0; iY < fullTileSize; iY++) {
				int tileY = (fullTileSize +1) * mTileY + iY;
				for (int iX = 0; iX < fullTileSize; iX++) {
					int tileX = (fullTileSize +1) * mTileX + iX;
					for (int l = 0; l < dbg_img.length; l++) {
						dbg_img[l][tileY * dbg_with + tileX] = macro_tiles[mtile][l][iY * fullTileSize + iX];
					}							
				}
			}
		}
		String [] dsrbg_titles = qthis.getDSRGGTitles(); //{"d", "s", "r", "b", "g"};
		(new ShowDoubleFloatArrays()).showArrays(
				dbg_img,
				dbg_with,
				dbg_height,
				true,
				title,
				dsrbg_titles);
		
	}

	
	
	/**
	 * Calculate and display comparison stack of reference and scene images 
	 * @param suffix Add this text to the end of image name
	 * blur_reference Process and blur the reference image same as the scene one 
	 * @param camera_xyz0 scene camera offset in world coordinates
	 * @param camera_atr0 scene camera orientation  in world coordinates
	 * @param reference_QuadCLT reference scene instance
	 * @param scene_QuadCLT scene instance
	 * @param iscale upsample for interpolation
	 */
	public void compareRefSceneTiles(
			String suffix,
			boolean blur_reference,
			double [] camera_xyz0,
			double [] camera_atr0,
			QuadCLT reference_QuadCLT,
			QuadCLT scene_QuadCLT,
			int iscale) // 8
	{
		String title =  reference_QuadCLT.getImageName()+"-"+scene_QuadCLT.image_name+suffix;
		double [][] dsrbg = transformCameraVew( // shifts previous image correctly (right)
				title,                   // final String    title,
				null, // final double [][] dsrbg_camera_in,
				camera_xyz0, // double [] camera_xyz, // camera center in world coordinates
				camera_atr0, //double [] camera_atr, // camera orientation relative to world frame
				scene_QuadCLT,       // QuadCLT   camera_QuadClt,
				reference_QuadCLT,       // reference
				iscale);
		double [][] dsrbg_ref;
		if (blur_reference) {
			dsrbg_ref = transformCameraVew( // shifts previous image correctly (right)
					title+"-reference",     // final String    title,
					null, // final double [][] dsrbg_camera_in,
					ZERO3, // camera_xyz0,  // double [] camera_xyz, // camera center in world coordinates
					ZERO3, // camera_atr0,  // double [] camera_atr, // camera orientation relative to world frame
					reference_QuadCLT,      // scene_QuadCLT,       // QuadCLT   camera_QuadClt,
					reference_QuadCLT,      // reference_QuadCLT,       // reference
					iscale);
		} else {
			dsrbg_ref= reference_QuadCLT.getDSRBG();
		}
		double [][][] pair = {dsrbg_ref,  dsrbg};
		TileProcessor tp = reference_QuadCLT.getTileProcessor();
        int tilesX = tp.getTilesX();
		int tilesY = tp.getTilesY();
		String [] dsrbg_titles = reference_QuadCLT.getDSRGGTitles(); // {"d", "s", "r", "b", "g"};

		// combine this scene with warped previous one
		String [] rtitles = new String[2* dsrbg_titles.length];
		double [][] dbg_rslt = new double [rtitles.length][];
		for (int i = 0; i < dsrbg_titles.length; i++) {
			rtitles[2*i] =    dsrbg_titles[i]+"0";
			rtitles[2*i+1] =  dsrbg_titles[i];
			dbg_rslt[2*i] =   pair[0][i];
			dbg_rslt[2*i+1] = pair[1][i];
		}
		(new ShowDoubleFloatArrays()).showArrays(
				dbg_rslt,
				tilesX,
				tilesY,
				true,
				title,
				rtitles);
		
	}
	
	public void compareRefSceneTiles(
			String suffix,
			boolean blur_reference,
			double [][][] scene_xyzatr, // includeS reference (last)
			double [][][] scene_ers_dt, // includeS reference (last)
			QuadCLT [] scenes,
			int iscale) // 8
	{
		int nscenes = scenes.length;
		int indx_ref = nscenes - 1;
		String title =  "previous_frames_matching"+suffix;
		double [][][] dsrbg = new double [nscenes][][];
		String [] time_stamps = new String[nscenes];
		// [0] - last scene before the reference one
		for (int i = 0; i < nscenes; i++) {
			int indx = dsrbg.length - i - 1;
			time_stamps[i] = scenes[indx].getImageName();
			System.out.println("\n"+i+ ": "+ time_stamps[i]+" compareRefSceneTiles()");
			if ((i == 0) && !blur_reference) {
				dsrbg[0]= scenes[indx_ref].getDSRBG();
			} else {
				ErsCorrection ers_scene = scenes[indx].getErsCorrection();
				double [] ers_scene_original_xyz_dt = ers_scene.getErsXYZ_dt();
				double [] ers_scene_original_atr_dt = ers_scene.getErsATR_dt();
				ers_scene.setErsDt(
						scene_ers_dt[indx][0], // double []    ers_xyz_dt,
						scene_ers_dt[indx][1]); // double []    ers_atr_dt)(ers_scene_original_xyz_dt);
				ers_scene.setupERS();
				dsrbg[i] = transformCameraVew(   // shifts previous image correctly (right) null pointer
						title,                   // final String    title,
						null, // final double [][] dsrbg_camera_in,
						scene_xyzatr[indx][0],   // double [] camera_xyz, // camera center in world coordinates
						scene_xyzatr[indx][1],   //double [] camera_atr, // camera orientation relative to world frame
						scenes[indx],            // QuadCLT   camera_QuadClt,
						scenes[indx_ref],        // reference
						iscale);
				ers_scene.setErsDt(
						ers_scene_original_xyz_dt, // double []    ers_xyz_dt,
						ers_scene_original_atr_dt); // double []    ers_atr_dt)(ers_scene_original_xyz_dt);
				ers_scene.setupERS();
			}
		}
		
		TileProcessor tp = scenes[indx_ref].getTileProcessor();
        int tilesX = tp.getTilesX();
		int tilesY = tp.getTilesY();
		String [] dsrbg_titles =  scenes[indx_ref].getDSRGGTitles(); // {"d", "s", "r", "b", "g"};{"d", "s", "r", "b", "g"};
		int nslices = dsrbg_titles.length;

		// combine this scene with warped previous one
		String [] rtitles = new String[nscenes * nslices];
		double [][] dbg_rslt = new double [rtitles.length][];
		for (int nslice = 0; nslice < nslices; nslice++) {
			for (int nscene = 0; nscene < nscenes; nscene++) {
				rtitles[nscenes * nslice + nscene] =    dsrbg_titles[nslice]+"-"+time_stamps[nscene];
				if (dsrbg[nscene] !=  null) {
					dbg_rslt[nscenes * nslice + nscene] =   dsrbg[nscene][nslice];
				}
			}
		}
		(new ShowDoubleFloatArrays()).showArrays(
				dbg_rslt,
				tilesX,
				tilesY,
				true,
				title,
				rtitles);

	}
	
	
	
	/**
	 * Show macrotiles in comparison, typically reference to scene ones
	 * @param title Image stack title. 
	 * @param source_tiles_sets typically {reference_tiles, scene tiles} pair, generated by prepareReferenceTiles()
	 *        and prepareSceneTiles(), respectively.
	 * @param qthis Scene instance to extract dimensions
	 * @param margin Extra margin around the tiles (not used currently, always 0)
	 */
	public void showCompareMacroTiles(
			String title,
			double [][][][] source_tiles_sets,
			final QuadCLT qthis,
			final int     margin) // extra margins over 16x16 tiles to accommodate distorted destination tiles
	{
		final TileProcessor tp =         qthis.getTileProcessor();
		final double [][] dsrbg =        qthis.getDSRBG();
		final int tilesX =               tp.getTilesX();
		final int tilesY =               tp.getTilesY();
		final int transform_size =       tp.getTileSize();
		final int macroTilesX =          tilesX/transform_size;
		final int macroTilesY =          tilesY/transform_size;
		final int fullTileSize =         2 * (transform_size + margin);
		
		// show debug image
		final int dbg_with =   macroTilesX * (fullTileSize +1) - 1;
		final int dbg_height = macroTilesY * (fullTileSize +1) - 1;
		final double [][] dbg_img = new double [dsrbg.length * source_tiles_sets.length][dbg_with * dbg_height];
		for (int l = 0; l < dbg_img.length; l++) {
			Arrays.fill(dbg_img[l],  Double.NaN);
		}
		String [] titles = qthis.getDSRGGTitles(); //{"d", "s", "r", "b", "g"};
		String [] dsrbg_titles = new String [titles.length * source_tiles_sets.length ]; 
		
		for (int iset = 0; iset < source_tiles_sets.length; iset ++) {
			for (int l = 0; l < titles.length; l++) {
				dsrbg_titles[l * source_tiles_sets.length + iset] = titles[l]+"-"+iset;
			}
			double [][][] source_tiles = source_tiles_sets[iset];
			for (int mtile = 0; mtile < source_tiles.length; mtile++) if (source_tiles[mtile] != null){
				int mTileY = mtile / macroTilesX;
				int mTileX = mtile % macroTilesX;
				for (int iY = 0; iY < fullTileSize; iY++) {
					int tileY = (fullTileSize +1) * mTileY + iY;
					for (int iX = 0; iX < fullTileSize; iX++) {
						int tileX = (fullTileSize +1) * mTileX + iX;
						for (int l = 0; l < titles.length; l++) {
							dbg_img[l*source_tiles_sets.length + iset][tileY * dbg_with + tileX] = source_tiles[mtile][l][iY * fullTileSize + iX];
						}							
					}
				}
			}
		}
//		String [] dsrbg_titles = {"d", "s", "r", "b", "g"};
		(new ShowDoubleFloatArrays()).showArrays(
				dbg_img,
				dbg_with,
				dbg_height,
				true,
				title,
				dsrbg_titles);
		
	}
	
	/**
	 * Show multiple sets of 2D correlation tiles as image stack
	 * @param title Image stack title. 
	 * @param corr_tiles An array of several 2D correlation tiles sets, each containing per-macrotile
	 *        2D correlation tiles (or null), typically 255-long.
	 * @param tilesX Number of macrotiles in a row (typically 40 for 5MPix images).
	 * @param tile_width Correlation tile width (typically 15)
	 * @param tile_height  Correlation tile height (typically 15)
	 */
	public void showCorrTiles(
			String title,
			double [][][] corr_tiles,
			int           tilesX,
			int           tile_width,
			int           tile_height)
	{
		// show debug image
		final int dbg_with =   tilesX * (tile_width +1) - 1;
		final int dbg_height = (corr_tiles[0].length / tilesX) * (tile_height +1) - 1;
		final double [][] dbg_img = new double [corr_tiles.length][dbg_with * dbg_height];
		for (int l = 0; l < dbg_img.length; l++) {
			Arrays.fill(dbg_img[l],  Double.NaN);
		}
		for (int slice = 0; slice < corr_tiles.length; slice++) {
			for (int mtile = 0; mtile < corr_tiles[slice].length; mtile++) if (corr_tiles[slice][mtile] != null){
				int mTileY = mtile / tilesX;
				int mTileX = mtile % tilesX;
				for (int iY = 0; iY < tile_height; iY++) {
					int tileY = (tile_height +1) * mTileY + iY;
					for (int iX = 0; iX < tile_width; iX++) {
						int tileX = (tile_width +1) * mTileX + iX;
						dbg_img[slice][tileY * dbg_with + tileX] = corr_tiles[slice][mtile][iY * tile_width + iX];
					}
				}
			}
		}
		(new ShowDoubleFloatArrays()).showArrays(
				dbg_img,
				dbg_with,
				dbg_height,
				true,
				title); //	dsrbg_titles);
	}
	
	public double [][] getSceneDisparityStrength(
			final boolean     to_ref_disparity, // false - return scene disparity, true - convert disparity back to the reference scene
			final double []   disparity_ref,   // invalid tiles - NaN in disparity
			final double []   disparity_scene, // invalid tiles - NaN in disparity (just for masking out invalid scene tiles)
			final double []   strength_scene,  // to calculate interpolated strength
			final double []   scene_xyz, // camera center in world coordinates
			final double []   scene_atr, // camera orientation relative to world frame
			final double []   scene_ers_xyz_dt, // camera ERS linear
			final double []   scene_ers_atr_dt, // camera ERS linear
			final QuadCLT     scene_QuadClt,
			final QuadCLT     reference_QuadClt,
			final int         margin,
			final double      tolerance_ref_absolute,
			final double      tolerance_ref_relative,
			final double      tolerance_scene_absolute, // of 4 bi-linear interpolation corners
			final double      tolerance_scene_relative)  // of 4 bi-linear interpolation corners
	{
		TileProcessor tp = reference_QuadClt.getTileProcessor();
		final int tilesX = tp.getTilesX();
		final int tilesY = tp.getTilesY();
		final int transform_size = tp.getTileSize();
		final int tiles = tilesX*tilesY;
		final double scene_disparity_cor = 0.0;
		final boolean [] valid_tiles = new boolean[tiles];
		scene_QuadClt.getErsCorrection().setErsDt(
				scene_ers_xyz_dt, // double []    ers_xyz_dt,
				scene_ers_atr_dt); // double []    ers_atr_dt)(ers_scene_original_xyz_dt);
		//setupERS() will be inside transformToScenePxPyD()
		// OK to use the same reference_QuadClt for both reference_QuadClt and scene_QuadClt
		double [][] scene_pXpYD = transformToScenePxPyD( // will be null for disparity == NaN
				null, // final Rectangle [] extra_woi,    // show larger than sensor WOI (or null)
				disparity_ref,      // final double []   disparity_ref, // invalid tiles - NaN in disparity (maybe it should not be masked by margins?)
				scene_xyz,          // final double []   scene_xyz, // camera center in world coordinates
				scene_atr,          // final double []   scene_atr, // camera orientation relative to world frame
				scene_QuadClt,      // final QuadCLT     scene_QuadClt,
				reference_QuadClt); // final QuadCLT     reference_QuadClt)
		
		TpTask[]  tp_tasks =  GpuQuad.setInterTasks( // just to calculate valid_tiles
				scene_QuadClt.getNumSensors(),
				scene_QuadClt.getGeometryCorrection().getSensorWH()[0],
				!scene_QuadClt.hasGPU(), // final boolean             calcPortsCoordinatesAndDerivatives, // GPU can calculate them centreXY
				scene_pXpYD, // final double [][]         pXpYD, // per-tile array of pX,pY,disparity triplets (or nulls)
				null,          // final boolean []          selection, // may be null, if not null do not  process unselected tiles
				scene_QuadClt.getGeometryCorrection(), // final GeometryCorrection  geometryCorrection,
				scene_disparity_cor, // final double              disparity_corr,
    			margin, // final int                 margin,      // do not use tiles if their centers are closer to the edges
    			valid_tiles, // final boolean []          valid_tiles,            
    			threadsMax); // final int                 threadsMax)  // maximal number of threads to launch
		//FIXME:  not clear here tp_tasks was supposed to go? no
		/*
		scene_QuadClt.getGPU().setInterTasks(
				scene_pXpYD, // final double [][]         pXpYD, // per-tile array of pX,pY,disparity triplets (or nulls)
       			scene_QuadClt.getGeometryCorrection(), // final GeometryCorrection  geometryCorrection,
       			scene_disparity_cor, // final double              disparity_corr,
    			margin, // final int                 margin,      // do not use tiles if their centers are closer to the edges
    			valid_tiles, // final boolean []          valid_tiles,            
    			threadsMax); // final int                 threadsMax)  // maximal number of threads to launch
		*/
		
		final double [][] disparity_strength = new double [2][tiles];
		Arrays.fill(disparity_strength[0], Double.NaN);
		final Thread[] threads = ImageDtt.newThreadArray(threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				public void run() {
					double [] corners =   new double [4];
					double [] strengths = new double [4];
					for (int nTile = ai.getAndIncrement(); nTile < tiles; nTile = ai.getAndIncrement()) if (valid_tiles[nTile]){
						double pX =            scene_pXpYD[nTile][0];
						double pY =            scene_pXpYD[nTile][1];
						double tX = (pX - transform_size/2)/transform_size;
						double tY = (pY - transform_size/2)/transform_size;
						double disp_from_ref = scene_pXpYD[nTile][2];
						int itX0 = (int) Math.floor(tX); 
						int itY0 = (int) Math.floor(tY);
						double kX = tX - itX0;
						double kY = tY - itY0;
						double [] corner_weights = {(1.0-kX)*(1.0-kY), kX*(1.0-kY),(1.0-kX)*kY, kX*kY};
						int itX1 = itX0 + 1;
						if (itX1 >= tilesX) {
							itX1 = tilesX - 1;
						}
						int itY1 = itY0 + 1;
						if (itY1 >= tilesY) {
							itY1 = tilesY - 1;
						}
						corners[0] =  disparity_scene[itX0 + itY0 * tilesX];
						corners[1] =  disparity_scene[itX1 + itY0 * tilesX];
						corners[2] =  disparity_scene[itX0 + itY1 * tilesX];
						corners[3] =  disparity_scene[itX1 + itY1 * tilesX];
						strengths[0] = strength_scene[itX0 + itY0 * tilesX];
						strengths[1] = strength_scene[itX1 + itY0 * tilesX];
						strengths[2] = strength_scene[itX0 + itY1 * tilesX];
						strengths[3] = strength_scene[itX1 + itY1 * tilesX];

						double disp_tol_ref = tolerance_ref_absolute + disp_from_ref * tolerance_ref_relative;
						double disp_min_ref =  disp_from_ref - disp_tol_ref; 
						double disp_max_ref =  disp_from_ref + disp_tol_ref;
						int idisp_min = 0, idisp_max = 0;
						int num_corners = 0;
						for (int i = 0; i < 4; i++) {
							if (!Double.isNaN(corners[i])) {
								if ((corners[i] < disp_min_ref) || (corners[i] > disp_max_ref)){
									corners[i] = Double.NaN;
									continue;
								}
								if (!(corners[i] >= corners[idisp_min])) {
									idisp_min = i; 
								}
								if (!(corners[i] <= corners[idisp_max])) {
									idisp_max = i; 
								}
								num_corners++;
							}
						}
						if (num_corners > 0) {
							double disp_half = 0.5* (corners[idisp_min] + corners[idisp_max]);
							double disp_range = 2 * (tolerance_scene_absolute + tolerance_ref_relative * disp_half);
							while ((num_corners > 0) && ((corners[idisp_max] - corners[idisp_min]) > disp_range )) {
								num_corners--;
								if (num_corners > 0) {
									if (disp_half > disp_from_ref) { // remove max
										corners[idisp_max] = Double.NaN;
										// find new max
										idisp_max = idisp_min;
										for (int i = 0; i < 4; i++) {
											if (!(corners[i] <= corners[idisp_max])) {
												idisp_max = i; 
											}
										}
									} else {
										corners[idisp_min] = Double.NaN;
										// find new min
										idisp_min = idisp_max;
										for (int i = 0; i < 4; i++) {
											if (!(corners[i] >= corners[idisp_min])) {
												idisp_min = i; 
											}
										}
									}
								}
							}
							if (num_corners > 0) {
								double sw=0.0, swd = 0.0, sws = 0.0;
								for (int i = 0; i < 4; i++) {
									if (!Double.isNaN(corners[i])) {
										sw += corner_weights[i];
										swd += corner_weights[i] * corners[i];
										sws += corner_weights[i] * strengths[i];
									}
								}
								disparity_strength[0][nTile] = swd / sw; 
								disparity_strength[1][nTile] = sws / sw; 
							}
						}
					}
				}
			};
		}		      
		ImageDtt.startAndJoin(threads);
		if (to_ref_disparity) {
			ai.set(0);
//double [][] scene_pXpYD
			for (int ithread = 0; ithread < threads.length; ithread++) {
				threads[ithread] = new Thread() {
					public void run() {
						for (int nTile = ai.getAndIncrement(); nTile < tiles; nTile = ai.getAndIncrement()) if (scene_pXpYD != null) {
							if (Double.isNaN(disparity_strength[0][nTile])) {
								scene_pXpYD[nTile] = null;
							} else { // replace scene disparity (strength will stay the same
								scene_pXpYD[nTile][2] = disparity_strength[0][nTile];
							}
						}
					}
				};
			}		      
			ImageDtt.startAndJoin(threads);
			double [][] toref_pXpYD = transformFromScenePxPyD( // does not look at identity scene_xyz, scene_atr 
					scene_pXpYD, // final double [][] pXpYD_scene, // tiles correspond to reference, pX,pY,D - for scene
					scene_xyz,   // final double []   scene_xyz,   // camera center in world coordinates
					scene_atr,   // final double []   scene_atr,   // camera orientation relative to world frame
					scene_QuadClt, // final QuadCLT     scene_QuadClt,
					reference_QuadClt); // final QuadCLT     reference_QuadClt)
			ai.set(0);
			for (int ithread = 0; ithread < threads.length; ithread++) {
				threads[ithread] = new Thread() {
					public void run() {
						for (int nTile = ai.getAndIncrement(); nTile < tiles; nTile = ai.getAndIncrement()) {
							if (toref_pXpYD[nTile] != null) {
								disparity_strength[0][nTile] = toref_pXpYD[nTile][2];
							} else {
								disparity_strength[0][nTile] = Double.NaN;
							}
						}
					}
				};
			}		      
			ImageDtt.startAndJoin(threads);
		}
		return disparity_strength;
	}

	/**
	 * Transform disparity to world XYZ to test distance in meters
	 * @param disparity_ref [tilesX * tilesY] per-tile disparity
	 * @param quadClt       scene
	 * @param threadsMax    
	 * @return per tile array (with possible nulls) of X,Y,Z triplets in meters (left, up, negative distance)
	 */
	public static double [][] transformToWorldXYZ(
			final double []   disparity_ref, // invalid tiles - NaN in disparity
			final QuadCLT     quadClt, // now - may be null - for testing if scene is rotated ref
			int               threadsMax)
	{
		TileProcessor tp = quadClt.getTileProcessor();
		final int tilesX = tp.getTilesX();
		final int tilesY = tp.getTilesY();
		final int tiles = tilesX*tilesY;
		final int transform_size = tp.getTileSize();
		final double [][] world_xyz=           new double [tiles][];
		final ErsCorrection ersCorrection = quadClt.getErsCorrection();
		ersCorrection.setupERS();
		final Thread[] threads = ImageDtt.newThreadArray(threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				public void run() {
					for (int nTile = ai.getAndIncrement(); nTile < tiles; nTile = ai.getAndIncrement()) if (!Double.isNaN(disparity_ref[nTile])) {
						double disparity = disparity_ref[nTile];
						int tileY = nTile / tilesX;  
						int tileX = nTile % tilesX;
						double centerX = tileX * transform_size + transform_size/2; //  - shiftX;
						double centerY = tileY * transform_size + transform_size/2; //  - shiftY;
						world_xyz[nTile] = ersCorrection.getWorldCoordinatesERS( // ersCorrection - reference
								centerX,        // double px,                // pixel coordinate X in the reference view
								centerY,        // double py,                // pixel coordinate Y in the reference view
								disparity,      // double disparity,         // reference disparity 
								true,           // boolean distortedView,    // This camera view is distorted (diff.rect), false - rectilinear
								ZERO3,          // double [] reference_xyz,  // this view position in world coordinates (typically ZERO3)
								ZERO3);         // double [] reference_atr,  // this view orientation relative to world frame  (typically ZERO3)
					}
				}
			};
		}		      
		ImageDtt.startAndJoin(threads);
		return world_xyz;
	}

	public double [][] transformToScenePxPyD(
			final Rectangle   full_woi_in,      // show larger than sensor WOI (or null) IN TILES
			final double []   disparity_ref, // invalid tiles - NaN in disparity
			final double []   scene_xyz, // camera center in world coordinates
			final double []   scene_atr, // camera orientation relative to world frame
			final QuadCLT     scene_QuadClt,
			final QuadCLT     reference_QuadClt) {
		return transformToScenePxPyD(
				full_woi_in,   // final Rectangle [] extra_woi,    // show larger than sensor WOI (or null)
				disparity_ref, // invalid tiles - NaN in disparity
				scene_xyz,     // camera center in world coordinates
				scene_atr,     // camera orientation relative to world frame
				scene_QuadClt,
				reference_QuadClt,
				this.threadsMax);
	}
	
	/**
	 * Calculate pX, pY, Disparity triplets for the rotated scene to match uniform grid of a virtual camera
	 * Supports reference window larger that the physical sensor to show more of the other frames with partial
	 * overlap.
	 * @param full_woi_in null or a larger reference window {width, height, left, top}
	 * @param disparity_ref_in disparity value - either full_woi size or a reference frame only (rest will be 0)  
	 * @param scene_xyz scene linear offset (in meters)
	 * @param scene_atr scene azimuth, tilt, roll offset
	 * @param scene_QuadClt 
	 * @param reference_QuadClt
	 * @param threadsMax
	 * @return pX, pY, Disparity of the other scene. pX, pY are measured from the sensor top left corner
	 */
	public static double [][] transformToScenePxPyD(
			final Rectangle   full_woi_in,      // show larger than sensor WOI (or null) IN TILES
			final double []   disparity_ref_in, // invalid tiles - NaN in disparity
			final double []   scene_xyz,        // camera center in world coordinates
			final double []   scene_atr,        // camera orientation relative to world frame
			final QuadCLT     scene_QuadClt,
			final QuadCLT     reference_QuadClt, // now - may be null - for testing if scene is rotated ref
			int               threadsMax)
	{
		
		TileProcessor tp = scene_QuadClt.getTileProcessor();
		final int tilesX = (full_woi_in==null) ? tp.getTilesX() : full_woi_in.width; // full width,includeing extra
		final int tilesY = (full_woi_in==null) ? tp.getTilesY() : full_woi_in.height;
		final int offsetX_ref = (full_woi_in==null) ? 0 : full_woi_in.x;
		final int offsetY_ref = (full_woi_in==null) ? 0 : full_woi_in.y;
		int ref_w = tp.getTilesX();
		int ref_h = tp.getTilesY();
		double [] dref = disparity_ref_in;
		if (full_woi_in!=null) {
			if ((ref_w + offsetX_ref) > tilesX) ref_w =  tilesX - offsetX_ref;
			if ((ref_h + offsetY_ref) > tilesY) ref_h =  tilesY - offsetY_ref;
			if (disparity_ref_in.length < (full_woi_in.width * full_woi_in.height)) {
				dref= new double[full_woi_in.width * full_woi_in.height];
				for (int i = 0; i < ref_h; i++) {
					System.arraycopy(
							disparity_ref_in,
							i * tp.getTilesX(), // not truncated
							dref,
							(i + offsetY_ref) * full_woi_in.width + offsetX_ref,
							ref_w); // may be truncated
				}
			}
		}
		final boolean ref_is_identity = 
				(scene_xyz[0]==0.0) && (scene_xyz[1]==0.0) && (scene_xyz[2]==0.0) &&
				(scene_atr[0]==0.0) && (scene_atr[1]==0.0) && (scene_atr[2]==0.0);
		final double []   disparity_ref = dref;
//		final int tilesX_ref = ref_w;
//		final int tilesY_ref = ref_h;
		final int tiles = tilesX*tilesY;
		final int transform_size = tp.getTileSize();
		final double [][] pXpYD=               new double [tiles][];
		final ErsCorrection ersSceneCorrection =     scene_QuadClt.getErsCorrection();
		final ErsCorrection ersReferenceCorrection = (reference_QuadClt!=null)? reference_QuadClt.getErsCorrection(): ersSceneCorrection;
		if (reference_QuadClt!=null) {
			ersReferenceCorrection.setupERS(); // just in case - setUP using instance parameters
		}
		ersSceneCorrection.setupERS();
		final Thread[] threads = ImageDtt.newThreadArray(threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				public void run() {
					for (int nTile = ai.getAndIncrement(); nTile < tiles; nTile = ai.getAndIncrement()) if (!Double.isNaN(disparity_ref[nTile])) {
						double disparity = disparity_ref[nTile];
						int tileY = nTile / tilesX;  
						int tileX = nTile % tilesX;
						double centerX = (tileX + 0.5 - offsetX_ref) * transform_size; //  - shiftX;
						double centerY = (tileY + 0.5 - offsetY_ref) * transform_size; //  - shiftY;
						if (disparity < 0) {
							disparity = 1.0* disparity; // 0.0;
						}
						if ((scene_QuadClt == reference_QuadClt) && (ref_is_identity)) {
							pXpYD[nTile] = new double [] {centerX, centerY, disparity};
						} else {
							pXpYD[nTile] = ersReferenceCorrection.getImageCoordinatesERS( // ersCorrection - reference
									scene_QuadClt,  // QuadCLT cameraQuadCLT, // camera station that got image to be to be matched 
									centerX,        // double px,                // pixel coordinate X in the reference view
									centerY,        // double py,                // pixel coordinate Y in the reference view
									disparity,      // double disparity,         // reference disparity 
									true,           // boolean distortedView,    // This camera view is distorted (diff.rect), false - rectilinear
									ZERO3,          // double [] reference_xyz,  // this view position in world coordinates (typically ZERO3)
									ZERO3,          // double [] reference_atr,  // this view orientation relative to world frame  (typically ZERO3)
									true,           // boolean distortedCamera,  // camera view is distorted (false - rectilinear)
									scene_xyz,      // double [] camera_xyz,     // camera center in world coordinates
									scene_atr,      // double [] camera_atr,     // camera orientation relative to world frame
									LINE_ERR);      // double    line_err)       // threshold error in scan lines (1.0)
							if (pXpYD[nTile] != null) {
								if (    (pXpYD[nTile][0] < 0.0) ||
										(pXpYD[nTile][1] < 0.0) ||
										(pXpYD[nTile][0] > ersSceneCorrection.getSensorWH()[0]) ||
										(pXpYD[nTile][1] > ersSceneCorrection.getSensorWH()[1])) {
									pXpYD[nTile] = null;
								}
							}
						}
					}
				}
			};
		}		      
		ImageDtt.startAndJoin(threads);
		return pXpYD;
	}
	
	@Deprecated
	public static double [][] transformToScenePxPyD(
			final double []   disparity_ref, // invalid tiles - NaN in disparity
			final double []   scene_xyz, // camera center in world coordinates
			final double []   scene_atr, // camera orientation relative to world frame
			final QuadCLT     scene_QuadClt,
			final QuadCLT     reference_QuadClt, // now - may be null - for testing if scene is rotated ref
			int               threadsMax)
	{
		TileProcessor tp = scene_QuadClt.getTileProcessor();
		final int tilesX = tp.getTilesX();
		final int tilesY = tp.getTilesY();
		final int tiles = tilesX*tilesY;
		final int transform_size = tp.getTileSize();
		final double [][] pXpYD=               new double [tiles][];
		final ErsCorrection ersSceneCorrection =     scene_QuadClt.getErsCorrection();
		final ErsCorrection ersReferenceCorrection = (reference_QuadClt!=null)? reference_QuadClt.getErsCorrection(): ersSceneCorrection;
		if (reference_QuadClt!=null) {
			ersReferenceCorrection.setupERS(); // just in case - setUP using instance paRAMETERS
		}
		ersSceneCorrection.setupERS();
		final Thread[] threads = ImageDtt.newThreadArray(threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				public void run() {
					for (int nTile = ai.getAndIncrement(); nTile < tiles; nTile = ai.getAndIncrement()) if (!Double.isNaN(disparity_ref[nTile])) {
						double disparity = disparity_ref[nTile];
						int tileY = nTile / tilesX;  
						int tileX = nTile % tilesX;
						double centerX = tileX * transform_size + transform_size/2; //  - shiftX;
						double centerY = tileY * transform_size + transform_size/2; //  - shiftY;
						if (disparity < 0) {
							disparity = 1.0* disparity; // 0.0;
						}
						if (scene_QuadClt == reference_QuadClt) {
							pXpYD[nTile] = new double [] {centerX, centerY, disparity};
						} else {
							pXpYD[nTile] = ersReferenceCorrection.getImageCoordinatesERS( // ersCorrection - reference
									scene_QuadClt,  // QuadCLT cameraQuadCLT, // camera station that got image to be to be matched 
									centerX,        // double px,                // pixel coordinate X in the reference view
									centerY,        // double py,                // pixel coordinate Y in the reference view
									disparity,      // double disparity,         // reference disparity 
									true,           // boolean distortedView,    // This camera view is distorted (diff.rect), false - rectilinear
									ZERO3,          // double [] reference_xyz,  // this view position in world coordinates (typically ZERO3)
									ZERO3,          // double [] reference_atr,  // this view orientation relative to world frame  (typically ZERO3)
									true,           // boolean distortedCamera,  // camera view is distorted (false - rectilinear)
									scene_xyz,      // double [] camera_xyz,     // camera center in world coordinates
									scene_atr,      // double [] camera_atr,     // camera orientation relative to world frame
									LINE_ERR);      // double    line_err)       // threshold error in scan lines (1.0)
							if (pXpYD[nTile] != null) {
								if (    (pXpYD[nTile][0] < 0.0) ||
										(pXpYD[nTile][1] < 0.0) ||
										(pXpYD[nTile][0] > ersSceneCorrection.getSensorWH()[0]) ||
										(pXpYD[nTile][1] > ersSceneCorrection.getSensorWH()[1])) {
									pXpYD[nTile] = null;
								}
							}
						}
					}
				}
			};
		}		      
		ImageDtt.startAndJoin(threads);
		return pXpYD;
	}


	
	//TODO: refine inter-scene pose to accommodate refined disparity map
	/**
	 * Removing BG tiles that are not visible because of the FG ones
	 * @param tp TileProcessor instance to get image dimensions
	 * @param pXpYD Array of pX, pY, Disparity triplets for the current camera calculated from the reference 3D model
	 * @param max_overlap maximal area overlap (for the full 16x16 image tiles) that allows the BG tile to be kept
	 * @param pXpYD_cam optional array of this camera disparity map to "cast shadows" from the objects that are not visible
	 * in the reference (accurate) 3D model TODO: pre-filter to remove those that should be visible in pXpYD? At least remove
	 * low-confidence triplets.
	 * @param min_adisp_cam minimal absolute disparity difference for pXpYD_cam to consider
	 * @param min_rdisp_cam minimal relative disparity difference for pXpYD_cam to consider
	 * @param debug_level debug level
	 * @return copy of pXpYD with occluded elements nulled
	 */
	
	public double [][] filterBG (
			final TileProcessor tp,
			final double [][] pXpYD,
			final double max_overlap,
//			final double [][] pXpYD_cam,
			final double []   disparity_cam,
//			final double min_str_cam,
			final double min_adisp_cam,
			final double min_rdisp_cam,
			final int    dbg_tileX,
			final int    dbg_tileY,
			final int    debug_level
			){
		final int tilesX = tp.getTilesX();
		final int tilesY = tp.getTilesY();
		final int dbg_nTile = dbg_tileY * tilesX + dbg_tileX;
		final int tileSize = tp.getTileSize();
//		final double tileSize2 = tileSize * 2;
		final Thread[] threads = ImageDtt.newThreadArray(threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		final int tiles = tilesX*tilesY;
		final TileNeibs tn =  new TileNeibs(tilesX, tilesY);
		ArrayList<List<Integer>> fg_bg_list = new ArrayList<List<Integer>> (tiles);
		for (int i = 0; i < tiles; i++) {
			fg_bg_list.add(Collections.synchronizedList(new ArrayList<Integer>()));
		}
		final int offs_range = 1; // (max_overlap < 0.5) ? 2 : 1; 
		final double[][] pXpYD_filtered = pXpYD.clone();
		final AtomicInteger ai_num_tiles = new AtomicInteger(0);
		final AtomicInteger ai_num_removed = new AtomicInteger(0);
		final int overlap_radius = 4; // 9x9
		final int overlap_diameter = 2 * overlap_radius + 1; // 9
		final int overlap_size =  overlap_diameter * overlap_diameter; // 81
		final double scale_dist = 2.00 * overlap_radius / tileSize; // 1.0 for overlap_radius==4

		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				public void run() {
					for (int indx = ai.getAndIncrement(); indx < pXpYD.length; indx = ai.getAndIncrement()) if (pXpYD[indx] != null) {
						int tx = (int)Math.round(pXpYD[indx][0]/tileSize);
						int ty = (int)Math.round(pXpYD[indx][1]/tileSize);
						if ((debug_level > 0) && (indx == dbg_nTile)) { //(tx == dbg_tileX) && (ty == dbg_tileY)) {
							System.out.println("filterBG(): tx = "+tx+", ty="+ty+", indx="+indx);
							System.out.print("");
						}
						if ((tx >=0) && (ty >=0) && (tx < tilesX) && (ty < tilesY)) {
							int nTile = ty * tilesX + tx;
							synchronized(fg_bg_list.get(nTile)) {
								fg_bg_list.get(nTile).add(indx);
							}
							ai_num_tiles.getAndIncrement();
						} else {
							pXpYD_filtered[indx] = null;
						}
					}
				}
			};
		}		      
		ImageDtt.startAndJoin(threads);
		ai.set(0);
		// filter by the reference model
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				public void run() {
					boolean [] overlap_staging = new boolean [overlap_size];
					for (int nTile = ai.getAndIncrement(); nTile < tiles; nTile = ai.getAndIncrement()) if (fg_bg_list.get(nTile).size() > 0) {
						for (int tindx_bg: fg_bg_list.get(nTile)) {
							double [] txyd_bg = pXpYD[tindx_bg];
// 		final int offs_range = (max_overlap < 0.5) ? 2 : 1;
							if ((debug_level > -1) && (tindx_bg == dbg_nTile)) {
								System.out.println("filterBG(): tindx_bg="+tindx_bg+", nTile = "+nTile+", txyd_bg[0]="+txyd_bg[0]+", txyd_bg[1]="+txyd_bg[1]+", txyd_bg[2]="+txyd_bg[2]);
								System.out.print("");
							}
							Arrays.fill(overlap_staging, false);
							boolean some_overlap = false;
							for (int dty = -offs_range; dty <= offs_range; dty++) {
								for (int dtx = -offs_range; dtx <= offs_range; dtx++) {
									int nTile_fg = tn.getNeibIndex(nTile, dtx, dty);
									if (nTile_fg >= 0) 	for (int tindx_fg: fg_bg_list.get(nTile_fg)) if (tindx_fg != tindx_bg){
										double [] txyd_fg = pXpYD[tindx_fg];
										// check if FG is closer than BG (here does not have to be significantly closer
										if (txyd_fg[2] > txyd_bg[2]) {
											// see if there is any overlap
											double x_fg_bg = txyd_fg[0] - txyd_bg[0];
											double y_fg_bg = txyd_fg[1] - txyd_bg[1];
											if ((Math.abs(x_fg_bg) < tileSize) && (Math.abs(y_fg_bg) < tileSize)) {
												applyOverlap(
														overlap_staging,          // boolean [] staging,
														overlap_diameter,         // int        overlap_diameter,
														scale_dist,               // double     scale_dist,
														x_fg_bg,  // double     dx,
														y_fg_bg); // double     dy
												some_overlap = true;
											}
										}
									}
								}
							}
							// apply camera disparity map
							if (disparity_cam != null) {
								double ddisp = min_adisp_cam +  txyd_bg[2] * min_rdisp_cam;
								int tx = (int) Math.round(txyd_bg[0]/tileSize);
								int ty = (int) Math.round(txyd_bg[1]/tileSize);
								// Limit to 0.. max?
								int nTile_fg_center = ty * tilesX + tx;
								for (int dty = -offs_range; dty <= offs_range; dty++) {
									for (int dtx = -offs_range; dtx <= offs_range; dtx++) {
										int nTile_fg = tn.getNeibIndex(nTile_fg_center, dtx, dty);
										if ((nTile_fg >= 0) && (disparity_cam[nTile_fg] - txyd_bg[2] > ddisp)){
											double x_fg_bg = (tx + dtx + 0.5) * tileSize - txyd_bg[0];
											double y_fg_bg = (ty + dty + 0.5) * tileSize - txyd_bg[1];
											if ((Math.abs(x_fg_bg) < tileSize) && (Math.abs(y_fg_bg) < tileSize)) {
												applyOverlap(
														overlap_staging,          // boolean [] staging,
														overlap_diameter,         // int        overlap_diameter,
														scale_dist,               // double     scale_dist,
														x_fg_bg,  // double     dx,
														y_fg_bg); // double     dy
												some_overlap = true;
											}
										}
									}
								}
							}
							if (some_overlap) { // count actual overlap
								int nuv_overlap = 0;
								for (boolean staging_point: overlap_staging) {
									if (staging_point) {
										nuv_overlap++;
									}
								}
								double frac_overlap = 1.0 * nuv_overlap / overlap_staging.length;
								if (frac_overlap> max_overlap) {
									ai_num_removed.getAndIncrement();
									pXpYD_filtered[tindx_bg] = null; // OK that it still remains in the lists
								}
							}


							/*
											double overlapX = Math.max(tileSize2 - Math.abs(txyd_fg[0] - txyd_bg[0]), 0)/tileSize2;
											double overlapY = Math.max(tileSize2 - Math.abs(txyd_fg[1] - txyd_bg[1]), 0)/tileSize2;
											if ((overlapX * overlapY) > max_overlap) { // remove BG tile
												pXpYD_filtered[tindx_bg] = null; // OK that it still remains in the lists
												ai_num_removed_ref.getAndIncrement();
												if ((debug_level > -1) && (nTile==dbg_nTile)) {
													System.out.println("+++++++++++++++ filterBG(): nTile = "+nTile+
															", txyd_bg[0]="+txyd_bg[0]+", txyd_bg[1]="+txyd_bg[1]+", txyd_bg[2]="+txyd_bg[2]+
															", txyd_fg[0]="+txyd_fg[0]+", txyd_fg[1]="+txyd_fg[1]+", txyd_fg[2]="+txyd_fg[2]);
													System.out.print("");
												}
											}
							 */

						}
					}
				}
			};
		}		      
		ImageDtt.startAndJoin(threads);
		/*
		// Maybe remove here from the list tiles that are removed from pXpYD_filtered?
		if (pXpYD_cam != null) {
			ai.set(0);
			// filter by the reference model
			for (int ithread = 0; ithread < threads.length; ithread++) {
				threads[ithread] = new Thread() {
					public void run() {
						for (int tindx_fg = ai.getAndIncrement(); tindx_fg < pXpYD_cam.length; tindx_fg = ai.getAndIncrement()) if (pXpYD_cam[tindx_fg] != null) {
							double [] txyd_fg = pXpYD[tindx_fg];
							double ddisp = min_adisp_cam +  txyd_fg[2] * min_rdisp_cam; 
							int tx = (int)Math.round(txyd_fg[0]/tileSize);
							int ty = (int)Math.round(txyd_fg[1]/tileSize);
							if ((tx >=0) && (ty >=0) && (tx < tilesX) && (ty < tilesY)) {
								int nTile_fg = ty * tilesX + tx;
								for (int dy = -offs_range; dy <= offs_range; dy++) {
									for (int dx = -offs_range; dx <= offs_range; dx++) {
										int nTile_bg = tn.getNeibIndex(nTile_fg, dx, dy);
										if (nTile_bg >= 0) 	for (int tindx_bg: fg_bg_list.get(nTile_bg)){
											double [] txyd_bg = pXpYD[tindx_bg];
											if ((txyd_fg[2] - txyd_bg[2]) > ddisp) { // FG is significantly closer than BG
												double overlapX = Math.max(tileSize2 - Math.abs(txyd_fg[0] - txyd_bg[0]), 0)/tileSize2;
												double overlapY = Math.max(tileSize2 - Math.abs(txyd_fg[1] - txyd_bg[1]), 0)/tileSize2;
												if ((overlapX * overlapY) > max_overlap) { // remove BG tile
													pXpYD_filtered[tindx_bg] = null; // OK that it still remains in the lists
													ai_num_removed_cam.getAndIncrement();
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
		}
		*/
		if (debug_level > -1){
			System.out.println("filterBG(): num_all_tiles = "+ai_num_tiles.get()+
					", num_removed="+ ai_num_removed.get()+
					", remaining tiles="+(ai_num_tiles.get() - ai_num_removed.get()));
			System.out.print("");
		}

		return pXpYD_filtered;
	}

	
	private void applyOverlap(
			boolean [] staging,
			int        overlap_diameter,
			double     scale_dist,
			double     dx,
			double     dy
			) { 
		int ix0 = (int) Math.round(dx * scale_dist);
		int ix1 = ix0 + overlap_diameter;
		int iy0 = (int) Math.round(dy * scale_dist);
		int iy1 = iy0 + overlap_diameter;
		if (ix0 < 0) ix0 = 0;
		if (ix1 > overlap_diameter) ix1 = overlap_diameter;
		if (iy0 < 0) iy0 = 0;
		if (iy1 > overlap_diameter) iy1 = overlap_diameter;
		for (int iy = iy0; iy < iy1; iy++) {
			int line_start = iy * overlap_diameter;
			Arrays.fill(staging,  line_start + ix0, line_start + ix1, true);
		}
	}
	
	
	/*
	ArrayList<Integer> neib_list = new ArrayList<Integer>(20);
	// here each tlist is accessed exclusively
	Collections.sort(fg_bg_list.get(nTile), new Comparator<Integer>() {
	    @Override
	    public int compare(Integer lhs, Integer rhs) { // descending // ascending
	        return pXpYD[lhs][2] > pXpYD[rhs][2]  ? -1 : (pXpYD[lhs][2]  < pXpYD[rhs][2] ) ? 1 : 0;
	    }
	});
	 * 
  List list = Collections.synchronizedList(new ArrayList());
      ...
  synchronized(list) {
      Iterator i = list.iterator(); // Must be in synchronized block
      while (i.hasNext())
          foo(i.next());
  }

	 */

	public double [][] transformFromScenePxPyD(
			final double [][] pXpYD_scene, // tiles correspond to reference, pX,pY,D - for scene
			final double []   scene_xyz,   // camera center in world coordinates
			final double []   scene_atr,   // camera orientation relative to world frame
			final QuadCLT     scene_QuadClt,
			final QuadCLT     reference_QuadClt)
	{
		TileProcessor tp = reference_QuadClt.getTileProcessor();
		final int tilesX = tp.getTilesX();
		final int tilesY = tp.getTilesY();
		final int tiles = tilesX*tilesY;
		final double [][] pXpYD=               new double [tiles][];
		final ErsCorrection ersReferenceCorrection = reference_QuadClt.getErsCorrection();
		final ErsCorrection ersSceneCorrection =     scene_QuadClt.getErsCorrection();
		ersReferenceCorrection.setupERS(); // just in case - setUP using instance paRAMETERS
		ersSceneCorrection.setupERS();
		final Thread[] threads = ImageDtt.newThreadArray(threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				public void run() {
					for (int nTile = ai.getAndIncrement(); nTile < tiles; nTile = ai.getAndIncrement()) if (pXpYD_scene[nTile] != null) {
						double centerX = pXpYD_scene[nTile][0]; //  - shiftX;
						double centerY = pXpYD_scene[nTile][1]; //  - shiftX;
						double disparity = pXpYD_scene[nTile][2];
						if (disparity < 0) {
							disparity = 0.0;
						}
						if (scene_QuadClt == reference_QuadClt) {
							pXpYD[nTile] = new double [] {centerX, centerY, disparity};
						} else {
							pXpYD[nTile] = ersSceneCorrection.getImageCoordinatesERS( // ersCorrection - reference
									reference_QuadClt, // QuadCLT cameraQuadCLT, // camera station that got image to be to be matched 
									centerX,           // double px,                // pixel coordinate X in the reference view
									centerY,           // double py,                // pixel coordinate Y in the reference view
									disparity,         // double disparity,         // reference disparity 
									true,              // boolean distortedView,    // This camera view is distorted (diff.rect), false - rectilinear
									scene_xyz,         // double [] reference_xyz,  // this view position in world coordinates (typically ZERO3)
									scene_atr,         // double [] reference_atr,  // this view orientation relative to world frame  (typically ZERO3)
									true,              // boolean distortedCamera,  // camera view is distorted (false - rectilinear)
									ZERO3,             // double [] camera_xyz,     // camera center in world coordinates
									ZERO3,             // double [] camera_atr,     // camera orientation relative to world frame
									LINE_ERR);         // double    line_err)       // threshold error in scan lines (1.0)
							if (pXpYD[nTile] != null) {
								if (    (pXpYD[nTile][0] < 0.0) ||
										(pXpYD[nTile][1] < 0.0) ||
										(pXpYD[nTile][0] > ersReferenceCorrection.getSensorWH()[0]) ||
										(pXpYD[nTile][1] > ersReferenceCorrection.getSensorWH()[1])) {
									pXpYD[nTile] = null;
								}
							}
						}
					}
				}
			};
		}		      
		ImageDtt.startAndJoin(threads);
		return pXpYD; // Only disparity matters, pX, py - just for testing
	}
	
	
	
	/**
	 * Transform scene view to visually match with a reference scene. It is not accurate as it uses resampling and
	 * related low pass filtering.
	 * @param title image title to print
	 * @param dsrbg_camera_in - null (old compatibility) or [variable_length][tiles] array of disparity, strength, ... for the camera tiles 
	 * @param scene_xyz Scene X (right),Y (up), Z (negative away form camera) in the reference camera coordinates
	 *        or null to use scene instance coordinates.
	 * @param scene_atr Scene azimuth, tilt and roll (or null to use scene instance).
	 * @param scene_QuadClt Scene QuadCLT instance.
	 * @param reference_QuadClt Reference QuadCLT instance.
	 * @param iscale interpolation scale (use finer grid), typically 8
	 * @return Per-tile array of resampled {disparity,strength,red,blue,green} values (or nulls).
	 */
	public double [][] transformCameraVew(
			final String    title,
			final double [][] dsrbg_camera_in,
			final double [] scene_xyz, // camera center in world coordinates
			final double [] scene_atr, // camera orientation relative to world frame
			final QuadCLT   scene_QuadClt,
			final QuadCLT   reference_QuadClt,
			final int       iscale)
	{
		boolean debug = (title != null) && (title.length() > 0);
		final double line_error = 0.5;
		TileProcessor tp = reference_QuadClt.getTileProcessor();
		final int tilesX = tp.getTilesX();
		final int tilesY = tp.getTilesY();
		final int tiles = tilesX*tilesY;
		final int transform_size = tp.getTileSize();
		final int rel_num_passes = 10;
		final int num_passes =    transform_size; // * 2;

		final int stilesX = iscale*tilesX; 
		final int stilesY = iscale*tilesY;
		final int stiles = stilesX*stilesY;
		final double sigma = 0.5 * iscale; // was 0.5
		final double scale =  1.0 * iscale/transform_size;
		final double [][] dsrbg_camera = (dsrbg_camera_in == null) ? scene_QuadClt.getDSRBG() : dsrbg_camera_in;
		if (dsrbg_camera == null) {
			return null;
		}
		final double [][] ds =        new double [dsrbg_camera.length][stiles]; // null pointer
		for (int i = 0; i <ds.length; i++) {
			for (int j = 0; j <ds[i].length; j++) {
				ds[i][j] = Double.NaN;
			}
		}
		
		final ErsCorrection ersReferenceCorrection = reference_QuadClt.getErsCorrection();
		final ErsCorrection ersSceneCorrection =     scene_QuadClt.getErsCorrection();
		ersReferenceCorrection.setupERS(); // just in case - setUP using instance paRAMETERS
		ersSceneCorrection.setupERS();
		if (debug) {
			System.out.println("\ntransformCameraVew(): transformCameraVew(): >> "+title +" <<");
			System.out.println("transformCameraVew(): Reference scene ("+reference_QuadClt.getImageName()+"):");
			ersReferenceCorrection.printVectors(null, null);
			System.out.println("transformCameraVew(): Target scene ("+scene_QuadClt.getImageName()+"):");
			ersSceneCorrection.printVectors    (scene_xyz, scene_atr);
		}
		final Thread[] threads = ImageDtt.newThreadArray(threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		DoubleAccumulator [] azbuffer = new DoubleAccumulator[tiles];
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				public void run() {
					for (int nTile = ai.getAndIncrement(); nTile < tiles; nTile = ai.getAndIncrement()) {
						azbuffer[nTile] = new DoubleAccumulator (Double::max, Double.NEGATIVE_INFINITY);
					}
				}
			};
		}		      
		ImageDtt.startAndJoin(threads);
		ai.set(0);
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				public void run() {
					for (int nTile = ai.getAndIncrement(); nTile < tiles; nTile = ai.getAndIncrement()) if (dsrbg_camera[QuadCLT.DSRBG_STRENGTH][nTile] > 0.0) {
						double disparity = dsrbg_camera[QuadCLT.DSRBG_DISPARITY][nTile];
						if (!Double.isNaN(disparity)) {
							int tileY = nTile / tilesX;  
							int tileX = nTile % tilesX;
							double centerX = tileX * transform_size + transform_size/2; //  - shiftX;
							double centerY = tileY * transform_size + transform_size/2; //  - shiftY;
							if (disparity < 0) {
								disparity = 0.0;
							}
							double [] pXpYD = ersSceneCorrection.getImageCoordinatesERS( // ersCorrection - reference
									reference_QuadClt, // QuadCLT cameraQuadCLT, // camera station that got image to be to be matched 
									centerX,        // double px,                // pixel coordinate X in the reference view
									centerY,        // double py,                // pixel coordinate Y in the reference view
									disparity,      // double disparity,         // reference disparity 
									true,           // boolean distortedView,    // This camera view is distorted (diff.rect), false - rectilinear
									scene_xyz,      // double [] reference_xyz,  // this view position in world coordinates (typically ZERO3)
									scene_atr,      // double [] reference_atr,  // this view orientation relative to world frame  (typically ZERO3)
									true,           // boolean distortedCamera,  // camera view is distorted (false - rectilinear)
									ZERO3,          // double [] camera_xyz,     // camera center in world coordinates
									ZERO3,          // double [] camera_atr,     // camera orientation relative to world frame
									line_error); // LINE_ERR);      // double    line_err)       // threshold error in scan lines (1.0)
							if (pXpYD != null) {
								int px = (int) Math.round(pXpYD[0]/transform_size);
								int py = (int) Math.round(pXpYD[1]/transform_size);
								int spx = (int) Math.round(pXpYD[0]*scale);
								int spy = (int) Math.round(pXpYD[1]*scale);
								if ((px >= 0) && (py >= 0) && (px < tilesX) & (py < tilesY)) {
									double d = pXpYD[2];
									azbuffer[nTile].accumulate(pXpYD[2]);
									//Z-buffer
									if (!(d < azbuffer[nTile].get())) {
										if ((spx >= 0) && (spy >= 0) && (spx < stilesX) & (spy < stilesY)) {
											int sTile = spx + spy* stilesX;
											ds[QuadCLT.DSRBG_DISPARITY][sTile] = d; // pXpYD[2]; //reduce*
											for (int i = QuadCLT.DSRBG_STRENGTH; i < dsrbg_camera.length; i++) {
												ds[i][sTile] = dsrbg_camera[i][nTile]; // reduce * 
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
		ai.set(0);
		if (debug) {
			(new ShowDoubleFloatArrays()).showArrays(
				ds,
				stilesX,
				stilesY,
				true,
				"ds-0",
				new String[] {"D","S"});
		}
		
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				public void run() {
					for (int i = ai.getAndIncrement(); i < ds.length; i = ai.getAndIncrement()) {
						ds[i] = (new DoubleGaussianBlur()).blurWithNaN(
								ds[i], // double[] pixels,
								null,  // double [] in_weight, // or null
								stilesX, // int width,
								stilesY, // int height,
								sigma, // double sigmaX,
								sigma, // double sigmaY,
								0.01); // double accuracy);
					}
				}
			};
		}		      
		ImageDtt.startAndJoin(threads);
		ai.set(0);
		if (debug) {
			(new ShowDoubleFloatArrays()).showArrays(
					ds,
					stilesX,
					stilesY,
					true,
					"ds-1",
					new String[] {"D","S"});
		}
		
		final double [][] dsrbg_out = new double [dsrbg_camera.length][tiles];
		final int [][] num_non_nan = new int [dsrbg_out.length] [tiles];
		
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				public void run() {
					for (int nTile = ai.getAndIncrement(); nTile < tiles; nTile = ai.getAndIncrement()) if (dsrbg_camera[QuadCLT.DSRBG_STRENGTH][nTile] > 0.0) {
						int tileY = nTile / tilesX;  
						int tileX = nTile % tilesX;
						int tile =  tileX +  tileY *  tilesX;
						int stileY0 = tileY * iscale;
						int stileY1 = stileY0 + iscale; 
						int stileX0 = tileX * iscale;
						int stileX1 = stileX0 + iscale; 
						for (int stileY = stileY0; stileY < stileY1; stileY++) {
							for (int stileX = stileX0; stileX < stileX1; stileX++) {
								int stile = stileX + stileY * stilesX;
								for (int i = 0; i < dsrbg_out.length; i++) {
									double d = ds[i][stile];
									if (!Double.isNaN(d)) {
										num_non_nan[i][tile] ++;
										dsrbg_out[i][tile] += d;
									}	
								}
							}							
						}
					}
				}
			};
		}		      
		ImageDtt.startAndJoin(threads);
		if (debug) {
			(new ShowDoubleFloatArrays()).showArrays(
					dsrbg_out,
					tilesX,
					tilesY,
					true,
					"dsrbg_out-0");
		}
		
		
		for (int i = 0; i < dsrbg_out.length; i++) {
			for (int j = 0; j < tiles; j++) {
				if (num_non_nan[i][j] == 0) {
					dsrbg_out[i][j] = Double.NaN;
				} else {
					dsrbg_out[i][j]/=num_non_nan[i][j];
				}
			}
		}
		if (debug) {
			(new ShowDoubleFloatArrays()).showArrays(
					dsrbg_out,
					tilesX,
					tilesY,
					true,
					"dsrbg_out-1");
		}
		/*
		if (num_passes > 0) {
			for (int i = 0; i < dsrbg_out.length; i++) {
				dsrbg_out[i] = tp.fillNaNs(
						dsrbg_out[i],                  // double [] data,
						tilesX,                        //int       width, 
						2 * num_passes,                // int       grow,
						0.5 * Math.sqrt(2.0),          // double    diagonal_weight, // relative to ortho
						num_passes * rel_num_passes,   // int       num_passes,
						threadsMax);                   // final int threadsMax) // maximal number of threads to launch                         
			}
		}
		if (debug) {
			(new ShowDoubleFloatArrays()).showArrays(
					dsrbg_out,
					tilesX,
					tilesY,
					true,
					"dsrbg_out-2");
		}
		*/
		return dsrbg_out;
	}

	@Deprecated
	public double [][] transformCameraVewSingle( // not used
			double [] scene_xyz, // camera center in world coordinates
			double [] scene_atr, // camera orientation relative to world frame
			QuadCLT   scene_QuadClt,
			QuadCLT   reference_QuadClt,
			int       iscale)
	{
		TileProcessor tp = reference_QuadClt.getTileProcessor();
		int tilesX = tp.getTilesX();
		int tilesY = tp.getTilesY();
		int tiles = tilesX*tilesY;
		int transform_size = tp.getTileSize();
		int rel_num_passes = 10;
		int num_passes =    transform_size; // * 2;

		int stilesX = iscale*tilesX; 
		int stilesY = iscale*tilesY;
		int stiles = stilesX*stilesY;
		double sigma = 0.5 * iscale;
		double scale =  1.0 * iscale/transform_size;
		double [][] dsrbg_camera =    scene_QuadClt.getDSRBG();
///		double [][] dsrbg_reference = reference_QuadClt.getDSRBG();
		double [][] ds =        new double [dsrbg_camera.length][stiles];
		for (int i = 0; i <ds.length; i++) {
			for (int j = 0; j <ds[i].length; j++) {
				ds[i][j] = Double.NaN;
			}
		}
		
		ErsCorrection ersReferenceCorrection = reference_QuadClt.getErsCorrection();
		ErsCorrection ersSceneCorrection =     scene_QuadClt.getErsCorrection();
		ersReferenceCorrection.setupERS(); // just in case - setUP using instance paRAMETERS
		ersSceneCorrection.setupERS();
		double [] zbuffer = new double [tiles];
		for (int tileY = 0; tileY < tilesY; tileY++) {
			for (int tileX = 0; tileX < tilesX; tileX++) {
				int nTile = tileX + tileY * tilesX;
				double centerX = tileX * transform_size + transform_size/2; //  - shiftX;
				double centerY = tileY * transform_size + transform_size/2; //  - shiftY;
				double disparity = dsrbg_camera[QuadCLT.DSRBG_DISPARITY][nTile];
				if (disparity < 0) {
					disparity = 0.0;
				}
				// found that there are tiles with strength == 0.0, while disparity is not NaN
				if (!Double.isNaN(disparity) && (dsrbg_camera[QuadCLT.DSRBG_STRENGTH][nTile] > 0.0)) {
					// swapping reference <-> scene
					double [] pXpYD = ersSceneCorrection.getImageCoordinatesERS( // ersCorrection - reference
							reference_QuadClt, // QuadCLT cameraQuadCLT, // camera station that got image to be to be matched 
							centerX,        // double px,                // pixel coordinate X in the reference view
							centerY,        // double py,                // pixel coordinate Y in the reference view
							disparity,      // double disparity,         // reference disparity 
							true,           // boolean distortedView,    // This camera view is distorted (diff.rect), false - rectilinear
							scene_xyz,      // double [] reference_xyz,  // this view position in world coordinates (typically ZERO3)
							scene_atr,      // double [] reference_atr,  // this view orientation relative to world frame  (typically ZERO3)
							true,           // boolean distortedCamera,  // camera view is distorted (false - rectilinear)
							ZERO3,          // double [] camera_xyz,     // camera center in world coordinates
							ZERO3,          // double [] camera_atr,     // camera orientation relative to world frame
							0.5); // LINE_ERR);      // double    line_err)       // threshold error in scan lines (1.0)

					if (pXpYD != null) {
						int px = (int) Math.round(pXpYD[0]/transform_size);
						int py = (int) Math.round(pXpYD[1]/transform_size);
						int spx = (int) Math.round(pXpYD[0]*scale);
						int spy = (int) Math.round(pXpYD[1]*scale);
						if ((px >= 0) && (py >= 0) && (px < tilesX) & (py < tilesY)) {
							//Z-buffer
							if (!(pXpYD[2] < zbuffer[px + py* tilesX])) {
								zbuffer[px + py* tilesX] = pXpYD[2];
								if ((spx >= 0) && (spy >= 0) && (spx < stilesX) & (spy < stilesY)) {
									int sTile = spx + spy* stilesX;
									ds[QuadCLT.DSRBG_DISPARITY][sTile] = pXpYD[2]; //reduce*
									for (int i = QuadCLT.DSRBG_STRENGTH; i < dsrbg_camera.length; i++) {
										ds[i][sTile] = dsrbg_camera[i][nTile]; // reduce * 
									}
								}								
							}
						}
					}
				}
			}
		}
		
		//dsrbg_out[DSRBG_DISPARITY]
		for (int i = 0; i < ds.length; i++) {
			ds[i] = (new DoubleGaussianBlur()).blurWithNaN(
					ds[i], // double[] pixels,
					null,  // double [] in_weight, // or null
					stilesX, // int width,
					stilesY, // int height,
					sigma, // double sigmaX,
					sigma, // double sigmaY,
					0.01); // double accuracy);
		}
		double [][] dsrbg_out = new double [dsrbg_camera.length][tiles];
		int [][] num_non_nan = new int [dsrbg_out.length] [tiles];
		
		for (int stileY = 0; stileY < stilesY; stileY++) {
			int tileY = stileY / iscale; 
			for (int stileX = 0; stileX < stilesX; stileX++) {
				int tileX = stileX / iscale;
				int stile = stileX + stileY * stilesX;
				int tile =  tileX +  tileY *  tilesX;
				for (int i = 0; i < dsrbg_out.length; i++) {
					double d = ds[i][stile];
					if (!Double.isNaN(d)) {
						num_non_nan[i][tile] ++;
						dsrbg_out[i][tile] += d;
					}	
				}
			}
		}
		for (int i = 0; i < dsrbg_out.length; i++) {
			for (int j = 0; j < tiles; j++) {
				if (num_non_nan[i][j] == 0) {
					dsrbg_out[i][j] = Double.NaN;
				} else {
					dsrbg_out[i][j]/=num_non_nan[i][j];
				}
			}
		}

		if (num_passes > 0) {
			for (int i = 0; i < dsrbg_out.length; i++) {
				dsrbg_out[i] = tp.fillNaNs(
						dsrbg_out[i],                  // double [] data,
						tilesX,                        //int       width, 
						2 * num_passes,                // int       grow,
						0.5 * Math.sqrt(2.0),          // double    diagonal_weight, // relative to ortho
						num_passes * rel_num_passes,   // int       num_passes,
						threadsMax);                   // final int threadsMax) // maximal number of threads to launch                         
			}
		}
		return dsrbg_out;
	}

	
	
	
	public void adjustPairsDualPass(
			CLTParameters  clt_parameters,			
			double k_prev, 
			QuadCLT [] scenes, // ordered by increasing timestamps
			int debug_level
			)
	{
		double scale_two_omegas = 2.00; // ers angular velocities contain double omegas
		
		// Maintain pose differences and velocities, they are saved after pass2
		double [][][] scenes_xyzatr = new double [scenes.length][][]; // previous scene relative to the next one  
		double [][][] ers_xyzatr =    new double [scenes.length][][]; // previous scene relative to the next one  
		
		// for scanning around
		int [][] offset_start_corner = {{-1,-1},{1,-1},{1,1},{-1,1}}; //{x,y}
		int [][] offset_move = {{1,0},{0,1},{-1,0},{0,-1}};
		
		boolean pattern_mode = clt_parameters.ofp.pattern_mode;  //  true; // not yet used
		double max_rms_maybe = clt_parameters.ofp.max_rms_maybe; //  7.0;
		double max_rms_sure =  clt_parameters.ofp.max_rms_sure;  //  2.1;
		int pix_step =         clt_parameters.ofp.pix_step;      // 20;
		int search_rad =       clt_parameters.ofp.search_rad;    //  2; // 0;
		boolean next_predict = true; // Do not use macrotiles for 3-rd and next scenes - assume linear motion 
		boolean use_second_pass = false; // Probably not needed - instead just set velocities from neighbors
		
		// pass I -  adjust pairs using angular and linear velocities from intrascene ERS.
		for (int i = 1; i < scenes.length; i++) {
			QuadCLT reference_QuadClt = scenes[i];
			QuadCLT scene_QuadClt = scenes[i - 1];
			double [][] pose = new double [2][3];
			if (!clt_parameters.ofp.ignore_ers) { // ignore pose (all 0 when ERS is not available)
				pose = 	getPoseFromErs(  
						k_prev,
						reference_QuadClt,
						scene_QuadClt,
						debug_level);
			}
			reference_QuadClt.getErsCorrection().setupERSfromExtrinsics();
			scene_QuadClt.getErsCorrection().setupERSfromExtrinsics();
			double [] rms2 = new double[2];
			if ((i <= 1) || !next_predict) {
				// for scanning around
				double angle_per_step = reference_QuadClt.getGeometryCorrection().getCorrVector().getTiltAzPerPixel() * pix_step;
				double [] rmses = new double [(2*search_rad+1)*(2*search_rad+1)];
				Arrays.fill(rmses, max_rms_maybe+ 1.0); // undefined - set larger than worst
				double [][] atrs = new double [rmses.length][];
				int rad = 0, dir=0, n=0;
				int ntry = 0;
				int num_good=0;
				try_around:
					for (rad = 0; rad <= search_rad; rad++) {
						for (dir = 0; dir < ((rad==0)?1:4); dir++) {
							int n_range = (rad > 0) ? (2* rad) : 1;
							for (n = 0; n < n_range; n++) {
								int ix = rad*offset_start_corner[dir][0] + n * offset_move[dir][0];
								int iy = rad*offset_start_corner[dir][1] + n * offset_move[dir][1];;
								double [] atr = {
										pose[1][0]+ix * angle_per_step,
										pose[1][1]+iy * angle_per_step,
										pose[1][2]};
								if (debug_level > -2) {
									System.out.println("interPairsLMA(): trying adjustPairsLMA() with initial offset azimuth: "+
											atr[0]+", tilt ="+atr[1]);
								}
								//						double [][] new_pose =
								scenes_xyzatr[i] = adjustPairsLMA(
										clt_parameters,     // CLTParameters  clt_parameters,			
										reference_QuadClt, // QuadCLT reference_QuadCLT,
										scene_QuadClt, // QuadCLT scene_QuadCLT,
										pose[0], // xyz
										atr, // atr
										clt_parameters.ilp.ilma_lma_select,             // final boolean[]   param_select,
										clt_parameters.ilp.ilma_regularization_weights, //  final double []   param_regweights,
										rms2, // 			double []      rms, // null or double [2]
										null, // double [][]    dbg_img,
										max_rms_maybe, // double         max_rms,
										clt_parameters.ofp.debug_level_optical); // int debug_level)
								if ((scenes_xyzatr[i] != null) & !Double.isNaN(rms2[0])) {
									atrs[ntry] = scenes_xyzatr[i][1].clone();
									rmses[ntry] =  rms2[0];
									num_good++;
									if (rms2[0] <= max_rms_sure) {
										break try_around;
									}
								}
								ntry++;
							}
						}
					}
				if (ntry >= rmses.length) { // otherwise all is done already, last atr is set to scene
					if (debug_level>-3) {
						System.out.println("interPairsLMA(): left "+ntry+" candidates");
					}
					// rerun with the best atr
					int best_indx = 0;
					for (int j = 1; j < rmses.length; j++) {
						if (rmses[j] < rmses[best_indx]) {
							best_indx = j;
						}
					}
					if (debug_level > -2) {
						System.out.println("---- interPairsLMA(): re-trying adjustPairsLMA() with initial offset azimuth: "+
								atrs[best_indx][0]+", tilt ="+atrs[best_indx][1]+" best_indx="+best_indx); /// null pointer
					}

					scenes_xyzatr[i] = adjustPairsLMA(
							clt_parameters, // CLTParameters  clt_parameters,
							// FIXME: *********** update getPoseFromErs to use QUADCLTCPU ! **********				
							reference_QuadClt, // QuadCLT reference_QuadCLT,
							scene_QuadClt, // QuadCLT scene_QuadCLT,
							pose[0], // xyz
							atrs[best_indx], // pose[1], // atr
							clt_parameters.ilp.ilma_lma_select,             // final boolean[]   param_select,
							clt_parameters.ilp.ilma_regularization_weights, //  final double []   param_regweights,
							rms2, // 			double []      rms, // null or double [2]
							null, // double [][]    dbg_img,
							max_rms_maybe, // double         max_rms,
							clt_parameters.ofp.debug_level_optical); // 1); // -1); // int debug_level);
					if (debug_level > -1) {
						System.out.println("Pass 1 MACRO scene "+i+" (of "+ scenes.length+") "+
								reference_QuadClt.getImageName() + "/" + scene_QuadClt.getImageName()+(" Done.\n"));
					}
					if (num_good == 0) {
						System.out.println("****** Error ! Could not find any good match in a pair of consecutive scenes ! *****");
					}
				}

			} else {
				double dt_prev = scenes[i - 1].getTimeStamp() - scenes[i - 2].getTimeStamp();
				double dt_this = scenes[i - 0].getTimeStamp() - scenes[i - 1].getTimeStamp();
				double tscale = dt_this / dt_prev;
				scenes_xyzatr[i] = new double [2][3];
				for (int t = 0; t < scenes_xyzatr[i].length; t++) {
					for (int j = 0; j < scenes_xyzatr[i][t].length; j++) {
						scenes_xyzatr[i][t][j] = tscale * scenes_xyzatr[i-1][t][j];
					}
				}
			}

			scenes_xyzatr[i] = adjustPairsLMAInterscene(
					clt_parameters,                                 // CLTParameters  clt_parameters,
					reference_QuadClt,                              // QuadCLT reference_QuadCLT,
					null, // double []        ref_disparity, // null or alternative reference disparity
					null, // boolean []     reliable_ref, // null or bitmask of reliable reference tiles
					scene_QuadClt,                                  // QuadCLT scene_QuadCLT,
					scenes_xyzatr[i][0],                            // xyz
					scenes_xyzatr[i][1],                            // atr
					clt_parameters.ilp.ilma_lma_select,             // final boolean[]   param_select,
					clt_parameters.ilp.ilma_regularization_weights, // final double []   param_regweights,
					rms2,                                           // double []      rms, // null or double [2]
					clt_parameters.imp.max_rms,                     // double         max_rms,
					clt_parameters.imp.debug_level);                // 1); // -1); // int debug_level);
			if (debug_level > -1) {
				System.out.println("Pass 1 scene "+i+" (of "+ scenes.length+") "+
						reference_QuadClt.getImageName() + "/" + scene_QuadClt.getImageName()+(" Done.\n"));
			}
		} // for (int i = 1; i < scenes.length; i++) {
		
		// TODO: Add setting velocities w/o second pass
		
		// updating ERS from delta pose
		for (int i = 0; i < scenes.length; i++) {
			int i_prev = i - ((i > 0) ? 1 : 0);
			int i_next = i + ((i < (scenes.length - 1)) ? 1 : 0);
			double dt =  scenes[i_next].getTimeStamp() - scenes[i_prev].getTimeStamp();
			ers_xyzatr[i] = new double[2][3];
			// both prev and next are differential, so they are added
			
			if (i>0) {
				for (int j = 0; j < 3; j++) {
					ers_xyzatr[i][0][j] += scenes_xyzatr[i][0][j];
					ers_xyzatr[i][1][j] += scenes_xyzatr[i][1][j];
				}
			}
			if (i < (scenes.length - 1)) {
				for (int j = 0; j < 3; j++) {
					ers_xyzatr[i][0][j] += scenes_xyzatr[i + 1][0][j]; // null
					ers_xyzatr[i][1][j] += scenes_xyzatr[i + 1][1][j];
				}
			}
			for (int j = 0; j < 3; j++) {
				ers_xyzatr[i][0][j] *= 1.0 / dt;
				ers_xyzatr[i][1][j] *= scale_two_omegas / dt;
			}
			ers_xyzatr[i][1][0] = -ers_xyzatr[i][1][0]; /// TESTING!
			ers_xyzatr[i][1][2] = -ers_xyzatr[i][1][2]; /// TESTING!
		}
		
		boolean show_results = clt_parameters.ofp.show_result_images; //  true;
		if (show_results) {
			int dbg_w = ers_xyzatr.length;
			int dbg_h = 6;
			double [] dbg_img = new double [dbg_w * dbg_h];
			for (int dh = 0; dh  < dbg_h; dh++) {
				for (int dw = 0; dw  < dbg_w; dw++) {
					dbg_img[dh*dbg_w + dw] = ers_xyzatr[dw][dh / 3][dh % 3];
				}
			}
			(new ShowDoubleFloatArrays()).showArrays(
					dbg_img,
					dbg_w,
					dbg_h,
					"ers_xyzatr_fixed"); //	dsrbg_titles);
		}
		if (clt_parameters.ofp.lpf_pairs > 0.0) {
			ers_xyzatr = LPFVelocities(
					ers_xyzatr,                    // double [][][]  ers_xyzatr_in,
					clt_parameters.ofp.lpf_pairs); // double         half_run_range
			if (show_results) {
				int dbg_w = ers_xyzatr.length;
				int dbg_h = 6;
				double [] dbg_img = new double [dbg_w * dbg_h];
				for (int dh = 0; dh  < dbg_h; dh++) {
					for (int dw = 0; dw  < dbg_w; dw++) {
						if (ers_xyzatr[dw] == null) {
							System.out.println("adjustPairsDualPass(): ers_xyzatr["+dw+"] == null");
							continue;
						}
						if (ers_xyzatr[dw][dh/3] == null) {
							System.out.println("adjustPairsDualPass(): ers_xyzatr["+dw+"]["+(dh/3)+"] == null");
							continue;
						}
						dbg_img[dh*dbg_w + dw] = ers_xyzatr[dw][dh / 3][dh % 3];
					}
				}
				(new ShowDoubleFloatArrays()).showArrays(
						dbg_img,
						dbg_w,
						dbg_h,
						"ers_xyzatr_lpf"+clt_parameters.ofp.lpf_pairs); //	dsrbg_titles);
			}
		}
	//{camera_xyz0, camera_atr0}	
		if (!use_second_pass) {
			for (int i = 1; i < scenes.length; i++) {
				QuadCLT reference_QuadClt = scenes[i];
				QuadCLT scene_QuadClt = scenes[i - 1];
				ErsCorrection ers_reference = reference_QuadClt.getErsCorrection();
//				ErsCorrection ers_scene =     scene_QuadClt.getErsCorrection();
				ers_reference.addScene(scene_QuadClt.getImageName(),
						scenes_xyzatr[i][0],
						scenes_xyzatr[i][1],
						ers_xyzatr[i][0],		
						ers_xyzatr[i][1]		
						);
			}		
		} else {
			// pass II - set scene velocities from offsets to 1 before and one after, freeze ERS parameters, and
			// adjust other ones.
			boolean[]   param_select2 =     clt_parameters.ilp.ilma_lma_select.clone();             // final boolean[]   param_select,
			double []   param_regweights2 = clt_parameters.ilp.ilma_regularization_weights; //  final double []   param_regweights,
			// freeze reference ERS, free scene ERS
			for (int j = 0; j <3; j++) {
				param_select2[ErsCorrection.DP_DVX  + j] = false;
				param_select2[ErsCorrection.DP_DVAZ + j] = false;
				//			param_select2[ErsCorrection.DP_DSVX  + j] = false; // disabling, may check with high rot speed
				//			param_select2[ErsCorrection.DP_DSVAZ + j] = true;  // so far ers correction noise is too high to compare
				param_regweights2[ErsCorrection.DP_DSVX +  j] = 0.0;
				param_regweights2[ErsCorrection.DP_DSVAZ + j] = 0.0;
			}
			double [][][] scenes_xyzatr1 = new double [scenes.length][][]; // previous scene relative to the next one  

			for (int i = 1; i < scenes.length; i++) {
				QuadCLT reference_QuadClt = scenes[i];
				QuadCLT scene_QuadClt = scenes[i - 1];
				double [][] pose = 	scenes_xyzatr[i];
				ErsCorrection ers_reference = reference_QuadClt.getErsCorrection();
				ErsCorrection ers_scene =     scene_QuadClt.getErsCorrection();
				ers_reference.ers_wxyz_center_dt = ers_xyzatr[i][0].clone();
				ers_reference.ers_watr_center_dt = ers_xyzatr[i][1].clone();
				int i_prev = i - ((i > 0) ? 1 : 0);
				ers_reference.setupERS(); // just in case - setUP using instance paRAMETERS
				ers_scene.ers_wxyz_center_dt = ers_xyzatr[i_prev][0].clone();
				ers_scene.ers_watr_center_dt = ers_xyzatr[i_prev][1].clone();
				ers_scene.setupERS();     // just in case - setUP using instance paRAMETERS
				scenes_xyzatr1[i] = adjustPairsLMA(
						clt_parameters,     // CLTParameters  clt_parameters,			
						reference_QuadClt,  // QuadCLT reference_QuadCLT,
						scene_QuadClt, // QuadCLT scene_QuadCLT,
						pose[0], // xyz
						pose[1], // atr
						param_select2,             // final boolean[]   param_select,
						param_regweights2, //  final double []   param_regweights,
						null, // 			double []      rms, // null or double [2]
						null, // double [][]    dbg_img,
						0.0, // double         max_rms,
						debug_level); // int debug_level)
				ers_reference.addScene(scene_QuadClt.getImageName(),
						scenes_xyzatr1[i][0],
						scenes_xyzatr1[i][1],
						ers_scene.getErsXYZ_dt(),		
						ers_scene.getErsATR_dt()		
						);
				if (debug_level > -1) {
					System.out.println("Pass 2 scene "+i+" (of "+ scenes.length+") "+
							reference_QuadClt.getImageName() + "/" + scene_QuadClt.getImageName()+" Done.");
				}
			}

			if (show_results) {
				double [][][] ers_xyzatr_dt = new double [scenes.length][][];
				for (int i = 0; i < scenes.length; i++) {
					ErsCorrection ers_scene =     scenes[i].getErsCorrection();
					ers_xyzatr_dt[i] = new double[][] {
						ers_scene.ers_wxyz_center_dt,
						ers_scene.ers_watr_center_dt};
				}
				int dbg_w = ers_xyzatr_dt.length;
				int dbg_h = 6;
				double [] dbg_img = new double [dbg_w * dbg_h];
				for (int dh = 0; dh  < dbg_h; dh++) {
					for (int dw = 0; dw  < dbg_w; dw++) {
						dbg_img[dh*dbg_w + dw] = ers_xyzatr_dt[dw][dh / 3][dh % 3];
					}
				}
				(new ShowDoubleFloatArrays()).showArrays(
						dbg_img,
						dbg_w,
						dbg_h,
						"ers_xyzatr_adjusted"); //	dsrbg_titles);
			}


			//TODO: Add update ers with lpf here?
			//		boolean show_results = true;
			if (show_results) {
				int dbg_w = scenes_xyzatr1.length;
				int dbg_h = 6;
				double [] dbg_img = new double [dbg_w * dbg_h];
				Arrays.fill(dbg_img, Double.NaN);
				for (int dh = 0; dh  < dbg_h; dh++) {
					for (int dw = 0; dw  < dbg_w; dw++) {
						if (scenes_xyzatr1[dw] != null) {
							dbg_img[dh*dbg_w + dw] = scenes_xyzatr1[dw][dh / 3][dh % 3];
						}
					}
				}
				(new ShowDoubleFloatArrays()).showArrays(
						dbg_img,
						dbg_w,
						dbg_h,
						"scenes_xyzatr1"); //	dsrbg_titles);
			}
		}


		
		// Update ERS velocities with running average?
		
		for (int i = 1; i < scenes.length; i++) {
			QuadCLT reference_QuadClt = scenes[i];
			reference_QuadClt.saveInterProperties( // save properties for interscene processing (extrinsics, ers, ...)
		            null, // String path,             // full name with extension or w/o path to use x3d directory
		            debug_level+1);
		}		
	}
	/*
	public class VideoPathParams{
		public String path;
		public int    width;
		public int    height;
		public int    gap;
		public double fps;
		public double real_fps;
	}
	*/
	/**
	 * Build series of poses from just a single (reference) scene.
	 * @param quadCLT_main
	 * @param ref_index
	 * @param ref_step
	 * @param clt_parameters
	 * @param debayerParameters
	 * @param colorProcParameters
	 * @param channelGainParameters
	 * @param rgbParameters
	 * @param equirectangularParameters
	 * @param properties
	 * @param reset_from_extrinsics
	 * @param threadsMax
	 * @param updateStatus
	 * @param debugLevel
	 * @return null on failure, or String path to the reference model directory 
	 * @throws Exception
	 */
    public String buildSeries(
    		boolean                                              batch_mode,
			QuadCLT                                              quadCLT_main, // tiles should be set
			int                                                  ref_index, // -1 - last
//			int                                                  start_index, 
			CLTParameters                                        clt_parameters,
			EyesisCorrectionParameters.DebayerParameters         debayerParameters,
			ColorProcParameters                                  colorProcParameters,
			CorrectionColorProc.ColorGainsParameters             channelGainParameters,
			EyesisCorrectionParameters.RGBParameters             rgbParameters,
			EyesisCorrectionParameters.EquirectangularParameters equirectangularParameters,
			Properties                                           properties,
			boolean                                              reset_from_extrinsics,
			String [][]                                          videos, // null or String[1][] list of generated avi or webm paths
			int [][]                                             stereo_widths, // null or int[1][] matching videos -
			int []                                               start_ref_pointers, // [0] - earliest valid scene, [1] ref_index
			                                                     // each element is 0 for non-stereo and full width for stereo 
			final int        threadsMax,  // maximal number of threads to launch
			final boolean    updateStatus,
			final int        debugLevel)  throws Exception
    {
    	int earliest_scene = 0;      // increase on failure
    	
    	boolean build_ref_dsi =              clt_parameters.imp.force_ref_dsi;
    	boolean force_initial_orientations = clt_parameters.imp.force_orientations ;
    	
//    	boolean build_interscene =   clt_parameters.imp.force_interscene;
    	
    	int min_num_interscene =             clt_parameters.imp.min_num_interscene; //  2; // make from parameters, should be >= 1
    	int min_num_orient =                 clt_parameters.imp.min_num_orient;    //  2; // make from parameters, should be >= 1
    	
    	boolean export_images =              clt_parameters.imp.export_images;
    	boolean export_dsi_image =           clt_parameters.imp.export_ranges;
    	boolean export_ml_files =             clt_parameters.imp.export_ml_files;
    	
    	boolean show_dsi_image =             clt_parameters.imp.show_ranges && !batch_mode;
    	boolean show_images =                clt_parameters.imp.show_images && !batch_mode;
//    	boolean show_color_nan =     clt_parameters.imp.show_color_nan;
//    	boolean show_mono_nan =      clt_parameters.imp.show_mono_nan;
    	int min_num_scenes =                 clt_parameters.imp.min_num_scenes; // abandon series if there are less than this number of scenes in it 
    	
    	boolean show_images_bgfg =           clt_parameters.imp.show_images_bgfg && !batch_mode;
    	boolean show_images_mono =           clt_parameters.imp.show_images_mono && !batch_mode;

		double  range_disparity_offset =     clt_parameters.imp.range_disparity_offset ; //   -0.08;
		double  range_min_strength =         clt_parameters.imp.range_min_strength ; // 0.5;
		double  range_max =                  clt_parameters.imp.range_max ; //  5000.0;
		
		boolean [] save_mapped_mono_color =  {clt_parameters.imp.save_mapped_mono, clt_parameters.imp.save_mapped_color};
		boolean [] gen_avi_mono_color =      {clt_parameters.imp.gen_avi_mono, clt_parameters.imp.gen_avi_color};
//		boolean show_mapped_color =  clt_parameters.imp.show_mapped_color && !batch_mode;
//		boolean show_mapped_mono =   clt_parameters.imp.show_mapped_mono && !batch_mode;
		boolean [] show_mono_color =        {
				clt_parameters.imp.show_mapped_mono && !batch_mode,
				clt_parameters.imp.show_mapped_color && !batch_mode};
		
		boolean [] gen_seq_mono_color =  {
				save_mapped_mono_color[0]  || gen_avi_mono_color[0] || show_mono_color[0],
				save_mapped_mono_color[1]  || gen_avi_mono_color[1] || show_mono_color[1]} ;
		// skip completely if no color or mono, tiff or video
		boolean generate_mapped =            clt_parameters.imp.generate_mapped &&
				(gen_seq_mono_color[0] || gen_seq_mono_color[1]);   // generate sequences - Tiff and/or video
		
		
		boolean [] annotate_mono_color = {clt_parameters.imp.annotate_mono,clt_parameters.imp.annotate_color}; 
		boolean annotate_transparent_mono = clt_parameters.imp.annotate_transparent_mono;
		
		boolean [] generate_modes3d = {
				clt_parameters.imp.generate_raw,
				clt_parameters.imp.generate_inf,
				clt_parameters.imp.generate_fg,
				clt_parameters.imp.generate_bg};
		boolean generate_stereo =            clt_parameters.imp.generate_stereo;
		
//		double [] stereo_bases =         clt_parameters.imp.stereo_bases; // {0.0, 200.0, 500.0, 1000.0}; 
		double [][] stereo_views =       clt_parameters.imp.stereo_views; // {0.0, 200.0, 500.0, 1000.0}; 
		boolean [] generate_stereo_var = clt_parameters.imp.generate_stereo_var;
		
		boolean stereo_merge =       clt_parameters.imp.stereo_merge;
		boolean anaglyth_en =        clt_parameters.imp.anaglyth_en;
		final Color anaglyph_left =  clt_parameters.imp.anaglyph_left;
		final Color anaglyph_right = clt_parameters.imp.anaglyph_right;
		
		int     stereo_gap =         clt_parameters.imp.stereo_gap;
//		double  stereo_intereye =    clt_parameters.imp.stereo_intereye;
//		double  stereo_phone_width = clt_parameters.imp.stereo_phone_width;
		
		int     extra_hor_tile =     clt_parameters.imp.extra_hor_tile;
		int     extra_vert_tile =    clt_parameters.imp.extra_vert_tile;
		boolean crop_3d =            clt_parameters.imp.crop_3d;
		int     sensor_mask =        clt_parameters.imp.sensor_mask; // -1 - all
		// video
		double  video_fps =          clt_parameters.imp.video_fps;
		int     mode_avi =           clt_parameters.imp.mode_avi;
		int     avi_JPEG_quality =   clt_parameters.imp.avi_JPEG_quality; // 90;
	    boolean run_ffmpeg =         clt_parameters.imp.run_ffmpeg;
		String  video_ext =          clt_parameters.imp.video_ext;
		String  video_codec =        clt_parameters.imp.video_codec.toLowerCase();
//		String  video_extra =        clt_parameters.imp.video_extra;
		int     video_crf =          clt_parameters.imp.video_crf;
		boolean remove_avi =         clt_parameters.imp.remove_avi;
		boolean um_mono =            clt_parameters.imp.um_mono;
		double  um_sigma =           clt_parameters.imp.um_sigma;
		double  um_weight =          clt_parameters.imp.um_weight;
		boolean mono_fixed =         clt_parameters.imp.mono_fixed;
		double  mono_range =         clt_parameters.imp.mono_range;
		
		boolean reuse_video =        clt_parameters.imp.reuse_video &&
				(gen_seq_mono_color[0] || gen_seq_mono_color[1]);   // generate sequences - Tiff and/or video

		
		final Color annotate_color_color = clt_parameters.imp.annotate_color_color;
		final Color annotate_color_mono =  clt_parameters.imp.annotate_color_mono;
		
//		boolean  readjust_orient =   clt_parameters.imp.readjust_orient;
		boolean  test_ers =          clt_parameters.imp.test_ers && !batch_mode;
		int     test_ers0 =          clt_parameters.imp.test_ers0; // try adjusting a pair of scenes with ERS. Reference scene index
		int     test_ers1 =          clt_parameters.imp.test_ers1; // try adjusting a pair of scenes with ERS. Other scene index
		test_ers &= (test_ers0 >= 0) && (test_ers1 >= 0);
		final int        debugLevelInner=clt_parameters.batch_run? -2: debugLevel; // copied from TQ
		
		double  min_ref_str =        clt_parameters.imp.min_ref_str;
		
		
		if (reuse_video) { // disable all other options
			generate_mapped = false;
			export_images = false;
			export_dsi_image = false;
			show_dsi_image = false;
			export_ml_files = false;
			test_ers = false;
		}
		
		ArrayList<String> video_list = new ArrayList<String>();
		ArrayList<Integer> stereo_widths_list = new ArrayList<Integer>();
		if ((quadCLT_main != null) && (quadCLT_main.getGPU() != null)) {
			quadCLT_main.getGPU().resetGeometryCorrection();
			quadCLT_main.gpuResetCorrVector(); // .getGPU().resetGeometryCorrectionVector();
		}
		// final boolean    batch_mode = clt_parameters.batch_run;
		this.startTime=System.nanoTime();
		String [] sourceFiles0=quadCLT_main.correctionsParameters.getSourcePaths();
		QuadCLT.SetChannels [] set_channels_main = quadCLT_main.setChannels(debugLevel);
		if ((set_channels_main == null) || (set_channels_main.length==0)) {
			System.out.println("buildSeriesTQ(): No files to process (of "+sourceFiles0.length+")");
			return null;
		}
		// set_channels will include all 99 scenes even as quadCLTs.length matches ref_index
		QuadCLT.SetChannels [] set_channels=quadCLT_main.setChannels(debugLevel); 
		if (ref_index < 0) {
			ref_index += set_channels.length; 
		}
		if (start_ref_pointers != null) {
			start_ref_pointers[0] = 0;
			start_ref_pointers[1] = ref_index;
		}
		
		QuadCLT [] quadCLTs = new QuadCLT [ref_index+1]; //  [set_channels.length];
		
		//start_index
		
		
		double [][][] scenes_xyzatr = new double [quadCLTs.length][][]; // previous scene relative to the next one
		scenes_xyzatr[ref_index] = new double[2][3]; // all zeros
		// See if build_ref_dsi is needed
		if (!build_ref_dsi) {
			quadCLTs[ref_index] = (QuadCLT) quadCLT_main.spawnNoModelQuadCLT( // will conditionImageSet
					set_channels[ref_index].set_name,
					clt_parameters,
					colorProcParameters, //
					threadsMax,
					debugLevel-2);
			if ((quadCLTs[ref_index] == null) || !quadCLTs[ref_index].dsiExists()) {
				if (debugLevel >-2) {
					System.out.println("DSI data for scene "+set_channels[ref_index].set_name+" does not exist, forcing its calculation.");
				}
				build_ref_dsi = true;
			}
		}
		
		// 1. Reference scene DSI
		if (build_ref_dsi) {
			TwoQuadCLT.copyJP4src( // actually there is no sense to process multiple image sets. Combine with other
					// processing?
					set_channels[ref_index].set_name, // String set_name
					quadCLT_main, // QuadCLT quadCLT_main,
					null, // QuadCLT quadCLT_aux,
					clt_parameters, // EyesisCorrectionParameters.DCTParameters dct_parameters,
					true, // boolean                                  skip_existing,
					true, // false, // boolean                                  search_KML,
					debugLevel);
			quadCLTs[ref_index] = (QuadCLT) quadCLT_main.spawnNoModelQuadCLT( // will conditionImageSet
					set_channels[ref_index].set_name,
					clt_parameters,
					colorProcParameters, //
					threadsMax,
					debugLevel-2);
			quadCLTs[ref_index].saveQuadClt(); // to re-load new set of Bayer images to the GPU (do nothing for CPU) and Geometry
			quadCLTs[ref_index].preExpandCLTQuad3d( // returns ImagePlus, but it already should be saved/shown
	                clt_parameters,
	                debayerParameters,
	                colorProcParameters,
	                rgbParameters,
	                threadsMax,  // maximal number of threads to launch
	                updateStatus,
	                debugLevelInner);
			double [][] aux_last_scan = TileProcessor.getDSLMA(
					quadCLTs[ref_index].tp.clt_3d_passes.get( quadCLTs[ref_index].tp.clt_3d_passes.size() -1),
					false); // boolean force_final);
			double [][] dsi = new double [TwoQuadCLT.DSI_SLICES.length][];
			dsi[TwoQuadCLT.DSI_DISPARITY_AUX] =     aux_last_scan[0];
			dsi[TwoQuadCLT.DSI_STRENGTH_AUX] =      aux_last_scan[1];
			dsi[TwoQuadCLT.DSI_DISPARITY_AUX_LMA] = aux_last_scan[2];
			
			if (quadCLT_main.correctionsParameters.clt_batch_dsi_cm_strength) {
				CLTPass3d scan = new CLTPass3d(quadCLTs[ref_index].tp);
				scan.setTileOpDisparity(aux_last_scan[0]); // measure w/o LMA, use just strength
				quadCLTs[ref_index].CLTMeas( // perform single pass according to prepared tiles operations and disparity
						clt_parameters, // EyesisCorrectionParameters.CLTParameters           clt_parameters,
						scan,           // final CLTPass3d   scan,
						false,          // final boolean     save_textures,
						false,          // final boolean       need_diffs,     // calculate diffs even if textures are not needed 
						0,              // final int         clust_radius,
						true,           // final boolean     save_corr,
						false,          // final boolean     run_lma, // =    true;
						threadsMax,     // final int         threadsMax,  // maximal number of threads to launch
						updateStatus,   // final boolean     updateStatus,
						debugLevel);    // final int         debugLevel);
				dsi[TwoQuadCLT.DSI_STRENGTH_AUX] =      scan.getStrength();
				if (debugLevel > 1) {
					quadCLTs[ref_index].tp.showScan(
							scan, // CLTPass3d   scan,
							"test-strength");
				}
			} else {
				dsi[TwoQuadCLT.DSI_STRENGTH_AUX] =      aux_last_scan[1];
			}
			quadCLTs[ref_index].saveDSIAll (
					"-DSI_MAIN", // String suffix, // "-DSI_MAIN"
					dsi);
			quadCLTs[ref_index].saveInterProperties( // save properties for interscene processing (extrinsics, ers, ...)
					null, // String path,             // full name with extension or w/o path to use x3d directory
					//							 	null, // Properties properties,   // if null - will only save extrinsics)
					debugLevel);
			// Copy source files to the model directory - they are needed for poses and interscene
			for (int scene_index =  ref_index - 1; scene_index >= 0 ; scene_index--) {
				TwoQuadCLT.copyJP4src( // actually there is no sense to process multiple image sets. Combine with other
						// processing?
						set_channels[scene_index].set_name, // String set_name
						quadCLT_main, // QuadCLT quadCLT_main,
						null, // QuadCLT quadCLT_aux,
						clt_parameters, // EyesisCorrectionParameters.DCTParameters dct_parameters,
						true, // boolean                                  skip_existing,					
						true, // false, // boolean                                  search_KML,
						debugLevel);
			} // split cycles to remove output clutter
			
			quadCLTs[ref_index].set_orient(0); // reset orientations
			quadCLTs[ref_index].set_accum(0);  // reset accumulations ("build_interscene") number
			earliest_scene = 0;  // reset failures, try to use again all scenes
			
		} // if (build_ref_dsi) {
		
		quadCLTs[ref_index] = (QuadCLT) quadCLT_main.spawnQuadCLT( // restores dsi from "DSI-MAIN"
				set_channels[ref_index].set_name,
				clt_parameters,
				colorProcParameters, //
				threadsMax,
				debugLevel);
		
		quadCLTs[ref_index].setDSRBG(
				clt_parameters, // CLTParameters  clt_parameters,
				threadsMax,     // int            threadsMax,  // maximal number of threads to launch
				updateStatus,   // boolean        updateStatus,
				debugLevel);    // int            debugLevel)
		ErsCorrection ers_reference = quadCLTs[ref_index].getErsCorrection();
		
		if (!force_initial_orientations) {
			if ((quadCLTs[ref_index].getNumOrient() == 0) ||(!quadCLTs[ref_index].propertiesContainString (".scenes_"))) {
				if (debugLevel >-2) {
					System.out.println("Egomotion data for scene "+set_channels[ref_index].set_name+" does not exist, forcing its calculation.");
				}
				force_initial_orientations = true;
			}
		}
		
		// Build initial orientations
		
		if (force_initial_orientations && !reuse_video) {
			double maximal_series_rms = 0.0;
			double [] lma_rms = new double[2];
			double [] use_atr = null;
			for (int scene_index =  ref_index - 1; scene_index >= 0 ; scene_index--) {
				quadCLTs[scene_index] = (QuadCLT) quadCLT_main.spawnNoModelQuadCLT(
						set_channels[scene_index].set_name,
						clt_parameters,
						colorProcParameters, //
						threadsMax,
						debugLevel-2);
			} // split cycles to remove output clutter
			int debug_scene = -15;
			
			boolean [] reliable_ref = null;
			if (min_ref_str > 0.0) {
				reliable_ref = quadCLTs[ref_index].getReliableTiles( // will be null if does not exist.
						min_ref_str,   // double min_strength,
						true);         // boolean needs_lma);
			}
			
			for (int scene_index =  ref_index - 1; scene_index >= earliest_scene ; scene_index--) {
				if (scene_index == debug_scene) {
					System.out.println("scene_index = "+scene_index);
					System.out.println("scene_index = "+scene_index);
				}
				QuadCLT scene_QuadClt = quadCLTs[scene_index];
				// get initial xyzatr:
				if (scene_index ==  ref_index - 1) { // search around for the best fit
					use_atr = spiralSearchATR(
							clt_parameters, // CLTParameters            clt_parameters,
							quadCLTs[ref_index], // QuadCLT             reference_QuadClt,
							scene_QuadClt, // QuadCLT                   scene_QuadClt,
							reliable_ref, // ********* boolean []                reliable_ref,
							debugLevel);
					if (use_atr == null) {
						earliest_scene = scene_index + 1;
						if (debugLevel > -3) {
							System.out.println("Pass multi scene "+scene_index+" (of "+ quadCLTs.length+") "+
									quadCLTs[ref_index].getImageName() + "/" + scene_QuadClt.getImageName()+
									" spiralSearchATR() FAILED. Setting earliest_scene to "+earliest_scene);
						}
						// set this and all previous to null
						for (; scene_index >= 0 ; scene_index--) {
							ers_reference.addScene(quadCLTs[scene_index].getImageName(), null);
						}
						break; // failed with even first before reference
					}
					scenes_xyzatr[scene_index] = new double [][] {new double[3], use_atr};
				} else { // assume linear motion
					double [][] last_diff = ErsCorrection.combineXYZATR(
							scenes_xyzatr[scene_index+1],
							ErsCorrection.invertXYZATR(scenes_xyzatr[scene_index+2]));
					scenes_xyzatr[scene_index] = ErsCorrection.combineXYZATR(
							scenes_xyzatr[scene_index+1],
							last_diff);
				}
				// Refine with LMA
				scenes_xyzatr[scene_index] = adjustPairsLMAInterscene(
						clt_parameters,      // CLTParameters  clt_parameters,
						quadCLTs[ref_index], // QuadCLT reference_QuadCLT,
						null,                // double []        ref_disparity, // null or alternative reference disparity
						reliable_ref,        // boolean []     reliable_ref, // null or bitmask of reliable reference tiles
						scene_QuadClt,                                  // QuadCLT scene_QuadCLT,
						scenes_xyzatr[scene_index][0],                                // xyz
						scenes_xyzatr[scene_index][1],                                // atr
						clt_parameters.ilp.ilma_lma_select,                                  // final boolean[]   param_select,
						clt_parameters.ilp.ilma_regularization_weights,                              //  final double []   param_regweights,
						lma_rms,                                        // double []      rms, // null or double [2]
						clt_parameters.imp.max_rms,                     // double         max_rms,
						clt_parameters.imp.debug_level);                // 1); // -1); // int debug_level);
				if (scenes_xyzatr[scene_index] == null) {
					earliest_scene = scene_index + 1;
					if (debugLevel > -3) {
						System.out.println("Pass multi scene "+scene_index+" (of "+ quadCLTs.length+") "+
								quadCLTs[ref_index].getImageName() + "/" + scene_QuadClt.getImageName()+
								" FAILED. Setting earliest_scene to "+earliest_scene);
					}
					// set this and all previous to null
					for (; scene_index >= 0 ; scene_index--) {
						ers_reference.addScene(quadCLTs[scene_index].getImageName(), null);
					}
					break;
				}
				
				// TODO: Maybe after testing high-ers scenes - restore velocities from before/after scenes
				ers_reference.addScene(scene_QuadClt.getImageName(),
						scenes_xyzatr[scene_index][0],
						scenes_xyzatr[scene_index][1],
						ZERO3, // ers_scene.getErsXYZ_dt(),		
						ZERO3 // ers_scene.getErsATR_dt()		
						);
			    if (lma_rms[0] > maximal_series_rms) {
			        maximal_series_rms = lma_rms[0];
			    }

				if (debugLevel > -3) {
					System.out.println("Pass multi scene "+scene_index+" (of "+ quadCLTs.length+") "+
							quadCLTs[ref_index].getImageName() + "/" + scene_QuadClt.getImageName()+
							" Done. RMS="+lma_rms[0]+", maximal so far was "+maximal_series_rms);
				}
			} // for (int scene_index =  ref_index - 1; scene_index >= 0 ; scene_index--)
			if ((ref_index - earliest_scene + 1) < min_num_scenes) {
				System.out.println("Total number of useful scenes = "+(ref_index - earliest_scene + 1)+
						" < "+min_num_scenes+". Scrapping this series.");
				if (start_ref_pointers != null) {
					start_ref_pointers[0] = earliest_scene;
				}
				return null;
			}
			if (earliest_scene > 0) {
				System.out.println("Not all scenes matched, earliest useful scene = "+earliest_scene+
						" (total number of scenes = "+(ref_index - earliest_scene + 1)+
						").Maximal RMSE was "+maximal_series_rms);
			} else if (debugLevel > -4) {
				System.out.println("All multi scene passes are Done. Maximal RMSE was "+maximal_series_rms);
			}
			quadCLTs[ref_index].set_orient(1); // first orientation
			quadCLTs[ref_index].set_accum(0);  // reset accumulations ("build_interscene") number
			quadCLTs[ref_index].saveInterProperties( // save properties for interscene processing (extrinsics, ers, ...)  // null pointer
					null, // String path,             // full name with extension or w/o path to use x3d directory
					debugLevel+1);
		} else {// if (build_orientations) {
			if (!reuse_video) { // reuse_video only uses reference scene
				for (int scene_index =  ref_index - 1; scene_index >= earliest_scene ; scene_index--) {
					quadCLTs[scene_index] = (QuadCLT) quadCLT_main.spawnNoModelQuadCLT( // restores image data
							set_channels[scene_index].set_name,
							clt_parameters,
							colorProcParameters, //
							threadsMax,
							debugLevel-2);
				}
			}
		}
		// just in case that orientations were calculated before:
		earliest_scene = getEarliestScene(quadCLTs);
		
		double [][] combo_dsn_final = null;
		while (!reuse_video && ((quadCLTs[ref_index].getNumOrient() < min_num_orient) || (quadCLTs[ref_index].getNumAccum() < min_num_interscene))) {
			if ((quadCLTs[ref_index].getNumAccum() < min_num_interscene) &&
					((quadCLTs[ref_index].getNumAccum() <  quadCLTs[ref_index].getNumOrient())||
							(quadCLTs[ref_index].getNumOrient() >= min_num_orient))) {
				// should skip scenes w/o orientation 06/29/2022
				
				combo_dsn_final = intersceneExport( // result indexed by COMBO_DSN_TITLES, COMBO_DSN_INDX_***
						clt_parameters,      // CLTParameters        clt_parameters,
						ers_reference,       // ErsCorrection        ers_reference,
						quadCLTs,            // QuadCLT []           scenes,
						colorProcParameters, // ColorProcParameters  colorProcParameters,
						debugLevel);         // int                  debug_level
				
				
				
				quadCLTs[ref_index].inc_accum();
				// save with updated num_accum
				quadCLTs[ref_index].saveInterProperties( // save properties for interscene processing (extrinsics, ers, ...)
						null, // String path,             // full name with extension or w/o path to use x3d directory
						debugLevel+1);
			}

			//if (readjust_orient) {
			if (quadCLTs[ref_index].getNumOrient() < min_num_orient) {
				/*
				if (!force_initial_orientations && !build_interscene) {
					for (int nscene = 0; nscene < (quadCLTs.length -1); nscene++) {
						quadCLTs[nscene] = (QuadCLT) quadCLT_main.spawnNoModelQuadCLT( // restores image data
								set_channels[nscene].set_name,
								clt_parameters,
								colorProcParameters, //
								threadsMax,
								debugLevel-2);
					}
				}
				*/
				boolean [] reliable_ref = null;
				if (min_ref_str > 0.0) {
					reliable_ref = quadCLTs[ref_index].getReliableTiles( // will be null if does not exist.
							min_ref_str,   // double min_strength,
							true);         // boolean needs_lma);
				}
				earliest_scene=reAdjustPairsLMAInterscene( // after combo dsi is available and preliminary poses are known
						clt_parameters,     // CLTParameters  clt_parameters,
						reliable_ref,       // boolean []     reliable_ref, // null or bitmask of reliable reference tiles
						quadCLTs,           // QuadCLT []     quadCLTs,
						debugLevel) ;       // int            debugLevel)
				// should update earliest_scene
				if ((ref_index - earliest_scene + 1) < min_num_scenes) {
					System.out.println("After reAdjustPairsLMAInterscene() total number of useful scenes = "+(ref_index - earliest_scene + 1)+
							" < "+min_num_scenes+". Scrapping this series.");
					if (start_ref_pointers != null) {
						start_ref_pointers[0] = earliest_scene;
					}
					return null;
				}
				if (earliest_scene > 0) {
					System.out.println("After reAdjustPairsLMAInterscene() not all scenes matched, earliest useful scene = "+earliest_scene+
							" (total number of scenes = "+(ref_index - earliest_scene + 1)+").");
				}
				quadCLTs[ref_index].inc_orient(); 
				quadCLTs[ref_index].saveInterProperties( // save properties for interscene processing (extrinsics, ers, ...)
						null, // String path,             // full name with extension or w/o path to use x3d directory
						debugLevel+1);
			}
		}

		if (test_ers) { // only debug feature
			test_ers0 = quadCLTs.length -1; // make it always == reference !
			
			testERS(
					clt_parameters, // CLTParameters clt_parameters,
					test_ers0, // int           indx0, // reference scene in a pair
					test_ers1, // int           indx1, // other scene in a pair
					//		    		double []     ref_disparity,			
					quadCLTs, // QuadCLT []    quadCLTs,
					debugLevel); // int           debugLevel)

			System.out.println("buildSeries(): ABORTED after test_ers"); //
			return quadCLTs[ref_index].getX3dTopDirectory();
		}
		
		// generates 3-d modes, colors, stereos, tiffs/videos
		
		if (generate_mapped || reuse_video) {
			int tilesX =  quadCLTs[ref_index].getTileProcessor().getTilesX();
	        int tilesY =  quadCLTs[ref_index].getTileProcessor().getTilesY();
	        double [] disparity_fg =  null;
	        double [] strength_fg =   null;
	        double [] disparity_bg =  null;
	        double [] strength_bg =   null;
	        double [] disparity_raw = null;
	        if (generate_mapped) {
	        	disparity_raw = new double [tilesX * tilesY];
	        	Arrays.fill(disparity_raw,clt_parameters.disparity);
	        	if (combo_dsn_final == null) {
	        		combo_dsn_final = quadCLTs[ref_index].readDoubleArrayFromModelDirectory(
	        				"-INTER-INTRA-LMA", // String      suffix,
	        				0, // int         num_slices, // (0 - all)
	        				null); // int []      wh);
	        	}
	        	double [][] dls = {
	        			combo_dsn_final[COMBO_DSN_INDX_DISP],
	        			combo_dsn_final[COMBO_DSN_INDX_LMA],
	        			combo_dsn_final[COMBO_DSN_INDX_STRENGTH]
	        	};
	        	double [][] ds = conditionInitialDS(
	        			true,                // boolean        use_conf,       // use configuration parameters, false - use following  
	        			clt_parameters,      // CLTParameters  clt_parameters,
	        			dls,                 // double [][]    dls
	        			quadCLTs[ref_index], // QuadCLT        scene,
	        			debugLevel);			

	        	disparity_fg = ds[0]; // combo_dsn_final[COMBO_DSN_INDX_DISP_FG];
	        	strength_fg =  ds[1];
	        	// BG mode
	        	double [] bg_lma = combo_dsn_final[COMBO_DSN_INDX_DISP_BG_ALL].clone();
	        	double [] bg_str = combo_dsn_final[COMBO_DSN_INDX_STRENGTH].clone();

	        	for (int i = 0; i < bg_lma.length; i++) {
	        		if (Double.isNaN(combo_dsn_final[COMBO_DSN_INDX_LMA][i])){
	        			bg_lma[i] = Double.NaN;
	        		}
	        		if (!Double.isNaN(combo_dsn_final[COMBO_DSN_INDX_DISP_BG][i])){
	        			bg_lma[i] = combo_dsn_final[COMBO_DSN_INDX_LMA_BG][i];
	        		}
	        		if (!Double.isNaN(combo_dsn_final[COMBO_DSN_INDX_STRENGTH_BG][i])){
	        			bg_str[i] = combo_dsn_final[COMBO_DSN_INDX_STRENGTH_BG][i];
	        		}
	        	}
	        	double [][] dls_bg = {
	        			combo_dsn_final[COMBO_DSN_INDX_DISP_BG_ALL],
	        			bg_lma,
	        			bg_str
	        	};
	        	double [][] ds_bg = conditionInitialDS(
	        			true, // boolean        use_conf,       // use configuration parameters, false - use following  
	        			clt_parameters,         // CLTParameters  clt_parameters,
	        			dls_bg,                 // double [][]    dls
	        			quadCLTs[ref_index],    // QuadCLT        scene,
	        			debugLevel);			
	        	disparity_bg = ds_bg[0]; // combo_dsn_final[COMBO_DSN_INDX_DISP_FG];
	        	strength_bg =  ds_bg[1];
	        	// for now using disparity for just standard size (90x64), later may use full size and at
	        	// minimum fill peripheral areas with Laplassian?
	        	double [][] dxyzatr_dt = new double [quadCLTs.length][];
	        	int [][] min_max_vel = getERSStats(
	        			clt_parameters,   //  CLTParameters clt_parameters,
	        			quadCLTs,         // QuadCLT []    quadCLTs,
	        			dxyzatr_dt);      // double [][] dxyzatr_dt);
	        	double [] vel_avg = new double [min_max_vel.length];
	        	int navg = 0;
	        	for (int n = 0; n < dxyzatr_dt.length; n++) if (dxyzatr_dt[n] != null){
	        		for (int i = 0; i < vel_avg.length; i++) {
	        			vel_avg[i] += dxyzatr_dt[n][i];
	        		}
	        		navg++;
	        	}
	        	for (int i = 0; i < vel_avg.length; i++) {
	        		vel_avg[i] /= navg;
	        	}

	        	String [] vel_names = {"Vx","Vy","Vz","Vaz","Vtl","Vrl"};
	        	String [] vel_units = {"m/s","m/s","m/s","mrad/s","mrad/s","mrad/s","mrad/s"};
	        	System.out.println("Minimal/maximal/average linear and angular velocities of the scenes (ref. scene: "+
	        			quadCLTs[quadCLTs.length - 1].getImageName()+")");
	        	for (int i = 0; i < min_max_vel.length; i++) {
	        		double v0 = dxyzatr_dt[min_max_vel[i][0]][i];
	        		double v1 = dxyzatr_dt[min_max_vel[i][1]][i];
	        		double va = vel_avg[i];
	        		if (i == 2) { // forward movement
	        			v0*=-1;
	        			v1*=-1;
	        			va*=-1;
	        		}
	        		if (i > 2) { // angular
	        			v0*=1000;
	        			v1*=1000;
	        			va*=1000;
	        		}
	        		System.out.println(String.format(
	        				"%3s: min: %3d (%8.3f%-6s), max: %3d (%8.3f%-6s), average: %8.3f%-6s",
	        				vel_names[i], min_max_vel[i][0], v0, vel_units[i],
	        				min_max_vel[i][1], v1, vel_units[i],
	        				va, vel_units[i]));
	        	}
	        }
	        
	        // for older compatibility mode3d = -1 for RAW, 0 - INF, 1 - FG, 2 BG
	        for (int mode3d = -1; mode3d < (generate_modes3d.length-1); mode3d++) if (generate_modes3d[mode3d+1]) {
	        	Rectangle fov_tiles = new Rectangle(
	        			extra_hor_tile,
	        			extra_vert_tile,
	        			tilesX + 2 * extra_hor_tile,
	        			tilesY + 2 * extra_vert_tile);
	        	if ((mode3d < 0) || (crop_3d && (mode3d > 0))) {
	        		fov_tiles = null; // use sensor dimensions
	        	}
	        	boolean is_3d = mode3d > 0;
	        	boolean gen_stereo = is_3d && generate_stereo;
//	        	double [] baselines = (gen_stereo)? stereo_bases : new double[] {0.0};
	        	double [][]  views = (gen_stereo)? stereo_views : new double [][] {{0.0,0.0,0.0}};
	        	
//	        	for (int ibase = 0; ibase < baselines.length; ibase++) if (!gen_stereo || generate_stereo_var[ibase]) {
	        	for (int ibase = 0; ibase < views.length; ibase++) if (!gen_stereo || generate_stereo_var[ibase]) {
	        		//    				double stereo_baseline = gen_stereo? stereo_bases[ibase] : 0.0;
	        		double stereo_baseline = gen_stereo? views[ibase][0] : 0.0;
	        		boolean is_stereo = gen_stereo && stereo_baseline > 0;
	        		double  stereo_baseline_meters =    0.001 * stereo_baseline;
	        		double  view_height_meters =        0.001 * views[ibase][1];
	        		double  view_back_meters =          0.001 * views[ibase][2];
	        		//    				double  stereo_back = 3.0; // 0; // -10.0; // meters
	        		// col_mode: 0 - mono, 1 - color
	        		for (int col_mode = 0; col_mode < 2; col_mode++) if (gen_seq_mono_color[col_mode]){ // skip if not needed
	        			double[] selected_disparity = (mode3d > 1)?disparity_bg:((mode3d > 0)?disparity_fg: disparity_raw); 
	        			double[] selected_strength =  (mode3d > 1)?strength_bg:((mode3d > 0)?strength_fg: null);
	        			if (selected_strength != null) { // for FG/BG only, fixing for transformCameraVew()
	        				for (int i = 0; i < selected_disparity.length; i++) {
	        					if (!Double.isNaN(selected_disparity[i]) && (selected_strength[i] == 0)) selected_strength[i] = 0.01; // transformCameraVew ignores strength= 0
	        				}
	        			}
	        			final boolean       toRGB = col_mode > 0;
	        			String scenes_suffix = quadCLTs[quadCLTs.length-1].getImageName()+
	        					"-SEQ-" + IntersceneMatchParameters.MODES3D[mode3d+1] + "-"+(toRGB?"COLOR":"MONO");
	        			if (!toRGB && um_mono) {
	        				if (mono_fixed) {
		        				scenes_suffix+=String.format("-UM%.1f_%.2f_%.0f",um_sigma,um_weight,mono_range);
	        				} else {
		        				scenes_suffix+=String.format("-UM%.1f_%.2f_A",um_sigma,um_weight);
	        				}
	        			}
	        			int num_stereo = (is_stereo && (mode3d > 0))? 2:1; // only for 3D views
	        			boolean combine_left_right = (num_stereo > 1) && (stereo_merge || (anaglyth_en && !toRGB));
	        			ImagePlus [] imp_scenes_pair = new ImagePlus[num_stereo];
	        			String scenes_suffix_pair = scenes_suffix;
	        			for (int nstereo = 0; nstereo < num_stereo; nstereo++) {
	        				double [] xyz_offset = {
	        						-stereo_baseline_meters * (nstereo - 0.5) * (num_stereo - 1), // x offset
	        						-view_height_meters,  // Y offset
	        						-view_back_meters};   // Z offset
	        				if (num_stereo > 1) {
	        					scenes_suffix = scenes_suffix_pair + ((nstereo > 0)?"-RIGHT":"-LEFT"); // check if opposite
	        					scenes_suffix += "-B"+String.format("%.0f",views[ibase][0]);
	        				}
	        				if (views[ibase][1] != 0) {
	        					scenes_suffix += "-Y"+String.format("%.0f",views[ibase][1]);
	        				}
	        				if (views[ibase][2] != 0) {
	        					scenes_suffix += "-Z"+String.format("%.0f",views[ibase][2]);
	        				}
	        				double [][] ds_vantage = new double[][] {selected_disparity,selected_strength};
	        				if ((views[ibase][0] != 0) || (views[ibase][1] != 0) || (views[ibase][2] != 0)) {
	        					ds_vantage = transformCameraVew(
	        							null, // (debug_ds_fg_virt?"transformCameraVew":null),     // final String    title,
	        							ds_vantage,                    // final double [][] dsrbg_camera_in,
	        							xyz_offset, // _inverse[0], // final double [] scene_xyz, // camera center in world coordinates
	        							ZERO3, // _inverse[1], // final double [] scene_atr, // camera orientation relative to world frame
	        							quadCLTs[ref_index],      // final QuadCLT   scene_QuadClt,
	        							quadCLTs[ref_index],      // final QuadCLT   reference_QuadClt,
	        							8); // iscale);                  // final int       iscale);
	        				}
	        				if (generate_mapped) {
	        					imp_scenes_pair[nstereo]= renderSceneSequence(
	        							clt_parameters,     // CLTParameters clt_parameters,
	        							fov_tiles,          // Rectangle     fov_tiles,
	        							mode3d,             // int           mode3d,
	        							toRGB,              // boolean       toRGB,
	        							xyz_offset,         // double []     stereo_offset, // offset reference camera {x,y,z}
	        							sensor_mask,        // int           sensor_mask,
	        							scenes_suffix,      // String        suffix,
	        							ds_vantage[0],      // selected_disparity, // double []     ref_disparity,			
	        							quadCLTs,           // QuadCLT []    quadCLTs,
	        							debugLevel);        // int           debugLevel);
	        					if (save_mapped_mono_color[col_mode]) {	        	
	        						quadCLTs[ref_index].saveImagePlusInModelDirectory(
	        								null, // "GPU-SHIFTED-D"+clt_parameters.disparity, // String      suffix,
	        								imp_scenes_pair[nstereo]); // imp_scenes); // ImagePlus   imp)
	        					}
	        				} else {
	        					if (fov_tiles==null) {
	        						fov_tiles = new Rectangle(0, 0, tilesX, tilesY);
	        					}
	        					FloatProcessor fp = new FloatProcessor(
	        							fov_tiles.width*quadCLTs[ref_index].getTileProcessor().getTileSize(),
	        							fov_tiles.height*quadCLTs[ref_index].getTileProcessor().getTileSize());

	        					boolean merge_all = clt_parameters.imp.merge_all;
	        					if (mode3d < 1) {
	        						merge_all = false;
	        					}
	        					imp_scenes_pair[nstereo]=new ImagePlus(scenes_suffix+((mode3d > 0)?(merge_all?"-MERGED":"-SINGLE"):""), fp);
	        				}
	        				// Save as AVI	        	
	        				if (gen_avi_mono_color[col_mode]) {
	        					if (annotate_mono_color[col_mode] && generate_mapped) {
	        						if (!toRGB) {
	        							// If it is mono, first convert to color
	        							if (mono_fixed && um_mono) {
	        								imp_scenes_pair[nstereo].getProcessor().setMinAndMax(-mono_range/2, mono_range/2);
	        							}
	        							ImageConverter imageConverter = new ImageConverter(imp_scenes_pair[nstereo]);
	        							imageConverter.convertToRGB(); // Did it convert imp_scenes ?
	        						}

	        						final Color fcolor = toRGB ? annotate_color_color: annotate_color_mono;
	        						final ImageStack fstack_scenes = imp_scenes_pair[nstereo].getImageStack();
	        						final int width =  imp_scenes_pair[nstereo].getWidth();
	        						final int height = imp_scenes_pair[nstereo].getHeight();
	        						final int posX= width - 119; // 521;
	        						final int posY= height + 1;  // 513;
	        						final Font font = new Font("Monospaced", Font.PLAIN, 12);
	        						final int nSlices = fstack_scenes.getSize();
	        						final Thread[] threads = ImageDtt.newThreadArray(QuadCLT.THREADS_MAX);
	        						final AtomicInteger ai = new AtomicInteger(0);
	        						for (int ithread = 0; ithread < threads.length; ithread++) {
	        							threads[ithread] = new Thread() {
	        								public void run() {
	        									for (int nSlice = ai.getAndIncrement(); nSlice < nSlices; nSlice = ai.getAndIncrement()) {
	        										String scene_title = fstack_scenes.getSliceLabel(nSlice+1);
	        										ImageProcessor ip = fstack_scenes.getProcessor(nSlice+1);
	        										ip.setColor(fcolor); // Color.BLUE);
	        										ip.setFont(font);
	        										if (toRGB || !annotate_transparent_mono) {
	        											ip.drawString(scene_title, posX, posY,Color.BLACK);
	        										} else {
	        											ip.drawString(scene_title, posX, posY); // transparent
	        										}
	        									}
	        								}
	        							};
	        						}		      
	        						ImageDtt.startAndJoin(threads);
	        					}
	        					if (combine_left_right && (nstereo == 0)) {
	        						continue;
	        					}
	        					// no_combine, stereo_2_images, stereo_anaglyth
	        					ImagePlus imp_video = imp_scenes_pair[nstereo];
	        					boolean [] combine_modes = {!combine_left_right, stereo_merge && combine_left_right, anaglyth_en && !toRGB && combine_left_right };
	        					for (int istereo_mode = 0; istereo_mode < combine_modes.length; istereo_mode++) if(combine_modes[istereo_mode]) {
//	        						if (combine_left_right) { // combine pairs multi-threaded
	        						if (istereo_mode == 1) { // combine pairs for "Google" VR 
	        							final int left_width = imp_scenes_pair[0].getWidth();
	        							final int right_width = imp_scenes_pair[1].getWidth();
	        							final int stereo_width =  left_width + right_width+stereo_gap;
	        							final int stereo_height = imp_scenes_pair[0].getHeight();
	        							final ImageStack stereo_stack = new ImageStack(stereo_width, stereo_height);
	        							final int nSlices = imp_scenes_pair[0].getStack().getSize();
	        							for (int i = 0; i < nSlices; i++) {
	        								stereo_stack.addSlice(
	        										imp_scenes_pair[0].getStack().getSliceLabel(i+1),
	        										new int[stereo_width * stereo_height]);
	        							}
	        							if (generate_mapped) {
	        								final Thread[] threads = ImageDtt.newThreadArray(QuadCLT.THREADS_MAX);
	        								final AtomicInteger ai = new AtomicInteger(0);
	        								for (int ithread = 0; ithread < threads.length; ithread++) {
	        									threads[ithread] = new Thread() {
	        										public void run() {
	        											for (int nSlice = ai.getAndIncrement(); nSlice < nSlices; nSlice = ai.getAndIncrement()) {
	        												int[] pixels_stereo = (int[]) stereo_stack.getPixels(nSlice+1);
	        												int[] pixels_left = (int[]) imp_scenes_pair[0].getStack().getPixels(nSlice+1);
	        												int[] pixels_right = (int[]) imp_scenes_pair[1].getStack().getPixels(nSlice+1);
	        												for (int row = 0; row < stereo_height; row++) {
	        													System.arraycopy(
	        															pixels_left,
	        															left_width * row,
	        															pixels_stereo,
	        															stereo_width * row,
	        															left_width);
	        													System.arraycopy(
	        															pixels_right,
	        															right_width * row,
	        															pixels_stereo,
	        															stereo_width * row + left_width + stereo_gap,
	        															right_width);
	        												}
	        											}
	        										}
	        									};
	        								}		      
	        								ImageDtt.startAndJoin(threads);
	        							}
	        							imp_video = new ImagePlus();
	        							imp_video.setImage(imp_scenes_pair[1]); // copy many attributes
	        							imp_video.setStack(stereo_stack);
	        							String title = imp_scenes_pair[1].getTitle();
	        							imp_video.setTitle(title.replace("-RIGHT","-STEREO"));

	        							// convert stereo_stack to imp_scenes_pair[1], keeping calibration and fps?
///	        							imp_scenes_pair[1].setStack(stereo_stack);
///	        							String title = imp_scenes_pair[1].getTitle();
///	        							imp_video = new ImagePlus(
///	        									imp_scenes_pair[1].getTitle().replace("-RIGHT","-STEREO"),
///	        									stereo_stack);
///	        							imp_scenes_pair[1].setTitle(title.replace("-RIGHT","-STEREO"));
	        						} else if (istereo_mode == 2) { // combine anaglyph
//	        							final Color anaglyph_left =      clt_parameters.imp.anaglyph_left;
//	        							final Color anaglyph_right =     clt_parameters.imp.anaglyph_right;
	        							final double [] left_rgb= {
	        									anaglyph_left.getRed()/255.0,
	        									anaglyph_left.getGreen()/255.0,
	        									anaglyph_left.getBlue()/255.0};  
	        							final double [] right_rgb= {
	        									anaglyph_right.getRed()/255.0,
	        									anaglyph_right.getGreen()/255.0,
	        									anaglyph_right.getBlue()/255.0};  
	        							final int left_width =  imp_scenes_pair[0].getWidth();
	        							final int left_height = imp_scenes_pair[0].getHeight();
	        							final int nSlices = imp_scenes_pair[0].getStack().getSize();
	        							final ImageStack stereo_stack = new ImageStack(left_width, left_height);
	        							for (int i = 0; i < nSlices; i++) {
	        								stereo_stack.addSlice(
	        										imp_scenes_pair[0].getStack().getSliceLabel(i+1),
	        										new int[left_width * left_height]);
	        							}
	        							
	        							if (generate_mapped) {
	        								final Thread[] threads = ImageDtt.newThreadArray(QuadCLT.THREADS_MAX);
	        								final AtomicInteger ai = new AtomicInteger(0);
	        								for (int ithread = 0; ithread < threads.length; ithread++) {
	        									threads[ithread] = new Thread() {
	        										public void run() {
    													int [] rgb = new int [3]; 
	        											for (int nSlice = ai.getAndIncrement(); nSlice < nSlices; nSlice = ai.getAndIncrement()) {
	        												int[] pixels_stereo = (int[]) stereo_stack.getPixels(nSlice+1);
	        												int[] pixels_left =   (int[]) imp_scenes_pair[0].getStack().getPixels(nSlice+1);
	        												int[] pixels_right =  (int[]) imp_scenes_pair[1].getStack().getPixels(nSlice+1);
	        												for (int pix = 0; pix < pixels_left.length; pix++) {
	        													int gl = ((pixels_left[pix] & 0xff00) >> 8);
	        													int gr = ((pixels_right[pix] & 0xff00) >> 8);
	        													rgb[0] = ((int) Math.min(gl * left_rgb[0] + gr * right_rgb[0], 255)) & 0xff;
	        													rgb[1] = ((int) Math.min(gl * left_rgb[1] + gr * right_rgb[1], 255)) & 0xff;
	        													rgb[2] = ((int) Math.min(gl * left_rgb[2] + gr * right_rgb[2], 255)) & 0xff;
	        													pixels_stereo[pix] = 0xff000000 + (rgb[0] << 16) + (rgb[1] << 8)  + rgb[2]; 
	        												}
	        											}
	        										}
	        									};
	        								}		      
	        								ImageDtt.startAndJoin(threads);
	        							}
	        							imp_video = new ImagePlus();
	        							imp_video.setImage(imp_scenes_pair[1]); // copy many attributes
	        							imp_video.setStack(stereo_stack);
	        							String title = imp_scenes_pair[1].getTitle();
	        							imp_video.setTitle(title.replace("-RIGHT","-ANAGLYPH"));
///	        							String title = imp_scenes_pair[1].getTitle();
///	        							imp_scenes_pair[1].setTitle(title.replace("-RIGHT","-ANAGLYPH"));
	        						} // if (istereo_mode == 1) {if (combine_left_right) { // combine pairs multi-threaded
	        						String avi_path=null;
	        						video:
	        						{
	        							try {
	        								avi_path=quadCLTs[ref_index].saveAVIInModelDirectory(
	        										!generate_mapped, // boolean     dry_run, 
	        										null,             // "GPU-SHIFTED-D"+clt_parameters.disparity, // String      suffix,
	        										mode_avi,         // int         avi_mode,
	        										avi_JPEG_quality, // int         avi_JPEG_quality,
	        										video_fps,        // double      fps,
	        										imp_video); // imp_scenes_pair[nstereo]);      // ImagePlus   imp)
	        							} catch (IOException e) {
	        								// TODO Auto-generated catch block
	        								e.printStackTrace();
	        								break video;

	        							}
	        							// Convert with ffmpeg?
	        							if (avi_path == null) {
	        								break video;
	        							}
//	        							int img_width=imp_scenes_pair[nstereo].getWidth();
	        							int img_width=imp_video.getWidth();
	        							int stereo_width = combine_left_right? img_width:0;
	        							stereo_widths_list.add(stereo_width);
	        							if (!run_ffmpeg) {
	        								video_list.add(avi_path);
	        								break video; // webm not requested
	        							}

	        							String webm_path = avi_path.substring(0, avi_path.length()-4)+video_ext;
	        							// added -y not to as "overwrite y/n?"
	        							String shellCommand = String.format("ffmpeg -y -i %s -c %s -b:v 0 -crf %d %s",
	        									avi_path, video_codec, video_crf, webm_path);
	        							Process p = null;
	        							if (generate_mapped) {
	        								int exit_code = -1;
	        								try {
	        									p = Runtime.getRuntime().exec(shellCommand);
	        								} catch (IOException e) {
	        									System.out.println("Failed shell command: \""+shellCommand+"\"");
	        								}
	        								if (p != null) {
	        									p.waitFor();
	        									exit_code = p.exitValue();
	        								}
	        								System.out.println("Ran shell command: \""+shellCommand+"\" -> "+exit_code);
	        								// Check if webm file exists
	        								if ((exit_code != 0) || !(new File(webm_path)).exists()) {
	        									System.out.println("Failed to create : \""+webm_path+"\"");
	        									video_list.add(avi_path);
	        									break video;
	        								}
	        							} else {
	        								System.out.println("Simulated shell command: \""+shellCommand);

	        							}
	        							video_list.add(webm_path);
	        							if (remove_avi && generate_mapped) {
	        								(new File(avi_path)).delete();
	        								System.out.println("Deleted AVI video file: \""+avi_path+"\"");
	        							}
	        						}

	        					} // for (int istereo_mode = 0; istereo_mode < stereo_modes.length; istereo_mode++) if(combine_modes[istereo_mode]) {
	        				} // if (gen_avi_mono_color[col_mode])
	        				if (show_mono_color[col_mode]  && generate_mapped) {
	        					imp_scenes_pair[nstereo].show();
	        				}
	        			} // for (int nstereo = 0; nstereo < num_stereo; nstereo++)
	        		} // for (int col_mode = 0; col_mode<2; col_mode++) {
	        	} // for (int ibase = 0; ibase < baselines.length; ibase++) {
	        }// for (int mode3d = 0; mode3d<generate_modes3d.length; mode3d++) if (generate_modes3d[mode3d]) {
		} // if (generate_mapped) {
		
		if (export_images) {
			if (combo_dsn_final == null) {
				combo_dsn_final = quadCLTs[ref_index].readDoubleArrayFromModelDirectory(
						"-INTER-INTRA-LMA", // String      suffix,
						0, // int         num_slices, // (0 - all)
						null); // int []      wh);
			}
			double [][] dls = {
					combo_dsn_final[COMBO_DSN_INDX_DISP],
					combo_dsn_final[COMBO_DSN_INDX_LMA],
					combo_dsn_final[COMBO_DSN_INDX_STRENGTH]
			};
			double [][] ds_fg = conditionInitialDS(
					true, // boolean        use_conf,       // use configuration parameters, false - use following  
					clt_parameters,      // CLTParameters  clt_parameters,
					dls,                 // double [][]    dls
					quadCLTs[ref_index], // QuadCLT        scene,
					debugLevel);			
			
			double [] fg_disparity = ds_fg[0]; // combo_dsn_final[COMBO_DSN_INDX_DISP_FG];
			
			double [] bg_lma = combo_dsn_final[COMBO_DSN_INDX_DISP_BG_ALL].clone();
			double [] bg_str = combo_dsn_final[COMBO_DSN_INDX_STRENGTH].clone();
			
			for (int i = 0; i < bg_lma.length; i++) {
				if (Double.isNaN(combo_dsn_final[COMBO_DSN_INDX_LMA][i])){
					bg_lma[i] = Double.NaN;
				}
				if (!Double.isNaN(combo_dsn_final[COMBO_DSN_INDX_DISP_BG][i])){
					bg_lma[i] = combo_dsn_final[COMBO_DSN_INDX_LMA_BG][i];
				}
				if (!Double.isNaN(combo_dsn_final[COMBO_DSN_INDX_STRENGTH_BG][i])){
					bg_str[i] = combo_dsn_final[COMBO_DSN_INDX_STRENGTH_BG][i];
				}
			}
			double [][] dls_bg = {
					combo_dsn_final[COMBO_DSN_INDX_DISP_BG_ALL],
					bg_lma,
					bg_str
			};
			double [][] ds_bg = conditionInitialDS(
					true, // boolean        use_conf,       // use configuration parameters, false - use following  
					clt_parameters,      // CLTParameters  clt_parameters,
					dls_bg,                 // double [][]    dls
					quadCLTs[ref_index], // QuadCLT        scene,
					debugLevel);			
			double [] bg_disparity = ds_bg[0]; // combo_dsn_final[COMBO_DSN_INDX_DISP_FG];


			double [] constant_disparity = new double [fg_disparity.length];
			Arrays.fill(constant_disparity,clt_parameters.disparity);
			//Rectangle testr = new Rectangle(10, 8, 100,80);
			ImagePlus imp_constant = QuadCLT.renderGPUFromDSI(
					-1,                  // final int         sensor_mask,
					false, // final boolean     merge_channels,
					null, // testr, // null,                // final Rectangle   full_woi_in,      // show larger than sensor WOI (or null)
					clt_parameters,      // CLTParameters     clt_parameters,
					constant_disparity,  // double []         disparity_ref,
					ZERO3,               // final double []   scene_xyz, // camera center in world coordinates
					ZERO3, // new double[] {.1,0.1,.1}, // ZERO3,               // final double []   scene_atr, // camera orientation relative to world frame
					quadCLTs[ref_index], // final QuadCLT     scene,
					quadCLTs[ref_index], // final QuadCLT     ref_scene, // now - may be null - for testing if scene is rotated ref
					true, // toRGB,               // final boolean     toRGB,
					"GPU-SHIFTED-D"+clt_parameters.disparity, // String            suffix,
					threadsMax,          // int               threadsMax,
					debugLevel);         // int         debugLevel)
			quadCLTs[ref_index].saveImagePlusInModelDirectory(
					null, // "GPU-SHIFTED-D"+clt_parameters.disparity, // String      suffix,
					imp_constant); // ImagePlus   imp)
			ImagePlus imp_constant_mono = QuadCLT.renderGPUFromDSI(
					-1,                  // final int         sensor_mask,
					false, // final boolean     merge_channels,
					null, // testr, // null,                // final Rectangle   full_woi_in,      // show larger than sensor WOI (or null)
					clt_parameters,      // CLTParameters     clt_parameters,
					constant_disparity,  // double []         disparity_ref,
					ZERO3,               // final double []   scene_xyz, // camera center in world coordinates
					ZERO3,               // final double []   scene_atr, // camera orientation relative to world frame
					quadCLTs[ref_index], // final QuadCLT     scene,
					quadCLTs[ref_index], // final QuadCLT     ref_scene, // now - may be null - for testing if scene is rotated ref
					false, // toRGB,               // final boolean     toRGB,
					"GPU-SHIFTED-D"+clt_parameters.disparity, // String            suffix,
					threadsMax,          // int               threadsMax,
					debugLevel);         // int         debugLevel)
			quadCLTs[ref_index].saveImagePlusInModelDirectory(
					null, // "GPU-SHIFTED-D"+clt_parameters.disparity, // String      suffix,
					imp_constant_mono); // ImagePlus   imp)
			if (show_images) {
				imp_constant.show();
				if (show_images_mono) {
					imp_constant_mono.show();
				}
			}
			
			boolean offset_fg_image = false; // true; // config later, generate FG image for all stereo views
			double [][] img_views = offset_fg_image ? stereo_views : (new double [][] {{0,0,0}});
			double min_str = 0.01; 
			for (int i = 0; i < ds_fg[0].length; i++) {
				if (!Double.isNaN(ds_fg[0][i]) && (ds_fg[1][i] == 0)) ds_fg[1][i] = min_str; // transformCameraVew ignores strength= 0
				if (!Double.isNaN(ds_bg[0][i]) && (ds_bg[1][i] == 0)) ds_bg[1][i] = min_str;
			}
			for (int ibase = 0; ibase < img_views.length; ibase++) if (!offset_fg_image || generate_stereo_var[ibase]) {
        		double  stereo_baseline_meters =    0.001 * img_views[ibase][0];
        		double  view_height_meters =        0.001 * img_views[ibase][1];
        		double  view_back_meters =          0.001 * img_views[ibase][2];
				double [] xyz_offset = {
						-stereo_baseline_meters, // x offset
						-view_height_meters,     // Y offset
						-view_back_meters};      // Z offset
				double [] atr_offset = ZERO3; 
        		String scenes_suffix = "";
				if (img_views[ibase][0] != 0) {
					scenes_suffix += "-B"+String.format("%.0f",img_views[ibase][0]);
					
				}
				if (img_views[ibase][1] != 0) {
					scenes_suffix += "-Y"+String.format("%.0f",img_views[ibase][1]);
					
				}
				if (img_views[ibase][2] != 0) {
					scenes_suffix += "-Z"+String.format("%.0f",img_views[ibase][2]);

				}
        		// calculate virtual view fg_ds_virt from the reference ds_fg;
				boolean debug_ds_fg_virt = false; // false;
				double [][] ds_fg_virt = ds_fg;
				if ((img_views[ibase][0] != 0) || (img_views[ibase][1] != 0) || (img_views[ibase][2] != 0)) {
					ds_fg_virt = transformCameraVew(
							(debug_ds_fg_virt?"transformCameraVew":null),     // final String    title,
							ds_fg,                    // final double [][] dsrbg_camera_in,
							xyz_offset, // _inverse[0], // final double [] scene_xyz, // camera center in world coordinates
							atr_offset, // _inverse[1], // final double [] scene_atr, // camera orientation relative to world frame
							quadCLTs[ref_index],      // final QuadCLT   scene_QuadClt,
							quadCLTs[ref_index],      // final QuadCLT   reference_QuadClt,
							8); // iscale);                  // final int       iscale);
					if (debug_ds_fg_virt){
						int dbgX =  quadCLTs[ref_index].getTileProcessor().getTilesX();
						int dbgY =  quadCLTs[ref_index].getTileProcessor().getTilesY();
						double [][] dbg_img = new double[][] {ds_fg[0],ds_fg_virt[0],ds_fg[1],ds_fg_virt[1]};
						(new ShowDoubleFloatArrays()).showArrays(
								dbg_img,
								dbgX,
								dbgY,
								true,
								"virtual-view-ds",
								new String[] {"d-ref","d-virt","s-ref","s-virt"}); //	dsrbg_titles);
					}
				}
				
				ImagePlus imp_fg = QuadCLT.renderGPUFromDSI(
						-1,                  // final int         sensor_mask,
						false,               // final boolean     merge_channels,
						null,                // final Rectangle   full_woi_in,      // show larger than sensor WOI (or null)
						clt_parameters,      // CLTParameters     clt_parameters,
						ds_fg_virt[0],       // fg_disparity,  // double []         disparity_ref,
						xyz_offset,          // ZERO3,               // final double []   scene_xyz, // camera center in world coordinates
						ZERO3,               // final double []   scene_atr, // camera orientation relative to world frame
						quadCLTs[ref_index], // final QuadCLT     scene,
						quadCLTs[ref_index], // final QuadCLT     ref_scene, // now - may be null - for testing if scene is rotated ref
						true, // toRGB,               // final boolean     toRGB,
						scenes_suffix+"GPU-SHIFTED-FOREGROUND", // String            suffix,
						threadsMax,          // int               threadsMax,
						debugLevel);         // int         debugLevel)
				quadCLTs[ref_index].saveImagePlusInModelDirectory(
						null, // "GPU-SHIFTED-FOREGROUND", // String      suffix,
						imp_fg); // ImagePlus   imp)
				ImagePlus imp_fg_mono = QuadCLT.renderGPUFromDSI(
						-1,                  // final int         sensor_mask,
						false, // final boolean     merge_channels,
						null,                // final Rectangle   full_woi_in,      // show larger than sensor WOI (or null)
						clt_parameters,      // CLTParameters     clt_parameters,
						ds_fg_virt[0], // fg_disparity,  // double []         disparity_ref,
						xyz_offset, // ZERO3,               // final double []   scene_xyz, // camera center in world coordinates
						ZERO3,               // final double []   scene_atr, // camera orientation relative to world frame
						quadCLTs[ref_index], // final QuadCLT     scene,
						quadCLTs[ref_index], // final QuadCLT     ref_scene, // now - may be null - for testing if scene is rotated ref
						false, // toRGB,               // final boolean     toRGB,
						scenes_suffix+"GPU-SHIFTED-FOREGROUND", // String            suffix,
						threadsMax,          // int               threadsMax,
						debugLevel);         // int         debugLevel)
				quadCLTs[ref_index].saveImagePlusInModelDirectory(
						null, // "GPU-SHIFTED-FOREGROUND", // String      suffix,
						imp_fg_mono); // ImagePlus   imp)
				if (show_images && show_images_bgfg) {
					imp_fg.show();
					if (show_images_mono) {
						imp_fg_mono.show();
					}
				}
			}

			
			ImagePlus imp_bg = QuadCLT.renderGPUFromDSI(
					-1,                  // final int         sensor_mask,
					false, // final boolean     merge_channels,
					null,                // final Rectangle   full_woi_in,      // show larger than sensor WOI (or null)
					clt_parameters,      // CLTParameters     clt_parameters,
					bg_disparity,        // double []         disparity_ref,
					ZERO3,               // final double []   scene_xyz, // camera center in world coordinates
					ZERO3,               // final double []   scene_atr, // camera orientation relative to world frame
					quadCLTs[ref_index], // final QuadCLT     scene,
					quadCLTs[ref_index], // final QuadCLT     ref_scene, // now - may be null - for testing if scene is rotated ref
					true,               // final boolean     toRGB,
					"GPU-SHIFTED-BACKGROUND", // String            suffix,
					threadsMax,          // int               threadsMax,
					debugLevel);         // int         debugLevel)
			quadCLTs[ref_index].saveImagePlusInModelDirectory(
					null, // "GPU-SHIFTED-BACKGROUND", // String      suffix,
					imp_bg); // ImagePlus   imp)
			ImagePlus imp_bg_mono = QuadCLT.renderGPUFromDSI(
					-1,                  // final int         sensor_mask,
					false, // final boolean     merge_channels,
					null,                // final Rectangle   full_woi_in,      // show larger than sensor WOI (or null)
					clt_parameters,      // CLTParameters     clt_parameters,
					bg_disparity,        // double []         disparity_ref,
					ZERO3,               // final double []   scene_xyz, // camera center in world coordinates
					ZERO3,               // final double []   scene_atr, // camera orientation relative to world frame
					quadCLTs[ref_index], // final QuadCLT     scene,
					quadCLTs[ref_index], // final QuadCLT     ref_scene, // now - may be null - for testing if scene is rotated ref
					false,               // final boolean     toRGB,
					"GPU-SHIFTED-BACKGROUND", // String            suffix,
					threadsMax,          // int               threadsMax,
					debugLevel);         // int         debugLevel)
			quadCLTs[ref_index].saveImagePlusInModelDirectory(
					null, // "GPU-SHIFTED-BACKGROUND", // String      suffix,
					imp_bg_mono); // ImagePlus   imp)
			if (show_images && show_images_bgfg) {
				imp_bg.show();
				if (show_images_mono) {
					imp_bg_mono.show();
				}
			}
		}
		
		if (export_dsi_image || show_dsi_image) {
			if (combo_dsn_final == null) {
				combo_dsn_final =quadCLTs[ref_index].readDoubleArrayFromModelDirectory(
						"-INTER-INTRA-LMA", // String      suffix,
						0, // int         num_slices, // (0 - all)
						null); // int []      wh);
			}
			
			// re-load , should create quadCLTs[ref_index].dsi
			double [][] dls = {
					combo_dsn_final[COMBO_DSN_INDX_DISP],
					combo_dsn_final[COMBO_DSN_INDX_LMA],
					combo_dsn_final[COMBO_DSN_INDX_STRENGTH]
			};
			double [][] ds = conditionInitialDS(
					true, // boolean        use_conf,       // use configuration parameters, false - use following  
					clt_parameters,      // CLTParameters  clt_parameters,
					dls,                 // double [][]    dls
					quadCLTs[ref_index], // QuadCLT        scene,
					debugLevel);         // int debug_level)
			double [] disparity = ds[0];
			double [] strength = ds[1];		
			double [][]szxy = getSceneSZXY(
					quadCLTs[ref_index],    // QuadCLT scene,
					range_disparity_offset, // double    disparity_offset,
					range_min_strength,     // double    min_strength,
					range_max,              // double    max_range,
					disparity,              // double [] disparity,
					strength);              // double [] strength)
			String [] szxy_titles = {"strength", "Z(m)", "X(m)", "Y(m)"};
			TileProcessor tp = quadCLTs[ref_index].getTileProcessor();
			int tilesX =         tp.getTilesX();
	        int tilesY =         tp.getTilesY();
	        ImagePlus impSZXY = (new ShowDoubleFloatArrays()).makeArrays(
        			szxy,
        			tilesX,
        			tilesY,
        			quadCLTs[ref_index].getImageName()+"_SZXY",
        			szxy_titles);
	        if (export_dsi_image) {
	        	quadCLTs[ref_index].saveImagePlusInModelDirectory(
	        			null, // "GPU-SHIFTED-BACKGROUND", // String      suffix,
	        			impSZXY); // ImagePlus   imp)
	        }
	        if (show_dsi_image) {
	        	impSZXY.show();

	        }
		}
		if (export_ml_files) { 
			if (combo_dsn_final == null) { // always re-read?
				combo_dsn_final =quadCLTs[ref_index].readDoubleArrayFromModelDirectory( // always re-read?
						"-INTER-INTRA-LMA", // String      suffix,
						0, // int         num_slices, // (0 - all)
						null); // int []      wh);
			}
			double [][] combo_dsn_final_filtered = 
					conditionComboDsnFinal(
							true,                // boolean        use_conf,       // use configuration parameters, false - use following  
							clt_parameters,      // CLTParameters  clt_parameters,
							combo_dsn_final,     // double [][]    combo_dsn_final, // dls,
							quadCLTs[ref_index], // QuadCLT        scene,
							debugLevel); // int            debugLevel);// > 0
			intersceneMlExport(
					clt_parameters,      // CLTParameters        clt_parameters,
					ers_reference,       // ErsCorrection        ers_reference,
					quadCLTs,            // QuadCLT []           scenes,
					colorProcParameters, // ColorProcParameters  colorProcParameters,
					combo_dsn_final_filtered, // double [][]          combo_dsn_final,
					debugLevel); // int                  debug_level
					
			
		}
		if (videos != null) {
			videos[0] = video_list.toArray(new String[0]);
		}
		if (stereo_widths != null) {
			stereo_widths[0] = new int [stereo_widths_list.size()];
			for (int i = 0; i < stereo_widths[0].length; i++) {
				stereo_widths[0][i] = stereo_widths_list.get(i);
			}
		}
		if (start_ref_pointers != null) {
			start_ref_pointers[0] = earliest_scene;
		}
		System.out.println("buildSeries(): DONE"); //
		return quadCLTs[ref_index].getX3dTopDirectory();
    }
    
    public void testERS(
    		CLTParameters clt_parameters,
    		int           indx0, // reference scene in a pair
    		int           indx1, // other scene in a pair
    		QuadCLT []    quadCLTs,
    		int           debugLevel) {
    	// First create a pair of images, similar to renderSceneSequence()
		boolean show_color =  clt_parameters.imp.show_mapped_color;
		boolean show_mono =   clt_parameters.imp.show_mapped_mono;
	    boolean use_combo_dsi =        clt_parameters.imp.use_combo_dsi;
	    boolean use_lma_dsi =          clt_parameters.imp.use_lma_dsi;
		
    	int ref_index = quadCLTs.length-1;
    	indx0 = ref_index; // disregard initial setting, set to reference
		int tilesX =  quadCLTs[ref_index].getTileProcessor().getTilesX();
        int tilesY =  quadCLTs[ref_index].getTileProcessor().getTilesY();
        double [] disparity_raw = new double [tilesX * tilesY];
        Arrays.fill(disparity_raw,clt_parameters.disparity);
        double [][] combo_dsn_final = quadCLTs[ref_index].readDoubleArrayFromModelDirectory(
        			"-INTER-INTRA-LMA", // String      suffix,
        			0, // int         num_slices, // (0 - all)
        			null); // int []      wh);
        double [][] dls = {
        		combo_dsn_final[COMBO_DSN_INDX_DISP],
        		combo_dsn_final[COMBO_DSN_INDX_LMA],
        		combo_dsn_final[COMBO_DSN_INDX_STRENGTH]
        };
        double [][] ds = conditionInitialDS(
        		true,                // boolean        use_conf,       // use configuration parameters, false - use following  
        		clt_parameters,      // CLTParameters  clt_parameters,
        		dls,                 // double [][]    dls
        		quadCLTs[ref_index], // QuadCLT        scene,
        		debugLevel);			
        double [] disparity_fg = ds[0]; // combo_dsn_final[COMBO_DSN_INDX_DISP_FG];
        double [] interscene_ref_disparity = null; // keep null to use old single-scene disparity for interscene matching
        if (use_combo_dsi) {
        	interscene_ref_disparity = ds[0].clone(); // use_lma_dsi ?
        	if (use_lma_dsi) {
        		for (int i = 0; i < interscene_ref_disparity.length; i++) {
        			if (Double.isNaN(dls[1][i])) {
        				interscene_ref_disparity[i] = Double.NaN;
        			}
        		}
        	}
        }
        
        int [] other_ref = {indx1, indx0};
		ErsCorrection ers_reference = quadCLTs[ref_index].getErsCorrection();
        ImageStack stack_scenes_color = null;
        ImageStack stack_scenes_mono = null;
        int sensor_num = 0;
        int sensor_mask = 1 << sensor_num;
        
        String suffix_color = quadCLTs[indx0].getImageName()+"-"+quadCLTs[indx1].getImageName()+"-ERS_TEST-COLOR";
        String suffix_mono = quadCLTs[indx0].getImageName()+"-"+quadCLTs[indx1].getImageName()+"-ERS_TEST-MONO";
        double [][] dxyzatr_dt = new double[quadCLTs.length][];
        for (int iscene = 0; iscene <2; iscene++) {
        	int nscene = other_ref[iscene];
			String ts = quadCLTs[nscene].getImageName();
        	int nscene0 = nscene - ((nscene >0)? 1:0);
        	int nscene1 = nscene + ((nscene < ref_index)? 1:0);
        	double dt = quadCLTs[nscene1].getTimeStamp() - quadCLTs[nscene0].getTimeStamp();
        	String ts0 = quadCLTs[nscene0].getImageName();
        	String ts1 = quadCLTs[nscene1].getImageName();
    		double [] scene_xyz0 = ers_reference.getSceneXYZ(ts0);
    		double [] scene_atr0 = ers_reference.getSceneATR(ts0);
    		double [] scene_xyz1 = (nscene1== ref_index)? ZERO3:ers_reference.getSceneXYZ(ts1);
    		double [] scene_atr1 = (nscene1== ref_index)? ZERO3:ers_reference.getSceneATR(ts1);
    		dxyzatr_dt[nscene] = new double[6];
    		for (int i = 0; i < 3; i++) {
    			dxyzatr_dt[nscene][i] = (scene_xyz1[i]-scene_xyz0[i])/dt;
    			dxyzatr_dt[nscene][i + 3] = (scene_atr1[i]-scene_atr0[i])/dt;
    		}
			double []   scene_xyz = ZERO3;
			double []   scene_atr = ZERO3;
			if (nscene != ref_index) {
				scene_xyz = ers_reference.getSceneXYZ(ts);
				scene_atr = ers_reference.getSceneATR(ts);
				double []   scene_ers_xyz_dt = ers_reference.getSceneErsXYZ_dt(ts);
				double []   scene_ers_atr_dt = ers_reference.getSceneErsATR_dt(ts);
				quadCLTs[nscene].getErsCorrection().setErsDt(
						scene_ers_xyz_dt, // double []    ers_xyz_dt,
						scene_ers_atr_dt); // double []    ers_atr_dt)(ers_scene_original_xyz_dt);
			}
			
			if (show_color) {
				ImagePlus imp_scene_color = QuadCLT.renderGPUFromDSI(
						sensor_mask,         // final int         sensor_mask,
						false,               // final boolean     merge_channels,
						null,                // final Rectangle   full_woi_in,      // show larger than sensor WOI (or null)
						clt_parameters,      // CLTParameters     clt_parameters,
						disparity_fg,        // double []         disparity_ref,
						scene_xyz,           // final double []   scene_xyz, // camera center in world coordinates
						scene_atr,           // final double []   scene_atr, // camera orientation relative to world frame
						quadCLTs[nscene],    // final QuadCLT     scene,
						quadCLTs[ref_index], // final QuadCLT     ref_scene, // now - may be null - for testing if scene is rotated ref
						true,                // final boolean     toRGB,
						"",                  // String            suffix, no suffix here
						threadsMax,          // int               threadsMax,
						debugLevel);         // int         debugLevel)
				if (stack_scenes_color == null) {
					stack_scenes_color = new ImageStack(imp_scene_color.getWidth(),imp_scene_color.getHeight());
				}
				stack_scenes_color.addSlice(
						nscene+":"+ts,
						imp_scene_color.getStack().getPixels(sensor_num+1));
			}
			if (show_mono) {
				ImagePlus imp_scene_mono = QuadCLT.renderGPUFromDSI(
						sensor_mask,         // final int         sensor_mask,
						false,               // final boolean     merge_channels,
						null,                // final Rectangle   full_woi_in,      // show larger than sensor WOI (or null)
						clt_parameters,      // CLTParameters     clt_parameters,
						disparity_fg,        // double []         disparity_ref,
						scene_xyz,           // final double []   scene_xyz, // camera center in world coordinates
						scene_atr,           // final double []   scene_atr, // camera orientation relative to world frame
						quadCLTs[nscene],    // final QuadCLT     scene,
						quadCLTs[ref_index], // final QuadCLT     ref_scene, // now - may be null - for testing if scene is rotated ref
						false,               // final boolean     toRGB,
						"",                  // String            suffix, no suffix here
						threadsMax,          // int               threadsMax,
						debugLevel);         // int         debugLevel)
				if (stack_scenes_mono == null) {
					stack_scenes_mono = new ImageStack(imp_scene_mono.getWidth(),imp_scene_mono.getHeight());
				}
				stack_scenes_mono.addSlice(
						nscene+":"+ts,
						imp_scene_mono.getStack().getPixels(sensor_num+1));
			}
        }
        if (show_color) {
        	ImagePlus imp_scenes_color = new ImagePlus(suffix_color, stack_scenes_color);
        	imp_scenes_color.getProcessor().resetMinAndMax();
        	imp_scenes_color.show();
        }
        if (show_mono) {
        	ImagePlus imp_scenes_mono = new ImagePlus(suffix_mono, stack_scenes_mono);
        	imp_scenes_mono.getProcessor().resetMinAndMax();
        	imp_scenes_mono.show();
        }
		
        String [] vel_names = {"Vx","Vy","Vz","Vaz","Vtl","Vrl"};
        String [] vel_units = {"m/s","m/s","m/s","rad/s","rad/s","rad/s","rad/s"};
        System.out.println("Pair: ref. scene "+indx0+": ("+quadCLTs[indx0].getImageName()+")"+
        ", other scene "+indx1+": ("+quadCLTs[indx1].getImageName()+")");

		for (int i = 0; i < vel_names.length; i++) {
        	System.out.println(String.format(
        			"%3s: other: %3d (%8.3f%-6s), ref: %3d (%8.3f%-6s)",
        			vel_names[i], other_ref[0], dxyzatr_dt[other_ref[0]][i], vel_units[i],
        			other_ref[1], dxyzatr_dt[other_ref[1]][i], vel_units[i]));
		}
		double []   scene_xyz_pre = ZERO3;
		double []   scene_atr_pre = ZERO3;
		String ts_other = quadCLTs[other_ref[0]].getImageName(); // other scene

		scene_xyz_pre = ers_reference.getSceneXYZ(ts_other);
		scene_atr_pre = ers_reference.getSceneATR(ts_other);
		double []   scene_ers_xyz_dt = ers_reference.getSceneErsXYZ_dt(ts_other);
		double []   scene_ers_atr_dt = ers_reference.getSceneErsATR_dt(ts_other);
		quadCLTs[other_ref[0]].getErsCorrection().setErsDt(
				scene_ers_xyz_dt, // double []    ers_xyz_dt,
				scene_ers_atr_dt); // double []    ers_atr_dt)(ers_scene_original_xyz_dt);
		double []      lma_rms = new double[2];

		double []ref_dt =   dxyzatr_dt[other_ref[1]]; 
		double []scene_dt = dxyzatr_dt[other_ref[0]]; 
		double []ref_xyz_dt = ZERO3; // {ref_dt[0],ref_dt[1],ref_dt[2]}; 
//		double []ref_xyz_dt = {ref_dt[0],ref_dt[1],ref_dt[2]}; 
//		double []ref_xyz_dt = {-ref_dt[0],-ref_dt[1],-ref_dt[2]}; 
		double k = 1.0;
		double []ref_atr_dt = {k*ref_dt[3],k*ref_dt[4],k*ref_dt[5]}; 
		double []scene_xyz_dt = ZERO3; // {scene_dt[0],scene_dt[1],scene_dt[2]}; 
//		double []scene_xyz_dt = {scene_dt[0],scene_dt[1],scene_dt[2]}; 
//		double []scene_xyz_dt = {-scene_dt[0],-scene_dt[1],-scene_dt[2]}; 
		double []scene_atr_dt = {k*scene_dt[3],k*scene_dt[4],k*scene_dt[5]};
		// first - try ATR only
		quadCLTs[other_ref[0]].getErsCorrection().setErsDt( // scene
				scene_xyz_dt, // double []    ers_xyz_dt,
				scene_atr_dt); // double []    ers_atr_dt)(ers_scene_original_xyz_dt);
		quadCLTs[other_ref[1]].getErsCorrection().setErsDt( // reference
				ref_xyz_dt,  // double []    ers_xyz_dt,
				ref_atr_dt); // double []    ers_atr_dt)(ers_scene_original_xyz_dt);
		
		double [][] adjusted_xyzatr = adjustPairsLMAInterscene(
				clt_parameters,                                 // CLTParameters  clt_parameters,
				quadCLTs[other_ref[1]],                         // QuadCLT reference_QuadCLT,
				null, // double []        ref_disparity, // null or alternative reference disparity
				null, // boolean []     reliable_ref, // null or bitmask of reliable reference tiles
				quadCLTs[other_ref[0]],                         // QuadCLT scene_QuadCLT,
				scene_xyz_pre,                                // xyz
				scene_atr_pre,                                // atr
				clt_parameters.ilp.ilma_lma_select,                                  // final boolean[]   param_select,
				clt_parameters.ilp.ilma_regularization_weights,                              //  final double []   param_regweights,
				lma_rms,                                        // double []      rms, // null or double [2]
				clt_parameters.imp.max_rms,                     // double         max_rms,
				clt_parameters.imp.debug_level);                // 1); // -1); // int debug_level);
//		System.out.println("lma_rms={"+lma_rms[0]+","+lma_rms[1]+"}");
		
        ImageStack stack_adjusted_color = null;
        ImageStack stack_adjusted_mono = null;
        String suffix_adjusted_color = quadCLTs[indx0].getImageName()+"-"+quadCLTs[indx1].getImageName()+"-ERS_TEST_ADJUSTED-COLOR";
        String suffix_adjusted_mono = quadCLTs[indx0].getImageName()+"-"+quadCLTs[indx1].getImageName()+"-ERS_TEST_ADJUSTED-MONO";
        for (int iscene = 0; iscene <2; iscene++) {
        	int nscene = other_ref[iscene];
			String ts = quadCLTs[nscene].getImageName();
			double []   scene_xyz = ZERO3;
			double []   scene_atr = ZERO3;
			if (nscene != ref_index) {
				// keep current/last adjusted
				scene_xyz = adjusted_xyzatr[0]; // ers_reference.getSceneXYZ(ts);
				scene_atr = adjusted_xyzatr[1]; // ers_reference.getSceneATR(ts);
				/*
				double []   scene_ers_xyz_dt1 = ers_reference.getSceneErsXYZ_dt(ts);
				double []   scene_ers_atr_dt1 = ers_reference.getSceneErsATR_dt(ts);
				quadCLTs[nscene].getErsCorrection().setErsDt(
						scene_ers_xyz_dt1, // double []    ers_xyz_dt,
						scene_ers_atr_dt1); // double []    ers_atr_dt)(ers_scene_original_xyz_dt);
				*/
			}
			if (show_color) {
				ImagePlus imp_scene_color = QuadCLT.renderGPUFromDSI(
						sensor_mask,         // final int         sensor_mask,
						false,               // final boolean     merge_channels,
						null,                // final Rectangle   full_woi_in,      // show larger than sensor WOI (or null)
						clt_parameters,      // CLTParameters     clt_parameters,
						disparity_fg,        // double []         disparity_ref,
						scene_xyz,           // final double []   scene_xyz, // camera center in world coordinates
						scene_atr,           // final double []   scene_atr, // camera orientation relative to world frame
						quadCLTs[nscene],    // final QuadCLT     scene,
						quadCLTs[ref_index], // final QuadCLT     ref_scene, // now - may be null - for testing if scene is rotated ref
						true,                // final boolean     toRGB,
						"",                  // String            suffix, no suffix here
						threadsMax,          // int               threadsMax,
						debugLevel);         // int         debugLevel)
				if (stack_adjusted_color == null) {
					stack_adjusted_color = new ImageStack(imp_scene_color.getWidth(),imp_scene_color.getHeight());
				}
				stack_adjusted_color.addSlice(
						nscene+":"+ts,
						imp_scene_color.getStack().getPixels(sensor_num+1));
			}
			if (show_mono) {
				ImagePlus imp_scene_mono = QuadCLT.renderGPUFromDSI(
						sensor_mask,         // final int         sensor_mask,
						false,               // final boolean     merge_channels,
						null,                // final Rectangle   full_woi_in,      // show larger than sensor WOI (or null)
						clt_parameters,      // CLTParameters     clt_parameters,
						disparity_fg,        // double []         disparity_ref,
						scene_xyz,           // final double []   scene_xyz, // camera center in world coordinates
						scene_atr,           // final double []   scene_atr, // camera orientation relative to world frame
						quadCLTs[nscene],    // final QuadCLT     scene,
						quadCLTs[ref_index], // final QuadCLT     ref_scene, // now - may be null - for testing if scene is rotated ref
						false,               // final boolean     toRGB,
						"",                  // String            suffix, no suffix here
						threadsMax,          // int               threadsMax,
						debugLevel);         // int         debugLevel)
				if (stack_adjusted_mono == null) {
					stack_adjusted_mono = new ImageStack(imp_scene_mono.getWidth(),imp_scene_mono.getHeight());
				}
				stack_adjusted_mono.addSlice(
						nscene+":"+ts,
						imp_scene_mono.getStack().getPixels(sensor_num+1));
			}
        }
        if (show_color) {
        	ImagePlus imp_adjusted_color = new ImagePlus(suffix_adjusted_color, stack_adjusted_color);
        	imp_adjusted_color.getProcessor().resetMinAndMax();
        	imp_adjusted_color.show();
        }
        if (show_mono) {
        	ImagePlus imp_adjusted_mono = new ImagePlus(suffix_adjusted_mono, stack_adjusted_mono);
        	imp_adjusted_mono.getProcessor().resetMinAndMax();
        	imp_adjusted_mono.show();
        }
		return;
    }
    
  
    public static int [][] getERSStats(
    		CLTParameters clt_parameters,
    		QuadCLT []    quadCLTs,
    		double [][] dxyzatr_dt) {
        int ref_index = quadCLTs.length -1;
		ErsCorrection ers_reference = quadCLTs[ref_index].getErsCorrection();
        int dbg_scene = -95;
        int [][] min_max_xyzatr = new int [6][2];
        for (int i = 0; i < min_max_xyzatr.length; i++) {
        	Arrays.fill(min_max_xyzatr[i],-1);
        }
		if (dxyzatr_dt == null) {
			dxyzatr_dt = new double[quadCLTs.length][];
		} else {
			Arrays.fill(dxyzatr_dt, null);
		}
		for (int nscene =  0; nscene < quadCLTs.length ; nscene++) if (quadCLTs[nscene] != null){
			if (nscene== dbg_scene) {
				System.out.println("renderSceneSequence(): nscene = "+nscene);
			}
			String ts = quadCLTs[nscene].getImageName();
			if ((ers_reference.getSceneXYZ(ts) != null) && (ers_reference.getSceneATR(ts) != null)) {
				int nscene0 = nscene - ((nscene >0)? 1:0);
				int nscene1 = nscene + ((nscene < ref_index)? 1:0);
				String ts0 = quadCLTs[nscene0].getImageName();
				if ((ers_reference.getSceneXYZ(ts0) == null) || (ers_reference.getSceneATR(ts0) == null)) {
					nscene0 = nscene;
					ts0 = quadCLTs[nscene0].getImageName();
				}
				double dt = quadCLTs[nscene1].getTimeStamp() - quadCLTs[nscene0].getTimeStamp();
//				String ts0 = quadCLTs[nscene0].getImageName();
				String ts1 = quadCLTs[nscene1].getImageName();
				double [] scene_xyz0 = ers_reference.getSceneXYZ(ts0);
				double [] scene_atr0 = ers_reference.getSceneATR(ts0);
				double [] scene_xyz1 = (nscene1== ref_index)? ZERO3:ers_reference.getSceneXYZ(ts1);
				double [] scene_atr1 = (nscene1== ref_index)? ZERO3:ers_reference.getSceneATR(ts1);
				dxyzatr_dt[nscene] = new double[6];
				for (int i = 0; i < 3; i++) {
					dxyzatr_dt[nscene][i] = (scene_xyz1[i]-scene_xyz0[i])/dt;
					dxyzatr_dt[nscene][i + 3] = (scene_atr1[i]-scene_atr0[i])/dt;
				}
				for (int i = 0; i < dxyzatr_dt[nscene].length; i++) {
					if ((min_max_xyzatr[i][0] < 0) || // update index of minimum
							(dxyzatr_dt[nscene][i] < dxyzatr_dt[min_max_xyzatr[i][0]][i]) ) {
						min_max_xyzatr[i][0] = nscene;
					}
					if ((min_max_xyzatr[i][1] < 0) || // update index of maximum
							(dxyzatr_dt[nscene][i] > dxyzatr_dt[min_max_xyzatr[i][1]][i]) ) {
						min_max_xyzatr[i][1] = nscene;
					}
				}
			}   
		}
		return min_max_xyzatr;
    }
    
    
    public static ImagePlus renderSceneSequence(
    		CLTParameters clt_parameters,
    		Rectangle     fov_tiles,
    		int           mode3d,
    		boolean       toRGB,
    		double []     stereo_xyz, // offset reference camera {x,y,z}
    		int           sensor_mask,
    		String        suffix_in,
    		double []     ref_disparity,			
    		QuadCLT []    quadCLTs,
    		int           debugLevel) {
    	double [] stereo_atr = ZERO3; // maybe later play with rotated camera
		boolean um_mono =            clt_parameters.imp.um_mono;
		double  um_sigma =           clt_parameters.imp.um_sigma;
		double  um_weight =          clt_parameters.imp.um_weight;
    	final float fum_weight = (float)  um_weight; 
    	
    	boolean merge_all = clt_parameters.imp.merge_all;
    	if (mode3d < 1) {
    		merge_all = false;
    	}
    	if (merge_all) {
    		sensor_mask = 1;
    	}
    	String        suffix = suffix_in+((mode3d > 0)?(merge_all?"-MERGED":"-SINGLE"):"");
        int ref_index = quadCLTs.length -1;
    	int num_sens = quadCLTs[ref_index].getNumSensors();
		ErsCorrection ers_reference = quadCLTs[ref_index].getErsCorrection();
        int num_used_sens = 0;
        for (int i = 0; i < num_sens; i++) if (((sensor_mask >> i) & 1) != 0) num_used_sens++;
        int [] channels = new int [num_used_sens];
        int nch = 0;
        for (int i = 0; i < num_sens; i++) if (((sensor_mask >> i) & 1) != 0) channels[nch++] = i;
        ImageStack stack_scenes = null;
        int dbg_scene = -95;
		for (int nscene =  0; nscene < quadCLTs.length ; nscene++) if (quadCLTs[nscene] != null){
			if (nscene== dbg_scene) {
				System.out.println("renderSceneSequence(): nscene = "+nscene);
			}
			String ts = quadCLTs[nscene].getImageName();
			double []   scene_xyz = ZERO3;
			double []   scene_atr = ZERO3;
//			if ((nscene != ref_index) && (mode3d >= 0)) {
			if (nscene != ref_index) { // Check even for raw, so video frames will match in all modes 
				scene_xyz = ers_reference.getSceneXYZ(ts);
				scene_atr = ers_reference.getSceneATR(ts);
				if ((scene_atr==null) || (scene_xyz == null)) {
					continue;
				}
				if (mode3d >= 0) {
					double []   scene_ers_xyz_dt = ers_reference.getSceneErsXYZ_dt(ts);
					double []   scene_ers_atr_dt = ers_reference.getSceneErsATR_dt(ts);
					quadCLTs[nscene].getErsCorrection().setErsDt(
							scene_ers_xyz_dt, // double []    ers_xyz_dt,
							scene_ers_atr_dt); // double []    ers_atr_dt)(ers_scene_original_xyz_dt);
				} else { // ugly, restore for raw mode that should not be rotated/shifted
					scene_xyz = ZERO3;
					scene_atr = ZERO3;
				}
			}
			if (stereo_xyz != null) { // offset all, including reference scene
				double [][] combo_xyzatr = ErsCorrection.combineXYZATR(
						stereo_xyz,  // double [] reference_xyz,
						stereo_atr,  // double [] reference_atr, 
						scene_xyz,   // double [] scene_xyz,
						scene_atr);  // double [] scene_atr) 
				scene_xyz = combo_xyzatr[0];
				scene_atr = combo_xyzatr[1];
			}
			int sm = merge_all? -1: sensor_mask;
			ImagePlus imp_scene = QuadCLT.renderGPUFromDSI(
					sm,                  // final int         sensor_mask,
					merge_all,           // final boolean     merge_channels,
					fov_tiles,           // testr, // null,                // final Rectangle   full_woi_in,      // show larger than sensor WOI (or null)
					clt_parameters,      // CLTParameters     clt_parameters,
					ref_disparity,       // double []         disparity_ref,
					// not used, just as null/not null now
					scene_xyz,           // final double []   scene_xyz, // camera center in world coordinates
					scene_atr,           // final double []   scene_atr, // camera orientation relative to world frame
					quadCLTs[nscene],    // final QuadCLT     scene,
					quadCLTs[ref_index], // final QuadCLT     ref_scene, // now - may be null - for testing if scene is rotated ref
					toRGB,               // final boolean     toRGB,
					"", // String            suffix, no suffix here
					QuadCLT.THREADS_MAX,          // int               threadsMax,
					debugLevel);         // int         debugLevel)
			if (stack_scenes == null) {
				stack_scenes = new ImageStack(imp_scene.getWidth(),imp_scene.getHeight());
			}
			for (int i = 0; i < channels.length; i++) {
				stack_scenes.addSlice(
						ts+"-"+channels[i],
						imp_scene.getStack().getPixels(i+1));
			}
		}
		// Apply unsharp mask here, in parallel
		if (um_mono && !toRGB) {
			final ImageStack fstack_scenes = stack_scenes;
			final int nSlices = fstack_scenes.getSize();
			final Thread[] threads = ImageDtt.newThreadArray(QuadCLT.THREADS_MAX);
			final AtomicInteger ai = new AtomicInteger(0);
			for (int ithread = 0; ithread < threads.length; ithread++) {
				threads[ithread] = new Thread() {
					public void run() {
						for (int nSlice = ai.getAndIncrement(); nSlice < nSlices; nSlice = ai.getAndIncrement()) {
							FloatProcessor fp = (FloatProcessor) fstack_scenes.getProcessor(nSlice+1);
							float [] fpixels = (float[]) fstack_scenes.getPixels(nSlice+1);
							float [] fpixels_orig = fpixels.clone();
							(new GaussianBlur()).blurFloat(
									fp,       // FloatProcessor ip,
									um_sigma, // double sigmaX,
									um_sigma, // double sigmaY,
									0.01);    // double accuracy)
							for (int i = 0; i < fpixels.length; i++) {
								fpixels[i] = fpixels_orig[i] - fum_weight * fpixels[i];
							}
						}
					}
				};
			}		      
			ImageDtt.startAndJoin(threads);
		}
		
		ImagePlus imp_scenes = new ImagePlus(suffix, stack_scenes);
		imp_scenes.getProcessor().resetMinAndMax();
    	return imp_scenes;
    }
    
    
    
    public double [][] getSceneSZXY(
    	QuadCLT   scene,
    	double    disparity_offset,
    	double    min_strength,
    	double    max_range,
    	double [] disparity0,
    	double [] strength){ // may be null
    	double [] disparity = disparity0.clone();
		for (int i = 0; i < disparity.length; i++) {
			disparity[i] -= disparity_offset;
		}
    	double [][] xyz = transformToWorldXYZ(
    			disparity, // [final double []   disparity_ref, // invalid tiles - NaN in disparity
    			scene, // final QuadCLT     quadClt, // now - may be null - for testing if scene is rotated ref
    			threadsMax); // int               threadsMax);
    	
    	double [][] szxy = new double [4][disparity.length];
    	szxy[0] = strength;
    	for (int i = 1; i < szxy.length; i++) {
    		Arrays.fill(szxy[i], Double.NaN);
    	}
    	for (int nTile = 0; nTile < xyz.length; nTile++) if (xyz[nTile] != null) {
    		if ((xyz[nTile][2] < 0) && (xyz[nTile][2] > -max_range) && ((strength == null) ||(strength[nTile] >= min_strength))) {
    			szxy[1][nTile] = -xyz[nTile][2];
    			szxy[2][nTile] = xyz[nTile][0];
    			szxy[3][nTile] = xyz[nTile][1];
    		}
    	}
    	return szxy;
    }
    

    public double [] spiralSearchATR(
    		CLTParameters             clt_parameters,
    		QuadCLT                   reference_QuadClt,
    		QuadCLT                   scene_QuadClt,
    		boolean []                reliable_ref,
    		int                       debugLevel
    		) {
    	int [][] offset_start_corner = {{-1,-1},{1,-1},{1,1},{-1,1}}; //{x,y}
    	int [][] offset_move = {{1,0},{0,1},{-1,0},{0,-1}};
    	int    margin =            clt_parameters.imp.margin;
    	int    sensor_mask_inter = clt_parameters.imp.sensor_mask_inter ; //-1;
    	float [][][] facc_2d_img = new float [1][][];
    	int    pix_step =          clt_parameters.imp.pix_step;      // 4;
    	int    search_rad =        clt_parameters.imp.search_rad;    // 10;
    	double maybe_sum =         clt_parameters.imp.maybe_sum;     // 8.0;
    	double sure_sum =          clt_parameters.imp.sure_sum;     // 30.0;
    	double maybe_avg =         clt_parameters.imp.maybe_avg;     // 0.005;
    	double sure_avg =          clt_parameters.imp.sure_avg;     // 0.015;
    	
    	double  max_search_rms =   clt_parameters.imp.max_search_rms; //             1.2;   // good - 0.34, so-so - 0.999
    	double  maybe_fom =        clt_parameters.imp.maybe_fom;      //             3.0;   // good - 30, second good - 5
    	double  sure_fom =         clt_parameters.imp.sure_fom;       //            10.0;   // good - 30, second good - 10
    	
		double [][] pose = new double [2][3];
    	double angle_per_step = reference_QuadClt.getGeometryCorrection().getCorrVector().getTiltAzPerPixel() * pix_step;
    	double [][] atrs = new double [(2*search_rad+1)*(2*search_rad+1)][];
    	boolean [] good_offset =    new boolean [atrs.length];
    	double [] confidences_sum = new double [atrs.length];
    	double [] confidences_avg = new double [atrs.length];
    	double [] confidences_fom = new double [atrs.length];
    	int min_defined =           75;
    	int rad = 0, dir=0, n=0;
    	int ntry = 0;
    	int num_good=0;
    	double [] use_atr = null;
    	try_around:
    		for (rad = 0; rad <= search_rad; rad++) {
    			for (dir = 0; dir < ((rad==0)?1:4); dir++) {
    				int n_range = (rad > 0) ? (2* rad) : 1;
    				for (n = 0; n < n_range; n++) {
    					int ix = rad*offset_start_corner[dir][0] + n * offset_move[dir][0];
    					int iy = rad*offset_start_corner[dir][1] + n * offset_move[dir][1];;
    					double [] atr = {
    							pose[1][0]+ix * angle_per_step,
    							pose[1][1]+iy * angle_per_step,
    							pose[1][2]};
    					if (debugLevel > -1) {
    						System.out.println("buildSeries(): trying adjustPairsLMA() with initial offset azimuth: "+
    								atr[0]+", tilt ="+atr[1]);
    					}
    					double [][][] coord_motion = interCorrPair( // new double [tilesY][tilesX][][];
    							clt_parameters,      // CLTParameters  clt_parameters,
    							reference_QuadClt,   // QuadCLT reference_QuadCLT,
    							null, // double []        ref_disparity, // null or alternative reference disparity  
    							scene_QuadClt,       // QuadCLT scene_QuadCLT,
    							pose[0],             // camera_xyz0,         // xyz
    							atr,                 // camera_atr0,         // pose[1], // atr
    							reliable_ref, // null,                // final boolean [] selection, // may be null, if not null do not  process unselected tiles
    							margin,              // final int        margin,
    							sensor_mask_inter,   // final int        sensor_mask_inter, // The bitmask - which sensors to correlate, -1 - all.
    							facc_2d_img,         // final float [][][]   accum_2d_corr, // if [1][][] - return accumulated 2d correlations (all pairs)final float [][][]   accum_2d_corr, // if [1][][] - return accumulated 2d correlations (all pairs)
    							null,                //	final float [][][] dbg_corr_fpn,
    							false,               // boolean            near_important, // do not reduce weight of the near tiles
    							true,                // boolean            all_fpn,        // do not lower thresholds for non-fpn (used during search)
    							clt_parameters.imp.debug_level, // int                imp_debug_level,
    							clt_parameters.imp.debug_level); // 1); // -1); // int debug_level);
    					int num_defined = 0;
    					double sum_strength = 0.0;
    					// calculate average vector and standard deviation for its components
    					double swx=0.0,swy = 0.0, swx2 = 0.0, swy2 = 0.0, sw2=0.0 ; 
    					
    					for (int i = 0; i < coord_motion[1].length; i++) if (coord_motion[1][i] != null){
    						swx +=  coord_motion[1][i][2] * coord_motion[1][i][0];
    						swx2 += coord_motion[1][i][2] * coord_motion[1][i][0] * coord_motion[1][i][0];
    						swy +=  coord_motion[1][i][2] * coord_motion[1][i][1];
    						swy2 += coord_motion[1][i][2] * coord_motion[1][i][1] * coord_motion[1][i][1];
    						sum_strength += coord_motion[1][i][2];
    						sw2+=coord_motion[1][i][2]*coord_motion[1][i][2];
    						num_defined++;
    					}
    					double vx_avg = swx/sum_strength;
    					double vy_avg = swy/sum_strength;
    					double l2x = swx2/sum_strength - vx_avg * vx_avg;
    					double l2y = swy2/sum_strength - vy_avg * vy_avg;
    					double rms_x = Math.sqrt(l2x);
    					double rms_y= Math.sqrt(l2y);
    					double rms = Math.sqrt(l2x+l2y);
    					double fom = sum_strength/rms/Math.sqrt(0.01 * num_defined); // reduce influence of number defined
    					double avg_strength = (num_defined > 0)?(sum_strength/num_defined) : 00;
    					atrs[ntry] = atr;
    					confidences_sum[ntry] = sum_strength;
    					confidences_avg[ntry] = avg_strength;
    					confidences_fom[ntry] = fom;
    					boolean good = (num_defined > min_defined) && (rms < max_search_rms) && (fom > maybe_fom) && (sum_strength > maybe_sum) && (avg_strength > maybe_avg);
    					good_offset[ntry] = good;
    					if (debugLevel > -2) {
    						System.out.println ((good?"* ":"  ")+"ntry = "+ntry+", num_defined = "+num_defined+" ("+min_defined+")"+
    								", FOM="+fom+
    								", vx_avg = "+vx_avg+", vy_avg="+vy_avg+
    								", sum_strength = "+sum_strength+", rms="+rms+
    								", str2="+Math.sqrt(sw2)+ ", str2_avg="+Math.sqrt(sw2/num_defined)+
    								", avg_strength = "+avg_strength+
    								", rms_x="+rms_x+", rms_y="+rms_y);
    					}
    					if (good) {
    						if ((fom > sure_fom) || (sum_strength > sure_sum) || (avg_strength > sure_avg)) {
    							use_atr = atr;
    							break try_around;
    						}
    						num_good++;
    					}
//    					if (ntry == 105) {
//    						System.out.println("spiralSearchATR(): ntry="+ntry);
//    					}
    					ntry++;
    				}
    			}
    		}
    	if ((use_atr == null) && (num_good == 0)) {
    		System.out.println("Failed to find a match between the reference scene ("+reference_QuadClt.getImageName() +
    				") and a previous one ("+scene_QuadClt.getImageName()+")");
    		return null;
    	}
    	if (use_atr == null) {
    		int best_ntry = -1;
    		for (int i = 0; i <  atrs.length; i++) if (good_offset[i]){
       			if ((best_ntry < 0) || (confidences_fom[i] > confidences_fom[best_ntry])) {
    				best_ntry = i;
    			}
    		}
    		use_atr = atrs[best_ntry];
    	}
    	return use_atr;
    }
	
	
	public void adjustSeries(
			CLTParameters  clt_parameters,			
			double         k_prev, 
			QuadCLT []     scenes, // ordered by increasing timestamps
			int            ref_index,
			int            debug_level
			)
	{
		boolean show_results = clt_parameters.ilp.ilma_debug_adjust_series;
		boolean pattern_mode = clt_parameters.ofp.pattern_mode;
		boolean high_res_motion_vectors = true; // use new 05/20/2022 mode
		if (pattern_mode) {
			if (clt_parameters.ofp.center_index < 0) {
				double [][][] atrs = new double [scenes.length][][];
				atrs[scenes.length - 1] = new double[2][3]; 
				for (int i =  scenes.length - 2; i >= 0 ; i--) {
					String scene_ts =            scenes[i].getImageName(); // it should be present in the scenes[i+1] scenes
					ErsCorrection ers_scene_last_known = scenes[i+1].getErsCorrection();
					double [][] last_known_atr = atrs[i+1];
					double [] new_from_last_xyz = ers_scene_last_known.getSceneXYZ(scene_ts);
					double [] new_from_last_atr = ers_scene_last_known.getSceneATR(scene_ts);
					// combine two rotations and two translations (translations will be zero for pattern_mode) 
					atrs[i]=ErsCorrection.combineXYZATR(
							last_known_atr[0],  // double [] reference_xyz,
							last_known_atr[1],  // double [] reference_atr,
							new_from_last_xyz,  // double [] scene_xyz,
							new_from_last_atr); // double [] scene_atr)
				}
				double [][] xyzatr_avg = new double[2][3];
				for (int i = 0; i < atrs.length; i++) {
					for (int k = 0; k < xyzatr_avg.length; k++) {
						for (int j = 0; j < 3; j++) {
							xyzatr_avg[k][j] += atrs[i][k][j] / atrs.length;
						}
					}
				}
				double [] wxyzatr = {0.0, 1.0}; // weight of xyz offset - not used now 
				int nearest = 0;
				for (int i = 1; i < atrs.length; i++) {
					double d2best = 0, d2this=0;
					for (int k = 0; k < xyzatr_avg.length; k++) {
						for (int j = 0; j < 3; j++) {
							double d= atrs[i][k][j] - xyzatr_avg[k][j];
							d2this += wxyzatr[k] * d * d;
							double d1= atrs[nearest][k][j] - xyzatr_avg[k][j];
							d2best += wxyzatr[k] * d1 * d1;
							if (d2this < d2best) {
								nearest = i;
							}
						}
					}
				}
				clt_parameters.ofp.center_index = nearest;
			}
			ref_index = clt_parameters.ofp.center_index;
		}
		
		if (ref_index < 0) {
			ref_index += scenes.length;
		}
		double [][][] scenes_xyzatr = new double [scenes.length][][]; // previous scene relative to the next one
		QuadCLT reference_QuadClt = scenes[ref_index]; // scenes.length-1]; // last acquired
		ErsCorrection ers_reference = reference_QuadClt.getErsCorrection();
		// modify LMA parameters to freeze reference ERS, remove pull on scene ERS
		boolean[]   param_select2 =     clt_parameters.ilp.ilma_lma_select.clone();             // final boolean[]   param_select,
		boolean[]   param_select3 =     clt_parameters.ilp.ilma_lma_select.clone();             // final boolean[]   param_select,
		double []   param_regweights2 = clt_parameters.ilp.ilma_regularization_weights; //  final double []   param_regweights,
		double []   param_regweights3 = clt_parameters.ilp.ilma_regularization_weights; //  final double []   param_regweights,
		boolean delete_scene_asap = false; //  (debug_level < 10); // to save memory
		// freeze reference ERS, free scene ERS
		for (int j = 0; j <3; j++) {
			param_select2[ErsCorrection.DP_DVX  + j] = false;
			param_select2[ErsCorrection.DP_DVAZ + j] = false;
			param_regweights2[ErsCorrection.DP_DSVX +  j] = 0.0;
			param_regweights2[ErsCorrection.DP_DSVAZ + j] = 0.0;
		}
		for (int j = 0; j <3; j++) {
			param_select3[ErsCorrection.DP_DVX  + j] = false;
			param_select3[ErsCorrection.DP_DVAZ + j] = false;
			param_select3[ErsCorrection.DP_DSVX  + j] = clt_parameters.ilp.ilma_ers_adj_lin; // disabling, may check with high rot speed
			param_select3[ErsCorrection.DP_DSVAZ + j] = clt_parameters.ilp.ilma_ers_adj_ang; // so far ers correction noise is too high to compare
			param_regweights3[ErsCorrection.DP_DSVX +  j] = 0.0;
			param_regweights3[ErsCorrection.DP_DSVAZ + j] = 0.0;
		}
		
		if (show_results) { // useless before references to the ref_scene
			double [][][] ers_current = new double [scenes.length][][];
			double [] scene_xyz_dt;
			double [] scene_atr_dt;
			for (int nscene = 0; nscene < scenes.length; nscene++) {
				if (nscene == ref_index) {
					scene_xyz_dt = ers_reference.getErsXYZ_dt();
					scene_atr_dt = ers_reference.getErsATR_dt();
				} else {
					String sts = scenes[nscene].getImageName();
					scene_xyz_dt = ers_reference.getSceneErsXYZ_dt(sts); 
					scene_atr_dt = ers_reference.getSceneErsATR_dt(sts); 
				}
				if ((scene_xyz_dt != null) &&(scene_atr_dt != null)) {
					ers_current[nscene] = new double[][] {scene_xyz_dt, scene_atr_dt};
				} else {
					System.out.println("adjustSeries(): null for nscene="+nscene);
				}
			}
			
			int dbg_w = ers_current.length;
			int dbg_h = 6;
			double [] dbg_img = new double [dbg_w * dbg_h];
			Arrays.fill(dbg_img, Double.NaN);
			for (int dh = 0; dh  < dbg_h; dh++) {
				for (int dw = 0; dw  < dbg_w; dw++) {
					if (ers_current[dw] != null) {
						dbg_img[dh*dbg_w + dw] = ers_current[dw][dh / 3][dh % 3];
					}
				}
			}
			(new ShowDoubleFloatArrays()).showArrays(
					dbg_img,
					dbg_w,
					dbg_h,
					"ers_pre_adjust"); //	dsrbg_titles);
			System.out.println("adjustSeries(): dbg");
		}
		
		// used only for debug
		TileProcessor tp = reference_QuadClt.getTileProcessor();
		int tilesX =         tp.getTilesX();
        int tilesY =         tp.getTilesY();
		int transform_size = tp.getTileSize();
		int macroTilesX =          tilesX/transform_size;
		int macroTilesY =          tilesY/transform_size;
		double [][] dbg_iterdata = show_results ? (new double [3][]): null;
		if (dbg_iterdata != null) {
			for (int ii = 0; ii < dbg_iterdata.length; ii++) {
				dbg_iterdata[ii] = new double [macroTilesY * scenes.length * macroTilesX * clt_parameters.ilp.ilma_num_corr];
			}
		}
		if (debug_level > -1000) {
			dbg_iterdata = null; // does not work now  (at 3964)
		}
		double []     dbg_rms_pre = new double [scenes.length];
		Arrays.fill(dbg_rms_pre, Double.NaN);
		// process scenes before reference
		int dbg_ref_index = 10; // wait for ref_index <= dbg_ref_index and print (manually update dbg_ref_index
		double maximal_series_rms = 0.0;
		double maximal_series_rms1 = 0.0;
		if (ref_index >  1) {
			for (int i =  ref_index - 2; i >= 0 ; i--) {
				QuadCLT scene_QuadClt =      scenes[i];
				String last_known_ts =           scenes[i+1].getImageName(); // it should be present in the reference scene scenes
				String scene_ts =                scenes[i].getImageName(); // it should be present in the scenes[i+1] scenes
				ErsCorrection ers_scene_last_known = scenes[i+1].getErsCorrection();
				ErsCorrection ers_scene =            scene_QuadClt.getErsCorrection();

				double [] last_known_xyz = ers_reference.getSceneXYZ(last_known_ts);
				double [] last_known_atr = ers_reference.getSceneATR(last_known_ts);

				double [] new_from_last_xyz = ers_scene_last_known.getSceneXYZ(scene_ts);
				double [] new_from_last_atr = ers_scene_last_known.getSceneATR(scene_ts);

				// combine two rotations and two translations 
				System.out.println("***** Processing scene "+i+": "+scene_QuadClt.getImageName()+" *****");
				double [][] combo_XYZATR = ErsCorrection.combineXYZATR(
						last_known_xyz,     // double [] reference_xyz,
						last_known_atr,     // double [] reference_atr, // null?
						new_from_last_xyz,  // double [] scene_xyz,
						new_from_last_atr); // double [] scene_atr)

				// before adjusting - save original ERS, restart afterwards
				double [] ers_scene_original_xyz_dt = ers_scene.getErsXYZ_dt();
				double [] ers_scene_original_atr_dt = ers_scene.getErsATR_dt();

				// ers should be correct for both
				double [] lma_rms = new double[2];
				double [][] dbg_img = (dbg_iterdata != null) ? (new double[3][]) : null;
				if(i <= dbg_ref_index) {
					System.out.println ("i = "+i+" <= dbg_ref_index");
				}
				if (high_res_motion_vectors) {
					scenes_xyzatr[i] = adjustPairsLMAInterscene(
							clt_parameters,                                 // CLTParameters  clt_parameters,
							reference_QuadClt,                              // QuadCLT reference_QuadCLT,
							null, // double []        ref_disparity, // null or alternative reference disparity
							null, // boolean []     reliable_ref, // null or bitmask of reliable reference tiles
							scene_QuadClt,                                  // QuadCLT scene_QuadCLT,
							combo_XYZATR[0],                                // xyz
							combo_XYZATR[1],                                // atr
							param_select2,                                  // final boolean[]   param_select,
							param_regweights2,                              //  final double []   param_regweights,
							lma_rms,                                        // double []      rms, // null or double [2]
							clt_parameters.imp.max_rms,                     // double         max_rms,
							clt_parameters.imp.debug_level);                // 1); // -1); // int debug_level);
				} else {
					scenes_xyzatr[i] = adjustPairsLMA(
							clt_parameters,     // CLTParameters  clt_parameters,			
							reference_QuadClt, // QuadCLT reference_QuadCLT,
							scene_QuadClt, // QuadCLT scene_QuadCLT,
							combo_XYZATR[0], // xyz
							combo_XYZATR[1], // atr
							param_select2,             // final boolean[]   param_select,
							param_regweights2, //  final double []   param_regweights,
							lma_rms, // 			double []      rms, // null or double [2]
							dbg_img, // double [][]    dbg_img,
							0.0, // double         max_rms,
							debug_level); // int debug_level)
				}
				if (dbg_iterdata != null) {
					int dbg_width = clt_parameters.ilp.ilma_num_corr * macroTilesX;
					for (int kk = 0; kk < dbg_iterdata.length; kk++) {
						for (int ii = 0; ii < macroTilesY; ii++) {
							System.arraycopy(dbg_img[kk], // null pointer
									ii * dbg_width,
									dbg_iterdata[kk],
									(i * macroTilesY + ii) * dbg_width,
									dbg_width);
						}
					}
				}
				if (lma_rms[0] > maximal_series_rms) {
					maximal_series_rms = lma_rms[0];
				}
				dbg_rms_pre[i] = lma_rms[0];
				ers_reference.addScene(scene_QuadClt.getImageName(),
						scenes_xyzatr[i][0],
						scenes_xyzatr[i][1],
						ers_scene.getErsXYZ_dt(),		
						ers_scene.getErsATR_dt()		
						);

				// restore original ers data
				ers_scene.setErsDt(
						ers_scene_original_xyz_dt, // double []    ers_xyz_dt,
						ers_scene_original_atr_dt); // double []    ers_atr_dt)(ers_scene_original_xyz_dt);
				ers_scene.setupERS();
				if (debug_level > -1) {
					System.out.println("Pass multi scene "+i+" (of "+ scenes.length+") "+
							reference_QuadClt.getImageName() + "/" + scene_QuadClt.getImageName()+
							" Done. RMS="+lma_rms[0]+", maximal so far was "+maximal_series_rms);
				}
				if (delete_scene_asap) {
					scenes[i+1] = null;
				}
			}
		}
		
		if (dbg_iterdata != null) {
			(new ShowDoubleFloatArrays()).showArrays(
					dbg_iterdata,
					clt_parameters.ilp.ilma_num_corr * macroTilesX,
					macroTilesY * scenes.length,
					true,
					"first_pass_"+scenes[scenes.length-1].getImageName(),
					new String[] {"pX","pY","disparity"}); //	dsrbg_titles);
		}
		// process scenes after reference (if it is not the last
		if (ref_index <  (scenes.length - 1)) {
			for (int i =  ref_index + 1; i < scenes.length ; i++) {
				QuadCLT scene_QuadClt =      scenes[i];
				String last_known_ts =           scenes[i-1].getImageName(); // it should be present in the reference scene scenes
				ErsCorrection ers_scene =            scene_QuadClt.getErsCorrection();

				double [] last_known_xyz = ers_reference.getSceneXYZ(last_known_ts);
				double [] last_known_atr = ers_reference.getSceneATR(last_known_ts);
				double [] last_from_new_xyz = ers_scene.getSceneXYZ(last_known_ts);
				double [] last_from_new_atr = ers_scene.getSceneATR(last_known_ts);
				double [][] new_from_last =   ErsCorrection.invertXYZATR(last_from_new_xyz, last_from_new_atr);

				double [] new_from_last_xyz = new_from_last[0]; // ers_scene_last_known.getSceneXYZ(scene_ts);
				double [] new_from_last_atr = new_from_last[1]; // ers_scene_last_known.getSceneATR(scene_ts);
				double [][] combo_XYZATR = new_from_last;
				System.out.println("Processing scene "+i+": "+scene_QuadClt.getImageName());
				if (i > ( ref_index + 1)) {
					// combine two rotations and two translations 
					combo_XYZATR = ErsCorrection.combineXYZATR(
							last_known_xyz,     // double [] reference_xyz,
							last_known_atr,     // double [] reference_atr, // null?
							new_from_last_xyz,  // double [] scene_xyz,
							new_from_last_atr); // double [] scene_atr)
				}
				// before adjusting - save original ERS, restore afterwards
				double [] ers_scene_original_xyz_dt = ers_scene.getErsXYZ_dt();
				double [] ers_scene_original_atr_dt = ers_scene.getErsATR_dt();

				// ers should be correct for both

				double [] lma_rms = new double[2];
				if (high_res_motion_vectors) {
					scenes_xyzatr[i] = adjustPairsLMAInterscene(
							clt_parameters,                                 // CLTParameters  clt_parameters,
							reference_QuadClt,                              // QuadCLT reference_QuadCLT,
							null, // double []        ref_disparity, // null or alternative reference disparity
							null, // boolean []     reliable_ref, // null or bitmask of reliable reference tiles
							scene_QuadClt,                                  // QuadCLT scene_QuadCLT,
							combo_XYZATR[0],                                // xyz
							combo_XYZATR[1],                                // atr
							param_select2,                                  // final boolean[]   param_select,
							param_regweights2,                              //  final double []   param_regweights,
							lma_rms,                                        // double []      rms, // null or double [2]
							clt_parameters.imp.max_rms,                     // double         max_rms,
							clt_parameters.imp.debug_level);                // 1); // -1); // int debug_level);
				} else {
					scenes_xyzatr[i] = adjustPairsLMA(
							clt_parameters,     // CLTParameters  clt_parameters,			
							reference_QuadClt, // QuadCLT reference_QuadCLT,
							scene_QuadClt, // QuadCLT scene_QuadCLT,
							combo_XYZATR[0], // xyz
							combo_XYZATR[1], // atr
							param_select2,             // final boolean[]   param_select,
							param_regweights2, //  final double []   param_regweights,
							lma_rms, // 			double []      rms, // null or double [2]
							null, // double [][]    dbg_img,
							0.0, // double         max_rms,
							debug_level); // int debug_level)
				}
			    if (lma_rms[0] > maximal_series_rms) {
			        maximal_series_rms = lma_rms[0];
			    }
				dbg_rms_pre[i] = lma_rms[0];
				ers_reference.addScene(scene_QuadClt.getImageName(),
						scenes_xyzatr[i][0],
						scenes_xyzatr[i][1],
						ers_scene.getErsXYZ_dt(),		
						ers_scene.getErsATR_dt()		
						);

				// restore original ers data
				ers_scene.setErsDt(
						ers_scene_original_xyz_dt, // double []    ers_xyz_dt,
						ers_scene_original_atr_dt); // double []    ers_atr_dt)(ers_scene_original_xyz_dt);
				ers_scene.setupERS();
				if (debug_level > -1) {
					System.out.println("Pass multi scene "+i+" (of "+ scenes.length+") "+
			                reference_QuadClt.getImageName() + "/" + scene_QuadClt.getImageName()+
			                " Done. RMS="+lma_rms[0]+", maximal so far was "+maximal_series_rms);
				}
				if (delete_scene_asap) {
					scenes[i-1] = null;
				}
			}		
		}
		if (debug_level > -3) {
			System.out.println("All multi scene passes are Done. Maximal RMSE was "+maximal_series_rms);
		}
		
//		boolean show_results = true;
		if (show_results) {
			int dbg_w = scenes_xyzatr.length;
			int dbg_h = 6;
			double [] dbg_img = new double [dbg_w * dbg_h];
			Arrays.fill(dbg_img, Double.NaN);
			for (int dh = 0; dh  < dbg_h; dh++) {
				for (int dw = 0; dw  < dbg_w; dw++) {
					if (scenes_xyzatr[dw] != null) {
						dbg_img[dh*dbg_w + dw] = scenes_xyzatr[dw][dh / 3][dh % 3];
					}
				}
			}
			(new ShowDoubleFloatArrays()).showArrays(
					dbg_img,
					dbg_w,
					dbg_h,
					"scenes_xyzatr"); //	dsrbg_titles);
		}
		
		double [][][] scenes_xyzatr_preadjust = new double [scenes.length][][]; // debug only
		double [][][] scenes_xyzatr_dt_preadjust = new double [scenes.length][][]; // debug only
		double [] dbg_rms = new double [scenes.length];
		if ((clt_parameters.ofp.lpf_series > 1.0) && pattern_mode) { // no velocities in pattern mode
			System.out.println ("adjustSeries(): No processing velosities (\"Second pass\" in pattern_mode");
		}
		if ((clt_parameters.ofp.lpf_series > 1.0) && !pattern_mode) { // no velocities in pattern mode
			boolean ers_invert = false; // true; // false;
			// get current ers
			double [][][]ers_xyzatr =  getVelocitiesFromScenes(
					 scenes, // QuadCLT []     scenes, // ordered by increasing timestamps
					 scenes_xyzatr,                    // double [][][]  scenes_xyzatr,
					 clt_parameters.ofp.lpf_series);   // double         half_run_range
			if (ers_invert) {
				for (int nscene = 0; nscene < ers_xyzatr.length; nscene++) if (ers_xyzatr[nscene] != null) {
					for (int m= 0; m < ers_xyzatr[nscene].length; m++) {
						for (int n= 0; n < ers_xyzatr[nscene][m].length; n++) {
							ers_xyzatr[nscene][m][n] = - ers_xyzatr[nscene][m][n]; 
						}
					}
				}
			}
			
			if (show_results) {
				int dbg_w = ers_xyzatr.length;
				int dbg_h = 6;
				double [] dbg_img = new double [dbg_w * dbg_h];
				Arrays.fill(dbg_img, Double.NaN);
				for (int dh = 0; dh  < dbg_h; dh++) {
					for (int dw = 0; dw  < dbg_w; dw++) {
						if (ers_xyzatr[dw] != null) {
							dbg_img[dh*dbg_w + dw] = ers_xyzatr[dw][dh / 3][dh % 3];
						}
					}
				}
				(new ShowDoubleFloatArrays()).showArrays(
						dbg_img,
						dbg_w,
						dbg_h,
						"scenes_xyzatr-ers_lpf"+clt_parameters.ofp.lpf_series); //	dsrbg_titles);
			}
			// Set reference frame ERS
			ers_reference.setErsDt_test(// scaled/sign
					(clt_parameters.ilp.ilma_ignore_ers ? (new double[3]) : (ers_xyzatr[ref_index][00])),  // double []    ers_xyz_dt,
					(clt_parameters.ilp.ilma_ignore_ers ? (new double[3]) : (ers_xyzatr[ref_index][1]))); // double []    ers_atr_dt)(ers_scene_original_xyz_dt);
			ers_reference.setupERS();
			Arrays.fill(dbg_rms, Double.NaN);
			for (int i =  00; i < scenes.length ; i++) if (i != ref_index){
				String scene_ts =          scenes[i].getImageName(); // it should be present in the scenes[i+1] scenes
				double [] scene_xyz =      ers_reference.getSceneXYZ(scene_ts);
				double [] scene_atr =      ers_reference.getSceneATR(scene_ts);
				ErsCorrection ers_scene =  scenes[i].getErsCorrection();
				ers_scene.setErsDt_test(// scaled/sign
						(clt_parameters.ilp.ilma_ignore_ers ? (new double[3]) : (ers_xyzatr[i][0])), // double []    ers_xyz_dt,
						(clt_parameters.ilp.ilma_ignore_ers ? (new double[3]) : (ers_xyzatr[i][1])));  // double []    ers_atr_dt)(ers_scene_original_xyz_dt);
				ers_scene.setupERS();

				// Debug only
				scenes_xyzatr_dt_preadjust[i] = new double[][] {ers_scene.getErsXYZ_dt().clone(),ers_scene.getErsATR_dt().clone()};
				scenes_xyzatr_preadjust[i] =    new double[][] {scene_xyz.clone(), scene_atr.clone()};
				double [] lma_rms = new double[2];
				
				if (high_res_motion_vectors) {
					scenes_xyzatr[i] = adjustPairsLMAInterscene(
							clt_parameters,                  // CLTParameters  clt_parameters,
							reference_QuadClt,               // QuadCLT reference_QuadCLT,
							null, // double []        ref_disparity, // null or alternative reference disparity 
							null, // boolean []     reliable_ref, // null or bitmask of reliable reference tiles
							scenes[i],                       // QuadCLT scene_QuadCLT,
							scene_xyz, // combo_XYZATR[0],   // xyz
							scene_atr, // combo_XYZATR[1],   // atr
							param_select3, // 3, // 2,       // final boolean[]   param_select,
							param_regweights2,               // final double []   param_regweights,
							lma_rms,                         //	double []      rms, // null or double [2]
							clt_parameters.imp.max_rms,      // double         max_rms,
							clt_parameters.imp.debug_level); // 1); // -1); // int debug_level);
				} else {
					scenes_xyzatr[i] = adjustPairsLMA(
							clt_parameters,                  // CLTParameters  clt_parameters,			
							reference_QuadClt,               // QuadCLT reference_QuadCLT,
							scenes[i],                       // QuadCLT scene_QuadCLT,
							scene_xyz,                       // combo_XYZATR[0], // xyz
							scene_atr,                       // combo_XYZATR[1], // atr
							param_select3,                   // final boolean[]   param_select,
							param_regweights2,               // final double []   param_regweights,
							lma_rms,                         //	double []      rms, // null or double [2]
							null,                            // double [][]    dbg_img,
							0.0,                             // double         max_rms,
							debug_level);                    // int debug_level)
				}
			    if (lma_rms[0] > maximal_series_rms1) {
			        maximal_series_rms1 = lma_rms[0];
			    }
				dbg_rms[i] = lma_rms[0];
				ers_reference.addScene(
						scene_ts,
						scenes_xyzatr[i][0],
						scenes_xyzatr[i][1],
						ers_scene.getErsXYZ_dt(),
						ers_scene.getErsATR_dt()		
						);
				
/* Was it just to undo LMA?
				// restore original ers data
				ers_scene.setErsDt(
						ers_scene_original_xyz_dt, // double []    ers_xyz_dt,
						ers_scene_original_atr_dt); // double []    ers_atr_dt)(ers_scene_original_xyz_dt);
				ers_scene.setupERS();
*/				
				if (debug_level > -1) {
					System.out.println("Pass 2 multi scene "+i+" (of "+ scenes.length+") "+
							reference_QuadClt.getImageName() + "/" + scene_ts+
							" Done. RMS="+lma_rms[0]+", maximal so far was "+maximal_series_rms);
				}
				if (delete_scene_asap) {
					scenes[i-1] = null;
				}
			}		
		}
		double rms_pre_mean = 0.0, rms_mean = 0.0;
		int    rms_pre_num =  0, rms_num =    0;
		for (int i = 0; (i < dbg_rms_pre.length) && (i < dbg_rms.length); i++) {
			if (!Double.isNaN(dbg_rms_pre[i]) && !Double.isNaN(dbg_rms[i])) {
				rms_pre_num++;
				rms_num++;
				rms_pre_mean +=dbg_rms_pre[i];
				rms_mean +=    dbg_rms[i];
			}
		}
		rms_pre_mean /= rms_pre_num;
		rms_mean /= rms_num;
		if (debug_level > -3) {
			System.out.println("adjustSeries() rms_pre_mean="+rms_pre_mean+", rms_mean="+
					rms_mean+", maximal RMSE="+maximal_series_rms1);
		}

		if (show_results) {
			int dbg_w = scenes_xyzatr_preadjust.length;
			int dbg_h = 6;
			double [] dbg_img = new double [dbg_w * (dbg_h + 1)];
			Arrays.fill(dbg_img, Double.NaN);
			for (int dh = 0; dh  < dbg_h; dh++) {
				for (int dw = 0; dw  < dbg_w; dw++) {
					if (scenes_xyzatr_preadjust[dw] != null) {
						dbg_img[dh*dbg_w + dw] = scenes_xyzatr_preadjust[dw][dh / 3][dh % 3];
					}
				}
			}
			for (int dw = 0; dw  < dbg_w; dw++) {
				dbg_img[dbg_h*dbg_w + dw] = dbg_rms_pre[dw];
			}
			(new ShowDoubleFloatArrays()).showArrays(
					dbg_img,
					dbg_w,
					dbg_h + 1,
					"scenes_xyzatr_preadjust"); //	dsrbg_titles);
			System.out.println("adjustSeries(): dbg-2.5");
		}
		
		if (show_results) {
			int dbg_w = scenes_xyzatr_dt_preadjust.length;
			int dbg_h = 6;
			double [] dbg_img = new double [dbg_w * dbg_h];
			Arrays.fill(dbg_img, Double.NaN);
			for (int dh = 0; dh  < dbg_h; dh++) {
				for (int dw = 0; dw  < dbg_w; dw++) {
					if (scenes_xyzatr_dt_preadjust[dw] != null) {
						dbg_img[dh*dbg_w + dw] = scenes_xyzatr_dt_preadjust[dw][dh / 3][dh % 3];
					}
				}
			}
			(new ShowDoubleFloatArrays()).showArrays(//
					dbg_img,
					dbg_w,
					dbg_h,
					"ers_preadjust-2"); //	dsrbg_titles);
			System.out.println("adjustSeries(): dbg-3");
		}

		if (show_results) {
			double [][][] xyzqtr_current = new double [scenes.length][][];
			double [] scene_xyz;
			double [] scene_atr;
			for (int nscene = 0; nscene < scenes.length; nscene++) {
				if (nscene == ref_index) {
					scene_xyz = ers_reference.getCameraXYZ();
					scene_atr = ers_reference.getCameraATR();
				} else {
					String sts = scenes[nscene].getImageName();
					scene_xyz = ers_reference.getSceneXYZ(sts); 
					scene_atr = ers_reference.getSceneATR(sts); 
				}
				if ((scene_xyz != null) &&(scene_atr != null)) {
					xyzqtr_current[nscene] = new double[][] {scene_xyz, scene_atr};
				} else {
					System.out.println("adjustSeries(): null for nscene="+nscene);
				}
			}
			
			int dbg_w = xyzqtr_current.length;
			int dbg_h = 6;
			double [] dbg_img = new double [dbg_w * (dbg_h + 1)];
			Arrays.fill(dbg_img, Double.NaN);
			for (int dh = 0; dh  < dbg_h; dh++) {
				for (int dw = 0; dw  < dbg_w; dw++) {
					if (xyzqtr_current[dw] != null) {
						dbg_img[dh*dbg_w + dw] = xyzqtr_current[dw][dh / 3][dh % 3];
					}
				}
			}
			for (int dw = 0; dw  < dbg_w; dw++) {
				dbg_img[dbg_h*dbg_w + dw] = dbg_rms[dw];
			}

			//dbg_rms
			(new ShowDoubleFloatArrays()).showArrays(
					dbg_img,
					dbg_w,
					dbg_h+1,
					"scenes_xyzatr_adjust-ers"); //	dsrbg_titles);
			System.out.println("adjustSeries(): dbg-2.1");
		}
		
		// get add show current ERS, including reference
		if (show_results) {
			double [][][] ers_current = new double [scenes.length][][];
			double [] scene_xyz_dt;
			double [] scene_atr_dt;
			for (int nscene = 0; nscene < scenes.length; nscene++) {
				if (nscene == ref_index) {
					scene_xyz_dt = ers_reference.getErsXYZ_dt();
					scene_atr_dt = ers_reference.getErsATR_dt();
				} else {
					String sts = scenes[nscene].getImageName();
					scene_xyz_dt = ers_reference.getSceneErsXYZ_dt(sts); 
					scene_atr_dt = ers_reference.getSceneErsATR_dt(sts); 
				}
				if ((scene_xyz_dt != null) &&(scene_atr_dt != null)) {
					ers_current[nscene] = new double[][] {scene_xyz_dt, scene_atr_dt};
				} else {
					System.out.println("adjustSeries(): null for nscene="+nscene);
				}
			}
			
			int dbg_w = ers_current.length;
			int dbg_h = 6;
			double [] dbg_img = new double [dbg_w * dbg_h];
			Arrays.fill(dbg_img, Double.NaN);
			for (int dh = 0; dh  < dbg_h; dh++) {
				for (int dw = 0; dw  < dbg_w; dw++) {
					if (ers_current[dw] != null) {
						dbg_img[dh*dbg_w + dw] = ers_current[dw][dh / 3][dh % 3];
					}
				}
			}
			(new ShowDoubleFloatArrays()).showArrays( //
					dbg_img,
					dbg_w,
					dbg_h,
					"ers_adjust-ers"); //	dsrbg_titles);
			System.out.println("adjustSeries(): dbg-2");
		}
		
// getVelocitiesFromScenes	
//		
		reference_QuadClt.saveInterProperties( // save properties for interscene processing (extrinsics, ers, ...)
	            null, // String path,             // full name with extension or w/o path to use x3d directory
	            debug_level+1);
		
		if (!delete_scene_asap && (debug_level > 1000)) { // Null pointer at compareRefSceneTiles ->ers_scene.setupERS();-> Quaternion quat_center1 = ...
			System.out.println("adjustSeries(): preparing image set...");
			int nscenes = scenes.length;
			int indx_ref = nscenes - 1; 
			double [][][] all_scenes_xyzatr = new double [scenes.length][][]; // includes reference (last)
			double [][][] all_scenes_ers_dt = new double [scenes.length][][]; // includes reference (last)
			all_scenes_xyzatr[indx_ref] = new double [][] {ZERO3,ZERO3};
			all_scenes_ers_dt[indx_ref] = new double [][] {
				ers_reference.getErsXYZ_dt(),
				ers_reference.getErsATR_dt()};
				for (int i = 0; i < nscenes; i++) if (i != indx_ref) {
					String ts = scenes[i].getImageName();
					all_scenes_xyzatr[i] = new double[][] {ers_reference.getSceneXYZ(ts),       ers_reference.getSceneATR(ts)}; 		
					all_scenes_ers_dt[i] = new double[][] {ers_reference.getSceneErsXYZ_dt(ts), ers_reference.getSceneErsATR_dt(ts)}; 		
				}
				compareRefSceneTiles(
						"" ,               // String suffix,
						false,             // boolean blur_reference,
						all_scenes_xyzatr, // double [][][] scene_xyzatr, // does not include reference
						all_scenes_ers_dt, // double [][][] scene_ers_dt, // does not include reference
						scenes,            // QuadCLT [] scenes,
						8);                // int iscale) // 8
		}		
		if (debug_level > -3) {
			System.out.println("adjustSeries() Done. Maximal RMSE in pass1 was "+maximal_series_rms+
					", in pass2 - "+maximal_series_rms1);
		}
	}
	
	/**
	 * Low-pass filter per-scene linear and angular velocities to be used for ERS 
	 * @param ers_xyzatr_in scene linear and angular velocities before filtering [scene][2][3]
	 * @param half_run_range number of scenes each way to average (0 - none,
	 *        1 - 3, 2 - 5) reduced when near the limits
	 * @return [scene][2][3] xyz/dt, atr/dt. atr[0] and atr[2] are inverted relative to input
	 */
	public static double [][][] LPFVelocities(
			double [][][]  ers_xyzatr_in,
			double         half_run_range
			){
		double [] weights = new double [(int) Math.floor(half_run_range)+1];
		weights[0] = 1;
		for (int i = 1; i < weights.length; i++) {
			weights[i] = 0.5* (Math.cos(i*Math.PI/half_run_range) + 1.0);
		}
		double [][][] ers_xyzatr = new double [ers_xyzatr_in.length][][];
		for (int nscene = 0; nscene < ers_xyzatr_in.length; nscene ++) {
			double sw = 0.0;
			double [][] swd =  new double[2][3]; 
			for (int ds = -weights.length + 1; ds < weights.length; ds++) {
				int ns = nscene + ds;
				if ((ns >= 0) && (ns < ers_xyzatr_in.length) && (ers_xyzatr_in[ns] != null)) {
					double w = (ds >= 0) ? weights[ds] : weights[-ds];
					sw += w;
					for (int m = 0; m < swd.length; m++) {
						for (int d = 0; d < swd[m].length; d++) {
							swd[m][d] += w * ers_xyzatr_in[ns][m][d];
						}
					}
				}
			}
			if (sw > 0) {
				ers_xyzatr[nscene] = new double[2][3];
				for (int m = 0; m < swd.length; m++) {
					for (int d = 0; d < swd[m].length; d++) {
						ers_xyzatr[nscene][m][d] = swd[m][d]/sw;
					}
				}
			}
		}
		return ers_xyzatr;
	}
	
	
	
	/**
	 * Calculate linear and angular velocities by running-average of the poses 
	 * @param scenes Scene objects (to get timestamps)
	 * @param scenes_xyzatr scene linear and angular coordinate [scene][2][3]
	 * @param half_run_range number of scenes each way to average (0 - none,
	 *        1 - 3, 2 - 5) reduced when near the limits
	 * @return [scene][2][3] xyz/dt, atr/dt. atr[0] and atr[2] are inverted relative to input
	 */
	public static double [][][] getVelocitiesFromScenes(
			QuadCLT []     scenes, // ordered by increasing timestamps
			double [][][]  scenes_xyzatr,
			double         half_run_range
			){
		double [] weights = new double [(int) Math.floor(half_run_range)+1];
		weights[0] = 1;
		for (int i = 1; i < weights.length; i++) {
			weights[i] = 0.5* (Math.cos(i*Math.PI/half_run_range) +1.0);
		}
		double [][][] ers_xyzatr = new double [scenes_xyzatr.length][][];
		for (int nscene = 0; nscene < scenes_xyzatr.length; nscene ++) if (scenes[nscene] != null){
			double tref = scenes[nscene].getTimeStamp();
			double s0 = 0.0, sx = 0.0, sx2 = 0.0;
			double [][] sy =  new double[2][3]; 
			double [][] sxy = new double[2][3];
			int nds = 0;
			for (int ds = -weights.length + 1; ds < weights.length; ds++) {
				int ns = nscene + ds;
				if ((ns >= 0) && (ns < scenes_xyzatr.length) && (scenes_xyzatr[ns] != null)) {
					double w = (ds >= 0) ? weights[ds] : weights[-ds];
					s0 += w;
					double dt = scenes[ns].getTimeStamp() - tref;
					double wx = w * dt;
					sx +=  wx;
					sx2 += wx * dt;
					for (int m = 0; m < sy.length; m++) {
						for (int d = 0; d < sy[m].length; d++) {
							double y = scenes_xyzatr[ns][m][d];
							sy [m][d] += w *  y;
							sxy[m][d] += wx * y;
						}
					}
					nds++; // number of different timestamps used
				}
			}
			if ((nds > 1) && (s0 > 0)) {
				ers_xyzatr[nscene] = new double[2][3];
				for (int m = 0; m < sy.length; m++) {
					for (int d = 0; d < sy[m].length; d++) {
						ers_xyzatr[nscene][m][d] = (s0*sxy[m][d]-sx*sy[m][d])/(s0*sx2 - sx*sx);
					}
				}
				ers_xyzatr[nscene][1][0] *= -1.0;
				ers_xyzatr[nscene][1][2] *= -1.0;
			}
		}
		return ers_xyzatr;
	}
	
	public void IntersceneAccumulate(
			CLTParameters        clt_parameters,
			ColorProcParameters  colorProcParameters,
			QuadCLT.SetChannels [] set_channels,
			QuadCLT              ref_scene, // ordered by increasing timestamps
			NoiseParameters noise_sigma_level, // only comes with no-noise here
			int                  debug_level
			)
	{
		System.out.println("IntersceneAccumulate(), scene timestamp="+ref_scene.getImageName());
		ErsCorrection ers_reference = ref_scene.getErsCorrection();
		String [] sts = ers_reference.getScenes(); // should have all scenes selected, last being reference
		ArrayList<String> sts_list = new ArrayList<String>();
		ArrayList<String> scenes_list = new ArrayList<String>();
		for (String s:sts) sts_list.add(s);
		for (int i = 0; i < set_channels.length; i++) {
			String scene_name = set_channels[i].name();
			if (sts_list.contains(scene_name)) {
				scenes_list.add(scene_name);
			}
		}
		String [] scene_names = scenes_list.toArray(new String [0]); 
		// get list of all other scenes
		int num_scenes = scene_names.length + 1;
		int indx_ref = num_scenes - 1; 
		QuadCLT [] scenes = new QuadCLT [num_scenes];
		scenes[indx_ref] = ref_scene;
		
		for (int i = 0; i < scene_names.length; i++) {
			scenes[i] = ref_scene.spawnQuadCLTWithNoise( // spawnQuadCLT(
					scene_names[i],
					clt_parameters,
					colorProcParameters, //
					noise_sigma_level,   // double []            noise_sigma_level,only comes with non-noise here, so noise_variant is not needed
					-1,                  // int                  noise_variant, // <0 - no-variants, compatible with old code
					ref_scene,           // QuadCLTCPU           ref_scene, // may be null if scale_fpn <= 0
					threadsMax,
					-1); // debug_level);
			scenes[i].setDSRBG(
					clt_parameters, // CLTParameters  clt_parameters,
					threadsMax,     // int            threadsMax,  // maximal number of threads to launch
					updateStatus,   // boolean        updateStatus,
					-1); // debug_level);    // int            debugLevel)
		}
		
		final double [][] combo_dsn =  prepareInitialComboDS(
				clt_parameters,   // final CLTParameters       clt_parameters,
				scenes,           // final QuadCLT []          scenes,
				indx_ref,         // final int                 indx_ref,
				debug_level-2);     // final int                 debug_level);
		// re-read ref's INTER if possible
		
		if (ref_scene.restoreDSI( // if there is already calculated with interscene - us it
				"-DSI_INTER",
				true // silent
				) >=0 ) {
			System.out.println("IntersceneAccumulate(): Using previously calculated interscene DSI (*-DSI_INTER) as initial DSI");
			combo_dsn[0] = ref_scene.dsi[ref_scene.is_aux?TwoQuadCLT.DSI_DISPARITY_AUX:TwoQuadCLT.DSI_DISPARITY_MAIN];
			combo_dsn[1] = ref_scene.dsi[ref_scene.is_aux?TwoQuadCLT.DSI_STRENGTH_AUX:TwoQuadCLT.DSI_STRENGTH_MAIN];
		}
		
		final double [][] combo_dsn_change = new double [combo_dsn.length+1][];
		for (int i = 0; i < combo_dsn.length; i++) {
			combo_dsn_change[i] = combo_dsn[i];
		}

		final int margin = 8;
		final int tilesX = ref_scene.getTileProcessor().getTilesX();
		final int tilesY = ref_scene.getTileProcessor().getTilesY();

		String [] combo_dsn_titles = {"disp", "strength", "num_valid","change"};
	
		if (debug_level >-3) {
			(new ShowDoubleFloatArrays()).showArrays(
					combo_dsn_change,
					tilesX,
					tilesY,
					true,
					"combo_dsn-initial"+ref_scene.getImageName(),
					combo_dsn_titles); //	dsrbg_titles);
		}
		
		final int max_refines = 10;
		final String [] iter_titles = {"disp", "diff", "strength"};
		final int [] iter_indices = {0,1,3};
		final int last_slices = combo_dsn_titles.length;
		final int last_initial_slices = last_slices + iter_titles.length;
		
		final double [][] refine_results = new double [last_slices + 3 * (max_refines + 1)][];
		String [] refine_titles = new String [refine_results.length];
		for (int i = 0; i < combo_dsn_titles.length; i++) {
			refine_results[i] = combo_dsn_change[i];
			refine_titles[i] = combo_dsn_titles[i]+"-last";
		}
		for (int i = 0; i < iter_titles.length; i++) {
			refine_titles[last_slices + i] = iter_titles[i]+"-initial";
			if (combo_dsn_change[iter_indices[i]] != null) {
				refine_results[last_slices + i] = combo_dsn_change[iter_indices[i]].clone();
			} else {
				refine_results[last_slices + i] = new double [tilesX * tilesY];
			}
		}
		for (int nrefine = 0; nrefine < max_refines; nrefine++) {
			for (int i = 0; i < iter_titles.length; i++) {
				refine_titles[last_initial_slices + i * max_refines + nrefine ] = combo_dsn_titles[i]+"-"+nrefine;
			}
		}		
		double [][] disparity_map = null;
		for (int nrefine = 0; nrefine < max_refines; nrefine++) {
	///		Runtime.getRuntime().gc();
	///		System.out.println("--- Free memory="+Runtime.getRuntime().freeMemory()+" (of "+Runtime.getRuntime().totalMemory()+")");
			int mcorr_sel = Correlation2d.corrSelEncode(clt_parameters.img_dtt,scenes[indx_ref].getNumSensors());
			
			disparity_map = correlateInterscene(
					clt_parameters, // final CLTParameters  clt_parameters,
					scenes,         // final QuadCLT []     scenes,
					indx_ref,       // final int            indx_ref,
					combo_dsn_change[0],   // final double []      disparity_ref,  // disparity in the reference view tiles (Double.NaN - invalid)
					null,           // final boolean []     selection, // may be null, if not null do not  process unselected tiles
					margin,         // final int            margin,
					nrefine,        // final int            nrefine, // just for debug title
					clt_parameters.inp.show_final_2d, // final boolean        show_2d_corr,
					mcorr_sel,      // final int            mcorr_sel, //  = Correlation2d.corrSelEncode(clt_parameters.img_dtt,scenes[nscene].getNumSensors());
					null,           // final float [][][]   accum_2d_corr, // if [1][][] - return accumulated 2d correlations (all pairs)
					false,          // final boolean        no_map, // do not generate disparity_map (time-consuming LMA)
					debug_level-5);   // final int            debug_level)

	///		Runtime.getRuntime().gc();
	///		System.out.println("--- Free memory="+Runtime.getRuntime().freeMemory()+" (of "+Runtime.getRuntime().totalMemory()+")");

			(new ShowDoubleFloatArrays()).showArrays(
					disparity_map,
					tilesX,
					tilesY,
					true,
					"accumulated_disparity_map-"+nrefine,
					ImageDtt.getDisparityTitles(ref_scene.getNumSensors(),ref_scene.isMonochrome()) // ImageDtt.DISPARITY_TITLES
					);
			// update disparities
			final int disparity_index = ImageDtt.DISPARITY_INDEX_CM; // 2
			final int strength_index =  ImageDtt.DISPARITY_STRENGTH_INDEX; // 10
			for (int nTile =0; nTile < combo_dsn_change[0].length; nTile++) {
				if (!Double.isNaN(combo_dsn_change[0][nTile]) && !Double.isNaN(disparity_map[disparity_index][nTile])) {
					combo_dsn_change[0][nTile] += disparity_map[disparity_index][nTile];
					combo_dsn_change[1][nTile]  = disparity_map[strength_index][nTile];
				}
			}
			combo_dsn_change[combo_dsn_change.length -1] = disparity_map[disparity_index]; 
			
			for (int i = 0; i < iter_titles.length; i++) {
				refine_results[last_initial_slices + (i * max_refines) + nrefine] = combo_dsn_change[iter_indices[i]].clone();
			}
			
			if (debug_level >-3) {
				(new ShowDoubleFloatArrays()).showArrays(
						combo_dsn_change,
						tilesX,
						tilesY,
						true,
						"combo_dsn-"+nrefine+"-"+ref_scene.getImageName(),
						combo_dsn_titles); //	dsrbg_titles);
			}
		}
		
		if (debug_level >-5) {
			(new ShowDoubleFloatArrays()).showArrays(
					refine_results,
					tilesX,
					tilesY,
					true,
					"combo-"+max_refines+"-"+ref_scene.getImageName(),
					refine_titles); //	dsrbg_titles);
		}
		ref_scene.saveDoubleArrayInModelDirectory(
				"-results-nonoise",  // String      suffix,
				refine_titles,       // null,          // String []   labels, // or null
				refine_results,      // dbg_data,         // double [][] data,
				tilesX,              // int         width,
				tilesY);             // int         height)
		
		// save _DSI_INTER - same format, as _DSI_MAIN, it will be used instead of _DSI_MAIN next time
		double [][] dsi = new double [TwoQuadCLT.DSI_SLICES.length][];
		dsi[ref_scene.is_aux?TwoQuadCLT.DSI_DISPARITY_AUX:TwoQuadCLT.DSI_DISPARITY_MAIN] = combo_dsn_change[0];
		dsi[ref_scene.is_aux?TwoQuadCLT.DSI_STRENGTH_AUX:TwoQuadCLT.DSI_STRENGTH_MAIN] =   combo_dsn_change[1];
		double [] disp_lma = disparity_map[ImageDtt.DISPARITY_INDEX_POLY];
		if (disp_lma != null) { 
			int indx_lma = ref_scene.is_aux?TwoQuadCLT.DSI_DISPARITY_AUX_LMA:TwoQuadCLT.DSI_DISPARITY_MAIN_LMA; 
			dsi[indx_lma] = combo_dsn_change[0].clone();
			for (int i = 0; i < disp_lma.length; i++) {
				if (Double.isNaN(disp_lma[i])) {
					dsi[indx_lma][i] = Double.NaN;		
				}
			}
		}
		ref_scene.saveDSIAll ( "-DSI_INTER",dsi); // delete/rename this file to start iterations from reference lma

		// save combo_dsn_change to model directory
		if (debug_level >-100) {
			return;
		}
		
		System.out.println("IntersceneAccumulate(), got previous scenes: "+sts.length);
		if (debug_level > 1) { // tested OK
			System.out.println("IntersceneAccumulate(): preparing image set...");
			int nscenes = scenes.length;
			//			int indx_ref = nscenes - 1; 
			double [][][] all_scenes_xyzatr = new double [scenes.length][][]; // includes reference (last)
			double [][][] all_scenes_ers_dt = new double [scenes.length][][]; // includes reference (last)
			all_scenes_xyzatr[indx_ref] = new double [][] {ZERO3,ZERO3};
			all_scenes_ers_dt[indx_ref] = new double [][] {
				ers_reference.getErsXYZ_dt(),
				ers_reference.getErsATR_dt()};
				
				for (int i = 0; i < nscenes; i++) if (i != indx_ref) {
					String ts = scenes[i].getImageName();
					all_scenes_xyzatr[i] = new double[][] {ers_reference.getSceneXYZ(ts),       ers_reference.getSceneATR(ts)}; 		
					all_scenes_ers_dt[i] = new double[][] {ers_reference.getSceneErsXYZ_dt(ts), ers_reference.getSceneErsATR_dt(ts)}; 		
				}
				compareRefSceneTiles(
						"" ,               // String suffix,
						true, // false,             // boolean blur_reference,
						all_scenes_xyzatr, // double [][][] scene_xyzatr, // does not include reference
						all_scenes_ers_dt, // double [][][] scene_ers_dt, // does not include reference
						scenes,            // QuadCLT [] scenes,
						8);                // int iscale) // 8
		}
		// create initial disparity map for the reference scene
		
	}
	
	public void intersceneNoiseDebug(
			CLTParameters        clt_parameters,
			boolean              ref_only, // process only reference frame (false - inter-scene)
			ColorProcParameters  colorProcParameters,
			QuadCLT              ref_scene, // ordered by increasing timestamps
//			double []            
			NoiseParameters		 noise_sigma_level,
			int                  debug_level
			)
	{
		System.out.println("IntersceneNoise(), scene timestamp="+ref_scene.getImageName());
		ErsCorrection ers_reference = ref_scene.getErsCorrection();
		String [] sts = ref_only ? (new String [0]) : ers_reference.getScenes();
		// get list of all other scenes
		int num_scenes = sts.length + 1;
		int indx_ref = num_scenes - 1; 
		QuadCLT [] scenes = new QuadCLT [num_scenes];
		scenes[indx_ref] = ref_scene;
		
		for (int i = 0; i < sts.length; i++) {
			scenes[i] = ref_scene.spawnQuadCLTWithNoise( // spawnQuadCLT(
					sts[i],
					clt_parameters,
					colorProcParameters, //
					noise_sigma_level,   // double []            noise_sigma_level,
					ref_scene,           // QuadCLTCPU           ref_scene, // may be null if scale_fpn <= 0
					threadsMax,
					-1); // debug_level);
			scenes[i].setDSRBG(
					clt_parameters, // CLTParameters  clt_parameters,
					threadsMax,     // int            threadsMax,  // maximal number of threads to launch
					updateStatus,   // boolean        updateStatus,
					-1); // debug_level);    // int            debugLevel)
		}
		
		double [][] combo_dsn =  null;
		if (noise_sigma_level == null) {
			combo_dsn = prepareInitialComboDS(
					clt_parameters,   // final CLTParameters       clt_parameters,
					scenes,           // final QuadCLT []          scenes,
					indx_ref,         // final int                 indx_ref,
					debug_level-2);     // final int                 debug_level);
		} else {
			combo_dsn = ref_scene. readDoubleArrayFromModelDirectory(
					"-results-nonoise", // String      suffix,
					3, // int         num_slices, // (0 - all)
					null); // int []      wh);

		}
		final double [][] combo_dsn_change = new double [combo_dsn.length+1][];
		for (int i = 0; i < combo_dsn.length; i++) {
			combo_dsn_change[i] = combo_dsn[i];
		}
		final int margin = 8;
		final int tilesX = ref_scene.getTileProcessor().getTilesX();
		final int tilesY = ref_scene.getTileProcessor().getTilesY();

		String [] combo_dsn_titles = {"disp", "strength", "num_valid","change"};
	
		if (debug_level > 0) {
			(new ShowDoubleFloatArrays()).showArrays(
					combo_dsn_change,
					tilesX,
					tilesY,
					true,
					"combo_dsn-initial"+ref_scene.getImageName(),
					combo_dsn_titles); //	dsrbg_titles);
		}
		
		final int max_refines = 10;
		final String [] iter_titles = {"disp", "diff", "strength"};
		final int [] iter_indices = {0,1,3};
		final int last_slices = combo_dsn_titles.length;
		final int last_initial_slices = last_slices + iter_titles.length;
		
		final double [][] refine_results = new double [last_slices + 3 * (max_refines + 1)][];
		String [] refine_titles = new String [refine_results.length];
		for (int i = 0; i < combo_dsn_titles.length; i++) {
			refine_results[i] = combo_dsn_change[i];
			refine_titles[i] = combo_dsn_titles[i]+"-last";
		}
		for (int i = 0; i < iter_titles.length; i++) {
			refine_titles[last_slices + i] = iter_titles[i]+"-initial";
			if (combo_dsn_change[iter_indices[i]] != null) {
				refine_results[last_slices + i] = combo_dsn_change[iter_indices[i]].clone();
			} else {
				refine_results[last_slices + i] = new double [tilesX * tilesY];
			}
		}
		for (int nrefine = 0; nrefine < max_refines; nrefine++) {
			for (int i = 0; i < iter_titles.length; i++) {
				refine_titles[last_initial_slices + i * max_refines + nrefine ] = combo_dsn_titles[iter_indices[i]]+"-"+nrefine;
			}
		}		
		if (noise_sigma_level != null) { // add initial offset to the expected disparity
			for (int i = 0; i < combo_dsn_change[0].length; i++) {
				combo_dsn_change[0][i] += noise_sigma_level.initial_offset; // [2]; 
			}
		}
		for (int nrefine = 0; nrefine < max_refines; nrefine++) {
///			Runtime.getRuntime().gc();
///			System.out.println("--- Free memory="+Runtime.getRuntime().freeMemory()+" (of "+Runtime.getRuntime().totalMemory()+")");

			double [][] disparity_map = 
					//		double [][][][][] clt_corr_partial =
					correlateIntersceneDebug(
							clt_parameters, // final CLTParameters  clt_parameters,
							scenes,         // final QuadCLT []     scenes,
							indx_ref,       // final int            indx_ref,
							combo_dsn_change[0],   // final double []      disparity_ref,  // disparity in the reference view tiles (Double.NaN - invalid)
							margin,         // final int            margin,
							nrefine,        // final int            nrefine, // just for debug title
							debug_level-5);   // final int            debug_level)

///			Runtime.getRuntime().gc();
///			System.out.println("--- Free memory="+Runtime.getRuntime().freeMemory()+" (of "+Runtime.getRuntime().totalMemory()+")");
			if (debug_level > 0) {
				(new ShowDoubleFloatArrays()).showArrays(
						disparity_map,
						tilesX,
						tilesY,
						true,
						"accumulated_disparity_map-"+nrefine,
						ImageDtt.getDisparityTitles(ref_scene.getNumSensors(),ref_scene.isMonochrome()) // ImageDtt.DISPARITY_TITLES
						);
			}
			// update disparities
			final int disparity_index = ImageDtt.DISPARITY_INDEX_CM; // 2
			final int strength_index =  ImageDtt.DISPARITY_STRENGTH_INDEX; // 10
			for (int nTile =0; nTile < combo_dsn_change[0].length; nTile++) {
				if (!Double.isNaN(combo_dsn_change[0][nTile]) && !Double.isNaN(disparity_map[disparity_index][nTile])) {
					combo_dsn_change[0][nTile] += disparity_map[disparity_index][nTile];
					combo_dsn_change[1][nTile]  = disparity_map[strength_index][nTile];
				}
			}
			combo_dsn_change[combo_dsn_change.length -1] = disparity_map[disparity_index]; 
			
			for (int i = 0; i < iter_titles.length; i++) {
				refine_results[last_initial_slices + (i * max_refines) + nrefine] = combo_dsn_change[iter_indices[i]].clone();
			}
			
			if (debug_level >0) {
				(new ShowDoubleFloatArrays()).showArrays(
						combo_dsn_change,
						tilesX,
						tilesY,
						true,
						"combo_dsn-"+nrefine+"-"+ref_scene.getImageName(),
						combo_dsn_titles); //	dsrbg_titles);
			}
		}
		
		if (debug_level > 0) {
			(new ShowDoubleFloatArrays()).showArrays(
					refine_results,
					tilesX,
					tilesY,
					true,
					"combo-"+max_refines+"-"+ref_scene.getImageName(),
					refine_titles); //	dsrbg_titles);
		}
		//noise_sigma_level
		String rslt_suffix = "-results-nonoise";
		if (noise_sigma_level != null) {
//			rslt_suffix = "-results-lev_"+noise_sigma_level[0]+"-sigma_"+noise_sigma_level[1]+"-offset"+noise_sigma_level[2];
			rslt_suffix =
					"-results-rnd_"+noise_sigma_level.scale_random+
					"-fpn_"+        noise_sigma_level.scale_fpn+
					"-sigma_"+      noise_sigma_level.sigma+ // [1]+
					"-offset"+      noise_sigma_level.initial_offset+ // [2];
					"-sensors"+     noise_sigma_level.used_sensors +
					(clt_parameters.correlate_lma?"-lma":"-nolma");

			if (ref_only) {
				rslt_suffix +="-nointer";
			} else {
				rslt_suffix +="-inter";
			}
//			rslt_suffix +="-mask"+clt_parameters.img_dtt.dbg_pair_mask;
		}
		ref_scene.saveDoubleArrayInModelDirectory(
				rslt_suffix,         // String      suffix,
				refine_titles,       // null,          // String []   labels, // or null
				refine_results,      // dbg_data,         // double [][] data,
				tilesX,              // int         width,
				tilesY);             // int         height)
		// save combo_dsn_change to model directory
		if (debug_level >-100) {
			return;
		}
		
		System.out.println("IntersceneAccumulate(), got previous scenes: "+sts.length);
		if (debug_level > 1) { // tested OK
			System.out.println("IntersceneAccumulate(): preparing image set...");
			int nscenes = scenes.length;
			//			int indx_ref = nscenes - 1; 
			double [][][] all_scenes_xyzatr = new double [scenes.length][][]; // includes reference (last)
			double [][][] all_scenes_ers_dt = new double [scenes.length][][]; // includes reference (last)
			all_scenes_xyzatr[indx_ref] = new double [][] {ZERO3,ZERO3};
			all_scenes_ers_dt[indx_ref] = new double [][] {
				ers_reference.getErsXYZ_dt(),
				ers_reference.getErsATR_dt()};
				
				for (int i = 0; i < nscenes; i++) if (i != indx_ref) {
					String ts = scenes[i].getImageName();
					all_scenes_xyzatr[i] = new double[][] {ers_reference.getSceneXYZ(ts),       ers_reference.getSceneATR(ts)}; 		
					all_scenes_ers_dt[i] = new double[][] {ers_reference.getSceneErsXYZ_dt(ts), ers_reference.getSceneErsATR_dt(ts)}; 		
				}
				compareRefSceneTiles(
						"" ,               // String suffix,
						true, // false,             // boolean blur_reference,
						all_scenes_xyzatr, // double [][][] scene_xyzatr, // does not include reference
						all_scenes_ers_dt, // double [][][] scene_ers_dt, // does not include reference
						scenes,            // QuadCLT [] scenes,
						8);                // int iscale) // 8
		}
		// create initial disparity map for the reference scene
		
	}
	
	public static ImagePlus generateSceneOutlines(
			QuadCLT    ref_scene, // ordered by increasing timestamps
			QuadCLT [] scenes,
			int        extra, // add around largest outline
			int        scale,
			int        line_width_outline,
			int        line_width_corners,
			Color      line_color_outline,
			Color      line_color_corners
			) {
		int step_outline = 5;
		int tilesX = ref_scene.getTileProcessor().getTilesX();
		int tilesY = ref_scene.getTileProcessor().getTilesY();
		int transform_size = ref_scene.getTileProcessor().getTileSize();
		ErsCorrection ers_reference = ref_scene.getErsCorrection();
		double [][] corners = {//new double [4][3]; // for each corner {px, py, d}
				{0.00,                         0.0,                         0.0},
				{tilesX * transform_size - 1, 0.0,                         0.0}, 
				{tilesX * transform_size - 1, tilesY * transform_size - 1, 0.0}, 
				{0.0,                         tilesY * transform_size - 1, 0.0}, 
		};
		double [][][] outlines = new double [scenes.length][corners.length][];
		double minx = tilesX * transform_size, maxx = 0.0, miny = tilesY * transform_size, maxy = 0.0; //
		for (int nscene = 0; nscene < scenes.length; nscene++ ) {
			String ts = scenes[nscene].getImageName();
			for (int ncorn = 0; ncorn < corners.length; ncorn++) {
				double []   ref_xyz = ers_reference.getSceneXYZ(ts);
				double []   ref_atr = ers_reference.getSceneATR(ts);

				double [] pXpYD = scenes[nscene].getErsCorrection().getImageCoordinatesERS(
						ref_scene,          // QuadCLT cameraQuadCLT,    // camera station that got image to be to be matched 
						corners[ncorn][0],  // double px,                // pixel coordinate X in the reference view
						corners[ncorn][1],  // double py,                // pixel coordinate Y in the reference view
						corners[ncorn][2],  // double disparity,         // this reference disparity 
						true,               // boolean distortedView,    // This camera view is distorted (diff.rect), false - rectilinear
						ref_xyz,            // double [] reference_xyz,  // this view position in world coordinates (typically zero3)
						ref_atr,            // double [] reference_atr,  // this view orientation relative to world frame  (typically zero3)
						true,               // boolean distortedCamera,  // camera view is distorted (false - rectilinear)
						new double [3], // double [] camera_xyz,     // camera center in world coordinates
						new double [3], // double [] camera_atr,     // camera orientation relative to world frame
						1.0); // double    line_err)       // threshold error in scan lines (1.0)
				outlines[nscene][ncorn] = pXpYD;
				if (minx > pXpYD[0]) minx = pXpYD[0];
				if (miny > pXpYD[1]) miny = pXpYD[1];
				if (maxx < pXpYD[0]) maxx = pXpYD[0];
				if (maxy < pXpYD[1]) maxy = pXpYD[1];
				System.out.println(String.format("%3d:%1d: px= %8.2f py= %8.2f disparity= %8.5f",nscene, ncorn, pXpYD[0], pXpYD[1], pXpYD[2]));
			}
			System.out.println();
		}
		int ix0 =    (int) Math.floor(minx) - extra;  // original pixels, subtract from data
		int width =  scale * (((int) Math.ceil (maxx) - ix0 ) + extra);  // scaled
		int iy0 =    (int) Math.floor(miny) - extra; 
		int height = scale * (((int) Math.ceil (maxy) - iy0 ) + extra);
		ImageStack stack = new 		ImageStack(width, height);
		int [] blank = new int [width*height];
		for (int nscene = 0; nscene < (scenes.length + 2); nscene++ ) {
			String slice_title = "";
			if (nscene < scenes.length) {
				slice_title = scenes[nscene].getImageName();
			} else if (nscene == scenes.length) {
				slice_title = "corners";
			} else {
				slice_title = "all";
			}
			stack.addSlice(slice_title, blank.clone());
		}
		// Draw outline on each individual slice
		for (int nscene = 0; nscene < scenes.length; nscene++ ) {
			ImageProcessor ip = stack.getProcessor(nscene+1);
			ip.setLineWidth(line_width_outline);
			ip.setColor    (line_color_outline); // does not work?
			for (int icorn = 0; icorn < outlines[nscene].length; icorn++) {
				int icorn_next = icorn+1;
				if (icorn_next >= outlines[nscene].length) {
					icorn_next = 0;
				}
				double [][] se = {
						{scale * (outlines[nscene][icorn][0]-ix0),     scale * (outlines[nscene][icorn][1]-iy0)},
						{scale * (outlines[nscene][icorn_next][0]-ix0),scale * (outlines[nscene][icorn_next][1]-iy0)}};
				for (int dy = -line_width_outline+1; dy < line_width_outline; dy+=2) {
					for (int dx = -line_width_outline+1; dx < line_width_outline; dx+=2) {
						Line line = new Line(se[0][0] + 0.5*dx, se[0][1] + 0.5*dy, se[1][0] + 0.5*dx, se[1][1] + 0.5*dy);
						line.drawPixels(ip);
					}
				}
			}
		}
		int indx_corners = scenes.length;
		int indx_all =    indx_corners+1;
		// Trace corners
		{
			ImageProcessor ip = stack.getProcessor(indx_corners+1);
			ip.setLineWidth(line_width_corners);
			ip.setColor    (line_color_corners); // does not work?
			for (int icorn = 0; icorn < outlines[0].length; icorn++) {
				for (int nscene = 0; nscene < (scenes.length - 1); nscene++ ) {
					double [][] se = {
					{scale * (outlines[nscene][icorn][0]-ix0),     scale * (outlines[nscene][icorn][1]-iy0)},
					{scale * (outlines[nscene + 1][icorn][0]-ix0),scale * (outlines[nscene+1][icorn][1]-iy0)}};
					for (int dy = -line_width_corners+1; dy < line_width_corners; dy+=2) {
						for (int dx = -line_width_corners+1; dx < line_width_corners; dx+=2) {
							Line line = new Line(se[0][0] + 0.5*dx, se[0][1] + 0.5*dy, se[1][0] + 0.5*dx, se[1][1] + 0.5*dy);
							line.drawPixels(ip);
						}
					}
				}
			}			
		}
		// Combine all in last slice
		{
			// outlines
			ImageProcessor ip = stack.getProcessor(indx_all+1);
			ip.setLineWidth(line_width_outline);
			ip.setColor    (line_color_outline); // does not work?
			for (int nscene = 0; nscene < scenes.length; nscene++ ) if (((scenes.length - nscene -1) % step_outline) == 0){
				for (int icorn = 0; icorn < outlines[nscene].length; icorn++) {
					int icorn_next = icorn+1;
					if (icorn_next >= outlines[nscene].length) {
						icorn_next = 0;
					}
					double [][] se = {
							{scale * (outlines[nscene][icorn][0]-ix0),     scale * (outlines[nscene][icorn][1]-iy0)},
							{scale * (outlines[nscene][icorn_next][0]-ix0),scale * (outlines[nscene][icorn_next][1]-iy0)}};
					for (int dy = -line_width_outline+1; dy < line_width_outline; dy+=2) {
						for (int dx = -line_width_outline+1; dx < line_width_outline; dx+=2) {
							Line line = new Line(se[0][0] + 0.5*dx, se[0][1] + 0.5*dy, se[1][0] + 0.5*dx, se[1][1] + 0.5*dy);
							line.drawPixels(ip);
						}
					}

				}
			}
			// corners
			ip.setLineWidth(line_width_corners);
			ip.setColor    (line_color_corners); // does not work?
			for (int icorn = 0; icorn < outlines[0].length; icorn++) {
				for (int nscene = 0; nscene < (scenes.length - 1); nscene++ ) {
					double [][] se = {
					{scale * (outlines[nscene][icorn][0]-ix0),     scale * (outlines[nscene][icorn][1]-iy0)},
					{scale * (outlines[nscene + 1][icorn][0]-ix0),scale * (outlines[nscene+1][icorn][1]-iy0)}};
					for (int dy = -line_width_corners+1; dy < line_width_corners; dy+=2) {
						for (int dx = -line_width_corners+1; dx < line_width_corners; dx+=2) {
							Line line = new Line(se[0][0] + 0.5*dx, se[0][1] + 0.5*dy, se[1][0] + 0.5*dx, se[1][1] + 0.5*dy);
							line.drawPixels(ip);
						}
					}
				}
			}			
			
		}
		
        ImagePlus imp_stack = new ImagePlus("scene_outlines", stack);
//        imp_stack.getProcessor().resetMinAndMax();
//       imp_stack.show();
		return imp_stack;
	}
	
	public double [][] intersceneExport(
			CLTParameters        clt_parameters,
			ErsCorrection        ers_reference,
			QuadCLT []           scenes,
			ColorProcParameters  colorProcParameters,
			int                  debug_level
			)
	{
		// empirical correction for both lma and non-lma step
	  	double corr_nonlma = 1.0; // 1.23;
	  	double corr_lma =    1.0; // 1.23;
		// reference scene is always added to the end, even is out of timestamp order
		int      indx_ref = scenes.length - 1; // Always added to the end even if out-of order
		
		QuadCLT  ref_scene = scenes[indx_ref]; // ordered by increasing timestamps
		boolean generate_outlines = false; // true; // TODO: move to configs
		System.out.println("intersceneExport(), scene timestamp="+ref_scene.getImageName());
		int num_scenes = scenes.length;
		String [] combo_dsn_titles_full = COMBO_DSN_TITLES.clone();
		String [] combo_dsn_titles = new String [COMBO_DSN_INDX_DISP_BG];
		for (int i = 0; i < combo_dsn_titles.length; i++) {
			combo_dsn_titles[i] = combo_dsn_titles_full[i]; 
		}
		if (clt_parameters.rig.mll_max_refines_bg <= 0) {
			combo_dsn_titles_full = combo_dsn_titles;
		}
		
		double min_disp_change = clt_parameters.rig.mll_min_disp_change_pre; // 0.001; // stop re-measure when difference is below
		final int max_refines =
				clt_parameters.rig.mll_max_refines_bg  +
				clt_parameters.rig.mll_max_refines_lma +
				clt_parameters.rig.mll_max_refines_pre;
		
		final int [] iter_indices = {
				COMBO_DSN_INDX_DISP,
				COMBO_DSN_INDX_STRENGTH,
				COMBO_DSN_INDX_LMA,
				COMBO_DSN_INDX_CHANGE}; // which to save for each iteration: {"disp", "strength","disp_lma","change"};
		final int [] initial_indices = {
				COMBO_DSN_INDX_DISP,
				COMBO_DSN_INDX_STRENGTH,
				COMBO_DSN_INDX_VALID}; // initial: "disp", "strength","num_valid"
		double [][] combo_dsn =  null;
		// USE reference scene if no individual DSI are available double [][] dsrbg = scene.getDSRBG(); Probably null for non-reference
		// make sure to use already integrated if available
        double [][] combo_dsn0 = scenes[indx_ref].readDoubleArrayFromModelDirectory(
    			"-INTER-INTRA-LMA", // String      suffix,
    			0, // int         num_slices, // (0 - all)
    			null); // int []      wh);
        if (combo_dsn0 == null) {
        	combo_dsn0 = prepareInitialComboDS( // 3
        			clt_parameters,   // final CLTParameters       clt_parameters,
        			scenes,           // final QuadCLT []          scenes,
        			indx_ref,         // final int                 indx_ref,
        			debug_level-2);     // final int                 debug_level);
        } else {
        	combo_dsn0 = conditionComboDsnFinal(
							true,                // boolean        use_conf,       // use configuration parameters, false - use following  
							clt_parameters,      // CLTParameters  clt_parameters,
							combo_dsn0,     // double [][]    combo_dsn_final, // dls,
							scenes[indx_ref], // QuadCLT        scene,
							debug_level); // int            debugLevel);
        }
        combo_dsn = new double[combo_dsn_titles.length - 1][];
        for (int i = 0; i < initial_indices.length; i++) {
        	combo_dsn[initial_indices[i]] = combo_dsn0[i]; // "disp", "strength", <null>, "num_valid"
        }
		
		
		final int margin = 8;
		final int tilesX = ref_scene.getTileProcessor().getTilesX();
		final int tilesY = ref_scene.getTileProcessor().getTilesY();
		final int tiles = tilesX * tilesY;
		// uses 2 GB - change format
		if (generate_outlines) { // debug_level > 100) { // add parameter?
			int        extra = 10; // pixels around largest outline
			int        scale = 4;

			int        line_width_outline = 1;
			int        line_width_corners = 3;
			Color      line_color_outline = new Color(0, 255, 0);   // green
			Color      line_color_corners = new Color(0, 255, 255); // cyan

			// generating scene outlines for results documentation
			ImagePlus imp_outlines = generateSceneOutlines(
					ref_scene,          // QuadCLT    ref_scene, // ordered by increasing timestamps
					scenes,             // QuadCLT [] scenes
					extra,              // int        extra // add around largest outline
					scale,              // int        scale,
					line_width_outline, // int        line_width_outline,
					line_width_corners, // int        line_width_corners,
					line_color_outline, // Color      line_color_outline,
					line_color_corners  // Color      line_color_corners
					);
			ref_scene.saveImagePlusInModelDirectory(
					"scene_outlines", // String      suffix,
					imp_outlines); // ImagePlus   imp)
//			imp_outlines.show();
		}
		final double [][] combo_dsn_change = new double [combo_dsn_titles.length] [tiles];
		for (int i = 0; i < combo_dsn.length; i++) { // 4 elements: "disp", "strength","disp_lma","num_valid"
			if (combo_dsn[i] != null) combo_dsn_change[i] = combo_dsn[i]; // all but change
		}
		if (debug_level > 0) {
			(new ShowDoubleFloatArrays()).showArrays(
					combo_dsn_change,
					tilesX,
					tilesY,
					true,
					"combo_dsn-initial"+ref_scene.getImageName(),
					combo_dsn_titles); //	dsrbg_titles);
		}
		
		final int last_slices = combo_dsn_titles.length;
		final int last_initial_slices = last_slices + initial_indices.length;
		final boolean [] defined_tiles = new boolean [tiles];
		for (int i = 0; i < defined_tiles.length; i++) {
			defined_tiles[i] = !Double.isNaN(combo_dsn_change[COMBO_DSN_INDX_DISP][i]); // initially defined tiles 
		}
		final double [][] refine_results = new double [last_slices + 4 * (max_refines + 1)][];
		String [] refine_titles = new String [refine_results.length];
		for (int i = 0; i < combo_dsn_titles.length; i++) {
			refine_results[i] = combo_dsn_change[i]; // first 5 - references to 5-element combo_dsn_change
			refine_titles[i] =  combo_dsn_titles[i]+"-last"; // "disp", "strength","disp_lma","num_valid","change"
		}
		for (int i = 0; i < initial_indices.length; i++) {
			refine_titles[last_slices + i] = combo_dsn_titles[initial_indices[i]]+"-initial"; // "disp", "strength","num_valid"
			if (combo_dsn_change[initial_indices[i]] != null) {
				refine_results[last_slices + i] = combo_dsn_change[initial_indices[i]].clone();
			} else {
				refine_results[last_slices + i] = new double [tiles];
			}
		}
		for (int nrefine = 0; nrefine < max_refines; nrefine++) {
			for (int i = 0; i < iter_indices.length; i++) {
				refine_titles[last_initial_slices + i * max_refines + nrefine ] = combo_dsn_titles[iter_indices[i]]+"-"+nrefine;
			}
		}
		double [] target_disparity = combo_dsn_change[COMBO_DSN_INDX_DISP].clone();
		double [] target_disparity_orig = target_disparity.clone(); // will just use NaN/not NaN to restore tasks before second pass with LMA
//		double [][] combo_dsn_final = new double [combo_dsn_titles.length][combo_dsn[0].length];
		double [][] combo_dsn_final = 
				new double [(clt_parameters.rig.mll_max_refines_bg > 0)? combo_dsn_titles_full.length:combo_dsn_titles.length][combo_dsn[0].length];
		combo_dsn_final[COMBO_DSN_INDX_DISP]= combo_dsn[COMBO_DSN_INDX_DISP].clone();
		for (int i = 1; i < combo_dsn_final.length; i++) {
			Arrays.fill(combo_dsn_final[i], Double.NaN);
		}
		// Save pair selection and minimize them for scanning, then restore;
		int num_sensors =scenes[indx_ref].getNumSensors();
		int save_pairs_selection = clt_parameters.img_dtt.getMcorr(num_sensors);
		int save_bimax_combine_mode = clt_parameters.img_dtt.bimax_combine_mode;
		boolean save_bimax_dual_only = clt_parameters.img_dtt.bimax_dual_only;
		if (clt_parameters.rig.mll_max_refines_bg > 0) {
			clt_parameters.img_dtt.bimax_combine_mode = Correlation2d.CAMEL_FG; // initially refine FG 
			clt_parameters.img_dtt.bimax_dual_only =    false; // not just dual only
		}
		
		// FIXME: uncomment next
		System.out.println ("++++ Uncomment next line (Done) ++++");
		clt_parameters.img_dtt.setMcorr(num_sensors, 0 ); // remove all
		clt_parameters.img_dtt.setMcorrNeib(num_sensors,true);
		clt_parameters.img_dtt.setMcorrSq  (num_sensors,true); // remove even more?
		clt_parameters.img_dtt.setMcorrDia (num_sensors,true); // remove even more?
		boolean save_run_lma = clt_parameters.correlate_lma;
		clt_parameters.correlate_lma = false;
		double []  disp_err =     new double [combo_dsn[0].length];  // last time measure disparity error from correlation (scaled for CM)
		double []  corr_scale =   new double [combo_dsn[0].length];  // applied correction scale
		boolean [] was_lma =      new boolean [combo_dsn[0].length]; // previous measurement used LMA (not to match LMA/non-LMA measurements)
		Arrays.fill(disp_err, Double.NaN);
		Arrays.fill(corr_scale, 1.0);
		double min_scale = 0.75, max_scale = 1.5; // reset correction scale to 1.0 if calculated is too big/small (jumping FG/BG)
		double [][] dbg_corr_scale = null;
		if (debug_level > 0) {
			dbg_corr_scale = new double[max_refines][];
		}
		boolean [] selection = new boolean [target_disparity.length];
		for (int i = 0; i < target_disparity.length; i++) {
			selection[i] = !Double.isNaN(target_disparity[i]);
		}
		boolean [] selection_orig = selection.clone();
		
		// Refine disparity map, each time recalculating each scene "projection" pixel X, pixel Y and 
		// disparity in each scene image set that correspond to the reference scene uniform grid tile.
		// First use reduced number of pairs and no LMA, then continue with all pairs and LMA 
		boolean bg_refine= false;
		for (int nrefine = 0; nrefine < max_refines; nrefine++) {
			if (nrefine == clt_parameters.rig.mll_max_refines_pre) {
				min_disp_change = clt_parameters.rig.mll_min_disp_change_lma;				
				clt_parameters.img_dtt.setMcorr(num_sensors, save_pairs_selection); // restore
				clt_parameters.correlate_lma = save_run_lma; // restore
				selection = selection_orig.clone();
				if (debug_level > -2) {
					int num_tomeas = 0;
					for (int nt = 0; nt < target_disparity.length; nt++) if (selection[nt]) { // (!Double.isNaN(target_disparity[nt])){
						num_tomeas++;
					}
					System.out.println ("a. nrefine pass = "+nrefine+", remaining "+num_tomeas+" tiles to re-measure");
				}
			} else if (nrefine == (clt_parameters.rig.mll_max_refines_pre + clt_parameters.rig.mll_max_refines_lma)) {
				clt_parameters.img_dtt.bimax_dual_only =    true; // Camel-max tiles only
				clt_parameters.img_dtt.bimax_combine_mode = Correlation2d.CAMEL_BG;
				bg_refine= true;
				selection = selection_orig.clone(); // re-select all original tiles
			} else if (nrefine == (clt_parameters.rig.mll_max_refines_pre + clt_parameters.rig.mll_max_refines_lma + 1)) {
				clt_parameters.img_dtt.bimax_dual_only =    false; // not just dual only
			}
			
			if (nrefine == (clt_parameters.rig.mll_max_refines_pre + clt_parameters.rig.mll_max_refines_lma + clt_parameters.rig.mll_max_bg_nearest)) {
				clt_parameters.img_dtt.bimax_combine_mode = Correlation2d.CAMEL_NEAREST;
			}
			
			int mcorr_sel = Correlation2d.corrSelEncode(clt_parameters.img_dtt,num_sensors);
			double [][] disparity_map = 
					correlateInterscene( // should skip scenes w/o orientation 06/29/2022
							clt_parameters, // final CLTParameters  clt_parameters,
							scenes,         // final QuadCLT []     scenes,
							indx_ref,       // final int            indx_ref,
							target_disparity, // combo_dsn_change[combo_dsn_indx_disp],   // final double []      disparity_ref,  // disparity in the reference view tiles (Double.NaN - invalid)
							selection,      // final boolean []     selection, // may be null, if not null do not  process unselected tiles
							margin,         // final int            margin,
							nrefine,        // final int            nrefine, // just for debug title
							false, // ( nrefine == (max_refines - 1)) && clt_parameters.inp.show_final_2d, // final boolean        show_2d_corr,
					        mcorr_sel,      // final int            mcorr_sel, //  = 							
							null,           // final float [][][]   accum_2d_corr, // if [1][][] - return accumulated 2d correlations (all pairs)
							false,          // final boolean        no_map, // do not generate disparity_map (time-consuming LMA)
							debug_level-8);   // final int            debug_level)
			
			if (debug_level > 0) { //-3) {
				(new ShowDoubleFloatArrays()).showArrays(
						disparity_map,
						tilesX,
						tilesY,
						true,
						"accumulated_disparity_map-"+nrefine,
						ImageDtt.getDisparityTitles(ref_scene.getNumSensors(),ref_scene.isMonochrome()) // ImageDtt.DISPARITY_TITLES
						);
			}
			// update disparities
			double [] map_disparity =     disparity_map[ImageDtt.DISPARITY_INDEX_CM]; // 2
			double [] map_strength =      disparity_map[ImageDtt.DISPARITY_STRENGTH_INDEX]; // 10
			double [] map_disparity_lma = disparity_map[ImageDtt.DISPARITY_INDEX_POLY]; // 8
			int num_tomeas = 0; // number of tiles to measure
			Arrays.fill(selection, false);
			for (int nTile =0; nTile < combo_dsn_change[0].length; nTile++) {
				if (defined_tiles[nTile]) { // originally defined, maybe not measured last time
					if (!Double.isNaN(map_disparity[nTile])) { // re-measured
						boolean is_lma = (map_disparity_lma != null) && !Double.isNaN(map_disparity_lma[nTile]);
						if (is_lma) {
							combo_dsn_change[COMBO_DSN_INDX_CHANGE][nTile] = map_disparity_lma[nTile] * corr_nonlma;
						} else if (!Double.isNaN(map_disparity[nTile])) {
							combo_dsn_change[COMBO_DSN_INDX_CHANGE][nTile] = map_disparity[nTile] / clt_parameters.ofp.magic_scale * corr_lma;
						}
						if (!Double.isNaN(combo_dsn_change[COMBO_DSN_INDX_CHANGE][nTile])) {
							// will keep previous new_corr_scale value when no previous data or switching non-LMA/LMA. Should be initialized to 1.0 
							if (!Double.isNaN(disp_err[nTile]) && (is_lma == was_lma[nTile])) {
								// current target_disparity: combo_dsn_change[combo_dsn_indx_disp][nTile]
								double this_target_disparity = combo_dsn_change[COMBO_DSN_INDX_DISP][nTile];
								double previous_target_disparity = this_target_disparity - corr_scale[nTile] * disp_err[nTile];
								double prev_disp_err = disp_err[nTile];
								double this_disp_err = combo_dsn_change[COMBO_DSN_INDX_CHANGE][nTile];
								double new_corr_scale = (previous_target_disparity - this_target_disparity) / (this_disp_err - prev_disp_err);
								if ((new_corr_scale <= max_scale) && (new_corr_scale >= min_scale)) {
									corr_scale[nTile] = new_corr_scale;
								}
							}
							
							combo_dsn_change[COMBO_DSN_INDX_DISP][nTile] +=     combo_dsn_change[COMBO_DSN_INDX_CHANGE][nTile] * corr_scale[nTile];
							was_lma[nTile] = is_lma;
							
							combo_dsn_change[COMBO_DSN_INDX_STRENGTH][nTile]  = map_strength[nTile]; // combine CM/LMA
							if (bg_refine) {
								combo_dsn_final[COMBO_DSN_INDX_DISP_BG][nTile] =    combo_dsn_change[COMBO_DSN_INDX_DISP][nTile];
								combo_dsn_final[COMBO_DSN_INDX_STRENGTH_BG][nTile] = combo_dsn_change[COMBO_DSN_INDX_STRENGTH][nTile];
								combo_dsn_final[COMBO_DSN_INDX_LMA_BG][nTile] = combo_dsn_change[COMBO_DSN_INDX_DISP][nTile];
								if (map_disparity_lma != null) {
									combo_dsn_final[COMBO_DSN_INDX_LMA_BG][nTile] = Double.isNaN(map_disparity_lma[nTile])? Double.NaN : combo_dsn_final[COMBO_DSN_INDX_DISP_BG][nTile];
								}
								combo_dsn_final[COMBO_DSN_INDX_CHANGE_BG][nTile] = combo_dsn_change[COMBO_DSN_INDX_CHANGE][nTile];
							} else {
								combo_dsn_final[COMBO_DSN_INDX_DISP][nTile] =    combo_dsn_change[COMBO_DSN_INDX_DISP][nTile];
								combo_dsn_final[COMBO_DSN_INDX_STRENGTH][nTile] = combo_dsn_change[COMBO_DSN_INDX_STRENGTH][nTile];
								combo_dsn_final[COMBO_DSN_INDX_LMA][nTile] = combo_dsn_change[COMBO_DSN_INDX_DISP][nTile];
								if (map_disparity_lma != null) {
									combo_dsn_final[COMBO_DSN_INDX_LMA][nTile] = Double.isNaN(map_disparity_lma[nTile])? Double.NaN : combo_dsn_final[COMBO_DSN_INDX_DISP][nTile];
								}
								combo_dsn_final[COMBO_DSN_INDX_VALID][nTile] = combo_dsn[COMBO_DSN_INDX_VALID][nTile]; // not much sense
								combo_dsn_final[COMBO_DSN_INDX_CHANGE][nTile] = combo_dsn_change[COMBO_DSN_INDX_CHANGE][nTile];
							}
						}
						disp_err[nTile] = combo_dsn_change[COMBO_DSN_INDX_CHANGE][nTile];
						
						if (Math.abs(combo_dsn_change[COMBO_DSN_INDX_CHANGE][nTile]) >= min_disp_change) {
							target_disparity[nTile] = combo_dsn_change[COMBO_DSN_INDX_DISP][nTile] ;
							selection[nTile] = true;
							num_tomeas ++;
						} else {
							num_tomeas+=0;
						}
						
					} else { // originally defined, but not re-measured
						num_tomeas+=0;
					}
					
				}
			}
			if (dbg_corr_scale != null) {
				dbg_corr_scale[nrefine] = corr_scale.clone();
			}
			// Copy disparity to disparity_lma and just mask out  tiles with no LMA data (keep all if LMA did not run at all)  
			System.arraycopy(combo_dsn_change[COMBO_DSN_INDX_DISP], 0, combo_dsn_change[COMBO_DSN_INDX_LMA], 0, combo_dsn_change[COMBO_DSN_INDX_DISP].length); // lma
			if (map_disparity_lma != null) { 
				for (int i = 0; i < map_disparity_lma.length; i++) {
					if (Double.isNaN(map_disparity_lma[i])) {
						combo_dsn_change[COMBO_DSN_INDX_LMA][i] = Double.NaN;		
					}
				}
			}
			
			for (int i = 0; i < iter_indices.length; i++) {
				refine_results[last_initial_slices + (i * max_refines) + nrefine] = combo_dsn_change[iter_indices[i]].clone();
			}
			//					clt_parameters.inp.show_final_2d, // final boolean        show_2d_corr,
			if (debug_level >0) {
				(new ShowDoubleFloatArrays()).showArrays(
						combo_dsn_change,
						tilesX,
						tilesY,
						true,
						"combo_dsn-"+nrefine+"-"+ref_scene.getImageName(),
						combo_dsn_titles); //	dsrbg_titles);
			}
			if (debug_level > -2) {
				System.out.println ("b. nrefine pass = "+nrefine+", remaining "+num_tomeas+" tiles to re-measure");
			}
			if (num_tomeas == 0) {
				break;
			}
			
		} //for (int nrefine = 0; nrefine < max_refines; nrefine++) {
		// Add duplicate of FG disparity and FG+BG disparity (FG where no BG) for visual comparison
		if (clt_parameters.rig.mll_max_refines_bg > 0) { 
			for (int nTile =0; nTile < combo_dsn_change[0].length; nTile++) {
				combo_dsn_final[COMBO_DSN_INDX_DISP_FG][nTile] = combo_dsn_final[COMBO_DSN_INDX_DISP][nTile];
				if (Double.isNaN(combo_dsn_final[COMBO_DSN_INDX_DISP_BG][nTile])) {
					combo_dsn_final[COMBO_DSN_INDX_DISP_BG_ALL][nTile] = combo_dsn_final[COMBO_DSN_INDX_DISP][nTile];
				} else {
					combo_dsn_final[COMBO_DSN_INDX_DISP_BG_ALL][nTile] = combo_dsn_final[COMBO_DSN_INDX_DISP_BG][nTile];
				}
			}
		}
		
		// restore modified parameters
		clt_parameters.img_dtt.bimax_combine_mode = save_bimax_combine_mode;
		clt_parameters.img_dtt.bimax_dual_only =    save_bimax_dual_only;
		if (dbg_corr_scale != null) {
			(new ShowDoubleFloatArrays()).showArrays(
					dbg_corr_scale,
					tilesX,
					tilesY,
					true,
					"Correction scales"
					);
		}
// Do above twice: with 40 pairs, no-lma and then with all pairs+LMA
		
		if (debug_level > 1) {
			(new ShowDoubleFloatArrays()).showArrays(
					combo_dsn_change,
					tilesX,
					tilesY,
					true,
					"combo_dsn_change-"+ref_scene.getImageName(),
					combo_dsn_titles); //	dsrbg_titles);
			(new ShowDoubleFloatArrays()).showArrays(
					combo_dsn,
					tilesX,
					tilesY,
					true,
					"combo_dsn-"+ref_scene.getImageName(),
					combo_dsn_titles); //	dsrbg_titles);
			
			(new ShowDoubleFloatArrays()).showArrays(
					combo_dsn_final,
					tilesX,
					tilesY,
					true,
					"combo_dsn-final-"+ref_scene.getImageName(),
					combo_dsn_titles_full); //	dsrbg_titles);
		}
		
		
		if (debug_level > 0) {
			(new ShowDoubleFloatArrays()).showArrays(
					refine_results,
					tilesX,
					tilesY,
					true,
					"combo-"+max_refines+"-"+ref_scene.getImageName(),
					refine_titles); //	dsrbg_titles);
		}
		//noise_sigma_level
		String rslt_suffix = "-INTER-INTRA-HISTORIC";
		rslt_suffix += (clt_parameters.correlate_lma?"-LMA":"-NOLMA");

		ref_scene.saveDoubleArrayInModelDirectory(
				rslt_suffix,         // String      suffix,
				refine_titles,       // null,          // String []   labels, // or null
				refine_results,      // dbg_data,         // double [][] data,
				tilesX,              // int         width,
				tilesY);             // int         height)
		rslt_suffix = "-INTER-INTRA";
		rslt_suffix += (clt_parameters.correlate_lma?"-LMA":"-NOLMA");

		ref_scene.saveDoubleArrayInModelDirectory( // error
				rslt_suffix,           // String      suffix,
				combo_dsn_titles_full, // null,          // String []   labels, // or null
				combo_dsn_final,       // dbg_data,         // double [][] data,
				tilesX,                // int         width,
				tilesY);               // int         height)
		
		
		// save combo_dsn_change to model directory
//		if (debug_level >-100) {
//			return;
//		}
//		System.out.println("IntersceneAccumulate(), got previous scenes: "+sts.length);
		if (debug_level > 1) { // tested OK
			System.out.println("IntersceneAccumulate(): preparing image set...");
			int nscenes = scenes.length;
			//			int indx_ref = nscenes - 1; 
			double [][][] all_scenes_xyzatr = new double [scenes.length][][]; // includes reference (last)
			double [][][] all_scenes_ers_dt = new double [scenes.length][][]; // includes reference (last)
			all_scenes_xyzatr[indx_ref] = new double [][] {ZERO3,ZERO3};
			all_scenes_ers_dt[indx_ref] = new double [][] {
				ers_reference.getErsXYZ_dt(),
				ers_reference.getErsATR_dt()};
				
				for (int i = 0; i < nscenes; i++) if (i != indx_ref) {
					String ts = scenes[i].getImageName();
					all_scenes_xyzatr[i] = new double[][] {ers_reference.getSceneXYZ(ts),       ers_reference.getSceneATR(ts)}; 		
					all_scenes_ers_dt[i] = new double[][] {ers_reference.getSceneErsXYZ_dt(ts), ers_reference.getSceneErsATR_dt(ts)}; 		
				}
				compareRefSceneTiles( // null pointer
						"" ,               // String suffix,
						true, // false,             // boolean blur_reference,
						all_scenes_xyzatr, // double [][][] scene_xyzatr, // does not include reference
						all_scenes_ers_dt, // double [][][] scene_ers_dt, // does not include reference
						scenes,            // QuadCLT [] scenes,
						8);                // int iscale) // 8
		}
		// create initial disparity map for the reference scene
		
		
		// Moved to intersceneMlExport()
		
		/*
		
		boolean add_combo =         clt_parameters.rig.mll_add_combo;         //true; add 121-st slice with combined pairs correlation
		boolean save_accum =        clt_parameters.rig.mll_save_accum;        //true;  // save accumulated 0-offset correlation
		boolean randomize_offsets = clt_parameters.rig.mll_randomize_offsets; // true; 
		double  disparity_low =     clt_parameters.rig.mll_disparity_low;     // -5.0;
		double  disparity_high =    clt_parameters.rig.mll_disparity_high;    // 5.0;
		double  disparity_pwr =     clt_parameters.rig.mll_disparity_pwr;     // 2.0;
		int     disparity_steps =   clt_parameters.rig.mll_disparity_steps;   // 20;
		double  tileMetaScale =     clt_parameters.rig.mll_tileMetaScale;     // 0.001;
		int     tileMetaSlice =     clt_parameters.rig.mll_tileMetaSlice;     // -1; // all slices
		int     tileStepX =         clt_parameters.rig.mll_tileStepX;         // 16;
		int     tileStepY =         clt_parameters.rig.mll_tileStepY;         // 16;
		String  suffix =            clt_parameters.rig.mll_suffix;            // "-ML";

		double  disp_ampl = Math.max(Math.abs(disparity_low),Math.abs(disparity_high));
		//		String suffix =ref_scene.correctionsParameters.mlDirectory; // now "ML32
		double  fat_zero_single = clt_parameters.getGpuFatZero(ref_scene.isMonochrome()); // for single scene
		
		ImageDtt image_dtt;
		image_dtt = new ImageDtt(
				numSens,
				clt_parameters.transform_size,
				clt_parameters.img_dtt,
				ref_scene.isAux(),
				ref_scene.isMonochrome(),
				ref_scene.isLwir(),
				clt_parameters.getScaleStrength(ref_scene.isAux()),
				ref_scene.getGPU());
		image_dtt.getCorrelation2d(); // initiate image_dtt.correlation2d, needed if disparity_map != null  
		
		if (save_accum) {
			int mcorr_sel = Correlation2d.corrSelEncodeAll(0); // all sensors
			        float [][][] facc_2d_img = new float [1][][];
					// FIXME: 							null,           // final boolean []     selection, // may be null, if not null do not  process unselected tiles
					correlateInterscene(
							clt_parameters, // final CLTParameters  clt_parameters,
							scenes,         // final QuadCLT []     scenes,
							indx_ref,       // final int            indx_ref,
							combo_dsn_change[0],   // final double []      disparity_ref,  // disparity in the reference view tiles (Double.NaN - invalid)
							null,           // final boolean []     selection, // may be null, if not null do not  process unselected tiles
							margin,         // final int            margin,
							-1,             // final int            nrefine, // just for debug title
							false,          // final boolean        show_2d_corr,
					        mcorr_sel,      // final int            mcorr_sel, //  = 							
					        facc_2d_img,    // final float [][][]   accum_2d_corr, // if [1][][] - return accumulated 2d correlations (all pairs)
							true,           // final boolean        no_map, // do not generate disparity_map (time-consuming LMA)
							debug_level-8); // final int            debug_level)
					float [][] corr_2d_img = facc_2d_img[0];					
//					double [] 
					target_disparity = combo_dsn_change[0].clone();
					double [][] payload = {
							target_disparity,
							combo_dsn_final[0], // GT disparity
							combo_dsn_final[1], // GT confidence
							combo_dsn_final[2], // disparity_lma
							combo_dsn_final[3], // frac_valid
							combo_dsn_final[4]  // last_diff
					};
					for (int i = 0; i < payload.length; i++) {
						add_tile_meta(
								corr_2d_img,   // final float [][] fimg,
								tilesX,        // final int tilesX,
								tilesY,        // final int tilesY,
								tileStepX,     // final int stepX,
								tileStepY,     // final int stepY,
								tileMetaScale, //final double payload_scale,
								payload[i],    // final double [] payload,
								tileMetaSlice, // final int slice,
								i,             // final int offsX,
								tileStepY-1);  // final int offsY)
					}
					String [] titles = new String [corr_2d_img.length]; // dcorr_tiles[0].length];
					int ind_length = image_dtt.getCorrelation2d().getCorrTitles().length;

					System.arraycopy(image_dtt.getCorrelation2d().getCorrTitles(), 0, titles, 0, ind_length);
					for (int i = ind_length; i < titles.length; i++) {
						titles[i] = "combo-"+(i - ind_length);
					}
					ImageStack ml_stack = (new ShowDoubleFloatArrays()).makeStack(
							corr_2d_img,                         // float[][] pixels,
							tilesX*(2*image_dtt.transform_size), // int width,
							tilesY*(2*image_dtt.transform_size), // int height,
							titles, // String [] titles,
							false); // boolean noNaN)
					String x3d_path = ref_scene.getX3dDirectory();
					String title = ref_scene.getImageName() + suffix+
							(ref_scene.isAux()?"-AUX":"-MAIN")+"-ACCUM";
					String mldir=ref_scene.correctionsParameters.mlDirectory;
					String aMldir=x3d_path + Prefs.getFileSeparator() + mldir;
					String file_path = aMldir + Prefs.getFileSeparator() + title + ".tiff";
					File dir = (new File(file_path)).getParentFile();
					if (!dir.exists()){
						dir.mkdirs();
					}
					ImagePlus imp_ml = new ImagePlus(title, ml_stack);
					imp_ml.setProperty("VERSION",             "2.0");
//					imp_ml.setProperty("tileWidth",    ""+ml_width);
					imp_ml.setProperty("numScenes",           ""+num_scenes);
					imp_ml.setProperty("indexReference",      ""+indx_ref);
					imp_ml.setProperty("fatZero",             ""+fat_zero_single);
					imp_ml.setProperty("dispOffset",          ""+0);
					imp_ml.setProperty("tileMetaScale",       ""+tileMetaScale);
					imp_ml.setProperty("tileMetaSlice",       ""+tileMetaSlice);
					imp_ml.setProperty("tileStepX",           ""+tileStepX);
					imp_ml.setProperty("tileStepY",           ""+tileStepY);
					imp_ml.setProperty("metaTargetDisparity", ""+0);
					imp_ml.setProperty("metaGTDisparity",     ""+1);
					imp_ml.setProperty("metaGTConfidence",    ""+2);
					imp_ml.setProperty("metaGTDisparityLMA",  ""+3);
					imp_ml.setProperty("metaFracValid",       ""+4);
					imp_ml.setProperty("metaLastDiff",        ""+5);
					(new JP46_Reader_camera(false)).encodeProperiesToInfo(imp_ml);			
					FileSaver fs=new FileSaver(imp_ml);
					fs.saveAsTiff(file_path);
					System.out.println("intersceneExport(): saved "+file_path);
		}
		if ( clt_parameters.ofp.pattern_mode) {
			return combo_dsn_final;
		}
		double [][] all_offsets = new double [disparity_steps][];
		String [] soffset_centers = new String [disparity_steps];
		for (int nstep = 0; nstep < disparity_steps; nstep++) {
			double [] disparity_offsets_rel = { // below, center, above
					(disparity_low + (disparity_high - disparity_low) * (nstep - 1) / (disparity_steps-1))/disp_ampl,
					(disparity_low + (disparity_high - disparity_low) * (nstep + 0) / (disparity_steps-1))/disp_ampl,
					(disparity_low + (disparity_high - disparity_low) * (nstep + 1) / (disparity_steps-1))/disp_ampl
			};
			double [] disparity_offset3 = new double [disparity_offsets_rel.length];
			for (int i = 0; i < disparity_offsets_rel.length; i++) {
				double dsgn = (disparity_offsets_rel[i] > 0.0) ? 1.0 : ((disparity_offsets_rel[i] < 0.0)? -1.0 : 0.0);
				disparity_offset3[i] = Math.pow(Math.abs(disparity_offsets_rel[i]), disparity_pwr) * dsgn * disp_ampl ;
			}
			if (disparity_offset3[0] > disparity_offset3[2]) { // can that happen?
				double d = disparity_offset3[0];
				disparity_offset3[0] =disparity_offset3[2];
				disparity_offset3[2] = d; 
			}
			double disparity_offset = disparity_offset3[1];
			if (debug_level > -2) {
				System.out.println("dispatity offset #"+(nstep + 1)+" (of "+disparity_steps+") = "+disparity_offset);
			}
			double [] disparity_offsets;
			if (randomize_offsets) {
				disparity_offsets = getLappedRandom(
						disparity_offset3[0], // double min_exclusive,
						disparity_offset3[2], // double max_exclusive,
						tiles); // int    nSamples);
			} else {
				disparity_offsets = new double[tiles];
				Arrays.fill(disparity_offsets, disparity_offset3[0]);
			}
			all_offsets[nstep] = disparity_offsets;
			
			float [][] corr_2d_img = generateOffset2DCorrelations(
					clt_parameters, // final CLTParameters  clt_parameters,
					scenes[indx_ref], // final QuadCLT        ref_scene,
					combo_dsn_change[0], // final double []      disparity_ref_in,  // disparity in the reference view tiles (Double.NaN - invalid)
					disparity_offsets, // disparity_offset, // final double         disparity_offset,
					margin, // final int            margin,
					add_combo, // final boolean        add_combo,
					debug_level-9); // final int            debug_level);
//			double [] 
			target_disparity = combo_dsn_change[0].clone();
			for (int i = 0; i < target_disparity.length; i++) {
				target_disparity[i]+= disparity_offsets[i];
			}
			double [][] payload = {
					target_disparity,
					combo_dsn_final[0], // GT disparity
					combo_dsn_final[1], // GT confidence - wrong
					combo_dsn_final[2], // disparity_lma - Wrong !
					combo_dsn_final[3], // frac_valid
					combo_dsn_final[4]  // last_diff
			};
			for (int i = 0; i < payload.length; i++) {
				add_tile_meta(
						corr_2d_img,   // final float [][] fimg,
						tilesX,        // final int tilesX,
						tilesY,        // final int tilesY,
						tileStepX,     // final int stepX,
						tileStepY,     // final int stepY,
						tileMetaScale, //final double payload_scale,
						payload[i],    // final double [] payload,
						tileMetaSlice, // final int slice,
						i,             // final int offsX,
						tileStepY-1);  // final int offsY)
			}
			
			String [] titles = new String [corr_2d_img.length]; // dcorr_tiles[0].length];
			int ind_length = image_dtt.getCorrelation2d().getCorrTitles().length;

			System.arraycopy(image_dtt.getCorrelation2d().getCorrTitles(), 0, titles, 0, ind_length);
			for (int i = ind_length; i < titles.length; i++) {
				titles[i] = "combo-"+(i - ind_length);
			}
			ImageStack ml_stack = (new ShowDoubleFloatArrays()).makeStack(
					corr_2d_img,                         // float[][] pixels,
					tilesX*(2*image_dtt.transform_size), // int width,
					tilesY*(2*image_dtt.transform_size), // int height,
					titles, // String [] titles,
					false); // boolean noNaN)
			String x3d_path = ref_scene.getX3dDirectory();
			String sdisparity_offset = String.format("%8.3f", disparity_offset).trim();
			soffset_centers[nstep] = sdisparity_offset;

			String title = ref_scene.getImageName() + suffix+
					(ref_scene.isAux()?"-AUX":"-MAIN");
			if (randomize_offsets) {
				title+="-RND";
			}			
			title+="-DOFFS"+ sdisparity_offset;
			String mldir=ref_scene.correctionsParameters.mlDirectory;
			String aMldir=x3d_path + Prefs.getFileSeparator() + mldir;
			String file_path = aMldir + Prefs.getFileSeparator() + title + ".tiff";
			File dir = (new File(file_path)).getParentFile();
			if (!dir.exists()){
				dir.mkdirs();
			}
			ImagePlus imp_ml = new ImagePlus(title, ml_stack);
			imp_ml.setProperty("VERSION",             "2.0");
			imp_ml.setProperty("numScenes",           ""+num_scenes);
			imp_ml.setProperty("indexReference",      ""+indx_ref);
			imp_ml.setProperty("fatZero",             ""+fat_zero_single);
			imp_ml.setProperty("dispOffset",          ""+disparity_offset3[1]);
			imp_ml.setProperty("randomize_offsets",   ""+randomize_offsets);
			if (randomize_offsets) {
				imp_ml.setProperty("dispOffsetLow",   ""+disparity_offset3[0]);
				imp_ml.setProperty("dispOffsetHigh",  ""+disparity_offset3[2]);
			}
			imp_ml.setProperty("disparity_low",       ""+disparity_low);
			imp_ml.setProperty("disparity_high",      ""+disparity_high);
			imp_ml.setProperty("disparity_pwr",       ""+disparity_pwr);
			imp_ml.setProperty("disparity_steps",     ""+disparity_steps);
			imp_ml.setProperty("tileMetaScale",       ""+tileMetaScale);
			imp_ml.setProperty("tileMetaSlice",       ""+tileMetaSlice);
			imp_ml.setProperty("tileStepX",           ""+tileStepX);
			imp_ml.setProperty("tileStepY",           ""+tileStepY);
			imp_ml.setProperty("metaTargetDisparity", ""+0);
			imp_ml.setProperty("metaGTDisparity",     ""+1);
			imp_ml.setProperty("metaGTConfidence",    ""+2);
			imp_ml.setProperty("metaFracValid",       ""+3);
			imp_ml.setProperty("metaLastDiff",        ""+4);
			(new JP46_Reader_camera(false)).encodeProperiesToInfo(imp_ml);			
			FileSaver fs=new FileSaver(imp_ml);
			fs.saveAsTiff(file_path);
			System.out.println("intersceneExport(): saved "+file_path);
		}
		if (disparity_steps > 0) {
			String offsets_suffix = "-DISP_OFFSETS";
			if (randomize_offsets) {
				offsets_suffix+="-RND";
			}			

			ref_scene.saveDoubleArrayInModelDirectory(
					offsets_suffix,      // String      suffix,
					soffset_centers,     // null,          // String []   labels, // or null
					all_offsets,         // dbg_data,         // double [][] data,
					tilesX,              // int         width,
					tilesY);             // int         height)
		}
		*/
		return combo_dsn_final;
	}
	
	public void intersceneMlExport(
			CLTParameters        clt_parameters,
			ErsCorrection        ers_reference,
			QuadCLT []           scenes,
			ColorProcParameters  colorProcParameters,
			double [][]          combo_dsn_final,
			int                  debug_level
			) {
		
		final int margin = 8;
		boolean add_combo =         clt_parameters.rig.mll_add_combo;         //true; add 121-st slice with combined pairs correlation
		boolean save_accum =        clt_parameters.rig.mll_save_accum;        //true;  // save accumulated 0-offset correlation
		boolean randomize_offsets = clt_parameters.rig.mll_randomize_offsets; // true; 
		double  disparity_low =     clt_parameters.rig.mll_disparity_low;     // -5.0;
		double  disparity_high =    clt_parameters.rig.mll_disparity_high;    // 5.0;
		double  disparity_pwr =     clt_parameters.rig.mll_disparity_pwr;     // 2.0;
		int     disparity_steps =   clt_parameters.rig.mll_disparity_steps;   // 20;
		double  tileMetaScale =     clt_parameters.rig.mll_tileMetaScale;     // 0.001;
		int     tileMetaSlice =     clt_parameters.rig.mll_tileMetaSlice;     // -1; // all slices
		int     tileStepX =         clt_parameters.rig.mll_tileStepX;         // 16;
		int     tileStepY =         clt_parameters.rig.mll_tileStepY;         // 16;
		String  suffix =            clt_parameters.rig.mll_suffix;            // "-ML";
		boolean show_input =        clt_parameters.img_dtt.lmamask_dbg; // true; // debug_level = -1;
		

		double  disp_ampl = Math.max(Math.abs(disparity_low),Math.abs(disparity_high));
		//		String suffix =ref_scene.correctionsParameters.mlDirectory; // now "ML32
		int      indx_ref = scenes.length - 1; // Always added to the end even if out-of order
		QuadCLT  ref_scene = scenes[indx_ref]; // ordered by increasing timestamps
//		int num_scenes =  scenes.length; // FIXME: skip unmatched
		final int num_scenes =  indx_ref - ref_scene.getEarliestScene(scenes) + 1;
		final int tilesX = ref_scene.getTileProcessor().getTilesX();
		final int tilesY = ref_scene.getTileProcessor().getTilesY();
		final int tiles = tilesX * tilesY;
		
		double  fat_zero_single = clt_parameters.getGpuFatZero(ref_scene.isMonochrome()); // for single scene
		ImageDtt image_dtt;
		image_dtt = new ImageDtt(
				numSens,
				clt_parameters.transform_size,
				clt_parameters.img_dtt,
				ref_scene.isAux(),
				ref_scene.isMonochrome(),
				ref_scene.isLwir(),
				clt_parameters.getScaleStrength(ref_scene.isAux()),
				ref_scene.getGPU());
		image_dtt.getCorrelation2d(); // initiate image_dtt.correlation2d, needed if disparity_map != null  
		
		if (show_input) {
			//COMBO_DSN_TITLES
			(new ShowDoubleFloatArrays()).showArrays(
					combo_dsn_final,
					tilesX,
					tilesY,
					true,
					"combo_dsn_final-"+ref_scene.getImageName(),
					COMBO_DSN_TITLES); //	dsrbg_titles);

		}
		
		if (save_accum) {
			int mcorr_sel = Correlation2d.corrSelEncodeAll(0); // all sensors
			float [][][] facc_2d_img = new float [1][][];
			// FIXME?: should work with non-matched
			boolean        no_map = true;
			boolean        show_2d_corr = false;
			if (clt_parameters.img_dtt.lmamask_dbg) {
				no_map = false;
				show_2d_corr = true;
			}
			
			double [][] disparity_map = correlateInterscene(
					clt_parameters, // final CLTParameters  clt_parameters,
					scenes,         // final QuadCLT []     scenes,
					indx_ref,       // final int            indx_ref,
					combo_dsn_final[COMBO_DSN_INDX_DISP],   // final double []      disparity_ref,  // disparity in the reference view tiles (Double.NaN - invalid)
					null,           // final boolean []     selection, // may be null, if not null do not  process unselected tiles
					margin,         // final int            margin,
					-1,             // final int            nrefine, // just for debug title
					show_2d_corr,   // final boolean        show_2d_corr,
					mcorr_sel,      // final int            mcorr_sel, //  = 							
					facc_2d_img,    // final float [][][]   accum_2d_corr, // if [1][][] - return accumulated 2d correlations (all pairs)
					no_map,         // final boolean        no_map, // do not generate disparity_map (time-consuming LMA)
					clt_parameters.img_dtt.lmamask_dbg? 1:(debug_level-8)); // final int            debug_level)
			if (disparity_map != null) {
				(new ShowDoubleFloatArrays()).showArrays(
						disparity_map,
						tilesX,
						tilesY,
						true,
						"disparity_map_debug-"+ref_scene.getImageName(),
						image_dtt.getDisparityTitles()); //	dsrbg_titles);

			}
			float [][] corr_2d_img = facc_2d_img[0];					
			double [] target_disparity = combo_dsn_final[COMBO_DSN_INDX_DISP].clone();
			double [][] payload = {
					target_disparity,
					combo_dsn_final[COMBO_DSN_INDX_DISP],       // GT disparity
					combo_dsn_final[COMBO_DSN_INDX_STRENGTH],   // GT confidence
					combo_dsn_final[COMBO_DSN_INDX_LMA],        // disparity_lma
					combo_dsn_final[COMBO_DSN_INDX_VALID],      // frac_valid
					combo_dsn_final[COMBO_DSN_INDX_CHANGE],     // last_diff
					combo_dsn_final[COMBO_DSN_INDX_DISP_BG],    // cumulative BG disparity (from CM or POLY)
					combo_dsn_final[COMBO_DSN_INDX_STRENGTH_BG],// background strength
					combo_dsn_final[COMBO_DSN_INDX_LMA_BG],     // masked copy from BG disparity 
					combo_dsn_final[COMBO_DSN_INDX_CHANGE_BG],  // increment, BG
					combo_dsn_final[COMBO_DSN_INDX_DISP_FG],    //cumulative disparity (from CM or POLY), FG == COMBO_DSN_INDX_DISP
					combo_dsn_final[COMBO_DSN_INDX_DISP_BG_ALL] // cumulative BG disparity (Use FG where no BG is available)
			};
			for (int i = 0; i < payload.length; i++) {
				add_tile_meta(
						corr_2d_img,   // final float [][] fimg,
						tilesX,        // final int tilesX,
						tilesY,        // final int tilesY,
						tileStepX,     // final int stepX,
						tileStepY,     // final int stepY,
						tileMetaScale, //final double payload_scale,
						payload[i],    // final double [] payload,
						tileMetaSlice, // final int slice,
						i,             // final int offsX,
						tileStepY-1);  // final int offsY)
			}
			String [] titles = new String [corr_2d_img.length]; // dcorr_tiles[0].length];
			int ind_length = image_dtt.getCorrelation2d().getCorrTitles().length;

			System.arraycopy(image_dtt.getCorrelation2d().getCorrTitles(), 0, titles, 0, ind_length);
			for (int i = ind_length; i < titles.length; i++) {
				titles[i] = "combo-"+(i - ind_length);
			}
			ImageStack ml_stack = (new ShowDoubleFloatArrays()).makeStack(
					corr_2d_img,                         // float[][] pixels,
					tilesX*(2*image_dtt.transform_size), // int width,
					tilesY*(2*image_dtt.transform_size), // int height,
					titles, // String [] titles,
					false); // boolean noNaN)
			String x3d_path = ref_scene.getX3dDirectory();
			String title = ref_scene.getImageName() + suffix+
					(ref_scene.isAux()?"-AUX":"-MAIN")+"-ACCUM";
			String mldir=ref_scene.correctionsParameters.mlDirectory;
			String aMldir=x3d_path + Prefs.getFileSeparator() + mldir;
			String file_path = aMldir + Prefs.getFileSeparator() + title + ".tiff";
			File dir = (new File(file_path)).getParentFile();
			if (!dir.exists()){
				dir.mkdirs();
			}
			ImagePlus imp_ml = new ImagePlus(title, ml_stack);
			imp_ml.setProperty("VERSION",             "2.0");
			//					imp_ml.setProperty("tileWidth",    ""+ml_width);
			imp_ml.setProperty("numScenes",           ""+num_scenes);
			imp_ml.setProperty("indexReference",      ""+indx_ref);
			imp_ml.setProperty("fatZero",             ""+fat_zero_single);
			imp_ml.setProperty("dispOffset",          ""+0);
			imp_ml.setProperty("tileMetaScale",       ""+tileMetaScale);
			imp_ml.setProperty("tileMetaSlice",       ""+tileMetaSlice);
			imp_ml.setProperty("tileStepX",           ""+tileStepX);
			imp_ml.setProperty("tileStepY",           ""+tileStepY);
			imp_ml.setProperty("metaTargetDisparity", ""+0);
			imp_ml.setProperty("metaGTDisparity",     ""+1);
			imp_ml.setProperty("metaGTConfidence",    ""+2);
			imp_ml.setProperty("metaGTDisparityLMA",  ""+3);
			imp_ml.setProperty("metaFracValid",       ""+4);
			imp_ml.setProperty("metaLastDiff",        ""+5);
			imp_ml.setProperty("metaBGDisparity",     ""+6);
			imp_ml.setProperty("metaBGConfidence",    ""+7);
			imp_ml.setProperty("metaBGDisparityLMA",  ""+8);
			imp_ml.setProperty("metaBGLastDiff",      ""+9);
			imp_ml.setProperty("metaFGDisparity",     ""+10); //=="metaGTDisparity"
			imp_ml.setProperty("metaBGDisparityAll",  ""+11);
			(new JP46_Reader_camera(false)).encodeProperiesToInfo(imp_ml);			
			FileSaver fs=new FileSaver(imp_ml);
			fs.saveAsTiff(file_path);
			System.out.println("intersceneExport(): saved "+file_path);
		}
		if ( clt_parameters.ofp.pattern_mode) {
			return; //  combo_dsn_final;
		}
		double [][] all_offsets = new double [disparity_steps][];
		String [] soffset_centers = new String [disparity_steps];
		for (int nstep = 0; nstep < disparity_steps; nstep++) {
			double [] disparity_offsets_rel = { // below, center, above
					(disparity_low + (disparity_high - disparity_low) * (nstep - 1) / (disparity_steps-1))/disp_ampl,
					(disparity_low + (disparity_high - disparity_low) * (nstep + 0) / (disparity_steps-1))/disp_ampl,
					(disparity_low + (disparity_high - disparity_low) * (nstep + 1) / (disparity_steps-1))/disp_ampl
			};
			double [] disparity_offset3 = new double [disparity_offsets_rel.length];
			for (int i = 0; i < disparity_offsets_rel.length; i++) {
				double dsgn = (disparity_offsets_rel[i] > 0.0) ? 1.0 : ((disparity_offsets_rel[i] < 0.0)? -1.0 : 0.0);
				disparity_offset3[i] = Math.pow(Math.abs(disparity_offsets_rel[i]), disparity_pwr) * dsgn * disp_ampl ;
			}
			if (disparity_offset3[0] > disparity_offset3[2]) { // can that happen?
				double d = disparity_offset3[0];
				disparity_offset3[0] =disparity_offset3[2];
				disparity_offset3[2] = d; 
			}
			double disparity_offset = disparity_offset3[1];
			if (debug_level > -2) {
				System.out.println("dispatity offset #"+(nstep + 1)+" (of "+disparity_steps+") = "+disparity_offset);
			}
			double [] disparity_offsets;
			if (randomize_offsets) {
				disparity_offsets = getLappedRandom(
						disparity_offset3[0], // double min_exclusive,
						disparity_offset3[2], // double max_exclusive,
						tiles); // int    nSamples);
			} else {
				disparity_offsets = new double[tiles];
				Arrays.fill(disparity_offsets, disparity_offset3[0]);
			}
			all_offsets[nstep] = disparity_offsets;
			
			float [][] corr_2d_img = generateOffset2DCorrelations(
					clt_parameters, // final CLTParameters  clt_parameters,
					scenes[indx_ref], // final QuadCLT        ref_scene,
					combo_dsn_final[COMBO_DSN_INDX_DISP], // final double []      disparity_ref_in,  // disparity in the reference view tiles (Double.NaN - invalid)
					disparity_offsets, // disparity_offset, // final double         disparity_offset,
					margin, // final int            margin,
					add_combo, // final boolean        add_combo,
					debug_level-9); // final int            debug_level);
			double [] target_disparity = combo_dsn_final[COMBO_DSN_INDX_DISP].clone();
			for (int i = 0; i < target_disparity.length; i++) {
				target_disparity[i]+= disparity_offsets[i];
			}
			double [][] payload = {
					target_disparity,
					combo_dsn_final[COMBO_DSN_INDX_DISP],       // GT disparity
					combo_dsn_final[COMBO_DSN_INDX_STRENGTH],   // GT confidence
					combo_dsn_final[COMBO_DSN_INDX_LMA],        // disparity_lma
					combo_dsn_final[COMBO_DSN_INDX_VALID],      // frac_valid
					combo_dsn_final[COMBO_DSN_INDX_CHANGE],     // last_diff
					combo_dsn_final[COMBO_DSN_INDX_DISP_BG],    // cumulative BG disparity (from CM or POLY)
					combo_dsn_final[COMBO_DSN_INDX_STRENGTH_BG],// background strength
					combo_dsn_final[COMBO_DSN_INDX_LMA_BG],     // masked copy from BG disparity 
					combo_dsn_final[COMBO_DSN_INDX_CHANGE_BG],  // increment, BG
					combo_dsn_final[COMBO_DSN_INDX_DISP_FG],    //cumulative disparity (from CM or POLY), FG == COMBO_DSN_INDX_DISP
					combo_dsn_final[COMBO_DSN_INDX_DISP_BG_ALL] // cumulative BG disparity (Use FG where no BG is available)
			};
			for (int i = 0; i < payload.length; i++) {
				add_tile_meta(
						corr_2d_img,   // final float [][] fimg,
						tilesX,        // final int tilesX,
						tilesY,        // final int tilesY,
						tileStepX,     // final int stepX,
						tileStepY,     // final int stepY,
						tileMetaScale, // final double payload_scale,
						payload[i],    // final double [] payload,
						tileMetaSlice, // final int slice,
						i,             // final int offsX,
						tileStepY-1);  // final int offsY)
			}
			
			String [] titles = new String [corr_2d_img.length]; // dcorr_tiles[0].length];
			int ind_length = image_dtt.getCorrelation2d().getCorrTitles().length;

			System.arraycopy(image_dtt.getCorrelation2d().getCorrTitles(), 0, titles, 0, ind_length);
			for (int i = ind_length; i < titles.length; i++) {
				titles[i] = "combo-"+(i - ind_length);
			}
			ImageStack ml_stack = (new ShowDoubleFloatArrays()).makeStack(
					corr_2d_img,                         // float[][] pixels,
					tilesX*(2*image_dtt.transform_size), // int width,
					tilesY*(2*image_dtt.transform_size), // int height,
					titles, // String [] titles,
					false); // boolean noNaN)
			String x3d_path = ref_scene.getX3dDirectory();
			String sdisparity_offset = String.format("%8.3f", disparity_offset).trim();
			soffset_centers[nstep] = sdisparity_offset;

			String title = ref_scene.getImageName() + suffix+
					(ref_scene.isAux()?"-AUX":"-MAIN");
			if (randomize_offsets) {
				title+="-RND";
			}			
			title+="-DOFFS"+ sdisparity_offset;
			String mldir=ref_scene.correctionsParameters.mlDirectory;
			String aMldir=x3d_path + Prefs.getFileSeparator() + mldir;
			String file_path = aMldir + Prefs.getFileSeparator() + title + ".tiff";
			File dir = (new File(file_path)).getParentFile();
			if (!dir.exists()){
				dir.mkdirs();
			}
			ImagePlus imp_ml = new ImagePlus(title, ml_stack);
			imp_ml.setProperty("VERSION",             "2.0");
			imp_ml.setProperty("numScenes",           ""+num_scenes);
			imp_ml.setProperty("indexReference",      ""+indx_ref);
			imp_ml.setProperty("fatZero",             ""+fat_zero_single);
			imp_ml.setProperty("dispOffset",          ""+disparity_offset3[1]);
			imp_ml.setProperty("randomize_offsets",   ""+randomize_offsets);
			if (randomize_offsets) {
				imp_ml.setProperty("dispOffsetLow",   ""+disparity_offset3[0]);
				imp_ml.setProperty("dispOffsetHigh",  ""+disparity_offset3[2]);
			}
			imp_ml.setProperty("disparity_low",       ""+disparity_low);
			imp_ml.setProperty("disparity_high",      ""+disparity_high);
			imp_ml.setProperty("disparity_pwr",       ""+disparity_pwr);
			imp_ml.setProperty("disparity_steps",     ""+disparity_steps);
			imp_ml.setProperty("tileMetaScale",       ""+tileMetaScale);
			imp_ml.setProperty("tileMetaSlice",       ""+tileMetaSlice);
			imp_ml.setProperty("tileStepX",           ""+tileStepX);
			imp_ml.setProperty("tileStepY",           ""+tileStepY);
			imp_ml.setProperty("metaTargetDisparity", ""+0);
			imp_ml.setProperty("metaGTDisparity",     ""+1);
			imp_ml.setProperty("metaGTConfidence",    ""+2);
			imp_ml.setProperty("metaGTDisparityLMA",  ""+3);
			imp_ml.setProperty("metaFracValid",       ""+4);
			imp_ml.setProperty("metaLastDiff",        ""+5);
			imp_ml.setProperty("metaBGDisparity",     ""+6);
			imp_ml.setProperty("metaBGConfidence",    ""+7);
			imp_ml.setProperty("metaBGDisparityLMA",  ""+8);
			imp_ml.setProperty("metaBGLastDiff",      ""+9);
			imp_ml.setProperty("metaFGDisparity",     ""+10); //=="metaGTDisparity"
			imp_ml.setProperty("metaBGDisparityAll",  ""+11);
			(new JP46_Reader_camera(false)).encodeProperiesToInfo(imp_ml);			
			FileSaver fs=new FileSaver(imp_ml);
			fs.saveAsTiff(file_path);
			System.out.println("intersceneExport(): saved "+file_path);
		}
		if (disparity_steps > 0) {
			String offsets_suffix = "-DISP_OFFSETS";
			if (randomize_offsets) {
				offsets_suffix+="-RND";
			}			

			ref_scene.saveDoubleArrayInModelDirectory(
					offsets_suffix,      // String      suffix,
					soffset_centers,     // null,          // String []   labels, // or null
					all_offsets,         // dbg_data,         // double [][] data,
					tilesX,              // int         width,
					tilesY);             // int         height)
		}
		
	}
	
	/**
	 * Generate a pseudo-random values in the range (min_exclusive, max_exclusive) - both exclusive
	 * with a histogram being a shifted cosine, so an overlap of such distributions shifted by half-range
	 * will result in even distribution over the full range
	 * @param min_exclusive minimal value (exclusive)
	 * @param max_exclusive maximal value (exclusive)
	 * @param nSamples number of samples to generate
	 * @return nSamples-long array of values
	 */
	public double [] getLappedRandom(
			double min_exclusive,
			double max_exclusive,
			int    nSamples) {
		return getLappedRandom(
				min_exclusive,
				max_exclusive,
				nSamples,
				1E-8,
				100);
	}

	public double [] getLappedRandom(
			double min_exclusive,
			double max_exclusive,
			int    nSamples,
			double e,
			int    num_iter)
	{
		final double [] rslt = new double [nSamples];
		final double offset =    0.5 * (max_exclusive + min_exclusive);
		final double amplitude = 0.5 * (max_exclusive - min_exclusive);
		final Thread[] threads = ImageDtt.newThreadArray(threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				public void run() {
					for (int nTile = ai.getAndIncrement(); nTile < nSamples; nTile = ai.getAndIncrement()) {
						double x = ThreadLocalRandom.current().nextDouble(); // only positives, will mirror later
						// use Newton method to invert x=y+sin(pi*y)/pi
						double y = 1 - Math.sqrt(1 - x); // initial approximation
						for (int n = 0; n < num_iter; n++) {
							double xi = y + Math.sin(Math.PI * y)/Math.PI;
							if (Math.abs (x - xi) < e) {
								break;
							}
							y += (x - xi) / (1.0 + Math.cos(Math.PI * y));
						}
						if (ThreadLocalRandom.current().nextBoolean()) { // negative values, now -1.0<y<1.0 
							y = -y; 
						}
						rslt[nTile] = offset + amplitude * y;
					}
				}
			};
		}		      
		ImageDtt.startAndJoin(threads);
		return rslt;
	}
	
	
	/**
	 * Add per-tile metada into the gaps between the image tiles
	 * @param fimg    image to add metadata [slices][(tilesY*stepY)*(tilesX*stepX)
	 * @param tilesX  number of image tiles in X direction
	 * @param tilesY  number of image tiles in Y direction
	 * @param stepX   width of the image tile in pixels
	 * @param stepY   height of the image tile in pixels 
	 * @param payload_scale Multiply payload data to reduce its contrast (0.001) 
	 * @param payload metadata [tilesY*tilesX] to be embedded, one double per tile
	 * @param slice   image slice number to embedd (first index in fimg). Negative - embed into all slices
	 * @param offsX   pixel horizontal (righth) offset to place metadata inside each tile. Typical 0,1,2,3 ...
	 * @param offsY   pixel vertical (down) offset to place metadata inside each tile Typical 15
	 */
	public void add_tile_meta(
			final float [][] fimg,
			final int tilesX,
			final int tilesY,
			final int stepX,
			final int stepY,
			final double payload_scale,
			final double [] payload,
			final int slice,
			final int offsX,
			final int offsY) {
		final int width = tilesX * stepX;
		final int offs = offsX + width*offsY;
		final int tiles = tilesX * tilesY; 
		final Thread[] threads = ImageDtt.newThreadArray(threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				public void run() {
					for (int nTile = ai.getAndIncrement(); nTile < tiles; nTile = ai.getAndIncrement()) {
						int tileY = nTile / tilesX;  
						int tileX = nTile % tilesX;
						float d = (float) (payload_scale* payload[nTile]);
						if (slice >=0) {
							fimg[slice][offs + tileY * width * stepY + tileX * stepX]  = d;
						} else {
							for (int s = 0; s < fimg.length; s ++) {
								fimg[s][offs + tileY * width * stepY + tileX * stepX]  = d;
							}
						}
					}
				}
			};
		}		      
		ImageDtt.startAndJoin(threads);
	}
	
	
	
	
	
	
	public void intersceneNoise(
			CLTParameters        clt_parameters,
			boolean              ref_only, // process only reference frame (false - inter-scene)
			ColorProcParameters  colorProcParameters,
			QuadCLT              ref_scene, // ordered by increasing timestamps
//			double []
			NoiseParameters		 noise_sigma_level,
			int                  noise_variant, // <0 - no-variants, compatible with old code			
			int                  debug_level
			)
	{
		System.out.println("IntersceneNoise(), scene timestamp="+ref_scene.getImageName());
		ErsCorrection ers_reference = ref_scene.getErsCorrection();
		String [] sts = ref_only ? (new String [0]) : ers_reference.getScenes();
		// get list of all other scenes
		int num_scenes = sts.length + 1;
		int indx_ref = num_scenes - 1; 
		QuadCLT [] scenes = new QuadCLT [num_scenes];
		scenes[indx_ref] = ref_scene;
		
		for (int i = 0; i < sts.length; i++) {
			scenes[i] = ref_scene.spawnQuadCLTWithNoise( // spawnQuadCLT(
					sts[i],
					clt_parameters,
					colorProcParameters, //
					noise_sigma_level,   // double []            noise_sigma_level,
					noise_variant,       // int                  noise_variant, // <0 - no-variants, compatible with old code
					ref_scene,          // QuadCLTCPU           ref_scene, // may be null if scale_fpn <= 0
					threadsMax,
					-1); // debug_level);
			scenes[i].setDSRBG(
					clt_parameters, // CLTParameters  clt_parameters,
					threadsMax,     // int            threadsMax,  // maximal number of threads to launch
					updateStatus,   // boolean        updateStatus,
					-1); // debug_level);    // int            debugLevel)
		}
//		String [] combo_dsn_titles = {"disp", "strength", "num_valid","change"};
		String [] combo_dsn_titles = {"disp", "strength","disp_lma","num_valid","change"};
		int combo_dsn_indx_disp =     0; // cumulative disparity (from CM or POLY)
		int combo_dsn_indx_strength = 1;
		int combo_dsn_indx_lma =      2; // masked copy from 0 - cumulative disparity
		int combo_dsn_indx_valid =    3; // initial only
		int combo_dsn_indx_change =   4; // increment
		boolean read_nonoise_lma =clt_parameters.correlate_lma || true; // read always
//		final String [] iter_titles = {"disp", "diff", "strength","disp_lma"};
		final int [] iter_indices = {
				combo_dsn_indx_disp,
				combo_dsn_indx_strength,
				combo_dsn_indx_lma,
				combo_dsn_indx_change}; // which to save for each iteration: {"disp", "strength","disp_lma","change"};
		final int [] initial_indices = {
				combo_dsn_indx_disp,
				combo_dsn_indx_strength,
				combo_dsn_indx_valid}; // initial: "disp", "strength","num_valid"
		double [][] combo_dsn =  null;
		if (noise_sigma_level == null) {
			double[][] combo_dsn0 = prepareInitialComboDS( // 3
					clt_parameters,   // final CLTParameters       clt_parameters,
					scenes,           // final QuadCLT []          scenes,
					indx_ref,         // final int                 indx_ref,
					debug_level-2);     // final int                 debug_level);
			combo_dsn = new double[combo_dsn_titles.length - 1][];
			for (int i = 0; i < combo_dsn0.length; i++) {
				combo_dsn[initial_indices[i]] = combo_dsn0[i]; // "disp", "strength", <null>, "num_valid"
			}
			
		} else {
			combo_dsn = ref_scene.readDoubleArrayFromModelDirectory( //"disp", "strength","disp_lma","num_valid"
					"-results-nonoise" + (read_nonoise_lma?"-lma":"-nolma"), // String      suffix,
					combo_dsn_titles.length - 1, // 4
					null); // int []      wh);

		}
		
		
		
//		final double [][] combo_dsn_change = new double [combo_dsn.length+1][];
		final int margin = 8;
		final int tilesX = ref_scene.getTileProcessor().getTilesX();
		final int tilesY = ref_scene.getTileProcessor().getTilesY();
		if (debug_level > 0) {
			int        extra = 10; // pixels around largest outline
			int        scale = 4;

			int        line_width_outline = 1;
			int        line_width_corners = 3;
			Color      line_color_outline = new Color(0, 255, 0);   // green
			Color      line_color_corners = new Color(0, 255, 255); // cyan

			// generating scene outlines for results documentation
			ImagePlus imp_outlines = generateSceneOutlines(
					ref_scene,          // QuadCLT    ref_scene, // ordered by increasing timestamps
					scenes,             // QuadCLT [] scenes
					extra,              // int        extra // add around largest outline
					scale,              // int        scale,
					line_width_outline, // int        line_width_outline,
					line_width_corners, // int        line_width_corners,
					line_color_outline, // Color      line_color_outline,
					line_color_corners  // Color      line_color_corners
					);		
			imp_outlines.show();
		}
		final double [][] combo_dsn_change = new double [combo_dsn_titles.length] [tilesX*tilesY];
		for (int i = 0; i < combo_dsn.length; i++) { // 4 elements: "disp", "strength","disp_lma","num_valid"
			if (combo_dsn[i] != null) combo_dsn_change[i] = combo_dsn[i]; // all but change
		}
		if (noise_sigma_level != null) { // add initial offset to the expected disparity
			for (int i = 0; i < combo_dsn_change[0].length; i++) {
				combo_dsn_change[combo_dsn_indx_disp][i] += noise_sigma_level.initial_offset; //initial offset
			}
		}
//		combo_dsn_change[combo_dsn_change.length - 1] = new double [tilesX*tilesY];
		if (debug_level > 0) {
			(new ShowDoubleFloatArrays()).showArrays(
					combo_dsn_change,
					tilesX,
					tilesY,
					true,
					"combo_dsn-initial"+ref_scene.getImageName(),
					combo_dsn_titles); //	dsrbg_titles);
		}
		
		final int max_refines = 10;
//		final String [] iter_titles = {"disp", "diff", "strength"};
		final int last_slices = combo_dsn_titles.length;
		final int last_initial_slices = last_slices + initial_indices.length;
		
//		final double [][] refine_results = new double [last_slices + 3 * (max_refines + 1)][];
		final double [][] refine_results = new double [last_slices + 4 * (max_refines + 1)][];
		String [] refine_titles = new String [refine_results.length];
		for (int i = 0; i < combo_dsn_titles.length; i++) {
			refine_results[i] = combo_dsn_change[i]; // first 5 - references to 5-element combo_dsn_change
			refine_titles[i] = combo_dsn_titles[i]+"-last"; // "disp", "strength","disp_lma","num_valid","change"
		}
		for (int i = 0; i < initial_indices.length; i++) {
			refine_titles[last_slices + i] = combo_dsn_titles[initial_indices[i]]+"-initial"; // "disp", "strength","num_valid"
			if (combo_dsn_change[initial_indices[i]] != null) {
				refine_results[last_slices + i] = combo_dsn_change[initial_indices[i]].clone();
			} else {
				refine_results[last_slices + i] = new double [tilesX * tilesY];
			}
		}
		for (int nrefine = 0; nrefine < max_refines; nrefine++) {
			for (int i = 0; i < iter_indices.length; i++) {
				refine_titles[last_initial_slices + i * max_refines + nrefine ] = combo_dsn_titles[iter_indices[i]]+"-"+nrefine;
			}
		}		
		for (int nrefine = 0; nrefine < max_refines; nrefine++) {
///			Runtime.getRuntime().gc();
///			System.out.println("--- Free memory="+Runtime.getRuntime().freeMemory()+" (of "+Runtime.getRuntime().totalMemory()+")");
			int mcorr_sel = Correlation2d.corrSelEncode(clt_parameters.img_dtt,scenes[indx_ref].getNumSensors());
			// FIXME: 							null,           // final boolean []     selection, // may be null, if not null do not  process unselected tiles
			double [][] disparity_map = 
					correlateInterscene(
							clt_parameters, // final CLTParameters  clt_parameters,
							scenes,         // final QuadCLT []     scenes,
							indx_ref,       // final int            indx_ref,
							combo_dsn_change[0],   // final double []      disparity_ref,  // disparity in the reference view tiles (Double.NaN - invalid)
							null,           // final boolean []     selection, // may be null, if not null do not  process unselected tiles
							margin,         // final int            margin,
							nrefine,        // final int            nrefine, // just for debug title
							( nrefine == (max_refines - 1)) && clt_parameters.inp.show_final_2d, // final boolean        show_2d_corr,
					        mcorr_sel,      // final int            mcorr_sel, //  = 							
							null,           // final float [][][]   accum_2d_corr, // if [1][][] - return accumulated 2d correlations (all pairs)
							false,          // final boolean        no_map, // do not generate disparity_map (time-consuming LMA)
							debug_level-5);   // final int            debug_level)

///			Runtime.getRuntime().gc();
///			System.out.println("--- Free memory="+Runtime.getRuntime().freeMemory()+" (of "+Runtime.getRuntime().totalMemory()+")");
			
			if (debug_level > 0) {
				(new ShowDoubleFloatArrays()).showArrays(
						disparity_map,
						tilesX,
						tilesY,
						true,
						"accumulated_disparity_map-"+nrefine,
						ImageDtt.getDisparityTitles(ref_scene.getNumSensors(),ref_scene.isMonochrome()) // ImageDtt.DISPARITY_TITLES
						);
			}
			// update disparities
			double [] map_disparity =     disparity_map[ImageDtt.DISPARITY_INDEX_CM]; // 2
			double [] map_strength =      disparity_map[ImageDtt.DISPARITY_STRENGTH_INDEX]; // 10
			double [] map_disparity_lma = disparity_map[ImageDtt.DISPARITY_INDEX_POLY]; // 8
			for (int nTile =0; nTile < combo_dsn_change[0].length; nTile++) {
				if (!Double.isNaN(combo_dsn_change[0][nTile])) {
					if ((map_disparity_lma != null) && !Double.isNaN(map_disparity_lma[nTile])) {
						combo_dsn_change[combo_dsn_indx_change][nTile] = map_disparity_lma[nTile];
					} else if (!Double.isNaN(map_disparity[nTile])) {
						combo_dsn_change[combo_dsn_indx_change][nTile] = map_disparity[nTile] / clt_parameters.ofp.magic_scale;
					}
					if (!Double.isNaN(combo_dsn_change[combo_dsn_indx_change][nTile])) {
						combo_dsn_change[combo_dsn_indx_disp][nTile] +=     combo_dsn_change[combo_dsn_indx_change][nTile]; 
						combo_dsn_change[combo_dsn_indx_strength][nTile]  = map_strength[nTile]; // combine CM/LMA
					}
				}
			}
			// Copy disparity to sisparity_lma and just mask out  tiles with no DMA data (keep all if LMA did not run at all)  
			System.arraycopy(combo_dsn_change[combo_dsn_indx_disp], 0, combo_dsn_change[combo_dsn_indx_lma], 0, combo_dsn_change[combo_dsn_indx_disp].length); // lma
			if (map_disparity_lma != null) { 
				for (int i = 0; i < map_disparity_lma.length; i++) {
					if (Double.isNaN(map_disparity_lma[i])) {
						combo_dsn_change[combo_dsn_indx_lma][i] = Double.NaN;		
					}
				}
			}
			
			for (int i = 0; i < iter_indices.length; i++) {
				refine_results[last_initial_slices + (i * max_refines) + nrefine] = combo_dsn_change[iter_indices[i]].clone();
			}
			//					clt_parameters.inp.show_final_2d, // final boolean        show_2d_corr,
			if (debug_level >0) {
				(new ShowDoubleFloatArrays()).showArrays(
						combo_dsn_change,
						tilesX,
						tilesY,
						true,
						"combo_dsn-"+nrefine+"-"+ref_scene.getImageName(),
						combo_dsn_titles); //	dsrbg_titles);
			}
		}
		
		if (debug_level > 0) {
			(new ShowDoubleFloatArrays()).showArrays(
					refine_results,
					tilesX,
					tilesY,
					true,
					"combo-"+max_refines+"-"+ref_scene.getImageName(),
					refine_titles); //	dsrbg_titles);
		}
		//noise_sigma_level
		String rslt_suffix = "-results-nonoise";
		if (noise_sigma_level != null) {
			rslt_suffix =
//					"-results-lev_"+noise_sigma_level[0]+
					"-results-rnd_"+noise_sigma_level.scale_random+
					"-fpn_"+        noise_sigma_level.scale_fpn+
					"-sigma_"+      noise_sigma_level.sigma+ // [1]+
					"-offset"+      noise_sigma_level.initial_offset+ // [2];
					"-sensors"+     noise_sigma_level.used_sensors;
			
			if (ref_only) {
				rslt_suffix +="-nointer";
			} else {
				rslt_suffix +="-inter";
			}
			//rslt_suffix +="-mask"+clt_parameters.img_dtt.dbg_pair_mask;
		}
		rslt_suffix += (clt_parameters.correlate_lma?"-lma":"-nolma");
		if (noise_variant >= 0) {
			rslt_suffix +="-variant"+noise_variant;
		}
		
		//			int                  noise_variant, // <0 - no-variants, compatible with old code			

		ref_scene.saveDoubleArrayInModelDirectory(
				rslt_suffix,         // String      suffix,
				refine_titles,       // null,          // String []   labels, // or null
				refine_results,      // dbg_data,         // double [][] data,
				tilesX,              // int         width,
				tilesY);             // int         height)
		// save combo_dsn_change to model directory
		if (debug_level >-100) {
			return;
		}
		
		System.out.println("IntersceneAccumulate(), got previous scenes: "+sts.length);
		if (debug_level > 1) { // tested OK
			System.out.println("IntersceneAccumulate(): preparing image set...");
			int nscenes = scenes.length;
			//			int indx_ref = nscenes - 1; 
			double [][][] all_scenes_xyzatr = new double [scenes.length][][]; // includes reference (last)
			double [][][] all_scenes_ers_dt = new double [scenes.length][][]; // includes reference (last)
			all_scenes_xyzatr[indx_ref] = new double [][] {ZERO3,ZERO3};
			all_scenes_ers_dt[indx_ref] = new double [][] {
				ers_reference.getErsXYZ_dt(),
				ers_reference.getErsATR_dt()};
				
				for (int i = 0; i < nscenes; i++) if (i != indx_ref) {
					String ts = scenes[i].getImageName();
					all_scenes_xyzatr[i] = new double[][] {ers_reference.getSceneXYZ(ts),       ers_reference.getSceneATR(ts)}; 		
					all_scenes_ers_dt[i] = new double[][] {ers_reference.getSceneErsXYZ_dt(ts), ers_reference.getSceneErsATR_dt(ts)}; 		
				}
				compareRefSceneTiles(
						"" ,               // String suffix,
						true, // false,             // boolean blur_reference,
						all_scenes_xyzatr, // double [][][] scene_xyzatr, // does not include reference
						all_scenes_ers_dt, // double [][][] scene_ers_dt, // does not include reference
						scenes,            // QuadCLT [] scenes,
						8);                // int iscale) // 8
		}
		// create initial disparity map for the reference scene
		
	}
	
	
	public double [][] prepareInitialComboDS(
			final CLTParameters       clt_parameters,
			final QuadCLT []          scenes,
			final int                 indx_ref,
			final int                 debug_level)
	{
		final QuadCLT ref_scene = scenes[indx_ref];
		final int num_scenes = scenes.length+0;
		final double [][][] initial_ds = new double[num_scenes][][];
		final ErsCorrection ers_reference = ref_scene.getErsCorrection();
		for (int i = 0; i < num_scenes; i++) {
			initial_ds[i] = conditionInitialDS(
					clt_parameters, // CLTParameters  clt_parameters,
					scenes[i],      // QuadCLT        scene,
					-1);            // int debug_level);
		}
		if (debug_level > -1) {
			String []   dbg_titles = new String [2*num_scenes];
			double [][] dbg_img =    new double [2*num_scenes][];
			for (int i = 0; i < num_scenes; i++) {
				dbg_titles[i + 0] =          "d-"+scenes[i].getImageName();
				dbg_titles[i + num_scenes] = "s-"+scenes[i].getImageName();
				if (initial_ds[i] != null) {
					dbg_img   [i + 0] =          initial_ds[i][0];
					dbg_img   [i + num_scenes] = initial_ds[i][1];
				}
			}
			(new ShowDoubleFloatArrays()).showArrays(
					dbg_img,
					ref_scene.getTileProcessor().getTilesX(),
					ref_scene.getTileProcessor().getTilesY(),
					true,
					"all_set_fill_dispatity_gaps"+ref_scene.getImageName(),
					dbg_titles); //	dsrbg_titles);
		}

		final int         margin = 8;
		final double      tolerance_ref_absolute = 0.2; // how much scene disparity may diverge from predicted from reference
		final double      tolerance_ref_relative = 0.02;
		final double      tolerance_scene_absolute = 0.2;   // how much 4 bi-linear interpolation corners may differ (remove extremes while not)
		final double      tolerance_scene_relative = 0.02;
		double [][][] initial_toref_ds = new double[num_scenes][][];
		initial_toref_ds[indx_ref] = initial_ds[indx_ref];
		// reference camera ERS should be already set
		for (int i = 0; i < num_scenes; i++) if ((i != indx_ref) && (initial_ds[i] != null)) {
			String ts = scenes[i].getImageName();
			double []   scene_xyz = ers_reference.getSceneXYZ(ts);
			double []   scene_atr = ers_reference.getSceneATR(ts);
			if ((scene_xyz != null) && (scene_atr != null)) { // skip mission scenes (now all but reference)
				double []   scene_ers_xyz_dt = ers_reference.getSceneErsXYZ_dt(ts);
				double []   scene_ers_atr_dt = ers_reference.getSceneErsATR_dt(ts);
				initial_toref_ds[i] = getSceneDisparityStrength(
						true,                      // final boolean     to_ref_disparity, // false - return scene disparity, true - convert disparity back to the reference scene
						initial_ds[indx_ref][0],   // final double []   disparity_ref,   // invalid tiles - NaN in disparity
						initial_ds[i][0],          // final double []   disparity_scene, // invalid tiles - NaN in disparity (just for masking out invalid scene tiles)
						initial_ds[i][1],          // initial_ds[i][1],          // final double []   strength_scene,  // to calculate interpolated strength
						scene_xyz,                 // final double []   scene_xyz, // camera center in world coordinates
						scene_atr,                 // final double []   scene_atr, // camera orientation relative to world frame
						scene_ers_xyz_dt,          // final double []   scene_ers_xyz_dt, // camera ERS linear
						scene_ers_atr_dt,          // final double []   scene_ers_atr_dt, // camera ERS linear
						scenes[i],                 // final QuadCLT     scene_QuadClt,
						scenes[indx_ref],          // final QuadCLT     reference_QuadClt,
						margin, // final int         margin,
						tolerance_ref_absolute,    // final double      tolerance_ref_absolute,
						tolerance_ref_relative,    // final double      tolerance_ref_relative,
						tolerance_scene_absolute,  // final double      tolerance_scene_absolute, // of 4 bi-linear interpolation corners
						tolerance_scene_relative); // final double      tolerance_scene_relative)  // of 4 bi-linear interpolation corners
			}
		}		
		if (debug_level > -1) { // **** Used to create images for report !
			String []   dbg_titles = new String [2*num_scenes];
			double [][] dbg_img =    new double [2*num_scenes][];
			for (int i = 0; i < num_scenes; i++) {
				dbg_titles[i + 0] =          "d-"+scenes[i].getImageName();
				dbg_titles[i + num_scenes] = "s-"+scenes[i].getImageName();
				if (initial_toref_ds[i] != null) {
					dbg_img   [i + 0] =          initial_toref_ds[i][0];
					dbg_img   [i + num_scenes] = initial_toref_ds[i][1];
				}
			}
			(new ShowDoubleFloatArrays()).showArrays(
					dbg_img,
					ref_scene.getTileProcessor().getTilesX(),
					ref_scene.getTileProcessor().getTilesY(),
					true,
					"mapped_scenes_disparity_strength"+ref_scene.getImageName(),
					dbg_titles); //	dsrbg_titles);
		}
		//ref_scene
		final int tilesX = ref_scene.getTileProcessor().getTilesX();
		final int tilesY = ref_scene.getTileProcessor().getTilesY();
		final int tiles =tilesX * tilesY;

		final double zero_strength = 0.05; // add where disparity exists and strength is 0.0
		final double [][] combo_dsn = new double[3][tiles];
		Arrays.fill(combo_dsn[0], Double.NaN);
		final Thread[] threads = ImageDtt.newThreadArray(threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				public void run() {
					for (int nTile = ai.getAndIncrement(); nTile < tiles; nTile = ai.getAndIncrement()) {
						int ndef = 0;
						double sw = 0.0, swd = 0.0;
						for (int i = 0; i < initial_toref_ds.length; i++) if (initial_toref_ds[i] != null){
							double d = initial_toref_ds[i][0][nTile];
							double w = initial_toref_ds[i][1][nTile]+zero_strength;
							if (!Double.isNaN(d)) {
								sw += w;
								swd += w * d;
								ndef++;
							}
							if (ndef > 0) {
								combo_dsn[0][nTile] = swd/sw;
								combo_dsn[1][nTile] = sw/ndef;
								combo_dsn[2][nTile] = 1.0*ndef/initial_toref_ds.length;
							}
						}
					}
				}
			};
		}		      
		ImageDtt.startAndJoin(threads);
		return combo_dsn;
	}
	
	
	@Deprecated
public double[][] correlateIntersceneDebug( // only uses GPU and quad
			final CLTParameters  clt_parameters,
			final QuadCLT []     scenes,
			final int            indx_ref,
			final double []      disparity_ref,  // disparity in the reference view tiles (Double.NaN - invalid)
			final int            margin,
			final int            nrefine, // just for debug title
			final int            debug_level
			)
	{
		// Debug images if debug_level > -2
		final boolean bayer_no_mix_hor_ver = true;
		final boolean debug_individual_dhv = true;
		final boolean debug_local_maps =     true;
		final boolean debug_local_offsets =  true;
		final boolean debug_map_cpu =        true;
		final double fat_zero_pre = (debug_level>-100)?-1.0:0.0; //100000.0; // 10000.0; // 1000.0; //
		final double output_amplitude = 30000; // 10000; // 1000; // 200; // 50.0;
		final int num_scenes = scenes.length;
		final QuadCLT ref_scene = scenes[indx_ref];
		final ErsCorrection ers_reference = ref_scene.getErsCorrection();
		final int tilesX = ref_scene.getTileProcessor().getTilesX();
		final int tilesY = ref_scene.getTileProcessor().getTilesY();
		final int tiles =tilesX * tilesY;
		float  [][][][][]     fcorrs_td =       new float[num_scenes][tilesY][tilesX][][];
		float  [][][][][]     fcorrs_combo_td = new float[num_scenes][4][tilesY][tilesX][];
		final double [][] debug_dhv =       debug_individual_dhv? (new double [3 *(num_scenes+1)][]):null;
		final double [][] debug_local_dhv = debug_local_maps?     (new double [3 *(num_scenes+1)][]):null;
		final double [][] debug_target =    debug_local_maps?     (new double [num_scenes+1][tilesX*tilesY]):null;
		final double debug_hist_max_disparity = 50.0; // for hor/vert vs. disparity
		final double debug_hist_disparity_step = 0.1;
		final double debug_hist_max_error =      1.5;
		final double debug_hist_error_step =     0.025;
		final int    debug_hist_gap =            3; // pixel rows
		final double [] debug_disparity_range = {15.0, 40.0}; // for hor/vert vs. hor/vert offsets
//		final double [] debug_disparity_range = {0.0,  50.0}; // for hor/vert vs. hor/vert offsets
		final double debug_disparity_bias =     0.0; // 2.0; // 0.0; //  0.1; //  0.0; // 0.1; // add to target_disparity
		if (debug_local_maps) {
			for (int i = 0; i < debug_target.length; i++) {
				Arrays.fill(debug_target[i], Double.NaN);
			}
		}
		final int num_offsets = 8;
		final double [][] debug_offsets =  debug_local_offsets? (new double [(num_scenes+1)*num_offsets][tilesX*tilesY]):null;
		if (debug_local_offsets) {
			for (int i = 0; i < debug_offsets.length; i++) {
				Arrays.fill(debug_offsets[i], Double.NaN);
			}
		}
		
		ImageDtt image_dtt = new ImageDtt(
				numSens,
				clt_parameters.transform_size,
				clt_parameters.img_dtt,
				ref_scene.isAux(),
				ref_scene.isMonochrome(),
				ref_scene.isLwir(),
				clt_parameters.getScaleStrength(ref_scene.isAux()),
				ref_scene.getGPU());
		double[][] disparity_map = new double [image_dtt.getDisparityTitles().length][];

		int disparity_modes = 
				ImageDtt.BITS_ALL_DISPARITIES |
				ImageDtt.BITS_ALL_DIFFS | // needs max_diff?
				ImageDtt.BITS_OVEREXPOSED; //  |
//		final Rectangle tile_woi = new Rectangle(127,2,40,60); // for visualizations
//		final Rectangle tile_woi = new Rectangle(81,77,64,97); // for visualizations
//		final Rectangle tile_woi = new Rectangle(30,70,50,110); // for visualizations
		final Rectangle tile_woi = new Rectangle(250,114,32,24); // for visualizations
		final int vis_gap = 2;
		final float [][] vis_corr_td = (debug_level > -2) ? (new float[num_scenes + 1][]) : null; // transform-domain visualization
		final float [][] vis_corr_pd = (debug_level > -2) ? (new float[num_scenes + 2][]) : null; // pixel-domain visualization
		final int [] wis_wh = new int [2];
		final float [][][][] fclt_corrs = new float [num_scenes+1][tilesX*tilesY][][]; // will only contain tile_woi tiles to save memory
		final double disparity_corr = debug_disparity_bias; // 0.0; // (z_correction == 0) ? 0.0 : geometryCorrection.getDisparityFromZ(1.0/z_correction);
		for (int nscene = 0; nscene < num_scenes; nscene++) {
			if (nscene == indx_ref) {
				System.out.println("Correlating reference scene, nrefine = "+nrefine);
			} else {
				System.out.println("Correlating scene "+nrefine+":"+nscene);
			}
			String ts = scenes[nscene].getImageName();
			double [][] scene_pXpYD;
			if (nscene == indx_ref) {
				// transform to self - maybe use a method that sets central points
				scene_pXpYD = transformToScenePxPyD(
						null, // final Rectangle [] extra_woi,    // show larger than sensor WOI (or null)
						disparity_ref,      // final double []   disparity_ref, // invalid tiles - NaN in disparity (maybe it should not be masked by margins?)
						ZERO3,              // final double []   scene_xyz, // camera center in world coordinates
						ZERO3,              // final double []   scene_atr, // camera orientation relative to world frame
						ref_scene,          // final QuadCLT     scene_QuadClt,
						ref_scene);         // final QuadCLT     reference_QuadClt)

			} else {
				double []   scene_xyz = ers_reference.getSceneXYZ(ts);
				double []   scene_atr = ers_reference.getSceneATR(ts);
				double []   scene_ers_xyz_dt = ers_reference.getSceneErsXYZ_dt(ts);
				double []   scene_ers_atr_dt = ers_reference.getSceneErsATR_dt(ts);
				scenes[nscene].getErsCorrection().setErsDt(
						scene_ers_xyz_dt, // double []    ers_xyz_dt,
						scene_ers_atr_dt); // double []    ers_atr_dt)(ers_scene_original_xyz_dt);
				//setupERS() will be inside transformToScenePxPyD()
				scene_pXpYD = transformToScenePxPyD( // will be null for disparity == NaN
						null, // final Rectangle [] extra_woi,    // show larger than sensor WOI (or null)
						disparity_ref,      // final double []   disparity_ref, // invalid tiles - NaN in disparity (maybe it should not be masked by margins?)
						scene_xyz,          // final double []   scene_xyz, // camera center in world coordinates
						scene_atr,          // final double []   scene_atr, // camera orientation relative to world frame
						scenes[nscene],      // final QuadCLT     scene_QuadClt,
						ref_scene); // final QuadCLT     reference_QuadClt)
			}
			if (scenes[nscene].getGPU().getQuadCLT() != scenes[nscene]) {
				scenes[nscene].getGPU().updateQuadCLT(scenes[nscene]); // to re-load new set of Bayer images to the GPU
			}
			final double gpu_sigma_corr =     clt_parameters.getGpuCorrSigma(scenes[nscene].isMonochrome());
			final double gpu_sigma_rb_corr =  scenes[nscene].isMonochrome()? 1.0 : clt_parameters.gpu_sigma_rb_corr;
			final double gpu_sigma_log_corr = clt_parameters.getGpuCorrLoGSigma(scenes[nscene].isMonochrome());
			final float  [][][]       fclt_corr = new float [tilesX * tilesY][][];
			image_dtt.quadCorrTD(
					clt_parameters.img_dtt,        // final ImageDttParameters  imgdtt_params,   // Now just extra correlation parameters, later will include, most others
					scene_pXpYD,                   // final double [][]         pXpYD,            // per-tile array of pX,pY,disparity triplets (or nulls)
					fcorrs_td[nscene],                  // final float  [][][][]     fcorr_td,        // [tilesY][tilesX][pair][4*64] transform domain representation of 6 corr pairs
					fcorrs_combo_td[nscene],            //  final float  [][][][]     fcorr_combo_td,  // [4][tilesY][tilesX][pair][4*64] TD of combo corrs: qud, cross, hor,vert
					scenes[nscene].getErsCorrection(),  // final GeometryCorrection  geometryCorrection,
					disparity_corr,                // final double              disparity_corr,   // disparity offset at infinity
					margin,                        // final int                 margin,           // do not use tiles if their centers are closer to the edges
					clt_parameters.gpu_sigma_r,    // 0.9, 1.1
					clt_parameters.gpu_sigma_b,    // 0.9, 1.1
					clt_parameters.gpu_sigma_g,    // 0.6, 0.7
					clt_parameters.gpu_sigma_m,    //  =       0.4; // 0.7;
					gpu_sigma_rb_corr,             // final double              gpu_sigma_rb_corr, //  = 0.5; // apply LPF after accumulating R and B correlation before G, monochrome ? 1.0 : gpu_sigma_rb_corr;
					gpu_sigma_corr,                //  =    0.9;gpu_sigma_corr_m
					gpu_sigma_log_corr,            // final double              gpu_sigma_log_corr,   // hpf to reduce dynamic range for correlations
					clt_parameters.corr_red,       // +used
					clt_parameters.corr_blue,      // +used
					threadsMax,                    // final int             threadsMax,       // maximal number of threads to launch
					debug_level);                  // final int                 globalDebugLevel)
			
			// experimental, not currently used
			if (fat_zero_pre >= 0.0) {
				// FIXME: will not work with combining pairs !!!
				int num_pairs = Correlation2d.getNumPairs(numSens);
				ImageDtt.corr_td_normalize(
						fcorrs_td[nscene], // final float [][][][] fcorr_td, // will be updated
						// if 0 - fcorr_combo_td = new float[4][tilesY][tilesX][];
						// if > 0 - fcorr_td =       new float[tilesY][tilesX][num_slices][];
						num_pairs, // GPUTileProcessor.NUM_PAIRS, // final int            num_slices,
						image_dtt.transform_size,   // final int            transform_size,
						fat_zero_pre, // final double         fat_zero_abs,
						output_amplitude, // final double         output_amplitude,
						threadsMax); // final int            threadsMax);     // maximal number of threads to launch
				ImageDtt.corr_td_normalize(
						fcorrs_combo_td[nscene], // final float [][][][] fcorr_td, // will be updated
						// if 0 - fcorr_combo_td = new float[4][tilesY][tilesX][];
						// if > 0 - fcorr_td =       new float[tilesY][tilesX][num_slices][];
						0, // final int            num_slices,
						image_dtt.transform_size,   // final int            transform_size,
						fat_zero_pre, // final double         fat_zero_abs,
						output_amplitude, //final double         output_amplitude,
						threadsMax); // final int            threadsMax);     // maximal number of threads to launch
			}
			if (vis_corr_td != null) {
				vis_corr_td[nscene] = ImageDtt.corr_td_wnd(
						fcorrs_td[nscene],             // final float [][][][] fcorr_td,       // float[tilesY][tilesX][num_slices][];
						fcorrs_combo_td[nscene],       // final float [][][][] fcorr_combo_td, // float[4][tilesY][tilesX][];
						tile_woi,                 // final Rectangle      woi,
						vis_gap,                  // final int            gap,
						wis_wh,                   // final int []         wh,
						image_dtt.transform_size, // final int            transform_size,
						threadsMax);              // final int            threadsMax)     // maximal number of threads to launch
///				Runtime.getRuntime().gc();
///				System.out.println("--- Free memory="+Runtime.getRuntime().freeMemory()+" (of "+Runtime.getRuntime().totalMemory()+")");
			}
			if ((vis_corr_pd != null) || (fclt_corrs != null)) { // calculate and extract correlation
				image_dtt.clt_process_tl_correlations_GPU(	// convert to pixel domain and process correlations already prepared in fcorr_td and/or fcorr_combo_td
						clt_parameters.img_dtt,				// final ImageDttParameters  imgdtt_params,   // Now just extra correlation parameters, later will include, most others
						// both arrays should have same non-null tiles
						fcorrs_td[nscene],				 	 		// final float  [][][][]     corr_td,         // [tilesY][tilesX][pair][4*64] transform domain representation of 6 corr pairs
						fcorrs_combo_td[nscene], 						// final float  [][][][]     corr_combo_td,   // [4][tilesY][tilesX][pair][4*64] TD of combo corrs: qud, cross, hor,vert
																// each of the top elements may be null to skip particular combo type
						null, // 	final double [][][][]     corr_tiles,      // [tilesY][tilesX][pair][] ([(2*gpu_corr_rad+1)*(2*gpu_corr_rad+1)]) or null
						null, //clt_corr_partial,           			// final double [][][][][]   clt_corr_partial,// [tilesY][tilesX][quad]color][(2*transform_size-1)*(2*transform_size-1)] // if null - will not calculate
						// [tilesY][tilesX] should be set by caller
						fclt_corr, // [tile][index][(2*transform_size-1)*(2*transform_size-1)] // if null - will not calculate
						// When clt_mismatch is non-zero, no far objects extraction will be attempted
						null, 								// final double [][]         clt_mismatch,    // [12][tilesY * tilesX] // ***** transpose unapplied ***** ?. null - do not calculate
						// values in the "main" directions have disparity (*_CM) subtracted, in the perpendicular - as is
						disparity_map,						// final double [][]         disparity_map,   // [8][tilesY][tilesX], only [6][] is needed on input or null - do not calculate
						// last 2 - contrast, avg/ "geometric average)
						disparity_modes,                      // final int                 disparity_modes, // bit mask of disparity_map slices to calculate/return
						clt_parameters.gpu_corr_scale,        //  gpu_corr_scale,  //  0.75; // reduce GPU-generated correlation values
						clt_parameters.getGpuFatZero(ref_scene.isMonochrome()), // final double              gpu_fat_zero,    // clt_parameters.getGpuFatZero(is_mono);absolute == 30.0\
						image_dtt.transform_size - 1,			// clt_parameters.gpu_corr_rad,   // = transform_size - 1 ?
						clt_parameters.max_corr_radius,		// final double              max_corr_radius, // 3.9;
						clt_parameters.clt_window,     		// final int                 window_type,     // GPU: will not be used
						clt_parameters.tileX,      		    // final int               debug_tileX,
						clt_parameters.tileY,          		// final int               debug_tileY,
						threadsMax,
						debug_level -1 );
				if (vis_corr_pd != null) {
					vis_corr_pd[nscene] = ImageDtt.corr_partial_wnd( // not used in lwir
							fclt_corr,                      // clt_corr_partial,               // final double [][][][][] corr_data,
							tilesX,                         // final int               tilesX,
							2*image_dtt.transform_size - 1, // final int               corr_size,
							tile_woi,                       // final Rectangle      woi,
							vis_gap,                        // final int            gap,
							wis_wh,                         // final int []         wh,
							threadsMax);                    // final int               threadsMax)
				}
				if (fclt_corrs != null) {
					fclt_corrs[nscene] = ImageDtt.extract_corr_woi(
							true, // final boolean      copy, // copy tiles stack, not reference
							fclt_corr,   // final float [][][] fcorr,
							tile_woi,    // final Rectangle    woi,
							tilesX,      // final int          tilesX,
							threadsMax); // final int          threadsMax)     // maximal number of threads to launch
							
				}
///				Runtime.getRuntime().gc();
///				System.out.println("--- Free memory="+Runtime.getRuntime().freeMemory()+" (of "+Runtime.getRuntime().totalMemory()+")");
				//nscene
				if (debug_dhv != null) {
					debug_dhv[3 * nscene + 0] = disparity_map[ImageDtt.DISPARITY_INDEX_CM];
					debug_dhv[3 * nscene + 1] = disparity_map[ImageDtt.DISPARITY_INDEX_HOR];
					debug_dhv[3 * nscene + 2] = disparity_map[ImageDtt.DISPARITY_INDEX_VERT];
				}
			}
			if (debug_local_dhv != null) {
				double [] disparity = scenes[nscene].getDSRBG()[0];
				double [][] local_map;
				if (debug_map_cpu) {
					local_map = scenes[nscene].CLTCorrDisparityMapCPUTasks(
							  "-"+String.format("%02d", (nscene+1)), // String            suffix
							  clt_parameters,        // CLTParameters     clt_parameters,
							  debug_disparity_range, // final double []   debug_disparity_range,
							  disparity,             // double []         disparity, // do not measure NaN
							  disparity_corr,        // double disparity_corr,
							  threadsMax,            // final int         threadsMax,  // maximal number of threads to launch
							  true,                  // final boolean     updateStatus,
							  debug_level -1 );      // final int         debugLevel);
				} else {
					local_map = scenes[nscene].CLTCorrDisparityMap(
							  clt_parameters,   // CLTParameters     clt_parameters,
							  disparity,        // double []         disparity, // do not measure NaN
							  disparity_corr,   // double disparity_corr,
							  threadsMax,       // final int         threadsMax,  // maximal number of threads to launch
							  true,             // final boolean     updateStatus,
							  debug_level -1 ); // final int         debugLevel);
				}
				debug_local_dhv[3 * nscene + 0] = local_map[ImageDtt.DISPARITY_INDEX_CM];
				debug_local_dhv[3 * nscene + 1] = local_map[ImageDtt.DISPARITY_INDEX_HOR];
				debug_local_dhv[3 * nscene + 2] = local_map[ImageDtt.DISPARITY_INDEX_VERT];
				debug_target[nscene] = disparity;
				if (debug_offsets != null) {
					double [][] scene_offsets = scenes[nscene].getPortsXY();
					for (int nt = 0; nt< scene_offsets.length; nt++) if (scene_offsets[nt] != null) {
						for (int i = 0; i < num_offsets; i++) {
							debug_offsets[nscene * num_offsets + i][nt] = scene_offsets[nt][i];
						}
					}
				}
			}
		}
		// Combine non-null correlations from all scenes (initially combined and individual for visualization and analysis)
		final float  [][][][]     fcorr_td =       new float[tilesY][tilesX][][];
		final float  [][][][]     fcorr_combo_td = new float[4][tilesY][tilesX][];
		final int first_combo = 0;
		final Thread[] threads = ImageDtt.newThreadArray(threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				public void run() {
					for (int nTile = ai.getAndIncrement(); nTile < tiles; nTile = ai.getAndIncrement()) {
						int tileY = nTile / tilesX;  
						int tileX = nTile % tilesX;
						int num_indiv = 0;
						int num_combo = 0;
						for (int i = 0; i < num_scenes; i++) {
							if (fcorrs_td[i][tileY][tileX] != null) {
								if (fcorr_td[tileY][tileX] == null) {
									fcorr_td[tileY][tileX] = new float [fcorrs_td[i][tileY][tileX].length][];
									for (int j = 0; j < fcorr_td[tileY][tileX].length; j++) {
										fcorr_td[tileY][tileX][j] = fcorrs_td[i][tileY][tileX][j].clone(); 
									}
								} else {
									for (int j = 0; j < fcorr_td[tileY][tileX].length; j++) {
										for (int k = 0; k < fcorr_td[tileY][tileX][j].length; k++) {
											fcorr_td[tileY][tileX][j][k] += fcorrs_td[i][tileY][tileX][j][k];
										}
									}
								}
								num_indiv++;
							}
							
							if (fcorrs_combo_td[i][first_combo][tileY][tileX] != null) {
								if (fcorr_combo_td[first_combo][tileY][tileX] == null) {
									for (int p = 0; p < fcorr_combo_td.length; p++) {
										if (fcorrs_combo_td[i][p][tileY][tileX]!=null) {
											fcorr_combo_td[p][tileY][tileX] = fcorrs_combo_td[i][p][tileY][tileX].clone(); 
										}
									}
								} else {
									for (int p = 0; p < fcorr_combo_td.length; p++) {
										for (int k = 0; k < fcorr_combo_td[p][tileY][tileX].length; k++) {
											fcorr_combo_td[p][tileY][tileX][k] += fcorrs_combo_td[i][p][tileY][tileX][k]; 
										}
									}
								}
								num_combo++;
							}
						}
						if (num_indiv > 0) {
							double s = 1.0/num_indiv;
							for (int j = 0; j < fcorr_td[tileY][tileX].length; j++) {
								for (int k = 0; k < fcorr_td[tileY][tileX][j].length; k++) {
									fcorr_td[tileY][tileX][j][k] *= s;
								}
							}
						}
						if (num_combo > 0) {
							double s = 1.0/num_combo;
							for (int p = 0; p < fcorr_combo_td.length; p++) {
								for (int k = 0; k < fcorr_combo_td[p][tileY][tileX].length; k++) {
									fcorr_combo_td[p][tileY][tileX][k] *= s; 
								}
							}
						}
					}
				}
			};
		}		      
		ImageDtt.startAndJoin(threads);
//		// TODO: Visualize each correlation - individual and combo
		// testing: overwrite with reference frame correlations
		if (debug_level < -100) { 
			for (int i = 0; i < fcorr_combo_td.length; i++) {
				fcorr_combo_td[i] = fcorrs_combo_td[indx_ref][i].clone();
			}
			for (int i = 0; i < fcorr_td.length; i++) {
				fcorr_td[i] = fcorrs_td[indx_ref][i].clone();
			}
		}
		for (int i = 0; i < num_scenes; i++) {
			fcorrs_td[i] = null;
			fcorrs_combo_td[i] = null;
		}
		if (debug_level > -2){ // -1
			int indx_corr = -1;
			indx_corr = -1;
			final float [][][][] fcorr_td_dbg =       (indx_corr < 0) ? fcorr_td : fcorrs_td[indx_corr];
			final float [][][][] fcorr_combo_td_dbg = (indx_corr < 0) ? fcorr_combo_td : fcorrs_combo_td[indx_corr];
			int [] wh = new int[2];
			// FIXME: will not work with combining pairs !!!
			int num_pairs = Correlation2d.getNumPairs(numSens);
			float [][] dbg_corr = ImageDtt.corr_td_dbg(
					fcorr_td_dbg,               // final float [][][][] fcorr_td,
					// if 0 - fcorr_combo_td = new float[4][tilesY][tilesX][];
					// if > 0 - fcorr_td =       new float[tilesY][tilesX][num_slices][];
					num_pairs,                  // final int            num_slices,
					image_dtt.transform_size,   // final int            transform_size,
					wh,                         // final int []         wh, // should be initialized as int[2];
					threadsMax); // final int            threadsMax)     // maximal number of threads to launch

			float [][] dbg_corr_combo = ImageDtt.corr_td_dbg(
					fcorr_combo_td_dbg,         // final float [][][][] fcorr_td,
					// if 0 - fcorr_combo_td = new float[4][tilesY][tilesX][];
					// if > 0 - fcorr_td =       new float[tilesY][tilesX][num_slices][];
					0,                          // final int            num_slices,
					image_dtt.transform_size,   // final int            transform_size,
					null,                         // final int []         wh, // should be initialized as int[2];
					threadsMax); // final int            threadsMax)     // maximal number of threads to launch
			
			float [][] dbg_img = new float[dbg_corr.length + dbg_corr_combo.length][];
			for (int i = 0; i < dbg_corr.length; i++) {
				dbg_img[i] = dbg_corr[i];
			}
			for (int i = 0; i < dbg_corr_combo.length; i++) {
				dbg_img[i + dbg_corr.length] = dbg_corr_combo[i];
			}
			
			// titles.length = 15, corr_rslt_partial.length=16!
			String [] dbg_titles = new String [dbg_img.length];
			System.arraycopy(ImageDtt.CORR_TITLES, 0, dbg_titles, 0, dbg_titles.length);
			(new ShowDoubleFloatArrays()).showArrays( // out of boundary 15
					dbg_img,
					wh[0],
					wh[1],
					true,
					ref_scene.getImageName()+"-TD-CORR-"+nrefine,
					dbg_titles);
		}
		
///		Runtime.getRuntime().gc();
///		System.out.println("--- Free memory="+Runtime.getRuntime().freeMemory()+" (of "+Runtime.getRuntime().totalMemory()+")");
		final float  [][][]       fclt_corr = new float [tilesX * tilesY][][];
		image_dtt.clt_process_tl_correlations_GPU(	// convert to pixel domain and process correlations already prepared in fcorr_td and/or fcorr_combo_td
				clt_parameters.img_dtt,				// final ImageDttParameters  imgdtt_params,   // Now just extra correlation parameters, later will include, most others
				// both arrays should have same non-null tiles
				fcorr_td,				 	 		// final float  [][][][]     corr_td,         // [tilesY][tilesX][pair][4*64] transform domain representation of 6 corr pairs
				fcorr_combo_td, 						// final float  [][][][]     corr_combo_td,   // [4][tilesY][tilesX][pair][4*64] TD of combo corrs: qud, cross, hor,vert
				// each of the top elements may be null to skip particular combo type
				null, // 	final double [][][][]     corr_tiles,      // [tilesY][tilesX][pair][] ([(2*gpu_corr_rad+1)*(2*gpu_corr_rad+1)]) or null
				null, // clt_corr_partial,           			// final double [][][][][]   clt_corr_partial,// [tilesY][tilesX][quad]color][(2*transform_size-1)*(2*transform_size-1)] // if null - will not calculate
				// [tilesY][tilesX] should be set by caller
				fclt_corr, // [tile][index][(2*transform_size-1)*(2*transform_size-1)] // if null - will not calculate
				// When clt_mismatch is non-zero, no far objects extraction will be attempted
				null, 								// final double [][]         clt_mismatch,    // [12][tilesY * tilesX] // ***** transpose unapplied ***** ?. null - do not calculate
				// values in the "main" directions have disparity (*_CM) subtracted, in the perpendicular - as is
				disparity_map,						// final double [][]         disparity_map,   // [8][tilesY][tilesX], only [6][] is needed on input or null - do not calculate
				// last 2 - contrast, avg/ "geometric average)
				disparity_modes,                      // final int                 disparity_modes, // bit mask of disparity_map slices to calculate/return
				clt_parameters.gpu_corr_scale,        //  gpu_corr_scale,  //  0.75; // reduce GPU-generated correlation values
				clt_parameters.getGpuFatZero(ref_scene.isMonochrome()), // final double              gpu_fat_zero,    // clt_parameters.getGpuFatZero(is_mono);absolute == 30.0\
				image_dtt.transform_size - 1,			// clt_parameters.gpu_corr_rad,   // = transform_size - 1 ?
				clt_parameters.max_corr_radius,		// final double              max_corr_radius, // 3.9;
				clt_parameters.clt_window,     		// final int                 window_type,     // GPU: will not be used
				clt_parameters.tileX,      		    // final int               debug_tileX,
				clt_parameters.tileY,          		// final int               debug_tileY,
				threadsMax,
				debug_level -1 );
		if (vis_corr_pd != null) { // add combined data as the last slice
			vis_corr_pd[num_scenes] = ImageDtt.corr_partial_wnd( // not used in lwir
					fclt_corr,                      // clt_corr_partial,               // final double [][][][][] corr_data,
					tilesX,                         // 
					2*image_dtt.transform_size - 1, // final int               corr_size,
					tile_woi,                       // final Rectangle      woi,
					vis_gap,                        // final int            gap,
					wis_wh,                         // final int []         wh,
					threadsMax);                    // final int               threadsMax)				
		}

		if (debug_dhv != null) {
			debug_dhv[3 * num_scenes + 0] = disparity_map[ImageDtt.DISPARITY_INDEX_CM];
			debug_dhv[3 * num_scenes + 1] = disparity_map[ImageDtt.DISPARITY_INDEX_HOR];
			debug_dhv[3 * num_scenes + 2] = disparity_map[ImageDtt.DISPARITY_INDEX_VERT];
			debug_target[num_scenes] = disparity_ref;

			String [] dhv_mode_titles= {"disparity","hor-disparity","vert-disparity"};
			String [] dbg_titles = new String[debug_dhv.length];
			String [] dbg_titles1 = new String[debug_dhv.length/dhv_mode_titles.length];
			double [][][]dhv_mode_data = new double [dhv_mode_titles.length][debug_dhv.length/dhv_mode_titles.length][];
			for (int nscene = 0; nscene <= num_scenes; nscene++) {
				String ts = (nscene == num_scenes)? "interscene-combo":scenes[nscene].getImageName();
				dbg_titles[3 * nscene + 0] = "disp_"+ts;
				dbg_titles[3 * nscene + 1] = "hor_"+ts;
				dbg_titles[3 * nscene + 2] = "vert_"+ts;
				dbg_titles1[nscene] = ts;
				for (int i = 0; i < dhv_mode_data.length; i++) {
					dhv_mode_data[i][nscene] = debug_dhv[dhv_mode_data.length * nscene + i];
				}
			}
			(new ShowDoubleFloatArrays()).showArrays(
					debug_dhv,
					tilesX,
					tilesY,
					true,
					"dispaity_hor_vert_"+nrefine+"-d"+disparity_corr,
					dbg_titles
					);
			for (int mode = 0; mode < dhv_mode_titles.length; mode++) {
				(new ShowDoubleFloatArrays()).showArrays(
						dhv_mode_data[mode],
						tilesX,
						tilesY,
						true,
						dhv_mode_titles[mode]+"_"+nrefine+"-d"+disparity_corr,
						dbg_titles1
						);
				
			}
			(new ShowDoubleFloatArrays()).showArrays(
					debug_target,
					tilesX,
					tilesY,
					true,
					"target_disparity_"+nrefine+"-d"+disparity_corr,
					dbg_titles1
					);
			
			
			double [][] debug_dtarget_dx = new double [debug_target.length][debug_target[0].length];
			double [][] debug_dtarget_dy = new double [debug_target.length][debug_target[0].length];
			for (int s = 0; s < debug_target.length; s++) {
				Arrays.fill(debug_dtarget_dx[s],Double.NaN);
				Arrays.fill(debug_dtarget_dy[s],Double.NaN);
				for (int ty = 1; ty < (tilesY - 1); ty++) {
					for (int tx = 1; tx < (tilesX - 1); tx++) {
						int nt = ty * tilesX + tx;
						debug_dtarget_dx[s][nt] = 0.5 * (debug_target[s][nt + 1]      - debug_target[s][nt - 1]);
						debug_dtarget_dy[s][nt] = 0.5 * (debug_target[s][nt + tilesX] - debug_target[s][nt - tilesX]);
					}
				}
			}
			/*
			(new ShowDoubleFloatArrays()).showArrays(
					debug_dtarget_dx,
					tilesX,
					tilesY,
					true,
					"d_target_disparity_dx_"+nrefine+"-d"+disparity_corr,
					dbg_titles1
					);
			(new ShowDoubleFloatArrays()).showArrays(
					debug_dtarget_dy,
					tilesX,
					tilesY,
					true,
					"d_target_disparity_dy_"+nrefine+"-d"+disparity_corr,
					dbg_titles1
					);
			*/
			
		}
		if (debug_local_dhv != null) {
			debug_local_dhv[3 * num_scenes + 0] = disparity_map[ImageDtt.DISPARITY_INDEX_CM];
			debug_local_dhv[3 * num_scenes + 1] = disparity_map[ImageDtt.DISPARITY_INDEX_HOR];
			debug_local_dhv[3 * num_scenes + 2] = disparity_map[ImageDtt.DISPARITY_INDEX_VERT];
			
			String [] dhv_mode_titles= {"disparity","hor-disparity","vert-disparity"};
			int dhv_modes = dhv_mode_titles.length;
			String [] dbg_titles = new String[debug_local_dhv.length];
			String [] dbg_titles1 = new String[debug_local_dhv.length/dhv_mode_titles.length];
			double [][][]dhv_mode_data = new double [dhv_mode_titles.length][debug_local_dhv.length/dhv_mode_titles.length][];
			for (int nscene = 0; nscene <= num_scenes; nscene++) {
				String ts = (nscene == num_scenes)? "interscene-combo":scenes[nscene].getImageName();
				dbg_titles[3 * nscene + 0] = "disp_"+ts;
				dbg_titles[3 * nscene + 1] = "hor_"+ts;
				dbg_titles[3 * nscene + 2] = "vert_"+ts;
				dbg_titles1[nscene] = ts;
				for (int i = 0; i < dhv_modes; i++) {
					dhv_mode_data[i][nscene] = debug_local_dhv[dhv_mode_data.length * nscene + i];
				}
			}
			(new ShowDoubleFloatArrays()).showArrays(
					debug_local_dhv,
					tilesX,
					tilesY,
					true,
					"dispaity_hor_vert_local"+nrefine+"-d"+disparity_corr,
					dbg_titles
					);
			for (int mode = 0; mode < dhv_mode_titles.length; mode++) {
				(new ShowDoubleFloatArrays()).showArrays(
						dhv_mode_data[mode],
						tilesX,
						tilesY,
						true,
						dhv_mode_titles[mode]+"_local"+nrefine+"-d"+disparity_corr,
						dbg_titles1
						);
			}
            if (debug_offsets != null) {
    			String [] dbg_titles2 = new String[debug_offsets.length];
    			for (int nscene = 0; nscene <= num_scenes; nscene++) {
    				String ts = (nscene == num_scenes)? "interscene-combo":scenes[nscene].getImageName();
    				for (int i = 0; i < num_offsets; i++) {
    					dbg_titles2[nscene * num_offsets + i] = (((i & 1) == 0)?"x":"y") + (i/2) + "-" + ts;
    				}
    			}
    			(new ShowDoubleFloatArrays()).showArrays(
    					debug_offsets,
    					tilesX,
    					tilesY,
    					true,
    					"xy-offsets-local"+nrefine+"-d"+disparity_corr,
    					dbg_titles2
    					);
    			
    			double [][] dbg_hor = new double [dbg_titles1.length][tilesX*tilesY];
    			double [][] dbg_vert = new double [dbg_titles1.length][tilesX*tilesY];
    			for (int nscene = 0; nscene <= num_scenes; nscene++) {
    				for (int nt = 0; nt < dbg_hor[0].length; nt++) {
    					int indx = nscene * num_offsets; 
    					dbg_hor[nscene][nt] = -0.25 * (
    							debug_offsets[indx + 2][nt] + debug_offsets[indx + 6][nt] - 
    							debug_offsets[indx + 0][nt] - debug_offsets[indx + 4][nt]);
    					dbg_vert[nscene][nt] = -0.25 * (
    							debug_offsets[indx + 5][nt] + debug_offsets[indx + 7][nt] - 
    							debug_offsets[indx + 1][nt] - debug_offsets[indx + 3][nt]);
    				}
    			}
    			
    			(new ShowDoubleFloatArrays()).showArrays(
    					dbg_hor,
    					tilesX,
    					tilesY,
    					true,
    					"hor-offsets-local"+nrefine+"-d"+disparity_corr,
    					dbg_titles1
    					);
    			(new ShowDoubleFloatArrays()).showArrays(
    					dbg_vert,
    					tilesX,
    					tilesY,
    					true,
    					"vert-offsets-local"+nrefine+"-d"+disparity_corr,
    					dbg_titles1
    					);
    			
    			String [] dbg_titles30= {
    					"x1-x0","y1-y0",
    					"x3-x2","y3-y2",
    					"x2-x0","y2-y0",
    					"x3-x1","y3-y1"
    					};
    			final double [][] debug_offsets_diff = new double [(num_scenes+1)*dbg_titles30.length][tiles];
    			int [][] diff_pairs= {
    					{2,0},{3,1},
    					{6,4},{7,5},
    					{4,0},{5,1},
    					{6,2},{7,3}};
    			String [] dbg_titles3 = new String[debug_offsets.length];
    			for (int i = 0; i < debug_offsets_diff.length; i++) {
    				Arrays.fill(debug_offsets_diff[i], Double.NaN);
    			}

    			for (int nscene = 0; nscene <= num_scenes; nscene++) {
    				String ts = (nscene == num_scenes)? "interscene-combo":scenes[nscene].getImageName();
    				for (int i = 0; i < dbg_titles30.length; i++) {
    					dbg_titles3[nscene * num_offsets + i] = dbg_titles30[i] + "-" + ts;
    				}
    			}
    			for (int nscene = 0; nscene <= num_scenes; nscene++) {
    				for (int i = 0; i < dbg_titles30.length; i++) {
    					int ip0 = nscene * num_offsets + diff_pairs[i][0];
    					int ip1 = nscene * num_offsets + diff_pairs[i][1];
    					int indx = nscene * dbg_titles30.length + i;
    					for (int nt = 0; nt< tiles; nt++) {
    						debug_offsets_diff[indx][nt] = debug_offsets[ip0][nt] -debug_offsets[ip1][nt]; 
    					}
    				}    				
    			}
    			
    			(new ShowDoubleFloatArrays()).showArrays(
    					debug_offsets_diff,
    					tilesX,
    					tilesY,
    					true,
    					"xy-diff-offsets-local"+nrefine+"-d"+disparity_corr,
    					dbg_titles3
    					);
    			
    			
// build 2D histograms: target disparity - horizontal, hor/vert disparity - vertical
// combo, horizontal - above, vertical - below
    			double [][] diff_ranges = {{-0.2, 0.2},{-0.15, 0.8}};
    			int num_diff_bins = 1; // 5;
    			double [][] dx_dy_range = new double[2][2];
    			for (int bin = 0; bin < num_diff_bins; bin++) {
    				for (int dir = 0; dir < dx_dy_range.length; dir++) {
    					double bin_range = (diff_ranges[dir][1] - diff_ranges[dir][0])/num_diff_bins;
    					dx_dy_range[dir][0] = diff_ranges[dir][0] + bin * bin_range;
    					dx_dy_range[dir][1] = dx_dy_range[dir][0] + bin_range;
    					if (bin == 0) {
    						dx_dy_range[dir][0] = Double.NEGATIVE_INFINITY;
    					} else if (bin == (num_diff_bins - 1)) {
    						dx_dy_range[dir][1] = Double.POSITIVE_INFINITY;
    					}
    				}
    			
// show Bayer differences (hor, vert, and all?)
    			double [][] bayer_green_hor =    new double [num_scenes+1][tiles];
    			double [][] bayer_green_vert =   new double [num_scenes+1][tiles];
    			double [][] bayer_green_square = new double [num_scenes+1][tiles];
    			double [][] bayer_green_all =    new double [num_scenes+1][tiles];
    			for (int ns = 0; ns < bayer_green_hor.length; ns++) {
    				Arrays.fill(bayer_green_hor[ns],   Double.NaN);
    				Arrays.fill(bayer_green_vert[ns],  Double.NaN);
    				Arrays.fill(bayer_green_square[ns],Double.NaN);
    				Arrays.fill(bayer_green_all[ns],   Double.NaN);
    				for (int nt = 0; nt < tiles; nt++) if (!Double.isNaN(debug_offsets[ns * num_offsets][nt])){
    					double [] bph = new double[4];
    					double [] bpv = new double[4];
    					double [] bpg = new double[4];
    					for (int chn = 0; chn < bph.length; chn++) {
    						bph[chn] = debug_offsets[ns * num_offsets + 2* chn + 0][nt];
    						bpv[chn] = debug_offsets[ns * num_offsets + 2* chn + 1][nt];
    						bph[chn] -= 2 * Math.floor(bph[chn]/2); // 0<=bph<2.0
    						bpv[chn] -= 2 * Math.floor(bpv[chn]/2); // 0<=bpv<2.0
    						bpg[chn] = bph[chn]+bpv[chn]; 
    					}
    					double htp, hbp, vlp, vrp, dmp, dop;
    					if (bayer_no_mix_hor_ver) {
    						htp = (bph[1] - bph[0]) + 1; 
    						hbp = (bph[3] - bph[2]) + 1;
    						vlp = (bpv[2] - bpv[0]) + 1; 
    						vrp = (bpv[3] - bpv[1]) + 1;
    						dmp = (bpg[3] - bpg[0]) + 1; 
    						dop = (bpg[1] - bpg[2]) + 1;
    					} else {
    						htp = (bpg[1] - bpg[0]) + 1; 
    						hbp = (bpg[3] - bpg[2]) + 1;
    						vlp = (bpg[2] - bpg[0]) + 1; 
    						vrp = (bpg[3] - bpg[1]) + 1;
    						dmp = (bpg[3] - bpg[0]) + 1; 
    						dop = (bpg[1] - bpg[2]) + 1;
    					}

    					
    					htp = htp - 2*Math.floor(htp/2) - 1;
    					hbp = hbp - 2*Math.floor(hbp/2) - 1;
    					
    					vlp = vlp - 2*Math.floor(vlp/2) - 1;
    					vrp = vrp - 2*Math.floor(vrp/2) - 1;
    					
    					dmp = dmp - 2*Math.floor(dmp/2) - 1;
    					dop = dop - 2*Math.floor(dop/2) - 1;
    					
    					double htps = Math.sin(Math.PI*htp);
    					double hbps = Math.sin(Math.PI*hbp);
    					double vlps = Math.sin(Math.PI*vlp);
    					double vrps = Math.sin(Math.PI*vrp);
    					double dmps = Math.sin(Math.PI*dmp);
    					double dops = Math.sin(Math.PI*dop);
    					/*
    					if ((nt == 36178) || (nt == 37804) || (nt == 39105)) {
    						System.out.println("ns="+ns+" nt="+nt);
        					for (int chn = 0; chn < bph.length; chn++) {
        						System.out.println(
        								"  bph["+chn+"]="+bph[chn]+
        								"  bpv["+chn+"]="+bpv[chn]+
        								"  bpg["+chn+"]="+bpg[chn]);
        					}
        					System.out.println(
    								" htp="+htp+
    								" hbp="+hbp+"\n"+
    								" vlp="+vlp+
    								" vrp="+vrp+"\n"+
    								" dmp="+dmp+
    								" dop="+dop+"\n");
    						System.out.println(
    								" htps="+htps+
    								" hbps="+hbps+"\n"+
    								" vlps="+vlps+
    								" vrps="+vrps+"\n"+
    								" dmps="+dmps+
    								" dops="+dops+"\n");
    					}
    					*/
    					
    					bayer_green_hor [ns][nt] =   0.5 * (htps + hbps);
    					bayer_green_vert[ns][nt] =   0.5 * (vlps + vrps);
    					bayer_green_square[ns][nt] = 0.5 * (bayer_green_hor [ns][nt] + bayer_green_vert [ns][nt]);
    					double diag = 0.5 * (dmps + dops);
    					bayer_green_all[ns][nt] = (2 * bayer_green_square [ns][nt] + diag) / 3.0;
    				}
    			}
    			String mix = bayer_no_mix_hor_ver ? "nomix": "mix";
    			(new ShowDoubleFloatArrays()).showArrays(
    					bayer_green_hor,
    					tilesX,
    					tilesY,
    					true,
    					"bayer_green_hor-"+mix+"-"+nrefine+"-d"+disparity_corr,
    					dbg_titles1
    					);
    			(new ShowDoubleFloatArrays()).showArrays(
    					bayer_green_vert,
    					tilesX,
    					tilesY,
    					true,
    					"bayer_green_vert-"+mix+"-"+nrefine+"-d"+disparity_corr,
    					dbg_titles1
    					);
    			(new ShowDoubleFloatArrays()).showArrays(
    					bayer_green_square,
    					tilesX,
    					tilesY,
    					true,
    					"bayer_green_square-"+mix+"-"+nrefine+"-d"+disparity_corr,
    					dbg_titles1
    					);
    			(new ShowDoubleFloatArrays()).showArrays(
    					bayer_green_all,
    					tilesX,
    					tilesY,
    					true,
    					"bayer_green_all-"+mix+"-"+nrefine+"-d"+disparity_corr,
    					dbg_titles1
    					);
				String          suffix = "bin"+bin+"-pass"+nrefine+"-d"+disparity_corr;
	    		String          suffix2 = "bin"+bin+"-pass"+nrefine+"-d"+disparity_corr;
				showVHDisparityHistograms(
						ref_scene,                  // QuadCLT         ref_scene,
						debug_hist_max_disparity,   // double          debug_hist_max_disparity,
						debug_hist_disparity_step,  // double          debug_hist_disparity_step,
						debug_hist_max_error,       // double          debug_hist_max_error,
						debug_hist_error_step,      // double          debug_hist_error_step,
						disparity_corr,             // double          disparity_corr,
						dhv_modes,                  // int             dhv_modes,
						debug_hist_gap,             // int             debug_hist_gap,
						suffix,                     // String          suffix, // +nrefine+"-d"+disparity_corr,
						suffix2,                    // String          suffix2, // +nrefine+"-d"+disparity_corr,
						dbg_titles1,                // String []       dbg_titles1,
						debug_target,               // double [][]     debug_target, // target disparity
						dhv_mode_data,              // double [][][]   dhv_mode_data,
						dbg_hor,                    // double [][]     dbg_hor,
						dbg_vert,                   // double [][]     dbg_vert,
						debug_disparity_range,      // double []       debug_disparity_range,
						dx_dy_range,                // double [][]     dx_dy_range
						bayer_green_hor,            // double [][]     bayer_green_hor,
						bayer_green_vert);          // double [][]     bayer_green_vert);

			}
   			
    			
    			
    			
// Now remove integer part (x-round(x) and show again
    			for (int i = 0; i < debug_offsets.length; i++) {
        			for (int j = 0; j < debug_offsets[i].length; j++) {
        				debug_offsets[i][j] -= Math.round(debug_offsets[i][j]); 
        			}
    			}
    			for (int i = 0; i < dbg_hor.length; i++) {
        			for (int j = 0; j < dbg_hor[i].length; j++) {
        				dbg_hor[i][j] -= Math.round(dbg_hor[i][j]); 
        			}
    			}
    			for (int i = 0; i < dbg_vert.length; i++) {
        			for (int j = 0; j < dbg_vert[i].length; j++) {
        				dbg_vert[i][j] -= Math.round(dbg_vert[i][j]); 
        			}
    			}
    			(new ShowDoubleFloatArrays()).showArrays(
    					debug_offsets,
    					tilesX,
    					tilesY,
    					true,
    					"xy-offsets-fract"+nrefine+"-d"+disparity_corr,
    					dbg_titles2
    					);
    			(new ShowDoubleFloatArrays()).showArrays(
    					dbg_hor,
    					tilesX,
    					tilesY,
    					true,
    					"hor-offsets-fract"+nrefine+"-d"+disparity_corr,
    					dbg_titles1
    					);
    			(new ShowDoubleFloatArrays()).showArrays(
    					dbg_vert,
    					tilesX,
    					tilesY,
    					true,
    					"vert-offsets-fract"+nrefine+"-d"+disparity_corr,
    					dbg_titles1
    					);

    			
            }			
			
			
		}		
		
		
///		Runtime.getRuntime().gc();
///		System.out.println("--- Free memory="+Runtime.getRuntime().freeMemory()+" (of "+Runtime.getRuntime().totalMemory()+")");

		if (vis_corr_td != null) { // add combined data as the last slice
			vis_corr_td[num_scenes] = ImageDtt.corr_td_wnd(
					fcorr_td,                 // final float [][][][] fcorr_td,       // float[tilesY][tilesX][num_slices][];
					fcorr_combo_td,           // final float [][][][] fcorr_combo_td, // float[4][tilesY][tilesX][];
					tile_woi,                 // final Rectangle      woi,
					vis_gap,                  // final int            gap,
					wis_wh,                   // final int []         wh,
					image_dtt.transform_size, // final int            transform_size,
					threadsMax);              // final int            threadsMax)     // maximal number of threads to launch
		}



		if (debug_level > -2){ // -1
			float [][] dbg_corr_rslt_partial = ImageDtt.corr_partial_dbg(
					fclt_corr, // final float  [][][]     fcorr_data,       // [tile][index][(2*transform_size-1)*(2*transform_size-1)] // if null - will not calculate
					tilesX,    //final int               tilesX,
					2*image_dtt.transform_size - 1,	//final int corr_size,
					10, // final int layers 4,	// final int pairs,
					clt_parameters.corr_border_contrast,
					threadsMax,
					debug_level);
			// titles.length = 15, corr_rslt_partial.length=16!
			(new ShowDoubleFloatArrays()).showArrays( // out of boundary 15
					dbg_corr_rslt_partial,
					tilesX*(2*image_dtt.transform_size),
					tilesY*(2*image_dtt.transform_size),
					true,
					ref_scene.getImageName()+"-PD-CORR-"+nrefine,
					ImageDtt.CORR_TITLES);
		}
		if (vis_corr_td != null) {
			String [] dbg_titles = new String[num_scenes+1];
			for (int i = 0; i < num_scenes; i++) {
				dbg_titles[i] = scenes[i].getImageName();
			}
			dbg_titles[num_scenes] = "combo";
			(new ShowDoubleFloatArrays()).showArrays( // out of boundary 15
					vis_corr_td,
					wis_wh[0],
					wis_wh[1],
					true,
					"TD-"+tile_woi.x+"_"+tile_woi.y+"-"+nrefine,
					dbg_titles);
		}
		float [][][] fcorr_extra = null; 
	    if (fclt_corrs != null) {
	        fclt_corrs[num_scenes] = ImageDtt.extract_corr_woi(
					true, // final boolean      copy, // copy tiles stack, not reference
	                fclt_corr,   // final float [][][] fcorr,
	                tile_woi,    // final Rectangle    woi,
	                tilesX,      // final int          tilesX,
	                threadsMax); // final int          threadsMax)     // maximal number of threads to launch
	        fcorr_extra = ImageDtt.corr_get_extra(
	        		fclt_corrs,    // final float [][][][] fcorrs,
	    			tilesX,        // final int            tilesX,
	    			num_scenes,    // final int            ncombo,
	    			10,            // final int            slices,
	    			threadsMax);   // final int            threadsMax)
	    }

		if (vis_corr_pd != null) {
			if (fcorr_extra != null) {
				vis_corr_pd[num_scenes + 1]=ImageDtt. corr_show_extra(
						fcorr_extra,              // final float [][][] fcorr_extra,  // float[tile][slices][extra];
						tilesX,                   // final int          tilesX,
						tile_woi,                 // final Rectangle    woi,
						vis_gap,                  // final int          gap,
						3,                        // final int          step,
						2,                        // final int          size,
						null,                     // final int []       wh,
						image_dtt.transform_size, // final int          transform_size,
						threadsMax);              //  final int          threadsMax)     // maximal number of threads to launch
			}
			String [] dbg_titles = new String[num_scenes+2];
			for (int i = 0; i < num_scenes; i++) {
				dbg_titles[i] = scenes[i].getImageName();
			}
			dbg_titles[num_scenes] =     "combo";
			dbg_titles[num_scenes + 1] = "lucky";
			(new ShowDoubleFloatArrays()).showArrays( // out of boundary 15
					vis_corr_pd,
					wis_wh[0],
					wis_wh[1],
					true,
					"PD-"+"FZ-"+(clt_parameters.getGpuFatZero(ref_scene.isMonochrome()))+"-"+tile_woi.x+"_"+tile_woi.y+"-"+nrefine,
					dbg_titles);
		}		  
		return disparity_map; // disparity_map
	}
	
	public static int getEarliestScene(QuadCLT []     scenes)
	{
		int ref_index = scenes.length - 1;
		ErsCorrection ers_reference = scenes[ref_index].getErsCorrection();
		for (int nscene = ref_index-1; nscene >= 00; nscene--) {
			String ts = scenes[nscene].getImageName();
			double []   scene_xyz = ers_reference.getSceneXYZ(ts);
			double []   scene_atr = ers_reference.getSceneATR(ts);
			if ((scene_xyz == null) || (scene_atr == null)){
				return nscene + 1; // scene is not matched
			}
		}
		return 0;
	}
	
// Cleaned up and optimized version to reduce memory usage (on-the-fly integration, not saving full correlation data)
	public double[][] correlateInterscene(
			final CLTParameters  clt_parameters,
			final QuadCLT []     scenes,
			final int            indx_ref,
			final double []      disparity_ref,  // disparity in the reference view tiles (Double.NaN - invalid)
			final boolean []     selection, // may be null, if not null do not  process unselected tiles
			final int            margin,
			final int            nrefine, // just for debug title
			final boolean        show_2d_corr,
			final int            mcorr_sel, //  = Correlation2d.corrSelEncode(clt_parameters.img_dtt,scenes[nscene].getNumSensors());
			final float [][][]   accum_2d_corr, // if [1][][] - return accumulated 2d correlations (all pairs)
			final boolean        no_map, // do not generate disparity_map (time-consuming LMA)
			final int            debug_level
			)
	{
//		final int num_scenes = scenes.length;
		final QuadCLT ref_scene = scenes[indx_ref];
//		final int num_scenes = ref_scene.getEarliestScene(scenes);
		final int num_scenes = indx_ref - ref_scene.getEarliestScene(scenes) + 1;
		
		final ErsCorrection ers_reference = ref_scene.getErsCorrection();
		final int tilesX = ref_scene.getTileProcessor().getTilesX();
		final int tilesY = ref_scene.getTileProcessor().getTilesY();
		final int num_pairs = Correlation2d.getNumPairs(ref_scene.getNumSensors());
		final double [][][][][] dcorr_td_acc  = new double[num_pairs][][][][];
		final float  [][][][]   fcorr_td_acc  = new float [tilesY][tilesX][][];
//		final int    [][][]     num_acc = new int [tilesY][tilesX][num_pairs];
		final float  [][][]     num_acc = new float [tilesY][tilesX][num_pairs];
		boolean show_accumulated_correlations = show_2d_corr || debug_level > -5;
		boolean show_reference_correlations =  show_2d_corr || debug_level > -5;
		final float  [][][]       fclt_corr = ((accum_2d_corr != null) || show_accumulated_correlations || show_reference_correlations) ?
				(new float [tilesX * tilesY][][]) : null;
		ImageDtt image_dtt;
		image_dtt = new ImageDtt(
				numSens,
				clt_parameters.transform_size,
				clt_parameters.img_dtt,
				ref_scene.isAux(),
				ref_scene.isMonochrome(),
				ref_scene.isLwir(),
				clt_parameters.getScaleStrength(ref_scene.isAux()),
				ref_scene.getGPU());
		if (ref_scene.getGPU() != null) {
			ref_scene.getGPU().setGpu_debug_level(debug_level);
		}
		image_dtt.getCorrelation2d(); // initiate image_dtt.correlation2d, needed if disparity_map != null  

		double[][] disparity_map = no_map ? null : new double [image_dtt.getDisparityTitles().length][];
		final double disparity_corr = 0.00; // (z_correction == 0) ? 0.0 : geometryCorrection.getDisparityFromZ(1.0/z_correction);
		TpTask[] tp_tasks_ref = null;
		for (int nscene = 0; nscene < scenes.length; nscene++) {
			String ts = scenes[nscene].getImageName();
			double [][] scene_pXpYD;
			if (nscene == indx_ref) {
				// transform to self - maybe use a method that sets central points
				scene_pXpYD = transformToScenePxPyD(
						null, // final Rectangle [] extra_woi,    // show larger than sensor WOI (or null)
						disparity_ref,      // final double []   disparity_ref, // invalid tiles - NaN in disparity (maybe it should not be masked by margins?)
						ZERO3,              // final double []   scene_xyz, // camera center in world coordinates
						ZERO3,              // final double []   scene_atr, // camera orientation relative to world frame
						ref_scene,          // final QuadCLT     scene_QuadClt,
						ref_scene);         // final QuadCLT     reference_QuadClt)

			} else {
				double []   scene_xyz = ers_reference.getSceneXYZ(ts);
				double []   scene_atr = ers_reference.getSceneATR(ts);
				if ((scene_xyz == null) || (scene_atr == null)){
					continue; // scene is not matched
				}
				double []   scene_ers_xyz_dt = ers_reference.getSceneErsXYZ_dt(ts);
				double []   scene_ers_atr_dt = ers_reference.getSceneErsATR_dt(ts);
				scenes[nscene].getErsCorrection().setErsDt(
						scene_ers_xyz_dt, // double []    ers_xyz_dt,
						scene_ers_atr_dt); // double []    ers_atr_dt)(ers_scene_original_xyz_dt);
				//setupERS() will be inside transformToScenePxPyD()
				double [][] scene_pXpYD_prefilter = transformToScenePxPyD( // will be null for disparity == NaN, total size - tilesX*tilesY
						null, // final Rectangle [] extra_woi,    // show larger than sensor WOI (or null)
						disparity_ref,      // final double []   disparity_ref, // invalid tiles - NaN in disparity (maybe it should not be masked by margins?)
						scene_xyz,          // final double []   scene_xyz, // camera center in world coordinates
						scene_atr,          // final double []   scene_atr, // camera orientation relative to world frame
						scenes[nscene],      // final QuadCLT     scene_QuadClt,
						ref_scene); // final QuadCLT     reference_QuadClt)
				
				double max_overlap = 0.6; 
//				double min_str_cam =    0.1;
				double min_adisp_cam =  0.2;
				double min_rdisp_cam =  0.03;
				double [][] scene_ds =conditionInitialDS(
						clt_parameters, // CLTParameters  clt_parameters,
						scenes[nscene], // QuadCLT        scene,
						-1); // int debug_level);
				if (scene_ds != null) {
					double [] disparity_cam = scene_ds[0]; // null; // for now
					scene_pXpYD = filterBG (
							scenes[indx_ref].getTileProcessor(), // final TileProcessor tp,
							scene_pXpYD_prefilter, // final double [][] pXpYD,
							max_overlap,           // final double max_overlap,
							null, // disparity_cam,         // final double [] disparity_cam,
							min_adisp_cam,         // final double min_adisp_cam,
							min_rdisp_cam,         // final double min_rdisp_cam,
							clt_parameters.tileX,  // final int    dbg_tileX,
							clt_parameters.tileY,  // final int    dbg_tileY,
							0); // 1); //debug_level);          // final int    debug_level);
				} else {
					scene_pXpYD = scene_pXpYD_prefilter;
				}
			}
			if (debug_level > -1) {
				if (nscene == indx_ref) {
					System.out.println("Correlating reference scene, nrefine = "+nrefine);
				} else {
					System.out.println("Correlating scene "+nrefine+":"+nscene);
				}
			}
			scenes[nscene].saveQuadClt(); // to re-load new set of Bayer images to the GPU (do nothing for CPU)
			final double gpu_sigma_corr =     clt_parameters.getGpuCorrSigma(scenes[nscene].isMonochrome());
			final double gpu_sigma_rb_corr =  scenes[nscene].isMonochrome()? 1.0 : clt_parameters.gpu_sigma_rb_corr;
			final double gpu_sigma_log_corr = clt_parameters.getGpuCorrLoGSigma(scenes[nscene].isMonochrome());
			//			final float  [][][]       fclt_corr = new float [tilesX * tilesY][][];

			TpTask[] tp_tasks =  GpuQuad.setInterTasks(
					scenes[nscene].getNumSensors(),
					scenes[nscene].getErsCorrection().getSensorWH()[0],
					!scenes[nscene].hasGPU(), // final boolean             calcPortsCoordinatesAndDerivatives, // GPU can calculate them centreXY
					scene_pXpYD,              // final double [][]         pXpYD, // per-tile array of pX,pY,disparity triplets (or nulls)
					selection,                // final boolean []          selection, // may be null, if not null do not  process unselected tiles
					scenes[nscene].getErsCorrection(), // final GeometryCorrection  geometryCorrection,
					disparity_corr,           // final double              disparity_corr,
					margin,                   // final int                 margin,      // do not use tiles if their centers are closer to the edges
					null,                     // final boolean []          valid_tiles,            
					threadsMax);              // final int                 threadsMax)  // maximal number of threads to launch
			if (nscene == indx_ref) {
				tp_tasks_ref = tp_tasks; // will use coordinates data for LMA ? disp_dist
			}
//			int mcorr_sel = Correlation2d.corrSelEncode(clt_parameters.img_dtt,scenes[nscene].getNumSensors());
			if (scenes[nscene].hasGPU()) {
				float  [][][][]     fcorr_td =       new float[tilesY][tilesX][][];
				image_dtt.quadCorrTD(
						clt_parameters.img_dtt,            // final ImageDttParameters imgdtt_params,    // Now just extra correlation parameters, later will include, most others
						tp_tasks, // *** will be updated inside from GPU-calculated geometry
						fcorr_td, // fcorrs_td[nscene],                 // [tilesY][tilesX][pair][4*64] transform domain representation of 6 corr pairs
//						scenes[nscene].getErsCorrection(), //
						clt_parameters.gpu_sigma_r,        // 0.9, 1.1
						clt_parameters.gpu_sigma_b,        // 0.9, 1.1
						clt_parameters.gpu_sigma_g,        // 0.6, 0.7
						clt_parameters.gpu_sigma_m,        //  =       0.4; // 0.7;
						gpu_sigma_rb_corr,                 // final double              gpu_sigma_rb_corr, //  = 0.5; // apply LPF after accumulating R and B correlation before G, monochrome ? 1.0 : gpu_sigma_rb_corr;
						gpu_sigma_corr,                    //  =    0.9;gpu_sigma_corr_m
						gpu_sigma_log_corr,                // final double              gpu_sigma_log_corr,   // hpf to reduce dynamic range for correlations
						clt_parameters.corr_red,           // +used
						clt_parameters.corr_blue,          // +used
						mcorr_sel,                         // final int                 mcorr_sel,    // Which pairs to correlate // +1 - all, +2 - dia, +4 - sq, +8 - neibs, +16 - hor + 32 - vert
						threadsMax,       // maximal number of threads to launch
						debug_level);
				if (image_dtt.getGPU().getGpu_debug_level() > -1) {
					System.out.println("==ooo=after image_dtt.quadCorrTD()");
				}
// Verify tasks are now updated
				accumulateCorrelations(
						num_acc,       // final int [][][]     num_acc,          // number of accumulated tiles [tilesY][tilesX][pair]
						fcorr_td,      // final float [][][][] fcorr_td,         // [tilesY][tilesX][pair][256] sparse transform domain representation of corr pairs 
						fcorr_td_acc); // final float [][][][] fcorr_td_acc      // [tilesY][tilesX][pair][256] sparse transform domain representation of corr pairs
				if (image_dtt.getGPU().getGpu_debug_level() > -1) {
					System.out.println("==ooo=accumulateCorrelations()");
				}
				
				if ((nscene == indx_ref) && show_reference_correlations) { // prepare 2d correlations for visualization, double/CPU mode
					//clt_parameters.getGpuFatZero(ref_scene.isMonochrome()),   // final double     gpu_fat_zero,
					//Use same for both CPU/GPU
					double [][][] dcorr_tiles = (fclt_corr != null)? (new double [tp_tasks_ref.length][][]):null;
					image_dtt.clt_process_tl_correlations( // convert to pixel domain and process correlations already prepared in fcorr_td and/or fcorr_combo_td
							clt_parameters.img_dtt,		   // final ImageDttParameters  imgdtt_params,   // Now just extra correlation parameters, later will include, most others
							fcorr_td,		 	 		   // final float  [][][][]     fcorr_td,        // [tilesY][tilesX][pair][4*64] transform domain representation of all selected corr pairs
							null, // num_acc,              // int [][][]                num_acc,         // number of accumulated tiles [tilesY][tilesX][pair] (or null)       
							null, // dcorr_weight,                  // double []                 dcorr_weight,    // alternative to num_acc, compatible with CPU processing (only one non-zero enough)
							clt_parameters.gpu_corr_scale, //  final double              gpu_corr_scale,  //  0.75; // reduce GPU-generated correlation values
							clt_parameters.getGpuFatZero(ref_scene.isMonochrome()),   // final double     gpu_fat_zero,    // clt_parameters.getGpuFatZero(is_mono);absolute == 30.0
							image_dtt.transform_size - 1,  // final int                 gpu_corr_rad,    // = transform_size - 1 ?
					        // The tp_tasks data should be decoded from GPU to get coordinates
							tp_tasks_ref,                  // final TpTask []           tp_tasks,        // data from the reference frame - will be applied to LMW for the integrated correlations
							ref_scene.getErsCorrection().getRXY(false), // final double [][]         rXY,             // from geometryCorrection
							// next both can be nulls
							null,                          // final double [][][][]     clt_corr_out,   // sparse (by the first index) [type][tilesY][tilesX][(2*transform_size-1)*(2*transform_size-1)] or null
						    // combo will be added as extra pair if mcorr_comb_width > 0 and clt_corr_out has a slot for it
							// to be converted to float
							dcorr_tiles,                   // final double  [][][]      dcorr_tiles,     // [tile][pair][(2*transform_size-1)*(2*transform_size-1)] // if null - will not calculate
							// When clt_mismatch is non-zero, no far objects extraction will be attempted
							//optional, may be null
							disparity_map,                 // final double [][]         disparity_map,   // [8][tilesY][tilesX], only [6][] is needed on input or null - do not calculate
							null,                          // final double [][]         ddnd,            // data for LY. SHould be either null or [num_sensors][]
							clt_parameters.correlate_lma,  // final boolean             run_lma,         // calculate LMA, false - CM only
				  		    // define combining of all 2D correlation pairs for CM (LMA does not use them)
							clt_parameters.img_dtt.mcorr_comb_width, //final int                 mcorr_comb_width,  // combined correlation tile width (set <=0 to skip combined correlations)
							clt_parameters.img_dtt.mcorr_comb_height,//final int                 mcorr_comb_height, // combined correlation tile full height
							clt_parameters.img_dtt.mcorr_comb_offset,//final int                 mcorr_comb_offset, // combined correlation tile height offset: 0 - centered (-height/2 to height/2), height/2 - only positive (0 to height)
							clt_parameters.img_dtt.mcorr_comb_disp,	 //final double              mcorr_comb_disp,   // Combined tile per-pixel disparity for baseline == side of a square
							clt_parameters.clt_window,     // final int                 window_type,     // GPU: will not be used
							clt_parameters.tileX,          // final int                 debug_tileX,
							clt_parameters.tileY,          // final int                 debug_tileY,
							threadsMax,                    // final int                 threadsMax,      // maximal number of threads to launch
							"reference",                   // final String              debug_suffix,
							debug_level + 2); // -1 );              // final int                 globalDebugLevel)
					ImageDtt.convertFcltCorr(
							dcorr_tiles, // double [][][] dcorr_tiles,// [tile][sparse, correlation pair][(2*transform_size-1)*(2*transform_size-1)] // if null - will not calculate
							fclt_corr);  // float  [][][] fclt_corr) //  new float [tilesX * tilesY][][] or null
				}
				
			} else { // CPU version
				//Correlation2d correlation2d = 
				image_dtt.getCorrelation2d();
				double [][][][]   dcorr_td = new double[tp_tasks.length][][][]; // [tile][pair][4][64] sparse by pair transform domain representation of corr pairs
				image_dtt.quadCorrTD(
						scenes[nscene].getImageData(),                      // final double [][][]       image_data,      // first index - number of image in a quad
						scenes[nscene].getErsCorrection().getSensorWH()[0], // final int                 width,
						tp_tasks,                                           // final TpTask []           tp_tasks,
						clt_parameters.img_dtt,                             // final ImageDttParameters  imgdtt_params,    // Now just extra correlation parameters, later will include, most others
						dcorr_td,                                           // final double [][][][][]   dcorr_td,        // [pair][tilesY][tilesX][4][64] sparse transform domain representation of corr pairs
						// no combo here - rotate, combine in pixel domain after interframe
						scenes[nscene].getCltKernels(),                     // final double [][][][][][] clt_kernels,     // [channel_in_quad][color][tileY][tileX][band][pixel] , size should match image (have 1 tile around)
//						clt_parameters.kernel_step,                         // final int                 kernel_step,
						clt_parameters.clt_window,                          // final int                 window_type,
						clt_parameters.corr_red,                            // final double              corr_red,
						clt_parameters.corr_blue,                           // final double              corr_blue,
						mcorr_sel,                                          // final int                 mcorr_sel,    // Which pairs to correlate // +1 - all, +2 - dia, +4 - sq, +8 - neibs, +16 - hor + 32 - vert
						clt_parameters.tileX,                               // final int                 debug_tileX,
						clt_parameters.tileY,                               // final int                 debug_tileY,
						threadsMax,                                         // final int                 threadsMax,       // maximal number of threads to launch
						debug_level);                                       // final int                 globalDebugLevel)
				
				accumulateCorrelations(
						tp_tasks,      // final TpTask []         tp_tasks,
						num_acc,       // final int [][][]        num_acc,     // number of accumulated tiles [tilesY][tilesX][pair]
						dcorr_td,      // final double [][][][][] dcorr_td,    // [tile][pair][4][64] sparse transform domain representation of corr pairs 
						dcorr_td_acc); // final double [][][][][] dcorr_td_acc // [pair][tilesY][tilesX][4][64] sparse transform domain representation of corr pairs 

				if ((nscene == indx_ref) && show_reference_correlations) { // prepare 2d correlations for visualization, double/CPU mode
					double [][][] dcorr_tiles = (fclt_corr != null)? (new double [tp_tasks_ref.length][][]):null;
					image_dtt.clt_process_tl_correlations( // convert to pixel domain and process correlations already prepared in fcorr_td and/or fcorr_combo_td
							clt_parameters.img_dtt,		   // final ImageDttParameters  imgdtt_params,   // Now just extra correlation parameters, later will include, most others
							tp_tasks_ref,                  // final TpTask []           tp_tasks,        // data from the reference frame - will be applied to LMW for the integrated correlations
							// only listed tiles will be processed
							ref_scene.getErsCorrection().getRXY(false), // final double [][]         rXY,             // from geometryCorrection
							tilesX,                        // final int                 tilesX,          // tp_tasks may lack maximal tileX, tileY  
							tilesY,                        // final int                 tilesY,
							// no fcorr_combo_td here both arrays should have same non-null tiles
							dcorr_td,                      // final double [][][][]     dcorr_td,        // [tile][pair][4][64] sparse by pair transform domain representation of corr pairs
							null,                          // final double []           dcorr_weight,    // [tile] weighted number of tiles averaged (divide squared fat zero by this)
							// next both can be nulls
							null,                          // final double [][][][]     clt_corr_out,   // sparse (by the first index) [type][tilesY][tilesX][(2*transform_size-1)*(2*transform_size-1)] or null
							// to be converted to float
							dcorr_tiles,                   // final double  [][][]      dcorr_tiles,     // [tile][pair][(2*transform_size-1)*(2*transform_size-1)] // if null - will not calculate
							// When clt_mismatch is non-zero, no far objects extraction will be attempted
							//optional, may be null
							disparity_map,                 // final double [][]         disparity_map,   // [8][tilesY][tilesX], only [6][] is needed on input or null - do not calculate
							null,                          // final double [][]         ddnd,            // data for LY. SHould be either null or [num_sensors][]
							clt_parameters.correlate_lma,  // final boolean             run_lma,         // calculate LMA, false - CM only
							// last 2 - contrast, avg/ "geometric average)
							clt_parameters.getGpuFatZero(ref_scene.isMonochrome()),   // clt_parameters.getGpuFatZero(ref_scene.isMonochrome()), // final double              afat_zero2,      // gpu_fat_zero ==30? clt_parameters.getGpuFatZero(is_mono); absolute fat zero, same units as components squared values
							clt_parameters.gpu_sigma_m,    // final double              corr_sigma,      //
							// define combining of all 2D correlation pairs for CM (LMA does not use them)

							clt_parameters.img_dtt.mcorr_comb_width, // final int                 mcorr_comb_width,  // combined correlation tile width
							clt_parameters.img_dtt.mcorr_comb_height,// final int                 mcorr_comb_height, // combined correlation tile full height
							clt_parameters.img_dtt.mcorr_comb_offset,// final int                 mcorr_comb_offset, // combined correlation tile height offset: 0 - centered (-height/2 to height/2), height/2 - only positive (0 to height)
							clt_parameters.img_dtt.mcorr_comb_disp,	 // final double              mcorr_comb_disp,   // Combined tile per-pixel disparity for baseline == side of a square

							clt_parameters.clt_window,     // final int                 window_type,     // GPU: will not be used
							clt_parameters.tileX,          // final int                 debug_tileX,
							clt_parameters.tileY,          // final int                 debug_tileY,
							threadsMax,                    // final int                 threadsMax,      // maximal number of threads to launch
							null,                          // final String              debug_suffix,
							debug_level + 2); // -1 );              // final int                 globalDebugLevel)
					ImageDtt.convertFcltCorr(
							dcorr_tiles, // double [][][] dcorr_tiles,// [tile][sparse, correlation pair][(2*transform_size-1)*(2*transform_size-1)] // if null - will not calculate
							fclt_corr);  // float  [][][] fclt_corr) //  new float [tilesX * tilesY][][] or null
				}
			} // GPU-version, CPU-version

			if ((nscene == indx_ref) && show_reference_correlations) { // visualize prepare ref_scene correlation data
				float [][] dbg_corr_rslt_partial = ImageDtt.corr_partial_dbg( // not used in lwir
						fclt_corr, // final float  [][][]     fcorr_data,       // [tile][pair][(2*transform_size-1)*(2*transform_size-1)] // if null - will not calculate
						tp_tasks_ref, // final TpTask []         tp_tasks,        //
						tilesX,    //final int               tilesX,
						tilesY,    //final int               tilesX,
						2*image_dtt.transform_size - 1,	// final int               corr_size,
						1000, // will be limited by available layersfinal int               layers0,
						clt_parameters.corr_border_contrast, // final double            border_contrast,
						threadsMax, // final int               threadsMax,     // maximal number of threads to launch
						debug_level); // final int               globalDebugLevel)
				
				String [] titles = new String [dbg_corr_rslt_partial.length]; // dcorr_tiles[0].length];
				int ind_length = image_dtt.getCorrelation2d().getCorrTitles().length;
				
				System.arraycopy(image_dtt.getCorrelation2d().getCorrTitles(), 0, titles, 0, ind_length);
				for (int i = ind_length; i < titles.length; i++) {
					titles[i] = "combo-"+(i - ind_length);
				}
				
				// titles.length = 15, corr_rslt_partial.length=16!
				(new ShowDoubleFloatArrays()).showArrays( // out of boundary 15
						dbg_corr_rslt_partial,
						tilesX*(2*image_dtt.transform_size),
						tilesY*(2*image_dtt.transform_size),
						true,
						ref_scene.getImageName()+"-CORR-REFSCENE-"+nrefine,
						titles); // image_dtt.getCorrelation2d().getCorrTitles()); //CORR_TITLES);

			}



		} // for (int nscene = 0; nscene < num_scenes; nscene++) {
///		Runtime.getRuntime().gc();
///		System.out.println("--- Free memory="+Runtime.getRuntime().freeMemory()+" (of "+Runtime.getRuntime().totalMemory()+")");
		// Normalize accumulated correlations
		if (ref_scene.hasGPU()) {
			 accumulateCorrelationsAcOnly(
					 num_acc,       // final int [][][]     num_acc,          // number of accumulated tiles [tilesY][tilesX][pair]
					 fcorr_td_acc); // final float [][][][] fcorr_td_acc      // [tilesY][tilesX][pair][256] sparse transform domain representation of corr pairs 
		} else {
			accumulateCorrelations(
					num_acc,       // final int [][][]        num_acc,     // number of accumulated tiles [tilesY][tilesX][pair]
					dcorr_td_acc); // final double [][][][][] dcorr_td_acc // [pair][tilesY][tilesX][4][64] sparse transform domain representation of corr pairs 
		}
		
		if (ref_scene.hasGPU()) {
			double [][][]     dcorr_tiles = (fclt_corr != null)? (new double [tp_tasks_ref.length][][]):null;
			image_dtt.clt_process_tl_correlations( // convert to pixel domain and process correlations already prepared in fcorr_td and/or fcorr_combo_td
					clt_parameters.img_dtt,		   // final ImageDttParameters  imgdtt_params,   // Now just extra correlation parameters, later will include, most others
					fcorr_td_acc,		 	     // final float  [][][][]     fcorr_td,        // [tilesY][tilesX][pair][4*64] transform domain representation of all selected corr pairs
					num_acc,                       // float [][][]                num_acc,         // number of accumulated tiles [tilesY][tilesX][pair] (or null)       
					null, // dcorr_weight,                  // double []                 dcorr_weight,    // alternative to num_acc, compatible with CPU processing (only one non-zero enough)
					clt_parameters.gpu_corr_scale, //  final double              gpu_corr_scale,  //  0.75; // reduce GPU-generated correlation values
					clt_parameters.getGpuFatZero(ref_scene.isMonochrome()),   // final double     gpu_fat_zero,    // clt_parameters.getGpuFatZero(is_mono);absolute == 30.0
					image_dtt.transform_size - 1,  // final int                 gpu_corr_rad,    // = transform_size - 1 ?
			        // The tp_tasks data should be decoded from GPU to get coordinates
					tp_tasks_ref,                  // final TpTask []           tp_tasks,        // data from the reference frame - will be applied to LMW for the integrated correlations
					ref_scene.getErsCorrection().getRXY(false), // final double [][]         rXY,             // from geometryCorrection
					// next both can be nulls
					null,                          // final double [][][][]     clt_corr_out,   // sparse (by the first index) [type][tilesY][tilesX][(2*transform_size-1)*(2*transform_size-1)] or null
				    // combo will be added as extra pair if mcorr_comb_width > 0 and clt_corr_out has a slot for it
					// to be converted to float
					dcorr_tiles,                   // final double  [][][]      dcorr_tiles,     // [tile][pair][(2*transform_size-1)*(2*transform_size-1)] // if null - will not calculate
					// When clt_mismatch is non-zero, no far objects extraction will be attempted
					//optional, may be null
					disparity_map,                 // final double [][]         disparity_map,   // [8][tilesY][tilesX], only [6][] is needed on input or null - do not calculate
					null,                          // final double [][]         ddnd,            // data for LY. SHould be either null or [num_sensors][]
					clt_parameters.correlate_lma,  // final boolean             run_lma,         // calculate LMA, false - CM only
		  		    // define combining of all 2D correlation pairs for CM (LMA does not use them)
					clt_parameters.img_dtt.mcorr_comb_width, //final int                 mcorr_comb_width,  // combined correlation tile width (set <=0 to skip combined correlations)
					clt_parameters.img_dtt.mcorr_comb_height,//final int                 mcorr_comb_height, // combined correlation tile full height
					clt_parameters.img_dtt.mcorr_comb_offset,//final int                 mcorr_comb_offset, // combined correlation tile height offset: 0 - centered (-height/2 to height/2), height/2 - only positive (0 to height)
					clt_parameters.img_dtt.mcorr_comb_disp,	 //final double              mcorr_comb_disp,   // Combined tile per-pixel disparity for baseline == side of a square
					clt_parameters.clt_window,     // final int                 window_type,     // GPU: will not be used
					clt_parameters.tileX,          // final int                 debug_tileX,
					clt_parameters.tileY,          // final int                 debug_tileY,
					threadsMax,                    // final int                 threadsMax,      // maximal number of threads to launch
					"accumulated",                 // final String              debug_suffix,
					debug_level + 2); // -1 );     // final int                 globalDebugLevel)
			ImageDtt.convertFcltCorr(
					dcorr_tiles, // double [][][] dcorr_tiles,// [tile][sparse, correlation pair][(2*transform_size-1)*(2*transform_size-1)] // if null - will not calculate
					fclt_corr);  // float  [][][] fclt_corr) //  new float [tilesX * tilesY][][] or null

		} else {
			double [][][][]   dcorr_td = new double[tp_tasks_ref.length][][][]; // [tile][pair][4][64] sparse transform domain representation of corr pairs
			double []         dcorr_weight = new double[tp_tasks_ref.length];
			for (int iTile = 0; iTile < dcorr_td.length; iTile++) {
				TpTask task = tp_tasks_ref[iTile];
				int tileY = task.getTileY(); // tilesX;  
				int tileX = task.getTileX(); // nTile % tilesX;
				dcorr_td[iTile] = new double [num_pairs][][];
				for (int npair = 0; npair < num_pairs; npair++) if (dcorr_td_acc[npair] != null){
					dcorr_td[iTile][npair] = dcorr_td_acc[npair][tileY][tileX];
					dcorr_weight[iTile] = num_acc[tileY][tileX][npair]; // number of accumulated tiles [tilesY][tilesX][pair]
				}
			}

			double [][][]     dcorr_tiles = (fclt_corr != null)? (new double [tp_tasks_ref.length][][]):null;
			image_dtt.clt_process_tl_correlations( // convert to pixel domain and process correlations already prepared in fcorr_td and/or fcorr_combo_td
					clt_parameters.img_dtt,		   // final ImageDttParameters  imgdtt_params,   // Now just extra correlation parameters, later will include, most others
					tp_tasks_ref,                  // final TpTask []           tp_tasks,        // data from the reference frame - will be applied to LMW for the integrated correlations
					// only listed tiles will be processed
					ref_scene.getErsCorrection().getRXY(false), // final double [][]         rXY,             // from geometryCorrection
					tilesX,                        // final int                 tilesX,          // tp_tasks may lack maximal tileX, tileY  
					tilesY,                        // final int                 tilesY,
					// no fcorr_combo_td here both arrays should have same non-null tiles
					dcorr_td,                      // final double [][][][]     dcorr_td,        // [tile][pair][4][64] sparse by pair transform domain representation of corr pairs
					dcorr_weight,                   // final double []           dcorr_weight,    // [tile] weighted number of tiles averaged (divide squared fat zero by this)

					// next both can be nulls
					null,                          // final double [][][][]     clt_corr_out,   // sparse (by the first index) [type][tilesY][tilesX][(2*transform_size-1)*(2*transform_size-1)] or null
					// to be converted to float
					dcorr_tiles,                   // final double  [][][]      dcorr_tiles,     // [tile][pair][(2*transform_size-1)*(2*transform_size-1)] // if null - will not calculate
					// When clt_mismatch is non-zero, no far objects extraction will be attempted
					//optional, may be null
					disparity_map,                 // final double [][]         disparity_map,   // [8][tilesY][tilesX], only [6][] is needed on input or null - do not calculate
					null,                          // final double [][]         ddnd,            // data for LY. SHould be either null or [num_sensors][]
					clt_parameters.correlate_lma, // true,                          // final boolean             run_lma,         // calculate LMA, false - CM only
					// last 2 - contrast, avg/ "geometric average)
					clt_parameters.getGpuFatZero(ref_scene.isMonochrome()),   // clt_parameters.getGpuFatZero(ref_scene.isMonochrome()), // final double              afat_zero2,      // gpu_fat_zero ==30? clt_parameters.getGpuFatZero(is_mono); absolute fat zero, same units as components squared values
					clt_parameters.gpu_sigma_m,    // final double              corr_sigma,      //
					// define combining of all 2D correlation pairs for CM (LMA does not use them)

					clt_parameters.img_dtt.mcorr_comb_width, // final int                 mcorr_comb_width,  // combined correlation tile width
					clt_parameters.img_dtt.mcorr_comb_height,// final int                 mcorr_comb_height, // combined correlation tile full height
					clt_parameters.img_dtt.mcorr_comb_offset,// final int                 mcorr_comb_offset, // combined correlation tile height offset: 0 - centered (-height/2 to height/2), height/2 - only positive (0 to height)
					clt_parameters.img_dtt.mcorr_comb_disp,	 // final double              mcorr_comb_disp,   // Combined tile per-pixel disparity for baseline == side of a square

					clt_parameters.clt_window,     // final int                 window_type,     // GPU: will not be used
					clt_parameters.tileX,          // final int                 debug_tileX,
					clt_parameters.tileY,          // final int                 debug_tileY,
					threadsMax,                    // final int                 threadsMax,      // maximal number of threads to launch
					null,                          // final String              debug_suffix,
					debug_level + (show_reference_correlations? 2 : -1));              // final int                 globalDebugLevel)
			ImageDttCPU.convertFcltCorr(
					dcorr_tiles, // double [][][] dcorr_tiles,// [tile][sparse, correlation pair][(2*transform_size-1)*(2*transform_size-1)] // if null - will not calculate
					fclt_corr);  // float  [][][] fclt_corr) //  new float [tilesX * tilesY][][] or null
			
		}

///		Runtime.getRuntime().gc();
///		System.out.println("--- Free memory="+Runtime.getRuntime().freeMemory()+" (of "+Runtime.getRuntime().totalMemory()+")");
		if (show_accumulated_correlations || (accum_2d_corr != null)){ // -1
			float [][] accum_2d_img = ImageDtt.corr_partial_dbg( // not used in lwir
					fclt_corr, // final float  [][][]     fcorr_data,       // [tile][pair][(2*transform_size-1)*(2*transform_size-1)] // if null - will not calculate
					tp_tasks_ref, // final TpTask []         tp_tasks,        //
					tilesX,    //final int               tilesX,
					tilesY,    //final int               tilesX,
					2*image_dtt.transform_size - 1,	// final int               corr_size,
					1000, // will be limited by available layersfinal int               layers0,
					clt_parameters.corr_border_contrast, // final double            border_contrast,
					threadsMax, // final int               threadsMax,     // maximal number of threads to launch
					debug_level); // final int               globalDebugLevel)
			String [] titles = new String [accum_2d_img.length]; // dcorr_tiles[0].length];
			int ind_length = image_dtt.getCorrelation2d().getCorrTitles().length;
			
			System.arraycopy(image_dtt.getCorrelation2d().getCorrTitles(), 0, titles, 0, ind_length);
			for (int i = ind_length; i < titles.length; i++) {
				titles[i] = "combo-"+(i - ind_length);
			}
			if ((accum_2d_corr != null)) {
				accum_2d_corr[0] = accum_2d_img;
			}
			if (show_accumulated_correlations) {
				(new ShowDoubleFloatArrays()).showArrays( // out of boundary 15
						accum_2d_img,
						tilesX*(2*image_dtt.transform_size),
						tilesY*(2*image_dtt.transform_size),
						true,
						ref_scene.getImageName()+"-CORR-ACCUM"+num_scenes+"-"+nrefine,
						titles);      // image_dtt.getCorrelation2d().getCorrTitles()); //CORR_TITLES);
			}
		}
		return disparity_map; // disparity_map
	}
	
	/**
	 * Generate (randomly) offset 2D correlations for DNN training
	 * @param clt_parameters
	 * @param ref_scene
	 * @param disparity_ref_in
	 * @param disparity_offsets
	 * @param margin
	 * @param add_combo
	 * @param debug_level
	 * @return
	 */
	public float[][] generateOffset2DCorrelations(
			final CLTParameters  clt_parameters,
			final QuadCLT        ref_scene,
			final double []      disparity_ref_in,  // disparity in the reference view tiles (Double.NaN - invalid)
			final double []      disparity_offsets,
			final int            margin,
			final boolean        add_combo,
			final int            debug_level)
	{
		final double [] disparity_ref = disparity_ref_in.clone();
		///		final ErsCorrection ers_reference = ref_scene.getErsCorrection();
		final int tilesX = ref_scene.getTileProcessor().getTilesX();
		final int tilesY = ref_scene.getTileProcessor().getTilesY();
		int mcorr_sel = Correlation2d.corrSelEncodeAll(0); // all sensors
		
		final float  [][][]  fclt_corr =  new float [tilesX * tilesY][][];
		ImageDtt image_dtt;
		image_dtt = new ImageDtt(
				numSens,
				clt_parameters.transform_size,
				clt_parameters.img_dtt,
				ref_scene.isAux(),
				ref_scene.isMonochrome(),
				ref_scene.isLwir(),
				clt_parameters.getScaleStrength(ref_scene.isAux()),
				ref_scene.getGPU());

		image_dtt.getCorrelation2d(); // initiate image_dtt.correlation2d, needed if disparity_map != null  

//		double[][] disparity_map = new double [image_dtt.getDisparityTitles().length][];
		final double disparity_corr = 0.0; // (z_correction == 0) ? 0.0 : geometryCorrection.getDisparityFromZ(1.0/z_correction);
		TpTask[] tp_tasks_ref = null;
//		String ts = ref_scene.getImageName();
		double [][] scene_pXpYD;
		// transform to self - maybe use a method that sets central points
		scene_pXpYD = transformToScenePxPyD(
				null, // final Rectangle [] extra_woi,    // show larger than sensor WOI (or null)
				disparity_ref,      // final double []   disparity_ref, // invalid tiles - NaN in disparity (maybe it should not be masked by margins?)
				ZERO3,              // final double []   scene_xyz, // camera center in world coordinates
				ZERO3,              // final double []   scene_atr, // camera orientation relative to world frame
				ref_scene,          // final QuadCLT     scene_QuadClt,
				ref_scene);         // final QuadCLT     reference_QuadClt)
	for (int i = 0; i < scene_pXpYD.length; i++) {
		if ((scene_pXpYD[i] != null) && !Double.isNaN(scene_pXpYD[i][2])) {
			scene_pXpYD[i][2] += disparity_offsets[i];
		}
	}

		ref_scene.saveQuadClt(); // to re-load new set of Bayer images to the GPU (do nothing for CPU)
		final double gpu_sigma_corr =     clt_parameters.getGpuCorrSigma(ref_scene.isMonochrome());
		final double gpu_sigma_rb_corr =  ref_scene.isMonochrome()? 1.0 : clt_parameters.gpu_sigma_rb_corr;
		final double gpu_sigma_log_corr = clt_parameters.getGpuCorrLoGSigma(ref_scene.isMonochrome());

		TpTask[] tp_tasks =  GpuQuad.setInterTasks(
				ref_scene.getNumSensors(),
				ref_scene.getErsCorrection().getSensorWH()[0],
				!ref_scene.hasGPU(), // final boolean             calcPortsCoordinatesAndDerivatives, // GPU can calculate them centreXY
				scene_pXpYD,              // final double [][]         pXpYD, // per-tile array of pX,pY,disparity triplets (or nulls)
				null,              // final boolean []          selection, // may be null, if not null do not  process unselected tiles
				ref_scene.getErsCorrection(), // final GeometryCorrection  geometryCorrection,
				disparity_corr,     // final double              disparity_corr,
				margin,             // final int                 margin,      // do not use tiles if their centers are closer to the edges
				null,               // final boolean []          valid_tiles,            
				threadsMax);        // final int                 threadsMax)  // maximal number of threads to launch
		tp_tasks_ref = tp_tasks; // will use coordinates data for LMA ? disp_dist
		if (ref_scene.hasGPU()) {
			float  [][][][]     fcorr_td =       new float[tilesY][tilesX][][];
			image_dtt.quadCorrTD(
					clt_parameters.img_dtt,            // final ImageDttParameters imgdtt_params,    // Now just extra correlation parameters, later will include, most others
					tp_tasks, // *** will be updated inside from GPU-calculated geometry
					fcorr_td, // fcorrs_td[nscene],                 // [tilesY][tilesX][pair][4*64] transform domain representation of 6 corr pairs
//					ref_scene.getErsCorrection(), //
					clt_parameters.gpu_sigma_r,        // 0.9, 1.1
					clt_parameters.gpu_sigma_b,        // 0.9, 1.1
					clt_parameters.gpu_sigma_g,        // 0.6, 0.7
					clt_parameters.gpu_sigma_m,        //  =       0.4; // 0.7;
					gpu_sigma_rb_corr,                 // final double              gpu_sigma_rb_corr, //  = 0.5; // apply LPF after accumulating R and B correlation before G, monochrome ? 1.0 : gpu_sigma_rb_corr;
					gpu_sigma_corr,                    //  =    0.9;gpu_sigma_corr_m
					gpu_sigma_log_corr,                // final double              gpu_sigma_log_corr,   // hpf to reduce dynamic range for correlations
					clt_parameters.corr_red,           // +used
					clt_parameters.corr_blue,          // +used
					mcorr_sel,                         // final int                 mcorr_sel,    // Which pairs to correlate // +1 - all, +2 - dia, +4 - sq, +8 - neibs, +16 - hor + 32 - vert
					threadsMax,       // maximal number of threads to launch
					debug_level);
			if (image_dtt.getGPU().getGpu_debug_level() > -1) {
				System.out.println("==ooo=after image_dtt.quadCorrTD()");
			}

			//Use same for both CPU/GPU
			double [][][] dcorr_tiles = (fclt_corr != null)? (new double [tp_tasks_ref.length][][]):null;
			image_dtt.clt_process_tl_correlations( // convert to pixel domain and process correlations already prepared in fcorr_td and/or fcorr_combo_td
					clt_parameters.img_dtt,		   // final ImageDttParameters  imgdtt_params,   // Now just extra correlation parameters, later will include, most others
					fcorr_td,		 	 		   // final float  [][][][]     fcorr_td,        // [tilesY][tilesX][pair][4*64] transform domain representation of all selected corr pairs
					null, // num_acc,              // int [][][]                num_acc,         // number of accumulated tiles [tilesY][tilesX][pair] (or null)       
					null, // dcorr_weight,                  // double []                 dcorr_weight,    // alternative to num_acc, compatible with CPU processing (only one non-zero enough)
					clt_parameters.gpu_corr_scale, //  final double              gpu_corr_scale,  //  0.75; // reduce GPU-generated correlation values
					clt_parameters.getGpuFatZero(ref_scene.isMonochrome()),   // final double     gpu_fat_zero,    // clt_parameters.getGpuFatZero(is_mono);absolute == 30.0
					image_dtt.transform_size - 1,  // final int                 gpu_corr_rad,    // = transform_size - 1 ?
					// The tp_tasks data should be decoded from GPU to get coordinates
					tp_tasks_ref,                  // final TpTask []           tp_tasks,        // data from the reference frame - will be applied to LMW for the integrated correlations
					ref_scene.getErsCorrection().getRXY(false), // final double [][]         rXY,             // from geometryCorrection
					// next both can be nulls
					null,                          // final double [][][][]     clt_corr_out,   // sparse (by the first index) [type][tilesY][tilesX][(2*transform_size-1)*(2*transform_size-1)] or null
					// combo will be added as extra pair if mcorr_comb_width > 0 and clt_corr_out has a slot for it
					// to be converted to float
					dcorr_tiles,                   // final double  [][][]      dcorr_tiles,     // [tile][pair][(2*transform_size-1)*(2*transform_size-1)] // if null - will not calculate
					// When clt_mismatch is non-zero, no far objects extraction will be attempted
					//optional, may be null
					null, // disparity_map,                 // final double [][]         disparity_map,   // [8][tilesY][tilesX], only [6][] is needed on input or null - do not calculate
					null,                          // final double [][]         ddnd,            // data for LY. SHould be either null or [num_sensors][]
					clt_parameters.correlate_lma,  // final boolean             run_lma,         // calculate LMA, false - CM only
					// define combining of all 2D correlation pairs for CM (LMA does not use them)
					(add_combo ? clt_parameters.img_dtt.mcorr_comb_width : 0), //final int                 mcorr_comb_width,  // combined correlation tile width (set <=0 to skip combined correlations)
					clt_parameters.img_dtt.mcorr_comb_height,//final int                 mcorr_comb_height, // combined correlation tile full height
					clt_parameters.img_dtt.mcorr_comb_offset,//final int                 mcorr_comb_offset, // combined correlation tile height offset: 0 - centered (-height/2 to height/2), height/2 - only positive (0 to height)
					clt_parameters.img_dtt.mcorr_comb_disp,	 //final double              mcorr_comb_disp,   // Combined tile per-pixel disparity for baseline == side of a square
					clt_parameters.clt_window,     // final int                 window_type,     // GPU: will not be used
					clt_parameters.tileX,          // final int                 debug_tileX,
					clt_parameters.tileY,          // final int                 debug_tileY,
					threadsMax,                    // final int                 threadsMax,      // maximal number of threads to launch
					null,                          // final String              debug_suffix,
					debug_level + 2); // -1 );              // final int                 globalDebugLevel)
			ImageDtt.convertFcltCorr(
					dcorr_tiles, // double [][][] dcorr_tiles,// [tile][sparse, correlation pair][(2*transform_size-1)*(2*transform_size-1)] // if null - will not calculate
					fclt_corr);  // float  [][][] fclt_corr) //  new float [tilesX * tilesY][][] or null

		} else { // CPU version
			//Correlation2d correlation2d = 
			image_dtt.getCorrelation2d();
			double [][][][]   dcorr_td = new double[tp_tasks.length][][][]; // [tile][pair][4][64] sparse by pair transform domain representation of corr pairs
			image_dtt.quadCorrTD(
					ref_scene.getImageData(),                      // final double [][][]       image_data,      // first index - number of image in a quad
					ref_scene.getErsCorrection().getSensorWH()[0], // final int                 width,
					tp_tasks,                                           // final TpTask []           tp_tasks,
					clt_parameters.img_dtt,                             // final ImageDttParameters  imgdtt_params,    // Now just extra correlation parameters, later will include, most others
					dcorr_td,                                           // final double [][][][][]   dcorr_td,        // [pair][tilesY][tilesX][4][64] sparse transform domain representation of corr pairs
					// no combo here - rotate, combine in pixel domain after interframe
					ref_scene.getCltKernels(),                     // final double [][][][][][] clt_kernels,     // [channel_in_quad][color][tileY][tileX][band][pixel] , size should match image (have 1 tile around)
//					clt_parameters.kernel_step,                         // final int                 kernel_step,
					clt_parameters.clt_window,                          // final int                 window_type,
					clt_parameters.corr_red,                            // final double              corr_red,
					clt_parameters.corr_blue,                           // final double              corr_blue,
					mcorr_sel,                                          // final int                 mcorr_sel,    // Which pairs to correlate // +1 - all, +2 - dia, +4 - sq, +8 - neibs, +16 - hor + 32 - vert
					clt_parameters.tileX,                               // final int                 debug_tileX,
					clt_parameters.tileY,                               // final int                 debug_tileY,
					threadsMax,                                         // final int                 threadsMax,       // maximal number of threads to launch
					debug_level);                                       // final int                 globalDebugLevel)

			double [][][] dcorr_tiles = (fclt_corr != null)? (new double [tp_tasks_ref.length][][]):null;
			image_dtt.clt_process_tl_correlations( // convert to pixel domain and process correlations already prepared in fcorr_td and/or fcorr_combo_td
					clt_parameters.img_dtt,		   // final ImageDttParameters  imgdtt_params,   // Now just extra correlation parameters, later will include, most others
					tp_tasks_ref,                  // final TpTask []           tp_tasks,        // data from the reference frame - will be applied to LMW for the integrated correlations
					// only listed tiles will be processed
					ref_scene.getErsCorrection().getRXY(false), // final double [][]         rXY,             // from geometryCorrection
					tilesX,                        // final int                 tilesX,          // tp_tasks may lack maximal tileX, tileY  
					tilesY,                        // final int                 tilesY,
					// no fcorr_combo_td here both arrays should have same non-null tiles
					dcorr_td,                      // final double [][][][]     dcorr_td,        // [tile][pair][4][64] sparse by pair transform domain representation of corr pairs
					null,                          // final double []           dcorr_weight,    // [tile] weighted number of tiles averaged (divide squared fat zero by this)
					// next both can be nulls
					null,                          // final double [][][][]     clt_corr_out,   // sparse (by the first index) [type][tilesY][tilesX][(2*transform_size-1)*(2*transform_size-1)] or null
					// to be converted to float
					dcorr_tiles,                   // final double  [][][]      dcorr_tiles,     // [tile][pair][(2*transform_size-1)*(2*transform_size-1)] // if null - will not calculate
					// When clt_mismatch is non-zero, no far objects extraction will be attempted
					//optional, may be null
					null, // disparity_map,        // final double [][]         disparity_map,   // [8][tilesY][tilesX], only [6][] is needed on input or null - do not calculate
					null,                          // final double [][]         ddnd,            // data for LY. SHould be either null or [num_sensors][]
					clt_parameters.correlate_lma,  // final boolean             run_lma,         // calculate LMA, false - CM only
					// last 2 - contrast, avg/ "geometric average)
					clt_parameters.getGpuFatZero(ref_scene.isMonochrome()),   // clt_parameters.getGpuFatZero(ref_scene.isMonochrome()), // final double              afat_zero2,      // gpu_fat_zero ==30? clt_parameters.getGpuFatZero(is_mono); absolute fat zero, same units as components squared values
					clt_parameters.gpu_sigma_m,    // final double              corr_sigma,      //
					// define combining of all 2D correlation pairs for CM (LMA does not use them)

					0, // clt_parameters.img_dtt.mcorr_comb_width, // final int                 mcorr_comb_width,  // combined correlation tile width
					clt_parameters.img_dtt.mcorr_comb_height,// final int                 mcorr_comb_height, // combined correlation tile full height
					clt_parameters.img_dtt.mcorr_comb_offset,// final int                 mcorr_comb_offset, // combined correlation tile height offset: 0 - centered (-height/2 to height/2), height/2 - only positive (0 to height)
					clt_parameters.img_dtt.mcorr_comb_disp,	 // final double              mcorr_comb_disp,   // Combined tile per-pixel disparity for baseline == side of a square

					clt_parameters.clt_window,     // final int                 window_type,     // GPU: will not be used
					clt_parameters.tileX,          // final int                 debug_tileX,
					clt_parameters.tileY,          // final int                 debug_tileY,
					threadsMax,                    // final int                 threadsMax,      // maximal number of threads to launch
					null,                          // final String              debug_suffix,
					debug_level + 2); // -1 );              // final int                 globalDebugLevel)
			image_dtt.convertFcltCorr(
					dcorr_tiles, // double [][][] dcorr_tiles,// [tile][sparse, correlation pair][(2*transform_size-1)*(2*transform_size-1)] // if null - will not calculate
					fclt_corr);  // float  [][][] fclt_corr) //  new float [tilesX * tilesY][][] or null
		} // GPU-version, CPU-version
		
		float [][] dbg_corr_rslt_partial = ImageDtt.corr_partial_dbg( // not used in lwir
				fclt_corr, // final float  [][][]     fcorr_data,       // [tile][pair][(2*transform_size-1)*(2*transform_size-1)] // if null - will not calculate
				tp_tasks_ref, // final TpTask []         tp_tasks,        //
				tilesX,    //final int               tilesX,
				tilesY,    //final int               tilesX,
				2*image_dtt.transform_size - 1,	// final int               corr_size,
				1000, // will be limited by available layersfinal int               layers0,
				clt_parameters.corr_border_contrast, // final double            border_contrast,
				threadsMax, // final int               threadsMax,     // maximal number of threads to launch
				debug_level); // final int               globalDebugLevel)

		return dbg_corr_rslt_partial;
	}
	
	
	@Deprecated
	public void accumulateCorrelations(
			final int [][][]        num_acc,     // number of accumulated tiles [tilesY][tilesX][pair]
			final double [][][][][] dcorr_td,    // [pair][tilesY][tilesX][4][64] sparse transform domain representation of corr pairs 
			final double [][][][][] dcorr_td_acc // [pair][tilesY][tilesX][4][64] sparse transform domain representation of corr pairs 
			) {
		int tX=-1, tY=-1; 
		for (int np = 0; np < dcorr_td.length; np++) if (dcorr_td[np] != null){
			for (int ity = 0; ity < dcorr_td[np].length; ity++) if (dcorr_td[np][ity] != null){
				tY = dcorr_td[np].length;
				tX = dcorr_td[np][ity].length;
				np = dcorr_td.length;
				break;
			}
		}
		final int tilesY = tY;
		final int tilesX = tX;
		final int tiles =  tilesY * tilesX;
		for (int pair = 0; pair < dcorr_td.length; pair++) {
			if ((dcorr_td[pair] != null) && (dcorr_td_acc[pair] == null)) {
				dcorr_td_acc[pair] = new double[tilesY][tilesX][][];
			}
		}
		final Thread[] threads = ImageDtt.newThreadArray(threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				public void run() {
					for (int nTile = ai.getAndIncrement(); nTile < tiles; nTile = ai.getAndIncrement()) {
						int tileY = nTile / tilesX;  
						int tileX = nTile % tilesX;
						for (int pair = 0; pair < dcorr_td.length; pair++) if ((dcorr_td[pair] != null) && (dcorr_td[pair][tileY][tileX] != null)){
//							if (dcorr_td_acc[pair] == null) {
//								dcorr_td_acc[pair] = new double[tilesY][tilesX][][];
//							}
							if (dcorr_td_acc[pair][tileY][tileX] == null) {
								dcorr_td_acc[pair][tileY][tileX] = new double [dcorr_td[pair][tileY][tileX].length][]; // 4
								for (int q = 0; q < dcorr_td_acc[pair][tileY][tileX].length; q++) {
									dcorr_td_acc[pair][tileY][tileX][q] = dcorr_td[pair][tileY][tileX][q].clone();
								}
							} else {
								for (int q = 0; q < dcorr_td_acc[pair][tileY][tileX].length; q++) {
									for (int i = 0; i < dcorr_td_acc[pair][tileY][tileX][q].length; i++) {
										dcorr_td_acc[pair][tileY][tileX][q][i] += dcorr_td[pair][tileY][tileX][q][i];
									}
								}
							}
							num_acc[tileY][tileX][pair]++;
						}
					}
				}
			};
		}		      
		ImageDtt.startAndJoin(threads);
	}

	public void accumulateCorrelations(
			final TpTask []         tp_tasks,
			final float [][][]      num_acc,     // number of accumulated tiles [tilesY][tilesX][pair] 
			final double [][][][]   dcorr_td,    // [tile][pair][4][64] sparse transform domain representation of corr pairs 
			final double [][][][][] dcorr_td_acc // [pair][tilesY][tilesX][4][64] sparse transform domain representation of corr pairs 
			) {
		accumulateCorrelations(
				tp_tasks,
				num_acc,     // number of accumulated tiles [tilesY][tilesX][pair] 
				dcorr_td,    // [tile][pair][4][64] sparse transform domain representation of corr pairs 
				dcorr_td_acc, // [pair][tilesY][tilesX][4][64] sparse transform domain representation of corr pairs
				threadsMax);
	}
	
	
	public static void accumulateCorrelations(
			final TpTask []         tp_tasks,
			final float [][][]      num_acc,     // number of accumulated tiles [tilesY][tilesX][pair] 
			final double [][][][]   dcorr_td,    // [tile][pair][4][64] sparse transform domain representation of corr pairs 
			final double [][][][][] dcorr_td_acc, // [pair][tilesY][tilesX][4][64] sparse transform domain representation of corr pairs
			final int               threadsMax
			) {
		final int tilesY = num_acc.length;
		final int tilesX = num_acc[0].length;
		final int num_pairs = num_acc[0][0].length;
		final AtomicBoolean [] acorrs = new AtomicBoolean[num_pairs];
		for (int i = 0; i < acorrs.length; i++) {
			acorrs[i] = new AtomicBoolean();
		}
		final Thread[] threads = ImageDtt.newThreadArray(threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				public void run() {
					for (int iTile = ai.getAndIncrement(); iTile < tp_tasks.length; iTile = ai.getAndIncrement()) if (dcorr_td[iTile] != null) {
						for (int pair = 0; pair < num_pairs; pair++) if (dcorr_td[iTile][pair] != null){
							acorrs[pair].set(true);
						}
					}
				}
			};
		}		      
		ImageDtt.startAndJoin(threads);
		for (int pair = 0; pair < acorrs.length; pair++) if ((dcorr_td_acc[pair] == null) && acorrs[pair].get()) {
			dcorr_td_acc[pair] = new double[tilesY][tilesX][][];
		}
		ai.set(0);
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				public void run() {
					for (int iTile = ai.getAndIncrement(); iTile < tp_tasks.length; iTile = ai.getAndIncrement()) if ((tp_tasks[iTile] != null) && (tp_tasks[iTile].getTask() != 0)) {
						int tileY = tp_tasks[iTile].getTileY(); // tilesX;  
						int tileX = tp_tasks[iTile].getTileX(); // nTile % tilesX;
						for (int pair = 0; pair < num_pairs; pair++) if ((dcorr_td[iTile] != null) && (dcorr_td[iTile][pair] != null)){
							if (dcorr_td_acc[pair][tileY][tileX] == null) {
								dcorr_td_acc[pair][tileY][tileX] = new double [dcorr_td[iTile][pair].length][]; // 4
								for (int q = 0; q < dcorr_td_acc[pair][tileY][tileX].length; q++) {
									dcorr_td_acc[pair][tileY][tileX][q] = dcorr_td[iTile][pair][q].clone();
								}
							} else {
								for (int q = 0; q < dcorr_td_acc[pair][tileY][tileX].length; q++) {
									for (int i = 0; i < dcorr_td_acc[pair][tileY][tileX][q].length; i++) {
										dcorr_td_acc[pair][tileY][tileX][q][i] += dcorr_td[iTile][pair][q][i];
									}
								}
							}
							num_acc[tileY][tileX][pair]++;
						}
					}
				}
			};
		}		      
		ImageDtt.startAndJoin(threads);
	}
	
	
// GPU (float) version
	public void accumulateCorrelations(
			final float [][][]   num_acc,          // number of accumulated tiles [tilesY][tilesX][pair]
			final float [][][][] fcorr_td,         // [tilesY][tilesX][pair][256] sparse transform domain representation of corr pairs 
			final float [][][][] fcorr_td_acc      // [tilesY][tilesX][pair][256] sparse transform domain representation of corr pairs
			) {
		accumulateCorrelations(
				num_acc,          // number of accumulated tiles [tilesY][tilesX][pair]
				fcorr_td,         // [tilesY][tilesX][pair][256] sparse transform domain representation of corr pairs 
				fcorr_td_acc,      // [tilesY][tilesX][pair][256] sparse transform domain representation of corr pairs
				threadsMax);
	}

	
	public static void accumulateCorrelations(
			final float [][][]   num_acc,          // number of accumulated tiles [tilesY][tilesX][pair]
			final float [][][][] fcorr_td,         // [tilesY][tilesX][pair][256] sparse transform domain representation of corr pairs 
			final float [][][][] fcorr_td_acc,      // [tilesY][tilesX][pair][256] sparse transform domain representation of corr pairs
			final int threadsMax
			) {
		int tX=-1; 
		for (int ity = 0; ity < fcorr_td.length; ity++) if (fcorr_td[ity] != null){
			tX = fcorr_td[ity].length;
			break;
		}
		final int tilesY = fcorr_td.length;
		final int tilesX = tX;
		final int tiles =  tilesY * tilesX;

		final Thread[] threads = ImageDtt.newThreadArray(threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				public void run() {
					for (int nTile = ai.getAndIncrement(); nTile < tiles; nTile = ai.getAndIncrement()) {
						int tileY = nTile / tilesX;  
						int tileX = nTile % tilesX;
						if (fcorr_td [tileY][tileX] != null) {
							if (fcorr_td_acc[tileY][tileX] == null) {
								fcorr_td_acc[tileY][tileX] = new float [fcorr_td [tileY][tileX].length][];
							}
							for (int pair = 0; pair < fcorr_td [tileY][tileX].length; pair++) if (fcorr_td [tileY][tileX][pair] != null){
								if (fcorr_td_acc[tileY][tileX][pair] == null) {
									fcorr_td_acc[tileY][tileX][pair] = fcorr_td[tileY][tileX][pair].clone();
								} else {
									for (int i = 0; i < fcorr_td [tileY][tileX][pair].length; i++) {
										fcorr_td_acc[tileY][tileX][pair][i] += fcorr_td[tileY][tileX][pair][i];
									}
								}
								num_acc[tileY][tileX][pair] += 1;
							}
						}
					}
				}
			};
		}		      
		ImageDtt.startAndJoin(threads);
	}

	public void accumulateCorrelations(  // normalize
			final float  [][][]     num_acc,      // number of accumulated tiles [tilesY][tilesX][pair]
			final double [][][][][] dcorr_td_acc) { // [pair][tilesY][tilesX][4][64] sparse transform domain representation of corr pairs
		accumulateCorrelations(  // normalize
				num_acc,     // number of accumulated tiles [tilesY][tilesX][pair]
				dcorr_td_acc, // [pair][tilesY][tilesX][4][64] sparse transform domain representation of corr pairs
				threadsMax);
	}
	
	public static void accumulateCorrelations(  // normalize
			final float  [][][]     num_acc,     // number of accumulated tiles [tilesY][tilesX][pair]
			final double [][][][][] dcorr_td_acc, // [pair][tilesY][tilesX][4][64] sparse transform domain representation of corr pairs
			final int               threadsMax
			) {
		int tX=-1, tY=-1; 
		for (int np = 0; np < dcorr_td_acc.length; np++) if (dcorr_td_acc[np] != null){
			for (int ity = 0; ity < dcorr_td_acc[np].length; ity++) if (dcorr_td_acc[np][ity] != null){
				tY = dcorr_td_acc[np].length;
				tX = dcorr_td_acc[np][ity].length;
				np = dcorr_td_acc.length;
				break;
			}
		}
		final int tilesY = tY;
		final int tilesX = tX;
		final int tiles =  tilesY * tilesX;

		final Thread[] threads = ImageDtt.newThreadArray(threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				public void run() {
					for (int nTile = ai.getAndIncrement(); nTile < tiles; nTile = ai.getAndIncrement()) {
						int tileY = nTile / tilesX;  
						int tileX = nTile % tilesX;
						for (int pair = 0; pair < dcorr_td_acc.length; pair++) if (num_acc [tileY][tileX][pair] > 0){
							double k = 1.0 / num_acc [tileY][tileX][pair];
							for (int q = 0; q < dcorr_td_acc[pair][tileY][tileX].length; q++) {
								for (int i = 0; i < dcorr_td_acc[pair][tileY][tileX][q].length; i++) {
									dcorr_td_acc[pair][tileY][tileX][q][i] *= k; 
								}
							}							
						}
					}
				}
			};
		}		      
		ImageDtt.startAndJoin(threads);
	}

	public void accumulateCorrelationsAcOnly( // normalize
			final float [][][]   num_acc,          // number of accumulated tiles [tilesY][tilesX][pair]
			final float [][][][] fcorr_td_acc      // [tilesY][tilesX][pair][256] sparse transform domain representation of corr pairs 
			) {
		accumulateCorrelationsAcOnly( // normalize
				num_acc,          // number of accumulated tiles [tilesY][tilesX][pair]
				fcorr_td_acc,     // [tilesY][tilesX][pair][256] sparse transform domain representation of corr pairs
				threadsMax);
	}

	public static void accumulateCorrelationsAcOnly( // normalize
			final float [][][]   num_acc,          // number of accumulated tiles [tilesY][tilesX][pair]
			final float [][][][] fcorr_td_acc,      // [tilesY][tilesX][pair][256] sparse transform domain representation of corr pairs
			final int            threadsMax
			) {
		int tX=-1; 
		for (int ity = 0; ity < fcorr_td_acc.length; ity++) if (fcorr_td_acc[ity] != null){
			tX = fcorr_td_acc[ity].length;
			break;
		}
		final int tilesY = fcorr_td_acc.length;
		final int tilesX = tX;
		final int tiles =  tilesY * tilesX;

		final Thread[] threads = ImageDtt.newThreadArray(threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				public void run() {
					for (int nTile = ai.getAndIncrement(); nTile < tiles; nTile = ai.getAndIncrement()) {
						int tileY = nTile / tilesX;  
						int tileX = nTile % tilesX;
						if (fcorr_td_acc [tileY][tileX] != null) {
							for (int pair = 0; pair < fcorr_td_acc [tileY][tileX].length; pair++) if (num_acc [tileY][tileX][pair] >0){
								double k = 1.0 / num_acc [tileY][tileX][pair];
								for (int i = 0; i < fcorr_td_acc[tileY][tileX][pair].length; i++) {
									fcorr_td_acc[tileY][tileX][pair][i] *= k; 
								}
							}
						}
					}
				}
			};
		}		      
		ImageDtt.startAndJoin(threads);
	}

	/**
	 * Merge FD correlations and weights for tile clusters, normalize, modify fcorr_td_acc and
	 * num_acc (to be used for fat zero). Put merged data to top left tile of each cluster, remove
	 * other ones.
	 * @param clust_size Cluster size (side of a square), rightmost and bottom may be partial (used ceil())
	 * @param num_acc weighs (currently just number) of accumulated pairs (individually for each pair). WILL BE MODIFIED 
	 * @param fcorr_td_acc Frequency-domain correlation per tile, per pair. . WILL BE MODIFIED 
	 * @param tile_clust_weights Per-tile weights to calculate per-cluster average coordinates. Sum for each cluster is 1.0
	 * @param threadsMax
	 */
	public static void  accumulateCorrelationsAcOnly( // normalize
			final int            clust_size,
			final float [][][]   num_acc,        // number of accumulated tiles [tilesY][tilesX][pair]
			final float [][][][] fcorr_td_acc,   // [tilesY][tilesX][pair][256] sparse transform domain representation of corr pairs
			final double [][]    tile_clust_weights,  // null or  [tilesY][tilesX]
			final int            threadsMax
			) {
		int tX=-1; 
		for (int ity = 0; ity < fcorr_td_acc.length; ity++) if (fcorr_td_acc[ity] != null){
			tX = fcorr_td_acc[ity].length;
			break;
		}
		final int tilesY =    fcorr_td_acc.length;
		final int tilesX =    tX;
		final int clustersX = (int) Math.ceil(1.0 * tilesX / clust_size);
		final int clustersY = (int) Math.ceil(1.0 * tilesY / clust_size);
		final int clusters =  clustersX * clustersY;
		final Thread[] threads = ImageDtt.newThreadArray(threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				public void run() {
					for (int nClust = ai.getAndIncrement(); nClust < clusters; nClust = ai.getAndIncrement()) {
						int clustX = nClust % clustersX;
						int clustY = nClust / clustersX;
						double [] total_weights = null; // for this cluster, per-pair
						double cluster_weight = 0.0;
						for (int ctY = 0; ctY < clust_size; ctY++) {
							int tileY = clustY * clust_size + ctY;
							if (tileY < tilesY) {
								for (int ctX = 0; ctX < clust_size; ctX++) {
									int tileX = clustX * clust_size + ctX;
									if (tileX < tilesX) {
										if (num_acc[tileY][tileX] != null) {
//											if (total_weights == null) {
//												total_weights = new double [num_acc[tileY][tileX].length];
//											}
											for (int pair = 0; pair < num_acc[tileY][tileX].length; pair++) {
												double w = num_acc[tileY][tileX][pair];
												if (w > 0.0) {
													if (total_weights == null) {
														total_weights = new double [num_acc[tileY][tileX].length];
													}
													total_weights[pair] += w; // for this cluster, per-pair
													tile_clust_weights[tileY][tileX] += w; // for this tile
													cluster_weight += w; // all pairs, this cluster
												}
											}
										}
									}
								}
							}
						}
						if (total_weights != null) {
							// normalize per-tile weights, so for each cluster sum is 1.0
							for (int ctY = 0; ctY < clust_size; ctY++) {
								int tileY = clustY * clust_size + ctY;
								if (tileY < tilesY) {
									for (int ctX = 0; ctX < clust_size; ctX++) {
										int tileX = clustX * clust_size + ctX;
										if (tileX < tilesX) {
//											for (int pair = 0; pair < total_weights.length; pair++) {
												tile_clust_weights[tileY][tileX] /= cluster_weight; // for this tile
//											}
										}
									}
								}
							}
							
							// add/normalize correlations from each tile in a cluster to the top-left tile
							int tileX0 = clustX * clust_size;
							int tileY0 = clustY * clust_size;
							if (num_acc[tileY0][tileX0] == null) {
								num_acc[tileY0][tileX0] = new float [total_weights.length];
							}
							if (fcorr_td_acc[tileY0][tileX0] == null) {
								fcorr_td_acc[tileY0][tileX0] = new float [total_weights.length][];
							}
							for (int pair = 0; pair < total_weights.length; pair++) { // if (total_weights[pair] > 0) {
								if (total_weights[pair] > 0) {
									num_acc[tileY0][tileX0][pair] = (float) total_weights[pair];
									for (int ctY = 0; ctY < clust_size; ctY++) {
										int tileY = clustY * clust_size + ctY;
										if (tileY < tilesY) {
											for (int ctX = 0; ctX < clust_size; ctX++) {
												int tileX = clustX * clust_size + ctX;
												if (tileX < tilesX) {
													if ((fcorr_td_acc[tileY][tileX] != null ) &&
															(fcorr_td_acc[tileY][tileX][pair] != null )) {
														if (fcorr_td_acc[tileY0][tileX0][pair] == null) {
															fcorr_td_acc[tileY0][tileX0][pair] = new float [fcorr_td_acc[tileY][tileX][pair].length];
														}
														
														for (int i = 0; i < fcorr_td_acc[tileY0][tileX0][pair].length; i++) {
															fcorr_td_acc[tileY0][tileX0][pair][i] += fcorr_td_acc[tileY][tileX][pair][i];
														}
														/*
													} else {
														if (num_acc[tileY][tileX] != null) { 
															num_acc[tileY][tileX][pair] =    0;
														}
														if (fcorr_td_acc[tileY][tileX] != null) {
															fcorr_td_acc[tileY][tileX][pair] = null;
														}
													}
													if ((ctY != 0) || (ctX != 0)){
														if (num_acc[tileY][tileX] != null) { 
															num_acc[tileY][tileX][pair] =    0;
														}
														if (fcorr_td_acc[tileY][tileX] != null) {
															fcorr_td_acc[tileY][tileX][pair] = null;
														}
														*/
													}
												}
											}
										}
									}
									double k = 1.0 / total_weights[pair];
									for (int i = 0; i < fcorr_td_acc[tileY0][tileX0][pair].length; i++) {
										fcorr_td_acc[tileY0][tileX0][pair][i] *= k;
									}
								}
							}
						}
						// Make null for all pairs (all tiles of empty clusters, all but top left for used ones)
						int aa=0;
						for (int ctY = 0; ctY < clust_size; ctY++) {
							int tileY = clustY * clust_size + ctY;
							if (tileY < tilesY) {
								for (int ctX = 0; ctX < clust_size; ctX++) {
									int tileX = clustX * clust_size + ctX;
									if (tileX < tilesX) {
										if ((ctY != 0) || (ctX != 0) || (total_weights == null)) {
											num_acc[tileY][tileX] =      null; // OK, if that tileX, tileY is not used in tp_tasks
											fcorr_td_acc[tileY][tileX] = null;
										}
									}
								}
							}
						}
						aa++;
					}
				}
			};
		}		      
		ImageDtt.startAndJoin(threads);
		return;
	}
	
	
	
	
	
	@Deprecated	
	public void showVHDisparityHistograms( // needs update
			QuadCLT         ref_scene,
			double          debug_hist_max_disparity,
			double          debug_hist_disparity_step,
			double          debug_hist_max_error,
			double          debug_hist_error_step,
			double          disparity_corr,
			int             dhv_modes,
			int             debug_hist_gap,
			String          suffix, // +nrefine+"-d"+disparity_corr,
			String          suffix2, // +nrefine+"-d"+disparity_corr,
			String []       dbg_titles1,
			double [][]     debug_target, // target disparity
			double [][][]   dhv_mode_data,
			double [][]     dbg_hor,
			double [][]     dbg_vert,
			double []       debug_disparity_range,
			double [][]     dx_dy_range,
			double [][]     bayer_green_hor,
			double [][]     bayer_green_vert
			)

	{
		final int tilesX = ref_scene.getTileProcessor().getTilesX();
		final int tilesY = ref_scene.getTileProcessor().getTilesY();
		final int tiles =tilesX * tilesY;
		final int num_scenes = dbg_titles1.length - 1;
		double [][] debug_dtarget_dx = new double [debug_target.length][debug_target[0].length];
		double [][] debug_dtarget_dy = new double [debug_target.length][debug_target[0].length];
		for (int s = 0; s < debug_target.length; s++) {
			Arrays.fill(debug_dtarget_dx[s],Double.NaN);
			Arrays.fill(debug_dtarget_dy[s],Double.NaN);
			for (int ty = 1; ty < (tilesY - 1); ty++) {
				for (int tx = 1; tx < (tilesX - 1); tx++) {
					int nt = ty * tilesX + tx;
					debug_dtarget_dx[s][nt] = 0.5 * (debug_target[s][nt + 1]      - debug_target[s][nt - 1]);
					debug_dtarget_dy[s][nt] = 0.5 * (debug_target[s][nt + tilesX] - debug_target[s][nt - tilesX]);
				}
			}
		}
		(new ShowDoubleFloatArrays()).showArrays(
				debug_dtarget_dx,
				tilesX,
				tilesY,
				true,
				"d_target_disparity_dx-"+suffix, // nrefine+"-d"+disparity_corr,
				dbg_titles1
				);
		(new ShowDoubleFloatArrays()).showArrays(
				debug_dtarget_dy,
				tilesX,
				tilesY,
				true,
				"d_target_disparity_dy-"+suffix, // +nrefine+"-d"+disparity_corr,
				dbg_titles1
				);
		
		// build 2D histograms: target disparity - horizontal, hor/vert disparity - vertical
		// combo, horizontal - above, vertical - below    			
		int disp_steps = (int) Math.ceil(debug_hist_max_disparity/debug_hist_disparity_step);
		int err_steps =  2*((int) Math.ceil(debug_hist_max_error/debug_hist_error_step)) - 1;
		double err_offs = -err_steps*debug_hist_error_step/2.0 - disparity_corr;

		int hist_width =  disp_steps;
		int hist_height = dhv_modes * (err_steps + debug_hist_gap) - debug_hist_gap;
		int hist_len = hist_height * hist_width;
		double [][] dbg_hist = new double [num_scenes+1][hist_len];
		for (int nscene = 0; nscene <= num_scenes; nscene++) {
			for (int mode = 0; mode < (dhv_modes-1); mode++) {
				Arrays.fill(dbg_hist[nscene],
						(mode* (err_steps +     debug_hist_gap) + err_steps                 ) * hist_width,
						(mode* (err_steps +     debug_hist_gap) + err_steps + debug_hist_gap) * hist_width,
						Double.NaN);
			}
			for (int nt = 0; nt < tiles; nt++ ) {
				if (!Double.isNaN(debug_target[nscene][nt])) {
					int x = (int) Math.floor(debug_target[nscene][nt]/debug_hist_disparity_step);
					if      (x < 0) x = 0;
					else if (x >= hist_width) x = hist_width - 1;
					for (int mode = 0; mode < dhv_modes; mode++) {
						if (mode == 1) { // limit by d/dx
							double d_dx = debug_dtarget_dx[nscene][nt];
							if ((d_dx < dx_dy_range[0][0]) || (d_dx > dx_dy_range[0][1])) {
								continue;
							}
							
						} else if (mode == 2) {// limit by d/dy
							double d_dy = debug_dtarget_dy[nscene][nt];
							if ((d_dy < dx_dy_range[1][0]) || (d_dy > dx_dy_range[1][1])) {
								continue;
							}
						}
						
						double d = dhv_mode_data[mode][nscene][nt]; //
						if (!Double.isNaN(d)) {
							int y = (int) Math.floor((d - err_offs)/debug_hist_error_step);
							if      (y < 0) y = 0;
							else if (y >= err_steps) y = err_steps - 1;
							dbg_hist[nscene][(mode*(err_steps + debug_hist_gap) + y) * hist_width+x] += 1.0;
						}
					}    					
				}
			}
		}

		(new ShowDoubleFloatArrays()).showArrays(
				dbg_hist,
				hist_width,
				hist_height,
				true,
				"histogram_disparity_hor-vert-"+suffix, // nrefine+"-d"+disparity_corr,
				dbg_titles1
				);
		// --- end of histogram 1 (by disparity) ---    			

		// Another histogram by fract disparity			
		double debug_hist_max_error2 = 0.6; //
		int err_steps2 =  2*((int) Math.ceil(debug_hist_max_error2/debug_hist_error_step)) - 1;
		double err_offs2 = -err_steps2 * debug_hist_error_step / 2.0;
		int hist_width2 =  3 * (err_steps2 + debug_hist_gap) - debug_hist_gap;
		int hist_height2 = err_steps;
		int mult = 1; // 2; // 1 - rounded to integer, 2 rounded to even integer
		double rmult = 1.0/mult;
		int hist_len2 = hist_height2 * hist_width2;
		double [][] dbg_hist2 = new double [num_scenes+1][hist_len2];
		double [][][] hor_vert_offs = {dbg_hor, dbg_vert};
		double [][][] hor_ver_err = {dhv_mode_data[1],dhv_mode_data[2]}; // hor, vert
		for (int nscene = 0; nscene <= num_scenes; nscene++) {
			for (int mode = 0; mode < 2; mode++) { // hor, vert only
				Arrays.fill(dbg_hist[nscene],
						(mode* (err_steps +     debug_hist_gap) + err_steps                 ) * hist_width,
						(mode* (err_steps +     debug_hist_gap) + err_steps + debug_hist_gap) * hist_width,
						Double.NaN);
			}
			for (int nt = 0; nt < tiles; nt++ ) {
				double td = debug_target[nscene][nt];
				if (!Double.isNaN(td) && (td >= debug_disparity_range[0]) && (td < debug_disparity_range[1])) {
					for (int mode = 0; mode < 3; mode++) { // hor, vert only
						double offs = rmult * ( hor_vert_offs[(mode < 2 )? mode : 0][nscene][nt]);
						offs -= Math.round(offs);
						if (!Double.isNaN(offs)) {
							int x = (int) Math.floor((offs - err_offs2)/debug_hist_error_step);
							if      (x < 0) x = 0;
							else if (x >= err_steps2) x = err_steps2 - 1;
							double d;
							int y;
							if (mode < 2) {
								if (mode == 0) { // limit by d/dx
									double d_dx = debug_dtarget_dx[nscene][nt];
									if ((d_dx < dx_dy_range[0][0]) || (d_dx > dx_dy_range[0][1])) {
										continue;
									}
								} else if (mode == 1) {// limit by d/dy
									double d_dy = debug_dtarget_dy[nscene][nt];
									if ((d_dy < dx_dy_range[1][0]) || (d_dy > dx_dy_range[1][1])) {
										continue;
									}
								}
								d = rmult * hor_ver_err[mode][nscene][nt];
								y = (int) Math.floor((d - err_offs)/debug_hist_error_step);
							} else {
								d = rmult * hor_vert_offs[1][nscene][nt];
								d -= Math.round(d);
								y = (int) Math.floor((d - err_offs2)/debug_hist_error_step);
							}
							if (!Double.isNaN(d)) {
								if      (y < 0) y = 0;
								else if (y >= err_steps) y = err_steps - 1;
								dbg_hist2[nscene][y * hist_width2 + mode*(err_steps2 + debug_hist_gap) + x] += 1.0;
							}
						}
					}    						
				}
			}
		}
		(new ShowDoubleFloatArrays()).showArrays(
				dbg_hist2,
				hist_width2,
				hist_height2,
				true,
				"histogram_offsets_hor-vert_"+debug_disparity_range[0]+"-"+debug_disparity_range[1]+suffix2+"-mult"+mult,
				//				"_"+nrefine+"-d"+disparity_corr+"-mult"+mult,
				dbg_titles1
				);

		// Another histogram by fract disparity			
		double debug_hist_max_bayer = 1.0; // debug_hist_max_error2
		int err_steps_bayer =  2*((int) Math.ceil(debug_hist_max_bayer/debug_hist_error_step)) - 1;// err_steps2
		double err_offs_bayer = -err_steps_bayer * debug_hist_error_step / 2.0; // err_offs2
		int hist_width_bayer =  3 * (err_steps_bayer + debug_hist_gap) - debug_hist_gap; //hist_width2
		int hist_height_bayer = err_steps; //hist_height2
		int hist_len_bayer = hist_height_bayer * hist_width_bayer;
		double [][] dbg_hist_bayer = new double [num_scenes+1][hist_len_bayer];
		double [][][] hor_vert_bayer = {bayer_green_hor, bayer_green_vert};
//		double [][][] hor_ver_err = {dhv_mode_data[1],dhv_mode_data[2]}; // hor, vert // reuse as is
		for (int nscene = 0; nscene <= num_scenes; nscene++) {
			for (int mode = 0; mode < 2; mode++) { // hor, vert only
				Arrays.fill(dbg_hist[nscene],
						(mode* (err_steps +     debug_hist_gap) + err_steps                 ) * hist_width,
						(mode* (err_steps +     debug_hist_gap) + err_steps + debug_hist_gap) * hist_width,
						Double.NaN);
			}
			for (int nt = 0; nt < tiles; nt++ ) {
				double td = debug_target[nscene][nt];
				if (!Double.isNaN(td) && (td >= debug_disparity_range[0]) && (td < debug_disparity_range[1])) {
					for (int mode = 0; mode < 3; mode++) { // hor, vert only
						double offs = rmult * ( hor_vert_bayer[(mode < 2 )? mode : 0][nscene][nt]);
						offs -= Math.round(offs);
						if (!Double.isNaN(offs)) {
							int x = (int) Math.floor((offs - err_offs_bayer)/debug_hist_error_step);
							if      (x < 0) x = 0;
							else if (x >= err_steps_bayer) x = err_steps_bayer - 1;
							double d;
							int y;
							if (mode < 2) {
								if (mode == 0) { // limit by d/dx
									double d_dx = debug_dtarget_dx[nscene][nt];
									if ((d_dx < dx_dy_range[0][0]) || (d_dx > dx_dy_range[0][1])) {
										continue;
									}
								} else if (mode == 1) {// limit by d/dy
									double d_dy = debug_dtarget_dy[nscene][nt];
									if ((d_dy < dx_dy_range[1][0]) || (d_dy > dx_dy_range[1][1])) {
										continue;
									}
								}
								d = rmult * hor_ver_err[mode][nscene][nt];
								y = (int) Math.floor((d - err_offs)/debug_hist_error_step);
							} else {
								d = rmult * hor_vert_bayer[1][nscene][nt];
								d -= Math.round(d);
								y = (int) Math.floor((d - err_offs_bayer)/debug_hist_error_step);
							}
							if (!Double.isNaN(d)) {
								if      (y < 0) y = 0;
								else if (y >= err_steps) y = err_steps - 1;
								dbg_hist_bayer[nscene][y * hist_width_bayer + mode*(err_steps_bayer + debug_hist_gap) + x] += 1.0;
							}
						}
					}    						
				}
			}
		}
		(new ShowDoubleFloatArrays()).showArrays(
				dbg_hist_bayer,
				hist_width_bayer,
				hist_height_bayer,
				true,
				"histogram_bayer_hor-vert_"+debug_disparity_range[0]+"-"+debug_disparity_range[1]+suffix2+"-mult"+mult,
				//				"_"+nrefine+"-d"+disparity_corr+"-mult"+mult,
				dbg_titles1
				);
		
	}
	
	// what does it do?
	
	private double [][] conditionInitialDS(
			CLTParameters  clt_parameters,
			QuadCLT        scene,
			int debug_level){
		double [][] dls = scene.getDLS();
		if (dls == null) {
			return null;
		}
		return conditionInitialDS(
				false, // boolean        use_conf,       // use configuration parameters, false - use following  
				clt_parameters,
				dls,
				scene,
				debug_level);
	}
	
	/**
	 * Separately filters foreground and background of the combo_dsn_final_in (see COMBO_DSN_TITLES)
	 * modifies outputs (source untouched) :COMBO_DSN_INDX_DISP, COMBO_DSN_INDX_STRENGTH, COMBO_DSN_INDX_LMA,
	 * COMBO_DSN_INDX_DISP_FG
	 * COMBO_DSN_INDX_DISP_BG, COMBO_DSN_INDX_STRENGTH_BG, COMBO_DSN_INDX_LMA_BG, COMBO_DSN_INDX_DISP_BG_ALL
	 * COMBO_DSN_INDX_STRENGTH_BG corresponds to COMBO_DSN_INDX_DISP_BG_ALL
	 * 
	 * @param use_conf
	 * @param clt_parameters
	 * @param combo_dsn_final_in
	 * @param scene
	 * @param debug_level
	 * @return
	 */
	
	private double [][] conditionComboDsnFinal(
			boolean        use_conf,       // use configuration parameters, false - use following  
			CLTParameters  clt_parameters,
			double [][]    combo_dsn_final, // dls,
			QuadCLT        scene,
			int            debugLevel) {
		double [][]    combo_dsn_final_out = combo_dsn_final.clone();
		for (int i = 0; i < combo_dsn_final.length; i++) {
			if (combo_dsn_final[i] != null) {
				combo_dsn_final_out[i] = combo_dsn_final[i].clone();
			}
		}
		
        double [][] dls = {
        		combo_dsn_final[COMBO_DSN_INDX_DISP],
        		combo_dsn_final[COMBO_DSN_INDX_LMA],
        		combo_dsn_final[COMBO_DSN_INDX_STRENGTH]
        };

        double [][] ds = conditionInitialDS(
        		true,                // boolean        use_conf,       // use configuration parameters, false - use following  
        		clt_parameters,      // CLTParameters  clt_parameters,
        		dls,                 // double [][]    dls
        		scene,               // QuadCLT        scene,
        		debugLevel);			
        combo_dsn_final_out[COMBO_DSN_INDX_DISP] =     ds[0];
        combo_dsn_final_out[COMBO_DSN_INDX_STRENGTH] = ds[1];
        combo_dsn_final_out[COMBO_DSN_INDX_DISP_FG] =  ds[0].clone();
        combo_dsn_final_out[COMBO_DSN_INDX_LMA] =      ds[0].clone();
        for (int i = 0; i < combo_dsn_final_out[COMBO_DSN_INDX_LMA].length; i++) {
        	if (Double.isNaN(combo_dsn_final[COMBO_DSN_INDX_LMA][i])) {
        		combo_dsn_final_out[COMBO_DSN_INDX_LMA][i] = Double.NaN; // mark them NaN if original LMA was NaN
        	}
        }
        // BG mode
        double [] bg_lma = combo_dsn_final[COMBO_DSN_INDX_DISP_BG_ALL].clone();
        double [] bg_str = combo_dsn_final[COMBO_DSN_INDX_STRENGTH].clone();
        for (int i = 0; i < bg_lma.length; i++) {
        	if (Double.isNaN(combo_dsn_final[COMBO_DSN_INDX_LMA][i])){
        		bg_lma[i] = Double.NaN;
        	}
        	if (!Double.isNaN(combo_dsn_final[COMBO_DSN_INDX_DISP_BG][i])){
        		bg_lma[i] = combo_dsn_final[COMBO_DSN_INDX_LMA_BG][i];
        	}
        	if (!Double.isNaN(combo_dsn_final[COMBO_DSN_INDX_STRENGTH_BG][i])){
        		bg_str[i] = combo_dsn_final[COMBO_DSN_INDX_STRENGTH_BG][i];
        	}
        }
        double [][] dls_bg = {
        		combo_dsn_final[COMBO_DSN_INDX_DISP_BG_ALL],
        		bg_lma,
        		bg_str};
        double [][] ds_bg = conditionInitialDS(
        		true, // boolean        use_conf,       // use configuration parameters, false - use following  
        		clt_parameters,         // CLTParameters  clt_parameters,
        		dls_bg,                 // double [][]    dls
        		scene,                  // QuadCLT        scene,
        		debugLevel);
        combo_dsn_final_out[COMBO_DSN_INDX_DISP_BG_ALL] = ds_bg[0];
        combo_dsn_final_out[COMBO_DSN_INDX_STRENGTH_BG] = ds_bg[1];
        combo_dsn_final_out[COMBO_DSN_INDX_DISP_BG] =     ds_bg[0].clone();
        combo_dsn_final_out[COMBO_DSN_INDX_LMA_BG] =      ds_bg[0].clone();
        for (int i = 0; i < combo_dsn_final_out[COMBO_DSN_INDX_LMA].length; i++) {
        	if (Double.isNaN(combo_dsn_final[COMBO_DSN_INDX_DISP_BG][i])) {
        		combo_dsn_final_out[COMBO_DSN_INDX_DISP_BG][i] = Double.NaN; // mark them NaN if original LMA was NaN
        		combo_dsn_final_out[COMBO_DSN_INDX_LMA_BG][i] =  Double.NaN; // mark them NaN if original LMA was NaN
        	} 
        	if (Double.isNaN(combo_dsn_final[COMBO_DSN_INDX_LMA_BG][i])) {
        		combo_dsn_final_out[COMBO_DSN_INDX_LMA_BG][i] = Double.NaN; // mark them NaN if original LMA was NaN
        	}
        }
		return combo_dsn_final_out;
	}
	
	
	private double [][] conditionInitialDS(
			boolean        use_conf,       // use configuration parameters, false - use following  
			CLTParameters  clt_parameters,
			double [][]    dls,
			QuadCLT        scene,
			int debug_level)
	{
		int         num_bottom =                      6; // average this number of lowest disparity neighbors (of 8)
		int         num_passes =                    100;
		double      max_change =                      1.0e-3;
		double      min_disparity =                  -0.2;
		double      max_sym_disparity =               0.2;
		double      min_strength_lma =                0.7; // weaker - treat as non-lma
		double      min_strength_replace =            0.05; ///  0.14; /// Before /// - LWIR, after - RGB
		double      min_strength_blur =               0.06; ///  0.2;
		double      sigma =                           2; /// 5;
		int         num_blur =                        1; // 3;
		double      disparity_corr =                  0.0;
		int         outliers_nth_fromextrem =         1; // second from min/max - removes dual-tile max/mins
		double      outliers_tolerance_absolute =     0.2;
		double      outliers_tolerance_relative =     0.02;
		int         outliers_max_iter =             100;
		double      outliers_max_strength2 =          1.0; // 0.5; any
		int         outliers_nth_fromextrem2 =        0; // second from min/max - removes dual-tile max/mins
		double      outliers_tolerance_absolute2 =    0.5;
		double      outliers_tolerance_relative2 =    0.1;
		double      outliers_lma_max_strength =       0.4; // 0.5;
		double      outliers_max_strength =           0.1; ///  0.25;
		double      outliers_from_lma_max_strength =  0.8;
		int         search_radius =                   1;       // Search farther if no LMA neighbor is found closer. Original value - 1 (8 neighbors)
		boolean     remove_no_lma_neib =              false;
		double      diff_from_lma_pos =             100.0;
		double      diff_from_lma_neg =               2.0;
		int         outliers_lma_nth_fromextrem =     0; // 1; 
		int         margin =                          8; // pixels
		double      weak_tolerance_absolute=          0.25;
		double      weak_tolerance_relative=          0.025;
		int         weak_min_neibs =                  5;
		double      strong_strength=                  0.5;
		double      weak_strength=                    0.2; // none is below 
		
		
		if (use_conf) {
			num_bottom                    = clt_parameters.imp.num_bottom;
			num_passes                    = clt_parameters.imp.num_passes;
			max_change                    = clt_parameters.imp.max_change;
			min_disparity                 = clt_parameters.imp.min_disparity;
			max_sym_disparity             = clt_parameters.imp.max_sym_disparity;
			min_strength_lma              = clt_parameters.imp.min_strength_lma;
			min_strength_replace          = clt_parameters.imp.min_strength_replace;
			min_strength_blur             = clt_parameters.imp.min_strength_blur;
			sigma                         = clt_parameters.imp.sigma;
			num_blur                      = clt_parameters.imp.num_blur;
			disparity_corr                = clt_parameters.imp.disparity_corr;
			outliers_nth_fromextrem       = clt_parameters.imp.outliers_nth_fromextrem;
			outliers_tolerance_absolute   = clt_parameters.imp.outliers_tolerance_absolute;
			outliers_tolerance_relative   = clt_parameters.imp.outliers_tolerance_relative;
			outliers_max_iter             = clt_parameters.imp.outliers_max_iter;
			outliers_max_strength2        = clt_parameters.imp.outliers_max_strength2;
			outliers_nth_fromextrem2      = clt_parameters.imp.outliers_nth_fromextrem2;
			outliers_tolerance_absolute2  = clt_parameters.imp.outliers_tolerance_absolute2;
			outliers_tolerance_relative2  = clt_parameters.imp.outliers_tolerance_relative2;
			outliers_lma_max_strength     = clt_parameters.imp.outliers_lma_max_strength; //.4
			outliers_max_strength         = clt_parameters.imp.outliers_max_strength;
			outliers_from_lma_max_strength= clt_parameters.imp.outliers_from_lma_max_strength;
			search_radius                 = clt_parameters.imp.search_radius;
			remove_no_lma_neib            = clt_parameters.imp.remove_no_lma_neib;
			diff_from_lma_pos             = clt_parameters.imp.diff_from_lma_pos;
			diff_from_lma_neg             = clt_parameters.imp.diff_from_lma_neg;
			outliers_lma_nth_fromextrem   = clt_parameters.imp.outliers_lma_nth_fromextrem; //0
			margin                        = clt_parameters.imp.filter_margin;
			weak_tolerance_absolute       = clt_parameters.imp.weak_tolerance_absolute; 
			weak_tolerance_relative       = clt_parameters.imp.weak_tolerance_relative; 
			weak_min_neibs                = clt_parameters.imp.weak_min_neibs;          
			strong_strength               = clt_parameters.imp.strong_strength;         
			weak_strength                 = clt_parameters.imp.weak_strength;           

		}
		
		return conditionInitialDS(
				clt_parameters,                 // CLTParameters  clt_parameters,
				dls,                            // double [][]    dls,
				scene,                          // QuadCLT        scene,
				debug_level,                    // int            debug_level,
				num_bottom,                     // final int         num_bottom,
				num_passes,                     // final int         num_passes,
				max_change,                     // final double      max_change,
				min_disparity,                  // final double      min_disparity,
				max_sym_disparity,              // final double      max_sym_disparity,
				min_strength_lma,               // final double      min_strength_lma,
				min_strength_replace,           // final double      min_strength_replace,
				min_strength_blur,              // final double      min_strength_blur,
				sigma,                          // final double      sigma,
				num_blur,                       // final int         num_blur,
				disparity_corr,                 // final double      disparity_corr,
				outliers_nth_fromextrem,        // final int         outliers_nth_fromextrem,
				outliers_tolerance_absolute,    // final double      outliers_tolerance_absolute,
				outliers_tolerance_relative,    // final double      outliers_tolerance_relative,
				outliers_max_iter,              // final int         outliers_max_iter,
				outliers_max_strength2,         // final double      outliers_max_strength2,
				outliers_nth_fromextrem2,       // final int         outliers_nth_fromextrem2,
				outliers_tolerance_absolute2,   // final double      outliers_tolerance_absolute2,
				outliers_tolerance_relative2,   // final double      outliers_tolerance_relative2,
				outliers_lma_max_strength,      // final double      outliers_lma_max_strength,
				outliers_max_strength,          // final double      outliers_max_strength,
				outliers_from_lma_max_strength, // final double      outliers_from_lma_max_strength,
				search_radius,       // final int         search_radius,       // Search farther if no LMA neighbor is found closer. Original value - 1 (8 neighbors)
				remove_no_lma_neib,  // final boolean     remove_no_lma_neib,  // remove without LMA neighbors
				diff_from_lma_pos,              // final double      diff_from_lma_pos,
				diff_from_lma_neg,              // final double      diff_from_lma_neg,
				outliers_lma_nth_fromextrem,    // final int         outliers_lma_nth_fromextrem, 
				margin,                         // final int         margin)
				weak_tolerance_absolute,	
				weak_tolerance_relative,
				weak_min_neibs,
				strong_strength,
				weak_strength);
    }		
	
	
	
	private double [][] conditionInitialDS(
			CLTParameters  clt_parameters,
			double [][]    dls,
			QuadCLT        scene,
			int            debug_level,
			final int         num_bottom,
			final int         num_passes,
			final double      max_change,
			final double      min_disparity,
			final double      max_sym_disparity,		
			final double      min_strength_lma, // remove weak LMA
			final double      min_strength_replace,
			final double      min_strength_blur,
			final double      sigma,
			final int         num_blur,
			final double      disparity_corr,
			final int         outliers_nth_fromextrem,
			final double      outliers_tolerance_absolute,
			final double      outliers_tolerance_relative,
			final int         outliers_max_iter,
			final double      outliers_max_strength2,
			final int         outliers_nth_fromextrem2,
			final double      outliers_tolerance_absolute2,
			final double      outliers_tolerance_relative2,
			final double      outliers_lma_max_strength,
			final double      outliers_max_strength,
			final double      outliers_from_lma_max_strength,
			final int         search_radius,       // Search farther if no LMA neighbor is found closer. Original value - 1 (8 neighbors)
			final boolean     remove_no_lma_neib,  // remove without LMA neighbors
			final double      diff_from_lma_pos,
			final double      diff_from_lma_neg,
			final int         outliers_lma_nth_fromextrem, 
			final int         margin,
			final double      weak_tolerance_absolute,
			final double      weak_tolerance_relative,
			final int         weak_min_neibs,
			final double      strong_strength,
			final double      weak_strength)
			
	{
		
		final int tilesX = scene.getTileProcessor().getTilesX();
		final int tilesY = scene.getTileProcessor().getTilesY();
		final int transform_size = scene.getTileProcessor().getTileSize();
		final int tiles =tilesX * tilesY;

		String [] dbg_titles = {"str", "lma", "clean-lma", "disp","-lma","by-lma","-nonlma", "few_weak", "old-disp","old-sngl","weak","filled"}; 
		double [][] dbg_img = new double [dbg_titles.length][];
		double [] clean_lma = dls[1].clone();
		for (int i = 0; i <clean_lma.length; i++) {
			if (dls[2][i] < min_strength_lma) {
				clean_lma[i] = Double.NaN;
			}
		}
		//Remove crazy LMA high-disparity tiles
		dbg_img[0] = dls[2].clone();
		dbg_img[1] = dls[1].clone();
		dbg_img[2] = clean_lma.clone();
		dbg_img[3] = dls[0].clone();
		double [] disp_outliers = QuadCLT.removeDisparityLMAOutliers( // nothing removed (trying to remove bad LMA)
				false,                       // final boolean     non_ma,
//				dls, //final double [][] dls,
				new double[][] {dls[1], dls[1], clean_lma}, //final double [][] dls,
				outliers_lma_max_strength,   // final double      max_strength,  // do not touch stronger
				outliers_lma_nth_fromextrem, // final int         nth_fromextrem, // 0 - compare to max/min. 1 - second max/min, ... 
				outliers_tolerance_absolute, // final double      tolerance_absolute,
				outliers_tolerance_relative, // final double      tolerance_relative,
				tilesX,                      //final int         width,               //tilesX
				outliers_max_iter,           // final int         max_iter,
				threadsMax,                  // final int         threadsMax,
				debug_level);                // final int         debug_level)
		dbg_img[4] = disp_outliers.clone();
		disp_outliers = QuadCLT.removeDisparityOutliersByLMA( // removed sky, keeps sky edge near strong objects
				new double[][] {disp_outliers, dls[1], clean_lma}, //final double [][] dls,
				outliers_from_lma_max_strength,  // final double      max_strength,  // do not touch stronger
				diff_from_lma_pos,           // final double      diff_from_lma_pos,   // Difference from farthest FG objects (OK to have large, e.g. 100)
				diff_from_lma_neg,           // final double      diff_from_lma_neg,   // Difference from nearest BG objects (small, as FG are usually more visible)
				search_radius,               // int               search_radius,       // Search farther if no LMA neighbor is found closer. Original value - 1 (8				
				remove_no_lma_neib,          // final boolean     remove_no_lma_neib,  // remove without LMA neighbors
				tilesX,                      // final int         width,               //tilesX
				threadsMax,                  // final int         threadsMax,
				debug_level);                // final int         debug_level)
		dbg_img[5] = disp_outliers.clone();
		// mostly filter infinity, clouds, sky
		disp_outliers = QuadCLT.removeDisparityLMAOutliers( // filter non-lma tiles // removed too few !!!
				true,                       // final boolean     non_ma,
				new double[][] {disp_outliers, dls[1], clean_lma}, //final double [][] dls,
				outliers_max_strength,   // final double      max_strength,  // do not touch stronger
				outliers_nth_fromextrem, // final int         nth_fromextrem, // 0 - compare to max/min. 1 - second max/min, ... 
				outliers_tolerance_absolute, // final double      tolerance_absolute,
				outliers_tolerance_relative, // final double      tolerance_relative,
				tilesX,                      //final int         width,               //tilesX
				outliers_max_iter,           // final int         max_iter,
				threadsMax,                  // final int         threadsMax,
				debug_level);                // final int         debug_level)
		dbg_img[6] = disp_outliers.clone();
		
		disp_outliers = QuadCLT.removeFewWeak( // filter non-lma tiles // removed too few !!!
				new double[][] {disp_outliers, dls[1], clean_lma}, //final double [][] dls,
				strong_strength,             // final double      strong,
				weak_strength,               // final double      weak,
				weak_min_neibs,              // final int         min_neibs, 
				weak_tolerance_absolute,     // final double      tolerance_absolute,
				weak_tolerance_relative,     // final double      tolerance_relative,
				tilesX,                      //final int         width,               //tilesX
				outliers_max_iter,           // final int         max_iter,
				threadsMax,                  // final int         threadsMax,
				debug_level);                // final int         debug_level)
		
		dbg_img[7] = disp_outliers.clone();
		// Pre- 2022 filters, some may be obsolete 
		disp_outliers = QuadCLT.removeDisparityOutliers(
				new double[][] {disp_outliers, clean_lma}, //final double [][] dls,
				outliers_max_strength,       // final double      max_strength,  // do not touch stronger
				outliers_nth_fromextrem,     // final int         nth_fromextrem, // 0 - compare to max/min. 1 - second max/min, ... 
				outliers_tolerance_absolute, // final double      tolerance_absolute,
				outliers_tolerance_relative, // final double      tolerance_relative,
				tilesX,                      // final int         width,
				outliers_max_iter,           // final int         max_iter,
				false,                       // final boolean     fit_completely, // do not add tolerance when replacing
				threadsMax,                  // final int         threadsMax,
				debug_level);                // final int         debug_level)
		dbg_img[8] = disp_outliers.clone();
		// remove extreme single-tile outliers (some may be strong - 0.404)
		disp_outliers = QuadCLT.removeDisparityOutliers(
				new double[][] {disp_outliers, clean_lma}, //final double [][] dls,
				outliers_max_strength2,       // final double      max_strength,  // do not touch stronger
				outliers_nth_fromextrem2,     // final int         nth_fromextrem, // 0 - compare to max/min. 1 - second max/min, ... 
				outliers_tolerance_absolute2, // final double      tolerance_absolute,
				outliers_tolerance_relative2, // final double      tolerance_relative,
				tilesX,                      // final int         width,
				outliers_max_iter,           // final int         max_iter,
				false,                       //  final boolean     fit_completely, // do not add tolerance when replacing
				threadsMax,                  // final int         threadsMax,
				debug_level);                // final int         debug_level)
		dbg_img[9] = disp_outliers.clone();
		double [] disp = QuadCLT.blurWeak(
				new double[][] {disp_outliers, clean_lma}, //final double [][] dls,
				min_strength_blur,    // double      min_strength_blur,
				min_strength_replace, // double      min_strength_replace,
				num_blur,             // int         n,
				tilesX,               // int         width,
				sigma);               // double      sigma);
		dbg_img[10] = disp.clone();
		double [][] ds = {disp,     clean_lma};
		final double [][] ds_filled = QuadCLT.fillDisparityStrength(
				ds,                // final double [][] ds0,
				min_disparity,     // final double      min_disparity,
				max_sym_disparity, // final double      max_sym_disparity, // lower disparity - num_bottom = 8; 
				num_bottom,        // final int         num_bottom, // average this number of lowest disparity neighbors (of 8)
				10000, // num_passes,        // final int         num_passes,
				max_change,        // final double      max_change,
				tilesX,            // final int         width,
				threadsMax,        // final int         threadsMax,
				debug_level); // final int         debug_level)
		dbg_img[11] = ds_filled[0].clone();
		
		if ((debug_level > 0)) {// && (clt_parameters.ofp.enable_debug_images)) {
			(new ShowDoubleFloatArrays()).showArrays(
					dbg_img,
					tilesX,
					tilesY,
					true,
					"filtered-"+scene.getImageName(),
					dbg_titles); //	dsrbg_titles);
		}
		
		// Mask bad tiles
		final boolean [] valid_tiles = new boolean [tiles]; // not used ?
		final double [][] pXpYD = new double [tiles][3];
		final Thread[] threads = ImageDtt.newThreadArray(threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				public void run() {
					for (int nTile = ai.getAndIncrement(); nTile < tiles; nTile = ai.getAndIncrement()) {
						int tileY = nTile / tilesX;
						int tileX = nTile % tilesX;
						pXpYD[nTile][0] =  tileX * transform_size + transform_size/2;
						pXpYD[nTile][1] =  tileY * transform_size + transform_size/2;
						pXpYD[nTile][2] =  ds_filled[0][nTile];
					}
				}
			};
		}		      
		ImageDtt.startAndJoin(threads);
///		TpTask[]  tp_tasks =  
				GpuQuad.setInterTasks(
				scene.getNumSensors(),
				scene.getGeometryCorrection().getSensorWH()[0],
				!scene.hasGPU(), // final boolean             calcPortsCoordinatesAndDerivatives, // GPU can calculate them centreXY
				pXpYD, // final double [][]         pXpYD, // per-tile array of pX,pY,disparity triplets (or nulls)
				null,          // final boolean []          selection, // may be null, if not null do not  process unselected tiles
    			scene.getGeometryCorrection(), // final GeometryCorrection  geometryCorrection,
    			disparity_corr, // final double              disparity_corr,
    			margin, // final int                 margin,      // do not use tiles if their centers are closer to the edges
    			valid_tiles, // final boolean []          valid_tiles,            
    			threadsMax); // final int                 threadsMax)  // maximal number of threads to launch
		
		//FIXME:  not clear where tp_tasks was supposed to go? They are not needed, but  valid_tiles - are

		/*
		scene.getGPU().setInterTasks(
				pXpYD, // final double [][]         pXpYD, // per-tile array of pX,pY,disparity triplets (or nulls)
    			scene.getGeometryCorrection(), // final GeometryCorrection  geometryCorrection,
    			disparity_corr, // final double              disparity_corr,
    			margin, // final int                 margin,      // do not use tiles if their centers are closer to the edges
    			valid_tiles, // final boolean []          valid_tiles,            
    			threadsMax); // final int                 threadsMax)  // maximal number of threads to launch
    			*/
		ai.set(0);
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				public void run() {
					for (int nTile = ai.getAndIncrement(); nTile < tiles; nTile = ai.getAndIncrement()) {
						if (!valid_tiles[nTile]) {
							ds_filled[0][nTile] = Double.NaN;
							ds_filled[1][nTile] = 0.0;
						}
					}
				}
			};
		}		      
		ImageDtt.startAndJoin(threads);
		
		
		if ((debug_level > 0) && (clt_parameters.ofp.enable_debug_images)) {
					dbg_titles = new String[]{"disp", "disp_outliers", "disp_blur", "disp_filled", "str", "str_filled"};
			dbg_img = new double[][]{dls[0], disp_outliers, ds[0], ds_filled[0], ds[1], ds_filled[1]};
			(new ShowDoubleFloatArrays()).showArrays(
					dbg_img,
					tilesX,
					tilesY,
					true,
					"ref_fill_dispatity_gaps-"+scene.getImageName(),
					dbg_titles); //	dsrbg_titles);
		}
		return ds_filled;
	}
	
	
	

	static double [][] getPoseFromErs(
			double k_prev,
			QuadCLT reference_QuadCLT,
			QuadCLT scene_QuadCLT,
			int debug_level)
	{
		ErsCorrection ersCorrection = reference_QuadCLT.getErsCorrection();
		String this_image_name = reference_QuadCLT.getImageName();
		if (debug_level > 0) {
			System.out.println("\n"+this_image_name+":\n"+ersCorrection.extrinsic_corr.toString());
			System.out.println(String.format("%s: ers_wxyz_center=     %f, %f, %f", this_image_name,
					ersCorrection.ers_wxyz_center[0], ersCorrection.ers_wxyz_center[1],ersCorrection.ers_wxyz_center[2] ));
			System.out.println(String.format("%s: ers_wxyz_center_dt=  %f, %f, %f",	this_image_name,
					ersCorrection.ers_wxyz_center_dt[0], ersCorrection.ers_wxyz_center_dt[1],ersCorrection.ers_wxyz_center_dt[2] ));
			System.out.println(String.format("%s: ers_wxyz_center_d2t= %f, %f, %f", this_image_name,
					ersCorrection.ers_wxyz_center_d2t[0], ersCorrection.ers_wxyz_center_d2t[1],ersCorrection.ers_wxyz_center_d2t[2] ));
			System.out.println(String.format("%s: ers_watr_center_dt=  %f, %f, %f", this_image_name,
					ersCorrection.ers_watr_center_dt[0], ersCorrection.ers_watr_center_dt[1],ersCorrection.ers_watr_center_dt[2] ));
			System.out.println(String.format("%s: ers_watr_center_d2t= %f, %f, %f", this_image_name,
					ersCorrection.ers_watr_center_d2t[0], ersCorrection.ers_watr_center_d2t[1],ersCorrection.ers_watr_center_d2t[2] ));
		}
		double dt = 0.0;
		if (scene_QuadCLT == null) {
			scene_QuadCLT = reference_QuadCLT;
		}
		ErsCorrection ersCorrectionPrev = (ErsCorrection) (scene_QuadCLT.geometryCorrection);
		dt = reference_QuadCLT.getTimeStamp() - scene_QuadCLT.getTimeStamp();
		if (dt < 0) {
			k_prev = (1.0-k_prev);
		}
		double frame_time_min = ersCorrectionPrev.line_time*ersCorrectionPrev.pixelCorrectionHeight;// (70fps)
		if (Math.abs(dt) > 2 * frame_time_min) { // 0.15) { // at least two frames TODO: use number of lines* line_time * ...? 
			k_prev = 0.5;
			System.out.println("Non-consecutive frames, dt = "+dt);
		}
		double [] wxyz_center_dt_prev =   ersCorrectionPrev.ers_wxyz_center_dt;
		double [] watr_center_dt_prev =   ersCorrectionPrev.ers_watr_center_dt; // is twice omega!
		double [] wxyz_delta = new double[3];
		double [] watr_delta = new double[3];
		for (int i = 0; i <3; i++) {
			wxyz_delta[i] =              dt * (k_prev * wxyz_center_dt_prev[i] + (1.0-k_prev) * ersCorrection.ers_wxyz_center_dt[i]);
			watr_delta[i] = 0.5 *        dt * (k_prev * watr_center_dt_prev[i] + (1.0-k_prev) * ersCorrection.ers_watr_center_dt[i]);
		}
		watr_delta[0] = -watr_delta[0]; /// TESTING!
		watr_delta[2] = -watr_delta[2]; /// TESTING!
		return new double [][] {wxyz_delta, watr_delta};
	}
	float [][] getInterCorrOffsetsDebug(
			final TpTask[] tp_tasks_ref,
			final TpTask[] tp_tasks,
			final int      tilesX,
			final int      tilesY
			){
		final int tiles = tilesX * tilesY;
		final float [][] offsets = new float [2*numSens+2][tiles];
		final int [] num_inp = new int[tiles];
		for (int i = 0; i < offsets.length; i++) {
			Arrays.fill(offsets[i], Float.NaN);
		}
		final Thread[] threads = ImageDtt.newThreadArray(threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				public void run() {
					for (int nTile = ai.getAndIncrement(); nTile <tp_tasks.length; nTile = ai.getAndIncrement()) {
						TpTask task = tp_tasks[nTile];
						int tile = task.ty * tilesX + task.tx;
						num_inp[tile] = 1;
						for (int i = 0; i < numSens; i++) {
							offsets[2*i+2][tile] = task.xy[i][0];
							offsets[2*i+3][tile] = task.xy[i][1];
						}
					}
				}
			};
		}		      
		ImageDtt.startAndJoin(threads);
		ai.set(0);
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				public void run() {
					for (int nTile = ai.getAndIncrement(); nTile <tp_tasks_ref.length; nTile = ai.getAndIncrement()) {
						TpTask task = tp_tasks_ref[nTile];
						int tile = task.ty * tilesX + task.tx;
						if (num_inp[tile] == 1) {
							num_inp[tile] = 2;
							float sx=0,sy=0;
							for (int i = 0; i < numSens; i++) {
								offsets[2*i+2][tile] -= task.xy[i][0];
								offsets[2*i+3][tile] -= task.xy[i][1];
								sx += offsets[2*i+2][tile];
								sy += offsets[2*i+3][tile];
							}
							offsets[0][tile] = sx/numSens;
							offsets[1][tile] = sy/numSens;
						}
					}
				}
			};
		}		      
		ImageDtt.startAndJoin(threads);
		return offsets;
	}
	/**
	 * Get average X,Y offsets between inter-scene correlated tiles
	 * to eliminate correlation between FPN of the same sensor.   
	 * @param max_offset maximal X, Y offset to keep (normally tile size == 8) 
	 * @param tp_tasks_ref reference tiasks (should be updated from GPU to have 
	 *        actual X,Y for each sensor
	 * @param tp_tasks scene correlated to the reference, should also have
	 *        per-sensor coordinates
	 * @param tilesX number of tile columns
	 * @param tilesY number of tile rows.
	 * @return array of X, Y pairs ([tilesX*tilesY][2]), each may be null if 
	 *         not defined or out of tile. Returns null if no tiles contain a valid offset
	 */
	double [][] getInterCorrOffsets(
			final double   max_offset,
			final TpTask[] tp_tasks_ref,
			final TpTask[] tp_tasks,
			final int      tilesX,
			final int      tilesY
			){
		final int tiles = tilesX * tilesY;
		final double [][] offset_pairs = new double [2*numSens][tiles];
		final int [] num_inp = new int[tiles];
		for (int i = 0; i < offset_pairs.length; i++) {
			Arrays.fill(offset_pairs[i], Float.NaN);
		}
		double [][] offsets = new double[tiles][]; 
		final Thread[] threads = ImageDtt.newThreadArray(threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				public void run() {
					for (int nTile = ai.getAndIncrement(); nTile <tp_tasks.length; nTile = ai.getAndIncrement()) {
						TpTask task = tp_tasks[nTile];
						int tile = task.ty * tilesX + task.tx;
						num_inp[tile] = 1;
						for (int i = 0; i < numSens; i++) {
							offset_pairs[2*i+0][tile] = task.xy[i][0];
							offset_pairs[2*i+1][tile] = task.xy[i][1];
						}
					}
				}
			};
		}		      
		ImageDtt.startAndJoin(threads);
		final AtomicInteger anum_pairs = new AtomicInteger(0);
		ai.set(0);
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				public void run() {
					for (int nTile = ai.getAndIncrement(); nTile <tp_tasks_ref.length; nTile = ai.getAndIncrement()) {
						TpTask task = tp_tasks_ref[nTile];
						int tile = task.ty * tilesX + task.tx;
						if (num_inp[tile] == 1) {
							num_inp[tile] = 2;
							double sx=0,sy=0;
							for (int i = 0; i < numSens; i++) {
								sx += offset_pairs[2*i+0][tile] - task.xy[i][0];
								sy += offset_pairs[2*i+1][tile] - task.xy[i][1];
							}
							sx /= numSens;
							sy /= numSens;
							if ((Math.abs(sx) <= max_offset) && (Math.abs(sy) <= max_offset )) {
								offsets[tile] = new double[] {sx,sy};
								anum_pairs.getAndIncrement();
							}
						}
					}
				}
			};
		}		      
		ImageDtt.startAndJoin(threads);
		if (anum_pairs.get() > 0) {
			return offsets;
		}
		return null;
	}
	
	public double [][][]  interCorrPair( // return [tilesX*telesY]{ref_pXpYD, dXdYS}
			CLTParameters      clt_parameters,			
			QuadCLT            ref_scene,
			double []          ref_disparity, // null or alternative reference disparity  
			QuadCLT            scene,
			double []          scene_xyz,
			double []          scene_atr,
			final boolean []   selection, // may be null, if not null do not  process unselected tiles
			final int          margin,
			final int          sensor_mask_inter, // The bitmask - which sensors to correlate, -1 - all.
			final float [][][] accum_2d_corr, // if [1][][] - return accumulated 2d correlations (all pairs)final float [][][]   accum_2d_corr, // if [1][][] - return accumulated 2d correlations (all pairs)
			final float [][]   dbg_corr_fpn,
			boolean            near_important, // do not reduce weight of the near tiles
			boolean            all_fpn,        // do not lower thresholds for non-fpn (used during search)
			int                imp_debug_level, // MANUALLY change to 2 for debug!
			int                debug_level)
	{
		TileProcessor tp = ref_scene.getTileProcessor();
		// Temporary reusing same ref scene ******
		boolean scene_is_ref_test =    clt_parameters.imp.scene_is_ref_test; // false; // true;
		boolean show_2d_correlations = clt_parameters.imp.show2dCorrelations(imp_debug_level); // true;
		boolean show_motion_vectors =  clt_parameters.imp.showMotionVectors(imp_debug_level); // true;
		boolean show_render_ref =      clt_parameters.imp.renderRef(imp_debug_level); // false; //true;
		boolean show_render_scene =    clt_parameters.imp.renderScene(imp_debug_level); // false; // true;
		boolean toRGB =                clt_parameters.imp.toRGB  ; // true;
		boolean show_coord_motion =    clt_parameters.imp.showCorrMotion(imp_debug_level); // mae its own
	    int     erase_clt = (toRGB? clt_parameters.imp.show_color_nan : clt_parameters.imp.show_mono_nan) ? 1:0;
	    boolean use_lma_dsi =          clt_parameters.imp.use_lma_dsi;
		boolean mov_en =               clt_parameters.imp.mov_en; // true;  // enable detection/removal of the moving objects during pose matching
	    boolean mov_debug_images =     clt_parameters.imp.showMovementDetection(imp_debug_level);
	    int mov_debug_level =          clt_parameters.imp.movDebugLevel(imp_debug_level);

	    boolean fpn_remove =           clt_parameters.imp.fpn_remove;
	    double  fpn_max_offset =       clt_parameters.imp.fpn_max_offset;
	    double  fpn_radius =           clt_parameters.imp.fpn_radius;
	    boolean fpn_ignore_border =    clt_parameters.imp.fpn_ignore_border; // only if fpn_mask != null - ignore tile if maximum touches fpn_mask
	    
		if (scene_is_ref_test) {
			scene_xyz = ZERO3.clone();
			scene_atr = ZERO3.clone();
			ref_scene = scene; // or opposite: 	scene = ref_scene;

		}
		int tilesX =         tp.getTilesX();
		int tilesY =         tp.getTilesY();
		/*
		if (dbg_corr_fpn != null) { // 2*16 or 2*17 (average, individual)
			for (int i = 0; i < dbg_corr_fpn.length; i++) {
				dbg_corr_fpn[i] = new float [tilesY * tilesX];
				Arrays.fill(dbg_corr_fpn[i], Float.NaN);
			}
		}
		*/
		double [][][] coord_motion = null; // new double [2][tilesX*tilesY][];
		final double [][][] motion_vectors = show_motion_vectors?new double [tilesY *tilesX][][]:null;
		final float  [][][] fclt_corr = ((accum_2d_corr != null) || show_2d_correlations) ?
				(new float [tilesX * tilesY][][]) : null;
		if (debug_level > 1) {
			System.out.println("interCorrPair():   "+IntersceneLma.printNameV3("ATR",scene_atr)+
					" "+IntersceneLma.printNameV3("XYZ",scene_xyz));
		}
		ImageDtt image_dtt;
		image_dtt = new ImageDtt(
				numSens,
				clt_parameters.transform_size,
				clt_parameters.img_dtt,
				ref_scene.isAux(),
				ref_scene.isMonochrome(),
				ref_scene.isLwir(),
				clt_parameters.getScaleStrength(ref_scene.isAux()),
				ref_scene.getGPU());
		if (ref_scene.getGPU() != null) {
			ref_scene.getGPU().setGpu_debug_level(debug_level - 4); // monitor GPU ops >=-1
		}
		final double disparity_corr = 0.0; // (z_correction == 0) ? 0.0 : geometryCorrection.getDisparityFromZ(1.0/z_correction);
//ref_disparity	
		if (ref_disparity == null) {
			ref_disparity = ref_scene.getDLS()[use_lma_dsi?1:0];
		}
		double [][] ref_pXpYD = transformToScenePxPyD( // full size - [tilesX*tilesY], some nulls
				null, // final Rectangle [] extra_woi,    // show larger than sensor WOI (or null)
				ref_disparity, // dls[0],  // final double []   disparity_ref, // invalid tiles - NaN in disparity (maybe it should not be masked by margins?)
				ZERO3,         // final double []   scene_xyz, // camera center in world coordinates
				ZERO3,         // final double []   scene_atr, // camera orientation relative to world frame
				ref_scene,     // final QuadCLT     scene_QuadClt,
				ref_scene);    // final QuadCLT     reference_QuadClt)
		final double gpu_sigma_corr =     clt_parameters.getGpuCorrSigma(ref_scene.isMonochrome());
		final double gpu_sigma_rb_corr =  ref_scene.isMonochrome()? 1.0 : clt_parameters.gpu_sigma_rb_corr;
		final double gpu_sigma_log_corr = clt_parameters.getGpuCorrLoGSigma(ref_scene.isMonochrome());
		TpTask[] tp_tasks_ref =  GpuQuad.setInterTasks(
				ref_scene.getNumSensors(),
				ref_scene.getErsCorrection().getSensorWH()[0],
				!ref_scene.hasGPU(),          // final boolean             calcPortsCoordinatesAndDerivatives, // GPU can calculate them centreXY
				ref_pXpYD,                    // final double [][]         pXpYD, // per-tile array of pX,pY,disparity triplets (or nulls)
				selection,                    // final boolean []          selection, // may be null, if not null do not  process unselected tiles
				ref_scene.getErsCorrection(), // final GeometryCorrection  geometryCorrection,
				disparity_corr,               // final double              disparity_corr,
				margin,                       // final int                 margin,      // do not use tiles if their centers are closer to the edges
				null,                         // final boolean []          valid_tiles,            
				threadsMax);                  // final int                 threadsMax)  // maximal number of threads to launch

		double [][] scene_pXpYD = transformToScenePxPyD( // will be null for disparity == NaN, total size - tilesX*tilesY
				null, // final Rectangle [] extra_woi,    // show larger than sensor WOI (or null)
				ref_disparity, // dls[0], // final double []   disparity_ref, // invalid tiles - NaN in disparity (maybe it should not be masked by margins?)
				scene_xyz,    // final double []   scene_xyz, // camera center in world coordinates
				scene_atr,    // final double []   scene_atr, // camera orientation relative to world frame
				scene,        // final QuadCLT     scene_QuadClt,
				scene_is_ref_test ? null: ref_scene);   // final QuadCLT     reference_QuadClt)
		TpTask[] tp_tasks =  GpuQuad.setInterTasks(
				scene.getNumSensors(),
				scene.getErsCorrection().getSensorWH()[0],
				!scene.hasGPU(),          // final boolean             calcPortsCoordinatesAndDerivatives, // GPU can calculate them centreXY
				scene_pXpYD,              // final double [][]         pXpYD, // per-tile array of pX,pY,disparity triplets (or nulls)
				selection,                // final boolean []          selection, // may be null, if not null do not  process unselected tiles
				scene.getErsCorrection(), // final GeometryCorrection  geometryCorrection,
				disparity_corr,           // final double              disparity_corr,
				margin,                   // final int                 margin,      // do not use tiles if their centers are closer to the edges
				null,                     // final boolean []          valid_tiles,            
				threadsMax);              // final int                 threadsMax)  // maximal number of threads to launch
		
		if (ref_scene.hasGPU()) { // direct convert to reference_clt
			float  [][][][]     fcorr_td =  null; // no accumulation, use data in GPU
			ref_scene.saveQuadClt(); // to re-load new set of Bayer images to the GPU (do nothing for CPU) and Geometry
			image_dtt.setReferenceTD( // tp_tasks_ref will be updated
					erase_clt,
					null,                       // final int []              wh,               // null (use sensor dimensions) or pair {width, height} in pixels
					clt_parameters.img_dtt,     // final ImageDttParameters  imgdtt_params,    // Now just extra correlation parameters, later will include, most others
					true, // final boolean             use_reference_buffer,
					tp_tasks_ref,               // final TpTask[]            tp_tasks,
					clt_parameters.gpu_sigma_r, // final double              gpu_sigma_r,     // 0.9, 1.1
					clt_parameters.gpu_sigma_b, // final double              gpu_sigma_b,     // 0.9, 1.1
					clt_parameters.gpu_sigma_g, // final double              gpu_sigma_g,     // 0.6, 0.7
					clt_parameters.gpu_sigma_m, // final double              gpu_sigma_m,     //  =       0.4; // 0.7;
					threadsMax,                 // final int                 threadsMax,       // maximal number of threads to launch
					debug_level);               // final int                 globalDebugLevel);
			if (show_render_ref) {
				ImagePlus imp_render_ref = ref_scene.renderFromTD (
						-1,                  // final int         sensor_mask,
						false,               // boolean             merge_channels,
						clt_parameters,                                 // CLTParameters clt_parameters,
						clt_parameters.getColorProcParameters(ref_scene.isAux()), //ColorProcParameters colorProcParameters,
						clt_parameters.getRGBParameters(),              //EyesisCorrectionParameters.RGBParameters rgbParameters,
						null, // int []  wh,
						toRGB, // boolean toRGB,
						true,  //boolean use_reference
						"GPU-SHIFTED-REFERENCE"); // String  suffix)
				imp_render_ref.show();
			}
			
			scene.saveQuadClt(); // to re-load new set of Bayer images to the GPU (do nothing for CPU) and Geometry
			image_dtt.interCorrTD(
					clt_parameters.img_dtt,     // final ImageDttParameters  imgdtt_params,    // Now just extra correlation parameters, later will include, most others
					tp_tasks,                   // final TpTask[]            tp_tasks,
					fcorr_td,                   // final float  [][][][]     fcorr_td,        // [tilesY][tilesX][pair][4*64] transform domain representation of 6 corr pairs
					clt_parameters.gpu_sigma_r, // final double              gpu_sigma_r,     // 0.9, 1.1
					clt_parameters.gpu_sigma_b, // final double              gpu_sigma_b,     // 0.9, 1.1
					clt_parameters.gpu_sigma_g, // final double              gpu_sigma_g,     // 0.6, 0.7
					clt_parameters.gpu_sigma_m, // final double              gpu_sigma_m,     //  =       0.4; // 0.7;
					gpu_sigma_rb_corr,          // final double              gpu_sigma_rb_corr,    //  = 0.5; // apply LPF after accumulating R and B correlation before G, monochrome ? 1.0 :
					gpu_sigma_corr,             // final double              gpu_sigma_corr,       //  =    0.9;gpu_sigma_corr_m
					gpu_sigma_log_corr,         // final double              gpu_sigma_log_corr,   // hpf to reduce dynamic range for correlations
					clt_parameters.corr_red,    // final double              corr_red, // +used
					clt_parameters.corr_blue,   // final double              corr_blue,// +used
					sensor_mask_inter,          // final int                 sensor_mask_inter, // The bitmask - which sensors to correlate, -1 - all.
					threadsMax,                 // final int                 threadsMax,       // maximal number of threads to launch
					debug_level);               // final int                 globalDebugLevel);
			if (show_render_scene) {
				ImagePlus imp_render_scene = scene.renderFromTD (
						-1,                  // final int         sensor_mask,
						false,               // boolean             merge_channels,
						clt_parameters,                                 // CLTParameters clt_parameters,
						clt_parameters.getColorProcParameters(ref_scene.isAux()), //ColorProcParameters colorProcParameters,
						clt_parameters.getRGBParameters(),              //EyesisCorrectionParameters.RGBParameters rgbParameters,
						null, // int []  wh,
						toRGB, // boolean toRGB,
						false, //boolean use_reference
						"GPU-SHIFTED-SCENE"); // String  suffix)

				imp_render_scene.show();
			}
			if (dbg_corr_fpn != null) { // 2*16 or 2*17 (average, individual)
				float [][] foffsets = getInterCorrOffsetsDebug(
						tp_tasks_ref, // final TpTask[] tp_tasks_ref,
						tp_tasks, // final TpTask[] tp_tasks,
						tilesX, // final int      tilesX,
						tilesY); // final int      tilesY
				for (int i = 0; (i < dbg_corr_fpn.length) && (i < foffsets.length); i++) {
					dbg_corr_fpn[i] = foffsets[i];
				}
			}
			
			double [][] fpn_offsets = null;
			if (fpn_remove) {
				fpn_offsets = getInterCorrOffsets(
						fpn_max_offset, // final double   max_offset,
						tp_tasks_ref,   // final TpTask[] tp_tasks_ref,
						tp_tasks,       // final TpTask[] tp_tasks,
						tilesX,         // final int      tilesX,
						tilesY);        // final int      tilesY);
			}
			
			double            half_disparity = near_important ? 0.0 : clt_parameters.imp.half_disparity;
	        double [][][]     dcorr_tiles = (fclt_corr != null)? (new double [tp_tasks.length][][]):null;
	        // will use num_acc with variable number of accumulations (e.g. clusters)
	        //all_fpn
	        double min_str =     all_fpn ? clt_parameters.imp.min_str_fpn : clt_parameters.imp.min_str;
	        double min_str_sum = all_fpn ? clt_parameters.imp.min_str_sum_fpn : clt_parameters.imp.min_str_sum;
	        coord_motion = image_dtt.clt_process_tl_interscene(       // convert to pixel domain and process correlations already prepared in fcorr_td and/or fcorr_combo_td
	        		clt_parameters.img_dtt,            // final ImageDttParameters  imgdtt_params,   // Now just extra correlation parameters, later will include, most others
	        		fcorr_td,                          // final float  [][][][]     fcorr_td,        // [tilesY][tilesX][pair][4*64] transform domain representation of all selected corr pairs
	        		null,                              // float [][][]              num_acc,         // number of accumulated tiles [tilesY][tilesX][pair] (or null). Can be inner null if not used in tp_tasks
	        		null,                              // double []                 dcorr_weight,    // alternative to num_acc, compatible with CPU processing (only one non-zero enough)
	        		clt_parameters.gpu_corr_scale,     // final double              gpu_corr_scale,  //  0.75; // reduce GPU-generated correlation values
	        		clt_parameters.getGpuFatZeroInter(ref_scene.isMonochrome()), // final double              gpu_fat_zero,    // clt_parameters.getGpuFatZero(is_mono);absolute == 30.0
	        		image_dtt.transform_size - 1,      // final int                 gpu_corr_rad,    // = transform_size - 1 ?
	        		// The tp_tasks data should be decoded from GPU to get coordinates
	        		tp_tasks,                          // final TpTask []           tp_tasks,        // data from the reference frame - will be applied to LMA for the integrated correlations
	        		// to be converted to float (may be null)
	        		dcorr_tiles,                       // final double  [][][]      dcorr_tiles,     // [tile][pair_abs, sparse][(2*transform_size-1)*(2*transform_size-1)] // if null - will not calculate
	        		ref_pXpYD,                         // final double [][]         pXpYD,           // pXpYD for the reference scene
	        		fpn_offsets,                       // final double [][]         fpn_offsets,     // null, or per-tile X,Y offset to be blanked
	        		fpn_radius,                        // final double              fpn_radius,      // radius to be blanked around FPN offset center
	        		fpn_ignore_border,                 // final boolean             fpn_ignore_border, // only if fpn_mask != null - ignore tile if maximum touches fpn_mask			
	        		motion_vectors,                    // final double [][][]       motion_vectors,  // [tilesY*tilesX][][] -> [][][num_sel_sensors+1][2]
	        		clt_parameters.imp.run_poly,       // final boolean             run_poly,        // polynomial max, if false - centroid
	        		clt_parameters.imp.use_partial,    // final boolean             use_partial,     // find motion vectors for individual pairs, false - for sum only
	        		clt_parameters.imp.centroid_radius,// final double              centroid_radius, // 0 - use all tile, >0 - cosine window around local max
	        		clt_parameters.imp.n_recenter,     // final int                 n_recenter,      // when cosine window, re-center window this many times
	        		clt_parameters.imp.td_weight,      // final double  td_weight,    // mix correlations accumulated in TD with 
	        		clt_parameters.imp.pd_weight,      // final double  pd_weight,    // correlations (post) accumulated in PD
	        		clt_parameters.imp.td_nopd_only,   // final boolean td_nopd_only, // only use TD accumulated data if no safe PD is available for the tile.
	        		min_str,                           // final double              min_str_nofpn,         //  = 0.25;
	        		min_str_sum,                       // final double              min_str_sum_nofpn,     // = 0.8; // 5;
	        		clt_parameters.imp.min_str_fpn,    // final double              min_str,         //  = 0.25;
	        		clt_parameters.imp.min_str_sum_fpn,// final double              min_str_sum,     // = 0.8; // 5;
	        		clt_parameters.imp.min_neibs,      // final int                 min_neibs,       //   2;	   // minimal number of strong neighbors (> min_str)
	        		clt_parameters.imp.weight_zero_neibs, // final double              weight_zero_neibs,//  0.2;   // Reduce weight for no-neib (1.0 for all 8)
	        		half_disparity,                    // final double              half_disparity,  //   5.0;   // Reduce weight twice for this disparity
	        		clt_parameters.imp.half_avg_diff,  // final double              half_avg_diff,   //   0.2;   // when L2 of x,y difference from average of neibs - reduce twice
	        		clt_parameters.tileX,              // final int                 debug_tileX,
	        		clt_parameters.tileY,      	       // final int                 debug_tileY,
	        		threadsMax,                        // final int                 threadsMax,       // maximal number of threads to launch
	        		debug_level);
	        // final int                 globalDebugLevel);
	        if (coord_motion == null) {
	        	System.out.println("clt_process_tl_interscene() returned null");
	        	return null;
	        }
	        if (mov_en) {
	        	String debug_image_name = mov_debug_images ? (scene.getImageName()+"-"+ref_scene.getImageName()+"-movements"): null;
	        	boolean [] move_mask = getMovementMask(
	        			clt_parameters, // CLTParameters clt_parameters,
	        			coord_motion[1], // double [][]   motion, // only x,y,w components
	        			tilesX, // int           tilesX,
	        			debug_image_name, // String        debug_image_name,
	        			mov_debug_level); // int           debug_level);
	        	if (move_mask != null) {
	        		for (int nTile=0; nTile < move_mask.length; nTile++) if (move_mask[nTile]) {
	        			coord_motion[1][nTile]= null;
	        		}
	        	}
	        }
	        
		    float [][][] fclt_corr1 = ImageDtt.convertFcltCorr( // partial length, matching corr_indices = gpuQuad.getCorrIndices(); // also sets num_corr_tiles
		            dcorr_tiles, // double [][][] dcorr_tiles,// [tile][sparse, correlation pair][(2*transform_size-1)*(2*transform_size-1)] // if null - will not calculate
		            fclt_corr);  // float  [][][] fclt_corr) //  new float [tilesX * tilesY][][] or null
		    if (show_2d_correlations) { // visualize prepare ref_scene correlation data
				float [][] dbg_corr_rslt_partial = ImageDtt.corr_partial_dbg( // not used in lwir
						fclt_corr1, // final float  [][][]     fcorr_data,       // [tile][pair][(2*transform_size-1)*(2*transform_size-1)] // if null - will not calculate
						image_dtt.getGPU().getCorrIndices(), // tp_tasks,  // final TpTask []         tp_tasks,        //
						tilesX,    //final int                tilesX,
						tilesY,    //final int                tilesX,
						2*image_dtt.transform_size - 1,	// final int               corr_size,
						1000, // will be limited by available layersfinal int               layers0,
						clt_parameters.corr_border_contrast, // final double            border_contrast,
						threadsMax, // final int               threadsMax,     // maximal number of threads to launch
						debug_level); // final int               globalDebugLevel)
				
				String [] titles = new String [dbg_corr_rslt_partial.length]; // dcorr_tiles[0].length];
				for (int i = 0; i < titles.length; i++) {
					if (i== (titles.length - 1)) {
						titles[i] = "sum";
					} else {
						titles[i] = "sens-"+i;
					}
				}
				
				// titles.length = 15, corr_rslt_partial.length=16!
				(new ShowDoubleFloatArrays()).showArrays( // out of boundary 15
						dbg_corr_rslt_partial,
						tilesX*(2*image_dtt.transform_size),
						tilesY*(2*image_dtt.transform_size),
						true,
						scene.getImageName()+"-"+ref_scene.getImageName()+"-interscene",
						titles);
			}
		} else {
			throw new IllegalArgumentException ("interCorrPair(): CPU mode not supported");
		}
		if (motion_vectors != null) {
			int num_sens = 0;
			int num_rslt = 0;
			find_num_layers:
				for (int nt = 0; nt <motion_vectors.length; nt++) {
					if (motion_vectors[nt] != null) {
						num_sens = motion_vectors[nt].length;
						for (int i = 0; i < num_sens; i++) {
							if (motion_vectors[nt][i] != null) {
								num_rslt = motion_vectors[nt][i].length;
								break find_num_layers;
							}
						}
					}
				}
			double [][] dbg_img = new double [num_sens * num_rslt][tilesX*tilesY];
			String [] titles = new String [dbg_img.length];
			for (int ns = 0; ns < num_sens; ns++) {
				String ss = (ns == (num_sens - 1))?"sum":(""+ns);
				for (int nr = 0; nr < num_rslt; nr++) {
					int indx = ns*num_rslt+nr;
					Arrays.fill(dbg_img[indx],Double.NaN);
					String sr="";
					switch (nr) {
					case 0:sr ="X";break;
					case 1:sr ="Y";break;
					case 2:sr ="S";break;
					default:
						sr = ""+nr;
					}
					titles[indx] = ss+":"+sr;
					for (int nt = 0; nt <motion_vectors.length; nt++) {
						if ((motion_vectors[nt] != null) && (motion_vectors[nt][ns] != null) ) {
							dbg_img[indx][nt] = motion_vectors[nt][ns][nr];
						}
					}
				}
			}
			(new ShowDoubleFloatArrays()).showArrays( // out of boundary 15
					dbg_img,
					tilesX,
					tilesY,
					true,
					scene.getImageName()+"-"+ref_scene.getImageName()+"-motion_vectors",
					titles);
			
		}
		if (show_coord_motion) {
			//coord_motion
			String [] mvTitles = {"dx", "dy", "conf", "pX", "pY","Disp","defined"}; // ,"blurX","blurY", "blur"};
			double [][] dbg_img = new double [mvTitles.length][tilesX*tilesY];
			for (int l = 0; l < dbg_img.length; l++) {
				Arrays.fill(dbg_img[l], Double.NaN);
			}
			for (int nTile = 0; nTile <  coord_motion[0].length; nTile++) {
				if (coord_motion[0][nTile] != null) {
					for (int i = 0; i <3; i++) {
						dbg_img[3+i][nTile] = coord_motion[0][nTile][i];
					}
				}
				if (coord_motion[1][nTile] != null) {
					for (int i = 0; i <3; i++) {
						dbg_img[0+i][nTile] = coord_motion[1][nTile][i];
					}
				}
				dbg_img[6][nTile] = ((coord_motion[0][nTile] != null)?1:0)+((coord_motion[0][nTile] != null)?2:0);
			}
			(new ShowDoubleFloatArrays()).showArrays( // out of boundary 15
					dbg_img,
					tilesX,
					tilesY,
					true,
					scene.getImageName()+"-"+ref_scene.getImageName()+"-coord_motion",
					mvTitles);
		}
		if (debug_level > 0){
			int num_defined = 0;
			double sum_strength = 0.0;
			for (int i = 0; i < coord_motion[1].length; i++) if (coord_motion[1][i] != null){
				sum_strength += coord_motion[1][i][2];
				num_defined++;
			}
			System.out.println ("interCorrPair(): num_defined = "+num_defined+
						", sum_strength = "+sum_strength+", avg_strength = "+(sum_strength/num_defined));
		}
		
		return coord_motion;
	}
	
	
	public static boolean [] getMovementMask(
			CLTParameters clt_parameters,
			double [][]   motion, // only x,y,w components
			int           tilesX,
			String        debug_image_name,
			int           debug_level)
	{
		boolean       show_debug_images = (debug_image_name != null);

		String [] mvTitles = {"dx", "dy", "conf", "blur", "blurX","blurY","clust","mask"};
		boolean mov_en =         clt_parameters.imp.mov_en; // true;  // enable detection/removal of the moving objects during pose matching
		if (!mov_en) {
			return null;
		}
		double  mov_sigma =      clt_parameters.imp.mov_sigma; // 1.5;   // pix - weighted-blur offsets before detection
		// next two to prevent motion detection while errors are too big
		double  mov_max_std =    clt_parameters.imp.mov_max_std; // 0.5;   // pix
		double  mov_thresh_rel = clt_parameters.imp.mov_thresh_rel; // 5.0;   // exceed average error
		double  mov_thresh_abs=  clt_parameters.imp.mov_thresh_abs; // 1.0;   // sqrt(dx^2+dy^2) in moving areas 
		double  mov_clust_max =  clt_parameters.imp.mov_clust_max; // 1.5;   // cluster maximum should exceed threshold this times
		int     mov_grow =       clt_parameters.imp.mov_grow; // 4;     // grow detected moving area
		double  mov_max_std2 =   mov_max_std * mov_max_std;
		double  mov_thresh_rel2 = mov_thresh_rel* mov_thresh_rel;
		double  mov_thresh_abs2 = mov_thresh_abs * mov_thresh_abs;
		int    tilesY = motion.length / tilesX;
		double [] wa = new double[3];
		for (int nTile = 0; nTile <  motion.length; nTile++) if (motion[nTile] != null) {
			double w = motion[nTile][2]; //?  1.0; //
			wa[0]+=w;
			wa[1]+=w*motion[nTile][0];
			wa[2]+=w*motion[nTile][1];
		}
		wa[1] /= wa[0];
		wa[2] /= wa[0];
		double [][] mov_obj = new double[show_debug_images?4:2][tilesX*tilesY];
		double sl2 = 0;
		for (int nTile = 0; nTile <  motion.length; nTile++) if (motion[nTile] != null) {
			double w = motion[nTile][2]; //?  1.0; // 
			mov_obj[0][nTile]= w;
			double x = motion[nTile][0]-wa[1];
			double y = motion[nTile][1]-wa[2];
			double l2 = x*x + y*y;
			mov_obj[1][nTile]= w*l2;
			if (show_debug_images) {
				mov_obj[2][nTile]= w*x;
				mov_obj[3][nTile]= w*y;
			}
			sl2+=w*l2;
		}
		sl2/=wa[0];
		if (debug_level>0) {
			System.out.println("getMovementMask(): std="+Math.sqrt(sl2));
		}
		double [][] dbg_img = null;
		if(show_debug_images) {		
			dbg_img = new double [mvTitles.length][tilesX*tilesY];
			for (int l = 0; l < dbg_img.length; l++) {
				Arrays.fill(dbg_img[l], Double.NaN);
			}
		}		
		boolean [] move_mask = null; 
		process_movements:
		{	
			if (sl2 > mov_max_std2) {
				if (debug_level>0) {
					System.out.println("getMovementMask(): Too high scene mismatch: std="+Math.sqrt(sl2)+" > "+mov_max_std+", aborting movement processing");
				}
				break process_movements;
			}
			DoubleGaussianBlur gb = new DoubleGaussianBlur();
			gb.blurDouble(mov_obj[0], tilesX, tilesY, mov_sigma, mov_sigma, 0.01);
			gb.blurDouble(mov_obj[1], tilesX, tilesY, mov_sigma, mov_sigma, 0.01);
			for (int nTile = 0; nTile <  motion.length; nTile++) {
				mov_obj[1][nTile] /= mov_obj[0][nTile];
			}		
			if (show_debug_images) {
				gb.blurDouble(mov_obj[2], tilesX, tilesY, mov_sigma, mov_sigma, 0.01);
				gb.blurDouble(mov_obj[3], tilesX, tilesY, mov_sigma, mov_sigma, 0.01);
				for (int nTile = 0; nTile <  motion.length; nTile++) {
					mov_obj[2][nTile] /= mov_obj[0][nTile];
					mov_obj[3][nTile] /= mov_obj[0][nTile];
				}
				for (int nTile = 0; nTile <  motion.length; nTile++) {
					 if (motion[nTile] != null) {
						 dbg_img[0][nTile] =  motion[nTile][0];
						 dbg_img[1][nTile] =  motion[nTile][1];
						 dbg_img[2][nTile] =  motion[nTile][2];
					 }
					 dbg_img[3][nTile] =  mov_obj[1][nTile];				
					 dbg_img[4][nTile] =  mov_obj[2][nTile];				
					 dbg_img[5][nTile] =  mov_obj[3][nTile];
				}				
			}
			double mov_thresh = sl2 * mov_thresh_rel2;
			if (mov_thresh < mov_thresh_abs2) {
				mov_thresh = mov_thresh_abs2;
			}
			int num_over = 0;
			boolean [] over_thresh = new boolean [mov_obj[0].length];
			for (int nTile = 0; nTile <  motion.length; nTile++) if (mov_obj[1][nTile] > mov_thresh){
				over_thresh[nTile] = true;
				num_over ++;
			}
			if (num_over == 0) {
				if (debug_level>0) {
					System.out.println("getMovementMask(): No moving tiles exceed ="+Math.sqrt(mov_thresh)+", aborting movement processing");
				}
				break process_movements;
			}
			// Find clusters
			TileNeibs tn = new TileNeibs(tilesX,tilesY);

			int [] clusters = tn.enumerateClusters(
					over_thresh, // boolean [] tiles,
					false); // boolean ordered)
			// check clusters contain super-threshold tile
			double max_cluster = mov_thresh*mov_clust_max; // each cluster should have larger tile value to survive
			int max_clust_num = 0;
			for (int nTile=0; nTile < clusters.length;nTile++) if (clusters[nTile] > max_clust_num) {
				max_clust_num = clusters[nTile];
			}
			double [] clust_max = new double [max_clust_num];
			for (int nTile=0; nTile < clusters.length;nTile++) if (clusters[nTile] > 0) {
				int i = clusters[nTile]-1;
				if (mov_obj[1][nTile] > clust_max[i]) {
					clust_max[i] = mov_obj[1][nTile];
				}
			}
			int num_good_clusters = 0;
			boolean [] good_clusters = new boolean[clust_max.length];
			for (int i = 0; i < clust_max.length; i++) if (clust_max[i] >= max_cluster){
				good_clusters[i] = true;
				num_good_clusters ++;
			}		
			if (debug_level>0) {
				System.out.println("getMovementMask(): Got "+max_clust_num+" clusters, of them good (>"+max_cluster+") - "+num_good_clusters );
				for (int i = 0; i < clust_max.length; i++) {
					System.out.println("getMovementMask(): Cluster "+i+(good_clusters[i]?"*":" ")+", max = "+clust_max[i]);
				}
			}
			if(show_debug_images) {		
				for (int nTile = 0; nTile <  motion.length; nTile++) if (clusters[nTile] > 0){
					 dbg_img[6][nTile] = clusters[nTile];
				}
			}
			
			if (num_good_clusters == 0) {
				if (debug_level>0) {
					System.out.println("getMovementMask(): No strong enough moving clusters (of total "+clust_max.length+"), aborting movement processing");
				}
				break process_movements;
			}
			move_mask = new boolean [clusters.length];
			for (int nTile=0; nTile < clusters.length; nTile++) if (clusters[nTile] > 0) {
				move_mask[nTile] = good_clusters[clusters[nTile] -1];
			}
			tn.growSelection(
					mov_grow,  // int        grow,           // grow tile selection by 1 over non-background tiles 1: 4 directions, 2 - 8 directions, 3 - 8 by 1, 4 by 1 more
					move_mask, // boolean [] tiles,
					null);     // boolean [] prohibit);
			if(show_debug_images) {		
				for (int nTile = 0; nTile <  motion.length; nTile++) if (move_mask[nTile]){
					 dbg_img[7][nTile] = 1.0;
				}
			}

		}
		if(show_debug_images) {	
			(new ShowDoubleFloatArrays()).showArrays( // out of boundary 15
					dbg_img,
					tilesX,
					tilesY,
					true,
					debug_image_name,
					mvTitles);
		}
		if (debug_level > -2) { // 0
			if (move_mask != null) {
				int num_mov = 0;
				for (int nTile=0; nTile<move_mask.length;nTile++) if (move_mask[nTile]) {
					num_mov++;
				}
				System.out.println("getMovementMask(): Created moving objects mask of "+num_mov+" tiles.");

			} else {
				System.out.println("getMovementMask(): No moving objects mask is created.");
			}
		}
		return move_mask;
	}
	
	
	public double[][]  adjustPairsLMA(
			CLTParameters  clt_parameters,			
			QuadCLT        reference_QuadCLT,
			QuadCLT        scene_QuadCLT,
			double []      camera_xyz0,
			double []      camera_atr0,
			boolean[]      param_select,
			double []      param_regweights,
			double []      rms_out, // null or double [2]
			double [][]    dbg_img, // null or [3][]
			double         max_rms,
			int            debug_level)
	{
		TileProcessor tp = reference_QuadCLT.getTileProcessor();
		final int iscale = 8;
		boolean blur_reference = false;
		
		int tilesX =         tp.getTilesX();
		int tilesY =         tp.getTilesY();
		int transform_size = tp.getTileSize();
		int macroTilesX =    tilesX/transform_size;
		int macroTilesY =    tilesY/transform_size;
		
		if (debug_level > -1) {
			System.out.println("adjustPairsLMA():  "+IntersceneLma.printNameV3("ATR",camera_atr0)+
					" "+IntersceneLma.printNameV3("XYZ",camera_xyz0));
		}
		if (dbg_img != null) {
			for (int i = 0; i < dbg_img.length; i++) {
				dbg_img[i] = new double [macroTilesX * macroTilesY * clt_parameters.ilp.ilma_num_corr];
				Arrays.fill(dbg_img[i],Double.NaN);
			}
		}
		if (clt_parameters.ofp.enable_debug_images && (debug_level > 0)) {
			compareRefSceneTiles(
					"-before_LMA",      // String suffix,
					blur_reference,    // boolean blur_reference,
					camera_xyz0,       // double [] camera_xyz0,
					camera_atr0,       // double [] camera_atr0,
					reference_QuadCLT, // QuadCLT reference_QuadCLT,
					scene_QuadCLT,     // QuadCLT scene_QuadCLT,
					iscale);           // int iscale) // 8
		}
		IntersceneLma intersceneLma = new IntersceneLma(
				this, // OpticalFlow opticalFlow
				clt_parameters.ilp.ilma_thread_invariant);
		int nlma = 0;
		int lmaResult = -1;
		// debug: had to clt_parameters.ofp.debug_level_iterate=-10 to stop images
		for (nlma = 0; nlma < clt_parameters.ilp.ilma_num_corr; nlma++) {
			boolean last_run = nlma == ( clt_parameters.ilp.ilma_num_corr - 1);
			int macroTiles = macroTilesX * macroTilesY; 
			double [][] flowXY = new double [macroTiles][2]; // zero pre-shifts
			//		double [][] flowXY_frac = new double [macroTiles][]; // Will contain fractional X/Y shift for CLT
			double [][] reference_tiles_macro = new double [macroTiles][];
			// hack below to pass nlma for debugging. Change  clt_parameters.ilp.ilma_num_corr to -11 after first pass
			int test_debug_level = clt_parameters.ilp.ilma_debug_invariant? (clt_parameters.ofp.debug_level_iterate-nlma): -11;
			double [][] vector_XYS = correlate2DIterate( // returns optical flow and confidence
					clt_parameters.img_dtt, // final ImageDttParameters  imgdtt_params,
					// for prepareSceneTiles()			
					camera_xyz0,                                 // final double []   scene_xyz,     // camera center in world coordinates
					camera_atr0,                                 // final double []   scene_atr,     // camera orientation relative to world frame
					scene_QuadCLT,                               // final QuadCLT     scene_QuadClt,
					reference_QuadCLT,                           // final QuadCLT     reference_QuadClt,
					reference_tiles_macro,                       //			final double [][] reference_tiles_macro,
					clt_parameters.ofp.center_occupancy_ref,         // final double      reference_center_occupancy,   // fraction of remaining  tiles in the center 8x8 area (<1.0)
					// flowXY should be initialized to all pairs of zeros (or deliberate pixel offset pairs if initial error is too high, will be modified with each iteration
					flowXY,                   // final double [][] flowXY, // per macro tile {mismatch in image pixels in X and Y directions // initialize to [reference_tiles.length][2]
					clt_parameters.ofp.tolerance_absolute_inter, // final double      tolerance_absolute, // absolute disparity half-range in each tile
					clt_parameters.ofp.tolerance_relative_inter, // final double      tolerance_relative, // relative disparity half-range in each tile
					clt_parameters.ofp.occupancy_inter,          // final double      occupancy,          // fraction of remaining  tiles (<1.0)
					clt_parameters.ofp.num_laplassian,           // final int         num_passes,
					clt_parameters.ofp.change_laplassian,        // final double      max_change,
					// for correlate2DSceneToReference ()
					clt_parameters.ofp.chn_weights,              // final double []   chn_weights, // absolute, starting from strength (strength,r,b,g)
					clt_parameters.ofp.corr_sigma,               // final double      corr_sigma,
					clt_parameters.ofp.fat_zero,                 //  final double      fat_zero,
					clt_parameters.ofp.late_normalize_iterate,   // final boolean     late_normalize,
					// for correlation2DToVectors_CM()
					clt_parameters.ofp.iradius_cm,               // final int         iradius,      // half-size of the square to process 
					clt_parameters.ofp.dradius_cm,               // final double      dradius,      // weight calculation (1/(r/dradius)^2 + 1)
					clt_parameters.ofp.refine_num_cm,            // final int         refine_num,   // number of iterations to apply weights around new center
					clt_parameters.ofp.num_refine_all,           // final int         num_run_all, // run all tiles for few iterations before filtering
					clt_parameters.ofp.max_refines,              // final int         max_tries,
					// for recalculateFlowXY()
					clt_parameters.ofp.magic_scale,              // final double      magic_scale, // 0.85 for CM
					clt_parameters.ofp.min_change,               // final double      min_change,
					clt_parameters.ofp.best_neibs_num,           // final int         best_num,
					clt_parameters.ofp.ref_stdev,                // final double      ref_stdev,
					test_debug_level, // clt_parameters.ofp.debug_level_iterate-nlma,      // final int         debug_level)
					clt_parameters.ofp.enable_debug_images);     //final boolean     enable_debug_images)
			if (dbg_img != null) {
				for (int iy = 0; iy < macroTilesY; iy++) {
					for (int ix = 0; ix < macroTilesX; ix++) {
						int ii = iy*macroTilesX + ix;
						if  (vector_XYS[ii] != null) {
							int ii1 = (iy * clt_parameters.ilp.ilma_num_corr + nlma) * macroTilesX + ix;
							for (int jj = 0; jj < dbg_img.length; jj++) {
								dbg_img[jj][ii1] = vector_XYS[ii][jj];
							}
						}
					}
				}
			}

			if (clt_parameters.ofp.enable_debug_images && (debug_level > 2)) {
				String dbg_title = "OpticalFlow-"+scene_QuadCLT.getImageName()+"-"+reference_QuadCLT.getImageName()+"-iteration_"+nlma;;
				showVectorXYConfidence(
						dbg_title, // String      title,
						vector_XYS, // double [][] flowXYS,
						macroTilesX); // int         width)	
			}
			int n = removeOutliers(
					clt_parameters.ofp.nsigma, // double nsigma, 1.5 - 2.0
					vector_XYS); // double [][] flowXYS)
			if (debug_level > 1) {
				System.out.println("Removed "+n+" outliers");
			}
			int n2 = removeOutliers(
					clt_parameters.ofp.nsigma2, // double nsigma, 1.5 - 2.0
					vector_XYS); // double [][] flowXYS)
			if (debug_level > 1) {
				System.out.println("Removed "+n2+" outliers in a second pass, total removed:"+(n+n2));
			}
			if (clt_parameters.ofp.enable_debug_images && (debug_level > 0)) {
				if ((debug_level > 1) || (nlma == 0)) { 
					String dbg_title = "OpticalFlowFiltered-"+scene_QuadCLT.getImageName()+"-"+reference_QuadCLT.getImageName()+"-iteration_"+nlma;
					showVectorXYConfidence(
							dbg_title, // String      title,
							vector_XYS, // double [][] flowXYS,
							macroTilesX); // int         width)
				}
			}

			if (clt_parameters.ilp.ilma_debug_invariant) { // dbg_img != null) {
				/*
				long [] long_xyz0 =  new long [camera_xyz0.length];
				long [] long_atr0 =  new long [camera_atr0.length];
				for (int i = 00; i < camera_xyz0.length;i++) {
					System.out.println("camera_xyz0["+i+"] = " + camera_xyz0[i] );
					long_xyz0[i] = Double.doubleToRawLongBits(camera_xyz0[i]);
				}
				for (int i = 0; i < camera_xyz0.length;i++) {
					System.out.println("camera_atr0["+i+"] = " + camera_atr0[i] );
					long_atr0[i] = Double.doubleToRawLongBits(camera_atr0[i]);
				}
				for (int i = 0; i < camera_xyz0.length;i++) {
					System.out.println("long_xyz0["+i+"] = " + long_xyz0[i] );
				}
				for (int i = 0; i < camera_xyz0.length;i++) {
					System.out.println("long_atr0["+i+"] = " + long_atr0[i] );
				}
				System.out.println("camera_xyz0.hashCode() = " +           camera_xyz0.hashCode());
				System.out.println("camera_atr0.hashCode() = " +           camera_atr0.hashCode());
				System.out.println("long_xyz0.hashCode() = " +             long_xyz0.hashCode());
				System.out.println("long_atr0.hashCode() = " +             long_atr0.hashCode());
				System.out.println("scene_QuadCLT.hashCode() = " +         scene_QuadCLT.hashCode());
				System.out.println("reference_QuadCLT.hashCode() = " +     reference_QuadCLT.hashCode());
				System.out.println("param_select.hashCode() = " +          param_select.hashCode());
				System.out.println("param_regweights.hashCode() = " +      param_regweights.hashCode());
				System.out.println("vector_XYS.hashCode() = " +            vector_XYS.hashCode());
				System.out.println("reference_tiles_macro.hashCode() = " + reference_tiles_macro.hashCode());
				*/
				try {
					System.out.println("getChecksum(camera_xyz0) = " +           IntersceneLma.getChecksum(camera_xyz0));
					System.out.println("getChecksum(camera_atr0) = " +           IntersceneLma.getChecksum(camera_atr0));
//					System.out.println("getChecksum(scene_QuadCLT) = " +         getChecksum(scene_QuadCLT));
//					System.out.println("getChecksum(reference_QuadCLT) = " +     getChecksum(reference_QuadCLT));
					System.out.println("getChecksum(param_select) = " +          IntersceneLma.getChecksum(param_select));
					System.out.println("getChecksum(param_regweights) = " +      IntersceneLma.getChecksum(param_regweights));
					System.out.println("getChecksum(vector_XYS) = " +            IntersceneLma.getChecksum(vector_XYS));
					System.out.println("getChecksum(reference_tiles_macro) = " + IntersceneLma.getChecksum(reference_tiles_macro));
				} catch (NoSuchAlgorithmException | IOException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
//getChecksum
			intersceneLma.prepareLMA(
					camera_xyz0,                                    // final double []   scene_xyz0,     // camera center in world coordinates (or null to use instance)
					camera_atr0,                                    // final double []   scene_atr0,     // camera orientation relative to world frame (or null to use instance)
					// reference atr, xyz are considered 0.0
					scene_QuadCLT,                                  // final QuadCLT     scene_QuadClt,
					reference_QuadCLT,                              // final QuadCLT     reference_QuadClt,
					param_select,                                   // final boolean[]   param_select,
					param_regweights,                               // final double []   param_regweights,
					vector_XYS,                                     // final double [][] vector_XYS, // optical flow X,Y, confidence obtained from the correlate2DIterate()
					reference_tiles_macro,                          // final double [][] centers,    // macrotile centers (in pixels and average disparities
					(nlma == 0),                                    // boolean           first_run,
					clt_parameters.ilp.ilma_debug_level);           // final int         debug_level)
			
			lmaResult = intersceneLma.runLma(
					clt_parameters.ilp.ilma_lambda, // double lambda,           // 0.1
					clt_parameters.ilp.ilma_lambda_scale_good, //  double lambda_scale_good,// 0.5
					clt_parameters.ilp.ilma_lambda_scale_bad,  // double lambda_scale_bad, // 8.0
					clt_parameters.ilp.ilma_lambda_max,        // double lambda_max,       // 100
					clt_parameters.ilp.ilma_rms_diff,          // double rms_diff,         // 0.001
					clt_parameters.ilp.ilma_num_iter,          // int    num_iter,         // 20
					last_run,                                  // boolean last_run,
					clt_parameters.ilp.ilma_debug_level);      // int    debug_level)
			// debugging:
//			String [] lines1 = intersceneLma.printOldNew(false); // boolean allvectors)
//			System.out.println("lmaResult="+lmaResult+", RMS="+intersceneLma.getLastRms()[0]);
//			for (String line : lines1) {
//				System.out.println(line);
//			}
			
			if (lmaResult < 0) {
				System.out.println("LMA failed, nlma="+nlma);
				break;
			}
			camera_xyz0 = intersceneLma.getSceneXYZ(false); // true for initial values
			camera_atr0 = intersceneLma.getSceneATR(false); // true for initial values
			
			if (clt_parameters.ofp.enable_debug_images && (debug_level > 1)) {
				compareRefSceneTiles(
						"iteration_"+nlma,      // String suffix,
						blur_reference,    // boolean blur_reference,
						camera_xyz0,       // double [] camera_xyz0,
						camera_atr0,       // double [] camera_atr0,
						reference_QuadCLT, // QuadCLT reference_QuadCLT,
						scene_QuadCLT,     // QuadCLT scene_QuadCLT,
						iscale);           // int iscale) // 8
			}
			if (lmaResult <= 1) {
				break;
			}
		}
		if (debug_level > 2) { // duplecate - already reported in intersceneLma.runLma()
			System.out.println("LMA: full RMS="+intersceneLma.getLastRms()[0]+", pure RMS="+intersceneLma.getLastRms()[1]);
			String [] lines = intersceneLma.printOldNew(false); // boolean allvectors)
			for (String line : lines) {
				System.out.println(line);
			}
		}

		if (clt_parameters.ofp.enable_debug_images && (debug_level == 1))  {
///		if (!clt_parameters.ofp.enable_debug_images || (clt_parameters.ofp.enable_debug_images && (debug_level == 1)))  {
			compareRefSceneTiles(
					"-after_lma",      // String suffix,
					blur_reference,    // boolean blur_reference,
					camera_xyz0,       // double [] camera_xyz0,
					camera_atr0,       // double [] camera_atr0,
					reference_QuadCLT, // QuadCLT reference_QuadCLT,
					scene_QuadCLT,     // QuadCLT scene_QuadCLT,
					iscale);           // int iscale) // 8
		}
		if (rms_out != null) {
			rms_out[0] = intersceneLma.getLastRms()[0];
			rms_out[1] = intersceneLma.getLastRms()[1];
			//if (lmaResult < 0) { last_rms[0]
		}
		if (max_rms > 0.0) {
			if (lmaResult < 0) { // = 0) {
				return null;
			}
			if (!(intersceneLma.getLastRms()[0] <= max_rms)) {
				System.out.println("RMS failed: "+intersceneLma.getLastRms()[0]+" >= " + max_rms);
				return null;
			}
		}
		return new double [][] {camera_xyz0, camera_atr0};
	}

	public int  reAdjustPairsLMAInterscene( // after combo dgi is available and preliminary poses are known
			CLTParameters  clt_parameters,
			boolean []     reliable_ref, // null or bitmask of reliable reference tiles			
    		QuadCLT []     quadCLTs,
			int            debugLevel)
	{   
		int earliest_scene = 0;
	    boolean use_combo_dsi =        clt_parameters.imp.use_combo_dsi;
	    boolean use_lma_dsi =          clt_parameters.imp.use_lma_dsi;
    	int ref_index = quadCLTs.length-1;
		int tilesX =  quadCLTs[ref_index].getTileProcessor().getTilesX();
        int tilesY =  quadCLTs[ref_index].getTileProcessor().getTilesY();
        double [] disparity_raw = new double [tilesX * tilesY];
        Arrays.fill(disparity_raw,clt_parameters.disparity);
        double [][] combo_dsn_final = quadCLTs[ref_index].readDoubleArrayFromModelDirectory(
        			"-INTER-INTRA-LMA", // String      suffix,
        			0, // int         num_slices, // (0 - all)
        			null); // int []      wh);
        double [][] dls = {
        		combo_dsn_final[COMBO_DSN_INDX_DISP],
        		combo_dsn_final[COMBO_DSN_INDX_LMA],
        		combo_dsn_final[COMBO_DSN_INDX_STRENGTH]
        };
        double [][] ds = conditionInitialDS(
        		true,                // boolean        use_conf,       // use configuration parameters, false - use following  
        		clt_parameters,      // CLTParameters  clt_parameters,
        		dls,                 // double [][]    dls
        		quadCLTs[ref_index], // QuadCLT        scene,
        		debugLevel);			
//        double [] disparity_fg = ds[0]; // combo_dsn_final[COMBO_DSN_INDX_DISP_FG];
        double [] interscene_ref_disparity = null; // keep null to use old single-scene disparity for interscene matching
        if (use_combo_dsi) {
        	interscene_ref_disparity = ds[0].clone(); // use_lma_dsi ?
        	if (use_lma_dsi) {
        		for (int i = 0; i < interscene_ref_disparity.length; i++) {
        			if (Double.isNaN(dls[1][i])) {
        				interscene_ref_disparity[i] = Double.NaN;
        			}
        		}
        	}
        }
		ErsCorrection ers_reference = quadCLTs[ref_index].getErsCorrection();
        double [][][] dxyzatr_dt = new double[quadCLTs.length][][];
		double [][][] scenes_xyzatr = new double [quadCLTs.length][][]; // previous scene relative to the next one
		scenes_xyzatr[ref_index] = new double[2][3]; // all zeros
		// should have at least next or previous non-null
		double maximal_series_rms = 0.00;
		
        for (int nscene = ref_index; nscene >= earliest_scene; nscene--) {
        	if ((quadCLTs[nscene] == null) ||
        			((nscene != ref_index) &&
        			((ers_reference.getSceneXYZ(quadCLTs[nscene].getImageName())== null) ||
        			(ers_reference.getSceneATR(quadCLTs[nscene].getImageName())== null)))) {
        		earliest_scene = nscene + 1;
        		break;
        	}
			String ts = quadCLTs[nscene].getImageName();
        	int nscene0 = nscene - 1;
        	if ((nscene0 < 0) ||
        			(quadCLTs[nscene0]== null)||
        			(ers_reference.getSceneXYZ(quadCLTs[nscene0].getImageName())== null) ||
        			(ers_reference.getSceneATR(quadCLTs[nscene0].getImageName())== null)) {
        		nscene0 = nscene;
        	}
        	int nscene1 = nscene + 1;
        	if ((nscene1 > ref_index) || (quadCLTs[nscene1]== null)) {
        		nscene1 = nscene;
        	}
        	if (nscene1 == nscene0) {
        		System.out.println("**** Isoloated scene!!! skippiung... now may only happen for a ref_scene****");
        		earliest_scene = nscene + 1;
        		break;
        	}
        	double dt = quadCLTs[nscene1].getTimeStamp() - quadCLTs[nscene0].getTimeStamp();
        	String ts0 = quadCLTs[nscene0].getImageName();
        	String ts1 = quadCLTs[nscene1].getImageName();
    		double [] scene_xyz0 = ers_reference.getSceneXYZ(ts0);
    		double [] scene_atr0 = ers_reference.getSceneATR(ts0);
    		if (scene_xyz0 == null) {
    			System.out.println ("BUG: No egomotion data for timestamp "+ts0);
    			System.out.println ("Need to re-run with Force egomotion calculation");
        		earliest_scene = nscene + 1;
        		break;
    		}
    		double [] scene_xyz1 = (nscene1== ref_index)? ZERO3:ers_reference.getSceneXYZ(ts1);
    		double [] scene_atr1 = (nscene1== ref_index)? ZERO3:ers_reference.getSceneATR(ts1);
    		dxyzatr_dt[nscene] = new double[2][3];
    		for (int i = 0; i < 3; i++) {
    			dxyzatr_dt[nscene][0][i] = 0.0; // (scene_xyz1[i]-scene_xyz0[i])/dt;
    			dxyzatr_dt[nscene][1][i] = (scene_atr1[i]-scene_atr0[i])/dt;
    		}
			double []   scene_xyz_pre = ZERO3;
			double []   scene_atr_pre = ZERO3;
			quadCLTs[nscene].getErsCorrection().setErsDt( // set for ref also (should be set before non-ref!)
					dxyzatr_dt[nscene][0], // double []    ers_xyz_dt,
					dxyzatr_dt[nscene][1]); // double []    ers_atr_dt)(ers_scene_original_xyz_dt);
			int debug_scene = -15;
			if (nscene != ref_index) {
				if (nscene == debug_scene) {
					System.out.println("nscene = "+nscene);
					System.out.println("nscene = "+nscene);
				}
				scene_xyz_pre = ers_reference.getSceneXYZ(ts);
				scene_atr_pre = ers_reference.getSceneATR(ts);
				double []      lma_rms = new double[2];
				scenes_xyzatr[nscene] = adjustPairsLMAInterscene(
						clt_parameters,                                 // CLTParameters  clt_parameters,
						quadCLTs[ref_index],                            // QuadCLT reference_QuadCLT,
						interscene_ref_disparity,                       // double []        ref_disparity, // null or alternative reference disparity
						reliable_ref,  // boolean []     reliable_ref, // null or bitmask of reliable reference tiles
						quadCLTs[nscene],                               // QuadCLT scene_QuadCLT,
						scene_xyz_pre,                                  // xyz
						scene_atr_pre,                                  // atr
						clt_parameters.ilp.ilma_lma_select,                                  // final boolean[]   param_select,
						clt_parameters.ilp.ilma_regularization_weights,                              //  final double []   param_regweights,
						lma_rms,                                        // double []      rms, // null or double [2]
						clt_parameters.imp.max_rms,                     // double         max_rms,
						clt_parameters.imp.debug_level);                // 1); // -1); // int debug_level);
				if (scenes_xyzatr[nscene] == null) {
					System.out.println("LMA failed at nscene = "+nscene+". Truncating series.");
	        		earliest_scene = nscene + 1;
					if (debugLevel > -3) {
						System.out.println("reAdjustPairsLMAInterscene "+nscene+" (of "+ quadCLTs.length+") "+
								quadCLTs[ref_index].getImageName() + "/" + ts+
								" FAILED. Setting earliest_scene to "+earliest_scene);
					}
					// set this and all previous to null
					for (; nscene >= 0 ; nscene--) {
						ers_reference.addScene(quadCLTs[nscene].getImageName(), null);
					}
	        		break;
				}
				
//				System.out.println("lma_rms={"+lma_rms[0]+","+lma_rms[1]+"}");
				ers_reference.addScene(ts,
						scenes_xyzatr[nscene][0],
						scenes_xyzatr[nscene][1],
						quadCLTs[nscene].getErsCorrection().getErsXYZ_dt(),	// same as dxyzatr_dt[nscene][0], just keep for future adjustments?	
						quadCLTs[nscene].getErsCorrection().getErsATR_dt()	// same as dxyzatr_dt[nscene][1], just keep for future adjustments?		
						);
			    if (lma_rms[0] > maximal_series_rms) {
			        maximal_series_rms = lma_rms[0];
			    }

				if (debugLevel > -3) {
					System.out.println("reAdjustPairsLMAInterscene "+nscene+" (of "+ quadCLTs.length+") "+
							quadCLTs[ref_index].getImageName() + "/" + ts+
							" Done. RMS="+lma_rms[0]+", maximal so far was "+maximal_series_rms);
				}
			}
        } // for (int nscene = ref_index; nscene > earliest_scene; nscene--) {
		
		if (debugLevel > -4) {
			System.out.println("All multi scene passes are Done. Maximal RMSE was "+maximal_series_rms);
		}

		quadCLTs[ref_index].saveInterProperties( // save properties for interscene processing (extrinsics, ers, ...)
				null, // String path,             // full name with extension or w/o path to use x3d directory
				debugLevel+1);
		return earliest_scene;
	}
	
	public double[][]  adjustPairsLMAInterscene(
			CLTParameters  clt_parameters,			
			QuadCLT        reference_QuadClt,
			double []      ref_disparity, // null or alternative reference disparity
			boolean []     reliable_ref, // null or bitmask of reliable reference tiles
			QuadCLT        scene_QuadClt,
			double []      camera_xyz,
			double []      camera_atr,
			boolean[]      param_select,
			double []      param_regweights,
			double []      rms_out, // null or double [2]
			double         max_rms,
			int            debug_level)
	{
		int        margin = clt_parameters.imp.margin;
		int        sensor_mask_inter = clt_parameters.imp.sensor_mask_inter ; //-1;
		float [][][] facc_2d_img = new float [1][][];
		IntersceneLma intersceneLma = new IntersceneLma(
				this, // OpticalFlow opticalFlow
				clt_parameters.ilp.ilma_thread_invariant);
		int lmaResult = -1;
		boolean last_run = false;
		double[] camera_xyz0 = camera_xyz.clone();
		double[] camera_atr0 = camera_atr.clone();
		double [][][] coord_motion = null;
		boolean show_corr_fpn = debug_level > -1; //  -3; *********** Change to debug FPN correleation *** 
		
		int nlma = 00;
		for (; nlma < clt_parameters.imp.max_cycles; nlma ++) {
			
			boolean near_important = nlma > 0;
			coord_motion = interCorrPair( // new double [tilesY][tilesX][][];
					clt_parameters,      // CLTParameters  clt_parameters,
					reference_QuadClt,   // QuadCLT reference_QuadCLT,
					ref_disparity,       // double []        ref_disparity, // null or alternative reference disparity  
					scene_QuadClt,       // QuadCLT scene_QuadCLT,
					camera_xyz0,         // xyz
					camera_atr0,         // pose[1], // atr
					reliable_ref,        // ****null,                // final boolean [] selection, // may be null, if not null do not  process unselected tiles
					margin,              // final int        margin,
					sensor_mask_inter,   // final int        sensor_mask_inter, // The bitmask - which sensors to correlate, -1 - all.
					facc_2d_img,         // final float [][][]   accum_2d_corr, // if [1][][] - return accumulated 2d correlations (all pairs)final float [][][]   accum_2d_corr, // if [1][][] - return accumulated 2d correlations (all pairs)
					null,                //	final float [][] dbg_corr_fpn,
					near_important,      // boolean            near_important, // do not reduce weight of the near tiles
					false,               // boolean            all_fpn,        // do not lower thresholds for non-fpn (used during search)
					clt_parameters.imp.debug_level, // int                imp_debug_level,
					debug_level); // 1); // -1); // int debug_level);
	        if (coord_motion == null) {
	        	System.out.println("adjustPairsLMAInterscene() returned null");
	        	return null;
	        }
			intersceneLma.prepareLMA(
					camera_xyz0,         // final double []   scene_xyz0,     // camera center in world coordinates (or null to use instance)
					camera_atr0,         // final double []   scene_atr0,     // camera orientation relative to world frame (or null to use instance)
					// reference atr, xyz are considered 0.0
					scene_QuadClt,       // final QuadCLT     scene_QuadClt,
					reference_QuadClt,   // final QuadCLT     reference_QuadClt,
					param_select,        // final boolean[]   param_select,
					param_regweights,    // final double []   param_regweights,
					coord_motion[1],     // final double [][] vector_XYS, // optical flow X,Y, confidence obtained from the correlate2DIterate()
					coord_motion[0],     // final double [][] centers,    // macrotile centers (in pixels and average disparities
					(nlma == 0),                                    // boolean           first_run,
					clt_parameters.imp.debug_level);           // final int         debug_level)
			lmaResult = intersceneLma.runLma(
					clt_parameters.ilp.ilma_lambda, // double lambda,           // 0.1
					clt_parameters.ilp.ilma_lambda_scale_good, //  double lambda_scale_good,// 0.5
					clt_parameters.ilp.ilma_lambda_scale_bad,  // double lambda_scale_bad, // 8.0
					clt_parameters.ilp.ilma_lambda_max,        // double lambda_max,       // 100
					clt_parameters.ilp.ilma_rms_diff,          // double rms_diff,         // 0.001
					clt_parameters.imp.max_LMA,                // int    num_iter,         // 20
					last_run,                                  // boolean last_run,
					debug_level);           // int    debug_level)
			if (lmaResult < 0) {
				System.out.println("Interscene adjustment failed, lmaResult="+lmaResult+" < 0");
				return null;
			}
			camera_xyz0 = intersceneLma.getSceneXYZ(false); // true for initial values
			camera_atr0 = intersceneLma.getSceneATR(false); // true for initial values
			double [] diffs_atr = intersceneLma.getV3Diff(ErsCorrection.DP_DSAZ);
			double [] diffs_xyz = intersceneLma.getV3Diff(ErsCorrection.DP_DSX);
			if ((diffs_atr[0] < clt_parameters.imp.exit_change_atr) &&
					(diffs_xyz[0] < clt_parameters.imp.exit_change_xyz)) {
				break;
			}
		}
		if (show_corr_fpn && (debug_level > -1)) { // now not needed, restore if needed
			String [] fpn_dbg_titles = new String[2 + scene_QuadClt.getNumSensors() * 2];
			fpn_dbg_titles[00] = "X_avg";
			fpn_dbg_titles[1] = "X_avg";
			for (int i = 0; i < scene_QuadClt.getNumSensors(); i++) {
				fpn_dbg_titles[2 + 2*i] = "X-"+i;
				fpn_dbg_titles[3 + 2*i] = "Y-"+i;
			}
			float [][] dbg_corr_fpn = new float [fpn_dbg_titles.length][];
			coord_motion = interCorrPair( // new double [tilesY][tilesX][][];
					clt_parameters,      // CLTParameters  clt_parameters,
					reference_QuadClt,   // QuadCLT reference_QuadCLT,
					ref_disparity,       // double []        ref_disparity, // null or alternative reference disparity  
					scene_QuadClt,       // QuadCLT scene_QuadCLT,
					camera_xyz0,         // xyz
					camera_atr0,         // pose[1], // atr
					null, // reliable_ref,        // final boolean [] selection, // may be null, if not null do not  process unselected tiles
					margin,              // final int        margin,
					sensor_mask_inter,   // final int        sensor_mask_inter, // The bitmask - which sensors to correlate, -1 - all.
					facc_2d_img,         // final float [][][]   accum_2d_corr, // if [1][][] - return accumulated 2d correlations (all pairs)final float [][][]   accum_2d_corr, // if [1][][] - return accumulated 2d correlations (all pairs)
					dbg_corr_fpn,        //	final float [][] dbg_corr_fpn,
					true,                // boolean            near_important, // do not reduce weight of the near tiles
					false,               // boolean            all_fpn,        // do not lower thresholds for non-fpn (used during search)
					2,                   // int                imp_debug_level,
					debug_level); // 1); // -1); // int debug_level);
			TileProcessor tp = reference_QuadClt.getTileProcessor();
			int tilesX = tp.getTilesX();
			int tilesY = tp.getTilesY();

			(new ShowDoubleFloatArrays()).showArrays(
					dbg_corr_fpn,
					tilesX,
					tilesY,
					true,
					scene_QuadClt.getImageName()+"-"+reference_QuadClt.getImageName()+"-CORR-FPN",
					fpn_dbg_titles);
		}
		if (show_corr_fpn) { // repeat after last adjustment to get images
			String [] fpn_dbg_titles = new String[2 + scene_QuadClt.getNumSensors() * 2];
			fpn_dbg_titles[00] = "X_avg";
			fpn_dbg_titles[1] = "X_avg";
			for (int i = 0; i < scene_QuadClt.getNumSensors(); i++) {
				fpn_dbg_titles[2 + 2*i] = "X-"+i;
				fpn_dbg_titles[3 + 2*i] = "Y-"+i;
			}
			float [][] dbg_corr_fpn = new float [fpn_dbg_titles.length][];
			coord_motion = interCorrPair( // new double [tilesY][tilesX][][];
					clt_parameters,      // CLTParameters  clt_parameters,
					reference_QuadClt,   // QuadCLT reference_QuadCLT,
					ref_disparity,       // double []        ref_disparity, // null or alternative reference disparity  
					scene_QuadClt,       // QuadCLT scene_QuadCLT,
					camera_xyz0,         // xyz
					camera_atr0,         // pose[1], // atr
					reliable_ref,        // final boolean [] selection, // may be null, if not null do not  process unselected tiles
					margin,              // final int        margin,
					sensor_mask_inter,   // final int        sensor_mask_inter, // The bitmask - which sensors to correlate, -1 - all.
					facc_2d_img,         // final float [][][]   accum_2d_corr, // if [1][][] - return accumulated 2d correlations (all pairs)final float [][][]   accum_2d_corr, // if [1][][] - return accumulated 2d correlations (all pairs)
					dbg_corr_fpn,        //	final float [][] dbg_corr_fpn,
					true,                // boolean            near_important, // do not reduce weight of the near tiles
					false,               // boolean            all_fpn,        // do not lower thresholds for non-fpn (used during search)
					2,                   // int                imp_debug_level,
					debug_level); // 1); // -1); // int debug_level);
			TileProcessor tp = reference_QuadClt.getTileProcessor();
			int tilesX = tp.getTilesX();
			int tilesY = tp.getTilesY();

			(new ShowDoubleFloatArrays()).showArrays(
					dbg_corr_fpn,
					tilesX,
					tilesY,
					true,
					scene_QuadClt.getImageName()+"-"+reference_QuadClt.getImageName()+"-CORR-FPN",
					fpn_dbg_titles);
		}
		
		
		
		if (clt_parameters.imp.debug_level > -2) {
			String [] lines1 = intersceneLma.printOldNew((clt_parameters.imp.debug_level > -1)); // false); // boolean allvectors)
			System.out.println("Adjusted interscene, iteration="+nlma+
					", last RMS = "+intersceneLma.getLastRms()[0]+
					" (pure RMS = "+intersceneLma.getLastRms()[1]+")"+
					", results:");
			for (String line : lines1) {
				System.out.println(line);
			}
//			if (intersceneLma.getLastRms()[0] > clt_parameters.imp.max_rms) {
//				System.out.println("Interscene adjustment failed, RMS="+
//						intersceneLma.getLastRms()[0]+" > "+ clt_parameters.imp.max_rms);
//				return null;
//			}
		}
		if ((rms_out != null) && (intersceneLma.getLastRms() != null)) {
			rms_out[0] = intersceneLma.getLastRms()[0];
			rms_out[1] = intersceneLma.getLastRms()[1];
			//if (lmaResult < 0) { last_rms[0]
		}
		if (max_rms > 0.0) {
			if (lmaResult < 0) { // = 0) {
				return null;
			}
			if (!(intersceneLma.getLastRms()[0] <= max_rms)) {
				System.out.println("RMS failed: "+intersceneLma.getLastRms()[0]+" >= " + max_rms);
				return null;
			}
		}
		return new double [][] {camera_xyz0, camera_atr0};
	}	

	
	/*
	public static String getChecksum(Serializable object) throws IOException, NoSuchAlgorithmException {
	    ByteArrayOutputStream baos = null;
	    ObjectOutputStream oos = null;
	    try {
	        baos = new ByteArrayOutputStream();
	        oos = new ObjectOutputStream(baos);
	        oos.writeObject(object);
	        MessageDigest md = MessageDigest.getInstance("MD5");
	        byte[] thedigest = md.digest(baos.toByteArray());
	        return DatatypeConverter.printHexBinary(thedigest);
	    } finally {
	        oos.close();
	        baos.close();
	    }
	}
	*/
	
	public double[][][]  test_LMA(
			CLTParameters  clt_parameters,			
			double k_prev,
			QuadCLT reference_QuadCLT,
			QuadCLT scene_QuadCLT,
			double corr_scale, //  = 0.75
			int debug_level)
	{
		TileProcessor tp = reference_QuadCLT.getTileProcessor();
		final int iscale = 8;
		boolean blur_reference = false;
		double ts =        reference_QuadCLT.getTimeStamp();
		double ts_prev =   ts;
		double [] camera_xyz0 = ZERO3.clone();
		double [] camera_atr0 = ZERO3.clone();
		
		ErsCorrection ersCorrection = reference_QuadCLT.getErsCorrection();
		String this_image_name = reference_QuadCLT.getImageName();
		
		System.out.println("\n"+this_image_name+":\n"+ersCorrection.extrinsic_corr.toString());
		System.out.println(String.format("%s: ers_wxyz_center=     %f, %f, %f", this_image_name,
				ersCorrection.ers_wxyz_center[0], ersCorrection.ers_wxyz_center[1],ersCorrection.ers_wxyz_center[2] ));
		System.out.println(String.format("%s: ers_wxyz_center_dt=  %f, %f, %f",	this_image_name,
				ersCorrection.ers_wxyz_center_dt[0], ersCorrection.ers_wxyz_center_dt[1],ersCorrection.ers_wxyz_center_dt[2] ));
		System.out.println(String.format("%s: ers_wxyz_center_d2t= %f, %f, %f", this_image_name,
				ersCorrection.ers_wxyz_center_d2t[0], ersCorrection.ers_wxyz_center_d2t[1],ersCorrection.ers_wxyz_center_d2t[2] ));
		System.out.println(String.format("%s: ers_watr_center_dt=  %f, %f, %f", this_image_name,
				ersCorrection.ers_watr_center_dt[0], ersCorrection.ers_watr_center_dt[1],ersCorrection.ers_watr_center_dt[2] ));
		System.out.println(String.format("%s: ers_watr_center_d2t= %f, %f, %f", this_image_name,
				ersCorrection.ers_watr_center_d2t[0], ersCorrection.ers_watr_center_d2t[1],ersCorrection.ers_watr_center_d2t[2] ));
		
		double dt = 0.0;
		if (scene_QuadCLT == null) {
			scene_QuadCLT = reference_QuadCLT;
		}
		if (scene_QuadCLT != null) {
			ts_prev = scene_QuadCLT.getTimeStamp();
			dt = ts-ts_prev;
			if (dt < 0) {
				k_prev = (1.0-k_prev);
			}
			if (Math.abs(dt) > 0.15) { // at least two frames TODO: use number of lines* line_time * ...? 
				k_prev = 0.5;
				System.out.println("Non-consecutive frames, dt = "+dt);
			}
			ErsCorrection ersCorrectionPrev = (ErsCorrection) (scene_QuadCLT.geometryCorrection);
			double [] wxyz_center_dt_prev =   ersCorrectionPrev.ers_wxyz_center_dt;
			double [] watr_center_dt_prev =   ersCorrectionPrev.ers_watr_center_dt;
			double [] wxyz_delta = new double[3];
			double [] watr_delta = new double[3];
			for (int i = 0; i <3; i++) {
				wxyz_delta[i] = corr_scale * dt * (k_prev * wxyz_center_dt_prev[i] + (1.0-k_prev) * ersCorrection.ers_wxyz_center_dt[i]);
				watr_delta[i] = corr_scale * dt * (k_prev * watr_center_dt_prev[i] + (1.0-k_prev) * ersCorrection.ers_watr_center_dt[i]);
			}
			camera_xyz0 = wxyz_delta;
			camera_atr0 = watr_delta;
		}
		
		int tilesX = tp.getTilesX();
		int tilesY = tp.getTilesY();
		if (debug_level > 0) {
			compareRefSceneTiles(
					"before_LMA",      // String suffix,
					blur_reference,    // boolean blur_reference,
					camera_xyz0,       // double [] camera_xyz0,
					camera_atr0,       // double [] camera_atr0,
					reference_QuadCLT, // QuadCLT reference_QuadCLT,
					scene_QuadCLT,     // QuadCLT scene_QuadCLT,
					iscale);           // int iscale) // 8
		}
		IntersceneLma intersceneLma = new IntersceneLma(
				this, // OpticalFlow opticalFlow
				clt_parameters.ilp.ilma_thread_invariant);
		for (int nlma = 0; nlma < clt_parameters.ilp.ilma_num_corr; nlma++) {
			boolean last_run = nlma == ( clt_parameters.ilp.ilma_num_corr - 1);
			int transform_size = tp.getTileSize();
			int macroTilesX =          tilesX/transform_size;
			int macroTilesY =          tilesY/transform_size;
			int macroTiles = macroTilesX * macroTilesY; 
			double [][] flowXY = new double [macroTiles][2]; // zero pre-shifts
			//		double [][] flowXY_frac = new double [macroTiles][]; // Will contain fractional X/Y shift for CLT
			double [][] reference_tiles_macro = new double [macroTiles][];
			double [][] vector_XYS = correlate2DIterate( // returns optical flow and confidence
					clt_parameters.img_dtt, // final ImageDttParameters  imgdtt_params,
					// for prepareSceneTiles()			
					camera_xyz0,                                 // final double []   scene_xyz,     // camera center in world coordinates
					camera_atr0,                                 // final double []   scene_atr,     // camera orientation relative to world frame
					scene_QuadCLT,                               // final QuadCLT     scene_QuadClt,
					reference_QuadCLT,                           // final QuadCLT     reference_QuadClt,
					reference_tiles_macro,                       //			final double [][] reference_tiles_macro,
					clt_parameters.ofp.center_occupancy_ref,         // final double      reference_center_occupancy,   // fraction of remaining  tiles in the center 8x8 area (<1.0)
					// flowXY should be initialized to all pairs of zeros (or deliberate pixel offset pairs if initial error is too high, will be modified with each iteration
					flowXY,                   // final double [][] flowXY, // per macro tile {mismatch in image pixels in X and Y directions // initialize to [reference_tiles.length][2]
					clt_parameters.ofp.tolerance_absolute_inter, // final double      tolerance_absolute, // absolute disparity half-range in each tile
					clt_parameters.ofp.tolerance_relative_inter, // final double      tolerance_relative, // relative disparity half-range in each tile
					clt_parameters.ofp.occupancy_inter,          // final double      occupancy,          // fraction of remaining  tiles (<1.0)
					clt_parameters.ofp.num_laplassian,           // final int         num_passes,
					clt_parameters.ofp.change_laplassian,        // final double      max_change,
					// for correlate2DSceneToReference ()
					clt_parameters.ofp.chn_weights,              // final double []   chn_weights, // absolute, starting from strength (strength,r,b,g)
					clt_parameters.ofp.corr_sigma,               // final double      corr_sigma,
					clt_parameters.ofp.fat_zero,                 //  final double      fat_zero,
					clt_parameters.ofp.late_normalize_iterate,   // final boolean     late_normalize,
					// for correlation2DToVectors_CM()
					clt_parameters.ofp.iradius_cm,               // final int         iradius,      // half-size of the square to process 
					clt_parameters.ofp.dradius_cm,               // final double      dradius,      // weight calculation (1/(r/dradius)^2 + 1)
					clt_parameters.ofp.refine_num_cm,            // final int         refine_num,   // number of iterations to apply weights around new center
					clt_parameters.ofp.num_refine_all,           // final int         num_run_all, // run all tiles for few iterations before filtering
					clt_parameters.ofp.max_refines,              // final int         max_tries,
					// for recalculateFlowXY()
					clt_parameters.ofp.magic_scale,              // final double      magic_scale, // 0.85 for CM
					clt_parameters.ofp.min_change,               // final double      min_change,
					clt_parameters.ofp.best_neibs_num,           // final int         best_num,
					clt_parameters.ofp.ref_stdev,                // final double      ref_stdev,
					clt_parameters.ofp.debug_level_iterate,      // final int         debug_level)
					clt_parameters.ofp.enable_debug_images);     //final boolean     enable_debug_images)

			if (debug_level > 2) {
				String dbg_title = "OpticalFlow-"+scene_QuadCLT.getImageName()+"-"+reference_QuadCLT.getImageName()+"-iteration_"+nlma;;
				showVectorXYConfidence(
						dbg_title, // String      title,
						vector_XYS, // double [][] flowXYS,
						macroTilesX); // int         width)	
			}
			int n = removeOutliers(
					clt_parameters.ofp.nsigma, // double nsigma, 1.5 - 2.0
					vector_XYS); // double [][] flowXYS)
			if (debug_level > -1) {
				System.out.println("Removed "+n+" outliers");
			}
			int n2 = removeOutliers(
					clt_parameters.ofp.nsigma2, // double nsigma, 1.5 - 2.0
					vector_XYS); // double [][] flowXYS)
			if (debug_level > -1) {
				System.out.println("Removed "+n2+" outliers in a second pass, total removed:"+(n+n2));
			}
			if (debug_level > 1) {
				String dbg_title = "OpticalFlowFiltered-"+scene_QuadCLT.getImageName()+"-"+reference_QuadCLT.getImageName()+"-iteration_"+nlma;
				showVectorXYConfidence(
						dbg_title, // String      title,
						vector_XYS, // double [][] flowXYS,
						macroTilesX); // int         width)	
			}

//			IntersceneLma intersceneLma = new IntersceneLma(
//					this); // OpticalFlow opticalFlow

			intersceneLma.prepareLMA(
					camera_xyz0,                                    // final double []   scene_xyz0,     // camera center in world coordinates (or null to use instance)
					camera_atr0,                                    // final double []   scene_atr0,     // camera orientation relative to world frame (or null to use instance)
					// reference atr, xyz are considered 0.0
					scene_QuadCLT,                                  // final QuadCLT     scene_QuadClt,
					reference_QuadCLT,                              //final QuadCLT     reference_QuadClt,
					clt_parameters.ilp.ilma_lma_select,             // final boolean[]   param_select,
					clt_parameters.ilp.ilma_regularization_weights, //  final double []   param_regweights,
					vector_XYS,                                     // final double [][] vector_XYS, // optical flow X,Y, confidence obtained from the correlate2DIterate()
					reference_tiles_macro,                          // final double [][] centers,    // macrotile centers (in pixels and average disparities
					(nlma == 0),                                    // boolean           first_run,
					debug_level);                                   // final int         debug_level)
			int lmaResult = intersceneLma.runLma(
					clt_parameters.ilp.ilma_lambda, // double lambda,           // 0.1
					clt_parameters.ilp.ilma_lambda_scale_good, //  double lambda_scale_good,// 0.5
					clt_parameters.ilp.ilma_lambda_scale_bad,  // double lambda_scale_bad, // 8.0
					clt_parameters.ilp.ilma_lambda_max,        // double lambda_max,       // 100
					clt_parameters.ilp.ilma_rms_diff,          // double rms_diff,         // 0.001
					clt_parameters.ilp.ilma_num_iter,          // int    num_iter,         // 20
					last_run,                                  // boolean last_run,
					clt_parameters.ilp.ilma_debug_level);      // int    debug_level)
			if (lmaResult < 0) {
				System.out.println("LMA failed");
				break;
			}
			camera_xyz0 = intersceneLma.getSceneXYZ(false); // true for initial values
			camera_atr0 = intersceneLma.getSceneATR(false); // true for initial values
			if (debug_level > 1) {
				compareRefSceneTiles(
						"iteration_"+nlma,      // String suffix,
						blur_reference,    // boolean blur_reference,
						camera_xyz0,       // double [] camera_xyz0,
						camera_atr0,       // double [] camera_atr0,
						reference_QuadCLT, // QuadCLT reference_QuadCLT,
						scene_QuadCLT,     // QuadCLT scene_QuadCLT,
						iscale);           // int iscale) // 8
			}
			if (lmaResult <= 1) {
				break;
			}
		}
		if (debug_level == 1)  {
			compareRefSceneTiles(
					"after_lma",      // String suffix,
					blur_reference,    // boolean blur_reference,
					camera_xyz0,       // double [] camera_xyz0,
					camera_atr0,       // double [] camera_atr0,
					reference_QuadCLT, // QuadCLT reference_QuadCLT,
					scene_QuadCLT,     // QuadCLT scene_QuadCLT,
					iscale);           // int iscale) // 8
		}
		reference_QuadCLT.getErsCorrection().addScene(scene_QuadCLT.getImageName(), camera_xyz0,camera_atr0);
		reference_QuadCLT.saveInterProperties( // save properties for interscene processing (extrinsics, ers, ...)
	            null, // String path,             // full name with extension or w/o path to use x3d directory
	            debug_level);
		return null;
	}
	
//		ErsCorrection ersCorrection = reference_QuadCLT.getErsCorrection();
	
	/**
	 * A top-level method for testing optical flow generating, currently includes temporary testing functionality 
	 * @param clt_parameters CLT parameters
	 * @param k_prev Coefficient of the previous (in time) frame weight to calculate initial estimation of the pose
	 *        differences from the single-scene ERS values determined from the Lazy Eye LMA adjustment. The ERS
	 *        parameters typically correspond to the second half of the image (top is usually inifinity/long range,
	 *        while the scene pose is calculated for the image center scanline.  Tested with k_prev = 0.75
	 * @param reference_QuadClt Reference QuadCLT instance.
	 * @param scene_QuadClt Scene QuadCLT instance.
	 * @param corr_scale Correction coefficient - still to find out the reason that the pose difference predicted
	 *        from the intrascene ERS should be reduced when calculating interscene pose difference. The heuristc
	 *        value is 0.75.
	 * @param debug_level Debug Level
	 * @return a pair of reference and interpolated scenes
	 */
	public double[][][]  get_pair(
			CLTParameters  clt_parameters,			
			double k_prev,
			QuadCLT reference_QuadCLT,
			QuadCLT scene_QuadCLT,
			double corr_scale, //  = 0.75 - REMOVE
			int debug_level)
	{
		TileProcessor tp = reference_QuadCLT.getTileProcessor();
		final int iscale = 8;
		double ts =        reference_QuadCLT.getTimeStamp();
		double ts_prev =   ts;
		double [] camera_xyz0 = ZERO3.clone();
		double [] camera_atr0 = ZERO3.clone();
		
		ErsCorrection ersCorrection = reference_QuadCLT.getErsCorrection();
		String this_image_name = reference_QuadCLT.getImageName();
		
		System.out.println("\n"+this_image_name+":\n"+ersCorrection.extrinsic_corr.toString());
		System.out.println(String.format("%s: ers_wxyz_center=     %f, %f, %f", this_image_name,
				ersCorrection.ers_wxyz_center[0], ersCorrection.ers_wxyz_center[1],ersCorrection.ers_wxyz_center[2] ));
		System.out.println(String.format("%s: ers_wxyz_center_dt=  %f, %f, %f",	this_image_name,
				ersCorrection.ers_wxyz_center_dt[0], ersCorrection.ers_wxyz_center_dt[1],ersCorrection.ers_wxyz_center_dt[2] ));
		System.out.println(String.format("%s: ers_wxyz_center_d2t= %f, %f, %f", this_image_name,
				ersCorrection.ers_wxyz_center_d2t[0], ersCorrection.ers_wxyz_center_d2t[1],ersCorrection.ers_wxyz_center_d2t[2] ));
		System.out.println(String.format("%s: ers_watr_center_dt=  %f, %f, %f", this_image_name,
				ersCorrection.ers_watr_center_dt[0], ersCorrection.ers_watr_center_dt[1],ersCorrection.ers_watr_center_dt[2] ));
		System.out.println(String.format("%s: ers_watr_center_d2t= %f, %f, %f", this_image_name,
				ersCorrection.ers_watr_center_d2t[0], ersCorrection.ers_watr_center_d2t[1],ersCorrection.ers_watr_center_d2t[2] ));
		
		double dt = 0.0;
		if (scene_QuadCLT == null) {
			scene_QuadCLT = reference_QuadCLT;
		}
		if (scene_QuadCLT != null) {
			ts_prev = scene_QuadCLT.getTimeStamp();
			dt = ts-ts_prev;
			if (dt < 0) {
				k_prev = (1.0-k_prev);
			}
			if (Math.abs(dt) > 0.15) { // at least two frames TODO: use number of lines* line_time * ...? 
				k_prev = 0.5;
				System.out.println("Non-consecutive frames, dt = "+dt);
			}
			ErsCorrection ersCorrectionPrev = (ErsCorrection) (scene_QuadCLT.geometryCorrection);
			double [] wxyz_center_dt_prev =   ersCorrectionPrev.ers_wxyz_center_dt;
			double [] watr_center_dt_prev =   ersCorrectionPrev.ers_watr_center_dt;
			double [] wxyz_delta = new double[3];
			double [] watr_delta = new double[3];
			for (int i = 0; i <3; i++) {
				wxyz_delta[i] = corr_scale * dt * (k_prev * wxyz_center_dt_prev[i] + (1.0-k_prev) * ersCorrection.ers_wxyz_center_dt[i]);
				watr_delta[i] = corr_scale * dt * (k_prev * watr_center_dt_prev[i] + (1.0-k_prev) * ersCorrection.ers_watr_center_dt[i]);
			}
			camera_xyz0 = wxyz_delta;
			camera_atr0 = watr_delta;
		}
		
		int tilesX = tp.getTilesX();
		int tilesY = tp.getTilesY();
		String [] dsrbg_titles = {"d", "s", "r", "b", "g"};
		String title = this_image_name+"-"+scene_QuadCLT.image_name+"-dt"+dt;
		double [][] dsrbg = transformCameraVew( // shifts previous image correctly (right)
				title,                   // final String    title,
				null, // final double [][] dsrbg_camera_in,
				camera_xyz0,             // double [] camera_xyz, // camera center in world coordinates
				camera_atr0,             //double [] camera_atr, // camera orientation relative to world frame
				scene_QuadCLT,           // QuadCLT   camera_QuadClt,
				reference_QuadCLT,       // reference
				iscale);
		double [][][] pair = {reference_QuadCLT.getDSRBG(),dsrbg};
		
		/*
		reference_QuadCLT.getErsCorrection().compareDSItoWorldDerivatives(
				reference_QuadCLT, // QuadCLT   scene_QuadClt,
				0.03,              // double    max_inf_disparity, // absolute value
				1);                // int       debug_level);
		*/
		reference_QuadCLT.getErsCorrection().comparePXYD_Derivatives(
				scene_QuadCLT,     // QuadCLT   scene_QuadClt,
				reference_QuadCLT, // QuadCLT   reference_QuadClt,
				0.03, // double    max_inf_disparity, // absolute value
				1); // int       debug_level
		
		
		if (debug_level > -100) {
			return pair;
		}
		
		// combine this scene with warped previous one
		if (debug_level > -2) {
			String [] rtitles = new String[2* dsrbg_titles.length];
			double [][] dbg_rslt = new double [rtitles.length][];
			for (int i = 0; i < dsrbg_titles.length; i++) {
				rtitles[2*i] =    dsrbg_titles[i]+"0";
				rtitles[2*i+1] =  dsrbg_titles[i];
				dbg_rslt[2*i] =   pair[0][i];
				dbg_rslt[2*i+1] = pair[1][i];
			}
//			String title = this_image_name+"-"+scene_QuadCLT.image_name+"-dt"+dt;
			(new ShowDoubleFloatArrays()).showArrays(
					dbg_rslt,
					tilesX,
					tilesY,
					true,
					title,
					rtitles);
		}
///		double      tolerance_absolute = 0.25; // absolute disparity half-range in each tile
///		double      tolerance_relative = 0.2; // relative disparity half-range in each tile
///		double      center_occupancy =   0.25;   // fraction of remaining  tiles in the center 8x8 area (<1.0)
///		int         num_passes = 100;
///		double      max_change = 0.005 ;

///		double      tolerance_absolute_inter = 0.25; // absolute disparity half-range in each tile
///		double      tolerance_relative_inter = 0.2; // relative disparity half-range in each tile
///		double      occupancy_inter =         0.25;   // fraction of remaining  tiles in the center 8x8 area (<1.0)

		// Add limitations on disparity ? To be used with multi-tile consolidation
		// Check with walking, not only rotating
		int transform_size = tp.getTileSize();
		int macroTilesX =          tilesX/transform_size;
		int macroTilesY =          tilesY/transform_size;
		int macroTiles = macroTilesX * macroTilesY; 
		double [][] flowXY = new double [macroTiles][2]; // zero pre-shifts
		double [][] flowXY_frac = new double [macroTiles][]; // Will contain fractional X/Y shift for CLT
///		double []   chn_weights = {1.0,1.0,1.0,1.0}; // strength, r,b,g
//		double []   chn_weights = {1.0,0.0,0.0,0.0}; // strength, r,b,g
//		double []   chn_weights = {0.0,1.0,1.0,1.0}; // strength, r,b,g
		// Apply DOG to colors, normalize by standard deviation?
///		double      corr_sigma = 0.5;
///		double      fat_zero =   0.05;
///		double      frac_radius = 0.9;  // add to integer radius for window calculation
///		double      tolerance_absolute_inter_macro = 0.25; // absolute disparity half-range to consolidate macro tiles
///		double      tolerance_relative_inter_macro = 0.2;  // relative disparity half-range to consolidate macro tiles
///		int         iradius =    3;      // half-size of the square to process 
///		double      dradius =    1.5;      // weight calculation (1/(r/dradius)^2 + 1)
///		int         refine_num_cm = 5;   // number of iterations to apply weights around new center
///		int         num_refine_all = 3;
///		int         max_refines =   50;
////		int         max_rad = 3;
		
///		int         best_num = 4; // use  4 best neighbors to calculate std deviation
///		double      ref_stdev = 5.0; // strength 0.5 if standard deviation of best neighbors to tile difference is this.
		
///		boolean     combine_empty_only = true; // false;
///		double      magic_scale = 0.85; // 2.0 * 0.85;

///		boolean late_normalize_iterate = true;
		
		// for recalculateFlowXY()
///		double      min_change = 0.1; // 01;//   sqrt (dx*dx + dy*dy) for correction (int tiles) in pixels
		
///		int         debug_level_iterate = -1; // 2;

		
		double [][] vector_XYS = correlate2DIterate( // returns optical flow and confidence
				clt_parameters.img_dtt, // final ImageDttParameters  imgdtt_params,
				// for prepareSceneTiles()			
				camera_xyz0,                                 // final double []   scene_xyz,     // camera center in world coordinates
				camera_atr0,                                 // final double []   scene_atr,     // camera orientation relative to world frame
				scene_QuadCLT,                               // final QuadCLT     scene_QuadClt,
				reference_QuadCLT,                           // final QuadCLT     reference_QuadClt,
				null,                                        // final double [][] reference_tiles_macro,
				clt_parameters.ofp.center_occupancy_ref,         // final double      reference_center_occupancy,   // fraction of remaining  tiles in the center 8x8 area (<1.0)
				// flowXY should be initialized to all pairs of zeros (or deliberate pixel offset pairs if initial error is too high, will be modified with each iteration
				flowXY,                   // final double [][] flowXY, // per macro tile {mismatch in image pixels in X and Y directions // initialize to [reference_tiles.length][2]
				clt_parameters.ofp.tolerance_absolute_inter, // final double      tolerance_absolute, // absolute disparity half-range in each tile
				clt_parameters.ofp.tolerance_relative_inter, // final double      tolerance_relative, // relative disparity half-range in each tile
				clt_parameters.ofp.occupancy_inter,          // final double      occupancy,          // fraction of remaining  tiles (<1.0)
				clt_parameters.ofp.num_laplassian,           // final int         num_passes,
				clt_parameters.ofp.change_laplassian,        // final double      max_change,
				// for correlate2DSceneToReference ()
				clt_parameters.ofp.chn_weights,              // final double []   chn_weights, // absolute, starting from strength (strength,r,b,g)
				clt_parameters.ofp.corr_sigma,               // final double      corr_sigma,
				clt_parameters.ofp.fat_zero,                 //  final double      fat_zero,
				clt_parameters.ofp.late_normalize_iterate,   // final boolean     late_normalize,
				// for correlation2DToVectors_CM()
				clt_parameters.ofp.iradius_cm,               // final int         iradius,      // half-size of the square to process 
				clt_parameters.ofp.dradius_cm,               // final double      dradius,      // weight calculation (1/(r/dradius)^2 + 1)
				clt_parameters.ofp.refine_num_cm,            // final int         refine_num,   // number of iterations to apply weights around new center
				clt_parameters.ofp.num_refine_all,           // final int         num_run_all, // run all tiles for few iterations before filtering
				clt_parameters.ofp.max_refines,              // final int         max_tries,
				// for recalculateFlowXY()
				clt_parameters.ofp.magic_scale,              // final double      magic_scale, // 0.85 for CM
				clt_parameters.ofp.min_change,               // final double      min_change,
				clt_parameters.ofp.best_neibs_num,           // final int         best_num,
				clt_parameters.ofp.ref_stdev,                // final double      ref_stdev,
				clt_parameters.ofp.debug_level_iterate,      // final int         debug_level)
				clt_parameters.ofp.enable_debug_images);     //final boolean     enable_debug_images)
		
		if (debug_level > -2) {
			String dbg_title = "OpticalFlow-"+scene_QuadCLT.getImageName()+"-"+reference_QuadCLT.getImageName();
			showVectorXYConfidence(
					dbg_title, // String      title,
					vector_XYS, // double [][] flowXYS,
					macroTilesX); // int         width)	
		}
		
		
		if (debug_level > 0) {

			double [][][] reference_tiles = prepareReferenceTiles(
					reference_QuadCLT,        // final QuadCLT     qthis,
					clt_parameters.ofp.tolerance_absolute_ref, // final double      tolerance_absolute, // absolute disparity half-range in each tile
					clt_parameters.ofp.tolerance_relative_ref, // final double      tolerance_relative, // relative disparity half-range in each tile
					clt_parameters.ofp.center_occupancy_ref,   // final double      center_occupancy,   // fraction of remaining  tiles in the center 8x8 area (<1.0)
					-1); // -1); // 2); // final int         debug_level)

			fillTilesNans(
					reference_tiles,          // final double [][][] nan_tiles,
					reference_QuadCLT,                 // final QuadCLT     qthis,
					clt_parameters.ofp.num_laplassian, //   num_passes,            // final int         num_passes,
					clt_parameters.ofp.change_laplassian, // max_change,            // final double      max_change,
					-1); //-1); // 2);                    // final int         debug_level)


			double [][][] scene_tiles = prepareSceneTiles(// to match to reference
					// null for {scene,reference}{xyz,atr} uses instances globals 
					camera_xyz0,              // final double []   scene_xyz,     // camera center in world coordinates
					camera_atr0,              // final double []   scene_atr,     // camera orientation relative to world frame
					scene_QuadCLT,                    // final QuadCLT     scene_QuadClt,
					reference_QuadCLT,                    // final QuadCLT     reference_QuadClt,
					reference_tiles,          // final double [][][] reference_tiles, // prepared with prepareReferenceTiles() + fillTilesNans();
					flowXY,                   // final double [][] flowXY, // per macro tile {mismatch in image pixels in X and Y directions
					flowXY_frac,              // final double [][] flowXY_frac, // should be initialized as [number of macro tiles][] - returns fractional shifts [-0.5, 0.5)
					clt_parameters.ofp.tolerance_absolute_inter, // final double      tolerance_absolute, // absolute disparity half-range in each tile
					clt_parameters.ofp.tolerance_relative_inter, // final double      tolerance_relative, // relative disparity half-range in each tile
					clt_parameters.ofp.occupancy_inter,          // final double      occupancy,          // fraction of remaining  tiles (<1.0)
					clt_parameters.ofp.num_laplassian, //   num_passes,            // final int         num_passes,
					clt_parameters.ofp.change_laplassian, // max_change,            // final double      max_change,
					-1); //-1); // 1); // 2);                       // final int         debug_level)

			String dbg_title = "flowXY_frac-"+scene_QuadCLT.getImageName()+"-"+reference_QuadCLT.getImageName();
			String [] dbg_titles = {"dpX", "dpY"};
			double [][] dbg_img = new double [dbg_titles.length][macroTilesX*macroTilesY];
				Arrays.fill(dbg_img[0], Double.NaN);
				Arrays.fill(dbg_img[1], Double.NaN);
				for (int i = 0; i < flowXY_frac.length; i++) if (flowXY_frac[i] != null){
					dbg_img[0][i] = flowXY_frac[i][0];
					dbg_img[1][i] = flowXY_frac[i][1];
				}
			(new ShowDoubleFloatArrays()).showArrays(
					dbg_img,
					macroTilesX,
					macroTilesY,
					true,
					dbg_title,
					dbg_titles);

			double [][][] corr2dscene_ref_multi = new double [clt_parameters.ofp.test_corr_rad_max + 2][][]; 
			corr2dscene_ref_multi[0] = correlate2DSceneToReference(// to match to reference
					clt_parameters.img_dtt, // final ImageDttParameters  imgdtt_params,
					scene_QuadCLT,          // final QuadCLT     scene_QuadClt,
					reference_QuadCLT,      // final QuadCLT     reference_QuadClt,
					scene_tiles,            // final double [][][] scene_tiles,     // prepared with prepareSceneTiles()
					reference_tiles,        // final double [][][] reference_tiles, // prepared with prepareReferenceTiles() + fillTilesNans(); - combine?
					flowXY_frac,            // final double [][] flowXY_frac, // X, YH fractional shift [-0.5,0.5) to implement with FD rotations
					clt_parameters.ofp.chn_weights,            // final double []   chn_weights, // absolute, starting from strength (strength,r,b,g)
					clt_parameters.ofp.corr_sigma,             // final double      corr_sigma,
					clt_parameters.ofp.fat_zero,               //  final double      fat_zero,
					false,                  // final boolean     late_normalize,
					clt_parameters.ofp.combine_empty_only,     //  final boolean     combine_empty_only, // only use neighbor correlations for empty tiles (false - any)
					0.0,     // final double      combine_dradius
					clt_parameters.ofp.tolerance_absolute_macro, // final double      tolerance_absolute, // absolute disparity half-range to consolidate tiles
					clt_parameters.ofp.tolerance_relative_macro, // final double      tolerance_relative, // relative disparity half-range to consolidate tiles
					-1); // 1); // final int         debug_level)

			for (int irad = 0; irad <= 3; irad++) {
				corr2dscene_ref_multi[irad+1]= correlate2DSceneToReference(// to match to reference
						clt_parameters.img_dtt, // final ImageDttParameters  imgdtt_params,
						scene_QuadCLT,          // final QuadCLT     scene_QuadClt,
						reference_QuadCLT,      // final QuadCLT     reference_QuadClt,
						scene_tiles,            // final double [][][] scene_tiles,     // prepared with prepareSceneTiles()
						reference_tiles,        // final double [][][] reference_tiles, // prepared with prepareReferenceTiles() + fillTilesNans(); - combine?
						flowXY_frac,            // final double [][] flowXY_frac, // X, YH fractional shift [-0.5,0.5) to implement with FD rotations
						clt_parameters.ofp.chn_weights,            // final double []   chn_weights, // absolute, starting from strength (strength,r,b,g)
						clt_parameters.ofp.corr_sigma,             // final double      corr_sigma,
						clt_parameters.ofp.fat_zero,               //  final double      fat_zero,
						true,                   // final boolean     late_normalize,
						clt_parameters.ofp.combine_empty_only,     //  final boolean     combine_empty_only, // only use neighbor correlations for empty tiles (false - any)
						irad + clt_parameters.ofp.frac_radius,     // final double      combine_dradius
						clt_parameters.ofp.tolerance_absolute_macro, // final double      tolerance_absolute, // absolute disparity half-range to consolidate tiles
						clt_parameters.ofp.tolerance_relative_macro, // final double      tolerance_relative, // relative disparity half-range to consolidate tiles
						-1); // 1); // final int         debug_level)
			}

			showCorrTiles(
					"scene:"+scene_QuadCLT.getImageName()+"-ref"+reference_QuadCLT.getImageName(),  //  String title,
					corr2dscene_ref_multi,           // double [][] source_tiles,
					tilesX/transform_size,     // int         tilesX,
					(2 * transform_size - 1),  // int         tile_width,
					(2 * transform_size - 1)); // int         tile_height) // extra margins over 16x16 tiles to accommodate distorted destination tiles
			//reference_tiles
			double [][][][] scene_to_ref = {reference_tiles, scene_tiles};
			if (debug_level > 0) {
				showCompareMacroTiles(
						"tiles_scene-"+scene_QuadCLT.getImageName()+"-ref"+reference_QuadCLT.getImageName(),// String title,
						scene_to_ref,      // double [][][][] source_tiles_sets,
						reference_QuadCLT, // final QuadCLT qthis,
						0);                // final int     margin) // extra margins over 16x16 tiles to accommodate distorted destination tiles
			}


			if (debug_level > 100) {
				String flowXYS_title =  (debug_level > 0)?("vectorXYS_"+scene_QuadCLT.getImageName()+"-ref"+reference_QuadCLT.getImageName()):null;
				//			double [][] vectorXYConfidence =  
				attachVectorConfidence(
						flowXY,         // final double [][] flowXY,
						macroTilesX,    // final int         width,
						clt_parameters.ofp.best_neibs_num,       // final int         best_num,
						clt_parameters.ofp.ref_stdev,      // final double      ref_stdev,
						flowXYS_title); // final String      debug_title);    

				double [][][] vectorsXYS = new double [corr2dscene_ref_multi.length][][];
				for (int i = 0; i < vectorsXYS.length; i++) {
					vectorsXYS[i] = correlation2DToVectors_CM(
							corr2dscene_ref_multi[i], // final double [][] corr2d_tiles, // per 2d calibration tiles (or nulls)
							transform_size, // final int         transform_size,
							clt_parameters.ofp.iradius_cm, // final int         iradius,      // half-size of the square to process 
							clt_parameters.ofp.dradius_cm,      // final double      dradius,      // weight calculation (1/(r/dradius)^2 + 1)
							clt_parameters.ofp.refine_num_cm,   // final int         refine_num,   // number of iterations to apply weights around new center
							1); //final int         debug_level)
				}
				if (debug_level > -1) {
					showVectorXYConfidence(
							"dXdYS-"+scene_QuadCLT.getImageName()+"-"+reference_QuadCLT.getImageName(), // String        title,
							vectorsXYS, // double [][][] flowXYSs,
							macroTilesX); // int           width)
				}
				
				int selected_index = 1; // single, post-norm
				double [][] flowXY1 = recalculateFlowXY(
						flowXY, // final double [][] currentFlowXY,
						vectorsXYS[selected_index], // final double [][] corr_vectorsXY,
						clt_parameters.ofp.magic_scale/transform_size); // final double      magic_scale) // 0.85 for CM

				double [][][] scene_tiles1 = prepareSceneTiles(// to match to reference
						// null for {scene,reference}{xyz,atr} uses instances globals 
						camera_xyz0,              // final double []   scene_xyz,     // camera center in world coordinates
						camera_atr0,              // final double []   scene_atr,     // camera orientation relative to world frame
						scene_QuadCLT,                    // final QuadCLT     scene_QuadClt,
						reference_QuadCLT,                    // final QuadCLT     reference_QuadClt,
						reference_tiles,          // final double [][][] reference_tiles, // prepared with prepareReferenceTiles() + fillTilesNans();
						flowXY1,                  // final double [][] flowXY, // per macro tile {mismatch in image pixels in X and Y directions
						flowXY_frac,              // final double [][] flowXY_frac, // should be initialized as [number of macro tiles][] - returns fractional shifts [-0.5, 0.5)
						clt_parameters.ofp.tolerance_absolute_inter, // final double      tolerance_absolute, // absolute disparity half-range in each tile
						clt_parameters.ofp.tolerance_relative_inter, // final double      tolerance_relative, // relative disparity half-range in each tile
						clt_parameters.ofp.occupancy_inter,          // final double      occupancy,          // fraction of remaining  tiles (<1.0)
						clt_parameters.ofp.num_laplassian, //   num_passes,            // final int         num_passes,
						clt_parameters.ofp.change_laplassian, // max_change,            // final double      max_change,
						1); //-1); // 1); // 2);                       // final int         debug_level)

				// single, late
				corr2dscene_ref_multi[0] = correlate2DSceneToReference(// to match to reference
						clt_parameters.img_dtt, // final ImageDttParameters  imgdtt_params,
						scene_QuadCLT,          // final QuadCLT     scene_QuadClt,
						reference_QuadCLT,      // final QuadCLT     reference_QuadClt,
						scene_tiles1,            // final double [][][] scene_tiles,     // prepared with prepareSceneTiles()
						reference_tiles,        // final double [][][] reference_tiles, // prepared with prepareReferenceTiles() + fillTilesNans(); - combine?
						flowXY_frac,            // final double [][] flowXY_frac, // X, YH fractional shift [-0.5,0.5) to implement with FD rotations
						clt_parameters.ofp.chn_weights,            // final double []   chn_weights, // absolute, starting from strength (strength,r,b,g)
						clt_parameters.ofp.corr_sigma,             // final double      corr_sigma,
						clt_parameters.ofp.fat_zero,               //  final double      fat_zero,
						false,                  // final boolean     late_normalize,
						clt_parameters.ofp.combine_empty_only,     //  final boolean     combine_empty_only, // only use neighbor correlations for empty tiles (false - any)
						0.0,                    // final double      combine_dradius
						clt_parameters.ofp.tolerance_absolute_macro, // final double      tolerance_absolute, // absolute disparity half-range to consolidate tiles
						clt_parameters.ofp.tolerance_relative_macro, // final double      tolerance_relative, // relative disparity half-range to consolidate tiles
						-1); // 1); // final int         debug_level)

				for (int irad = 0; irad <= 3; irad++) {
					corr2dscene_ref_multi[irad+1]= correlate2DSceneToReference(// to match to reference
							clt_parameters.img_dtt, // final ImageDttParameters  imgdtt_params,
							scene_QuadCLT,          // final QuadCLT     scene_QuadClt,
							reference_QuadCLT,      // final QuadCLT     reference_QuadClt,
							scene_tiles1,           // final double [][][] scene_tiles,     // prepared with prepareSceneTiles()
							reference_tiles,        // final double [][][] reference_tiles, // prepared with prepareReferenceTiles() + fillTilesNans(); - combine?
							flowXY_frac,            // final double [][] flowXY_frac, // X, YH fractional shift [-0.5,0.5) to implement with FD rotations
							clt_parameters.ofp.chn_weights,            // final double []   chn_weights, // absolute, starting from strength (strength,r,b,g)
							clt_parameters.ofp.corr_sigma,             // final double      corr_sigma,
							clt_parameters.ofp.fat_zero,               //  final double      fat_zero,
							true, // final boolean     late_normalize,
							clt_parameters.ofp.combine_empty_only,     //  final boolean     combine_empty_only, // only use neighbor correlations for empty tiles (false - any)
							irad + clt_parameters.ofp.frac_radius,     // final double      combine_dradius
							clt_parameters.ofp.tolerance_absolute_macro, // final double      tolerance_absolute, // absolute disparity half-range to consolidate tiles
							clt_parameters.ofp.tolerance_relative_macro, // final double      tolerance_relative, // relative disparity half-range to consolidate tiles
							-1); // 1); // final int         debug_level)
				}


				showCorrTiles(
						"scene:"+scene_QuadCLT.getImageName()+"-ref"+reference_QuadCLT.getImageName(),  //  String title,
						corr2dscene_ref_multi,           // double [][] source_tiles,
						tilesX/transform_size,     // int         tilesX,
						(2 * transform_size - 1),  // int         tile_width,
						(2 * transform_size - 1)); // int         tile_height) // extra margins over 16x16 tiles to accommodate distorted destination tiles

				double [][][] vectorsXYS1 = new double [corr2dscene_ref_multi.length][][];
				for (int i = 0; i < vectorsXYS1.length; i++) {
					vectorsXYS1[i] = correlation2DToVectors_CM(
							corr2dscene_ref_multi[i], // final double [][] corr2d_tiles, // per 2d calibration tiles (or nulls)
							transform_size, // final int         transform_size,
							clt_parameters.ofp.iradius_cm, // final int         iradius,      // half-size of the square to process 
							clt_parameters.ofp.dradius_cm,      // final double      dradius,      // weight calculation (1/(r/dradius)^2 + 1)
							clt_parameters.ofp.refine_num_cm,   // final int         refine_num,   // number of iterations to apply weights around new center
							1); //final int         debug_level)
				}
				if (debug_level > -1) {
					showVectorXYConfidence(
							"dXdYS1-"+scene_QuadCLT.getImageName()+"-"+reference_QuadCLT.getImageName(), // String        title,
							vectorsXYS1, // double [][][] flowXYSs,
							macroTilesX); // int           width)
				}
				double [][] flowXY2 = recalculateFlowXY(
						flowXY1, // final double [][] currentFlowXY,
						vectorsXYS1[selected_index], // final double [][] corr_vectorsXY,
						clt_parameters.ofp.magic_scale/transform_size); // final double      magic_scale) // 0.85 for CM

				//			double [][][] scene_tiles2 = 
				prepareSceneTiles(// to match to reference
						// null for {scene,reference}{xyz,atr} uses instances globals 
						camera_xyz0,              // final double []   scene_xyz,     // camera center in world coordinates
						camera_atr0,              // final double []   scene_atr,     // camera orientation relative to world frame
						scene_QuadCLT,                    // final QuadCLT     scene_QuadClt,
						reference_QuadCLT,                    // final QuadCLT     reference_QuadClt,
						reference_tiles,          // final double [][][] reference_tiles, // prepared with prepareReferenceTiles() + fillTilesNans();
						flowXY2,                  // final double [][] flowXY, // per macro tile {mismatch in image pixels in X and Y directions
						flowXY_frac,              // final double [][] flowXY_frac, // should be initialized as [number of macro tiles][] - returns fractional shifts [-0.5, 0.5)
						clt_parameters.ofp.tolerance_absolute_inter, // final double      tolerance_absolute, // absolute disparity half-range in each tile
						clt_parameters.ofp.tolerance_relative_inter, // final double      tolerance_relative, // relative disparity half-range in each tile
						clt_parameters.ofp.occupancy_inter,          // final double      occupancy,          // fraction of remaining  tiles (<1.0)
						clt_parameters.ofp.num_laplassian, //   num_passes,            // final int         num_passes,
						clt_parameters.ofp.change_laplassian, // max_change,            // final double      max_change,
						1); //-1); // 1); // 2);                       // final int         debug_level)
			}

		}	
		return pair;
	}
	

}
