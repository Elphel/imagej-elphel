package com.elphel.imagej.tileprocessor;
/**
 **
 ** QuadCLT - Process images with CLT-based methods (code specific to ImageJ plugin)
 ** Using CPU+GPU
 **
 ** Copyright (C) 2017-2020 Elphel, Inc.
 **
 ** -----------------------------------------------------------------------------**
 **
 **  QuadCLT.java is free software: you can redistribute it and/or modify
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

//import static jcuda.driver.JCudaDriver.cuMemcpyDtoH;

import java.awt.Rectangle;
import java.io.DataOutputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.channels.Channels;
import java.nio.channels.WritableByteChannel;
import java.util.Arrays;
import java.util.Properties;
import java.util.Random;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.atomic.DoubleAccumulator;

import com.elphel.imagej.cameras.CLTParameters;
import com.elphel.imagej.cameras.ColorProcParameters;
import com.elphel.imagej.cameras.EyesisCorrectionParameters;
import com.elphel.imagej.common.DoubleGaussianBlur;
import com.elphel.imagej.common.PolynomialApproximation;
import com.elphel.imagej.common.ShowDoubleFloatArrays;
import com.elphel.imagej.correction.CorrectionColorProc;
import com.elphel.imagej.correction.EyesisCorrections;
import com.elphel.imagej.gpu.GPUTileProcessor;
import com.elphel.imagej.gpu.GpuQuad;
import com.elphel.imagej.gpu.TpTask;
import com.elphel.imagej.tileprocessor.QuadCLTCPU.SetChannels;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.Prefs;
import ij.io.FileSaver;
public class QuadCLT extends QuadCLTCPU {
	int dbg_lev = 1;
	private GpuQuad gpuQuad =              null; 	// use updateQuadCLT() to update after switching to a different
	// QuadCLT instance, this.getGPU().getQuadCLT() should be equal to this
	/*
			if (getGPU().getQuadCLT() != this) {
				getGPU().updateQuadCLT(this); // to re-load new set of Bayer images to the GPU
			}

	 */
	
	public QuadCLT( // creates AUX, does not know numSensors?
			String                                          prefix,
			Properties                                      properties,
			EyesisCorrections                               eyesisCorrections,
			EyesisCorrectionParameters.CorrectionParameters correctionsParameters
			){
		super (prefix,
				properties,
				eyesisCorrections,
				correctionsParameters
				);
	}
	
	public GpuQuad getGPUQuad() {
		return gpuQuad;
	}
	
	public boolean hasGPU() {
		return (gpuQuad != null);
	}
	
	
	public QuadCLT saveQuadClt() {
		if (gpuQuad == null) {
			return null;
		}
		QuadCLT savedQuadClt =  gpuQuad.getQuadCLT();
		if (savedQuadClt != this) {
			gpuQuad.updateQuadCLT(this); // to re-load new set of Bayer images to the GPU
		} else {
			savedQuadClt = null;
		}
		return savedQuadClt;
	}
	void restoreQuadClt(QuadCLT savedQuadClt) {
		if (savedQuadClt !=  null) {
			gpuQuad.updateQuadCLT(savedQuadClt);
		}
	}

	public void setQuadClt() {
		if (gpuQuad != null) {
			gpuQuad.updateQuadCLT(this); // to re-load new set of Bayer images to the GPU
		}
	}
	
	public QuadCLT(QuadCLTCPU pq, String name) {
		super (pq, name);
		if (pq instanceof QuadCLT) {
			this.gpuQuad = ((QuadCLT) pq).gpuQuad; //  careful when switching - reset Geometry, vectors, bayer images. Kernels should be the same
		}
	}
	
	/**
	 * Remove weak non-LMA tiles if they do not have any LMA or strong neighbors and
	 * too few weak neighbors. Single strong neighbor within range is enough, strong/LMA
	 * that are not within range are still counted towards total "good" neighbors. For
	 * example, weak tiles of the sky just above strong sky-line will survive if the
	 * number of similar weak neighbor tiles representing clouds plus number of sky-line
	 * tiles is sufficient. Range may be somewhat wider than that for stronger tiles.
	 * Strong neibs with lower disparity count as weak ones (towards number of good neibs).
	 * There should be at least one "near".
	 *  
	 * @param dls {disparity, lma_disparity, strength} - not to be modified
	 * @param strong strength to be considered non weak 
	 * @param weak   strength to be considered weak, below - just delete (may be 0 if not used) 
	 * @param min_neibs minimal number of neighbors (of 8) to survive
	 * @param tolerance_absolute - absolute tolerance
	 * @param tolerance_relative
	 * @param width
	 * @param max_iter
	 * @param threadsMax
	 * @param debug_level
	 * @return updated disparity with some tiles replaced with Double.NaN
	 */
	
	public static double [] removeFewWeak(
			final double [][] dls,
			final double      strong,
			final double      weak,
			final int         min_neibs, 
			final double      tolerance_absolute,
			final double      tolerance_relative,
			final int         width,
			final int         max_iter,
			final int         threadsMax,
			final int         debug_level)
	{
		final int tiles = dls[0].length;
		final double [] disparity =     dls[0].clone();
		final double [] disparity_out = dls[0].clone();
		final double [] disparity_lma = dls[1].clone();
		final double [] strength =      dls[2]; // will not be updated
		final TileNeibs tn =  new TileNeibs(width, tiles/width);
		final Thread[] threads = ImageDtt.newThreadArray(threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		final AtomicInteger anum_updated = new AtomicInteger(0);
		final int dbg_tile = 2512;
		for (int iter = 0; iter < max_iter; iter++) {
			ai.set(0);
			anum_updated.set(0);
			for (int ithread = 0; ithread < threads.length; ithread++) {
				threads[ithread] = new Thread() {
					public void run() {
						for (int nTile = ai.getAndIncrement(); nTile < tiles; nTile = ai.getAndIncrement()) {
							if ((debug_level >0) && (nTile == dbg_tile)) {
								System.out.println("removeFewWeak():removeDisparityOutliers() nTile="+nTile);
							}
							if (    !Double.isNaN(disparity[nTile]) &&
									Double.isNaN(disparity_lma[nTile]) && // is not lma 
									(strength[nTile] < strong)) {         // weak
								if (strength[nTile] < weak) {
									disparity_out[nTile] = Double.NaN;
									anum_updated.getAndIncrement();
									continue;
								}
								double tolerance = tolerance_absolute +  ((disparity[nTile]>0)? disparity[nTile]*tolerance_relative:0);
								double lim_max = disparity[nTile] + tolerance; 
								double lim_min = disparity[nTile] - tolerance;
								int num_goog_neibs = 0; // strong or close
								int num_near = 0;
								boolean keep = false;
								for (int dir = 0; dir < 8; dir++) {
									int ineib = tn.getNeibIndex(nTile, dir);
									if ((ineib >= 0) && !Double.isNaN(disparity[ineib])) {
										boolean near = (disparity[ineib] >= lim_min) && (disparity[ineib] <= lim_max);
										boolean is_strong = !Double.isNaN(disparity_lma[ineib]) ||
												(strength[ineib] >= strong);
										// Strong and near only counts for strong in FG (nearer), far ones same as weak
										if (near && is_strong && (disparity[ineib] >= disparity[nTile])) {
											keep = true;
											break;
										}
										if (near || is_strong) {
											num_goog_neibs++;
											if (near) num_near++;
											if ((num_goog_neibs >= min_neibs) && (num_near > 0)) {
												keep = true;
												break;
											}
										}
									}
								}
								if (!keep) {
									disparity_out[nTile] = Double.NaN;
									anum_updated.getAndIncrement();
								}
							}
						}
					}
				};
			}		      
			ImageDtt.startAndJoin(threads);
			if (anum_updated.get() ==0) {
				break;
			}
			System.arraycopy(disparity_out,0,disparity,0,tiles);
		}
		return disparity_out;
	}
	
	
	
	
	
	
	

	public static double [] removeDisparityLMAOutliers( // just LMA FG
			final boolean     non_lma,
			final double [][] dls,
			final double      max_strength,  // do not touch stronger
			final int         nth_fromextrem, // 0 - compare to max/min. 1 - second max/min, ... 
			final double      tolerance_absolute,
			final double      tolerance_relative,
			final int         width,
			final int         max_iter,
			final int         threadsMax,
			final int         debug_level)
	{
		final int tiles = dls[0].length;
		final double [] disparity =     dls[0].clone();
		final double [] disparity_out = dls[0].clone();
		final double [] disparity_lma = dls[1].clone();
		final double [] strength =      dls[2]; // will not be updated
		final TileNeibs tn =  new TileNeibs(width, tiles/width);
		final Thread[] threads = ImageDtt.newThreadArray(threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		final AtomicInteger anum_updated = new AtomicInteger(0);
		final int dbg_tile = 2997;
		for (int iter = 0; iter < max_iter; iter++) {
			ai.set(0);
			anum_updated.set(0);
			for (int ithread = 0; ithread < threads.length; ithread++) {
				threads[ithread] = new Thread() {
					public void run() {
						double [] neibs = new double[8];
						for (int nTile = ai.getAndIncrement(); nTile < tiles; nTile = ai.getAndIncrement()) {
							if ((debug_level >0) && (nTile == dbg_tile)) {
								System.out.println("removeDisparityOutliers() nTile="+nTile);
							}
							if (    !Double.isNaN(disparity[nTile]) &&
									(!Double.isNaN(disparity_lma[nTile]) ^ non_lma) && // is_lma 
									(strength[nTile] < max_strength)) { // weak LMA
								Arrays.fill(neibs, Double.NaN);
								for (int dir = 0; dir < 8; dir++) {
									int ineib = tn.getNeibIndex(nTile, dir);
									if (ineib >= 0) {
										neibs[dir] = non_lma ? disparity[ineib] : disparity_lma[ineib];
										if (Double.isNaN(disparity[ineib])) {
											neibs[dir] =Double.NaN;
										}
									}
								}
								Arrays.sort(neibs); // increasing, NaNs - in the end
								int num_defined = neibs.length;
								for (int i = 0; i < neibs.length; i++) {
									if (Double.isNaN(neibs[i])) {
										num_defined = i;
										break;
									}
								}
								if (num_defined > 0) {
									int nth_min = nth_fromextrem;
									if (nth_min >= num_defined) {
										nth_min = num_defined - 1;
									}
									int nth_max = num_defined - 1 - nth_min;
									double d_max = neibs[nth_max] + tolerance_absolute + neibs[nth_max]*tolerance_relative;
									if ((non_lma ? disparity[nTile] : disparity_lma[nTile]) > d_max) {
										disparity_out[nTile] = Double.NaN;
									}
									if (disparity_out[nTile] != disparity[nTile]) {
										anum_updated.getAndIncrement();
									}
								}
							}
						}
					}
				};
			}		      
			ImageDtt.startAndJoin(threads);
			if (anum_updated.get() ==0) {
				break;
			}
			System.arraycopy(disparity_out,0,disparity,0,tiles);
		}
		return disparity_out;
	}
	
	/**
	 * Remove no-LMA tiles if they do not have a close neighbors that have LMA
	 * result. As an option - do not touch tiles that do not have LMA neighbors
	 * @param dls disparity, lma-disparity, strength (disparity and LMA disparity
	 *            have the same values, just some of LMA have NaN. 
	 * @param max_strength Do not touch tiles stronger than that
	 * @param diff_from_lma_pos Difference from farthest FG objects (OK to have large, e.g. 100)
	 * @param diff_from_lma_neg diff_from_lma_neg Difference from nearest BG objects (small, as FG are usually more visible)
	 * @param remove_no_lma_neib remove tiles that do not have LMA neighbors
	 * @param width       DSI width
	 * @param threadsMax
	 * @param debug_level
	 * @return new disparity array, removed tiles are Double.NaN
	 */
	public static double [] removeDisparityOutliersByLMA(
			final double [][] dls,
			final double      max_strength,        // do not touch stronger
			final double      diff_from_lma_pos,   // Difference from farthest FG objects (OK to have large, e.g. 100)
			final double      diff_from_lma_neg,   // Difference from nearest BG objects (small, as FG are usually more visible)
			final int         search_radius,       // Search farther if no LMA neighbor is found closer. Original value - 1 (8 neighbors)
			final boolean     remove_no_lma_neib,  // remove without LMA neighbors
			final int         width,               //tilesX
			final int         threadsMax,
			final int         debug_level)
	{
		final int tiles = dls[0].length;
		final double [] disparity =     dls[0].clone();
		final double [] disparity_out = dls[0].clone();
		final double [] disparity_lma = dls[1].clone();
		final double [] strength =      dls[2]; // will not be updated
		final TileNeibs tn =  new TileNeibs(width, tiles/width);
		final Thread[] threads = ImageDtt.newThreadArray(threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		final AtomicInteger anum_updated = new AtomicInteger(0);
		final int dbg_tile = 1944;
		anum_updated.set(0);
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				public void run() {
					for (int nTile = ai.getAndIncrement(); nTile < tiles; nTile = ai.getAndIncrement()) {
						if ((debug_level >0) && (nTile == dbg_tile)) {
							System.out.println("removeDisparityOutliers() nTile="+nTile);
						}
						if (Double.isNaN(disparity_lma[nTile]) && !Double.isNaN(disparity[nTile]) && (strength[nTile] < max_strength)) {
							double best_fit_pos = Double.NaN; // Closest higher disparity than this
							double best_fit_neg = Double.NaN; // Closest lower disparity than this
							for (int rad = 1; rad <= search_radius; rad++) {
								int numdir = TileNeibs.getNumDirs(rad);
								for (int dir = 0; dir < numdir; dir++) {
									int ineib = tn.getNeibIndexRadius(nTile, dir, rad);
									if (    (ineib >= 0) &&
											!Double.isNaN(disparity_lma[ineib]) &&
											!Double.isNaN(disparity[ineib])) {
										double d = disparity[nTile] - disparity_lma[ineib];
										if (d > 0) {
											if (!(d >= best_fit_neg)) {
												best_fit_neg = d;
											}
										} else {
											if (!(-d >= best_fit_pos)) {
												best_fit_neg = -d;
											}
										}
									}
								}
								if (!(Double.isNaN(best_fit_pos) && Double.isNaN(best_fit_neg))) {
									break;
								}
							}
							if (    (best_fit_neg > diff_from_lma_neg) ||
									(best_fit_pos > diff_from_lma_pos) ||
									(Double.isNaN(best_fit_pos) && Double.isNaN(best_fit_neg) && remove_no_lma_neib)) {
								disparity_out[nTile] = Double.NaN;
							}
						}
					}
				}
			};
		}		      
		ImageDtt.startAndJoin(threads);
		return disparity_out;
	}
	
	
	
	public static double [] removeDisparityOutliers(
			final double [][] ds0,
			final double      max_strength,  // do not touch stronger
			final int         nth_fromextrem, // 0 - compare to max/min. 1 - second max/min, ... 
			final double      tolerance_absolute,
			final double      tolerance_relative,
			final int         width,
			final int         max_iter,
			final boolean     fit_completely, // do not add tolerance when replacing
			final int         threadsMax,
			final int         debug_level)
	{
		final int tiles = ds0[0].length;
		final double [] disparity =     ds0[0].clone();
		final double [] disparity_out = ds0[0].clone();
		final double [] strength =      ds0[1]; // will not be updated
		final TileNeibs tn =  new TileNeibs(width, tiles/width);
		final Thread[] threads = ImageDtt.newThreadArray(threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		final AtomicInteger anum_updated = new AtomicInteger(0);
		final int dbg_tile = 8309;
		for (int iter = 0; iter < max_iter; iter++) {
			ai.set(0);
			anum_updated.set(0);
			for (int ithread = 0; ithread < threads.length; ithread++) {
				threads[ithread] = new Thread() {
					public void run() {
						double [] neibs = new double[8];
						for (int nTile = ai.getAndIncrement(); nTile < tiles; nTile = ai.getAndIncrement()) {
							if ((debug_level >0) && (nTile == dbg_tile)) {
								System.out.println("removeDisparityOutliers() nTile="+nTile);
							}
							if ((strength[nTile] < max_strength) && !Double.isNaN(disparity[nTile])) {
								Arrays.fill(neibs, Double.NaN);
								for (int dir = 0; dir < 8; dir++) {
									int ineib = tn.getNeibIndex(nTile, dir);
									if (ineib >= 0) {
										neibs[dir] = disparity[ineib];
									}
								}
								Arrays.sort(neibs);
								int num_defined = neibs.length;
								for (int i = 0; i < neibs.length; i++) {
									if (Double.isNaN(neibs[i])) {
										num_defined = i;
										break;
									}
								}
								if (num_defined > 0) {
									int nth_min = nth_fromextrem;
									if (nth_min >= num_defined) {
										nth_min = num_defined - 1;
									}
									int nth_max = num_defined - 1 - nth_min;
									double d_min = neibs[nth_min] - tolerance_absolute - neibs[nth_min]*tolerance_relative;
									double d_max = neibs[nth_max] + tolerance_absolute + neibs[nth_max]*tolerance_relative;
									if  (disparity[nTile] < d_min) {
										disparity_out[nTile] = fit_completely? neibs[nth_min]: d_min;
									} else if (disparity[nTile] > d_max) {
										disparity_out[nTile] = fit_completely? neibs[nth_max]: d_max;
									}
									if (disparity_out[nTile] != disparity[nTile]) {
										anum_updated.getAndIncrement();
									}
								}
							}
						}
					}
				};
			}		      
			ImageDtt.startAndJoin(threads);
			if (anum_updated.get() ==0) {
				break;
			}
			System.arraycopy(disparity_out,0,disparity,0,tiles);
		}
		return disparity_out;
	}

	// trying new version
	public static double [][] fillDisparityStrength(
			final double [][] ds0,
			final double      min_disparity,
			final double      max_sym_disparity, // lower disparity - num_bottom = 8; 
			final int         num_bottom, // average this number of lowest disparity neighbors (of 8)
			final int         num_passes,
			final double      max_change,
			final int         width,
			final int         threadsMax,
			final int         debug_level)
	{
		final double max_disp_diff =     1.00; // no pull from pixels that are this high than average of lowest
		final double max_disp_diff_rel = 0.2; // no pull from pixels that are this high than average of lowest
		final double anchor_up_abs =     0.2; // if anchor is more than that larger, pull as if it is that latger
		final double anchor_up_rel =     0.2; // add to anchor_up_abs multiplied by anchor disparity
		final int    min_float_neibs =   2; // minimal number of float neibs to override fixed;
		
        final double   diagonal_weight = 0.5 * Math.sqrt(2.0); // relative to ortho
		double wdiag = 0.25 *diagonal_weight / (diagonal_weight + 1.0);
		double wortho = 0.25 / (diagonal_weight + 1.0);
		final double [] neibw = {wortho, wdiag, wortho, wdiag, wortho, wdiag, wortho, wdiag}; 

		final double max_change2 = max_change * max_change;
		final double [][] ds = new double [][] {ds0[0].clone(),ds0[1].clone()};
		final int tiles = ds[0].length;
		final TileNeibs tn =  new TileNeibs(width, tiles/width);
		
		final boolean [] floating =      new boolean[tiles]; // which tiles will change
		final double [] anchors =        new double [tiles]; // average of the known fixed neighbors
		final double [] anchor_weights = new double [tiles]; // weights of anchors
		Arrays.fill(anchors,Double.NaN);
		
		final int [] tile_indices = new int [tiles];
		final Thread[] threads = ImageDtt.newThreadArray(threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		final AtomicInteger anum_gaps = new AtomicInteger(0);
		final int dbg_tile = -3379;
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				public void run() {
					for (int nTile = ai.getAndIncrement(); nTile < tiles; nTile = ai.getAndIncrement()) {
						if (ds[0][nTile] < min_disparity) {
							ds[0][nTile] = Double.NaN; // min_disparity;
						}
						if (Double.isNaN(ds[0][nTile]) || (ds[1][nTile] <= 0)) { // blue_sky is made small >0 weight
							ds[0][nTile] = Double.NaN;
							ds[1][nTile] = 0.0;
							int indx = anum_gaps.getAndIncrement();
							tile_indices[indx] = nTile;
							floating[nTile] = true;
						}
					}
				}
			};
		}		      
		ImageDtt.startAndJoin(threads);
		ai.set(0);
		final int num_gaps = anum_gaps.get(); 
		if (num_gaps == 0) {
			return ds; // no gaps already
		}
		
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				public void run() {
					double [] neibs = new double[8];
					for (int indx = ai.getAndIncrement(); indx < num_gaps; indx = ai.getAndIncrement()) {
						int nTile = tile_indices[indx];
						if ((debug_level >0) && (nTile == dbg_tile)) {
							System.out.println("fillDisparityStrength() nTile="+nTile);
						}
						Arrays.fill(neibs, Double.NaN);
						double min_def_disp = Double.NaN;
						for (int dir = 0; dir < 8; dir++) {
							int nt_neib = tn.getNeibIndex(nTile, dir);
							if ((nt_neib >= 0) && !floating[nt_neib]) {
								neibs[dir] = ds[0][nt_neib];
								if (!(neibs[dir] >= min_def_disp)) { // true if was NaN or new is smaller
									min_def_disp = neibs[dir];
								}
							}
						}
						if (!Double.isNaN(min_def_disp)) { // no fixed neighbors
							double max_to_vag = min_def_disp * (1.0 + max_disp_diff_rel) + max_disp_diff;
							double swd = 0.0, sw = 0.0;
							for (int dir = 0; dir < 8; dir++) {
								if (neibs[dir] < max_to_vag) {// NaN OK, will be false
									sw +=  neibw[dir];
									swd += neibw[dir] * neibs[dir]; 
								}
							}
							anchors[nTile] = swd / sw; // here = not zero as there will be at least one tile for min_def_disp
							anchor_weights[nTile] = sw;
						}
					}
				}
			};
		}		      
		ImageDtt.startAndJoin(threads);
		
		ai.set(0);
		final boolean [] fill_all = {false};
		final double [] disp_new = ds[0].clone();
		DoubleAccumulator amax_diff =  new DoubleAccumulator (Double::max, Double.NEGATIVE_INFINITY);
		for (int npass = 0; npass < num_passes; npass+= fill_all[0]? 1:0 ) { // do not limit initial passes
			anum_gaps.set(0);
			amax_diff.reset();
			for (int ithread = 0; ithread < threads.length; ithread++) {
				threads[ithread] = new Thread() {
					public void run() {
						double [] neibs = new double[8];
						for (int indx = ai.getAndIncrement(); indx < num_gaps; indx = ai.getAndIncrement()) {
							int nTile = tile_indices[indx];
							if ((debug_level >0) && (nTile == dbg_tile)) {
								System.out.println("fillDisparityStrength() nTile="+nTile);
							}
							if (!fill_all[0] && !Double.isNaN(ds[0][nTile])) {
								continue; // fill only new
							}
							Arrays.fill(neibs, Double.NaN);
							// average only floating, use anchors for fixed
							double swd = 0.0, sw = 0.0;
							double new_val; //  = Double.NaN;
							int num_floats = 0;
							for (int dir = 0; dir < 8; dir++) {
								int nt_neib = tn.getNeibIndex(nTile, dir);
								if ((nt_neib >= 0) && floating[nt_neib]) {
									neibs[dir] = ds[0][nt_neib];
									if (! Double.isNaN(neibs[dir])) {
										sw +=  neibw[dir];
										swd += neibw[dir] * neibs[dir];
										num_floats++;
									}
								}
							}
							if ((sw > 0) || (anchor_weights[nTile] > 0)) {
								// (num_floats >= min_float_neibs) to disallow floats BG pull down through just diagonal 
								if ((sw > 0) && (anchor_weights[nTile] > 0) && (num_floats >= min_float_neibs)) {
									double avg_flt = swd/sw;
									double max_fix = avg_flt + anchor_up_abs+ avg_flt*anchor_up_rel;
									double avg_fix = anchors[nTile];
									if (avg_fix > max_fix) {
										avg_fix = max_fix;
									}
									new_val = (avg_fix * anchor_weights[nTile] + swd)/(sw + anchor_weights[nTile]);
								} else if (anchor_weights[nTile] > 0) { // (sw == 0) || (num_floats < min_float_neibs)
									new_val = anchors[nTile];
								} else { // OK even if number of float neighbors is <  min_float_neibs - there are no fixed competitors 
									new_val = swd/sw;
								}
								if (fill_all[0]) {
									double d2 = max_change2 * 2;
									if (!Double.isNaN(disp_new[nTile])){
										double d = new_val -  disp_new[nTile];
										d2 = d * d;
										if ((debug_level > 0) &&  (d2 > 15)) {
											System.out.println(
													"fillDisparityStrength(): nTile="+nTile+
													" tileX="+(nTile%width)+ " tileY="+(nTile/width)+
													", d="+d);
										}
									}
									amax_diff.accumulate(d2);
								} else {
									anum_gaps.getAndIncrement();
								}
								disp_new[nTile] = new_val;
							}
						}
					}
				};
			}		      
			ImageDtt.startAndJoin(threads);
			ai.set(0);
			System.arraycopy(disp_new, 0, ds[0], 0, tiles);
			if ((debug_level > 0) && fill_all[0]) {
				System.out.println("fillDisparityStrength() num_gaps="+num_gaps+", npass="+npass+", change="+Math.sqrt(amax_diff.get())+" ("+max_change+")");
			}
			if (fill_all[0] && (amax_diff.get() < max_change2)) {
				break;
			}
			if (anum_gaps.get() == 0) { // no new tiles filled
				fill_all[0] = true; 
			}
			if ((debug_level>0) && (npass == (num_passes-1))){
				System.out.println("fillDisparityStrength() LAST PASS ! npass="+npass+", change="+Math.sqrt(amax_diff.get())+" ("+max_change+")");
				System.out.println("fillDisparityStrength() LAST PASS ! npass="+npass+", change="+Math.sqrt(amax_diff.get())+" ("+max_change+")");
				System.out.println("fillDisparityStrength() LAST PASS ! npass="+npass+", change="+Math.sqrt(amax_diff.get())+" ("+max_change+")");
			}
		} // for (int npass = 0; npass < num_passes; npass+= fill_all[0]? 1:0 )
		return ds;
	}
	
	public static double[] blurWeak(
			double [][] ds,
			double      min_strength_blur,
			double      min_strength_replace,
			int         n,
			int         width,
			double      sigma) {
		double [] disparity = new double[ds[0].length];
		Arrays.fill(disparity, Double.NaN);
		for (int i = 0; i < disparity.length; i++) {
			if (Double.isNaN(ds[0][i]) || (ds[1][i] < min_strength_blur)) {
				disparity[i] = ds[0][i];
			}
		}
		for (int i = 0; i < n; i++) { 
			disparity =   (new DoubleGaussianBlur()).blurWithNaN(
					disparity,              // double[] pixels,
					null,                   // double [] in_weight, // or null
					width,                  // int width,
					disparity.length/width, // int height,
					sigma,                  // double sigmaX,
					sigma,                  // double sigmaY,
					0.00001);               // double accuracy);
		}
		for (int i = 0; i < disparity.length; i++) {
			if (ds[1][i] > min_strength_replace) {
				disparity[i] = ds[0][i];
			}
		}
		return disparity;
	}
	
	public void setGPU(GpuQuad gpuQuad) {
		this.gpuQuad = gpuQuad;
	}
	public GpuQuad getGPU() {
		return this.gpuQuad;
	}
	
	
	public void genSave4sliceImage(
			CLTParameters                                   clt_parameters,
			String                                          suffix,
			EyesisCorrectionParameters.DebayerParameters    debayerParameters,
			ColorProcParameters                             colorProcParameters,
			CorrectionColorProc.ColorGainsParameters        channelGainParameters,
			EyesisCorrectionParameters.RGBParameters        rgbParameters,
			final int                                       threadsMax,  // maximal number of threads to launch
			final int                                       debugLevel){
///		String x3d_path = getX3dDirectory();
///		String file_name = image_name + suffix;
///		String file_path = x3d_path + Prefs.getFileSeparator() + file_name + ".tiff";
		if ((getGPU() != null) && (getGPU().getQuadCLT() != this)) {
			getGPU().updateQuadCLT(this); // to re-load new set of Bayer images to the GPU
		}

///		ImagePlus img_noise = 
				processCLTQuadCorrGPU(
				null,                  // ImagePlus []                        imp_quad, //null will be OK
				null,                  // boolean [][]                                    saturation_imp, // (near) saturated pixels or null // Not needed use this.saturation_imp
				clt_parameters,        // CLTParameters                                   clt_parameters,
				debayerParameters,     // EyesisCorrectionParameters.DebayerParameters    debayerParameters,
				colorProcParameters,   // ColorProcParameters                             colorProcParameters,
				channelGainParameters, // CorrectionColorProc.ColorGainsParameters        channelGainParameters,
				rgbParameters,         // EyesisCorrectionParameters.RGBParameters        rgbParameters,
				null,                  // double []	                                      scaleExposures, // probably not needed here - restores brightness of the final image
				true,                  // boolean                                         only4slice,
				threadsMax,            // final int                                       threadsMax,  // maximal number of threads to launch
				false,                 // final boolean                                   updateStatus,
				debugLevel);           // final int                                       debugLevel);
//		FileSaver fs=new FileSaver(img_noise); // is null, will be saved inside to /home/elphel/lwir16-proc/proc1/results_cuda/1626032208_613623-AUX-SHIFTED-D0.0
//		fs.saveAsTiff(file_path);
	}	
	
	
	// added 12/21
	public void processCLTQuadCorrs( // not used in lwir
			CLTParameters           clt_parameters,
			EyesisCorrectionParameters.DebayerParameters     debayerParameters,
			ColorProcParameters colorProcParameters,
			CorrectionColorProc.ColorGainsParameters     channelGainParameters,
			EyesisCorrectionParameters.RGBParameters             rgbParameters,
			//			  int          convolveFFTSize, // 128 - fft size, kernel size should be size/2
			final boolean    apply_corr,    // calculate and apply additional fine geometry correction
			final boolean    infinity_corr, // calculate and apply geometry correction at infinity
			final int        threadsMax,    // maximal number of threads to launch
			final boolean    updateStatus,
			final int        debugLevel)
	{
		if (infinity_corr && (clt_parameters.z_correction != 0.0)){
			System.out.println(
					"****************************************\n"+
							"* Resetting manual infinity correction *\n"+
					"****************************************\n");
			clt_parameters.z_correction = 0.0;
		}

		this.startTime=System.nanoTime();
		String [] sourceFiles=correctionsParameters.getSourcePaths();
		SetChannels [] set_channels=setChannels(debugLevel);
		if ((set_channels == null) || (set_channels.length==0)) {
			System.out.println("No files to process (of "+sourceFiles.length+")");
			return;
		}
		// multiply each image by this and divide by individual (if not NaN)
		double [] referenceExposures = null;
		//		  if (!colorProcParameters.lwir_islwir) {
		if (!isLwir()) {
			referenceExposures=eyesisCorrections.calcReferenceExposures(debugLevel);
		}
		for (int nSet = 0; nSet < set_channels.length; nSet++){
			int [] channelFiles = set_channels[nSet].fileNumber();
			boolean [][] saturation_imp = (clt_parameters.sat_level > 0.0)? new boolean[channelFiles.length][] : null;
			double [] scaleExposures = new double[channelFiles.length];

			//			  ImagePlus [] imp_srcs = 
			conditionImageSet(
					clt_parameters,             // EyesisCorrectionParameters.CLTParameters  clt_parameters,
					colorProcParameters,
					sourceFiles,                // String []                                 sourceFiles,
					set_channels[nSet].name(),  // String                                    set_name,
					referenceExposures,         // double []                                 referenceExposures,
					channelFiles,               // int []                                    channelFiles,
					scaleExposures,   //output  // double [] scaleExposures
					saturation_imp,   //output  // boolean [][]                              saturation_imp,
					threadsMax,                 // int                                       threadsMax,
					debugLevel); // int                                       debugLevel);


			// once per quad here
			processCLTQuadCorrGPU(  // returns ImagePlus, but it already should be saved/shown
					null,           // imp_srcs, // [srcChannel], // should have properties "name"(base for saving results), "channel","path"
					saturation_imp, // boolean [][] saturation_imp, // (near) saturated pixels or null
					clt_parameters,
					debayerParameters,
					colorProcParameters,
					channelGainParameters,
					rgbParameters,
					scaleExposures,
					false, // true,            // boolean                                         only4slice,
					threadsMax,                // final int                                       threadsMax,  // maximal number of threads to launch
					false,                     // final boolean                                   updateStatus,
					debugLevel);               // final int                                       debugLevel);
			if (debugLevel >-1) System.out.println("Processing set "+(nSet+1)+" (of "+set_channels.length+") finished at "+
					IJ.d2s(0.000000001*(System.nanoTime()-this.startTime),3)+" sec, --- Free memory1="+Runtime.getRuntime().freeMemory()+" (of "+Runtime.getRuntime().totalMemory()+")");
			if (eyesisCorrections.stopRequested.get()>0) {
				System.out.println("User requested stop");
				System.out.println("Processing "+(nSet + 1)+" file sets (of "+set_channels.length+") finished at "+
						IJ.d2s(0.000000001*(System.nanoTime()-this.startTime),3)+" sec, --- Free memory2="+Runtime.getRuntime().freeMemory()+" (of "+Runtime.getRuntime().totalMemory()+")");
				return;
			}
		}
		System.out.println("processCLTQuadCorrs(): processing "+getTotalFiles(set_channels)+" files ("+set_channels.length+" file sets) finished at "+
				IJ.d2s(0.000000001*(System.nanoTime()-this.startTime),3)+" sec, --- Free memory3="+Runtime.getRuntime().freeMemory()+" (of "+Runtime.getRuntime().totalMemory()+")");
	}

	public static boolean calibratePhotometric2( // quadratic
			CLTParameters     clt_parameters,
			final QuadCLT     ref_scene, // now - may be null - for testing if scene is rotated ref
			final double      min_strength,
			final double      max_diff,     //  30.0
			int               photo_order,  //   0 - offset only, 1 - linear, 2 - quadratic
			final double      photo_std_1,  //  50.0;   // Minimal standard deviation of the filtered values for poly order 1
			final double      photo_std_2,	// 200.0;   // Minimal standard deviation of the filtered values for poly order 2
			final int         photo_offs_set, // 0;     // 0 - keep weighted offset average, 1 - balance result image, 2 - set weighted average to specific value
			final double      photo_offs,     // 21946; // weighted average offset target value, if photo_offs_set (and not photo_offs_balance)
			final int         num_refines, // 2
			final int         min_good,     // minimal number of "good" pixels
			final double [][] combo_dsn_final,     // double [][]    combo_dsn_final, // dls,
			final boolean[]   blue_sky,
			int               threadsMax,
			final boolean     debug)
	{
		
		// filter disparity by LMA only, same bg and fg
			double []         disparity_ref= combo_dsn_final[OpticalFlow.COMBO_DSN_INDX_DISP].clone();
			for (int i = 0; i < disparity_ref.length; i++) {
				if (combo_dsn_final[OpticalFlow.COMBO_DSN_INDX_STRENGTH][i] < min_strength) {
					disparity_ref[i] = Double.NaN;
				} else if (combo_dsn_final[OpticalFlow.COMBO_DSN_INDX_DISP_BG_ALL][i] !=
						combo_dsn_final[OpticalFlow.COMBO_DSN_INDX_DISP][i]) {
					disparity_ref[i] = Double.NaN;
				}
			}
			if (blue_sky != null) {
				for (int i = 0; i < blue_sky.length; i++) if (blue_sky[i]) {
						disparity_ref[i] = 0.0;
				}
			}
			
			ImagePlus img_ref = renderGPUFromDSI(
					-1,                 // final int         sensor_mask,
					false,              // final boolean     merge_channels,
					null,               // final Rectangle   full_woi_in,      // show larger than sensor WOI in tiles (or null)
					 clt_parameters,    // CLTParameters     clt_parameters,
					 disparity_ref,     // double []         disparity_ref,
					 OpticalFlow.ZERO3, // final double []   scene_xyz, // camera center in world coordinates
					 OpticalFlow.ZERO3, // final double []   scene_atr, // camera orientation relative to world frame
					 ref_scene,         // final QuadCLT     scene,
					 ref_scene,         // final QuadCLT     ref_scene, // now - may be null - for testing if scene is rotated ref
					false,              // final boolean     toRGB,
					true,               // final boolean     show_nan,
					"PHOTOMETRIC",      // String            suffix,
					threadsMax,         // int               threadsMax,
					-2);        // final int         debugLevel);
			if (debug) {
				img_ref.show();
			}
			ImageStack imageStack = img_ref.getStack();
			int num_sens=imageStack.getSize();
			float [] fpixels;
			double [][] dpixels = new double [num_sens][];
			for (int n = 0; n < num_sens; n++) {
				fpixels = (float[]) imageStack.getPixels(n + 1);
				dpixels[n] = new double [fpixels.length];
				for (int i = 0; i < fpixels.length; i++) {
					dpixels[n][i] = fpixels[i];
				}
			}
//			boolean quadratic = true;
			int width =  img_ref.getWidth();
			int height = img_ref.getHeight();
			int len = width* height;
			double [] avg_pix = new double [len];
			boolean [] good_pix = new boolean[len];
			double [][] pa_coeff = new double[num_sens][];
			double min_abs_a = 1E-9;
			Arrays.fill(good_pix, true);
			
			double [] offs_old =    ref_scene.getLwirOffsets();
			double [] scales_old =  ref_scene.getLwirScales();
			double [] scales2_old = ref_scene.getLwirScales2();
			// Calculate avg_pix - average (for all 16 channels)  pixel value and initial good_pix
			// next passes avg_pix will be modified and more bad pixels added
			
			double avg_img = 0.0;
			int num_good = 0;
			for (int i = 0; i < len; i++) {
				for (int n = 0; n < num_sens; n++) {
					avg_pix[i]+=dpixels[n][i];
					good_pix[i] &= !Double.isNaN(dpixels[n][i]);
				}
				avg_pix[i] /= dpixels.length;
				if (good_pix[i]) {
					avg_img += avg_pix[i];
					num_good++;
				}
			}
			if (num_good < min_good) {
				System.out.println("calibratePhotometric2(): Too few good pixels for calibration: "+
			num_good+" < "+min_good+", abandoning photometric calibration.");
				return false;
			}
			avg_img /= num_good;
			double wavg_offs = photo_offs;
			boolean set_offs = true;
			if (photo_offs_set == 1) { // balance
				for (int i = 0; i < len; i++) {
					avg_pix[i] -= avg_img;
				}
				set_offs = false;
			} else if (photo_offs_set == 0) {
				double s = 1.0;
				wavg_offs = 0.0;
				for (int n = 0; n < num_sens; n++) {
					s *= scales_old[n];
					wavg_offs += offs_old[n]*scales_old[n];
				}
				s = Math.pow(s, 1.0/num_sens);
				wavg_offs /= num_sens * s;
			}
			
			// recover original "raw" pixel values undoing old corrections
			double [][] raw = new double [num_sens][len];
			for (int nsens = 0; nsens < num_sens; nsens++) {
				Arrays.fill(raw[nsens],Double.NaN); // debug only
				for (int i = 0; i < len; i++) if (good_pix[i]){
					double A0 = scales2_old[nsens];
					double B0 = scales_old[nsens];
					double C0 = offs_old[nsens];
					double a = A0;
					double b = (B0 - 2 * C0 * A0);
					double c = (C0*C0*A0 -C0 * B0 - dpixels[nsens][i]);
					double p = -c/b;
					if (Math.abs(a) >= min_abs_a) {
						double d2 = b*b - 4 * a * c;
						if (d2 < 0) {
							System.out.println("calibratePhotometric2() 0: Failed quadratic: a = "+a+", b = "+b+", c="+c);
							return false;
						}
						p= (-b + Math.sqrt(d2))/(2 * a);
					}
					raw[nsens][i] = p;
				}
            }
			
			double [] offs_new =    new double [num_sens];
			double [] scales_new =  new double [num_sens];
			double [] scales2_new = new double [num_sens];
			double [] avg_chn =     new double [num_sens];
			double [] std_chn =     new double [num_sens];
			double std = 0.0;
			for (int nref = 0; nref < num_refines; nref++) {
				Arrays.fill(avg_chn, 0.0);
				Arrays.fill(std_chn, 0.0);
				num_good = 0;
				for (int i = 0; i < len; i++) {
					if (good_pix[i]) {
						num_good++;
						for (int n = 0; n < num_sens; n++) {
							avg_chn[n] += dpixels[n][i];
							std_chn[n] += dpixels[n][i] * dpixels[n][i];
						}
					}
				}
				std = 0.0;
				for (int n = 0; n < num_sens; n++) {
					avg_chn[n] /= num_good;
					std_chn[n] /= num_good;
					std_chn[n] = Math.sqrt(std_chn[n] - avg_chn[n] * avg_chn[n]);
					std += std_chn[n] * std_chn[n];
				}
				std = Math.sqrt(std/num_sens);
				if ((photo_order > 0) && (std < photo_std_1)) {
					photo_order = 0;
					System.out.println ("Standard deviation = "+std+" < "+photo_std_1+", adjusting only offsets. nref = "+nref);
					nref = -1;
					continue; // restart
				} else if ((photo_order > 1) && (std < photo_std_2)) {
					photo_order = 1;
					System.out.println ("Standard deviation = "+std+" < "+photo_std_2+", adjusting only offsets and scales. nref = "+nref);
					nref = -1;
					continue; // restart
				} else if (debug) {
					System.out.println ("Standard deviation = "+std+", nref = "+nref);
				}
				
				double [][][] pa_data = new double [num_sens][num_good][2];
				for (int nsens = 0; nsens < num_sens; nsens++) {
					int indx = 0;
					for (int i = 0; i < len; i++) if (good_pix[i]){
						pa_data[nsens][indx][0] = raw[nsens][i]; // dpixels[nsens][i]; 
						pa_data[nsens][indx][1] = avg_pix[i];
						indx++;
					}
					double c, b, a;
					int poly_debug = 0;
					if (photo_order < 1) { // zero out scales2, keep old scales (even if not normalized)
						a = 0;
						b = scales_old[nsens];
						double s = 0.0;
						for (int i = 0; i < num_good; i++) {
							s += pa_data[nsens][i][0] - pa_data[nsens][i][1]/b; 
						}
						c = -s/num_good*b;
						if (Double.isNaN(c)) {
							System.out.println("calibratePhotometric2().1 c=NaN");
							return false;
						}
					} else {
					// need to balance quadratic so their average is 0	
					// linear or quadratic
						pa_coeff[nsens] =(new PolynomialApproximation(poly_debug)).polynomialApproximation1d(pa_data[nsens], photo_order);
						c = pa_coeff[nsens][0];
						b = pa_coeff[nsens][1];
						a = (pa_coeff[nsens].length > 2) ? pa_coeff[nsens][2] : 0.0;
						if (Double.isNaN(c)) {
							System.out.println("calibratePhotometric2().2 c=NaN");
							return false;
						}
					}
					double A = a;
					double C = -c/b;
					if (Double.isNaN(C)) {
						System.out.println("calibratePhotometric2().3 C=NaN");
						return false;
					}

					if (Math.abs(a) >= min_abs_a) {
						double d2 = b*b - 4*a*c;
						if (d2 < 0) {
							System.out.println("calibratePhotometric2() 1: Failed quadratic: a = "+a+", b = "+b+", c="+c);
							return false;
						}
						C = (-b + Math.sqrt(d2))/(2 * a);
						if (Double.isNaN(C)) {
							System.out.println("calibratePhotometric2().4 C=NaN");
							return false;
						}
					}
					double B = 2 * C * a + b;
					scales2_new[nsens] = A;
					scales_new [nsens] = B;
					offs_new[nsens] =    C;
				}
				// Make scales_new to be (geometric) average 1.0 (except photo_order == 0), update pa_data to match
				if (photo_order > 0) {
					double scales_avg = 1.0;
					for (int nsens = 0; nsens < num_sens; nsens++) {
						scales_avg *= scales_new [nsens];
					}
					scales_avg = Math.pow(scales_avg, 1.0/num_sens);
					for (int nsens = 0; nsens < num_sens; nsens++) {
						scales_new [nsens] /= scales_avg;
						scales2_new [nsens] /= scales_avg;
						for (int i = 0; i < num_good; i++) {
							pa_data[nsens][i][1] /= scales_avg;
						}
					}
					for (int i = 0; i < avg_pix.length; i++) { // for the next iteration
						avg_pix[i] /= scales_avg;
					}
					//avg_pix[i]
					if (debug) {
						System.out.println ("scales_avg = "+scales_avg+", nref = "+nref);
					}

				}				
				// make scales2 be average zero, re-run scales
				if (photo_order > 1) {
					double scales2_offset = 0;
					for (int nsens = 0; nsens < num_sens; nsens++) {
						scales2_offset += scales2_new[nsens];
					}
					scales2_offset /= num_sens;
					for (int nsens = 0; nsens < num_sens; nsens++) {
						scales2_new[nsens] -= scales2_offset;
					}
					// modify pa_data, re-run linear
					
					for (int nsens = 0; nsens < num_sens; nsens++) {
						int indx = 0;
//						for (int i = 0; i < num_good; i++) {
						for (int i = 0; i < avg_pix.length; i++) if (good_pix[i]){
							double poffs = pa_data[nsens][indx][0] - offs_new[nsens];
							pa_data[nsens][indx][1] -= poffs*poffs*scales2_new[nsens];
//							avg_pix[i] -= poffs*poffs*scales2_offset/num_sens; // for the next iteration
							indx++;
						}
						pa_coeff[nsens] =(new PolynomialApproximation(0)).polynomialApproximation1d(
								pa_data[nsens], 1);
						double c = pa_coeff[nsens][0];
						double b = pa_coeff[nsens][1];
						scales_new [nsens] = b;
						offs_new[nsens] =   -c/b;
						if (Double.isNaN(offs_new[nsens])) {
							System.out.println("calibratePhotometric2().5 offs_new[nsens]=NaN");
							return false;
						}

					}
					// re-normalize scales_new
					double scales_avg = 1.0;
					for (int nsens = 0; nsens < num_sens; nsens++) {
						scales_avg *= scales_new [nsens];
					}
					scales_avg = Math.pow(scales_avg, 1.0/num_sens);
					for (int nsens = 0; nsens < num_sens; nsens++) {
						scales_new [nsens] /= scales_avg;
						scales2_new [nsens] /= scales_avg;
						for (int i = 0; i < num_good; i++) {
							pa_data[nsens][i][1] /= scales_avg;
						}
					}
					for (int i = 0; i < avg_pix.length; i++) { // for the next iteration
						avg_pix[i] /= scales_avg;
					}
					//avg_pix[i]
					if (debug) {
						System.out.println ("normalization after quadratic: scales_avg = "+scales_avg+", nref = "+nref);
					}
				}
				
				// update offsets if set_offs
				if (set_offs) {
					double wa = 0.0;
					for (int nsens = 0; nsens < num_sens; nsens++) {
						wa += offs_new[nsens] * scales_new [nsens];
					}
					wa /= num_sens;
					double a = wavg_offs - wa;
					for (int nsens = 0; nsens < num_sens; nsens++) {
						offs_new[nsens] += a/scales_new [nsens];
					}
					for (int i = 0; i < avg_pix.length; i++) if (good_pix[i]){
						avg_pix[i] -= a;
					}
					if (debug) {
						System.out.println ("calibratePhotometric2(): wa = "+wa+", a="+a);
					}
				}
				
				
				
				if (true) { // debug) { during debug,
					System.out.println("DEBUG: calibratePhotometric() nref="+nref);
					System.out.println(String.format("%3s %10s %8s %8s %10s %8s %8s",
						//  "chn","   c   ","   b   ", " 1e6*a "));          
							"chn","  offs0  ","scale0 ", "scale20","  offs "," scale ", " scale2"));
					for (int n = 0; n < num_sens; n++) {
						System.out.println(String.format("%2d: %10.4f %8.5f %8.5f %10.4f %8.5f %8.5f",
								n,offs_old[n],scales_old[n], 1e6*scales2_old[n],
								offs_new[n],scales_new[n], 1e6*scales2_new[n]));
					}
					System.out.println();
				}
				if (debug) {
					double [][] diffs = new double [num_sens][len];
					for (int n = 0; n < num_sens; n++) {
						Arrays.fill(diffs[n], Double.NaN);
					}
					for (int i = 0; i < len; i++)  if (good_pix[i]) { // if (!Double.isNaN(avg_pix[i])){
						for (int n = 0; n < num_sens; n++) {
							double poffs = raw[n][i] - offs_new[n];
							diffs[n][i] = -(avg_pix[i] -
									(poffs * scales_new[n] + poffs*poffs*scales2_new[n]));  
						}
					}
					ShowDoubleFloatArrays.showArrays( // out of boundary 15
							diffs,
							width,
							height,
							true,
							"photometric-quad-err"+nref);
					double [][] corrected = new double [num_sens][len];
					for (int n = 0; n < num_sens; n++) {
						Arrays.fill(corrected[n], Double.NaN);
						for (int i = 0; i < len; i++)  if (!Double.isNaN(dpixels[n][i])){
							double poffs = raw[n][i] - offs_new[n];
							corrected[n][i] =	(poffs * scales_new[n] + poffs*poffs*scales2_new[n]);  
						}
					}

					ShowDoubleFloatArrays.showArrays( // out of boundary 15
							corrected,
							width,
							height,
							true,
							"photometric-quad-corr"+nref);
				}				
				//max_diff
				if (nref < (num_refines-1)) {
					for (int i = 0; i < len; i++)  if (good_pix[i]) { // !Double.isNaN(avg_pix[i])){
						for (int n = 0; n < num_sens; n++) {
							double poffs = raw[n][i] - offs_new[n];
							double diff = Math.abs(avg_pix[i] -
									(poffs * scales_new[n] + poffs*poffs*scales2_new[n]));
							if (diff > max_diff) {
								good_pix[i] = false;
								break;
							}
						}
					}
				}
			}
			System.out.println("calibratePhotometric() Updated calibration:");
			for (int n = 0; n < num_sens; n++) {
				System.out.println(String.format("%2d: %10.4f %8.6f %8.6f", n,offs_new[n],scales_new[n], 1E6*scales2_new[n]));
			}
			System.out.println();
			ref_scene.setLwirOffsets (offs_new);
			ref_scene.setLwirScales  (scales_new);
			ref_scene.setLwirScales2 (scales2_new);
			ref_scene.setPhotometricScene(); // use ref_scene name
			return true;
	}
	
	
	public static ImagePlus renderGPUFromDSI(
			final int         sensor_mask,
			final boolean     merge_channels,
			final Rectangle   full_woi_in,      // show larger than sensor WOI in tiles (or null)
			CLTParameters     clt_parameters,
			double []         disparity_ref,
			final double []   scene_xyz, // camera center in world coordinates
			final double []   scene_atr, // camera orientation relative to world frame
			final QuadCLT     scene,
			final QuadCLT     ref_scene, // now - may be null - for testing if scene is rotated ref
			final boolean     toRGB,
			final boolean     show_nan,
			String            suffix,
			int               threadsMax,
			final int         debugLevel){
		return renderGPUFromDSI(
				sensor_mask,
				merge_channels,
				full_woi_in,      // show larger than sensor WOI in tiles (or null)
				clt_parameters,
				disparity_ref,
				// motion blur compensation 
				0.0, // double            mb_tau,      // 0.008; // time constant, sec
				0.0, // mb_max_gain, // 5.0;   // motion blur maximal gain (if more - move second point more than a pixel
				null, // double [][]       mb_vectors,  //
				
				scene_xyz, // camera center in world coordinates
				scene_atr, // camera orientation relative to world frame
				scene,
				ref_scene, // now - may be null - for testing if scene is rotated ref
				toRGB,
				show_nan,
				suffix,
				threadsMax,
				debugLevel);
	}	
	
	/**
	 * Making scene_xyz==null or scene_atr==null to use uniform grid (for rendering raw images)
	 * @param sensor_mask
	 * @param merge_channels
	 * @param full_woi_in
	 * @param clt_parameters
	 * @param disparity_ref
	 * @param mb_tau
	 * @param mb_max_gain
	 * @param mb_vectors
	 * @param scene_xyz
	 * @param scene_atr
	 * @param scene
	 * @param ref_scene
	 * @param toRGB
	 * @param show_nan
	 * @param suffix
	 * @param threadsMax
	 * @param debugLevel
	 * @return
	 */
	
	public static ImagePlus renderGPUFromDSI(
			final int         sensor_mask,
			final boolean     merge_channels,
			final Rectangle   full_woi_in,      // show larger than sensor WOI in tiles (or null)
			CLTParameters     clt_parameters,
			double []         disparity_ref,
			// motion blur compensation 
			double            mb_tau,      // 0.008; // time constant, sec
			double            mb_max_gain, // 5.0;   // motion blur maximal gain (if more - move second point more than a pixel
			double [][]       mb_vectors,  // now [2][ntiles];
			
			double []         scene_xyz, // camera center in world coordinates. If null - no shift, no ers
			double []         scene_atr, // camera orientation relative to world frame
			final QuadCLT     scene,
			final QuadCLT     ref_scene, // now - may be null - for testing if scene is rotated ref
			final boolean     toRGB,
			final boolean     show_nan,
			String            suffix,
			int               threadsMax,
			final int         debugLevel){
		double [][] pXpYD;
		if ((scene_xyz == null) || (scene_atr == null)) {
			scene_xyz = new double[3];
			scene_atr = new double[3];
			pXpYD=OpticalFlow.transformToScenePxPyD( // now should work with offset ref_scene
					full_woi_in,   // final Rectangle [] extra_woi,    // show larger than sensor WOI (or null)
					disparity_ref, // final double []   disparity_ref, // invalid tiles - NaN in disparity
					scene_xyz,     // final double []   scene_xyz, // camera center in world coordinates
					scene_atr,     // final double []   scene_atr, // camera orientation relative to world frame
					ref_scene,     // final QuadCLT     scene_QuadClt,
					ref_scene,     // final QuadCLT     reference_QuadClt, // now - may be null - for testing if scene is rotated ref
					threadsMax);   // int               threadsMax)
		} else {
			pXpYD=OpticalFlow.transformToScenePxPyD( // now should work with offset ref_scene
					full_woi_in,   // final Rectangle [] extra_woi,    // show larger than sensor WOI (or null)
					disparity_ref, // final double []   disparity_ref, // invalid tiles - NaN in disparity
					scene_xyz,     // final double []   scene_xyz, // camera center in world coordinates
					scene_atr,     // final double []   scene_atr, // camera orientation relative to world frame
					scene,         // final QuadCLT     scene_QuadClt,
					ref_scene,     // final QuadCLT     reference_QuadClt, // now - may be null - for testing if scene is rotated ref
					threadsMax);   // int               threadsMax)
		}
		int rendered_width = scene.getErsCorrection().getSensorWH()[0];
		if (full_woi_in != null) {
			rendered_width = full_woi_in.width * GPUTileProcessor.DTT_SIZE;
		}
		boolean showPxPyD = false;
		if (showPxPyD) {
			int dbg_width = rendered_width/GPUTileProcessor.DTT_SIZE;
			int dbg_height = pXpYD.length/dbg_width;
			double [][] dbg_img = new double [3 + ((mb_vectors!=null)? 2:0)][pXpYD.length];
			String [] dbg_titles = (mb_vectors!=null)?
					(new String[] {"pX","pY","Disparity","mb_X","mb_Y"}):
						(new String[] {"pX","pY","Disparity"});
			for (int i = 0; i < dbg_img.length; i++) {
				Arrays.fill(dbg_img[i], Double.NaN);
			}
			for (int nTile = 0; nTile < pXpYD.length; nTile++){
				if (pXpYD[nTile] != null) {
					for (int i = 0; i < pXpYD[nTile].length; i++) {
						dbg_img[i][nTile] = pXpYD[nTile][i];
					}
				}
				if (mb_vectors!=null) {
					for (int i = 0; i <2; i++) {
						dbg_img[3 + i][nTile] =  mb_tau * mb_vectors[i][nTile];
					}
				}
			}
			ShowDoubleFloatArrays.showArrays( // out of boundary 15
					dbg_img,
					dbg_width,
					dbg_height,
					true,
					scene.getImageName()+"-pXpYD",
					dbg_titles);
		}
		TpTask[][] tp_tasks;
		if (mb_vectors!=null) {
			tp_tasks = GpuQuad.setInterTasksMotionBlur( // "true" reference, with stereo actual reference will be offset
					scene.getNumSensors(),
					rendered_width,           // should match output size, pXpYD.length
					!scene.hasGPU(),          // final boolean             calcPortsCoordinatesAndDerivatives, // GPU can calculate them centreXY
					pXpYD,                    // final double [][]         pXpYD, // per-tile array of pX,pY,disparity triplets (or nulls)
					null,                     // final boolean []          selection, // may be null, if not null do not  process unselected tiles
					// motion blur compensation 
					mb_tau,                   // final double              mb_tau,      // 0.008; // time constant, sec
					mb_max_gain,              // final double              mb_max_gain, // 5.0;   // motion blur maximal gain (if more - move second point more than a pixel
					mb_vectors,               //final double [][]         mb_vectors,  //
					scene.getErsCorrection(), // final GeometryCorrection  geometryCorrection,
					0.0,                      // final double              disparity_corr,
					-1, // 0, // margin,      // final int                 margin,      // do not use tiles if their centers are closer to the edges
					null,                     // final boolean []          valid_tiles,            
					threadsMax);              // final int                 threadsMax)  // maximal number of threads to launch
		} else {
			tp_tasks = new TpTask[1][];
			tp_tasks[0] =  GpuQuad.setInterTasks( // "true" reference, with stereo actual reference will be offset
					scene.getNumSensors(),
					rendered_width,           // should match output size, pXpYD.length
					!scene.hasGPU(),          // final boolean             calcPortsCoordinatesAndDerivatives, // GPU can calculate them centreXY
					pXpYD,                    // final double [][]         pXpYD, // per-tile array of pX,pY,disparity triplets (or nulls)
					null,                     // final boolean []          selection, // may be null, if not null do not  process unselected tiles
					scene.getErsCorrection(), // final GeometryCorrection  geometryCorrection,
					0.0,                      // final double              disparity_corr,
					-1, // 0, // margin,      // final int                 margin,      // do not use tiles if their centers are closer to the edges
					null,                     // final boolean []          valid_tiles,            
					threadsMax);              // final int                 threadsMax)  // maximal number of threads to launch
		}
	    scene.saveQuadClt(); // to re-load new set of Bayer images to the GPU (do nothing for CPU) and Geometry
	    ImageDtt image_dtt = new ImageDtt(
	    		scene.getNumSensors(),
	    		clt_parameters.transform_size,
	    		clt_parameters.img_dtt,
	    		scene.isAux(),
	    		scene.isMonochrome(),
	    		scene.isLwir(),
	    		clt_parameters.getScaleStrength(scene.isAux()),
	    		scene.getGPU());
	    boolean use_reference = false;
	    int [] wh = (full_woi_in == null)? null: new int[]{
	    		full_woi_in.width * GPUTileProcessor.DTT_SIZE,
	    		full_woi_in.height * GPUTileProcessor.DTT_SIZE};
	    int                 erase_clt = show_nan ? 1:0;
//	    boolean test1 = true;
	    if (mb_vectors!=null) {// && test1) {
	    	image_dtt.setReferenceTDMotionBlur( // change to main?
	    			erase_clt, //final int                 erase_clt,
	    			wh, // null,                       // final int []              wh,               // null (use sensor dimensions) or pair {width, height} in pixels
	    			clt_parameters.img_dtt,     // final ImageDttParameters  imgdtt_params,    // Now just extra correlation parameters, later will include, most others
	    			use_reference, // true, // final boolean             use_reference_buffer,
	    			tp_tasks,               // final TpTask[]            tp_tasks,
	    			clt_parameters.gpu_sigma_r, // final double              gpu_sigma_r,     // 0.9, 1.1
	    			clt_parameters.gpu_sigma_b, // final double              gpu_sigma_b,     // 0.9, 1.1
	    			clt_parameters.gpu_sigma_g, // final double              gpu_sigma_g,     // 0.6, 0.7
	    			clt_parameters.gpu_sigma_m, // final double              gpu_sigma_m,     //  =       0.4; // 0.7;
	    			threadsMax,                 // final int                 threadsMax,       // maximal number of threads to launch
	    			debugLevel);                // final int                 globalDebugLevel);
	    } else {
	    	image_dtt.setReferenceTD( // change to main?
	    			erase_clt, //final int                 erase_clt,
	    			wh, // null,                       // final int []              wh,               // null (use sensor dimensions) or pair {width, height} in pixels
	    			clt_parameters.img_dtt,     // final ImageDttParameters  imgdtt_params,    // Now just extra correlation parameters, later will include, most others
	    			use_reference, // true, // final boolean             use_reference_buffer,
	    			tp_tasks[0],               // final TpTask[]            tp_tasks,
	    			clt_parameters.gpu_sigma_r, // final double              gpu_sigma_r,     // 0.9, 1.1
	    			clt_parameters.gpu_sigma_b, // final double              gpu_sigma_b,     // 0.9, 1.1
	    			clt_parameters.gpu_sigma_g, // final double              gpu_sigma_g,     // 0.6, 0.7
	    			clt_parameters.gpu_sigma_m, // final double              gpu_sigma_m,     //  =       0.4; // 0.7;
	    			threadsMax,                 // final int                 threadsMax,       // maximal number of threads to launch
	    			debugLevel);                // final int                 globalDebugLevel);
	    }
        ImagePlus imp_render = scene.renderFromTD (
        		sensor_mask,                  // final int         sensor_mask,
        		merge_channels,               // boolean             merge_channels,
                clt_parameters,                                 // CLTParameters clt_parameters,
                clt_parameters.getColorProcParameters(scene.isAux()), //ColorProcParameters colorProcParameters,
                clt_parameters.getRGBParameters(),              //EyesisCorrectionParameters.RGBParameters rgbParameters,\
                wh,                                             // null, // int []  wh,
                toRGB,                                          // boolean toRGB,
                use_reference,                                  // boolean use_reference
                suffix);                                        // String  suffix)
		return imp_render;
	}
	
	/**
	 * Prepare 16x16 texture tiles using GPU from disparity reference. Includes motion blur correction
	 * Does not use scene.saveQuadClt() to re-load new set of Bayer images to the GPU -- should be done by caller
	 * @param clt_parameters processing parameters
	 * @param disparity_ref disparity scene reference
	 * @param mb_tau     	thermal constant
	 * @param mb_max_gain   maximal gain when correcting MB
	 * @param mb_vectors    motion blur vectors {vx[tiles], vy[tiles]}
	 * @param scene_xyz     scene offset X,Y,Z in meters
	 * @param scene_atr     scene azimuth, tilt, roll in radians
	 * @param scene         scene QuadCLT instance
	 * @param ref_scene     reference scene QuadCLT instance
	 * @param filter_bg     remove BG tiles (possibly occluded)
	 * @param max_distortion maximal neighbor tiles offset as a fraction of tile size (8) 
	 * @param cluster_index which cluster the tile belongs - to verify tile distortions
	 * @param border        border tiles, may have neighbors (other border tiles) over discontinuity
	 * @param keep_channels output per-channel Y(G) after YA (RBGA) 
	 * @param debugLevel
	 * @return
	 */
	public static double [][][][] texturesGPUFromDSI(
			CLTParameters     clt_parameters,
			double []         disparity_ref,
			// motion blur compensation 
			double            mb_tau,      // 0.008; // time constant, sec
			double            mb_max_gain, // 5.0;   // motion blur maximal gain (if more - move second point more than a pixel
			double [][]       mb_vectors,  // now [2][ntiles];
			final double []   scene_xyz, // camera center in world coordinates
			final double []   scene_atr, // camera orientation relative to world frame
			final QuadCLT     scene,
			final QuadCLT     ref_scene, // now - may be null - for testing if scene is rotated ref
			final boolean     filter_bg, // remove bg tiles (possibly occluded)
			final double      max_distortion, // maximal neighbor tiles offset as a fraction of tile size (8) 
			final int []      cluster_index,  // 
			final boolean []  border, // border tiles
			final boolean     keep_channels,
			final int         debugLevel){
		// FIXME: Move to clt_parameters;
		final double max_overlap =    0.6; 
		final double min_adisp_cam =  0.2;
		final double min_rdisp_cam =  0.03;


		double [][] pXpYD_prefilter =OpticalFlow.transformToScenePxPyD( // now should work with offset ref_scene
				null,          // full_woi_in,   // final Rectangle [] extra_woi,    // show larger than sensor WOI (or null)
				disparity_ref, // final double []   disparity_ref, // invalid tiles - NaN in disparity
				scene_xyz,     // final double []   scene_xyz, // camera center in world coordinates
				scene_atr,     // final double []   scene_atr, // camera orientation relative to world frame
				scene,         // final QuadCLT     scene_QuadClt,
				ref_scene,     // final QuadCLT     reference_QuadClt, // now - may be null - for testing if scene is rotated ref
				OpticalFlow.THREADS_MAX);   // int               threadsMax)
		double [][] pXpYD;
		if (filter_bg) {
			double [][] scene_ds = OpticalFlow.conditionInitialDS(
					clt_parameters, // CLTParameters  clt_parameters,
					scene, // QuadCLT        scene,
					-1); // int debug_level);
			if (scene_ds != null) {
				double [] disparity_cam = scene_ds[0]; // null; // for now
				pXpYD = OpticalFlow.filterBG (
						ref_scene.getTileProcessor(), // final TileProcessor tp,
						pXpYD_prefilter,       // final double [][] pXpYD,
						max_overlap,           // final double max_overlap,
						null, // disparity_cam,         // final double [] disparity_cam,
						min_adisp_cam,         // final double min_adisp_cam,
						min_rdisp_cam,         // final double min_rdisp_cam,
						clt_parameters.tileX,  // final int    dbg_tileX,
						clt_parameters.tileY,  // final int    dbg_tileY,
						0); // 1); //debug_level);          // final int    debug_level);
			} else {
				pXpYD = pXpYD_prefilter;
			}
		} else {
			pXpYD = pXpYD_prefilter;
		}

		int rendered_width = scene.getErsCorrection().getSensorWH()[0];
		boolean showPxPyD = false;
		if (showPxPyD) {
			int dbg_width = rendered_width/GPUTileProcessor.DTT_SIZE;
			int dbg_height = pXpYD.length/dbg_width;
			double [][] dbg_img = new double [3 + ((mb_vectors!=null)? 2:0)][pXpYD.length];
			String [] dbg_titles = (mb_vectors!=null)?
					(new String[] {"pX","pY","Disparity","mb_X","mb_Y"}):
						(new String[] {"pX","pY","Disparity"});
			for (int i = 0; i < dbg_img.length; i++) {
				Arrays.fill(dbg_img[i], Double.NaN);
			}
			for (int nTile = 0; nTile < pXpYD.length; nTile++){
				if (pXpYD[nTile] != null) {
					for (int i = 0; i < pXpYD[nTile].length; i++) {
						dbg_img[i][nTile] = pXpYD[nTile][i];
					}
				}
				if (mb_vectors!=null) {
					for (int i = 0; i <2; i++) {
						dbg_img[3 + i][nTile] =  mb_tau * mb_vectors[i][nTile];
					}
				}
			}
			ShowDoubleFloatArrays.showArrays( // out of boundary 15
					dbg_img,
					dbg_width,
					dbg_height,
					true,
					scene.getImageName()+"-pXpYD",
					dbg_titles);
		}
		TpTask[][] tp_tasks;
		if (mb_vectors!=null) {
			tp_tasks = GpuQuad.setInterTasksMotionBlur( // "true" reference, with stereo actual reference will be offset
					scene.getNumSensors(),
					rendered_width,           // should match output size, pXpYD.length
					!scene.hasGPU(),          // final boolean             calcPortsCoordinatesAndDerivatives, // GPU can calculate them centreXY
					pXpYD,                    // final double [][]         pXpYD, // per-tile array of pX,pY,disparity triplets (or nulls)
					null,                     // final boolean []          selection, // may be null, if not null do not  process unselected tiles
					// motion blur compensation 
					mb_tau,                   // final double              mb_tau,      // 0.008; // time constant, sec
					mb_max_gain,              // final double              mb_max_gain, // 5.0;   // motion blur maximal gain (if more - move second point more than a pixel
					mb_vectors,               //final double [][]         mb_vectors,  //
					scene.getErsCorrection(), // final GeometryCorrection  geometryCorrection,
					0.0,                      // final double              disparity_corr,
					-1, // 0, // margin,      // final int                 margin,      // do not use tiles if their centers are closer to the edges
					null,                     // final boolean []          valid_tiles,            
					OpticalFlow.THREADS_MAX); // final int                 threadsMax)  // maximal number of threads to launch
		} else {
			tp_tasks = new TpTask[1][];
			tp_tasks[0] =  GpuQuad.setInterTasks( // "true" reference, with stereo actual reference will be offset
					scene.getNumSensors(),
					rendered_width,           // should match output size, pXpYD.length
					!scene.hasGPU(),          // final boolean             calcPortsCoordinatesAndDerivatives, // GPU can calculate them centreXY
					pXpYD,                    // final double [][]         pXpYD, // per-tile array of pX,pY,disparity triplets (or nulls)
					null,                     // final boolean []          selection, // may be null, if not null do not  process unselected tiles
					scene.getErsCorrection(), // final GeometryCorrection  geometryCorrection,
					0.0,                      // final double              disparity_corr,
					-1, // 0, // margin,      // final int                 margin,      // do not use tiles if their centers are closer to the edges
					null,                     // final boolean []          valid_tiles,            
					OpticalFlow.THREADS_MAX);              // final int                 threadsMax)  // maximal number of threads to launch
		}
		if (tp_tasks[0].length == 0) {
			if (debugLevel > -1) {
				System.out.println("texturesGPUFromDSI(): no tiles to process");
				
			}
			return null;
		}
		
		///		scene.saveQuadClt(); // to re-load new set of Bayer images to the GPU (do nothing for CPU) and Geometry
		ImageDtt image_dtt = new ImageDtt(
				scene.getNumSensors(),
				clt_parameters.transform_size,
				clt_parameters.img_dtt,
				scene.isAux(),
				scene.isMonochrome(),
				scene.isLwir(),
				clt_parameters.getScaleStrength(scene.isAux()),
				scene.getGPU());
		boolean use_reference = false;
		int                 erase_clt = 0; // show_nan ? 1:0;
		if (mb_vectors!=null) {// && test1) {
			image_dtt.setReferenceTDMotionBlur( // change to main?
					erase_clt, //final int                 erase_clt,
					null, // wh, // null,                       // final int []              wh,               // null (use sensor dimensions) or pair {width, height} in pixels
					clt_parameters.img_dtt,     // final ImageDttParameters  imgdtt_params,    // Now just extra correlation parameters, later will include, most others
					use_reference, // true, // final boolean             use_reference_buffer,
					tp_tasks,               // final TpTask[]            tp_tasks,
					clt_parameters.gpu_sigma_r, // final double              gpu_sigma_r,     // 0.9, 1.1
					clt_parameters.gpu_sigma_b, // final double              gpu_sigma_b,     // 0.9, 1.1
					clt_parameters.gpu_sigma_g, // final double              gpu_sigma_g,     // 0.6, 0.7
					clt_parameters.gpu_sigma_m, // final double              gpu_sigma_m,     //  =       0.4; // 0.7;
					OpticalFlow.THREADS_MAX,    // final int                 threadsMax,       // maximal number of threads to launch
					debugLevel);                // final int                 globalDebugLevel);
		} else {
			image_dtt.setReferenceTD( // change to main?
					erase_clt, //final int                 erase_clt,
					null, // wh, // null,                       // final int []              wh,               // null (use sensor dimensions) or pair {width, height} in pixels
					clt_parameters.img_dtt,     // final ImageDttParameters  imgdtt_params,    // Now just extra correlation parameters, later will include, most others
					use_reference, // true, // final boolean             use_reference_buffer,
					tp_tasks[0],               // final TpTask[]            tp_tasks,
					clt_parameters.gpu_sigma_r, // final double              gpu_sigma_r,     // 0.9, 1.1
					clt_parameters.gpu_sigma_b, // final double              gpu_sigma_b,     // 0.9, 1.1
					clt_parameters.gpu_sigma_g, // final double              gpu_sigma_g,     // 0.6, 0.7
					clt_parameters.gpu_sigma_m, // final double              gpu_sigma_m,     //  =       0.4; // 0.7;
					OpticalFlow.THREADS_MAX,    // final int                 threadsMax,       // maximal number of threads to launch
					debugLevel);                // final int                 globalDebugLevel);
		}
		// now obeys keep_weights 
		final double [][][][] texture_tiles =  image_dtt.process_texture_tiles( // from transform domain all sensors to combined texture
				clt_parameters.corr_red,       // double    corr_red,
				clt_parameters.corr_blue,      // double    corr_blue,
				clt_parameters.min_shot,       // double    min_shot,           // 10.0
				clt_parameters.scale_shot,     // double    scale_shot,         // 3.0
				clt_parameters.diff_sigma,     // double    diff_sigma,         // pixel value/pixel change Used much larger sigma = 10.0 instead of 1.5
				clt_parameters.diff_threshold, // double    diff_threshold,     // pixel value/pixel change
				clt_parameters.min_agree,      // double    min_agree,          // minimal number of channels to agree on a point (real number to work with fuzzy averages)
				clt_parameters.dust_remove,    // boolean   dust_remove
				0); // int       keep_weights       // 2 bits now, move to parameters


		if (max_distortion > 0) {//  remove distorted tiles
			double max_distortion2 = max_distortion * max_distortion; 
			final int tilesX =           ref_scene.getTileProcessor().getTilesX();
			final int tilesY =           ref_scene.getTileProcessor().getTilesY();
			final int tiles =            tilesX*tilesY;
			int tile_size =              ref_scene.getTileProcessor().getTileSize();
			final TileNeibs tn =         new TileNeibs(tilesX, tilesY);
			final boolean [] distorted = new boolean [tiles];
			final double [][][] tileXY = new double [tilesY][tilesX][];
			final TpTask[][] ftp_tasks = tp_tasks;
			final Thread[] threads =     ImageDtt.newThreadArray(OpticalFlow.THREADS_MAX);
			final AtomicInteger ai =     new AtomicInteger(0);
			for (int ithread = 0; ithread < threads.length; ithread++) {
				threads[ithread] = new Thread() {
					public void run() {
						for (int itile = ai.getAndIncrement(); itile < tp_tasks[0].length; itile = ai.getAndIncrement()) {
							TpTask task = ftp_tasks[0][itile];
							if (task != null) {
								tileXY[task.getTileY()][task.getTileX()] = task.getDoubleCenterXY();
							}
						}
					}
				};
			}		      
			ImageDtt.startAndJoin(threads);
			ai.set(0);
			final double [] dbg_distort = (debugLevel>2)? (new double [tiles]):null;
			for (int ithread = 0; ithread < threads.length; ithread++) {
				threads[ithread] = new Thread() {
					public void run() {
						for (int nTile = ai.getAndIncrement(); nTile < tiles; nTile = ai.getAndIncrement()) {
							int tileX = nTile % tilesX;
							int tileY = nTile / tilesX;
							if (texture_tiles[tileY][tileX] != null) {
								double [] centerXY =tileXY[tileY][tileX]; 
								if ((centerXY != null) && (cluster_index[nTile] >= 0)){
									// see if any tile
									for (int dir = 0; dir < TileNeibs.DIRS; dir++) {
										int tile1 = tn.getNeibIndex(nTile, dir);
										if (tile1 >= 0) {
											if (cluster_index[tile1] == cluster_index[nTile]) {
												int tileX1 = tn.getX(tile1);
												int tileY1 = tn.getY(tile1);
												double [] thisXY = tileXY[tileY1][tileX1];
												// border-border neighbor may have discontinuity even belonging to the same cluster
												// (coming close AGAIN)
												if ((thisXY != null) && !(border[nTile] && border[tile1])) {
													double dx = thisXY[0] - centerXY[0] - tile_size * TileNeibs.getDX(dir);
													double dy = thisXY[1] - centerXY[1] - tile_size * TileNeibs.getDY(dir);
													double e2 = dx*dx+dy*dy;
													if (dbg_distort != null) {
														if (e2 > dbg_distort[nTile]) {
															dbg_distort[nTile] = Math.sqrt(e2);
														}
													} else {
														if (e2 > max_distortion2) {
															distorted[nTile] = true;
															break;
														}
													}
												} else {
													continue; // check why task is null here for >=0 cluster index 
												}
											}
										}
									}
								} else {
									if (debugLevel > -3) {
										System.out.println("Non-null texture for no-cluster, nTile="+nTile+
												", tileX="+tileX+", tileY="+tileY+", cluster_index["+nTile+"]="+cluster_index[nTile]);
									}
								}
							}
						}
					}
				};
			}		      
			ImageDtt.startAndJoin(threads);
			ai.set(0);
			if (dbg_distort != null) {
				ShowDoubleFloatArrays.showArrays( // out of boundary 15
						dbg_distort,
						tilesX,
						tilesY,
						"test-distort");

			}


			for (int ithread = 0; ithread < threads.length; ithread++) {
				threads[ithread] = new Thread() {
					public void run() {
						for (int nTile = ai.getAndIncrement(); nTile < distorted.length; nTile = ai.getAndIncrement()) if (distorted[nTile]){
							int tileX = nTile % tilesX;
							int tileY = nTile / tilesX;
							if (texture_tiles[tileY] != null) {
								texture_tiles[tileY][tileX] = null;
							}
						}
					}
				};
			}		      
			ImageDtt.startAndJoin(threads);
		}
		return texture_tiles;
	}
	
	public static double [][][][] texturesNoOverlapGPUFromDSI(
			CLTParameters     clt_parameters,
			double []         disparity_ref,
			// motion blur compensation 
			double            mb_tau,      // 0.008; // time constant, sec
			double            mb_max_gain, // 5.0;   // motion blur maximal gain (if more - move second point more than a pixel
			double [][]       mb_vectors,  // now [2][ntiles];
			final double []   scene_xyz, // camera center in world coordinates
			final double []   scene_atr, // camera orientation relative to world frame
			final QuadCLT     scene,
			final QuadCLT     ref_scene, // now - may be null - for testing if scene is rotated ref
			final boolean     filter_bg, // remove bg tiles (possibly occluded)
			final double      max_distortion, // maximal neighbor tiles offset as a fraction of tile size (8) 
			final int []      cluster_index,  // [tilesX*tilesY]
			final boolean []  border, // border tiles
			final boolean     keep_channels,
			final int         debugLevel){
		// FIXME: Move to clt_parameters;
		final double max_overlap =    0.6; 
		final double min_adisp_cam =  0.2;
		final double min_rdisp_cam =  0.03;
		int  keep_weights = 2; // (clt_parameters.replace_weights? 2 : 0); // just 0/2


		double [][] pXpYD_prefilter =OpticalFlow.transformToScenePxPyD( // now should work with offset ref_scene
				null,          // full_woi_in,   // final Rectangle [] extra_woi,    // show larger than sensor WOI (or null)
				disparity_ref, // final double []   disparity_ref, // invalid tiles - NaN in disparity
				scene_xyz,     // final double []   scene_xyz, // camera center in world coordinates
				scene_atr,     // final double []   scene_atr, // camera orientation relative to world frame
				scene,         // final QuadCLT     scene_QuadClt,
				ref_scene,     // final QuadCLT     reference_QuadClt, // now - may be null - for testing if scene is rotated ref
				OpticalFlow.THREADS_MAX);   // int               threadsMax)
		double [][] pXpYD;
		if (filter_bg) {
			double [][] scene_ds = OpticalFlow.conditionInitialDS(
					clt_parameters, // CLTParameters  clt_parameters,
					scene, // QuadCLT        scene,
					-1); // int debug_level);
			if (scene_ds != null) {
				double [] disparity_cam = scene_ds[0]; // null; // for now
				pXpYD = OpticalFlow.filterBG (
						ref_scene.getTileProcessor(), // final TileProcessor tp,
						pXpYD_prefilter,       // final double [][] pXpYD,
						max_overlap,           // final double max_overlap,
						null, // disparity_cam,         // final double [] disparity_cam,
						min_adisp_cam,         // final double min_adisp_cam,
						min_rdisp_cam,         // final double min_rdisp_cam,
						clt_parameters.tileX,  // final int    dbg_tileX,
						clt_parameters.tileY,  // final int    dbg_tileY,
						0); // 1); //debug_level);          // final int    debug_level);
			} else {
				pXpYD = pXpYD_prefilter;
			}
		} else {
			pXpYD = pXpYD_prefilter;
		}

		int rendered_width = scene.getErsCorrection().getSensorWH()[0];
		boolean showPxPyD = false;
		if (showPxPyD) {
			int dbg_width = rendered_width/GPUTileProcessor.DTT_SIZE;
			int dbg_height = pXpYD.length/dbg_width;
			double [][] dbg_img = new double [3 + ((mb_vectors!=null)? 2:0)][pXpYD.length];
			String [] dbg_titles = (mb_vectors!=null)?
					(new String[] {"pX","pY","Disparity","mb_X","mb_Y"}):
						(new String[] {"pX","pY","Disparity"});
			for (int i = 0; i < dbg_img.length; i++) {
				Arrays.fill(dbg_img[i], Double.NaN);
			}
			for (int nTile = 0; nTile < pXpYD.length; nTile++){
				if (pXpYD[nTile] != null) {
					for (int i = 0; i < pXpYD[nTile].length; i++) {
						dbg_img[i][nTile] = pXpYD[nTile][i];
					}
				}
				if (mb_vectors!=null) {
					for (int i = 0; i <2; i++) {
						dbg_img[3 + i][nTile] =  mb_tau * mb_vectors[i][nTile];
					}
				}
			}
			ShowDoubleFloatArrays.showArrays( // out of boundary 15
					dbg_img,
					dbg_width,
					dbg_height,
					true,
					scene.getImageName()+"-pXpYD",
					dbg_titles);
		}
		TpTask[][] tp_tasks;
		if (mb_vectors!=null) {
			tp_tasks = GpuQuad.setInterTasksMotionBlur( // "true" reference, with stereo actual reference will be offset
					scene.getNumSensors(),
					rendered_width,           // should match output size, pXpYD.length
					!scene.hasGPU(),          // final boolean             calcPortsCoordinatesAndDerivatives, // GPU can calculate them centreXY
					pXpYD,                    // final double [][]         pXpYD, // per-tile array of pX,pY,disparity triplets (or nulls)
					null,                     // final boolean []          selection, // may be null, if not null do not  process unselected tiles
					// motion blur compensation 
					mb_tau,                   // final double              mb_tau,      // 0.008; // time constant, sec
					mb_max_gain,              // final double              mb_max_gain, // 5.0;   // motion blur maximal gain (if more - move second point more than a pixel
					mb_vectors,               //final double [][]         mb_vectors,  //
					scene.getErsCorrection(), // final GeometryCorrection  geometryCorrection,
					0.0,                      // final double              disparity_corr,
					-1, // 0, // margin,      // final int                 margin,      // do not use tiles if their centers are closer to the edges
					null,                     // final boolean []          valid_tiles,            
					OpticalFlow.THREADS_MAX); // final int                 threadsMax)  // maximal number of threads to launch
		} else {
			tp_tasks = new TpTask[1][];
			tp_tasks[0] =  GpuQuad.setInterTasks( // "true" reference, with stereo actual reference will be offset
					scene.getNumSensors(),
					rendered_width,           // should match output size, pXpYD.length
					!scene.hasGPU(),          // final boolean             calcPortsCoordinatesAndDerivatives, // GPU can calculate them centreXY
					pXpYD,                    // final double [][]         pXpYD, // per-tile array of pX,pY,disparity triplets (or nulls)
					null,                     // final boolean []          selection, // may be null, if not null do not  process unselected tiles
					scene.getErsCorrection(), // final GeometryCorrection  geometryCorrection,
					0.0,                      // final double              disparity_corr,
					-1, // 0, // margin,      // final int                 margin,      // do not use tiles if their centers are closer to the edges
					null,                     // final boolean []          valid_tiles,            
					OpticalFlow.THREADS_MAX);              // final int                 threadsMax)  // maximal number of threads to launch
		}
		if (tp_tasks[0].length == 0) {
			if (debugLevel > -1) {
				System.out.println("texturesGPUFromDSI(): no tiles to process");
				
			}
			return null;
		}
		
		///		scene.saveQuadClt(); // to re-load new set of Bayer images to the GPU (do nothing for CPU) and Geometry
		ImageDtt image_dtt = new ImageDtt(
				scene.getNumSensors(),
				clt_parameters.transform_size,
				clt_parameters.img_dtt,
				scene.isAux(),
				scene.isMonochrome(),
				scene.isLwir(),
				clt_parameters.getScaleStrength(scene.isAux()),
				scene.getGPU());
		boolean use_reference = false;
		int                 erase_clt = 0; // show_nan ? 1:0;
		if (mb_vectors!=null) {// && test1) {
			image_dtt.setReferenceTDMotionBlur( // change to main?
					erase_clt, //final int                 erase_clt,
					null, // wh, // null,                       // final int []              wh,               // null (use sensor dimensions) or pair {width, height} in pixels
					clt_parameters.img_dtt,     // final ImageDttParameters  imgdtt_params,    // Now just extra correlation parameters, later will include, most others
					use_reference, // true, // final boolean             use_reference_buffer,
					tp_tasks,               // final TpTask[]            tp_tasks,
					clt_parameters.gpu_sigma_r, // final double              gpu_sigma_r,     // 0.9, 1.1
					clt_parameters.gpu_sigma_b, // final double              gpu_sigma_b,     // 0.9, 1.1
					clt_parameters.gpu_sigma_g, // final double              gpu_sigma_g,     // 0.6, 0.7
					clt_parameters.gpu_sigma_m, // final double              gpu_sigma_m,     //  =       0.4; // 0.7;
					OpticalFlow.THREADS_MAX,    // final int                 threadsMax,       // maximal number of threads to launch
					debugLevel);                // final int                 globalDebugLevel);
		} else {
			image_dtt.setReferenceTD( // change to main?
					erase_clt, //final int                 erase_clt,
					null, // wh, // null,                       // final int []              wh,               // null (use sensor dimensions) or pair {width, height} in pixels
					clt_parameters.img_dtt,     // final ImageDttParameters  imgdtt_params,    // Now just extra correlation parameters, later will include, most others
					use_reference, // true, // final boolean             use_reference_buffer,
					tp_tasks[0],               // final TpTask[]            tp_tasks,
					clt_parameters.gpu_sigma_r, // final double              gpu_sigma_r,     // 0.9, 1.1
					clt_parameters.gpu_sigma_b, // final double              gpu_sigma_b,     // 0.9, 1.1
					clt_parameters.gpu_sigma_g, // final double              gpu_sigma_g,     // 0.6, 0.7
					clt_parameters.gpu_sigma_m, // final double              gpu_sigma_m,     //  =       0.4; // 0.7;
					OpticalFlow.THREADS_MAX,    // final int                 threadsMax,       // maximal number of threads to launch
					debugLevel);                // final int                 globalDebugLevel);
		}
		// now obeys keep_weights
		Rectangle woi = new Rectangle(); // will be filled out to match actual available image
//		int       keep_weights = (clt_parameters.keep_weights? 1 : 0) + (clt_parameters.replace_weights? 2 : 0);
//		int  keep_weights = (clt_parameters.replace_weights? 2 : 0); // just 0/2
		int  text_colors = (scene.isMonochrome() ? 1 : 3);
		int  texture_layers = (text_colors + 1)+((keep_weights != 0)?(text_colors * scene.getNumSensors()):0);   
		double [] col_weights = new double[text_colors];
		if (text_colors < 3) {
			col_weights[0] = 1.0;
		} else {
			col_weights[2] = 1.0/(1.0 +  clt_parameters.corr_red + clt_parameters.corr_blue);    // green color
			col_weights[0] = clt_parameters.corr_red *  col_weights[2];
			col_weights[1] = clt_parameters.corr_blue * col_weights[2];
		}
		scene.gpuQuad.execRBGA(
				col_weights,                   // double [] color_weights,
				scene.isLwir(),                // boolean   is_lwir,
				clt_parameters.min_shot,       // double    min_shot,           // 10.0
				clt_parameters.scale_shot,     // double    scale_shot,         // 3.0
				clt_parameters.diff_sigma,     // double    diff_sigma,         // pixel value/pixel change Used much larger sigma = 10.0 instead of 1.5
				clt_parameters.diff_threshold, // double    diff_threshold,     // pixel value/pixel change
				clt_parameters.min_agree,      // double    min_agree,          // minimal number of channels to agree on a point (real number to work with fuzzy averages)
				clt_parameters.dust_remove,    // boolean   dust_remove,
				keep_weights);                            // int       keep_weights)
		float [][] rbga = scene.gpuQuad.getRBGA(
				texture_layers - 1, // (isMonochrome() ? 1 : 3), // int     num_colors,
				woi); // woi in pixels
		boolean show_rbga= false;
		if (show_rbga) {
			ShowDoubleFloatArrays.showArrays( // show slices RBGA (colors - 256, A - 1.0)
					rbga,
					woi.width,
					woi.height,
					true,
					scene.getImageName()+"-RGBA-STACK-D"+clt_parameters.disparity+
					":"+clt_parameters.gpu_woi_tx+":"+clt_parameters.gpu_woi_ty+
					":"+clt_parameters.gpu_woi_twidth+":"+clt_parameters.gpu_woi_theight+
					":"+(clt_parameters.gpu_woi_round?"C":"R")
					//,new String[] {"R","B","G","A"}
					);
			show_rbga = false;
		}
		
//		if (debugLevel > -1000) {
//			return null;
//		}
		final int tilesX =    scene.getTileProcessor().getTilesX();
		final int tilesY =    scene.getTileProcessor().getTilesY();
		final int tiles =     tilesX*tilesY;
		int tile_size =       scene.getTileProcessor().getTileSize();
		int tile_len =        tile_size*tile_size;
		final double [][][][] texture_tiles88 = new double [tilesY][tilesX][][]; // here - non-overlapped!
		// copy from windowed output to 
		final TpTask[] ftp_tasks =    tp_tasks[0];
		final Rectangle fwoi =        woi;
		final Thread[] threads =      ImageDtt.newThreadArray(OpticalFlow.THREADS_MAX);
		final AtomicInteger ai =      new AtomicInteger(0);
		final double [][][] tileXY =  new double [tilesY][tilesX][];
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				public void run() {
					for (int itile = ai.getAndIncrement(); itile < ftp_tasks.length; itile = ai.getAndIncrement()) {
						TpTask task = ftp_tasks[itile];
						if (task != null) {
							int x0 = tile_size * task.tx - fwoi.x;
							int y0 = tile_size * task.ty - fwoi.y;
							if ((x0 >=0) && (y0>=0) && (x0 < fwoi.width) && (y0 < fwoi.height)) {
								int sindx0 =  x0 + y0 * fwoi.width;
								texture_tiles88[task.ty][task.tx] = new double [rbga.length][tile_len];
								for (int nchn = 0; nchn < rbga.length; nchn++) {
									int indx=0;
									double [] tt = texture_tiles88[task.ty][task.tx][nchn];
									for (int row = 0; row <tile_size; row++) {
										int sindx = sindx0 + fwoi.width * row;
										for (int col = 0; col <tile_size; col++) {
											tt[indx++] = rbga[nchn][sindx++];
										}
									}
								}
							} else {
								System.out.println("texturesNoneoverlapGPUFromDSI() task tile outside of texture:\n"+
										"task.tx="+task.tx+", task.ty="+task.ty+", woi.x="+fwoi.x+", woi.y="+fwoi.y+
												", woi.width="+fwoi.width+", woi.height="+fwoi.height);
							}
							tileXY[task.ty][task.tx] = task.getDoubleCenterXY();
						}
					}
				}
			};
		}		      
		ImageDtt.startAndJoin(threads);

		if (max_distortion > 0) {//  remove distorted tiles
			double max_distortion2 = max_distortion * max_distortion; 
			final TileNeibs tn =         new TileNeibs(tilesX, tilesY);
			final boolean [] distorted = new boolean [tiles];
			ai.set(0);
			final double [] dbg_distort = (debugLevel>2)? (new double [tiles]):null;
			for (int ithread = 0; ithread < threads.length; ithread++) {
				threads[ithread] = new Thread() {
					public void run() {
						for (int nTile = ai.getAndIncrement(); nTile < tiles; nTile = ai.getAndIncrement()) {
							int tileX = nTile % tilesX;
							int tileY = nTile / tilesX;
							if (tileXY[tileY][tileX] != null) {
								double [] centerXY =tileXY[tileY][tileX]; 
								if ((centerXY != null) && (cluster_index[nTile] >= 0)){
									// see if any tile
									for (int dir = 0; dir < TileNeibs.DIRS; dir++) {
										int tile1 = tn.getNeibIndex(nTile, dir);
										if (tile1 >= 0) {
											if (cluster_index[tile1] == cluster_index[nTile]) {
												int tileX1 = tn.getX(tile1);
												int tileY1 = tn.getY(tile1);
												double [] thisXY = tileXY[tileY1][tileX1];
												// border-border neighbor may have discontinuity even belonging to the same cluster
												// (coming close AGAIN)
												if ((thisXY != null) && !(border[nTile] && border[tile1])) {
													double dx = thisXY[0] - centerXY[0] - tile_size * TileNeibs.getDX(dir);
													double dy = thisXY[1] - centerXY[1] - tile_size * TileNeibs.getDY(dir);
													double e2 = dx*dx+dy*dy;
													if (dbg_distort != null) {
														if (e2 > dbg_distort[nTile]) {
															dbg_distort[nTile] = Math.sqrt(e2);
														}
													} else {
														if (e2 > max_distortion2) {
															distorted[nTile] = true;
															break;
														}
													}
												} else {
													continue; // check why task is null here for >=0 cluster index 
												}
											}
										}
									}
								} else {
									if (debugLevel > -3) {
										System.out.println("Non-null texture for no-cluster, nTile="+nTile+
												", tileX="+tileX+", tileY="+tileY+", cluster_index["+nTile+"]="+cluster_index[nTile]);
									}
								}
							}
						}
					}
				};
			}		      
			ImageDtt.startAndJoin(threads);
			ai.set(0);
			if (dbg_distort != null) {
				ShowDoubleFloatArrays.showArrays( // out of boundary 15
						dbg_distort,
						tilesX,
						tilesY,
						"test-distort");

			}
			for (int ithread = 0; ithread < threads.length; ithread++) {
				threads[ithread] = new Thread() {
					public void run() {
						for (int nTile = ai.getAndIncrement(); nTile < distorted.length; nTile = ai.getAndIncrement()) if (distorted[nTile]){
							int tileX = nTile % tilesX;
							int tileY = nTile / tilesX;
							if (texture_tiles88[tileY] != null) {
								texture_tiles88[tileY][tileX] = null;
							}
						}
					}
				};
			}		      
			ImageDtt.startAndJoin(threads);
		}
		return texture_tiles88;
	}
	
	
	
	
	// add similar for textures
	public ImagePlus renderFromTD (
			int                 sensor_mask,
			boolean             merge_channels,
			CLTParameters       clt_parameters,
			ColorProcParameters colorProcParameters,
			EyesisCorrectionParameters.RGBParameters rgbParameters,
			int []              wh,
			boolean             toRGB,
			boolean             use_reference,
			String              suffix
			) {
		suffix+=(toRGB?"-COLOR":"-MONO");
        gpuQuad.execImcltRbgAll(
        		isMonochrome(),
        		use_reference,
				wh); //int [] wh
		// get data back from GPU
		final float [][][] iclt_fimg = new float [getNumSensors()][][];
		int nchn = 0;
		int ncol = 0;
		int nTiles = 0;
		for (int ncam = 0; ncam < iclt_fimg.length; ncam++) if (((1 << ncam) & sensor_mask) != 0){
			iclt_fimg[ncam] = gpuQuad.getRBG(ncam); // updated window
			ncol = iclt_fimg[ncam].length;
			nTiles = iclt_fimg[ncam][0].length;
			nchn++;
		}
		if (merge_channels) { // [0] - combo, other sensor channels - nulls
			final double scale = 1.0 / nchn;
			final float [][] iclt_fimg_combo = new float [ncol][nTiles];
			final Thread[] threads = ImageDtt.newThreadArray(THREADS_MAX);
			final AtomicInteger ai = new AtomicInteger(0);
			for (int ithread = 0; ithread < threads.length; ithread++) {
				threads[ithread] = new Thread() {
					public void run() {
						for (int nTile = ai.getAndIncrement(); nTile < iclt_fimg_combo[0].length; nTile = ai.getAndIncrement()) {
							for (int ncol = 0; ncol < iclt_fimg_combo.length; ncol++) {
								double d = 0;
								for (int i = 0; i < iclt_fimg.length; i++) if (iclt_fimg[i] != null) {
									d+=iclt_fimg[i][ncol][nTile];
								}
								iclt_fimg_combo[ncol][nTile] = (float) (d * scale);
							}
						}
					}
				};
			}		      
			ImageDtt.startAndJoin(threads);
			iclt_fimg[0] = iclt_fimg_combo;
			for (int i = 1; i < iclt_fimg.length; i++) {
				iclt_fimg[i] = null;
			}
		}
		
		// 2022/06/15 - handles variable window size
		int out_width =  gpuQuad.getImageWidth();//   + gpuQuad.getDttSize(); // 2022/05/12 removed margins from gpuQuad.getRBG(ncam);
		int out_height = gpuQuad.getImageHeight(); // + gpuQuad.getDttSize(); // 2022/05/12 removed margins from gpuQuad.getRBG(ncam);
		if (isLwir() && colorProcParameters.lwir_autorange) {
			double rel_low =  colorProcParameters.lwir_low;
			double rel_high = colorProcParameters.lwir_high;
			if (!Double.isNaN(getLwirOffset())) {
				rel_low -=  getLwirOffset();
				rel_high -= getLwirOffset();
			}
			double [] cold_hot =  autorange(
					iclt_fimg,                         // iclt_data, // double [][][] iclt_data, //  [iQuad][ncol][i] - normally only [][2][] is non-null
					rel_low,                           // double hard_cold,// matches data, DC (this.lwir_offset)  subtracted
					rel_high,                          // double hard_hot,   // matches data, DC (this.lwir_offset)  subtracted
					colorProcParameters.lwir_too_cold, // double too_cold, // pixels per image
					colorProcParameters.lwir_too_hot,  // double too_hot,  // pixels per image
					1024); // int num_bins)
			if (cold_hot != null) {
				if (!Double.isNaN(getLwirOffset())) {
					cold_hot[0] += getLwirOffset();
					cold_hot[1] += getLwirOffset();
				}
			}
			setColdHot(cold_hot); // will be used for shifted images and for texture tiles
		}
		/* Prepare 4-channel images*/
		ImagePlus [] imps_RGB = new ImagePlus[iclt_fimg.length];
		for (int ncam = 0; ncam < iclt_fimg.length; ncam++) if (iclt_fimg[ncam] != null){
            String title=String.format("%s%s-%02d",image_name, sAux(), ncam);
			imps_RGB[ncam] = linearStackToColor( // probably no need to separate and process the second half with quadCLT_aux (!)
					clt_parameters,
					colorProcParameters,
					rgbParameters,
					title, // String name,
					"-D"+clt_parameters.disparity, //String suffix, // such as disparity=...
					toRGB, // does not work here?
					!correctionsParameters.jpeg, // boolean bpp16, // 16-bit per channel color mode for result
					false, // true, // boolean saveShowIntermediate, // save/show if set globally
					false, // boolean saveShowFinal,        // save/show result (color image?)
					iclt_fimg[ncam],
					out_width,
					out_height,
					1.0, // scaleExposures[iAux][iSubCam], // double scaleExposure, // is it needed?
					-1); // debugLevel );
		}

		// combine to a sliced color image
		int [] slice_seq = {0,1,3,2}; //clockwise
		if (imps_RGB.length > 4) {
			slice_seq = new int [imps_RGB.length];
			for (int i = 0; i < slice_seq.length; i++) {
				slice_seq[i] = i;
			}
		}
		int width = imps_RGB[0].getWidth();
		int height = imps_RGB[0].getHeight();
		ImageStack array_stack=new ImageStack(width,height);
		for (int i = 0; i<slice_seq.length; i++) if (imps_RGB[slice_seq[i]] != null){
			array_stack.addSlice("port_"+slice_seq[i], imps_RGB[slice_seq[i]].getProcessor().getPixels());
		}
		ImagePlus imp_stack = new ImagePlus(image_name+sAux()+suffix, array_stack);
		imp_stack.getProcessor().resetMinAndMax();
		return imp_stack;
	}
	
	
	
	
	public ImagePlus processCLTQuadCorrGPU(
			ImagePlus []                                    imp_quad, //null will be OK
			boolean [][]                                    saturation_imp, // (near) saturated pixels or null // Not needed use this.saturation_imp
			CLTParameters                                   clt_parameters,
			EyesisCorrectionParameters.DebayerParameters    debayerParameters,
			ColorProcParameters                             colorProcParameters,
			CorrectionColorProc.ColorGainsParameters        channelGainParameters,
			EyesisCorrectionParameters.RGBParameters        rgbParameters,
			double []	                                    scaleExposures, // probably not needed here - restores brightness of the final image
			boolean                                         only4slice,
			final int                                       threadsMax,  // maximal number of threads to launch
			final boolean                                   updateStatus,
			final int                                       debugLevel){
		if (gpuQuad == null) {
			System.out.println("GPU instance is not initialized, using CPU mode");
			processCLTQuadCorrCPU(
					saturation_imp,        // boolean [][]                  saturation_imp, // (near) saturated pixels or null // Not needed use this.saturation_imp
					clt_parameters,        // CLTParameters                 clt_parameters,
					debayerParameters,     // EyesisCorrectionParameters.DebayerParameters     debayerParameters,
					colorProcParameters,   // ColorProcParameters                              colorProcParameters,
					channelGainParameters, // CorrectionColorProc.ColorGainsParameters         channelGainParameters,
					rgbParameters, // EyesisCorrectionParameters.RGBParameters             rgbParameters,
					scaleExposures, // double []	       scaleExposures, // probably not needed here
					false, // final boolean    apply_corr, // calculate and apply additional fine geometry correction
					false, // final boolean    infinity_corr, // calculate and apply geometry correction at infinity
					threadsMax,  // final int        threadsMax,  // maximal number of threads to launch
					updateStatus, // final boolean    updateStatus,
					debugLevel); // final int        debugLevel)			
			
			return null;
		}
		saveQuadClt(); // to re-load new set of Bayer images to the GPU (do nothing for CPU) and Geometry
// GPU-specific
		boolean showCoord = debugLevel > 1;
		boolean is_mono = isMonochrome();
		boolean is_lwir = isLwir();
		final boolean      batch_mode = clt_parameters.batch_run; //disable any debug images
		boolean test_execCorr2D =    batch_mode?false: false; // true;
		boolean try_lores =          batch_mode?false: true;
		boolean show_textures_rgba = batch_mode?false: clt_parameters.show_rgba_color;
		boolean try_textures =       batch_mode?false: true;
		
		ImageDtt image_dtt = new ImageDtt(
				getNumSensors(),
				clt_parameters.transform_size,
				clt_parameters.img_dtt,
				isAux(),
				is_mono,
				is_lwir,
				clt_parameters.getScaleStrength(isAux())); // 1.0);
		image_dtt.getCorrelation2d(); // initiate image_dtt.correlation2d, needed if disparity_map != null  
		float [][] lpf_rgb = new float[][] {
			image_dtt.floatGetCltLpfFd(clt_parameters.gpu_sigma_r),
			image_dtt.floatGetCltLpfFd(clt_parameters.gpu_sigma_b),
			image_dtt.floatGetCltLpfFd(clt_parameters.gpu_sigma_g),
			image_dtt.floatGetCltLpfFd(clt_parameters.gpu_sigma_m)
		};
		gpuQuad.setLpfRbg( // constants memory - same for all cameras
				lpf_rgb,
				true); // !batch_mode);

		float [] lpf_flat = image_dtt.floatGetCltLpfFd(clt_parameters.getGpuCorrSigma(is_mono));

		gpuQuad.setLpfCorr(// constants memory - same for all cameras
				"lpf_corr", // String const_name, // "lpf_corr"
				lpf_flat,
				!batch_mode);

		float [] lpf_rb_flat = image_dtt.floatGetCltLpfFd(clt_parameters.getGpuCorrRBSigma(is_mono));
		gpuQuad.setLpfCorr(// constants memory - same for all cameras
				"lpf_rb_corr", // String const_name, // "lpf_corr"
				lpf_rb_flat,
				!batch_mode);
		final boolean use_aux = false; // currently GPU is configured for a single quad camera
		
//		String image_name = image_name; // correctionsParameters.getModelName((String) imp_quad_main[0].getProperty("name"));
//		String image_path = image_path; //  (String) imp_quad_main[0].getProperty("path"); // Only for debug output
		// now set to only 4 !
		if (debugLevel>1) System.out.println("processing: "+image_path);

		
// end of GPU-specific		

		boolean advanced=  correctionsParameters.zcorrect || correctionsParameters.equirectangular;
		boolean toRGB=     advanced? true: correctionsParameters.toRGB;
		if (isLwir()) {
			toRGB = colorProcParameters.lwir_pseudocolor;
		}
		
		ImagePlus [] results = null;
		if (imp_quad !=  null) {
			results = new ImagePlus[imp_quad.length];
			for (int i = 0; i < results.length; i++) {
				results[i] = imp_quad[i];
				results[i].setTitle(results[i].getTitle()+"RAW");
			}
		}
		if (debugLevel>1) System.out.println("processing: "+gpuQuad.quadCLT);

		double z_correction =  clt_parameters.z_correction;
		if (clt_parameters.z_corr_map.containsKey(image_name)){
			z_correction +=clt_parameters.z_corr_map.get(image_name);// not used in lwir
		}
		final double disparity_corr = (z_correction == 0) ? 0.0 : geometryCorrection.getDisparityFromZ(1.0/z_correction);
		// the following will be replaced by setFullFrameImages() with target disparities
		TpTask [] tp_tasks  = gpuQuad.setFullFrameImages( // when disparities array is known - use different arguments of setFullFrameImages
				null,                                 // Rectangle                 woi,
				clt_parameters.gpu_woi_round,         // boolean                   round_woi,
	    		(float) clt_parameters.disparity,     // float                     target_disparity, // apply same disparity to all tiles
	    		(float) disparity_corr,
	    		0xf,                                  // int                       out_image, // from which tiles to generate image (currently 0/1)
	    		0x3f,                                 // int                       corr_mask,  // which correlation pairs to generate (maybe later - reduce size from 15x15)
	    		debugLevel);                          // final int                 debugLevel) - not yet used
		if (tp_tasks.length == 0) {
			System.out.println("Empty tasks - nothing to do");
			return null;
		}

		gpuQuad.setTasks( // copy tp_tasks to the GPU memory
				tp_tasks, // TpTask [] tile_tasks,
				use_aux, // boolean use_aux)
				clt_parameters.img_dtt.gpu_verify); // boolean verify

		gpuQuad.execSetTilesOffsets(true); // prepare tiles offsets in GPU memory // calculate tile centers

		if (showCoord) {
			double max_r_by_rdist_err = gpuQuad.maxRbyRDistErr();
			System.out.println("Maximal R-by-RDist error = "+max_r_by_rdist_err);
			showCentersPortsXY("Tile_Virtual_Centers");
//			showCentersXY("Tile_Virtual_Centers");
//			showPortsXY("Tile_Physical_Sensor_Centers");
		}
		
		
		gpuQuad.execConvertDirect(-1); // boolean erase_clt
		int mcorr_sel = Correlation2d.corrSelEncode(clt_parameters.img_dtt, getNumSensors());

		if (test_execCorr2D) {
			int [] i_mcorr_sel = Correlation2d.intCorrPairs(
					mcorr_sel,
					getNumSensors(),
					4); // int num_out); should be 4 int
			gpuQuad.setCorrMask(i_mcorr_sel);
			
			double fat_zero = 30.0; // 2000.0; // 30.00;
			int gpu_corr_rad =  7;
			int corr_size = 2* gpu_corr_rad + 1;

			final int numcol = isMonochrome()?1:3;
			final double [] col_weights= new double [numcol]; // colors are RBG
			if (isMonochrome()) {
				col_weights[0] = 1.00; // Was 0 ! 0;
			} else {
				col_weights[2] = 1.0/(1.0 + clt_parameters.corr_red + clt_parameters.corr_blue);    // green color
				col_weights[0] = clt_parameters.corr_red *  col_weights[2];
				col_weights[1] = clt_parameters.corr_blue * col_weights[2];
			}

			//Generate 2D phase correlations from the CLT representation
			// Try to zero out memory before calculating?	
			gpuQuad.eraseGpuCorrs(); // no-diff
			gpuQuad.execCorr2D(
					gpuQuad.getCorrMask(),    // boolean [] pair_select,
					col_weights,              // double [] scales,
					fat_zero,                 // double fat_zero);
					gpu_corr_rad);            // int corr_radius
			// Should be done before execCorr2D_TD as corr_indices are shared to save memory
			int [] corr_indices = gpuQuad.getCorrIndices();
			// the following is not yet shared
			float [][] corr2D =  gpuQuad.getCorr2D(
					gpu_corr_rad); //  int corr_rad);
			final int tilesX=gpuQuad.getTilesX(); // width/transform_size;
			final int tilesY=gpuQuad.getTilesY(); // final int tilesY=height/transform_size;
//			int num_tiles = tilesX * tilesY;
			int sq = 16;
			int num_pairs = gpuQuad.getNumUsedPairs();
			float [][] corr_img = new float [num_pairs][tilesY * sq * tilesX * sq];
			for (int pair = 0; pair < num_pairs; pair++) {
				Arrays.fill(corr_img[pair], Float.NaN);
			}
			for (int ict = 0; ict < corr_indices.length; ict++){
				//    	int ct = cpu_corr_indices[ict];
				int ctt = ( corr_indices[ict] >>  GPUTileProcessor.CORR_NTILE_SHIFT);
				int cpair = corr_indices[ict] & ((1 << GPUTileProcessor.CORR_NTILE_SHIFT) - 1);
				int ty = ctt / tilesX;
				int tx = ctt % tilesX;
				int dst_offs0 = (ty * sq * tilesX * sq) + (tx * sq);
				for (int iy = 0; iy < corr_size; iy++){
					int dst_offs = dst_offs0 + iy * (tilesX * sq);
					System.arraycopy(corr2D[ict], iy * corr_size, corr_img[cpair], dst_offs, corr_size);
				}
			}
			ShowDoubleFloatArrays.showArrays( // out of boundary 15
					corr_img,
					tilesX * sq,
					tilesY * sq,
					true,
					"test-corr1");
		}
		
		
		gpuQuad.execImcltRbgAll(is_mono); 

		// get data back from GPU
//		float [][][] iclt_fimg = new float [gpuQuad.getNumCams()][][];
		float [][][] iclt_fimg = new float [getNumSensors()][][];

		for (int ncam = 0; ncam < iclt_fimg.length; ncam++) {
			iclt_fimg[ncam] = gpuQuad.getRBG(ncam);
		}

		int out_width =  gpuQuad.getImageWidth();//   + gpuQuad.getDttSize(); // 2022/05/12 removed margins from gpuQuad.getRBG(ncam);
		int out_height = gpuQuad.getImageHeight(); // + gpuQuad.getDttSize(); // 2022/05/12 removed margins from gpuQuad.getRBG(ncam);
		if (isLwir() && colorProcParameters.lwir_autorange) {
			double rel_low =  colorProcParameters.lwir_low;
			double rel_high = colorProcParameters.lwir_high;
			if (!Double.isNaN(getLwirOffset())) {
				rel_low -=  getLwirOffset();
				rel_high -= getLwirOffset();
			}
			double [] cold_hot =  autorange(
					iclt_fimg,                         // iclt_data, // double [][][] iclt_data, //  [iQuad][ncol][i] - normally only [][2][] is non-null
					rel_low,                           // double hard_cold,// matches data, DC (this.lwir_offset)  subtracted
					rel_high,                          // double hard_hot,   // matches data, DC (this.lwir_offset)  subtracted
					colorProcParameters.lwir_too_cold, // double too_cold, // pixels per image
					colorProcParameters.lwir_too_hot,  // double too_hot,  // pixels per image
					1024); // int num_bins)
			if (cold_hot != null) {
				if (!Double.isNaN(getLwirOffset())) {
					cold_hot[0] += getLwirOffset();
					cold_hot[1] += getLwirOffset();
				}
			}
			setColdHot(cold_hot); // will be used for shifted images and for texture tiles
		}
		
		/* Prepare 4-channel images*/
		ImagePlus [] imps_RGB = new ImagePlus[iclt_fimg.length];
		for (int ncam = 0; ncam < iclt_fimg.length; ncam++) {
            String title=String.format("%s%s-%02d",image_name, sAux(), ncam);
			imps_RGB[ncam] = linearStackToColor( // probably no need to separate and process the second half with quadCLT_aux (!)
					clt_parameters,
					colorProcParameters,
					rgbParameters,
					title, // String name,
					"-D"+clt_parameters.disparity, //String suffix, // such as disparity=...
					toRGB,
					!correctionsParameters.jpeg, // boolean bpp16, // 16-bit per channel color mode for result
					!batch_mode, // true, // boolean saveShowIntermediate, // save/show if set globally
					false, // boolean saveShowFinal,        // save/show result (color image?)
					iclt_fimg[ncam],
					out_width,
					out_height,
					1.0, // scaleExposures[iAux][iSubCam], // double scaleExposure, // is it needed?
					debugLevel );
		}
		
		if (clt_parameters.gen_chn_img || only4slice) { // save and show 4-slice image
			// combine to a sliced color image
			int [] slice_seq = {0,1,3,2}; //clockwise
			if (imps_RGB.length > 4) {
				slice_seq = new int [imps_RGB.length];
				for (int i = 0; i < slice_seq.length; i++) {
					slice_seq[i] = i;
				}
			}
			int width = imps_RGB[0].getWidth();
			int height = imps_RGB[0].getHeight();
			ImageStack array_stack=new ImageStack(width,height);
			for (int i = 0; i<slice_seq.length; i++){
				///				if (imps_RGB[slice_seq[i]] != null) {
				array_stack.addSlice("port_"+slice_seq[i], imps_RGB[slice_seq[i]].getProcessor().getPixels());
				///				} else {
				///					array_stack.addSlice("port_"+slice_seq[i], results[slice_seq[i]].getProcessor().getPixels());
				///				}
			}
			ImagePlus imp_stack = new ImagePlus(image_name+sAux()+"GPU-SHIFTED-D"+clt_parameters.disparity, array_stack);
			imp_stack.getProcessor().resetMinAndMax();
			if (only4slice) {
				return imp_stack;
			} else {
				if (!batch_mode) {
					imp_stack.updateAndDraw();
				}
				eyesisCorrections.saveAndShowEnable(
						imp_stack,  // ImagePlus             imp,
						correctionsParameters, // EyesisCorrectionParameters.CorrectionParameters  correctionsParameters,
						true, // boolean               enableSave,
						!batch_mode) ;// boolean               enableShow);
			}
		}

		if (clt_parameters.gen_4_img) { // save 4 JPEG images
			// Save as individual JPEG images in the model directory
			String x3d_path = getX3dDirectory();
			for (int sub_img = 0; sub_img < imps_RGB.length; sub_img++){
				EyesisCorrections.saveAndShow(
						imps_RGB[sub_img],
						x3d_path,
						correctionsParameters.png && !clt_parameters.black_back,
						!batch_mode && clt_parameters.show_textures,
						correctionsParameters.JPEG_quality, // jpegQuality); // jpegQuality){//  <0 - keep current, 0 - force Tiff, >0 use for JPEG
						(debugLevel > 0) ? debugLevel : 1); // int debugLevel (print what it saves)
			}
			String model_path= correctionsParameters.selectX3dDirectory(
					image_name, // quad timestamp. Will be ignored if correctionsParameters.use_x3d_subdirs is false
					null,
					true,  // smart,
					true);  //newAllowed, // save

			createThumbNailImage(
					imps_RGB[0],
					model_path,
					"thumb"+sAux(),
					debugLevel);

		}
// try textures here
		if (show_textures_rgba) {
			float [][] texture_img = new float [isMonochrome()?2:4][];
//			double [] col_weights = new double[3];
			double [] col_weights = new double[isMonochrome() ? 1 :3];
			if (isMonochrome()) {
				col_weights[0] = 1.0;
//			if (isMonochrome()) {
//				col_weights[0] = 1.0;
//				col_weights[1] = 0.0;
//				col_weights[2] = 0.0;// green color/mono
			} else {
				col_weights[2] = 1.0/(1.0 +  clt_parameters.corr_red + clt_parameters.corr_blue);    // green color
				col_weights[0] = clt_parameters.corr_red *  col_weights[2];
				col_weights[1] = clt_parameters.corr_blue * col_weights[2];
			}
			Rectangle woi = new Rectangle(); // will be filled out to match actual available image
//			int       keep_weights = (clt_parameters.keep_weights? 1 : 0) + (clt_parameters.replace_weights? 2 : 0);
			int  keep_weights = (clt_parameters.replace_weights? 2 : 0); // just 0/2
			int text_colors = (isMonochrome() ? 1 : 3);
			int  texture_layers = (text_colors + 1)+((keep_weights != 0)?(text_colors * getNumSensors()):0);   
			gpuQuad.execRBGA(
					col_weights,                   // double [] color_weights,
					isLwir(),                      // boolean   is_lwir,
					clt_parameters.min_shot,       // double    min_shot,           // 10.0
					clt_parameters.scale_shot,     // double    scale_shot,         // 3.0
					clt_parameters.diff_sigma,     // double    diff_sigma,         // pixel value/pixel change Used much larger sigma = 10.0 instead of 1.5
					clt_parameters.diff_threshold, // double    diff_threshold,     // pixel value/pixel change
					clt_parameters.min_agree,      // double    min_agree,          // minimal number of channels to agree on a point (real number to work with fuzzy averages)
					clt_parameters.dust_remove,    // boolean   dust_remove,
					keep_weights);                            // int       keep_weights)
			float [][] rbga = gpuQuad.getRBGA(
					texture_layers - 1, // (isMonochrome() ? 1 : 3), // int     num_colors,
					woi);
			for (int ncol = 0; ncol < texture_img.length; ncol++) if (ncol < rbga.length) {
				texture_img[ncol] = rbga[ncol];
			}
			// first try just multi-layer, then with palette
			
			ShowDoubleFloatArrays.showArrays( // show slices RBGA (colors - 256, A - 1.0)
					rbga,
					woi.width,
					woi.height,
					true,
					getImageName()+"-RGBA-STACK-D"+clt_parameters.disparity+
					":"+clt_parameters.gpu_woi_tx+":"+clt_parameters.gpu_woi_ty+
					":"+clt_parameters.gpu_woi_twidth+":"+clt_parameters.gpu_woi_theight+
					":"+(clt_parameters.gpu_woi_round?"C":"R")
					//,new String[] {"R","B","G","A"}
					);
			float [][] rgb_text; //  = new float [(isMonochrome() ? 1 : 3)][]; // {rbga[0],rbga[1],rbga[2]};
			float [][] rgba_text; //  = new float [(isMonochrome() ? 2 : 4)][]; // {rbga[0],rbga[1],rbga[2],rbga[3]};
			if (isMonochrome()) {
				rgb_text = new float [][]{rbga[0]};
				rgba_text = new float [][]{rbga[0]};
			} else {
				rgb_text = new float [][]{rbga[0],rbga[2],rbga[1]};
				rgba_text = new float [][]{rbga[0],rbga[2],rbga[1],rbga[3]};
			}
			//isAux()
			ImagePlus imp_rgba = linearStackToColor(
					clt_parameters,
					colorProcParameters,
					rgbParameters,
					getImageName()+"-texture", // String name,
					"-D"+clt_parameters.disparity+"-"+(isAux()?"AUX":"MAIN")+ "GPU", //String suffix, // such as disparity=...
					toRGB,
					!correctionsParameters.jpeg, // boolean bpp16, // 16-bit per channel color mode for result
					false, // true, // boolean saveShowIntermediate, // save/show if set globally
					false, // true, // boolean saveShowFinal,        // save/show result (color image?)
					((clt_parameters.alpha1 > 0)? rgba_text: rgb_text),
					woi.width, // clt_parameters.gpu_woi_twidth * image_dtt.transform_size, //  tilesX *  image_dtt.transform_size,
					woi.height, // clt_parameters.gpu_woi_theight *image_dtt.transform_size, //  tilesY *  image_dtt.transform_size,
					1.0,         // double scaleExposure, // is it needed?
					debugLevel );
			imp_rgba.getProcessor().resetMinAndMax();
			imp_rgba.show();
		}
		if (try_lores) {
			//Generate non-overlapping (16x16) texture tiles, prepare
//			double [] col_weights = new double[3];
			double [] col_weights = new double[isMonochrome() ? 1 :3];
			if (isMonochrome()) {
				col_weights[0] = 1.0;
			} else {
				col_weights[2] = 1.0/(1.0 +  clt_parameters.corr_red + clt_parameters.corr_blue);    // green color
				col_weights[0] = clt_parameters.corr_red *  col_weights[2];
				col_weights[1] = clt_parameters.corr_blue * col_weights[2];
			}
			gpuQuad.execTextures(
				col_weights,                   // double [] color_weights,
				isLwir(),                      // boolean   is_lwir,
				clt_parameters.min_shot,       // double    min_shot,           // 10.0
				clt_parameters.scale_shot,     // double    scale_shot,         // 3.0
				clt_parameters.diff_sigma,     // double    diff_sigma,         // pixel value/pixel change Used much larger sigma = 10.0 instead of 1.5
				clt_parameters.diff_threshold, // double    diff_threshold,     // pixel value/pixel change
				clt_parameters.min_agree,      // double    min_agree,          // minimal number of channels to agree on a point (real number to work with fuzzy averages)
				clt_parameters.dust_remove,    // boolean   dust_remove,
				0,                             // int       keep_weights,       // 2 bits now, move to parameters
				false,                         // boolean   calc_textures,
				true,                          // boolean   calc_extra
				false);                        // boolean   linescan_order) // TODO: use true to avoid reordering of the low-res output 
			float [][] extra = gpuQuad.getExtra(); // now 4*numSensors
			int num_cams = getNumSensors();
			
			ShowDoubleFloatArrays.showArrays( // show slices RBGA (colors - 256, A - 1.0)
					extra,
					gpuQuad.img_width / GPUTileProcessor.DTT_SIZE,
					gpuQuad.img_height / GPUTileProcessor.DTT_SIZE,
					true,
					getImageName()+"-LOW-RES"
					//,new String[] {"R","B","G","A"}
					);
			for (int nc = 00; nc < (extra.length - num_cams); nc++) {
				int sindx = nc + num_cams;
				
			}
			if (try_textures) {
				int       keep_weights = (clt_parameters.keep_weights? 1 : 0) + (clt_parameters.replace_weights? 2 : 0);
				//Generate non-overlapping (16x16) texture tiles, prepare 
				gpuQuad.execTextures(
					col_weights,                   // double [] color_weights,
					isLwir(),         // boolean   is_lwir,
					clt_parameters.min_shot,       // double    min_shot,           // 10.0
					clt_parameters.scale_shot,     // double    scale_shot,         // 3.0
					clt_parameters.diff_sigma,     // double    diff_sigma,         // pixel value/pixel change Used much larger sigma = 10.0 instead of 1.5
					clt_parameters.diff_threshold, // double    diff_threshold,     // pixel value/pixel change
					clt_parameters.min_agree,      // double    min_agree,          // minimal number of channels to agree on a point (real number to work with fuzzy averages)
					clt_parameters.dust_remove,    // boolean   dust_remove,
					keep_weights,                  // int       keep_weights,       // 2 bits now, move to parameters
					true,                          // boolean   calc_textures,
					false,                         // boolean   calc_extra
					false);                        // boolean   linescan_order) 

				int [] texture_indices = gpuQuad.getTextureIndices();
				int numcol = isMonochrome()? 1 : 3;
//				int    num_src_slices = numcol + 1; //  + (clt_parameters.keep_weights?(ports + numcol + 1):0); // 12 ; // calculate
				int    num_src_slices = numcol + 1  + (clt_parameters.keep_weights?(getNumSensors() + numcol + 1):0); // 12 ; // calculate
				float [] flat_textures =  gpuQuad.getFlatTextures( // fatal error has been detected by the Java Runtime Environment:
						texture_indices.length,
						numcol, // int     num_colors,
			    		(keep_weights != 0)); // false); // clt_parameters.keep_weights); // boolean keep_weights);
				if (keep_weights != 0) {
					
				}
				// numcol + 1  + (clt_parameters.keep_weights?(getNumSensors() + numcol + 1):0); // 12 ; // calculate
				int tilesX = gpuQuad.img_width / GPUTileProcessor.DTT_SIZE;
				int tilesY = gpuQuad.img_height / GPUTileProcessor.DTT_SIZE;
				double [][][][] texture_tiles = new double [tilesY][tilesX][][];
				GpuQuad.doubleTextures(
			    		new Rectangle(0, 0, tilesX, tilesY), // Rectangle    woi,
			    		texture_tiles,                       // double [][][][] texture_tiles, // null or [tilesY][tilesX]
			    		texture_indices,                     // int []       indices,
			    		flat_textures,                       // float [][][] ftextures,
			    		tilesX,                              // int          full_width,
			    		num_src_slices, // isMonochrome()? 2: 4,                // rbga only /int          num_slices
			    		num_src_slices                       // int          num_src_slices
			    		);
				int num_out_slices = 0;
				for (int nt = 0; nt < tilesY * tilesX; nt++) {
					if (texture_tiles[nt / tilesX][nt % tilesX] != null) {
						num_out_slices = texture_tiles[nt / tilesX][nt % tilesX].length;
						break;
					}
				}
				if (num_out_slices > 0) {
					String [] dbg_titles = new String[num_out_slices];
					fill_titles: 
					{
						int indx = 0;
						for (int i = 0; i < numcol; i++) {
							if (indx >= num_out_slices)	break fill_titles;
							else dbg_titles[indx++] = "c"+i;
						}
						if (indx >= num_out_slices)	break fill_titles;
						else dbg_titles[indx++] = "alpha";
						for (int i = 0; i < getNumSensors(); i++) {
							if (indx >= num_out_slices)	break fill_titles;
							else dbg_titles[indx++] = (((keep_weights & 2) != 0)?"t":"w")+i;
						}
						for (int i = 0; i < numcol; i++) {
							if (indx >= num_out_slices)	break fill_titles;
							else dbg_titles[indx++] = "crms"+i;
						}
						if (indx >= num_out_slices)	break fill_titles;
						else dbg_titles[indx++] = "crms";
					}
							//crms
					
					int ssize = 2*GPUTileProcessor.DTT_SIZE;
					int width =  tilesX * ssize;
					int height = tilesY * ssize;
					double [][] dbg_nonoverlap = new double[num_out_slices][width * height];
					for (int slice = 0; slice < num_out_slices; slice++) {
						Arrays.fill(dbg_nonoverlap[slice], Double.NaN);
					}
					for (int ty = 0; ty < tilesY; ty++) {
						for (int tx = 0; tx < tilesX; tx++) {
							if (texture_tiles[ty][tx] != null) {
								for (int slice = 0; slice < num_out_slices; slice++) {
									for (int row = 0; row < ssize; row++) {
										System.arraycopy(
												texture_tiles[ty][tx][slice],
												row * ssize,
												dbg_nonoverlap[slice],
												(ty * ssize + row) * width + (tx * ssize),
												ssize);
									}
								}
							}
						}
						
					}
					ShowDoubleFloatArrays.showArrays( // show slices RBGA (colors - 256, A - 1.0)
							dbg_nonoverlap,
							width,
							height,
							true,
							getImageName()+"-textures-"+keep_weights,
							dbg_titles);
				}
				System.out.println("try_textures DONE");
			}
		}
// try low-res and non-overlap textures
		return null;	  
	}
	
	
	
	
	
	public static ImagePlus [] processCLTQuadCorrPairGpu(
			QuadCLT                                         quadCLT_main, // use gpuQuad_main.quadCLT
			QuadCLT                                         quadCLT_aux, // use gpuQuad_aux.quadCLT
			ImagePlus []                                    imp_quad_main,
			ImagePlus []                                    imp_quad_aux,
			boolean [][]                                    saturation_main, // (near) saturated pixels or null
			boolean [][]                                    saturation_aux, // (near) saturated pixels or null
			CLTParameters                                   clt_parameters,
			EyesisCorrectionParameters.CorrectionParameters ecp,
			EyesisCorrectionParameters.DebayerParameters    debayerParameters,
			ColorProcParameters                             colorProcParameters,
			ColorProcParameters                             colorProcParameters_aux,
			EyesisCorrectionParameters.RGBParameters        rgbParameters,
			double []	                                    scaleExposures_main, // probably not needed here - restores brightness of the final image
			double []	                                    scaleExposures_aux, // probably not needed here - restores brightness of the final image
			boolean                                         notch_mode, // not used here: use pole-detection mode for inter-camera correlation
			final int                                       lt_rad,     // not used here: low texture mode - inter-correlation is averaged between the neighbors before argmax-ing, using
			final int        threadsMax,  // maximal number of threads to launch
			final boolean    updateStatus,
			final int        debugLevel){
		boolean resetGC = false;
		boolean resetEV = false;
		
		if (resetGC) {
			if (quadCLT_main.getGPU() != null) quadCLT_main.getGPU().resetGeometryCorrection();
			if (quadCLT_aux.getGPU() != null)  quadCLT_aux.getGPU().resetGeometryCorrection();
		}

		if (resetEV) {
			if (quadCLT_main.getGPU() != null) quadCLT_main.gpuResetCorrVector(); // getGPU().resetGeometryCorrectionVector();
			if (quadCLT_aux.getGPU() != null)  quadCLT_aux.gpuResetCorrVector();   // .getGPU().resetGeometryCorrectionVector();
		}
		
// get fat_zero (absolute) and color scales
		boolean is_mono = quadCLT_main.isMonochrome();
		boolean is_lwir = quadCLT_main.isLwir();
		boolean is_aux =  quadCLT_main.isAux();


		double    fat_zero = clt_parameters.getGpuFatZero(is_mono); //   30.0;
/*		
		double [] scales = (is_mono) ? (new double [] {1.0}) :(new double [] {
				clt_parameters.gpu_weight_r, // 0.25
				clt_parameters.gpu_weight_b, // 0.25
				1.0 - clt_parameters.gpu_weight_r - clt_parameters.gpu_weight_b}); // 0.5
*/
		double corr_green = 1.0/(1.0 + clt_parameters.gpu_weight_r + clt_parameters.gpu_weight_b);    // green color
		// the same, as col_weights, for separate debugging
		double [] scales= (is_mono) ? (new double [] {1.0}) :(new double [] {
				clt_parameters.gpu_weight_r *  corr_green,
				clt_parameters.gpu_weight_b * corr_green,
				corr_green});
		
		
		double cwgreen = 1.0/(1.0 + clt_parameters.corr_red + clt_parameters.corr_blue);    // green color
		double [] col_weights= (is_mono) ? (new double [] {1.0}) :(new double [] {
				clt_parameters.corr_red *  cwgreen,
				clt_parameters.corr_blue * cwgreen,
				cwgreen});

		ImageDtt image_dtt = new ImageDtt(
				quadCLT_main.getNumSensors(),
				clt_parameters.transform_size,
				clt_parameters.img_dtt,
				is_aux,
				is_mono,
				is_lwir,
				1.0);
		image_dtt.getCorrelation2d(); // initiate image_dtt.correlation2d, needed if disparity_map != null  
		float [][] lpf_rgb = new float[][] {
			image_dtt.floatGetCltLpfFd(clt_parameters.gpu_sigma_r),
			image_dtt.floatGetCltLpfFd(clt_parameters.gpu_sigma_b),
			image_dtt.floatGetCltLpfFd(clt_parameters.gpu_sigma_g),
			image_dtt.floatGetCltLpfFd(clt_parameters.gpu_sigma_m)
		};
		quadCLT_main.getGPU().setLpfRbg( // constants memory - same for all cameras
				lpf_rgb,
				debugLevel > 2); // -3

		float [] lpf_flat = image_dtt.floatGetCltLpfFd(clt_parameters.getGpuCorrSigma(is_mono));

		quadCLT_main.getGPU().setLpfCorr(// constants memory - same for all cameras
				"lpf_corr", // String const_name, // "lpf_corr"
				lpf_flat,
				debugLevel > 2); // -3

		float [] lpf_rb_flat = image_dtt.floatGetCltLpfFd(clt_parameters.getGpuCorrRBSigma(is_mono));
		quadCLT_main.getGPU().setLpfCorr(// constants memory - same for all cameras
				"lpf_rb_corr", // String const_name, // "lpf_corr"
				lpf_rb_flat,
				debugLevel > 2); // -3

		final boolean use_aux = false; // currently GPU is configured for a single quad camera

		final boolean      batch_mode = clt_parameters.batch_run; //disable any debug images
		boolean toRGB=     quadCLT_main.correctionsParameters.toRGB;

		// may use this.StartTime to report intermediate steps execution times
		String name = quadCLT_main.image_name; // correctionsParameters.getModelName((String) imp_quad_main[0].getProperty("name"));
		String path = quadCLT_main.image_path; //  (String) imp_quad_main[0].getProperty("path"); // Only for debug output
		// now set to only 4 !
		ImagePlus [] results = new ImagePlus[imp_quad_main.length]; // + imp_quad_aux.length];
		for (int i = 0; i < results.length; i++) {
			if (i< imp_quad_main.length) {
				results[i] = imp_quad_main[i];
			} else {
				results[i] = imp_quad_aux[i-imp_quad_main.length];
			}
			results[i].setTitle(results[i].getTitle()+"RAW");
		}
		if (debugLevel>1) System.out.println("processing: "+path);

		// Set task clt_parameters.disparity
		Rectangle twoi = clt_parameters.gpu_woi? (new Rectangle ( // normally - full window (0,0,324,242)
				clt_parameters.gpu_woi_tx,
				clt_parameters.gpu_woi_ty,
				clt_parameters.gpu_woi_twidth,
				clt_parameters.gpu_woi_theight)): null;
		
		// the following will be replaced by setFullFrameImages() with target disparities
		
		double z_correction =  clt_parameters.z_correction;
		if (clt_parameters.z_corr_map.containsKey(quadCLT_main.image_name)){
			z_correction +=clt_parameters.z_corr_map.get(quadCLT_main.image_name);// not used in lwir
		}
		final double disparity_corr = (z_correction == 0) ? 0.0 : quadCLT_main.geometryCorrection.getDisparityFromZ(1.0/z_correction);
		
		
		TpTask [] tp_tasks  = quadCLT_main.getGPU().setFullFrameImages( // when disparities array is known - use different arguments of setFullFrameImages
				twoi,                                 // Rectangle                 woi,
				clt_parameters.gpu_woi_round,         // boolean                   round_woi,
	    		(float) clt_parameters.disparity,     // float                     target_disparity, // apply same disparity to all tiles
	    		(float) disparity_corr,               // float                     disparity_corr, // add to disparity (at infinity)
	    		0xf,                                  // int                       out_image, // from which tiles to generate image (currently 0/1)
	    		0x3f,                                 // int                       corr_mask,  // which correlation pairs to generate (maybe later - reduce size from 15x15)
	    		debugLevel);                          // final int                 debugLevel) - not yet used
		if (tp_tasks.length == 0) {
			System.out.println("--- Empty tasks - nothing to do ---");
			return null;
		}

		quadCLT_main.getGPU().setTasks( // copy tp_tasks to the GPU memory
				tp_tasks, // TpTask [] tile_tasks,
				use_aux, // boolean use_aux)
				clt_parameters.img_dtt.gpu_verify); // boolean verify
				
/*
		quadCLT_main.getGPU().setGeometryCorrection( // copy Geometry correction data to the GPU memory (once per camera group?), allocate GPU memory for the rotation/deriv. matrices
				quadCLT_main.getGeometryCorrection(),
				false); // boolean use_java_rByRDist) { // false - use newer GPU execCalcReverseDistortions); // once
		
		quadCLT_main.getGPU().setExtrinsicsVector(quadCLT_main.getGeometryCorrection().getCorrVector()); // for each new image - copy geometry correction vector to the GPU
		
		// Optionally save offsets here?
		if (clt_parameters.gpu_save_ports_xy) {
			savePortsXY(
					ecp,        // EyesisCorrectionParameters.CorrectionParameters ecp,
					tp_tasks);  // GPUTileProcessor.TpTask [] tp_tasks
		}
*/		
		// All set, run kernel (correct and convert)
		int NREPEAT = 1; // 00;
		System.out.println("\n------------ Running GPU "+NREPEAT+" times ----------------");
		long startGPU=System.nanoTime();
		for (int i = 0; i < NREPEAT; i++ ) {
			// Calculate reverse radial distortion table from gpu_geometry_correction
			// Needed once during initialization of GeometryCorrection
/// Following will be executed when/if needed			
///			quadCLT_main.getGPU().execCalcReverseDistortions();       

		}
		long startRotDerivs=System.nanoTime();
		for (int i = 0; i < NREPEAT; i++ ) {
			// Calculate rotation matrices and their derivatives
			// Needed after extrinsics changed
/// Following will be executed when/if needed			
///			quadCLT_main.getGPU().execRotDerivs();
		}

		long startTasksSetup=System.nanoTime();
		for (int i = 0; i < NREPEAT; i++ ) {
			// Calculate tiles offsets (before each direct conversion run)
			quadCLT_main.getGPU().execSetTilesOffsets(true); // prepare tiles offsets in GPU memory // calculate tile centers
		}

		long startDirectConvert=System.nanoTime();

		for (int i = 0; i < NREPEAT; i++ ) {
			// Direct CLT conversion and aberration correction
			quadCLT_main.getGPU().execConvertDirect(-1); // boolean erase_clt
		}

		long startIMCLT=System.nanoTime();
		for (int i = 0; i < NREPEAT; i++ ) {
			// Generate corrected image(s) from the CLT representation created by execConvertDirect()
			quadCLT_main.getGPU().execImcltRbgAll(quadCLT_main.isMonochrome());
		}
		long endImcltTime = System.nanoTime();

		// FIXME: Temporary setting all pairs to correlate
		int num_pairs = quadCLT_main.getGPU().getNumPairs();
		boolean [] sel_pairs = new boolean [num_pairs];
		Arrays.fill(sel_pairs, true);
		
		long startCorr2d=System.nanoTime();   // System.nanoTime();
		
		for (int i = 0; i < NREPEAT; i++ ) {
			//Generate 2D phase correlations from the CLT representation
			quadCLT_main.getGPU().execCorr2D(
					sel_pairs,                    // boolean [] pair_select,
					scales,                       // double [] scales,
					fat_zero,                     // double fat_zero);
					clt_parameters.gpu_corr_rad); // int corr_radius
		}
		long endCorr2d = System.nanoTime();
		
		// Should be done before execCorr2D_TD as corr_indices are shared to save memory
		int [] corr_indices = quadCLT_main.getGPU().getCorrIndices();
		// the following is not yet shared
		float [][] corr2D = quadCLT_main.getGPU().getCorr2D(
				clt_parameters.gpu_corr_rad); //  int corr_rad);
		
// calculate correlations, keep TD
		int mcorr_sel = Correlation2d.corrSelEncode(clt_parameters.img_dtt,quadCLT_main.getNumSensors());
		int [] i_mcorr_sel = Correlation2d.intCorrPairs(
				mcorr_sel,
				quadCLT_main.getNumSensors(),
				4); // int num_out); should be 4 int

		quadCLT_main.getGPU().execCorr2D_TD(
	    		scales,i_mcorr_sel);// double [] scales,
		
		quadCLT_main.getGPU().execCorr2D_combine( // calculate cross pairs
		        true, // boolean init_corr,    // initialize output tiles (false - add to current)
		        6,    // int     num_pairs_in, // typically 6 - number of pairs per tile (tile task should have same number per each tile
		        0x0f); // int     pairs_mask    // selected pairs (0x3 - horizontal, 0xc - vertical, 0xf - quad, 0x30 - cross)
		
		quadCLT_main.getGPU().execCorr2D_normalize(
        		true, // boolean combo, // normalize combo correlations (false - per-pair ones) 
	    		fat_zero, // double fat_zero);
				null, // float [] fcorr_weights, // null or one per correlation tile (num_corr_tiles) to divide fat zero2
	    		clt_parameters.gpu_corr_rad); // int corr_radius
// run textures
		long startTextures = System.nanoTime();   // System.nanoTime();
		boolean   calc_textures = clt_parameters.gpu_show_jtextures; //  true;
		boolean   calc_extra =    clt_parameters.gpu_show_extra; //  true;
		for (int i = 0; i < NREPEAT; i++ )
			//Generate non-overlapping (16x16) texture tiles, prepare 
			quadCLT_main.getGPU().execTextures(
				col_weights,                   // double [] color_weights,
				quadCLT_main.isLwir(),         // boolean   is_lwir,
				clt_parameters.min_shot,       // double    min_shot,           // 10.0
				clt_parameters.scale_shot,     // double    scale_shot,         // 3.0
				clt_parameters.diff_sigma,     // double    diff_sigma,         // pixel value/pixel change
				clt_parameters.diff_threshold, // double    diff_threshold,     // pixel value/pixel change - never used in GPU ?
				clt_parameters.min_agree,      // double    min_agree,          // minimal number of channels to agree on a point (real number to work with fuzzy averages)
				clt_parameters.dust_remove,    // boolean   dust_remove,        // Do not reduce average weight when only one image differs much from the average
				0,                             // int       keep_weights,       // 2 bits now, move to parameters
				calc_textures,                 // boolean   calc_textures,
				calc_extra,                    // boolean   calc_extra)
				false);                        // boolean   linescan_order) // TODO: use true to avoid reordering of the low-res output 
				

		long endTextures = System.nanoTime();
// run texturesRBGA
		long startTexturesRBGA = System.nanoTime();   // System.nanoTime();

		for (int i = 0; i < NREPEAT; i++ )
// Generate combined (overlapping) texture
			quadCLT_main.getGPU().execRBGA(
				col_weights,                   // double [] color_weights,
				quadCLT_main.isLwir(),         // boolean   is_lwir,
				clt_parameters.min_shot,       // double    min_shot,           // 10.0
				clt_parameters.scale_shot,     // double    scale_shot,         // 3.0
				clt_parameters.diff_sigma,     // double    diff_sigma,         // pixel value/pixel change
				clt_parameters.diff_threshold, // double    diff_threshold,     // pixel value/pixel change
				clt_parameters.min_agree,      // double    min_agree,          // minimal number of channels to agree on a point (real number to work with fuzzy averages)
				clt_parameters.dust_remove,   // boolean   dust_remove,
				0);             // int       keep_weights)

		long endTexturesRBGA = System.nanoTime();
		long endGPUTime = System.nanoTime();

		long calcReverseTime=      (startRotDerivs-     startGPU)           /NREPEAT;
		long rotDerivsTime=        (startTasksSetup-    startRotDerivs)     /NREPEAT;
		long tasksSetupTime=       (startDirectConvert- startTasksSetup)    /NREPEAT;
		long firstGPUTime=         (startIMCLT-         startDirectConvert) /NREPEAT;
		long runImcltTime =        (endImcltTime -      startIMCLT)         /NREPEAT;
		long runCorr2DTime =       (endCorr2d -         startCorr2d)        /NREPEAT;
		long runTexturesTime =     (endTextures -       startTextures)      /NREPEAT;
		long runTexturesRBGATime = (endTexturesRBGA -   startTexturesRBGA)  /NREPEAT;
		long runGPUTime =          (endGPUTime -        startGPU)           /NREPEAT;
		// run corr2d
//RotDerivs
		System.out.println("\n------------ End of running GPU "+NREPEAT+" times ----------------");
		System.out.println("GPU run time ="+        (runGPUTime * 1.0e-6)+"ms");
		System.out.println(" - calc reverse dist.: "+(calcReverseTime*1.0e-6)+"ms");
		System.out.println(" - rot/derivs:         "+(rotDerivsTime*1.0e-6)+"ms");
		System.out.println(" - tasks setup:        "+(tasksSetupTime*1.0e-6)+"ms");
		System.out.println(" - direct conversion:  "+(firstGPUTime*1.0e-6)+"ms");
		System.out.println(" - imclt:              "+(runImcltTime*1.0e-6)+"ms");
		System.out.println(" - corr2D:             "+(runCorr2DTime*1.0e-6)+"ms");
		System.out.println(" - textures:           "+(runTexturesTime*1.0e-6)+"ms");
		System.out.println(" - RGBA:               "+(runTexturesRBGATime*1.0e-6)+"ms");
		// get data back from GPU
//		float [][][] iclt_fimg = new float [quadCLT_main.getGPU().getNumCams()][][];
		float [][][] iclt_fimg = new float [quadCLT_main.getNumSensors()][][];
		
		for (int ncam = 0; ncam < iclt_fimg.length; ncam++) {
			iclt_fimg[ncam] = quadCLT_main.getGPU().getRBG(ncam);
		}

		int out_width =  quadCLT_main.getGPU().getImageWidth(); //  + quadCLT_main.getGPU().getDttSize(); // 2022/05/12 removed margins from gpuQuad.getRBG(ncam);
		int out_height = quadCLT_main.getGPU().getImageHeight(); // + quadCLT_main.getGPU().getDttSize();// 2022/05/12 removed margins from gpuQuad.getRBG(ncam);
		int tilesX =     quadCLT_main.getGPU().getImageWidth()  / quadCLT_main.getGPU().getDttSize();
		int tilesY =     quadCLT_main.getGPU().getImageHeight() / quadCLT_main.getGPU().getDttSize();
		
		if (clt_parameters.gpu_show_geometry) {
			//			GPUTileProcessor.TpTask []
			tp_tasks = quadCLT_main.getGPU().getTasks (false); // boolean use_aux)
			double [][] geom_dbg = new double [ImageDtt.GEOM_TITLES_DBG.length][tilesX*tilesY];
			int num_cams = quadCLT_main.getNumSensors(); // GPUTileProcessor.numSensors; // NUM_CAMS;
			for (int nt = 0; nt < tp_tasks.length; nt++) {
				for (int i = 0; i < num_cams; i++) {
					TpTask task = tp_tasks[nt];
					int nTile = task.ty * tilesX + task.tx;
					geom_dbg[2 * i + 0][nTile] = task.xy[i][0]; // x
					geom_dbg[2 * i + 1][nTile] = task.xy[i][1]; // y
					for (int j = 0; j < 4; j++) {
						geom_dbg[2 * num_cams + 4 * i + j][nTile] = task.disp_dist[i][j];
					}
				}
			}
		    ShowDoubleFloatArrays.showArrays(
		            geom_dbg,
		            tilesX,
		            tilesY,
		            true,
		            name+"-GEOM-DBG-D"+clt_parameters.disparity,
		            ImageDtt.GEOM_TITLES_DBG);
			
		}

		/*
        public TpTask [] getTasks (boolean use_aux)
        {
        	float [] ftasks = new float [TPTASK_SIZE * num_task_tiles];
        	cuMemcpyDtoH(Pointer.to(ftasks), gpu_tasks, TPTASK_SIZE * num_task_tiles * Sizeof.FLOAT);
        	TpTask [] tile_tasks = new TpTask[num_task_tiles];
        	for (int i = 0; i < num_task_tiles; i++) {
        		tile_tasks[i] = new TpTask(ftasks, i* TPTASK_SIZE, use_aux);
        	}
        	return tile_tasks;
        }

		 */
		
		
		// Read extra data for macro generation: 4 DIFFs, 4 of R,  4 of B, 4 of G
		// Available after gpu_diff_rgb_combo is generated in execTextures
		if (calc_extra) {
			String [] extra_group_titles = {"DIFF","Red","Blue","Green"};
			int numSensors = quadCLT_main.getNumSensors();
			String [] extra_titles = new String [extra_group_titles.length* numSensors];
			for (int g = 0; g < extra_group_titles.length;g++) {
				for (int ncam=0; ncam < numSensors;ncam++) {
					extra_titles[g * numSensors + ncam]= extra_group_titles[g]+"-"+ncam;
				}
			}
			float [][] extra = quadCLT_main.getGPU().getExtra();
			ShowDoubleFloatArrays.showArrays(
					extra,
					tilesX,
					tilesY,
					true,
					name+"-EXTRA-D"+clt_parameters.disparity,
					extra_titles);
		}
		/* Prepare 4-channel images*/
		ImagePlus [] imps_RGB = new ImagePlus[iclt_fimg.length];
		for (int ncam = 0; ncam < iclt_fimg.length; ncam++) {
			String title=name+"-"+String.format("%02d", ncam);
			imps_RGB[ncam] = quadCLT_main.linearStackToColor( // probably no need to separate and process the second half with quadCLT_aux
					clt_parameters,
					colorProcParameters,
					rgbParameters,
					title, // String name,
					"-D"+clt_parameters.disparity, //String suffix, // such as disparity=...
					toRGB,
					!quadCLT_main.correctionsParameters.jpeg, // boolean bpp16, // 16-bit per channel color mode for result
					!batch_mode, // true, // boolean saveShowIntermediate, // save/show if set globally
					false, // boolean saveShowFinal,        // save/show result (color image?)
					iclt_fimg[ncam],
					out_width,
					out_height,
					1.0, // scaleExposures[iAux][iSubCam], // double scaleExposure, // is it needed?
					debugLevel );
		}
		
		
		//Show 2D correlations
		int [] wh = new int[2];
		if (clt_parameters.show_corr) {
			int [] corr_quad_indices = quadCLT_main.getGPU().getCorrComboIndices(); // get quad
			float [][] corr2D_quad = quadCLT_main.getGPU().getCorr2DCombo(clt_parameters.gpu_corr_rad);
// calculate and get cross here!			
			quadCLT_main.getGPU().execCorr2D_combine( // calculate cross pairs
			        true, // boolean init_corr,    // initialize output tiles (false - add to current)
			        6,    // int     num_pairs_in, // typically 6 - number of pairs per tile (tile task should have same number per each tile
			        0x30); // int     pairs_mask    // selected pairs (0x3 - horizontal, 0xc - vertical, 0xf - quad, 0x30 - cross)
			
			quadCLT_main.getGPU().execCorr2D_normalize(
	        		true, // boolean combo, // normalize combo correlations (false - per-pair ones) 
		    		fat_zero, // double fat_zero);
					null, // float [] fcorr_weights, // null or one per correlation tile (num_corr_tiles) to divide fat zero2
		    		clt_parameters.gpu_corr_rad); // int corr_radius
			
			int [] corr_cross_indices = quadCLT_main.getGPU().getCorrComboIndices(); // get quad
			float [][] corr2D_cross = quadCLT_main.getGPU().getCorr2DCombo(clt_parameters.gpu_corr_rad);
			
			double [][] dbg_corr_pairs = GPUTileProcessor.getCorr2DView(
					quadCLT_main.getNumSensors(),
					tilesX,
					tilesY,
					corr_indices,
					corr2D,
					wh);
			double [][] dbg_corr_quad = GPUTileProcessor.getCorr2DView(
					quadCLT_main.getNumSensors(),
					tilesX,
					tilesY,
					corr_quad_indices,
					corr2D_quad,
					wh);

			double [][] dbg_corr_cross = GPUTileProcessor.getCorr2DView(
					quadCLT_main.getNumSensors(),
					tilesX,
					tilesY,
					corr_cross_indices,
					corr2D_cross,
					wh);
			
			double [][] dbg_corr = {
					dbg_corr_pairs[0],
					dbg_corr_pairs[1],
					dbg_corr_pairs[2],
					dbg_corr_pairs[3],
					dbg_corr_pairs[4],
					dbg_corr_pairs[5],
					dbg_corr_quad[15],
					dbg_corr_cross[48]
			};

			ShowDoubleFloatArrays.showArrays(
					dbg_corr,
					wh[0],
					wh[1],
					true,
					name+"-CORR2D-D"+clt_parameters.disparity,
					GPUTileProcessor.getCorrTitles());
		}
		
// convert to overlapping and show
		if (clt_parameters.gen_chn_img) { // save and show 4-slice image
			// combine to a sliced color image
			// assuming total number of images to be multiple of 4
			//			  int [] slice_seq = {0,1,3,2}; //clockwise
			int [] slice_seq = new int[results.length];
			for (int i = 0; i < slice_seq.length; i++) {
				slice_seq[i] = i ^ ((i >> 1) & 1); // 0,1,3,2,4,5,7,6, ...
			}
			int width = imps_RGB[0].getWidth();
			int height = imps_RGB[0].getHeight();
			ImageStack array_stack=new ImageStack(width,height);
			for (int i = 0; i<slice_seq.length; i++){
				if (imps_RGB[slice_seq[i]] != null) {
					array_stack.addSlice("port_"+slice_seq[i], imps_RGB[slice_seq[i]].getProcessor().getPixels());
				} else {
					array_stack.addSlice("port_"+slice_seq[i], results[slice_seq[i]].getProcessor().getPixels());
				}
			}
			ImagePlus imp_stack = new ImagePlus(name+"GPU-SHIFTED-D"+clt_parameters.disparity, array_stack);
			imp_stack.getProcessor().resetMinAndMax();
			if (!batch_mode) {
				imp_stack.updateAndDraw();
			}
			quadCLT_main.eyesisCorrections.saveAndShowEnable(
					imp_stack,  // ImagePlus             imp,
					quadCLT_main.correctionsParameters, // EyesisCorrectionParameters.CorrectionParameters  correctionsParameters,
					true, // boolean               enableSave,
					!batch_mode) ;// boolean               enableShow);
		}

		if (clt_parameters.gen_4_img) { // save 4 JPEG images
			// Save as individual JPEG images in the model directory
			String x3d_path= quadCLT_main.getX3dDirectory();
			for (int sub_img = 0; sub_img < imps_RGB.length; sub_img++){
				EyesisCorrections.saveAndShow(
						imps_RGB[sub_img],
						x3d_path,
						quadCLT_main.correctionsParameters.png && !clt_parameters.black_back,
						!batch_mode && clt_parameters.show_textures,
						quadCLT_main.correctionsParameters.JPEG_quality, // jpegQuality); // jpegQuality){//  <0 - keep current, 0 - force Tiff, >0 use for JPEG
						(debugLevel > 0) ? debugLevel : 1); // int debugLevel (print what it saves)
			}
			String model_path= quadCLT_main.correctionsParameters.selectX3dDirectory(
					name, // quad timestamp. Will be ignored if correctionsParameters.use_x3d_subdirs is false
					null,
					true,  // smart,
					true);  //newAllowed, // save

			quadCLT_main.createThumbNailImage(
					imps_RGB[0],
					model_path,
					"thumb",
					debugLevel);

		}
		// Use GPU prepared RBGA
		if (clt_parameters.show_rgba_color) {
			Rectangle woi = new Rectangle();
			float [][] rbga = quadCLT_main.getGPU().getRBGA(
					(is_mono?1:3), // int     num_colors,
					woi);
			ShowDoubleFloatArrays.showArrays( // show slices RBGA (colors - 256, A - 1.0)
					rbga,
					woi.width,
					woi.height,
					true,
					name+"-RGBA-STACK-D"+clt_parameters.disparity+
					":"+clt_parameters.gpu_woi_tx+":"+clt_parameters.gpu_woi_ty+
					":"+clt_parameters.gpu_woi_twidth+":"+clt_parameters.gpu_woi_theight+
					":"+(clt_parameters.gpu_woi_round?"C":"R"),
					new String[] {"R","B","G","A"}
					);



			// for now - use just RGB. Later add option for RGBA
			float [][] rgb_main = {rbga[0],rbga[1],rbga[2]};
			float [][] rgba_main = {rbga[0],rbga[1],rbga[2],rbga[3]};
			ImagePlus imp_rgba_main = quadCLT_main.linearStackToColor(
					clt_parameters,
					colorProcParameters,
					rgbParameters,
					name+"-texture", // String name,
					"-D"+clt_parameters.disparity+"-MAINGPU", //String suffix, // such as disparity=...
					toRGB,
					!quadCLT_main.correctionsParameters.jpeg, // boolean bpp16, // 16-bit per channel color mode for result
					false, // true, // boolean saveShowIntermediate, // save/show if set globally
					false, // true, // boolean saveShowFinal,        // save/show result (color image?)
					((clt_parameters.alpha1 > 0)? rgba_main: rgb_main),
					woi.width, // clt_parameters.gpu_woi_twidth * image_dtt.transform_size, //  tilesX *  image_dtt.transform_size,
					woi.height, // clt_parameters.gpu_woi_theight *image_dtt.transform_size, //  tilesY *  image_dtt.transform_size,
					1.0,         // double scaleExposure, // is it needed?
					debugLevel );

			int width = imp_rgba_main.getWidth();
			int height =imp_rgba_main.getHeight();
			ImageStack texture_stack=new ImageStack(width,height);
			texture_stack.addSlice("main",      imp_rgba_main.getProcessor().getPixels()); // single slice
			ImagePlus imp_texture_stack = new ImagePlus(
					name+"-RGBA-D"+clt_parameters.disparity+
					":"+clt_parameters.gpu_woi_tx+":"+clt_parameters.gpu_woi_ty+
					":"+clt_parameters.gpu_woi_twidth+":"+clt_parameters.gpu_woi_theight+
					":"+(clt_parameters.gpu_woi_round?"C":"R"),
					texture_stack);
			imp_texture_stack.getProcessor().resetMinAndMax();
			String results_path= quadCLT_main.correctionsParameters.selectResultsDirectory( // selectX3dDirectory(
					true,  // smart,
					true);  //newAllowed, // save
			EyesisCorrections.saveAndShow( // save and show color RGBA texture
					imp_texture_stack,
					results_path,
					true, // quadCLT_main.correctionsParameters.png && !clt_parameters.black_back,
					true, // !batch_mode && clt_parameters.show_textures,
					0, // quadCLT_main.correctionsParameters.JPEG_quality, // jpegQuality); // jpegQuality){//  <0 - keep current, 0 - force Tiff, >0 use for JPEG
					(debugLevel > 0) ? debugLevel : 1); // int debugLevel (print what it saves)
		}
		

		// convert textures to RGBA in Java
//		if (clt_parameters.show_rgba_color && (debugLevel > 100)) { // disabling
		if (clt_parameters.show_rgba_color && clt_parameters.gpu_show_jtextures) {  // seem to have wrong colors
			
			int numcol = quadCLT_main.isMonochrome()?1:3;
			int ports = imp_quad_main.length;
			int [] texture_indices = quadCLT_main.getGPU().getTextureIndices();
			int          num_src_slices = numcol + 1 + (clt_parameters.keep_weights?(ports + numcol + 1):0); // 12 ; // calculate
			float [] flat_textures =  quadCLT_main.getGPU().getFlatTextures( // fatal error has been detected by the Java Runtime Environment:
					texture_indices.length,
		    		(is_mono?1:3), // int     num_colors,
		    		clt_parameters.keep_weights); // boolean keep_weights);
	    	int texture_slice_size = (2 * quadCLT_main.getGPU().getDttSize())* (2 * quadCLT_main.getGPU().getDttSize());
	    	int texture_tile_size = texture_slice_size * num_src_slices ;

			if (debugLevel > -1) {
		    	for (int indx = 0; indx < texture_indices.length; indx++) if ((texture_indices[indx] & (1 << GPUTileProcessor.LIST_TEXTURE_BIT)) != 0){
		    		int tile = texture_indices[indx] >> GPUTileProcessor.TEXT_NTILE_SHIFT;
		    		int tileX = tile % tilesX;
		    		int tileY = tile / tilesX;
		    		if ((tileY == clt_parameters.tileY) && (tileX == clt_parameters.tileX)) {

		    			System.out.println("=== tileX= "+tileX+" tileY= "+tileY+" tile="+tile+" ===");

		    			for (int slice =0; slice < num_src_slices; slice++) {
		    				System.out.println("=== Slice="+slice+" ===");
		    				for (int i = 0; i < 2 * quadCLT_main.getGPU().getDttSize(); i++) {
		    					for (int j = 0; j < 2 * quadCLT_main.getGPU().getDttSize(); j++) {
		    						System.out.print(String.format("%10.4f ",
		    								flat_textures[indx*texture_tile_size + slice* texture_slice_size + 2 * quadCLT_main.getGPU().getDttSize() * i + j]));
		    					}
		    					System.out.println();
		    				}
		    			}
		    		}
		    	}
			}
			double [][][][] texture_tiles = new double [tilesY][tilesX][][];
			quadCLT_main.getGPU();
			GpuQuad.doubleTextures(
		    		new Rectangle(0, 0, tilesX, tilesY), // Rectangle    woi,
		    		texture_tiles,                       // double [][][][] texture_tiles, // null or [tilesY][tilesX]
		    		texture_indices,                     // int []       indices,
		    		flat_textures,                       // float [][][] ftextures,
		    		tilesX,                              // int          full_width,
		    		4,                                   // rbga only /int          num_slices
		    		num_src_slices                       // int          num_src_slices
		    		);

			if ((debugLevel > -1) && (clt_parameters.tileX >= 0) && (clt_parameters.tileY >= 0) && (clt_parameters.tileX < tilesX) && (clt_parameters.tileY < tilesY)) {
				String [] rgba_titles = {"red","blue","green","alpha"};
				String [] rgba_weights_titles = {"red","blue","green","alpha","port0","port1","port2","port3","r-rms","b-rms","g-rms","w-rms"};
				double [][] texture_tile = texture_tiles[clt_parameters.tileY][clt_parameters.tileX];
				int tile = +clt_parameters.tileY * tilesX  +clt_parameters.tileX;
    			System.out.println("=== tileX= "+clt_parameters.tileX+" tileY= "+clt_parameters.tileY+" tile="+tile+" ===");

    			for (int slice =0; slice < texture_tile.length; slice++) {
    				System.out.println("\n=== Slice="+slice+" ===");
    				for (int i = 0; i < 2 * quadCLT_main.getGPU().getDttSize(); i++) {
    					for (int j = 0; j < 2 * quadCLT_main.getGPU().getDttSize(); j++) {
    						System.out.print(String.format("%10.4f ",
    								texture_tile[slice][2 * quadCLT_main.getGPU().getDttSize() * i + j]));
    					}
    					System.out.println();
    				}
    			}
    			ShowDoubleFloatArrays.showArrays(
						texture_tile,
						2 * image_dtt.transform_size,
						2 * image_dtt.transform_size,
						true,
						name + "-TXTNOL-GPU-D"+clt_parameters.disparity+"-X"+clt_parameters.tileX+"-Y"+clt_parameters.tileY,
						(clt_parameters.keep_weights?rgba_weights_titles:rgba_titles));

			}


			int alpha_index = 3;
			// in monochrome mode only MONO_CHN == GREEN_CHN is used, R and B are null
			double [][] texture_overlap_main = image_dtt.combineRBGATiles(
					texture_tiles,                 // array [tp.tilesY][tp.tilesX][4][4*transform_size] or [tp.tilesY][tp.tilesX]{null}
//					image_dtt.transform_size,
					true,                         // when false - output each tile as 16x16, true - overlap to make 8x8
					clt_parameters.sharp_alpha,    // combining mode for alpha channel: false - treat as RGB, true - apply center 8x8 only
					threadsMax,                    // maximal number of threads to launch
					debugLevel);
			if (clt_parameters.alpha1 > 0){ // negative or 0 - keep alpha as it was
				double scale = (clt_parameters.alpha1 > clt_parameters.alpha0) ? (1.0/(clt_parameters.alpha1 - clt_parameters.alpha0)) : 0.0;
				for (int i = 0; i < texture_overlap_main[alpha_index].length; i++){
					double d = texture_overlap_main[alpha_index][i];
					if      (d >=clt_parameters.alpha1) d = 1.0;
					else if (d <=clt_parameters.alpha0) d = 0.0;
					else d = scale * (d- clt_parameters.alpha0);
					texture_overlap_main[alpha_index][i] = d;
				}
			}

			// for now - use just RGB. Later add option for RGBA
			double [][] texture_rgb_main = {texture_overlap_main[0],texture_overlap_main[1],texture_overlap_main[2]};
			double [][] texture_rgba_main = {texture_overlap_main[0],texture_overlap_main[1],texture_overlap_main[2],texture_overlap_main[3]};
			ImagePlus imp_texture_main = quadCLT_main.linearStackToColor(
					clt_parameters,
					colorProcParameters,
					rgbParameters,
					name+"-texture", // String name,
					"-D"+clt_parameters.disparity+"-MAINGPU", //String suffix, // such as disparity=...
					toRGB,
					!quadCLT_main.correctionsParameters.jpeg, // boolean bpp16, // 16-bit per channel color mode for result
					false, // true, // boolean saveShowIntermediate, // save/show if set globally
					false, // true, // boolean saveShowFinal,        // save/show result (color image?)
					((clt_parameters.alpha1 > 0)? texture_rgba_main: texture_rgb_main),
					tilesX *  image_dtt.transform_size,
					tilesY *  image_dtt.transform_size,
					1.0,         // double scaleExposure, // is it needed?
					debugLevel );

			int width = imp_texture_main.getWidth();
			int height =imp_texture_main.getHeight();
			ImageStack texture_stack=new ImageStack(width,height);
			texture_stack.addSlice("main",      imp_texture_main.getProcessor().getPixels()); // single slice
			ImagePlus imp_texture_stack = new ImagePlus(name+"-TEXTURES-D"+clt_parameters.disparity, texture_stack);
			imp_texture_stack.getProcessor().resetMinAndMax();
			imp_texture_stack.show();
		}

		return results;
	}
	
	
	public static boolean savePortsXY(
			EyesisCorrectionParameters.CorrectionParameters ecp,
			TpTask [] tp_tasks
			) {
		if ((ecp.tile_processor_gpu != null) && !ecp.tile_processor_gpu.isEmpty()) {
			int quad = 4;
			String file_prefix = ecp.tile_processor_gpu +"clt/main";
			for (int chn = 0; chn < quad; chn++) {
				String img_path =  file_prefix+"_chn"+chn+".portsxy";
				FileOutputStream fos;
				try {
					fos = new FileOutputStream(img_path);
				} catch (FileNotFoundException e) {
					// TODO Auto-generated catch block
					System.out.println("Could not write to "+img_path+" (file not found) port offsets");
					break;
				}
				DataOutputStream dos = new DataOutputStream(fos);
				WritableByteChannel channel = Channels.newChannel(dos);
				ByteBuffer bb = ByteBuffer.allocate(tp_tasks.length * 2 * 4);
				bb.order(ByteOrder.LITTLE_ENDIAN);
				bb.clear();
				for (int i = 0; i <  tp_tasks.length; i++) {
					bb.putFloat((tp_tasks[i].xy[chn][0])); // x-offset
					bb.putFloat((tp_tasks[i].xy[chn][1])); // y-offset
				}
				bb.flip();
				try {
					channel.write(bb);
				} catch (IOException e) {
					System.out.println("Could not write to "+img_path+" port offsets");
					break;
				}
				try {
					dos.close();
				} catch (IOException e) {
					System.out.println("Could not close DataOutputStream for "+img_path+" port offsets");
				}
				System.out.println("Wrote port offsets to "+img_path+".");
			}
			return true;
		} else {
			return false;
		}
		
	}

	public ImagePlus getBackgroundImage(
			boolean []                                bgnd_tiles, 
			CLTPass3d                                 bgnd_data,
			CLTParameters                             clt_parameters,
			ColorProcParameters                       colorProcParameters,
			EyesisCorrectionParameters.RGBParameters  rgbParameters,
			String     name,
			//			  int        disparity_index, // index of disparity value in disparity_map == 2 (0,2 or 4)
			int        threadsMax,  // maximal number of threads to launch
			boolean    updateStatus,
			int        debugLevel
			) {
		if ((gpuQuad == null) || !(isAux()?clt_parameters.gpu_use_aux : clt_parameters.gpu_use_main)) {
			return super.getBackgroundImage( // measure background // USED in lwir
					bgnd_tiles, 
					bgnd_data,
					clt_parameters,
					colorProcParameters,
					rgbParameters,
					name,
					//					  disparity_index, // index of disparity value in disparity_map == 2 (0,2 or 4)
					threadsMax,  // maximal number of threads to launch
					updateStatus,
					debugLevel);
		}

		final int tilesX = tp.getTilesX();
		final int tilesY = tp.getTilesY();
		//		  CLTPass3d bgnd_data = tp.clt_3d_passes.get(0);
		boolean [] bgnd_tiles_grown2 = bgnd_data.getSelected().clone(); // only these have non 0 alpha
		tp.growTiles(
				2,      // grow tile selection by 1 over non-background tiles 1: 4 directions, 2 - 8 directions, 3 - 8 by 1, 4 by 1 more
				bgnd_tiles_grown2,
				null); // prohibit
		//				  bgnd_data.getTextureSelection()
		bgnd_data.setTextureSelection(bgnd_tiles_grown2);
		Rectangle  texture_woi_pix = new Rectangle(); // in pixels
		float [][] texture_img = GetTextureGPU( // returns texture CUDA error at com.elphel.imagej.gpu.GpuQuad.execRBGA_noDP(GpuQuad.java:2180)
				clt_parameters,      // CLTParameters       clt_parameters,
				texture_woi_pix,     // Rectangle  texture_woi, // = new Rectangle();
				bgnd_data.disparity, // double [][]         disparity_array, // [tilesY][tilesX]
				bgnd_tiles_grown2,   // bgnd_tiles_grown2,   // boolean []          selection,
				threadsMax,          // final int           threadsMax,  // maximal number of threads to launch
				updateStatus,        // final boolean       updateStatus,
				debugLevel);         // final int           debugLevel);
		int out_width =  tilesX *  clt_parameters.transform_size;
		int out_height = tilesY *  clt_parameters.transform_size;
		ImagePlus imp_texture_full = null;
		if (texture_img == null) {
			System.out.println("No background visible");
			texture_woi_pix = new Rectangle(0,0,1,1);
			texture_img = new float[][] {{0.0f},{0.0f}};
		}

		// for now - use just RGB. Later add option for RGBA
		float [][] texture_rgb = isMonochrome() ? (new float [][] {texture_img[0]}): (new float [][] {texture_img[0],texture_img[1],texture_img[2]});
		float [][] texture_rgba = isMonochrome() ? (new float [][] {texture_img[0],texture_img[1]}) : (new float [][] {texture_img[0],texture_img[1],texture_img[2],texture_img[3]});

		boolean toRGB = !isLwir() || colorProcParameters.lwir_pseudocolor;	
		ImagePlus imp_texture_bgnd = linearStackToColor(
				clt_parameters,
				colorProcParameters,
				rgbParameters,
				name+"-texture-bgnd", // String name,
				"", //String suffix, // such as disparity=...
				toRGB, // toRGB,
				!this.correctionsParameters.jpeg, // boolean bpp16, // 16-bit per channel color mode for result
				true, // boolean saveShowIntermediate, // save/show if set globally
				false, //true, // boolean saveShowFinal,        // save/show result (color image?)
				((clt_parameters.alpha1 > 0)? texture_rgba: texture_rgb),
				texture_woi_pix.width, // tilesX *  clt_parameters.transform_size,
				texture_woi_pix.height, // tilesY *  clt_parameters.transform_size,
				1.0,         // double scaleExposure, // is it needed?
				debugLevel);
		//		  imp_texture_bgnd.show();
		// resize for backdrop here!
		//	public double getFOVPix(){ // get ratio of 1 pixel X/Y to Z (distance to object)
		imp_texture_full = resizeToFull(
				out_width, // int       out_width,
				out_height, // int       out_height,
				texture_woi_pix.x, // int       x0, // image offset-x 
				texture_woi_pix.y, // int       y0, // image offset-y
				imp_texture_bgnd, // ImagePlus imp_texture_bgnd,
				false, // boolean   fillBlack,
				false, // boolean noalpha, // only with fillBlack, otherwise ignored
				debugLevel);
		//		  imp_texture_full.show();
		return imp_texture_full;

	}

	  /**
	   * Get GPU-calculated tiles centers coordinates for each sesnor
	   * @return per-tile null 2*num_sesnors array of alternating X and Y 
	   */
	  public double [][] getPortsXY() // run after GPU 
	  {
		  boolean use_aux = false;
		  final int tilesX = tp.getTilesX();
		  final int tilesY = tp.getTilesY();
		  double [][] ports_xy = new double [tilesY*tilesX][];
		  TpTask [] tp_tasks = gpuQuad.getTasks(use_aux); // use_aux);
		  for (TpTask tp_task:tp_tasks) {
			  float [][] txy =  tp_task.getXY(); // use_aux);
			  double [] tile_ports_xy = new double [txy.length * txy[0].length];
			  int indx = 0;
			  for (int i = 0; i < txy.length; i++) {
				  for (int j = 0; j < txy[i].length; j++) {
					  tile_ports_xy[indx++] = txy[i][j];
				  }
			  }
			  ports_xy[tp_task.getTileY()*tilesX + tp_task.getTileX()] = tile_ports_xy;
		  }
		  return ports_xy;
	  }

	  /**
	   * Get GPU-calculated tiles centerXY (tile centers for the virtual camera)
	   * @return per-tile null or pair of X,Y
	   */
	  public double [][] getCenterXY() // run after GPU 
	  {
		  boolean use_aux = false;
		  final int tilesX = tp.getTilesX();
		  final int tilesY = tp.getTilesY();

		  double [][] centers_xy = new double [tilesY*tilesX][];
		  TpTask [] tp_tasks = gpuQuad.getTasks(use_aux); // use_aux);
		  for (TpTask tp_task:tp_tasks) {
			  centers_xy[tp_task.getTileY()*tilesX + tp_task.getTileX()] = tp_task.getDoubleCenterXY(); // use_aux);
		  }
		  return centers_xy;
	  }

	  public void showPortsXY(String title) {
		  int tilesX = tp.getTilesX();
		  int tilesY = tp.getTilesY();
		  int num_sensors = getNumSensors();
		  String [] titles = new String[2*num_sensors];
		  for (int n = 0; n < num_sensors; n++) {
			  titles[n]= "x"+n;
			  titles[n+num_sensors]= "y"+n;
		  }
		  double [][] ports_xy =  getPortsXY();
		  double [][] data = new double [2*num_sensors][tilesX*tilesY];
		  for (int l = 0; l < data.length; l++)  Arrays.fill(data[l], Double.NaN);
		  for (int i = 0; i < ports_xy.length; i++) if (ports_xy[i]!= null) {
			  for (int j = 0; j < data.length; j++) {
				  data[j][i] = ports_xy[i][2 * (j % num_sensors)  + (j / num_sensors)];
			  }
		  }
		  ShowDoubleFloatArrays.showArrays( 
				  data,
				  tilesX,
				  tilesY,
				  true,
				  title,
				  titles);
	  }

	  public void showCentersXY(String title) {
		  int tilesX = tp.getTilesX();
		  int tilesY = tp.getTilesY();
		  String [] titles = {"centerX","centerY"};
		  double [][] centers_xy =  getCenterXY();
		  double [][] data = new double [2][tilesX*tilesY];
		  for (int l = 0; l < data.length; l++)  Arrays.fill(data[l], Double.NaN);
		  for (int i = 0; i < centers_xy.length; i++) if (centers_xy[i]!= null) {
			  data[0][i] = centers_xy[i][0];
			  data[1][i] = centers_xy[i][1];
		  }
		  ShowDoubleFloatArrays.showArrays( 
				  data,
				  tilesX,
				  tilesY,
				  true,
				  title,
				  titles);
	  }
	  
	  public void showCentersPortsXY(String title) {
		  int tilesX = tp.getTilesX();
		  int tilesY = tp.getTilesY();
		  int num_sensors = getNumSensors();
		  String [] titles = new String[2*num_sensors + 2];
		  titles[0]= "centerX" ;
		  titles[num_sensors + 1]= "centerY" ;
		  for (int n = 0; n < num_sensors; n++) {
			  titles[n + 1]= "x"+n ;
			  titles[n + 1 + (num_sensors + 1)]= "y"+n;
		  }
		  double [][] ports_xy =  getPortsXY();
		  double [][] centers_xy =  getCenterXY();

		  double [][] data = new double [2*(num_sensors + 1)][tilesX*tilesY];
		  for (int l = 0; l < data.length; l++)  Arrays.fill(data[l], Double.NaN);
		  for (int i = 0; i < ports_xy.length; i++) if (ports_xy[i]!= null) {
			  for (int j = 0; j < data.length; j++) {
				  int dir = j / (num_sensors + 1); // 0 - X, 1 - Y
				  int n = j - dir * (num_sensors + 1);
				  if (n == 0) {
					  data[j][i] = centers_xy[i][dir];
				  } else {
					  n --;
					  data[j][i] = ports_xy[i][2 * n  + dir];
				  }
			  }
		  }
		  ShowDoubleFloatArrays.showArrays( 
				  data,
				  tilesX,
				  tilesY,
				  true,
				  title,
				  titles);
	  }

	  
	  
	  
	  public double [][] CLTCorrDisparityMap(
			  CLTParameters     clt_parameters,
			  double []         disparity, // do not measure NaN
			  double            disparity_corr,
			  final int         threadsMax,  // maximal number of threads to launch
			  final boolean     updateStatus,
			  final int         debugLevel)
	  {
		  if ((gpuQuad == null) || !(isAux()?clt_parameters.gpu_use_aux : clt_parameters.gpu_use_main)) {
			  return null; // no non-GPU 
		  }
		  final int tilesX = tp.getTilesX();
		  final int tilesY = tp.getTilesY();
		  int op = ImageDtt.setImgMask(0, 0xf);
		  op =     ImageDtt.setPairMask(op,0xf);
		  op =     ImageDtt.setForcedDisparity(op,true);
		  int [][] tile_op = new int [tilesY][tilesX];
		  double [][]  disparity_array = new double[tilesY][tilesX];
		  for (int ty = 0; ty < tilesY; ty++) {
			  for (int tx = 0; tx < tilesX; tx++) {
				  int indx = ty*tilesX+tx;
				  if (!Double.isNaN(disparity[indx])) {
					  tile_op        [ty][tx] = op; // 511; // d;
					  disparity_array[ty][tx] = disparity[indx];
					  
				  }
			  }
		  }
//		  double [][] disparity_map = new double [ImageDtt.DISPARITY_TITLES.length][];
		  double [][] disparity_map = new double [ImageDtt.getDisparityTitles(getNumSensors(), isMonochrome()).length][];

		  ImageDtt image_dtt = new ImageDtt(
				  getNumSensors(),
				  clt_parameters.transform_size,
				  clt_parameters.img_dtt,
				  isAux(),
				  isMonochrome(),
				  isLwir(),
				  clt_parameters.getScaleStrength(isAux()),
				  gpuQuad);
		  image_dtt.getCorrelation2d(); // initiate image_dtt.correlation2d, needed if disparity_map != null  
//		  double z_correction =  clt_parameters.z_correction;
//		  if (clt_parameters.z_corr_map.containsKey(image_name)){ // not used in lwir
//			  z_correction +=clt_parameters.z_corr_map.get(image_name);
//		  }
		  // final double disparity_corr = (z_correction == 0) ? 0.0 : geometryCorrection.getDisparityFromZ(1.0/z_correction);
		  // TODO: find where BITS_ALL_DIFFS is not needed and remove - it takes much longer than correlation in the GPU
		  int disparity_modes = 
				  ImageDtt.BITS_ALL_DISPARITIES |
				  ImageDtt.BITS_ALL_DIFFS | // needs max_diff?
				  ImageDtt.BITS_OVEREXPOSED; //  |
//				  ImageDtt.BITS_TONE_RGB;

		  image_dtt.clt_aberrations_quad_corr_GPU( // USED in LWIR
				  clt_parameters.img_dtt,        // final ImageDttParameters  imgdtt_params,   // Now just extra correlation parameters, later will include, most others
				  1,                             // final int  macro_scale, // to correlate tile data instead of the pixel data: 1 - pixels, 8 - tiles
				  tile_op,                       // per-tile operation bit codes
				  disparity_array,               // clt_parameters.disparity,     // final double            disparity,
				  null,							 // final float  [][][][]     corr_td,         // [tilesY][tilesX][pair][4*64] transform domain representation of 6 corr pairs
				  null, 						 //	final float  [][][][]     corr_combo_td,   // [4][tilesY][tilesX][pair][4*64] TD of combo corrs: qud, cross, hor,vert
					                             // each of the top elements may be null to skip particular combo type
				  null,                          // final float  [][][][]     fdisp_dist,      // [tilesY][tilesX][cams][4], // disparity derivatives vectors or null
				  null,                          //	final float  [][][][]     fpxpy,           // [tilesY][tilesX][cams][2], tile {pX,pY}
				  //// Uses quadCLT from gpuQuad			
				                                 // correlation results - final and partial
				  null,                          // [type][tilesY][tilesX][(2*transform_size-1)*(2*transform_size-1)] // if null - will not calculate
				                                 // [type][tilesY][tilesX] should be set by caller
				                                 // types: 0 - selected correlation (product+offset), 1 - sum
				  null,                          // clt_corr_partial,// [tilesY][tilesX][quad]color][(2*transform_size-1)*(2*transform_size-1)] // if null - will not calculate
				                                 // [tilesY][tilesX] should be set by caller
				                                 // When clt_mismatch is non-zero, no far objects extraction will be attempted
				  null,                          // clt_mismatch,    // [12][tilesY * tilesX] // ***** transpose unapplied ***** ?. null - do not calculate
				                                 // values in the "main" directions have disparity (*_CM) subtracted, in the perpendicular - as is
				  disparity_map,                 // disparity_map,   // [8][tilesY][tilesX], only [6][] is needed on input or null - do not calculate
				                                 // last 2 - contrast, avg/ "geometric average)
				  disparity_modes,               // disparity_modes, // bit mask of disparity_map slices to calculate/return
				  null,                          // 	final double [][][][]     texture_tiles,   // compatible with the CPU ones      
				  null, // texture_img,          // texture_img,     // null or [3][] (RGB) or [4][] RGBA
				  null, // texture_woi,          // texture_woi,     // null or generated texture location/size
				  null,                          //  iclt_fimg,       // will return quad images or null to skip, use quadCLT.linearStackToColor 
				  clt_parameters.gpu_corr_scale, //  gpu_corr_scale,  //  0.75; // reduce GPU-generated correlation values
				  clt_parameters.getGpuFatZero(isMonochrome()), // final double              gpu_fat_zero,    // clt_parameters.getGpuFatZero(is_mono);absolute == 30.0\
				  clt_parameters.gpu_sigma_r,    // 0.9, 1.1
				  clt_parameters.gpu_sigma_b,    // 0.9, 1.1
				  clt_parameters.gpu_sigma_g,    // 0.6, 0.7
				  clt_parameters.gpu_sigma_m,    //  =       0.4; // 0.7;
				  (isMonochrome()? 1.0 : clt_parameters.gpu_sigma_rb_corr), // final double              gpu_sigma_rb_corr, //  = 0.5; // apply LPF after accumulating R and B correlation before G, monochrome ? 1.0 : gpu_sigma_rb_corr;
				  clt_parameters.gpu_sigma_corr, //  =    0.9;gpu_sigma_corr_m
				  image_dtt.transform_size - 1, // clt_parameters.gpu_corr_rad,   // = transform_size - 1 ?
				  clt_parameters.corr_red,       // +used
				  clt_parameters.corr_blue,      // +used
				  clt_parameters.max_corr_radius,// 3.9;
				  clt_parameters.corr_mode,      // Correlation mode: 0 - integer max, 1 - center of mass, 2 - polynomial
				  clt_parameters.min_shot,       // 10.0;  // Do not adjust for shot noise if lower than
				  clt_parameters.scale_shot,     // 3.0;   // scale when dividing by sqrt ( <0 - disable correction)
				  clt_parameters.diff_sigma,     // 5.0;//RMS difference from average to reduce weights (~ 1.0 - 1/255 full scale image)
				  clt_parameters.diff_threshold, // 5.0;   // RMS difference from average to discard channel (~ 1.0 - 1/255 full scale image)
				  clt_parameters.diff_gauss,     // true;  // when averaging images, use gaussian around average as weight (false - sharp all/nothing)
				  clt_parameters.min_agree,      // 3.0;   // minimal number of channels to agree on a point (real number to work with fuzzy averages)
				  clt_parameters.dust_remove,    // Do not reduce average weight when only one image differes much from the average
				  geometryCorrection,            // final GeometryCorrection  geometryCorrection,
				  null,                          // final GeometryCorrection  geometryCorrection_main, // if not null correct this camera (aux) to the coordinates of the main
				  clt_parameters.clt_window,
				  disparity_corr,                // final double              disparity_corr, // disparity at infinity
				  clt_parameters.tileX,          // final int               debug_tileX,
				  clt_parameters.tileY,          // final int               debug_tileY,
				  threadsMax,
				  debugLevel);
		  return disparity_map;
	  }
	  
	  
	  public String getPassImage( // get image from a single pass, return relative path for x3d // USED in lwir
			  CLTParameters           clt_parameters,
			  ColorProcParameters colorProcParameters,
			  EyesisCorrectionParameters.RGBParameters             rgbParameters,
			  String     name,
			  int        scanIndex,
			  int        threadsMax,  // maximal number of threads to launch
			  boolean    updateStatus,
			  int        debugLevel)
	  {
		  if ((gpuQuad == null) || !(isAux()?clt_parameters.gpu_use_aux : clt_parameters.gpu_use_main)) {
			  return super.getPassImage( // get image from a single pass, return relative path for x3d // USED in lwir
					  clt_parameters,
					  colorProcParameters,
					  rgbParameters,
					  name,
					  scanIndex,
					  threadsMax,  // maximal number of threads to launch
					  updateStatus,
					  debugLevel);
		  }

		  //		  final int tilesX = tp.getTilesX();
		  //		  final int tilesY = tp.getTilesY();
		  final int transform_size =clt_parameters.transform_size;
		  CLTPass3d scan = tp.clt_3d_passes.get(scanIndex);
		  Rectangle  texture_woi_pix = new Rectangle(); // in pixels
		  float [][] texture_img = GetTextureGPU( // returns texture
				  clt_parameters,             // CLTParameters       clt_parameters,
				  texture_woi_pix,                // Rectangle  texture_woi, // = new Rectangle();
				  scan.disparity,             // double [][]         disparity_array, // [tilesY][tilesX]
				  scan.getTextureSelection(), // bgnd_tiles_grown2,   // boolean []          selection,
				  threadsMax,                 // final int           threadsMax,  // maximal number of threads to launch
				  updateStatus,               // final boolean       updateStatus,
				  debugLevel);                // final int           debugLevel);

		  // for now - use just RGB. Later add option for RGBA
		  float [][] texture_rgb = {texture_img[0],texture_img[1],texture_img[2]};
		  float [][] texture_rgba = {texture_img[0],texture_img[1],texture_img[2],texture_img[3]};
		  float [][] texture_rgbx = ((clt_parameters.alpha1 > 0)? texture_rgba: texture_rgb);
		  scan.texture_bounds = new Rectangle(
				  texture_woi_pix.x/transform_size,
				  texture_woi_pix.y/transform_size,
				  texture_woi_pix.width/transform_size,
				  texture_woi_pix.height/transform_size);
		  scan.setSelected(scan.getTextureSelection().clone()); // null
		  boolean toRGB = !isLwir() || colorProcParameters.lwir_pseudocolor;	
		  ImagePlus imp_texture_cluster = linearStackToColor(
				  clt_parameters,
				  colorProcParameters,
				  rgbParameters,
				  name+"-texture", // String name,
				  "", //String suffix, // such as disparity=...
				  toRGB, // true, // toRGB,
				  !this.correctionsParameters.jpeg, // boolean bpp16, // 16-bit per channel color mode for result
				  true, // boolean saveShowIntermediate, // save/show if set globally
				  false, //true, // boolean saveShowFinal,        // save/show result (color image?)
				  texture_rgbx,
				  scan.texture_bounds.width * transform_size, // tilesX *  clt_parameters.transform_size,
				  scan.texture_bounds.height * transform_size, // tilesY *  clt_parameters.transform_size,
				  1.0,         // double scaleExposure, // is it needed?
				  debugLevel);

		  String path= correctionsParameters.selectX3dDirectory(
				  //TODO: Which one to use - name or this.image_name ?
				  correctionsParameters.getModelName(this.image_name), // quad timestamp. Will be ignored if correctionsParameters.use_x3d_subdirs is false
				  correctionsParameters.x3dModelVersion,
				  true,  // smart,
				  true);  //newAllowed, // save
		  // only show/save original size if debug or debug_filters)
		  eyesisCorrections.saveAndShow(
				  imp_texture_cluster,
				  path,
				  correctionsParameters.png,
				  clt_parameters.show_textures,
				  -1); // jpegQuality){//  <0 - keep current, 0 - force Tiff, >0 use for JPEG
		  
		  scan.is_measured =   true; // should it be here?
		  scan.is_combo =      false;
		  return imp_texture_cluster.getTitle()+".png"; // imp_texture_cluster;
	  }	

	  public void setPassAvgRBGA( // get image from a single pass, return relative path for x3d // USED in lwir
			  CLTParameters           clt_parameters,
			  CLTPass3d  scan,			  
			  int        threadsMax,  // maximal number of threads to launch
			  boolean    updateStatus,
			  int        debugLevel)
	  {
		  if ((gpuQuad != null) && (isAux()?clt_parameters.gpu_use_aux : clt_parameters.gpu_use_main)) {
			  final int tilesX = scan.getTileProcessor().getTilesX();
			  final int tilesY = scan.getTileProcessor().getTilesY();
			  double [] disparity = scan.getDisparity();
			  double [] strength =  scan.getStrength();
			  boolean [] selection = null; // scan.getSelected();
			  if (selection == null) {
				  selection = new boolean[tilesX*tilesY];
				  for (int nTile = 0; nTile < selection.length; nTile++) {
					  selection[nTile] = !Double.isNaN(disparity[nTile]) && (strength[nTile] > 0.0);
				  }
				  scan.setSelected(selection);
			  }
			  if ((scan.disparity == null) || (scan.tile_op == null)) {
				  int d = ImageDtt.setImgMask(0, 0xf); // no correlations
				  d =     ImageDtt.setForcedDisparity(d,true);
				  // New as for 11/18/2022 - no CPU support yet 
				  d |= (1 << GPUTileProcessor.TASK_TEXT_EN);
				  scan.setTileOpDisparity(
						  d,
						  scan.getSelected(), // boolean [] selection,
						  scan.getDisparity()); // double []  disparity)
			  }
			  
			  scan.texture_tiles = getTextureTilesGPU( // returns texture tiles, as with CPU
					  clt_parameters,     // CLTParameters       clt_parameters,
					  scan.tile_op,       // int [][]            tile_op,
					  scan.disparity,             // double [][]         disparity_array, // [tilesY][tilesX]
					  threadsMax,                 // final int           threadsMax,  // maximal number of threads to launch
					  updateStatus,               // final boolean       updateStatus,
					  debugLevel);                // final int           debugLevel);
		  } else {
			  CLTMeasureTextures( // perform single pass according to prepared tiles operations and disparity // not used in lwir
					  clt_parameters,   // final CLTParameters           clt_parameters,
					  scan,             // final CLTPass3d   scan,
					  threadsMax,
					  updateStatus,
					  debugLevel);
			  
		  }
		  super.setPassAvgRBGA( // get image from a single pass, return relative path for x3d // USED in lwir
				  clt_parameters,
				  scan,         // scanIndex,
				  threadsMax,   // maximal number of threads to launch
				  updateStatus,
				  debugLevel);
	  }
	  
			  /*
		  final int tilesX = tp.getTilesX();
		  final int tilesY = tp.getTilesY();
		  final int transform_size =clt_parameters.transform_size;
		  CLTPass3d scan = tp.clt_3d_passes.get(scanIndex);
		  Rectangle  texture_woi = new Rectangle();
// need	to set tasks first!	  
		  float [][] texture_img = GetTextureGPU( // returns texture 
				  clt_parameters,             // CLTParameters       clt_parameters,
				  texture_woi,                // Rectangle  texture_woi, // = new Rectangle();
				  scan.disparity,             // double [][]         disparity_array, // [tilesY][tilesX]
				  scan.getTextureSelection(), // bgnd_tiles_grown2,   // boolean []          selection,
				  threadsMax,                 // final int           threadsMax,  // maximal number of threads to launch
				  updateStatus,               // final boolean       updateStatus,
				  debugLevel);                // final int           debugLevel);

		  // for now - use just RGB. Later add option for RGBA
		  float [][] texture_rgb = {texture_img[0],texture_img[1],texture_img[2]};
		  float [][] texture_rgba = {texture_img[0],texture_img[1],texture_img[2],texture_img[3]};
		  float [][] texture_rgbx = ((clt_parameters.alpha1 > 0)? texture_rgba: texture_rgb);
		  scan.texture_bounds = new Rectangle(
				  texture_woi.x * transform_size,
				  texture_woi.y * transform_size,
				  texture_woi.width * transform_size,
				  texture_woi.height * transform_size);
		  scan.selected = scan.getTextureSelection().clone();
		  int numTiles = tilesX * tilesY; 
		  int num_layers = texture_rgbx.length; // here only 3/4
		  double [][] tileTones = new double [num_layers][numTiles];
		  double [] scales = new double [num_layers];
		  for (int n = 0; n < num_layers; n++){
			  if       (n < 3)  scales[n] = 1.0/255.0; // R,B,G
			  else if  (n == 3) scales[n] = 1.0; //alpha
			  else if  (n < 8)  scales[n] = 1.0; // ports 0..3
			  else              scales[n] = 1.0/255.0; // RBG rms, in 1/255 units, but small
		  }
		  int tx0 = scan.texture_bounds.x;
		  int ty0 = scan.texture_bounds.y;
		  int tw =  scan.texture_bounds.width;
		  int th =  scan.texture_bounds.height;
		  int text_step = transform_size * tw;
		  for (int ty = ty0; ty < (ty0 + th); ty++) {
			  for (int tx = tx0; tx < (tx0 + tw); tx++) {
				  int text_start = ((ty - ty0) * transform_size * tw  + (tx-tx0)) * transform_size;
				  double s = 0.0;
				  for (int n = 0; n < num_layers; n++) {
					  for (int i = 0; i < transform_size; i++) {
						  for (int j = 0; j < transform_size; j++) {
							 s += texture_rgbx[n][text_start + i * text_step + j]; 
						  }
					  }
					  tileTones[n][ty * tilesX + tx] = s * scales[n] /(transform_size * transform_size);
				  }
			  }
		  }
		  scan.setTilesRBGA(tileTones); 
		  ShowDoubleFloatArrays.showArrays( // out of boundary 15
				  tileTones,
				  tilesX,
				  tilesY,
				  true,
				  image_name+sAux()+"-TileTones");
	  }			  
*/


	  
	  public float [][] GetTextureGPU( // returns texture
			  CLTParameters       clt_parameters,
			  Rectangle           texture_woi_pix, // = new Rectangle(); in pixels!
			  double [][]         disparity_array, // [tilesY][tilesX]
			  boolean []          selection,
			  final int           threadsMax,  // maximal number of threads to launch
			  final boolean       updateStatus,
			  final int           debugLevel)
	  {
			  
		  final int tilesX = tp.getTilesX();
		  final int tilesY = tp.getTilesY();
		  int d = ImageDtt.setImgMask(0, 0xf);
		  d =     ImageDtt.setForcedDisparity(d,true);
		  int [][] tile_op = new int [tilesY][tilesX];
		  int indx = 0;
		  if (selection != null) {
			  for (int ty = 0; ty < tilesY; ty++) {
				  for (int tx = 0; tx < tilesX; tx++) {
					  if (selection[indx++]) {
						  tile_op[ty][tx] = d;
					  }
				  }
			  }
		  } else {
			  for (int ty = 0; ty < tilesY; ty++) {
				  for (int tx = 0; tx < tilesX; tx++) {
					  tile_op[ty][tx] = d;
				  }
			  }
		  }

		  /*
		  double [][] shiftXY = new double [4][2];
		  // not used
		  if (!clt_parameters.fine_corr_ignore) {
			  double [][] shiftXY0 = {
					  {clt_parameters.fine_corr_x_0,clt_parameters.fine_corr_y_0},
					  {clt_parameters.fine_corr_x_1,clt_parameters.fine_corr_y_1},
					  {clt_parameters.fine_corr_x_2,clt_parameters.fine_corr_y_2},
					  {clt_parameters.fine_corr_x_3,clt_parameters.fine_corr_y_3}};
			  shiftXY = shiftXY0;
		  }
          */
		  ImageDtt image_dtt = new ImageDtt(
				  getNumSensors(),
				  clt_parameters.transform_size,
				  clt_parameters.img_dtt,
				  isAux(),
				  isMonochrome(),
				  isLwir(),
				  clt_parameters.getScaleStrength(isAux()),
				  gpuQuad);
		  image_dtt.getCorrelation2d(); // initiate image_dtt.correlation2d, needed if disparity_map != null  
		  double z_correction =  clt_parameters.z_correction;
		  if (clt_parameters.z_corr_map.containsKey(image_name)){ // not used in lwir
			  z_correction +=clt_parameters.z_corr_map.get(image_name);
		  }
		  final double disparity_corr = (z_correction == 0) ? 0.0 : geometryCorrection.getDisparityFromZ(1.0/z_correction);
		  float [][] texture_img = new float [isMonochrome()?2:4][];
		  image_dtt.clt_aberrations_quad_corr_GPU( // USED in LWIR
				  clt_parameters.img_dtt,        // final ImageDttParameters  imgdtt_params,   // Now just extra correlation parameters, later will include, most others
				  1,                             // final int  macro_scale, // to correlate tile data instead of the pixel data: 1 - pixels, 8 - tiles
				  tile_op,                       // per-tile operation bit codes
				  disparity_array,               // clt_parameters.disparity,     // final double            disparity,
				  null,							 // final float  [][][][]     corr_td,         // [tilesY][tilesX][pair][4*64] transform domain representation of 6 corr pairs
				  null, 						 //	final float  [][][][]     corr_combo_td,   // [4][tilesY][tilesX][pair][4*64] TD of combo corrs: qud, cross, hor,vert
					                             // each of the top elements may be null to skip particular combo type
				  null,                          // final float  [][][][]    fdisp_dist,      // [tilesY][tilesX][cams][4], // disparity derivatives vectors or null
				  null,                          //	final float  [][][][]     fpxpy,           // [tilesY][tilesX][cams][2], tile {pX,pY}
				  //// Uses quadCLT from gpuQuad			
				                                 // correlation results - final and partial
				  null,                          // [type][tilesY][tilesX][(2*transform_size-1)*(2*transform_size-1)] // if null - will not calculate
				                                 // [type][tilesY][tilesX] should be set by caller
				                                 // types: 0 - selected correlation (product+offset), 1 - sum
				  null,                          // clt_corr_partial,// [tilesY][tilesX][quad]color][(2*transform_size-1)*(2*transform_size-1)] // if null - will not calculate
				                                 // [tilesY][tilesX] should be set by caller
				                                 // When clt_mismatch is non-zero, no far objects extraction will be attempted
				  null,                          // clt_mismatch,    // [12][tilesY * tilesX] // ***** transpose unapplied ***** ?. null - do not calculate
				                                 // values in the "main" directions have disparity (*_CM) subtracted, in the perpendicular - as is
				  null,                          // disparity_map,   // [8][tilesY][tilesX], only [6][] is needed on input or null - do not calculate
				                                 // last 2 - contrast, avg/ "geometric average)
				  0,                             // disparity_modes, // bit mask of disparity_map slices to calculate/return
				  null,                          // 	final double [][][][]     texture_tiles,   // compatible with the CPU ones      
				  texture_img,                   // texture_img,     // null or [3][] (RGB) or [4][] RGBA
				  texture_woi_pix,               // texture_woi_pix,     // null or generated texture location/size
				  null,                          //  iclt_fimg,       // will return quad images or null to skip, use quadCLT.linearStackToColor 
					// new parameters, will replace some other?
//				  clt_parameters.getFatZero(isMonochrome()),      // final double              gpu_fat_zero,    // clt_parameters.getGpuFatZero(is_mono);absolute == 30.0\
				  clt_parameters.gpu_corr_scale, //  gpu_corr_scale,  //  0.75; // reduce GPU-generated correlation values
				  clt_parameters.getGpuFatZero(isMonochrome()), // final double              gpu_fat_zero,    // clt_parameters.getGpuFatZero(is_mono);absolute == 30.0\
				  clt_parameters.gpu_sigma_r,    // 0.9, 1.1
				  clt_parameters.gpu_sigma_b,    // 0.9, 1.1
				  clt_parameters.gpu_sigma_g,    // 0.6, 0.7
				  clt_parameters.gpu_sigma_m,    //  =       0.4; // 0.7;
				  (isMonochrome()? 1.0 : clt_parameters.gpu_sigma_rb_corr), // final double              gpu_sigma_rb_corr, //  = 0.5; // apply LPF after accumulating R and B correlation before G, monochrome ? 1.0 : gpu_sigma_rb_corr;
				  clt_parameters.gpu_sigma_corr, //  =    0.9;gpu_sigma_corr_m
				  image_dtt.transform_size - 1, // clt_parameters.gpu_corr_rad,   // = transform_size - 1 ?
				  clt_parameters.corr_red,       // +used
				  clt_parameters.corr_blue,      // +used
				  clt_parameters.max_corr_radius,// 3.9;
				  clt_parameters.corr_mode,      // Correlation mode: 0 - integer max, 1 - center of mass, 2 - polynomial
				  clt_parameters.min_shot,       // 10.0;  // Do not adjust for shot noise if lower than
				  clt_parameters.scale_shot,     // 3.0;   // scale when dividing by sqrt ( <0 - disable correction)
				  clt_parameters.diff_sigma,     // 5.0;//RMS difference from average to reduce weights (~ 1.0 - 1/255 full scale image)
				  clt_parameters.diff_threshold, // 5.0;   // RMS difference from average to discard channel (~ 1.0 - 1/255 full scale image)
				  clt_parameters.diff_gauss,     // true;  // when averaging images, use gaussian around average as weight (false - sharp all/nothing)
				  clt_parameters.min_agree,      // 3.0;   // minimal number of channels to agree on a point (real number to work with fuzzy averages)
				  clt_parameters.dust_remove,    // Do not reduce average weight when only one image differes much from the average
				  geometryCorrection,            // final GeometryCorrection  geometryCorrection,
				  null,                          // final GeometryCorrection  geometryCorrection_main, // if not null correct this camera (aux) to the coordinates of the main
				  clt_parameters.clt_window,
				  disparity_corr,                // final double              disparity_corr, // disparity at infinity
				  clt_parameters.tileX,          // final int               debug_tileX,
				  clt_parameters.tileY,          // final int               debug_tileY,
				  threadsMax,
				  debugLevel);
		  if ((texture_img[0] == null) || (texture_img[1]== null)) {
			  return null;
		  }
		  return texture_img;
	  }
	  
	  public double [][][][] getTextureTilesGPU( // returns texture tiles, as with CPU
			  CLTParameters       clt_parameters,
			  int [][]            tile_op,
			  double [][]         disparity_array, // [tilesY][tilesX]
//			  boolean []          selection,       // null - all
			  final int           threadsMax,      // maximal number of threads to launch
			  final boolean       updateStatus,
			  final int           debugLevel)
	  {
			  
		  final int tilesX = tp.getTilesX(); // may be different from last GPU run !
		  final int tilesY = tp.getTilesY(); // may be different from last GPU run !
		  /*
		  int d = ImageDtt.setImgMask(0, 0xf);
		  d =     ImageDtt.setForcedDisparity(d,true);
		  int [][] tile_op = new int [tilesY][tilesX];
		  int indx = 0;
		  if (selection != null) {
			  for (int ty = 0; ty < tilesY; ty++) {
				  for (int tx = 0; tx < tilesX; tx++) {
					  if (selection[indx++]) {
						  tile_op[ty][tx] = d;
					  }
				  }
			  }
		  } else {
			  for (int ty = 0; ty < tilesY; ty++) {
				  for (int tx = 0; tx < tilesX; tx++) {
					  tile_op[ty][tx] = d;
				  }
			  }
		  }
         */
		  /*
		  double [][] shiftXY = new double [4][2];
		  // not used
		  if (!clt_parameters.fine_corr_ignore) {
			  double [][] shiftXY0 = {
					  {clt_parameters.fine_corr_x_0,clt_parameters.fine_corr_y_0},
					  {clt_parameters.fine_corr_x_1,clt_parameters.fine_corr_y_1},
					  {clt_parameters.fine_corr_x_2,clt_parameters.fine_corr_y_2},
					  {clt_parameters.fine_corr_x_3,clt_parameters.fine_corr_y_3}};
			  shiftXY = shiftXY0;
		  }
          */
		  ImageDtt image_dtt = new ImageDtt(
				  getNumSensors(),
				  clt_parameters.transform_size,
				  clt_parameters.img_dtt,
				  isAux(),
				  isMonochrome(),
				  isLwir(),
				  clt_parameters.getScaleStrength(isAux()),
				  gpuQuad);
		  image_dtt.getCorrelation2d(); // initiate image_dtt.correlation2d, needed if disparity_map != null  
		  double z_correction =  clt_parameters.z_correction;
		  if (clt_parameters.z_corr_map.containsKey(image_name)){ // not used in lwir
			  z_correction +=clt_parameters.z_corr_map.get(image_name);
		  }
		  final double disparity_corr = (z_correction == 0) ? 0.0 : geometryCorrection.getDisparityFromZ(1.0/z_correction);
		  double [][][][] texture_tiles = new double [tilesY][tilesX][][];
		  image_dtt.clt_aberrations_quad_corr_GPU( // USED in LWIR
				  clt_parameters.img_dtt,        // final ImageDttParameters  imgdtt_params,   // Now just extra correlation parameters, later will include, most others
				  1,                             // final int  macro_scale, // to correlate tile data instead of the pixel data: 1 - pixels, 8 - tiles
				  tile_op,                       // per-tile operation bit codes
				  disparity_array,               // clt_parameters.disparity,     // final double            disparity,
				  null,							 // final float  [][][][]     corr_td,         // [tilesY][tilesX][pair][4*64] transform domain representation of 6 corr pairs
				  null, 						 //	final float  [][][][]     corr_combo_td,   // [4][tilesY][tilesX][pair][4*64] TD of combo corrs: qud, cross, hor,vert
					                             // each of the top elements may be null to skip particular combo type
				  null,                          // final float  [][][][]     fdisp_dist,      // [tilesY][tilesX][cams][4], // disparity derivatives vectors or null
				  null,                          //	final float  [][][][]     fpxpy,           // [tilesY][tilesX][cams][2], tile {pX,pY}
				  //// Uses quadCLT from gpuQuad			
				                                 // correlation results - final and partial
				  null,                          // [type][tilesY][tilesX][(2*transform_size-1)*(2*transform_size-1)] // if null - will not calculate
				                                 // [type][tilesY][tilesX] should be set by caller
				                                 // types: 0 - selected correlation (product+offset), 1 - sum
				  null,                          // clt_corr_partial,// [tilesY][tilesX][quad]color][(2*transform_size-1)*(2*transform_size-1)] // if null - will not calculate
				                                 // [tilesY][tilesX] should be set by caller
				                                 // When clt_mismatch is non-zero, no far objects extraction will be attempted
				  null,                          // clt_mismatch,    // [12][tilesY * tilesX] // ***** transpose unapplied ***** ?. null - do not calculate
				                                 // values in the "main" directions have disparity (*_CM) subtracted, in the perpendicular - as is
				  null,                          // disparity_map,   // [8][tilesY][tilesX], only [6][] is needed on input or null - do not calculate
				                                 // last 2 - contrast, avg/ "geometric average)
				  0,                             // disparity_modes, // bit mask of disparity_map slices to calculate/return
				  texture_tiles,                 // 	final double [][][][]     texture_tiles,   // compatible with the CPU ones      
				  null,                          // texture_img,     // null or [3][] (RGB) or [4][] RGBA
				  null,                          // texture_woi,     // null or generated texture location/size
				  null,                          //  iclt_fimg,       // will return quad images or null to skip, use quadCLT.linearStackToColor 
					// new parameters, will replace some other?
//				  clt_parameters.getFatZero(isMonochrome()),      // final double              gpu_fat_zero,    // clt_parameters.getGpuFatZero(is_mono);absolute == 30.0\
				  clt_parameters.gpu_corr_scale, //  gpu_corr_scale,  //  0.75; // reduce GPU-generated correlation values
				  clt_parameters.getGpuFatZero(isMonochrome()), // final double              gpu_fat_zero,    // clt_parameters.getGpuFatZero(is_mono);absolute == 30.0\
				  clt_parameters.gpu_sigma_r,    // 0.9, 1.1
				  clt_parameters.gpu_sigma_b,    // 0.9, 1.1
				  clt_parameters.gpu_sigma_g,    // 0.6, 0.7
				  clt_parameters.gpu_sigma_m,    //  =       0.4; // 0.7;
				  (isMonochrome()? 1.0 : clt_parameters.gpu_sigma_rb_corr), // final double              gpu_sigma_rb_corr, //  = 0.5; // apply LPF after accumulating R and B correlation before G, monochrome ? 1.0 : gpu_sigma_rb_corr;
				  clt_parameters.gpu_sigma_corr, //  =    0.9;gpu_sigma_corr_m
				  image_dtt.transform_size - 1, // clt_parameters.gpu_corr_rad,   // = transform_size - 1 ?
				  clt_parameters.corr_red,       // +used
				  clt_parameters.corr_blue,      // +used
				  clt_parameters.max_corr_radius,// 3.9;
				  clt_parameters.corr_mode,      // Correlation mode: 0 - integer max, 1 - center of mass, 2 - polynomial
				  clt_parameters.min_shot,       // 10.0;  // Do not adjust for shot noise if lower than
				  clt_parameters.scale_shot,     // 3.0;   // scale when dividing by sqrt ( <0 - disable correction)
				  clt_parameters.diff_sigma,     // 5.0;//RMS difference from average to reduce weights (~ 1.0 - 1/255 full scale image)
				  clt_parameters.diff_threshold, // 5.0;   // RMS difference from average to discard channel (~ 1.0 - 1/255 full scale image)
				  clt_parameters.diff_gauss,     // true;  // when averaging images, use gaussian around average as weight (false - sharp all/nothing)
				  clt_parameters.min_agree,      // 3.0;   // minimal number of channels to agree on a point (real number to work with fuzzy averages)
				  clt_parameters.dust_remove,    // Do not reduce average weight when only one image differes much from the average
				  geometryCorrection,            // final GeometryCorrection  geometryCorrection,
				  null,                          // final GeometryCorrection  geometryCorrection_main, // if not null correct this camera (aux) to the coordinates of the main
				  clt_parameters.clt_window,
				  disparity_corr,                // final double              disparity_corr, // disparity at infinity
				  clt_parameters.tileX,          // final int               debug_tileX,
				  clt_parameters.tileY,          // final int               debug_tileY,
				  threadsMax,
				  debugLevel);
		  
		  return texture_tiles;
	  }

	  
//	  public CLTPass3d  CLTMeasureTextures( // perform single pass according to prepared tiles operations and disparity // not used in lwir
	  public void CLTMeasureTextures( // perform single pass according to prepared tiles operations and disparity // not used in lwir
			  final CLTParameters           clt_parameters,
			  final CLTPass3d          scan,
			  final int         threadsMax,  // maximal number of threads to launch
			  final boolean     updateStatus,
			  final int         debugLevel)
	  {
		  if ((gpuQuad == null) || !(isAux()?clt_parameters.gpu_use_aux : clt_parameters.gpu_use_main)) {
			  super.CLTMeasureTextures( // measure background // USED in lwir
					  clt_parameters,
					  scan,
					  threadsMax,  // maximal number of threads to launch
					  updateStatus,
					  debugLevel);
			  return;
		  }
		  final int tilesX = tp.getTilesX();
		  final int tilesY = tp.getTilesY();

//		  CLTPass3d scan = tp.clt_3d_passes.get(scanIndex);
		  // maybe just rely on tile_op ?
		  boolean [] text_sel = new boolean [tilesX * tilesY];
		  for (int ty = 0; ty < tilesY; ty++) {
			  for (int tx = 0; tx < tilesX; tx++) {
				  text_sel[ty*tilesX + tx] =  scan.tile_op[ty][tx] != 0;
			  }
		  }
//		  scan.setTextureSelection(scan.selected); // .clone(); GPU will measure when the texture will be requested
		  scan.setTextureSelection(text_sel); // .clone(); GPU will measure when the texture will be requested
//		  scan.is_measured =   true; // but no disparity map/textures
		  scan.is_combo =      false;
		  if (debugLevel > 10) {
			  tp.showScan(
					  scan,   // CLTPass3d   scan,
					  "CLTMeasureTextures->");
		  }
//		  return scan; // what about processed and other attributes?
		  
	  }
	
	  public void CLTMeasureLY( // perform single pass according to prepared tiles operations and disparity // USED in lwir
			  final CLTParameters clt_parameters,
			  final CLTPass3d     scan,
			  final int           bgIndex, // combine, if >=0
			  final int           threadsMax,  // maximal number of threads to launch
			  final boolean       updateStatus,
			  int           debugLevel)
	  {
		  final int dbg_x = -295-debugLevel;
		  final int dbg_y = -160-debugLevel;
		  final int tilesX = tp.getTilesX();
		  final int tilesY = tp.getTilesY();
		  final int cluster_size =clt_parameters.tileStep;
		  final int clustersX= (tilesX + cluster_size - 1) / cluster_size;
		  final int clustersY= (tilesY + cluster_size - 1) / cluster_size;

//		  CLTPass3d scan = tp.clt_3d_passes.get(scanIndex);
		  scan.setLazyEyeClusterSize(cluster_size);
		  boolean [] force_disparity= new boolean[clustersX * clustersY];
		  //		  scan.setLazyEyeForceDisparity(force_disparity);
		  if (bgIndex >= 0) {
			  CLTPass3d    bg_scan = tp.clt_3d_passes.get(bgIndex);
			  // if at least one tile in a cluster is BG, use BG for the whole cluster and set lazy_eye_force_disparity
			  for (int cY = 0; cY < clustersY; cY ++) {
				  for (int cX = 0; cX < clustersX; cX ++) {
					  boolean has_bg = false;
					  for (int cty = 0; (cty < cluster_size) && !has_bg; cty++) {
						  int ty = cY * cluster_size + cty;
						  if (ty < tilesY) for (int ctx = 0; ctx < cluster_size; ctx++) {
							  int tx = cX * cluster_size + ctx;
							  if ((tx < tilesX ) && (bg_scan.tile_op[ty][tx] > 0)) {
								  has_bg = true;
								  break;
							  }
						  }
					  }
					  if (has_bg) {
						  for (int cty = 0; cty < cluster_size; cty++) {
							  int ty = cY * cluster_size + cty;
							  if (ty < tilesY) for (int ctx = 0; ctx < cluster_size; ctx++) {
								  int tx = cX * cluster_size + ctx;
								  if (tx < tilesX ) {
									  scan.tile_op[ty][tx] =   bg_scan.tile_op[ty][tx];
									  scan.disparity[ty][tx] = bg_scan.disparity[ty][tx];
								  }
							  }
						  }
						  force_disparity[cY * clustersX + cX] = true;
					  }
				  }
			  }
			  scan.setLazyEyeForceDisparity(force_disparity);
		  }
		  int [][]     tile_op =         scan.tile_op;
		  double [][]  disparity_array = scan.disparity;

		  // Should not happen !
		  if (scan.disparity == null) { // not used in lwir
			  System.out.println ("** BUG: should not happen - scan.disparity == null ! **");
			  System.out.println ("Trying to recover");
			  double [] backup_disparity = scan.getDisparity(0);
			  if (backup_disparity == null) {
				  System.out.println ("** BUG: no disparity at all !");
				  backup_disparity = new double[tilesX*tilesY];
			  }
			  scan.disparity = new double[tilesY][tilesX];
			  for (int ty = 0; ty < tilesY; ty++) {
				  for (int tx = 0; tx < tilesX; tx++) {
					  scan.disparity[ty][tx] = backup_disparity[ty*tilesX + tx];
					  if (Double.isNaN(scan.disparity[ty][tx])) {
						  scan.disparity[ty][tx] = 0;
						  tile_op[ty][tx] = 0;
					  }
				  }
			  }
			  disparity_array = scan.disparity;
		  }

		  if (debugLevel > -1){
			  int numTiles = 0;
			  for (int ty = 0; ty < tile_op.length; ty ++) for (int tx = 0; tx < tile_op[ty].length; tx ++){
				  if (tile_op[ty][tx] != 0) numTiles ++;
			  }
//			  System.out.println("CLTMeasure("+scanIndex+"): numTiles = "+numTiles);
			  System.out.println("CLTMeasure(): numTiles = "+numTiles);
			  if ((dbg_y >= 0) && (dbg_x >= 0) && (tile_op[dbg_y][dbg_x] != 0)){
//				  System.out.println("CLTMeasure("+scanIndex+"): tile_op["+dbg_y+"]["+dbg_x+"] = "+tile_op[dbg_y][dbg_x]);
				  System.out.println("CLTMeasure(): tile_op["+dbg_y+"]["+dbg_x+"] = "+tile_op[dbg_y][dbg_x]);
			  }
		  }
		  double min_corr_selected = clt_parameters.min_corr;
		  double [][] shiftXY = new double [getNumSensors()][2];
		  if (!clt_parameters.fine_corr_ignore) {// invalid for AUX!
			  double [][] shiftXY0 = {
					  {clt_parameters.fine_corr_x_0,clt_parameters.fine_corr_y_0},
					  {clt_parameters.fine_corr_x_1,clt_parameters.fine_corr_y_1},
					  {clt_parameters.fine_corr_x_2,clt_parameters.fine_corr_y_2},
					  {clt_parameters.fine_corr_x_3,clt_parameters.fine_corr_y_3}};
			  // FIXME - only first 4 sensors have correction. And is it the same for aux and main? 
			  for (int i = 0; i < shiftXY0.length;i++) {
				  shiftXY[i] = shiftXY0[i];
			  }
		  }

		  ImageDtt image_dtt = new ImageDtt(
				  getNumSensors(),
				  clt_parameters.transform_size,
				  clt_parameters.img_dtt,
				  isAux(),
				  isMonochrome(),
				  isLwir(),
				  clt_parameters.getScaleStrength(isAux()),
				  gpuQuad);
		  image_dtt.getCorrelation2d(); // initiate image_dtt.correlation2d, needed if disparity_map != null  
		  
		  double z_correction =  clt_parameters.z_correction;
		  if (clt_parameters.z_corr_map.containsKey(image_name)){ // not used in lwir
			  z_correction +=clt_parameters.z_corr_map.get(image_name);
		  }
		  final double disparity_corr = (z_correction == 0) ? 0.0 : geometryCorrection.getDisparityFromZ(1.0/z_correction);

		  if (debugLevel > -3){ // now -3
			  tp.showScan(
					  scan,   // CLTPass3d   scan,
					  "LY-combo_scan-"+scan+"_post"); //String title)
		  }
		  // use new, LMA-based mismatch calculation
		  double [][] lazy_eye_data = null;
		  if ((gpuQuad == null) || !(isAux()?clt_parameters.gpu_use_aux_adjust : clt_parameters.gpu_use_main_adjust)) { // CPU
			  lazy_eye_data = image_dtt.cltMeasureLazyEye ( // returns d,s lazy eye parameters 
					  clt_parameters.img_dtt,       // final ImageDttParameters  imgdtt_params,   // Now just extra correlation parameters, later will include, most others
					  tile_op,                      // final int [][]            tile_op,         // [tilesY][tilesX] - what to do - 0 - nothing for this tile
					  disparity_array,              // final double [][]         disparity_array, // [tilesY][tilesX] - individual per-tile expected disparity
					  image_data,                   // final double [][][]      imade_data, // first index - number of image in a quad
					  saturation_imp,               // final boolean [][]        saturation_imp, // (near) saturated pixels or null
					  tilesX * image_dtt.transform_size, // 	final int                 width,
					  clt_parameters.getFatZero(isMonochrome()),      // final double              corr_fat_zero,    // add to denominator to modify phase correlation (same units as data1, data2). <0 - pure sum
					  clt_parameters.corr_red,      // final double              corr_red,
					  clt_parameters.corr_blue,     // final double              corr_blue,
					  clt_parameters.getCorrSigma(image_dtt.isMonochrome()), // final double              corr_sigma,
					  min_corr_selected, // 0.0001; //final double              min_corr,        // 0.02; // minimal correlation value to consider valid
					  geometryCorrection,            // final GeometryCorrection  geometryCorrection,
					  null,                          // final GeometryCorrection  geometryCorrection_main, // if not null correct this camera (aux) to the coordinates of the main
					  clt_kernels,                   // final double [][][][][][] clt_kernels, // [channel_in_quad][color][tileY][tileX][band][pixel] , size should match image (have 1 tile around)
//					  clt_parameters.kernel_step,    // final int                 kernel_step,
					  clt_parameters.clt_window,     // final int                 window_type,
					  shiftXY,                       // final double [][]         shiftXY, // [port]{shiftX,shiftY}
					  disparity_corr,                // final double              disparity_corr, // disparity at infinity
					  clt_parameters.shift_x,        // final double              shiftX, // shift image horizontally (positive - right) - just for testing
					  clt_parameters.shift_y,        // final double              shiftY, // shift image vertically (positive - down)
					  clt_parameters.tileStep,       // final int                 tileStep, // process tileStep x tileStep cluster of tiles when adjusting lazy eye parameters
					  clt_parameters.img_dtt.getMcorrSelLY(getNumSensors()), //    final int                 mcorr_sel, // +1 - all, +2 - dia, +4 - sq, +8 - neibs, +16 - hor + 32 - vert

					  clt_parameters.img_dtt.mcorr_comb_width,					// final int                 mcorr_comb_width,  // combined correlation tile width
					  clt_parameters.img_dtt.mcorr_comb_height,					// final int                 mcorr_comb_height, // combined correlation tile full height
					  clt_parameters.img_dtt.mcorr_comb_offset,					// final int                 mcorr_comb_offset, // combined correlation tile height offset: 0 - centered (-height/2 to height/2), height/2 - only positive (0 to height)
					  clt_parameters.img_dtt.mcorr_comb_disp,					// final double              mcorr_comb_disp,   // Combined tile per-pixel disparity for baseline == side of a square
					  
					  clt_parameters.tileX,        // final int                 debug_tileX,
					  clt_parameters.tileY,         // final int                 debug_tileY,
					  threadsMax, // final int                 threadsMax,  // maximal number of threads to launch
					  debugLevel - 2); // final int                 globalDebugLevel)
			  
		  } else {// GPU
			  // First - measure and get 	fcorr_td, fdisp_dist, fpxpy
			  float  [][][][]     fcorr_td =   new float[tilesY][tilesX][][];
			  float  [][][][]     fdisp_dist = new float[tilesY][tilesX][][];
			  float  [][][][]     fpxpy =      new float[tilesY][tilesX][][];
			  final double        gpu_fat_zero = clt_parameters.getGpuFatZero(isMonochrome()); 
			  
			  image_dtt.clt_aberrations_quad_corr_GPU( // measures 2D correlations
					  clt_parameters.img_dtt,        // final ImageDttParameters  imgdtt_params,   // Now just extra correlation parameters, later will include, most others
					  1,                             // final int  macro_scale, // to correlate tile data instead of the pixel data: 1 - pixels, 8 - tiles
					  tile_op,                       // per-tile operation bit codes
					  disparity_array,               // clt_parameters.disparity,     // final double            disparity,
					  fcorr_td,				 		 // final float  [][][][]     corr_td,         // [tilesY][tilesX][pair][4*64] transform domain representation of 6 corr pairs
					  null,			  				 //	final float  [][][][]     corr_combo_td,   // [4][tilesY][tilesX][pair][4*64] TD of combo corrs: qud, cross, hor,vert
					  // each of the top elements may be null to skip particular combo type
					  fdisp_dist,                    // final float  [][][][]     fdisp_dist,      // [tilesY][tilesX][cams][4], // disparity derivatives vectors or null
					  fpxpy,	                     //	final float  [][][][]     fpxpy,           // [tilesY][tilesX][cams][2], tile {pX,pY}
					  //// Uses quadCLT from gpuQuad			
					  // correlation results - final and partial
					  null,                          // [type][tilesY][tilesX][(2*transform_size-1)*(2*transform_size-1)] // if null - will not calculate
					  // [type][tilesY][tilesX] should be set by caller
					  // types: 0 - selected correlation (product+offset), 1 - sum
					  null,     			         // clt_corr_partial,// [tilesY][tilesX][quad]color][(2*transform_size-1)*(2*transform_size-1)] // if null - will not calculate
					  // [tilesY][tilesX] should be set by caller
					  // When clt_mismatch is non-zero, no far objects extraction will be attempted
					  null,                          // clt_mismatch,    // [12][tilesY * tilesX] // ***** transpose unapplied ***** ?. null - do not calculate
					  // values in the "main" directions have disparity (*_CM) subtracted, in the perpendicular - as is
					  null,              			 // disparity_map,   // [8][tilesY][tilesX], only [6][] is needed on input or null - do not calculate
					  // last 2 - contrast, avg/ "geometric average)
					  0, // disparity_modes,         // disparity_modes, // bit mask of disparity_map slices to calculate/return
					  null,                          // 	final double [][][][]     texture_tiles,   // compatible with the CPU ones      
					  null, // texture_img,          // texture_img,     // null or [3][] (RGB) or [4][] RGBA
					  null, // texture_woi,          // texture_woi,     // null or generated texture location/size
					  null,                          //  iclt_fimg,       // will return quad images or null to skip, use quadCLT.linearStackToColor 
					  clt_parameters.gpu_corr_scale, //.gpu_corr_scale, //  gpu_corr_scale,  //  0.75; // reduce GPU-generated correlation values
					  gpu_fat_zero, // final double              gpu_fat_zero,    // clt_parameters.getGpuFatZero(is_mono);absolute == 30.0\
					  clt_parameters.gpu_sigma_r,    // 0.9, 1.1
					  clt_parameters.gpu_sigma_b,    // 0.9, 1.1
					  clt_parameters.gpu_sigma_g,    // 0.6, 0.7
					  clt_parameters.gpu_sigma_m,    //  =       0.4; // 0.7;
					  (isMonochrome()? 1.0 : clt_parameters.gpu_sigma_rb_corr), // final double              gpu_sigma_rb_corr, //  = 0.5; // apply LPF after accumulating R and B correlation before G, monochrome ? 1.0 : gpu_sigma_rb_corr;
					  clt_parameters.gpu_sigma_corr, //  =    0.9;gpu_sigma_corr_m
					  image_dtt.transform_size - 1, // clt_parameters.gpu_corr_rad,   // = transform_size - 1 ?
					  clt_parameters.corr_red,       // +used
					  clt_parameters.corr_blue,      // +used
					  clt_parameters.max_corr_radius,// 3.9;
					  clt_parameters.corr_mode,      // Correlation mode: 0 - integer max, 1 - center of mass, 2 - polynomial
					  clt_parameters.min_shot,       // 10.0;  // Do not adjust for shot noise if lower than
					  clt_parameters.scale_shot,     // 3.0;   // scale when dividing by sqrt ( <0 - disable correction)
					  clt_parameters.diff_sigma,     // 5.0;//RMS difference from average to reduce weights (~ 1.0 - 1/255 full scale image)
					  clt_parameters.diff_threshold, // 5.0;   // RMS difference from average to discard channel (~ 1.0 - 1/255 full scale image)
					  clt_parameters.diff_gauss,     // true;  // when averaging images, use gaussian around average as weight (false - sharp all/nothing)
					  clt_parameters.min_agree,      // 3.0;   // minimal number of channels to agree on a point (real number to work with fuzzy averages)
					  clt_parameters.dust_remove,    // Do not reduce average weight when only one image differes much from the average
					  geometryCorrection,            // final GeometryCorrection  geometryCorrection,
					  null,                          // final GeometryCorrection  geometryCorrection_main, // if not null correct this camera (aux) to the coordinates of the main
					  clt_parameters.clt_window,
					  disparity_corr,                // final double              disparity_corr, // disparity at infinity
					  clt_parameters.tileX,          // final int               debug_tileX,
					  clt_parameters.tileY,          // final int               debug_tileY,
					  threadsMax,
					  debugLevel);
if (debugLevel < -100) {
			  double [][] lazy_eye_data_dbg = image_dtt.cltMeasureLazyEyeGPU ( // processes 2D correlations, returns d,s lazy eye parameters
					  clt_parameters.img_dtt,        // final ImageDttParameters  imgdtt_params,   // Now just extra correlation parameters, later will include, most others
					  disparity_array,               // final double [][]         disparity_array, // [tilesY][tilesX] - individual per-tile expected disparity
					  fcorr_td,                      // final float  [][][][]     fcorr_td,        // [tilesY][tilesX][pair][4*64] transform domain representation of 6 corr
					  fdisp_dist,                    // final float  [][][][]     fdisp_dist,      // [tilesY][tilesX][cams][4], // disparity derivatives vectors or null
					  fpxpy,                         // 						  final float  [][][][]     fpxpy,           // [tilesY][tilesX][cams][2], tile {pX,pY}
					  clt_parameters.gpu_corr_scale, // gpu_corr_scale,  //  final double              gpu_corr_scale,  //  0.75; // reduce GPU-generated correlation values
					  gpu_fat_zero,                  // final double              gpu_fat_zero,    // clt_parameters.getGpuFatZero(is_mono);absolute == 30.0\
					  geometryCorrection,            // final GeometryCorrection  geometryCorrection,
					  null,                          // geometryCorrection_main, // final GeometryCorrection  geometryCorrection_main, // if not null correct this camera (aux) to the coordinates of the main
					  disparity_corr,                // final double              disparity_corr, // disparity at infinity
					  1, // clt_parameters.tileStep,       // final int                 tileStep, // process tileStep x tileStep cluster of tiles when adjusting lazy eye parameters
					  0, // clt_parameters.super_radius,   // final int                 super_radius, // 0 - none, 1 - 3x3, 2 - 5x5,

					  clt_parameters.tileX,          // final int               debug_tileX,
					  clt_parameters.tileY,          // final int               debug_tileY,
					  threadsMax,
					  debugLevel);
			  ExtrinsicAdjustment ea_dbg = new ExtrinsicAdjustment (
					  geometryCorrection, // GeometryCorrection gc,
					  1, // int         clusterSize,
					  tilesX, // int         clustersX,
					  tilesY); // int         clustersY)
			  
//			  	double [][] dbg_img = new double[ExtrinsicAdjustment.DATA_TITLES.length][lazy_eye_data_dbg.length];
			  	double [][] dbg_img = new double[ea_dbg.data_titles.length][lazy_eye_data_dbg.length];
			  	for (int i = 0; i < dbg_img.length; i++) {
			  		Arrays.fill(dbg_img[i], Double.NaN);
			  	}
			  	for (int nt = 0; nt < lazy_eye_data_dbg.length; nt++) if (lazy_eye_data_dbg[nt] != null){
				  	for (int i = 0; i < dbg_img.length; i++) {
				  		dbg_img[i][nt] = lazy_eye_data_dbg[nt][i];
				  	}
			  	}
			  
				ShowDoubleFloatArrays.showArrays(
						dbg_img,
						tp.getTilesX(),
						tp.getTilesY(),
						true,
						"LY-everytile",
						ea_dbg.data_titles); //	ExtrinsicAdjustment.DATA_TITLES);
			  
}			  
			  
			  lazy_eye_data = image_dtt.cltMeasureLazyEyeGPU ( // processes 2D correlations, returns d,s lazy eye parameters
					  clt_parameters.img_dtt,        // final ImageDttParameters  imgdtt_params,   // Now just extra correlation parameters, later will include, most others
					  disparity_array,               // final double [][]         disparity_array, // [tilesY][tilesX] - individual per-tile expected disparity
					  fcorr_td,                      // final float  [][][][]     fcorr_td,        // [tilesY][tilesX][pair][4*64] transform domain representation of 6 corr
					  fdisp_dist,                    // final float  [][][][]     fdisp_dist,      // [tilesY][tilesX][cams][4], // disparity derivatives vectors or null
					  fpxpy,                         // 						  final float  [][][][]     fpxpy,           // [tilesY][tilesX][cams][2], tile {pX,pY}
					  clt_parameters.gpu_corr_scale, // gpu_corr_scale,  //  final double              gpu_corr_scale,  //  0.75; // reduce GPU-generated correlation values
					  gpu_fat_zero,                  // final double              gpu_fat_zero,    // clt_parameters.getGpuFatZero(is_mono);absolute == 30.0\
					  geometryCorrection,            // final GeometryCorrection  geometryCorrection,
					  null,                          // geometryCorrection_main, // final GeometryCorrection  geometryCorrection_main, // if not null correct this camera (aux) to the coordinates of the main
					  disparity_corr,                // final double              disparity_corr, // disparity at infinity
					  clt_parameters.tileStep,       // final int                 tileStep, // process tileStep x tileStep cluster of tiles when adjusting lazy eye parameters
					  clt_parameters.super_radius,   // final int                 super_radius, // 0 - none, 1 - 3x3, 2 - 5x5,

					  clt_parameters.tileX,          // final int               debug_tileX,
					  clt_parameters.tileY,          // final int               debug_tileY,
					  threadsMax,
					  debugLevel);
			  
			  
		  }

		  scan.setLazyEyeData(lazy_eye_data);
		  scan.is_measured =   true; // but no disparity map/textures
		  scan.is_combo =      false;
		  scan.resetProcessed();
//		  return scan;
	  }
	  
	  public void gpuResetCorrVector() {
		  if (getGPU() != null) {
			  getGPU().resetGeometryCorrectionVector();
			  getGPU().resetGeometryCorrection(); // testing, just in case, should be removed
		  }
	  }
	  public void gpuResetGeometryCorrection() {
		  if (getGPU() != null) {
			  getGPU().resetGeometryCorrection();
		  }
	  }
	  
	  public void resetBayer() {
		  if (getGPU() != null) {
			  getGPU().resetBayer();
		  }
	  }
	  
		public QuadCLT spawnQuadCLT(
				String              set_name,
				CLTParameters       clt_parameters,
				ColorProcParameters colorProcParameters, //
				int                 threadsMax,
				int                 debugLevel)
		{
			QuadCLT quadCLT = new QuadCLT(this, set_name); //null
			
			quadCLT.restoreFromModel(
					clt_parameters,
					colorProcParameters,
					null,                 // double []    noise_sigma_level,
					-1,                   // noise_variant, // <0 - no-variants, compatible with old code
					null,                 // final QuadCLTCPU     ref_scene, // may be null if scale_fpn <= 0
					threadsMax,
					debugLevel);
			
			return quadCLT;
		}
		public QuadCLT spawnNoModelQuadCLT(
				String              set_name,
				CLTParameters       clt_parameters,
				ColorProcParameters colorProcParameters, //
				int                 threadsMax,
				int                 debugLevel)
		{
			QuadCLT quadCLT = new QuadCLT(this, set_name); //null
			
			return (QuadCLT) quadCLT.restoreNoModel(
					clt_parameters,
					colorProcParameters,
					null,                 // double []    noise_sigma_level,
					-1,                   // noise_variant, // <0 - no-variants, compatible with old code
					null,                 // final QuadCLTCPU     ref_scene, // may be null if scale_fpn <= 0
					threadsMax,
					debugLevel);
		}

	  
}
