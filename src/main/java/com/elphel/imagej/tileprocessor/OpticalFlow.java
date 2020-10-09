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

import java.util.Arrays;
import java.util.Comparator;
import java.util.concurrent.atomic.AtomicInteger;

import com.elphel.imagej.cameras.CLTParameters;
import com.elphel.imagej.common.DoubleGaussianBlur;
import com.elphel.imagej.common.ShowDoubleFloatArrays;

public class OpticalFlow {
	public static double [] ZERO3 = {0.0,0.0,0.0};
	public static double  LINE_ERR = 0.1;
	public int            threadsMax = 100;  // maximal number of threads to launch
	public boolean        updateStatus = true;

	public OpticalFlow (
			int            threadsMax,  // maximal number of threads to launch
			boolean        updateStatus) {
		this.threadsMax = threadsMax;
		this.updateStatus = updateStatus;
		
	}

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
///		final int fullTileLen =          fullTileSize * fullTileSize;
		double wdiag = 0.25 *diagonal_weight / (diagonal_weight + 1.0);
		double wortho = 0.25 / (diagonal_weight + 1.0);
		final double [] neibw = {wortho, wdiag, wortho, wdiag, wortho, wdiag, wortho, wdiag}; 
		
		final Thread[] threads = ImageDtt.newThreadArray(threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		final int dbg_mtile = 0; // 203;
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				public void run() {
					for (int iMTile = ai.getAndIncrement(); iMTile < nan_tiles.length; iMTile = ai.getAndIncrement()) {
						if (iMTile == dbg_mtile) {
//							System.out.println("iMTile = "+iMTile);
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
		
		if (debug_level > 0) {
			// show debug image
			String title = qthis.getImageName()+"-NO-NaN";
			showMacroTiles(
					title,        // String title,
					nan_tiles, // double [][][] source_tiles,
					qthis,        // final QuadCLT qthis,
					margin);      // final int     margin); // extra margins over 16x16 tiles to accommodate distorted destination tiles
		}
//		System.out.println("fillTilesNans() DONE.");
	}

	/**
	 * Calculatde confidence for the interscene X,Y correlation
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
			// show debug image
			String [] titles ={"dX","dY","Strength"};
			final double [][] dbg_img = new double [titles.length][width * height];
			for (int l = 0; l < dbg_img.length; l++) {
				Arrays.fill(dbg_img[l],  Double.NaN);
			}
			for (int mtile = 0; mtile < flowXYS.length; mtile++) if (flowXYS[mtile] != null){
				for (int l = 0; l < dbg_img.length; l++) {
					dbg_img[l][mtile] = flowXYS[mtile][l];
				}
			}
			(new ShowDoubleFloatArrays()).showArrays(
					dbg_img,
					width,
					height,
					true,
					debug_title,
					titles);
		}
		//		System.out.println("fillTilesNans() DONE.");
		return flowXYS;
	}
	
	
	
	public double [][] correlate2DIterate( // returns optical flow and confidence
			final ImageDttParameters  imgdtt_params,   // Now just extra correlation parameters, later will include, most others
			// for prepareSceneTiles()			
			final double []   scene_xyz,     // camera center in world coordinates
			final double []   scene_atr,     // camera orientation relative to world frame
			final QuadCLT     scene_QuadClt,
			final QuadCLT     reference_QuadClt,
			final double [][][] reference_tiles, // prepared with prepareReferenceTiles() + fillTilesNans();
			// flowXY should be initialized to all pairs of zeros (or deliberate pixel offset pairs if initial error is too high, will be modified with each iteration
			final double [][] flowXY, // per macro tile {mismatch in image pixels in X and Y directions // initialize to [reference_tiles.length][2]
			final double      tolerance_absolute, // absolute disparity half-range in each tile
			final double      tolerance_relative, // relative disparity half-range in each tile
			final double      occupancy,          // fraction of remaining  tiles (<1.0)
			final int         num_passes,
			final double      max_change,
			// for correlate2DSceneToReference ()
			final double []   chn_weights, // absolute, starting from strength (strength,r,b,g)
			final double      corr_sigma,
			final double      fat_zero,
			final boolean     late_normalize,
			// for correlation2DToVectors_CM()
			final int         iradius,      // half-size of the square to process 
			final double      dradius,      // weight calculation (1/(r/dradius)^2 + 1)
			final int         refine_num,   // number of iterations to apply weights around new center
			
			final int         num_run_all, // run all tiles for few iterations before filtering
			final int         max_tries,
			
			// for recalculateFlowXY()
			final double      magic_scale, // 0.85 for CM
			final double      min_change,
			
			final int         debug_level)
	{
		final TileProcessor tp =         reference_QuadClt.getTileProcessor();
		final int transform_size =       tp.getTileSize();
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
					num_passes,               // final int         num_passes,
					max_change,               // final double      max_change,
					-1); //-1); // 1); // 2);                       // final int         debug_level)
			// undefine tiles in flowXY that are never used
			if (ntry == 0) {
				for (int i = 0; i <flowXY.length; i++) {
//					if (flowXY_run[i] == null) {
					if ((scene_tiles[i] == null) || (reference_tiles[i] == null)) {
						flowXY[i] = null;	
					}
				}
			}
			
			//
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
			boolean     ignore_worsening = ntry < num_run_all; // (num_run_all + 10);
			if (debug_level > 0) {
				System.out.println("======== NTRY "+ntry +" ========");
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
					debug_level);                         // final int         debug_level);
			if (flowXY_run == null) { // nothing to do left
				break; 
			}
		}
		return flowXY; // it is also in input arguments
	}
	
	
	
	
	double [][] recalculateFlowXY(
			final double [][] currentFlowXY,
			final double [][] corr_vectorsXY,
			final double      magic_scale) // 0.85 for CM
	{
		
		final Thread[] threads = ImageDtt.newThreadArray(threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		final double [][]   flowXY =  new double [currentFlowXY.length][];
		final int dbg_mtile = 620; // 453; // 500;
		final double rmagic_scale = 1.0/magic_scale;

		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				public void run() {
					for (int iMTile = ai.getAndIncrement(); iMTile < currentFlowXY.length; iMTile = ai.getAndIncrement())
						if ((currentFlowXY[iMTile] != null) && (corr_vectorsXY[iMTile] != null)){
  						if (iMTile == dbg_mtile) {
//							System.out.println("iMTile = "+iMTile);
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
		final int dbg_mtile = 994; // 473; // 295; // 15/7 620; // 453; // 500;
		final double rmagic_scale = 1.0/magic_scale;
		final AtomicInteger aupdate = new AtomicInteger(0); //number of tiles to recalculate
		final double reduce_step = 0.5; //multiply step if calculated difference is larger thart the previous
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				public void run() {
					for (int iMTile = ai.getAndIncrement(); iMTile < flowXY.length; iMTile = ai.getAndIncrement()) {
  						if (iMTile == dbg_mtile) {
//							System.out.println("iMTile = "+iMTile);
						}
  						if (flowXY[iMTile] != null){
  							if (corr_vectorsXY[iMTile] == null) {
  								if (min_change <= 0.0) {
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
  								// apply correction in any case
//  								flowXY[iMTile][0] += dx;
//  								flowXY[iMTile][1] += dy;
  								double new_diff = Math.sqrt(dx*dx + dy*dy);
//  		  						if ((debug_level > 0) && (iMTile == dbg_mtile)) {
//  									System.out.print(String.format("iMTile = %4d (%2d / %2d) flowXY = [%8.6f/%8.6f] dx = %8.6f dy =  %8.6f  abs= %8.6f ",
//  											iMTile,  (iMTile %40), (iMTile / 40), flowXY[iMTile][0], flowXY[iMTile][1], dx,dy,new_diff));
//  								}
  								
  								double last_change = abs_change[iMTile]; // may be NaN;
  								abs_change[iMTile] = new_diff;
  								
  								if ((debug_level >2) && (new_diff > last_change) && (min_change > 0.0)) {
  									System.out.println("iMTile="+iMTile+", new_diff="+ new_diff+", last_change="+last_change);
  								}
//  								if ((new_diff > min_change) || (min_change <= 0.0)) { // not worse (any value are not worse than NaN) and larger than threshold
  								if (new_diff < min_change) {
  									if (flowXY_prev[iMTile] != null) {
  										flowXY_prev[iMTile][0] = flowXY[iMTile][0];
  										flowXY_prev[iMTile][1] = flowXY[iMTile][1];
  									} else {
  										flowXY_prev[iMTile] = flowXY[iMTile].clone();
  									}
  									flowXY[iMTile][0] += dx;
  									flowXY[iMTile][1] += dy;
  								} else {
  									if (ignore_worsening || !(new_diff >= last_change)) { // better or ignore - continue iterations
  										//
  										if ((debug_level > 1) && (iMTile == dbg_mtile))  {
  											System.out.println(String.format("iMTile = %4d (%2d / %2d) flowXY = [%8.6f/%8.6f] step_scale = %8.6f dx = %8.6f dy =  %8.6f  abs= %8.6f previous = %8.6f CONTINUE",
  													iMTile, (iMTile %40), (iMTile / 40), flowXY[iMTile][0], flowXY[iMTile][1], step_scale[iMTile], dx,dy,new_diff, last_change));
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
  									} else if ((new_diff >= last_change) && (min_change > 0)) { // worse - reduce step, but still apply
   										if (debug_level > 2) {
   											System.out.println(String.format("iMTile = %4d (%2d / %2d) flowXY = [%8.6f/%8.6f] step_scale = %8.6f dx = %8.6f dy =  %8.6f  abs= %8.6f previous = %8.6f REDUCED STEP",
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
   										
  	  									step_scale[iMTile] *= reduce_step;
   										dx *= reduce_step;
   										dx *= reduce_step;
   										
   										flowXY[iMTile][0] += dx;
   										flowXY[iMTile][1] += dy;
   										flowXY_task[iMTile] = flowXY[iMTile]; // set to measure
   										abs_change[iMTile] = last_change; // restore previous step
   										aupdate.getAndIncrement();
//   									} else if (debug_level > 1) {
//	  									System.out.println(String.format("iMTile = %4d (%2d / %2d) dx = %8.6f dy =  %8.6f  abs= %8.6f previous = %8.6f",
//  	  											iMTile, (iMTile %40), (iMTile / 40), dx,dy,new_diff, last_change));
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
	
	
	
	
	
	public double [][] correlation2DToVectors_CM(
			final double [][] corr2d_tiles, // per 2d calibration tiles (or nulls)
			final int         transform_size,
			final int         iradius,      // half-size of the square to process 
			final double      dradius,      // weight calculation (1/(r/dradius)^2 + 1)
			final int         refine_num,   // number of iterations to apply weights around new center
			final int         debug_level)
	{
		final Thread[] threads = ImageDtt.newThreadArray(threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		final double [][]   vectors_xys =    new double [corr2d_tiles.length][];
		final int dbg_mtile = 620; // 453; // 500;

		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				public void run() {
					for (int iMTile = ai.getAndIncrement(); iMTile < corr2d_tiles.length; iMTile = ai.getAndIncrement()) if (corr2d_tiles[iMTile] != null) {
  						if (iMTile == dbg_mtile) {
//							System.out.println("iMTile = "+iMTile);
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
	
	public double [] getCorrCenterXYS_CM(
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
		if (iYMmax >= transform_size) {
			iYMmax = transform_size -1;
		}
		if (iYMmax >= transform_size) {
			iXMmax = transform_size -1;
		}
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
	
	
	public double [][] correlate2DSceneToReference(// to match to reference
			final ImageDttParameters  imgdtt_params,   // Now just extra correlation parameters, later will include, most others
			final QuadCLT     scene_QuadClt,
			final QuadCLT     reference_QuadClt,
			final double [][][] scene_tiles,     // prepared with prepareSceneTiles()
			final double [][][] reference_tiles, // prepared with prepareReferenceTiles() + fillTilesNans(); - combine?
			final double [][] flowXY_frac, // X, YH fractional shift [-0.5,0.5) to implement with FD rotations
			final double []   chn_weights, // absolute, starting from strength (strength,r,b,g)
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
		
		final int dbg_mtile = 0; // 203;
		final int num_channels = chn_weights.length;
		final int chn_offset = QuadCLT.DSRBG_STRENGTH; // start scene_tiles, reference tiles with this 2-nd index
		
//		final int corr_size = transform_size * 2 -1;
		final ImageDtt image_dtt = new ImageDtt(
				transform_size,
				false,
				false,
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
							imgdtt_params,              // ImageDttParameters  imgdtt_params,
							transform_size,             // int transform_size,
							2.0,                        //  double wndx_scale, // (wndy scale is always 1.0)
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
//								System.out.println("iMTile = "+iMTile);
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
/*
								corr2d.normalize_TD(
										corr_tiles_TD[iMTile],         // double [][] td,
										filter,          // double []   lpf, // or null
										fat_zero);       // double      fat_zero);

								//						    	double [] corrs_ortho = corr2d.convertCorrToPD(
								corr_tiles[iMTile] = corr2d.convertCorrToPD(
										corr_tiles_TD[iMTile]);     // double [][] td);
*/										
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
//						DttRad2 dtt = new DttRad2(transform_size);
//						dtt.set_window(1);
						final TileNeibs tn =  new TileNeibs(macroTilesX, macroTilesY);

						double [][] corr_tile_2D = new double [4][tile_length];
						Correlation2d corr2d = new Correlation2d(
								imgdtt_params,              // ImageDttParameters  imgdtt_params,
								transform_size,             // int transform_size,
								2.0,                        //  double wndx_scale, // (wndy scale is always 1.0)
								false,                      // boolean monochrome,
								false);                     //   boolean debug)
						// reference tile should not be null, scene = may be
						for (int iMTile = ai.getAndIncrement(); iMTile < macroTiles; iMTile = ai.getAndIncrement()) if (reference_tiles[iMTile] != null){ 
//							if (((scene_tiles[iMTile] != null) && (reference_tiles[iMTile] != null)) || (combine_radius > 0)) {
//							if ((scene_tiles[iMTile] != null) && (reference_tiles[iMTile] != null)) {
							if (true) { // !combine_empty_only  || (corr_tiles_TD[iMTile] == null)) {
								if (iMTile == dbg_mtile) {
//									System.out.println("iMTile = "+iMTile);
								}
								if ((combine_radius > 0) && (!combine_empty_only  || (corr_tiles_TD[iMTile] == null))) { // 
									for (int q = 0; q< 4; q++) {
										Arrays.fill(corr_tile_2D[q], 0.0);
									}
									double disp_tol = tolerance_absolute + tolerance_relative * avg_disparity_ref[iMTile];
									double sw = 0;
//									int num_tiles = 0;
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

								//						    	double [] corrs_ortho = corr2d.convertCorrToPD(
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
//		final int min_remain_center = (int) Math.round(center_occupancy * transform_size * transform_size);
		final Thread[] threads = ImageDtt.newThreadArray(threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
//		final double [][] pXpYD = new double [tilesX * tilesY][]; // to hold scene pX, pY, disparity, where pX, pY include flowXY (image pre-shift)
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
		final int dbg_mtile = 54; // 453; // 500;
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				public void run() {
					//					double [][] pXpYD = new double [fullTileLen][];
					double [][] tXtYD = new double [fullTileLen][]; // measured in tiles, not pixels (disparity - still pixels)
					for (int iMTile = ai.getAndIncrement(); iMTile < reference_tiles.length; iMTile = ai.getAndIncrement())
						if ((reference_tiles[iMTile] != null) && (flowXY[iMTile] != null)){
							if (iMTile == dbg_mtile) {
//								System.out.println("iMTile = "+iMTile);
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
												if (xyd != null) {
													tXtYD[iTile] = new double [] {
															((xyd[0] + flowXY[iMTile][0])/transform_size),
															((xyd[1] + flowXY[iMTile][1])/transform_size),
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
							double [][] scene_slices = new double [dsrbg_scene.length][iwidth*iheight];
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
								String [] dbg_titles= {"d","s","r","b","g"};
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
								String [] dbg_titles= {"d","s","r","b","g"};
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
								String [] dbg_titles= {"d","s","r","b","g"};
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
								String [] dbg_titles= {"d","s","r","b","g"};
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
									iwidth);      // final int         width

							if ((debug_level>1) && (iMTile == dbg_mtile)) {
								String [] dbg_titles= {"d","s","r","b","g"};
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
								String [] dbg_titles= {"d","s","r","b","g"};
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
									String [] dbg_titles= {"d","s","r","b","g"};
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
			String title = scene_QuadClt.getImageName()+"-scene_tiles";
			showMacroTiles(
					title,        // String title,
					scene_tiles, // double [][][] source_tiles,
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
			if (debug_level > 1) {
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
	
	
	// helper to be called from thread
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
//			if (!Double.isNaN(slices[QuadCLT.DSRBG_DISPARITY][i]) &&  (strength[i] > 0.0)) {
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
	
	
/*
	public double [] getImageCoordinatesERS(
			QuadCLT cameraQuadCLT, // camera station that got image to be to be matched 
			double px,                // pixel coordinate X in this camera view
			double py,                // pixel coordinate Y in this camera view
			double disparity,         // this view disparity 
			boolean distortedView,    // This camera view is distorted (diff.rect), false - rectilinear
			double [] reference_xyz,  // this view position in world coordinates (typically zero3)
			double [] reference_atr,  // this view orientation relative to world frame  (typically zero3)
			boolean distortedCamera,  // camera view is distorted (false - rectilinear)
			double [] camera_xyz,     // camera center in world coordinates
			double [] camera_atr,     // camera orientation relative to world frame
			double    line_err)       // threshold error in scan lines (1.0)
	
 */
	
	public double [][][] prepareReferenceTiles(
			final QuadCLT     qthis,
//			final int         margin, // extra margins over 16x16 tiles to accommodate distorted destination tiles
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
	
	public void showMacroTiles(
			String title,
			double [][][] source_tiles,
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
		for (int mtile = 0; mtile < source_tiles.length; mtile++) if (source_tiles[mtile] != null){
			int mTileY = mtile / macroTilesX;
			int mTileX = mtile % macroTilesX;
			for (int iY = 0; iY < fullTileSize; iY++) {
				int tileY = (fullTileSize +1) * mTileY + iY;
				for (int iX = 0; iX < fullTileSize; iX++) {
					int tileX = (fullTileSize +1) * mTileX + iX;
					for (int l = 0; l < dbg_img.length; l++) {
						dbg_img[l][tileY * dbg_with + tileX] = source_tiles[mtile][l][iY * fullTileSize + iX];
					}							
				}
			}
		}
		String [] dsrbg_titles = {"d", "s", "r", "b", "g"};
		(new ShowDoubleFloatArrays()).showArrays(
				dbg_img,
				dbg_with,
				dbg_height,
				true,
				title,
				dsrbg_titles);
		
	}

	public void showMacroTiles(
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
		String [] titles = {"d", "s", "r", "b", "g"};
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
	
	
	
	
	public void showCorrTiles(
			String title,
			double [][][] source_tiles,
			int           tilesX,
			int           tile_width,
			int           tile_height) // extra margins over 16x16 tiles to accommodate distorted destination tiles
	{
		// show debug image
		final int dbg_with =   tilesX * (tile_width +1) - 1;
		final int dbg_height = (source_tiles[0].length / tilesX) * (tile_height +1) - 1;
		final double [][] dbg_img = new double [source_tiles.length][dbg_with * dbg_height];
		for (int l = 0; l < dbg_img.length; l++) {
			Arrays.fill(dbg_img[l],  Double.NaN);
		}
		for (int slice = 0; slice < source_tiles.length; slice++) {
			for (int mtile = 0; mtile < source_tiles[slice].length; mtile++) if (source_tiles[slice][mtile] != null){
				int mTileY = mtile / tilesX;
				int mTileX = mtile % tilesX;
				for (int iY = 0; iY < tile_height; iY++) {
					int tileY = (tile_height +1) * mTileY + iY;
					for (int iX = 0; iX < tile_width; iX++) {
						int tileX = (tile_width +1) * mTileX + iX;
						dbg_img[slice][tileY * dbg_with + tileX] = source_tiles[slice][mtile][iY * tile_width + iX];
					}
				}
			}
		}
//		String [] dsrbg_titles = {"d", "s", "r", "b", "g"};
		(new ShowDoubleFloatArrays()).showArrays(
				dbg_img,
				dbg_with+0,
				dbg_height,
				true,
				title); //	dsrbg_titles);
	}
	
	
	
	public double[][][]  get_pair(
			CLTParameters  clt_parameters,			
			double k_prev,
			QuadCLT reference_QuadCLT,
			QuadCLT scene_QuadCLT,
			double corr_scale, //  = 0.75
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
		
		double [][] dsrbg = transformCameraVew( // shifts previous image correctly (right)
				reference_QuadCLT,       // reference
				scene_QuadCLT,       // QuadCLT   camera_QuadClt,
				camera_xyz0, // double [] camera_xyz, // camera center in world coordinates
				camera_atr0, //double [] camera_atr, // camera orientation relative to world frame
				iscale);
		double [][][] pair = {reference_QuadCLT.getDSRBG(),dsrbg};
		
		// combine this scene with warped previous one
		if (debug_level > 0) {
			String [] rtitles = new String[2* dsrbg_titles.length];
			double [][] dbg_rslt = new double [rtitles.length][];
			for (int i = 0; i < dsrbg_titles.length; i++) {
				rtitles[2*i] =    dsrbg_titles[i]+"0";
				rtitles[2*i+1] =  dsrbg_titles[i];
				dbg_rslt[2*i] =   pair[0][i];
				dbg_rslt[2*i+1] = pair[1][i];
			}
			String title = this_image_name+"-"+scene_QuadCLT.image_name+"-dt"+dt;
			(new ShowDoubleFloatArrays()).showArrays(
					dbg_rslt,
					tilesX,
					tilesY,
					true,
					title,
					rtitles);
		}
		/* */
		double      tolerance_absolute = 0.25; // absolute disparity half-range in each tile
		double      tolerance_relative = 0.2; // relative disparity half-range in each tile
		double      center_occupancy =   0.25;   // fraction of remaining  tiles in the center 8x8 area (<1.0)
		int         num_passes = 100;
		double      max_change = 0.005 ;

		double      tolerance_absolute_inter = 0.25; // absolute disparity half-range in each tile
		double      tolerance_relative_inter = 0.2; // relative disparity half-range in each tile
		double      occupancy_inter =         0.25;   // fraction of remaining  tiles in the center 8x8 area (<1.0)

		// Add limitations on disparity ? To be used with multi-tile consolidation
		// Check with walkingh, not only rotating
		double [][][] reference_tiles = prepareReferenceTiles(
				reference_QuadCLT,        // final QuadCLT     qthis,
				tolerance_absolute, // final double      tolerance_absolute, // absolute disparity half-range in each tile
				tolerance_relative, // final double      tolerance_relative, // relative disparity half-range in each tile
				center_occupancy,   // final double      center_occupancy,   // fraction of remaining  tiles in the center 8x8 area (<1.0)
				-1); // -1); // 2); // final int         debug_level)
		
		fillTilesNans(
				reference_tiles,          // final double [][][] nan_tiles,
				reference_QuadCLT,                 // final QuadCLT     qthis,
				num_passes,            // final int         num_passes,
				max_change,            // final double      max_change,
				-1); //-1); // 2);                    // final int         debug_level)
		
		double [][] flowXY = new double [reference_tiles.length][2]; // zero pre-shifts
		double [][] flowXY_frac = new double [reference_tiles.length][]; // Will contain fractional X/Y shift for CLT
		double []   chn_weights = {1.0,1.0,1.0,1.0}; // strength, r,b,g
//		double []   chn_weights = {1.0,0.0,0.0,0.0}; // strength, r,b,g
//		double []   chn_weights = {0.0,1.0,1.0,1.0}; // strength, r,b,g
		// Apply DOG to colors, normalize by standard deviation?
		double      corr_sigma = 0.5;
		double      fat_zero =   0.05;
		double      frac_radius = 0.9;  // add to integer radius for window calculation
		double      tolerance_absolute_inter_macro = 0.25; // absolute disparity half-range to consolidate macro tiles
		double      tolerance_relative_inter_macro = 0.2;  // relative disparity half-range to consolidate macro tiles
		int         iradius =    3;      // half-size of the square to process 
		double      dradius =    1.5;      // weight calculation (1/(r/dradius)^2 + 1)
		int         refine_num = 5;   // number of iterations to apply weights around new center
		int max_rad = 3;
		
		boolean combine_empty_only = true; // false;
		double magic_scale = 0.85; // 2.0 * 0.85;

		boolean late_normalize_iterate = true;
		
		int          num_run_all = 3; // 5; // 5; // run all tiles for few iterations before filtering
		int          max_tries =   50; // 100;
		
		// for recalculateFlowXY()
		double      min_change = 0.1; // 01;//   sqrt (dx*dx + dy*dy) for correction (int tiles) in pixels
		
//		double []   abs_change = new double [reference_tiles.length]; // updated 
		int         debug_level_iterate = -1; // 2;

		int transform_size = tp.getTileSize();
		int macroTilesX =          tilesX/transform_size;
		int macroTilesY =          tilesY/transform_size;
		
//		public double [][] 
		correlate2DIterate( // returns optical flow and confidence
				clt_parameters.img_dtt, // final ImageDttParameters  imgdtt_params,   // Now just extra correlation parameters, later will include, most others
				// for prepareSceneTiles()			
				camera_xyz0,              // final double []   scene_xyz,     // camera center in world coordinates
				camera_atr0,              // final double []   scene_atr,     // camera orientation relative to world frame
				scene_QuadCLT,                    // final QuadCLT     scene_QuadClt,
				reference_QuadCLT,                    // final QuadCLT     reference_QuadClt,
				reference_tiles,          // final double [][][] reference_tiles, // prepared with prepareReferenceTiles() + fillTilesNans();
				// flowXY should be initialized to all pairs of zeros (or deliberate pixel offset pairs if initial error is too high, will be modified with each iteration
				flowXY,                   // final double [][] flowXY, // per macro tile {mismatch in image pixels in X and Y directions // initialize to [reference_tiles.length][2]
				tolerance_absolute_inter, // final double      tolerance_absolute, // absolute disparity half-range in each tile
				tolerance_relative_inter, // final double      tolerance_relative, // relative disparity half-range in each tile
				occupancy_inter,          // final double      occupancy,          // fraction of remaining  tiles (<1.0)
				num_passes,               // final int         num_passes,
				max_change,               // final double      max_change,
				// for correlate2DSceneToReference ()
				chn_weights,              // final double []   chn_weights, // absolute, starting from strength (strength,r,b,g)
				corr_sigma,               // final double      corr_sigma,
				fat_zero,                 //  final double      fat_zero,
				late_normalize_iterate,   // final boolean     late_normalize,
				// for correlation2DToVectors_CM()
				iradius,                  // final int         iradius,      // half-size of the square to process 
				dradius,                  // final double      dradius,      // weight calculation (1/(r/dradius)^2 + 1)
				refine_num,               // final int         refine_num,   // number of iterations to apply weights around new center
				num_run_all,              // final int         num_run_all, // run all tiles for few iterations before filtering
				max_tries,                // final int         max_tries,
				// for recalculateFlowXY()
//				abs_change,               // final double []   abs_change, // updated 
				magic_scale,              // final double      magic_scale, // 0.85 for CM
				min_change,               // final double      min_change,
				debug_level_iterate);     // final int         debug_level)
		
		if (debug_level > 0) {
			String dbg_title = "dXdY-"+scene_QuadCLT.getImageName()+"-"+reference_QuadCLT.getImageName();
			String [] dbg_titles = {"dpX", "dpY"};
			double [][] dbg_img = new double [dbg_titles.length][macroTilesX*macroTilesY];
				Arrays.fill(dbg_img[0], Double.NaN);
				Arrays.fill(dbg_img[1], Double.NaN);
				for (int i = 0; i < flowXY.length; i++) if (flowXY[i] != null){
					dbg_img[0][i] = flowXY[i][0];
					dbg_img[1][i] = flowXY[i][1];
				}
			(new ShowDoubleFloatArrays()).showArrays(
					dbg_img,
					macroTilesX,
					macroTilesY,
					true,
					dbg_title,
					dbg_titles);
		}
		
		
		
		double [][][] scene_tiles = prepareSceneTiles(// to match to reference
				// null for {scene,reference}{xyz,atr} uses instances globals 
				camera_xyz0,              // final double []   scene_xyz,     // camera center in world coordinates
				camera_atr0,              // final double []   scene_atr,     // camera orientation relative to world frame
				scene_QuadCLT,                    // final QuadCLT     scene_QuadClt,
				reference_QuadCLT,                    // final QuadCLT     reference_QuadClt,
				reference_tiles,          // final double [][][] reference_tiles, // prepared with prepareReferenceTiles() + fillTilesNans();
				flowXY,                   // final double [][] flowXY, // per macro tile {mismatch in image pixels in X and Y directions
				flowXY_frac,              // final double [][] flowXY_frac, // should be initialized as [number of macro tiles][] - returns fractional shifts [-0.5, 0.5)
				tolerance_absolute_inter, // final double      tolerance_absolute, // absolute disparity half-range in each tile
				tolerance_relative_inter, // final double      tolerance_relative, // relative disparity half-range in each tile
				occupancy_inter,          // final double      occupancy,          // fraction of remaining  tiles (<1.0)
				num_passes,               // final int         num_passes,
				max_change,               // final double      max_change,
				-1); //-1); // 1); // 2);                       // final int         debug_level)

		if (debug_level > 0) {
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
		}

		
		/* */
		double [][][] corr2dscene_ref_multi = new double [max_rad + 2][][]; 
		
		corr2dscene_ref_multi[0] = correlate2DSceneToReference(// to match to reference
				clt_parameters.img_dtt, // final ImageDttParameters  imgdtt_params,   // Now just extra correlation parameters, later will include, most others
				scene_QuadCLT,          // final QuadCLT     scene_QuadClt,
				reference_QuadCLT,      // final QuadCLT     reference_QuadClt,
				scene_tiles,            // final double [][][] scene_tiles,     // prepared with prepareSceneTiles()
				reference_tiles,        // final double [][][] reference_tiles, // prepared with prepareReferenceTiles() + fillTilesNans(); - combine?
//				final double [][] flowXY,      // per macro tile {mismatch in IMAGE_PIXELS in X and Y directions (pre-shift)
				flowXY_frac,            // final double [][] flowXY_frac, // X, YH fractional shift [-0.5,0.5) to implement with FD rotations
				chn_weights,            // final double []   chn_weights, // absolute, starting from strength (strength,r,b,g)
				corr_sigma,             // final double      corr_sigma,
				fat_zero,               //  final double      fat_zero,
				false,                  // final boolean     late_normalize,
				combine_empty_only,     //  final boolean     combine_empty_only, // only use neighbor correlations for empty tiles (false - any)
				0.0,     // final double      combine_dradius
				tolerance_absolute_inter_macro, // final double      tolerance_absolute, // absolute disparity half-range to consolidate tiles
				tolerance_relative_inter_macro, // final double      tolerance_relative, // relative disparity half-range to consolidate tiles
				-1); // 1); // final int         debug_level)

		for (int irad = 0; irad <= 3; irad++) {
			corr2dscene_ref_multi[irad+1]= correlate2DSceneToReference(// to match to reference
				clt_parameters.img_dtt, // final ImageDttParameters  imgdtt_params,   // Now just extra correlation parameters, later will include, most others
				scene_QuadCLT,          // final QuadCLT     scene_QuadClt,
				reference_QuadCLT,      // final QuadCLT     reference_QuadClt,
				scene_tiles,            // final double [][][] scene_tiles,     // prepared with prepareSceneTiles()
				reference_tiles,        // final double [][][] reference_tiles, // prepared with prepareReferenceTiles() + fillTilesNans(); - combine?
				flowXY_frac,            // final double [][] flowXY_frac, // X, YH fractional shift [-0.5,0.5) to implement with FD rotations
				chn_weights,            // final double []   chn_weights, // absolute, starting from strength (strength,r,b,g)
				corr_sigma,             // final double      corr_sigma,
				fat_zero,               //  final double      fat_zero,
				true,                   // final boolean     late_normalize,
				combine_empty_only,     //  final boolean     combine_empty_only, // only use neighbor correlations for empty tiles (false - any)
				irad + frac_radius,     // final double      combine_dradius
				tolerance_absolute_inter_macro, // final double      tolerance_absolute, // absolute disparity half-range to consolidate tiles
				tolerance_relative_inter_macro, // final double      tolerance_relative, // relative disparity half-range to consolidate tiles
				-1); // 1); // final int         debug_level)
		}
		
		if (debug_level > 0) {
			showCorrTiles(
					"scene:"+scene_QuadCLT.getImageName()+"-ref"+reference_QuadCLT.getImageName(),  //  String title,
					corr2dscene_ref_multi,           // double [][] source_tiles,
					tilesX/transform_size,     // int         tilesX,
					(2 * transform_size - 1),  // int         tile_width,
					(2 * transform_size - 1)); // int         tile_height) // extra margins over 16x16 tiles to accommodate distorted destination tiles
		}
		//reference_tiles
		//double [][][] scene_tiles			
		double [][][][] scene_to_ref = {reference_tiles, scene_tiles};
		if (debug_level > 0) {
			showMacroTiles(
					"tiles_scene-"+scene_QuadCLT.getImageName()+"-ref"+reference_QuadCLT.getImageName(),// String title,
					scene_to_ref,      // double [][][][] source_tiles_sets,
					reference_QuadCLT, // final QuadCLT qthis,
					0);                // final int     margin) // extra margins over 16x16 tiles to accommodate distorted destination tiles
		}
		String flowXYS_title =  (debug_level > 0)?("vectorXYS_"+scene_QuadCLT.getImageName()+"-ref"+reference_QuadCLT.getImageName()):null;
		int         best_num = 4; // use  4 best neighbors to calculate std deviation
		double      ref_stdev = 5.0; // strength 0.5 if standard deviation of best neighbors to tile difference is this.
		
		double [][] vectorXYConfidence =  attachVectorConfidence(
				flowXY,         // final double [][] flowXY,
				macroTilesX,    // final int         width,
				best_num,       // final int         best_num,
				ref_stdev,      // final double      ref_stdev,
				flowXYS_title); // final String      debug_title);    
		
		
		
/*		
		public double [][] attachVectorConfidence(
				final double [][] flowXY,
				final int         width,
				final int         best_num,
				final double      ref_stdev,
				final String      debug_title)
		*/
		
		
		
		
		if (debug_level > 100) {
			double [][][] vectorsXYS = new double [corr2dscene_ref_multi.length][][];
			for (int i = 0; i < vectorsXYS.length; i++) {
				vectorsXYS[i] = correlation2DToVectors_CM(
						corr2dscene_ref_multi[i], // final double [][] corr2d_tiles, // per 2d calibration tiles (or nulls)
						transform_size, // final int         transform_size,
						iradius, // final int         iradius,      // half-size of the square to process 
						dradius,      // final double      dradius,      // weight calculation (1/(r/dradius)^2 + 1)
						refine_num,   // final int         refine_num,   // number of iterations to apply weights around new center
						1); //final int         debug_level)
			}
			if (debug_level > -1) {
				String dbg_title = "dXdYS-"+scene_QuadCLT.getImageName()+"-"+reference_QuadCLT.getImageName();
				int slices = vectorsXYS.length;
				String [] dbg_titles = new String[slices * 3];
				double [][] dbg_img = new double [dbg_titles.length][macroTilesX*macroTilesY];
				for (int slice = 0; slice < vectorsXYS.length; slice++) {
					dbg_titles[slice + 0 * slices] = "dX"+slice;
					dbg_titles[slice + 1 * slices] = "dY"+slice;
					dbg_titles[slice + 2 * slices] = "S"+slice;
					Arrays.fill(dbg_img[slice + slices * 0], Double.NaN);
					Arrays.fill(dbg_img[slice + slices * 1], Double.NaN);
					Arrays.fill(dbg_img[slice + slices * 2], Double.NaN);
					for (int i = 0; i < vectorsXYS[slice].length; i++) if (vectorsXYS[slice][i] != null){
						dbg_img[slice + 0 * slices][i] = vectorsXYS[slice][i][0];
						dbg_img[slice + 1 * slices][i] = vectorsXYS[slice][i][1];
						dbg_img[slice + 2 * slices][i] = vectorsXYS[slice][i][2];
					}
				}
				(new ShowDoubleFloatArrays()).showArrays(
						dbg_img,
						macroTilesX,
						macroTilesY,
						true,
						dbg_title,
						dbg_titles);


			}
			int selected_index = 1; // single, post-norm
			double [][] flowXY1 = recalculateFlowXY(
					flowXY, // final double [][] currentFlowXY,
					vectorsXYS[selected_index], // final double [][] corr_vectorsXY,
					magic_scale/transform_size); // final double      magic_scale) // 0.85 for CM

			double [][][] scene_tiles1 = prepareSceneTiles(// to match to reference
					// null for {scene,reference}{xyz,atr} uses instances globals 
					camera_xyz0,              // final double []   scene_xyz,     // camera center in world coordinates
					camera_atr0,              // final double []   scene_atr,     // camera orientation relative to world frame
					scene_QuadCLT,                    // final QuadCLT     scene_QuadClt,
					reference_QuadCLT,                    // final QuadCLT     reference_QuadClt,
					reference_tiles,          // final double [][][] reference_tiles, // prepared with prepareReferenceTiles() + fillTilesNans();
					flowXY1,                  // final double [][] flowXY, // per macro tile {mismatch in image pixels in X and Y directions
					flowXY_frac,              // final double [][] flowXY_frac, // should be initialized as [number of macro tiles][] - returns fractional shifts [-0.5, 0.5)
					tolerance_absolute_inter, // final double      tolerance_absolute, // absolute disparity half-range in each tile
					tolerance_relative_inter, // final double      tolerance_relative, // relative disparity half-range in each tile
					occupancy_inter,          // final double      occupancy,          // fraction of remaining  tiles (<1.0)
					num_passes,               // final int         num_passes,
					max_change,               // final double      max_change,
					1); //-1); // 1); // 2);                       // final int         debug_level)

			// single, late
			corr2dscene_ref_multi[0] = correlate2DSceneToReference(// to match to reference
					clt_parameters.img_dtt, // final ImageDttParameters  imgdtt_params,   // Now just extra correlation parameters, later will include, most others
					scene_QuadCLT,          // final QuadCLT     scene_QuadClt,
					reference_QuadCLT,      // final QuadCLT     reference_QuadClt,
					scene_tiles1,            // final double [][][] scene_tiles,     // prepared with prepareSceneTiles()
					reference_tiles,        // final double [][][] reference_tiles, // prepared with prepareReferenceTiles() + fillTilesNans(); - combine?
					flowXY_frac,            // final double [][] flowXY_frac, // X, YH fractional shift [-0.5,0.5) to implement with FD rotations
					chn_weights,            // final double []   chn_weights, // absolute, starting from strength (strength,r,b,g)
					corr_sigma,             // final double      corr_sigma,
					fat_zero,               //  final double      fat_zero,
					false,                  // final boolean     late_normalize,
					combine_empty_only,     //  final boolean     combine_empty_only, // only use neighbor correlations for empty tiles (false - any)
					0.0,                    // final double      combine_dradius
					tolerance_absolute_inter_macro, // final double      tolerance_absolute, // absolute disparity half-range to consolidate tiles
					tolerance_relative_inter_macro, // final double      tolerance_relative, // relative disparity half-range to consolidate tiles
					-1); // 1); // final int         debug_level)

			for (int irad = 0; irad <= 3; irad++) {
				corr2dscene_ref_multi[irad+1]= correlate2DSceneToReference(// to match to reference
						clt_parameters.img_dtt, // final ImageDttParameters  imgdtt_params,   // Now just extra correlation parameters, later will include, most others
						scene_QuadCLT,          // final QuadCLT     scene_QuadClt,
						reference_QuadCLT,      // final QuadCLT     reference_QuadClt,
						scene_tiles1,           // final double [][][] scene_tiles,     // prepared with prepareSceneTiles()
						reference_tiles,        // final double [][][] reference_tiles, // prepared with prepareReferenceTiles() + fillTilesNans(); - combine?
						flowXY_frac,            // final double [][] flowXY_frac, // X, YH fractional shift [-0.5,0.5) to implement with FD rotations
						chn_weights,            // final double []   chn_weights, // absolute, starting from strength (strength,r,b,g)
						corr_sigma,             // final double      corr_sigma,
						fat_zero,               //  final double      fat_zero,
						true, // final boolean     late_normalize,
						combine_empty_only,     //  final boolean     combine_empty_only, // only use neighbor correlations for empty tiles (false - any)
						irad + frac_radius,     // final double      combine_dradius
						tolerance_absolute_inter_macro, // final double      tolerance_absolute, // absolute disparity half-range to consolidate tiles
						tolerance_relative_inter_macro, // final double      tolerance_relative, // relative disparity half-range to consolidate tiles
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
						iradius, // final int         iradius,      // half-size of the square to process 
						dradius,      // final double      dradius,      // weight calculation (1/(r/dradius)^2 + 1)
						refine_num,   // final int         refine_num,   // number of iterations to apply weights around new center
						1); //final int         debug_level)
			}
			if (debug_level > -1) {
				String dbg_title = "dXdYS1-"+scene_QuadCLT.getImageName()+"-"+reference_QuadCLT.getImageName();
				int slices = vectorsXYS1.length;
				String [] dbg_titles = new String[slices * 3];
				double [][] dbg_img = new double [dbg_titles.length][macroTilesX*macroTilesY];
				for (int slice = 0; slice < vectorsXYS1.length; slice++) {
					dbg_titles[slice + 0 * slices] = "dX"+slice;
					dbg_titles[slice + 1 * slices] = "dY"+slice;
					dbg_titles[slice + 2 * slices] = "S"+slice;
					Arrays.fill(dbg_img[slice + slices * 0], Double.NaN);
					Arrays.fill(dbg_img[slice + slices * 1], Double.NaN);
					Arrays.fill(dbg_img[slice + slices * 2], Double.NaN);
					for (int i = 0; i < vectorsXYS1[slice].length; i++) if (vectorsXYS1[slice][i] != null){
						dbg_img[slice + 0 * slices][i] = vectorsXYS1[slice][i][0];
						dbg_img[slice + 1 * slices][i] = vectorsXYS1[slice][i][1];
						dbg_img[slice + 2 * slices][i] = vectorsXYS1[slice][i][2];
					}
				}
				(new ShowDoubleFloatArrays()).showArrays(
						dbg_img,
						macroTilesX,
						macroTilesY,
						true,
						dbg_title,
						dbg_titles);


			}
			double [][] flowXY2 = recalculateFlowXY(
					flowXY1, // final double [][] currentFlowXY,
					vectorsXYS1[selected_index], // final double [][] corr_vectorsXY,
					magic_scale/transform_size); // final double      magic_scale) // 0.85 for CM

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
					tolerance_absolute_inter, // final double      tolerance_absolute, // absolute disparity half-range in each tile
					tolerance_relative_inter, // final double      tolerance_relative, // relative disparity half-range in each tile
					occupancy_inter,          // final double      occupancy,          // fraction of remaining  tiles (<1.0)
					num_passes,               // final int         num_passes,
					max_change,               // final double      max_change,
					1); //-1); // 1); // 2);                       // final int         debug_level)
		}
		
		
		
		
		return pair;
	}

	
	
	public double [][] transformCameraVew(
			QuadCLT   reference_QuadClt,
			QuadCLT   camera_QuadClt,
			double [] camera_xyz, // camera center in world coordinates
			double [] camera_atr, // camera orientation relative to world frame
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
		double [][] dsrbg_camera =    camera_QuadClt.getDSRBG();
///		double [][] dsrbg_reference = reference_QuadClt.getDSRBG();
		double [][] ds =        new double [dsrbg_camera.length][stiles];
		for (int i = 0; i <ds.length; i++) {
			for (int j = 0; j <ds[i].length; j++) {
				ds[i][j] = Double.NaN;
			}
		}
		
		ErsCorrection ersReferenceCorrection = reference_QuadClt.getErsCorrection();
		ersReferenceCorrection.setupERS(); // just in case - setUP using instance paRAMETERS
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
					double [] pXpYD = ersReferenceCorrection.getImageCoordinatesReferenceERS( // ersCorrection - reference
							
							camera_QuadClt, // QuadCLT cameraQuadCLT, // camera station that got image to be to be matched 
							centerX,        // double px,                // pixel coordinate X in the reference view
							centerY,        // double py,                // pixel coordinate Y in the reference view
							disparity,      // double disparity,         // reference disparity 
							true,           // boolean distortedView,    // This camera view is distorted (diff.rect), false - rectilinear
							ZERO3,          // double [] reference_xyz,  // this view position in world coordinates (typically ZERO3)
							ZERO3,          // double [] reference_atr,  // this view orientation relative to world frame  (typically ZERO3)
							true,           // boolean distortedCamera,  // camera view is distorted (false - rectilinear)
							camera_xyz,     // double [] camera_xyz,     // camera center in world coordinates
							camera_atr,     // double [] camera_atr,     // camera orientation relative to world frame
							LINE_ERR);       // double    line_err)       // threshold error in scan lines (1.0)
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
//				dsrbg_out[i] = reference_QuadClt.fillNaNGaps(dsrbg_out[i], num_passes, rel_num_passes, 100); // threadsMax);
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
	
	

}
