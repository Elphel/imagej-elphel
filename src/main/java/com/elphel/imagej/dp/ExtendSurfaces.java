package com.elphel.imagej.dp;
/**
 **
 ** ExtendSurfaces predict disparity in the unknown areas
 **
 ** Copyright (C) 2017 Elphel, Inc.
 **
 ** -----------------------------------------------------------------------------**
 **  
 **  ExtendSurfaces.java is free software: you can redistribute it and/or modify
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
import java.util.concurrent.atomic.AtomicInteger;

import com.elphel.imagej.common.ShowDoubleFloatArrays;


public class ExtendSurfaces {
	TileProcessor tp;
	public ExtendSurfaces(TileProcessor tp){
		this.tp = tp;
	}

	public int smoothNew( // return maximal correction value 
			// should first be initialized with all not known = Double.NaN
			final double []  values,     // will be modified, Double.NaN is not yet assigned
//			final double []  strengths,  // cell weights
			final boolean [] known,      // cells with known values (will not be modified)
			final boolean [] new_cells,  // cells that need to  be calculated
			final int        smpl_size, // == 5
			final int        min_points, // == 3
			final boolean    use_wnd,   // use window function fro the neighbors 
			final double     tilt_cost,
			final double     split_threshold, // if full range of the values around the cell higher, need separate fg, bg
			final double     same_range,      // modify 
			final double     diff_continue,   // maximal difference from the old value (for previously defined tiles
			final double     max_abs_tilt, //   = 2.0; // pix per tile
			final double     max_rel_tilt, //   = 0.2; // (pix / disparity) per tile
			final boolean    en_normal, // expand sngle-range cells (that have all neighbors within split_threshold)
			final boolean    en_far,    // expand background - use neighbors cells within same_range from the lowest
			// en_near is only valid if !en_far 
			final boolean    en_near,   // expand foreground - use neighbors cells within same_range from the highest
			final int        max_tries, // maximal number of smoothing steps
			final double     final_diff, // maximal change to finish iterations
			final boolean    new_only, // expand only Double.NaN cells (not yet assigned), use for first passes
			final int        dbg_x,
			final int        dbg_y,
			final int        debugLevel)
	{
		for (int nTile = 0; nTile < values.length; nTile++){
			if (!known[nTile]) values [nTile] = Double.NaN;
		}
		for (int num_try = 0; num_try < max_tries; num_try ++){
			double max_diff = smoothNewStep( // return maximal correction value 
					// should first be initialized with all not known = Double.NaN
					values,          // final double []  values,     // will be modified, Double.NaN is not yet assigned
					new_cells,       // final boolean [] new_cells,  // cells that need to  be calculated
					smpl_size,       // final int        smpl_size,  // == 5
					min_points,      // final int        min_points, // == 3
					use_wnd,         // final boolean    use_wnd,    // use window function fro the neighbors 
					tilt_cost,       // final double     tilt_cost,
					split_threshold, // final double     split_threshold, // if full range of the values around the cell higher, need separate fg, bg
					same_range,      // final double     same_range,      // modify 
					diff_continue,   // final double     diff_continue,   // maximal difference from the old value (for previously defined tiles
					max_abs_tilt,    // final double     max_abs_tilt, //   = 2.0; // pix per tile
					max_rel_tilt,    // final double     max_rel_tilt, //   = 0.2; // (pix / disparity) per tile
					en_normal,       // final boolean    en_normal,  // expand sngle-range cells (that have all neighbors within split_threshold)
					en_far,          // final boolean    en_far,     // expand background - use neighbors cells within same_range from the lowest
					// en_near is only valid if !en_far 
					en_near,         // final boolean    en_near,    // expand foreground - use neighbors cells within same_range from the highest
					new_only,        // final boolean    new_only,   // expand only Double.NaN cells (not yet assigned), use for first passes
					true,            // final boolean    no_zero,  // do not use values[i] == 0.0, exact 0.0 is only background
					dbg_x,           // final int        dbg_x,
					dbg_y,           // final int        dbg_y,
					debugLevel);     // final int        debugLevel);

			if (max_diff <= final_diff) break; 
		}
		
		
		
		int num_new = 0;
		for (int nTile = 0; nTile < values.length; nTile++){
			if (new_cells[nTile] && !Double.isNaN(values[nTile])) num_new++;
		}
		return num_new;
	}
	
	private int [] getNeibIndices(
			int size)
	{
		final int [] neib_indices = new int [8]; // indices of the sample array immediately around the center;
		int ni = 0;
		for (int sy = -1; sy <= 1; sy++) {
			for (int sx = -1; sx <= 1; sx++) {
				if ((sy != 0) || (sx != 0)){
					neib_indices[ni++] = (size/2 + sx) + (size/2 + sy) * size; 
				}
			}		
		}		
		return neib_indices;
	}
	
	
	private double smoothNewStep( // return maximal correction value 
			// should first be initialized with all not known = Double.NaN
			final double []  values,     // will be modified, Double.NaN is not yet assigned
			final boolean [] new_cells,  // cells that need to  be calculated
			final int        smpl_size, // == 5
			final int        min_points, // == 3
			final boolean    use_wnd,   // use window function for the neighbors 
			final double     tilt_cost,
			final double     split_threshold, // if full range of the values around the cell higher, need separate fg, bg
			final double     same_range,      // modify for old - +/- 0.5?
			final double     diff_continue,   // maximal difference from the old value (for previously defined tiles
			final double     max_abs_tilt, //   = 2.0; // pix per tile
			final double     max_rel_tilt, //   = 0.2; // (pix / disparity) per tile
			final boolean    en_normal, // expand single-range cells (that have all neighbors within split_threshold)
			final boolean    en_far,    // expand background - use neighbors cells within same_range from the lowest
			// en_near is only valid if !en_far 
			final boolean    en_near,   // expand foreground - use neighbors cells within same_range from the highest
			final boolean    new_only, // expand only Double.NaN cells (not yet assigned), use for first passes
			final boolean    no_zero,  // do not use values[i] == 0.0, exact 0.0 is only background
			final int        dbg_x,
			final int        dbg_y,
			final int        debugLevel)
	{
		final int tilesX = tp.getTilesX();
		final int tilesY = tp.getTilesY();
		final int smpl_half = smpl_size / 2; // 2
		// probably should not be needed
		final double thresholdLin = 1.0E-20;  // threshold ratio of matrix determinant to norm for linear approximation (det too low - fail)
		final double [] damping = {tilt_cost,tilt_cost,0.0}; // 0.0 will be applied to average value, tilt_cost - to both tilts
		final int center_index = smpl_half * (smpl_size + 1); 
//		final int smpl_high = smpl_size  - smpl_low; // 3
		final int smpl_len = smpl_size * smpl_size;
		final int dbg_tile = dbg_x + tilesX * dbg_y;
		final double [] smpl_weights = MeasuredLayers.getSampleWindow(smpl_size, !use_wnd);
		final Thread[] threads = ImageDtt.newThreadArray(tp.threadsMax);
		final double []  max_corrs = new double [threads.length];
		final AtomicInteger ai_thread = new AtomicInteger(0);
		final AtomicInteger ai = new AtomicInteger(0);
		final double [] new_vals = values.clone();
		final int [] neib_indices = getNeibIndices(smpl_size);
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				public void run() {
					PolynomialApproximation pa = new PolynomialApproximation();
					int this_thread = ai_thread.getAndIncrement();
					for (int nTile = ai.getAndIncrement(); nTile < new_vals.length; nTile = ai.getAndIncrement()) {
						if (new_cells[nTile] && (!new_only || Double.isNaN(values[nTile]))) {
							int dl = (nTile == dbg_tile) ? debugLevel : -1;
							int tileX = nTile % tilesX;
							int tileY = nTile / tilesX;

//							boolean [] smpl_new = new boolean [smpl_len];
							boolean [] smpl_sel = new boolean [smpl_len];
							double [] smpl_d =  new double [smpl_len];
//							double [] smpl_p =  new double [smpl_len];
//							double [] smpl_w =  new double [smpl_len];
							int num_in_sample = 0;
							for (int sy = 0; sy < smpl_size; sy++) {
								int y = tileY + sy - smpl_half;
								if ((y >= 0) && (y < tilesY)) {
									for (int sx = 0; sx < smpl_size; sx++) {
										int x = tileX + sx - smpl_half;
										if ((x >= 0) && (x < tilesX)) {
											int indx = x + y * tilesX;
//											if (!Double.isNaN(values[indx]) && (indx != center_index)){
											if (!Double.isNaN(values[indx]) && (!no_zero || (values[indx] != 0.0))){
												int indxs = sx + sy * smpl_size;
//												double w = smpl_weights[indxs];
												smpl_sel[indxs] = true;
												smpl_d[indxs] = values[indx];
												num_in_sample++;
											}
										}
									}
								}
							}
							if (dl > 1){
								System.out.println("smoothNewStep(): tileX="+tileX+", tileY="+tileY+", nTile="+nTile+" num_in_sample="+num_in_sample);
								for (int sy = 0; sy < smpl_size; sy++) {
									int y = tileY + sy - smpl_half;
									if ((y >= 0) && (y < tilesY)) {
										for (int sx = 0; sx < smpl_size; sx++) {
											int x = tileX + sx - smpl_half;
											if ((x >= 0) && (x < tilesX)) {
												int indx = x + y * tilesX;
												if (!Double.isNaN(values[indx])){
													int indxs = sx + sy * smpl_size;
													System.out.print(String.format("%8.3f ", smpl_d[indxs]));
												} else {
													System.out.print(String.format("%8s ", "--NaN--"));
												}
											}
										}
									}
									System.out.println();
								}
								
							}
							if (num_in_sample >= min_points) {
								// if center was NaN - diff is the absolute value.
								
								// Two branches - for existing and new cells.
								// New cells may stick to the far/near, already defined - only to similar,
								// otherwise they can all switch from some distant cell during many iteration steps
								boolean process_sample = false;
								if (smpl_sel[center_index]) {
									// Keep only cells within +/- diff_continue from the center value
									for (int indxs = 0; indxs < smpl_len; indxs ++ ) if (smpl_sel[indxs]){
										if (Math.abs(smpl_d[indxs] - smpl_d[center_index]) > diff_continue){
											smpl_sel[indxs] = false;
											num_in_sample--;
										}
									}
									process_sample = (num_in_sample > 1);
								} else {
									// add starting from the known cells. When using 5x5 samples, having immediate neighbor
									// to the center usually means that likely there is enough data to determine the tilt.
									boolean has_neib = false;
									for (int i = 0; i < neib_indices.length; i++) if (smpl_sel[neib_indices[i]]) {
										has_neib = true;
										break;
									}
//									process_sample = false;
									if (has_neib) { // all conditions met, find linear 2-d damped approximation for the center
										int imin = -1, imax = -1;
										for (int indxs = 0; indxs < smpl_len; indxs ++ ) if (smpl_sel[indxs]){
											if ((imin < 0) || (smpl_d[indxs] < smpl_d[imin])) imin = indxs ;
											if ((imax < 0) || (smpl_d[indxs] > smpl_d[imax])) imax = indxs ;
										}
										// does it all fit
										if ((smpl_d[imax] - smpl_d[imin]) <= split_threshold){
											process_sample = en_normal;
										} else {
											if (en_far) {
												double threshold = smpl_d[imin] + same_range;  
												for (int indxs = 0; indxs < smpl_len; indxs ++ ) if (smpl_sel[indxs]){
													if (smpl_d[indxs] > threshold) {
														smpl_sel[indxs] = false;
														num_in_sample--;
													}
												}										
												process_sample = true;
											} else if (en_near){
												double threshold = smpl_d[imax] - same_range;  
												for (int indxs = 0; indxs < smpl_len; indxs ++ ) if (smpl_sel[indxs]){
													if (smpl_d[indxs] < threshold) {
														smpl_sel[indxs] = false;
														num_in_sample--;
													}
												}										
												process_sample = true;
											}
											// removed some points - recheck neighbors
											if (process_sample){
												has_neib = false;
												for (int i = 0; i < neib_indices.length; i++) if (smpl_sel[neib_indices[i]]) {
													has_neib = true;
													break;
												}
												if (!has_neib) process_sample = false; // neighbor was removed, will process
												if (dl > 1){
													System.out.println("smoothNewStep(): tileX="+tileX+", tileY="+tileY+", nTile="+nTile+
															" num_in_sample="+num_in_sample+" has_neib="+has_neib+ " process_sample="+process_sample);
												}

												// with other direction
											}
										}
									}
								}
								if (dl > 1){
									System.out.println("smoothNewStep(): tileX="+tileX+", tileY="+tileY+", nTile="+nTile+
											" num_in_sample="+num_in_sample+" process_sample="+process_sample);
								}
								if (num_in_sample < min_points) {
									process_sample = false;
								}
								if (process_sample) {
									// build linear 2d approximation, damped as there could be co-linear cells
									// or even a single cell
									
									double sw = 0.0;

									double [][][] mdata = new double [num_in_sample][3][];
									int mindx = 0;
									for (int sy = 0; sy < smpl_size; sy++){
										for (int sx = 0; sx < smpl_size; sx++){
											int indxs = sy * smpl_size + sx;
											if (smpl_sel[indxs]) {
												mdata[mindx][0] = new double [2];
												mdata[mindx][0][0] =  sx - smpl_half; // should match integer center
												mdata[mindx][0][1] =  sy - smpl_half;
												mdata[mindx][1] = new double [1];
												mdata[mindx][1][0] = smpl_d[indxs];
												mdata[mindx][2] = new double [1];
												mdata[mindx][2][0] =  smpl_weights[indxs];
												if (mdata[mindx][1][0] == 0.0){
													System.out.println("zero!!!");
												}
												mindx ++;
											}
										}
									}
									double[][] approx2d = pa.quadraticApproximation(
											mdata,
											true,          // boolean forceLinear,  // use linear approximation
											damping,       // double [] damping,
											thresholdLin,  // threshold ratio of matrix determinant to norm for linear approximation (det too low - fail)
											0.0,           // threshold ratio of matrix determinant to norm for quadratic approximation (det too low - fail)
											dl-1);           // int debugLevel;
									if (approx2d == null) {
										System.out.println("A BUG in smoothNewStep() - failed quadraticApproximation()");
										continue;
									}
									double tilt2 = approx2d[0][0]*approx2d[0][0] + approx2d[0][1]*approx2d[0][1]; 
									if (dl > 1){
										System.out.println("smoothNewStep(): tileX="+tileX+" coeffs: D="+
												approx2d[0][0]+" E="+approx2d[0][1]+" F="+approx2d[0][2]+
												", tilt = "+Math.sqrt(tilt2));
									}
									double max_tilt = Math.min(max_abs_tilt, max_rel_tilt * approx2d[0][2]);
									if (tilt2 > (max_tilt * max_tilt)){
										continue;
									}
									double new_center_val = approx2d[0][2]; // Last coefficient, F (of 6)
									double this_corr;
									if (Double.isNaN(values[nTile])) this_corr = Math.abs(new_center_val);
									else  this_corr = Math.abs(new_center_val - values[nTile]);
									new_vals[nTile] = new_center_val;
									max_corrs[this_thread] = Math.max(max_corrs[this_thread], this_corr);
								}
							}
						}
					}	
				}
			};
		}		      
		ImageDtt.startAndJoin(threads);
		System.arraycopy(new_vals, 0, values, 0, new_vals.length);
		double max_corr = 0;
		for (int i = 0; i < max_corrs.length; i++) max_corr = Math.max(max_corr,max_corrs[i]);
		return max_corr;
		
	}
	
	private boolean [] startDiscontinuity(
			final double []  values,     // will be modified, Double.NaN is not yet assigned
			final boolean [] known,      // cells with new values (will be modified)
			final int        smpl_size, // == 5
			final double     split_threshold, // if full range of the values around the cell higher, need separate fg, bg
			final boolean    extend_far,    // extend far (false - extend near)
			final boolean    no_zero,  // do not use values[i] == 0.0, exact 0.0 is only background			
			final int        dbg_x,
			final int        dbg_y,
			final int        debugLevel)
	{
		final int tilesX = tp.getTilesX();
		final int tilesY = tp.getTilesY();
		final int smpl_half = smpl_size / 2; // 2
		final boolean [] discont_cells = new boolean [tilesX * tilesY]; 
		final Thread[] threads = ImageDtt.newThreadArray(tp.threadsMax);
		final int []  num_new = new int [threads.length];
		final AtomicInteger ai_thread = new AtomicInteger(0);
		final AtomicInteger ai = new AtomicInteger(0);
		final int dbg_tile = dbg_x + tilesX * dbg_y;
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				public void run() {
					int this_thread = ai_thread.getAndIncrement();
					for (int nTile = ai.getAndIncrement(); nTile < known.length; nTile = ai.getAndIncrement()) if (known[nTile]){
						double threshold = extend_far ? (values[nTile] - split_threshold): (values[nTile] + split_threshold);
						int dl = (nTile == dbg_tile) ? debugLevel : -1;

//						int dl = (nTile == dbg_tile) ? debugLevel : -1;
						int tileX = nTile % tilesX;
						int tileY = nTile / tilesX;
						label_tile:
						{	
							for (int sy = 0; sy < smpl_size; sy++) {
								int y = tileY + sy - smpl_half;
								if ((y >= 0) && (y < tilesY)) {
									for (int sx = 0; sx < smpl_size; sx++) {
										int x = tileX + sx - smpl_half;
										if ((x >= 0) && (x < tilesX)) {
											int indx = x + y * tilesX;
											if (known[indx] &&  (!no_zero || (values[indx] != 0.0))) {
												if (	(extend_far && (values[indx] < threshold)) ||
														(!extend_far && (values[indx] > threshold))){
													discont_cells[nTile] = true;
													num_new[this_thread] ++;
													break label_tile;
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
		ai.set(0);
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				public void run() {
					int this_thread = ai_thread.getAndIncrement();
					for (int nTile = ai.getAndIncrement(); nTile < known.length; nTile = ai.getAndIncrement()){
						if (discont_cells[nTile]){
//							values[nTile] = Double.NaN;
							known[nTile] = false;
						}
						if (!known[nTile]){
							values[nTile] = Double.NaN;
						}
						
					}	
				}
			};
		}
		ImageDtt.startAndJoin(threads);

		for (int i = 0; i < num_new.length; i++) if (num_new[i] > 0){
			return discont_cells;
		}
		return null; // nothing found
	}

	private boolean [] startExpanding(
			final double []  values,     // will be modified, Double.NaN is not yet assigned
			final boolean [] known,      // cells with new values (will not be modified)
			final int        smpl_size, // == 5
			final boolean    no_zero,  // do not use values[i] == 0.0, exact 0.0 is only background			
			final int        dbg_x,
			final int        dbg_y,
			final int        debugLevel)
	{
		final int tilesX = tp.getTilesX();
		final int tilesY = tp.getTilesY();
		final boolean [] expanding = new boolean [tilesX * tilesY]; 
		final Thread[] threads = ImageDtt.newThreadArray(tp.threadsMax);
		final int []  num_new = new int [threads.length];
		final AtomicInteger ai_thread = new AtomicInteger(0);
		final AtomicInteger ai = new AtomicInteger(0);
		final int dbg_tile = dbg_x + tilesX * dbg_y;
		final TileNeibs tnSurface = new TileNeibs(tilesX, tilesY);
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				public void run() {
					int this_thread = ai_thread.getAndIncrement();
					for (int nTile = ai.getAndIncrement(); nTile < known.length; nTile = ai.getAndIncrement()) if (!known[nTile]){
						int dl = (nTile == dbg_tile) ? debugLevel : -1;
						boolean has_exp_neib = false;
						for (int dir = 0; dir <8; dir++){
							int nTile1 = tnSurface.getNeibIndex(nTile, dir);
							if ((nTile1 >= 0) && known[nTile1] && !Double.isNaN(values[nTile1]) && (!no_zero || (values[nTile1] != 0.0))){
								has_exp_neib = true;
								break;
							}
						}
						if (has_exp_neib) {
							expanding[nTile] = true;
							num_new[this_thread] ++;
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
					for (int nTile = ai.getAndIncrement(); nTile < known.length; nTile = ai.getAndIncrement()){
						if (expanding[nTile] || !known[nTile]){
							values[nTile] = Double.NaN;
						}
					}	
				}
			};
		}
		ImageDtt.startAndJoin(threads);
		if (debugLevel > 2){
			String [] titles = {"values","known", "expanding"};
			double [][] dbg_img = new double [titles.length][];
			dbg_img[0] = values;
			dbg_img[1] = new double [tilesX*tilesY]; 
			dbg_img[2] = new double [tilesX*tilesY];
			for (int i = 0; i < dbg_img[1].length; i++){
				dbg_img[1][i] = known[i] ?     50.0: 0.0;
				dbg_img[2][i] = expanding[i] ? 50.0: 0.0;
			}
			
			(new ShowDoubleFloatArrays()).showArrays(dbg_img,  tilesX, tilesY, true, "startExpanding",titles);

		}
		
		for (int i = 0; i < num_new.length; i++) if (num_new[i] > 0){
			return expanding;
		}
		return null; // nothing found
	}
	
	
	private int expandDiscontinuityStep(
			final double []  values,     // will be modified, Double.NaN is not yet assigned
			final boolean [] known,      // cells with known values (will be modified)
			final boolean [] prohibited, // cells that can not be used (tried before with similar disparity)
			final boolean [] expanding,  // cells with new values (will be modified)
			final int        smpl_size, // == 5
			final double     split_threshold, // if full range of the values around the cell higher, need separate fg, bg
			final boolean    extend_far,    // extend far (false - extend near)
			final boolean    unknown_only,  // pure expanding, not over previously found values
			final int        dbg_x,
			final int        dbg_y,
			final int        debugLevel)
	{
		final int tilesX = tp.getTilesX();
		final int tilesY = tp.getTilesY();
		final int smpl_half = smpl_size / 2; // 2
		final boolean [] expanding_prev = expanding.clone();
		final Thread[] threads = ImageDtt.newThreadArray(tp.threadsMax);
		final int []  num_new = new int [threads.length];
		final AtomicInteger ai_thread = new AtomicInteger(0);
		final AtomicInteger ai = new AtomicInteger(0);
		final int dbg_tile = dbg_x + tilesX * dbg_y;
		final TileNeibs tnSurface = new TileNeibs(tilesX, tilesY);
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				public void run() {
					int this_thread = ai_thread.getAndIncrement();
					for (int nTile = ai.getAndIncrement(); nTile < known.length; nTile = ai.getAndIncrement()){
						if (!expanding_prev[nTile] && ((prohibited == null) || !prohibited[nTile])) {
							int dl = (nTile == dbg_tile) ? debugLevel : -1;
							boolean has_exp_neib = false;
							for (int dir = 0; dir <8; dir++){
								int nTile1 = tnSurface.getNeibIndex(nTile, dir);
								if ((nTile1 >= 0) && expanding_prev[nTile1] && !Double.isNaN(values[nTile1])){
									has_exp_neib = true;
									break;
								}
							}
							if (has_exp_neib) {
								if (!known[nTile]){
									expanding[nTile] = true;
									num_new[this_thread] ++;
								} else if (!unknown_only){
									double threshold = extend_far ? (values[nTile] - split_threshold): (values[nTile] + split_threshold);
									int tileX = nTile % tilesX;
									int tileY = nTile / tilesX;
									label_tile:
									{	
										for (int sy = 0; sy < smpl_size; sy++) {
											int y = tileY + sy - smpl_half;
											if ((y >= 0) && (y < tilesY)) {
												for (int sx = 0; sx < smpl_size; sx++) {
													int x = tileX + sx - smpl_half;
													if ((x >= 0) && (x < tilesX)) {
														int indx = x + y * tilesX;
														if ((known[indx] || (expanding_prev[indx]) && !Double.isNaN(values[indx]))) {
															if (	(extend_far && (values[indx] < threshold)) ||
																	(!extend_far && (values[indx] > threshold))){
																expanding[nTile] = true;
																num_new[this_thread] ++;
																break label_tile;
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
					}	
				}
			};
		}		      
		ImageDtt.startAndJoin(threads);
		ai.set(0);
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				public void run() {
//					int this_thread = ai_thread.getAndIncrement();
					for (int nTile = ai.getAndIncrement(); nTile < known.length; nTile = ai.getAndIncrement()){
						if (expanding[nTile] && !expanding_prev[nTile]){
							values[nTile] = Double.NaN;
							known[nTile] = false;
						}
					}	
				}
			};
		}
		ImageDtt.startAndJoin(threads);
		
		int total_new = 0;
		for (int i = 0; i < num_new.length; i++){
			total_new +=num_new[i];
		}
		return total_new;
	}
	
	public boolean [] expandDiscontinuity( // return maximal correction value 
			final ArrayList <CLTPass3d> passes,
			final int                   firstPass,
			final int                   lastPassPlus1,
			final double []  values,     // will be modified, Double.NaN is not yet assigned
			final boolean [] known,      // cells with known values (will not be modified)
			final int        num_steps,  // how far to extend
			final int        smpl_size, // == 5
			final int        min_points, // == 3
			final boolean    use_wnd,   // use window function for the neighbors 
			final double     tilt_cost,
			final double     split_threshold, // if full range of the values around the cell higher, need separate fg, bg
			final double     same_range,      // modify
			final double     diff_continue,   // maximal difference from the old value (for previously defined tiles
			final double     max_abs_tilt, //   = 2.0; // pix per tile
			final double     max_rel_tilt, //   = 0.2; // (pix / disparity) per tile
			final boolean    extend_far,    // extend far (false - extend near)
			final int        max_tries, // maximal number of smoothing steps
			final double     final_diff, // maximal change to finish iterations 
		    final double     grow_disp_min,
			final double     grow_disp_max,
			final double     grow_disp_step,
			final double     unique_pre_tolerance, // usually larger than clt_parameters.unique_tolerance
			final boolean    try_next, // try next disparity range if the current one was already tried
			final int        dbg_x,
			final int        dbg_y,
			final int        debugLevel)
	{
		boolean [] expanding = startDiscontinuity(
				values,           // final double []  values,     // will be modified, Double.NaN is not yet assigned
				known,            // final boolean [] known,      // cells with known values (will not be modified)
				smpl_size,        // final int        smpl_size, // == 5
				split_threshold,  // final double     split_threshold, // if full range of the values around the cell higher, need separate fg, bg
				extend_far,       // final boolean    extend_far,    // extend far (false - extend near)
				true,             // final boolean    no_zero,  // do not use values[i] == 0.0, exact 0.0 is only background			
				dbg_x,           // final int        dbg_x,
				dbg_y,           // final int        dbg_y,
				debugLevel);     // final int        debugLevel)
		if (expanding == null){
			return null; // no discontinuity seeds
		}

		boolean [] prohibited = tp.getTriedBefore( 
				passes,         // final ArrayList <CLTPass3d> passes,
				firstPass,      // final int                   firstPass,
				lastPassPlus1,  // final int                   lastPassPlus1,
				expanding,      // final boolean []            selected,
				null,           // final boolean []            prohibited,
				values,         // final double []             disparity,
				true,           // final boolean               mod_selected,
				true,           // final boolean               mod_disparity,
				false,          // final boolean               en_near,
				false,          // final boolean               en_far,
				grow_disp_min,  // final double                grow_disp_min,
				grow_disp_max,  // final double                grow_disp_max,
				grow_disp_step, // final double                grow_disp_step,
				try_next ? 0.0: unique_pre_tolerance); // final double unique_pre_tolerance) // disable for try_next
		
		for (int nstep = 0; nstep < num_steps; nstep++){
			int num_new_cells = 0;
			if (nstep > 0){
				num_new_cells = expandDiscontinuityStep(
						values,          // final double []  values,     // will be modified, Double.NaN is not yet assigned
						known,           // final boolean [] known,      // cells with known values (will not be modified)
						prohibited,      // final boolean [] prohibited, // cells that can not be used (tried before with similar disparity)
						expanding,       // final boolean [] expanding,  // cells with known values (will not be modified)
						smpl_size,       //final int        smpl_size, // == 5
						split_threshold, // final double     split_threshold, // if full range of the values around the cell higher, need separate fg, bg
						extend_far,      // final boolean    extend_far,    // extend far (false - extend near)
						false,           // inal boolean    unknown_only,  // pure expanding, not over previously found values
						dbg_x,           // final int        dbg_x,
						dbg_y,           // final int        dbg_y,
						debugLevel);     // final int        debugLevel)

				tp.getTriedBefore(
						passes,          // final ArrayList <CLTPass3d> passes,
						firstPass,       // final int                   firstPass,
						lastPassPlus1,   // final int                   lastPassPlus1,
						expanding,       // final boolean []            selected,
						prohibited,      //  final boolean []            prohibited,
						values,          // final double []             disparity,
						true,            //  final boolean               mod_selected,
						true,            //  final boolean               mod_disparity,
						false,          // final boolean               en_near,
						false,          // final boolean               en_far,
						grow_disp_min,  // final double                grow_disp_min,
						grow_disp_max,  // final double                grow_disp_max,
						grow_disp_step, // final double                grow_disp_step,
						try_next ? 0.0: unique_pre_tolerance); // final double unique_pre_tolerance) // disable for try_next
				
			}
			int mt = 2 * (nstep + 1); // max_tries; // or use k * nstep ???
			// increase tries last iteration?
			if (((num_new_cells == 0) && (nstep > 0)) || (nstep == (num_steps -1))) {
				mt = max_tries; // last time smooth more 
			}
			for (int num_try = 0; num_try < mt; num_try ++){
				if (debugLevel > 1){
					System.out.print("expandDiscontinuity(): nstep="+nstep+" num_try="+num_try+" ");
					if (num_try == (mt - 1)){
						System.out.print("expandDiscontinuity(): LAST try ");
					}
				}
				double max_diff = smoothNewStep( // return maximal correction value 
						// should first be initialized with all not known = Double.NaN
						values,          // final double []  values,     // will be modified, Double.NaN is not yet assigned
						expanding,       // final boolean [] new_cells,  // cells that need to  be calculated
						smpl_size,       // final int        smpl_size, // == 5
						min_points,      // final int        min_points, // == 3
						use_wnd,         // final boolean    use_wnd,   // use window function fro the neighbors 
						tilt_cost,       // final double     tilt_cost,
						split_threshold, // final double     split_threshold, // if full range of the values around the cell higher, need separate fg, bg
						same_range,      // final double     same_range,      // modify
						diff_continue,   // final double     diff_continue,   // maximal difference from the old value (for previously defined tiles
						max_abs_tilt,    // final double     max_abs_tilt, //   = 2.0; // pix per tile
						max_rel_tilt,    // final double     max_rel_tilt, //   = 0.2; // (pix / disparity) per tile
						true,            // final boolean    en_normal, // expand sngle-range cells (that have all neighbors within split_threshold)
						extend_far,      // final boolean    en_far,    // expand background - use neighbors cells within same_range from the lowest
						// en_near is only valid if !en_far 
						!extend_far,     // final boolean    en_near,   // expand foreground - use neighbors cells within same_range from the highest
						(num_try == 0),  // new_only,       // first run - set new cells, later - smooth all
						true,            // final boolean    no_zero,  // do not use values[i] == 0.0, exact 0.0 is only background
						dbg_x,           // final int        dbg_x,
						dbg_y,           // final int        dbg_y,
						debugLevel);     // final int        debugLevel);
				if (debugLevel > 1){
					System.out.println("expandDiscontinuity() -> max_diff="+max_diff);
				}
				if (max_diff <= final_diff) break; 
			}
			if (debugLevel > 1){
				System.out.println("expandDiscontinuity(): nstep="+nstep+" mt="+mt);
			}
		}
		if (try_next){
			tp.getTriedBefore(
				passes,          // final ArrayList <CLTPass3d> passes,
				firstPass,       // final int                   firstPass,
				lastPassPlus1,   // final int                   lastPassPlus1,
				expanding,       // final boolean []            selected,
				prohibited,      //  final boolean []            prohibited,
				values,          // final double []             disparity,
				true,            //  final boolean               mod_selected,
				true,            //  final boolean               mod_disparity,
				// TODO: some version to enable both en_near and en_far? or make en_far assume en_near and make it first?
				!extend_far,     // final boolean               en_near,
				extend_far,      // final boolean               en_far,
				grow_disp_min,   // final double                grow_disp_min,
				grow_disp_max,   // final double                grow_disp_max,
				grow_disp_step,  // final double                grow_disp_step,
				unique_pre_tolerance); // final double unique_pre_tolerance)
		}		
		return expanding;
	}

	
	public boolean[] expandKnown( // return maximal correction value 
			// should first be initialized with all not known = Double.NaN
			final ArrayList <CLTPass3d> passes,
			final int                   firstPass,
			final int                   lastPassPlus1,
			final double []  values,     // will be modified, Double.NaN is not yet assigned
			final boolean [] known,      // cells with known values (will not be modified)
			final int        num_steps,  // how far to extend
			final int        smpl_size, // == 5
			final int        min_points, // == 3
			final boolean    use_wnd,   // use window function fro the neighbors 
			final double     tilt_cost,
			final double     split_threshold, // if full range of the values around the cell higher, need separate fg, bg
			final double     same_range,      // modify
			final double     diff_continue,   // maximal difference from the old value (for previously defined tiles
			final double     max_abs_tilt, //   = 2.0; // pix per tile
			final double     max_rel_tilt, //   = 0.2; // (pix / disparity) per tile
			final boolean    smooth_only,   // Do not expand ambiguous cells (with discontinuity)
			final boolean    extend_far,    // extend far (false - extend near)
			final int        max_tries, // maximal number of smoothing steps
			final double     final_diff, // maximal change to finish iterations
		    final double     grow_disp_min,
			final double     grow_disp_max,
			final double     grow_disp_step,
			final double     unique_pre_tolerance, // usually larger than clt_parameters.unique_tolerance
			final boolean    try_next, // try next disparity range if the current one was already tried
			final int        dbg_x,
			final int        dbg_y,
			final int        debugLevel)
	{
		boolean [] expanding = startExpanding(
				values,           // final double []  values,     // will be modified, Double.NaN is not yet assigned
				known,            // final boolean [] known,      // cells with known values (will not be modified)
				smpl_size,        // final int        smpl_size, // == 5
				true,             // final boolean    no_zero,  // do not use values[i] == 0.0, exact 0.0 is only background			
				dbg_x,            // final int        dbg_x,
				dbg_y,            // final int        dbg_y,
				debugLevel+0);     // final int        debugLevel) // +2 for image
		
		if (expanding == null){
			return null; // no discontinuity seeds
		}
		
		boolean [] prohibited = tp.getTriedBefore( 
				passes,         // final ArrayList <CLTPass3d> passes,
				firstPass,      // final int                   firstPass,
				lastPassPlus1,  // final int                   lastPassPlus1,
				expanding,      // final boolean []            selected,
				null,           // final boolean []            prohibited,
				values,         // final double []             disparity,
				true,           // final boolean               mod_selected,
				true,           // final boolean               mod_disparity,
				false,          // final boolean               en_near,
				false,          // final boolean               en_far,
				grow_disp_min,  // final double                grow_disp_min,
				grow_disp_max,  // final double                grow_disp_max,
				grow_disp_step, // final double                grow_disp_step,
				try_next ? 0.0: unique_pre_tolerance); // final double unique_pre_tolerance) // disable for try_next

		for (int nstep = 0; nstep < num_steps; nstep++){
			int num_new_cells = 0;
			if (nstep > 0){
				num_new_cells = expandDiscontinuityStep(
						values,          // final double []  values,     // will be modified, Double.NaN is not yet assigned
						known,           // final boolean [] known,      // cells with known values (will not be modified)
						prohibited,      // final boolean [] prohibited, // cells that can not be used (tried before with similar disparity)
						expanding,       // final boolean [] expanding,  // cells with known values (will not be modified)
						smpl_size,       //final int        smpl_size, // == 5
						split_threshold, // final double     split_threshold, // if full range of the values around the cell higher, need separate fg, bg
						extend_far,      // final boolean    extend_far,    // extend far (false - extend near)
						true,            // final boolean    unknown_only,  // pure expanding, not over previously found values
						dbg_x,           // final int        dbg_x,
						dbg_y,           // final int        dbg_y,
						debugLevel);     // final int        debugLevel)
				
				tp.getTriedBefore(
						passes,          // final ArrayList <CLTPass3d> passes,
						firstPass,       // final int                   firstPass,
						lastPassPlus1,   // final int                   lastPassPlus1,
						expanding,       // final boolean []            selected,
						prohibited,      //  final boolean []            prohibited,
						values,          // final double []             disparity,
						true,            //  final boolean               mod_selected,
						true,            //  final boolean               mod_disparity,
						false,          // final boolean               en_near,
						false,          // final boolean               en_far,
						grow_disp_min,  // final double                grow_disp_min,
						grow_disp_max,  // final double                grow_disp_max,
						grow_disp_step, // final double                grow_disp_step,
						try_next ? 0.0: unique_pre_tolerance); // final double unique_pre_tolerance) // disable for try_next
			}
			int mt = 2 * (nstep + 1); // max_tries; // or use k * nstep ???
			// increase tries last iteration?
			if (((num_new_cells == 0) && (nstep > 0)) || (nstep == (num_steps -1))) {
				mt = max_tries; // last time smooth more 
			}
			for (int num_try = 0; num_try < mt; num_try ++){
				if (debugLevel >1){
					System.out.print("expandKnown(): nstep="+nstep+" num_try="+num_try+" ");
					if (num_try == (mt - 1)){
						System.out.print("expandKnown(): LAST try ");
					}
				}
				double max_diff = smoothNewStep( // return maximal correction value 
						// should first be initialized with all not known = Double.NaN
						values,          // final double []  values,     // will be modified, Double.NaN is not yet assigned
						expanding,       // final boolean [] new_cells,  // cells that need to  be calculated
						smpl_size,       // final int        smpl_size, // == 5
						min_points,      // final int        min_points, // == 3
						use_wnd,         // final boolean    use_wnd,   // use window function fro the neighbors 
						tilt_cost,       // final double     tilt_cost,
						split_threshold, // final double     split_threshold, // if full range of the values around the cell higher, need separate fg, bg
						same_range,      // final double     same_range,      // modify 
						diff_continue,   // final double     diff_continue,   // maximal difference from the old value (for previously defined tiles
						max_abs_tilt,    // final double     max_abs_tilt, //   = 2.0; // pix per tile
						max_rel_tilt,    // final double     max_rel_tilt, //   = 0.2; // (pix / disparity) per tile
						true,            // final boolean    en_normal, // expand sngle-range cells (that have all neighbors within split_threshold)
						!smooth_only && extend_far,      // final boolean    en_far,    // expand background - use neighbors cells within same_range from the lowest
						// en_near is only valid if !en_far 
						!smooth_only && !extend_far,     // final boolean    en_near,   // expand foreground - use neighbors cells within same_range from the highest
						(num_try == 0),  // new_only,       // first run - set new cells, later - smooth all 
						true,            // final boolean    no_zero,  // do not use values[i] == 0.0, exact 0.0 is only background
						dbg_x,           // final int        dbg_x,
						dbg_y,           // final int        dbg_y,
						debugLevel);     // final int        debugLevel);
				if (debugLevel > 1){
					System.out.println("expandKnown() -> max_diff="+max_diff);
				}

				if (max_diff <= final_diff) break; 
			}
			if (debugLevel > 1){
				System.out.println("expandKnown(): nstep="+nstep+" mt="+mt);
			}
		}
		if (try_next){
			tp.getTriedBefore(
				passes,          // final ArrayList <CLTPass3d> passes,
				firstPass,       // final int                   firstPass,
				lastPassPlus1,   // final int                   lastPassPlus1,
				expanding,       // final boolean []            selected,
				prohibited,      //  final boolean []            prohibited,
				values,          // final double []             disparity,
				true,            //  final boolean               mod_selected,
				true,            //  final boolean               mod_disparity,
				// TODO: some version to enable both en_near and en_far? or make en_far assume en_near and make it first?
				!extend_far,     // final boolean               en_near,
				extend_far,      // final boolean               en_far,
				grow_disp_min,   // final double                grow_disp_min,
				grow_disp_max,   // final double                grow_disp_max,
				grow_disp_step,  // final double                grow_disp_step,
				unique_pre_tolerance); // final double unique_pre_tolerance)
		}		
		
		return expanding;
	}
	
}
