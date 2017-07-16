import java.util.concurrent.atomic.AtomicInteger;

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
			final boolean    use_wnd,   // use window function fro the neighbors 
			final double     tilt_cost,
			final double     split_threshold, // if full range of the values around the cell higher, need separate fg, bg
			final double     same_range,      // modify 
			final boolean    en_normal, // expand sngle-range cells (that have all neighbors within split_threshold)
			final boolean    en_far,    // expand background - use neighbors cells within same_range from the lowest
			// en_near is only valid if !en_far 
			final boolean    en_near,   // expand foreground - use neighbors cells within same_range from the highest
			final int        max_tries, // maximal number of smoothing steps
			final double     final_diff, // maximal change to finish iterations 
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
					smpl_size,       // final int        smpl_size, // == 5
					use_wnd,         // final boolean    use_wnd,   // use window function fro the neighbors 
					tilt_cost,       // final double     tilt_cost,
					split_threshold, // final double     split_threshold, // if full range of the values around the cell higher, need separate fg, bg
					same_range,      // final double     same_range,      // modify 
					en_normal,       // final boolean    en_normal, // expand sngle-range cells (that have all neighbors within split_threshold)
					en_far,          // final boolean    en_far,    // expand background - use neighbors cells within same_range from the lowest
					// en_near is only valid if !en_far 
					en_near,         // final boolean    en_near,   // expand foreground - use neighbors cells within same_range from the highest
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
	
	
	
	
	public double smoothNewStep( // return maximal correction value 
			// should first be initialized with all not known = Double.NaN
			final double []  values,     // will be modified, Double.NaN is not yet assigned
			final boolean [] new_cells,  // cells that need to  be calculated
			final int        smpl_size, // == 5
			final boolean    use_wnd,   // use window function fro the neighbors 
			final double     tilt_cost,
			final double     split_threshold, // if full range of the values around the cell higher, need separate fg, bg
			final double     same_range,      // modify 
			final boolean    en_normal, // expand sngle-range cells (that have all neighbors within split_threshold)
			final boolean    en_far,    // expand background - use neighbors cells within same_range from the lowest
			// en_near is only valid if !en_far 
			final boolean    en_near,   // expand foreground - use neighbors cells within same_range from the highest
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
		final int [] neib_indices = new int [8]; // indices of the sample array immediately around the center;
		int ni = 0;
		for (int sy = -1; sy < 1; sy++) {
			for (int sx = -1; sx < 1; sx++) {
				if ((sy != 0) || (sx != 0)){
					neib_indices[ni++] = (smpl_half + sx) + (smpl_half + sy) * smpl_size; 
				}
			}		
		}		
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				public void run() {
					PolynomialApproximation pa = new PolynomialApproximation();
					int this_thread = ai_thread.getAndIncrement();
					for (int nTile = ai.getAndIncrement(); nTile < new_vals.length; nTile = ai.getAndIncrement()) {
						if (new_cells[nTile]) {
							int dl = (nTile == dbg_tile) ? debugLevel : -1;
							int tileX = nTile % tilesX;
							int tileY = nTile / tilesX;
//							boolean [] smpl_new = new boolean [smpl_len];
							boolean [] smpl_sel = new boolean [smpl_len];
							double [] smpl_d =  new double [smpl_len];
//							double [] smpl_p =  new double [smpl_len];
							double [] smpl_w =  new double [smpl_len];
							int num_in_sample = 0;
							for (int sy = 0; sy < smpl_size; sy++) {
								int y = tileY + sy - smpl_half;
								if ((y >= 0) && (y < tilesY)) {
									for (int sx = 0; sx < smpl_size; sx++) {
										int x = tileX + sx - smpl_half;
										if ((x >= 0) && (x < tilesX)) {
											int indx = x + y * tilesX;
											if (!Double.isNaN(values[indx]) && (indx != center_index)){
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
							if (num_in_sample > 0) {
								// add starting from the known cells. When using 5x5 samples, having immediate neighbor
								// to the center usually means that likely there is enough data to determine the tilt.
								boolean has_neib = false;
								for (int i = 0; i < neib_indices.length; i++) if (smpl_sel[neib_indices[i]]) {
									has_neib = true;
									break;
								}
								boolean process_sample = false;
								if (has_neib) { // all conditions met, find linear 2-d damped approximation for the center
									// if center was NaN - diff is the absolute value.
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

											double threshold = smpl_d[imax] - smpl_d[imax];  
											for (int indxs = 0; indxs < smpl_len; indxs ++ ) if (smpl_sel[indxs]){
												if (smpl_d[indxs] < threshold) {
													smpl_sel[indxs] = false;
													num_in_sample--;
												}
											}										
											process_sample = true;
										} else if (en_near){
											double threshold = smpl_d[imin] + smpl_d[imax];  
											for (int indxs = 0; indxs < smpl_len; indxs ++ ) if (smpl_sel[indxs]){
												if (smpl_d[indxs] > threshold) {
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
											// with other direction
										}
									}
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
												mdata[mindx][2][0] =  smpl_w[indxs];
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
											dl);           // int debugLevel;
									if (approx2d == null) {
										System.out.println("A BUG in smoothNewStep() - failed quadraticApproximation()");
										continue;
									}
									double new_center_val = approx2d[0][0]; // Last coefficient, F (of 6)
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
		
		System.arraycopy(new_vals, 0, values, 0, new_vals.length);
		double max_corr = 0;
		for (int i = 0; i < max_corrs.length; i++) max_corr = Math.max(max_corr,max_corrs[i]);
		return max_corr;
		
	}

}
/*
final Thread[] threads = ImageDtt.newThreadArray(threadsMax);
final AtomicInteger ai = new AtomicInteger(0);
for (int ithread = 0; ithread < threads.length; ithread++) {
	threads[ithread] = new Thread() {
		public void run() {
			for (int nTile = ai.getAndIncrement(); nTile < tiles; nTile = ai.getAndIncrement()) {
//				int dl = (nTile == dbg_tile) ? debugLevel : -1;
				int tileX = nTile % tilesX;
				int tileY = nTile / tilesX;
				measured[nTile] = measured_scan.tile_op[tileY][tileX] != 0;
				if (!measured[nTile]) disp_strength[1][nTile] = -1.0;
			}	
		}
	};
}		      
ImageDtt.startAndJoin(threads);
*/