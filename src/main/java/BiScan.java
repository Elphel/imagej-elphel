/**
 ** BiCamScan - calss to represent bultiple bi-quad camera measurements
 **
 ** Copyright (C) 2018 Elphel, Inc.
 **
 ** -----------------------------------------------------------------------------**
 **
 **  BiCamScan.java is free software: you can redistribute it and/or modify
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
import java.util.concurrent.atomic.AtomicInteger;

public class BiScan {
	final static double THRESHOLD_LIN = 1.0E-20;  // threshold ratio of matrix determinant to norm for linear approximation (det too low - fail)
	final static double THRESHOLD_QUAD = 1.0E-30; // threshold ratio of matrix determinant to norm for quadratic approximation (det too low - fail)

	double []  disparity_measured;
	double []  target_disparity;
	double []  strength_measured;
	boolean [] strong_trusted; // sufficient strength without neighbors
	boolean [] trusted;
	boolean [] cond_trusted;
	boolean [] disabled_measurement; // should disable source also
	int     [] src_index;       // index of the source scan which measured data is used here (applies to disparity_measured, strength_measured, disabled_measurement
	int        list_index = -1;
	BiCamDSI biCamDSI;
//	public BiScan(BiCamDSI biCamDSI) {
//		this.biCamDSI = biCamDSI;
//		int num_tiles = biCamDSI.tnImage.getSizeX()*biCamDSI.tnImage.getSizeY();
//		disparity= new double[num_tiles];
//		strength=  new double[num_tiles];
//		trusted=   new boolean[num_tiles];
//	}
	public BiScan(
			BiCamDSI   biCamDSI,
			int        indx,
			double []  disparity,
			double []  strength,
			boolean [] trusted,
			boolean [] disabled) {
		this.biCamDSI = biCamDSI;
		this.list_index = indx;
		int num_tiles = biCamDSI.tnImage.getSizeX()*biCamDSI.tnImage.getSizeY();
		if (disparity == null) {
			disparity= new double[num_tiles];
		} else {
			this.disparity_measured = disparity.clone();
			if (strength == null)  {
				strength=  new double[num_tiles];
				for (int nTile = 0; nTile < num_tiles; nTile++) strength[nTile] = Double.isNaN(disparity[nTile])?0.0:1.0;
			} else {
				this.strength_measured = strength.clone();
				for (int nTile = 0; nTile < num_tiles; nTile++) {
					if (Double.isNaN(disparity[nTile])) this.strength_measured [nTile] = 0.0;
					if (strength[nTile] == 0.0)         this.disparity_measured[nTile] = Double.NaN;
				}
			}
		}
		if (trusted == null)   trusted=   new boolean[num_tiles];
		if (disabled == null)  disabled=  new boolean[num_tiles];
		this.trusted = trusted;
		this.disabled_measurement = disabled;
		this.strong_trusted = new boolean[num_tiles];
		this.cond_trusted =   new boolean[num_tiles];
		src_index = new int[num_tiles];
		// set new measurement index to this, other to -1
		for (int i = 0; i < num_tiles; i++) {
			src_index[i] = (strength[i] > 0.0)? list_index:-1;
		}
	}

	public double []  getDisparityMeasured()   { return this.disparity_measured;} // FIXME!
	public double []  getStrengthMeasured()    { return this.strength_measured;} // FIXME
	public boolean [] getTrusted()             { return this.trusted;}
	public boolean [] getDisabledMeasurement() { return this.disabled_measurement;}

	public void disableTile (int nTile) {
		trusted[nTile] =        false;
		strong_trusted[nTile] = false;
		cond_trusted[nTile] =   false;
		//    	disabled[nTile] =       true;
		//    	if ((src_index[nTile] >= 0) && (src_index[nTile] != list_index)) {
		if (src_index[nTile] >= 0) {
			biCamDSI.getBiScan(src_index[nTile]).disabled_measurement[nTile] = false; // may be source tile or this tile
		}
	}
    /**
     * Get disparity and strength from the scan, mask by boolean attributes
     * @param only_strong keep only trusted strong tiles
     * @param only_trusted keep any trusted tiles
     * @param only_enabled keep all but disabled tiles
     * @return array of two arrays {disparity, strength}
     */

    public double [][] getDisparityStrength( // FIXME
    		final boolean only_strong,
    		final boolean only_trusted,
    		final boolean only_enabled){
    	final int num_tiles = disparity_measured.length;
    	final double [][] ds = {disparity_measured.clone(), strength_measured.clone()}; // just to start with
		final Thread[] threads = ImageDtt.newThreadArray(biCamDSI.threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		final BiScan this_scan = this;
		// find definitely trusted and conditionally trusted tiles
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				@Override
				public void run() {
					for (int nTile = ai.getAndIncrement(); nTile < num_tiles; nTile = ai.getAndIncrement()) {
						if ((only_strong && !strong_trusted[nTile]) ||
								(only_trusted && !trusted[nTile])) {
							ds[0][nTile] = Double.NaN;
							ds[1][nTile] = 0.0;
						} else { //if (src_index[nTile] != list_index){ // only one level of indirection?
							int src = src_index[nTile]; // same tile or different
							BiScan scan = this_scan;
							if (src < 0)  {
								src = list_index;
							} else {
								scan = biCamDSI.biScans.get(src);
							}
							boolean dsbl = scan.disabled_measurement[nTile];
							if (dsbl && only_enabled) { // src <0 only for the first scan where no data is available
								ds[0][nTile] = Double.NaN;
								ds[1][nTile] = 0.0;
							} else {
								ds[0][nTile] = scan.disparity_measured[nTile];
								ds[1][nTile] = scan.strength_measured[nTile];
							} // if (dsbl && only_enabled) -- else
						}
					}
				}
			};
		}
		ImageDtt.startAndJoin(threads);
    	return ds;
    }

	public void showScan(String title) { // FIXME!
		String [] titles= {
				"all",               // 0
				"enabled",           // 1
				"cond_trusted",      // 2
				"weak trusted",      // 3
				"strong trusted",    // 4
				"measured",          // 5
				"suggested",         // 6
				"measured strength", // 7
				"strength" };        // 8
		double [][] ds_all = getDisparityStrength( // FIXME
	    		false,   // final boolean only_strong,
	    		false,   // final boolean only_trusted,
	    		false) ; // final boolean only_enabled);
		double [][] ds = getDisparityStrength( // FIXME
	    		false,   // final boolean only_strong,
	    		false,   // final boolean only_trusted,
	    		true) ; // final boolean only_enabled);

		double [][] dbg_data = new double[titles.length][];
		dbg_data[5] = this.disparity_measured;
		dbg_data[0] = ds_all[0];
		dbg_data[7] = this.strength_measured;
		dbg_data[8] = ds_all[1];
		dbg_data[1] = ds[0];

		if (this.cond_trusted != null) {
			dbg_data[2] = ds[0].clone();
			for (int i = 0; i < this.cond_trusted.length; i++) if (!cond_trusted[i]) dbg_data[2][i] = Double.NaN;
		}
		if (this.trusted != null) {
			dbg_data[3] = ds[0].clone();
			for (int i = 0; i < this.trusted.length; i++) if (!trusted[i]) dbg_data[3][i] = Double.NaN;
		}
		if (this.strong_trusted != null) {
			dbg_data[4] =  ds[0].clone();
			for (int i = 0; i < this.strong_trusted.length; i++) if (!strong_trusted[i]) dbg_data[4][i] = Double.NaN;
		}
		if (this.target_disparity != null) {
			dbg_data[6] = this.target_disparity.clone();
		}
		(new showDoubleFloatArrays()).showArrays(
				dbg_data,
				biCamDSI.tnImage.getSizeX(),
				biCamDSI.tnImage.getSizeY(),
				true,
				title,
				titles);
	}

	/**
	 * Reduce outliers on DSI when multiple "refined" disparity values exist for the same tile and the strongest does not seem to be the best
	 * Each disparity solution is compared to the weighted average of the neighbors and the strength is divided by the difference from that
	 * average value, so the closest to the neighbors gets strength boost.
	 * @param str_floor absolute strength to subtract from the measured
	 * @param pf_disp_afloor offset disparity to add to the disparity difference to avoid division by 0 or small numbers
	 * @param pf_disp_rfloor realtive to the disparity portion of the offset
	 * @return number of replaced tiles
	 */
	public int copyFittestEnabled(
			final double  str_floor,      // absolute strength floor
			final double  pf_disp_afloor, // =            0.1;    // When selecting the best fit from the alternative disparities, divide by difference increased by this
			final double  pf_disp_rfloor) //  =            0.02;   // Increase pf_disp_afloor for large disparities
		{
		final TileNeibs  tnImage = biCamDSI.tnImage;
		final int num_tiles = biCamDSI.tnImage.getSizeX()*biCamDSI.tnImage.getSizeY();

		final double [][] ds = getDisparityStrength( // FIXME
	    		false,   // final boolean only_strong,
	    		false,   // final boolean only_trusted,
	    		true) ; // final boolean only_enabled);

		final Thread[] threads = ImageDtt.newThreadArray(biCamDSI.threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		final int [] new_src = new int[num_tiles];
		final AtomicInteger num_changes = new AtomicInteger(0);
		int dbg_x = 157;
		int dbg_y = 212;
		int debugLevel = -1;
		final int dbg_tile = (debugLevel>-2)?(dbg_x + tnImage.sizeX*dbg_y):-1;

		ai.set(0);
		// find definitely trusted and conditionally trusted tiles
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				@Override
				public void run() {
					for (int nTile = ai.getAndIncrement(); nTile < num_tiles; nTile = ai.getAndIncrement()) {
						if (nTile == dbg_tile) {
							System.out.println("copyFittestEnabled(): nTile="+nTile);
						}
						if (ds[1][nTile] > 0.0 ) { // && (src_index[nTile] != list_index)){ // should be already found
							int num_neib = 0;
							double sw = 0.0, swd = 0.0;
							double wdiag = 0.7;
							for (int dir = 0; dir <8; dir++) {
								int nTile1 = tnImage.getNeibIndex(nTile, dir);
								if ((nTile1 >= 0) && (ds[1][nTile1] > str_floor)) {
									double w = ds[1][nTile1] - str_floor;
									if ((dir & 1) != 0) w *= wdiag;
									sw += w;
									swd += ds[0][nTile1] * w;
									num_neib++;
								}
							}
							if (num_neib > 0) {
								double disp_mean = swd/sw;
								double disp_floor = pf_disp_afloor + disp_mean * pf_disp_rfloor;
								double disp_floor2 = disp_floor*disp_floor;

								int best_indx=-1;
								double best_strength =0.0;
								for (int indx = list_index; indx >= 0; indx--) { // include the latest measurement
									BiScan scan = biCamDSI.getBiScan(indx);
									if (scan.disabled_measurement[nTile]) { //  || (scan.src_index[nTile] != indx)){ // skip all but enabled
										continue;
									}
									// calculate effective strength
									double strength = scan.strength_measured[nTile] - str_floor;
									if (strength <= 0) {
										continue;
									}

									double diff = scan.disparity_measured[nTile] - disp_mean;
									double eff_strength = strength/Math.sqrt(diff*diff + disp_floor2);

									if (eff_strength > best_strength) {
										best_strength = eff_strength;
										best_indx = indx;
									}
								}
								if ((best_indx >= 0) &&  (best_indx != src_index[nTile])) { // not the same as already set
									new_src[nTile] = best_indx+1; // +1 so initial 0 will be "not set"
									num_changes.getAndIncrement();
								}
							}
						}
					}
				}
			};
		}
		ImageDtt.startAndJoin(threads);
		if (num_changes.get() > 0) {
			ai.set(0);
			// find definitely trusted and conditionally trusted tiles
			for (int ithread = 0; ithread < threads.length; ithread++) {
				threads[ithread] = new Thread() {
					@Override
					public void run() {
						for (int nTile = ai.getAndIncrement(); nTile < num_tiles; nTile = ai.getAndIncrement()) {
							if (nTile == dbg_tile) {
								System.out.println("copyFittestEnabled() 2 : nTile="+nTile);
							}
							if (new_src[nTile] > 0){
								int best_indx = new_src[nTile]-1;
								src_index[nTile] = best_indx;
							}
						}
					}
				};
			}
			ImageDtt.startAndJoin(threads);
		}
		// will need trusted* recalculated
		return num_changes.get();
	}

	/**
	 * Copy data (if the current was not measured) from one of the previous scans - strongest that is not disabled. If last_priority is true
	 * the latest not disabled scan will be used, even if it is not the strongest
	 * @param last_priority
	 */

	public void copyLastStrongestEnabled(
			final boolean last_priority) // use last if exists and is not disabled
		{
		final int num_tiles = biCamDSI.tnImage.getSizeX()*biCamDSI.tnImage.getSizeY();
		final Thread[] threads = ImageDtt.newThreadArray(biCamDSI.threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
//		final boolean [] cond_trusted = new boolean [num_tiles];
		strong_trusted = new boolean [num_tiles];
		trusted =        new boolean [num_tiles];
		cond_trusted =   new boolean [num_tiles];

		ai.set(0);
		// find definitely trusted and conditionally trusted tiles
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				@Override
				public void run() {
					for (int nTile = ai.getAndIncrement(); nTile < num_tiles; nTile = ai.getAndIncrement()) if (strength_measured[nTile] <= 0.0 ){ // keep last measurement for a while, even if it is not the best
						int best_indx=-1;
						boolean latest = true;
						double best_strength =0.0;
						for (int indx = list_index-1; indx >= 0; indx--) { // no need to try the latest - it is empty (strength_measured[nTile] <= 0.0 )
							BiScan scan = biCamDSI.getBiScan(indx);
							if (scan.disabled_measurement[nTile] || (scan.src_index[nTile] != indx)){ // skip all but enabled sources
								continue;
							}
							if (scan.strength_measured[nTile] > best_strength) {
								best_strength = scan.strength_measured[nTile];
								best_indx = indx;
							}
							if (last_priority && latest) {
								break; // best_indx should be set correctly, as strength > 0.0
							}
							latest = false; // first not disabled with strength>0 gets here
						}
						if (best_indx >= 0) {
							src_index[nTile] = best_indx;
						}
					}
				}
			};
		}
		ImageDtt.startAndJoin(threads);
	}




	public void showScans(String title) {
	}
	/**
	 * Suggest disparities to try for the tiles in poorly textured areas by fitting planes in DSI
	 * calcTrusted should be called before to set up trusted/cond_trusted tiles
	 * suggested tiles will be compared against and made sure they differ by more than a specified margin
	 * 1) current measured (refined) disparity value
	 * 2) target disaprity that lead to the current measurement after refinement
	 * 3) any other disable measurement
	 * 4) any target disparity that lead to the disabled measurement
	 * @param trusted_strength strength to trust unconditionally
	 * @param strength_rfloor strength floor to subrtact as a fraction of the trusted strength
	 * @param discard_cond if true may suggest new disparities for conditionally trusted tiles
	 * @param discard_weak if true may suggest new disparities over trusted weak tiles
	 * @param discard_stron if true may suggest new disparities over any tile
	 * @param strength_pow raise strength to thyis power (normally just 1.0)
	 * @param smpl_radius sample "radius", square side is  2 * smpl_radius + 1
	 * @param smpl_num minimal absolute number of samples required to try fit a plane and validate a tile
	 * @param smpl_fract minimal fraction number of the neighbor samples that fit the rms filter required to try fit a plane and validate a tile
	 * @param smpl_num_narrow minimal absolute number of samples for preliminary fitting plane to trhe center area
	 * @param max_adiff maximal absolute difference from the center tile for initial neighbors selection
	 * @param max_rdiff maximal (additional) relative (to tile disparity) difference from the center tile for initial neighbors selection
	 * @param max_atilt maximal absolute tilt (pix/tile) for the tilted planes to fit
	 * @param max_rtilt maximal relative tilt (pix/tile per disparity pixel). min(max_rtilt*disparity, max_atilt) will be used
	 * @param smpl_arms maximal absolute rms of the weighted remaining samples for the successful plane fitting
	 * @param smpl_rrms maximal relative (additional)rms of the weighted remaining samples for the successful plane fitting
	 * @param damp_tilt regularization value to handle planes if the remaining samples are co-linear (or just a single tile)
	 * @param rwsigma weight Gaussian sigma to reduce influence of far tiles relative to smpl_radius
	 * @param rwsigma_narrow Gaussian sigma for the preliminary plain fitting using the closesttiles ~= 1/smpl_radius
	 * @param new_diff minimal difference between the new suggested and the already tried/measured one
	 * @param remove_all_tried remove from suggested - not only disabled, but all tried
	 * @param dbg_x tileX to debug
	 * @param dbg_y tileY to debug
	 * @param debugLevel debug level
	 * @return number of new tiles to measure in the  array of suggested disparities - Double.NaN - nothing suggested
	 *  for the tile. May need additional filtering to avoid suggested already tried disparities
	 */

	int  suggestNewScan(
		    final double     trusted_strength, // trusted correlation strength
			final double     strength_rfloor,   // strength floor - relative to trusted
		    final boolean    discard_cond,      // consider conditionally trusted tiles (not promoted to trusted) as empty
		    final boolean    discard_weak,      // consider weak trusted tiles (not promoted to trusted) as empty
		    final boolean    discard_strong,    // suggest new disparitieas even for strong tiles
			final double     strength_pow,      // raise strength-floor to this power
			final int        smpl_radius,
			final int        smpl_num,   //         = 3;      // Number after removing worst (should be >1)
			final double     smpl_fract, // Number of friends among all neighbors
			final int        smpl_num_narrow,   //         = 3;      // Number after removing worst (should be >1)
			final double     max_adiff,  // Maximal absolute difference betweenthe center tile and friends
			final double     max_rdiff, //  Maximal relative difference between the center tile and friends
			final double     max_atilt, //  = 2.0; // pix per tile
			final double     max_rtilt, //  = 0.2; // (pix / disparity) per tile
			final double     smpl_arms, //         = 0.1;    // Maximal RMS of the remaining tiles in a sample
			final double     smpl_rrms,        //      = 0.005;  // Maximal RMS/disparity in addition to smplRms
			final double     damp_tilt, //   =     0.001; // Tilt cost for damping insufficient plane data
			final double     rwsigma,           //  = 0.7; // influence of far neighbors diminish as a Gaussian with this sigma
			final double     rwsigma_narrow,    //  = used to determine initial tilt
			final double     new_diff,          // minimal difference between the new suggested and the already tried/measured one
			final boolean    remove_all_tried,  // remove from suggested - not only disabled, but all tried
			final int        dbg_x,
			final int        dbg_y,
			final int        debugLevel
			) {
		final TileNeibs  tnImage = biCamDSI.tnImage;
		final double     wsigma = rwsigma*smpl_radius;
		final double     wsigma_narrow = rwsigma_narrow*smpl_radius;
		// prepare window
		final double [][] weights =        new double [smpl_radius + 1][smpl_radius + 1];
		final double [][] weights_narrow = new double [smpl_radius + 1][smpl_radius + 1];
		for (int i = 0; i <weights.length; i++) {
			for (int j = i; j <weights[i].length; j++) {
				weights[i][j] = (wsigma >0.0) ?Math.exp(-(i*i+j*j)/(2*wsigma*wsigma)):1.0;
				weights[j][i] = weights[i][j];
				weights_narrow[i][j] = (wsigma_narrow >0.0) ?Math.exp(-(i*i+j*j)/(2*wsigma_narrow*wsigma_narrow)):1.0;
				weights_narrow[j][i] = weights_narrow[i][j];
			}
		}
		final double [][] ds = getDisparityStrength( // FIXME
	    		false,   // final boolean only_strong,
	    		false,   // final boolean only_trusted,
	    		true) ; // final boolean only_enabled);
		final int num_tiles = tnImage.getSizeX()*tnImage.getSizeY();
		final double     strength_floor = trusted_strength * strength_rfloor;
		final int        smpl_side = 2 * smpl_radius + 1; // Sample size (side of a square)
		final int        smpl_len = smpl_side*smpl_side;

		final boolean [] trusted_sw = discard_weak ? (this.strong_trusted) : (discard_cond ? this.trusted: this.cond_trusted);

		final Thread[] threads = ImageDtt.newThreadArray(biCamDSI.threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
//		final double [] suggested_disparity = new double [num_tiles];
		target_disparity = new double [num_tiles];
		for (int i = 0; i < num_tiles; i++) target_disparity[i] = Double.NaN;
//		cond_trusted and trusted should be set;
		ai.set(0);
//		final BiScan this_scan = this;
		final AtomicInteger num_new = new AtomicInteger(0);
		// find definitely trusted and conditionally trusted tiles
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				@Override
				public void run() {
					int smpl_center = (smpl_side + 1) * smpl_radius;
					for (int nTile = ai.getAndIncrement(); nTile < num_tiles; nTile = ai.getAndIncrement()) if (discard_strong || !trusted_sw[nTile]){
						// Select all neighbors, then filter
						double [] smpl_d =           new double  [smpl_len];
						double [] smpl_w =           new double  [smpl_len];
						double [] smpl_w_narrow =    new double  [smpl_len];
						double [] smpl_p =           new double  [smpl_len]; // plane disparity,
						int nall = 0;
						double sw = 0, swd = 0;
						for (int dy = -smpl_radius; dy <= smpl_radius; dy++) {
							int ady = (dy > 0)? dy:(-dy);
							for (int dx = -smpl_radius; dx <= smpl_radius; dx++) if ((dx != 0) || (dy != 0)){
								int nTile1 =  tnImage.getNeibIndex(nTile, dx, dy);
								if ((nTile1 >= 0) && trusted[nTile1]) { // weak trusted OK
									nall++;
									int adx = (dx > 0)? dx:(-dx);
									int smpl_indx = smpl_center + dy*smpl_side + dx;
									double w = ds[1][nTile1]-strength_floor;
									smpl_d[smpl_indx] =       ds[0][nTile1];
									smpl_w[smpl_indx] =        w * weights[ady][adx];
									smpl_w_narrow[smpl_indx] = w * weights_narrow[ady][adx];
									if (strength_pow != 1.0) {
										smpl_w[smpl_indx] = Math.pow(smpl_w[smpl_indx], strength_pow);
										smpl_w_narrow[smpl_indx] = Math.pow(smpl_w_narrow[smpl_indx], strength_pow);
									}
									sw += smpl_w_narrow[smpl_indx];
									swd += smpl_w_narrow[smpl_indx]* smpl_d[smpl_indx];
								}
							}
						}

						if (sw == 0.0) {
							continue; //
						}
						double disp_mean = swd/sw; // preliminary reference disparity
						double max_tilt = max_rtilt * disp_mean;
						if (max_tilt > max_atilt) {
							max_tilt = max_atilt;
						}

						int fin_samples= (int) ( nall * smpl_fract);
						if (fin_samples < smpl_num) fin_samples = smpl_num;
//						int fin_samples_narrow= (int) ( nall * smpl_fract_narrow);
//						if (fin_samples_narrow < smpl_num) fin_samples_narrow = smpl_num_narrow;

						// fit plane to mostly centertiles
						int nsmpls = nall;
						if (nsmpls < smpl_num_narrow) { // no tiles even to start
							continue; //
						}
						double [] fit_rslt = fitPlaneRemoveOutliers(
								smpl_radius,        // int                     smpl_radius,
								max_tilt,           // double                  max_tilt,
								damp_tilt,          // double                  damp_tilt, //   =     0.001; // Tilt cost for damping insufficient plane data
								true,               // boolean                 full_plane,
								smpl_d,             // double []               smpl_d,
								smpl_w_narrow,      // double []               smpl_w, // will be modified,
								smpl_p,             // double []               smpl_p, // will be set if provided
								smpl_num_narrow,    // int                     fin_samples, // remove until this number remain
								debugLevel);        // int                     debugLevel)
						if ( (fit_rslt == null) ||  (fit_rslt[0] > (smpl_arms + smpl_rrms * fit_rslt[1]))){
							continue; // narrow selection - too high rms
						}
						disp_mean = fit_rslt[1]; // smpl_p[smpl_center]; // center of the fitted plane
						// re-select tiles to fit the plane and use wide weights
						double max_diff = max_adiff + max_rdiff * disp_mean; // no provisions for tilt (or add a fraction)?
						nsmpls = 0;
						for (int indxs = 0; indxs < smpl_len; indxs++) if (smpl_w[indxs]>0){
							if (Math.abs(smpl_d[indxs] - smpl_p[indxs]) < max_diff) {
								nsmpls++;
							} else {
								smpl_w[indxs] = 0.0;
							}
						}
						if (nsmpls < fin_samples) { // no tiles even to satrt
							continue; //
						}
						fit_rslt = fitPlaneRemoveOutliers(
								smpl_radius,        // int                     smpl_radius,
								max_tilt,           // double                  max_tilt,
								damp_tilt,          // double                  damp_tilt, //   =     0.001; // Tilt cost for damping insufficient plane data
								true,               // boolean                 full_plane,
								smpl_d,             // double []               smpl_d,
								smpl_w,             // double []               smpl_w, // will be modified,
								smpl_p,             // double []               smpl_p, // will be set if provided
								fin_samples,        // int                     fin_samples, // remove until this number remain
								debugLevel);        // int                     debugLevel)
						if ( (fit_rslt == null) ||  (fit_rslt[0] > (smpl_arms + smpl_rrms * fit_rslt[1]))){
							continue; // narrow selection - too high rms
						}
						if (Math.abs(fit_rslt[1] - ds[0][nTile]) < new_diff) { // suggested is too close to already measured
							continue; // already measured for this tile
						}
						// compare to the previous suggestion
						int previous_indx = ((src_index[nTile] < 0)? list_index:src_index[nTile]) -1;
						double previous_target = (previous_indx >= 0)? biCamDSI.getBiScan(previous_indx).target_disparity[nTile]:Double.NaN; // Nothing is known about the target of the 0-scan
						if (Math.abs(fit_rslt[1] - previous_target) < new_diff) { // suggested is too close to already suggested and result - disabled
							continue; // already measured for this tile
						}
						// see if close one was already disabled
						// do not compare with the scans that were not disabled - re-try them?remove_all_tried
						boolean valid_suggestion = true;
						for (BiScan other_scan:biCamDSI.biScans) if (other_scan.disabled_measurement[nTile] || remove_all_tried) {
//							int other_indx = (other_scan.src_index[nTile] < 0)? other_scan.list_index:other_scan.src_index[nTile];
//							double other_disparity = biCamDSI.getBiScan(other_indx).disparity_measured[nTile];
							double other_disparity = other_scan.disparity_measured[nTile];
							if (Math.abs(fit_rslt[1] - other_disparity) < new_diff) { // suggested is too close to already measured and disabled
								valid_suggestion = false;
								break; // already measured for this tile
							}
							int other_indx = other_scan.list_index;
							double other_target = (other_indx > 0)? biCamDSI.getBiScan(other_indx - 1).target_disparity[nTile]:Double.NaN; // Nothing is known about the target of the 0-scan
							if (Math.abs(fit_rslt[1] - other_target) < new_diff) { // suggested is too close to already suggested and result - disabled
								valid_suggestion = false;
								break; // already measured for this tile
							}
						}
						if (valid_suggestion) {
							target_disparity[nTile] = fit_rslt[1];
							num_new.getAndIncrement();
						}
					}
				}
			};
		}
		ImageDtt.startAndJoin(threads);
		// remove duplicates from what was suggested or measured before
		return num_new.get();
	}


// Simple version for non-flat strong areas - try duplicating neighbor
	int  suggestNewScan(
			final int []     dxy,               //up,down,right,left
		    final boolean    discard_cond,      // consider conditionally trusted tiles (not promoted to trusted) as empty
		    final boolean    discard_weak,      // consider weak trusted tiles (not promoted to trusted) as empty
		    final boolean    discard_strong,    // suggest new disparitieas even for strong tiles
			final double     new_diff,          // minimal difference between the new suggested and the already tried/measured one
			final boolean    remove_all_tried,  // remove from suggested - not only disabled, but all tried
			final int        dbg_x,
			final int        dbg_y,
			final int        debugLevel
			) {
		final TileNeibs  tnImage = biCamDSI.tnImage;
		final double [][] ds = getDisparityStrength( // FIXME
	    		false,   // final boolean only_strong,
	    		false,   // final boolean only_trusted,
	    		true) ;  // final boolean only_enabled);

		final int num_tiles = tnImage.getSizeX()*tnImage.getSizeY();

		final boolean [] trusted_sw = discard_weak ? (this.strong_trusted) : (discard_cond ? this.trusted: this.cond_trusted);
		final boolean [] trusted_weak = this.cond_trusted; // or even any?

		final Thread[] threads = ImageDtt.newThreadArray(biCamDSI.threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		target_disparity = new double [num_tiles];
		for (int i = 0; i < num_tiles; i++) target_disparity[i] = Double.NaN;
		ai.set(0);
		final AtomicInteger num_new = new AtomicInteger(0);
		// find definitely trusted and conditionally trusted tiles
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				@Override
				public void run() {
					for (int nTile = ai.getAndIncrement(); nTile < num_tiles; nTile = ai.getAndIncrement()) if (discard_strong || !trusted_sw[nTile]){
						int nTile1 =  tnImage.getNeibIndex(nTile, dxy[0], dxy[1]);
						if ((nTile1 >= 0) && trusted_weak[nTile1]) { // weak trusted OK, maybe even any measured
							double new_disp = ds[0][nTile1];
							if (Math.abs(new_disp - ds[0][nTile]) < new_diff) { // suggested is too close to already measured
								continue; // already measured for this tile
							}
							// compare to the previous suggestion
							int previous_indx = ((src_index[nTile] < 0)? list_index:src_index[nTile]) -1;
							double previous_target = (previous_indx >= 0)? biCamDSI.getBiScan(previous_indx).target_disparity[nTile]:Double.NaN; // Nothing is known about the target of the 0-scan
							if (Math.abs(new_disp - previous_target) < new_diff) { // suggested is too close to already suggested and result - disabled
								continue; // already measured for this tile
							}
							// see if close one was already disabled
							boolean valid_suggestion = true;
							for (BiScan other_scan:biCamDSI.biScans) if (other_scan.disabled_measurement[nTile] || remove_all_tried) {
//								int other_indx = (other_scan.src_index[nTile] < 0)? other_scan.list_index:other_scan.src_index[nTile];
//								double other_disparity = biCamDSI.getBiScan(other_indx).disparity_measured[nTile];
								double other_disparity = other_scan.disparity_measured[nTile];
								if (Math.abs(new_disp - other_disparity) < new_diff) { // suggested is too close to already measured and disabled
									valid_suggestion = false;
									break; // already measured for this tile
								}
								int other_indx = other_scan.list_index;
								double other_target = (other_indx > 0)? biCamDSI.getBiScan(other_indx - 1).target_disparity[nTile]:Double.NaN; // Nothing is known about the target of the 0-scan
								if (Math.abs(new_disp - other_target) < new_diff) { // suggested is too close to already suggested and result - disabled
									valid_suggestion = false;
									break; // already measured for this tile
								}
							}
							if (valid_suggestion) {
								target_disparity[nTile] = new_disp;
								num_new.getAndIncrement();
							}
						}
					}
				}
			};
		}
		ImageDtt.startAndJoin(threads);
		return num_new.get();
	}






	/**
	 * Calculate trusted tiles from the strength and disparity. Trusted are tiles that are really strong
	 * (then they are trusted regardless of neighbors) or are somewhat strong and have sufficient neighbors
	 * that (together with this tile) make a good (tilted) plane
	 * @param trusted_strength strength to trust unconditionally
	 * @param strength_rfloor strength floor to subrtact as a fraction of the trusted strength
	 * @param cond_rtrusted fraction of the trusted strength (after subtracting str4ength_floor) that is sufficient
	 * to participate in plane fitting, if successful - make a tile trusted
	 * @param strength_pow raise strength to thyis power (normally just 1.0)
	 * @param smpl_radius sample "radius", square side is  2 * smpl_radius + 1
	 * @param smpl_num minimal absolute number of samples required to try fit a plane and validate a tile
	 * @param smpl_fract minimal fraction number of the neighbor samples that fit the rms filter required to try fit a plane and validate a tile
	 * @param max_adiff maximal absolute difference from the center tile for initial neighbors selection
	 * @param max_rdiff maximal (additional) relative (to tile disparity) difference from the center tile for initial neighbors selection
	 * @param max_atilt maximal absolute tilt (pix/tile) for the tilted planes to fit
	 * @param max_rtilt maximal relative tilt (pix/tile per disparity pixel). min(max_rtilt*disparity, max_atilt) will be used
	 * @param smpl_arms maximal absolute rms of the weighted remaining samples for the successful plane fitting
	 * @param smpl_rrms maximal relative (additional)rms of the weighted remaining samples for the successful plane fitting
	 * @param damp_tilt regularization value to handle planes if the remaining samples are co-linear (or just a single tile)
	 * @param rwsigma weight Gaussina sigma to reduce influence of far tiles relative to smpl_radius
	 * @param dbg_x tileX to debug
	 * @param dbg_y tileY to debug
	 * @param debugLevel debug level
	 * @return array of 3 numers: number of trusted strong tiles, number of additional trusted by plane fitting, and number of all
	 * somewhat strong tiles
	 */


	int [] calcTrusted(
		    final double     trusted_strength, // trusted correlation strength
			final double     strength_rfloor,   // strength floor - relative to trusted
		    final double     cond_rtrusted,     // minimal strength to consider - fraction of trusted
			final double     strength_pow,      // raise strength-floor to this power
			final int        smpl_radius,
			final int        smpl_num,   //         = 3;      // Number after removing worst (should be >1)
			final double     smpl_fract, // Number of friends among all neighbors
			final double     max_adiff,  // Maximal absolute difference between the center tile and friends
			final double     max_rdiff, //  Maximal relative difference between the center tile and friends
			final double     max_atilt, //  = 2.0; // pix per tile
			final double     max_rtilt, //  = 0.2; // (pix / disparity) per tile
			final double     smpl_arms, //         = 0.1;    // Maximal RMS of the remaining tiles in a sample
			final double     smpl_rrms,        //      = 0.005;  // Maximal RMS/disparity in addition to smplRms
			final double     damp_tilt, //   =     0.001; // Tilt cost for damping insufficient plane data
			final double     rwsigma,           //  = 0.7; // influence of far neighbors diminish as a Gaussian with this sigma
			final int        dbg_x,
			final int        dbg_y,
			final int        debugLevel

			) {
		final TileNeibs  tnImage = biCamDSI.tnImage;
		final double     wsigma = rwsigma*smpl_radius;
		final double [][] ds = getDisparityStrength( // already has disabled zeroed
	    		false,   // final boolean only_strong,
	    		false,   // final boolean only_trusted,
	    		true) ;  // final boolean only_enabled);
		final int num_tiles = tnImage.getSizeX()*tnImage.getSizeY();
		final double     strength_floor = trusted_strength * strength_rfloor;
		final double     min_strength = strength_floor + (trusted_strength - strength_floor) *  cond_rtrusted;
		final int        smpl_side = 2 * smpl_radius + 1; // Sample size (side of a square)
		final int        smpl_len = smpl_side*smpl_side;
		// prepare window
		final double [][] weights = new double [smpl_radius + 1][smpl_radius + 1];
		for (int i = 0; i <weights.length; i++) {
			for (int j = i; j <weights[i].length; j++) {
				weights[i][j] = (wsigma >0.0) ?Math.exp(-(i*i+j*j)/(2*wsigma*wsigma)):1.0;
				weights[j][i] = weights[i][j];
			}
		}
		final Thread[] threads = ImageDtt.newThreadArray(biCamDSI.threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		final AtomicInteger num_non_weak = new AtomicInteger(0);
		final AtomicInteger num_trusted_strong = new AtomicInteger(0);
		final AtomicInteger num_trusted_plane = new AtomicInteger(0);
//		final boolean [] cond_trusted = new boolean [num_tiles];
		strong_trusted = new boolean [num_tiles];
		trusted =        new boolean [num_tiles];
		cond_trusted =   new boolean [num_tiles];

		ai.set(0);
		// find definitely trusted and conditionally trusted tiles
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				@Override
				public void run() {
					for (int nTile = ai.getAndIncrement(); nTile < num_tiles; nTile = ai.getAndIncrement()){
						if (!Double.isNaN(ds[0][nTile]) && (ds[1][nTile] >= min_strength) || (ds[1][nTile] > 0)) {
							num_non_weak.getAndIncrement();
							cond_trusted[nTile] = true;
							if (ds[1][nTile] > trusted_strength) {
								strong_trusted[nTile] = true;
								trusted[nTile] =        true;
								num_trusted_strong.getAndIncrement();
							}
						}
					}
				}
			};
		}
		ImageDtt.startAndJoin(threads);

		ai.set(0);
		// find definitely trusted and conditionally trusted tiles
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				@Override
				public void run() {
					int smpl_center = (smpl_side + 1) * smpl_radius;
					for (int nTile = ai.getAndIncrement(); nTile < num_tiles; nTile = ai.getAndIncrement()) if (cond_trusted[nTile] && !trusted[nTile]){
						double [] smpl_d =    new double  [smpl_len];
						double [] smpl_w =    new double  [smpl_len];
						double [] smpl_p =    new double  [smpl_len]; // plane disparity,
						// copy neighbor tiles
						double disp_center = ds[0][nTile];
						double max_diff = max_adiff + max_rdiff * disp_center;
						double max_tilt = max_rtilt * disp_center;
						if (max_tilt > max_atilt) {
							max_tilt = max_atilt;
						}

						int nsmpls = 0;
						int nall = 0;

						for (int dy = -smpl_radius; dy <= smpl_radius; dy++) {
							int ady = (dy > 0)? dy:(-dy);
							for (int dx = -smpl_radius; dx <= smpl_radius; dx++) {
								int nTile1 =  tnImage.getNeibIndex(nTile, dx, dy);
								if ((nTile1 >= 0) && cond_trusted[nTile1]) {
									nall++;
									int adx = (dx > 0)? dx:(-dx);
									double max_fdiff = max_diff + (ady+adx) * max_tilt;
									int smpl_indx = smpl_center + dy*smpl_side + dx;
									if (Math.abs(ds[0][nTile1] - disp_center) <= max_fdiff) {
//										smpl_sel[smpl_indx] = true;
										smpl_d[smpl_indx] = ds[0][nTile1];
										smpl_w[smpl_indx] = (ds[1][nTile1]-strength_floor) * weights[ady][adx];
										if (strength_pow != 1.0) {
											smpl_w[smpl_indx] = Math.pow(smpl_w[smpl_indx], strength_pow);
										}
										nsmpls ++;
									}
								}
							}
						}
						int fin_samples= (int) ( nall * smpl_fract);
						if (fin_samples < smpl_num) fin_samples = smpl_num;

						if (nsmpls >= fin_samples) {
							double [] fit_rslt = fitPlaneRemoveOutliers(
									smpl_radius, // int                     smpl_radius,
									max_tilt,    // double                  max_tilt,
									damp_tilt,   // double                  damp_tilt, //   =     0.001; // Tilt cost for damping insufficient plane data
									false,       // boolean                 full_plane,
									smpl_d,      // double []               smpl_d,
									smpl_w,      // double []               smpl_w, // will be modified,
									smpl_p,      // double []               smpl_p, // will be set if provided
									fin_samples, // int                     fin_samples, // remove until this number remain
									debugLevel); // int                     debugLevel)
							if ( (fit_rslt != null) && (fit_rslt[0] < (smpl_arms + smpl_rrms * fit_rslt[1]))){
								// Trusted tile, save it
								trusted[nTile] = true;
								num_trusted_plane.getAndIncrement();
							}
						}
					}
				}
			};
		}
		ImageDtt.startAndJoin(threads);
		int [] numTrustedAll = {num_trusted_strong.get(), num_trusted_plane.get(), num_non_weak.get()};
		return numTrustedAll;
	}

	// FG edge should be strong
	// Trimming(disabling) weak (trusted but not strong_trusted) tiles if on any one side:
	// a) there are no same plane or closer tiles
	// b) there are no in-plane or closer strong tiles, but there are some (strong or any?) farther tiles
	// repeat while more are trimmed
	// maybe, if there are both strong in-plane and far - see which are closer

	/**
	 * Disable low-textured tiles are not between strong tiles, but on one side of it.
	 * This method relies on the assumption that FG edge should bave strong correlation, so it tries multiple directions
	 * from the weak (not trusted strong) tiles and trims tiles that eithre do not have anything in that direction or have
	 * farther tiles.
	 * Trimming(disabling) weak (trusted but not strong_trusted) tiles if on any one side:
	 *   a) there are no same plane or closer tiles
	 *   b) there are no in-plane or closer strong tiles, but there are some (strong or any?) farther tiles
	 * repeat while more are trimmed
	 * maybe, if there are both strong in-plane and far - see which are closer
	 * @param trusted_strength strength to trust unconditionally
	 * @param strength_rfloor strength floor to subrtact as a fraction of the trusted strength
	 * @param cond_rtrusted fraction of the trusted strength (after subtracting str4ength_floor) that is sufficient
	 * to participate in plane fitting, if successful - make a tile trusted
	 * @param strength_pow raise strength to thyis power (normally just 1.0)
	 * @param smpl_radius sample "radius", square side is  2 * smpl_radius + 1
	 * @param smpl_num minimal absolute number of samples required to try fit a plane and validate a tile
	 * @param smpl_fract minimal fraction number of the neighbor samples that fit the rms filter required to try fit a plane and validate a tile
	 * @param max_adiff maximal absolute difference from the center tile for initial neighbors selection
	 * @param max_rdiff maximal (additional) relative (to tile disparity) difference from the center tile for initial neighbors selection
	 * @param max_atilt maximal absolute tilt (pix/tile) for the tilted planes to fit
	 * @param max_rtilt maximal relative tilt (pix/tile per disparity pixel). min(max_rtilt*disparity, max_atilt) will be used
	 * @param smpl_arms maximal absolute rms of the weighted remaining samples for the successful plane fitting
	 * @param smpl_rrms maximal relative (additional)rms of the weighted remaining samples for the successful plane fitting
	 * @param damp_tilt regularization value to handle planes if the remaining samples are co-linear (or just a single tile)
	 * @param rwsigma weight Gaussina sigma to reduce influence of far tiles relative to smpl_radius
	 * @param atolerance absolute disparity tolerance to what to consider "far"
	 * @param rtolerance relative to disparity tolerance to what to consider "far"
	 * @param num_dirs number of directions per 2*PI to try
	 * @param blind_dist when trying directions require distance in that direction to exceed this value to count
	 * @param strong_only_far true: require far tiles to be trusted_strong, false - just trusted is OK
	 * @param num_strong_far number of strong trusted far tiles to make this FG tile hanging
	 * @param num_weak_far number of strong trusted far tiles to make this FG tile hanging
	 * @param dbg_x tileX to debug
	 * @param dbg_y tileY to debug
	 * @param debugLevel debug level
	 * @return
	 */

	int  trimWeakFG(
		    final double     trusted_strength, // trusted correlation strength
			final double     strength_rfloor,   // strength floor - relative to trusted
		    final double     cond_rtrusted,     // minimal strength to consider - fraction of trusted
			final double     strength_pow,      // raise strength-floor to this power
			final int        smpl_radius,
			final int        smpl_num,   //         = 3;      // Number after removing worst (should be >1)
			final double     smpl_fract, // Number of friends among all neighbors
			final double     max_adiff,  // Maximal absolute difference betweenthe center tile and friends
			final double     max_rdiff, //  Maximal relative difference between the center tile and friends
			final double     max_atilt, //  = 2.0; // pix per tile
			final double     max_rtilt, //  = 0.2; // (pix / disparity) per tile
			final double     smpl_arms, //         = 0.1;    // Maximal RMS of the remaining tiles in a sample
			final double     smpl_rrms,        //      = 0.005;  // Maximal RMS/disparity in addition to smplRms
			final double     damp_tilt, //   =     0.001; // Tilt cost for damping insufficient plane data
			final double     rwsigma,           //  = 0.7; // influence of far neighbors diminish as a Gaussian with this sigma

			final double     atolerance,  // When deciding closer/farther
			final double     rtolerance,  // same, scaled with disparity
			final int        num_dirs,    // number of directions to try
			final double     blind_dist,  // analyze only tiles farther than this in the selected direction
			final boolean    strong_only_far, // in variant b) only compare with strong far
			final int        num_strong_far,    // number of directions to try
			final int        num_weak_far,     // number of directions to try

			final int        dbg_x,
			final int        dbg_y,
			final int        debugLevel

			) {
		final TileNeibs  tnImage = biCamDSI.tnImage;
		final int dbg_tile = (debugLevel>-2)?(dbg_x + tnImage.sizeX*dbg_y):-1;
		final double     wsigma = rwsigma*smpl_radius;
		final double [][] ds = getDisparityStrength( // already has disabled zeroed
	    		false,   // final boolean only_strong,
	    		false,   // final boolean only_trusted,
	    		true) ;  // final boolean only_enabled);

		final int num_tiles = tnImage.getSizeX()*tnImage.getSizeY();
		final double     strength_floor = trusted_strength * strength_rfloor;
//		final double     min_strength = strength_floor + (trusted_strength - strength_floor) *  cond_rtrusted;
		final int        smpl_side = 2 * smpl_radius + 1; // Sample size (side of a square)
		final int        smpl_len = smpl_side*smpl_side;
		// prepare window
		final double [][] weights = new double [smpl_radius + 1][smpl_radius + 1];
		for (int i = 0; i <weights.length; i++) {
			for (int j = i; j <weights[i].length; j++) {
				weights[i][j] = (wsigma >0.0) ?Math.exp(-(i*i+j*j)/(2*wsigma*wsigma)):1.0;
				weights[j][i] = weights[i][j];
			}
		}
		final Thread[] threads = ImageDtt.newThreadArray(biCamDSI.threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		final AtomicInteger num_disabled = new AtomicInteger(0);
		int total_disabled = 0;
        int npass=0;
		while (true) {
			num_disabled.set(0);
			final boolean [] new_disabled = new boolean[num_tiles]; // this.disabled.clone();
			ai.set(0);
			// find definitely trusted and conditionally trusted tiles
			for (int ithread = 0; ithread < threads.length; ithread++) {
				threads[ithread] = new Thread() {
					@Override
					public void run() {
						int smpl_center = (smpl_side + 1) * smpl_radius;
						for (int nTile = ai.getAndIncrement(); nTile < num_tiles; nTile = ai.getAndIncrement()) if (trusted[nTile] && !strong_trusted[nTile]){
							boolean dbg_this = (dbg_tile == nTile);
							if (dbg_this) {
								System.out.println("trimWeakFG(): debugging tile"+nTile);
							}
							double [] smpl_d =       new double  [smpl_len];
							double [] smpl_w =       new double  [smpl_len];
							double [] smpl_w_all =   new double  [smpl_len];
							double [] smpl_p =       new double  [smpl_len]; // plane disparity,
							boolean [] smpl_strong = new boolean[smpl_len];
							// copy neighbor tiles
							double disp_center = ds[0][nTile];
							double max_diff = max_adiff + max_rdiff * disp_center;
							double max_tilt = max_rtilt * disp_center;
							if (max_tilt > max_atilt) {
								max_tilt = max_atilt;
							}

							int nsmpls = 0;
							int nall = 0;

							for (int dy = -smpl_radius; dy <= smpl_radius; dy++) {
								int ady = (dy > 0)? dy:(-dy);
								for (int dx = -smpl_radius; dx <= smpl_radius; dx++) {
									int nTile1 =  tnImage.getNeibIndex(nTile, dx, dy);
									//								if ((nTile1 >= 0) && cond_trusted[nTile1]) {
									if ((nTile1 >= 0) && trusted[nTile1]) {

										nall++;
										int adx = (dx > 0)? dx:(-dx);
										double max_fdiff = max_diff + (ady+adx) * max_tilt;
										int smpl_indx = smpl_center + dy*smpl_side + dx;
										smpl_d[smpl_indx] = ds[0][nTile1];
										double w = (ds[1][nTile1]-strength_floor) * weights[ady][adx];
										if (smpl_d[smpl_indx] < 0.0) { // discard stray negative disparity
											w = 0.0;
										}
										if (strength_pow != 1.0) {
											w = Math.pow(w, strength_pow);
										}
										smpl_w_all[smpl_indx] = w;
										if (Math.abs(ds[0][nTile1] - disp_center) <= max_fdiff) {
											smpl_w[smpl_indx] = w;
											nsmpls ++;
										}
										smpl_strong[smpl_indx]= strong_trusted[nTile1];
									}
								}
							}
							int fin_samples= (int) ( nall * smpl_fract);
							if (fin_samples < smpl_num) fin_samples = smpl_num;

							if (nsmpls >= fin_samples) {
								double [] fit_rslt = fitPlaneRemoveOutliers(
										smpl_radius, // int                     smpl_radius,
										max_tilt,    // double                  max_tilt,
										damp_tilt,   // double                  damp_tilt, //   =     0.001; // Tilt cost for damping insufficient plane data
										true,       // boolean                 full_plane,
										smpl_d,      // double []               smpl_d,
										smpl_w,      // double []               smpl_w, // will be modified,
										smpl_p,      // double []               smpl_p, // will be set if provided
										fin_samples, // int                     fin_samples, // remove until this number remain
										debugLevel); // int                     debugLevel)
								if ( (fit_rslt != null) && (fit_rslt[0] < (smpl_arms + smpl_rrms * fit_rslt[1]))){
									// check directions here, use smpl_w_all, smpl_d and smpl_p
									//								int [] num_this_closer =        new int [num_dirs];
									//								int [] num_this_closer_strong = new int [num_dirs];
									//								int [] num_farther =            new int [num_dirs];



									double thershold_disp = -(atolerance + rtolerance*fit_rslt[1]);
									double [][] dbg_data = dbg_this?(new double [5+num_dirs][smpl_len]):null;
									String [] dbg_titles = dbg_this?(new String [5+num_dirs]):null;
									if (dbg_this) {
										for (int indxs = 0; indxs < smpl_len; indxs++) {
											dbg_data[0][indxs] = smpl_d[indxs];
											dbg_data[1][indxs] = smpl_strong[indxs]? smpl_d[indxs]: Double.NaN;
											dbg_data[2][indxs] = smpl_p[indxs];
											dbg_data[3][indxs] = smpl_w_all[indxs];
											dbg_data[4][indxs] = smpl_d[indxs]-smpl_p[indxs];
										}
										dbg_titles[0]= "smpl_d";
										dbg_titles[1]= "smpl_d-strong";
										dbg_titles[2]= "smpl_p";
										dbg_titles[3]= "smpl_w_all";
										dbg_titles[4]= "smpl_d-smpl_p";
										for (int idir = 0; idir<num_dirs; idir++) {
											dbg_titles[5 + idir]="dir-"+idir;
											for (int indxs = 0; indxs < smpl_len; indxs++) {
												dbg_data[5+idir][indxs] = Double.NaN;
											}
										}
									}

									for (int idir = 0; idir<num_dirs; idir++) {
										double a = 2*Math.PI*idir/num_dirs;
										double ca = Math.cos(a);
										double sa = Math.sin(a);
										int num_this_closer =        0;
										int num_this_closer_strong = 0;
										int num_farther =            0;
										int num_farther_strong =     0;


										for (int indxs = 0; indxs < smpl_len; indxs++) if (smpl_w_all[indxs] > 0.0) {
											int dy = indxs/smpl_side - smpl_radius;
											int dx = indxs%smpl_side - smpl_radius;
											double d = dx*ca + dy*sa;
											if (dbg_this) {
												dbg_data[5+idir][indxs] = d - blind_dist;
											}
											if (d >= blind_dist) {
												boolean farther = (smpl_d[indxs]-smpl_p[indxs]) < thershold_disp;
												if (farther) {
													num_farther ++;
//													if (!strong_only_far || smpl_strong[indxs]) num_farther++;
													if (smpl_strong[indxs]) num_farther_strong++;
												} else {
													num_this_closer++;
													if (smpl_strong[indxs]) num_this_closer_strong++;
												}
											}

										}
										if (dbg_this) {
											System.out.println("trimWeakFG(): idir="+idir+" num_this_closer="+num_this_closer+
													", num_this_closer_strong="+num_this_closer_strong+", num_farther="+num_farther);
										}
										boolean far_exist = (num_farther_strong > num_strong_far) || (!strong_only_far && (num_farther > num_weak_far));

										if ((num_this_closer ==0) || far_exist) {
											new_disabled[nTile] =    true;
											num_disabled.getAndIncrement();
											if (dbg_this) {
												System.out.println("trimWeakFG(): DISABLED idir="+idir+" num_this_closer="+num_this_closer+
														", num_this_closer_strong="+num_this_closer_strong+", num_farther="+num_farther+", num_farther_strong="+num_farther_strong);
											}
											if (!dbg_this) {
												break; // already disabled
											}

										}
									}
									if (dbg_this) {
									(new showDoubleFloatArrays()).showArrays(
											dbg_data,
											smpl_side,
											smpl_side,
											true,
											"trimWeakFG_T"+nTile,
											dbg_titles);
									}


								} else { // may get even after initial filtering (by calcTrusted()) here because some neighbors became disabled? loosen flatness? Keep non-flat
									if (dbg_this) {
										System.out.println("trimWeakFG(): Too high RMS");
									}

								}
							} else { // may get even after initial filtering (by calcTrusted()) here because some neighbors became disabled?
								if (dbg_this) {
									System.out.println("trimWeakFG(): Insufficient points");
								}
							}
						}
					}
				};
			}
			ImageDtt.startAndJoin(threads);
			npass++;
			if (debugLevel > -2) {
				System.out.println("trimWeakFG(): npass="+npass+" removed="+num_disabled.get());
			}
			if (num_disabled.get() == 0) {
				break;
			}
			ai.set(0);
			// find definitely trusted and conditionally trusted tiles
			for (int ithread = 0; ithread < threads.length; ithread++) {
				threads[ithread] = new Thread() {
					@Override
					public void run() {
						for (int nTile = ai.getAndIncrement(); nTile < num_tiles; nTile = ai.getAndIncrement())	if (new_disabled[nTile]) {
							disableTile(nTile); // disabled source and this trusted*
						}
					}
				};
			}
			ImageDtt.startAndJoin(threads);
			if (debugLevel > -2) {
				biCamDSI.getLastBiScan().showScan(
						"trimWeakFG_"+npass);
			}

			total_disabled += num_disabled.get();
		}
		return total_disabled; // numTrustedAll;
	}



//					smpl_p[indxs] = approx2d[0][0] * (sx - smpl_radius) + approx2d[0][1] * (sy - smpl_radius) + approx2d[0][2];

	/**
	 * Remove worst outliers when fitting planes until the specified number of samples are left, re-fitting
	 * plane after each iteration
	 * @param smpl_radius sample "radius", square side is  2 * smpl_radius + 1
	 * @param max_tilt maximal DSI plane tilt (pixels/tile)
	 * @param damp_tilt tilt cost for damping insufficient plane data
	 * @param full_plane generate full tilted plane, not only for non-zero tiles
	 * @param smpl_d  array of disparity values to approximate
	 * @param smpl_w  weights - only >0.0 are processed
	 * @param smpl_p approximated disparity values, amy be null. When provided (not null), will have calculated disparity approximation
	 * @param fin_samples
	 * @param debugLevel debug level
	 * @return array of 4 values: rmas of the remaining samples,  average (center) disparity value, tiltX, tiltY
	 */
	public double [] fitPlaneRemoveOutliers(
			int                     smpl_radius,
			double                  max_tilt,
			double                  damp_tilt, //   =     0.001; // Tilt cost for damping insufficient plane data
			boolean                 full_plane,
			double []               smpl_d,
			double []               smpl_w, // will be modified,
			double []               smpl_p, // will be set if provided
			int                     fin_samples, // remove until this number remain
			int                     debugLevel)
	{
		double [] disp_tilts = null;
		final int        smpl_side = 2 * smpl_radius + 1; // Sample size (side of a square)
		int num_in_sample = 0;
		for (int i = 0; i < smpl_w.length; i++) if (smpl_w[i]> 0.0) num_in_sample++;
		int smpl_len = smpl_w.length;
		if (num_in_sample <fin_samples) {
			return null;
		}
		double rms = Double.NaN;
		if (smpl_p == null) smpl_p = new double [smpl_len];
		while (num_in_sample >= fin_samples) {
			// fit plane to the selected and remove outlayers after each iteration
			disp_tilts = fitPlane(
					smpl_radius,  //int        smpl_radius,
					max_tilt,     //double     max_tilt,
					damp_tilt,    //double     damp_tilt, //   =     0.001; // Tilt cost for damping insufficient plane data
					full_plane,   // boolean                 full_plane,
					smpl_d,       // double []  smpl_d,
					smpl_w,       // double []  smpl_w,
					smpl_p,       //double []  smpl_p, // will be set if provided
					debugLevel);  // int        debugLevel)
			// calculate fitting quality
			if (disp_tilts == null) {
				return null;
			}
//			double d_center = disp_tilts[0];
			double sw = 0.0, sd2 = 0.0;
			for (int indxs = 0; indxs < smpl_len;indxs++) if (smpl_w[indxs] > 0) {
				double d = smpl_d[indxs] - smpl_p[indxs];
				double dw = d * smpl_w[indxs];
				sd2 += dw * smpl_d[indxs];
				sw +=       smpl_w[indxs];
			}
			rms = (sw > 0.0)? Math.sqrt(sd2/sw): Double.NaN;

			// remove worst - it should not make remaining set
			if (num_in_sample > fin_samples) { // remove worst if it is not the last run where only calculations are needed
				int iworst = -1;
				double dworst2 = 0.0;
				for (int indxs = 0; indxs < smpl_len; indxs++) if (smpl_w[indxs] > 0.0) {
					double d2 = smpl_d[indxs] - smpl_p[indxs];
					d2 *=d2;
					if (d2 > dworst2) {
						if ((damp_tilt !=0.0) || notColinearWithout (
								indxs, // int        indx,
								smpl_w, // boolean [] sel,
								smpl_side)) { // int side))
							iworst = indxs;
							dworst2 = d2;
						}
					}
				}
				if (iworst < 0){
					if (debugLevel > 0) {
						System.out.println("**** this may be BUG in fitPlaneRemoveOutliers() can not find the worst sample  - all tiles fit perfectly ****");
					}
					// this can happen if some samples are the same and all the pixels fit exactly - use all of them
					break;
				}
				// remove worst sample
				smpl_w[iworst] = 0.0;
				num_in_sample --;
			} else {
				break;
			}
		} // while (nsmpls >= min_sample) {
		double [] rslt = {rms, disp_tilts[0], disp_tilts[1], disp_tilts[1]};
		return rslt;
	}

	/**
	 * Verify that selected points are not all on the same line, even if the specified one is removed
	 * @param indx index of the point to be removed
	 * @param smpl_w 2-d sample weights in linescan order (used zero/non-zero only)
	 * @param side square samples side
	 * @return true if there are enough samples for plane extraction, false otherwise
	 */

	public boolean notColinearWithout (
			int        indx,
			double [] smpl_w,
			int side)
	{
		if (smpl_w[indx] == 0.0){
			throw new IllegalArgumentException ("notCoplanarWithout(): specified is the non existing index");
		}
		double w = smpl_w[indx];
		smpl_w[indx] = 0.0;
		boolean rslt = notColinear ( smpl_w, side);
		smpl_w[indx] = w; // restore value
		return rslt;
	}
	/**
	 * Verify that selected points are not all on the same line
	 * @param smpl_w 2-d sample weights in linescan order (used zero/non-zero only)
	 * @param side square samples side
	 * @return true if there are enough samples for plane extraction, false otherwise
	 */
	public boolean notColinear (
			double [] smpl_w,
			int side)
	{
		int indx0, indx1;
		for (indx0 = 0; indx0 < smpl_w.length; indx0++){
			if (smpl_w[indx0] > 0.0) break;
		}
		for (indx1 = indx0+1; indx1 < smpl_w.length; indx1++){
			if (smpl_w[indx0] > 0.0) break;
		}
		if (indx1 >= smpl_w.length) return false; // too few points;
		int sx0 = indx0 % side;
		int sy0 = indx0 / side;
		int sx1 = indx1 % side;
		int sy1 = indx1 / side;
		for (int indx = indx1 +1; indx < smpl_w.length; indx++){
			int sx = indx % side;
			int sy = indx / side;
			if ((sx - sx0) * (sy - sy1) != (sx - sx1) * (sy - sy0)){
				return true;
			}
		}
		return false;
	}



	/**
	 * Approximate tile disparity by a tilted DSI plane
	 * @param smpl_radius sample "radius", square side is  2 * smpl_radius + 1
	 * @param max_tilt maximal DSI plane tilt (pixels/tile)
	 * @param damp_tilt tilt cost for damping insufficient plane data
	 * @param full_plane generate full tilted plane, not only for non-zero tiles
	 * @param smpl_d  array of disparity values to approximate
	 * @param smpl_w  weights - only >0.0 are processed
	 * @param smpl_p approximated disparity values, amy be null. When provided (not null), will have calculated disparity approximation
	 * @param debugLevel debug level
	 * @return array of 3 values: average (center) disparity value, tiltX, tiltY
	 */
	public double [] fitPlane(
			int                     smpl_radius,
			double                  max_tilt,
			double                  damp_tilt, //   =     0.001; // Tilt cost for damping insufficient plane data
			boolean                 full_plane,
			double []               smpl_d,
			double []               smpl_w,
			double []               smpl_p, // will be set if provided
			int                     debugLevel)
	{
		PolynomialApproximation pa = new PolynomialApproximation();
		final double [] damping = {damp_tilt, damp_tilt, 0.0}; // 0.0 will be applied to average value, tilt_cost - to both tilts
		final int        smpl_side = 2 * smpl_radius + 1; // Sample size (side of a square)
		final int        smpl_len = smpl_side*smpl_side;
		if (smpl_p == null) smpl_p = new double [smpl_len];
		int num_in_sample = 0;
		for (int i = 0; i < smpl_w.length; i++) if (smpl_w[i] > 0.0) num_in_sample++;

		double [][][] mdata = new double [num_in_sample][3][];
		int mindx = 0;
		for (int sy = 0; sy < smpl_side; sy++){
			for (int sx = 0; sx < smpl_side; sx++){
				int indxs = sy * smpl_side + sx;
				if (smpl_w[indxs] > 0.0) {
					mdata[mindx][0] = new double [2];
					mdata[mindx][0][0] =  sx - smpl_radius;
					mdata[mindx][0][1] =  sy - smpl_radius;
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
				THRESHOLD_LIN,  // threshold ratio of matrix determinant to norm for linear approximation (det too low - fail)
				THRESHOLD_QUAD, // threshold ratio of matrix determinant to norm for quadratic approximation (det too low - fail)
				debugLevel);
		if (approx2d == null){
			if (debugLevel > -1){
				System.out.println("getDisparityStrengthML(): can not find linear approximation");
			}
			return null;
		}
		// limit tilt to be within range
		//											double     max_abs_tilt, //  = 2.0; // pix per tile
		//											double     max_rel_tilt, //  = 0.2; // (pix / disparity) per tile
//		double max_tilt = Math.min(mlfp.max_abs_tilt, mlfp.max_rel_tilt * approx2d[0][2]);
		boolean overlimit = (Math.abs(approx2d[0][0]) > max_tilt) || (Math.abs(approx2d[0][1]) > max_tilt);
		if (overlimit) {
			approx2d[0][0] = Math.min(approx2d[0][0],  max_tilt);
			approx2d[0][1] = Math.min(approx2d[0][1],  max_tilt);
			approx2d[0][0] = Math.max(approx2d[0][0], -max_tilt);
			approx2d[0][1] = Math.max(approx2d[0][1], -max_tilt);
		}
		// subtract tilt from disparity
		for (int sy = 0; sy < smpl_side; sy++){
			for (int sx = 0; sx < smpl_side; sx++){
				int indxs = sy * smpl_side + sx;
				if ((smpl_w[indxs] > 0.0) || full_plane) {
					smpl_p[indxs] = approx2d[0][0] * (sx - smpl_radius) + approx2d[0][1] * (sy - smpl_radius) + approx2d[0][2];
				}
			}
		}

		if (overlimit){ // re-calculate disparity average (in the center)
			double sw = 0.0;
			double sd=0.0;
			for (int indxs = 0; indxs < smpl_len;indxs++) if (smpl_w[indxs] > 0.0) {
				double d = smpl_d[indxs] - smpl_p[indxs];
				double dw = d * smpl_w[indxs];
				sd += dw;
				sw += smpl_w[indxs];
			}
			sd /= sw;
			for (int indxs = 0; indxs < smpl_len;indxs++) if ((smpl_w[indxs] > 0.0) || full_plane) {
				smpl_p[indxs] += sd;
			}
			approx2d[0][2] += sd;
		}
		double [] rslt = {approx2d[0][2], approx2d[0][0],approx2d[0][1]}; // {center, tiltX, tiltY
		return rslt;
	}



}
