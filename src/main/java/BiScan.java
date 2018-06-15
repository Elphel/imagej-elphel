/**
 ** BiCamScan - class to represent multiple bi-quad camera measurements
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
	final static int    BISCAN_ANY =       -1;
	final static int    BISCAN_SINGLECORR = 0;
	final static int    BISCAN_AVGCORR =    1; // combined with low-texture averaging correlation
	final static int    BISCAN_POLE =       2; // combined with low-texture averaging correlation

	double []  disparity_measured;
	double []  target_disparity;
	double []  strength_measured;
	boolean [] strong_trusted; // sufficient strength without neighbors
	boolean [] trusted;
	boolean [] cond_trusted;
	boolean [] disabled_measurement; // should disable source, not this!
	int     [] src_index;       // index of the source scan which measured data is used here (applies to disparity_measured, strength_measured, disabled_measurement
	int        list_index = -1;
	int        scan_type = -1;
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
			boolean [] disabled,
			int scan_type) {
		this.biCamDSI = biCamDSI;
		this.list_index = indx;
		this.scan_type = scan_type;
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
			biCamDSI.getBiScan(src_index[nTile]).disabled_measurement[nTile] = true; // false; // may be source tile or this tile
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

    // trusted should be set, copied and replaced as needed
    public double [][] getFilteredDisparityStrength( // FIXME
			final boolean [] area_of_interest,
			final double [][] disparityStrength,
			final double     min_disparity,    // keep original disparity far tiles
		    final double     trusted_strength, // trusted correlation strength
			final double     strength_rfloor,   // strength floor - relative to trusted
		    final boolean    discard_unreliable,// replace v
		    final boolean    discard_weak,      // consider weak trusted tiles (not promoted to trusted) as empty
		    final boolean    discard_strong,    // suggest new disparitieas even for strong tiles
			final double     strength_pow,      // raise strength-floor to this power
			final double []  smpl_radius_array, // space-variant radius
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
			final double     center_weight,     // use center tile too (0.0 - do not use)
			final boolean    use_alt,           // use tiles from other scans if they fit better
			final double     goal_fraction_rms, // Try to make rms to be this fraction of maximal acceptable by removing outliers
			final double     boost_low_density, // 0 - strength is proportional to 1/density, 1.0 - same as remaining tiles
			final int        fourq_min,         // each of the 4 corners should have at least this number of tiles.
			final int        fourq_gap,         // symmetrical vertical and horizontal center areas that do not belong to any corner
			final int        dbg_x,
			final int        dbg_y,
			final int        debugLevel
    		){
		final int num_tiles = biCamDSI.tnImage.getSizeX()*biCamDSI.tnImage.getSizeY();
		double [][] ds0 = 	getDisparityStrength( // FIXME
	    		false,   // final boolean only_strong,
	    		false,   // final boolean only_trusted,
	    		true) ;  // final boolean only_enabled);


    	double [][] ds = new double[2][num_tiles];
    	for (int i = 0; i < num_tiles; i++) ds[0][i] = Double.NaN;
//    	double boost_low_density = 0.8; // 1.0; //0.2;
		  suggestNewScan(
				  area_of_interest,  // final boolean [] area_of_interest,
				  disparityStrength, // final double [][] disparityStrength,
				  trusted_strength,  // final double     trusted_strength, // trusted correlation strength
				  strength_rfloor,   // final double     strength_rfloor,   // strength floor - relative to trusted
				  true,              // final boolean    discard_cond,      // consider conditionally trusted tiles (not promoted to trusted) as empty
				  true,              // final boolean    discard_weak,      // consider conditionally trusted tiles (not promoted to trusted) as empty
				  true,              // final boolean    discard_strong,      // consider conditionally trusted tiles (not promoted to trusted) as empty
				  strength_pow,      // final double     strength_pow,      // raise strength-floor to this power
				  smpl_radius,       // final int        smpl_radius,
				  smpl_num,          // final int        smpl_num,   //         = 3;      // Number after removing worst (should be >1)
				  smpl_fract,        // final double     smpl_fract, // Number of friends among all neighbors
				  smpl_num_narrow,   // final int        smpl_num_narrow,   //         = 3;      // Number after removing worst (should be >1)
				  max_adiff,         // final double     max_adiff,  // Maximal absolute difference betweenthe center tile and friends
				  max_rdiff,         // final double     max_rdiff, //  Maximal relative difference between the center tile and friends
				  max_atilt,         // final double     max_atilt, //  = 2.0; // pix per tile
				  max_rtilt,         // final double     max_rtilt, //  = 0.2; // (pix / disparity) per tile
				  smpl_arms,         // final double     smpl_arms, //         = 0.1;    // Maximal RMS of the remaining tiles in a sample
				  smpl_rrms,         // final double     smpl_rrms,        //      = 0.005;  // Maximal RMS/disparity in addition to smplRms
				  damp_tilt,         // final double     damp_tilt, //   =     0.001; // Tilt cost for damping insufficient plane data
				  rwsigma,           // final double     rwsigma,           //  = 0.7; // influence of far neighbors diminish as a Gaussian with this sigma
				  rwsigma_narrow,    // final double     rwsigma_narrow,    //  = used to determine initial tilt
				  0.0,               // final double     new_diff,            // minimal difference between the new suggested and the already tried/measured one
				  false,             // final boolean    remove_all_tried,  // remove from suggested - not only disabled, but all tried
				  center_weight,     // final double     center_weight,     // use center tile too (0.0 - do not use)
				  use_alt,           // final boolean    use_alt,           // use tiles from other scans if they fit better
				  boost_low_density, // final double     boost_low_density, // 0 - strength is proportional to 1/density, 1.0 - same as remaining tiles
				  goal_fraction_rms, // final double     goal_fraction_rms, // Try to make rms to be this fraction of maximal acceptable by removing outliers
				  ds,                // final double [][]  smooth_ds,   // optionally fill strength array when used for smoothing DSI
				  fourq_min,         // final int        fourq_min,         // each of the 4 corners should have at least this number of tiles.
				  fourq_gap,         // final int        fourq_gap,         // symmetrical vertical and horizontal center areas that do not belong to any corner
				  dbg_x,             // final int        dbg_x,
				  dbg_y,             // final int        dbg_y,
				  debugLevel);       // final int        debugLevel);
		for (int nTile = 0; nTile < num_tiles; nTile++) {
			if ((ds0[1][nTile] > 0.0) && (
					(ds[0][nTile] <= 0.0) ||
					!discard_unreliable ||
					(strong_trusted[nTile] && !discard_strong) ||
					(trusted[nTile] && !discard_weak) ||
					(ds[0][nTile] <  min_disparity))) {
				ds[0][nTile] = ds0[0][nTile];
				ds[1][nTile] = ds0[1][nTile];
			}
		}
    	return ds;
    }
	public void showScan(String title) {
		showScan(title,null);
	}

	public void showScan(String title, double [][] ext_ds) {
		String [] titles= {
				"ext_disp",          //    0
				"all",               // 0  1
				"enabled",           // 1  2
				"cond_trusted",      // 2  3
				"weak trusted",      // 3  4
				"strong trusted",    // 4  5
				"measured",          // 5  6
				"suggested",         // 6  7
				"ext strength",      //    8
				"measured strength", // 7  9
				"strength" };        // 8 10
		double [][] ds_all = getDisparityStrength(
	    		false,   // final boolean only_strong,
	    		false,   // final boolean only_trusted,
	    		false) ; // final boolean only_enabled);
		double [][] ds = getDisparityStrength(
	    		false,   // final boolean only_strong,
	    		false,   // final boolean only_trusted,
	    		true) ; // final boolean only_enabled);

		double [][] dbg_data = new double[titles.length][];
		dbg_data[ 6] = this.disparity_measured;
		dbg_data[ 1] = ds_all[0];
		dbg_data[ 9] = this.strength_measured;
		dbg_data[10] = ds_all[1];
		dbg_data[ 2] = ds[0];
		if (ext_ds != null) {
			dbg_data[ 0] = ext_ds[0];
			dbg_data[ 8] = ext_ds[1];
		}

		if (this.cond_trusted != null) {
			dbg_data[3] = ds[0].clone();
			for (int i = 0; i < this.cond_trusted.length; i++) if (!cond_trusted[i]) dbg_data[3][i] = Double.NaN;
		}
		if (this.trusted != null) {
			dbg_data[4] = ds[0].clone();
			for (int i = 0; i < this.trusted.length; i++) if (!trusted[i]) dbg_data[4][i] = Double.NaN;
		}
		if (this.strong_trusted != null) {
			dbg_data[5] =  ds[0].clone();
			for (int i = 0; i < this.strong_trusted.length; i++) if (!strong_trusted[i]) dbg_data[5][i] = Double.NaN;
		}
		if (this.target_disparity != null) {
			dbg_data[7] = this.target_disparity.clone();
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
		int debugLevel = -10;
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
	 * 2) target disparity that lead to the current measurement after refinement
	 * 3) any other disable measurement
	 * 4) any target disparity that lead to the disabled measurement
	 * @param area_of_interest - limit results to these tiles (if provided)
	 * @param disparityStrength - a pair of array or null. If null, will calculate fro the current scan
	 *        if not null - use as is
	 * @param trusted_strength strength to trust unconditionally
	 * @param strength_rfloor strength floor to subtract as a fraction of the trusted strength
	 * @param discard_cond if true may suggest new disparities for conditionally trusted tiles
	 * @param discard_weak if true may suggest new disparities over trusted weak tiles
	 * @param discard_stron if true may suggest new disparities over any tile
	 * @param strength_pow raise strength to this power (normally just 1.0)
	 * @param smpl_radius sample "radius", square side is  2 * smpl_radius + 1
	 * @param smpl_num minimal absolute number of samples required to try fit a plane and validate a tile
	 * If smpl_num == 0, faster calculation (single pass) using only *_narrow settings
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
	 * @param rwsigma_narrow Gaussian sigma for the preliminary plain fitting using the closest tiles ~= 1/smpl_radius
	 * @param new_diff minimal difference between the new suggested and the already tried/measured one
	 * @param remove_all_tried remove from suggested - not only disabled, but all tried
	 * @param center_weight weight of the tile itself (0.0 - do not use). Should be set to 0.0 for suggesting, >0 - for "smoothing"
	 * @param use_alt use tiles from other scans if they fit better
	 * @param goal_fraction_rms try to make rms to be this fraction of maximal acceptable by removing outliers
	 * @param boost_low_density 0.0 - strength is proportional to 1/density, 1.0 - same as remaining tiles
	 * @param smooth_ds optionally fill disparity/strength  instead of the target_disparity
	 * @param fourq_min each of the 4 corners should have at least this number of tiles.
	 * @param fourq_gap symmetrical vertical and horizontal center areas that do not belong to any corner
	 * @param dbg_x tileX to debug
	 * @param dbg_y tileY to debug
	 * @param debugLevel debug level
	 * @return number of new tiles to measure in the  array of suggested disparities - Double.NaN - nothing suggested
	 *  for the tile. May need additional filtering to avoid suggested already tried disparities
	 */

	int  suggestNewScan(
			final boolean [] area_of_interest,
			final double [][] disparityStrength,
		    final double     trusted_strength, // trusted correlation strength
			final double     strength_rfloor,   // strength floor - relative to trusted
		    final boolean    discard_cond,      // consider conditionally trusted tiles (not promoted to trusted) as empty
		    final boolean    discard_weak,      // consider weak trusted tiles (not promoted to trusted) as empty
		    final boolean    discard_strong,    // suggest new disparities even for strong tiles
			final double     strength_pow,      // raise strength-floor to this power
			final int        smpl_radius,
			final int        smpl_num,   //         = 3;      // Number after removing worst (should be >1)
			final double     smpl_fract, // Number of friends among all neighbors
			final int        smpl_num_narrow,   //         = 3;      // Number after removing worst (should be >1)
			final double     max_adiff,  // Maximal absolute difference between the center tile and friends
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
			final double     center_weight,     // use center tile too (0.0 - do not use)
			final boolean    use_alt,           // use tiles from other scans if they fit better
			final double     goal_fraction_rms, // Try to make rms to be this fraction of maximal acceptable by removing outliers
			final double     boost_low_density, // 0 - strength is proportional to 1/density, 1.0 - same as remaining tiles
			final double [][] smooth_ds,        // optionally fill disparity/strength  instead of the target_disparity
			final int        fourq_min,         // each of the 4 corners should have at least this number of tiles.
			final int        fourq_gap,         // symmetrical vertical and horizontal center areas that do not belong to any corner
			final int        dbg_x,
			final int        dbg_y,
			final int        debugLevel
			) {
		final TileNeibs  tnImage = biCamDSI.tnImage;
		final int dbg_tile = (debugLevel > 0)?(dbg_y * tnImage.sizeX + dbg_x): -1;
		final double     wsigma = rwsigma*smpl_radius;
		final double     wsigma_narrow = rwsigma_narrow*smpl_radius;
		final int        smpl_side = 2 * smpl_radius + 1; // Sample size (side of a square)
		final int        smpl_len = smpl_side*smpl_side;
		final int        smpl_center = (smpl_side + 1) * smpl_radius;
		// prepare window
		final double [][] weights =        new double [smpl_radius + 1][smpl_radius + 1];
		final double [][] weights_narrow = new double [smpl_radius + 1][smpl_radius + 1];
		{ // normalize
			for (int i = 0; i <= smpl_radius; i++) {
				for (int j = i; j <= smpl_radius; j++) {
					weights[i][j] = (wsigma >0.0) ?Math.exp(-(i*i+j*j)/(2*wsigma*wsigma)):1.0;
					weights[j][i] = weights[i][j];
					weights_narrow[i][j] = (wsigma_narrow >0.0) ?Math.exp(-(i*i+j*j)/(2*wsigma_narrow*wsigma_narrow)):1.0;
					weights_narrow[j][i] = weights_narrow[i][j];
				}
			}
			weights[0][0] *=        center_weight;
			weights_narrow[0][0] *= center_weight;
			double sw_full = 0.0, sw_narrow = 0.0;
			for (int i = 0; i <= smpl_radius; i++) {
				for (int j = i; j <= smpl_radius; j++) {
					sw_full +=   weights[i][j];
					sw_narrow += weights_narrow[i][j];
				}
			}
			double k_full = 1.0/sw_full;
			double k_narrow = 1.0/sw_narrow;
			for (int i = 0; i <= smpl_radius; i++) {
				for (int j = i; j <= smpl_radius; j++) {
					weights[i][j] *=k_full;
					weights_narrow[i][j]*=k_narrow;
				}
			}
		}
		final int [] fourq_corner = new int[smpl_len];
		for (int i = 0; i < smpl_side; i++) {
			for (int j = 0; j < smpl_side; j++) {
				int indx = i* smpl_side + j;
				if     (((i > (smpl_radius - fourq_gap)) && (i < (smpl_radius + fourq_gap))) ||
						((j > (smpl_radius - fourq_gap)) && (j < (smpl_radius + fourq_gap)))){
					fourq_corner[indx] = 4; // will not be used

				} else {
					fourq_corner[indx] = ((i > smpl_radius) ? 2 : 0) + ((j > smpl_radius) ? 1 : 0); //122
				}
			}
		}


		final double [][] ds = (disparityStrength != null) ? disparityStrength: getDisparityStrength(
	    		false,   // final boolean only_strong,
	    		false,   // final boolean only_trusted,
	    		true) ; // final boolean only_enabled);
		final int num_tiles = tnImage.getSizeX()*tnImage.getSizeY();
		final double     strength_floor = trusted_strength * strength_rfloor;

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
					for (int nTile = ai.getAndIncrement(); nTile < num_tiles; nTile = ai.getAndIncrement()) if (
							((area_of_interest == null) || area_of_interest[nTile]) && (discard_strong || !trusted_sw[nTile])){
						boolean debug = nTile == dbg_tile;
						if (debug) {
							System.out.println("suggestNewScan(): debbugging nTile="+nTile);
							System.out.println("suggestNewScan(): debbugging nTile="+nTile);
						}
						// Select all neighbors, then filter
						double [] smpl_d =           new double  [smpl_len];
						double [] smpl_w =           new double  [smpl_len];
						double [] smpl_w_narrow =    new double  [smpl_len];
						double [] smpl_p =           new double  [smpl_len]; // plane disparity,
						double [] smpl_wsw= smpl_w_narrow;
						double [][] weights_sw = weights_narrow;

						int nall = 0;
						double sw = 0, swd = 0;
						for (int dy = -smpl_radius; dy <= smpl_radius; dy++) {
							int ady = (dy > 0)? dy:(-dy);
							for (int dx = -smpl_radius; dx <= smpl_radius; dx++) { // if ((dx != 0) || (dy != 0)){
								int nTile1 =  tnImage.getNeibIndex(nTile, dx, dy);
								if ((nTile1 >= 0) && trusted[nTile1]) { // weak trusted OK
									nall++;
									int adx = (dx > 0)? dx:(-dx);
									int smpl_indx = smpl_center + dy*smpl_side + dx;
									double w = ds[1][nTile1]-strength_floor;
									if ( w > 0) {
										if (strength_pow != 1.0) {
											w = Math.pow(w, strength_pow);
										}

										smpl_d[smpl_indx] =       ds[0][nTile1];
										smpl_w[smpl_indx] =        w * weights[ady][adx];
										smpl_w_narrow[smpl_indx] = w * weights_narrow[ady][adx];
										sw += smpl_w_narrow[smpl_indx];
										swd += smpl_w_narrow[smpl_indx]* smpl_d[smpl_indx];
									}
								}
							}
						}
						if (debug) {
							System.out.println("suggestNewScan(): sw="+sw);
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
						// smpl_num == 0 - special (fast) case do not use wide selection at all
						if ((smpl_num == 0) || (fin_samples < smpl_num)) fin_samples = smpl_num;

						// fit plane to mostly centertiles
						int nsmpls = nall;
						if (nsmpls < smpl_num_narrow) { // no tiles even to start
							continue; //
						}
						double max_rms = smpl_arms + smpl_rrms * disp_mean; // do not need to wait fro the final disparity for this estimate
						double [] fit_rslt = fitPlaneRemoveOutliers(
								smpl_radius,        // int                     smpl_radius,
								max_tilt,           // double                  max_tilt,
								damp_tilt,          // double                  damp_tilt, //   =     0.001; // Tilt cost for damping insufficient plane data
								true,               // boolean                 full_plane,
								smpl_d,             // double []               smpl_d,
								smpl_w_narrow,      // double []               smpl_w, // will be modified,
								smpl_p,             // double []               smpl_p, // will be set if provided
								smpl_num_narrow,    // int                     fin_samples, // remove until this number remain
								goal_fraction_rms*max_rms, // double                  fin_rms,
								false,              // boolean                 keep_center, // do not remove center tile - it is the tile that should be verified by neighbors
								fourq_min,          // int                     fourq_min,         // each of the 4 corners should have at least this number of tiles.
								fourq_corner,       // int [] 	               fourq_corner, //  array specifying corner number (0..3), -1 - gap. null when not used
								debugLevel);        // int                     debugLevel)
//						if ( (fit_rslt == null) ||  (fit_rslt[0] > (smpl_arms + smpl_rrms * fit_rslt[1]))){
						if ( (fit_rslt == null) ||  (fit_rslt[0] > max_rms)){
							continue; // narrow selection - too high rms
						}
						max_rms = smpl_arms + smpl_rrms * fit_rslt[1]; // updated center disparity and so rms
						disp_mean = fit_rslt[1]; // smpl_p[smpl_center]; // center of the fitted plane
						// try to use alternatives for the discarded (not only) tiles
						if (use_alt) {
							int nnew = findBetterFitToPlane( // ; // new assignments
									smpl_radius,     // int       smpl_radius,
									nTile,           // int       nTile,
									strength_floor,  // double    strength_floor,
									strength_pow,    // double    strength_pow,
									weights,         // double [][] weights,
									weights_narrow,  // double [][] weights_narrow,
									smpl_w,          // double [] smpl_w,
									smpl_w_narrow,   // double [] smpl_w_narrow,
									smpl_d,          // double [] smpl_d,
									smpl_p);         // double [] smpl_p);

							if (nnew > 0) { // there were some changes, recalculate narrow fit plane
								fit_rslt = fitPlaneRemoveOutliers(
										smpl_radius,        // int                     smpl_radius,
										max_tilt,           // double                  max_tilt,
										damp_tilt,          // double                  damp_tilt, //   =     0.001; // Tilt cost for damping insufficient plane data
										true,               // boolean                 full_plane,
										smpl_d,             // double []               smpl_d,
										smpl_w_narrow,      // double []               smpl_w, // will be modified,
										smpl_p,             // double []               smpl_p, // will be set if provided
										smpl_num_narrow,    // int                     fin_samples, // remove until this number remain
										goal_fraction_rms*max_rms,            // double                  fin_rms,
										false,              // boolean                 keep_center, // do not remove center tile - it is the tile that should be verified by neighbors
										fourq_min,          // int                     fourq_min,         // each of the 4 corners should have at least this number of tiles.
										fourq_corner,       // int [] 	               fourq_corner, //  array specifying corner number (0..3), -1 - gap. null when not used
										debugLevel);        // int                     debugLevel)
//								if ( (fit_rslt == null) ||  (fit_rslt[0] > (smpl_arms + smpl_rrms * fit_rslt[1]))){
								if ( (fit_rslt == null) ||  (fit_rslt[0] > max_rms)){
									continue; // narrow selection - too high rms
								}
								disp_mean = fit_rslt[1]; // smpl_p[smpl_center]; // center of the fitted plane
							}
// plane may have changed, look for the best fit in history again
							 nnew = findBetterFitToPlane( // ; // new assignments
										smpl_radius,     // int       smpl_radius,
										nTile,           // int       nTile,
										strength_floor,  // double    strength_floor,
										strength_pow,    // double    strength_pow,
										weights,         // double [][] weights,
										weights_narrow,  // double [][] weights_narrow,
										smpl_w,          // double [] smpl_w,
										smpl_w_narrow,   // double [] smpl_w_narrow,
										smpl_d,          // double [] smpl_d,
										smpl_p);         // double [] smpl_p);
						} //if (use_alt) {

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
						if (fin_samples > 0) {

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
									goal_fraction_rms*max_rms,            // double                  fin_rms,
									false,              // boolean                 keep_center, // do not remove center tile - it is the tile that should be verified by neighbors
									fourq_min,          // int                     fourq_min,         // each of the 4 corners should have at least this number of tiles.
									fourq_corner,       // int [] 	               fourq_corner, //  array specifying corner number (0..3), -1 - gap. null when not used
									debugLevel);        // int                     debugLevel)
							smpl_wsw = smpl_w;
							weights_sw = weights;
							//						if ( (fit_rslt == null) ||  (fit_rslt[0] > (smpl_arms + smpl_rrms * fit_rslt[1]))){
							if ( (fit_rslt == null) ||  (fit_rslt[0] > max_rms)){
								continue; // narrow selection - too high rms
							}
						}
						boolean valid_suggestion = true;
						if ((disparityStrength == null) && (new_diff > 0.0)) {
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
							// If disparityStrength[][] is provided, do not use history. If difference is <= 0 - also no sense to look there
							//							for (BiScan other_scan:biCamDSI.biScans) if (other_scan.disabled_measurement[nTile] || remove_all_tried) {
							for (BiScan other_scan:biCamDSI.biScans) if ( (other_scan.list_index <= list_index) && (other_scan.disabled_measurement[nTile] || remove_all_tried)) {

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
						}
						if (valid_suggestion) {
							num_new.getAndIncrement();
							if (smooth_ds == null) { // fill target disparity
								target_disparity[nTile] = fit_rslt[1];
							} else { // this method is used to provide a filtered DSI - no changes to the target_disparity
								smooth_ds[0][nTile] = fit_rslt[1];
								double s0 = 0, s1 = 0;
								for (int dy = -smpl_radius; dy <= smpl_radius; dy++) {
									int ady = (dy > 0)? dy:(-dy);
									for (int dx = -smpl_radius; dx <= smpl_radius; dx++) { // if ((dx != 0) || (dy != 0)){
										int adx = (dx > 0)? dx:(-dx);
										int smpl_indx = smpl_center + dy*smpl_side + dx;
										if (smpl_wsw[smpl_indx] > 0.0) { // either smpl_w_narrow or smpl_w
											s0+= weights_sw[ady][adx];
//											s1+= weights[ady][adx] * smpl_w[smpl_indx];
											s1+= smpl_wsw[smpl_indx]; // already was multiplied by weights[ady][adx]
										}
									}
								}
								//boost_low_density, // false weight of low density tiles is reduced, true - boosted
								double w = s1/s0;
//								if (boost_low_density > 0.0) {
//									w/= Math.pow(s0, boost_low_density);
//								}
								if (strength_pow != 1.0) {
									w = Math.pow(w, 1.0 / strength_pow);
								}
								w += strength_floor;
								smooth_ds[1][nTile] = w;
								if (Double.isNaN(w)) {
									System.out.println("suggestNewScan(): nTile="+nTile+" w="+w);
									System.out.println("suggestNewScan(): nTile="+nTile+" w="+w);
								}
//								if (nTile == 66945) {
//									System.out.println("suggestNewScan(): nTile="+nTile+" w="+w);
//									System.out.println("suggestNewScan(): nTile="+nTile+" w="+w);
//								}

							}
						}
					}
				}
			};
		}
		ImageDtt.startAndJoin(threads);
		// remove duplicates from what was suggested or measured before
		return num_new.get();
	}



    private int findBetterFitToPlane(
    		int       smpl_radius,
    		int       nTile,
    		double    strength_floor,
    		double    strength_pow,
    		double [][] weights,
    		double [][] weights_narrow,
    		double [] smpl_w,
    		double [] smpl_w_narrow,
    		double [] smpl_d,
    		double [] smpl_p
    		)
    {
		final TileNeibs  tnImage = biCamDSI.tnImage;
		final int        smpl_side = 2 * smpl_radius + 1; // Sample size (side of a square)
//		final int        smpl_len = smpl_side*smpl_side;
		final int        smpl_center = (smpl_side + 1) * smpl_radius;
		int nnew = 0; // new assignments
		for (int dy = -smpl_radius; dy <= smpl_radius; dy++) {
			int ady = (dy > 0)? dy:(-dy);
			for (int dx = -smpl_radius; dx <= smpl_radius; dx++){
				int adx = (dx > 0)? dx:(-dx);
				int smpl_indx = smpl_center + dy*smpl_side + dx;
				if  (smpl_w [smpl_indx] > 0.0) { // is in "untouched" weights
					int nTile1 =  tnImage.getNeibIndex(nTile, dx, dy);
//					if ((nTile1 >= 0) && trusted[nTile1]) { // weak trusted OK
					if (nTile1 >= 0) { // will use any strength > 0
						// see if there is a measurement that fits better to the plane
						int best_indx=-1;
						double best_diff = Math.abs(smpl_d[smpl_indx] - smpl_p[smpl_indx]);
						for (int scan_indx = list_index; scan_indx >= 0; scan_indx--) { // include the latest measurement
							BiScan scan = biCamDSI.getBiScan(scan_indx);
							if (scan.disabled_measurement[nTile]) { //  || (scan.src_index[nTile] != indx)){ // skip all but enabled
								continue;
							}
							double disp = scan.disparity_measured[nTile1];
							if (Double.isNaN(disp)) {
								continue;
							}
							if (scan.strength_measured[nTile1] <= strength_floor) {
								continue;
							}
							double diff = Math.abs(disp - smpl_p[smpl_indx]);
							if ((best_indx < 0) || (diff < best_diff)){
								best_indx = scan_indx;
								best_diff = diff;
							}
						}
						if (best_indx >= 0) {
							BiScan scan = biCamDSI.getBiScan(best_indx);
							smpl_d[smpl_indx] = scan.disparity_measured[nTile1];
							double w =          scan.strength_measured[nTile1] -strength_floor ;
							if (w > 0.0) { // should be anyway
								if (strength_pow != 1.0) {
									w = Math.pow(w, strength_pow);
								}
								smpl_w[smpl_indx] =        w * weights[ady][adx];
								smpl_w_narrow[smpl_indx] = w * weights_narrow[ady][adx];
								nnew++;
							}
						}
					}
				}
			}
		}
		return nnew;
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
							for (BiScan other_scan:biCamDSI.biScans) if ( (other_scan.list_index <= list_index) && (other_scan.disabled_measurement[nTile] || remove_all_tried)) {
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
		final double goal_fraction_rms = 0.5;
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
		final boolean [][] halves = new boolean [8][smpl_len];
		{
			for (int i = 0; i < smpl_side; i++) {
				for (int j = 0; j < smpl_side; j++) {
					int indxs = i*smpl_side+ j;
					halves[0][indxs] = (i <= smpl_radius);
					halves[4][indxs] = (i >= smpl_radius);
					halves[2][indxs] = (j >= smpl_radius);
					halves[6][indxs] = (j <= smpl_radius);
					halves[1][indxs] = (j >= i);
					halves[5][indxs] = (j <= i);
					halves[3][indxs] = ((i + j) >= (2 * smpl_radius));
					halves[7][indxs] = ((i + j) <= (2 * smpl_radius));
				}
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
						boolean [] smpl_trusted = new boolean [smpl_len];
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
									int adx = (dx > 0)? dx:(-dx);
									double max_fdiff = max_diff + (ady+adx) * max_tilt;
									int smpl_indx = smpl_center + dy*smpl_side + dx;
									smpl_trusted[smpl_indx] = true;
									nall++;
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
						double max_rms = smpl_arms + smpl_rrms * disp_center;
						if (nsmpls >= fin_samples) {
							double [] smpl_w_persistent = smpl_w.clone();
							double [] fit_rslt = fitPlaneRemoveOutliers(
									smpl_radius, // int                     smpl_radius,
									max_tilt,    // double                  max_tilt,
									damp_tilt,   // double                  damp_tilt, //   =     0.001; // Tilt cost for damping insufficient plane data
									false,       // boolean                 full_plane,
									smpl_d,      // double []               smpl_d,
									smpl_w,      // double []               smpl_w, // will be modified,
									smpl_p,      // double []               smpl_p, // will be set if provided
									fin_samples, // int                     fin_samples, // remove until this number remain
									goal_fraction_rms*max_rms,    // double                  fin_rms,
									true,        // boolean                 keep_center, // do not remove center tile - it is the tile that should be verified by neighbors
									0,           // int                     fourq_min,         // each of the 4 corners should have at least this number of tiles.
									null,        // int [] 	                fourq_corner, //  array specifying corner number (0..3), -1 - gap. null when not used
									debugLevel); // int                     debugLevel)
//							if ( (fit_rslt != null) && (fit_rslt[0] < (smpl_arms + smpl_rrms * fit_rslt[1]))){
							if ( (fit_rslt != null) && (fit_rslt[0] < max_rms)){
								// Trusted tile, save it
								trusted[nTile] = true;
								num_trusted_plane.getAndIncrement();
								continue;
							}
							// try 8 "halves" around the tile - it may be weak background close to the strong foreground (or stray FG tile)
							for (int dir = 0; dir < halves.length; dir++) {
								nsmpls = 0;
								nall = 0;
								smpl_w = new double [smpl_len];
								for (int indxs = 0; indxs < smpl_len; indxs++) if (halves[dir][indxs]) {
									if (smpl_trusted[indxs]) nall++;
									if (smpl_w_persistent[indxs] > 0.0) {
										smpl_w[indxs] = smpl_w_persistent[indxs];
										nsmpls++;
									}
								}
								fin_samples= (int) ( nall * smpl_fract);
								if (fin_samples < smpl_num) fin_samples = smpl_num;
								if (nsmpls >= fin_samples) {
									fit_rslt = fitPlaneRemoveOutliers(
											smpl_radius, // int                     smpl_radius,
											max_tilt,    // double                  max_tilt,
											damp_tilt,   // double                  damp_tilt, //   =     0.001; // Tilt cost for damping insufficient plane data
											false,       // boolean                 full_plane,
											smpl_d,      // double []               smpl_d,
											smpl_w,      // double []               smpl_w, // will be modified,
											smpl_p,      // double []               smpl_p, // will be set if provided
											fin_samples, // int                     fin_samples, // remove until this number remain
											goal_fraction_rms*max_rms,    // double                  fin_rms,
											true,        // boolean                 keep_center, // do not remove center tile - it is the tile that should be verified by neighbors
											0,           // int                     fourq_min,         // each of the 4 corners should have at least this number of tiles.
											null,        // int [] 	                fourq_corner, //  array specifying corner number (0..3), -1 - gap. null when not used
											debugLevel); // int                     debugLevel)
									if ((fit_rslt != null) && (fit_rslt[0] < max_rms)){
										// Trusted tile, save it
										trusted[nTile] = true;
										num_trusted_plane.getAndIncrement();
										break; // from for (int dir = 0; dir < halves.length; dir++) {
									}
								}
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

	/**
	 * Remove stray tiles that are closer than closest neighbor and weaker than it or
	 *  trusted_strength * min_rstrength
	 * @param trusted_strength absolute raw strength to trust
	 * @param min_rstrength minimal strength to allow lone FG, relative to trusted_strength
	 * @param max_disp_inc maximal disparity difference between this tile and the nearest neighbor
	 * @param dbg_x
	 * @param dbg_y
	 * @param debugLevel
	 * @return number of disabled tiles
	 */

	int  trimWeakLoneFG(
		    final double     trusted_strength, // trusted correlation strength
			final double     min_rstrength,   // strength floor - relative to trusted
			final double     max_disp_inc,
			final int        dbg_x,
			final int        dbg_y,
			final int        debugLevel
			) {
		final AtomicInteger ai_trimmed = new AtomicInteger(0);
		final double min_strength = trusted_strength * min_rstrength;
		final TileNeibs  tnImage = biCamDSI.tnImage;
		final int num_tiles = tnImage.sizeX * tnImage.sizeY;
		final Thread[] threads = ImageDtt.newThreadArray(biCamDSI.threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		final double [][] ds = getDisparityStrength( // already has disabled zeroed
	    		false,   // final boolean only_strong,
	    		false,   // final boolean only_trusted,
	    		true) ;  // final boolean only_enabled);

		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				@Override
				public void run() {
					for (int nTile = ai.getAndIncrement(); nTile < num_tiles; nTile = ai.getAndIncrement()) if (!Double.isNaN(ds[0][nTile])){
						double max_disp = 0;
						double max_disp_w = 0;
						double d_lim = ds[0][nTile] - max_disp_inc;
						double w =     ds[1][nTile];
						for (int dir = 0; dir < 8; dir++) {
							int nTile1 = tnImage.getNeibIndex(nTile, dir);
							if ((nTile1 >=0) && (ds[0][nTile1] > max_disp)){
								 max_disp =   ds[0][nTile1];
								 max_disp_w = ds[1][nTile1];
							}
							if (max_disp > d_lim) {
								break;
							}
						}
						if ((max_disp <= d_lim) && ((w < max_disp_w) || (w < min_strength))) {
							disableTile(nTile);
							ai_trimmed.getAndIncrement();
							if (debugLevel > -4) {
								System.out.println("trimWeakLoneFG: removing tile "+nTile+" ("+(nTile%tnImage.sizeX)+":"+(nTile/tnImage.sizeX));
							}
						}
					}
				}
			};
		}
		ImageDtt.startAndJoin(threads);
		return ai_trimmed.get();
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
		final double goal_fraction_rms = 0.5;
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
							double max_rms = smpl_arms + smpl_rrms * disp_center;
//		final double goal_fraction_rms = 0.5;

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
										goal_fraction_rms*max_rms,     // double                  fin_rms,
										true,        // boolean                 keep_center, // do not remove center tile - it is the tile that should be verified by neighbors
										0,           // int                     fourq_min,         // each of the 4 corners should have at least this number of tiles.
										null,        // int [] 	               fourq_corner, //  array specifying corner number (0..3), -1 - gap. null when not used
										debugLevel); // int                     debugLevel)
//								if ( (fit_rslt != null) && (fit_rslt[0] < (smpl_arms + smpl_rrms * fit_rslt[1]))){
								if ( (fit_rslt != null) && (fit_rslt[0] < max_rms)){
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
				biCamDSI.getLastBiScan(BISCAN_ANY).showScan(
						"trimWeakFG_"+npass);
			}

			total_disabled += num_disabled.get();
		}
		return total_disabled; // numTrustedAll;
	}



//					smpl_p[indxs] = approx2d[0][0] * (sx - smpl_radius) + approx2d[0][1] * (sy - smpl_radius) + approx2d[0][2];

	/**
	 * Remove worst outliers when fitting planes until the specified number of samples are left, re-fitting
	 * plane after each iteration.
	 * Additional mode for filling gaps (not extending) - "four quadrants" is activated when  fourq_min > 0
	 * In this mode each of the four corners of the sample square (after removing center 2*fourq_gap -1 rows
	 * and columns should have at least fourq_min remaining tiles
	 * @param smpl_radius sample "radius", square side is  2 * smpl_radius + 1
	 * @param max_tilt maximal DSI plane tilt (pixels/tile)
	 * @param damp_tilt tilt cost for damping insufficient plane data
	 * @param full_plane generate full tilted plane, not only for non-zero tiles
	 * @param smpl_d  array of disparity values to approximate
	 * @param smpl_w  weights - only >0.0 are processed
	 * @param smpl_p approximated disparity values, may be null. When provided (not null), will have calculated disparity approximation
	 * @param fin_samples remove until this number of tiles remain
	 * @param fin_rms - OK RMS - exit if it is below (set to 0.0 if unknown yet)
	 * @param keep_center do not remove center tile - it is the tile that should be verified by neighbors
	 * @param fourq_min each of the 4 corners should have at least this number of tiles.
	 * @param fourq_corner array specifying corner number (0..3), -1 - gap. null when not used
	 * @param debugLevel debug level
	 * @return array of 4 values: rms of the remaining samples,  average (center) disparity value, tiltX, tiltY
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
			double                  fin_rms,
			boolean                 keep_center, // do not remove center tile - it is the tile that should be verified by neighbors
			int                     fourq_min,
			int []                  fourq_corner,
			int                     debugLevel)
	{
		int [] num_in_corners = ((fourq_corner == null) || (fourq_min == 0)) ? null: new int [5]; // last is not used
		final int keep_tile = keep_center? (2*smpl_radius*(smpl_radius+1)): -1; // do not remove tile with this number
		double [] disp_tilts = null;
		final int        smpl_side = 2 * smpl_radius + 1; // Sample size (side of a square)
		int num_in_sample = 0;
		for (int i = 0; i < smpl_w.length; i++) if (smpl_w[i]> 0.0) {
			num_in_sample++;
			if (num_in_corners != null) num_in_corners[fourq_corner[i]]++;
		}

		int smpl_len = smpl_w.length;
		if (num_in_sample < fin_samples) {
			return null;
		}
		if (num_in_corners != null) {
			for (int i = 0; i < 4; i++) {
				if (num_in_corners[i] < fourq_min) {
					return null;
				}
			}
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
			if ((num_in_sample > fin_samples) && !(rms < fin_rms)) { // remove worst if it is not the last run where only calculations are needed
				int iworst = -1;
				double dworst2 = 0.0;
				for (int indxs = 0; indxs < smpl_len; indxs++) if (smpl_w[indxs] > 0.0) {
					if (indxs == keep_tile) { // it was the center tile disabled from removing
						continue;
					}
					// verify that the worst tile will not leave the corner too empty
					int num_corner = 4; // unused, gap
					if (num_in_corners != null) {
						num_corner =  fourq_corner[indxs];
						if ((num_corner < 4) && (num_in_corners[num_corner] <= fourq_min)) {
							continue; // can not remove this tile
						}
					}

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
				if (iworst < 0){ // can happen if damp_tilt == 0.0 and will become colinear
					if (debugLevel > 0) {
						System.out.println("**** this may be BUG in fitPlaneRemoveOutliers() can not find the worst sample  - all tiles fit perfectly ****");
					}
					// this can happen if some samples are the same and all the pixels fit exactly - use all of them
					break;
				}
				// remove worst sample
				if (num_in_corners != null) {
					num_in_corners[fourq_corner[iworst]]--; // corner counter - can remove from the gap, but nobody count that
				}
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

	public double [] getDensity(
			final boolean strong_only,
			final int need_tiles,
			final int max_radius,
			final int dbg_x,
			final int dbg_y,
			final int debugLevel )
	{
		final int [][][] incr = {
				{{ 1,-1},{ 0, 1}},  // top right corner, going down
				{{ 1, 1},{-1, 0}},  // bottom right corner, going left
				{{-1, 1},{ 0,-1}},  // bottom left corner, going up
				{{-1,-1},{ 1, 0}}}; // top left corner, going right

		final TileNeibs  tnImage = biCamDSI.tnImage;
		final int dbg_tile = (debugLevel > -2)? (dbg_y * tnImage.sizeX + dbg_x):-1;
		final boolean [] en = strong_only ? strong_trusted: trusted;
    	final int num_tiles = en.length;
    	final double [] density = new double[num_tiles];
		final Thread[] threads = ImageDtt.newThreadArray(biCamDSI.threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				@Override
				public void run() {
					for (int nTile = ai.getAndIncrement(); nTile < num_tiles; nTile = ai.getAndIncrement()) {
						if (nTile == dbg_tile) {
							System.out.println("getDensity(): nTile="+nTile);
						}
//						int num_tried = 0, num_found =0;
						int sd = 0, sdr2 = 0, s0=0, sr2 = 0;
						int dx = 0,dy = 0;
						int r = 0;
						label_snail: {
							for (; r <= max_radius; r++) {
								int nd = (r == 0) ? 1: 4;
								int nl = (r == 0) ? 1: (2 * r);
								for (int dir = 0; dir < nd; dir ++) {
									for (int l = 0; l < nl; l++) {
										dx = r * incr[dir][0][0] + l * incr[dir][1][0];
										dy = r * incr[dir][0][1] + l * incr[dir][1][1];
										int nTile1 = tnImage.getNeibIndex(nTile, dx,dy);
										if (nTile1 >= 0) {
											int r2 = dx*dx + dy*dy;
											s0++;
											sr2 += r2;
											if (en[nTile1]) {
												sd++;
												sdr2 += r2;
											}
//											num_found++;
											if (sd >= need_tiles) {
												break label_snail;
											}
										}
									}
								}
							}
						}
						int r02 = 2 * r * r;
						int num = r02 * sd - sdr2;
						int denom = r02 * s0 - sr2;
//						int num = sd;
//						int denom = s0;

						if (denom > 0) {
							density[nTile] = (1.0* num) / denom;
						}
					}
				}
			};
		}
		ImageDtt.startAndJoin(threads);
		return density;
	}

	/**
	 * Select low-textured tiles for averaging measurements
	 * @param min_disparity minimal disparity to accept
	 * @param max_density maximal trusted tile density (density varies from 0.0 to 1.0)
	 * @param grow how many layers of tiles should be added after filling gaps and removing small clusters
	 * @param max_gap_radius maximal radius of a void to be filled
	 * @param min_clust_radius minimal original cluster radius to survive
	 * @param density per-tile values of the density of trusted tiles around it.
	 * @param src_disparity - source disparity array. If null will only use density (that should be > 0)
	 * @return selection of the low-textured tiles to be processed with averaging correlation (3x3 or 5x5 tiles)
	 */
	public boolean [] selectLowTextures(
			double    min_disparity,
			double    max_density,
			int       grow,
			int       max_gap_radius,
			int       min_clust_radius,
			double [] density,
			double [] src_disparity)
	{
		boolean [] selection = new boolean [density.length];
		if (src_disparity == null) {
			for (int nTile = 0; nTile < selection.length; nTile++) {
				if ((density[nTile] <= max_density) && (density[nTile] <= max_density)) { // disparity has NaN-s, they will fail comparisons
					selection[nTile] = true;
				}
			}
		} else {
			for (int nTile = 0; nTile < selection.length; nTile++) {
				if ((src_disparity[nTile] >= min_disparity) && (density[nTile] <= max_density)) { // disparity has NaN-s, they will fail comparisons
					selection[nTile] = true;
				}
			}
		}
		final TileNeibs  tnImage = biCamDSI.tnImage;

		tnImage.growSelection(
				2* max_gap_radius, // int        grow,           // grow tile selection by 1 over non-background tiles 1: 4 directions, 2 - 8 directions, 3 - 8 by 1, 4 by 1 more
				selection, // boolean [] tiles,
				null); // boolean [] prohibit)

		tnImage.shrinkSelection(
				2*(max_gap_radius + min_clust_radius), // int        grow,           // grow tile selection by 1 over non-background tiles 1: 4 directions, 2 - 8 directions, 3 - 8 by 1, 4 by 1 more
				selection, // boolean [] tiles,
				null); // boolean [] prohibit)
		tnImage.growSelection(
				2 * (min_clust_radius + grow), // int        grow,           // grow tile selection by 1 over non-background tiles 1: 4 directions, 2 - 8 directions, 3 - 8 by 1, 4 by 1 more
				selection, // boolean [] tiles,
				null); // boolean [] prohibit)
		return selection;
	}


	/**
	 * Fill the gaps (where provided src_disparity is NaN) and optionally smooth it using pull to the weighted average of 8 neighbors
	 * @param src_disparity source disparity values, some may be undefined (Double.NaN)
	 * @param src_strength optional strengths of the initial values (should be with floor subtracted)
	 * @param selection tile selection to process
	 * @param neib_pull pull of the average neighbor weight relative to the original disparity value of a tile. If 0 - only gaps
	 *  (Double.NaN) values are filled, all the defined disparities remain as they were provided
	 * @param max_iterations Maximal number of the iteration steps
	 * @param min_change exit iterations when the maximal disparity change to a tile is less than this value
	 * @param dbg_x debug tile X coordinate
	 * @param dbg_y debug tile Y coordinate
	 * @param debugLevel debug level
	 * @return processed disparity and optional strength array. Normally only unselected tiles should remain Double.NaN, all other should be interpolated
	 */
	public double [][] fillAndSmooth(
			final double [] src_disparity,
			final double [] src_strength, // if not null will be used for weighted pull
			final boolean [] selection,
			final double     neib_pull, // pull to weighted average relative to pull to the original disparity value. If 0.0 - will only update former NaN-s
			final int max_iterations,
			final double min_change,
			final int dbg_x,
			final int dbg_y,
			final int debugLevel) {
		final double [] weights = {1.0, 0.7}; // {ortho, corners},
		final boolean adjust_all = (neib_pull > 0.0);
//		final double fneib_pull = adjust_all ? neib_pull: 1.0;
		final TileNeibs  tnImage = biCamDSI.tnImage;
//		final int dbg_tile = (debugLevel > -2)? (dbg_y * tnImage.sizeX + dbg_x):-1;
		final int dbg_tile = (debugLevel > 0)? (dbg_y * tnImage.sizeX + dbg_x):-1;
    	final int num_tiles = src_disparity.length;
		final Thread[] threads = ImageDtt.newThreadArray(biCamDSI.threadsMax);
		// max_changes may have Double.NaN value (here meaning larger than any)
		final double [] max_changes = new double [biCamDSI.threadsMax]; // Each thread provides its own maximal change, then they are combined
		final double [] disparity = src_disparity.clone();
		final double [] new_disparity = src_disparity.clone();
		final double [] strength = (src_strength != null)? src_strength: null; // new double [src_disparity.length];
		final double [] new_strength = (src_strength != null)? src_disparity.clone(): null;
//		if (src_strength == null) for (int i = 0; i < strength.length; i++) if (selection[i] && !Double.isNaN(src_disparity[i])) strength[i] = 1.0;
		final AtomicInteger ai = new AtomicInteger(0);
		final AtomicInteger ai_numThread = new AtomicInteger(0);
		for (int num_iter = 0; num_iter < max_iterations; num_iter++) {
			ai.set(0);
			ai_numThread.set(0);
			final AtomicInteger ai_count=new AtomicInteger(0);
			for (int i = 0; i < max_changes.length; i++) max_changes[i] = 0.0;
			final int fnum_iter = num_iter;
			for (int ithread = 0; ithread < threads.length; ithread++) {
				threads[ithread] = new Thread() {
					@Override
					public void run() {
						int numThread = ai_numThread.getAndIncrement(); // unique number of thread to write to max_changes[numThread]
//						max_changes[numThread]
						for (int nTile = ai.getAndIncrement(); nTile < num_tiles; nTile = ai.getAndIncrement()) if (selection[nTile]){
							if (nTile == dbg_tile) {
								System.out.println("fillAndSmooth(): iteration "+fnum_iter+" nTile="+nTile);
							}
							if (!adjust_all && !Double.isNaN(src_disparity[nTile])) {
// nothing to do - new_disparity and new_strength are already same as new_*
							} else {
								double sw=0.0, swd =0.0, sww = 0.0;


								for (int dir = 0; dir < 8; dir++) {
									int nTile1 = tnImage.getNeibIndex(nTile, dir);
									if (nTile1 >= 0) {
										double w = weights[dir & 1];
										double s = (strength == null) ? 1.0 : (Math.max(strength[nTile1], 0.0));
										double d = disparity[nTile1];
										double ww = w * s;
										if ((ww > 0.0) && !Double.isNaN(d)) {
											sw += w;
											sww += ww;
											swd += ww * d;
										}
									}
								}
								if (sww > 0.0) { // then sw
									double d_mean = swd/sww;
									double w_mean = sww/sw;
									if (Double.isNaN(src_disparity[nTile])) {
										new_disparity[nTile] =  d_mean;
										if (strength != null) {
											new_strength[nTile] = w_mean;
										}
									} else {
										double pull_origin = (strength==null)? 1.0: strength[nTile];
										new_disparity[nTile] =  (d_mean *(neib_pull * w_mean) + src_disparity[nTile] * pull_origin)/((neib_pull * w_mean) + pull_origin);
										if (strength != null) {
											new_strength[nTile] =  (w_mean * neib_pull + src_strength[nTile])/(neib_pull + 1);
										}
									}
									ai_count.getAndIncrement();
									double adiff = Math.abs(new_disparity[nTile] - disparity[nTile]); // disparity[nTile] may be NaN then adiff will be NaN as intended
									if (!(adiff < max_changes[numThread])) {
										max_changes[numThread] = adiff; // NaN will be copied
									}
								}
							}
						}
					}
				};
			}
			ImageDtt.startAndJoin(threads);
			System.arraycopy(new_disparity, 0, disparity, 0, num_tiles);
			if (strength != null) {
				System.arraycopy(new_strength, 0, strength, 0, num_tiles);
			}
			double change = 0.0;
			for (int i = 0; i < ai_numThread.get(); i++) {
				if (Double.isNaN(max_changes[i]) || (change < max_changes[i])) { // max_changes[i] may be NaN
					change = max_changes[i];
					if (!(change <= min_change)) { // change may be NaN
						if (debugLevel < -2) {
							// may exit here, continue to get debug info
							break;
						}
					}
				}
			}
			if (debugLevel > 0) {
				System.out.println("fillAndSmooth(): iteration "+fnum_iter+" change="+change+" (min_change="+min_change+")+ tiles updated="+ai_count.get());
			}
			if (change <= min_change) { // change may be NaN
				break; // from the main loop
			}
		}
		double [][] ds = {disparity,strength}; // strength may be null
		return ds;
	}
	/**
	 * Extend low texture areas horizontally if both ends are OK: either belong to the low texture areas selection (lt_sel)
	 * or are trusted and closer or the same (to the provided tolerance)
	 * @param tolerance strong limit should not have disparity lower by more than tolerance than low textured area
	 * @param ds_lt disparity/strength for the low textured area
	 * @param d_single disparity measured for the single-tile correlation
	 * @param lt_sel low -textured selection
	 * @param exp_sel expanded selection (does not intersect with lt_sel
	 * @param trusted trusted tiles selection
	 * @return extended disparity/strength data
	 */
	  public double [][] getLTExpanded(
			  final double      tolerance,
			  final double [][] ds_lt,
			  final double []   d_single,
			  final boolean []  lt_sel, // lt_sel and exp_sel do not intersect
			  final boolean []  exp_sel,
			  final boolean []  trusted)
	  {
//		final int dbg_tile = (debugLevel>-2)?(dbg_x + tnImage.sizeX*dbg_y):-1;
		  final int num_tiles = exp_sel.length;
		  final double [][] ds = new double [2][num_tiles];
		  for (int i = 0; i < num_tiles; i++) ds[0][i] = Double.NaN;
		  final Thread[] threads = ImageDtt.newThreadArray(biCamDSI.threadsMax);
		  final AtomicInteger ai = new AtomicInteger(0);
		  final TileNeibs         tnImage = biCamDSI.tnImage;
		  for (int ithread = 0; ithread < threads.length; ithread++) {
			  threads[ithread] = new Thread() {
				  @Override
				  public void run() {
					  //						max_changes[numThread]
					  for (int nTile = ai.getAndIncrement(); nTile < num_tiles; nTile = ai.getAndIncrement())  if (lt_sel[nTile]) {
						  // is low texture = - just copy
						  ds[0][nTile] = ds_lt[0][nTile];
						  ds[1][nTile] = ds_lt[1][nTile];
					  } else if (exp_sel[nTile]) {
						  int dbg_tileX = nTile%tnImage.sizeX;
						  int dbg_tileY = nTile/tnImage.sizeX;
//						  if ((dbg_tileY == 156) || (dbg_tileY == 157)) {
//							  System.out.println("getLTExpanded(): tileX="+dbg_tileX+", tileY="+dbg_tileY);
//							  System.out.println("getLTExpanded(): tileX="+dbg_tileX+", tileY="+dbg_tileY);
//						  }
						  int nTile0= tnImage.getNeibIndex(nTile,-1,0);
						  if ((nTile0 < 0) || !exp_sel[nTile0]){
							  boolean OK0 = (nTile0 < 0) || lt_sel[nTile0] || (trusted[nTile0] && (d_single[nTile0] >= (ds_lt[0][nTile] - tolerance)));
							  if (OK0) {
								  int nTile1= tnImage.getNeibIndex(nTile,1,0);
								  int l = 1;
								  while ((nTile1 >= 0) && exp_sel[nTile1]) {
									  nTile1= tnImage.getNeibIndex(nTile1,1,0);
									  l++;
								  }

								  boolean OK1 = (nTile1 < 0) || lt_sel[nTile1] || (trusted[nTile1] && (d_single[nTile1] >= (ds_lt[0][nTile1-1] - tolerance)));
								  if (OK1) {
									  for (int i = 0; i < l; i++) {
										  int nt = nTile+ i;
										  ds[0][nt] = ds_lt[0][nt];
										  ds[1][nt] = ds_lt[1][nt];
									  }
								  }
							  }
						  }
					  }
				  }
			  };
		  }
		  ImageDtt.startAndJoin(threads);
		  return ds;

	  }


}
