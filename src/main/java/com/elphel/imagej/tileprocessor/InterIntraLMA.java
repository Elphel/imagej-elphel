package com.elphel.imagej.tileprocessor;

import java.util.Arrays;
import java.util.concurrent.atomic.AtomicInteger;

import com.elphel.imagej.common.ShowDoubleFloatArrays;

import Jama.Matrix;

/**
 **
 ** InterIntraLMA - Comparing contrast for different interscene and intrascene
 ** connfifurations
 **
 ** Copyright (C) 2021 Elphel, Inc.
 **
 ** -----------------------------------------------------------------------------**
 **
 **  InterIntraLMA.java is free software: you can redistribute it and/or modify
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

public class InterIntraLMA {
	public static int threadsMax = 100; // maximal number of threads to use
	/**
	 * Estimate noise threshold for each combination of inter/intra and number of
	 * sensors.
	 * @param noise_file synthetic noise level per file (either rnd or fpn)
	 * @param senosor_mode_file sensor mode (0 - 16, 1 - 8, 2 - 4, 3 - 2 sensors) per file
	 * @param inter_file true - interscene , false - intrascene only
	 * @param good_file_tile true if tile of this file is considered "good"
	 * @param min_modes minimal number of modes the tile should be defined for (undefine if less than)
	 * @return per tile per mode average between highest noise of the good tile and lowest noise
	 * of the bad tile. Double.NaN if noise threshold can not be determined. null if the tile is
	 * undefined for all modes   
	 */
	public static double [][] getNoiseThreshold(
			double []    noise_file, //  = new double [noise_files.length];
			int []       sensor_mode_file,
			boolean []   inter_file,
			boolean [][] good_file_tile,
			double       min_inter16_noise_level,
			int          min_modes)
	{
		int dbg_tile = 1222;
//		boolean remove_non_monotonic = false; // true;
		boolean zero_all_bad = true;    // set noise_level to zero if all noise levels result in bad tiles
		boolean all_inter =    true;    // tile has to be defined for all inter
		boolean need_same_inter = true; // do not use intra sample if same inter is bad for all noise levels  
		int num_sensor_modes = 0;
		int num_tiles = good_file_tile[0].length;
		for (int i = 0; i < sensor_mode_file.length; i++) {
			if (sensor_mode_file[i] > num_sensor_modes) {
				num_sensor_modes = sensor_mode_file[i];
			}
		}
		num_sensor_modes ++;
		int num_modes = 2 * num_sensor_modes;
		double [][] rslt = new double [num_tiles][]; // number of tiles  
		double [][][] noise_interval = new double[num_modes][num_tiles][2]; // [modes][tiles] {max_good, min_bad}
		for (int i = 0; i < num_modes; i++) {
			for (int j= 0; j < num_tiles; j++) {
				noise_interval[i][j][0] =Double.NaN;
				noise_interval[i][j][1] =Double.NaN;
			}
		}
		for (int nf = 0; nf < noise_file.length; nf++) {
			double noise = noise_file[nf];
			int mode = sensor_mode_file[nf] + (inter_file[nf] ? 0: num_sensor_modes); 
			for (int ntile = 0; ntile < num_tiles; ntile++) {
				if (ntile == dbg_tile) {
					System.out.println("ntile = "+ntile+", nf ="+nf);
				}
				if (good_file_tile[nf][ntile]) { // good tile
					if (!(noise <= noise_interval[mode][ntile][0])){ // including Double.isNaN(noise_interval[mode][ntile][0]
						noise_interval[mode][ntile][0] = noise;
					}
				} else { // bad tile
					if (!(noise >= noise_interval[mode][ntile][1])){ // including Double.isNaN(noise_interval[mode][ntile][1]
						noise_interval[mode][ntile][1] = noise;
					}
				}
			}
		}
		
		for (int ntile = 0; ntile < num_tiles; ntile++){ 
			if (ntile == dbg_tile) {
				System.out.println("ntile = "+ntile);
			}
			
			int num_defined = 0;
			int num_defined_inter = 0;
			for (int mode = 0; mode < num_modes; mode++) {
				if (!Double.isNaN(noise_interval[mode][ntile][0]) && !Double.isNaN(noise_interval[mode][ntile][1])) {
					num_defined++;
					if (mode < 4) {
						num_defined_inter++;
					}
				}
			}
			//all_inter
			if ((num_defined >= min_modes) && (!all_inter || (num_defined_inter >= 4))) {
				rslt[ntile] = new double [num_modes];
				for (int mode = 0; mode < num_modes; mode++) {
//					if (need_same_inter && (mode >= 4) && Double.isNaN(noise_interval[mode & 3][ntile][0])) {  // no good for same sensors inter
					if (need_same_inter && Double.isNaN(noise_interval[mode & 3][ntile][0])) { // no good for same sensors inter
						rslt[ntile][mode] = Double.NaN;
					} else 	if (!Double.isNaN(noise_interval[mode][ntile][0]) && !Double.isNaN(noise_interval[mode][ntile][1])) {
						/*
						rslt[ntile][mode] = 0.5 * (noise_interval[mode][ntile][0] + noise_interval[mode][ntile][1]);
						if (remove_non_monotonic && (noise_interval[mode][ntile][0] > noise_interval[mode][ntile][1])) {
							rslt[ntile][mode] = Double.NaN;
						}
						*/
						// use the lowest failed noise level assuming that false positive may happen even for much higher noise level 
						rslt[ntile][mode] = noise_interval[mode][ntile][1]; // lowest noise for bad
						
//					} else if (zero_all_bad && Double.isNaN(noise_interval[mode][ntile][1])) {
					} else if (zero_all_bad && Double.isNaN(noise_interval[mode][ntile][0])) {
						rslt[ntile][mode] = 0.0;
					} else {
						rslt[ntile][mode] = Double.NaN;
					}
				}
			}
			if ((rslt[ntile] != null) && (min_inter16_noise_level >0)){ // filter by to weak inter-16 (mode 0)
				if (!(rslt[ntile][0] >= min_inter16_noise_level)){
					rslt[ntile] = null;
				}
			}
			if (rslt[ntile] != null) {
				boolean all_nan = true;
				boolean has_nan = false;
				boolean has_inter_nan = false;
				for (int mode = 0; mode < rslt[ntile].length; mode++) {
					if (Double.isNaN(rslt[ntile][mode])) {
						has_nan = true;
						if (mode < 4) {
							has_inter_nan = true;
						}
					} else {
						all_nan = false;
					}
				}
				if (all_nan) {
					System.out.println("All NaN for tile = "+ntile);
				}
				if (has_nan) {
					System.out.println("Has NaN for tile = "+ntile);
				}
				if (has_inter_nan) {
					System.out.println("Has has_inter_nan for tile = "+ntile);
				}
			}
		}
		return rslt;
	}
	//Monotonic function
	public int       debug_level = 0;
	public double    offset;
	public double    lambda = 1.0;
	public double    lambda_good = 0.5;
	public double    lambda_bad =  8.0;
//	public double    rms =  Double.NaN;
//	public double    rms0 = Double.NaN;
	public double [] vector; // N0, g[1]... [g7], St[i]
	public int [][]  sample_indx; // pairs of {tile_index, mode}
	public double [] gi;
	
	public double [][] last_jt;// 			this.last_jt = new double [num_pars][num_points];

	public double    N0;
	public double [] Y;
	public double [] K; // scale noise levels to make them near-relative
	public int    [] tile_index;
	public double [] St;
	public double [] weights;
	public double [] fx;
	public double    last_rms = Double.NaN;
	public double    initial_rms = Double.NaN;
	public double    good_or_bad_rms = Double.NaN;
	double []        last_ymfx;
	
	boolean adjust_N0 = true;
	boolean adjust_Gi = true;
	boolean adjust_St = true;
	
	int dbgTilesX = 80;
	int dbgTilesY = 64;
	
	public InterIntraLMA(
			double [][] noise_thresh,
			double      offset, // initial value for N0
			int         tilesX, // debug images only
			int         debug_level)
	{
		boolean debug_img =  (debug_level > -1);
//		this.gi =     g0.clone();
		this.offset = offset;
		this.debug_level = debug_level;
		int num_samples = 0;
		int num_tiles = 0;
		int num_modes = 0;
		for (double [] sample : noise_thresh) if (sample != null) {
			for (double d:sample) if (!Double.isNaN(d)) {
				num_samples++;
			}
			num_tiles++;
			num_modes = sample.length;
		}
		this.gi = new double[num_modes];
		this.gi[0] = 1.0; // all, inter - ga1n= 1.0
		sample_indx = new int [num_samples][2];
		tile_index = new int[num_tiles];
		N0 = 0.03; // offset; // .01; // offset;
		Y = new double[num_samples];
		K = new double[num_samples];
		weights = new double[num_samples];
		St = new double [num_tiles];
		Arrays.fill(St, Double.NaN);
		last_rms = Double.NaN;

		int indx = 0;
		int t_indx = 0;
		for (int ntile = 0; ntile < noise_thresh.length; ntile++) if (noise_thresh[ntile] != null) {
			tile_index[t_indx] = ntile;
			for (int mode = 0; mode < noise_thresh[ntile].length; mode++) {
				if (!Double.isNaN(noise_thresh[ntile][mode])) {
					sample_indx[indx][0] = t_indx; // ntile;
					sample_indx[indx][1] = mode;
					indx++;
				}
			}
			t_indx++;
		}
		
		
		// create Y, K and weights vectors
		for (int nsample = 0; nsample < num_samples; nsample++) {
			int tile = tile_index[sample_indx[nsample][0]];
			double d = noise_thresh[tile][sample_indx[nsample][1]];
			K[nsample] = 1.0/(d + offset);
			Y[nsample] = d * K[nsample];
			// may be modified, but sum (weights) should be == 1.0;
			weights[nsample] = 1.0/num_samples;
		}
		// initial approximation
		double N0 =  offset;
		double N02 = N0*N0;
		// set St for tiles that are defined for mode==0 (inter16) 
		for (int nsample = 0; nsample < num_samples; nsample++) if (sample_indx[nsample][1] == 0){
			double sqrt2 = Y[nsample]/K[nsample];
			sqrt2 *= sqrt2;
			double st2 = sqrt2 - N02;
			int ntile = sample_indx[nsample][0]; // used tile index
			if (st2 > 0.0) {
				St[ntile] = Math.sqrt(st2);	
			}
		}
		// iteratively find gi (other than 0), then refine St
		int num_initial_refines = 10; //5;
		double [][] dbg_img = debug_img ? (new double [num_initial_refines+1][]) : null;
		double [][] dbg_ratio = debug_img ? (new double [num_modes][noise_thresh.length]) : null;
		
		if (dbg_img != null) {
			double [] dbg_st = new double [noise_thresh.length];
			Arrays.fill(dbg_st,Double.NaN);
			for (int ntile = 0; ntile < tile_index.length; ntile++) {
				dbg_st[tile_index[ntile]] = St[ntile];
			}
			dbg_img[0] = dbg_st;
			for (int i = 0; i < dbg_ratio.length; i++) {
				Arrays.fill(dbg_ratio[i], Double.NaN);
			}
			for (int ntile = 0; ntile < noise_thresh.length; ntile++) if (noise_thresh[ntile] != null) {
				double [] noise_per_mode = noise_thresh[ntile];
				for (int mode = 0; mode < num_modes; mode++) {
					int mode_ref = (mode > 4) ? 4 : 0;
					if (!Double.isNaN(noise_per_mode[mode]) && !Double.isNaN(noise_per_mode[mode_ref])) {
						dbg_ratio[mode][ntile] = noise_per_mode[mode]/noise_per_mode[mode_ref];
					}
				}
			}
			String [] dbg_ratio_titles = {"0/0","1/0","2/0", "3/0", "4/0","5/4","6/4","7/4"};
			(new ShowDoubleFloatArrays()).showArrays(
					dbg_ratio,
					tilesX,
					noise_thresh.length/ tilesX,
					true,
					"dbg_ratio",
					dbg_ratio_titles);
		}
		dbgTilesX = tilesX;
		dbgTilesY = noise_thresh.length/ tilesX;
		
		for (int nrefine = 0; nrefine < num_initial_refines; nrefine++) {
			// finding g[i]
			gi = new double[num_modes];
//			gi[0] = 1.0; // all, inter - gain= 1.0
			double [] wgi = new double [num_modes]; // weights
			double scale_zero = 0.1; // weight of zero samples
			for (int nsample = 0; nsample < num_samples; nsample++) { // if (sample_indx[nsample][1] != 0){
				int ntile = sample_indx[nsample][0]; // used tile index
				double st = St[ntile]; 
				if (!Double.isNaN(st)) {
					int mode =  sample_indx[nsample][1];
					double sqrt2 = Y[nsample]/K[nsample];
					sqrt2 *= sqrt2;
					double gi2 = (sqrt2 - N02) / (st * st);
					if (gi2 > 0.0) {
						wgi[mode] += weights[nsample];
						gi[mode] +=  weights[nsample] * Math.sqrt(gi2);
					} else {
						wgi[mode] += weights[nsample] * scale_zero;
					}
					
				}
				
			}
//			for (int mode = 1; mode < num_modes; mode++) {
			for (int mode = 0; mode < num_modes; mode++) {
				gi[mode] /= wgi[mode];
			}
			if (debug_level > -1) {
				System.out.println("Iteration "+nrefine+":");
				for (int i = 0; i < gi.length; i++) {
					System.out.println("g["+i+"] = "+gi[i]);
				}
			}
			gi[0] = 1.0; // all, inter - gain= 1.0
			// Recalculate St, using current gi
			St = new double [num_tiles];
			double [] stw = new double [num_tiles];
			for (int nsample = 0; nsample < num_samples; nsample++) { //  if (sample_indx[nsample][1] != 0){
				int ntile = sample_indx[nsample][0]; // used tile index
				int mode =  sample_indx[nsample][1];
				if (mode < 4) { // trying - only use inter data to adjust St
					double w =  weights[nsample];
					double sqrt2 = Y[nsample]/K[nsample];
					sqrt2 *= sqrt2;
					double st2 = (sqrt2 - N02) / (gi[mode] * gi[mode]);
					if (st2 > 0.0) {
						stw[ntile] += w;
						St[ntile] +=  w * Math.sqrt(st2);	
					}
				}
			}
			for (int ntile = 0; ntile < num_tiles; ntile++) {
				if (stw[ntile] > 0.0) {
					St[ntile] /= stw[ntile]; 
				} else {
					St[ntile] = Double.NaN;
				}
			}
			if (dbg_img != null) {
				double [] dbg_st = new double [noise_thresh.length];
				Arrays.fill(dbg_st,Double.NaN);
				for (int ntile = 0; ntile < tile_index.length; ntile++) {
					dbg_st[tile_index[ntile]] = St[ntile];
				}
				dbg_img[nrefine+1] = dbg_st;
			}
			
		}
		if (dbg_img  != null) {
			(new ShowDoubleFloatArrays()).showArrays(
					dbg_img,
					tilesX,
					noise_thresh.length/ tilesX,
					true,
					"pre_lma_st");
		}
		return;
	}
	private double [] getVector(
			boolean adjust_N0,
			boolean adjust_Gi,
			boolean adjust_St,
			double n0,
			double [] gi,
			double [] st) {
		this.adjust_N0 = adjust_N0;
		this.adjust_Gi = adjust_Gi;
		this.adjust_St = adjust_St;
		double [] v = new double [(adjust_N0 ? 1 : 0) + (adjust_Gi ? (gi.length - 1) : 0) + (adjust_St ? st.length : 0)];
		if (adjust_N0) {
			v[0] = n0;
		}
		if (adjust_Gi) {
			System.arraycopy(gi, 1, v, (adjust_N0 ? 1 : 0), gi.length-1);
		}
		if (adjust_St) {
			System.arraycopy(st, 0, v, (adjust_N0 ? 1 : 0) + (adjust_Gi ? (gi.length - 1) : 0), st.length);
		}
		return v;
	}
	
	private double getN0(double [] v) {
		return adjust_N0 ? v[0] : N0;
	}
	
	private double [] getGi(double [] v) {
		double [] gi = new double [this.gi.length];
		
		gi[0] = 1.0;
		if (adjust_Gi) {
			System.arraycopy(v, (adjust_N0 ? 1 : 0), gi, 1, gi.length-1);
		} else {
			System.arraycopy(this.gi, 1, gi, 1, gi.length-1);
		}
		return gi;
	}
	
	private double [] getSt(double [] v) {
		double [] st; // = new double [this.St.length]; // .length - this.gi.length];
		if (adjust_St) {
			st = new double [this.St.length]; // .length - this.gi.length];
			System.arraycopy(v, (adjust_N0 ? 1 : 0) + (adjust_Gi ? (gi.length - 1) : 0), st, 0, st.length);
			return st;
		} else {
			return St.clone();
		}
	}
	private double [] getYminusFxWeighted(double [] fx)
	{
		double [] y_minus_fx_weighted = new double [fx.length];
		for (int i = 0; i < fx.length; i++) {
			y_minus_fx_weighted[i] = weights[i]*(Y[i] - fx[i]);
			if (Double.isNaN(y_minus_fx_weighted[i])) {
				System.out.println("y_minus_fx_weighted["+i+"]= NaN");
				y_minus_fx_weighted[i] = 0.0;
			}
		}
		return y_minus_fx_weighted;
	}

	
	
	private double [] getFxJt(
			double []   vector,
			double [][] jt) { // should be either [vector.length][samples.size()] or null - then only fx is calculated
		double n0= getN0(vector);
		double [] gi = getGi(vector);
		double [] st = getSt(vector);
		double [] fx = new double [sample_indx.length];
		if (jt != null) { // calculate Jacobian
			for (int i = 0; i < jt.length; i++) {
				Arrays.fill(jt[i],0.0);
			}
		}

		for (int i = 0; i < fx.length; i++) {
			int itile = sample_indx[i][0];
			int mode =  sample_indx[i][1];
			double nv2 = st[itile]*gi[mode];
			nv2 *= nv2;
			nv2 -= n0*n0;
			if (nv2 > 0) {  // if <=0 - keep 0.0
				double sqrt = Math.sqrt(nv2);
				fx[i] = K[i] * sqrt;
				if (jt != null) {
					double Amti = K[i]/sqrt;
					int indx = 0;
					if (adjust_N0) {
						jt[indx++][i] = - Amti * n0;
					}
					double asg = Amti*st[itile]*gi[mode];
					if (adjust_Gi && (mode > 0)) {
						jt[indx + mode -1][i] = asg *st[itile];
						indx += gi.length -1;
					}
					if (adjust_St) {
						jt[indx + itile][i] = asg * gi[mode];
					}
				}
			}
		}
		return fx;
	}
	public double [][] getYDbg() {
		double [][] dbg_Y = new double [gi.length][dbgTilesX*dbgTilesY];
		for (int mode = 0; mode < dbg_Y.length; mode++) {
			Arrays.fill(dbg_Y[mode], Double.NaN);
		}
		for (int i = 0; i < Y.length; i++) {
			int itile = sample_indx[i][0];
			int mode =  sample_indx[i][1];
			int tile = tile_index[itile];
			dbg_Y[mode][tile] = Y[i];
		}
		return dbg_Y;
	}

	public double [][] getFxDbg() {
		double [][] dbg_Fx = new double [gi.length][dbgTilesX*dbgTilesY];
		for (int mode = 0; mode < dbg_Fx.length; mode++) {
			Arrays.fill(dbg_Fx[mode], Double.NaN);
		}
		double [] fx = getFxJt(
				vector, // double []   vector,
				null);  // double [][] jt)
		for (int i = 0; i < Y.length; i++) {
			int itile = sample_indx[i][0];
			int mode =  sample_indx[i][1];
			int tile =  tile_index[itile];
			dbg_Fx[mode][tile] = fx[i];
		}
		return dbg_Fx;
	}
	
	
	private double [] getWYMinusFx(
			double [] vector,
			double [] weights) {
		double [] fx = getFxJt(vector, null);
		double [] wymfx = new double [fx.length];
		for (int i = 0; i < wymfx.length; i++) {
			wymfx[i] = Y[i] - fx[i];
		}
		if (weights != null) {
			for (int i = 0; i < wymfx.length; i++) {
				wymfx[i] *= weights[i];
			}
		}
		return wymfx;
	}

	private double getRms(
			double [] fx) // weights should be normalized to sum==1.0
	{ 
		double swd2 = 0.0;
		for (int i = 0; i < fx.length; i++) {
			double d = Y[i] - fx[i];
			if (Double.isNaN(d)) {
				System.out.println("getRms(): Y["+i+"] - fx["+i+"] = NaN");
				d=0.0;
			}
			swd2 += weights[i]*d *d;
		}
		return Math.sqrt(swd2);
	}
	
	/*
	public double [][] getWJtJlambda_single(
			double      lambda,
			double [][] jt){
		int num_pars = jt.length;
		int nup_points = (num_pars > 0)? jt[0].length : 0;
		double [][] wjtjl = new double [num_pars][num_pars];
		for (int i = 0; i < num_pars; i++) {
			for (int j = i; j < num_pars; j++) {
				double d = 0.0;
				for (int k = 0; k < nup_points; k++) {
					d += this.weights[k]*jt[i][k]*jt[j][k];
				}
				wjtjl[i][j] = d;
				if (i == j) {
					wjtjl[i][j] += d * lambda;
				} else {
					wjtjl[j][i] = d;
				}
			}
		}
		return wjtjl;
	}
	*/
	
//public static int threadsMax = 100; // maximal number of threads to use
	public double [][] getWJtJlambda(
			final double      lambda,
			final double [][] jt){
		final int num_pars = jt.length;
		final int num_points = (num_pars > 0)? jt[0].length : 0;
		final Thread[] threads = ImageDtt.newThreadArray(threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		final double [][] wjtjl = new double [num_pars][num_pars];
		final int num_cells = num_pars * num_pars;
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				@Override
				public void run() {
					for (int iCell = ai.getAndIncrement(); iCell < num_cells; iCell = ai.getAndIncrement()) {
						int i = iCell / num_pars;
						int j = iCell % num_pars;
						if (j >= i) { // just skip unneeded
							double d = 0.0;
							for (int k = 0; k < num_points; k++) {
								d += weights[k]*jt[i][k]*jt[j][k];
							}
							wjtjl[i][j] = d;
							if (i == j) {
								wjtjl[i][j] += d * lambda;
								if (d == 0) {
									System.out.println("Diagonal ZERO for i=j="+i+" absolute tile = "+tile_index[i-8]); // assuming N0, gi[1]...gi[7]
									wjtjl[i][j] = 1.0; // Jt * (y-fx) will anyway be 0, so any value here should work.
								}
							} else {
								wjtjl[j][i] = d;
							}
						}
					}
				}
			};
		}
		ImageDtt.startAndJoin(threads);
		return wjtjl;
	}
	
	
	//last_rms = Double.NaN;
	public boolean runLma(
			boolean adjust_N0,
			boolean adjust_Gi,
			boolean adjust_St,
			double lambda,           // 0.1
			double lambda_scale_good,// 0.5
			double lambda_scale_bad, // 8.0
			double lambda_max,       // 100
			double rms_diff,         // 0.001
			int    num_iter,         // 20
			int    debug_level)
	{
		vector =  getVector(
				adjust_N0, // boolean adjust_N0,
				adjust_Gi, // boolean adjust_Gi,
				adjust_St, // boolean adjust_St,
				N0, // double n0,
				gi, // double [] gi,
				St); // double [] st) {
			
		boolean [] rslt = {false,false};
		this.last_rms = Double.NaN;
		int iter = 0;
		for (iter = 0; iter < num_iter; iter++) {
			rslt =  lmaStep(
					lambda,
					rms_diff,
					debug_level);
			if (rslt == null) {
				return false; // need to check
			}
			if (debug_level > 1) {
				System.out.println("LMA step "+iter+": {"+rslt[0]+","+rslt[1]+"} full RMS= "+good_or_bad_rms+
						" ("+initial_rms+"), lambda="+lambda);
			}
			if (rslt[1]) {
				break;
			}
			if (rslt[0]) { // good
				lambda *= lambda_scale_good;
			} else {
				lambda *= lambda_scale_bad;
				if (lambda > lambda_max) {
					break; // not used in lwir
				}
			}
		}
		if (rslt[0]) { // better, but num tries exceeded
			if (iter >= num_iter) {
				if (debug_level > 0) System.out.println("Step "+iter+": Improved, but number of steps exceeded maximal");
			} else {
				if (debug_level > 0) System.out.println("Step "+iter+": LMA: Success");
			}

		} else { // improved over initial ?
			if (last_rms < initial_rms) {
				rslt[0] = true;
				if (debug_level > 0) System.out.println("Step "+iter+": Failed to converge, but result improved over initial");
			} else {
				if (debug_level > 0) System.out.println("Step "+iter+": Failed to converge");
			}
		}
		if (debug_level > 0) {
			System.out.println("LMA: full RMS="+last_rms+" ("+initial_rms+"), lambda="+lambda);
		}
		if (rslt[0]) { // success
			if (adjust_N0) {
				N0 = getN0(vector);
			}
			if (adjust_Gi) {
				gi = getGi(vector);
			}
			if (adjust_St) {
				St = getSt(vector);
			}
		}
		
		return rslt[0];
	}
	
	
	
	public boolean [] lmaStep(
			double lambda,
			double rms_diff,
			int debug_level) {
		int num_points = this.weights.length; // includes 2 extra for regularization
		int num_pars = vector.length;
		boolean [] rslt = {false,false};
		boolean dbg_img = debug_level > 2;
		if (Double.isNaN(last_rms)) { //first time, need to calculate all (vector is valid)
			last_jt = new double [num_pars][num_points];
			double [] fx = getFxJt(
					vector, // double []   vector,
					last_jt); // double [][] jt) { // should be either [vector.length][samples.size()] or null - then only fx is calculated
			if (fx == null) {
				return null; // need to re-init/restart LMA
			}
			last_ymfx = getYminusFxWeighted (fx);
			last_rms =  getRms(fx);
			initial_rms = last_rms;
			good_or_bad_rms = last_rms;
			if (dbg_img) {
				double [][] dbg_Y =  getYDbg();
				(new ShowDoubleFloatArrays()).showArrays(
						dbg_Y,
						dbgTilesX,
						dbgTilesY,
						true,
						"dbg_Y");
				double [][] dbg_Fx = getFxDbg();
				(new ShowDoubleFloatArrays()).showArrays(
						dbg_Fx,
						dbgTilesX,
						dbgTilesY,
						true,
						"dbg_Fx");
			}
			/*
			if (debug_level > 3) {
				debugJt(
						0.000001, // double      delta, // 0.2, //
						this.vector); // double []   vector);
			}
			*/
		}
		Matrix y_minus_fx_weighted = new Matrix(last_ymfx, last_ymfx.length);

		Matrix wjtjlambda = new Matrix(getWJtJlambda(
				lambda, // *10, // temporary
				last_jt)); // double [][] jt)
		
		if (debug_level>2) {
			System.out.println("JtJ + lambda*diag(JtJ");
			wjtjlambda.print(18, 6);
		}
		Matrix jtjl_inv = null;
		try {
			jtjl_inv = wjtjlambda.inverse(); // check for errors
		} catch (RuntimeException e) {
			rslt[1] = true;
			if (debug_level > 0) {
				System.out.println("Singular Matrix!");
			}
			return rslt;
		}
		if (debug_level>2) {
			System.out.println("(JtJ + lambda*diag(JtJ)).inv()");
			jtjl_inv.print(18, 6);
		}
//last_jt has NaNs
		Matrix jty = (new Matrix(this.last_jt)).times(y_minus_fx_weighted);
		if (debug_level>2) {
			System.out.println("Jt * (y-fx)");
			jty.print(18, 6);
		}

		Matrix mdelta = jtjl_inv.times(jty);
		if (debug_level>2) {
			System.out.println("mdelta");
			mdelta.print(18, 6);
		}

		double []  delta = mdelta.getColumnPackedCopy();
		double [] new_vector = this.vector.clone();
		for (int i = 0; i < num_pars; i++) new_vector[i]+= delta[i];
		// being optimistic, modify jt and last_ymfx in place, restore if failed
		double [] fx = getFxJt( // 
				new_vector, // double []   vector,
				last_jt); // double [][] jt) { // should be either [vector.length][samples.size()] or null - then only fx is calculated
		if (fx == null) {
			return null; // need to re-init/restart LMA
		}
		double [] new_ymfx = getYminusFxWeighted (fx); // weighted
		double rms =  getRms(fx); // new rms
		good_or_bad_rms = rms; 
		if (rms < last_rms) { // improved
			rslt[0] = true;
			rslt[1] = rms >=(last_rms * (1.0 - rms_diff));
			last_rms = rms;               // update
			vector = new_vector.clone();  // update
			last_ymfx = new_ymfx;         // update
			if (debug_level > 1) {
				System.out.print("New vector: ");
				for (int np = 0; np < vector.length; np++) {
					System.out.print(this.vector[np]+" ");
				}
				System.out.println();
			}

		} else { // worsened
			rslt[0] = false;
			rslt[1] = false; // do not know, caller will decide
			// restore state
			/*
			fx = getFxJt( // 
					vector, // it was not updated // new_vector, // double []   vector,
					last_jt); // double [][] jt) { // should be either [vector.length][samples.size()] or null - then only fx is calculated
			if (fx == null) {
				return null; // need to re-init/restart LMA
			}
			last_ymfx = getYminusFxWeighted (fx);
			last_rms =  getRms(fx);
			*/
			/*
			if (debug_level > 2) {
				debugJt(
						0.000001, // double      delta,
						this.vector); // double []   vector);
			}
            */
		}
		return rslt;
	}
}
