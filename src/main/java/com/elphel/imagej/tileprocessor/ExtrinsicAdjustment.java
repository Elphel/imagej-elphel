package com.elphel.imagej.tileprocessor;

import com.elphel.imagej.common.ShowDoubleFloatArrays;

import Jama.Matrix;

/**
 **
 ** ExtrinsicAdjustment - Adjust cameras extrinsic parameters and ERS correction
 **
 ** Copyright (C) 2019 Elphel, Inc.
 **
 ** -----------------------------------------------------------------------------**
 **
 **  ExtrinsicAdjustment.java is free software: you can redistribute it and/or modify
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




public class ExtrinsicAdjustment {
	static final int NUM_SENSORS =   4;

	static final int INDX_DISP =     0; // composite
	static final int INDX_STRENGTH = 1;
	static final int INDX_X0 =       2;
	static final int INDX_TARGET =  10; // target disparity
	static final int INDX_DIFF =    11; // differential disparity (now w/o magic composite =target + diff)
	static final int INDX_DYDDISP0 =12; // derivative of pixel y over disparity (for ERS)
	static final int INDX_PX =      16;
	static final int INDX_DD0 =     18;
	static final int INDX_ND0 =     22;
	static final int INDX_LENGTH =  26;
	static final int POINTS_SAMPLE = 2 * NUM_SENSORS +1; // points per sample residual disparity, 4*dd, *nd
	static final String [] DATA_TITLES = {
			"Disparity", "Strength",
			"DX-0","DY-0","DX-1","DY-1","DX-2","DY-2","DX-3","DY-3",
			"Target Disparity","Diff. Disparity",
			"dY_dD-0","dY_dD-1","dY_dD-2","dY_dD-3",
			"pX","pY",
			"DD-0", "DD-1","DD-2","DD-3","ND-0", "ND-1","ND-2","ND-3"};

	// next values are only updated after success
	private double []         last_rms =        null; // {rms, rms_pure}, matching this.vector
	private double []         good_or_bad_rms = null; // just for diagnostics, to read last (failed) rms
	private double []         initial_rms =     null; // {rms, rms_pure}, first-calcualted rms
	private double []         last_ymfx =       null;
	private double [][]       last_jt =         null;
	private double []         weights; // normalized so sum is 1.0 for all - samples and extra regularization terms
	private boolean [] 		  force_disparity = null;                     // boolean [] force_disparity, // same dimension as dsdn, true if disparity should be controlled
	private double            pure_weight; // weight of samples only
//	private double []         values;
	private GeometryCorrection.CorrVector corr_vector = null;
	private boolean []        par_mask =        null;
	private boolean           use_rig_offsets = false;
	private double [][]       measured_dsxy =   null;
	private double [][]       dy_ddisparity =   null; // conveniently extracted from  dsdn
	private double [][]       x0y0 =            null; //
	private double [][]       world_xyz =       null;
	private double []         weight_window =   null; // center area is more reliable

	public GeometryCorrection geometryCorrection = null;
	public int clusterSize;
	public int clustersX;
	public int clustersY;

	public double [] getOldNewRMS() {
		double [] on_rms = new double[2];
		if (initial_rms != null) {
			on_rms[0] = initial_rms[0];
		} else {
			on_rms[0] = Double.NaN;
		}
		if (last_rms != null) {
			on_rms[1] = last_rms[1];
		} else {
			on_rms[1] = Double.NaN;
		}
		return on_rms;
	}

	public ExtrinsicAdjustment (
			GeometryCorrection gc,
			int         clusterSize,
	  		int         clustersX,
	  		int         clustersY) {
		geometryCorrection = gc;
		this.clusterSize = clusterSize;
		this.clustersX = clustersX;
		this.clustersY = clustersY;
	}

	public void showInput(double[][] data, String title) {
		int clusters = clustersX * clustersY;
		double [][] pixels = new double [ExtrinsicAdjustment.INDX_LENGTH][clusters];
		for (int cluster = 0; cluster < clusters; cluster++) {
			if (data[cluster] != null) {
				for (int c = 0; c < data[cluster].length; c++) {
					pixels[c][cluster] = data[cluster][c];
				}
			} else {
				for (int c = 0; c < pixels.length; c++) {
					pixels[c][cluster] = (c == ExtrinsicAdjustment.INDX_STRENGTH)? 0.0: Double.NaN;
				}
			}
		}
		 (new ShowDoubleFloatArrays()).showArrays(
				 pixels,
				 clustersX,
				 clustersY,
				 true,
				 title,
				 ExtrinsicAdjustment.DATA_TITLES);
	}

	private void showX0Y0(double [][] xy0, String title) {
		String [] titles = {"xnd-0","ynd-0","xnd-1","ynd-1","xnd-2","ynd-2","xnd-3","ynd-3"};
		int clusters = clustersX * clustersY;
		double [][] pixels = new double [titles.length][clusters];
		for (int cluster = 0; cluster < clusters; cluster++) {
			if (xy0[cluster] != null) {
				for (int c = 0; c < xy0[cluster].length; c++) {
					pixels[c][cluster] = xy0[cluster][c];
				}
			} else {
				for (int c = 0; c < pixels.length; c++) {
					pixels[c][cluster] = Double.NaN;
				}
			}
		}
		 (new ShowDoubleFloatArrays()).showArrays(
				 pixels,
				 clustersX,
				 clustersY,
				 true,
				 title,
				 titles);
	}

	public GeometryCorrection.CorrVector  solveCorr (
			double      marg_fract,        // part of half-width, and half-height to reduce weights
			boolean     use_disparity,     // adjust disparity-related extrinsics
			boolean     use_aztilts,       // Adjust azimuths and tilts excluding disparity
			boolean     use_diff_rolls,    // Adjust differential rolls (3 of 4 angles)
//			boolean     force_convergence, // if true try to adjust convergence (disparity, symmetrical parameter 0) even with no disparity
//			                               // data, using just radial distortions
			int         min_num_forced,    // minimal number of clusters with forced disparity to use it
	  		boolean     common_roll,       // Enable common roll (valid for high disparity range only)
			boolean     corr_focalLength,  // Correct scales (focal length temperature? variations)
			boolean     ers_rot,           // Enable ERS correction of the camera rotation
			boolean     ers_forw,          // Enable ERS correction of the camera linear movement in z direction
			boolean     ers_side,          // Enable ERS correction of the camera linear movement in x direction
			boolean     ers_vert,          // Enable ERS correction of the camera linear movement in y direction
			// add balancing-related here?
	  		int         manual_par_sel,    // Manually select the parameter mask bit 0 - sym0, bit1 - sym1, ... (0 - use boolean flags, != 0 - ignore boolean flags)
			double      weight_disparity,
			double      weight_lazyeye,
	  		double [][] measured_dsxy_in,     //
	  		boolean [] force_disparity_in,    // boolean [] force_disparity,
//			GeometryCorrection geometryCorrection,
	  		boolean     use_main,          // corr_rots_aux != null;
			GeometryCorrection.CorrVector corr_vector_meas,
			double []   old_new_rms,       // should be double[2]
			int         debugLevel)
	{

		this.corr_vector = corr_vector_meas.clone();
		this.use_rig_offsets = false;
		this.measured_dsxy = measured_dsxy_in;
		this.force_disparity = force_disparity_in;
		boolean dbg_images = debugLevel > 1; // 2; // -3; // 2; // 1;

		weight_window = getWeightWindow(marg_fract);

		if (dbg_images) {
			 (new ShowDoubleFloatArrays()).showArrays(
					 weight_window,
					 clustersX,
					 clustersY,
					 "weight_window");

			showInput(
					measured_dsxy, // double[][] data,
					"input data");// String title);
		}
		world_xyz =  getWorldXYZ();

		x0y0 = getXYNondistorted(
				corr_vector,
				true); // boolean set_dydisp)

		if (dbg_images) {
			showX0Y0(
					x0y0, // double[][] data,
					"nondistorted X0Y0");// String title);
		}

		this.par_mask = geometryCorrection.getParMask(
				use_disparity, // has_disparity, // boolean use_disparity,
				use_aztilts,       // Adjust azimuths and tilts excluding disparity
				use_diff_rolls,    // Adjust differential rolls (3 of 4 angles)
				common_roll,// boolean common_roll,
				corr_focalLength,  // boolean corr_focalLength);
				ers_rot,           // boolean ers_rot,           // Enable ERS correction of the camera rotation
				ers_forw,          // Enable ERS correction of the camera linear movement in z direction
				ers_side,          // Enable ERS correction of the camera linear movement in x direction
				ers_vert,          // Enable ERS correction of the camera linear movement in y direction
		  		manual_par_sel);   // Manually select the parameter mask bit 0 - sym0, bit1 - sym1, ... (0 - use boolean flags, != 0 - ignore boolean flags)

		this.weights = getWeights(
				measured_dsxy,     // double  [][] measured_dsxy,
				force_disparity,   // boolean [] force_disparity, // same dimension as dsdn, true if disparity should be controlled
				min_num_forced,    //				int min_num_forced,
				weight_disparity,  // double weight_disparity,
				weight_lazyeye);   // double weight_lazyeye);

		double lambda = 0.1;
		double lambda_scale_good = 0.5;
		double lambda_scale_bad =  8.0;
		double lambda_max =      100;
		double rms_diff =          0.001;
		int    num_iter =          20;

		boolean lma_OK = runLma(
				lambda,            // double lambda,           // 0.1
				lambda_scale_good, // double lambda_scale_good,// 0.5
				lambda_scale_bad,  // double lambda_scale_bad, // 8.0
				lambda_max,        // double lambda_max,       // 100
				rms_diff,          // double rms_diff,         // 0.001
				num_iter,          // int    num_iter,         // 20
				debugLevel);       // int    debug_level)
		if (old_new_rms != null) {
			double [] on_rms = getOldNewRMS();
			old_new_rms[0] = on_rms[0];
			old_new_rms[1] = on_rms[1];
		}
		return lma_OK? corr_vector : null;
	}

	private double [][] getXYNondistorted(
			GeometryCorrection.CorrVector corr_vector,
			boolean set_dydisp){
		int clusters =clustersX * clustersY;
		Matrix []   corr_rots =  corr_vector.getRotMatrices(); // get array of per-sensor rotation matrices
		Matrix [][] deriv_rots = corr_vector.getRotDeriveMatrices();
		double []   imu =        corr_vector.getIMU(); // i)

		double [][] xyND =   new double[clusters][];
		if (set_dydisp) {
			dy_ddisparity = new double[clusters][];
		}
		for (int cluster = 0; cluster < clusters; cluster++) {
			if (measured_dsxy[cluster] != null) {
				if (set_dydisp) {
					dy_ddisparity[cluster] = new double[NUM_SENSORS];
					for (int i = 0; i < NUM_SENSORS; i++) {
						dy_ddisparity[cluster][i] = measured_dsxy[cluster][INDX_DYDDISP0 + i];
					}
				}
				xyND[cluster] = geometryCorrection.getPortsNonDistortedCoordinatesAndDerivatives( // USED in lwir
						geometryCorrection,                  // GeometryCorrection gc_main,
						use_rig_offsets,                     // boolean     use_rig_offsets,
						corr_rots,                           // Matrix []   rots,
						deriv_rots,                          // Matrix [][] deriv_rots,
						null,                                // double [][] pXYNDderiv, // if not null, should be double[8][]
						dy_ddisparity[cluster],              // dy_ddisparity,   // double [][] disp_dist, //disp_dist[i][2] or null
						imu,                                 // double []   imu,
						world_xyz[cluster],                  // double []   xyz, // world XYZ for ERS correction
						measured_dsxy[cluster][INDX_PX + 0], // double px,
						measured_dsxy[cluster][INDX_PX + 1], // double py,
						measured_dsxy[cluster][INDX_TARGET]); // double disparity);
			}
		}
		return xyND;
	}

	private double [][] getWorldXYZ(){
		int clusters =clustersX * clustersY;
		double [][] world_xyz =   new double[clusters][];
		for (int cluster = 0; cluster < clusters; cluster++) {
			if (measured_dsxy[cluster] != null) {
				double disparity = measured_dsxy[cluster][INDX_TARGET];
				if (disparity > 0.0) {
				world_xyz[cluster] = geometryCorrection.getWorldCoordinates( // USED in lwir
						measured_dsxy[cluster][INDX_PX + 0], // double px,
						measured_dsxy[cluster][INDX_PX + 1], // double py,
						disparity,                           // double disparity,
						true);                               // boolean correctDistortions)
				}
			}
		}
		return world_xyz;
	}




/*
	private double [] getYminusFx(
			GeometryCorrection.CorrVector corr_vector)
	{
		int clusters = clustersX * clustersY;
		Matrix []   corr_rots =  corr_vector.getRotMatrices(); // get array of per-sensor rotation matrices
		Matrix [][] deriv_rots = corr_vector.getRotDeriveMatrices();
		double [] y_minus_fx = new double  [clusters * POINTS_SAMPLE];
		double [] imu = corr_vector.getIMU(); // i)
		for (int cluster = 0; cluster < clusters;  cluster++) if (measured_dsxy[cluster] != null){
			double [] ddnd = geometryCorrection.getPortsDDNDAndDerivatives( // USED in lwir
					geometryCorrection,     // GeometryCorrection gc_main,
					use_rig_offsets,        // boolean     use_rig_offsets,
					corr_rots,              // Matrix []   rots,
					deriv_rots,             // Matrix [][] deriv_rots,
					null,                   // double [][] DDNDderiv,     // if not null, should be double[8][]
					dy_ddisparity[cluster], // double []   dy_ddisparity,   // double [][] disp_dist, //disp_dist[i][2] or null
					imu,                    // double []   imu,
					x0y0[cluster],          // double []   pXYND0,        // per-port non-distorted coordinates corresponding to the correlation measurements
					world_xyz[cluster],     // double []   xyz, // world XYZ for ERS correction
					measured_dsxy[cluster][ExtrinsicAdjustment.INDX_PX + 0],  // double      px,
					measured_dsxy[cluster][ExtrinsicAdjustment.INDX_PX + 1],  // double      py,
					measured_dsxy[cluster][ExtrinsicAdjustment.INDX_TARGET]); // double      disparity);
			//arraycopy(Object src, int srcPos, Object dest, int destPos, int length)
			ddnd[0] = -ddnd[0];
			if ((force_disparity != null) && force_disparity[cluster]) {
				ddnd[0] -= measured_dsxy[cluster][ExtrinsicAdjustment.INDX_DIFF];
			}
///			ddnd[0] = -measured_dsxy[cluster][ExtrinsicAdjustment.INDX_DIFF] - ddnd[0];
			for (int i = 0; i < NUM_SENSORS; i++) {
				ddnd[i + 1] = -measured_dsxy[cluster][ExtrinsicAdjustment.INDX_DD0 + i] - ddnd[i + 1];
				ddnd[i + 5] = -measured_dsxy[cluster][ExtrinsicAdjustment.INDX_ND0 + i] - ddnd[i + 5];
			}
			System.arraycopy(ddnd, 0, y_minus_fx, cluster * POINTS_SAMPLE, POINTS_SAMPLE);
		}
		return y_minus_fx;
	}
*/
	private double [] getWYmFxRms( // USED in lwir
			double []   fx) {
		int clusters = clustersX * clustersY;
		double rms = 0, rms_pure = 0;
		for (int cluster = 0; cluster < clusters;  cluster++) if (measured_dsxy[cluster] != null){
			int indx0 = POINTS_SAMPLE * cluster + 0;
			// force_disparity - compensate measured_dsxy[cluster][ExtrinsicAdjustment.INDX_DIFF],
			// false - keep (so force fx==0
			double d = - fx[indx0];
				if ((force_disparity != null) && force_disparity[cluster]) {
					d -= measured_dsxy[cluster][ExtrinsicAdjustment.INDX_DIFF];
				}
			fx[indx0] = this.weights[indx0] * d;
			rms += fx[indx0]*d; // sum of weights
			for (int cam = 0; cam < NUM_SENSORS; cam++) {
				int indx = indx0 + cam + 1;
				d = (-measured_dsxy[cluster][ExtrinsicAdjustment.INDX_DD0 + cam] - fx[indx]);
				fx[indx] = this.weights[indx] * d;
				rms += fx[indx] * d; // sum of weights
				indx = indx0 + cam + 5; // nd
				d = (-measured_dsxy[cluster][ExtrinsicAdjustment.INDX_ND0 + cam] - fx[indx]);
				fx[indx] = this.weights[indx] * d;
				rms += fx[indx] * d; // sum of weights
			}
		}
		rms_pure = Math.sqrt(rms)/this.pure_weight;
		// Calculate other regularization terms here if needed
		rms = Math.sqrt(rms);
		double [] rslt = {rms, rms_pure};
		return rslt;
	}

	private double [] getWeightWindow(double marg_fraction) { // 0.0 - no margins, 1.0 - pure cosine
		double mf_hor =  marg_fraction;
		double mf_vert = marg_fraction;
		double [] wx = new double [clustersX];
		double [] wy = new double [clustersY];
		double [] w = new double [clustersX * clustersY];
		int [] boost_wnd = {33,15,40,35};
		double boost_scale = 1.0; // 100.0;

		double center = 0.5 * (clustersX - 1);
		double marg = center * mf_hor;
		for (int x = 0; x <= clustersX / 2; x++) {
			if (x < marg) {
				wx[x] = Math.sin(Math.PI * x / 2.0 / marg);
				wx[x] *= wx[x];
			} else {
				wx[x] = 1.0;
			}
			wx[clustersX - 1 -x ] = wx[x];
		}
		center = 0.5 * (clustersY - 1);
		marg = center * mf_vert;
		for (int y = 0; y <= clustersY / 2; y++) {
			if (y < marg) {
				wy[y] = Math.sin(Math.PI * y / 2.0 / marg);
				wy[y] *= wx[y];
			} else {
				wy[y] = 1.0;
			}
			wy[clustersY - 1 - y ] = wy[y];
		}
		for (int y = 0; y < clustersY; y++) {
			for (int x = 0; x < clustersX; x++) {
				w[y * clustersX + x] = wx[x]*wy[y];
				if (boost_scale > 1.0) {
					if ((x >= boost_wnd[0]) && (x < boost_wnd[2]) && (y >= boost_wnd[1]) && (y < boost_wnd[3])) {
						w[y * clustersX + x] *= boost_scale;
					}
				}
			}
		}
		return w;

	}

	private double [] getWeights(
			double  [][] measured_dsxy,
			boolean [] force_disparity, // same dimension as dsdn, true if disparity should be controlled
			int min_num_forced,  // if number of forced samples exceeds this, zero out weights of non-forced

			double weight_disparity,
			double weight_lazyeye)
	{
		int clusters = clustersX * clustersY;
		double [] weights = new double  [clusters * POINTS_SAMPLE];
		double sw = 0.0;
		int num_forced = 0;
		if (force_disparity != null) for (int cluster = 0; cluster < clusters;  cluster++) if (force_disparity[cluster])num_forced ++;
		boolean use_forced = num_forced >= min_num_forced;
		for (int cluster = 0; cluster < clusters;  cluster++) if (measured_dsxy[cluster] != null){
			double s = measured_dsxy[cluster][ExtrinsicAdjustment.INDX_STRENGTH] * weight_window[cluster];
			double w = s * weight_disparity;
			if (use_forced && !force_disparity[cluster]) {
				w = 0.0;
			}
			weights[cluster * POINTS_SAMPLE + 0] =  w;
			sw += w;

			w = s * weight_lazyeye;
			for (int i = 1; i < POINTS_SAMPLE; i++) {
				weights[cluster * POINTS_SAMPLE + i] = w;
			}
			sw += w * (POINTS_SAMPLE - 1);
		}
		if (sw <= 0.0) {
			return null;
		}
		double k = 1.0/sw;
		for (int cluster = 0; cluster < clusters;  cluster++) if (measured_dsxy[cluster] != null){
			weights[cluster * POINTS_SAMPLE + 0] *= k;
			double w = weights[cluster * POINTS_SAMPLE + 1] * k;
			for (int i = 1; i < POINTS_SAMPLE; i++) {
				weights[cluster * POINTS_SAMPLE + i] = w;
			}
		}
		this.pure_weight = 1.0;
		return weights;
	}

	private double [] getFx(
			GeometryCorrection.CorrVector corr_vector)
	{
		int clusters = clustersX * clustersY;
		Matrix []   corr_rots =  corr_vector.getRotMatrices(); // get array of per-sensor rotation matrices
		Matrix [][] deriv_rots = corr_vector.getRotDeriveMatrices();
		double [] imu = corr_vector.getIMU(); // i)
		double [] y_minus_fx = new double  [clusters * POINTS_SAMPLE];
		for (int cluster = 0; cluster < clusters;  cluster++) if (measured_dsxy[cluster] != null){
			double [] ddnd = geometryCorrection.getPortsDDNDAndDerivatives( // USED in lwir
					geometryCorrection,     // GeometryCorrection gc_main,
					use_rig_offsets,        // boolean     use_rig_offsets,
					corr_rots,              // Matrix []   rots,
					deriv_rots,             // Matrix [][] deriv_rots,
					null,                   // double [][] DDNDderiv,     // if not null, should be double[8][]
					dy_ddisparity[cluster], // double []   dy_ddisparity,   // double [][] disp_dist, //disp_dist[i][2] or null
					imu,                    // double []   imu,
					x0y0[cluster],          // double []   pXYND0,        // per-port non-distorted coordinates corresponding to the correlation measurements
					world_xyz[cluster],     // double []   xyz, // world XYZ for ERS correction
					measured_dsxy[cluster][ExtrinsicAdjustment.INDX_PX + 0],  // double      px,
					measured_dsxy[cluster][ExtrinsicAdjustment.INDX_PX + 1],  // double      py,
					measured_dsxy[cluster][ExtrinsicAdjustment.INDX_TARGET]); // double      disparity);
			//arraycopy(Object src, int srcPos, Object dest, int destPos, int length)
            //		    System.arraycopy(src_pixels, 0, dst_pixels, 0, src_pixels.length); /* for the borders closer to 1/2 kernel size*/
			ddnd[0] = ddnd[0];
			for (int i = 0; i < NUM_SENSORS; i++) {
				ddnd[i + 1] = ddnd[i + 1];
				ddnd[i + 5] = ddnd[i + 5];
			}
			System.arraycopy(ddnd, 0, y_minus_fx, cluster*POINTS_SAMPLE, POINTS_SAMPLE);
		}
		return y_minus_fx;
	}

	private double [][] getJacobianTransposed(
			GeometryCorrection.CorrVector corr_vector,
			double delta){
		int clusters = clustersX * clustersY;
		int num_pars = getNumPars();
		double [][] jt = new double  [num_pars][clusters * POINTS_SAMPLE ];
		double rdelta = 1.0/delta;
		for (int par = 0; par < num_pars; par++) {
			double [] pars = new double[num_pars];
			pars[par] =  delta;
			GeometryCorrection.CorrVector corr_delta = geometryCorrection.getCorrVector(pars, par_mask);
			GeometryCorrection.CorrVector corr_vectorp = corr_vector.clone();
			GeometryCorrection.CorrVector corr_vectorm = corr_vector.clone();
			corr_vectorp.incrementVector(corr_delta,  0.5);
			corr_vectorm.incrementVector(corr_delta, -0.5);
			double [] fx_p = getFx(corr_vectorp);
			double [] fx_m = getFx(corr_vectorm);
			for (int i = 0; i < fx_p.length; i++) {
				jt[par][i] = (fx_p[i] - fx_m[i])*rdelta;
			}
		}
		return jt;
	}

	private String [] getSymNames() {
		int num_pars = getNumPars();
		String [] names = new String[num_pars];
		int indx = 0;
		for (int i = 0; i < par_mask.length; i++) if (par_mask[i]){
			names[indx++] = "S"+i;
		}
		return names;
	}

	private double dbgJacobians(
			GeometryCorrection.CorrVector corr_vector,
			double delta,
			boolean graphic) {
		int gap = 10; //pix
		int clusters = clustersX * clustersY;
		int num_pars = getNumPars();
		String [] titles = getSymNames();
		double [][] jt =        getJacobianTransposed(corr_vector);
		double [][] jt_delta =  getJacobianTransposed(corr_vector, delta);
		double tot_error = 0;
		double [][] err =   new double [num_pars][POINTS_SAMPLE];
		double [] err_par = new double [num_pars];
		for (int par = 0; par < num_pars; par++) {
			for (int cluster = 0; cluster < clusters;  cluster++) if (measured_dsxy[cluster] != null){
				for (int i = 0; i < POINTS_SAMPLE; i++) {
					int indx = cluster *  POINTS_SAMPLE+ i;
					if (Math.abs(jt[par][indx] - jt_delta[par][indx]) > err[par][i]) {
						err[par][i] = Math.abs(jt[par][indx] - jt_delta[par][indx]);
					}
				}
			}
			for (int i = 0; i < POINTS_SAMPLE; i++) {
				if (err[par][i] > err_par[par]) {
					err_par[par] = err[par][i];
				}
			}
			if (err_par[par] > tot_error) {
				tot_error = err_par[par];
			}
		}
		System.out.println("Maximal error for all parameters = "+tot_error);
		System.out.println(String.format("%4s %10.6s %10s %10s %10s %10s %10s %10s %10s %10s %10s",
				"","max","disp", "dd0", "dd1", "dd2", "dd3", "nd0", "nd1", "nd2", "nd3"));

		for (int par = 0; par < num_pars; par++) {
			System.out.println(String.format("%4s %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f",
					titles[par],err_par[par],err[par][0],
					err[par][1],err[par][2],err[par][3],err[par][4],err[par][5],err[par][6],err[par][7],err[par][8]));

		}
		if (graphic) {
			String [] titles3 = new String[num_pars * 3];

			int width  = 3 * clustersX + 2 * gap;
			int height = 3 * clustersY + 2 * gap;
			double [][] dbg_img = new double [num_pars * 3][width*height];
			for (int par = 0; par < num_pars; par++) {
				titles3[3 * par + 0] = titles[par]+"";
				titles3[3 * par + 1] = titles[par]+"_delta";
				titles3[3 * par + 2] = titles[par]+"_diff";
				for (int mode = 0; mode < POINTS_SAMPLE; mode++) {
					int x0 = (mode % 3) * (clustersX + gap);
					int y0 = (mode / 3) * (clustersY + gap);
					for (int cluster = 0; cluster < clusters;  cluster++) {
						int x = x0 + (cluster % clustersX);
						int y = y0 + (cluster / clustersX);
						int pix = x + y * width;
						int indx = cluster * POINTS_SAMPLE + mode;

						if (measured_dsxy[cluster] != null){
							dbg_img[3 * par + 0][pix] =  jt[par][indx];
							dbg_img[3 * par + 1][pix] =  jt_delta[par][indx];
							dbg_img[3 * par + 2][pix] =  jt[par][indx] - jt_delta[par][indx];
						} else {
							dbg_img[3 * par + 0][pix] =  Double.NaN;
							dbg_img[3 * par + 1][pix] =  Double.NaN;
							dbg_img[3 * par + 2][pix] =  Double.NaN;
						}
					}
				}
			}
			 (new ShowDoubleFloatArrays()).showArrays(
					 dbg_img,
					 width,
					 height,
					 true,
					 "Debug_Jacobians",
					 titles3);
		}
		return tot_error;
	}

	private void dbgYminusFx(
			double []   fx,
			String title) {
		int gap = 10; //pix
		int clusters = clustersX * clustersY;
		String [] titles = {"Y", "-fX", "Y+fx"};
		int width  = 3 * clustersX + 2 * gap;
		int height = 3 * clustersY + 2 * gap;
		double [][] dbg_img = new double [3][width*height];
		for (int mode = 0; mode < POINTS_SAMPLE; mode++) {
			int x0 = (mode % 3) * (clustersX + gap);
			int y0 = (mode / 3) * (clustersY + gap);
			for (int cluster = 0; cluster < clusters;  cluster++) {
				int x = x0 + (cluster % clustersX);
				int y = y0 + (cluster / clustersX);
				int pix = x + y * width;
				int indx = cluster * POINTS_SAMPLE + mode;
				if (measured_dsxy[cluster] != null){
					if (mode ==0) {
						dbg_img[0][pix] =  measured_dsxy[cluster][ExtrinsicAdjustment.INDX_DIFF];
						dbg_img[1][pix] = -fx[indx];
						dbg_img[2][pix] =  measured_dsxy[cluster][ExtrinsicAdjustment.INDX_DIFF] + fx[indx];
					} else {
						dbg_img[0][pix] =  measured_dsxy[cluster][ExtrinsicAdjustment.INDX_DD0 + mode - 1];
						dbg_img[1][pix] = -fx[indx];
						dbg_img[2][pix] =  measured_dsxy[cluster][ExtrinsicAdjustment.INDX_DD0 + mode - 1] + fx[indx];
					}
				} else {
					dbg_img[ 0][pix] =  Double.NaN;
					dbg_img[ 1][pix] =  Double.NaN;
					dbg_img[ 2][pix] =  Double.NaN;
				}
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

	private void dbgXY(
			GeometryCorrection.CorrVector corr_vector,
			String title) {
		int gap = 10; //pix
		int clusters = clustersX * clustersY;
		String [] titles = {"meas", "correction", "diff"};
		int width  = 3 * clustersX + 2 * gap;
		int height = 3 * clustersY + 2 * gap;
		double [][] xy = getXYNondistorted(corr_vector, false);

		double [][] dbg_img = new double [3][width*height];
		for (int mode = 0; mode < 2 * NUM_SENSORS; mode++) {
			int x0 = (mode % 3) * (clustersX + gap);
			int y0 = (mode / 3) * (clustersY + gap);
			for (int cluster = 0; cluster < clusters;  cluster++) {
				int x = x0 + (cluster % clustersX);
				int y = y0 + (cluster / clustersX);
				int pix = x + y * width;
//				int indx = cluster * POINTS_SAMPLE + mode;
				if (measured_dsxy[cluster] != null){
						dbg_img[0][pix] =  measured_dsxy[cluster][ExtrinsicAdjustment.INDX_X0+mode];
						dbg_img[1][pix] =  xy[cluster][mode] - x0y0[cluster][mode];
						dbg_img[2][pix] =  xy[cluster][mode] - x0y0[cluster][mode] - measured_dsxy[cluster][ExtrinsicAdjustment.INDX_X0+mode];
				} else {
					dbg_img[ 0][pix] =  Double.NaN;
					dbg_img[ 1][pix] =  Double.NaN;
					dbg_img[ 2][pix] =  Double.NaN;
				}
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





	private double [][] getJacobianTransposed(
			GeometryCorrection.CorrVector corr_vector)
	{
		int clusters = clustersX * clustersY;
		int num_pars = getNumPars();
		double [][] jt = new double  [num_pars][clusters * POINTS_SAMPLE ];
		Matrix []   corr_rots =  corr_vector.getRotMatrices(); // get array of per-sensor rotation matrices
		Matrix [][] deriv_rots = corr_vector.getRotDeriveMatrices();
		double [] imu = corr_vector.getIMU(); // i)
		for (int cluster = 0; cluster < clusters;  cluster++) if (measured_dsxy[cluster] != null){
//			Mismatch mm = mismatch_list.get(indx);
//			double [] pXY = mm.getPXY();
			// will calculate 9 rows (disparity, dd0, dd1,cdd2, dd3, nd0, nd1, nd2, nd3}, columns - parameters
			double [][] deriv = geometryCorrection.getPortsDDNDDerivatives( // USED in lwir
					geometryCorrection,     // GeometryCorrection gc_main,
					use_rig_offsets,        // boolean     use_rig_offsets,
					corr_rots,              // Matrix []   rots,
					deriv_rots,             // Matrix [][] deriv_rots,
					dy_ddisparity[cluster], // double []   dy_ddisparity,   // double [][] disp_dist, //disp_dist[i][2] or null
					imu,                    // double []   imu,
					world_xyz[cluster],                  // double []   xyz, // world XYZ for ERS correction
					measured_dsxy[cluster][ExtrinsicAdjustment.INDX_PX + 0], // double      px,
					measured_dsxy[cluster][ExtrinsicAdjustment.INDX_PX + 1], // double      py,
					measured_dsxy[cluster][ExtrinsicAdjustment.INDX_TARGET]); // double      disparity);
///			int dbg_index = cluster; // dbg_index (pXY, dbg_decimate);
			// convert to symmetrical coordinates
			// derivatives of each port coordinates (in pixels) for each of selected symmetric all parameters (sym0 is convergence for disparity)
			 double [][] jt_partial = corr_vector.getJtPartial(
						deriv, // double [][] port_coord_deriv,
						par_mask); // boolean [] par_mask

			// put partial (for 1 cluster 9  (disparity, dd0, dd1,cdd2, dd3, nd0, nd1, nd2, nd3}, transposed jacobian into full (all tiles) transposed Jacobian
			for (int npar = 0; npar < jt.length; npar++){
				for (int n = 0; n < jt_partial[npar].length; n++){
						jt[npar][jt_partial[npar].length * cluster + n] = jt_partial[npar][n];
						if (Double.isNaN(jt_partial[npar][n])) {
							System.out.println("getJacobianTransposed(): npar="+npar+", cluster="+cluster+", n="+n);
						}
				}
			}
		}
		return jt;
	}

	private double [][] getWJtJlambda( // USED in lwir
			double      lambda,
			double [][] jt){
		int num_pars = jt.length;
		int nup_points = jt[0].length;
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


	private boolean runLma(
			double lambda,           // 0.1
			double lambda_scale_good,// 0.5
			double lambda_scale_bad, // 8.0
			double lambda_max,       // 100
			double rms_diff,         // 0.001
			int    num_iter,         // 20
			int    debug_level)
	{
		boolean [] rslt = {false,false};
		this.last_rms = null;
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
				System.out.println("LMA step "+iter+": {"+rslt[0]+","+rslt[1]+"} full RMS= "+good_or_bad_rms[0]+
						" ("+initial_rms[0]+"), pure RMS="+good_or_bad_rms[1]+" ("+initial_rms[1]+") + lambda="+lambda);
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
		if (rslt[0]) { // better
			if (iter >= num_iter) { // better, but num tries exceeded
				if (debug_level > 0) System.out.println("Step "+iter+": Improved, but number of steps exceeded maximal");
			} else {
				if (debug_level > 0) System.out.println("Step "+iter+": LMA: Success");
			}

		} else { // improved over initial ?
			if (last_rms[0] < initial_rms[0]) {
				rslt[0] = true;
				if (debug_level > 0) System.out.println("Step "+iter+": Failed to converge, but result improved over initial");
			} else {
				if (debug_level > 0) System.out.println("Step "+iter+": Failed to converge");
			}
		}
		if (debug_level > 0) {
			System.out.println("LMA: full RMS="+last_rms[0]+" ("+initial_rms[0]+"), pure RMS="+last_rms[1]+" ("+initial_rms[1]+") + lambda="+lambda);
		}

		return rslt[0];
	}

	private int getNumPars() {
		int num_pars = 0;
		for (int i = 0; i < par_mask.length; i++) if (par_mask[i]) num_pars ++;
		return num_pars;
	}

	// returns {success, done}
	private boolean [] lmaStep(
			double lambda,
			double rms_diff,
			int debug_level) {
//		int num_points = this.weights.length; // includes 2 extra for regularization
//		int num_pars = getNumPars();
		boolean [] rslt = {false,false};
		if (this.last_rms == null) { //first time, need to calculate all (vector is valid)
//			this.last_ymfx = getFxJt(
//					this.vector, // double []   vector,
//					this.last_jt); // double [][] jt) { // should be either [vector.length][samples.size()] or null - then only fx is calculated
			this.last_jt =  getJacobianTransposed(corr_vector); // new double [num_pars][num_points];
			this.last_ymfx = getFx(corr_vector);

			if (debug_level > 2) {
				dbgYminusFx(this.last_ymfx, "Initial y-fX");
			}


			if (last_ymfx == null) {
				return null; // need to re-init/restart LMA
			}
			this.last_rms = getWYmFxRms(
					this.last_ymfx); // modifies this.last_ymfx (weights and subtracts fx
			this.initial_rms = this.last_rms.clone();
			this.good_or_bad_rms = this.last_rms.clone();
			// TODO: Restore/implement
			if (debug_level > 3) {
				 dbgJacobians(
							corr_vector, // GeometryCorrection.CorrVector corr_vector,
							1E-5, // double delta,
							true); //boolean graphic)
			}
		}
		Matrix y_minus_fx_weighted = new Matrix(this.last_ymfx, this.last_ymfx.length);

		Matrix wjtjlambda = new Matrix(getWJtJlambda(
				lambda, // *10, // temporary
				this.last_jt)); // double [][] jt)
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
			System.out.println("(JtJ + lambda*diag(JtJ).inv()");
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
		GeometryCorrection.CorrVector corr_delta = geometryCorrection.getCorrVector(delta, par_mask);

///		double [] new_vector = this.vector.clone();
		GeometryCorrection.CorrVector new_vector = this.corr_vector.clone();
		double scale = 1.0;
//		boolean ok =
		new_vector.incrementVector(corr_delta, scale); // ok = false if there are nay NaN-s

///		for (int i = 0; i < num_pars; i++) new_vector[i]+= delta[i];
		// being optimistic, modify jt and last_ymfx in place, restore if failed


///		this.last_ymfx = getFxJt(
///				new_vector, // double []   vector,
///				this.last_jt); // double [][] jt) { // should be either [vector.length][samples.size()] or null - then only fx is calculated
		this.last_jt =  getJacobianTransposed(new_vector); // new double [num_pars][num_points];
		this.last_ymfx = getFx(new_vector);
		if (debug_level > 2) {
			dbgYminusFx(this.last_ymfx, "next y-fX");
			dbgXY(new_vector, "XY-correction");
		}

		if (last_ymfx == null) {
			return null; // need to re-init/restart LMA
		}

		double [] rms = getWYmFxRms(this.last_ymfx); // modifies this.last_ymfx
		this.good_or_bad_rms = rms.clone();
		if (rms[0] < this.last_rms[0]) { // improved
			rslt[0] = true;
			rslt[1] = rms[0] >=(this.last_rms[0] * (1.0 - rms_diff));
			this.last_rms = rms.clone();

			this.corr_vector = new_vector.clone();
			if (debug_level > 2) {
				System.out.print("New vector: "+new_vector.toString());
///				for (int np = 0; np < vector.length; np++) {
///					System.out.print(this.vector[np]+" ");
///				}
				System.out.println();
			}
		} else { // worsened
			rslt[0] = false;
			rslt[1] = false; // do not know, caller will decide
			// restore state
///			this.last_ymfx = getFxJt( // recalculate fx
///					this.vector, // double []   vector,
///					this.last_jt); // double [][] jt) { // should be either [vector.length][samples.size()] or null - then only fx is calculated
			this.last_jt =  getJacobianTransposed(corr_vector); // new double [num_pars][num_points];
			this.last_ymfx = getFx(corr_vector);

			if (last_ymfx == null) {
				return null; // need to re-init/restart LMA
			}
			this.last_rms = getWYmFxRms(this.last_ymfx); // modifies this.last_ymfx
			if (debug_level > 2) {
				 dbgJacobians(
							corr_vector, // GeometryCorrection.CorrVector corr_vector,
							1E-5, // double delta,
							true); //boolean graphic)
			}
		}
		return rslt;
	}
}
