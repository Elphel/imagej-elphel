package com.elphel.imagej.tileprocessor;

import com.elphel.imagej.common.DoubleGaussianBlur;
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
	/*
	static final int NUM_SENSORS =   4;
	static final int INDX_DISP =     0; // composite
	static final int INDX_STRENGTH = 1;
	static final int INDX_X0 =       2; // 2 * NUM_SENSORS long
	static final int INDX_TARGET =  10; // target disparity
	static final int INDX_DIFF =    11; // differential disparity (now w/o magic composite =target + diff)
	static final int INDX_DYDDISP0 =12; // derivative of pixel y over disparity (for ERS) NUM_SENSORS long
	static final int INDX_PX =      16; // 2 ?
	static final int INDX_DD0 =     18; // NUM_SENSORS long
	static final int INDX_ND0 =     22; // NUM_SENSORS long
	static final int INDX_PYDIST =  26; // average pY to calculate ERS time difference
	static final int INDX_LENGTH =  26 + 4; // total length
	static final int POINTS_SAMPLE = 2 * NUM_SENSORS +1; // points per sample residual disparity, 4*dd, *nd
	static final String [] DATA_TITLES = {
			"Disparity", "Strength",
			"DX-0","DY-0","DX-1","DY-1","DX-2","DY-2","DX-3","DY-3",
			"Target Disparity","Diff. Disparity",
			"dY_dD-0","dY_dD-1","dY_dD-2","dY_dD-3",
			"pX","pY",
			"DD-0", "DD-1","DD-2","DD-3",
			"ND-0", "ND-1","ND-2","ND-3",
			"pY0",  "pY1", "pY2", "pY3"};
    */
	// changing values to consolidate fixed values (not dependent on the number of sensors to use static)
//    static final int NUM_SENSORS =   4;
	static final int INDX_DISP =     0; // composite
	static final int INDX_STRENGTH = 1;
	static final int INDX_TARGET =   2; // target disparity
	static final int INDX_DIFF =     3; // differential disparity (now w/o magic composite =target + diff)
	static final int INDX_PX =       4; // 4,5
	static final int INDX_X0 =       6; // 2 * NUM_SENSORS long
	// keep static for above, replace non-static 
//	static final int INDX_DYDDISP0 =14; // derivative of pixel y over disparity (for ERS) NUM_SENSORS long
//	static final int INDX_DD0 =     18; // NUM_SENSORS long
//	static final int INDX_ND0 =     22; // NUM_SENSORS long
//	static final int INDX_PYDIST =  26; // average pY to calculate ERS time difference
//	static final int INDX_LENGTH =  26 + 4; // total length
//	static final int POINTS_SAMPLE = 2 * NUM_SENSORS +1; // points per sample residual disparity, 4*dd, *nd
	
	final int indx_dyddisp0; //INDX_X0 + 2 * num_sensors  =14; // derivative of pixel y over disparity (for ERS) NUM_SENSORS long
	final int indx_dd0; //  =     18; // NUM_SENSORS long
	final int indx_nd0; //  =     22; // NUM_SENSORS long
	final int indx_pydist; //  =  26; // average pY to calculate ERS time difference
	final int indx_length; //  =  26 + 4; // total length
	final int points_sample; //  = 2 * NUM_SENSORS +1; // points per sample residual disparity, 4*dd, *nd
	final String [] data_titles; 
	/*
	static final String [] DATA_TITLES = {
			"Disparity", "Strength",
			"DX-0","DY-0","DX-1","DY-1","DX-2","DY-2","DX-3","DY-3",
			"Target Disparity","Diff. Disparity",
			"dY_dD-0","dY_dD-1","dY_dD-2","dY_dD-3",
			"pX","pY",
			"DD-0", "DD-1","DD-2","DD-3",
			"ND-0", "ND-1","ND-2","ND-3",
			"pY0",  "pY1", "pY2", "pY3"};
	*/
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
	private CorrVector corr_vector = null;
	private boolean []        par_mask =        null;
	private boolean           use_rig_offsets = false;
	private double [][]       measured_dsxy =   null; // only set in solveCorr()
//	private double [][]       dy_ddisparity =   null; // conveniently extracted from  dsdn
	private double [][]       pY_offset =       null; // conveniently extracted from  dsdn - per-sensor average pY to calculate ERS offset (set in getXYNondistorted())
	private double [][]       x0y0 =            null; //
	private double[][]        world_xyz =       null; // only set in solveCorr()
	private double []         weight_window =   null; // center area is more reliable

	public final GeometryCorrection geometryCorrection;
	public int clusterSize;
	public int clustersX;
	public int clustersY;
	
	public final int num_sensors;
	
	public double dbg_delta = 0; // 1.0E-5; // if not 0 - use delta instead of the derivatives in getJacobianTransposed

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

	/*
	final int indx_dyddisp0; //INDX_X0 + 2 * num_sensors  =14; // derivative of pixel y over disparity (for ERS) NUM_SENSORS long
	final int indx_dd0; //  =     18; // NUM_SENSORS long
	final int indx_nd0; //  =     22; // NUM_SENSORS long
	final int indx_pydist; //  =  26; // average pY to calculate ERS time difference
	final int indx_length; //  =  26 + 4; // total length
	final int points_sample; //  = 2 * NUM_SENSORS +1; // points per sample residual disparity, 4*dd, *nd
//	static final int INDX_DYDDISP0 =14; // derivative of pixel y over disparity (for ERS) NUM_SENSORS long
//	static final int INDX_DD0 =     18; // NUM_SENSORS long
//	static final int INDX_ND0 =     22; // NUM_SENSORS long
//	static final int INDX_PYDIST =  26; // average pY to calculate ERS time difference
//	static final int INDX_LENGTH =  26 + 4; // total length
	static final int POINTS_SAMPLE = 2 * NUM_SENSORS +1; // points per sample residual disparity, 4*dd, *nd
	
	 */
	public static int get_INDX_DYDDISP0(int num_sensors) {return INDX_X0 + 2 * num_sensors;}
	public static int get_INDX_DD0     (int num_sensors) {return INDX_X0 + 3 * num_sensors;}
	public static int get_INDX_ND0     (int num_sensors) {return INDX_X0 + 4 * num_sensors;}
	public static int get_INDX_PYDIST  (int num_sensors) {return INDX_X0 + 5 * num_sensors;}
	public static int get_INDX_LENGTH  (int num_sensors) {return INDX_X0 + 6 * num_sensors;}
	public static int get_POINTS_SAMPLE(int num_sensors) {return 2 * num_sensors + 1;}
	
	public ExtrinsicAdjustment (
			GeometryCorrection gc,
			int         clusterSize,
	  		int         clustersX,
	  		int         clustersY) {
		this.geometryCorrection = gc;
		this.num_sensors = gc.getNumSensors();
		this.clusterSize = clusterSize;
		this.clustersX =   clustersX;
		this.clustersY =   clustersY;
		
		indx_dyddisp0 = get_INDX_DYDDISP0(num_sensors); // = 14; // derivative of pixel y over disparity (for ERS) NUM_SENSORS long
		indx_dd0 =      get_INDX_DD0     (num_sensors); // = 18; // NUM_SENSORS long
		indx_nd0 =      get_INDX_ND0     (num_sensors); // = 22; // NUM_SENSORS long
		indx_pydist =   get_INDX_PYDIST  (num_sensors); // = 26; // average pY to calculate ERS time difference
		indx_length =   get_INDX_LENGTH  (num_sensors); // = 26 + 4; // total length
		points_sample = get_POINTS_SAMPLE(num_sensors); // = 2 * NUM_SENSORS +1; //points per sample residual disparity, 4*dd, *nd
		data_titles = new String[indx_length];
		String [] titles_static = {"Disparity", "Strength",
				"Target Disparity","Diff. Disparity",
				"pX","pY"};
		int indx = 0;
		for (int i = 0; i < titles_static.length; i++) {
			data_titles[indx++] = titles_static[i];
		}
		for (int i = 0; i < num_sensors; i++) {
			data_titles[indx++] = "DX-"+i;
			data_titles[indx++] = "DY-"+i;
		}
		for (int i = 0; i < num_sensors; i++) {
			data_titles[indx++] = "dY_dD-"+i;
		}
		for (int i = 0; i < num_sensors; i++) {
			data_titles[indx++] = "DD-"+i;
		}
		for (int i = 0; i < num_sensors; i++) {
			data_titles[indx++] = "ND-"+i;
		}
		for (int i = 0; i < num_sensors; i++) {
			data_titles[indx++] = "pY"+i;
		}
	}

	public void showInput(double[][] data, String title) {
		int clusters = clustersX * clustersY;
//		double [][] pixels = new double [ExtrinsicAdjustment.INDX_LENGTH+4][clusters];
//		String [] titles = new String[ExtrinsicAdjustment.INDX_LENGTH+4];
//		for (int i = 0; i < ExtrinsicAdjustment.INDX_LENGTH; i++) {
		double [][] pixels = new double [indx_length + 4][clusters];
		String [] titles =  new String[indx_length + 4];
			for (int i = 0; i < indx_length; i++) {
			//			titles[i] = ExtrinsicAdjustment.DATA_TITLES[i];
			titles[i] = data_titles[i];
		}
//		titles[ExtrinsicAdjustment.INDX_LENGTH+0]="Force_disparity";
//		titles[ExtrinsicAdjustment.INDX_LENGTH+1]="dx-sum";
//		titles[ExtrinsicAdjustment.INDX_LENGTH+2]="dy_sum";
//		titles[ExtrinsicAdjustment.INDX_LENGTH+3]="dd_sum";
		titles[indx_length+0]="Force_disparity";
		titles[indx_length+1]="dx-sum";
		titles[indx_length+2]="dy_sum";
		titles[indx_length+3]="dd_sum";
		
		for (int cluster = 0; cluster < clusters; cluster++) {
			if (data[cluster] != null) {
				for (int c = 0; c < data[cluster].length; c++) {
					pixels[c][cluster] = data[cluster][c];
				}
				for (int i = 0;i <4; i++) {
//					pixels[ExtrinsicAdjustment.INDX_LENGTH+1][cluster] += 0.25 * data[cluster][ExtrinsicAdjustment.INDX_X0 + 2 * i];
//					pixels[ExtrinsicAdjustment.INDX_LENGTH+2][cluster] += 0.25 * data[cluster][ExtrinsicAdjustment.INDX_X0 + 2 * i + 1];
//					pixels[ExtrinsicAdjustment.INDX_LENGTH+3][cluster] += 0.25 * data[cluster][ExtrinsicAdjustment.INDX_DD0 + i];
					pixels[indx_length+1][cluster] += 0.25 * data[cluster][INDX_X0 + 2 * i];
					pixels[indx_length+2][cluster] += 0.25 * data[cluster][INDX_X0 + 2 * i + 1];
					pixels[indx_length+3][cluster] += 0.25 * data[cluster][indx_dd0 + i];
				}
			} else {
				for (int c = 0; c < pixels.length; c++) {
					//					pixels[c][cluster] = (c == ExtrinsicAdjustment.INDX_STRENGTH)? 0.0: Double.NaN;
					pixels[c][cluster] = Double.NaN;
				}
			}
			if (force_disparity!=null) {
//				pixels[ExtrinsicAdjustment.INDX_LENGTH][cluster] = force_disparity[cluster]?1.0:0.0;
				pixels[indx_length][cluster] = force_disparity[cluster]?1.0:0.0;
			}
		}
		 (new ShowDoubleFloatArrays()).showArrays(
				 pixels,
				 clustersX,
				 clustersY,
				 true,
				 title,
				 titles); //ExtrinsicAdjustment.DATA_TITLES);
	}
	/*
	public static void showLYInput(
			double[][] data,
			String     title,
			int        clustersX,
			int        clustersY) {
		int clusters = clustersX * clustersY;
		double [][] pixels = new double [ExtrinsicAdjustment.INDX_LENGTH+3][clusters];
		String [] titles = new String[ExtrinsicAdjustment.INDX_LENGTH+3];
		for (int i = 0; i < ExtrinsicAdjustment.INDX_LENGTH; i++) {
			titles[i] = ExtrinsicAdjustment.DATA_TITLES[i];
//			titles[i] = data_titles[i];

		}
		titles[ExtrinsicAdjustment.INDX_LENGTH+0]="dx-sum";
		titles[ExtrinsicAdjustment.INDX_LENGTH+1]="dy_sum";
		titles[ExtrinsicAdjustment.INDX_LENGTH+2]="dd_sum";
		
		for (int cluster = 0; cluster < clusters; cluster++) {
			if (data[cluster] != null) {
				for (int c = 0; c < data[cluster].length; c++) {
					pixels[c][cluster] = data[cluster][c];
				}
				for (int i = 0;i <4; i++) {
					pixels[ExtrinsicAdjustment.INDX_LENGTH+0][cluster] += 0.25 * data[cluster][ExtrinsicAdjustment.INDX_X0 + 2 * i];
					pixels[ExtrinsicAdjustment.INDX_LENGTH+1][cluster] += 0.25 * data[cluster][ExtrinsicAdjustment.INDX_X0 + 2 * i + 1];
					pixels[ExtrinsicAdjustment.INDX_LENGTH+2][cluster] += 0.25 * data[cluster][ExtrinsicAdjustment.INDX_DD0 + i];
				}
			} else {
				for (int c = 0; c < pixels.length; c++) {
					//					pixels[c][cluster] = (c == ExtrinsicAdjustment.INDX_STRENGTH)? 0.0: Double.NaN;
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
				 titles); //ExtrinsicAdjustment.DATA_TITLES);
	}
	*/
	
	private void showX0Y0(double [][] xy0, String title) {
		int num_sensors = this.geometryCorrection.getNumSensors();
//		String [] titles = {"xnd-0","ynd-0","xnd-1","ynd-1","xnd-2","ynd-2","xnd-3","ynd-3"};
		String [] titles = new String [ 2 * num_sensors];
		for (int i = 0; i < num_sensors; i++) {
			titles[2*i + 0] = "xnd-"+i;
			titles[2*i + 1] = "ynd-"+i;
		}
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

	public CorrVector  solveCorr (
			double      marg_fract,        // part of half-width, and half-height to reduce weights
			boolean     use_disparity,     // adjust disparity-related extrinsics
			double      inf_min_disparity, // minimal disparity for infinity (from average, possibly tilted) 
			double      inf_max_disparity, // maximal disparity for infinity (from average, possibly tilted) 
			double      inf_min_disp_abs,  // minimal disparity for infinity (absolute) 
			double      inf_max_disp_abs,  // maximal disparity for infinity (absolute) 
			boolean     en_infinity_tilt,  // select infinity tiles form right/left tilted (false - from average)  
			boolean     infinity_right_left, // balance weights between right and left halves of infinity
			boolean     use_aztilts,       // Adjust azimuths and tilts excluding disparity
			boolean     use_diff_rolls,    // Adjust differential rolls (3 of 4 angles)
			int         min_num_forced,    // minimal number of clusters with forced disparity to use it
	  		boolean     common_roll,       // Enable common roll (valid for high disparity range only)
			boolean     corr_focalLength,  // Correct scales (focal length temperature? variations)
			boolean     ers_rot,           // Enable ERS correction of the camera rotation
			boolean     ers_forw,          // Enable ERS correction of the camera linear movement in z direction
			boolean     ers_side,          // Enable ERS correction of the camera linear movement in x direction
			boolean     ers_vert,          // Enable ERS correction of the camera linear movement in y direction
			// add balancing-related here?
	  		int         manual_par_sel,      // Manually select the parameter mask bit 0 - sym0, bit1 - sym1, ... (0 - use boolean flags, != 0 - ignore boolean flags)
			double      weight_infinity,     // 0.3, total weight of infinity tiles fraction (0.0 - 1.0) 
			double      weight_disparity,    // 0.0 disparity weight relative to the sum of 8 lazy eye values of the same tile 
			double      weight_disparity_inf,// 0.5 disparity weight relative to the sum of 8 lazy eye values of the same tile for infinity 
			double      max_disparity_far,   // 5.0 reduce weights of near tiles proportional to sqrt(max_disparity_far/disparity) 
			double      max_disparity_use,   // 5.0 (default 1000)disable near objects completely - use to avoid ERS
			double      min_dfe, // = 1.75;
			double      max_dfe, // = 5.0; // <=0 - disable feature
			// moving objects filtering
			boolean     moving_en,         // enable filtering areas with potentially moving objects 
			boolean     moving_apply,      // apply filtering areas with potentially moving objects 
			double      moving_sigma,      // blurring sigma for moving objects = 1.0;
			double      max_mov_disparity, // disparity limit for moving objects detection = 75.0;
			double      rad_to_hdiag_mov,  // radius to half-diagonal ratio to remove high-distortion corners = 0.7 ; // 0.8
			double      max_mov_average,   // do not attempt to detect moving objects if ERS is not accurate for terrain = .25;
			double      mov_min_L2,        // threshold for moving objects = 0.75;
			
	  		double [][] measured_dsxy_in,     //
	  		boolean [] force_disparity_in,    // boolean [] force_disparity,
	  		boolean     use_main,          // corr_rots_aux != null;
			CorrVector corr_vector_meas,
			double []   old_new_rms,       // should be double[2]
			int         debugLevel)
	{

		this.corr_vector = corr_vector_meas.clone(); // current correction vector (before adjustment)
		this.use_rig_offsets = false;
		this.measured_dsxy = measured_dsxy_in;
		// FIXME - had to invert to work (still need to check for use_disparity = false)
		for (int i = 0; i < this.measured_dsxy.length; i++) if (this.measured_dsxy[i] != null) {
			this.measured_dsxy[i][INDX_DIFF] = -this.measured_dsxy[i][INDX_DIFF]; // invert disparity sign 
		}
		
		this.force_disparity = force_disparity_in;
		boolean dbg_images = debugLevel > 0; // 2; // -3; // 2; // 1; manually set here
		double rad_to_hdiag = 0.8;
		weight_window = getWeightWindow(
				marg_fract,
				rad_to_hdiag); // limit corners 

		if (dbg_images  && (debugLevel > 0)) { // disabled
			 (new ShowDoubleFloatArrays()).showArrays(
					 weight_window,
					 clustersX,
					 clustersY,
					 "weight_window");

			showInput(
					measured_dsxy, // double[][] data,
					"input data");// String title);
		}
		world_xyz =  getWorldXYZ(); // freeze world coordinates for measured pX,pY and disparity
		// calculate x,y non-distorted offsets for current correction vectors (to subtract from the new (modified) values
		// 0/0 in the center (in the optical center), in pixels 
		x0y0 = getXYNondistorted( // Sets pY_offset[], needed for .getPortsDDNDDerivativesNew() (in getJacobianTransposed)
				corr_vector,
				true); // boolean set_dydisp)

		if (dbg_images && (debugLevel > 0)) { // disabled
			showX0Y0(
					x0y0, // double[][] data,
					"nondistorted X0Y0");// String title);
		}
/*
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
*/		  		
		boolean [] filtered_infinity = null;
		
		double [] dfe = null;
		if (max_dfe > 0) {
			dfe = distanceFromEdge ( // all 0? too little, non-continuous
					force_disparity,
					measured_dsxy,
					min_dfe, // 1.75
					max_dfe, // 5.0
					dbg_images  && (debugLevel > 0));  // disabled
		}
//		boolean en_tilt =            true;
//		double min_infinity_abs =    -1.0;
//		double max_infinity_abs =     0.5;
		if (use_disparity) {
			filtered_infinity = filterInfinity(
					measured_dsxy,      // double  [][] measured_dsxy,
					force_disparity,    // boolean [] force_disparity, // same dimension as dsdn, true if disparity should be controlled
					dfe,                // double [] distance_from_edge,
					en_infinity_tilt,   // boolean en_tilt,  // consider right/left infinity tilt (will be disabled if more than *abs difference over width)
					min_num_forced/4,   // int min_in_half,
					inf_min_disp_abs,   // double min_infinity_abs,
					inf_max_disp_abs,   // double max_infinity_abs,
					inf_min_disparity,  // double min_infinity,
					inf_max_disparity); // double max_infinity
		 }
/*
		this.weights = getWeights( // will ignore window for infinity (already used for selection)
				 measured_dsxy,        // double  [][] measured_dsxy,
				 (use_disparity? force_disparity: null),   // boolean [] force_disparity, // same dimension as dsdn, true if disparity should be controlled
				 filtered_infinity,    // boolean [] filtered_infinity,
				 dfe, // 	double  []  distance_from_edge,// to reduce weight of the mountain ridge, increase clouds (or null)
				 min_num_forced,       //				int min_num_forced,
				 infinity_right_left, //	boolean infinity_right_left,  // each halve should have > min_num_forced, will calculate separate average
				 weight_infinity,      // double weight_infinity,       // total weight of infinity tiles fraction (0.0 - 1.0) 
				 weight_disparity,     // double weight_disparity,      // disparity weight relative to the sum of 8 lazy eye values of the same tile 
				 weight_disparity_inf, // double weight_disparity_inf,  // disparity weight relative to the sum of 8 lazy eye values of the same tile for infinity 
				 max_disparity_far,   // 					double max_disparity_far)     // reduce weights of near tiles proportional to sqrt(max_disparity_far/disparity)
				 max_disparity_use);
*/		
		int [] inf_stat = setWeights( // number right, number left
				 measured_dsxy,        // double  [][] measured_dsxy,
				 (use_disparity? force_disparity: null),   // boolean [] force_disparity, // same dimension as dsdn, true if disparity should be controlled
				 filtered_infinity,    // boolean [] filtered_infinity,
				 dfe, // 	double  []  distance_from_edge,// to reduce weight of the mountain ridge, increase clouds (or null)
				 min_num_forced,       //				int min_num_forced,
				 infinity_right_left, // boolean infinity_right_left,  // each halve should have > min_num_forced, will calculate separate average
				 weight_infinity,      // double weight_infinity,       // total weight of infinity tiles fraction (0.0 - 1.0) 
				 weight_disparity,     // double weight_disparity,      // disparity weight relative to the sum of 8 lazy eye values of the same tile 
				 weight_disparity_inf, // double weight_disparity_inf,  // disparity weight relative to the sum of 8 lazy eye values of the same tile for infinity 
				 max_disparity_far,   // 					double max_disparity_far)     // reduce weights of near tiles proportional to sqrt(max_disparity_far/disparity)
				 max_disparity_use);
		if (use_disparity) {
			if (inf_stat[0] < min_num_forced) { // 4
				System.out.println("Too few infinity tiles ("+inf_stat[0]+"<"+(min_num_forced)+") to adjust disparity");
				// disable parameters...all extrinsic
				use_disparity =    false;
				use_aztilts =      false;
				use_diff_rolls =   false;
				common_roll =      false;
				corr_focalLength = false;
			} else if ((inf_stat[1] < min_num_forced/4) || (inf_stat[1] < 1)) { 
				System.out.println("Too few infinity tiles ("+inf_stat[1]+"<"+(min_num_forced/2)+") in one half to balance right/left");
				// disable parameters: all extrinsic but disparity (S0)
				use_aztilts =      false;
				use_diff_rolls =   false;
				common_roll =      false;
				corr_focalLength = false;
			}
		}
		// verify some parameters still remain?

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
		 if (this.par_mask == null) {
				System.out.println("All parameters disabled, nothing to adjust");
				return null;
		 }
		
		 if (moving_en) { // debugLevel > -10) { // temporary
			 String title_moving =       (debugLevel > -3) ? "Moving_objects": null;
			 String title_moving_extra = (debugLevel > -2) ? "Moving_objects_filtering" : null;
			 boolean debug_text =        (debugLevel > -3);
			 
			 boolean [] moving_maybe = selectMovingObjects(
					 moving_sigma,              // double sigma,             //  = 1.0;
					 max_mov_disparity,  // double max_mov_disparity, //  = 75.0;
					 rad_to_hdiag_mov,   // double rad_to_hdiag_mov,  //  = 0.7 ; // 0.8
					 max_mov_average,    // double max_mov_average,   //  = .25;
					 mov_min_L2,             // double min_l2,            //  = 0.75;
					 title_moving,       // String title,             // "Moving_objects"
					 title_moving_extra, // String title_extra,       // "Moving_objects_filtering"
					 debug_text          // boolean debug_text
					);
			 
			 if (moving_apply) {
				 blockSelectedClusters( // will do nothing for null
						 measured_dsxy, // double  [][] measured_dsxy,
						 this.weights, // double []    weights, // will be updated
						 moving_maybe); // 	boolean []   prohibit
			 }
		 }		

		 if (dbg_images) { // temporary
			 dbgYminusFxWeight(
					 this.last_ymfx,
					 this.weights,
					 "Initial_y-fX_after_moving_objects");
		 }
		 
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
			CorrVector corr_vector,
			boolean set_dydisp){
		int clusters =clustersX * clustersY;
		Matrix []   corr_rots =  corr_vector.getRotMatrices(); // get array of per-sensor rotation matrices
		Matrix [][] deriv_rots = corr_vector.getRotDeriveMatrices();
		double []   imu =        corr_vector.getIMU(); // i)

		double [][] xyND =   new double[clusters][];
		double [][] pXY0 = geometryCorrection.getPXY0();
		int [] top_woi = geometryCorrection.getWOITops();
/*		
		double [] py_time0 = new double [NUM_SENSORS];
		for (int i = 0; i < NUM_SENSORS; i++) {
			py_time0[i] = pXY0[i][1] - top_woi[i];
		}
*/
		if (set_dydisp) {
///			dy_ddisparity = new double[clusters][];
			pY_offset =     new double[clusters][];
			
		}
		for (int cluster = 0; cluster < clusters; cluster++) {
			if ((cluster == 1735 ) || (cluster==1736)){
				System.out.print("");
			}
			if (measured_dsxy[cluster] != null) {
				if (set_dydisp) {
					pY_offset [cluster] = new double[num_sensors];
					double lines_avg = 0;
					for (int i = 0; i < num_sensors; i++) {
						// find average time
						pY_offset[cluster][i] =  measured_dsxy[cluster][indx_pydist + i] - top_woi[i]; // time (scanlines) since frame start
						lines_avg += pY_offset[cluster][i]; 
					}
					lines_avg /= num_sensors;
					for (int i = 0; i < num_sensors; i++) {
						pY_offset[cluster][i] -=  lines_avg;
					}
					
				}
				xyND[cluster] = geometryCorrection.getPortsNonDistortedCoordinatesAndDerivativesNew( // USED in lwir
						geometryCorrection,                  // GeometryCorrection gc_main,
						use_rig_offsets,                     // boolean     use_rig_offsets,
						corr_rots,                           // Matrix []   rots,
						deriv_rots,                          // Matrix [][] deriv_rots,
						null,                                // double [][] pXYNDderiv, // if not null, should be double[8][]
						pY_offset[cluster],                  // double []   py_offset,  // array of per-port average pY offset from the center (to correct ERS) or null (for no ERS)
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


	private double [] getWYmFxRms( // USED in lwir
			double []   fx) {
		int clusters = clustersX * clustersY;
		double rms = 0, rms_pure = 0;
		for (int cluster = 0; cluster < clusters;  cluster++) if (measured_dsxy[cluster] != null){
			int indx0 = points_sample * cluster + 0;
			// force_disparity - compensate measured_dsxy[cluster][ExtrinsicAdjustment.INDX_DIFF],
			// false - keep (so force fx==0
			double d = - fx[indx0];
				if ((force_disparity != null) && force_disparity[cluster]) {
					d -= measured_dsxy[cluster][ExtrinsicAdjustment.INDX_DIFF];
				}
			fx[indx0] = this.weights[indx0] * d;
			rms += fx[indx0]*d; // sum of weights
			for (int cam = 0; cam < num_sensors; cam++) {
				int indx = indx0 + cam + 1;
//				d = (-measured_dsxy[cluster][ExtrinsicAdjustment.INDX_DD0 + cam] - fx[indx]);
				d = (-measured_dsxy[cluster][indx_dd0 + cam] - fx[indx]);
				fx[indx] = this.weights[indx] * d;
				rms += fx[indx] * d; // sum of weights
				indx = indx0 + cam + 5; // nd
//				d = (-measured_dsxy[cluster][ExtrinsicAdjustment.INDX_ND0 + cam] - fx[indx]);
				d = (-measured_dsxy[cluster][indx_nd0 + cam] - fx[indx]);
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
	/*
	@Deprecated
	private double [] getWeightWindowOld(double marg_fraction) { // 0.0 - no margins, 1.0 - pure cosine
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
	*/

	private double [] getWeightWindow(
			double marg_fraction, // 0.0 - no margins, 1.0 - pure cosine
			double rad_to_hdiag) // 0.8 limit by radius too (1.0 - radius is half diagonal, same margins, use min)
	{			
		double mf_hor =  marg_fraction;
		double mf_vert = marg_fraction;
		double [] wx = new double [clustersX];
		double [] wy = new double [clustersY];
		double [] w = new double [clustersX * clustersY];
//		int [] boost_wnd = {33,15,40,35};
//		double boost_scale = 1.0; // 100.0;

		double centerX = 0.5 * (clustersX - 1);
		double marg = centerX * mf_hor;
		for (int x = 0; x <= clustersX / 2; x++) {
			if (x < marg) {
				wx[x] = Math.sin(Math.PI * x / 2.0 / marg);
				wx[x] *= wx[x];
			} else {
				wx[x] = 1.0;
			}
			wx[clustersX - 1 -x ] = wx[x];
		}
		double centerY = 0.5 * (clustersY - 1);
		marg = centerY * mf_vert;
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
				/*
				if (boost_scale > 1.0) {
					if ((x >= boost_wnd[0]) && (x < boost_wnd[2]) && (y >= boost_wnd[1]) && (y < boost_wnd[3])) {
						w[y * clustersX + x] *= boost_scale;
					}
				}
				*/
			}
		}
		double r_outer = Math.sqrt(centerX * centerX + centerY*centerY)* rad_to_hdiag; // relative to horizontal with
		double r_marg = r_outer * marg_fraction;
		for (int y = 0; y < clustersY; y++) {
			for (int x = 0; x < clustersX; x++) {
				double r = Math.sqrt((x - centerX)*(x - centerX) + (y - centerY)*(y - centerY));
				double k = 1.0;
				if (r > r_outer) {
					k = 0.0;
				} else if (r > (r_outer - r_marg)) {
					k = Math.sin(Math.PI * (r_outer - r) / 2.0 / r_marg);
				}
				if (k < w[y * clustersX + x]) {
					w[y * clustersX + x] = k;
				}
			}
		}
		return w;

	}
	
	
	// Make total weight of disparity (forced) samples - weight_disparity, weight of all others (8 per cluster) weight_lazyeye;
	private double []  distanceFromEdge (
			boolean []  force_disparity,
			double [][] measured_dsxy,
			double      min_dfe,
			double      max_dfe,
			boolean     debug) {
		TileNeibs tn = new TileNeibs(clustersX, clustersY);
		int [] idfe = tn.distanceFromEdge(force_disparity);
		double [] dfe = new double [idfe.length];
		int nft = 0;
		double aw = 0.0;
		for (int i = 0; i < idfe.length; i++) if ((idfe[i] > min_dfe) && (measured_dsxy[i]!= null)){
			nft++;
			double d = Math.min(idfe[i] - min_dfe, max_dfe); // max_dfe
			aw += d;
		}
		if (nft > 0) {
			aw /= nft;
			for (int i = 0; i < idfe.length; i++) if ((idfe[i] > min_dfe) && (measured_dsxy[i]!= null)){
				double d = Math.min(idfe[i] - min_dfe, max_dfe); // max_dfe
				dfe[i] = d / aw;
			}
			if (debug) {
			 (new ShowDoubleFloatArrays()).showArrays(
					 dfe,
					 clustersX,
					 clustersY,
					 "dist_from_edge_"+min_dfe+"-"+max_dfe);
			}
			
		}
		return dfe;
	}
	
	
	
	private boolean [] filterInfinity(
			double  [][] measured_dsxy,
			boolean [] force_disparity, // same dimension as dsdn, true if disparity should be controlled
			double [] distance_from_edge, // all 0?
			boolean en_tilt,  // consider right/left infinity tilt (will be disabled if more than *abs difference over width)
			int min_in_half,
			double min_infinity_abs,
			double max_infinity_abs,
			double min_infinity,
			double max_infinity
			) {
		if (min_in_half < 1) {
			min_in_half = 1;
		}
		int clusters = clustersX * clustersY;
	    double half_width = 0.5 * clustersX;
		boolean [] true_infinity = new boolean[clusters];
		double s0 = 0.0;
		double sx = 0.0;
		double sx2 = 0.0;
		double sy = 0.0;
		double sxy= 0.0;
		int nright = 0;
		int nleft = 0;
		for (int cluster = 0; cluster < clusters;  cluster++) if ((measured_dsxy[cluster] != null) && force_disparity[cluster]){
			double s = measured_dsxy[cluster][ExtrinsicAdjustment.INDX_STRENGTH] * weight_window[cluster]; // cluster weight
			if (distance_from_edge != null) {
				s *= distance_from_edge[cluster];
			}
			if (s > 0.0) {
				double d = measured_dsxy[cluster][ExtrinsicAdjustment.INDX_DIFF];
				if (( d >= min_infinity_abs) && ( d <= max_infinity_abs)) {
					double x = cluster % clustersX;
					s0 +=  s;
					sx +=  s * x;
					sx2 += s * x * x;
					sy +=  s * d;
					sxy += s *x * d;
					if (x >= half_width)  nright ++;
					else nleft ++;
				}
			}
		}
		if (s0 > 0.0) {
			double dnm = s0 * sx2 - sx*sx;
			double a = 0.0;
			double b = sy/s0;
			if (en_tilt && (nright >= min_in_half) && (nleft >= min_in_half)) {
				a = (sxy * s0 - sy * sx) / dnm;
				b = (sy * sx2 - sxy * sx) / dnm;
				if ((Math.abs(a) * clustersX) > (max_infinity_abs - min_infinity_abs)/2) {
					System.out.println("Infinity tilt too high ("+a+" > "+((max_infinity_abs - min_infinity_abs)/2)+
							"), disabling tilt");
					a = 0.0;
					b = sy/s0;
				}
			}
			for (int cluster = 0; cluster < clusters;  cluster++) if ((measured_dsxy[cluster] != null) && force_disparity[cluster]){
				if ((distance_from_edge == null) || (distance_from_edge[cluster] > 0.0)) {
					double x = cluster % clustersX;
					double d = measured_dsxy[cluster][ExtrinsicAdjustment.INDX_DIFF] - (a * x + b) ;
					true_infinity[cluster] =  (d >= min_infinity) && (d <= max_infinity);
				}
			}		
		}
		return true_infinity;
	}
	
	// find near foreground tiles with disparity exceeding average for the horizontal row
	
	private boolean [] selectERS(
			double  [][] measured_dsxy,
			boolean [] selection, // or null
			double min_fg_disp,
			double min_rel_over) { // 0.25 typical, always
		int clusters = clustersX * clustersY;
		boolean [] fg = new boolean[clusters];
		for (int row = 0; row < clustersY; row++) {
			double sw = 0.0;
			double swd =0.0;
			for (int col = 0; col < clustersX; col++) {
				int cluster = row*clustersX + col;
				if ((measured_dsxy[cluster] != null) && ((selection==null) || selection[cluster])){
					double s = measured_dsxy[cluster][ExtrinsicAdjustment.INDX_STRENGTH] * weight_window[cluster]; // cluster weight
					sw += s;
					swd += s * measured_dsxy[cluster][ExtrinsicAdjustment.INDX_DISP];
				}
			}
			if (sw > 0.0) {
				swd /= sw;
				double d_min = Math.max(min_fg_disp, swd*(1.0 + min_rel_over));
				for (int col = 0; col < clustersX; col++) {
					int cluster = row*clustersX + col;
					if ((measured_dsxy[cluster] != null) && ((selection==null) || selection[cluster])){
						fg[cluster] = measured_dsxy[cluster][ExtrinsicAdjustment.INDX_DISP] >= d_min;
					}
				}
			}
			
		}
		return fg;
	}

	private void boostERS(
			double []    weights,
			boolean []   fg,
			int          min_num_fg,
			double       fg_boost_fraction) {
		int num_fg = 0;
		double sw =  0.0;
		double swf = 0.0;
		int clusters = clustersX * clustersY;
		for (int cluster = 0; cluster < clusters; cluster++) if (measured_dsxy[cluster] != null){
			for (int i = 0; i < points_sample; i++) {
				int indx = cluster * points_sample + i;
				sw += weights[indx];
				if (fg[cluster]) {
					swf += weights[indx];
				}
			}
			if (fg[cluster]) {
				num_fg++;
			}
		}
		swf /= sw; // normally sw == 1.0;
		if ((num_fg >= min_num_fg) && (swf > 0.0)) {
			double k_fg = fg_boost_fraction/swf;
			double k_bg = (1.0 - fg_boost_fraction)/(1.0 - swf);
			for (int cluster = 0; cluster < clusters; cluster++) if (measured_dsxy[cluster] != null){
				for (int i = 0; i < points_sample; i++) {
					int indx = cluster * points_sample + i;
					weights[indx] *= fg[cluster] ? k_fg : k_bg;
				}
			}
		}
		return;
	}
	
	private void limitDisparity(
			double  [][] measured_dsxy,
			double []    weights,
			double max_disparity
			) {
		double sw =  0.0;
		double sw_near = 0.0;
		int clusters = clustersX * clustersY;
		for (int cluster = 0; cluster < clusters; cluster++) if (measured_dsxy[cluster] != null){
			double disparity = measured_dsxy[cluster][ExtrinsicAdjustment.INDX_TARGET]; // INDX_DISP];
			if (disparity <= max_disparity) {
				for (int i = 0; i < points_sample; i++) {
					sw_near += weights[cluster * points_sample + i];
				}
			}
			for (int i = 0; i < points_sample; i++) {
				sw += weights[cluster * points_sample + i];
			}
		}
		
		if (sw_near > 0.0) {
			double wboost = sw/sw_near;
			for (int cluster = 0; cluster < clusters; cluster++) if (measured_dsxy[cluster] != null){
				double disparity = measured_dsxy[cluster][ExtrinsicAdjustment.INDX_TARGET]; // INDX_DISP];
				if (disparity <= max_disparity) {
					for (int i = 0; i < points_sample; i++) {
						weights[cluster * points_sample + i] *= wboost;
					}
				} else {
					for (int i = 0; i < points_sample; i++) {
						weights[cluster * points_sample + i] = 0;
					}
				}
			}
		}
	}
	
	
	// Create mask so that moving objects are only allowed inside permitted area
	private boolean [] getMovingObjectsSel(
			double [][] measured_dsxy, // to use target_disparity
			double []   cluster_weights,
			double      max_disparity, // disable higher disparity
			double      rad_to_hdiag) // 0.8 limit by radius too (1.0 - radius is half diagonal, same margins, use min)
	{
		double centerX = 0.5 * (clustersX - 1);
		double centerY = 0.5 * (clustersY - 1);
		int clusters = clustersX * clustersY;
		double [] weights = cluster_weights.clone();
		double [] target_disparity = new double [clusters];
		boolean [] selection = new boolean [clusters];
		for (int cluster = 0; cluster < clusters; cluster++) {
			if (measured_dsxy[cluster] != null){
				target_disparity[cluster] = measured_dsxy[cluster][ExtrinsicAdjustment.INDX_TARGET];
			} else {
				weights[cluster] = 0.0;
			}
		}
		target_disparity = (new DoubleGaussianBlur()).blurWithNaN(
				target_disparity, // double[] pixels,
				weights, // null, // double [] in_weight, // or null
				clustersX,// int width,
				clustersY,// int height,
				0.5, // double sigmaX, slightly blur, mostly - just fill small gaps
				0.5, // double sigmaY,
				0.01); // double accuracy
		double radius = Math.sqrt(centerX * centerX + centerY*centerY)* rad_to_hdiag; // relative to horizontal with
		for (int y = 0; y < clustersY; y++) {
			for (int x = 0; x < clustersX; x++) {
				double r = Math.sqrt((x - centerX)*(x - centerX) + (y - centerY)*(y - centerY));
				if (r < radius) {
					int cluster = y * clustersX + x;
					if (target_disparity[cluster] <= max_disparity) {
						selection[cluster] = true;
					}
				}
			}
		}
		return selection;
	}

	
	private void blockSelectedClusters(
			double  [][] measured_dsxy,
			double  []   weights,
			boolean []   prohibit
			) {
		if (prohibit == null) {
			return;
		}
		double sw =  0.0;
		double sw_en = 0.0;
		int clusters = clustersX * clustersY;
		for (int cluster = 0; cluster < clusters; cluster++) if (measured_dsxy[cluster] != null){
			if (!prohibit[cluster]) {
				for (int i = 0; i < points_sample; i++) {
					sw_en += weights[cluster * points_sample + i];
				}
			}
			for (int i = 0; i < points_sample; i++) {
				sw += weights[cluster * points_sample + i];
			}
		}
		
		if (sw_en > 0.0) {
			double wboost = sw/sw_en;
			for (int cluster = 0; cluster < clusters; cluster++) if (measured_dsxy[cluster] != null){
				if (!prohibit[cluster]) {
					for (int i = 0; i < points_sample; i++) {
						weights[cluster * points_sample + i] *= wboost;
					}
				} else {
					for (int i = 0; i < points_sample; i++) {
						weights[cluster * points_sample + i] = 0;
					}
				}
			}
		}

	}
	
	@Deprecated
	private double [] getWeightsOld(
			double  [][] measured_dsxy,
			boolean [] force_disparity, // same dimension as dsdn, true if disparity should be controlled
			boolean [] filtered_infinity,
			int min_num_forced,  // if number of forced samples exceeds this, zero out weights of non-forced
			double weight_disparity, // now 0.0 - 1.0 fraction of disparity in all samples
			double weight_lazyeye)   // relative weight of disparity to 1/points_sample 
	{
//		weight_disparity =.5; // FIXME: Fix!
		int clusters = clustersX * clustersY;
		double [] weights = new double  [clusters * points_sample];
		boolean [] disable = new boolean [clusters];
		if (filtered_infinity != null) {
			for (int cluster = 0; cluster < clusters;  cluster++)	{
				disable[cluster] = force_disparity[cluster] && !filtered_infinity[cluster];
				}
			}
		
		double sw = 0.0;
		double swf = 0.0;
		int num_forced = 0;
//		if (force_disparity != null) for (int cluster = 0; cluster < clusters;  cluster++) if (force_disparity[cluster])num_forced ++;
//		boolean use_forced = num_forced >= min_num_forced;
		for (int cluster = 0; cluster < clusters;  cluster++) if ((measured_dsxy[cluster] != null) && !disable[cluster]){
			double s = measured_dsxy[cluster][ExtrinsicAdjustment.INDX_STRENGTH] * weight_window[cluster]; // cluster weight
			if ((force_disparity != null) && (force_disparity[cluster])) {
				swf += s;
				num_forced++;
			}
			sw += s;
		}
		if (sw <= 0.0) {
			return null;
		}
		boolean use_forced = num_forced >= min_num_forced;
		double k_inf_disp, k_inf_ly, k_other_disp, k_other_ly;
		if (use_forced) {
			// disparity weight non-zero only for infinity, non-disparity - equal scale for all
			k_inf_disp = weight_disparity / swf;
			k_inf_ly = (1.0 - weight_disparity) / sw / (points_sample - 1);
			k_other_disp = 0.0;
			k_other_ly = k_inf_ly;
		} else {
			k_other_disp = weight_lazyeye / sw / points_sample;
			k_other_ly = (1.0 - weight_lazyeye/points_sample) / (points_sample - 1) / sw;
			k_inf_disp = k_other_disp; 
			k_inf_ly = k_other_ly;
		}
		
		for (int cluster = 0; cluster < clusters;  cluster++) if ((measured_dsxy[cluster] != null) && !disable[cluster]){
			double s = measured_dsxy[cluster][ExtrinsicAdjustment.INDX_STRENGTH] * weight_window[cluster]; // cluster weight
			if ((force_disparity != null) && force_disparity[cluster]) {
				weights[cluster * points_sample + 0] = s * k_inf_disp;
				for (int i = 1; i < points_sample; i++) {
					weights[cluster * points_sample + i] = s * k_inf_ly;
				}
			} else {
				weights[cluster * points_sample + 0] = s * k_other_disp;
				for (int i = 1; i < points_sample; i++) {
					weights[cluster * points_sample + i] = s * k_other_ly;
				}
			}
		}		
		this.pure_weight = 1.0;
		return weights;
	}
	
	private double [][] detectClouds(
			double  [][] measured_dsxy,
			boolean [] force_disparity)   // same dimension as dsdn, true if disparity should be controlled
	{
		int clusters = clustersX * clustersY;
		if (force_disparity == null) {
			force_disparity = new boolean[clusters];
			for (int cluster = 0; cluster < clusters; cluster++) {
				if ((measured_dsxy[cluster] != null) && (measured_dsxy[cluster][ExtrinsicAdjustment.INDX_TARGET] == 0.0)) {
					force_disparity[cluster] = true;
				}
			}
		}
		
		double [][] clouds = new double [3][clusters];
		for (int clusterX = 0; clusterX < clustersX; clusterX++) {
			double sw = 0.0;
			double sdw = 0.0;
			for (int clusterY = 0; clusterY < clustersY; clusterY++) {
				int cluster = clusterX + clustersX * clusterY;
				if (force_disparity[cluster] && (measured_dsxy[cluster] != null)) {
					double s = measured_dsxy[cluster][ExtrinsicAdjustment.INDX_STRENGTH];
					double d = -measured_dsxy[cluster][ExtrinsicAdjustment.INDX_DIFF];
					sw += s;
					sdw += s * d;
				}
				if (sw >= 0) {
					clouds[0][cluster] = sdw/sw;
				} else {
					clouds[0][cluster] = Double.NaN;
				}
				if (measured_dsxy[cluster] != null) {
					clouds[1][cluster] = -measured_dsxy[cluster][ExtrinsicAdjustment.INDX_DIFF];
					clouds[2][cluster] = measured_dsxy[cluster][ExtrinsicAdjustment.INDX_STRENGTH];
				}
			}
			for (int cluster = 0; cluster < clusters; cluster++) {
				if (!force_disparity[cluster]) {
					clouds[0][cluster] = Double.NaN;
				}
			}
		}
		return clouds;
	}
	@Deprecated
	private double [] getWeights(
			double  [][] measured_dsxy,
			boolean [] force_disparity,   // same dimension as dsdn, true if disparity should be controlled
			boolean [] filtered_infinity, // tiles known to be infinity
			double  []  distance_from_edge,// to reduce weight of the mountain ridge, increase clouds (or null)
			int min_num_forced,           // if number of forced samples exceeds this, zero out weights of non-forced
			boolean infinity_right_left,  // each halve should have > min_num_forced, will calculate separate average
			double weight_infinity,       // total weight of infinity tiles fraction (0.0 - 1.0) 
			double weight_disparity,      // disparity weight relative to the sum of 8 lazy eye values of the same tile 
			double weight_disparity_inf,  // disparity weight relative to the sum of 8 lazy eye values of the same tile for infinity 
			double max_disparity_far,     // reduce weights of near tiles proportional to sqrt(max_disparity_far/disparity)
			double max_disparity_use) // do not process near objects to avoid ERS
	{
		int clusters = clustersX * clustersY;
		double [] weights = new double  [clusters * points_sample];
		double [][] strength = new double [2][clusters]; // 0 - for disparity, 1 - for ly
		boolean [] disable = new boolean [clusters];
		if (filtered_infinity != null) {
			for (int cluster = 0; cluster < clusters;  cluster++)	{
				disable[cluster] = force_disparity[cluster] && !filtered_infinity[cluster];
			}
		}
		double sw = 0.0;
		double swf = 0.0;
		int num_forced = 0;
		for (int cluster = 0; cluster < clusters;  cluster++) if ((measured_dsxy[cluster] != null) && !disable[cluster]){
			double [] sdl =getClustWeight(
					cluster,
					measured_dsxy,
					force_disparity,   // same dimension as dsdn, true if disparity should be controlled
					distance_from_edge,
					max_disparity_far,     // reduce weights of near tiles proportional to sqrt(max_disparity_far/disparity)
					max_disparity_use);
			if ((force_disparity != null) && (force_disparity[cluster])) {
				strength[0][cluster] = weight_disparity_inf * sdl[0];
				strength[1][cluster] = (1.0 - weight_disparity_inf) * sdl[1];
				swf += strength[0][cluster] + strength[1][cluster];
				num_forced++;
			} else {
				strength[0][cluster] = weight_disparity * sdl[0];
				strength[1][cluster] = (1.0 - weight_disparity) * sdl[1];
			}
			
			sw += strength[0][cluster] + strength[1][cluster];
		}
		if (sw <= 0.0) {
			return null;
		}
		boolean use_forced = num_forced >= min_num_forced;
		double k_inf_disp, k_inf_ly, k_other_disp, k_other_ly;
		if (sw <= swf) {
			use_forced = false;
		}
		if (use_forced) {
			k_inf_disp =   weight_infinity / swf;
			k_inf_ly =     weight_infinity / swf / (points_sample - 1);
			k_other_disp = (1.0 - weight_infinity) / (sw - swf);
			k_other_ly =   (1.0 - weight_infinity) / (sw - swf) / (points_sample - 1);
			
		} else {
			k_other_disp = 1.0  / sw;
			k_other_ly =   1.0 / sw / (points_sample - 1);
			k_inf_disp =   k_other_disp; 
			k_inf_ly =     k_other_ly;
		}

		for (int cluster = 0; cluster < clusters;  cluster++) if ((measured_dsxy[cluster] != null) && !disable[cluster]){
			double sd = strength[0][cluster];
			double sly =strength[1][cluster];
			if ((force_disparity != null) && force_disparity[cluster]) {
				weights[cluster * points_sample + 0] = sd * k_inf_disp;
				for (int i = 1; i < points_sample; i++) {
					weights[cluster * points_sample + i] = sly * k_inf_ly;
				}
			} else {
				weights[cluster * points_sample + 0] = sd * k_other_disp;
				for (int i = 1; i < points_sample; i++) {
					weights[cluster * points_sample + i] = sly * k_other_ly;
				}
			}
		}		
		this.pure_weight = 1.0;
		return weights;
	}
	
	private int [] setWeights( // number right, number left
			double  [][] measured_dsxy,
			boolean [] force_disparity,   // same dimension as dsdn, true if disparity should be controlled
			boolean [] filtered_infinity, // tiles known to be infinity
			double  []  distance_from_edge,// to reduce weight of the mountain ridge, increase clouds (or null)
			int min_num_forced,           // if number of forced samples exceeds this, zero out weights of non-forced
			boolean infinity_right_left,  // each halve should have > min_num_forced, will calculate separate average
			double weight_infinity,       // total weight of infinity tiles fraction (0.0 - 1.0) 
			double weight_disparity,      // disparity weight relative to the sum of 8 lazy eye values of the same tile 
			double weight_disparity_inf,  // disparity weight relative to the sum of 8 lazy eye values of the same tile for infinity 
			double max_disparity_far,     // reduce weights of near tiles proportional to sqrt(max_disparity_far/disparity)
			double max_disparity_use) // do not process near objects to avoid ERS
	{
		int clusters = clustersX * clustersY;
	    double half_width = 0.5 * clustersX;
		double [] weights = new double  [clusters * points_sample];
		double [][] strength = new double [2][clusters]; // 0 - for disparity, 1 - for ly
		boolean [] disable = new boolean [clusters];
		if (filtered_infinity != null) {
			for (int cluster = 0; cluster < clusters;  cluster++)	{
				disable[cluster] = force_disparity[cluster] && !filtered_infinity[cluster];
			}
		}
		double sw = 0.0;
		double swf = 0.0;
		int num_forced_right = 0;
		int num_forced_left = 0;
		for (int cluster = 0; cluster < clusters;  cluster++) if ((measured_dsxy[cluster] != null) && !disable[cluster]){
			double [] sdl =getClustWeight(
					cluster,
					measured_dsxy,
					force_disparity,   // same dimension as dsdn, true if disparity should be controlled
					distance_from_edge,
					max_disparity_far,     // reduce weights of near tiles proportional to sqrt(max_disparity_far/disparity)
					max_disparity_use);
			if ((force_disparity != null) && (force_disparity[cluster])) {
				strength[0][cluster] = weight_disparity_inf * sdl[0];
				strength[1][cluster] = (1.0 - weight_disparity_inf) * sdl[1];
				swf += strength[0][cluster] + strength[1][cluster];
				double x = cluster % clustersX;
				if (x >= half_width)  num_forced_right ++;
				else                  num_forced_left ++;
			} else {
				strength[0][cluster] = weight_disparity * sdl[0];
				strength[1][cluster] = (1.0 - weight_disparity) * sdl[1];
			}
			
			sw += strength[0][cluster] + strength[1][cluster];
		}
		int num_forced =     num_forced_right + num_forced_left;
		int num_forced_min = (num_forced_right > num_forced_left) ? num_forced_left : num_forced_right;
		if (sw <= 0.0) {
			return null;
		}
		boolean use_forced = num_forced >= min_num_forced;
		infinity_right_left &= (num_forced_min > (min_num_forced/2)); // do not try to balance if at least one half has too few
		double k_inf_disp, k_inf_ly, k_other_disp, k_other_ly;
		if (sw <= swf) {
			use_forced = false;
		}
		if (use_forced) {
			k_inf_disp =   weight_infinity / swf;
			k_inf_ly =     weight_infinity / swf / (points_sample - 1);
			k_other_disp = (1.0 - weight_infinity) / (sw - swf);
			k_other_ly =   (1.0 - weight_infinity) / (sw - swf) / (points_sample - 1);
			
		} else {
			k_other_disp = 1.0  / sw;
			k_other_ly =   1.0 / sw / (points_sample - 1);
			k_inf_disp =   k_other_disp; 
			k_inf_ly =     k_other_ly;
		}

		for (int cluster = 0; cluster < clusters;  cluster++) if ((measured_dsxy[cluster] != null) && !disable[cluster]){
			double sd = strength[0][cluster];
			double sly =strength[1][cluster];
			if ((force_disparity != null) && force_disparity[cluster]) {
				weights[cluster * points_sample + 0] = sd * k_inf_disp;
				for (int i = 1; i < points_sample; i++) {
					weights[cluster * points_sample + i] = sly * k_inf_ly;
				}
			} else {
				weights[cluster * points_sample + 0] = sd * k_other_disp;
				for (int i = 1; i < points_sample; i++) {
					weights[cluster * points_sample + i] = sly * k_other_ly;
				}
			}
		}
		// optionally balance right/left for the infinity
		if (infinity_right_left && (force_disparity != null)) {
			double swdr = 0.0;
			double swdl = 0.0;
			double swlyr = 0.0;
			double swlyl = 0.0;
			for (int cluster = 0; cluster < clusters;  cluster++) {
				if ((measured_dsxy[cluster] != null) && force_disparity[cluster] && !disable[cluster]){
					double x = cluster % clustersX;
					double wd = weights[cluster * points_sample + 0];
					double wly = 0;
					for (int i = 1; i < points_sample; i++) {
						wly += weights[cluster * points_sample + i];
					}
					if (x >= half_width) {
						swdr += wd;
						swlyr += wly;
					} else {
						swdl += wd;
						swlyl += wly;
					}
				}
			}
			double kwdr = (swdr+swdl)/(2*swdr);
			double kwdl = (swdr+swdl)/(2*swdl);
			double kwlyr = (swlyr+swlyl)/(2*swlyr);
			double kwlyl = (swlyr+swlyl)/(2*swlyl);
			for (int cluster = 0; cluster < clusters;  cluster++) {
				if ((measured_dsxy[cluster] != null) && force_disparity[cluster] && !disable[cluster]){
					double x = cluster % clustersX;
					if (x >= half_width) {
						weights[cluster * points_sample + 0] *= kwdr;
						for (int i = 1; i < points_sample; i++) {
							weights[cluster * points_sample + i] *= kwlyr;
						}
					} else {
						weights[cluster * points_sample + 0] *= kwdl;
						for (int i = 1; i < points_sample; i++) {
							weights[cluster * points_sample + i] *= kwlyl;
						}
					}
				}
			}
		}
		this.pure_weight = 1.0;
		this.weights = weights;
		return new int [] {num_forced, num_forced_min};
	}

	
	
	
	private double []  getClustWeight( // {weight for disparity, weight for ly}
			int cluster,
			double  [][] measured_dsxy,
			boolean [] force_disparity,   // same dimension as dsdn, true if disparity should be controlled
			double [] distance_from_edge,
			double max_disparity_far,     // reduce weights of near tiles proportional to sqrt(max_disparity_far/disparity)
			double max_disparity_use) // do not process near objects to avoid ERS
		{
		double s = measured_dsxy[cluster][ExtrinsicAdjustment.INDX_STRENGTH];
		double s1 = Double.NaN;
		if (distance_from_edge != null) { // use vignetting for infinity too
			s *= weight_window[cluster];
			if ((force_disparity != null) && force_disparity[cluster]) {
				s1 = s;
				s *= distance_from_edge[cluster];
			} 
		} else { // do not use vignetting for infinity ?
			if ((force_disparity == null) || (!force_disparity[cluster])) {
				s *= weight_window[cluster];
			} 
		}
		double d = measured_dsxy[cluster][ExtrinsicAdjustment.INDX_TARGET];
		if (d > max_disparity_far) {
			s *= Math.sqrt(max_disparity_far/d);
		}
		if (d > max_disparity_use) {
			s = 0.0;
		}
		if (Double.isNaN(s1)) {
			s1 = s;
		}
		return new double[] {s,s1};
		
	}

	private double [] getFx(
			CorrVector corr_vector)
	{
		int clusters = clustersX * clustersY;
		Matrix []   corr_rots =  corr_vector.getRotMatrices(); // get array of per-sensor rotation matrices
		Matrix [][] deriv_rots = corr_vector.getRotDeriveMatrices();
		double [] imu = corr_vector.getIMU(); // i)
		double [] y_minus_fx = new double  [clusters * points_sample];
		for (int cluster = 0; cluster < clusters;  cluster++) {
			if (measured_dsxy[cluster] != null){
				double [] ddnd = geometryCorrection.getPortsDDNDAndDerivativesNew( // USED in lwir
						geometryCorrection,     // GeometryCorrection gc_main,
						use_rig_offsets,        // boolean     use_rig_offsets,
						corr_rots,              // Matrix []   rots,
						deriv_rots,             // Matrix [][] deriv_rots,
						null,                   // double [][] DDNDderiv,     // if not null, should be double[8][]
						pY_offset[cluster],     // double []   py_offset,  // array of per-port average pY offset from the center (to correct ERS) or null (for no ERS)
						imu,                    // double []   imu,
						x0y0[cluster],          // double []   pXYND0,        // per-port non-distorted coordinates corresponding to the correlation measurements
						world_xyz[cluster],     // double []   xyz, // world XYZ for ERS correction
						measured_dsxy[cluster][ExtrinsicAdjustment.INDX_PX + 0],  // double      px,
						measured_dsxy[cluster][ExtrinsicAdjustment.INDX_PX + 1],  // double      py,
						measured_dsxy[cluster][ExtrinsicAdjustment.INDX_TARGET]); // double      disparity);
				//arraycopy(Object src, int srcPos, Object dest, int destPos, int length)
				//		    System.arraycopy(src_pixels, 0, dst_pixels, 0, src_pixels.length); /* for the borders closer to 1/2 kernel size*/
				ddnd[0] = ddnd[0]; // ?
				for (int i = 0; i < num_sensors; i++) {
					ddnd[i + 1] = ddnd[i + 1];
					ddnd[i + 5] = ddnd[i + 5];
				}
				System.arraycopy(ddnd, 0, y_minus_fx, cluster*points_sample, points_sample);
			}
		}
		return y_minus_fx;
	}

	private double [][] getJacobianTransposed(
			CorrVector corr_vector,
			double delta){
		int clusters = clustersX * clustersY;
		int num_pars = getNumPars();
		double [][] jt = new double  [num_pars][clusters * points_sample ];
		double rdelta = 1.0/delta;
		for (int par = 0; par < num_pars; par++) {
			double [] pars = new double[num_pars];
			pars[par] =  delta;
			CorrVector corr_delta = geometryCorrection.getCorrVector(pars, par_mask);
			CorrVector corr_vectorp = corr_vector.clone();
			CorrVector corr_vectorm = corr_vector.clone();
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

	public String [] getSymNames() {
		int num_pars = getNumPars();
		String [] names = new String[num_pars];
		int indx = 0;
		for (int i = 0; i < par_mask.length; i++) if (par_mask[i]){
			names[indx++] = "S"+i;
		}
		return names;
	}

	// get full (s0..s18) index from jacobian parameter number 
	public int getFullParIndex(int indx) {
		for (int i = 0; i < par_mask.length; i++) if (par_mask[i]){
			if (indx == 0) {
				return i;
			}
			indx --;
		}
		return -1;
	}
	public int getParIndex(int full_index) {
		int indx = 0;
		for (int i = 0; i < par_mask.length; i++) if (par_mask[i]){
			if (i == full_index) {
				return indx;
			}
			indx ++;
		}
		return -1;
	}
	
	

	private double dbgJacobians(
			CorrVector corr_vector,
			double delta,
			boolean graphic) {
//		delta *= 0.1;
		int gap = 10; //pix
		int clusters = clustersX * clustersY;
		int num_pars = getNumPars();
		String [] titles = getSymNames();
		double [][] jt =        getJacobianTransposed(corr_vector);
		double [][] jt_delta =  getJacobianTransposed(corr_vector, delta);
		double tot_error = 0;
		double [][] err =   new double [num_pars][points_sample];
		double [] err_par = new double [num_pars];
		for (int par = 0; par < num_pars; par++) {
			for (int cluster = 0; cluster < clusters;  cluster++) if (measured_dsxy[cluster] != null){
				for (int i = 0; i < points_sample; i++) {
					int indx = cluster *  points_sample+ i;
					if (Math.abs(jt[par][indx] - jt_delta[par][indx]) > err[par][i]) {
						err[par][i] = Math.abs(jt[par][indx] - jt_delta[par][indx]);
					}
				}
			}
			for (int i = 0; i < points_sample; i++) {
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
			int rows = getRowsCols()[0];
			int cols = getRowsCols()[1];
			int width  = cols * (clustersX + gap) - gap;
			int height = rows * (clustersY + gap) - gap;
			double [][] dbg_img = new double [num_pars * 3][width*height];
			for (int par = 0; par < num_pars; par++) {
				titles3[3 * par + 0] = titles[par]+"";
				titles3[3 * par + 1] = titles[par]+"_delta";
				titles3[3 * par + 2] = titles[par]+"_diff";
				for (int mode = 0; mode < points_sample; mode++) {
					int x0 = (mode % cols) * (clustersX + gap);
					int y0 = (mode / cols) * (clustersY + gap);
					for (int cluster = 0; cluster < clusters;  cluster++) {
						int x = x0 + (cluster % clustersX);
						int y = y0 + (cluster / clustersX);
						int pix = x + y * width;
						int indx = cluster * points_sample + mode;

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
			 
				dbg_img = new double [num_pars][width*height];
				for (int par = 0; par < num_pars; par++) {
					for (int mode = 0; mode < points_sample; mode++) {
						int x0 = (mode % cols) * (clustersX + gap);
						int y0 = (mode / cols) * (clustersY + gap);
						for (int cluster = 0; cluster < clusters;  cluster++) {
							int x = x0 + (cluster % clustersX);
							int y = y0 + (cluster / clustersX);
							int pix = x + y * width;
							int indx = cluster * points_sample + mode;

							if (measured_dsxy[cluster] != null){
								dbg_img[par][pix] =  jt[par][indx];
							} else {
								dbg_img[par][pix] =  Double.NaN;
							}
						}
					}
				}
				 (new ShowDoubleFloatArrays()).showArrays(
						 dbg_img,
						 width,
						 height,
						 true,
						 "Jacobians",
						 titles);
			 
			 
		}
		return tot_error;
	}

	private void dbgYminusFx(
			double []   fx,
			String title) {
		int gap = 10; //pix
		int clusters = clustersX * clustersY;
		String [] titles = {"Y", "-fX", "Y+fx"};
		int rows = getRowsCols()[0];
		int cols = getRowsCols()[1];
		int width  = cols * (clustersX + gap) - gap;
		int height = rows * (clustersY + gap) - gap;
		double [][] dbg_img = new double [3][width*height];
		for (int mode = 0; mode < points_sample; mode++) {
			int x0 = (mode % cols) * (clustersX + gap);
			int y0 = (mode / cols) * (clustersY + gap);
			for (int cluster = 0; cluster < clusters;  cluster++) {
				int x = x0 + (cluster % clustersX);
				int y = y0 + (cluster / clustersX);
				int pix = x + y * width;
				int indx = cluster * points_sample + mode;
				if (measured_dsxy[cluster] != null){
					if (mode ==0) {
						dbg_img[0][pix] =  measured_dsxy[cluster][ExtrinsicAdjustment.INDX_DIFF];
						dbg_img[1][pix] = -fx[indx];
						dbg_img[2][pix] =  measured_dsxy[cluster][ExtrinsicAdjustment.INDX_DIFF] + fx[indx];
					} else {
//						dbg_img[0][pix] =  measured_dsxy[cluster][ExtrinsicAdjustment.INDX_DD0 + mode - 1];
						dbg_img[0][pix] =  measured_dsxy[cluster][indx_dd0 + mode - 1];
						dbg_img[1][pix] = -fx[indx];
//						dbg_img[2][pix] =  measured_dsxy[cluster][ExtrinsicAdjustment.INDX_DD0 + mode - 1] + fx[indx];
						dbg_img[2][pix] =  measured_dsxy[cluster][indx_dd0 + mode - 1] + fx[indx];
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

	private boolean [] selectMovingObjects(
			 double sigma,             //  = 1.0;
			 double max_mov_disparity, //  = 75.0;
			 double rad_to_hdiag_mov,  //  = 0.7 ; // 0.8
			 double max_mov_average,   //  = .25;
			 double min_l2,            //  = 0.75;
			 String title,             // "Moving_objects"
			 String title_extra,       // "Moving_objects_filtering"
			 boolean debug_text
			) {
		int clusters = clustersX * clustersY;
		double max_reasonable = 100.0;
		
		 double [][] ers_tilt_az = getMovingObjects( // will output null if ERS tilt or roll are disabled
				 corr_vector, // CorrVector corr_vector,
				 null ); // double [] fx
		 boolean [] moving_maybe = null;
		 if (ers_tilt_az != null) {
			 double [] l2_moving = new double [clusters];
			 for (int cluster = 0; cluster < clusters; cluster++) {
				 if (ers_tilt_az[cluster] != null) {
					 l2_moving[cluster] = Math.sqrt(0.5*(ers_tilt_az[cluster][0]*ers_tilt_az[cluster][0] + ers_tilt_az[cluster][1]*ers_tilt_az[cluster][1]));
					 if (l2_moving[cluster] > max_reasonable) {
						 l2_moving[cluster] = Double.NaN;
					 }
				 } else {
					 l2_moving[cluster] = Double.NaN;
				 }
			 }
			 double [] ers_blured = l2_moving.clone();
			 double [] cluster_weights = getClusterWeight();
			 if (sigma > 0.0){
				 ers_blured = (new DoubleGaussianBlur()).blurWithNaN(
						 ers_blured, // double[] pixels,
						 cluster_weights, // null, // double [] in_weight, // or null
						 clustersX,// int width,
						 clustersY,// int height,
						 sigma, // double sigmaX,
						 sigma, // double sigmaY,
						 0.01); // double accuracy
			 }
			 // get area permitted for moving objects (not too close to the edges)
			 boolean [] moving_en= getMovingObjectsSel(
					 measured_dsxy, // to use target_disparity
					 cluster_weights,
					 max_mov_disparity, // disable higher disparity
					 rad_to_hdiag_mov); // 0.8 limit by radius too (1.0 - radius is half diagonal, same margins, use min)
			 double sw = 0.0;
			 double mov_avg = 0.0;
			 for (int cluster = 0; cluster < clusters; cluster++) {
				 if (moving_en[cluster] && !Double.isNaN(ers_blured[cluster])) {
					 mov_avg += ers_blured[cluster] * cluster_weights[cluster];
					 sw += cluster_weights[cluster];
				 }
			 }
			 mov_avg /= sw; // clusters;
			 if (debug_text ) {
				 System.out.println("Average moving object detection value = "+mov_avg+
						 ",  max_mov_average ="+ max_mov_average);
			 }
			 if (mov_avg <= max_mov_average) {
				 //			 double min_l2 = 1.0 * mov_avg; 
				 //				 double min_l2 = 5.0 * mov_avg; //4.0 * mov_avg; //3.50 * mov_avg; // 2.50 * mov_avg; // 1.5 * mov_avg; // make absolute?
				 moving_maybe =  extractMovingObjects(
						 title_extra, //  String    title,
						 ers_blured,  //double [] ers_blured, // has NaNs
						 moving_en,   // boolean [] moving_en
						 min_l2,      // double min_l2); // minimal value (can use average or fraction of it?
						 4,           // int       shrink, // prevents growing outer border (==4)
						 2);          // int       grow) { // located moving areas

				 int num_moving = 0;
				 for (int cluster = 0; cluster < clusters; cluster++) {
					 if (moving_maybe[cluster]) num_moving++;
				 }
				 if (debug_text ) {
					 System.out.println("Number of clusters to remove from LMA fitting = "+num_moving);
				 }
			 }
			 if (title != null) {
				 String [] titles = {"vert", "hor", "amplitude", "blured", "mask","maybe"};
				 double [][] dbg_img = new double [titles.length][] ; // [clusters];
				 dbg_img[0] = new double [clusters];
				 dbg_img[1] = new double [clusters];
				 dbg_img[2] = new double [clusters];
				 for (int cluster = 0; cluster < clusters; cluster++) {
					 if (ers_tilt_az[cluster] != null) {
							dbg_img[0][cluster] = ers_tilt_az[cluster][0];
							dbg_img[1][cluster] = ers_tilt_az[cluster][1];
							dbg_img[2][cluster] = l2_moving[cluster];
					 } else {
							dbg_img[0][cluster] = Double.NaN;
							dbg_img[1][cluster] = Double.NaN;
							dbg_img[2][cluster] = Double.NaN;
					 }
				 }
				 dbg_img[3] = ers_blured;
				 dbg_img[4] = new double [clusters];
				 for (int cluster = 0; cluster < clusters; cluster++) {
					 dbg_img[4][cluster] = moving_en[cluster] ? 1.0:0.0;
				 }
				 if (moving_maybe != null) {
					 dbg_img[5] = new double [clusters];
					 for (int cluster = 0; cluster < clusters; cluster++) {
						 dbg_img[5][cluster] = moving_maybe[cluster] ? 1.0:0.0;
					 }
				 }
				 (new ShowDoubleFloatArrays()).showArrays(
						 dbg_img,
						 clustersX,
						 clustersY,
						 true,
						 title,
						 titles);
			 }
		 }
		 return moving_maybe;
	}
	
	
	private double [][] getMovingObjects(
			CorrVector corr_vector,
			double [] fx // or null;
			) {
		int ers_tilt_index =    getParIndex(corr_vector. getIMUIndex());
		int ers_azimuth_index = getParIndex(corr_vector. getIMUIndex() + 1);
		if ((ers_tilt_index < 0) || (ers_azimuth_index < 0)) {
			return null;
		}
		int clusters = clustersX * clustersY;
		double [][] ers_tilt_az = new double [clusters][];
		double [][] jt = getJacobianTransposed(corr_vector);
		for (int cluster = 0; cluster < clusters; cluster++) if ( measured_dsxy[cluster] != null) {
			ers_tilt_az[cluster] = new double [2];
			for (int i = 0; i < (points_sample-1); i++) {
				int indx = cluster * points_sample + i + 1; // skipping disparity
//				double d = measured_dsxy[cluster][ExtrinsicAdjustment.INDX_DD0 + i];
				double d = measured_dsxy[cluster][indx_dd0 + i];
				if (fx != null) {
					d += fx[indx];
				}
				ers_tilt_az[cluster][0] += d * jt[ers_tilt_index][indx];
				ers_tilt_az[cluster][1] += d * jt[ers_azimuth_index][indx];
			}
		}
		return ers_tilt_az;
	}


	private double [] getClusterWeight() { // sum == 1.0;
		int clusters = clustersX * clustersY;
		double [] cluster_weights = new double [clusters];
		double sw = 0.0;
		for (int cluster = 0; cluster < clusters; cluster++) {
			for (int i = 0; i < points_sample; i++) {
				cluster_weights[cluster] += this.weights[cluster * points_sample + i];
				sw += cluster_weights[cluster];
			}
		}
		if (sw > 0.0){ // normally should be == 1.0
			double k = 1.0/sw;
			for (int cluster = 0; cluster < clusters; cluster++) {
				cluster_weights[cluster] *= k;
			}			
		}
		return cluster_weights;
		
	}
	
	private void showMovingObjects(
			String      title,
			double [][] ers_tilt_az) {
		showMovingObjects(
				title,
				ers_tilt_az,
				0, // double sigma,
				0,
				true);
	}

	private double [] showMovingObjects(
			String      title,
			double [][] ers_tilt_az,
			double sigma,
			double accuracy,
			boolean debug
			)
	{
		int clusters = clustersX * clustersY;
		
		String [] titles = {"vert", "hor", "amplitude"};
		double [][] dbg_img = new double [titles.length][clusters];
		double max_reasonable = 100.0;
		for (int cluster = 0; cluster < clusters; cluster++) {
			if (ers_tilt_az[cluster] != null) {
				dbg_img[0][cluster] = ers_tilt_az[cluster][0];
				dbg_img[1][cluster] = ers_tilt_az[cluster][1];
				dbg_img[2][cluster] = Math.sqrt(0.5*(ers_tilt_az[cluster][0]*ers_tilt_az[cluster][0] + ers_tilt_az[cluster][1]*ers_tilt_az[cluster][1]));
				if (dbg_img[2][cluster] > max_reasonable) {
					dbg_img[0][cluster] = Double.NaN;
					dbg_img[1][cluster] = Double.NaN;
					dbg_img[2][cluster] = Double.NaN;
				}
			} else {
				dbg_img[0][cluster] = Double.NaN;
				dbg_img[1][cluster] = Double.NaN;
				dbg_img[2][cluster] = Double.NaN;
			}
		}
		
		if (sigma > 0.0){
			double [] cluster_weights = getClusterWeight();
			for (int i = 0; i < dbg_img.length; i++) {
				dbg_img[i] = (new DoubleGaussianBlur()).blurWithNaN(
						dbg_img[i], // double[] pixels,
						cluster_weights, // null, // double [] in_weight, // or null
						clustersX,// int width,
						clustersY,// int height,
						sigma, // double sigmaX,
						sigma, // double sigmaY,
						0.01); // double accuracy
			}
		}
		
		(new ShowDoubleFloatArrays()).showArrays(
				dbg_img,
				clustersX,
				clustersY,
				true,
				title,
				titles);
		return dbg_img[2]; 

	}
	
	public boolean [] extractMovingObjects(
			String    title,
			double [] ers_blured, // has NaNs
			boolean [] moving_en,
			double    min_l2, // minimal value (can use average or fraction of it?
			int       shrink, // prevents growing outer border (==4)
			int       grow) { // located moving areas
		int clusters = clustersX * clustersY;
		TileNeibs tn = new TileNeibs(clustersX, clustersY);
		boolean [] terrain = new boolean[clusters];
		for (int i = 0; i < clusters; i++) {
			terrain[i] = ers_blured[i] <= min_l2; // NaNs not included
//			if (moving_en != null) {
//				terrain[i] &= moving_en[i]; 
//			}
		}
		// find largest terrain area
		int [] iterrains = tn. enumerateClusters(
				terrain, // boolean [] tiles,
				true); //boolean ordered)
		boolean [] terrain_main = new boolean[clusters];
		for (int i = 0; i < clusters; i++) {
			terrain_main[i] = (iterrains[i] == 1); // largest cluster
		}
		
		boolean [] filled = tn.fillConcave(terrain_main);
		
		tn.shrinkSelection(
				clustersX, // shrink, // int        shrink, remove all open outside of moving_en
				filled, // boolean [] tiles,
				moving_en); // null); // boolean [] prohibit)
		boolean [] moving = filled.clone();
		for (int i = 0; i < clusters; i++) {
			moving[i] &= !terrain_main[i]; 
		}
		boolean [] grown = moving.clone();
		tn.growSelection(
				grow, // int        grow,           // grow tile selection by 1 over non-background tiles 1: 4 directions, 2 - 8 directions, 3 - 8 by 1, 4 by 1 more
				grown, // boolean [] tiles,
				null); // boolean [] prohibit)
		if (title != null) {
			String [] titles = {"terrain","clusters","main_terrain","filled","moving","grown"};
			double [][] dbg_img = new double [titles.length][clusters];
			for (int cluster = 0; cluster < clusters; cluster++) {
				dbg_img[0][cluster] = terrain[cluster]? 1.0 : 0.0;
				dbg_img[1][cluster] = iterrains[cluster];
				dbg_img[2][cluster] = terrain_main[cluster]? 1.0 : 0.0;
				dbg_img[3][cluster] = filled[cluster]? 1.0 : 0.0;
				dbg_img[4][cluster] = moving[cluster]? 1.0 : 0.0;
				dbg_img[5][cluster] = grown[cluster]? 1.0 : 0.0;
			}
			(new ShowDoubleFloatArrays()).showArrays(
					dbg_img,
					clustersX,
					clustersY,
					true,
					title,
					titles);
		}
		return grown;
	}
	
	
	
	
	public int [] getRowsCols() {
		int rows = (int) Math.round(Math.sqrt(points_sample));
		int cols = points_sample/rows;
		if (rows*cols < points_sample) {
			cols++;
		}
		return new int [] {rows,cols};
	}

	private void dbgYminusFxWeight(
			double []   fx,
			double []   weights,
			String title) {
		if (fx == null) {
			fx = new double [weights.length];
		}
		int gap = 10; //pix
		int clusters = clustersX * clustersY;
		String [] titles = {"Y", "-fX", "Y+fx", "Weight", "W*(Y+fx)", "Masked Y+fx"};
		int rows = getRowsCols()[0];
		int cols = getRowsCols()[1];
		int width  = cols * (clustersX + gap) - gap;
		int height = rows * (clustersY + gap) - gap;
		double [][] dbg_img = new double [titles.length][width*height];
		for (int mode = 0; mode < points_sample; mode++) {
			int x0 = (mode % cols) * (clustersX + gap);
			int y0 = (mode / cols) * (clustersY + gap);
			for (int cluster = 0; cluster < clusters;  cluster++) {
				int x = x0 + (cluster % clustersX);
				int y = y0 + (cluster / clustersX);
				int pix = x + y * width;
				int indx = cluster * points_sample + mode;
				if (measured_dsxy[cluster] != null){
					if (mode ==0) {
						dbg_img[0][pix] =  measured_dsxy[cluster][ExtrinsicAdjustment.INDX_DIFF];
						dbg_img[1][pix] = -fx[indx];
						dbg_img[2][pix] =  measured_dsxy[cluster][ExtrinsicAdjustment.INDX_DIFF] + fx[indx];
						if (weights[indx] > 0.0) {
							dbg_img[3][pix] =  weights[indx]*clusters;
							dbg_img[4][pix] =  weights[indx]*clusters*(measured_dsxy[cluster][ExtrinsicAdjustment.INDX_DIFF] + fx[indx]);
							dbg_img[5][pix] =  measured_dsxy[cluster][ExtrinsicAdjustment.INDX_DIFF] + fx[indx];
						} else {
							dbg_img[3][pix] =  Double.NaN;
							dbg_img[4][pix] =  Double.NaN;
							dbg_img[5][pix] =  Double.NaN;
						}
					} else {
//						dbg_img[0][pix] =  measured_dsxy[cluster][ExtrinsicAdjustment.INDX_DD0 + mode - 1];
						dbg_img[0][pix] =  measured_dsxy[cluster][indx_dd0 + mode - 1];
						dbg_img[1][pix] = -fx[indx];
//						dbg_img[2][pix] =  measured_dsxy[cluster][ExtrinsicAdjustment.INDX_DD0 + mode - 1] + fx[indx];
						dbg_img[2][pix] =  measured_dsxy[cluster][indx_dd0 + mode - 1] + fx[indx];
						if (weights[indx] > 0.0) {
							dbg_img[3][pix] =  weights[indx]*clusters;
//							dbg_img[4][pix] =  weights[indx]*clusters*(measured_dsxy[cluster][ExtrinsicAdjustment.INDX_DD0 + mode - 1] + fx[indx]);
//							dbg_img[5][pix] =  measured_dsxy[cluster][ExtrinsicAdjustment.INDX_DD0 + mode - 1] + fx[indx];
							dbg_img[4][pix] =  weights[indx]*clusters*(measured_dsxy[cluster][indx_dd0 + mode - 1] + fx[indx]);
							dbg_img[5][pix] =  measured_dsxy[cluster][indx_dd0 + mode - 1] + fx[indx];
						} else {
							dbg_img[3][pix] =  Double.NaN;
							dbg_img[4][pix] =  Double.NaN;
							dbg_img[5][pix] =  Double.NaN;
						}
					}
				} else {
					if (pix >= dbg_img[0].length) {
						System.out.println("pix="+pix+" >= "+dbg_img[0].length);
						System.out.println("pix="+pix+" >= "+dbg_img[0].length);
						System.out.println("pix="+pix+" >= "+dbg_img[0].length);
					}
					dbg_img[0][pix] =  Double.NaN; //Index 16240 out of bounds for length 16240
					dbg_img[1][pix] =  Double.NaN;
					dbg_img[2][pix] =  Double.NaN;
					dbg_img[3][pix] =  Double.NaN;
					dbg_img[4][pix] =  Double.NaN;
					dbg_img[5][pix] =  Double.NaN;
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
			CorrVector corr_vector,
			String title) {
		int gap = 10; //pix
		int clusters = clustersX * clustersY;
		String [] titles = {"meas", "correction", "diff"};
		int rows = getRowsCols()[0];
		int cols = getRowsCols()[1];
		int width  = cols * (clustersX + gap) - gap;
		int height = rows * (clustersY + gap) - gap;
		double [][] xy = getXYNondistorted(corr_vector, false);

		double [][] dbg_img = new double [3][width*height];
		double [][] moving_objects = new double [3][clusters];
		for (int mode = 0; mode < 2 * num_sensors + 1; mode++) {
			int x0 = (mode % cols) * (clustersX + gap);
			int y0 = (mode / cols) * (clustersY + gap);
			for (int cluster = 0; cluster < clusters;  cluster++) {
				int x = x0 + (cluster % clustersX);
				int y = y0 + (cluster / clustersX);
				int pix = x + y * width;
				if (measured_dsxy[cluster] != null){
					if (mode < (2 * num_sensors)) {
						dbg_img[0][pix] =  measured_dsxy[cluster][ExtrinsicAdjustment.INDX_X0+mode];
						dbg_img[1][pix] =  xy[cluster][mode] - x0y0[cluster][mode];
						dbg_img[2][pix] =  measured_dsxy[cluster][ExtrinsicAdjustment.INDX_X0+mode] - (xy[cluster][mode] - x0y0[cluster][mode]);
						moving_objects[0][cluster] += dbg_img[0][pix]*dbg_img[0][pix];
						moving_objects[1][cluster] += dbg_img[1][pix]*dbg_img[0][pix];
						moving_objects[2][cluster] += dbg_img[2][pix]*dbg_img[0][pix];
					} else {
						dbg_img[0][pix] = Math.sqrt(moving_objects[0][cluster]/(2 * num_sensors));
						dbg_img[1][pix] = Math.sqrt(moving_objects[0][cluster]/(2 * num_sensors));
						dbg_img[2][pix] = Math.sqrt(moving_objects[0][cluster]/(2 * num_sensors));
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
		System.out.println("dbgXY() DONE");
	}





	private double [][] getJacobianTransposed(
			CorrVector corr_vector)
	{
		if (dbg_delta > 0.0) {
			return getJacobianTransposed(corr_vector, dbg_delta); // running LMA with delta
		}
		int clusters = clustersX * clustersY;
		int num_pars = getNumPars();
		double [][] jt = new double  [num_pars][clusters * points_sample ];
		Matrix []   corr_rots =  corr_vector.getRotMatrices(); // get array of per-sensor rotation matrices
		Matrix [][] deriv_rots = corr_vector.getRotDeriveMatrices();
		double [] imu = corr_vector.getIMU(); // i)
		for (int cluster = 0; cluster < clusters;  cluster++) if (measured_dsxy[cluster] != null){
//			Mismatch mm = mismatch_list.get(indx);
//			double [] pXY = mm.getPXY();
			// will calculate 9 rows (disparity, dd0, dd1,cdd2, dd3, nd0, nd1, nd2, nd3}, columns - parameters
			// Now 2*num_sensors+1
			double [][] deriv = geometryCorrection.getPortsDDNDDerivativesNew( // USED in lwir
					geometryCorrection,     // GeometryCorrection gc_main,
					use_rig_offsets,        // boolean     use_rig_offsets,
					corr_rots,              // Matrix []   rots,
					deriv_rots,             // Matrix [][] deriv_rots,
					pY_offset[cluster],     // double []   py_offset,  // array of per-port average pY offset from the center (to correct ERS) or null (for no ERS)
					imu,                    // double []   imu,
					world_xyz[cluster],                  // double []   xyz, // world XYZ for ERS correction
					measured_dsxy[cluster][ExtrinsicAdjustment.INDX_PX + 0], // double      px,
					measured_dsxy[cluster][ExtrinsicAdjustment.INDX_PX + 1], // double      py,
					measured_dsxy[cluster][ExtrinsicAdjustment.INDX_TARGET]); // double      disparity);
///			int dbg_index = cluster; // dbg_index (pXY, dbg_decimate);
			// convert to symmetrical coordinates
			// derivatives of each port coordinates (in pixels) for each of selected symmetric all parameters (sym0 is convergence for disparity)
			 double [][] jt_partial = corr_vector.getJtPartial( // extract common?
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

	private double [] getWJtDiagSqrt(
			double [][] jt){
		int num_pars = jt.length;
		int nup_points = jt[0].length;
		double [] diag = new double [num_pars];
		for (int i = 0; i < num_pars; i++) {
				double d = 0.0;
				for (int k = 0; k < nup_points; k++) {
					d += this.weights[k]*jt[i][k]*jt[i][k];
				}
				diag[i] = Math.sqrt(d);
		}
		return diag;
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
		boolean [] rslt = {false,false};
		if (this.last_rms == null) { //first time, need to calculate all (vector is valid)
			this.last_jt =  getJacobianTransposed(corr_vector); // new double [num_pars][num_points];
			this.last_ymfx = getFx(corr_vector);
			if (debug_level > -2) { // temporary
				dbgYminusFxWeight(
						this.last_ymfx,
						this.weights,
						"Initial_y-fX_after_moving_objects");
			}

			if (last_ymfx == null) {
				return null; // need to re-init/restart LMA
			}
			this.last_rms = getWYmFxRms(
					this.last_ymfx); // modifies this.last_ymfx (weights and subtracts fx)
			this.initial_rms = this.last_rms.clone();
			this.good_or_bad_rms = this.last_rms.clone();
			if (debug_level > -1) {
				showDerivatives(0);
				showDerivatives(1);
				showDerivatives(2);
				showDerivatives(3);
			}
			// TODO: Restore/implement
			if (debug_level > 3) { // compares true jacobians and calculated with delta (same for all parameters)
				 dbgJacobians(
							corr_vector, // CorrVector corr_vector,
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
			double [] diag_sqrt = getWJtDiagSqrt(this.last_jt);
			System.out.print("diag_sqrt={");
			for (int i = 0; i < diag_sqrt.length; i++) {
				System.out.print(diag_sqrt[i]);
				if (i < (diag_sqrt.length-1)) System.out.print(", ");
			}
			System.out.println("}");
			System.out.print("diag_inverted_sqrt={");
			for (int i = 0; i < diag_sqrt.length; i++) {
				System.out.print(1.0/diag_sqrt[i]);
				if (i < (diag_sqrt.length-1)) System.out.print(", ");
			}
			System.out.println("}");
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
		CorrVector corr_delta = geometryCorrection.getCorrVector(delta, par_mask);

		CorrVector new_vector = this.corr_vector.clone();
		double scale = 1.0;
		new_vector.incrementVector(corr_delta, scale); // ok = false if there are nay NaN-s

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
				System.out.print("delta: "+corr_delta.toString()+"\n");
				System.out.print("New vector: "+new_vector.toString()+"\n");
				System.out.println();
			}
		} else { // worsened
			rslt[0] = false;
			rslt[1] = false; // do not know, caller will decide
			// restore state
			this.last_jt =  getJacobianTransposed(corr_vector); // new double [num_pars][num_points];
			this.last_ymfx = getFx(corr_vector);

			if (last_ymfx == null) {
				return null; // need to re-init/restart LMA
			}
			this.last_rms = getWYmFxRms(this.last_ymfx); // modifies this.last_ymfx
			if (debug_level > 2) {
				 dbgJacobians(
							corr_vector, // CorrVector corr_vector,
							1E-5, // double delta,
							true); //boolean graphic)
			}
		}
		return rslt;
	}
	
	public void showDerivatives(int typ4) {
		// typ == 0 -> DDND, 1 - XT
		int typ = typ4 % 2;
		boolean use_sym = typ4 > 1;
		int gap = 10;
		int clusters = clustersX * clustersY;
		//		int num_pars = getNumPars(); // only used
		int num_pars = corr_vector.toArray().length; // all parameters
		String [] titles = new String [num_pars]; //ea.getSymNames();
		for (int i = 0; i < num_pars; i++) {
			titles[i] = "S"+i;
		}
		Matrix []   corr_rots =  corr_vector.getRotMatrices(); // get array of per-sensor rotation matrices
		Matrix [][] deriv_rots = corr_vector.getRotDeriveMatrices();
		Matrix from_sym = corr_vector.symmVectorsSet.from_sym;
		double [] imu = corr_vector.getIMU(); //
		int rows = getRowsCols()[0];
		int cols = getRowsCols()[1];
		int width  = cols * (clustersX + gap) - gap;
		int height = rows * (clustersY + gap) - gap;
		double [][] dbg_img = new double [num_pars][width*height];
		//		double [][][] derivs = new double [points_sample][][];
		double [][][] derivs = new double [clusters][][];
		for (int cluster = 0; cluster < clusters;  cluster++) if (measured_dsxy[cluster] != null){
			//			Mismatch mm = mismatch_list.get(indx);
			//			double [] pXY = mm.getPXY();
			// will calculate 9 rows (disparity, dd0, dd1,cdd2, dd3, nd0, nd1, nd2, nd3}, columns - parameters
			// Now 2*num_sensors+1
			if (typ == 0) {
				derivs[cluster] = geometryCorrection.getPortsDDNDDerivativesNew( // USED in lwir
						this.geometryCorrection,     // GeometryCorrection gc_main,
						this.use_rig_offsets,        // boolean     use_rig_offsets,
						corr_rots,              // Matrix []   rots,
						deriv_rots,             // Matrix [][] deriv_rots,
						this.pY_offset[cluster],     // double []   py_offset,  // array of per-port average pY offset from the center (to correct ERS) or null (for no ERS)
						imu,                    // double []   imu,
						this.world_xyz[cluster],                  // double []   xyz, // world XYZ for ERS correction
						this.measured_dsxy[cluster][ExtrinsicAdjustment.INDX_PX + 0], // double      px,
						this.measured_dsxy[cluster][ExtrinsicAdjustment.INDX_PX + 1], // double      py,
						this.measured_dsxy[cluster][ExtrinsicAdjustment.INDX_TARGET]); // double      disparity);
			} else {
				derivs[cluster] = geometryCorrection.getPortsXYDerivativesNew( // USED in lwir
						this.geometryCorrection,     // GeometryCorrection gc_main,
						this.use_rig_offsets,        // boolean     use_rig_offsets,
						corr_rots,              // Matrix []   rots,
						deriv_rots,             // Matrix [][] deriv_rots,
						this.pY_offset[cluster],     // double []   py_offset,  // array of per-port average pY offset from the center (to correct ERS) or null (for no ERS)
						imu,                    // double []   imu,
						this.world_xyz[cluster],                  // double []   xyz, // world XYZ for ERS correction
						this.measured_dsxy[cluster][ExtrinsicAdjustment.INDX_PX + 0], // double      px,
						this.measured_dsxy[cluster][ExtrinsicAdjustment.INDX_PX + 1], // double      py,
						this.measured_dsxy[cluster][ExtrinsicAdjustment.INDX_TARGET]); // double      disparity);
				
			}
			// optionally convert to use sym vectors
			if (use_sym) {
				Matrix derivs_tarz = new Matrix(derivs[cluster]);
				Matrix derivs_sym = derivs_tarz.times(from_sym);
				derivs[cluster] = derivs_sym.getArray();
			}
		}
		for (int par = 0; par < num_pars; par++) {
			for (int mode = 0; mode < points_sample; mode++) {
				int x0 = (mode % cols) * (clustersX + gap);
				int y0 = (mode / cols) * (clustersY + gap);
				for (int cluster = 0; cluster < clusters;  cluster++) {
					int x = x0 + (cluster % clustersX);
					int y = y0 + (cluster / clustersX);
					int pix = x + y * width;
//					int indx = (mode == 0) ? INDX_DIFF : (indx_dd0 + mode - 1);
					if (typ == 0) {
						double d = (derivs[cluster] == null)? Double.NaN : derivs[cluster][mode][par];
						if (mode == 0) {
							dbg_img[par][pix] = -d;
						} else {
							dbg_img[par][pix] = d;
						}
					} else {
						if (mode == 0) {
							dbg_img[par][pix] = Double.NaN; // skip
						} else {
							dbg_img[par][pix] = (derivs[cluster] == null)? Double.NaN : derivs[cluster][mode-1][par];
						}
					}
				}
			}
		}
		String title = ((typ == 0)?"dDND_dpar":"dXY_dpar") + (use_sym? "_sym":"_tarz");
		(new ShowDoubleFloatArrays()).showArrays(
				dbg_img,
				width,
				height,
				true,
				title, // +(update_disparity?"U":""),
				(use_sym? titles: corr_vector.getCorrNames())); // titles);
	}
}
