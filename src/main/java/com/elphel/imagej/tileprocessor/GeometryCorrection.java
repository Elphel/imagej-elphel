package com.elphel.imagej.tileprocessor;
import java.util.ArrayList;
import java.util.Properties;

import com.elphel.imagej.common.GenericJTabbedDialog;
import com.elphel.imagej.common.ShowDoubleFloatArrays;

import Jama.Matrix;
import ij.IJ;

/**
 **
 ** GeometryCorrection - geometry correction for multiple sensors sharing the same
 ** lens radial distortion model
 **
 ** Copyright (C) 2017 Elphel, Inc.
 **
 ** -----------------------------------------------------------------------------**
 **
 **  GeometryCorrection.java is free software: you can redistribute it and/or modify
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

public class GeometryCorrection {
	static final double FOCAL_LENGTH = 4.5; // nominal focal length - used as default and to convert editable parameters to pixels
	static final double DISTORTION_RADIUS = 2.8512; // nominal distortion radius - half width of the sensor
	static final double PIXEL_SIZE = 2.2; //um
	static final String[] RIG_PAR_NAMES = {"azimuth", "tilt", "roll", "zoom", "angle", "baseline"};
	public static String RIG_PREFIX =     "rig-";
	static double SCENE_UNITS_SCALE = 0.001;
	static String SCENE_UNITS_NAME = "m";
	static final String [] CORR_NAMES = {"tilt0","tilt1","tilt2","azimuth0","azimuth1","azimuth2","roll0","roll1","roll2","roll3","zoom0","zoom1","zoom2"};

	public int    debugLevel = 0;
	public int    pixelCorrectionWidth=2592;   // virtual camera center is at (pixelCorrectionWidth/2, pixelCorrectionHeight/2)
	public int    pixelCorrectionHeight=1936;
	public double focalLength=FOCAL_LENGTH;
	public double pixelSize=  PIXEL_SIZE; //um
	public double distortionRadius=  DISTORTION_RADIUS; // mm - half width of the sensor
	public double distortionA8=0.0; //r^8 (normalized to focal length or to sensor half width?)
	public double distortionA7=0.0; //r^7 (normalized to focal length or to sensor half width?)
	public double distortionA6=0.0; //r^6 (normalized to focal length or to sensor half width?)
	public double distortionA5=0.0; //r^5 (normalized to focal length or to sensor half width?)
	public double distortionA=0.0; // r^4 (normalized to focal length or to sensor half width?)
	public double distortionB=0.0; // r^3
	public double distortionC=0.0; // r^2

	// parameters, common for all sensors
	public double    elevation = 0.0; // degrees, up - positive;
	public double    heading  =  0.0;  // degrees, CW (from top) - positive

	public int       numSensors = 4;
	private double [] forward =   null;
	private double [] right =     null;
	private double [] height =    null;
	private double [] roll  =     null;  // degrees, CW (to target) - positive
	public  double [][] pXY0 =    null;  // sensor center XY in pixels

	private double common_right;    // mm right, camera center
	private double common_forward;  // mm forward (to target), camera center
	private double common_height;   // mm up, camera center
	private double common_roll;     // degrees CW (to target) camera as a whole
	private double [][] XYZ_he;     // all cameras coordinates transformed to eliminate heading and elevation (rolls preserved)
	private double [][] XYZ_her = null; // XYZ of the lenses in a corrected CCS (adjusted for to elevation, heading,  common_roll)
	private double [][] rXY =     null; // XY pairs of the in a normal plane, relative to disparityRadius
	private double [][] rXY_ideal = {{-0.5, -0.5}, {0.5,-0.5}, {-0.5, 0.5}, {0.5,0.5}};

	public double cameraRadius=0; // average distance from the "mass center" of the sensors to the sensors
	public double disparityRadius=150.0; // distance between cameras to normalize disparity units to. sqrt(2)*disparityRadius for quad camera (~=150mm)?

	private double [] rByRDist=null;
	private double    stepR=0.001;
	private double    maxR=2.0; // calculate up to this*distortionRadius

	public  CorrVector extrinsic_corr;

	public RigOffset   rigOffset =    null;
	public int []             woi_tops; // used to calculate scanline timing

	public int [] getWOITops() {
		return woi_tops;
	}


	public int [] getSensorWH() {
		int [] wh = {this.pixelCorrectionWidth, this.pixelCorrectionHeight};
		return wh;
	}

	public GeometryCorrection(double [] extrinsic_corr)
	{
		this.extrinsic_corr = 	new CorrVector(extrinsic_corr);
	}

	public boolean isInitialized() {
		return roll != null;
	}

	public double [][] getRXY(boolean use_rig){
		return (use_rig && (rigOffset != null)) ? rigOffset.rXY_aux: rXY ;
	}

//	public double [] getAuxOffset(boolean use_rig){
//		double [] main_offset = {0.0,0.0};
//		return (use_rig && (rigOffset != null)) ? rigOffset.getAuxOffset(): main_offset ;
//	}

	public Matrix getRotMatrix(boolean use_rig){
		return (use_rig && (rigOffset != null)) ? rigOffset.getRotMatrix(): null ;
	}

	public double getDisparityRadius() {
		return disparityRadius;
	}
	public double getBaseline() {
		return (rigOffset==null)?Double.NaN:rigOffset.baseline;
	}

	public double [][] getAuxOffsetAndDerivatives(
			GeometryCorrection gc_main) {
		if (rigOffset == null) return null;
		return rigOffset.getAuxOffsetAndDerivatives(gc_main);
	}

	public Matrix getAuxRotMatrix() {
		if (rigOffset == null) return null;
		return rigOffset.getRotMatrix();
	}

	public Matrix [] getAuxRotDeriveMatrices() {
		if (rigOffset == null) return null;
		return rigOffset.getRotDeriveMatrices();
	}

	public RigOffset rigOffsetClone() {
		if (rigOffset == null) return null;
		return rigOffset.clone();
	}

	public void rigOffestSetParNorm(int index, double value) {
		rigOffset.setParNorm(index, value);
	}
	public void rigOffestSetParNorm(RigOffset ro, int index, double value) {
		ro. setParNorm(index, value);
	}

	public double rigOffestGetParNorm(int index) {
		return rigOffset.getParNorm(index);
	}
	public double rigOffestGetParNorm(RigOffset ro, int index) {
		return ro.getParNorm(index);
	}


	public double [] getRigCorrection(
			double             infinity_importance, // of all measurements
			double             dx_max, //  = 0.3;
			double             dx_pow, //  = 1.0;
			boolean            adjust_orientation,
			boolean            adjust_roll,
			boolean            adjust_zoom,
			boolean            adjust_angle,
			boolean            adjust_distance,
			boolean            adjust_forward, // not used
			double             scale_correction,
			double             infinity_disparity,
			ArrayList<Integer> tile_list,
			QuadCLT            qc_main,
			double []          strength,
			double []          diff_x, // used only with target_disparity == 0
			double []          diff_y,
			double []          target_disparity,
			int debugLevel) {
		if (rigOffset == null) return null;
		return rigOffset.getRigCorrection(
				infinity_importance, // of all measurements
				dx_max, // double             dx_max, //  = 0.3;
				dx_pow, // double             dx_pow, //  = 1.0;
				adjust_orientation, // boolean            adjust_orientation,
				adjust_roll,        // boolean            adjust_roll,
				adjust_zoom,        // boolean            adjust_zoom,
				adjust_angle,       // boolean            adjust_angle,
				adjust_distance,    // boolean            adjust_distance,
				adjust_forward,     // boolean            adjust_forward, // not used
				scale_correction,   // double             scale_correction,
				infinity_disparity, // double             infinity_disparity,
				tile_list,          // ArrayList<Integer> tile_list,
				qc_main,            // QuadCLT            qc_main,
				strength,           // double []          strength,
				diff_x,             // double []          diff_x, // used only with target_disparity == 0
				diff_y,             // double []          diff_y,
				target_disparity,   // double []          target_disparity,
				debugLevel);        // int debugLevel)
	}

	// correction of cameras mis-alignment
	public CorrVector getCorrVector(double [] vector){
		return new CorrVector(vector);
	}
	public CorrVector getCorrVector(
			double [] vector,
			boolean [] par_mask){
		return new CorrVector(vector, par_mask);
	}

	public CorrVector getCorrVector(){
		return extrinsic_corr;
	}

	public void setCorrVector(CorrVector vector){
		if (vector == null){
			vector = new CorrVector();
		}
		extrinsic_corr = vector;
	}

	public boolean [] getParMask(
//			boolean disparity_only,
//			boolean use_disparity,
			boolean use_disparity,
//			boolean use_other_extr,
			boolean use_aztilts,       // Adjust azimuths and tilts excluding disparity
			boolean use_diff_rolls,    // Adjust differential rolls (3 of 4 angles)


			boolean common_roll,
			boolean corr_focalLength,
	  		int     manual_par_sel)    // Manually select the parameter mask bit 0 - sym0, bit1 - sym1, ... (0 - use boolean flags, != 0 - ignore boolean flags)

	{
		return (new CorrVector()).getParMask(
				use_disparity, // disparity_only,
//				use_other_extr, // boolean use_other_extr,
				use_aztilts,       // Adjust azimuths and tilts excluding disparity
				use_diff_rolls,    // Adjust differential rolls (3 of 4 angles)
				common_roll,
				corr_focalLength,
		  		manual_par_sel);    // Manually select the parameter mask bit 0 - sym0, bit1 - sym1, ... (0 - use boolean flags, != 0 - ignore boolean flags)

	}
	/**
	 * Position of the auxiliary camera relative to the main one (uses main camera CS)
	 */
	public class RigOffset{
		static final int VECTOR_LENGTH =      6;
		static final int AUX_AZIMUTH_INDEX =  0;
		static final int AUX_TILT_INDEX =     1;
		static final int AUX_ROLL_INDEX =     2;
		static final int AUX_ZOOM_INDEX =     3;
		static final int AUX_ANGLE_INDEX =    4;
		static final int AUX_BASELINE_INDEX = 5;
		static final double ROT_AZ_SGN = -1.0; // sign of first sin for azimuth rotation
		static final double ROT_TL_SGN =  1.0; // sign of first sin for tilt rotation
		static final double ROT_RL_SGN =  1.0; // sign of first sin for roll rotation
		static final double BASELINE =  1256.0; // default baseline value
		public double baseline =   BASELINE; // mm, distance between camera centers
		public double aux_angle =     0.0; // radians, 0 - aux camera to the right of the main, pi/2 - up
		// consider aux_z small enough for now, will need for SfM
		public double aux_z     =     0.0; // mm auxiliary camera distance from the plane of the main one (positive towards the scene)
		public double aux_azimuth =   0.0; // radians, azimuth of the auxiliary camera (positive - looks to the right)
		public double aux_tilt =      0.0; // radians, tilt of the auxiliary camera (positive - looks up)
		public double aux_roll =      0.0; // radians, roll of the auxiliary camera (positive - looks clockwise)
		public double aux_zoom =      0.0; // relative global zoom of the aux camera relative to the main one, difference from 1.0
		public double  [][] rXY_aux =  null; // XY pairs of the in a normal plane, relative to disparityRadius
		public double  [] vector =     null;
		public boolean [] par_select = null;
		double [] w_vector = null;
		double [] y_vector = null;
		double [] xy_vector = null;
		double [] d_vector = null;
		int    [] tile_vector = null; // just for debugging
		double [] par_scales = null; // scale vector components to have similar derivatives values
		int    [] full_par_index;

		public RigOffset () {
			System.out.println("created RigOffset");
			par_scales = new double [VECTOR_LENGTH];
			par_scales[AUX_AZIMUTH_INDEX] =  1000.0*FOCAL_LENGTH/PIXEL_SIZE;
			par_scales[AUX_TILT_INDEX] =     1000.0*FOCAL_LENGTH/PIXEL_SIZE;
			par_scales[AUX_ROLL_INDEX] =     1000.0*DISTORTION_RADIUS/PIXEL_SIZE;
			par_scales[AUX_ZOOM_INDEX] =     1000.0*DISTORTION_RADIUS/PIXEL_SIZE;
			par_scales[AUX_ANGLE_INDEX] =    1.0; // 1000.0*BASELINE/pixelSize;
			par_scales[AUX_BASELINE_INDEX] = 1.0/DISTORTION_RADIUS; // pixels per disparity pixel

		}

		@Override
		public RigOffset clone() {
			RigOffset ro =      new RigOffset();
			ro.baseline =       this.baseline;
			ro.aux_angle =      this.aux_angle;
			ro.aux_z =          this.aux_z;
			ro.aux_azimuth =    this.aux_azimuth;
			ro.aux_tilt =       this.aux_tilt;
			ro.aux_roll =       this.aux_roll;
			ro.aux_zoom =       this.aux_zoom;
			ro.full_par_index = this.full_par_index.clone();
			ro.par_scales =     this.par_scales.clone();
			return ro;
		}

		public void setParNorm(int index, double value) {
			value /= par_scales[index];
			switch (index) {
			case AUX_AZIMUTH_INDEX:  aux_azimuth = value; break;
			case AUX_TILT_INDEX:     aux_tilt =    value; break;
			case AUX_ROLL_INDEX:     aux_roll =    value; break;
			case AUX_ZOOM_INDEX:     aux_zoom =    value; break;
			case AUX_ANGLE_INDEX:    aux_angle =   value; break;
			case AUX_BASELINE_INDEX: baseline =    value; break;
			}
		}
		public double getParNorm(int index) {
			switch (index) {
			case AUX_AZIMUTH_INDEX:  return aux_azimuth * par_scales[index];
			case AUX_TILT_INDEX:     return aux_tilt *    par_scales[index];
			case AUX_ROLL_INDEX:     return aux_roll *    par_scales[index];
			case AUX_ZOOM_INDEX:     return aux_zoom *    par_scales[index];
			case AUX_ANGLE_INDEX:    return aux_angle *   par_scales[index];
			case AUX_BASELINE_INDEX: return baseline *    par_scales[index];
			}
			return Double.NaN;
		}


		public void setVector(
				boolean adjust_orientation,
				boolean adjust_roll,
				boolean adjust_zoom,
				boolean adjust_angle,
				boolean adjust_distance,
				boolean adjust_forward) // not used
		{
			par_select = new boolean [VECTOR_LENGTH];
			par_select[AUX_AZIMUTH_INDEX] =  adjust_orientation;
			par_select[AUX_TILT_INDEX] =     adjust_orientation;
			par_select[AUX_ROLL_INDEX] =     adjust_roll;
			par_select[AUX_ZOOM_INDEX] =     adjust_zoom;
			par_select[AUX_ANGLE_INDEX] =    adjust_angle;
			par_select[AUX_BASELINE_INDEX] = adjust_distance;
			setVector();
		}
		public void setVector() {
			int num_pars = 0;
			for (int i = 0; i < par_select.length; i++) if (par_select[i]) num_pars++;
			vector = new double[num_pars];
			full_par_index = new int [num_pars];

			int par_index = 0;
			for (int i = 0; i < par_select.length; i++) if (par_select[i]) {
				switch(i) {
				case AUX_AZIMUTH_INDEX:  vector[par_index] = aux_azimuth * par_scales[i]; break;
				case AUX_TILT_INDEX:     vector[par_index] = aux_tilt    * par_scales[i]; break;
				case AUX_ROLL_INDEX:     vector[par_index] = aux_roll    * par_scales[i]; break;
				case AUX_ZOOM_INDEX:     vector[par_index] = aux_zoom    * par_scales[i]; break;
				case AUX_ANGLE_INDEX:    vector[par_index] = aux_angle   * par_scales[i]; break;
				case AUX_BASELINE_INDEX: vector[par_index] = baseline    * par_scales[i]; break;
				}
				full_par_index[par_index] = i;
				par_index++;
			}
			//full_par_index
		}
		public void commitVector(double [] v) {
			vector = v;
			int par_index = 0;
			for (int i = 0; i < par_select.length; i++) if (par_select[i]) {
				switch(i) {
				case AUX_AZIMUTH_INDEX:  aux_azimuth = vector[par_index]/par_scales[i]; break;
				case AUX_TILT_INDEX:     aux_tilt =    vector[par_index]/par_scales[i]; break;
				case AUX_ROLL_INDEX:     aux_roll =    vector[par_index]/par_scales[i]; break;
				case AUX_ZOOM_INDEX:     aux_zoom =    vector[par_index]/par_scales[i]; break;
				case AUX_ANGLE_INDEX:    aux_angle =   vector[par_index]/par_scales[i]; break;
				case AUX_BASELINE_INDEX: baseline =    vector[par_index]/par_scales[i]; break;
				}
				par_index++;
			}
			recalcRXY();
		}

		double [][] getJacobianTransposed(
				GeometryCorrection gc_main,
				int debugLevel){
			double [][] jt = new double[vector.length][xy_vector.length]; // npe
			double [][] pXYderiv =          new double [2][];
			Matrix      aux_rot =           getAuxRotMatrix();
			Matrix []   aux_rot_derivs =    getAuxRotDeriveMatrices();
			double [][] aux_offset_derivs = getAuxOffsetAndDerivatives(gc_main);

			for (int i = 0; i < xy_vector.length; i += 2) {
//				double [] xy =
				getRigAuxCoordinatesAndDerivatives(
						gc_main,           // GeometryCorrection gc_main,
						aux_rot,           // Matrix      aux_rot,
						aux_rot_derivs,    // Matrix []   aux_rot_derivs,
						aux_offset_derivs, // double [][] aux_offset_derivs,
						pXYderiv,          // double [][] pXYderiv, // if not null, should be double[6][]
						xy_vector[i + 0],  // double px,
						xy_vector[i + 1],  // double py,
						d_vector[i >> 1]); // double disparity);
				int npar = 0;
				for (int j = 0; j < VECTOR_LENGTH; j++) if (par_select[j]) {
					jt[npar][i + 0] = pXYderiv[0][j]/par_scales[full_par_index[npar]];
					jt[npar][i + 1] = pXYderiv[1][j]/par_scales[full_par_index[npar]];
					npar++;
				}
			}
			return jt;
		}

		// dbug method;
		double [][] getJacobianTransposed(
				double delta,
				GeometryCorrection gc_main,
				int debugLevel){
			double [][] jt = new double[vector.length][xy_vector.length];

			Matrix []   aux_rot_p =             new Matrix[RigOffset.VECTOR_LENGTH];
			double [][][] aux_offset_derivs_p = new double[RigOffset.VECTOR_LENGTH][][];
			Matrix []   aux_rot_m =             new Matrix[RigOffset.VECTOR_LENGTH];
			double [][][] aux_offset_derivs_m = new double[RigOffset.VECTOR_LENGTH][][];
			for (int npar = 0; npar < RigOffset.VECTOR_LENGTH; npar++ ) {
				RigOffset ro = rigOffsetClone();
				rigOffestSetParNorm(ro, npar, rigOffestGetParNorm(ro, npar) +delta);
				aux_offset_derivs_p[npar] = ro.getAuxOffsetAndDerivatives(gc_main);
				aux_rot_p[npar]=            ro.getRotMatrix();
				rigOffestSetParNorm(ro, npar, rigOffestGetParNorm(ro, npar) - 2 * delta);
				aux_offset_derivs_m[npar] = ro.getAuxOffsetAndDerivatives(gc_main);
				aux_rot_m[npar]=            ro.getRotMatrix();
			}

			for (int i = 0; i < xy_vector.length; i += 2) {
				double [][] pXYderiv =          new double [2][RigOffset.VECTOR_LENGTH];
				for (int npar = 0; npar < RigOffset.VECTOR_LENGTH; npar++ ) {
					double [] xy_p = getRigAuxCoordinatesAndDerivatives(
							gc_main,                   // GeometryCorrection gc_main,
							aux_rot_p[npar],           // Matrix      aux_rot,
							null,                      // Matrix []   aux_rot_derivs,
							aux_offset_derivs_p[npar], // double [][] aux_offset_derivs,
							null,                      // double [][] pXYderiv, // if not null, should be double[6][]
							xy_vector[i + 0],  // double px,
							xy_vector[i + 1],  // double py,
							d_vector[i >> 1]); // double disparity);
					double [] xy_m = getRigAuxCoordinatesAndDerivatives(
							gc_main,                   // GeometryCorrection gc_main,
							aux_rot_m[npar],           // Matrix      aux_rot,
							null,                      // Matrix []   aux_rot_derivs,
							aux_offset_derivs_m[npar], // double [][] aux_offset_derivs,
							null,                      // double [][] pXYderiv, // if not null, should be double[6][]
							xy_vector[i + 0],  // double px,
							xy_vector[i + 1],  // double py,
							d_vector[i >> 1]); // double disparity);
					pXYderiv[0][npar] = (xy_p[0]-xy_m[0])/ (2*delta);
					pXYderiv[1][npar] = (xy_p[1]-xy_m[1])/ (2*delta);
				}
				int npar = 0;
				for (int j = 0; j < VECTOR_LENGTH; j++) if (par_select[j]) {
					jt[npar][i + 0] = pXYderiv[0][j];
					jt[npar][i + 1] = pXYderiv[1][j];
					npar++;
				}
			}
			return jt;
		}





		double [][] getJTJWeighted(
				double [][] jt)
		{
			double [][] jtj = new double [jt.length][jt.length];
			for (int i = 0; i < jt.length; i++){
				for (int j = 0; j < i; j++){
					jtj[i][j] = jtj[j][i];
				}
				for (int j = i; j < jt.length; j++){
					for (int k = 0; k < jt[0].length; k++){
						jtj[i][j] += jt[i][k] * jt[j][k] * w_vector[k];
					}
				}
			}
			return jtj;
		}

		double [] getJTYWeighted(double [][] jt) {
			double [] jtyw = new double [jt.length];
			for (int i = 0; i < jt.length; i++){
				for (int k=0; k < jt[i].length; k++){
					jtyw[i] += jt[i][k] * y_vector[k] * w_vector[k];
				}
			}
			return jtyw;
		}

		public double [] getRigCorrection(
				double             infinity_importance, // of all measurements
				double             dx_max, //  = 0.3;
				double             dx_pow, //  = 1.0;
				boolean            adjust_orientation,
				boolean            adjust_roll,
				boolean            adjust_zoom,
				boolean            adjust_angle,
				boolean            adjust_distance,
				boolean            adjust_forward, // not used
				double             scale_correction,
				double             infinity_disparity,
				ArrayList<Integer> tile_list,
				QuadCLT            qc_main,
				double []          strength,
				double []          diff_x, // used only with target_disparity == 0
				double []          diff_y,
				double []          target_disparity,
				int debugLevel) {

			setVector(adjust_orientation,
					adjust_roll,
					adjust_zoom,
					adjust_angle,
					adjust_distance,
					adjust_forward); // not used

			double rms = setupYW(
					infinity_importance, // of all measurements
					dx_max, // double             dx_max, //  = 0.3;
					dx_pow, // double             dx_pow, //  = 1.0;
					infinity_disparity, // double             infinity_disparity,
					tile_list,
					qc_main,
					strength,
					diff_x, // used only with target_disparity == 0
					diff_y,
					target_disparity);
			if (Double.isNaN(rms)) {
				System.out.println("rms= NaN");
			}

			if (debugLevel>-4) {
				System.out.println("getRigCorrection(): Current RMS = "+rms+ "(debugLevel= "+debugLevel+")");
			};
			double [][] jt = getJacobianTransposed( // npe
					qc_main.geometryCorrection, // GeometryCorrection gc_main,
					debugLevel);                // int debugLevel)
			if (debugLevel > 2) {
				double [][] jt_delta = getJacobianTransposed(
						1.0E-6,
						qc_main.geometryCorrection, // GeometryCorrection gc_main,
						debugLevel);                // int debugLevel)
				int pars = jt_delta.length;
				int tilesX = qc_main.tp.getTilesX();
				int tilesY = qc_main.tp.getTilesX();
				String [] titles = new String[6 * pars];
				double [][] dbg_data = new double [6 * pars][tilesY * tilesX];
				for (int i = 0; i < dbg_data.length; i++) for (int j = 0; j < dbg_data[i].length; j++) dbg_data[i][j]=Double.NaN;
				int npar = 0;
				for (int indx = 0; indx < VECTOR_LENGTH; indx++) if (par_select[indx]) {
					titles[npar + 0*pars] = "dx/d_"+RIG_PAR_NAMES[npar];
					titles[npar + 1*pars] = "dy/d_"+RIG_PAR_NAMES[npar];
					titles[npar + 2*pars] = "delta_x/d_"+RIG_PAR_NAMES[npar];
					titles[npar + 3*pars] = "delta_y/d_"+RIG_PAR_NAMES[npar];
					titles[npar + 4*pars] = "diff_x/d_"+RIG_PAR_NAMES[npar];
					titles[npar + 5*pars] = "delt_ay/d_"+RIG_PAR_NAMES[npar];
					for (int i = 0; i < jt[npar].length; i+=2) {
						int nTile = tile_vector[i >> 1];
						dbg_data[npar + 0*pars][nTile] = jt[npar][i + 0];
						dbg_data[npar + 1*pars][nTile] = jt[npar][i + 1];
						dbg_data[npar + 2*pars][nTile] = jt_delta[npar][i + 0];
						dbg_data[npar + 3*pars][nTile] = jt_delta[npar][i + 1];
						dbg_data[npar + 4*pars][nTile] = jt[npar][i + 0] - jt_delta[npar][i + 0];
						dbg_data[npar + 5*pars][nTile] = jt[npar][i + 1] - jt_delta[npar][i + 1];
					}
					npar++;
				}
				(new ShowDoubleFloatArrays()).showArrays(
						dbg_data,
						tilesX,
						tilesY,
						true,
						"Rig Jacobians",
						titles);
			}

			Matrix jtj = new Matrix(getJTJWeighted(jt));
			Matrix jty = new Matrix( getJTYWeighted(jt),jt.length);
			if (debugLevel>-1) {
				System.out.println("getRigCorrection(): jtj= ");
				jtj.print(18, 6);
				System.out.println("getRigCorrection(): jty= ");
				jty.print(18, 6);
			}
			Matrix jtj_inv = jtj.inverse();
			if (debugLevel>-1) {
				System.out.println("getRigCorrection(): jtj_inv= ");
				jtj_inv.print(18, 6);
			}
			Matrix mrslt = jtj_inv.times(jty);
			if (debugLevel>-1) {
				System.out.println("getRigCorrection(): mrslt= ");
				mrslt.print(18, 6);
			}
			double []  drslt = mrslt.getColumnPackedCopy();
			// wrong sign?
			for (int i = 0; i < drslt.length; i++){
				drslt[i] *= -scale_correction;
			}

			for (int i = 0; i < drslt.length; i++){
				vector[i] += drslt[i];
			}
			if (debugLevel>-1) {
				System.out.println("getRigCorrection(): vector = ");
				for (int i = 0; i < vector.length; i++) {
					System.out.println(i+": "+vector[i]);
				}
			}
			commitVector(vector);
			return vector;
		}

		public double setupYW(
				double             infinity_importance, // of all measurements
				double             dx_max, //  = 0.3;
				double             dx_pow, //  = 1.0;
				double             infinity_disparity,
				ArrayList<Integer> tile_list,
				QuadCLT            qc,
				double []          strength,
				double []          diff_x, // used only with target_disparity == 0
				double []          diff_y,
				double []          target_disparity) {
			y_vector =  new double[2*tile_list.size()];
			w_vector =  new double[2*tile_list.size()];
			xy_vector = new double[2*tile_list.size()];
			d_vector =  new double[tile_list.size()];
			tile_vector=  new int[tile_list.size()];
			boolean [] is_inf = new boolean[tile_list.size()];
			double sumw_inf = 0.0,sumw_near=0.0;
			double sum2_inf = 0.0,sum2_near = 0.0;
			int tilesX = qc.tp.getTilesX();
			double tileSize = qc.tp.getTileSize();
//			double dx_max = 0.2;
//			double dx_pow = 1.0;

			for (int i = 0; i < tile_list.size(); i++) {
				int nTile = tile_list.get(i);
				int tileY = nTile / tilesX;
				int tileX = nTile % tilesX;
				double w = strength[nTile];
				if (Double.isNaN(diff_y[nTile])) {
					System.out.println("==== setupYW(): diff_y["+nTile+"]= NaN (can happen if the same tile is now measured differently)");
					w = 0.0;
					diff_y[nTile] = 0.0;
				}
				double inf_w_corr = 1.0;
				if ((dx_max > 0) && (dx_pow > 0)){
					inf_w_corr = (dx_max - diff_x[nTile])/dx_max;
					if (inf_w_corr < 0.0) inf_w_corr = 0.0; // favor negative (more infinity)
					if (dx_pow != 1.0) {
						inf_w_corr = Math.pow(inf_w_corr,dx_pow);
					}
				}
				if (target_disparity[nTile] == infinity_disparity) { // only for infinity tiles
					w *= inf_w_corr;
					if (Double.isNaN(diff_x[nTile])) {
						System.out.println("==== setupYW(): diff_x["+nTile+"]= NaN (can happen if the same tile is now measured differently)");
						w = 0.0;
						diff_x[nTile] = 0.0;
					}
					y_vector[2*i + 0] = diff_x[nTile];
					w_vector[2*i + 0] = w;
					y_vector[2*i + 1] = diff_y[nTile];
					w_vector[2*i + 1] = w;
					sumw_inf += 2*w;
					sum2_inf += w*(diff_x[nTile]*diff_x[nTile]+diff_y[nTile]*diff_y[nTile]);
					is_inf[i] = true;
				} else {
					y_vector[2*i + 1] = diff_y[nTile];
					w_vector[2*i + 1] = w;
					sumw_near +=        w;
					sum2_near += w*diff_y[nTile]*diff_y[nTile];
				}
				xy_vector[2*i + 0] = (tileX + 0.5) * tileSize;
				xy_vector[2*i + 1] = (tileY + 0.5) * tileSize;
				d_vector[i] = target_disparity[nTile];
				tile_vector[i] = nTile;
			}
			if (infinity_importance > 1.0) infinity_importance = 1.0;
			else if (infinity_importance < 0.0) infinity_importance =0.0;
			double k_inf = 0.0, k_near = 0.0;
			if ((sumw_inf  > 0.0)  && (sumw_near > 0.0)){
				k_inf =  infinity_importance/sumw_inf;
				k_near =  (1.0 - infinity_importance)/sumw_near;
			} else if (sumw_inf  > 0.0){
				infinity_importance = 1.0;
				k_inf =  infinity_importance/sumw_inf;
			} else if (sumw_near > 0.0) {
				infinity_importance = 0.0;
				k_near = (1.0 - infinity_importance)/sumw_near;
			}
			System.out.println("setupYW(): k_inf="+k_inf+" k_near="+k_near+" (sum2_inf="+sum2_inf+", sum2_near="+sum2_near+")");

			double sum2 = k_inf*sum2_inf+k_near*sum2_near;

			for (int i = 0; i < is_inf.length; i++) {
				if (is_inf[i]) {
					w_vector[2 * i + 0] *= k_inf;
					w_vector[2 * i + 1] *= k_inf;
				} else {
					w_vector[2 * i + 1] *= k_near;
				}
			}
			return Math.sqrt(sum2); // RMS
		}

		public void recalcRXY() {
			if (rXY != null) {
				//			rXY_aux = rXY; // FIXME: put real stuff !!!
				double xc_pix = baseline * Math.cos(aux_angle)/getDisparityRadius();
				double yc_pix = baseline * Math.sin(aux_angle)/getDisparityRadius();
				rXY_aux = new double [rXY.length][2];
				double ssr = 1.0;
				for (int i = 0; i <rXY.length;i++) {
					rXY_aux[i][0] = xc_pix +     Math.cos(aux_roll)*rXY[i][0] + ssr*Math.sin(aux_roll)*rXY[i][1];
					rXY_aux[i][1] = yc_pix - ssr*Math.sin(aux_roll)*rXY[i][0] +     Math.cos(aux_roll)*rXY[i][1];
				}
				if (debugLevel > 0) {
					System.out.println("Auxiliary camera offsets per 1 nominal disparity pixel");
					for (int i = 0; i <rXY_aux.length;i++) {
						System.out.println(String.format("Camera %1d x = %8f y = %8f",i,rXY_aux[i][0],rXY_aux[i][1]));
					}
				}
			}
		}
		/**
		 * Get pixel offset of the idealized aux camera center in main camera coordinates per
		 * 1 pixel of the main camera disparity
		 * @param gc_main Instance of the main camera GeometryCorrection class
		 * @return {{xc, yc},{dxc/dAngle,dyc/dAngle},{dxc/dBaseline,dyc/dBaseline}}
		 */
		public double [][] getAuxOffsetAndDerivatives(
				GeometryCorrection gc_main) {
			double blp = baseline /gc_main.getDisparityRadius();
			double xc_pix = blp * Math.cos(aux_angle);
			double yc_pix = blp * Math.sin(aux_angle);
			double dxc_dangle = -yc_pix;
			double dyc_dangle =  xc_pix;
			double dxc_baseline = xc_pix/baseline;
			double dyc_baseline = xc_pix/baseline;

			double [][] rslt = {
					{ xc_pix,      yc_pix},
					{dxc_dangle,   dyc_dangle},
					{dxc_baseline, dyc_baseline}};
			return rslt;
		}

		public Matrix getRotMatrix()
		{
			//			Matrix [] rots = new Matrix [4];
			//			double [] azimuths = getAzimuths();
			//			double [] tilts =    getTilts();
			//			double [] rolls =    getFullRolls();
			//			double [] zooms =    getZooms();
			double ca = Math.cos(aux_azimuth);
			double sa = Math.sin(aux_azimuth);
			double ct = Math.cos(aux_tilt);
			double st = Math.sin(aux_tilt);
			double zoom = (1.0 + aux_zoom);
			double cr = Math.cos(aux_roll) * zoom;
			double sr = Math.sin(aux_roll) * zoom;
			double [][] a_az = { // inverted - OK
					{ ca,               0.0,  sa * ROT_AZ_SGN },
					{ 0.0,              1.0,  0.0},
					{ -sa* ROT_AZ_SGN,  0.0,  ca}};

			double [][] a_t =  { // inverted - OK
					{ 1.0,  0.0, 0.0},
					{ 0.0,  ct,               st * ROT_TL_SGN},
					{ 0.0, -st * ROT_TL_SGN,  ct}};

			double [][] a_r =  { // inverted OK
					{ cr,                sr * ROT_RL_SGN,  0.0},
					{ -sr * ROT_RL_SGN,  cr,               0.0},
					{ 0.0,               0.0,              1.0}};

			Matrix rot  = (new Matrix(a_r).times(new Matrix(a_t).times(new Matrix(a_az))));
			return rot;
		}

		/**
		 * Get derivatives of the auxiliary camera rotation matrix, per axis (azimuth, tilt, roll, zoom)
		 * d/dx and d/dy should be normalized by z-component of the vector (not derivative)
		 * @return 2-d array array of derivatives matrices
		 */
//TODO: UPDATE to include scales
		public Matrix [] getRotDeriveMatrices()
		{
			Matrix [] rot_derivs = new Matrix [4]; // channel, azimuth-tilt-roll-zoom

				double ca = Math.cos(aux_azimuth);
				double sa = Math.sin(aux_azimuth);
				double ct = Math.cos(aux_tilt);
				double st = Math.sin(aux_tilt);
				double zoom = (1.0 + aux_zoom);
				double cr = Math.cos(aux_roll);
				double sr = Math.sin(aux_roll);
				double [][] a_az = { // inverted - OK
						{  ca,                     0.0,                     sa * ROT_AZ_SGN },
						{  0.0,                    1.0,                     0.0},
						{ -sa* ROT_AZ_SGN,         0.0,                     ca}};

				double [][] a_t =  { // inverted - OK
						{  1.0,                    0.0,                     0.0},
						{  0.0,                    ct,                      st * ROT_TL_SGN},
						{  0.0,                   -st * ROT_TL_SGN,         ct}};

				double [][] a_r =  { // inverted OK
						{  cr * zoom,              sr * zoom * ROT_RL_SGN,  0.0},
						{ -sr * zoom * ROT_RL_SGN, cr * zoom,               0.0},
						{  0.0,                    0.0,                     1.0}};

				double [][] a_daz = { // inverted - OK
						{ -sa,                     0.0,                     ca * ROT_AZ_SGN },
						{  0.0,                    0.0,                     0.0},
						{ -ca* ROT_AZ_SGN,         0.0,                    -sa}};

				double [][] a_dt =  { // inverted - OK
						{ 0.0,  0.0,               0.0},
						{ 0.0, -st,                ct * ROT_TL_SGN},
						{ 0.0, -ct * ROT_TL_SGN,  -st}};

				double [][] a_dr =  { // inverted OK
						{ -sr * zoom,              cr * zoom * ROT_RL_SGN,  0.0},
						{ -cr * zoom *ROT_RL_SGN, -sr * zoom,               0.0},
						{ 0.0,                     0.0,                     0.0}};

				double [][] a_dzoom =  { // inverted OK
						{ cr,                      sr * ROT_RL_SGN,         0.0},
						{ -sr * ROT_RL_SGN,        cr,                      0.0},
						{ 0.0,                     0.0,                     0.0}};


				// d/d_az
				rot_derivs[0] = (new Matrix(a_r ).times(new Matrix(a_t ).times(new Matrix(a_daz))));
				rot_derivs[1] = (new Matrix(a_r ).times(new Matrix(a_dt).times(new Matrix(a_az ))));
				rot_derivs[2] = (new Matrix(a_dr).times(new Matrix(a_t ).times(new Matrix(a_az ))));
				rot_derivs[3] = (new Matrix(a_dzoom).times(new Matrix(a_t ).times(new Matrix(a_az ))));
			return rot_derivs;
		}





		public void setProperties(String parent_prefix,Properties properties){
			String prefix = parent_prefix + RIG_PREFIX;
			properties.setProperty(prefix+"baseline",      this.baseline+"");
			properties.setProperty(prefix+"aux_angle",     this.aux_angle+"");
			properties.setProperty(prefix+"aux_z",         this.aux_z+"");
			properties.setProperty(prefix+"aux_azimuth",   this.aux_azimuth+"");
			properties.setProperty(prefix+"aux_tilt",      this.aux_tilt+"");
			properties.setProperty(prefix+"aux_roll",      this.aux_roll+"");
			properties.setProperty(prefix+"aux_zoom",      this.aux_zoom+"");
		}
		public boolean getProperties(String parent_prefix,Properties properties){
			String prefix = parent_prefix + RIG_PREFIX;
			boolean got_data = false;
			if (properties.getProperty(prefix+"baseline")!=null)    {this.baseline=Double.parseDouble(properties.getProperty(prefix+"baseline")); got_data=true;}
			if (properties.getProperty(prefix+"aux_angle")!=null)   {this.aux_angle=Double.parseDouble(properties.getProperty(prefix+"aux_angle"));got_data=true;}
			if (properties.getProperty(prefix+"aux_z")!=null)       {this.aux_z=Double.parseDouble(properties.getProperty(prefix+"aux_z"));got_data=true;}
			if (properties.getProperty(prefix+"aux_azimuth")!=null) {this.aux_azimuth=Double.parseDouble(properties.getProperty(prefix+"aux_azimuth"));got_data=true;}
			if (properties.getProperty(prefix+"aux_tilt")!=null)    {this.aux_tilt=Double.parseDouble(properties.getProperty(prefix+"aux_tilt"));got_data=true;}
			if (properties.getProperty(prefix+"aux_roll")!=null)    {this.aux_roll=Double.parseDouble(properties.getProperty(prefix+"aux_roll"));got_data=true;}
			if (properties.getProperty(prefix+"aux_zoom")!=null)    {this.aux_zoom=Double.parseDouble(properties.getProperty(prefix+"aux_zoom"));got_data=true;}
			recalcRXY();
			return got_data;
		}
		// 9:%8.5f° 10: %8.5f‰
		public boolean editOffsetsDegrees() {
  			GenericJTabbedDialog gd = new GenericJTabbedDialog("Set CLT parameters",800,900);
			gd.addNumericField("Baseline",                                                            this.baseline,  1,6,"mm",
					"Distance between quad camera centers");
			gd.addNumericField("Angle to the aux camera from the main",                               180.0/Math.PI*this.aux_angle,  4,9,"°",
					"Directly to the right - 0°, directly up - 90°, ...");
			gd.addNumericField("Auxilliary camera forward from the plane of the main one (not used)", this.aux_z,  3,6,"mm",
					"Distance from the plane perpendicualr to the main camera axis to the auxiliary camera (positive for aux moved forward)");
			gd.addNumericField("Auxilliary camera azimuth  (positive - to the right)",                180.0/Math.PI*this.aux_azimuth,  3,6,"°",
					"Relative to the main camera axis");
			gd.addNumericField("Auxilliary camera tilt (positive - looking up)",                      180.0/Math.PI*this.aux_tilt,  3,6,"°",
					"Relative to the main camera");
			gd.addNumericField("Auxilliary camera roll (positive - clockwise)",                       180.0/Math.PI*this.aux_roll,  3,6,"°",
					"Roll of a camera as a whole relative to the main camera");
			gd.addNumericField("Relative zoom",                                                       1000.0*this.aux_zoom,  3,6,"‰",
					"Zoom ratio minus 1.0 multiplied by 1000.0");
  			gd.showDialog();
			if (gd.wasCanceled()) return false;
			this.baseline=     gd.getNextNumber();
			this.aux_angle=    gd.getNextNumber() * Math.PI/180;
			this.aux_z=        gd.getNextNumber();
			this.aux_azimuth=  gd.getNextNumber() * Math.PI/180;
			this.aux_tilt=     gd.getNextNumber() * Math.PI/180;
			this.aux_roll=     gd.getNextNumber() * Math.PI/180;
			this.aux_zoom=     gd.getNextNumber()/1000.0;
			recalcRXY();
			return true;
		}
		public boolean editOffsetsPixels() {
  			GenericJTabbedDialog gd = new GenericJTabbedDialog("Set dual camera rig parameters (auxiliary camera relative to the main one)",800,300);
			gd.addNumericField("Baseline",                                                            this.baseline,  1,6,"mm",
					"Distance between quad camera centers");
			gd.addNumericField("Angle to the aux camera from the main",                               180.0/Math.PI*this.aux_angle,  4,9,"°",
					"Directly to the right - 0°, directly up - 90°, ...");
			gd.addNumericField("Auxilliary camera forward from the plane of the main one (not used)", this.aux_z,  3,6,"mm",
					"Distance from the plane perpendicualr to the main camera axis to the auxiliary camera (positive for aux moved forward)");
			gd.addNumericField("Auxilliary camera azimuth  (positive - image to the right)",                par_scales[AUX_AZIMUTH_INDEX] * this.aux_azimuth,  3,6,"pix",
					"Positive - converge more: aux (right) camera to the left, image - to the right");
			gd.addNumericField("Auxilliary camera tilt (positive - move aux image up)",                      par_scales[AUX_TILT_INDEX] * this.aux_tilt,  3,6,"pix",
					"Positive: aux looks down, image moves up");
			gd.addNumericField("Auxilliary camera roll (positive - clockwise)",                       par_scales[AUX_ROLL_INDEX] * this.aux_roll,  3,6,"pix",
					"Roll of a camera as a whole relative to the main camera, shift at the image half-width from the center");
			gd.addNumericField("Relative zoom - difference from 1.0 in parts parts per 1/1000",       par_scales[AUX_ZOOM_INDEX] * this.aux_zoom,  3,6,"pix",
					"Zoom ratio, shift at the image half-width from the center");
  			gd.showDialog();
			if (gd.wasCanceled()) return false;
			this.baseline=     gd.getNextNumber();
			this.aux_angle=    gd.getNextNumber() * Math.PI/180;
			this.aux_z=        gd.getNextNumber();
			this.aux_azimuth=  gd.getNextNumber()/par_scales[AUX_AZIMUTH_INDEX] ;
			this.aux_tilt=     gd.getNextNumber()/par_scales[AUX_TILT_INDEX];
			this.aux_roll=     gd.getNextNumber()/par_scales[AUX_ROLL_INDEX];
			this.aux_zoom=     gd.getNextNumber()/par_scales[AUX_ZOOM_INDEX];
			recalcRXY();
			return true;
		}
		public void showRigOffsets()
		{
			System.out.println("=== Inter-camera adjustments ===");
			System.out.println("                                                           Baseline "+ this.baseline +"mm");
			System.out.println("                              Angle to the aux camera from the main "+ (180.0/Math.PI*this.aux_angle)+"°");
			System.out.println("Auxilliary camera forward from the plane of the main one (not used) "+ this.aux_z +"mm");
			System.out.println("Auxilliary camera azimuth  (positive: camera - left, image - right) "+ (par_scales[AUX_AZIMUTH_INDEX] * this.aux_azimuth) + "pix");
			System.out.println("   Auxilliary camera tilt (positive: aux camera - down, image - up) "+ (par_scales[AUX_TILT_INDEX] * this.aux_tilt)+"pix");
			System.out.println("                      Auxilliary camera roll (positive - clockwise) "+ (par_scales[AUX_ROLL_INDEX] * this.aux_roll)+"pix");
			System.out.println("      Relative zoom - difference from 1.0 in parts parts per 1/1000 "+ (par_scales[AUX_ZOOM_INDEX] * this.aux_zoom) +"pix");
		}

		public double getRigOffsetParameter(int indx, boolean inPix)
		{
			if ((indx <0.0) || (indx >=par_scales.length)){
				return Double.NaN;
			}
			double k = inPix? par_scales[indx] : 1.0;
			switch (indx) {
			case AUX_AZIMUTH_INDEX:  return k * this.aux_azimuth;
			case AUX_TILT_INDEX:     return k * this.aux_tilt;
			case AUX_ROLL_INDEX:     return k * this.aux_roll;
			case AUX_ZOOM_INDEX:     return k * this.aux_zoom;
			case AUX_ANGLE_INDEX:    return k * this.aux_angle;
			case AUX_BASELINE_INDEX: return k * this.baseline; // if inPix - show rig baseline/main baseline
			}
			return Double.NaN;
		}

	}

	public boolean editRig() {
		if (this.rigOffset == null) {
			this.rigOffset = new RigOffset();
		}
		return this.rigOffset.editOffsetsPixels();
	}

	public void showRig() {
		this.rigOffset.showRigOffsets();
	}

	public double getRigOffsetParameter(int indx, boolean inPix) {
		return this.rigOffset.getRigOffsetParameter(indx,inPix);
	}

	public boolean setRigOffsetFromProperies(String parent_prefix,Properties properties) {
		RigOffset rigOffset = new RigOffset();
		boolean gotit = rigOffset.getProperties(parent_prefix, properties);
		if (gotit) {
			this.rigOffset = rigOffset;
		}
		return gotit;
	}


	public class CorrVector{
		static final int LENGTH =       13; //  10;
		static final int LENGTH_ANGLES =10;
		static final int TILT_INDEX =    0;
		static final int AZIMUTH_INDEX = 3;
		static final int ROLL_INDEX =    6;
		static final int ZOOM_INDEX =    10;
		static final double ROT_AZ_SGN = -1.0; // sign of first sin for azimuth rotation
		static final double ROT_TL_SGN =  1.0; // sign of first sin for tilt rotation
		static final double ROT_RL_SGN =  1.0; // sign of first sin for roll rotation
		double [] vector;

		public Matrix [] getRotMatrices(Matrix rigMatrix)
		{
			Matrix [] rots = getRotMatrices();
			if (rigMatrix != null) {
				for (int chn = 0; chn < rots.length; chn++) {
//					rots[chn] = rigMatrix.times(rots[chn]);
					rots[chn] = rots[chn].times(rigMatrix); // rots[chn] left matrix is rotation (closest to the camera), it should remain leftmost
				}
			}
			return rots;
		}

		// not yet used
		public Matrix [][] getRotDeriveMatrices(Matrix rigMatrix)
		{
			Matrix [][] derivs = getRotDeriveMatrices();
			if (rigMatrix != null) {
				for (int chn = 0; chn < derivs.length; chn++) {
					for (int deriv_index=0; deriv_index < derivs[chn].length; deriv_index++) {
//						derivs[chn][deriv_index] = rigMatrix.times(derivs[chn][deriv_index]);
						derivs[chn][deriv_index] = derivs[chn][deriv_index].times(rigMatrix); // left matrix is rotation (closest to the camera), it should remain leftmost
					}
				}
			}
			return derivs;
		}

		public Matrix [] getRotMatrices()
		{
			Matrix [] rots = new Matrix [4];
			double [] azimuths = getAzimuths();
			double [] tilts =    getTilts();
			double [] rolls =    getFullRolls();
			double [] zooms =    getZooms();
			for (int chn = 0; chn < rots.length; chn++) {
				double ca = Math.cos(azimuths[chn]);
				double sa = Math.sin(azimuths[chn]);
				double ct = Math.cos(tilts[chn]);
				double st = Math.sin(tilts[chn]);
				double zoom = (1.0 + zooms[chn]);
				double cr = Math.cos(rolls[chn]) * zoom;
				double sr = Math.sin(rolls[chn]) * zoom;
				double [][] a_az = { // inverted - OK
						{ ca,               0.0,  sa * ROT_AZ_SGN },
						{ 0.0,              1.0,  0.0},
						{ -sa* ROT_AZ_SGN,  0.0,  ca}};

				double [][] a_t =  { // inverted - OK
						{ 1.0,  0.0, 0.0},
						{ 0.0,  ct,               st * ROT_TL_SGN},
						{ 0.0, -st * ROT_TL_SGN,  ct}};

				double [][] a_r =  { // inverted OK
						{ cr,                sr * ROT_RL_SGN,  0.0},
						{ -sr * ROT_RL_SGN,  cr,               0.0},
						{ 0.0,               0.0,              1.0}};

				rots[chn] = (new Matrix(a_r).times(new Matrix(a_t).times(new Matrix(a_az))));
			}
			return rots;
		}
		/**
		 * Get derivatives of the rotation matrices, per sensor per axis (azimuth, tilt, roll, zoom)
		 * d/dx and d/dy should be normalized by z-component of the vector (not derivative)
		 * @return 2-d array array of derivatives matrices
		 */
//TODO: UPDATE to include scales
		public Matrix [][] getRotDeriveMatrices()
		{
			Matrix [][] rot_derivs = new Matrix [4][4]; // channel, azimuth-tilt-roll-zoom
			double [] azimuths = getAzimuths();
			double [] tilts =    getTilts();
			double [] rolls =    getFullRolls();
			double [] zooms =    getZooms();

			for (int chn = 0; chn < rot_derivs.length; chn++) {
				double ca = Math.cos(azimuths[chn]);
				double sa = Math.sin(azimuths[chn]);
				double ct = Math.cos(tilts[chn]);
				double st = Math.sin(tilts[chn]);
				double zoom = (1.0 + zooms[chn]);
				double cr = Math.cos(rolls[chn]);
				double sr = Math.sin(rolls[chn]);
				double [][] a_az = { // inverted - OK
						{  ca,                     0.0,                     sa * ROT_AZ_SGN },
						{  0.0,                    1.0,                     0.0},
						{ -sa* ROT_AZ_SGN,         0.0,                     ca}};

				double [][] a_t =  { // inverted - OK
						{  1.0,                    0.0,                     0.0},
						{  0.0,                    ct,                      st * ROT_TL_SGN},
						{  0.0,                   -st * ROT_TL_SGN,         ct}};

/*
				double [][] a_r =  { // inverted OK
						{  cr,                     sr * ROT_RL_SGN,         0.0},
						{ -sr * ROT_RL_SGN,        cr,                      0.0},
						{  0.0,                    0.0,                     1.0}};
 */
				double [][] a_r =  { // inverted OK
						{  cr * zoom,              sr * zoom * ROT_RL_SGN,  0.0},
						{ -sr * zoom * ROT_RL_SGN, cr * zoom,               0.0},
						{  0.0,                    0.0,                     1.0}};

				double [][] a_daz = { // inverted - OK
						{ -sa,                     0.0,                     ca * ROT_AZ_SGN },
						{  0.0,                    0.0,                     0.0},
						{ -ca* ROT_AZ_SGN,         0.0,                    -sa}};

				double [][] a_dt =  { // inverted - OK
						{ 0.0,  0.0,               0.0},
						{ 0.0, -st,                ct * ROT_TL_SGN},
						{ 0.0, -ct * ROT_TL_SGN,  -st}};

				double [][] a_dr =  { // inverted OK
						{ -sr * zoom,              cr * zoom * ROT_RL_SGN,  0.0},
						{ -cr * zoom *ROT_RL_SGN, -sr * zoom,               0.0},
						{ 0.0,                     0.0,                     0.0}};

				double [][] a_dzoom =  { // inverted OK
						{ cr,                      sr * ROT_RL_SGN,         0.0},
						{ -sr * ROT_RL_SGN,        cr,                      0.0},
						{ 0.0,                     0.0,                     0.0}};

				// d/d_az
				rot_derivs[chn][0] = (new Matrix(a_r ).times(new Matrix(a_t ).times(new Matrix(a_daz))));    // d/daz   - wrong, needs *zoom
				rot_derivs[chn][1] = (new Matrix(a_r ).times(new Matrix(a_dt).times(new Matrix(a_az ))));    // d/dt    - wrong, needs *zoom
				rot_derivs[chn][2] = (new Matrix(a_dr).times(new Matrix(a_t ).times(new Matrix(a_az ))));    // d/dr    - correct, has zoom
				rot_derivs[chn][3] = (new Matrix(a_dzoom).times(new Matrix(a_t ).times(new Matrix(a_az )))); // d/dzoom - correct
			}
			return rot_derivs;
		}




		public CorrVector ()
		{
			this.vector = new double[LENGTH];
		}

		public CorrVector (
				double [] sym_vector,
				boolean [] par_mask)
		{
			this.vector = toTarArray(sym_vector, par_mask);

		}

		public CorrVector (
				double tilt0,    double tilt1,    double tilt2,
				double azimuth0, double azimuth1, double azimuth2,
				double roll0,    double roll1,    double roll2, double roll3,
				double zoom0,   double zoom1,   double zoom2)
		{
			double [] v = {
					tilt0,    tilt1,    tilt2,
					azimuth0, azimuth1, azimuth2,
					roll0,    roll1,    roll2, roll3,
					zoom0,   zoom1,   zoom2};
			this.vector = v;
		}


		public CorrVector (double [] vector)
		{

			if (vector != null) {
				if (vector.length != LENGTH) {
					throw new IllegalArgumentException("vector.length = "+vector.length+" != "+LENGTH);
				}
				this.vector = vector;
			}
		}
		/**
		 * Set subcamera corrections from a single array
		 * @param tilt    for subcameras 0..2, radians, positive - up (subcamera 3 so sum == 0)
		 * @param azimuth for subcameras 0..2, radians, positive - right (subcamera 3 so sum == 0)
		 * @param roll    for subcameras 0..3, radians, positive - CW looking to the target
		 * @param zoom    for subcameras 0..2, difference from 1.0 . Positive - image is too small, needs to be zoomed in by (1.0 + scale)
		 */
		public CorrVector (double [] tilt, double [] azimuth, double [] roll, double [] zoom)
		{
			double [] vector = {
					tilt[0],    tilt[1],    tilt[2],
					azimuth[0], azimuth[1], azimuth[2],
					roll[0], roll[1], roll[2], roll[3],
					zoom[0], zoom[1], zoom[2]};
			this.vector = vector;
		}

		public CorrVector getCorrVector(double [] vector){
			return new CorrVector(vector);
		}



		public double [] toArray()
		{
			return vector;
		}
		public double [] getTilts()
		{
			double [] tilts =     {vector[0], vector[1], vector[2], - (vector[0] + vector[1] +vector[2])};
			return tilts;
		}
		public double getTilt(int indx)
		{
			if (indx == 3) return - (vector[0] + vector[1] +vector[2]);
			else           return vector[0 + indx];
		}
		public double [] getAzimuths()
		{
			double [] azimuths =  {vector[3], vector[4], vector[5], -(vector[3] + vector[4] + vector[5])};
			return azimuths;
		}
		public double getAzimuth(int indx)
		{
			if (indx == 3) return - (vector[3] + vector[4] +vector[5]);
			else           return vector[3 + indx];
		}
		public double [] getRolls()
		{
			double [] rolls =     {vector[6],vector[7],vector[8], vector[9]};
			return rolls;
		}

		public double getRoll(int indx)
		{
			return vector[6 + indx];
		}

		public double [] getZooms()
		{
			double [] zooms =     {vector[10], vector[11], vector[12], - (vector[10] + vector[11] +vector[12])};
			return zooms;
		}
		public double getZoom(int indx)
		{
			if (indx == 3) return - (vector[10] + vector[11] +vector[12]);
			else           return vector[10 + indx];
		}


		// Include factory calibration rolls
		public double [] getFullRolls()
		{
			double d2r= Math.PI/180.0;
			double [] rolls =     {
					vector[6] + d2r * roll[0],
					vector[7] + d2r * roll[1],
					vector[8] + d2r * roll[2],
					vector[9] + d2r * roll[3]};
			return rolls;
		}
		public double getFullRoll(int indx)
		{
			return vector[6 + indx] + roll[indx] * Math.PI/180.0;
		}
		/**
		 * Return parameter value for reports
		 * @param indx parameter index (use CorrVector.XXX static integers)
		 * @param inPix show result in pixels , false - in radians (even for zooms)
		 * @return parameter value
		 */
		public double getExtrinsicParameterValue(int indx, boolean inPix) {
			if (indx <0) return Double.NaN;
			if (indx <    ROLL_INDEX) return vector[indx]* (inPix? (1000.0*focalLength     /pixelSize): 1.0); // tilt and azimuth
			if (indx < LENGTH_ANGLES) return vector[indx]* (inPix? (1000.0*distortionRadius/pixelSize): 1.0); //  rolls
			if (indx <        LENGTH) return vector[indx]* (inPix? (1000.0*distortionRadius/pixelSize): 1.0); //  zooms
			return Double.NaN;
		}
		public double getExtrinsicSymParameterValue(int indx, boolean inPix) {
			double [] sym_vect = toSymArray(null);
			if (indx <0) return Double.NaN;
			if (indx <    ROLL_INDEX) return sym_vect[indx]* (inPix? (1000.0*focalLength     /pixelSize): 1.0); // tilt and azimuth
			if (indx < LENGTH_ANGLES) return sym_vect[indx]* (inPix? (1000.0*distortionRadius/pixelSize): 1.0); //  rolls
			if (indx <        LENGTH) return sym_vect[indx]* (inPix? (1000.0*distortionRadius/pixelSize): 1.0); //  zooms
			return Double.NaN;
		}

		@Override
		public String toString()
		{
			String s;
			double [] sym_vect = toSymArray(null);
			double [] v = new double [vector.length];
			double [] sv = new double [vector.length];
			for (int i = 0; i < ROLL_INDEX; i++){
				v[i] = vector[i]*1000.0*focalLength/pixelSize; // tilt and azimuth
				sv[i] = sym_vect[i]*1000.0*focalLength/pixelSize; // all 3 angles
			}
			for (int i = ROLL_INDEX; i < LENGTH_ANGLES; i++){
				v[i] = vector[i]*1000.0*distortionRadius/pixelSize; // rolls
				sv[i] = sym_vect[i]*1000.0*distortionRadius/pixelSize; // combined rolls
			}
			for (int i = LENGTH_ANGLES; i < LENGTH; i++){
				v[i] = vector[i]*1000.0*distortionRadius/pixelSize; // zooms
				sv[i] = sym_vect[i]*1000.0*distortionRadius/pixelSize; // zooms
			}

			s  = String.format("tilt    (up):    %8.5fpx %8.5fpx %8.5fpx %8.5fpx (shift of he image center)\n" , v[0], v[1], v[2], -(v[0] + v[1] + v[2]) );
			s += String.format("azimuth (right): %8.5fpx %8.5fpx %8.5fpx %8.5fpx (shift of he image center)\n" , v[3], v[4], v[5], -(v[3] + v[4] + v[5]) );
			s += String.format("roll    (CW):    %8.5fpx %8.5fpx %8.5fpx %8.5fpx (shift at the image half-width from the center)\n" , v[6], v[7], v[8], v[9] );
			s += String.format("diff zoom (in):  %8.5fpx %8.5fpx %8.5fpx %8.5fpx (shift at the image half-width from the center)\n" , v[10], v[11],  v[12], -(v[10] + v[11] + v[12]) );
			s += "Symmetrical vector:\n";
			s += "    |↘ ↙|     |↘ ↗|     |↗ ↘|     |↙ ↘|      |↙ ↗|     |↖  ↘| 6: common roll 7:(r0-r3)/2,           |- +|     |- -|     |- +|\n";
			s += " 0: |↗ ↖|  1: |↙ ↖|  2: |↖ ↙|  3: |↖ ↗|  4:  |↗ ↙|  5: |↘  ↖| 8:(r1-r2)/2    9:(r0+r3-r1-r2)/4  10: |- +| 11: |+ +| 12: |+ -|\n";

			s += String.format(" 0:%9.6fpx 1:%8.5fpx 2:%8.5fpx 3:%8.5fpx 4:%8.5fpx 5:%8.5fpx 6:%8.5fpx 7:%8.5fpx 8:%8.5fpx 9:%8.5fpx 10: %8.5fpx 11:%8.5fpx 12:%8.5fpx\n" ,
					sv[0], sv[1], sv[2], sv[3], sv[4], sv[5], sv[6], sv[7], sv[8], sv[9], sv[10], sv[11], sv[12] );
			return s;
		}

		public String toStringDegrees()
		{
			String s;
			double [] sym_vect = toSymArray(null);
			double [] v = new double [vector.length];
			double [] sv = new double [vector.length];
			for (int i = 0; i < LENGTH; i++){
				if (i < LENGTH_ANGLES) {
					v[i] =  vector[i]*180/Math.PI;
					sv[i] = sym_vect[i]*180/Math.PI;
				} else {
					v[i] =  vector[i];
					sv[i] = sym_vect[i];
				}
			}

			s  = String.format("tilt    (up):    %8.5f° %8.5f° %8.5f° %8.5f°\n" , v[0], v[1], v[2], -(v[0] + v[1] + v[2]) );
			s += String.format("azimuth (right): %8.5f° %8.5f° %8.5f° %8.5f°\n" , v[3], v[4], v[5], -(v[3] + v[4] + v[5]) );
			s += String.format("roll    (CW):    %8.5f° %8.5f° %8.5f° %8.5f°\n" , v[6], v[7], v[8], v[9] );
			s += String.format("diff zoom (in):  %8.5f‰ %8.5f‰ %8.5f‰ %8.5f‰\n" , 1000*v[10],1000*v[11],1000*v[12], -1000*(v[10] + v[11] + v[12]) );
			s += "Symmetrical vector:\n";
			s += "    |↘ ↙|     |↘ ↗|     |↗ ↘|     |↙ ↘|      |↙ ↗|     |↖  ↘| 6: common roll 7:(r0-r3)/2,           |- +|     |- -|     |- +|\n";
			s += " 0: |↗ ↖|  1: |↙ ↖|  2: |↖ ↙|  3: |↖ ↗|  4:  |↗ ↙|  5: |↘  ↖| 8:(r1-r2)/2    9:(r0+r3-r1-r2)/4  10: |- +| 11: |+ +| 12: |+ -|\n";

			s += String.format(" 0:%9.6f° 1:%8.5f° 2:%8.5f° 3:%8.5f° 4:%8.5f° 5:%8.5f° 6:%8.5f° 7:%8.5f° 8:%8.5f° 9:%8.5f° 10: %8.5f‰ 11:%8.5f‰ 12:%8.5f‰\n" ,
					sv[0], sv[1], sv[2], sv[3], sv[4], sv[5], sv[6], sv[7], sv[8], sv[9], 1000*sv[10], 1000*sv[11], 1000*sv[12] );

			return s;
		}




		public void incrementVector(double [] incr,
				double scale)
		{
			for (int i = 0; i < incr.length; i++){
				vector[i]+= incr[i] * scale;
			}
		}

		public void incrementVector(CorrVector incr, double scale)
		{
			incrementVector(incr.toArray(), scale);
		}

		@Override
		public CorrVector clone(){
			return new CorrVector(this.vector.clone());
		}

		/**
		 * Convert manual pixel shift between images to azimuth/tilt rotations (distortions ignored)
		 * and apply (add) them to the current vector (normally should be all 0.0)
		 * @param pXY_shift manula XY pixel corrections (shiftXY made of clt_parameters.fine_corr_[xy]_[0123])
		 */
		public void applyPixelShift(double [][] pXY_shift){
			double [] pXY_avg = {0.0,0.0};
			for (int i = 0; i < numSensors; i++){
				for (int j = 0; j < 2; j++) {
					pXY_avg[j] += pXY_shift[i][j]/numSensors;
				}
			}
			double pix_to_rads = (0.001*pixelSize)/focalLength;
			for (int i = 0; i < (numSensors -1); i++){
				vector[TILT_INDEX + i]    -= pix_to_rads * (pXY_shift[i][1] - pXY_avg[1]);
				vector[AZIMUTH_INDEX + i] += pix_to_rads * (pXY_shift[i][0] - pXY_avg[0]);
			}
		}

		/**
		 * Get conversion matrix from symmetrical coordinates
		 * @return
		 *
		 *     |↘ ↙|     |↘ ↗|     |↗ ↘|     |↙ ↘|      |↙ ↗|     |↖  ↘|
		 *  0: |↗ ↖|  1: |↙ ↖|  2: |↖ ↙|  3: |↖ ↗|  4:  |↗ ↙|  5: |↘  ↖|
		 *
		 */
		public double [][] dSym_j_dTar_i()
		{
			double [][] tar_to_sym = {
					{-2.0, -2.0,  2.0, -2.0,  0.0,  0.0, 0.0,   0.0,  0.0,  0.0,   0.0,  0.0,  0.0}, // t0
					{-2.0,  0.0,  0.0, -2.0,  2.0, -2.0, 0.0,   0.0,  0.0,  0.0,   0.0,  0.0,  0.0}, // t1
					{ 0.0, -2.0,  2.0,  0.0,  2.0, -2.0, 0.0,   0.0,  0.0,  0.0,   0.0,  0.0,  0.0}, // t2

					{ 2.0,  2.0,  2.0, -2.0,  0.0,  0.0, 0.0,   0.0,  0.0,  0.0,   0.0,  0.0,  0.0}, // a0
					{ 0.0,  2.0,  2.0,  0.0,  2.0,  2.0, 0.0,   0.0,  0.0,  0.0,   0.0,  0.0,  0.0}, // a1
					{ 2.0,  0.0,  0.0, -2.0,  2.0,  2.0, 0.0,   0.0,  0.0,  0.0,   0.0,  0.0,  0.0}, // a2

					{ 0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.25,  0.5,  0.0,  0.25,  0.0,  0.0,  0.0}, // roll 0
					{ 0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.25,  0.0,  0.5, -0.25,  0.0,  0.0,  0.0}, // roll 1
					{ 0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.25,  0.0, -0.5, -0.25,  0.0,  0.0,  0.0}, // roll 2
					{ 0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.25, -0.5,  0.0,  0.25,  0.0,  0.0,  0.0}, // roll 3

					{ 0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.0,   0.0,  0.0,  0.0,  -1.0, -1.0,  0.0}, // scale 0
					{ 0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.0,   0.0,  0.0,  0.0,   0.0, -1.0,  1.0}, // scale 1
					{ 0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.0,   0.0,  0.0,  0.0,  -1.0,  0.0,  1.0}  // scale 2
			};

			return tar_to_sym;
		}
		public double [][] dTar_j_dSym_i(){
			//			Matrix sym_to_tar = new Matrix(symToTARArray());
			//			Matrix tar_to_sym = sym_to_tar.inverse();
			//			return tar_to_sym.getArray();
			/*
a=[
[-2.0, -2.0,  2.0, -2.0,  0.0,  0.0, 0.0,   0.0,  0.0,  0.0,  0.0,  0.0,  0.0],
[-2.0,  0.0,  0.0, -2.0,  2.0, -2.0, 0.0,   0.0,  0.0,  0.0,  0.0,  0.0,  0.0],
[ 0.0, -2.0,  2.0,  0.0,  2.0, -2.0, 0.0,   0.0,  0.0,  0.0,  0.0,  0.0,  0.0],
[ 2.0,  2.0,  2.0, -2.0,  0.0,  0.0, 0.0,   0.0,  0.0,  0.0,  0.0,  0.0,  0.0],
[ 0.0,  2.0,  2.0,  0.0,  2.0,  2.0, 0.0,   0.0,  0.0,  0.0,  0.0,  0.0,  0.0],
[ 2.0,  0.0,  0.0, -2.0,  2.0,  2.0, 0.0,   0.0,  0.0,  0.0,  0.0,  0.0,  0.0],
[ 0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.25,  0.5,  0.0,  0.25, 0.0,  0.0,  0.0],
[ 0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.25,  0.0,  0.5, -0.25, 0.0,  0.0,  0.0],
[ 0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.25,  0.0, -0.5, -0.25, 0.0,  0.0,  0.0],
[ 0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.25, -0.5,  0.0,  0.25, 0.0,  0.0,  0.0],
[ 0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.0,   0.0,  0.0,  0.0, -1.0, -1.0,  0.0],
[ 0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.0,   0.0,  0.0,  0.0,  0.0, -1.0,  1.0],
[ 0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.0,   0.0,  0.0,  0.0, -1.0,  0.0,  1.0]]


matrix([[-0.125, -0.125,  0.125,  0.125, -0.125,  0.125, -0.   , -0.   ,   -0.   , -0.   , -0.   , -0.   , -0.   ],
        [-0.125,  0.125, -0.125,  0.125,  0.125, -0.125,  0.   ,  0.   ,    0.   ,  0.   ,  0.   ,  0.   ,  0.   ],
        [ 0.125, -0.125,  0.125,  0.125,  0.125, -0.125,  0.   ,  0.   ,    0.   ,  0.   ,  0.   ,  0.   ,  0.   ],
        [-0.125, -0.125,  0.125, -0.125,  0.125, -0.125,  0.   ,  0.   ,    0.   ,  0.   ,  0.   ,  0.   ,  0.   ],
        [-0.125,  0.125,  0.125, -0.125,  0.125,  0.125,  0.   ,  0.   ,    0.   ,  0.   ,  0.   ,  0.   ,  0.   ],
        [ 0.125, -0.125, -0.125, -0.125,  0.125,  0.125,  0.   ,  0.   ,    0.   ,  0.   ,  0.   ,  0.   ,  0.   ],
        [ 0.   ,  0.   ,  0.   ,  0.   ,  0.   ,  0.   ,  1.   ,  1.   ,    1.   ,  1.   ,  0.   ,  0.   ,  0.   ],
        [-0.   , -0.   , -0.   , -0.   , -0.   , -0.   ,  1.   , -0.   ,   -0.   , -1.   , -0.   , -0.   , -0.   ],
        [-0.   , -0.   , -0.   , -0.   , -0.   , -0.   , -0.   ,  1.   ,   -1.   , -0.   , -0.   , -0.   , -0.   ],
        [-0.   , -0.   , -0.   , -0.   , -0.   , -0.   ,  1.   , -1.   ,   -1.   ,  1.   , -0.   , -0.   , -0.   ],
        [-0.   , -0.   , -0.   , -0.   , -0.   , -0.   , -0.   , -0.   ,   -0.   , -0.   , -0.5  ,  0.5  , -0.5  ],
        [-0.   , -0.   , -0.   , -0.   , -0.   , -0.   , -0.   , -0.   ,   -0.   , -0.   , -0.5  , -0.5  ,  0.5  ],
        [ 0.   ,  0.   ,  0.   ,  0.   ,  0.   ,  0.   ,  0.   ,  0.   ,    0.   ,  0.   , -0.5  ,  0.5  ,  0.5  ]])
*/
			double [][] sym_to_tar=	{
					// t0     t1     t2     a0     a1     a2     r0    r1    r2    r3    s0    s1    s2
					{-0.125,-0.125, 0.125, 0.125,-0.125, 0.125, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0 },  // sym0
			        {-0.125, 0.125,-0.125, 0.125, 0.125,-0.125, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0 },  // sym1
			        { 0.125,-0.125, 0.125, 0.125, 0.125,-0.125, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0 },  // sym2
			        {-0.125,-0.125, 0.125,-0.125, 0.125,-0.125, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0 },  // sym3
			        {-0.125, 0.125, 0.125,-0.125, 0.125, 0.125, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0 },  // sym4
			        { 0.125,-0.125,-0.125,-0.125, 0.125, 0.125, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0 },  // sym5
			        { 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 1.0,  1.0,  1.0,  1.0,  0.0,  0.0,  0.0 },  // sym6  = (r0+r1+r2+r3)/4
			        { 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 1.0,  0.0,  0.0, -1.0,  0.0,  0.0,  0.0 },  // sym7  = (r0-r3)/2
			        { 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0,  1.0, -1.0,  0.0,  0.0,  0.0,  0.0 },  // sym8  = (r1-r2)/2
			        { 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 1.0, -1.0, -1.0,  1.0,  0.0,  0.0,  0.0 },  // sym9  = (r0+r3-r1-r2)/4
			        { 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0,  0.0, -1.0,  1.0, -0.5,  0.5, -0.5 },  // sym10 = -s0 - s2
			        { 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0,  0.0, -1.0,  1.0, -0.5, -0.5,  0.5 },  // sym11 = -s0 - s1
			        { 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0,  0.0, -1.0,  1.0, -0.5,  0.5,  0.5 }   // sym12 =  s1 + s2
			        };
			return sym_to_tar;
		}

		public boolean [] getParMask(
				boolean use_disparity,
//				boolean use_other_extr,
				boolean use_aztilts,       // Adjust azimuths and tilts excluding disparity
				boolean use_diff_rolls,    // Adjust differential rolls (3 of 4 angles)
//				boolean disparity_only,
//				boolean use_disparity,
				boolean common_roll,
				boolean corr_focalLength,
		  		int     manual_par_sel)    // Manually select the parameter mask bit 0 - sym0, bit1 - sym1, ... (0 - use boolean flags, != 0 - ignore boolean flags)

		{
//			common_roll &=      !disparity_only;
//			corr_focalLength &= !disparity_only;
//			use_disparity |=     disparity_only;
			boolean [] par_mask = {
					use_disparity,    //sym0
					use_aztilts,  //sym1
					use_aztilts,  //sym2
					use_aztilts,  //sym3
					use_aztilts,  //sym4
					use_aztilts,  //sym5
					common_roll,      //sym6 // common roll
					use_diff_rolls,  //sym7
					use_diff_rolls,  //sym8
					use_diff_rolls,  //sym9
					corr_focalLength, //sym10
					corr_focalLength, //sym11
					corr_focalLength  //sym12
			};
			if (manual_par_sel != 0) {
				for (int i = 0; i < par_mask.length; i++) {
					par_mask[i] = ((manual_par_sel >> i) & 1) != 0;
				}
				System.out.println("*** Using manual parameter mask, overwriting boolean flags:");
				for (int i = 0; i < par_mask.length; i++) {
					System.out.println("Sym"+i+": "+par_mask[i]);
				}
			}
			return par_mask;
		}


		/**
		 * Get partial transposed Jacobian as 2d array (for one measurement set) from partial Jacobian for each sample
		 * with derivatives of port coordinates (all 4) by 3 tilts (ports 0..2), 3 azimuths (ports 0..2) and all 4 rolls
		 * Tilt and azimuth for port 3 is calculated so center would not move. Tilt is positive up, azimuth - right and
		 * roll - clockwise. Added zooms (difference from 1.0) for sensors 0..2
		 *
		 * Result is transposed Jacobian (rows (9 , 10,12 or 13) - parameters, columns - port coordinate components (8). Parameters
		 * here are symmetrical, 0 is disparity-related (all to the center), remaining 9 preserve disparity and
		 * are only responsible for the "lazy eye" (last 4  are were roll, now differential zooms are added):
		 *
		 *     |↘ ↙|     |↘ ↗|     |↗ ↘|     |↙ ↘|      |↙ ↗|     |↖  ↘|
		 *  0: |↗ ↖|  1: |↙ ↖|  2: |↖ ↙|  3: |↖ ↗|  4:  |↗ ↙|  5: |↘  ↖|
		 *
		 * @param port_coord_deriv result of getPortsCoordinatesAndDerivatives: first index: port0x, port0y... port3y,
		 *                         second index - parameters [tilt0, tilt1, ..., roll3]
		 * @param par_mask array of 10 elements - which parameters to use (normally all true or all but first
		 * @return
		 */

		public double [][] getJtPartial(
				double [][] port_coord_deriv,
				boolean [] par_mask)
		{
			int num_pars = 0;
			for (int npar = 0; npar < LENGTH; npar++) if ((par_mask==null) || par_mask[npar]) {
				num_pars++;
			}
			double [][] sym_to_tar=	dTar_j_dSym_i();
			double [][] jt_part = new double[num_pars][port_coord_deriv.length];
			int opar = 0;
			for (int npar = 0; npar < LENGTH; npar++) if ((par_mask==null) || par_mask[npar]) {
				for (int i = 0; i < port_coord_deriv.length; i++){ // pxy index (0..7)
					for (int k = 0; k < sym_to_tar[npar].length; k++){
						jt_part[opar][i] += sym_to_tar[npar][k]* port_coord_deriv[i][k]; // transposing port_coord_deriv
					}
				}
				opar++;
			}
			return jt_part;
		}
		// convert tilt0,... roll3 array to symmetrical coordinates [0] - to the center (disparity)
		public double [] toSymArray(boolean [] par_mask)
		{
			return toSymArray(this.vector, par_mask);
		}

		public double [] toSymArray(
				double [] tar_array,
				boolean [] par_mask)
		{
			double [][] tar_to_sym = dSym_j_dTar_i();
			int num_pars = 0;
			for (int npar = 0; npar < vector.length; npar++) if ((par_mask==null) || par_mask[npar]) {
				num_pars++;
			}
			double [] sym_array = new double  [num_pars];
			int opar = 0;
			for (int npar = 0; npar < LENGTH; npar++) if ((par_mask==null) || par_mask[npar]) {
				for (int i = 0; i < tar_array.length; i++){
					sym_array[opar] += tar_array[i] * tar_to_sym[i][npar];
				}
				opar++;
			}
			return sym_array;
		}

		public double [] toTarArray(
				double [] sym_array,
				boolean [] par_mask)
		{
			double [][] sym_to_tar = dTar_j_dSym_i();
			double [] tar_array = new double  [sym_to_tar[0].length];
			for (int npar = 0; npar < LENGTH; npar++)  {
				int spar = 0;
				for (int i = 0; i < LENGTH; i++) if ((par_mask==null) || par_mask[i]){
					tar_array[npar] += sym_array[spar] * sym_to_tar[i][npar];
					spar++;
				}
			}
			return tar_array;
		}
	}




	public void setDistortion(
			double focalLength,
			double distortionC,
			double distortionB,
			double distortionA,
			double distortionA5,
			double distortionA6,
			double distortionA7,
			double distortionA8,
			double distortionRadius,
			int    pixelCorrectionWidth,   // virtual camera center is at (pixelCorrectionWidth/2, pixelCorrectionHeight/2)
			int    pixelCorrectionHeight,
			double pixelSize

			)	{
		if (!Double.isNaN(focalLength))      this.focalLength = focalLength;
		if (!Double.isNaN(distortionC))      this.distortionC = distortionC;
		if (!Double.isNaN(distortionB))      this.distortionB = distortionB;
		if (!Double.isNaN(distortionA))      this.distortionA = distortionA;
		if (!Double.isNaN(distortionA5))     this.distortionA5 = distortionA5;
		if (!Double.isNaN(distortionA6))     this.distortionA6 = distortionA6;
		if (!Double.isNaN(distortionA7))     this.distortionA7 = distortionA7;
		if (!Double.isNaN(distortionA8))     this.distortionA8 = distortionA8;
		if (!Double.isNaN(distortionRadius)) this.distortionRadius = distortionRadius;
		if (pixelCorrectionWidth >= 0)       this.pixelCorrectionWidth = pixelCorrectionWidth;
		if (pixelCorrectionHeight >= 0)      this.pixelCorrectionHeight = pixelCorrectionHeight;
		if (!Double.isNaN(pixelSize))        this.pixelSize = pixelSize;
		//	imp.setProperty("distortion_formula",  "(normalized by distortionRadius in mm) Rdist/R=A8*R^7+A7*R^6+A6*R^5+A5*R^4+A*R^3+B*R^2+C*R+(1-A6-A7-A6-A5-A-B-C)");
		//	imp.setProperty("distortionRadius", ""+subCam.distortionRadius);
	}

	public void setSensors(
			int         numSensors, // <=0 - keep current
			double      elevation,  // NaN - keep
			double      heading,    // NaN - keep
			double []   forward,    // null - keep all, NaN - keep individual
			double []   right,      // null - keep all, NaN - keep individual
			double []   height,      // null - keep all, NaN - keep individual
			double []   roll,       // null - keep all, NaN - keep individual
			double [][] pXY0){      // null - keep all, [] null - keep individual

		if (numSensors > 0) this.numSensors = numSensors;
		if (!Double.isNaN(elevation))    this.elevation = elevation;
		if (!Double.isNaN(heading))      this.heading = heading;
		if (forward != null){
			if (forward.length != numSensors){
				throw new IllegalArgumentException ("forward.length ("+forward.length+") != numSensors ("+numSensors+")");
			}
			if ((this.forward == null) || (this.forward.length != numSensors)) this.forward = new double [numSensors];
			for (int i = 0; i < numSensors; i++) if (!Double.isNaN(forward[i]))      this.forward[i] = forward[i];
		}

		if (right != null){
			if (right.length != numSensors){
				throw new IllegalArgumentException ("right.length ("+right.length+") != numSensors ("+numSensors+")");
			}
			if ((this.right == null) || (this.right.length != numSensors)) this.right = new double [numSensors];
			for (int i = 0; i < numSensors; i++) if (!Double.isNaN(right[i]))      this.right[i] = right[i];
		}

		if (height != null){
			if (height.length != numSensors){
				throw new IllegalArgumentException ("height.length ("+height.length+") != numSensors ("+numSensors+")");
			}
			if ((this.height == null) || (this.height.length != numSensors)) this.height = new double [numSensors];
			for (int i = 0; i < numSensors; i++) if (!Double.isNaN(height[i]))      this.height[i] = height[i];
		}

		if (roll != null){
			if (roll.length != numSensors){
				throw new IllegalArgumentException ("roll.length ("+roll.length+") != numSensors ("+numSensors+")");
			}
			if ((this.roll == null) || (this.roll.length != numSensors)) this.roll = new double [numSensors];
			for (int i = 0; i < numSensors; i++) if (!Double.isNaN(roll[i]))      this.roll[i] = roll[i];
		}
		if (pXY0 != null){
			if (pXY0.length != numSensors){
				throw new IllegalArgumentException ("pXY0.length ("+pXY0.length+") != numSensors ("+numSensors+")");
			}
			if ((this.pXY0 == null) || (this.pXY0.length != numSensors)) this.pXY0 = new double [numSensors][];
			for (int i = 0; i < numSensors; i++) if (pXY0[i] != null)      this.pXY0[i] = pXY0[i].clone();
		}
	}

	public void planeProjectLenses(){ // calculate XYZ_he (any number of sensors)
		// get center of the adjusted camera
		common_right = 0;
		common_forward = 0;
		common_height = 0;
		for (int i = 0; i < numSensors; i++){
			common_right += right[i];
			common_forward += forward[i];
			common_height += height[i];
		}
		common_right   /= numSensors;
		common_forward /= numSensors;
		common_height  /= numSensors;
		//    	double [][]
		this.XYZ_he = new double [numSensors][3]; // after heading, then elevation rotation
		/*
    	rotate by phi around C2Y:Vc3= R3*Vc2
    	| cos(phi)   0   -sin(phi)   |   |X|
    	|     0      1         0     | * |Y|
    	| sin(phi)   0    cos(phi)   |   |Z|
		 */
		double c_head=  Math.cos(heading*Math.PI/180);
		double s_head=  Math.sin(heading*Math.PI/180);
		double [][] aR_head={{c_head,0.0,-s_head},{0.0,1.0,0.0},{s_head,0.0,c_head}};
		Matrix R_head=new Matrix(aR_head);
		/*
    	rotate by theta around C1X:Vc2= R2*Vc1
    	|    1        0         0        |   |X|
    	|    0   cos(theta)  -sin(theta) | * |Y|
    	|    0   sin(theta)   cos(theta) |   |Z|
		 */
		double c_elev=  Math.cos(elevation*Math.PI/180);
		double s_elev=  Math.sin(elevation*Math.PI/180);
		double [][] aR_elev={{1.0,0.0,0.0},{0.0,c_elev, -s_elev},{0.0, s_elev, c_elev}};
		Matrix R_elev=new Matrix(aR_elev);
		Matrix R_head_elev = R_elev.times(R_head);

		for (int i = 0; i<numSensors; i++){
			double [][] aXYZi_ccs = {
					{  right[i] -   common_right},
					{- (height[i] - common_height)},
					{  forward[i] - common_forward}};
			Matrix XYZi_ccs = new Matrix(aXYZi_ccs);
			Matrix mXYZ_he = R_head_elev.times(XYZi_ccs);
			for (int j = 0; j<3;j++) this.XYZ_he[i][j] = mXYZ_he.get(j, 0);
		}
		// Calculate average radius
		cameraRadius = 0;
		for (int i = 0; i < numSensors; i++){
			cameraRadius += this.XYZ_he[i][0] * this.XYZ_he[i][0] + this.XYZ_he[i][1] * this.XYZ_he[i][1];
		}
		cameraRadius = Math.sqrt(cameraRadius/numSensors);
	}

	// cameras should be Z-numbered (looking to the target, X - right, Y - down)
	public void adustSquare(){ // rotate heading/elevation aligned cameras around the Z-axis to make it more "square"
		if (numSensors != 4 ){
			throw new IllegalArgumentException ("adjustSquare() is valid only for quad-cameras, numSensors="+numSensors);
		}
		this.disparityRadius = Math.sqrt(2.0) * this.cameraRadius;
		double Sx = - XYZ_he[0][1] + XYZ_he[1][0] - XYZ_he[2][0] + XYZ_he[3][1];
		double Sy = - XYZ_he[0][0] - XYZ_he[1][1] + XYZ_he[2][1] + XYZ_he[3][0];
		double psi = 0.25*Math.PI - Math.atan2(Sy, Sx);
		common_roll = psi*180/Math.PI;
		/*
    	Converting from the sub-camera coordinates to the target coordinates
    	rotate by -psi around CZ
    	| cos(psi)  sin(psi)    0  |   |Xc0|
    	|-sin(psi)  cos(psi)    0  | * |Yc0|
    	|    0         0        1  |   |Zc0|
		 */
		double c_roll=  Math.cos(psi*Math.PI/180);
		double s_roll=  Math.sin(psi*Math.PI/180);

		double [][] aR_roll={
				{ c_roll, s_roll, 0.0},
				{-s_roll, c_roll, 0.0},
				{ 0.0,    0.0,    1.0}};
		Matrix R_roll = new Matrix(aR_roll);
		this.XYZ_her =  new double [numSensors][3];
		this.rXY =      new double [numSensors][2]; // XY pairs of the in a normal plane, relative to disparityRadius
		for (int i = 0; i<numSensors; i++){
			double [][] aXYZi_he = {
					{this.XYZ_he[i][0]},
					{this.XYZ_he[i][1]},
					{this.XYZ_he[i][2]}};
			Matrix mXYZi_he = new Matrix(aXYZi_he);
			Matrix mXYZ_her = R_roll.times(mXYZi_he);
			for (int j = 0; j<3;j++) this.XYZ_her[i][j] = mXYZ_her.get(j, 0);
			for (int j = 0; j<2;j++) this.rXY[i][j] = this.XYZ_her[i][j]/this.disparityRadius;
		}
		if (rigOffset != null) {
			rigOffset.recalcRXY();
		}
	}

	public void listGeometryCorrection(boolean showAll){
		System.out.println("'=== Constant parameters ===");
		System.out.println("pixelCorrectionWidth =\t"+  pixelCorrectionWidth+"\tpix");
		System.out.println("pixelCorrectionHeight =\t"+ pixelCorrectionHeight+"\tpix");
		System.out.println("pixelSize =\t"+             pixelSize+"\tum");
		System.out.println("distortionRadius =\t"+      distortionRadius+"\tmm");
		System.out.println("'=== Common input parameters ===");
		System.out.println("focalLength =\t"+  focalLength + " mm");
		System.out.println("distortionA8 =\t"+ distortionA8);
		System.out.println("distortionA7 =\t"+ distortionA7);
		System.out.println("distortionA6 =\t"+ distortionA6);
		System.out.println("distortionA5 =\t"+ distortionA5);
		System.out.println("distortionA =\t"+  distortionA);
		System.out.println("distortionB =\t"+  distortionB);
		System.out.println("distortionC =\t"+  distortionC);
		System.out.println("elevation =\t"+    elevation+"\tdegrees");
		System.out.println("heading =\t"+      heading+"\tdegrees");
		System.out.println("numSensors =\t"+   numSensors);
		System.out.println("'=== Individual input parameters ===");
		System.out.print  ("forward = ");for (int i = 0; i < numSensors;i++) System.out.print("\t"+forward[i]); System.out.println("\tmm");
		System.out.print  ("right = ");  for (int i = 0; i < numSensors;i++) System.out.print("\t"+right[i]);   System.out.println("\tmm");
		System.out.print  ("height = "); for (int i = 0; i < numSensors;i++) System.out.print("\t"+height[i]);  System.out.println("\tmm");
		System.out.print  ("roll = ");   for (int i = 0; i < numSensors;i++) System.out.print("\t"+roll[i]);      System.out.println("\tdegrees");
		System.out.print  ("px0 = ");    for (int i = 0; i < numSensors;i++) System.out.print("\t"+pXY0[i][0]);  System.out.println("\tpix");
		System.out.print  ("py0 = ");    for (int i = 0; i < numSensors;i++) System.out.print("\t"+pXY0[i][1]);  System.out.println("\tpix");

		System.out.println("'=== Common calculated parameters ===");
		System.out.println("common_right =\t"+common_right + "\tmm");
		System.out.println("common_forward =\t"+common_forward + "\tmm");
		System.out.println("common_height =\t"+common_height + "\tmm");
		System.out.println("common_roll =\t"+common_roll + "\tdegrees");
		System.out.println("cameraRadius =\t"+cameraRadius + "\tmm");
		System.out.println("disparityRadius =\t"+disparityRadius + "\tmm");

		if (showAll){
			System.out.println("'=== Intermediate data: coordinates corrected for common elevation and heading ===");
			System.out.print  ("X_he =");    for (int i = 0; i < numSensors;i++) System.out.print("\t"+XYZ_he[i][0]);  System.out.println("\tmm");
			System.out.print  ("Y_he =");    for (int i = 0; i < numSensors;i++) System.out.print("\t"+XYZ_he[i][1]);  System.out.println("\tmm");
			System.out.print  ("Z_he =");    for (int i = 0; i < numSensors;i++) System.out.print("\t"+XYZ_he[i][2]);  System.out.println("\tmm");
			System.out.println("'=== Intermediate data: coordinates corrected for common elevation, heading and roll ===");
			System.out.print  ("X_her =");    for (int i = 0; i < numSensors;i++) System.out.print("\t"+XYZ_her[i][0]);  System.out.println("\tmm");
			System.out.print  ("Y_her =");    for (int i = 0; i < numSensors;i++) System.out.print("\t"+XYZ_her[i][1]);  System.out.println("\tmm");
			System.out.print  ("Z_her =");    for (int i = 0; i < numSensors;i++) System.out.print("\t"+XYZ_her[i][2]);  System.out.println("\tmm");
		}

		System.out.println("'=== Individual calculated parameters ===");
		System.out.print  ("residual_roll = ");   for (int i = 0; i < numSensors;i++) System.out.print("\t"+(roll[i]-common_roll));System.out.println("\tdegrees");
		System.out.print  ("X_rel =");    for (int i = 0; i < numSensors;i++) System.out.print("\t"+rXY[i][0]);  System.out.println("\trelative to disparityRadius");
		System.out.print  ("Y_rel =");    for (int i = 0; i < numSensors;i++) System.out.print("\t"+rXY[i][1]);  System.out.println("\trelative to disparityRadius");
	}

	// return distance from disparity (in pixel units) for the current camera geometry
	public double getZFromDisparity(double disparity){
		return SCENE_UNITS_SCALE * this.focalLength * this.disparityRadius / (disparity * 0.001*this.pixelSize);
	}

	public double getDisparityFromZ(double z){
		return (1000.0 * SCENE_UNITS_SCALE * this.focalLength * this.disparityRadius / this.pixelSize) / z;
	}


	public double getFOVPix(){ // get ratio of 1 pixel X/Y to Z (distance to object)
		return 0.001 * this.pixelSize / this.focalLength;
	}

	public double getFOVWidth(){ // get FOV ratio: width to distance
		return this.pixelCorrectionWidth * 0.001 * this.pixelSize / this.focalLength;
	}
	public double getFOVHeight(){ // get FOV ratio: width to distance
		return this.pixelCorrectionHeight * 0.001 * this.pixelSize / this.focalLength;
	}

	public double getScaleDzDx()
	{
		return ( 0.001 * this.pixelSize) / this.focalLength;
	}


	/**
	 * Get real world coordinates from pixel coordinates and nominal disparity
	 * @param px horizontal pixel coordinate (right)
	 * @param py vertical pixel coordinate (down)
	 * @param disparity nominal disparity (pixels)
	 * @param correctDistortions true: correct lens distortions, false - no lens distortions
	 * @return {x, y, z} in meters
	 */
	public double [] getWorldCoordinates(
			double px,
			double py,
			double disparity,
			boolean correctDistortions) // correct distortion (will need corrected background too !)
	{
		double pXcd = px - 0.5 * this.pixelCorrectionWidth;
		double pYcd = py - 0.5 * this.pixelCorrectionHeight;
		double rD = Math.sqrt(pXcd*pXcd + pYcd*pYcd)*0.001*this.pixelSize; // distorted radius in a virtual center camera
		double rND2R = correctDistortions?(getRByRDist(rD/this.distortionRadius, false)): 1.0;
		double pXc = pXcd * rND2R; // non-distorted coordinates relative to the (0.5 * this.pixelCorrectionWidth, 0.5 * this.pixelCorrectionHeight)
		double pYc = pYcd * rND2R; // in pixels
		double z = -SCENE_UNITS_SCALE * this.focalLength * this.disparityRadius / (disparity * 0.001*this.pixelSize); // "+" - near, "-" far
		double x =  SCENE_UNITS_SCALE * pXc * this.disparityRadius / disparity;
		double y = -SCENE_UNITS_SCALE * pYc * this.disparityRadius / disparity;
		double [] xyz = {x,y,z};
		return xyz;
	}

	/**
	 * Find disparity for the intersection of the view ray (px, py) and a real-world plane orthogonal through the end of the
	 * vector norm_xyz
	 * @param norm_xyz vector from the origin (camera) orthogonal to the plane, length is a distance to the plane
	 * @param px pixel coordinate horizontal
	 * @param py pixel coordinate vertical
	 * @param correctDistortions true for lens distortion correction, false otherwise
	 * @return disparity for the point on the plane specified by norm_xyz and known view coordinates px, py
	 */
	public double getPlaneDisparity(
			double [] norm_xyz,
			double px,
			double py,
			boolean correctDistortions) // correct distortion (will need corrected background too !)
	{
		double pXcd = px - 0.5 * this.pixelCorrectionWidth;
		double pYcd = py - 0.5 * this.pixelCorrectionHeight;
		double rD = Math.sqrt(pXcd*pXcd + pYcd*pYcd)*0.001*this.pixelSize; // distorted radius in a virtual center camera
		double rND2R = correctDistortions?(getRByRDist(rD/this.distortionRadius, false)): 1.0;
		double pXc = pXcd * rND2R; // non-distorted coordinates relative to the (0.5 * this.pixelCorrectionWidth, 0.5 * this.pixelCorrectionHeight)
		double pYc = pYcd * rND2R; // in pixels
		// point for the unity disparity
		double x =  SCENE_UNITS_SCALE * pXc * this.disparityRadius;
		double y = -SCENE_UNITS_SCALE * pYc * this.disparityRadius;
		double z = -SCENE_UNITS_SCALE * this.focalLength * this.disparityRadius / (0.001*this.pixelSize); // "+" - near, "-" far
		double vect_dot_norm = (x * norm_xyz[0]) +           (y * norm_xyz[1]) +           (z * norm_xyz[2]);
		double norm_dot_norm = (norm_xyz[0] * norm_xyz[0]) + (norm_xyz[1] * norm_xyz[1]) + (norm_xyz[2] * norm_xyz[2]);
		return vect_dot_norm / norm_dot_norm;
	}

	/* Just for testing using delta instead of d */
	public double [][] getWorldJacobian(
			double px,
			double py,
			double disparity,
			boolean correctDistortions, // correct distortion (will need corrected background too !)
			double delta)
	{
		double [] dpxpy0 = {disparity,px,py};
		double [][] jacobian = new double [3][3];
		double [] xyz0 = getWorldCoordinates(dpxpy0[1],dpxpy0[2],dpxpy0[0],correctDistortions);
		for (int i = 0; i< 3; i++){
			double [] dpxpy = dpxpy0.clone();
			dpxpy[i] += delta;
			double [] xyz = getWorldCoordinates(dpxpy[1],dpxpy[2],dpxpy[0],correctDistortions);
			for (int j = 0; j<3; j++){
				jacobian[j][i] = (xyz[j] - xyz0[j])/delta;
			}
		}
		return jacobian;
	}


	/**
	 * Get jacobian matrix - derivatives of real world [x,y,z] (in meters, right, up, to camera) by [disparity, px,py] (in pixels disparity, right, down)
	 * @param px horizontal pixel coordinate (right)
	 * @param py vertical pixel coordinate (down)
	 * @param disparity nominal disparity (pixels)
	 * @param correctDistortions true: correct lens distortions, false - no lens distortions
	 * @return {{dx/ddisparity, dx/dpx, dx/dpy},{dy/ddisparity, dy/dpx, dy/dpy},{dz/ddisparity, dz/dpx, dz/dpy}}
	 */
	public double [][] getWorldJacobian(
			double px,
			double py,
			double disparity,
			boolean correctDistortions, // correct distortion (will need corrected background too !)
			boolean debug)
	{
		double pXcd = px - 0.5 * this.pixelCorrectionWidth;
		double pYcd = py - 0.5 * this.pixelCorrectionHeight;

		double d_pXcd_d_px = 1.0;
		double d_pYcd_d_py = 1.0;

		double pXcd2_pYcd2 = pXcd*pXcd + pYcd*pYcd;
		double rD = Math.sqrt(pXcd2_pYcd2)*0.001*this.pixelSize; // distorted radius in a virtual center camera

		double rrD = rD / this.distortionRadius;
		double d_rRD_d_px = pXcd * rD / pXcd2_pYcd2  / this.distortionRadius;
		double d_rRD_d_py = pYcd * rD / pXcd2_pYcd2  / this.distortionRadius;



		double rND2R = correctDistortions?(getRByRDist(rrD, false)): 1.0;



		double rrND = rND2R * rrD; // relative to distortion radius non-distorted radius


		// k = rD/r
		double d_k_d_rrND = correctDistortions?getDerivRDistFromR(rrND):0.0;
		double d_rND2R_d_rrD = - rND2R * rND2R * d_k_d_rrND / ( d_k_d_rrND * rrND + 1.0/ rND2R); // rrND);

		double d_rND2R_d_px = d_rND2R_d_rrD * d_rRD_d_px;
		double d_rND2R_d_py = d_rND2R_d_rrD * d_rRD_d_py;

		double pXc = pXcd * rND2R; // non-distorted coordinates relative to the (0.5 * this.pixelCorrectionWidth, 0.5 * this.pixelCorrectionHeight)
		double pYc = pYcd * rND2R; // in pixels

		double d_pXc_d_px = pXcd * d_rND2R_d_px + d_pXcd_d_px * rND2R; // incorrect
		double d_pYc_d_py = pYcd * d_rND2R_d_py + d_pYcd_d_py * rND2R; // incorrect

		double d_pXc_d_py = pXcd * d_rND2R_d_py; // incorrect (too small)?
		double d_pYc_d_px = pYcd * d_rND2R_d_px; // incorrect (too small)?

		double z =  -SCENE_UNITS_SCALE * this.focalLength * this.disparityRadius / (disparity * 0.001*this.pixelSize); // "+" - near, "-" far
		double kx =  SCENE_UNITS_SCALE * this.disparityRadius / disparity;
		double ky = -SCENE_UNITS_SCALE * this.disparityRadius / disparity;
		double x =  kx * pXc;
		double y =  ky * pYc;

		double d_z_d_disparity = -z / disparity;
		double d_x_d_disparity = -x / disparity;
		double d_y_d_disparity = -y / disparity;
		// d_z_d_px == d_z_d_py ==0
		double [][] jacobian =
			{		{d_x_d_disparity, kx * d_pXc_d_px, kx * d_pXc_d_py},
					{d_y_d_disparity, ky * d_pYc_d_px, ky * d_pYc_d_py},
					{d_z_d_disparity, 0.0,             0.0}};
		return jacobian; // xyz;
	}

	/**
	 * Get pixel disparity and coordinates from the real world coordinates (in meters)
	 * @param xyz real world coordinates {x, y, z} in meters (right up, towards camera)
	 * @param correctDistortions true: correct lens distortions, false - no lens distortions
	 * @return {disparity, px, py} (right, down)
	 */
	public double [] getImageCoordinates(
			double [] xyz,
			boolean correctDistortions) // correct distortion (will need corrected background too !)
	{
		double x = xyz[0];
		double y = xyz[1];
		double z = xyz[2];
		double disparity = -SCENE_UNITS_SCALE * this.focalLength * this.disparityRadius / (z * 0.001*this.pixelSize);
		// non-distorted coordinates relative to the (0.5 * this.pixelCorrectionWidth, 0.5 * this.pixelCorrectionHeight)in mm
		double pXc = x * disparity / (SCENE_UNITS_SCALE * this.disparityRadius); // pixels
		double pYc =-y * disparity / (SCENE_UNITS_SCALE * this.disparityRadius); // pixels
		double rND = Math.sqrt(pXc*pXc + pYc*pYc)*0.001*this.pixelSize; // mm
		double rD2RND = correctDistortions?getRDistByR(rND/this.distortionRadius):1.0;
		double px = pXc * rD2RND + 0.5 * this.pixelCorrectionWidth;  // distorted coordinates relative to the (0.5 * this.pixelCorrectionWidth, 0.5 * this.pixelCorrectionHeight)
		double py = pYc * rD2RND + 0.5 * this.pixelCorrectionHeight; // in pixels
		double [] dxy = {disparity, px, py};
		return dxy;
	}

	/* Just for testing using delta instead of d */
	public double [][] getImageJacobian(
			double [] xyz0,
			boolean correctDistortions, // correct distortion (will need corrected background too !)
			double delta)
	{
		double [][] jacobian = new double [3][3];
		double [] dpxpy0 = getImageCoordinates(xyz0,correctDistortions);
		for (int i = 0; i< 3; i++){
			double [] xyz = xyz0.clone();
			xyz[i] += delta;
			double [] dpxpy = getImageCoordinates(xyz,correctDistortions);
			for (int j = 0; j<3; j++){
				jacobian[j][i] = (dpxpy[j] - dpxpy0[j])/delta;
			}
		}
		return jacobian;
	}


	/**
	 * Get jacobian matrix - derivatives of [disparity, px,py] (in pixels disparity, right, down) by real world [x,y,z] (in meters, right, up, to camera)
	 * @param xyz real world coordinates {x, y, z} in meters (right up, towards camera)
	 * @param correctDistortions true: correct lens distortions, false - no lens distortions
	 * @return {{dx/ddisparity, dx/dpx, dx/dpy},{dy/ddisparity, dy/dpx, dy/dpy},{dz/ddisparity, dz/dpx, dz/dpy}}
	 */
	public double [][] getImageJacobian(
			double [] xyz,
			boolean correctDistortions,
			int debugLevel)
	{
		if (debugLevel > 1){
			System.out.println("getImageJacobian():");
		}
		double x = xyz[0];
		double y = xyz[1];
		double z = xyz[2];
		double disparity = -SCENE_UNITS_SCALE * this.focalLength * this.disparityRadius / (z * 0.001*this.pixelSize);

		double d_disparity_d_z = -disparity / z; // no dependence on x,y
		// non-distorted coordinates relative to the (0.5 * this.pixelCorrectionWidth, 0.5 * this.pixelCorrectionHeight)in mm
		double pXc = x * disparity / (SCENE_UNITS_SCALE * this.disparityRadius);
		double pYc =-y * disparity / (SCENE_UNITS_SCALE * this.disparityRadius);

		double d_pXc_d_z = -pXc / z; // pXc/disparity * d_disparity_d_z = pXc/disparity * (-disparity / z);
		double d_pYc_d_z = -pYc / z; // pYc/disparity * d_disparity_d_z = pYc/disparity * (-disparity / z);

		double d_pXc_d_x =  disparity / (SCENE_UNITS_SCALE * this.disparityRadius);
		double d_pYc_d_y = -disparity / (SCENE_UNITS_SCALE * this.disparityRadius);


		double rND2 = pXc*pXc + pYc*pYc;
		double rrND = Math.sqrt(rND2)*0.001*this.pixelSize/this.distortionRadius; // relative to distortion radius;


		double d_rrND_d_pXc =  pXc*rrND/rND2;
		double d_rrND_d_pYc =  pYc*rrND/rND2;


		double rD2RND = correctDistortions?getRDistByR(rrND):1.0;
		double d_rD2RND_d_rrND = correctDistortions?getDerivRDistFromR(rrND) : 0.0;
/*
		double d_rD2RND_d_rrND0 = correctDistortions?getDerivRDistFromR(rrND,0.000001) : 0.0;
		if (debugLevel > 0){
			System.out.println("getImageJacobian():, d_rD2RND_d_rrND0 = "+d_rD2RND_d_rrND0+" ("+0.000001+"), d_rD2RND_d_rrND="+d_rD2RND_d_rrND);
		}
*/
		double p_px_d_x = d_pXc_d_x * (pXc * (d_rD2RND_d_rrND * d_rrND_d_pXc) + rD2RND); // / (0.001 * this.pixelSize);
		double p_py_d_y = d_pYc_d_y * (pYc * (d_rD2RND_d_rrND * d_rrND_d_pYc) + rD2RND); // / (0.001 * this.pixelSize);


		double p_px_d_y = pXc * (d_rD2RND_d_rrND * d_rrND_d_pYc * d_pYc_d_y); // / (0.001 * this.pixelSize);
		double p_py_d_x = pYc * (d_rD2RND_d_rrND * d_rrND_d_pXc * d_pXc_d_x); // / (0.001 * this.pixelSize);

		double d_rrND_d_z = d_rrND_d_pXc * d_pXc_d_z  + d_rrND_d_pYc * d_pYc_d_z;

		double p_px_d_z = (pXc * (d_rD2RND_d_rrND * d_rrND_d_z) + (d_pXc_d_z * rD2RND)); // / (0.001 * this.pixelSize);
		double p_py_d_z = (pYc * (d_rD2RND_d_rrND * d_rrND_d_z) + (d_pYc_d_z * rD2RND)); // / (0.001 * this.pixelSize);

		double [][] jacobian = {
				{     0.0,     0.0, d_disparity_d_z},
				{p_px_d_x,p_px_d_y, p_px_d_z},
				{p_py_d_x,p_py_d_y, p_py_d_z}};

		return jacobian;
	}

	/**
	 * Convert pixel coordinates to relative (from -1.0 to +1.0 each) to apply polynomial corrections
	 * @param pXY pair of pixel  X, pixel Y image coordinates
	 * @return pair of relative X, Y coordinates -n -1.0 ..+1.0 range
	 */
	public double [] getRelativeCoords(double [] pXY){
		double [] relXY ={
				2.0 * (pXY[0]/this.pixelCorrectionWidth - 0.5),
				2.0 * (pXY[1]/this.pixelCorrectionWidth - 0.5),
		};
		return relXY;
	}

	/**
	 * Calculate pixel coordinates for each of numSensors images, for a given (px,py) of the idealized "center" (still distorted) image
	 * and generic disparity, measured in pixels
	 * @param gc_main - GeometryCorrection instance for the main camera, for which px,py are specified
	 * @param use_rig_offsets - for the auxiliary camera - use offsets from the main one
	 * @param rots misalignment correction (now includes zoom in addition to rotations
	 * @param deriv_rots derivatives by d_az, f_elev, d_rot, d_zoom
	 * @param px pixel X coordinate
	 * @param py pixel Y coordinate
	 * @param disparity disparity
	 * @return array of per port pairs of pixel shifts
	 */
	public double [][] getPortsCoordinatesAndDerivatives(
			GeometryCorrection gc_main,
			boolean     use_rig_offsets,
			Matrix []   rots,
			Matrix [][] deriv_rots,
			double [][] pXYderiv, // if not null, should be double[8][]
			double px,
			double py,
			double disparity)
	{
//		String dbg_s = corr_vector.toString();
/* Starting with required tile center X, Y and nominal distortion, for each sensor port:
 * 1) unapply common distortion (maybe for different - master camera)
 * 2) apply disparity
 * 3) apply rotations and zoom
 * 4) re-apply distortion
 * 5) return port center X and Y
 */
		double [][] rXY = getRXY(use_rig_offsets); // may include rig offsets

		double [][] pXY = new double [numSensors][2];

		double pXcd = px - 0.5 * gc_main.pixelCorrectionWidth;
		double pYcd = py - 0.5 * gc_main.pixelCorrectionHeight;
		double rD = Math.sqrt(pXcd*pXcd + pYcd*pYcd)*0.001*gc_main.pixelSize; // distorted radius in a virtual center camera
		double rND2R=gc_main.getRByRDist(rD/gc_main.distortionRadius, (debugLevel > -1));
		double pXc = pXcd * rND2R; // non-distorted coordinates relative to the (0.5 * this.pixelCorrectionWidth, 0.5 * this.pixelCorrectionHeight)
		double pYc = pYcd * rND2R; // in pixels
		// next radial distortion coefficients are for this, not master camera (may be the same)
		double [] rad_coeff={this.distortionC,this.distortionB,this.distortionA,this.distortionA5,this.distortionA6,this.distortionA7,this.distortionA8};
		double fl_pix = focalLength/(0.001*pixelSize); // focal length in pixels - this camera
		double  ri_scale = 0.001 * this.pixelSize / this.distortionRadius;

		for (int i = 0; i < numSensors; i++){
			// non-distorted XY of the shifted location of the individual sensor
			double pXci0 = pXc - disparity *  rXY[i][0]; // in pixels
			double pYci0 = pYc - disparity *  rXY[i][1];
			// rectilinear, end of dealing with possibly other (master) camera, below all is for this camera distortions


			// Convert a 2-d non-distorted vector to 3d at fl_pix distance in z direction
			double [][] avi = {{pXci0}, {pYci0},{fl_pix}};
			Matrix vi = new Matrix(avi); // non-distorted sensor channel view vector in pixels (z -along the common axis)

			// Apply port-individual combined rotation/zoom matrix

			Matrix rvi = rots[i].times(vi);

			// get back to the projection plane by normalizing vector
			double norm_z = fl_pix/rvi.get(2, 0);
			double pXci =  rvi.get(0, 0) * norm_z;
			double pYci =  rvi.get(1, 0) * norm_z;

			// Re-apply distortion
			double rNDi = Math.sqrt(pXci*pXci + pYci*pYci); // in pixels
			//		Rdist/R=A8*R^7+A7*R^6+A6*R^5+A5*R^4+A*R^3+B*R^2+C*R+(1-A6-A7-A6-A5-A-B-C)");
			double ri = rNDi* ri_scale; // relative to distortion radius
			//    		double rD2rND = (1.0 - distortionA8 - distortionA7 - distortionA6 - distortionA5 - distortionA - distortionB - distortionC);
			double rD2rND = 1.0;
			double rri = 1.0;
			for (int j = 0; j < rad_coeff.length; j++){
				rri *= ri;
				rD2rND += rad_coeff[j]*(rri - 1.0); // Fixed
			}

			// Get port pixel coordiantes by scaling the 2d vector with Rdistorted/Dnondistorted coefficient)
			double pXid = pXci * rD2rND;
			double pYid = pYci * rD2rND;



			pXY[i][0] =  pXid + this.pXY0[i][0];
			pXY[i][1] =  pYid + this.pXY0[i][1];



			if (pXYderiv != null) {
				pXYderiv[2 * i] =   new double [CorrVector.LENGTH];
				pXYderiv[2 * i+1] = new double [CorrVector.LENGTH];
				Matrix drvi_daz = deriv_rots[i][0].times(vi);
				Matrix drvi_dtl = deriv_rots[i][1].times(vi);
				Matrix drvi_drl = deriv_rots[i][2].times(vi);
				Matrix drvi_dzm = deriv_rots[i][3].times(vi);

				double dpXci_dazimuth = drvi_daz.get(0, 0) * norm_z - pXci * drvi_daz.get(2, 0) / rvi.get(2, 0);
				double dpYci_dazimuth = drvi_daz.get(1, 0) * norm_z - pYci * drvi_daz.get(2, 0) / rvi.get(2, 0);

				double dpXci_dtilt =    drvi_dtl.get(0, 0) * norm_z - pXci * drvi_dtl.get(2, 0) / rvi.get(2, 0);
				double dpYci_dtilt =    drvi_dtl.get(1, 0) * norm_z - pYci * drvi_dtl.get(2, 0) / rvi.get(2, 0);

				double dpXci_droll =    drvi_drl.get(0, 0) * norm_z - pXci * drvi_drl.get(2, 0) / rvi.get(2, 0);
				double dpYci_droll =    drvi_drl.get(1, 0) * norm_z - pYci * drvi_drl.get(2, 0) / rvi.get(2, 0);

				double dpXci_dzoom =    drvi_dzm.get(0, 0) * norm_z - pXci * drvi_dzm.get(2, 0) / rvi.get(2, 0);
				double dpYci_dzoom =    drvi_dzm.get(1, 0) * norm_z - pYci * drvi_dzm.get(2, 0) / rvi.get(2, 0);

				double dri_dazimuth =  ri_scale / rNDi* (pXci * dpXci_dazimuth +  pYci * dpYci_dazimuth);
				double dri_dtilt =     ri_scale / rNDi* (pXci * dpXci_dtilt +     pYci * dpYci_dtilt);

				double dri_dzoom =     ri_scale / rNDi* (pXci * dpXci_dzoom +     pYci * dpYci_dzoom);

/*
				double dri_droll =     ri_scale / rNDi* (pXci * dpXci_droll +     pYci * dpYci_droll); // Not used anywhere ?
// TODO: verify dri_droll == 0 and remove
*/
				double drD2rND_dri = 0.0;
				rri = 1.0;
				for (int j = 0; j < rad_coeff.length; j++){
					drD2rND_dri += rad_coeff[j] * (j+1) * rri;
					rri *= ri;
				}

				double drD2rND_dazimuth = drD2rND_dri * dri_dazimuth;
				double drD2rND_dtilt =    drD2rND_dri * dri_dtilt;

				double drD2rND_dzoom =    drD2rND_dri * dri_dzoom; // new

				double dpXid_dazimuth = dpXci_dazimuth * rD2rND + pXci * drD2rND_dazimuth;
				double dpYid_dazimuth = dpYci_dazimuth * rD2rND + pYci * drD2rND_dazimuth;
				double dpXid_dtilt =    dpXci_dtilt    * rD2rND + pXci * drD2rND_dtilt;
				double dpYid_dtilt =    dpYci_dtilt    * rD2rND + pYci * drD2rND_dtilt;


				double dpXid_droll =    dpXci_droll    * rD2rND;
				double dpYid_droll =    dpYci_droll    * rD2rND;

				double dpXid_dzoom =    dpXci_dzoom    * rD2rND + pXci * drD2rND_dzoom; // new second term
				double dpYid_dzoom =    dpYci_dzoom    * rD2rND + pYci * drD2rND_dzoom; // new second term


				// verify that d/dsym are well, symmetrical
				if (i < (numSensors - 1)){
					pXYderiv[2 * i + 0][CorrVector.TILT_INDEX+i] =         dpXid_dtilt;
					pXYderiv[2 * i + 1][CorrVector.TILT_INDEX+i] =         dpYid_dtilt;
					pXYderiv[2 * i + 0][CorrVector.AZIMUTH_INDEX+i] =      dpXid_dazimuth;
					pXYderiv[2 * i + 1][CorrVector.AZIMUTH_INDEX+i] =      dpYid_dazimuth;

					pXYderiv[2 * i + 0][CorrVector.ZOOM_INDEX+i] =         dpXid_dzoom;
					pXYderiv[2 * i + 1][CorrVector.ZOOM_INDEX+i] =         dpYid_dzoom;
				} else {
					for (int j = 0; j < (numSensors - 1); j++){
						pXYderiv[2 * i + 0][CorrVector.TILT_INDEX+j] =    -dpXid_dtilt;
						pXYderiv[2 * i + 1][CorrVector.TILT_INDEX+j] =    -dpYid_dtilt;
						pXYderiv[2 * i + 0][CorrVector.AZIMUTH_INDEX+j] = -dpXid_dazimuth;
						pXYderiv[2 * i + 1][CorrVector.AZIMUTH_INDEX+j] = -dpYid_dazimuth;

						pXYderiv[2 * i + 0][CorrVector.ZOOM_INDEX+j] =    -dpXid_dzoom;
						pXYderiv[2 * i + 1][CorrVector.ZOOM_INDEX+j] =    -dpYid_dzoom;
					}
				}
				pXYderiv[2 * i + 0][CorrVector.ROLL_INDEX+i] = dpXid_droll;
				pXYderiv[2 * i + 1][CorrVector.ROLL_INDEX+i] = dpYid_droll;
			}
		}
		return pXY;
	}

	/**
	 * Calculate pixel common "idealized" coordinates of the auxiliary camera image tile matching  the specified tile (px,py) of the idealized
	 * "center" (still distorted) image of the main camera and generic disparity, measured in pixels
	 * @param gc_main - GeometryCorrection instance for the main camera, for which px,py are specified
	 * @param aux_rot  misalignment correction matrix // getAuxRotMatrix()
	 * @param aux_rot_derivs derivatives by d_az, f_elev, d_rot, d_zoom // getAuxRotDeriveMatrices()
	 * @param aux_offset_derivs aux camera x,y offset per one pixel disparity of the main camera and derivatives by angle and baseline
	 * @param pXYderivs x and y (aux camera idealized coordinates) derivatives by azimuth, tilt, roll, zoom, angle and baseline
	 * If it is not null it will be used to return derivatives. should be initialized as double[2][]
	 * @param px pixel X coordinate (master camera idealized)
	 * @param py pixel Y coordinate (master camera idealized)
	 * @param disparity disparity  (master camera pixels)
	 * @return a pair of x,y coordinates of the matching image tile of the aux camera image
	 */

	public double [] getRigAuxCoordinatesAndDerivatives(
			GeometryCorrection gc_main,
			Matrix      aux_rot,
			Matrix []   aux_rot_derivs,
			double [][] aux_offset_derivs,
			double [][] pXYderiv, // if not null, should be double[2][]
			double px,
			double py,
			double disparity)
	{
/* Starting with required tile center X, Y and nominal distortion, for each sensor port:
 * 1) unapply common distortion of the main camera
 * 2) apply disparity (shift)
 * 3) apply rotations and zoom
 * 4) re-apply distortion (of the main camera!)
 * 5) return aux center X and Y in main camera coordinates
 */

		double pXcd = px - 0.5 * gc_main.pixelCorrectionWidth;
		double pYcd = py - 0.5 * gc_main.pixelCorrectionHeight;
		double rD = Math.sqrt(pXcd*pXcd + pYcd*pYcd)*0.001*gc_main.pixelSize; // distorted radius in a virtual center camera
		double rND2R=gc_main.getRByRDist(rD/gc_main.distortionRadius, (debugLevel > -1));
		double pXc = pXcd * rND2R; // non-distorted coordinates relative to the (0.5 * this.pixelCorrectionWidth, 0.5 * this.pixelCorrectionHeight)
		double pYc = pYcd * rND2R; // in pixels

		// for re-applying distortion
		double [] rad_coeff={gc_main.distortionC,gc_main.distortionB,gc_main.distortionA,gc_main.distortionA5,gc_main.distortionA6,gc_main.distortionA7,gc_main.distortionA8};
		double fl_pix = gc_main.focalLength/(0.001*gc_main.pixelSize); // focal length in pixels - this camera
		double  ri_scale = 0.001 * gc_main.pixelSize / gc_main.distortionRadius;

		// non-distorted XY relative to the auxiliary camera center if it was parallel to the main one
		// Allow aux_offset_derivs null only if disparity == 0.0;
		double pXci0 = pXc; //  - disparity *  aux_offset_derivs[0][0]; // in pixels
		double pYci0 = pYc; // - disparity *  aux_offset_derivs[0][1]; // in pixels
		if (disparity != 0.0) { // aux_offset_derivs != null)
			pXci0 -= disparity *  aux_offset_derivs[0][0]; // in pixels
			pYci0 -= disparity *  aux_offset_derivs[0][1]; // in pixels
		}
		// rectilinear here

		// Convert a 2-d non-distorted vector to 3d at fl_pix distance in z direction
		double [][] avi = {{pXci0}, {pYci0},{fl_pix}};
		Matrix vi = new Matrix(avi); // non-distorted sensor channel view vector in pixels (z -along the common axis)

		// Apply port-individual combined rotation/zoom matrix

		Matrix rvi = (aux_rot == null) ? vi: aux_rot.times(vi);


		// get back to the projection plane by normalizing vector
		double norm_z = fl_pix/rvi.get(2, 0);
		double pXci =  rvi.get(0, 0) * norm_z;
		double pYci =  rvi.get(1, 0) * norm_z;

		// non-distorted XY relative to the auxiliary camera center in the auxiliary camera coordinate system
		// if converted back to distorted main camera it should match original px, py

		// Re-apply distortion
		double rNDi = Math.sqrt(pXci*pXci + pYci*pYci); // in pixels
		//		Rdist/R=A8*R^7+A7*R^6+A6*R^5+A5*R^4+A*R^3+B*R^2+C*R+(1-A6-A7-A6-A5-A-B-C)");
		double ri = rNDi* ri_scale; // relative to distortion radius
		//    		double rD2rND = (1.0 - distortionA8 - distortionA7 - distortionA6 - distortionA5 - distortionA - distortionB - distortionC);
		double rD2rND = 1.0;
		double rri = 1.0;
		for (int j = 0; j < rad_coeff.length; j++){
			rri *= ri;
			rD2rND += rad_coeff[j]*(rri - 1.0); // Fixed
		}

		// Get aux camera center in master camera pixel coordinates by scaling the 2d vector with Rdistorted/Dnondistorted coefficient)
		double pXid = pXci * rD2rND;
		double pYid = pYci * rD2rND;

        double [] pXY = {pXid, pYid};


		if (pXYderiv != null) {
			pXYderiv[0] = new double [RigOffset.VECTOR_LENGTH];
			pXYderiv[1] = new double [RigOffset.VECTOR_LENGTH];

			double [][] advi_dangle =   {{ -disparity *  aux_offset_derivs[1][0]}, {-disparity *  aux_offset_derivs[1][1]},{0}};
			double [][] advi_baseline = {{ -disparity *  aux_offset_derivs[2][0]}, {-disparity *  aux_offset_derivs[2][1]},{0}};
			Matrix dvi_dangle =    new Matrix(advi_dangle);   // non-distorted pXci, pYci derivatives by angle of the aux camera position in the main camera principal plane
			Matrix dvi_dbaseline = new Matrix(advi_baseline); // non-distorted pXci, pYci derivatives by the rig baseline

			Matrix drvi_dangle =    aux_rot.times(dvi_dangle);
			Matrix drvi_dbaseline = aux_rot.times(dvi_dbaseline);

			Matrix drvi_daz = aux_rot_derivs[0].times(vi);
			Matrix drvi_dtl = aux_rot_derivs[1].times(vi);
			Matrix drvi_drl = aux_rot_derivs[2].times(vi);
			Matrix drvi_dzm = aux_rot_derivs[3].times(vi);

			double dpXci_dangle = drvi_dangle.get(0, 0) * norm_z - pXci * drvi_dangle.get(2, 0) / rvi.get(2, 0);
			double dpYci_dangle = drvi_dangle.get(1, 0) * norm_z - pYci * drvi_dangle.get(2, 0) / rvi.get(2, 0);

			double dpXci_dbaseline = drvi_dbaseline.get(0, 0) * norm_z - pXci * drvi_dbaseline.get(2, 0) / rvi.get(2, 0);
			double dpYci_dbaseline = drvi_dbaseline.get(1, 0) * norm_z - pYci * drvi_dbaseline.get(2, 0) / rvi.get(2, 0);

			double dpXci_dazimuth = drvi_daz.get(0, 0) * norm_z - pXci * drvi_daz.get(2, 0) / rvi.get(2, 0);
			double dpYci_dazimuth = drvi_daz.get(1, 0) * norm_z - pYci * drvi_daz.get(2, 0) / rvi.get(2, 0);

			double dpXci_dtilt =    drvi_dtl.get(0, 0) * norm_z - pXci * drvi_dtl.get(2, 0) / rvi.get(2, 0);
			double dpYci_dtilt =    drvi_dtl.get(1, 0) * norm_z - pYci * drvi_dtl.get(2, 0) / rvi.get(2, 0);

			double dpXci_droll =    drvi_drl.get(0, 0) * norm_z - pXci * drvi_drl.get(2, 0) / rvi.get(2, 0);
			double dpYci_droll =    drvi_drl.get(1, 0) * norm_z - pYci * drvi_drl.get(2, 0) / rvi.get(2, 0);

			double dpXci_dzoom =    drvi_dzm.get(0, 0) * norm_z - pXci * drvi_dzm.get(2, 0) / rvi.get(2, 0);
			double dpYci_dzoom =    drvi_dzm.get(1, 0) * norm_z - pYci * drvi_dzm.get(2, 0) / rvi.get(2, 0);

			double dri_dangle =    ri_scale / rNDi* (pXci * dpXci_dangle +    pYci * dpYci_dangle);
			double dri_dbaseline = ri_scale / rNDi* (pXci * dpXci_dbaseline + pYci * dpYci_dbaseline);
			double dri_dazimuth =  ri_scale / rNDi* (pXci * dpXci_dazimuth +  pYci * dpYci_dazimuth);
			double dri_dtilt =     ri_scale / rNDi* (pXci * dpXci_dtilt +     pYci * dpYci_dtilt);
			double dri_droll =     ri_scale / rNDi* (pXci * dpXci_droll +     pYci * dpYci_droll); // Not used anywhere ?
			double dri_dzoom =     ri_scale / rNDi* (pXci * dpXci_dzoom +     pYci * dpYci_dzoom);

			// TODO: verify dri_droll == 0 and remove

			double drD2rND_dri = 0.0;
			rri = 1.0;
			for (int j = 0; j < rad_coeff.length; j++){
				drD2rND_dri += rad_coeff[j] * (j+1) * rri;
				rri *= ri;
			}

			double drD2rND_dangle =    drD2rND_dri * dri_dangle;
			double drD2rND_dbaseline = drD2rND_dri * dri_dbaseline;

			double drD2rND_dazimuth = drD2rND_dri * dri_dazimuth;
			double drD2rND_dtilt =    drD2rND_dri * dri_dtilt;
			double drD2rND_droll =    drD2rND_dri * dri_droll; // new
			double drD2rND_dzoom =    drD2rND_dri * dri_dzoom; // new


			double dpXid_dangle =     dpXci_dangle * rD2rND + pXci * drD2rND_dangle;
			double dpYid_dangle =     dpYci_dangle * rD2rND + pYci * drD2rND_dangle;
			double dpXid_dbaseline =  dpXci_dbaseline * rD2rND + pXci * drD2rND_dbaseline;
			double dpYid_dbaseline =  dpYci_dbaseline * rD2rND + pYci * drD2rND_dbaseline;

			double dpXid_dazimuth =   dpXci_dazimuth * rD2rND + pXci * drD2rND_dazimuth;
			double dpYid_dazimuth =   dpYci_dazimuth * rD2rND + pYci * drD2rND_dazimuth;
			double dpXid_dtilt =      dpXci_dtilt    * rD2rND + pXci * drD2rND_dtilt;
			double dpYid_dtilt =      dpYci_dtilt    * rD2rND + pYci * drD2rND_dtilt;

			double dpXid_droll =      dpXci_droll    * rD2rND + pXci * drD2rND_droll; // probably should be zero anyway
			double dpYid_droll =      dpYci_droll    * rD2rND + pYci * drD2rND_droll; // probably should be zero anyway

			double dpXid_dzoom =      dpXci_dzoom    * rD2rND + pXci * drD2rND_dzoom; // new second term
			double dpYid_dzoom =      dpYci_dzoom    * rD2rND + pYci * drD2rND_dzoom; // new second term

			pXYderiv[0][RigOffset.AUX_AZIMUTH_INDEX] = dpXid_dazimuth;
			pXYderiv[1][RigOffset.AUX_AZIMUTH_INDEX] = dpYid_dazimuth;

			pXYderiv[0][RigOffset.AUX_TILT_INDEX] = dpXid_dtilt;
			pXYderiv[1][RigOffset.AUX_TILT_INDEX] = dpYid_dtilt;

			pXYderiv[0][RigOffset.AUX_ROLL_INDEX] = dpXid_droll;
			pXYderiv[1][RigOffset.AUX_ROLL_INDEX] = dpYid_droll;

			pXYderiv[0][RigOffset.AUX_ZOOM_INDEX] = dpXid_dzoom;
			pXYderiv[1][RigOffset.AUX_ZOOM_INDEX] = dpYid_dzoom;

			pXYderiv[0][RigOffset.AUX_ANGLE_INDEX] = dpXid_dangle;
			pXYderiv[1][RigOffset.AUX_ANGLE_INDEX] = dpYid_dangle;

			pXYderiv[0][RigOffset.AUX_BASELINE_INDEX] = dpXid_dbaseline;
			pXYderiv[1][RigOffset.AUX_BASELINE_INDEX] = dpYid_dbaseline;
		}
		return pXY;
	}

	public double [][] getPortsCoordinatesAndDerivatives( // uses rotations - used in AlignmentCorrection class
			boolean use_rig_offsets,
			double [] dbg_a_vector, // replace actual radial distortion coefficients (not currently used)
			double delta, // 1e-6
			CorrVector corr_vector,
			double [][] pXYderiv, // if not null, should be double[8][]
			double px,
			double py,
			double disparity)
	{
		// slower, will re-calculate matrices for each tile, but for debug - that is OK
		Matrix []   corr_rots =  corr_vector.getRotMatrices(); // get array of per-sensor rotation matrices
		double [][] rslt = getPortsCoordinatesAndDerivatives(
				this,
				use_rig_offsets,
				corr_rots, // Matrix []   rots,
				null, // deriv_rots, // Matrix [][] deriv_rots,
				null, // pXYderiv0, // null, // false, // boolean calc_deriv,
				px, // double px,
				py, // double py,
				disparity // double disparity
				);

		for (int i = 0; i < pXYderiv.length; i++){
			pXYderiv[i] = new double [CorrVector.LENGTH];
		}
		for (int nPar = 0; nPar < CorrVector.LENGTH; nPar++){
			CorrVector cv_delta_p = corr_vector.clone();
			CorrVector cv_delta_m = corr_vector.clone();
			cv_delta_p.vector[nPar] += 0.5 *delta;
			cv_delta_m.vector[nPar] -= 0.5 *delta;
			Matrix []   corr_rots_p =  cv_delta_p.getRotMatrices(); // get array of per-sensor rotation matrices
			Matrix []   corr_rots_m =  cv_delta_m.getRotMatrices(); // get array of per-sensor rotation matrices

			double [][] rslt_p = getPortsCoordinatesAndDerivatives(
					this,
					use_rig_offsets,
					corr_rots_p, // Matrix []   rots,
					null, // Matrix [][] deriv_rots,
					null, // boolean calc_deriv,
					px, // double px,
					py, // double py,
					disparity // double disparity
					);
			double [][] rslt_m = getPortsCoordinatesAndDerivatives(
					this,
					use_rig_offsets,
					corr_rots_m, // Matrix []   rots,
					null, // Matrix [][] deriv_rots,
					null, // boolean calc_deriv,
					px, // double px,
					py, // double py,
					disparity // double disparity
					);
			for (int sens = 0; sens < numSensors; sens++){
				for (int dir = 0; dir < 2; dir++){
					pXYderiv[2 * sens +dir][nPar] = (rslt_p[sens][dir] - rslt_m[sens][dir]) / delta;
				}
			}
		}
		return rslt;
	}


	// should return same as input if disparity==0
	public double [][] getPortsCoordinatesIdeal( // used in macro mode
			double px,
			double py,
			double disparity)
	{
		// reverse getPortsCoordinates
		double c_roll = 1.0; // Math.cos(( - this.common_roll) * Math.PI/180.0);
		double s_roll = 0.0; // Math.sin(( - this.common_roll) * Math.PI/180.0);
		double pXcd0 = px -  0.5 * this.pixelCorrectionWidth;
		double pYcd0 = py -  0.5 * this.pixelCorrectionWidth;
		double pXcd = c_roll *  pXcd0 - s_roll* pYcd0;
		double pYcd = s_roll *  pXcd0 + c_roll* pYcd0;
		double rD = Math.sqrt(pXcd*pXcd + pYcd*pYcd)*0.001*this.pixelSize; // distorted radius in a virtual center camera
		double rND2R=getRByRDist(rD/this.distortionRadius, (debugLevel > -1));
		double pXc = pXcd * rND2R; // non-distorted coordinates relative to the (0.5 * this.pixelCorrectionWidth, 0.5 * this.pixelCorrectionHeight)
		double pYc = pYcd * rND2R; // in pixels
		double [] a={this.distortionC,this.distortionB,this.distortionA,this.distortionA5,this.distortionA6,this.distortionA7,this.distortionA8};
		double [][] pXY = new double [numSensors][2];
		for (int i = 0; i < numSensors; i++){
			// non-distorted XY of the shifted location of the individual sensor
			double pXci = pXc - disparity *  this.rXY_ideal[i][0]; // in pixels
			double pYci = pYc - disparity *  this.rXY_ideal[i][1];

			// calculate back to distorted
			double rNDi = Math.sqrt(pXci*pXci + pYci*pYci); // in pixels
			//		Rdist/R=A8*R^7+A7*R^6+A6*R^5+A5*R^4+A*R^3+B*R^2+C*R+(1-A8-A7-A6-A5-A-B-C)");
			double ri = rNDi* 0.001 * this.pixelSize / this.distortionRadius; // relative to distortion radius
			//    		double rD2rND = (1.0 - distortionA8 - distortionA7 - distortionA6 - distortionA5 - distortionA - distortionB - distortionC);
			double rD2rND = 1.0;
			double rri = 1.0;
			for (int j = 0; j < a.length; j++){
				rri *= ri;
				rD2rND += a[j]*(rri - 1.0);
			}
			double pXid = pXci * rD2rND;
			double pYid = pYci * rD2rND;
			pXY[i][0] =  c_roll *  pXid + s_roll* pYid + 0.5 * this.pixelCorrectionWidth; // this.pXY0[i][0];
			pXY[i][1] = -s_roll *  pXid + c_roll* pYid + 0.5 * this.pixelCorrectionWidth; // this.pXY0[i][1];
		}
		return pXY;
	}

	/*
	public double [] getAuxCoordinatesRigIdeal( // used in macro mode
			GeometryCorrection gc_main,
			Matrix rots,
			double px,
			double py,
			double disparity)
	{
		// reverse getPortsCoordinates
		double c_roll = 1.0; // Math.cos(( - this.common_roll) * Math.PI/180.0);
		double s_roll = 0.0; // Math.sin(( - this.common_roll) * Math.PI/180.0);
		double pXcd0 = px -  0.5 * this.pixelCorrectionWidth;
		double pYcd0 = py -  0.5 * this.pixelCorrectionWidth;
		double pXcd = c_roll *  pXcd0 - s_roll* pYcd0;
		double pYcd = s_roll *  pXcd0 + c_roll* pYcd0;
		double rD = Math.sqrt(pXcd*pXcd + pYcd*pYcd)*0.001*this.pixelSize; // distorted radius in a virtual center camera
		double rND2R=getRByRDist(rD/this.distortionRadius, (debugLevel > -1));
		double pXc = pXcd * rND2R; // non-distorted coordinates relative to the (0.5 * this.pixelCorrectionWidth, 0.5 * this.pixelCorrectionHeight)
		double pYc = pYcd * rND2R; // in pixels
		double [] a={this.distortionC,this.distortionB,this.distortionA,this.distortionA5,this.distortionA6,this.distortionA7,this.distortionA8};
		double [] pXY = new double[2];
		// calculate for aux (this) camera
		double [][] aux_offset = getAuxOffsetAndDerivatives(gc_main);
			// non-distorted XY of the shifted location of the individual sensor
			double pXci = pXc - disparity *  aux_offset[0][0]; // in pixels
			double pYci = pYc - disparity *  aux_offset[0][1];

			// calculate back to distorted
			double rNDi = Math.sqrt(pXci*pXci + pYci*pYci); // in pixels
			//		Rdist/R=A8*R^7+A7*R^6+A6*R^5+A5*R^4+A*R^3+B*R^2+C*R+(1-A8-A7-A6-A5-A-B-C)");
			double ri = rNDi* 0.001 * this.pixelSize / this.distortionRadius; // relative to distortion radius
			//    		double rD2rND = (1.0 - distortionA8 - distortionA7 - distortionA6 - distortionA5 - distortionA - distortionB - distortionC);
			double rD2rND = 1.0;
			double rri = 1.0;
			for (int j = 0; j < a.length; j++){
				rri *= ri;
				rD2rND += a[j]*(rri - 1.0);
			}
			double pXid = pXci * rD2rND;
			double pYid = pYci * rD2rND;
			pXY[i][0] =  c_roll *  pXid + s_roll* pYid + 0.5 * this.pixelCorrectionWidth; // this.pXY0[i][0];
			pXY[i][1] = -s_roll *  pXid + c_roll* pYid + 0.5 * this.pixelCorrectionWidth; // this.pXY0[i][1];
		return pXY;
	}
*/


	public double [][] getPortsCoordinatesIdeal(
			int    macro_scale, // 1 for pixels, 8 - for tiles when correlating tiles instead of the pixels
			double px,
			double py,
			double disparity)
	{
		double [][] coords = getPortsCoordinatesIdeal(
				px * macro_scale,
				py * macro_scale,
				disparity);
		for (int i = 0; i < coords.length; i++){
			for (int j = 0; j< coords[i].length; j++){
				coords[i][j] /= macro_scale;
			}
		}
		return coords;
	}

	public double [] getRigAuxCoordinatesIdeal(
			int                macro_scale, // 1 for pixels, 8 - for tiles when correlating tiles instead of the pixels
			GeometryCorrection gc_main,
			Matrix             aux_rot,
			double             px,
			double             py,
			double             disparity)
	{
		double [] xy = getRigAuxCoordinatesAndDerivatives(
				gc_main,           // GeometryCorrection gc_main,
				aux_rot,           // Matrix      aux_rot,
				null,              // Matrix []   aux_rot_derivs,
				null,              // double [][] aux_offset_derivs,
				null,              // double [][] pXYderiv, // if not null, should be double[6][]
				px * macro_scale,  // double px,
				py * macro_scale,  // double py,
				disparity);        // double disparity);
		double [] coords = {xy[0]/macro_scale,xy[1]/macro_scale};
		return coords;
	}





	// Copied from PixelMapping
	/**
	 * Calculate reverse distortion table - from pixel radius to non-distorted radius
	 * Rdist/R=A5*R^4+A*R^3+B*R^2+C*R+(1-A5-A-B-C)
	 * @return false if distortion is too high
	 */
	public boolean calcReverseDistortionTable(){
		boolean debugThis=false; //true;
		double delta=1E-8;
		double minDerivative=0.1;
		int numIterations=1000;
		double drDistDr=1.0;
		//	public double distortionA5=0.0; //r^5 (normalized to focal length or to sensor half width?)
		//	public double distortionA=0.0; // r^4 (normalized to focal length or to sensor half width?)
		//	public double distortionB=0.0; // r^3
		//	public double distortionC=0.0; // r^2
		boolean use8=(this.distortionA8!=0.0) || (this.distortionA7!=0.0) || (this.distortionA6!=0.0);
		double d=1.0-this.distortionA8-this.distortionA7-this.distortionA6-this.distortionA5-this.distortionA-this.distortionB-this.distortionC;
		double rPrev=0.0;
		this.rByRDist=new double [(int) Math.ceil(this.maxR/this.stepR)+1];
		for (int j=1;j<this.rByRDist.length;j++) this.rByRDist[j]=Double.NaN;
		this.rByRDist[0]=1.0/d;
		boolean bailOut=false;
		if (debugThis)	System.out.println("calcReverseDistortionTable()");

		for (int i=1;i<this.rByRDist.length;i++) {
			double rDist=this.stepR*i;
			double r=rPrev+this.stepR/drDistDr;
			//		if (debugThis)	System.out.println("calcReverseDistortionTable() i="+i+" rDist="+rDist+" r="+r+" rPrev="+rPrev);

			for (int iteration=0;iteration<numIterations;iteration++){
				double k;
				if (use8){
					k=(((((((this.distortionA8)*r+this.distortionA7)*r+this.distortionA6)*r+this.distortionA5)*r + this.distortionA)*r+this.distortionB)*r+this.distortionC)*r+d;
					drDistDr=(((((((8*this.distortionA8)*r + 7*this.distortionA7)*r + 6*this.distortionA6)*r + 5*this.distortionA5)*r + 4*this.distortionA)*r+3*this.distortionB)*r+2*this.distortionC)*r+d;
				} else {
					k=(((this.distortionA5*r + this.distortionA)*r+this.distortionB)*r+this.distortionC)*r+d;
					drDistDr=(((5*this.distortionA5*r + 4*this.distortionA)*r+3*this.distortionB)*r+2*this.distortionC)*r+d;
				}
				double rD=r*k;
				if (drDistDr<minDerivative) {
					bailOut=true;
					break; // too high distortion
				}
				if (Math.abs(rD-rDist)<delta) break; // success
				r+=(rDist-rD)/drDistDr;
			}
			if (bailOut) {
				if (debugThis)	System.out.println("calcReverseDistortionTable() i="+i+" Bailing out, drDistDr="+drDistDr);
				return false;
			}
			rPrev=r;
			this.rByRDist[i]=r/rDist;
			if (debugThis)	System.out.println("calcReverseDistortionTable() i="+i+" rDist="+rDist+" r="+r+" rPrev="+rPrev+" this.rByRDist[i]="+this.rByRDist[i]);
		}
		return true;
	}

	/**
	 * Get relative (to distortion radius) distorted radius ratio to non-distorted from non-distorted
	 * @param r non-distorted radius (1.0 is 2.8512mm)
	 * @return ratio of distorted to non-distorted radius
	 */
	public double getRDistByR(double r) // relative to distortion radius
	{
		boolean use8=(this.distortionA8!=0.0) || (this.distortionA7!=0.0) || (this.distortionA6!=0.0);
		double d=1.0-this.distortionA8-this.distortionA7-this.distortionA6-this.distortionA5-this.distortionA-this.distortionB-this.distortionC;

		double k;
		if (use8){
			k=(((((((this.distortionA8)*r+this.distortionA7)*r+this.distortionA6)*r+this.distortionA5)*r + this.distortionA)*r+this.distortionB)*r+this.distortionC)*r+d;
		} else {
			k=(((this.distortionA5*r + this.distortionA)*r+this.distortionB)*r+this.distortionC)*r+d;
		}
		return k;
	}

	/**
	 * Get derivative of relative (to distortion radius) d_Rdist/d_R from relative (to distortion radius) non-distorted radius
	 * @param r non-distorted relative radius
	 * @return derivative d_Rdist/d_R from (relative to relative)
	 */
	public double getDerivRDistFromR(double r) // relative to distortion radius
	{
		boolean use8=(this.distortionA8!=0.0) || (this.distortionA7!=0.0) || (this.distortionA6!=0.0);
		double drDistDr;
		if (use8){
			drDistDr=(((((((7*this.distortionA8)*r + 6*this.distortionA7)*r + 5*this.distortionA6)*r + 4*this.distortionA5)*r + 3*this.distortionA)*r+2*this.distortionB)*r+1*this.distortionC); // +d;
		} else {
			drDistDr=((4*this.distortionA5*r + 3*this.distortionA)*r+2*this.distortionB)*r+1*this.distortionC;
		}
		return drDistDr;
	}

	public double getDerivRDistFromR(double r, double delta) // relative to distortion radius
	{
		return (getRDistByR(r+delta) -getRDistByR(r))/delta;

	}


	public double getRByRDist(double rDist, boolean debug){
		// add exceptions;
		if (this.rByRDist==null) {
			calcReverseDistortionTable();
			if (debug)System.out.println("getRByRDist("+IJ.d2s(rDist,3)+"): this.rByRDist==null");
			//		return Double.NaN;
		}
		if (rDist<0) {
			if (debug)System.out.println("getRByRDist("+IJ.d2s(rDist,3)+"): rDist<0");
			return Double.NaN;
		}
		int index=(int) Math.floor(rDist/this.stepR);
		if (index>=(this.rByRDist.length-1)) {
			if (debug) System.out.println("getRByRDist("+IJ.d2s(rDist,3)+"): index="+index+">="+(this.rByRDist.length-1));
			return Double.NaN;
		}
		double result=this.rByRDist[index]+(this.rByRDist[index+1]-this.rByRDist[index])*(rDist/this.stepR-index);
		if (Double.isNaN(result)){
			if (debug) System.out.println("this.rByRDist["+index+"]="+this.rByRDist[index]);
			if (debug) System.out.println("this.rByRDist["+(index+1)+"]="+this.rByRDist[index+1]);
			if (debug) System.out.println("rDist="+rDist);
			if (debug) System.out.println("(rDist/this.stepR="+(rDist/this.stepR));

		}
		return result;
	}

}
