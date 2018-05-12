import java.util.Properties;

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
	public static String RIG_PREFIX =     "rig-";
	static double SCENE_UNITS_SCALE = 0.001;
	static String SCENE_UNITS_NAME = "m";
	static final String [] CORR_NAMES = {"tilt0","tilt1","tilt2","azimuth0","azimuth1","azimuth2","roll0","roll1","roll2","roll3","zoom0","zoom1","zoom2"};

	public int    debugLevel = 0;
	public int    pixelCorrectionWidth=2592;   // virtual camera center is at (pixelCorrectionWidth/2, pixelCorrectionHeight/2)
	public int    pixelCorrectionHeight=1936;
	public double focalLength=4.5;
	public double pixelSize=  2.2; //um
	public double distortionRadius=  2.8512; // mm - half width of the sensor
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

	public Matrix getRotMatrix(boolean use_rig){
		return (use_rig && (rigOffset != null)) ? rigOffset.getRotMatrix(): null ;
	}

	public double getDisparityRadius() {
		return disparityRadius;
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
		static final double ROT_AZ_SGN = -1.0; // sign of first sin for azimuth rotation
		static final double ROT_TL_SGN =  1.0; // sign of first sin for tilt rotation
		static final double ROT_RL_SGN =  1.0; // sign of first sin for roll rotation
		public double baseline =   1256.0; // mm, distance between camera centers
		public double aux_angle =     0.0; // radians, 0 - aux camera to the right of the main, pi/2 - up
		// consider aux_z small enough for now, will need for SfM
		public double aux_z     =     0.0; // mm auxiliary camera distance from the plane of the main one (positive towards the scene)
		public double aux_azimuth =   0.0; // radians, azimuth of the auxiliary camera (positive - looks to the right)
		public double aux_tilt =      0.0; // radians, tilt of the auxiliary camera (positive - looks up)
		public double aux_roll =      0.0; // radians, roll of the auxiliary camera (positive - looks clockwise)
		public double aux_zoom =      0.0; // relative global zoom of the aux camera relative to the main one, difference from 1.0
		public double [][] rXY_aux = null; // XY pairs of the in a normal plane, relative to disparityRadius

/*
	private double [][] rXY =     null; // XY pairs of the in a normal plane, relative to disparityRadius
	private double [][] rXY_ideal = {{-0.5, -0.5}, {0.5,-0.5}, {-0.5, 0.5}, {0.5,0.5}};

 */
		public RigOffset () {
			System.out.println("created RigOffset");
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
				if (debugLevel > -2) {
					System.out.println("Auxiliary camera offsets per 1 nominal disparity pixel");
					for (int i = 0; i <rXY_aux.length;i++) {
						System.out.println(String.format("Camera %1d x = %8f y = %8f",i,rXY_aux[i][0],rXY_aux[i][1]));
					}
				}
			}
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
						{ ca,               0.0,  sa * ROT_AZ_SGN },
						{ 0.0,              1.0,  0.0},
						{ -sa* ROT_AZ_SGN,  0.0,  ca}};

				double [][] a_t =  { // inverted - OK
						{ 1.0,  0.0,              0.0},
						{ 0.0,  ct,               st * ROT_TL_SGN},
						{ 0.0, -st * ROT_TL_SGN,  ct}};

				double [][] a_r =  { // inverted OK
						{ cr,                sr * ROT_RL_SGN,  0.0},
						{ -sr * ROT_RL_SGN,  cr,               0.0},
						{ 0.0,               0.0,              1.0}};

				double [][] a_daz = { // inverted - OK
						{ -sa,              0.0,   ca * ROT_AZ_SGN },
						{  0.0,             0.0,   0.0},
						{ -ca* ROT_AZ_SGN,  0.0,  -sa}};

				double [][] a_dt =  { // inverted - OK
						{ 0.0,  0.0,               0.0},
						{ 0.0, -st,                ct * ROT_TL_SGN},
						{ 0.0, -ct * ROT_TL_SGN,  -st}};

				double [][] a_dr =  { // inverted OK
						{ -sr * zoom,               cr * zoom * ROT_RL_SGN,  0.0},
						{ -cr * zoom *ROT_RL_SGN,  -sr * zoom,               0.0},
						{ 0.0,               0.0,              0.0}};

				double [][] a_dzoom =  { // inverted OK
						{ cr,                sr * ROT_RL_SGN,  0.0},
						{ -sr * ROT_RL_SGN,  cr,               0.0},
						{ 0.0,               0.0,              0.0}};


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
			gd.addNumericField("Auxilliary camera azimuth  (positive - to the right)",                1000.0*focalLength/pixelSize * this.aux_azimuth,  3,6,"pix",
					"Relative to the main camera axis, shift of the center of the image in pixels");
			gd.addNumericField("Auxilliary camera tilt (positive - looking up)",                      1000.0*focalLength/pixelSize * this.aux_tilt,  3,6,"pix",
					"Relative to the main camera, shift of the center of the image in pixels");
			gd.addNumericField("Auxilliary camera roll (positive - clockwise)",                       1000.0*distortionRadius/pixelSize * this.aux_roll,  3,6,"pix",
					"Roll of a camera as a whole relative to the main camera, shift at the image half-width from the center");
			gd.addNumericField("Relative zoom - difference from 1.0 in parts parts per 1/1000",       1000.0*distortionRadius/pixelSize * this.aux_zoom,  3,6,"pix",
					"Zoom ratio, shift at the image half-width from the center");
  			gd.showDialog();
			if (gd.wasCanceled()) return false;
			this.baseline=     gd.getNextNumber();
			this.aux_angle=    gd.getNextNumber() * Math.PI/180;
			this.aux_z=        gd.getNextNumber();
			this.aux_azimuth=  gd.getNextNumber()/(1000.0*focalLength/pixelSize) ;
			this.aux_tilt=     gd.getNextNumber()/(1000.0*focalLength/pixelSize);
			this.aux_roll=     gd.getNextNumber()/(1000.0*distortionRadius/pixelSize);
			this.aux_zoom=     gd.getNextNumber()/(1000.0*distortionRadius/pixelSize);
			recalcRXY();
			return true;
		}
	}

	public boolean editRig() {
		if (this.rigOffset == null) {
			this.rigOffset = new RigOffset();
		}
		return this.rigOffset.editOffsetsPixels();
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
					rots[chn] = rigMatrix.times(rots[chn]);
				}
			}
			return rots;
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
						{ ca,               0.0,  sa * ROT_AZ_SGN },
						{ 0.0,              1.0,  0.0},
						{ -sa* ROT_AZ_SGN,  0.0,  ca}};

				double [][] a_t =  { // inverted - OK
						{ 1.0,  0.0,              0.0},
						{ 0.0,  ct,               st * ROT_TL_SGN},
						{ 0.0, -st * ROT_TL_SGN,  ct}};

				double [][] a_r =  { // inverted OK
						{ cr,                sr * ROT_RL_SGN,  0.0},
						{ -sr * ROT_RL_SGN,  cr,               0.0},
						{ 0.0,               0.0,              1.0}};

				double [][] a_daz = { // inverted - OK
						{ -sa,              0.0,   ca * ROT_AZ_SGN },
						{  0.0,             0.0,   0.0},
						{ -ca* ROT_AZ_SGN,  0.0,  -sa}};

				double [][] a_dt =  { // inverted - OK
						{ 0.0,  0.0,               0.0},
						{ 0.0, -st,                ct * ROT_TL_SGN},
						{ 0.0, -ct * ROT_TL_SGN,  -st}};

				double [][] a_dr =  { // inverted OK
						{ -sr * zoom,               cr * zoom * ROT_RL_SGN,  0.0},
						{ -cr * zoom *ROT_RL_SGN,  -sr * zoom,               0.0},
						{ 0.0,               0.0,              0.0}};

				double [][] a_dzoom =  { // inverted OK
						{ cr,                sr * ROT_RL_SGN,  0.0},
						{ -sr * ROT_RL_SGN,  cr,               0.0},
						{ 0.0,               0.0,              0.0}};


				// d/d_az
				rot_derivs[chn][0] = (new Matrix(a_r ).times(new Matrix(a_t ).times(new Matrix(a_daz))));
				rot_derivs[chn][1] = (new Matrix(a_r ).times(new Matrix(a_dt).times(new Matrix(a_az ))));
				rot_derivs[chn][2] = (new Matrix(a_dr).times(new Matrix(a_t ).times(new Matrix(a_az ))));
				rot_derivs[chn][3] = (new Matrix(a_dzoom).times(new Matrix(a_t ).times(new Matrix(a_az ))));
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
	 * @param use_rig_offsets - for the auxiliary camera - use offsets from the main one
	 * @param rots misalignment correction (now includes zoom in addition to rotations
	 * @param deriv_rots derivatives by d_az, f_elev, d_rot, d_zoom
	 * @param px pixel X coordinate
	 * @param py pixel Y coordinate
	 * @param disparity disparity
	 * @return array of per port pairs of pixel shifts
	 */

	public double [][] getPortsCoordinatesAndDerivatives(
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
 * 1) unapply common distortion
 * 2) apply disparity
 * 3) apply rotations and zoom
 * 4) re-apply distortion
 * 5) return port center X and Y
 */
		double [][] rXY = getRXY(use_rig_offsets); // may include rig offsets

		double [][] pXY = new double [numSensors][2];

		double pXcd = px - 0.5 * this.pixelCorrectionWidth;
		double pYcd = py - 0.5 * this.pixelCorrectionHeight;
		double rD = Math.sqrt(pXcd*pXcd + pYcd*pYcd)*0.001*this.pixelSize; // distorted radius in a virtual center camera
		double rND2R=getRByRDist(rD/this.distortionRadius, (debugLevel > -1));
		double pXc = pXcd * rND2R; // non-distorted coordinates relative to the (0.5 * this.pixelCorrectionWidth, 0.5 * this.pixelCorrectionHeight)
		double pYc = pYcd * rND2R; // in pixels
		double [] a={this.distortionC,this.distortionB,this.distortionA,this.distortionA5,this.distortionA6,this.distortionA7,this.distortionA8};
		double fl_pix = focalLength/(0.001*pixelSize); // focal length in pixels
		double  ri_scale = 0.001 * this.pixelSize / this.distortionRadius;

		for (int i = 0; i < numSensors; i++){
			// non-distorted XY of the shifted location of the individual sensor
			double pXci0 = pXc - disparity *  rXY[i][0]; // in pixels
			double pYci0 = pYc - disparity *  rXY[i][1];

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
			for (int j = 0; j < a.length; j++){
				rri *= ri;
				rD2rND += a[j]*(rri - 1.0); // Fixed
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

				double dri_droll =     ri_scale / rNDi* (pXci * dpXci_droll +     pYci * dpYci_droll); // Not used anywhere ?
// TODO: verify dri_droll == 0 and remove

				double drD2rND_dri = 0.0;
				rri = 1.0;
				for (int j = 0; j < a.length; j++){
					drD2rND_dri += a[j] * (j+1) * rri;
					rri *= ri;
				}

				double drD2rND_dazimuth = drD2rND_dri * dri_dazimuth;
				double drD2rND_dtilt =    drD2rND_dri * dri_dtilt;

				double dpXid_dazimuth = dpXci_dazimuth * rD2rND + pXci * drD2rND_dazimuth;
				double dpYid_dazimuth = dpYci_dazimuth * rD2rND + pYci * drD2rND_dazimuth;
				double dpXid_dtilt =    dpXci_dtilt    * rD2rND + pXci * drD2rND_dtilt;
				double dpYid_dtilt =    dpYci_dtilt    * rD2rND + pYci * drD2rND_dtilt;


				double dpXid_droll =    dpXci_droll    * rD2rND;
				double dpYid_droll =    dpYci_droll    * rD2rND;

				double dpXid_dzoom =    dpXci_dzoom    * rD2rND;
				double dpYid_dzoom =    dpYci_dzoom    * rD2rND;


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
					use_rig_offsets,
					corr_rots_p, // Matrix []   rots,
					null, // Matrix [][] deriv_rots,
					null, // boolean calc_deriv,
					px, // double px,
					py, // double py,
					disparity // double disparity
					);
			double [][] rslt_m = getPortsCoordinatesAndDerivatives(
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
