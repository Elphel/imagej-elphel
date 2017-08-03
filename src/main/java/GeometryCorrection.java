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
	static double SCENE_UNITS_SCALE = 0.001;
	static String SCENE_UNITS_NAME = "m";
	static final String [] CORR_NAMES = {"tilt0","tilt1","tilt2","azimuth0","azimuth1","azimuth2","roll0","roll1","roll2","roll3"};

	public int    debugLevel = 0;
	//	public double azimuth; // azimuth of the lens entrance pupil center, degrees, clockwise looking from top
	//	public double radius;  // mm, distance from the rotation axis
	//	public double height;       // mm, up - from the origin point
	//	public double phi;     // degrees, optical axis from azimuth/r vector, clockwise heading
	//	public double theta;   // degrees, optical axis from the eyesis horizon, positive - up elevation
	//	public double psi;     // degrees, rotation (of the sensor) around the optical axis. Positive if camera is rotated clockwise looking to the target roll
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
	//	public double px0=1296.0;          // center of the lens on the sensor, pixels
	//	public double py0=968.0;           // center of the lens on the sensor, pixels


	// parameters, common for all sensors	
	public double    elevation = 0.0; // degrees, up - positive;
	public double    heading  =  0.0;  // degrees, CW (from top) - positive
	//	public double    roll_common  =     0.0;  // degrees, CW (to target) - positive


	public int       numSensors = 4;
	public double [] forward =   null;
	public double [] right =     null;
	public double [] height =    null;
	public double [] roll  =     null;  // degrees, CW (to target) - positive
	public double [][] pXY0 =    null;  // sensor center XY in pixels 

	public double common_right;    // mm right, camera center
	public double common_forward;  // mm forward (to target), camera center
	public double common_height;   // mm up, camera center
	public double common_roll;     // degrees CW (to target) camera as a whole
	public double [][] XYZ_he;     // all cameras coordinates transformed to eliminate heading and elevation (rolls preserved) 
	public double [][] XYZ_her = null; // XYZ of the lenses in a corrected CCS (adjusted for to elevation, heading,  common_roll)
	public double [][] rXY =     null; // XY pairs of the in a normal plane, relative to disparityRadius
	public double [][] rXY_ideal = {{-0.5, -0.5}, {0.5,-0.5}, {-0.5, 0.5}, {0.5,0.5}}; 

	public double cameraRadius=0; // average distance from the "mass center" of the sensors to the sensors
	public double disparityRadius=0; // distance between cameras to normalize disparity units to. sqrt(2)*disparityRadius for quad camera (~=150mm)?

	private double [] rByRDist=null;
	private double    stepR=0.001;
	private double    maxR=2.0; // calculate up to this*distortionRadius

	public  CorrVector extrinsic_corr;

	public GeometryCorrection(double [] extrinsic_corr)
	{
		this.extrinsic_corr = 	new CorrVector(extrinsic_corr);
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

	
	
	public class CorrVector{
		static final int LENGTH =       10;
		static final int TILT_INDEX =    0;
		static final int AZIMUTH_INDEX = 3;
		static final int ROLL_INDEX =    6;
		double [] vector;
		public CorrVector ()
		{
			this.vector = new double[10];
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
				double roll0,    double roll1,    double roll2, double roll3)
		{}
		public CorrVector (double [] vector)
		{
			this.vector = vector;		
		}
		/**
		 * Set subcamera corrections from a single array
		 * @param tilt    for subcameras 0..2, radians, positive - up (subcamera 3 so sum == 0)
		 * @param azimuth for subcameras 0..2, radians, positive - right (subcamera 3 so sum == 0)
		 * @param roll    for subcameras 0..3, radians, positive - CW looking to the target
		 */
		public CorrVector (double [] tilt, double [] azimuth, double [] roll)
		{
			double [] vector = {
					tilt[0],    tilt[1],    tilt[2],
					azimuth[0], azimuth[1], azimuth[2],
					roll[0], roll[1], roll[2], roll[3]};
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
		public String toString()
		{
			String s;
			double [] sym_vect = toSymArray(null);
			double [] v = new double [vector.length];
			double [] sv = new double [vector.length];
			for (int i = 0; i < vector.length; i++){
				v[i] =  vector[i]*180/Math.PI;
				sv[i] = sym_vect[i]*180/Math.PI;
			}
			
			s  = String.format("tilt    (up):    %8.5f° %8.5f° %8.5f° %8.5f°\n" , v[0], v[1], v[2], -(v[0] + v[1] + v[2]) );
			s += String.format("azimuth (right): %8.5f° %8.5f° %8.5f° %8.5f°\n" , v[3], v[4], v[5], -(v[3] + v[4] + v[5]) );
			s += String.format("roll    (CW):    %8.5f° %8.5f° %8.5f° %8.5f°\n" , v[6], v[7], v[8], v[9] );
			s += "Symmetrical vector first 5 as last 4 (rolls) are the same:\n";
			s += "   |↘ ↙|     |↘ ↗|     |↗ ↘|     |↙ ↘|      |↙ ↗|     |↖  ↘|\n";
			s += "0: |↗ ↖|  1: |↙ ↖|  2: |↖ ↙|  3: |↖ ↗|  4:  |↗ ↙|  5: |↘  ↖|\n";
			s += String.format("0: %8.5f° 1:%8.5f° 2:%8.5f° 3:%8.5f° 4:%8.5f° 5:%8.5f°\n" , sv[0], sv[1], sv[2], sv[3], sv[4], sv[5]);
			
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
		
		
		public CorrVector clone(){
			return new CorrVector(this.vector.clone());
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
			//Previous does not consider movement of the 4-th sensor, so t3 = -t0+t1+t2), a3 = -(a0+a1+a2) 					
			double [][] tar_to_sym = {
					{-2.0, -2.0,  2.0, -2.0,  0.0,  0.0, 0.0, 0.0, 0.0, 0.0}, // t0
					{-2.0,  0.0,  0.0, -2.0,  2.0, -2.0, 0.0, 0.0, 0.0, 0.0}, // t1
					{ 0.0, -2.0,  2.0,  0.0,  2.0, -2.0, 0.0, 0.0, 0.0, 0.0}, // t2
					
					{ 2.0,  2.0,  2.0, -2.0,  0.0,  0.0, 0.0, 0.0, 0.0, 0.0}, // a0
					{ 0.0,  2.0,  2.0,  0.0,  2.0,  2.0, 0.0, 0.0, 0.0, 0.0}, // a1
					{ 2.0,  0.0,  0.0, -2.0,  2.0,  2.0, 0.0, 0.0, 0.0, 0.0}, // a2
					
					{ 0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 1.0, 0.0, 0.0, 0.0}, // roll 0
					{ 0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.0, 1.0, 0.0, 0.0}, // roll 1
					{ 0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.0, 0.0, 1.0, 0.0}, // roll 2
					{ 0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.0, 0.0, 0.0, 1.0}};// roll 3

			return tar_to_sym;
		}
		public double [][] dTar_j_dSym_i(){
			//			Matrix sym_to_tar = new Matrix(symToTARArray());
			//			Matrix tar_to_sym = sym_to_tar.inverse();
			//			return tar_to_sym.getArray();
			/*
>>> a=[[-2, -2, 2, -2, 0, 0, 0, 0, 0, 0], [-2, 0, 0, -2, 2, -2, 0, 0, 0, 0], [0, -2, 2, 0, 2, -2, 0, 0, 0, 0], [2, 2, 2, -2, 0, 0, 0, 0, 0, 0], [0, 2, 2, 0, 2, 2, 0, 0, 0, 0],
... [2, 0, 0, -2, 2, 2, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 1, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 1, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 1, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 1]]
>>> m = numpy.matrix(a)
>>> m
matrix([[-2, -2,  2, -2,  0,  0,  0,  0,  0,  0],
        [-2,  0,  0, -2,  2, -2,  0,  0,  0,  0],
        [ 0, -2,  2,  0,  2, -2,  0,  0,  0,  0],
        [ 2,  2,  2, -2,  0,  0,  0,  0,  0,  0],
        [ 0,  2,  2,  0,  2,  2,  0,  0,  0,  0],
        [ 2,  0,  0, -2,  2,  2,  0,  0,  0,  0],
        [ 0,  0,  0,  0,  0,  0,  1,  0,  0,  0],
        [ 0,  0,  0,  0,  0,  0,  0,  1,  0,  0],
        [ 0,  0,  0,  0,  0,  0,  0,  0,  1,  0],
        [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  1]])
>>> numpy.linalg.inv(m)
matrix([[-0.125, -0.125,  0.125,  0.125, -0.125,  0.125,  0.   , -0.   ,  -0.   , -0.   ],
        [-0.125,  0.125, -0.125,  0.125,  0.125, -0.125,  0.   ,  0.   ,   0.   ,  0.   ],
        [ 0.125, -0.125,  0.125,  0.125,  0.125, -0.125,  0.   ,  0.   ,   0.   ,  0.   ],
        [-0.125, -0.125,  0.125, -0.125,  0.125, -0.125,  0.   ,  0.   ,   0.   ,  0.   ],
        [-0.125,  0.125,  0.125, -0.125,  0.125,  0.125,  0.   ,  0.   ,   0.   ,  0.   ],
        [ 0.125, -0.125, -0.125, -0.125,  0.125,  0.125,  0.   ,  0.   ,   0.   ,  0.   ],
        [ 0.   ,  0.   ,  0.   ,  0.   ,  0.   ,  0.   ,  1.   ,  0.   ,   0.   ,  0.   ],
        [ 0.   ,  0.   ,  0.   ,  0.   ,  0.   ,  0.   ,  0.   ,  1.   ,   0.   ,  0.   ],
        [ 0.   ,  0.   ,  0.   ,  0.   ,  0.   ,  0.   ,  0.   ,  0.   ,   1.   ,  0.   ],
        [ 0.   ,  0.   ,  0.   ,  0.   ,  0.   ,  0.   ,  0.   ,  0.   ,   0.   ,  1.   ]])
			 */
			double [][] sym_to_tar=	{
					// t0       t1      t2     a0       a1     a2     r0       r1      r2       r3    
					{-0.125, -0.125,  0.125,  0.125, -0.125,  0.125,  0.0  ,  0.0  ,   0.0  ,  0.0  },  // sym0
			        {-0.125,  0.125, -0.125,  0.125,  0.125, -0.125,  0.0  ,  0.0  ,   0.0  ,  0.0  },  // sym1
			        { 0.125, -0.125,  0.125,  0.125,  0.125, -0.125,  0.0  ,  0.0  ,   0.0  ,  0.0  },  // sym2
			        {-0.125, -0.125,  0.125, -0.125,  0.125, -0.125,  0.0  ,  0.0  ,   0.0  ,  0.0  },  // sym3
			        {-0.125,  0.125,  0.125, -0.125,  0.125,  0.125,  0.0  ,  0.0  ,   0.0  ,  0.0  },  // sym4
			        { 0.125, -0.125, -0.125, -0.125,  0.125,  0.125,  0.0  ,  0.0  ,   0.0  ,  0.0  },  // sym5
			        { 0.0  ,  0.0  ,  0.0  ,  0.0  ,  0.0  ,  0.0  ,  1.0  ,  0.0  ,   0.0  ,  0.0  },  // sym6 = r0
			        { 0.0  ,  0.0  ,  0.0  ,  0.0  ,  0.0  ,  0.0  ,  0.0  ,  1.0  ,   0.0  ,  0.0  },  // sym7 = r1
			        { 0.0  ,  0.0  ,  0.0  ,  0.0  ,  0.0  ,  0.0  ,  0.0  ,  0.0  ,   1.0  ,  0.0  },  // sym8 = r2
			        { 0.0  ,  0.0  ,  0.0  ,  0.0  ,  0.0  ,  0.0  ,  0.0  ,  0.0  ,   0.0  ,  1.0  }}; // sym9 = r3
			return sym_to_tar;
		}
		/**
		 * Get partial transposed Jacobian as 2d array (for one measurement set) from partial Jacobian for each sample
		 * with derivatives of port coordinates (all 4) by 3 tilts (ports 0..2), 3 azimuths (ports 0..2) and all 4 rolls
		 * Tilt and azimuth for port 3 is calculated so center would not move. Tilt is positive up, azimuth - right and
		 * roll - clockwise
		 * 
		 * Result is transposed Jacobian (rows (9 or 10) - parameters, columns - port coordinate components (8). Parameters
		 * here are symmetrical, 0 is disparity-related (all to the center), remaining 9 preserve disparity and 
		 * are only responsible for the "lazy eye" (last 4  are for roll):
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
	 * Calculate pixel coordinates for each of numSensors images, for a given (px,py) of the idealized "center" (still distorted) image
	 * and generic disparity, measured in pixels 
	 * @param px pixel X coordinate
	 * @param py pixel Y coordinate
	 * @param disparity disparity
	 * @return array of per port pairs of pixel shifts
	 */
	
	public double [][] getPortsCoordinates(
			double px,
			double py,
			double disparity)
	{
//		boolean new_ports = true; // false;
//		if (new_ports) {
		return getPortsCoordinatesAndDerivatives(
				getCorrVector(), // CorrVector corr_vector,
				null,           // boolean calc_deriv,
				px, // double px,
				py, // double py,
				disparity); // double disparity)
//		} else {
//			return  getPortsCoordinates_nocorr( // old, tested
//					px,
//					py,
//					disparity);
//		}
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

	public double [][] getPortsCoordinates_nocorr( // old, tested
			double px,
			double py,
			double disparity)
	{
		double [][] pXY = new double [numSensors][2];
		double pXcd = px - 0.5 * this.pixelCorrectionWidth;
		double pYcd = py - 0.5 * this.pixelCorrectionHeight;
		double rD = Math.sqrt(pXcd*pXcd + pYcd*pYcd)*0.001*this.pixelSize; // distorted radius in a virtual center camera
		double rND2R=getRByRDist(rD/this.distortionRadius, (debugLevel > -1));
		double pXc = pXcd * rND2R; // non-distorted coordinates relative to the (0.5 * this.pixelCorrectionWidth, 0.5 * this.pixelCorrectionHeight)
		double pYc = pYcd * rND2R; // in pixels
		double [] a={this.distortionC,this.distortionB,this.distortionA,this.distortionA5,this.distortionA6,this.distortionA7,this.distortionA8};
		for (int i = 0; i < numSensors; i++){
			// non-distorted XY of the shifted location of the individual sensor
			double pXci = pXc - disparity *  this.rXY[i][0]; // in pixels
			double pYci = pYc - disparity *  this.rXY[i][1];
			// calculate back to distorted
			double rNDi = Math.sqrt(pXci*pXci + pYci*pYci); // in pixels
			//		Rdist/R=A8*R^7+A7*R^6+A6*R^5+A5*R^4+A*R^3+B*R^2+C*R+(1-A6-A7-A6-A5-A-B-C)");
			double ri = rNDi* 0.001 * this.pixelSize / this.distortionRadius; // relative to distortion radius
			//    		double rD2rND = (1.0 - distortionA8 - distortionA7 - distortionA6 - distortionA5 - distortionA - distortionB - distortionC);
			double rD2rND = 1.0;
			double rri = 1.0;
			for (int j = 0; j < a.length; j++){
				rri *= ri;
//				rD2rND += a[j]*(rri - a[j]); // BUG here !!!! - fix later, will need to re-adjust all fine corr
				rD2rND += a[j]*(rri - 1.0); // Fixed
			}
			double pXid = pXci * rD2rND;  
			double pYid = pYci * rD2rND;
			// individual rotate (check sign)
//			double a_roll = (this.roll[i] - this.common_roll) * Math.PI/180.0 ;
			double a_roll = (this.roll[i]) * Math.PI/180.0 ;
			double c_roll = Math.cos(a_roll);
			double s_roll = Math.sin(a_roll);
//			double c_roll = Math.cos((this.roll[i] - this.common_roll) * Math.PI/180.0);
//			double s_roll = Math.sin((this.roll[i] - this.common_roll) * Math.PI/180.0);
			pXY[i][0] =  c_roll *  pXid + s_roll* pYid + this.pXY0[i][0];
			pXY[i][1] = -s_roll *  pXid + c_roll* pYid + this.pXY0[i][1];
		}
		return pXY;
	}
	
	
	/**
	 * Calculate pixel coordinates for each of numSensors images, for a given (px,py) of the idealized "center" (still distorted) image
	 * and generic disparity, measured in pixels 
	 * @param corr_vector misalignment correction 
	 * @param px pixel X coordinate
	 * @param py pixel Y coordinate
	 * @param disparity disparity
	 * @return array of per port pairs of pixel shifts
	 */
	
	public double [][] getPortsCoordinatesAndDerivatives(
			CorrVector corr_vector,
			double [][] pXYderiv, // if not null, should be double[8][]
			double px,
			double py,
			double disparity)
	{
		String dbg_s = corr_vector.toString();
		
		double [][] pXY = new double [numSensors][2];
		
		double pXcd = px - 0.5 * this.pixelCorrectionWidth;
		double pYcd = py - 0.5 * this.pixelCorrectionHeight;
		double rD = Math.sqrt(pXcd*pXcd + pYcd*pYcd)*0.001*this.pixelSize; // distorted radius in a virtual center camera
		double rND2R=getRByRDist(rD/this.distortionRadius, (debugLevel > -1));
		double pXc = pXcd * rND2R; // non-distorted coordinates relative to the (0.5 * this.pixelCorrectionWidth, 0.5 * this.pixelCorrectionHeight)
		double pYc = pYcd * rND2R; // in pixels
		double [] a={this.distortionC,this.distortionB,this.distortionA,this.distortionA5,this.distortionA6,this.distortionA7,this.distortionA8};
		double corr_scale = focalLength/(0.001*pixelSize);
		double  ri_scale = 0.001 * this.pixelSize / this.distortionRadius;
		for (int i = 0; i < numSensors; i++){
			// non-distorted XY of the shifted location of the individual sensor
			double pXci = pXc - disparity *  this.rXY[i][0]; // in pixels
			double pYci = pYc - disparity *  this.rXY[i][1];
			// add misalignment, for small angles express in non-distorted pixels 
				pXci +=  -corr_vector.getAzimuth(i) * corr_scale;
				pYci +=   corr_vector.getTilt(i)    * corr_scale;
			// calculate back to distorted
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
			double pXid = pXci * rD2rND;  
			double pYid = pYci * rD2rND;
			// individual rotate (check sign)
//			double a_roll = (this.roll[i] - this.common_roll) * Math.PI/180.0 + corr_vector.getRoll(i) ;
			double a_roll = (this.roll[i]) * Math.PI/180.0 + corr_vector.getRoll(i) ;
			double c_roll = Math.cos(a_roll);
			double s_roll = Math.sin(a_roll);
			pXY[i][0] =  c_roll *  pXid + s_roll * pYid + this.pXY0[i][0];
			pXY[i][1] = -s_roll *  pXid + c_roll * pYid + this.pXY0[i][1];
			//		double [][] pXYderiv = calc_deriv? new double [numSensors*2][CorrVector.LENGTH] : null;

			if (pXYderiv != null) {
				pXYderiv[2 * i] =   new double [CorrVector.LENGTH];
				pXYderiv[2 * i+1] = new double [CorrVector.LENGTH];
				double dpXci_dazimuth = -corr_scale;
				double dpYci_dtilt =     corr_scale;
				double dri_dazimuth = pXci/rNDi * dpXci_dazimuth * ri_scale;
				double dri_dtilt =    pYci/rNDi * dpYci_dtilt * ri_scale;
				double drD2rND_dri = 0.0;
				rri = 1.0;
				for (int j = 0; j < a.length; j++){
					drD2rND_dri += a[j] * (j+1) * rri;
					rri *= ri;
				}
				
				double drD2rND_dazimuth = drD2rND_dri * dri_dazimuth;  
				double drD2rND_dtilt =    drD2rND_dri * dri_dtilt;  
//				double pXid = pXci * rD2rND;  
//				double pYid = pYci * rD2rND;
				double dpXid_dazimuth = dpXci_dazimuth * rD2rND + pXid * drD2rND_dazimuth;  
				double dpYid_dazimuth =                           pYid * drD2rND_dazimuth;  
				double dpXid_dtilt =                              pXid * drD2rND_dtilt;  
				double dpYid_dtilt =   dpYci_dtilt     * rD2rND + pYid * drD2rND_dtilt;  
				
//  			pXY[i][0] =  c_roll *  pXid + s_roll * pYid + this.pXY0[i][0];
//				pXY[i][1] = -s_roll *  pXid + c_roll * pYid + this.pXY0[i][1];

				double dxi_dazimuth =  c_roll *  dpXid_dazimuth + s_roll * dpYid_dazimuth;
				double dyi_dazimuth = -s_roll *  dpXid_dazimuth + c_roll * dpYid_dazimuth;
				
				double dxi_dtilt =     c_roll *  dpXid_dtilt + s_roll * dpYid_dtilt;
				double dyi_dtilt =    -s_roll *  dpXid_dtilt + c_roll * dpYid_dtilt;
				// verify that d/dsym are well, symmetrical
				if (i < (numSensors - 1)){
					pXYderiv[2 * i + 0][CorrVector.TILT_INDEX+i] =      dxi_dtilt; 
					pXYderiv[2 * i + 1][CorrVector.TILT_INDEX+i] =      dyi_dtilt; 
					pXYderiv[2 * i + 0][CorrVector.AZIMUTH_INDEX+i] =   dxi_dazimuth; 
					pXYderiv[2 * i + 1][CorrVector.AZIMUTH_INDEX+i] =   dyi_dazimuth; 
				} else {
					for (int j = 0; j < (numSensors - 1); j++){
						pXYderiv[2 * i + 0][CorrVector.TILT_INDEX+j] = -dxi_dtilt; 
						pXYderiv[2 * i + 1][CorrVector.TILT_INDEX+j] = -dyi_dtilt; 
						pXYderiv[2 * i + 0][CorrVector.AZIMUTH_INDEX+j] = -dxi_dazimuth; 
						pXYderiv[2 * i + 1][CorrVector.AZIMUTH_INDEX+j] = -dyi_dazimuth; 
					}
				}
				double dxi_droll = -s_roll * pXid + c_roll * pYid;   
				double dyi_droll = -c_roll * pXid - s_roll * pYid;
				pXYderiv[2 * i + 0][CorrVector.ROLL_INDEX+i] = dxi_droll; 
				pXYderiv[2 * i + 1][CorrVector.ROLL_INDEX+i] = dyi_droll; 
			}
		}
		return pXY;
	}

	public double [][] getPortsCoordinatesAndDerivatives(
			double [] dbg_a_vector, // replace actual radial distortion coefficients
			CorrVector corr_vector,
			double [][] pXYderiv, // if not null, should be double[8][]
			double px,
			double py,
			double disparity)
	{
//		String dbg_s = corr_vector.toString();
		
		double [][] pXY = new double [numSensors][2];
		
		double pXcd = px - 0.5 * this.pixelCorrectionWidth;
		double pYcd = py - 0.5 * this.pixelCorrectionHeight;
		double rD = Math.sqrt(pXcd*pXcd + pYcd*pYcd)*0.001*this.pixelSize; // distorted radius in a virtual center camera
		double rND2R=getRByRDist(rD/this.distortionRadius, (debugLevel > -1));
		double pXc = pXcd * rND2R; // non-distorted coordinates relative to the (0.5 * this.pixelCorrectionWidth, 0.5 * this.pixelCorrectionHeight)
		double pYc = pYcd * rND2R; // in pixels
		double [] a={this.distortionC,this.distortionB,this.distortionA,this.distortionA5,this.distortionA6,this.distortionA7,this.distortionA8};
		if (dbg_a_vector != null){
			a = dbg_a_vector;
		}
		double corr_scale = focalLength/(0.001*pixelSize);
		double  ri_scale = 0.001 * this.pixelSize / this.distortionRadius;
		for (int i = 0; i < numSensors; i++){
			// non-distorted XY of the shifted location of the individual sensor
			double pXci0 = pXc - disparity *  this.rXY[i][0]; // in pixels
			double pYci0 = pYc - disparity *  this.rXY[i][1];
			// add misalignment, for small angles express in non-distorted pixels 
			double pXci = pXci0  - corr_vector.getAzimuth(i) * corr_scale;
			double pYci = pYci0  + corr_vector.getTilt(i)    * corr_scale;
			// calculate back to distorted
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
			double pXid = pXci * rD2rND;  
			double pYid = pYci * rD2rND;
			// individual rotate (check sign)
//			double a_roll = (this.roll[i] - this.common_roll) * Math.PI/180.0 + corr_vector.getRoll(i) ;
			double a_roll = (this.roll[i]) * Math.PI/180.0 + corr_vector.getRoll(i) ;
			double c_roll = Math.cos(a_roll);
			double s_roll = Math.sin(a_roll);
			pXY[i][0] =  c_roll *  pXid + s_roll * pYid + this.pXY0[i][0];
			pXY[i][1] = -s_roll *  pXid + c_roll * pYid + this.pXY0[i][1];
			//		double [][] pXYderiv = calc_deriv? new double [numSensors*2][CorrVector.LENGTH] : null;

			if (pXYderiv != null) {
				pXYderiv[2 * i] =   new double [CorrVector.LENGTH];
				pXYderiv[2 * i+1] = new double [CorrVector.LENGTH];
				double dpXci_dazimuth = -corr_scale;
				double dpYci_dtilt =     corr_scale;
				double dri_dazimuth = pXci/rNDi * dpXci_dazimuth * ri_scale;
				double dri_dtilt =    pYci/rNDi * dpYci_dtilt *    ri_scale;
				double drD2rND_dri = 0.0;
				rri = 1.0;
				for (int j = 0; j < a.length; j++){
					drD2rND_dri += a[j] * (j+1) * rri;
					rri *= ri;
				}
				
				double drD2rND_dazimuth = drD2rND_dri * dri_dazimuth;  
				double drD2rND_dtilt =    drD2rND_dri * dri_dtilt;
				
				// TODO: understand! hack - testing seems to work, but why? 
				drD2rND_dazimuth /=  rD2rND;  
				drD2rND_dtilt    /=  rD2rND;
				
				
//				double pXid = pXci * rD2rND;  
//				double pYid = pYci * rD2rND;
				
				double dpXid_dazimuth = dpXci_dazimuth * rD2rND + pXid * drD2rND_dazimuth;  
				double dpYid_dazimuth =                           pYid * drD2rND_dazimuth;  
				double dpXid_dtilt =                              pXid * drD2rND_dtilt;  
				double dpYid_dtilt =   dpYci_dtilt     * rD2rND + pYid * drD2rND_dtilt;  
				
//  			pXY[i][0] =  c_roll *  pXid + s_roll * pYid + this.pXY0[i][0];
//				pXY[i][1] = -s_roll *  pXid + c_roll * pYid + this.pXY0[i][1];

				double dxi_dazimuth =  c_roll *  dpXid_dazimuth + s_roll * dpYid_dazimuth;
				double dyi_dazimuth = -s_roll *  dpXid_dazimuth + c_roll * dpYid_dazimuth;
				
				double dxi_dtilt =     c_roll *  dpXid_dtilt + s_roll * dpYid_dtilt;
				double dyi_dtilt =    -s_roll *  dpXid_dtilt + c_roll * dpYid_dtilt;
				// verify that d/dsym are well, symmetrical
				if (i < (numSensors - 1)){
					pXYderiv[2 * i + 0][CorrVector.TILT_INDEX+i] =      dxi_dtilt; 
					pXYderiv[2 * i + 1][CorrVector.TILT_INDEX+i] =      dyi_dtilt; 
					pXYderiv[2 * i + 0][CorrVector.AZIMUTH_INDEX+i] =   dxi_dazimuth; 
					pXYderiv[2 * i + 1][CorrVector.AZIMUTH_INDEX+i] =   dyi_dazimuth; 
				} else {
					for (int j = 0; j < (numSensors - 1); j++){
						pXYderiv[2 * i + 0][CorrVector.TILT_INDEX+j] = -dxi_dtilt; 
						pXYderiv[2 * i + 1][CorrVector.TILT_INDEX+j] = -dyi_dtilt; 
						pXYderiv[2 * i + 0][CorrVector.AZIMUTH_INDEX+j] = -dxi_dazimuth; 
						pXYderiv[2 * i + 1][CorrVector.AZIMUTH_INDEX+j] = -dyi_dazimuth; 
					}
				}
				double dxi_droll = -s_roll * pXid + c_roll * pYid;   
				double dyi_droll = -c_roll * pXid - s_roll * pYid;
				pXYderiv[2 * i + 0][CorrVector.ROLL_INDEX+i] = dxi_droll; 
				pXYderiv[2 * i + 1][CorrVector.ROLL_INDEX+i] = dyi_droll; 
			}
		}
		return pXY;
	}
	
	
	
	// debug version, calculates derivatives as differences
	public double [][] getPortsCoordinatesAndDerivatives(
			double [] dbg_a_vector, // replace actual radial distortion coefficients
			double delta, // 1e-6
			CorrVector corr_vector,
			double [][] pXYderiv, // if not null, should be double[8][]
			double px,
			double py,
			double disparity)
	{
		
		double [][] pXYderiv0 = new double[8][];
		double [][] rslt = getPortsCoordinatesAndDerivatives(
				dbg_a_vector, // double [] dbg_a_vector, // replace actual radial distortion coefficients
				corr_vector, // CorrVector corr_vector,
				pXYderiv0, // null, // false, // boolean calc_deriv,
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
			double [][] rslt_p = getPortsCoordinatesAndDerivatives(
					dbg_a_vector, // double [] dbg_a_vector, // replace actual radial distortion coefficients
					cv_delta_p, // CorrVector corr_vector,
					null, // boolean calc_deriv,
					px, // double px,
					py, // double py,
					disparity // double disparity
					);
			double [][] rslt_m = getPortsCoordinatesAndDerivatives(
					dbg_a_vector, // double [] dbg_a_vector, // replace actual radial distortion coefficients
					cv_delta_m, // CorrVector corr_vector,
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
	public double [][] getPortsCoordinatesIdeal(
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
		if (macro_scale == 1){
			return getPortsCoordinates(
					px * macro_scale,
					py * macro_scale,
					disparity);
		}
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
//		double d=1.0-this.distortionA8-this.distortionA7-this.distortionA6-this.distortionA5-this.distortionA-this.distortionB-this.distortionC;
		double drDistDr;
		if (use8){
			drDistDr=(((((((7*this.distortionA8)*r + 6*this.distortionA7)*r + 5*this.distortionA6)*r + 4*this.distortionA5)*r + 3*this.distortionA)*r+2*this.distortionB)*r+1*this.distortionC); // +d;
		} else {
//			drDistDr=(((4*this.distortionA5*r + 3*this.distortionA)*r+2*this.distortionB)*r+1*this.distortionC)*r;
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
	
	/*
	public double getDerivRByRDist(double rDist, boolean debug, double delta){
		return (getRByRDist(rDist+delta, false) - getRByRDist(rDist, false))/delta;
	}

	public double getDerivRByRDist(double rDist, boolean debug){
		// add exceptions;
		if (this.rByRDist==null) {
			calcReverseDistortionTable();
			if (debug)System.out.println("getDerivRByRDist("+IJ.d2s(rDist,3)+"): this.rByRDist==null");
			//		return Double.NaN;
		}
		if (rDist<0) {
			if (debug)System.out.println("getDerivRByRDist("+IJ.d2s(rDist,3)+"): rDist<0");
			return Double.NaN;
		}
		int index=(int) Math.floor(rDist/this.stepR);
		if (index>=(this.rByRDist.length-1)) {
			if (debug) System.out.println("getDerivRByRDist("+IJ.d2s(rDist,3)+"): index="+index+">="+(this.rByRDist.length-1));
			return Double.NaN;
		}
		double result=this.rByRDist[index+1]-this.rByRDist[index];
		if ((index > 0) || (index < (this.rByRDist.length-2))){ // interpolate
			double a = rDist/this.stepR-index; // fractional part, 0..1.0
			if ( a <= 0.5){
				result= 2  * (0.5 - a) * (this.rByRDist[index+1] - this.rByRDist[index]) +
					                a *  (this.rByRDist[index+1] - this.rByRDist[index-1]);
			} else {
				result= 2  * (1.0 - a) * (this.rByRDist[index+1] - this.rByRDist[index]) +
						     (a - 0.5) * (this.rByRDist[index+2] - this.rByRDist[index]);
				
			}
			
		}
		if (Double.isNaN(result)){
			if (debug) System.out.println("this.rByRDist["+index+"]="+this.rByRDist[index]);
			if (debug) System.out.println("this.rByRDist["+(index+1)+"]="+this.rByRDist[index+1]);
			if (debug) System.out.println("rDist="+rDist);
			if (debug) System.out.println("(rDist/this.stepR="+(rDist/this.stepR));

		}
		return result/this.stepR;
	}
	*/


}
