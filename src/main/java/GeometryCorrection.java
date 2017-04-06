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
	public double cameraRadius=0; // average distance from the "mass center" of the sensors to the sensors
	public double disparityRadius=0; // distance between cameras to normalize disparity units to. sqrt(2)*disparityRadius for quad camera (~=150mm)?

	private double [] rByRDist=null;
	private double    stepR=0.001;
	private double    maxR=2.0; // calculate up to this*distortionRadius

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
	
	/*
	 * Calculate pixel coordinates for each of numSensors images, for a given (px,py) of the idealized "center" (still distorted) image
	 * and generic disparity, measured in pixels 
	 */

	public double [][] getPortsCoordinates(
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
				rD2rND += a[j]*(rri - a[j]);
			}
			double pXid = pXci * rD2rND;  
			double pYid = pYci * rD2rND;
			// individual rotate (check sign)
			double c_roll = Math.cos((this.roll[i] - this.common_roll) * Math.PI/180.0);
			double s_roll = Math.sin((this.roll[i] - this.common_roll) * Math.PI/180.0);
			pXY[i][0] =  c_roll *  pXid + s_roll* pYid + this.pXY0[i][0];
			pXY[i][1] = -s_roll *  pXid + c_roll* pYid + this.pXY0[i][1];
		}
		return pXY;
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
