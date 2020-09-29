/**
 **
 ** ErsCorrection - inter-frame ERS correction extension for GeometryCorrection
 **
 ** Copyright (C) 2020 Elphel, Inc.
 **
 ** -----------------------------------------------------------------------------**
 **
 **  ErsCorrection.java is free software: you can redistribute it and/or modify
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

package com.elphel.imagej.tileprocessor;

import java.util.ArrayList;
import java.util.HashMap;

import org.apache.commons.math3.complex.Quaternion;
import org.apache.commons.math3.geometry.euclidean.threed.Rotation;
import org.apache.commons.math3.geometry.euclidean.threed.RotationConvention;
import org.apache.commons.math3.geometry.euclidean.threed.RotationOrder;
import org.apache.commons.math3.geometry.euclidean.threed.Vector3D;

import Jama.Matrix;

public class ErsCorrection extends GeometryCorrection {
	public static RotationConvention ROT_CONV = RotationConvention.FRAME_TRANSFORM;
	public static double THRESHOLD = 1E-10;

//	Rotation rotation; // (double[][] m, double THRESHOLD)
	double [][]  ers_xyz;           // per scan line
	double [][]  ers_xyz_dt;        // linear velocitiesper scan line
	Quaternion[] ers_quaternion;    // per scan line
	Quaternion[] ers_quaternion_dt; // per scan line
	double [][]  ers_atr;           // azimuth-tilt-roll per scan line
	double [][]  ers_atr_dt;        // angular velocities per scan line
	/**
	 * Position+orientation (world XYZ, Azimuth, Tilt, Roll) of other scenes relative to the position of this camera.
	 * Positions/orientations are sampled during scanning of the center line 
	 */
	HashMap <String, XyzAtr> scenes_poses = null;
	public class XyzAtr {
		double [] xyz;
		double [] atr;
		public XyzAtr() {
			xyz = new double[3];
			atr = new double[3];
		}
		public XyzAtr(String s) {
			ArrayList<Double> lxyzatr = new ArrayList<Double>();
			for (String snumber : s.split(","))
				lxyzatr.add(Double.parseDouble(snumber));
			Double [] xyzatr = new Double[6];
			xyzatr = lxyzatr.toArray(xyzatr);
			xyz = new double [] {xyzatr[0],xyzatr[1],xyzatr[2]};
			atr = new double [] {xyzatr[3],xyzatr[4],xyzatr[5]};
		}
		
		public String toString() {
			return String.format("%f,  %f, %f, %f, %f, %f",xyz[0],xyz[1],xyz[2],atr[0],atr[1],atr[2]);
		}
		public double [] getXYZ() {
			return xyz;
		}
		public double [] getATR() {
			return atr;
		}
		public void setXYZ(double [] d) {
			xyz[0] = d[0];
			xyz[1] = d[1];
			xyz[2] = d[2];
		}
		public void setATR(double [] d) {
			atr[0] = d[0];
			atr[1] = d[1];
			atr[2] = d[2];
		}
		
	}
	
	
	public ErsCorrection(GeometryCorrection gc) {
		debugLevel =           gc.debugLevel;
		line_time =            gc.line_time;            // 26.5E-6; // duration of sensor scan line (for ERS)
		pixelCorrectionWidth=  gc.pixelCorrectionWidth; // 2592;   // virtual camera center is at (pixelCorrectionWidth/2, pixelCorrectionHeight/2)
		pixelCorrectionHeight= gc.pixelCorrectionHeight; // 1936;
		focalLength =          gc.focalLength; // =FOCAL_LENGTH;
		pixelSize =            gc.pixelSize; // =  PIXEL_SIZE; //um
		distortionRadius =     gc.distortionRadius; // =  DISTORTION_RADIUS; // mm - half width of the sensor
		distortionA8=		   gc.distortionA8; //=0.0; //r^8 (normalized to focal length or to sensor half width?)
		distortionA7 =         gc.distortionA7; //=0.0; //r^7 (normalized to focal length or to sensor half width?)
		distortionA6 =         gc.distortionA6; //=0.0; //r^6 (normalized to focal length or to sensor half width?)
		distortionA5 =         gc.distortionA5; //=0.0; //r^5 (normalized to focal length or to sensor half width?)
		distortionA =          gc.distortionA; //=0.0; // r^4 (normalized to focal length or to sensor half width?)
		distortionB =          gc.distortionB; //=0.0; // r^3
		distortionC =          gc.distortionC; //=0.0; // r^2
		// parameters, common for all sensors
		elevation =            gc.elevation; // 0.0; // degrees, up - positive;
		heading  =             gc.heading; // 0.0;  // degrees, CW (from top) - positive
		numSensors =           gc.numSensors; // 4;
		forward =              gc.forward; // null;
		right =                gc.right; //    null;
		height =               gc.height; //   null;
		roll  =                gc.roll; //    null;  // degrees, CW (to target) - positive
		pXY0 =                 gc.pXY0; //   null;  // sensor center XY in pixels
		common_right =         gc.common_right;    // mm right, camera center
		common_forward =       gc.common_forward;  // mm forward (to target), camera center
		common_height =        gc.common_height;   // mm up, camera center
		common_roll =          gc.common_roll;     // degrees CW (to target) camera as a whole
		XYZ_he =               gc.XYZ_he;     // all cameras coordinates transformed to eliminate heading and elevation (rolls preserved)
		XYZ_her =              gc.XYZ_her; // = null; // XYZ of the lenses in a corrected CCS (adjusted for to elevation, heading,  common_roll)
		rXY =                  gc.rXY; // =     null; // XY pairs of the in a normal plane, relative to disparityRadius
		rXY_ideal =            gc.rXY_ideal; // = {{-0.5, -0.5}, {0.5,-0.5}, {-0.5, 0.5}, {0.5,0.5}};
		cameraRadius =         gc.cameraRadius; // =0; // average distance from the "mass center" of the sensors to the sensors
		disparityRadius =      gc.disparityRadius; //  150.0; // distance between cameras to normalize disparity units to. sqrt(2)*disparityRadius for quad camera (~=150mm)?
		rByRDist =             gc.rByRDist; // =null;
		stepR =                gc.stepR; //  =0.0004; // 0004 - double, 0.0002 - float to fit into GPU shared memory (was 0.001);
		maxR =                 gc.maxR; // =2.0; // calculate up to this*distortionRadius
		m_balance_xy =         gc.m_balance_xy; //  = null; // [2*numSensors][2*numSensors] 8x8 matrix to make XY ports correction to have average == 0
		m_balance_dd =         gc.m_balance_dd; //  = null; // [2*numSensors+1)][2*numSensors] 9x8 matrix to extract disparity from dd
		extrinsic_corr =       gc.extrinsic_corr; // ;
		rigOffset =            gc.rigOffset; //  =    null;
		woi_tops =             gc.woi_tops; //  =     null; // used to calculate scanline timing
	}
	
	public void setupERS(
			double [] wxyz_center,     // world camera XYZ (meters) for the frame center
			double [] wxyz_center_dt,  // world camera Vx, Vy, Vz (m/s)
			double [] wxyz_center_d2t, // world camera Vx, Vy, Vz (m/s^2)
//			double [] watr_center,     // camera orientation (az, tilt, roll in radians, corresponding to the frame center)
			double [] watr_center_dt,  // camera rotaions (az, tilt, roll in radians/s, corresponding to the frame center)
			double [] watr_center_d2t) // camera rotaions (az, tilt, roll in radians/s, corresponding to the frame center)
	{
		ers_xyz=            new double [pixelCorrectionHeight][3];
		ers_xyz_dt=         new double [pixelCorrectionHeight][3];
		ers_quaternion =    new Quaternion [pixelCorrectionHeight];
		ers_quaternion_dt = new Quaternion [pixelCorrectionHeight];
		ers_atr=            new double [pixelCorrectionHeight][3];
		ers_atr_dt=         new double [pixelCorrectionHeight][3];
		int cent_h = pixelCorrectionHeight/2;
		//Rotation r1= new Rotation(RotationOrder.YXZ, RC, ro.aux_azimuth, ro.aux_tilt, ro.aux_roll)
		//Rotation r1= new Rotation(RotationOrder.YXZ, RC, ro.aux_azimuth, ro.aux_tilt, ro.aux_roll)
//		Rotation rcenter0= new Rotation(RotationOrder.YXZ, ROT_CONV, watr_center[0],     watr_center[1],     watr_center[2]);
		Rotation rcenter0= new Rotation(RotationOrder.YXZ, ROT_CONV, 0.0, 0.0, 0.0);
		
		Quaternion quat_center0 = new Quaternion (rcenter0.getQ0(),rcenter0.getQ1(),rcenter0.getQ2(),rcenter0.getQ3());
		Quaternion quat_center1 = new Quaternion (0.0,watr_center_dt[1],   watr_center_dt[0],   watr_center_dt[2]); // angular velocity 1/s :tilt, az, roll
		Quaternion quat_center2 = new Quaternion (0.0,watr_center_d2t[1],  watr_center_d2t[0],  watr_center_d2t[2]); // angular velocity 1/s :tilt, az, roll
		
		// integration to the bottom of the image
		double dt = line_time; 
		double [] wxy0 = wxyz_center.clone();
		double [] wxy1 = wxyz_center_dt.clone();
		double [] wxy2 = wxyz_center_d2t.clone();
		// bottom half rotations
		dt =  line_time; 
		Quaternion q0 = quat_center0.multiply(1.0); // clone() orientation
		Quaternion q1 = quat_center1.multiply(1.0); // clone() angular velocity (pure)
		Quaternion q2 = quat_center2.multiply(1.0); // clone() angular accelerations (pure)
		for (int h = cent_h; h < pixelCorrectionHeight; h ++) {
			ers_quaternion[h] = q0; 
			ers_quaternion_dt[h] = q1;
			Quaternion q1_next = q1.add(q2.multiply(dt));
			Quaternion q_step = q1.add(q1_next).multiply(0.25*dt); 
			q0 = q0.add(q0.multiply(q_step)).normalize(); //
			q1 = q1_next;
		}
		// top half-frame rotations
		dt =  -line_time; 
		q0 = quat_center0.multiply(1.0); // clone() orientation
		q1 = quat_center1.multiply(1.0); // clone() angular velocity (pure)
		q2 = quat_center2.multiply(1.0); // clone() angular accelerations (pure)
		for (int h = cent_h; h >= 0; h--) {
			ers_quaternion[h] = q0; 
			ers_quaternion_dt[h] = q1;
			Quaternion q1_next = q1.add(q2.multiply(dt));
			Quaternion q_step = q1.add(q1_next).multiply(0.25*dt); 
			q0 = q0.add(q0.multiply(q_step)).normalize(); //
			q1 = q1_next;
		}
		
		// convert quaternions for orientations and angular velocities to A,T,R
		for (int h = 0; h < pixelCorrectionHeight; h ++) {
			double [] angles = (new Rotation(
					ers_quaternion[h].getQ0(),
					ers_quaternion[h].getQ1(),
					ers_quaternion[h].getQ2(),
					ers_quaternion[h].getQ3(), true)).getAngles(RotationOrder.YXZ,ROT_CONV); // boolean needsNormalization)
			ers_atr[h][0] = angles[0]; 
			ers_atr[h][1] = angles[1]; 
			ers_atr[h][2] = angles[2];
			ers_atr_dt[h][0] = ers_quaternion_dt[h].getQ2(); // az
			ers_atr_dt[h][1] = ers_quaternion_dt[h].getQ1(); // tl
			ers_atr_dt[h][2] = ers_quaternion_dt[h].getQ3(); // roll
		}
		// TODO: Make linear integration along the camera local axes, coordinates relative to the 
		// Integrate linear velocities/accelerations
		for (int h = cent_h; h < pixelCorrectionHeight; h ++) {
			for (int i = 0; i < 3; i++) {
				ers_xyz[h][i] =    wxy0[i];
				ers_xyz_dt[h][i] = wxy1[i];
				double wxy1_next = wxy1[i] +  wxy2[i] * dt;
				wxy0[i] +=  0.5*(wxy1[i] + wxy1_next) * dt;
				wxy1[i] = wxy1_next;
			}
		}
		dt = -line_time; 
		wxy0 = wxyz_center.clone();
		wxy1 = wxyz_center_dt.clone();
		for (int h = cent_h; h >= 0; h--) {
			for (int i = 0; i < 3; i++) {
				ers_xyz[h][i] =    wxy0[i];
				ers_xyz_dt[h][i] = wxy1[i];
				double wxy1_next = wxy1[i] +  wxy2[i] * dt;
				wxy0[i] +=  0.5*(wxy1[i] + wxy1_next) * dt;
				wxy1[i] = wxy1_next;
			}
		}
		
	}

	/**
	 * Get real world coordinates from pixel coordinates and nominal disparity
	 * @param px horizontal pixel coordinate (right)
	 * @param py vertical pixel coordinate (down)
	 * @param disparity nominal disparity (pixels)
	 * @param correctDistortions true: correct lens distortions, false - no lens distortions
	 * @param camera_xyz camera lens position during centerline acquisition in world coordinates
	 * @param camera_atr camera orientation during centerline acquisition in world frame
	 * @return a vector {x, y, z, 1.0 } in meters. For infinity : {x, y, z, 0} 
	 */
	public double [] getWorldCoordinatesERS( // USED in lwir
			double px,
			double py,
			double disparity,
			boolean correctDistortions,// correct distortion (will need corrected background too !)
			double [] camera_xyz, // camera center in world coordinates
			double [] camera_atr)  // camera orientation relative to world frame
	{
		double pXcd = px - 0.5 * this.pixelCorrectionWidth;
		double pYcd = py - 0.5 * this.pixelCorrectionHeight;
		double rD = Math.sqrt(pXcd*pXcd + pYcd*pYcd)*0.001*this.pixelSize; // distorted radius in a virtual center camera
		double rND2R = correctDistortions?(getRByRDist(rD/this.distortionRadius, false)): 1.0;
		double pXc = pXcd * rND2R; // non-distorted coordinates relative to the (0.5 * this.pixelCorrectionWidth, 0.5 * this.pixelCorrectionHeight)
		double pYc = pYcd * rND2R; // in pixels
		double zh = -SCENE_UNITS_SCALE * this.focalLength * this.disparityRadius / (0.001*this.pixelSize); // "+" - near, "-" far
		double xh =  SCENE_UNITS_SCALE * pXc * this.disparityRadius;
		double yh = -SCENE_UNITS_SCALE * pYc * this.disparityRadius;
		int line = (int) Math.round(py);
		if (line < 0) {
			line = 0;
		} else if (line >= pixelCorrectionHeight) {
			line = pixelCorrectionHeight - 1;
		}
///		double dt = (line - 0.5 * this.pixelCorrectionHeight) * line_time;
		double [] xyz = {xh, yh, zh};
		
		Vector3D world_xyz;
		Rotation   cam_orient_center= new Rotation(RotationOrder.YXZ, ROT_CONV, camera_atr[0],camera_atr[1],camera_atr[2]);
		// current camera offset in the centerline camera frame 
		Vector3D cam_now_local = new Vector3D(ers_xyz[line]); 
		// camera orientation during pixel acquisition :
		Quaternion qpix = ers_quaternion[line];
		Rotation cam_orient_now_local = new Rotation(qpix.getQ0(), qpix.getQ1(), qpix.getQ2(),qpix.getQ3(), true); // boolean needsNormalization)
		boolean is_infinity = Math.abs(disparity) < THRESHOLD;
		if (!is_infinity) {
			for (int i = 0; i < 3; i++) {
				xyz[i] /= disparity;
			}
		}
		Vector3D v3= new Vector3D(xyz);
		// convert to frame parallel to the camera during center line
		Vector3D cam_center_now_local = cam_orient_now_local.applyInverseTo(v3);
		// get real world xyz relative to the camera acquiring a center line
		Vector3D cam_center_local = (is_infinity) ? cam_center_now_local :  cam_center_now_local.add(cam_now_local); // skip translation for infinity
		// undo camera rotation during acquisition of the center line. 
		Vector3D cam_center_world = cam_orient_center.applyInverseTo(cam_center_local);
		// convert to the real world coordinates
		world_xyz =   (is_infinity) ? cam_center_world : cam_center_world.add(new Vector3D(camera_xyz));
		double [] wxyz = world_xyz.toArray();
		double [] wxyz4 = {wxyz[0],wxyz[1],wxyz[2], 1.0};
		if (is_infinity) {
			wxyz4[3] = 0.0;
		}
		return wxyz4;
	}
	
	/**
	 * Get pixel disparity and coordinates from the real world coordinates (in meters)
	 * @param xyzw real world coordinates {x, y, z, w} in meters (right up, towards camera)
	 * @param correctDistortions true: correct lens distortions, false - no lens distortions
	 * @param camera_xyz camera lens position during centerline acquisition in world coordinates
	 * @param camera_atr camera orientation during centerline acquisition in world frame
	 * @param line_err iterate until the line (pY) correction is below this value
	 * @return {disparity, px, py} (right, down)
	 */
	public double [] getImageCoordinatesERS( // USED in lwir
			double [] xyzw,
			boolean correctDistortions, // correct distortion (will need corrected background too !)
			double [] camera_xyz, // camera center in world coordinates
			double [] camera_atr,  // camera orientation relative to world frame
			double    line_err) // threshold error in scan lines (1.0)
	{
		boolean is_infinity = xyzw[3] == 0;
		Vector3D world_xyz = new Vector3D(xyzw[0],xyzw[1],xyzw[2]);
		if (!is_infinity) {
			world_xyz.scalarMultiply(1.0/xyzw[3]);
		}
		// convert to camera-centered, world-parallel coordinates 
		Vector3D cam_center_world = (is_infinity) ? world_xyz : world_xyz.subtract(new Vector3D(camera_xyz));
		// rotate to match camera coordinates when scanning the center line
		Rotation   cam_orient_center= new Rotation(RotationOrder.YXZ, ROT_CONV, camera_atr[0],camera_atr[1],camera_atr[2]);
		Vector3D cam_center_local = cam_orient_center.applyTo(cam_center_world);
		int line = pixelCorrectionHeight / 2;
		double err = pixelCorrectionHeight / 2;
		double [] dxy = null;
		// multiple iterations starting with no ERS distortions
		while (err > line_err) {
			// current camera offset in the centerline camera frame 
			Vector3D cam_now_local = new Vector3D(ers_xyz[line]);  
			Vector3D cam_center_now_local = (is_infinity) ? cam_center_local :  cam_center_local.subtract(cam_now_local); // skip translation for infinity
			Quaternion qpix = ers_quaternion[line];
			Rotation cam_orient_now_local = new Rotation(qpix.getQ0(), qpix.getQ1(), qpix.getQ2(),qpix.getQ3(), true); // boolean 			
			Vector3D v3 = cam_orient_now_local.applyTo(cam_center_now_local);
			double [] xyz = v3.toArray();
			if (Math.abs(xyz[2]) < THRESHOLD) {
				return null; // object too close to the lens
			}
			double pXc =       -(1000.0*focalLength / pixelSize) * xyz[0] / xyz[2];
			double pYc =        (1000.0*focalLength / pixelSize) * xyz[1] / xyz[2];
			double disparity = is_infinity ? 0.0 : (-(1000.0*focalLength / pixelSize)          / xyz[2] * SCENE_UNITS_SCALE * disparityRadius);
			double rND = Math.sqrt(pXc*pXc + pYc*pYc)*0.001*this.pixelSize; // mm

			double rD2RND = correctDistortions?getRDistByR(rND/this.distortionRadius):1.0;
			double px = pXc * rD2RND + 0.5 * this.pixelCorrectionWidth;  // distorted coordinates relative to the (0.5 * this.pixelCorrectionWidth, 0.5 * this.pixelCorrectionHeight)
			double py = pYc * rD2RND + 0.5 * this.pixelCorrectionHeight; // in pixels
			dxy = new double [] {disparity, px, py};
			int line1 = (int) Math.round(py);
			if (line1 < 0) {
				line1 = 0;
			} else if (line1 >= pixelCorrectionHeight) {
				line1 = pixelCorrectionHeight - 1;
			}
			err = Math.abs(line1 - line);
			line = line1;
		}
		return dxy;
	}
	
	
	
	public static void test_rotations() {
		GeometryCorrection.RigOffset ro = (new GeometryCorrection()). new RigOffset();
		ro.aux_azimuth =   Math.PI/2; // 0.1; // 0.1; // radians, azimuth of the auxiliary camera (positive - looks to the right)
		ro.aux_tilt =      0.1; // 2; //.2; // radians, tilt of the auxiliary camera (positive - looks up)
		ro.aux_roll =      0.0; // 3; //.3; // radians, roll of the auxiliary camera (positive - looks clockwise)
		
		double vx = 1.0;
		double vy = 0.0;
		double vz = 0.0;		
/*		
		Matrix mro = ro.getRotMatrix(); // -az -> +tl -> + roll //RotationOrder.ZXY
		mro.print(8, 5);
		double [][] ar = mro.getArray();
		
		Rotation r = new Rotation(ar, THRESHOLD);
		
		(new Matrix(r.getMatrix())).print(8, 5);
*/		
		/*
		double [] angles;
		RotationConvention RC = RotationConvention.FRAME_TRANSFORM;
		System.out.println("from ro.getRotMatrix():");
//		RotationConvention RC = RotationConvention.VECTOR_OPERATOR;
//		angles = r.getAngles(RotationOrder.XYZ, RC); printAngle("XYZ",angles);
//		angles = r.getAngles(RotationOrder.XZY, RC); printAngle("XZY",angles);
		angles = r.getAngles(RotationOrder.YXZ, RC); printAngle("YXZ",angles);
//		angles = r.getAngles(RotationOrder.YZX, RC); printAngle("YZX",angles);
//		angles = r.getAngles(RotationOrder.ZXY, RC); printAngle("ZXY",angles);
//		angles = r.getAngles(RotationOrder.ZYX, RC); printAngle("ZYX",angles);
///		System.out.println( );
		System.out.println("from new Rotation(RotationOrder.YXZ, RC, ro.aux_azimuth, ro.aux_tilt, ro.aux_roll):");
		Rotation r1= new Rotation(RotationOrder.YXZ, RC, ro.aux_azimuth, ro.aux_tilt, ro.aux_roll);
//		Rotation r1= new Rotation(RotationOrder.ZXY, RC, ro.aux_azimuth, ro.aux_tilt, ro.aux_roll);
		
///		(new Matrix(r1.getMatrix())).print(8, 5);
		
//		angles = r1.getAngles(RotationOrder.XYZ, RC); printAngle("XYZ",angles);
//		angles = r1.getAngles(RotationOrder.XZY, RC); printAngle("XZY",angles);
		angles = r1.getAngles(RotationOrder.YXZ, RC); printAngle("YXZ",angles);
//		angles = r1.getAngles(RotationOrder.YZX, RC); printAngle("YZX",angles);
//		angles = r1.getAngles(RotationOrder.ZXY, RC); printAngle("ZXY",angles);
//		angles = r1.getAngles(RotationOrder.ZYX, RC); printAngle("ZYX",angles);
		
		Quaternion q1 = new Quaternion (r1.getQ0(),r1.getQ1(),r1.getQ2(),r1.getQ3());

		System.out.println("Quaternion from rotation above q1 = "+q1.toString());
		System.out.println("Normalized q1 =                     "+q1.normalize().toString());
//		System.out.println();

		Rotation r2 = new Rotation(q1.getQ0(), q1.getQ1(), q1.getQ2(),q1.getQ3(), true); // boolean needsNormalization)
		// rig.getRotMatrix() matches (with signs)	new Rotation(RotationOrder.YXZ, RotationConvention.FRAME_TRANSFORM, ro.aux_azimuth, ro.aux_tilt, ro.aux_roll);
		printAngle("From q1 : YXZ",r2.getAngles(RotationOrder.YXZ, RC));
		
		// angular velocities
		double w_az = 3; // 3; 
		double w_tl = 2; // 2; 
		double w_rl = 1;
		double [] w = {w_tl, w_az, w_rl};
		Quaternion qw = new Quaternion(w);
		System.out.println("Quaternion from rotation velocity = "+qw.toString());
		double t = 0.1; // full time
		int n = 100;
		double dt = t/n;
		Quaternion q_step =    qw.multiply(0.5*dt);
		Quaternion q_add =q1.multiply(1.0);// Quaternion.IDENTITY; //
		
		for (int i = 0; i < n; i ++) {
			q_add=(q_add.add(q_add.multiply(q_step))).normalize(); // local space (rotation relative to the camera)
//			q_add=(q_add.add(q_step.multiply(q_add))).normalize(); // local space (rotation relative to world)
		}
		System.out.println("Integrated quaternion rq_add = "+q_add.toString());
		System.out.println("Normalized rq_add =            "+q_add.normalize().toString());
		
		Rotation rq_add = new Rotation(q_add.getQ0(), q_add.getQ1(), q_add.getQ2(),q_add.getQ3(), true); // boolean needsNormalization)
		printAngle("rq_add -> YXZ",rq_add.getAngles(RotationOrder.YXZ, RC));
*/
		System.out.println("\n");
//		double az = Math.PI/2;
//		double tl = 0.0;
//		double rl = 0.1;
		
		double az = ro.aux_azimuth;
		double tl = ro.aux_tilt; //  =      0.0; // 2; //.2; // radians, tilt of the auxiliary camera (positive - looks up)
		double rl = ro.aux_roll; //  =      0.0; // 3; //.3; // radians, roll of the auxiliary camera (positive - looks clockwise)

		Rotation rr= new Rotation(RotationOrder.YXZ, ROT_CONV, az,tl,rl);
		printAngle("From ATR : YXZ",rr.getAngles(RotationOrder.YXZ, ROT_CONV));
//		double vx = 0.0;
//		double vy = 1.0;
//		double vz = 0.0;		
		double [] v1 = {vx, vy, vz};
		double [] v_out = new double[3];
		double [] v2 = new double[3];
		rr.applyTo(v1,v_out);
		rr.applyInverseTo(v_out,v2);
		
		printAngle(" v_in", v1);
		printAngle("v_out", v_out);
		printAngle("   v2", v2);
//		System.out.println("v_in=["+v1[0]+","+v1[1]+","+v1[2]+"]");
//		System.out.println("v_out["+v_out[0]+","+v_out[1]+","+v_out[2]+"]");
		
		
		
/*		
		System.out.println("\nfrom new rotation 0.3, 0.2, 0.1:");
		Rotation r2a= new Rotation(RotationOrder.YXZ, RC, 0.03, 0.02, 0.01);
		Quaternion q2 = new Quaternion (r2a.getQ0(),r2a.getQ1(),r2a.getQ2(),r2a.getQ3());// q1.multiply(q1);
		System.out.println("Quaternion q2 = "+q2.toString());
		System.out.println("Normalized q2 = "+q2.normalize().toString());
		
		Rotation r3 = new Rotation(q2.getQ0(), q2.getQ1(), q2.getQ2(),q2.getQ3(), true); // boolean needsNormalization)
		printAngle("From q2: YXZ",r3.getAngles(RotationOrder.YXZ, RC));
		
		Quaternion q12 = q1.multiply(q2);
		System.out.println("Quaternion q12=q1*q2 = "+q12.toString());
		System.out.println("Normalized q12=q1*q2 = "+q12.normalize().toString());
		Rotation r12 = new Rotation(q12.getQ0(), q12.getQ1(), q12.getQ2(),q12.getQ3(), true); // boolean needsNormalization)
		printAngle("q12=q1*q2->YXZ",r12.getAngles(RotationOrder.YXZ, RC));
		System.out.println();
		
		

		int n = 1000;
		Quaternion qd = q2.multiply(1.0/n);
		System.out.println("Quaternion qd=q2/"+n+" = "+qd.toString());
		
		Quaternion q3 = Quaternion.IDENTITY; //q1.multiply(1.0);
		for (int i = 0; i < n; i ++) {
			q3=(q3.add(q3.multiply(qd))).normalize(); // no difference
		}
		System.out.println("Integrated quaternion q3 = "+q3.toString());
		System.out.println("Normalized q3 =            "+q3.normalize().toString());
		
		Rotation r4 = new Rotation(q3.getQ0(), q3.getQ1(), q3.getQ2(),q3.getQ3(), true); // boolean needsNormalization)
		printAngle("Q3 -> YXZ",r4.getAngles(RotationOrder.YXZ, RC));
		
		Vector3D axis =  r2a.getAxis(RC);
		double angle = r2a.getAngle();
		double dangle = angle/n;
//		System.out.println("Quaternion q3 = "+q3.toString());
		Rotation dr = new Rotation(axis, dangle, RC);
		Quaternion qdm = new Quaternion (dr.getQ0(),dr.getQ1(),dr.getQ2(),dr.getQ3());
		System.out.println("Quaternion qdm = "+qdm.toString());

		Quaternion q3m = Quaternion.IDENTITY; // q1.multiply(1.0);
		for (int i = 0; i < n; i ++) {
			q3m=q3m.multiply(qdm); //.normalize();
		}
		System.out.println("Multiplied quaternion q3m = "+q3m.toString());
		System.out.println("Multiplied normalized q3m = "+q3m.normalize().toString());
		Rotation r3m = new Rotation(q3m.getQ0(), q3m.getQ1(), q3m.getQ2(),q3m.getQ3(), true); // boolean needsNormalization)
		printAngle("Multiplied YXZ",r3m.getAngles(RotationOrder.YXZ, RC));

		Quaternion qda = qdm.subtract(Quaternion.IDENTITY);
		System.out.println("Mult-diff quaternion qda = "+qda.toString());
		Quaternion q3a = Quaternion.IDENTITY; //q1.multiply(1.0); // Quaternion.IDENTITY; //
		for (int i = 0; i < n; i ++) {
			q3a=q3a.add(q3a.multiply(qda)); //.normalize();
		}
		System.out.println("Integrated mult-diff quaternion q3a = "+q3a.toString());
		System.out.println("Normalized q3a =                      "+q3a.normalize().toString());
		
		Rotation r4a = new Rotation(q3a.getQ0(), q3a.getQ1(), q3a.getQ2(),q3a.getQ3(), true); // boolean needsNormalization)
		printAngle("YXZ",r4a.getAngles(RotationOrder.YXZ, RC));
*/		
		
		
		
	}
	public static void printAngle (String name, double [] v) {
		System.out.println(String.format("%s: [%9f, %9f, %9f]", name, v[0],v[1],v[2]));
	}
	
	
}





