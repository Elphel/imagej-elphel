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
import java.util.Arrays;
import java.util.Enumeration;
import java.util.HashMap;
import java.util.Properties;

import org.apache.commons.math3.complex.Quaternion;
import org.apache.commons.math3.geometry.euclidean.threed.Rotation;
import org.apache.commons.math3.geometry.euclidean.threed.RotationConvention;
import org.apache.commons.math3.geometry.euclidean.threed.RotationOrder;
import org.apache.commons.math3.geometry.euclidean.threed.Vector3D;

import com.elphel.imagej.common.ShowDoubleFloatArrays;

import Jama.Matrix;

public class ErsCorrection extends GeometryCorrection {
	static final double INTER_ROT_AZ_SGN = -1.0; // sign of first sin for azimuth rotation
	static final double INTER_ROT_TL_SGN =  1.0; // sign of first sin for tilt rotation
	static final double INTER_ROT_RL_SGN =  1.0; // sign of first sin for roll rotation

	
	static final String XYZ_PREFIX =    "xyz";
	static final String ATR_PREFIX =    "atr";
	static final String ERS_PREFIX =    "ers";
	static final String SCENES_PREFIX = "scenes";
	
	static final String ERS_XYZ_PREFIX =     ERS_PREFIX + "_xyz";
	static final String ERS_XYZ_DT_PREFIX =  ERS_PREFIX + "_xyz_dt";
	static final String ERS_XYZ_D2T_PREFIX = ERS_PREFIX + "_xyz_d2t";
	static final String ERS_ATR_DT_PREFIX =  ERS_PREFIX + "_atr_dt";
	static final String ERS_ATR_D2T_PREFIX = ERS_PREFIX + "_atr_d2t";
	
	static final String [] DW_DERIV_NAMES = {
			"dw_dpX",  //  (pix)     0
			"dw_dpY",  // (pix)     1
			"dw_dd",   // (pix)     2
			"dw_dvaz", // (rad/sec) 3
			"dw_dvtl", // (rad/sec) 4
			"dw_dvrl", // (rad/sec) 5
			"dw_dvx",  // (m/s)     6
			"dw_dvy",  // (m/s)     7
			"dw_dvz",  // (m/s)     8
			"dw_daz",  // (rad)     9
			"dw_dtl",  // (rad)    10
			"dw_drl",  // (rad)    11
			"dw_dx",   // (m)      12
			"dw_dy",   // (m)	   13
			"dw_dz"};  // (m)      14

	// returned arrays have the zero element with coordinates, not derivatives
	static final int DW_DPX =  0; // dw_dpX,  (pix)
	static final int DW_DPY =  1; // dw_dpY   (pix)
	static final int DW_DD =   2; // dw_dd,   (pix)
	static final int DW_DVAZ = 3; // dw_dvaz, (rad/sec)
	static final int DW_DVTL = 4; // dw_dvtl, (rad/sec)
	static final int DW_DVRL = 5; // dw_dvrl, (rad/sec)
	static final int DW_DVX =  6; // dw_dvx,  (m/s)
	static final int DW_DVY =  7; // dw_dvy,  (m/s)
	static final int DW_DVZ =  8; // dw_dvz,  (m/s)
	static final int DW_DAZ =  9; // dw_daz,  (rad)
	static final int DW_DTL = 10; // dw_dtl,  (rad)
	static final int DW_DRL = 11; // dw_drl,  (rad)
	static final int DW_DX =  12; // dw_dx,   (m)
	static final int DW_DY =  13; // dw_dy,   (m)
	static final int DW_DZ =  14; // dw_dz};  (m)

	static final String [] DP_DERIV_NAMES = {
			"pX",            // (pix)     0
			"pY",            // (pix)     1
			"disp",          // (pix)     2
			"ers_vaz_ref",   // (rad/sec) 3
			"ers_vtl_ref",   // (rad/sec) 4
			"ers_vrl_ref",   // (rad/sec) 5
			"ers_dvx_ref",   // (m/s)     6
			"ers_dvy_ref",   // (m/s)     7
			"ers_dvz_ref",   // (m/s)     8
			"azimuth_ref",   // (rad)     9
			"tilt_ref",      // (rad)    10
			"roll_ref",      // (rad)    11
			"X_ref",         // (m)      12
			"Y_ref",         // (m)	    13
			"Z_ref",         // (m)      14
			"ers_vaz_scene", // (rad/sec)15
			"ers_vtl_scene", // (rad/sec)16
			"ers_vrl_scene", // (rad/sec)17
			"ers_dvx_scene", // (m/s)    18
			"ers_dvy_scene", // (m/s)    19
			"ers_dvz_scene", // (m/s)    20
			"azimuth_scene", // (rad)    21
			"tilt_scene",    // (rad)    22
			"Roll_scene",    // (rad)    23
			"X_scene",       // (m)      24
			"Y_scene",       // (m)	     25
			"Z_scene"};      // (m)      26

	static final String [] DP_VECTORS_NAMES = {
			"pXpYD",          // (pix)     0
			null,             // (pix)     1
			null,             // (pix)     2
			"ERS_VATR_REF",   // (rad/sec) 3
			null,             // (rad/sec) 4
			null,             // (rad/sec) 5
			"ERS_VXYZ_REF",   // (m/s)     6
			null,             // (m/s)     7
			null,             // (m/s)     8
			"ATR_REF",        // (rad)     9
			null,             // (rad)    10
			null,             // (rad)    11
			"XYZ_REF",        // (m)      12
			null,             // (m)	    13
			null,             // (m)      14
			"ERS_VATR_SCENE", // (rad/sec)15
			null,             // (rad/sec)16
			null,             // (rad/sec)17
			"ERS_VXYZ_SCENE", // (m/s)    18
			null,             // (m/s)    19
			null,             // (m/s)    20
			"ATR_SCENE",      // (rad)    21
			null,             // (rad)    22
			null,             // (rad)    23
			"XYZ_SCENE",      // (m)      24
			null,             // (m)	     25
			null};            // (m)      26
	
	// returned arrays have the zero element with coordinates, not derivatives
	// Reference parameters
	static final int DP_DPX =   0; // dw_dpX,  (pix)
	static final int DP_DPY =   1; // dw_dpY   (pix)
	static final int DP_DD =    2; // dw_dd,   (pix)
	static final int DP_DVAZ =  3; // dw_dvaz, (rad/sec)
	static final int DP_DVTL =  4; // dw_dvtl, (rad/sec)
	static final int DP_DVRL =  5; // dw_dvrl, (rad/sec)
	static final int DP_DVX =   6; // dw_dvx,  (m/s)
	static final int DP_DVY =   7; // dw_dvy,  (m/s)
	static final int DP_DVZ =   8; // dw_dvz,  (m/s)
	static final int DP_DAZ =   9; // dw_daz,  (rad)
	static final int DP_DTL =  10; // dw_dtl,  (rad)
	static final int DP_DRL =  11; // dw_drl,  (rad)
	static final int DP_DX =   12; // dw_dx,   (m)
	static final int DP_DY =   13; // dw_dy,   (m)
	static final int DP_DZ =   14; // dw_dz};  (m)
	// Scene parameters
	static final int DP_DSVAZ =15; // dw_dvaz, (rad/sec)
	static final int DP_DSVTL =16; // dw_dvtl, (rad/sec)
	static final int DP_DSVRL =17; // dw_dvrl, (rad/sec)
	static final int DP_DSVX = 18; // dw_dvx,  (m/s)
	static final int DP_DSVY = 19; // dw_dvy,  (m/s)
	static final int DP_DSVZ = 20; // dw_dvz,  (m/s)
	static final int DP_DSAZ = 21; // dw_daz,  (rad)
	static final int DP_DSTL = 22; // dw_dtl,  (rad)
	static final int DP_DSRL = 23; // dw_drl,  (rad)
	static final int DP_DSX =  24; // dw_dx,   (m)
	static final int DP_DSY =  25; // dw_dy,   (m)
	static final int DP_DSZ =  26; // dw_dz};  (m)
	static final int DP_NUM_PARS = DP_DSZ+1; 
	
	static final RotationConvention ROT_CONV = RotationConvention.FRAME_TRANSFORM;
	static final double THRESHOLD = 1E-10;
	static final double LINE_ERR = 0.001; // line accuracy for ERS when converting from world to pixels.
	// parameters for the ERS distortion calculation
	public double [] ers_wxyz_center;     // world camera XYZ (meters) for the lens center (in camera coordinates, typically 0)
	public double [] ers_wxyz_center_dt;  // world camera Vx, Vy, Vz (m/s)
	public double [] ers_wxyz_center_d2t; // world camera Vx, Vy, Vz (m/s^2)
	public double [] ers_watr_center_dt;  // camera rotaions (az, tilt, roll in radians/s, corresponding to the frame center)
	public double [] ers_watr_center_d2t; // camera rotaions (az, tilt, roll in radians/s, corresponding to the frame center)
	
// absolute position and orientation of this camera in World coordinates	
	public double []    camera_xyz = new double[3];  // camera center in world coordinates
	public double []    camera_atr = new double[3];  // camera orientation relative to world frame

	public HashMap <String, XyzAtr> scenes_poses = new HashMap <String, XyzAtr>(); // scene timestamp as a key
	
	// Save above data through properties, separately scenes_poses(relative to this), this camera pose camera pose 
	
	
//	Rotation rotation; // (double[][] m, double THRESHOLD)
	private double [][]  ers_xyz;           // per scan line
	private double [][]  ers_xyz_dt;        // linear velocities per scan line
	private Quaternion[] ers_quaternion;    // per scan line
	private Quaternion[] ers_quaternion_dt; // per scan line
	private double [][]  ers_atr;           // azimuth-tilt-roll per scan line
	private double [][]  ers_atr_dt;        // angular velocities per scan line. It is now actually 2*omega!
	

	/*
	static final double ERS_MIN_DISPARITY = 1.0E-6; // to avoid disparity == 0.0 in division
	private double limitedReciprocal(double d) {
		if (d >  ERS_MIN_DISPARITY) return 1/0/d;
		if (d < -ERS_MIN_DISPARITY) return 1/0/d;
		if (d > 0) return 1.0/ERS_MIN_DISPARITY;
		return -1.0/ERS_MIN_DISPARITY;
	}
	*/
	
	public void setPose(
			double []    camera_xyz,
			double []    camera_atr) {
		this.camera_xyz = camera_xyz;
		this.camera_atr = camera_atr;
	}
	public double [] getCameraXYZ() {
		return camera_xyz;
	}
	public double [] getCameraATR() {
		return camera_atr;
	}
	public double [] getErsXYZ_dt() {
		return ers_wxyz_center_dt;
	}
	public double [] getErsATR_dt() {
		return ers_watr_center_dt;
	}
	public double [] getErsXYZ_d2t() {
		return ers_wxyz_center_d2t;
	}
	public double [] getErsATR_d2t() {
		return ers_watr_center_d2t;
	}

	public void setErsDt(
			double []    ers_xyz_dt,
			double []    ers_atr_dt) {
		this.ers_wxyz_center_dt = ers_xyz_dt;
		this.ers_watr_center_dt = ers_atr_dt;
	}
	public void setErsDt_test(
			double []    ers_xyz_dt,
			double []    ers_atr_dt) {
		double k = 1.0; // 0.5;
		this.ers_wxyz_center_dt = new double[] {-ers_xyz_dt[0],-ers_xyz_dt[1],-ers_xyz_dt[2]};
		this.ers_watr_center_dt = new double[] {k* ers_atr_dt[0],-k*ers_atr_dt[1],k*ers_atr_dt[2]};//ers_atr_dt;
	}
	
	
	public void setErsD2t(
			double []    ers_xyz_d2t,
			double []    ers_atr_d2t) {
		this.ers_wxyz_center_d2t = ers_xyz_d2t;
		this.ers_watr_center_d2t = ers_atr_d2t;
	}

	
	public void setPropertiesLineTime(String prefix, Properties properties){
		properties.setProperty(prefix+"_line_time", line_time+"");
	}

	public boolean getPropertiesLineTime(String prefix, Properties properties){
		if (properties.getProperty(prefix+"_line_time")!=null)  {
			line_time =     Double.parseDouble(properties.getProperty(prefix+"_line_time"));
			return true;
		} else {
			line_time = OLD_LINE_TIME;
			return false;
		}
	}
	
	
	public void setPropertiesPose(String prefix, Properties properties){
		properties.setProperty(prefix+XYZ_PREFIX, String.format("%f, %f, %f", camera_xyz[0], camera_xyz[1], camera_xyz[2]));
		properties.setProperty(prefix+ATR_PREFIX, String.format("%f, %f, %f", camera_atr[0], camera_atr[1], camera_atr[2]));
	}

	public boolean getPropertiesPose(String prefix,Properties properties){
		boolean got_data = false;
		if (properties.getProperty(prefix+XYZ_PREFIX)!=null) {camera_xyz = parseDoublesCSV(properties.getProperty(prefix+XYZ_PREFIX)); got_data=true;}
		if (properties.getProperty(prefix+ATR_PREFIX)!=null) {camera_atr = parseDoublesCSV(properties.getProperty(prefix+ATR_PREFIX)); got_data=true;}
		return got_data;
	}

	public void setPropertiesERS(String prefix, Properties properties){
		properties.setProperty(prefix+ERS_XYZ_PREFIX,     String.format("%f, %f, %f", ers_wxyz_center[0],     ers_wxyz_center[1],     ers_wxyz_center[2]));
		properties.setProperty(prefix+ERS_XYZ_DT_PREFIX,  String.format("%f, %f, %f", ers_wxyz_center_dt[0],  ers_wxyz_center_dt[1],  ers_wxyz_center_dt[2]));
		properties.setProperty(prefix+ERS_XYZ_D2T_PREFIX, String.format("%f, %f, %f", ers_wxyz_center_d2t[0], ers_wxyz_center_d2t[1], ers_wxyz_center_d2t[2]));
		properties.setProperty(prefix+ERS_ATR_DT_PREFIX,  String.format("%f, %f, %f", ers_watr_center_dt[0],  ers_watr_center_dt[1],  ers_watr_center_dt[2]));
		properties.setProperty(prefix+ERS_ATR_D2T_PREFIX, String.format("%f, %f, %f", ers_watr_center_d2t[0], ers_watr_center_d2t[1], ers_watr_center_d2t[2]));
	}

	public boolean getPropertiesERS(String prefix,Properties properties){
		boolean got_data = false;
		if (properties.getProperty(prefix+ERS_XYZ_PREFIX)!=null)     {ers_wxyz_center =     parseDoublesCSV(properties.getProperty(prefix+ERS_XYZ_PREFIX)); got_data=true;}
		if (properties.getProperty(prefix+ERS_XYZ_DT_PREFIX)!=null)  {ers_wxyz_center_dt =  parseDoublesCSV(properties.getProperty(prefix+ERS_XYZ_DT_PREFIX)); got_data=true;}
		if (properties.getProperty(prefix+ERS_XYZ_D2T_PREFIX)!=null) {ers_wxyz_center_d2t = parseDoublesCSV(properties.getProperty(prefix+ERS_XYZ_D2T_PREFIX)); got_data=true;}
		if (properties.getProperty(prefix+ERS_ATR_DT_PREFIX)!=null)  {ers_watr_center_dt =  parseDoublesCSV(properties.getProperty(prefix+ERS_ATR_DT_PREFIX)); got_data=true;}
		if (properties.getProperty(prefix+ERS_ATR_D2T_PREFIX)!=null) {ers_watr_center_d2t = parseDoublesCSV(properties.getProperty(prefix+ERS_ATR_D2T_PREFIX)); got_data=true;}
		if (got_data) {
			setupERS(); // calculate arrays
		}
		return got_data;
	}
	public void setPropertiesScenes(String prefix, Properties properties){
		String [] timestamps = getScenes();
		for (String k : timestamps) {
			if (getScene(k) != null) {
				String [] s_scenes = getScene(k).toStrings(); // null pointer
				properties.setProperty(prefix+SCENES_PREFIX+"_"+k,        s_scenes[0]);
				properties.setProperty(prefix+SCENES_PREFIX+"_"+k+"_dt",  s_scenes[1]);
				properties.setProperty(prefix+SCENES_PREFIX+"_"+k+"_d2t", s_scenes[2]);
				//			properties.setProperty(prefix+SCENES_PREFIX+"_"+k, getScene(k).toString());
			}
		}
	}
	
	// do not forget to reset scenes when switching to a new "this" scene.
	public boolean getPropertiesScenes(String parent_prefix,Properties properties){
		boolean got_data = false;
		ArrayList<String> timestamps = new ArrayList<String>();
		String prefix = parent_prefix+SCENES_PREFIX+"_"; 
		for (Enumeration<?> e = properties.propertyNames(); e.hasMoreElements();) {
			String key = (String) e.nextElement();
			if (key.startsWith(prefix) && !key.endsWith("t")) { // _dt, _d2t
				timestamps.add(key.substring(prefix.length()));
			}
		}
		if (!timestamps.isEmpty()) {
			got_data = true;
			for (String ts:timestamps) {
				String pose =    properties.getProperty(prefix+ts);
				String ers_dt =  properties.getProperty(prefix+ts+"_dt");
				String ers_d2t = properties.getProperty(prefix+ts+"_d2t");
//				addScene(ts, new XyzAtr(properties.getProperty(prefix+ts)));
				addScene(ts, new XyzAtr(new String [] {pose, ers_dt, ers_d2t}));
			}
		}
		return got_data;
	}

	

	//propertyNames()
	
	public double [] parseDoublesCSV(String s) {
		String[] snumbers = s.split(",");
		double [] data = new double [snumbers.length];
		for (int i = 0; i < data.length; i++) {
			data[i] = Double.parseDouble(snumbers[i]);
		}
		return data;
	}
	
	public void resetScenes() {
		scenes_poses = new HashMap <String, XyzAtr>();
	}

	public void removeScene(String timestamp) {
		scenes_poses.remove(timestamp);
	}

	public String [] getScenes() {
		int num_scenes = scenes_poses.size();
		String [] scenes = new String[num_scenes];
		int i = 0;
		for (String ts:scenes_poses.keySet()) {
			scenes[i++] = ts;
		}
		Arrays.sort(scenes);
		return scenes;
	}

	public void addScene(String timestamp, XyzAtr scene) {
		scenes_poses.put(timestamp, scene);
	}
	public void addScene(String timestamp, double [] xyz, double [] atr) {
		scenes_poses.put(timestamp, new XyzAtr(xyz, atr));
	}
	
	public void addScene(String timestamp, double [] xyz, double [] atr, double [] ers_xyz_dt, double [] ers_atr_dt) {
		scenes_poses.put(timestamp, new XyzAtr(xyz, atr, ers_xyz_dt, ers_atr_dt));
	}
	//not used
	public void addScene(String timestamp, double [] xyz, double [] atr, double [] ers_xyz_dt, double [] ers_atr_dt, double [] ers_xyz_d2t, double [] ers_atr_d2t) {
		scenes_poses.put(timestamp, new XyzAtr(xyz, atr, ers_xyz_dt, ers_atr_dt, ers_xyz_d2t, ers_atr_d2t));
	}
	
	public XyzAtr getScene(String timestamp) { // null if not found
		return scenes_poses.get(timestamp);
	}

	public double[] getSceneXYZ(String timestamp) {
		XyzAtr scene = scenes_poses.get(timestamp);
		if (scene == null) return null;
		return scene.getXYZ();
	}
	public double[] getSceneATR(String timestamp) {
		XyzAtr scene = scenes_poses.get(timestamp);
		if (scene == null) return null;
		return scene.getATR();
	}
	public double[] getSceneErsXYZ_dt(String timestamp) {
		XyzAtr scene = scenes_poses.get(timestamp);
		if (scene == null) return null;
		return scene.getErsXYZ_dt();
	}
	public double[] getSceneErsATR_dt(String timestamp) {
		XyzAtr scene = scenes_poses.get(timestamp);
		if (scene == null) return null;
		return scene.getErsATR_dt();
	}
	public double[] getSceneErsXYZ_d2t(String timestamp) {
		XyzAtr scene = scenes_poses.get(timestamp);
		if (scene == null) return null;
		return scene.getErsXYZ_d2t();
	}
	public double[] getSceneErsATR_d2t(String timestamp) {
		XyzAtr scene = scenes_poses.get(timestamp);
		if (scene == null) return null;
		return scene.getErsATR_d2t();
	}
	
	
	
	/**
	 * Position+orientation (world XYZ, Azimuth, Tilt, Roll) of other scenes relative to the position of this camera.
	 * Positions/orientations are sampled during scanning of the center line 
	 */
	public class XyzAtr {
		double [] xyz;
		double [] atr;
		
		double [] ers_xyz_dt;  // world camera Vx, Vy, Vz (m/s)
		double [] ers_atr_dt;  // camera rotaions (az, tilt, roll in radians/s, corresponding to the frame center)
		double [] ers_atr_d2t; // camera rotaions (az, tilt, roll in radians/s, corresponding to the frame center)
		double [] ers_xyz_d2t; // world camera Vx, Vy, Vz (m/s^2)
		
		public XyzAtr() {
			xyz = new double[3];
			atr = new double[3];
		}
		
		public XyzAtr(double [] xyz, double [] atr) {
			this.xyz = xyz;
			this.atr = atr;
			this.ers_xyz_dt =  new double[3];
			this.ers_atr_dt =  new double[3];
			this.ers_xyz_d2t = new double[3];
			this.ers_atr_d2t = new double[3];
		}

		public XyzAtr(
				double [] xyz,
				double [] atr,
				double [] ers_xyz_dt,
				double [] ers_atr_dt) {
			this.xyz = xyz;
			this.atr = atr;
			this.ers_xyz_dt = ers_xyz_dt;
			this.ers_atr_dt = ers_atr_dt;
			this.ers_xyz_d2t = new double[3];
			this.ers_atr_d2t = new double[3];
		}

		public XyzAtr(
				double [] xyz,
				double [] atr,
				double [] ers_xyz_dt,
				double [] ers_atr_dt,
				double [] ers_xyz_d2t,
				double [] ers_atr_d2t
				) {
			this.xyz = xyz;
			this.atr = atr;
			this.ers_xyz_dt = ers_xyz_dt;
			this.ers_atr_dt = ers_xyz_dt;
			this.ers_xyz_d2t = ers_xyz_d2t;
			this.ers_atr_d2t = ers_xyz_d2t;
		}
		
		
		
		
		public XyzAtr(String s) {
			double [] d = parseDoublesCSV(s);
			xyz = new double [] {d[0], d[1], d[2]};
			atr = new double [] {d[3], d[4], d[5]};
		}

		public XyzAtr(String [] ss) {
			if (ss[0] != null) {
				double [] d = parseDoublesCSV(ss[0]);
				xyz = new double [] {d[0], d[1], d[2]};
				atr = new double [] {d[3], d[4], d[5]};
				this.ers_xyz_dt =  new double[3];
				this.ers_atr_dt =  new double[3];
				this.ers_xyz_d2t = new double[3];
				this.ers_atr_d2t = new double[3];
				if (ss.length > 1) {
					if (ss[1] != null) {
						d = parseDoublesCSV(ss[1]);
						ers_xyz_dt = new double [] {d[0], d[1], d[2]};
						ers_atr_dt = new double [] {d[3], d[4], d[5]};
					}
					if (ss.length > 2) {
						if (ss[2] != null) {
							d = parseDoublesCSV(ss[2]);
							ers_xyz_d2t = new double [] {d[0], d[1], d[2]};
							ers_atr_d2t = new double [] {d[3], d[4], d[5]};
						}
					}
				} 
			}
		}
		
		public String toString() {
			return String.format("%f,  %f, %f, %f, %f, %f",xyz[0],xyz[1],xyz[2],atr[0],atr[1],atr[2]);
		}
		
		public String [] toStrings() {
			return new String[] {
					String.format("%f,  %f, %f, %f, %f, %f",xyz[0],         xyz[1],         xyz[2],         atr[0],         atr[1],         atr[2]),
					String.format("%f,  %f, %f, %f, %f, %f",ers_xyz_dt[0],  ers_xyz_dt[1],  ers_xyz_dt[2],  ers_atr_dt[0],  ers_atr_dt[1],  ers_atr_dt[2]),
					String.format("%f,  %f, %f, %f, %f, %f",ers_xyz_d2t[0], ers_xyz_d2t[1], ers_xyz_d2t[2], ers_atr_d2t[0], ers_atr_d2t[1], ers_atr_d2t[2])};
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
		
		public double [] getErsXYZ_dt() {
			return ers_xyz_dt;
		}
		public double [] getErsATR_dt() {
			return ers_atr_dt;
		}
		public double [] getErsXYZ_d2t() {
			return ers_xyz_d2t;
		}
		public double [] getErsATR_d2t() {
			return ers_atr_d2t;
		}
		public void setErsXYZ_dt(double [] d) {
			ers_xyz_dt[0] = d[0];
			ers_xyz_dt[1] = d[1];
			ers_xyz_dt[2] = d[2];
		}
		public void setErsATR_dt(double [] d) {
			ers_atr_dt[0] = d[0];
			ers_atr_dt[1] = d[1];
			ers_atr_dt[2] = d[2];
		}
		public void setErsXYZ_d2t(double [] d) {
			ers_xyz_d2t[0] = d[0];
			ers_xyz_d2t[1] = d[1];
			ers_xyz_d2t[2] = d[2];
		}
		public void setErsATR_d2t(double [] d) {
			ers_atr_d2t[0] = d[0];
			ers_atr_d2t[1] = d[1];
			ers_atr_d2t[2] = d[2];
		}
		
		
		
	}
	
	/**
	 * Calculate rotational matrix and its derivatives by az, tl and rl
	 * @param atr
	 * @return
	 */
	public static Matrix [] getInterRotDeriveMatrices(
			double [] atr,
			boolean invert){
		double ca = Math.cos(atr[0]); // azimuth
		double sa = Math.sin(atr[0]);
		double ct = Math.cos(atr[1]);
		double st = Math.sin(atr[1]);
		double cr = Math.cos(atr[2]);
		double sr = Math.sin(atr[2]);
		double rot_az_sign = invert ? -INTER_ROT_AZ_SGN : INTER_ROT_AZ_SGN; 
		double rot_tl_sign = invert ? -INTER_ROT_TL_SGN : INTER_ROT_TL_SGN; 
		double rot_rl_sign = invert ? -INTER_ROT_RL_SGN : INTER_ROT_RL_SGN; 
		Matrix m_az = new Matrix(new double[][]{
			{  ca,                     0.0,                     sa * rot_az_sign },
			{  0.0,                    1.0,                     0.0},
			{ -sa* rot_az_sign,        0.0,                     ca}});

		Matrix m_tl = new Matrix(new double[][]  {
			{  1.0,                    0.0,                     0.0},
			{  0.0,                    ct,                      st * rot_tl_sign},
			{  0.0,                   -st * rot_tl_sign,   ct}});

		Matrix m_rl = new Matrix(new double[][]  {
			{  cr,                     sr * rot_rl_sign,        0.0},
			{ -sr * rot_rl_sign,       cr,                      0.0},
			{  0.0,                    0.0,                     1.0}});

		Matrix m_daz = new Matrix(new double[][] {
			{ -sa,                     0.0,                     ca * rot_az_sign },
			{  0.0,                    0.0,                     0.0},
			{ -ca* rot_az_sign,        0.0,                    -sa}});

		Matrix m_dtl = new Matrix(new double[][] { // inverted - OK
			{ 0.0,                     0.0,                     0.0},
			{ 0.0,                    -st,                      ct * rot_tl_sign},
			{ 0.0,                    -ct * rot_tl_sign,       -st}});

		Matrix m_drl = new Matrix(new double[][] { // inverted OK
			{ -sr,                     cr * rot_rl_sign,        0.0},
			{ -cr *rot_rl_sign,       -sr,                      0.0},
			{ 0.0,                     0.0,                     0.0}});
		if (invert) {
			return new Matrix[] {
					m_az.times (m_tl.times (m_rl)),
					m_daz.times(m_tl.times (m_rl)),
					m_az.times (m_dtl.times(m_rl)),
					m_az.times (m_tl.times (m_drl))};
		} else {
			return new Matrix[] {
					m_rl.times (m_tl.times (m_az)),
					m_rl.times (m_tl.times (m_daz)),
					m_rl.times (m_dtl.times(m_az)),
					m_drl.times(m_tl.times (m_az))};
		}
	}
	
	public static Vector3D matrixTimesVector(Matrix m, Vector3D v3) {
		return new Vector3D(m.times(new Matrix(v3.toArray(),3)).getColumnPackedCopy());
	}
	
	// use deep=true for the independent instance (clone), false - to "upgrade" GeometryCorrection to ErsCorrection
	public ErsCorrection(GeometryCorrection gc, boolean deep) {
		debugLevel =           gc.debugLevel;
		line_time =            gc.line_time;            // 36.38! //26.5E-6; // duration of sensor scan line (for ERS)
		monochrome =           gc.monochrome;
		lwir =                 gc.lwir;
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
		
//		rXY_ideal =            gc.rXY_ideal; // = {{-0.5, -0.5}, {0.5,-0.5}, {-0.5, 0.5}, {0.5,0.5}};
		set_rXY_ideal(gc.get_rXY_ideal()); //)
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
		camera_heights =       gc.camera_heights; //  =     null; // used to calculate scanline timing
		if (deep) {
			forward =   clone1d(forward);
			right =     clone1d(right);
			height =    clone1d(height);
			roll =      clone1d(roll);
			pXY0 =      clone2d(pXY0);
			XYZ_he =    clone2d(XYZ_he);
			XYZ_her =   clone2d(XYZ_her);
			rXY =       clone2d(rXY);
//			rXY_ideal = clone2d(rXY_ideal);
			set_rXY_ideal(clone2d(gc.get_rXY_ideal())); //)
			rByRDist =  clone1d(rByRDist); // probably it is not needed
			extrinsic_corr = extrinsic_corr.clone(); 
			if (rigOffset!=null) rigOffset = rigOffset.clone();
			woi_tops =  clone1d(woi_tops);
			camera_heights = clone1d(camera_heights);
		}
		resetScenes(); // no scenes yet
		// generate initial ers velocity and roll
		setupERSfromExtrinsics();
	}
	
	
	
	public static double [] clone1d(double [] din){
		if (din == null) return null;
		return din.clone(); 
	}
	public static int [] clone1d(int [] din){
		if (din == null) return null;
		return din.clone(); 
	}
	public static boolean [] clone1d(boolean [] din){
		if (din == null) return null;
		return din.clone(); 
	}
	
	public static double [][] clone2d(double [][] din){
		if (din == null) return null;
		double [][] dout = new double [din.length][];
		for (int i = 0; i < dout.length; i++) {
			dout[i] = clone1d(din[i]);
		}
		return dout;
	}
	public static int [][] clone2d(int [][] din){
		if (din == null) return null;
		int [][] dout = new int [din.length][];
		for (int i = 0; i < dout.length; i++) {
			dout[i] = clone1d(din[i]);
		}
		return dout;
	}

	public static boolean [][] clone2d(boolean [][] din){
		if (din == null) return null;
		boolean [][] dout = new boolean [din.length][];
		for (int i = 0; i < dout.length; i++) {
			dout[i] = clone1d(din[i]);
		}
		return dout;
	}
	
	
	public static double [][][] clone3d(double [][][] din){
		if (din == null) return null;
		double [][][] dout = new double [din.length][][];
		for (int i = 0; i < dout.length; i++) {
			dout[i] = clone2d(din[i]);
		}
		return dout;
	}

	// setup from extrinsics vector
	public void setupERSfromExtrinsics() 
	{
		double BUGS = 1.0;
		if (line_time == OLD_LINE_TIME) {
			BUGS = 26.5/36.38; // temporary correction for the wrong scan line time!
		}
		
		double [] ersv = getCorrVector().getIMU();
		setupERS(
				new double [3],                            // double [] wxyz_center,     // world camera XYZ (meters) for the frame center
				new double [] {BUGS*ersv[3], BUGS*ersv[4], BUGS*ersv[5]}, // double [] wxyz_center_dt,  // world camera Vx, Vy, Vz (m/s)
				new double [3],                            // double [] wxyz_center_d2t, // world camera Vx, Vy, Vz (m/s^2)
//				double [] watr_center,     // camera orientation (az, tilt, roll in radians, corresponding to the frame center)
//				new double [] {ersv[1], ersv[0], ersv[2]}, // double [] watr_center_dt,  // camera rotaions (az, tilt, roll in radians/s, corresponding to the frame center)
				// REVERSING tilt sign !
///				new double [] {ersv[1], -ersv[0], ersv[2]}, // double [] watr_center_dt,  // camera rotaions (az, tilt, roll in radians/s, corresponding to the frame center)
////				new double [] {-ersv[1], -ersv[0], ersv[2]}, // double [] watr_center_dt,  // camera rotaions (az, tilt, roll in radians/s, corresponding to the frame center)
				
				// 2* - because watr_center_dt is now 2 * omega!
//				new double [] {-2*BUGS*ersv[1], -2*BUGS*ersv[0], 2*BUGS*ersv[2]}, // double [] watr_center_dt,  // camera rotaions (az, tilt, roll in radians/s, corresponding to the frame center)
				new double [] {-2*BUGS*ersv[1], -2*BUGS*ersv[0], -2*BUGS*ersv[2]}, // double [] watr_center_dt,  // camera rotaions (az, tilt, roll in radians/s, corresponding to the frame center)
				new double [3]);                           // double [] watr_center_d2t) // camera rotaions (az, tilt, roll in radians/s, corresponding to the frame center)
	
	}
	
	public void setupERS(
			double [] wxyz_center,     // world camera XYZ (meters) for the frame center
			double [] wxyz_center_dt,  // world camera Vx, Vy, Vz (m/s)
			double [] wxyz_center_d2t, // world camera Vx, Vy, Vz (m/s^2)
//			double [] watr_center,     // camera orientation (az, tilt, roll in radians, corresponding to the frame center)
			double [] watr_center_dt,  // camera rotaions (az, tilt, roll in radians/s, corresponding to the frame center)
			double [] watr_center_d2t) // camera rotaions (az, tilt, roll in radians/s, corresponding to the frame center)
	{
		this.ers_wxyz_center =     wxyz_center;     // world camera XYZ (meters) for the lens center (in camera coordinates, typically 0)
		this.ers_wxyz_center_dt =  wxyz_center_dt;  // world camera Vx, Vy, Vz (m/s)
		this.ers_wxyz_center_d2t = wxyz_center_d2t; // world camera Vx, Vy, Vz (m/s^2)
		this.ers_watr_center_dt =  watr_center_dt;  // camera rotaions (az, tilt, roll in radians/s, corresponding to the frame center)
		this.ers_watr_center_d2t = watr_center_d2t; // camera rotaions (az, tilt, roll in radians/s, corresponding to the frame center)
		setupERS();
	}
	
	public void printVectors(
			double [] xyz,
			double [] atr)
	{
		if (xyz != null) 	printAngle("       XYZ",xyz);
		if (atr != null) 	printAngle("       ATR",atr);
		printAngle(                    "ERS XYZ_dt",ers_wxyz_center_dt);
		printAngle(                    "ERS ATR_dt",ers_watr_center_dt);
	}
	
	public void setupERS()
	{
		double ers_sign =      1.0; // -1.0; // invert all corrections to opposite? 
		
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
		Quaternion quat_center1 = new Quaternion (0.0,ers_watr_center_dt[1],   ers_watr_center_dt[0],   ers_watr_center_dt[2]); // angular velocity 1/s :tilt, az, roll
		Quaternion quat_center2 = new Quaternion (0.0,ers_watr_center_d2t[1],  ers_watr_center_d2t[0],  ers_watr_center_d2t[2]); // angular velocity 1/s :tilt, az, roll
		
//		Quaternion dbg_quat_center0 = quat_center0.normalize();
//		Quaternion dbg_quat_center1 = quat_center1.normalize();
//		Quaternion dbg_quat_center2 = quat_center2.normalize(); // not a rotation, can not be normalized
		// integration to the bottom of the image
		double dt = ers_sign*line_time; 
		double [] wxy0 = ers_wxyz_center.clone();
		double [] wxy1 = ers_wxyz_center_dt.clone();
		double [] wxy2 = ers_wxyz_center_d2t.clone();
		// bottom half rotations
		dt =  ers_sign*line_time; 
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
		dt =  -ers_sign*line_time; 
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
		dt = ers_sign*line_time; 
		for (int h = cent_h; h < pixelCorrectionHeight; h ++) {
			for (int i = 0; i < 3; i++) {
				ers_xyz[h][i] =    wxy0[i];
				ers_xyz_dt[h][i] = wxy1[i];
				double wxy1_next = wxy1[i] +  wxy2[i] * dt;
				wxy0[i] +=  0.5*(wxy1[i] + wxy1_next) * dt;
				wxy1[i] = wxy1_next;
			}
		}
		dt = -ers_sign*line_time; 
		wxy0 = ers_wxyz_center.clone();
		wxy1 = ers_wxyz_center_dt.clone();
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
	 * Match other camera px, py, disparity to the reference one
	 * @param px horizontal pixel coordinate (right) of the reference camera view
	 * @param py vertical pixel coordinate (down) of the reference camera view
	 * @param disparity nominal disparity (pixels)  of the reference camera view
	 * @param distortedView true: Radially-distorted reference view, false - rectilinear
	 * @param reference_xyz reference view lens position during centerline acquisition in world coordinates (typically zero), null - use instance global
	 * @param reference_atr reference view orientation during centerline acquisition in world frame (typically zero), null - use instance global
	 * @param distortedCamera true: Radially-distorted camera view, false - rectilinear
	 * @param camera_xyz camera lens position during centerline acquisition in world coordinates, null - use instance global
	 * @param camera_atr camera orientation during centerline acquisition in world frame, null - use instance global
	 * @param line_err iterate until the line (pY) correction is below this value
	 * @return {px, py, disparity } (right, down) or null if behind the camera of the other camera
	 */
	
	public double [] getImageCoordinatesERS(
			QuadCLT cameraQuadCLT, // camera station that got image to be to be matched 
			double px,                // pixel coordinate X in the reference view
			double py,                // pixel coordinate Y in the reference view
			double disparity,         // this reference disparity 
			boolean distortedView,    // This camera view is distorted (diff.rect), false - rectilinear
			double [] reference_xyz,  // this view position in world coordinates (typically zero3)
			double [] reference_atr,  // this view orientation relative to world frame  (typically zero3)
			boolean distortedCamera,  // camera view is distorted (false - rectilinear)
			double [] camera_xyz,     // camera center in world coordinates
			double [] camera_atr,     // camera orientation relative to world frame
			double    line_err)       // threshold error in scan lines (1.0)
	{
		if (reference_xyz == null)	reference_xyz = this.camera_xyz;
		if (reference_atr == null)	reference_atr = this.camera_atr;

		// Find world coordinates of the reference pixel
		double [] xyzw = getWorldCoordinatesERS( // {x - left,y - up, z (0 at camera, negative away), 1} for real, {x,y,z,0} - for infinity
				px,
				py,
				disparity,
				distortedView,   // correct distortion (will need corrected background too !)
				reference_xyz,   // camera center in world coordinates
				reference_atr);  // camera orientation relative to world frame
		if (xyzw == null) {
			return null;
		}
		if (xyzw[2] > 0) {
			xyzw[2] = xyzw[2];
		}
		ErsCorrection ers_camera = this;
		if (cameraQuadCLT != null) {
			ers_camera = cameraQuadCLT.getErsCorrection();
		}
		if (camera_xyz == null)	camera_xyz = ers_camera.camera_xyz;
		if (camera_atr == null)	camera_atr = ers_camera.camera_atr;

		double [] pXpYD = ers_camera.getImageCoordinatesERS( // USED in lwir
				xyzw,
				distortedCamera,
				camera_xyz,  // camera center in world coordinates
				camera_atr,  // camera orientation relative to world frame
				line_err);             // threshold error in scan lines (1.0)

		return pXpYD;
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
	public double [] getWorldCoordinatesERS(
			double px,
			double py,
			double disparity,
			boolean correctDistortions)// correct distortion (will need corrected background too !)
	{
		return getWorldCoordinatesERS(
				px,
				py,
				disparity,
				correctDistortions,// correct distortion (will need corrected background too !)
				camera_xyz, // camera center in world coordinates
				camera_atr);
	}
	
	public double [] getWorldCoordinatesERS(
			double px,
			double py,
			double disparity,
			boolean correctDistortions,// correct distortion (will need corrected background too !)
			double [] camera_xyz, // camera center in world coordinates
			double [] camera_atr)  // camera orientation relative to world frame
	{
		if (Double.isNaN(disparity)) {
			return null;
		}
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
///		boolean is_infinity = Math.abs(disparity) < THRESHOLD; // Maybe all negative - too?
///		if (!is_infinity) {
		if (Math.abs(disparity) >= THRESHOLD) {
			for (int i = 0; i < 3; i++) {
				xyz[i] /= disparity;
			}
		} else {
			if (disparity >= 0) {
				for (int i = 0; i < 3; i++) {
					xyz[i] /= THRESHOLD;
				}
			} else {
				for (int i = 0; i < 3; i++) {
					xyz[i] /= -THRESHOLD;
				}
			}
		}
		Vector3D v3= new Vector3D(xyz);
		// convert to frame parallel to the camera during center line
		Vector3D cam_center_now_local = cam_orient_now_local.applyInverseTo(v3);
		// get real world xyz relative to the camera acquiring a center line
///		Vector3D cam_center_local = (is_infinity) ? cam_center_now_local :  cam_center_now_local.add(cam_now_local); // skip translation for infinity
		Vector3D cam_center_local = cam_center_now_local.add(cam_now_local); // skip translation for infinity
		// undo camera rotation during acquisition of the center line. 
		Vector3D cam_center_world = cam_orient_center.applyInverseTo(cam_center_local);
		// convert to the real world coordinates
///		world_xyz =   (is_infinity) ? cam_center_world : cam_center_world.add(new Vector3D(camera_xyz));
		world_xyz =   cam_center_world.add(new Vector3D(camera_xyz));
		double [] wxyz = world_xyz.toArray();
		double [] wxyz4 = {wxyz[0],wxyz[1],wxyz[2], 1.0};
		if (Double.isNaN(wxyz4[0])) {
			wxyz4[0] = Double.NaN;
		}

///		if (is_infinity) {
///			wxyz4[3] = 0.0;
///		}
		return wxyz4;
	}
	
	/**
	 * Get pixel disparity and coordinates from the real world coordinates (in meters)
	 * @param xyzw real world coordinates {x, y, z, w} in meters (right up, towards camera)
	 * @param correctDistortions true: correct lens distortions, false - no lens distortions
	 * @param camera_xyz camera lens position during centerline acquisition in world coordinates
	 * @param camera_atr camera orientation during centerline acquisition in world frame
	 * @param line_err iterate until the line (pY) correction is below this value
	 * @return {px, py, disparity } (right, down)
	 */
	@Deprecated
	public double [] getImageCoordinatesERS( // not used?
			double [] xyzw,
			boolean correctDistortions, // correct distortion (will need corrected background too !)
			double    line_err) // threshold error in scan lines (1.0)
	{
		return getImageCoordinatesERS( // USED in lwir
				xyzw,
				correctDistortions, // correct distortion (will need corrected background too !)
				camera_xyz, // camera center in world coordinates
				camera_atr,  // camera orientation relative to world frame
				line_err); // threshold error in scan lines (1.0)
	}
	
	public void comparePXYD_Derivatives(
			QuadCLT   scene_QuadClt,
			QuadCLT   reference_QuadClt,
			double    max_inf_disparity, // absolute value
			int       debug_level
			)
	{
		double [] deltas0 = {
				2.0,    // dw_dpX,  (pix)     0
				2.0,    // dw_dpY   (pix)     1
				0.1,    // dw_dd,   (pix)     2
				0.01,   // dw_dvaz, (rad/sec) 3 was 0.2
				0.01,   // dw_dvtl, (rad/sec) 4 was 0.2
				0.02,   // dw_dvrl, (rad/sec) 5 was 0.2
				0.05,   // dw_dvx,  (m/s)     6 was 0.1
				0.05,   // dw_dvy,  (m/s)     7 was 0.1
				0.05,   // dw_dvz,  (m/s)     8 was 0.1
				0.0002, // dw_daz,  (rad)     9 was 0.01
				0.0002, // dw_dtl,  (rad)    10 was 0.01
				0.0002, // dw_drl,  (rad)    11 was 0.01
				0.002,  // dw_dx,   (m)      12 was 0.1
				0.002,  // dw_dy,   (m)	   13 was 0.1
				0.002}; // dw_dz};  (m)      14 was 0.1
		double scale_delta = 1.0; // 0.1; // 1.0; // 0.1; // 0.5;
//		double [] deltas = deltas0.clone();
		double [] deltas = new double [DP_NUM_PARS];
		System.arraycopy(deltas0, 0, deltas, 0,              deltas0.length);
		System.arraycopy(deltas0, 3, deltas, deltas0.length, deltas0.length - 3);
		for (int i = 0; i < deltas.length; i++) deltas[i] *= scale_delta;
	    int dbg_tile = 7508; // 56:23 // 16629;
	    ErsCorrection scene_ers = scene_QuadClt.getErsCorrection();
		TileProcessor tp = reference_QuadClt.getTileProcessor();
		int tilesX = tp.getTilesX();
		int tilesY = tp.getTilesY();
		int tiles = tilesX*tilesY;
		int transform_size = tp.getTileSize();
		double [][] dsrbg_reference = reference_QuadClt.getDSRBG();
		double [][][] derivs_delta =  new double [deltas.length][tiles][];
		double [][][] derivs_true =   new double [deltas.length][tiles][];
		double [][]  pXpYD_original = new double [tiles][];
		double [][]  pXpYD_test =     new double [tiles][];
		double [][]  wxyz =           new double[tiles][];
		boolean []   is_infinity =    new boolean [tiles];
		boolean zero_ers =            false; // true; // true;// false; // true
		boolean fake_linear_ers =     true; // false; //true; // false; // true;
		boolean undefine_outside =    true; // undefine data outside window
		boolean remove_infinity =     false; //true;
		boolean only_infinity =       false; //true; // false;
		
		// ref should have Z more negative (stepping back) to avoid disparity going through infinity
//////	double [] reference_xyz = new double [] {0.0, 0.0,  5.0}; // reference center in world coordinates
/////		double [] reference_xyz = new double [] {-1.0, -0.3,  5.0}; // reference center in world coordinates
////		double [] reference_xyz = new double [] {1.0, 0.3,  -5.0}; // reference center in world coordinates
///		double [] reference_xyz = new double [] {0.0, 0.0,  0.0}; // reference center in world coordinates
//		double [] reference_atr = new double [] {0.0, 0.0,  0.0}; // reference orientation relative to world frame
//		double [] scene_xyz =     new double [] {1.0, 0.3, -5.0}; // scene center in world coordinates
//		double [] scene_atr =     new double [] {.1, .15,   0.2}; // scene orientation relative to world frame
///		double [] scene_xyz =     new double [] {1.0, 0.3,  -5.0}; // scene center in world coordinates
////	double [] scene_xyz =     new double [] {0.0, 0.0,  0.0}; // scene center in world coordinates
//		double [] scene_atr =     new double [] {0.0, 0.0,  0.0}; // scene orientation relative to world frame
//		double [] scene_xyz =     new double [] {-.3,-0.2,  2.0}; // scene center in world coordinates
//		double [] scene_atr =     new double [] {0.05, -0.03,  0.05}; // scene orientation relative to world frame

		double [] reference_xyz = new double [] {0.5, 0.2,-3.0}; // reference center in world coordinates
//		double [] reference_xyz = new double [3]; // reference center in world coordinates
		double [] reference_atr = new double [] {0.1, 0.05,   0.15}; // reference orientation relative to world frame
		double [] scene_xyz =     new double [] {0.2, -0.1,  -2.0}; // reference_xyz.clone(); //new double[3]; // 
		double [] scene_atr =     new double [] {0.05, 0.15, -0.05}; // new double[3]; // reference_atr.clone();
		System.out.println(String.format("reference_xyz= {%8f, %8f, %8f}",reference_xyz[0],reference_xyz[1],reference_xyz[2]));
		System.out.println(String.format("reference_atr= {%8f, %8f, %8f}",reference_atr[0],reference_atr[1],reference_atr[2]));
		System.out.println(String.format("scene_xyz=     {%8f, %8f, %8f}",scene_xyz[0],scene_xyz[1],scene_xyz[2]));
		System.out.println(String.format("scene_atr=     {%8f, %8f, %8f}",scene_atr[0],scene_atr[1],scene_atr[2]));

		System.out.println(String.format("Reference ers_wxyz_center_dt= {%8f, %8f, %8f}",this.ers_wxyz_center_dt[0],this.ers_wxyz_center_dt[1],this.ers_wxyz_center_dt[2]));
		System.out.println(String.format("Reference ers_watr_center_dt= {%8f, %8f, %8f}",this.ers_watr_center_dt[0],this.ers_watr_center_dt[1],this.ers_watr_center_dt[2]));
		System.out.println(String.format("Scene     ers_wxyz_center_dt= {%8f, %8f, %8f}",scene_ers.ers_wxyz_center_dt[0],scene_ers.ers_wxyz_center_dt[1],scene_ers.ers_wxyz_center_dt[2]));
		System.out.println(String.format("Scene     ers_watr_center_dt= {%8f, %8f, %8f}",scene_ers.ers_watr_center_dt[0],scene_ers.ers_watr_center_dt[1],scene_ers.ers_watr_center_dt[2]));
		if (zero_ers) {
			this.ers_wxyz_center_dt=new double[3];
			this.ers_watr_center_dt=new double[3];
			scene_ers.ers_wxyz_center_dt=new double[3];
			scene_ers.ers_watr_center_dt=new double[3];
			System.out.println("Updated:");
			System.out.println(String.format("Reference ers_wxyz_center_dt= {%8f, %8f, %8f}",this.ers_wxyz_center_dt[0],this.ers_wxyz_center_dt[1],this.ers_wxyz_center_dt[2]));
			System.out.println(String.format("Reference ers_watr_center_dt= {%8f, %8f, %8f}",this.ers_watr_center_dt[0],this.ers_watr_center_dt[1],this.ers_watr_center_dt[2]));
			System.out.println(String.format("Scene     ers_wxyz_center_dt= {%8f, %8f, %8f}",scene_ers.ers_wxyz_center_dt[0],scene_ers.ers_wxyz_center_dt[1],scene_ers.ers_wxyz_center_dt[2]));
			System.out.println(String.format("Scene     ers_watr_center_dt= {%8f, %8f, %8f}",scene_ers.ers_watr_center_dt[0],scene_ers.ers_watr_center_dt[1],scene_ers.ers_watr_center_dt[2]));
		}
		if (fake_linear_ers) {
			this.ers_wxyz_center_dt=      new double[] { 1.0, 0.25, -2.0}; // this.ers_watr_center_dt.clone(); // new double[3];
			scene_ers.ers_wxyz_center_dt= new double[] { 0.5, 0.3,  -1.5}; // new double[3];
//			this.ers_wxyz_center_dt=      new double[3];// { 1.0, 0.25, -2.0}; // this.ers_watr_center_dt.clone(); // new double[3];
//			scene_ers.ers_wxyz_center_dt= new double[3];// { 0.5, 0.3,  -1.5}; // new double[3];
			System.out.println("Updated2:");
			System.out.println(String.format("Reference ers_wxyz_center_dt= {%8f, %8f, %8f}",this.ers_wxyz_center_dt[0],this.ers_wxyz_center_dt[1],this.ers_wxyz_center_dt[2]));
			System.out.println(String.format("Reference ers_watr_center_dt= {%8f, %8f, %8f}",this.ers_watr_center_dt[0],this.ers_watr_center_dt[1],this.ers_watr_center_dt[2]));
			System.out.println(String.format("Scene     ers_wxyz_center_dt= {%8f, %8f, %8f}",scene_ers.ers_wxyz_center_dt[0],scene_ers.ers_wxyz_center_dt[1],scene_ers.ers_wxyz_center_dt[2]));
			System.out.println(String.format("Scene     ers_watr_center_dt= {%8f, %8f, %8f}",scene_ers.ers_watr_center_dt[0],scene_ers.ers_watr_center_dt[1],scene_ers.ers_watr_center_dt[2]));
		}

		
		Matrix [] reference_matrices_inverse = getInterRotDeriveMatrices(
				reference_atr, // double [] atr);
				true);         // boolean invert));

		Matrix [] scene_matrices_inverse = getInterRotDeriveMatrices(
				scene_atr, // double [] atr);
				true);         // boolean invert));

		Matrix  scene_rot_matrix = getInterRotDeriveMatrices(
				scene_atr, // double [] atr);
				false)[0]; // boolean invert));
		
		this.setupERS(); // just in case - setup using instance parameters
		scene_ers.setupERS(); // just in case - setup using instance parameters
		
		// measure derivatives by px,py,disparity
		for (int tileY = 0; tileY < tilesY; tileY++) {
			for (int tileX = 0; tileX < tilesX; tileX++) {
				int nTile = tileX + tileY * tilesX;
				double centerX = tileX * transform_size + transform_size/2; //  - shiftX;
				double centerY = tileY * transform_size + transform_size/2; //  - shiftY;
				double disparity = dsrbg_reference[QuadCLT.DSRBG_DISPARITY][nTile];
				
				//DEBUG!!!
//				if (nTile != 8156) {
					if (disparity < max_inf_disparity) { // or only for abs() ?
						disparity = 0.0;
						is_infinity[nTile] = true;
					} 
//				}
				pXpYD_original[nTile] = new double[] {centerX, centerY,disparity};
				// find derivatives by pX, pY, pZ
				if (nTile == dbg_tile) {
					System.out.println("comparePXYD_Derivatives()_0: nTile = "+nTile+" ("+tileX+"/"+tileY+")");
				}

				Vector3D[] wxyz_derivs = getDWorldDPixels(
						true,                          // boolean correctDistortions,
						is_infinity[nTile],            // boolean is_infinity,
						new double[] {centerX, centerY, disparity}, //  double [] pXpYD)
						reference_xyz,                 // double [] camera_xyz, // camera center in world coordinates
						reference_matrices_inverse,    // Matrix [] rot_deriv_matrices_inverse, Array of 4 matrices - 
						//                             inverse transformation matrix and 3 their derivatives by az, tl, rl 
						(nTile == dbg_tile)? 2:0);
				// store world xyz to be used for all scne parameters
				if ((wxyz_derivs != null) && (wxyz_derivs[0] != null)) {
					wxyz[nTile] = wxyz_derivs[0].toArray();
				}
				
				double [][][] xyz_pm = {
						{	getWorldCoordinatesERS(centerX + deltas[DP_DPX], centerY, disparity, true, reference_xyz, reference_atr),
							getWorldCoordinatesERS(centerX - deltas[DP_DPX], centerY, disparity, true, reference_xyz, reference_atr)
						},{ getWorldCoordinatesERS(centerX, centerY + deltas[DP_DPY], disparity, true, reference_xyz, reference_atr),
							getWorldCoordinatesERS(centerX,	centerY - deltas[DP_DPY], disparity, true, reference_xyz, reference_atr)
						},{ getWorldCoordinatesERS(centerX, centerY, disparity + deltas[DP_DD],	 true, reference_xyz, reference_atr),
							getWorldCoordinatesERS(centerX, centerY, disparity - deltas[DP_DD],	 true, reference_xyz, reference_atr)}};
				for (int n = 0; n < 3; n++) {
					int par_index =DP_DPX + n; // DP_DPX, DP_DPY, DP_DD
					if ((xyz_pm[n][0] != null) && (xyz_pm[n][1] != null)) {
						double [][] pm = new double [2][];
						for (int d = 0; d <2; d++){
							double [] xyz3 = xyz_pm[n][d];
							if (xyz3 != null) {
								double [] wxyz4 = {xyz3[0], xyz3[1], xyz3[2], is_infinity[nTile]? 0.0: 1.0};
								pm[d] = scene_ers.getImageCoordinatesERS(
										wxyz4,             // double [] xyzw,
										true,              // boolean correctDistortions, // correct distortion (will need corrected background too !)
										scene_xyz,         // double [] camera_xyz, // camera center in world coordinates
										scene_rot_matrix,  // Matrix    rot_matrix,  // 1-st of 4 matrices 
										LINE_ERR);         // double    line_err) // threshold error in scan lines (1.0)
							}
						}
						if ((xyz_pm[n][0] != null) && (xyz_pm[n][1] != null)) {
							derivs_delta[par_index][nTile] = new double[3];
							for (int i = 0; i < 3; i++) {
								derivs_delta[par_index][nTile][i] = (pm[0][i] - pm[1][i]) / (2.0 * deltas[par_index]); 
							}
						}
					}
				}
				// find derivatives by 	reference_xyz,	reference_atr; that do not require re-programming of the ERS arrays
				double [][][] atrxyz_pm = new double [6][2][];
				for (int n = 0; n < 3; n++ ){
					int par_index =DP_DAZ + n;
					double [] cam_atr = reference_atr.clone();
					cam_atr[n] = reference_atr[n] + deltas[par_index];
					atrxyz_pm[n][0] = getWorldCoordinatesERS(centerX, centerY, disparity, true, reference_xyz, cam_atr);
					cam_atr[n] = reference_atr[n] - deltas[par_index];
					atrxyz_pm[n][1] = getWorldCoordinatesERS(centerX, centerY, disparity, true, reference_xyz, cam_atr);
					par_index =DP_DX + n;
					double [] cam_xyz = reference_xyz.clone();
					cam_xyz[n] = reference_xyz[n] + deltas[par_index];
					atrxyz_pm[n+3][0] = getWorldCoordinatesERS(centerX, centerY, disparity, true, cam_xyz, reference_atr);
					cam_xyz[n] = reference_xyz[n] - deltas[par_index];
					atrxyz_pm[n+3][1] = getWorldCoordinatesERS(centerX, centerY, disparity, true, cam_xyz, reference_atr);
				}
				for (int n = 0; n < 6; n++ ){
					int par_index =DP_DAZ + n;
					if ((atrxyz_pm[n][0] != null) && (atrxyz_pm[n][1] != null)) {
						double [][] pm = new double [2][];
						for (int d = 0; d < 2; d++) {
							double [] xyz3 = atrxyz_pm[n][d];
							if (xyz3 != null) {
								double [] wxyz4 = {xyz3[0], xyz3[1], xyz3[2], is_infinity[nTile]? 0.0: 1.0};
								pm[d] = scene_ers.getImageCoordinatesERS(
										wxyz4,             // double [] xyzw,
										true,              // boolean correctDistortions, // correct distortion (will need corrected background too !)
										scene_xyz,         // double [] camera_xyz, // camera center in world coordinates
										scene_rot_matrix,  // Matrix    rot_matrix,  // 1-st of 4 matrices 
										LINE_ERR);         // double    line_err) // threshold error in scan lines (1.0)
							}
						}
						derivs_delta[par_index][nTile] = new double[3];
						for (int i = 0; i < 3; i++) {
							derivs_delta[par_index][nTile][i] = (pm[0][i] - pm[1][i]) / (2.0 * deltas[par_index]); 
						}
					}
				}					
			}
		}
	
		double [][][][] ers_ref_pm = new double [6][2][tiles][];
		double [] ers_ref_watr_center_dt_original = this.ers_watr_center_dt.clone();
		double [] ers_ref_wxyz_center_dt_original = this.ers_wxyz_center_dt.clone();
		for (int nf = 0; nf < 12; nf++) {
			int n = nf/2;
			int sign =   ((nf & 1) != 0) ? -1 : 1;
			int sign01 = ((nf & 1) != 0) ?  1 : 0;
			boolean lin_not_ang = n >= 3;
			int component = n % 3;
			int par_index =DP_DVAZ + n;
			this.ers_watr_center_dt = ers_ref_watr_center_dt_original.clone();
			this.ers_wxyz_center_dt = ers_ref_wxyz_center_dt_original.clone();
			if (lin_not_ang) {
				this.ers_wxyz_center_dt[component] = ers_ref_wxyz_center_dt_original[component] + sign * deltas[par_index];
			} else {
				this.ers_watr_center_dt[component] = ers_ref_watr_center_dt_original[component] + sign * deltas[par_index];
			}
			this.setupERS();
			for (int tileY = 0; tileY < tilesY; tileY++) {
				for (int tileX = 0; tileX < tilesX; tileX++) {
					int nTile = tileX + tileY * tilesX;
					double centerX = tileX * transform_size + transform_size/2; //  - shiftX;
					double centerY = tileY * transform_size + transform_size/2; //  - shiftY;
					double disparity = dsrbg_reference[QuadCLT.DSRBG_DISPARITY][nTile];
					if (is_infinity[nTile]) { // or only for abs() ?
						disparity = 0.0;
					}
					double [] xyz3 =this.getWorldCoordinatesERS(centerX, centerY, disparity, true, reference_xyz, reference_atr);
					// convert to scene pixels
					if (xyz3!= null) {
						double [] wxyz4 = {xyz3[0], xyz3[1], xyz3[2], is_infinity[nTile]? 0.0: 1.0};
						ers_ref_pm[n][sign01][nTile]= scene_ers.getImageCoordinatesERS(
								wxyz4,             // double [] xyzw,
								true,              // boolean correctDistortions, // correct distortion (will need corrected background too !)
								scene_xyz,         // double [] camera_xyz, // camera center in world coordinates
								scene_rot_matrix,  // Matrix    rot_matrix,  // 1-st of 4 matrices 
								LINE_ERR);         // double    line_err) // threshold error in scan lines (1.0)
					}
				}
			}
		}
		// restore original ERS data
		this.ers_watr_center_dt = ers_ref_watr_center_dt_original.clone();
		this.ers_wxyz_center_dt = ers_ref_wxyz_center_dt_original.clone();
		this.setupERS();
		for (int n = 0; n < 6; n++) {
			int par_index =DP_DVAZ + n;
			for (int nTile = 0; nTile < tiles; nTile++) {
				if ((ers_ref_pm[n][0][nTile] != null) && (ers_ref_pm[n][1][nTile] != null)) {
					derivs_delta[par_index][nTile] = new double[3];
					for (int i = 0; i <3; i++) {
						derivs_delta[par_index][nTile][i] = (ers_ref_pm[n][0][nTile][i] - ers_ref_pm[n][1][nTile][i]) / (2.0 * deltas[par_index]);
					}
				}
			}
		}
//Use wxyz for all the scene parameters 		
/// Now derivatives by scene parameters
		// measure derivatives by scene_xyz,	scene_atr
		for (int tileY = 0; tileY < tilesY; tileY++) {
			for (int tileX = 0; tileX < tilesX; tileX++) {
				int nTile = tileX + tileY * tilesX;
				// find derivatives by 	reference_xyz,	reference_atr; that do not require re-programming of the ERS arrays
				double [][][] atrxyz_pm = new double [6][2][];
				for (int n = 0; n < 3; n++ ){
					double [] xyz3 = wxyz[nTile];
					if (xyz3 != null) {
						double [] wxyz4 = {xyz3[0], xyz3[1], xyz3[2], is_infinity[nTile]? 0.0: 1.0};
						int par_index =DP_DSAZ + n;
						double [] cam_atr = scene_atr.clone();
						cam_atr[n] = scene_atr[n] + deltas[par_index];
						atrxyz_pm[n][0] = scene_ers.getImageCoordinatesERS(
								wxyz4,             // double [] xyzw,
								true,              // boolean correctDistortions, // correct distortion (will need corrected background too !)
								scene_xyz,         // double [] camera_xyz, // camera center in world coordinates
								cam_atr,           // double [] camera_atr,  // camera orientation relative to world frame 
								LINE_ERR);         // double    line_err) // threshold error in scan lines (1.0)
						cam_atr[n] = scene_atr[n] - deltas[par_index];
						atrxyz_pm[n][1] = scene_ers.getImageCoordinatesERS(
								wxyz4,             // double [] xyzw,
								true,              // boolean correctDistortions, // correct distortion (will need corrected background too !)
								scene_xyz,         // double [] camera_xyz, // camera center in world coordinates
								cam_atr,           // double [] camera_atr,  // camera orientation relative to world frame 
								LINE_ERR);         // double    line_err) // threshold error in scan lines (1.0)
						par_index =DP_DSX + n;
						double [] cam_xyz = scene_xyz.clone();
						cam_xyz[n] = scene_xyz[n] + deltas[par_index];
						atrxyz_pm[n+3][0] = scene_ers.getImageCoordinatesERS(
								wxyz4,             // double [] xyzw,
								true,              // boolean correctDistortions, // correct distortion (will need corrected background too !)
								cam_xyz,         // double [] camera_xyz, // camera center in world coordinates
								scene_atr,           // double [] camera_atr,  // camera orientation relative to world frame 
								LINE_ERR);         // double    line_err) // threshold error in scan lines (1.0)
						cam_xyz[n] = scene_xyz[n] - deltas[par_index];
						atrxyz_pm[n+3][1] = scene_ers.getImageCoordinatesERS(
								wxyz4,             // double [] xyzw,
								true,              // boolean correctDistortions, // correct distortion (will need corrected background too !)
								cam_xyz,         // double [] camera_xyz, // camera center in world coordinates
								scene_atr,           // double [] camera_atr,  // camera orientation relative to world frame 
								LINE_ERR);         // double    line_err) // threshold error in scan lines (1.0)
					}
				}
				for (int n = 0; n < 6; n++ ){
					int par_index =DP_DSAZ + n;
					if ((atrxyz_pm[n][0] != null) && (atrxyz_pm[n][1] != null)) {
						derivs_delta[par_index][nTile] = new double[3];
						for (int i = 0; i < 3; i++) {
							derivs_delta[par_index][nTile][i] = (atrxyz_pm[n][0][i] - atrxyz_pm[n][1][i]) / (2.0 * deltas[par_index]); 
						}
					}
				}					
			}
		}
		/// Now derivatives by scene ERS parameters
		double [][][][] ers_scene_pm = new double [6][2][tiles][];
		double [] ers_scene_watr_center_dt_original = scene_ers.ers_watr_center_dt.clone();
		double [] ers_scene_wxyz_center_dt_original = scene_ers.ers_wxyz_center_dt.clone();
		for (int nf = 0; nf < 12; nf++) {
			int n = nf/2;
			int sign =   ((nf & 1) != 0) ? -1 : 1;
			int sign01 = ((nf & 1) != 0) ?  1 : 0;
			boolean lin_not_ang = n >= 3;
			int component = n % 3;
			int par_index =DP_DSVAZ + n;
			scene_ers.ers_watr_center_dt = ers_scene_watr_center_dt_original.clone();
			scene_ers.ers_wxyz_center_dt = ers_scene_wxyz_center_dt_original.clone();
			if (lin_not_ang) {
				scene_ers.ers_wxyz_center_dt[component] = ers_scene_wxyz_center_dt_original[component] + sign * deltas[par_index];
			} else {
				scene_ers.ers_watr_center_dt[component] = ers_scene_watr_center_dt_original[component] + sign * deltas[par_index];
			}
			scene_ers.setupERS();
			for (int tileY = 0; tileY < tilesY; tileY++) {
				for (int tileX = 0; tileX < tilesX; tileX++) {
					int nTile = tileX + tileY * tilesX;
					double [] xyz3 = wxyz[nTile];
					// convert to scene pixels
					if (xyz3!= null) {
						double [] wxyz4 = {xyz3[0], xyz3[1], xyz3[2], is_infinity[nTile]? 0.0: 1.0};
						ers_scene_pm[n][sign01][nTile]= scene_ers.getImageCoordinatesERS(
								wxyz4,             // double [] xyzw,
								true,              // boolean correctDistortions, // correct distortion (will need corrected background too !)
								scene_xyz,         // double [] camera_xyz, // camera center in world coordinates
								scene_rot_matrix,  // Matrix    rot_matrix,  // 1-st of 4 matrices 
								LINE_ERR);         // double    line_err) // threshold error in scan lines (1.0)
					}
				}
			}
		}
		// restore original ERS data
		scene_ers.ers_watr_center_dt = ers_scene_watr_center_dt_original.clone();
		scene_ers.ers_wxyz_center_dt = ers_scene_wxyz_center_dt_original.clone();
		scene_ers.setupERS();
		for (int n = 0; n < 6; n++) {
			int par_index =DP_DSVAZ + n;
			for (int nTile = 0; nTile < tiles; nTile++) {
				if ((ers_scene_pm[n][0][nTile] != null) && (ers_scene_pm[n][1][nTile] != null)) {
					derivs_delta[par_index][nTile] = new double[3];
					for (int i = 0; i <3; i++) {
						derivs_delta[par_index][nTile][i] = (ers_scene_pm[n][0][nTile][i] - ers_scene_pm[n][1][nTile][i]) / (2.0 * deltas[par_index]);
					}
				}
			}
		}
// get real derivatives		
// Display results		
		for (int tileY = 0; tileY < tilesY; tileY++) {
			for (int tileX = 0; tileX < tilesX; tileX++) {
				int nTile = tileX + tileY * tilesX;
				double centerX = tileX * transform_size + transform_size/2; //  - shiftX;
				double centerY = tileY * transform_size + transform_size/2; //  - shiftY;
				double disparity = dsrbg_reference[QuadCLT.DSRBG_DISPARITY][nTile];
				if (is_infinity[nTile]) { // or only for abs() ?
					disparity = 0.0;
				}
				if (nTile == dbg_tile) {
					System.out.println("comparePXYD_Derivatives(): nTile = "+nTile+" ("+tileX+"/"+tileY+"), is_infinity="+is_infinity[nTile]);
//					double [][] deriv_params = 
					getDPxSceneDParameters(
							scene_ers,                                  // ErsCorrection    ers_scene,
							true,                                       // boolean correctDistortions,
							!is_infinity [nTile],                        // boolean is_infinity,
							new double[] {centerX, centerY, disparity}, // double [] pXpYD_reference,
							reference_xyz,                              // double [] reference_xyz, // reference (this) camera center in world coordinates
							scene_xyz,                                  // double [] scene_xyz,    // (other, ers_scene) camera center in world coordinates
							reference_matrices_inverse,                 // Matrix [] reference_matrices_inverse,
							scene_matrices_inverse,                     // Matrix [] scene_matrices_inverse,
							scene_rot_matrix,                           // Matrix    scene_rot_matrix,        // single rotation (direct) matrix for the scene
							(nTile == dbg_tile)? 2:0);                  // int debug_level);					
					System.out.println("comparePXYD_Derivatives(): nTile = "+nTile+" ("+tileX+"/"+tileY+"), is_infinity="+is_infinity[nTile]);
				}
//				double [] xyz3 = wxyz[nTile];
				double [][] deriv_params = 	getDPxSceneDParameters(
						scene_ers,                                  // ErsCorrection    ers_scene,
						true,                                       // boolean correctDistortions,
						is_infinity [nTile],                        // boolean is_infinity,
						new double[] {centerX, centerY, disparity}, // double [] pXpYD_reference,
						reference_xyz,                              // double [] reference_xyz, // reference (this) camera center in world coordinates
						scene_xyz,                                  // double [] scene_xyz,    // (other, ers_scene) camera center in world coordinates
						reference_matrices_inverse,                 // Matrix [] reference_matrices_inverse,
						scene_matrices_inverse,                     // Matrix [] scene_matrices_inverse,
						scene_rot_matrix,                           // Matrix    scene_rot_matrix,        // single rotation (direct) matrix for the scene
						(nTile == dbg_tile)? 2:0);                  // int debug_level);
				if (deriv_params != null) {
					pXpYD_test[nTile] = deriv_params[0]; // it should not match if scene != reference, but be close
					for (int i = 0; i < derivs_true.length; i++) {
						if (deriv_params[i+1] != null) {
							derivs_true[i][nTile] = deriv_params[i+1];
						}
					}
				}
			}
		}
		if (undefine_outside) {
			System.out.println("Removing all tiles that fall outside of the FoV of the scene");
			for (int nTile = 0; nTile < tiles; nTile++) {
				boolean bad_tile = true;
				if (pXpYD_test[nTile] != null){
					double pX = pXpYD_test[nTile][0];
					double pY = pXpYD_test[nTile][1];
					if ((pX >=0.0) && (pY >= 0.0) && (pX <= pixelCorrectionWidth) && (pY <= pixelCorrectionHeight)) {
						bad_tile = false;
					}
				}
				if (is_infinity[nTile] && remove_infinity) {
					bad_tile = true;
				}
				if (!is_infinity[nTile] && only_infinity) {
					bad_tile = true;
				}
				if (bad_tile ) {
					pXpYD_test[nTile] = null;
					for (int n = 0; n < derivs_true.length; n++) {
						derivs_true[n][nTile] = null;
						derivs_delta[n][nTile] = null;
					}
				}				
			}
		}
		
		
// display results		
		int show_gap_h = 6; // 324+6 = 330 
		int show_gap_v = 8; // 242 + 8 = 250
		int tilesX1 = tilesX + show_gap_h; 
		int tilesY1 = tilesY + show_gap_v;
		int tilesX2 = 2*tilesX + show_gap_h; 
		int tilesY2 = 2*tilesY + show_gap_v;
		int tiles2 = tilesX2 * tilesY2; 
		if (debug_level > 0) {
			String title3 =     scene_QuadClt.getImageName()+"-pXpYD_ref_scene";
			String [] titles3 = {"reference","scene","ref-scene_diff"};
			double [][] dbg_img3 = new double [titles3.length][tiles2];
			for (int i = 0; i < dbg_img3.length; i++) {
				Arrays.fill(dbg_img3[i], Double.NaN);
			}

			for (int tileY = 0; tileY < tilesY; tileY++) {
				for (int tileX = 0; tileX < tilesX; tileX++) {
					int nTile = tileX + tileY * tilesX;
					if (pXpYD_original[nTile] != null) {
						dbg_img3[0][tileX + 0 * tilesX1 + tilesX2 *(tileY + 0 * tilesY1)] = pXpYD_original[nTile][0];
						dbg_img3[0][tileX + 1 * tilesX1 + tilesX2 *(tileY + 0 * tilesY1)] = pXpYD_original[nTile][1];
						dbg_img3[0][tileX + 0 * tilesX1 + tilesX2 *(tileY + 1 * tilesY1)] = pXpYD_original[nTile][2];
						if (pXpYD_test[nTile] != null) {
							dbg_img3[1][tileX + 0 * tilesX1 + tilesX2 *(tileY + 0 * tilesY1)] = pXpYD_test[nTile][0];
							dbg_img3[1][tileX + 1 * tilesX1 + tilesX2 *(tileY + 0 * tilesY1)] = pXpYD_test[nTile][1];
							dbg_img3[1][tileX + 0 * tilesX1 + tilesX2 *(tileY + 1 * tilesY1)] = pXpYD_test[nTile][2];

							dbg_img3[2][tileX + 0 * tilesX1 + tilesX2 *(tileY + 0 * tilesY1)] = pXpYD_test[nTile][0] - pXpYD_original[nTile][0];
							dbg_img3[2][tileX + 1 * tilesX1 + tilesX2 *(tileY + 0 * tilesY1)] = pXpYD_test[nTile][1] - pXpYD_original[nTile][1];
							dbg_img3[2][tileX + 0 * tilesX1 + tilesX2 *(tileY + 1 * tilesY1)] = pXpYD_test[nTile][2] - pXpYD_original[nTile][2];
						}
					}
				}
			}
			(new ShowDoubleFloatArrays()).showArrays(
					dbg_img3,
					tilesX2,
					tilesY2,
					true,
					title3,
					titles3);
		}
		
		if (debug_level > 0) {
			String title_deriv =     scene_QuadClt.getImageName()+"-scene-ref-derivatives_test";
			String [] titles_deriv = new String [3 * DP_DERIV_NAMES.length];
			for (int i = 0; i < DP_DERIV_NAMES.length; i++) {
				titles_deriv[3 * i + 0] = DP_DERIV_NAMES[i]+"-deriv";
				titles_deriv[3 * i + 1] = DP_DERIV_NAMES[i]+"-delta";
				titles_deriv[3 * i + 2] = DP_DERIV_NAMES[i]+"-diff";
			}
			double [][] dbg_img = new double [titles_deriv.length][tiles2];
			for (int i = 0; i < dbg_img.length; i++) {
				Arrays.fill(dbg_img[i], Double.NaN);
			}
			for (int tileY = 0; tileY < tilesY; tileY++) {
				for (int tileX = 0; tileX < tilesX; tileX++) {
					int nTile = tileX + tileY * tilesX;
					for (int n = 0; n < DP_DERIV_NAMES.length; n++) {
						if (derivs_true[n][nTile] != null) {
							dbg_img[3 * n + 0][tileX + 0 * tilesX1 + tilesX2 *(tileY + 0 * tilesY1)] = derivs_true[n][nTile][0];
							dbg_img[3 * n + 0][tileX + 1 * tilesX1 + tilesX2 *(tileY + 0 * tilesY1)] = derivs_true[n][nTile][1];
							dbg_img[3 * n + 0][tileX + 0 * tilesX1 + tilesX2 *(tileY + 1 * tilesY1)] = derivs_true[n][nTile][2];
						}
						if (derivs_delta[n][nTile] != null) {
							dbg_img[3 * n + 1][tileX + 0 * tilesX1 + tilesX2 *(tileY + 0 * tilesY1)] = derivs_delta[n][nTile][0];
							dbg_img[3 * n + 1][tileX + 1 * tilesX1 + tilesX2 *(tileY + 0 * tilesY1)] = derivs_delta[n][nTile][1];
							dbg_img[3 * n + 1][tileX + 0 * tilesX1 + tilesX2 *(tileY + 1 * tilesY1)] = derivs_delta[n][nTile][2];
						}
						if ((derivs_true[n][nTile] != null) && (derivs_delta[n][nTile] != null)) {
							dbg_img[3 * n + 2][tileX + 0 * tilesX1 + tilesX2 *(tileY + 0 * tilesY1)] = derivs_true[n][nTile][0] - derivs_delta[n][nTile][0];
							dbg_img[3 * n + 2][tileX + 1 * tilesX1 + tilesX2 *(tileY + 0 * tilesY1)] = derivs_true[n][nTile][1] - derivs_delta[n][nTile][1];
							dbg_img[3 * n + 2][tileX + 0 * tilesX1 + tilesX2 *(tileY + 1 * tilesY1)] = derivs_true[n][nTile][2] - derivs_delta[n][nTile][2];
						}
					}
				}
			}
			(new ShowDoubleFloatArrays()).showArrays(
					dbg_img,
					tilesX2,
					tilesY2,
					true,
					title_deriv,
					titles_deriv);
			System.out.println("comparePXYD_Derivatives() DONE.");
			
		}
	}
	
	public void compareDSItoWorldDerivatives(
			QuadCLT   scene_QuadClt,
			double    max_inf_disparity, // absolute value
			int       debug_level
			)
	{
		double [] deltas0 = {
				2.0,  // dw_dpX,  (pix)     0
				2.0,  // dw_dpY   (pix)     1
				0.1,  // dw_dd,   (pix)     2
				0.2,  // dw_dvaz, (rad/sec) 3
				0.2,  // dw_dvtl, (rad/sec) 4
				0.2,  // dw_dvrl, (rad/sec) 5
				0.1,  // dw_dvx,  (m/s)     6
				0.1,  // dw_dvy,  (m/s)     7
				0.1,  // dw_dvz,  (m/s)     8
				0.01, // dw_daz,  (rad)     9
				0.01, // dw_dtl,  (rad)    10
				0.01, // dw_drl,  (rad)    11
				0.1,  // dw_dx,   (m)      12
				0.1,  // dw_dy,   (m)	   13
				0.1}; // dw_dz};  (m)      14
		double scale_delta = 1.0; // 0.1; // 0.5;
		double [] deltas = deltas0.clone();
		for (int i = 0; i < deltas.length; i++) deltas[i] *= scale_delta;
		
		
	    int dbg_tile = 16629;
		TileProcessor tp = scene_QuadClt.getTileProcessor();
		int tilesX = tp.getTilesX();
		int tilesY = tp.getTilesY();
		int tiles = tilesX*tilesY;
		int transform_size = tp.getTileSize();
		double [][] dsrbg_camera =    scene_QuadClt.getDSRBG();
		double [][][] derivs_delta = new double [deltas.length][tiles][];
		double [][][] derivs_true = new double [deltas.length][tiles][];
		double [][]  pXpYD_original = new double [tiles][];
		double [][]  pXpYD_test =     new double [tiles][];
		double [][]  wxyz = new double[tiles][];
		boolean zero_ers = false; // true
		boolean fake_linear_ers = true;
		
		double [] camera_xyz = new double [] {1.0, 0.3, -5.0}; // camera center in world coordinates
		double [] camera_atr = new double [] {.1, .15,   0.2};  // camera orientation relative to world frame
//		double [] camera_atr = new double [] {.2, 0.0,   0.0};  // camera orientation relative to world frame
//		double [] camera_atr = new double [] {0.0, 0.2,  0.0};  // camera orientation relative to world frame
//		double [] camera_atr = new double [] {0.0, 0.0,  0.2};  // camera orientation relative to world frame

		System.out.println(String.format("ers_wxyz_center_dt= {%8f, %8f, %8f}",ers_wxyz_center_dt[0],ers_wxyz_center_dt[1],ers_wxyz_center_dt[2]));
		System.out.println(String.format("ers_watr_center_dt= {%8f, %8f, %8f}",ers_watr_center_dt[0],ers_watr_center_dt[1],ers_watr_center_dt[2]));
		if (zero_ers) {
			ers_wxyz_center_dt=new double[3];
			ers_watr_center_dt=new double[3];
			System.out.println("Updated:");
			System.out.println(String.format("ers_wxyz_center_dt= {%8f, %8f, %8f}",ers_wxyz_center_dt[0],ers_wxyz_center_dt[1],ers_wxyz_center_dt[2]));
			System.out.println(String.format("ers_watr_center_dt= {%8f, %8f, %8f}",ers_watr_center_dt[0],ers_watr_center_dt[1],ers_watr_center_dt[2]));
		}
		if (fake_linear_ers) {
			ers_wxyz_center_dt=new double[] {1.0, 0.25, -2.0}; ers_watr_center_dt.clone(); // new double[3];
//			ers_watr_center_dt=new double[3];
			System.out.println("Updated:");
			System.out.println(String.format("ers_wxyz_center_dt= {%8f, %8f, %8f}",ers_wxyz_center_dt[0],ers_wxyz_center_dt[1],ers_wxyz_center_dt[2]));
			System.out.println(String.format("ers_watr_center_dt= {%8f, %8f, %8f}",ers_watr_center_dt[0],ers_watr_center_dt[1],ers_watr_center_dt[2]));
		}

		
		Matrix [] rot_deriv_matrices_inverse = getInterRotDeriveMatrices(
				camera_atr, // double [] atr);
				true); // boolean invert));

		Matrix  rot_matrix = getInterRotDeriveMatrices(
				camera_atr, // double [] atr);
				false)[0]; // boolean invert));
		
		setupERS(); // just in case - setup using instance parameters
		// measure derivatives by px,py,disparity
		for (int tileY = 0; tileY < tilesY; tileY++) {
			for (int tileX = 0; tileX < tilesX; tileX++) {
				int nTile = tileX + tileY * tilesX;
				double centerX = tileX * transform_size + transform_size/2; //  - shiftX;
				double centerY = tileY * transform_size + transform_size/2; //  - shiftY;
				double disparity = dsrbg_camera[QuadCLT.DSRBG_DISPARITY][nTile];
				if (disparity < max_inf_disparity) { // or only for abs() ?
					disparity = 0.0;
				}
				
				pXpYD_original[nTile] = new double[] {centerX, centerY,disparity};
				// find derivatives by pX, pY, pZ
				if (nTile == dbg_tile) {
					System.out.println("compareDSItoWorldDerivatives()_0: nTile = "+nTile+" ("+tileX+"/"+tileY+")");
				}

				double [][][] xyz_pm = {
						{	getWorldCoordinatesERS(centerX + deltas[DW_DPX], centerY, disparity, true, camera_xyz, camera_atr),
							getWorldCoordinatesERS(centerX - deltas[DW_DPX], centerY, disparity, true, camera_xyz, camera_atr)
						},{ getWorldCoordinatesERS(centerX, centerY + deltas[DW_DPY], disparity, true, camera_xyz, camera_atr),
							getWorldCoordinatesERS(centerX,	centerY - deltas[DW_DPY], disparity, true, camera_xyz, camera_atr)
						},{ getWorldCoordinatesERS(centerX, centerY, disparity + deltas[DW_DD],	 true, camera_xyz, camera_atr),
							getWorldCoordinatesERS(centerX, centerY, disparity - deltas[DW_DD],	 true, camera_xyz, camera_atr)}};
				for (int n = 0; n < 3; n++) {
					int par_index =DW_DPX + n;
					if ((xyz_pm[n][0] != null) && (xyz_pm[n][1] != null)) {
						derivs_delta[par_index][nTile] = new double[3];
						for (int i = 0; i < 3; i++) {
							derivs_delta[par_index][nTile][i] = (xyz_pm[n][0][i] - xyz_pm[n][1][i]) / (2.0 * deltas[par_index]); 
						}
					}
				}
				// find derivatives by 	camera_xyz,	camera_atr); that do not require re-programming of the ERS arrays
				double [][][] atrxyz_pm = new double [6][2][];
				for (int n = 0; n < 3; n++ ){
					int par_index =DW_DAZ + n;
					double [] cam_atr = camera_atr.clone();
					cam_atr[n] = camera_atr[n] + deltas[par_index];
					atrxyz_pm[n][0] = getWorldCoordinatesERS(centerX, centerY, disparity, true, camera_xyz, cam_atr);
					cam_atr[n] = camera_atr[n] - deltas[par_index];
					atrxyz_pm[n][1] = getWorldCoordinatesERS(centerX, centerY, disparity, true, camera_xyz, cam_atr);
					par_index =DW_DX + n;
					double [] cam_xyz = camera_xyz.clone();
					cam_xyz[n] = camera_xyz[n] + deltas[par_index];
					atrxyz_pm[n+3][0] = getWorldCoordinatesERS(centerX, centerY, disparity, true, cam_xyz, camera_atr);
					cam_xyz[n] = camera_xyz[n] - deltas[par_index];
					atrxyz_pm[n+3][1] = getWorldCoordinatesERS(centerX, centerY, disparity, true, cam_xyz, camera_atr);
				}
				for (int n = 0; n < 6; n++ ){
					int par_index =DW_DAZ + n;
					if ((atrxyz_pm[n][0] != null) && (atrxyz_pm[n][1] != null)) {
						derivs_delta[par_index][nTile] = new double[3];
						for (int i = 0; i < 3; i++) {
							derivs_delta[par_index][nTile][i] = (atrxyz_pm[n][0][i] - atrxyz_pm[n][1][i]) / (2.0 * deltas[par_index]); 
						}
					}
				}					
			}
		}
		double [][][][] ers_pm = new double [6][2][tiles][];
		double [] ers_watr_center_dt_original = ers_watr_center_dt.clone();
		double [] ers_wxyz_center_dt_original = ers_wxyz_center_dt.clone();
		for (int nf = 0; nf < 12; nf++) {
			int n = nf/2;
			int sign =   ((nf & 1) != 0) ? -1 : 1;
			int sign01 = ((nf & 1) != 0) ?  1 : 0;
			boolean lin_not_ang = n >= 3;
			int component = n % 3;
			int par_index =DW_DVAZ + n;
			ers_watr_center_dt = ers_watr_center_dt_original.clone();
			ers_wxyz_center_dt = ers_wxyz_center_dt_original.clone();
			if (lin_not_ang) {
				ers_wxyz_center_dt[component] = ers_wxyz_center_dt_original[component] + sign * deltas[par_index];
			} else {
				ers_watr_center_dt[component] = ers_watr_center_dt_original[component] + sign * deltas[par_index];
			}
			setupERS();
			for (int tileY = 0; tileY < tilesY; tileY++) {
				for (int tileX = 0; tileX < tilesX; tileX++) {
					int nTile = tileX + tileY * tilesX;
					double centerX = tileX * transform_size + transform_size/2; //  - shiftX;
					double centerY = tileY * transform_size + transform_size/2; //  - shiftY;
					double disparity = dsrbg_camera[QuadCLT.DSRBG_DISPARITY][nTile];
					if (disparity < max_inf_disparity) { // or only for abs() ?
						disparity = 0.0;
					}
					ers_pm[n][sign01][nTile]=  getWorldCoordinatesERS(centerX, centerY, disparity, true, camera_xyz, camera_atr);
				}
			}
		}
		// restore original ERS data
		ers_watr_center_dt = ers_watr_center_dt_original.clone();
		ers_wxyz_center_dt = ers_wxyz_center_dt_original.clone();
		setupERS();
		
		for (int n = 0; n < 6; n++) {
			int par_index =DW_DVAZ + n;
			for (int nTile = 0; nTile < tiles; nTile++) {
				if ((ers_pm[n][0][nTile] != null) && (ers_pm[n][1][nTile] != null)) {
					derivs_delta[par_index][nTile] = new double[3];
					for (int i = 0; i <3; i++) {
						derivs_delta[par_index][nTile][i] = (ers_pm[n][0][nTile][i] - ers_pm[n][1][nTile][i]) / (2.0 * deltas[par_index]);
					}
				}
			}
		}

		
		for (int tileY = 0; tileY < tilesY; tileY++) {
			for (int tileX = 0; tileX < tilesX; tileX++) {
				int nTile = tileX + tileY * tilesX;
				double centerX = tileX * transform_size + transform_size/2; //  - shiftX;
				double centerY = tileY * transform_size + transform_size/2; //  - shiftY;
				double disparity = dsrbg_camera[QuadCLT.DSRBG_DISPARITY][nTile];
				if (disparity < max_inf_disparity) { // or only for abs() ?
					disparity = 0.0;
				}
				boolean is_infinity = (disparity == 0.0);
				if (nTile == dbg_tile) {
					System.out.println("compareDSItoWorldDerivatives(): nTile = "+nTile+" ("+tileX+"/"+tileY+")");
				}
				Vector3D[] wxyz_derivs = getDWorldDPixels(
						true,                       // boolean correctDistortions,
						is_infinity,                // boolean is_infinity,
						new double[] {centerX, centerY, disparity}, //  double [] pXpYD)
						camera_xyz,                 // double [] camera_xyz, // camera center in world coordinates
//						camera_atr,                 // double [] camera_atr);  // camera orientation relative to world frame
						rot_deriv_matrices_inverse, // Matrix [] rot_deriv_matrices_inverse, Array of 4 matrices - 
						//                             inverse transformation matrix and 3 their derivatives by az, tl, rl 
						(nTile == dbg_tile)? 2:0);
				if (wxyz_derivs != null) {
					if (wxyz_derivs[0] != null) {
						wxyz[nTile] = wxyz_derivs[0].toArray();
						double [] wxyz4 = {wxyz[nTile][0],wxyz[nTile][1],wxyz[nTile][2], is_infinity? 0.0: 1.0};
/*						
						pXpYD_test[nTile] = getImageCoordinatesERS(
								wxyz4, // double [] xyzw,
								true, // boolean correctDistortions, // correct distortion (will need corrected background too !)
								camera_xyz,     // double [] camera_xyz, // camera center in world coordinates
								camera_atr,    // double [] camera_atr);  // camera orientation relative to world frame
								LINE_ERR); // double    line_err) // threshold error in scan lines (1.0)
*/
						pXpYD_test[nTile] = getImageCoordinatesERS(
								wxyz4,          // double [] xyzw,
								true,           // boolean correctDistortions, // correct distortion (will need corrected background too !)
								camera_xyz,     // double [] camera_xyz, // camera center in world coordinates
								rot_matrix,     // Matrix    rot_matrix,  // 1-st of 4 matrices 
								LINE_ERR);           // double    line_err) // threshold error in scan lines (1.0)

					}
					for (int i = 0; i < derivs_true.length; i++) {
						if (wxyz_derivs[i+1] != null) {
							derivs_true[i][nTile] = wxyz_derivs[i+1].toArray();
						}
					}
				}
			}
		}
		int show_gap_h = 6; // 324+6 = 330 
		int show_gap_v = 8; // 242 + 8 = 250
		int tilesX1 = tilesX + show_gap_h; 
		int tilesY1 = tilesY + show_gap_v;
		int tilesX2 = 2*tilesX + show_gap_h; 
		int tilesY2 = 2*tilesY + show_gap_v;
		int tiles2 = tilesX2 * tilesY2; 
		
		if (debug_level > 0) {
			String title3 =     scene_QuadClt.getImageName()+"-pXpYD_comprison";
			String [] titles3 = {"orig","restored","diff"};
			double [][] dbg_img3 = new double [titles3.length][tiles2];
			for (int i = 0; i < dbg_img3.length; i++) {
				Arrays.fill(dbg_img3[i], Double.NaN);
			}

			for (int tileY = 0; tileY < tilesY; tileY++) {
				for (int tileX = 0; tileX < tilesX; tileX++) {
					int nTile = tileX + tileY * tilesX;
					if (pXpYD_original[nTile] != null) {
						dbg_img3[0][tileX + 0 * tilesX1 + tilesX2 *(tileY + 0 * tilesY1)] = pXpYD_original[nTile][0];
						dbg_img3[0][tileX + 1 * tilesX1 + tilesX2 *(tileY + 0 * tilesY1)] = pXpYD_original[nTile][1];
						dbg_img3[0][tileX + 0 * tilesX1 + tilesX2 *(tileY + 1 * tilesY1)] = pXpYD_original[nTile][2];
						if (pXpYD_test[nTile] != null) {
							dbg_img3[1][tileX + 0 * tilesX1 + tilesX2 *(tileY + 0 * tilesY1)] = pXpYD_test[nTile][0];
							dbg_img3[1][tileX + 1 * tilesX1 + tilesX2 *(tileY + 0 * tilesY1)] = pXpYD_test[nTile][1];
							dbg_img3[1][tileX + 0 * tilesX1 + tilesX2 *(tileY + 1 * tilesY1)] = pXpYD_test[nTile][2];

							dbg_img3[2][tileX + 0 * tilesX1 + tilesX2 *(tileY + 0 * tilesY1)] = pXpYD_test[nTile][0] - pXpYD_original[nTile][0];
							dbg_img3[2][tileX + 1 * tilesX1 + tilesX2 *(tileY + 0 * tilesY1)] = pXpYD_test[nTile][1] - pXpYD_original[nTile][1];
							dbg_img3[2][tileX + 0 * tilesX1 + tilesX2 *(tileY + 1 * tilesY1)] = pXpYD_test[nTile][2] - pXpYD_original[nTile][2];
						}
					}
				}
			}
			(new ShowDoubleFloatArrays()).showArrays(
					dbg_img3,
					tilesX2,
					tilesY2,
					true,
					title3,
					titles3);
		}
		if (debug_level > 0) {
			String title_deriv =     scene_QuadClt.getImageName()+"-derivatives_comprison";
			String [] titles_deriv = new String [3 * DW_DERIV_NAMES.length];
			for (int i = 0; i < DW_DERIV_NAMES.length; i++) {
				titles_deriv[3 * i + 0] = DW_DERIV_NAMES[i]+"-deriv";
				titles_deriv[3 * i + 1] = DW_DERIV_NAMES[i]+"-delta";
				titles_deriv[3 * i + 2] = DW_DERIV_NAMES[i]+"-diff";
			}
			double [][] dbg_img = new double [titles_deriv.length][tiles2];
			for (int i = 0; i < dbg_img.length; i++) {
				Arrays.fill(dbg_img[i], Double.NaN);
			}
			for (int tileY = 0; tileY < tilesY; tileY++) {
				for (int tileX = 0; tileX < tilesX; tileX++) {
					int nTile = tileX + tileY * tilesX;
					for (int n = 0; n < DW_DERIV_NAMES.length; n++) {
						if (derivs_true[n][nTile] != null) {
							dbg_img[3 * n + 0][tileX + 0 * tilesX1 + tilesX2 *(tileY + 0 * tilesY1)] = derivs_true[n][nTile][0];
							dbg_img[3 * n + 0][tileX + 1 * tilesX1 + tilesX2 *(tileY + 0 * tilesY1)] = derivs_true[n][nTile][1];
							dbg_img[3 * n + 0][tileX + 0 * tilesX1 + tilesX2 *(tileY + 1 * tilesY1)] = derivs_true[n][nTile][2];
						}
						if (derivs_delta[n][nTile] != null) {
							dbg_img[3 * n + 1][tileX + 0 * tilesX1 + tilesX2 *(tileY + 0 * tilesY1)] = derivs_delta[n][nTile][0];
							dbg_img[3 * n + 1][tileX + 1 * tilesX1 + tilesX2 *(tileY + 0 * tilesY1)] = derivs_delta[n][nTile][1];
							dbg_img[3 * n + 1][tileX + 0 * tilesX1 + tilesX2 *(tileY + 1 * tilesY1)] = derivs_delta[n][nTile][2];
						}
						if ((derivs_true[n][nTile] != null) && (derivs_delta[n][nTile] != null)) {
							dbg_img[3 * n + 2][tileX + 0 * tilesX1 + tilesX2 *(tileY + 0 * tilesY1)] = derivs_true[n][nTile][0] - derivs_delta[n][nTile][0];
							dbg_img[3 * n + 2][tileX + 1 * tilesX1 + tilesX2 *(tileY + 0 * tilesY1)] = derivs_true[n][nTile][1] - derivs_delta[n][nTile][1];
							dbg_img[3 * n + 2][tileX + 0 * tilesX1 + tilesX2 *(tileY + 1 * tilesY1)] = derivs_true[n][nTile][2] - derivs_delta[n][nTile][2];
						}
					}
				}
			}
			(new ShowDoubleFloatArrays()).showArrays(
					dbg_img,
					tilesX2,
					tilesY2,
					true,
					title_deriv,
					titles_deriv);
			System.out.println("compareDSItoWorldDerivatives() DONE.");
			
		}
	}
	
	public static double [][] combineXYZATR(
			double [][] reference_xyzatr,
			double [][] scene_xyzatr){
		return combineXYZATR(
				reference_xyzatr[0],
				reference_xyzatr[1],
				scene_xyzatr[0],
				scene_xyzatr[1]);
	}
		
	public static double [][] invertXYZATR(
			double [][] source_xyzatr){
		return  invertXYZATR(
				source_xyzatr[0],
				source_xyzatr[1]);
	}
	
	public static double [][] combineXYZATR(
			double [] reference_xyz,
			double [] reference_atr,
			double [] scene_xyz,
			double [] scene_atr)
	{
		Rotation   ref_rotation=   new Rotation(RotationOrder.YXZ, ROT_CONV, reference_atr[0],reference_atr[1],reference_atr[2]); // null
		Rotation   scene_rotation= new Rotation(RotationOrder.YXZ, ROT_CONV, scene_atr[0],    scene_atr[1],    scene_atr[2]);
		Vector3D   ref_offset =    new Vector3D(reference_xyz);
		Vector3D   scene_offset =  new Vector3D(scene_xyz);
		Vector3D   offset = ref_offset.add(ref_rotation.applyTo(scene_offset));
		Rotation   rotation = ref_rotation.applyTo(scene_rotation);
		double []  angles = rotation.getAngles(RotationOrder.YXZ, ROT_CONV);
		double []  xyz = offset.toArray();
		return new double[][] {xyz,angles};
		
	}

	public static double [][] invertXYZATR(
			double [] source_xyz,
			double [] source_atr)
	{
		Rotation   scene_rotation= new Rotation(RotationOrder.YXZ, ROT_CONV, source_atr[0],    source_atr[1],    source_atr[2]);
		Vector3D   scene_offset =  new Vector3D(-source_xyz[0], -source_xyz[1], -source_xyz[2]);
		Rotation   rotation =      scene_rotation.revert(); // 
		Vector3D   offset = rotation.applyTo(scene_offset); // junk
		double []  angles = rotation.getAngles(RotationOrder.YXZ, ROT_CONV);
		double []  xyz = offset.toArray();
		return new double[][] {xyz,angles};
	}
	
	
	
	public double [] getImageCoordinatesERS(
			double [] xyzw,
			boolean correctDistortions, // correct distortion (will need corrected background too !)
			double [] camera_xyz, // camera center in world coordinates
			double [] camera_atr,  // camera orientation relative to world frame
			double    line_err) // threshold error in scan lines (1.0)
	{
///		boolean is_infinity = xyzw[3] == 0;
		double rxyzw3 = 0.0; //  = 1.0/xyzw[3]
		if (Math.abs(xyzw[3]) >= THRESHOLD) {
			rxyzw3 =  1.0/xyzw[3];
		} else if (xyzw[3] > 0.0) {
			rxyzw3 =  1.0/THRESHOLD;
		} else {
			rxyzw3 =  -1.0/THRESHOLD;
		}
		Vector3D world_xyz = new Vector3D(xyzw[0],xyzw[1],xyzw[2]);
		world_xyz.scalarMultiply(rxyzw3);
		/*
		if (!is_infinity) {
			world_xyz.scalarMultiply(1.0/xyzw[3]);
		}
		*/
		// convert to camera-centered, world-parallel coordinates 
///		Vector3D cam_center_world = (is_infinity) ? world_xyz : world_xyz.subtract(new Vector3D(camera_xyz));
		Vector3D cam_center_world = world_xyz.subtract(new Vector3D(camera_xyz));
		
		// rotate to match camera coordinates when scanning the center line
		Rotation   cam_orient_center= new Rotation(RotationOrder.YXZ, ROT_CONV, camera_atr[0],camera_atr[1],camera_atr[2]);
		Vector3D cam_center_local = cam_orient_center.applyTo(cam_center_world);
		double line = pixelCorrectionHeight / 2;
		double err = pixelCorrectionHeight / 2;
		double [] dxy = null;
		// multiple iterations starting with no ERS distortions
		for (int ntry = 0; (ntry < 100) && (err > line_err); ntry++) {
			// current camera offset in the centerline camera frame
			int iline0 = (int) Math.floor(line);
			double fline = line-iline0;
			int iline1 = iline0+1;
			if (iline1 >= pixelCorrectionHeight) {
				iline1 = pixelCorrectionHeight-1;
			}
			Vector3D cam_now_local0 = new Vector3D(ers_xyz[iline0]);
			Vector3D cam_now_local1 = new Vector3D(ers_xyz[iline1]);
			Vector3D cam_now_local = cam_now_local0.scalarMultiply(1.0 - fline).add(fline,cam_now_local1); 
			
///			Vector3D cam_center_now_local = (is_infinity) ? cam_center_local :  cam_center_local.subtract(cam_now_local); // skip translation for infinity
			Vector3D cam_center_now_local = cam_center_local.subtract(cam_now_local); // skip translation for infinity
			Quaternion qpix0 = ers_quaternion[iline0];
			Quaternion qpix1 = ers_quaternion[iline1];
			Quaternion qpix= (qpix0.multiply(1.0-fline)).add(qpix1.multiply(fline));
			Rotation cam_orient_now_local = new Rotation(qpix.getQ0(), qpix.getQ1(), qpix.getQ2(),qpix.getQ3(), true); // boolean 			
			Vector3D v3 = cam_orient_now_local.applyTo(cam_center_now_local);
			
			double [] xyz = v3.toArray();
///			if (Math.abs(xyz[2]) < THRESHOLD) {
///				return null; // object too close to the lens
///			}
			double recip_disp; // = xyz[2];
			if (Math.abs(xyz[2]) >= THRESHOLD) {
				recip_disp = 1.0/xyz[2];
				// object too close to the lens (positive or negative)				
			} else if (xyz[2] > 0) {
				recip_disp = 1.0/THRESHOLD;
			} else {
				recip_disp = -1.0/THRESHOLD;
			}
			double pXc =       -(1000.0*focalLength / pixelSize) * xyz[0] * recip_disp; // / xyz[2];
			double pYc =        (1000.0*focalLength / pixelSize) * xyz[1] * recip_disp; // / xyz[2];
//			double disparity = is_infinity ? 0.0 : (-(1000.0*focalLength / pixelSize) / xyz[2] * SCENE_UNITS_SCALE * disparityRadius);
			double disparity = -(1000.0*focalLength / pixelSize) * recip_disp * SCENE_UNITS_SCALE * disparityRadius;
			
			double rND = Math.sqrt(pXc*pXc + pYc*pYc)*0.001*this.pixelSize; // mm

			double rD2RND = correctDistortions?getRDistByR(rND/this.distortionRadius):1.0;
			double px = pXc * rD2RND + 0.5 * this.pixelCorrectionWidth;  // distorted coordinates relative to the (0.5 * this.pixelCorrectionWidth, 0.5 * this.pixelCorrectionHeight)
			double py = pYc * rD2RND + 0.5 * this.pixelCorrectionHeight; // in pixels
			dxy = new double [] {px, py, disparity};
			double line1 = py;
			if (line1 < 0.0) {
				line1 = 0.0;
			} else if (line1 > (pixelCorrectionHeight-1)) {
				line1 = pixelCorrectionHeight -1;
			}
			err = Math.abs(line1 - line);
			line = line1;
		}
		return dxy;
	}

	public double [] getImageCoordinatesERS(
			double [] xyzw,
			boolean correctDistortions, // correct distortion (will need corrected background too !)
			double [] camera_xyz,       // camera center in world coordinates
			//double [] camera_atr,     // camera orientation relative to world frame
			Matrix    rot_matrix,      // 1-st of 4 matricesArray of 4 matrices - 
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
//		Rotation   cam_orient_center= new Rotation(RotationOrder.YXZ, ROT_CONV, camera_atr[0],camera_atr[1],camera_atr[2]);
//		Vector3D cam_center_local = cam_orient_center.applyTo(cam_center_world);
		Vector3D cam_center_local = matrixTimesVector(rot_matrix, cam_center_world);
		
		double line = pixelCorrectionHeight / 2;
		double err = pixelCorrectionHeight / 2;
		double [] dxy = null;
		// multiple iterations starting with no ERS distortions
		for (int ntry = 0; (ntry < 100) && (err > line_err); ntry++) {
			// current camera offset in the centerline camera frame
			int iline0 = (int) Math.floor(line);
			double fline = line-iline0;
			int iline1 = iline0+1;
			if (iline1 >= pixelCorrectionHeight) {
				iline1 = pixelCorrectionHeight-1;
			}
			Vector3D cam_now_local0 = new Vector3D(ers_xyz[iline0]);
			Vector3D cam_now_local1 = new Vector3D(ers_xyz[iline1]);
			Vector3D cam_now_local = cam_now_local0.scalarMultiply(1.0 - fline).add(fline,cam_now_local1); 
			
			Vector3D cam_center_now_local = (is_infinity) ? cam_center_local :  cam_center_local.subtract(cam_now_local); // skip translation for infinity
			Quaternion qpix0 = ers_quaternion[iline0];
			Quaternion qpix1 = ers_quaternion[iline1];
			Quaternion qpix= (qpix0.multiply(1.0-fline)).add(qpix1.multiply(fline));
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
			dxy = new double [] {px, py, disparity};
			double line1 = py;
			if (line1 < 0.0) {
				line1 = 0.0;
			} else if (line1 > (pixelCorrectionHeight-1)) {
				line1 = pixelCorrectionHeight -1;
			}
			err = Math.abs(line1 - line);
			line = line1;
		}
		return dxy;
	}
	
	
	// Derivatives section
	/**
	 * Calculate scene pixel coordinates plus disparity from the reference scene pixel coordinates and 
	 * partial derivative of the scene pixel coordinates for reference scene pixel coordinates and
	 * ERS and pose parameters of both reference and the target scene. First index 0 - new coordinates,
	 * other 27 correspond to DP_DPX..DP_DSZ 
	 * @param ers_scene ErsCorrection instance of the scene (this is ErsCorrection instance for the reference)
	 * @param correctDistortions Apply radial distortion correction (not tested when false)
	 * @param is_infinity True for the points at infinity
	 * @param pXpYD_reference pX, pY and disparity of the reference scene.
	 * @param reference_xyz X,Y,Z location of the reference scene lens center during middle line acquisition
	 *                      (typically {0,0,0}).
	 * @param scene_xyz World location of the scene camera (lens center) in the reference scene coordinates.
	 * @param reference_matrices_inverse Reference scene rotation matrix and its derivatives by Azimuth, Tilt, Roll
	 *                                   inverted, to transform local to world.
	 * @param scene_matrices_inverse Scene rotation matrix and its derivatives by Azimuth, Tilt, Roll inverted,
	 *                                   used to calculate derivatives of world XYZ by scene pX,pY Disparity
	 *                                   and ERS and pose parameters.
	 * @param scene_rot_matrix Scene rotation matrix direct, to calculate scene pX, pY, Disparity from world XYZ.
	 * @param debug_level Debug level
	 * @return Array of 3-element arrays. [0] - {pX, pY, disparity} of the scene matching pXpYD_reference. Next 27
	 *         vectors correspond to derivatives, ordered according to by DP_DPX...DP_DSZ (plus 1).
	 */
	
	
	double[][] getDPxSceneDParameters(
			ErsCorrection    ers_scene,
			boolean correctDistortions,
			boolean is_infinity,
			double [] pXpYD_reference,
			double [] reference_xyz, // reference (this) camera center in world coordinates
			double [] scene_xyz,    // (other, ers_scene) camera center in world coordinates
			//			double [] camera_atr,  // camera orientation relative to world frame
			Matrix [] reference_matrices_inverse,
			Matrix [] scene_matrices_inverse,
			Matrix    scene_rot_matrix,        // single rotation (direct) matrix for the scene
			int debug_level)
	{
		if (Double.isNaN(pXpYD_reference[2])) {
			return null;
		}
		Vector3D[]  reference_vectors = getDWorldDPixels(
				correctDistortions,
				is_infinity,
				pXpYD_reference,
				reference_xyz, // camera center in world coordinates
				reference_matrices_inverse,
				debug_level);
		double [] xyz3 =  reference_vectors[0].toArray();
		double [] xyz4 = {xyz3[0], xyz3[1], xyz3[2],is_infinity?0.0:1.0};

		double [] pXpYD_scene = ers_scene.getImageCoordinatesERS(
				xyz4,               // double [] xyzw,
				correctDistortions, // boolean correctDistortions, // correct distortion (will need corrected background too !)
				scene_xyz,          // double [] camera_xyz,       // camera center in world coordinates
				scene_rot_matrix,   // Matrix    rot_matrix,      // 1-st of 4 matricesArray of 4 matrices - 
				LINE_ERR);          // double    line_err) // threshold error in scan lines (1.0)

		Vector3D[]  scene_vectors = ers_scene.getDWorldDPixels( // [0] X,Y,Z, other ones [DW_D* + 1] 
				correctDistortions,
				is_infinity,
				pXpYD_scene,
				scene_xyz,
				scene_matrices_inverse,
				debug_level);
		Matrix dx_dpscene = new Matrix(new double[][] {
			scene_vectors[DW_DPX + 1].toArray(),
			scene_vectors[DW_DPY + 1].toArray(),
			scene_vectors[DW_DD +  1].toArray()}).transpose();
		// next variables are {x,y,z} * disparity
		double xref = xyz4[0];
		double yref = xyz4[1];
		double zref = xyz4[2];
		double z2ref = zref * zref;
		double xscene = scene_vectors[0].getX();
		double yscene = scene_vectors[0].getY();
		double zscene = scene_vectors[0].getZ();
		double z2scene = zscene * zscene;
		
		if (is_infinity) {
			// dividing X, Y by Z, getting 2x2 Matrix
			/*
			 * U=x/z, V = y/z
			 * dU/dx = 1/z, dU/dy = 0,   dU/dz = -x/z^2
			 * dV/dx = 0,   dV/dy = 1/z, dV/dz = -y/z^2
			 * dU/dpX = 1/z * (dx/dpX) - x/z^2 * (dz/dpX)
			 * dU/dpY = 1/z * (dx/dpY) - x/z^2 * (dz/dpY)
			 * dV/dpX = 1/z * (dy/dpX) - y/z^2 * (dz/dpX)
			 * dV/dpY = 1/z * (dy/dpY) - y/z^2 * (dz/dpY)  
			 */
			double dU_dpX = dx_dpscene.get(0, 0) / zscene - xscene * dx_dpscene.get(2, 0) / z2scene;
			double dU_dpY = dx_dpscene.get(0, 1) / zscene - xscene * dx_dpscene.get(2, 1) / z2scene;
			double dV_dpX = dx_dpscene.get(1, 0) / zscene - yscene * dx_dpscene.get(2, 0) / z2scene;
			double dV_dpY = dx_dpscene.get(1, 1) / zscene - yscene * dx_dpscene.get(2, 1) / z2scene;  
			dx_dpscene = new Matrix (new double[][] {
				{ dU_dpX, dU_dpY, 0.0},
				{ dV_dpX, dV_dpY, 0.0},
				{ 0.0,    0.0,    1.0}});
			xyz4[0] /= xyz4[2];
			xyz4[1] /= xyz4[2];
			xyz4[2] =  1.0;
		}
		
		Matrix dpscene_dxyz = dx_dpscene.inverse();
		Matrix dpscene_dxyz_minus = dpscene_dxyz.times(-1.0); // negated to calculate /d{pX,pY,D) for the scene parameters
		double[][] derivatives = new double[DP_NUM_PARS+1][]; // includes [0] - pXpYD vector
		// scene pX, pY, Disparity
		derivatives[0] = pXpYD_scene;
		if (Double.isNaN(pXpYD_scene[0]) || Double.isNaN(pXpYD_scene[1]) || Double.isNaN(pXpYD_scene[2])) {
			return null;
		}
		// derivatives by the reference parameters, starting with /dpX, /dpY, /dd
		for (int indx = DW_DPX; indx <= DW_DZ; indx++) {
			int vindx = indx+1;
			if (is_infinity) {
				Vector3D dw_dp = new Vector3D(
						reference_vectors[vindx].getX()/zref - xref * reference_vectors[vindx].getZ() / z2ref,
						reference_vectors[vindx].getY()/zref - yref * reference_vectors[vindx].getZ() / z2ref,
						0.0);
				derivatives[vindx] = matrixTimesVector(dpscene_dxyz, dw_dp).toArray();
			} else {
				derivatives[vindx] = matrixTimesVector(dpscene_dxyz, reference_vectors[vindx]).toArray();
			}
			if (Double.isNaN(derivatives[vindx][0]) || Double.isNaN(derivatives[vindx][1]) || Double.isNaN(derivatives[vindx][2])) {
				return null;
			}
		}
		for (int indx = DP_DSVAZ; indx < DP_NUM_PARS; indx++) { // 15,16, ...
			int indx_out = indx+1; // 16, 17,
			int indx_in =  indx_out - DP_DSVAZ + DW_DVAZ; // 4, 5, ...
			if (is_infinity) {
				Vector3D dw_dp = new Vector3D(
						scene_vectors[indx_in].getX()/zscene - xscene * scene_vectors[indx_in].getZ() / z2scene,
						scene_vectors[indx_in].getY()/zscene - yscene * scene_vectors[indx_in].getZ() / z2scene,
						0.0);
				derivatives[indx_out] = matrixTimesVector(dpscene_dxyz_minus, dw_dp).toArray();
			} else {
				derivatives[indx_out] = matrixTimesVector(dpscene_dxyz_minus, scene_vectors[indx_in]).toArray();
			}
		}
		// derivatives by the scene parameters, starting with /dvaz - angular rotaion velocity around vertical axis (ERS distortion) 
		
		return derivatives;
	}
	/**
	 * Get Derivatives of World X (right), Y (up) and Z (to the camera) by image pX, pY, disparity
	 * Both image coordinates (pX, pY, disparity) and world ones (homogeneous 4-vector) should be known
	 * e.g. calculated with getWorldCoordinatesERS() or getImageCoordinatesERS() methods
	 * @param pXpYD Pixel coordinates (pX, pY, disparity)
	 * @param camera_xyz Camera center offset in world coordinates
	 * @param rot_deriv_matrices_inverse Array of 4 matrices - inverse transformation matrix and 3 their derivatives by az, tl, rl 
	 * @return Vector3D[] for XYZ and derivatives
	 */
	Vector3D[] getDWorldDPixels(
			boolean correctDistortions,
			boolean is_infinity,
			double [] pXpYD,
			double [] camera_xyz, // camera center in world coordinates
//			double [] camera_atr,  // camera orientation relative to world frame
			Matrix [] rot_deriv_matrices_inverse,
			int debug_level)
			{
		if (Double.isNaN(pXpYD[2])) {
			return null;
		}
		int line = (int) Math.round(pXpYD[1]);
		if (line < 0) {
			line = 0;
		} else if (line >= pixelCorrectionHeight) {
			line = pixelCorrectionHeight - 1;
		}
		double pXcd = pXpYD[0] - 0.5 * this.pixelCorrectionWidth;
		double pYcd = pXpYD[1] - 0.5 * this.pixelCorrectionHeight;
		double r_scale = 0.001 * this.pixelSize / this.distortionRadius;
		double rDn =    Math.sqrt(pXcd*pXcd + pYcd*pYcd) * r_scale; // relative to distortion radius
		double rND2RD = correctDistortions?(getRByRDist(rDn, false)): 1.0; // NaN
		// rD - in mm
//		double rD = rDn * this.distortionRadius; // Math.sqrt(pXcd*pXcd + pYcd*pYcd)*0.001*this.pixelSize; // distorted radius in a virtual center camera
		// pXc, pYc - in pixels
		double pXc = pXcd * rND2RD; // non-distorted coordinates relative to the (0.5 * this.pixelCorrectionWidth, 0.5 * this.pixelCorrectionHeight)
		double pYc = pYcd * rND2RD; // in pixels
		double rn =  rDn *  rND2RD; // normalized to distortion radius
		// Calculating derivative drD2rND_drn
		double []  rad_coeff={this.distortionC,this.distortionB,this.distortionA,this.distortionA5,this.distortionA6,this.distortionA7,this.distortionA8};
		double drD2rND_drn = 0.0; // derivative of (rD/rND) over d_rn; where rn = rND * ri_scale; // relative to distortion radius
	    double rr = 1.0;
	    for (int j = 0; j < rad_coeff.length; j++){
	        drD2rND_drn += rad_coeff[j] * (j+1) * rr;
	        rr *= rn;
	    }
	    // inverting, for rNDn = rND2RD * rDn calculating d_rND2RD_d_rDn
	    double rD2rND = 1.0/rND2RD;
	    double rNDn = rDn * rND2RD; 
	    double d_rND2rD_d_rDn = -drD2rND_drn / ((rD2rND * rD2rND) * (drD2rND_drn * rNDn + rD2rND));    
	    
	    double d_rDn_d_pX = pXcd * (r_scale * r_scale) / rDn;
		double d_rDn_d_pY = pYcd * (r_scale * r_scale) / rDn;
		
 		double d_pXc_d_pX = rND2RD + pXcd * d_rND2rD_d_rDn * d_rDn_d_pX; 	    
 		double d_pXc_d_pY =          pXcd * d_rND2rD_d_rDn * d_rDn_d_pY; 	    
 		double d_pYc_d_pX =          pYcd * d_rND2rD_d_rDn * d_rDn_d_pX; 	    
 		double d_pYc_d_pY = rND2RD + pYcd * d_rND2rD_d_rDn * d_rDn_d_pY; 

		double zh = -SCENE_UNITS_SCALE * this.focalLength * this.disparityRadius / (0.001*this.pixelSize); // "+" - near, "-" far
		double xh =  SCENE_UNITS_SCALE * pXc * this.disparityRadius;
		double yh = -SCENE_UNITS_SCALE * pYc * this.disparityRadius;
		
		double d_xh_d_px =  SCENE_UNITS_SCALE * d_pXc_d_pX * this.disparityRadius;
		double d_xh_d_py =  SCENE_UNITS_SCALE * d_pXc_d_pY * this.disparityRadius;
		double d_yh_d_px = -SCENE_UNITS_SCALE * d_pYc_d_pX * this.disparityRadius;
		double d_yh_d_py = -SCENE_UNITS_SCALE * d_pYc_d_pY * this.disparityRadius;
		double [] xyz = {xh, yh, zh};
		
		Vector3D [] dWdP= new Vector3D[3];
		if (is_infinity) {
			dWdP[0] = new Vector3D(new double [] {d_xh_d_px, d_yh_d_px, 0.0});
			dWdP[1] = new Vector3D(new double [] {d_xh_d_py, d_yh_d_py, 0.0});
			dWdP[2] = new Vector3D(new double [] {0.0,       0.0,       0.0});
		} else { // not infinity
			for (int i = 0; i < 3; i++) {
				xyz[i] /= pXpYD[2];
			}
			dWdP[0] = new Vector3D(new double [] {d_xh_d_px / pXpYD[2], d_yh_d_px / pXpYD[2], 0.0});
			dWdP[1] = new Vector3D(new double [] {d_xh_d_py / pXpYD[2], d_yh_d_py / pXpYD[2], 0.0});
			dWdP[2] = new Vector3D(new double [] {-xyz[0] /   pXpYD[2],-xyz[1] / pXpYD[2],   -xyz[2] / pXpYD[2]});
		}
		
		
		// 1) Apply rotations to 3 derivative vectors d/dpx, d/dpy, d/ddisparity
		Vector3D cam_now_local = new Vector3D(ers_xyz[line]); 
		// camera orientation during pixel acquisition :
		Quaternion qpix = ers_quaternion[line];
		Rotation cam_orient_now_local = new Rotation(qpix.getQ0(), qpix.getQ1(), qpix.getQ2(),qpix.getQ3(), true); // boolean needsNormalization)
		Vector3D dw_dpX = cam_orient_now_local.applyInverseTo(dWdP[0]);
		Vector3D dw_dpY = cam_orient_now_local.applyInverseTo(dWdP[1]);
		Vector3D dw_dd =  cam_orient_now_local.applyInverseTo(dWdP[2]);

		Vector3D v3= new Vector3D(xyz);
		Vector3D cam_center_now_local = cam_orient_now_local.applyInverseTo(v3); // undone rotation relative to center

		// 2) add rotation and displacement by dt proportional to dpY to the second vector

		double d_t_d_pY = line_time; // difference in pY modifies world coordinates using linear and angular velocities
		Quaternion qpix_dt = ers_quaternion_dt[line]; // .multiply(d_t_d_pY); // rotation // ers_quaternion_dt is actually omega
		Vector3D local_omega = new Vector3D(qpix_dt.getVectorPart());
		Vector3D dw_dt = Vector3D.crossProduct(local_omega,cam_center_now_local);
		// maybe temporarily disable this term
		
		/** */
		dw_dpY = dw_dpY.add(d_t_d_pY, dw_dt); // add scaled
		if (!is_infinity) {
			dw_dpY = dw_dpY.add(d_t_d_pY, new Vector3D(ers_xyz_dt[line]));
		}
		/** */
		Vector3D cam_center_local = (is_infinity) ? cam_center_now_local :  cam_center_now_local.add(cam_now_local); // skip translation for infinity

		// undo camera rotation during acquisition of the center line. 
////	Rotation cam_orient_center= new Rotation(RotationOrder.YXZ, ROT_CONV, camera_atr[0],camera_atr[1],camera_atr[2]);
		Vector3D cam_center_world = matrixTimesVector(rot_deriv_matrices_inverse[0], cam_center_local);
// Apply same rotation to the d_d* vectors
		dw_dpX = matrixTimesVector(rot_deriv_matrices_inverse[0], dw_dpX);
		dw_dpY = matrixTimesVector(rot_deriv_matrices_inverse[0], dw_dpY);
		dw_dd =  matrixTimesVector(rot_deriv_matrices_inverse[0], dw_dd);
/*		
		Vector3D dw_dpX0 = cam_orient_center.applyInverseTo(dw_dpX);
		Vector3D dw_dpY0 = cam_orient_center.applyInverseTo(dw_dpY);
		Vector3D dw_dd0 =  cam_orient_center.applyInverseTo(dw_dd);
*/		
		// convert coordinates to the real world coordinates
		Vector3D world_xyz = (is_infinity) ? cam_center_world : cam_center_world.add(new Vector3D(camera_xyz));
		double [] wxyz = world_xyz.toArray();
		double [] wxyz4 = {wxyz[0],wxyz[1],wxyz[2], 1.0};
		if (Double.isNaN(wxyz4[0])) {
			wxyz4[0] = Double.NaN;
		}
		if (is_infinity) {
			wxyz4[3] = 0.0;
		}

// calculate other derivatives: dw_d_xyz, dw_dvxyz, dw_datr, dw_dvatr
// dw_dvatr
		double delta_t =     line_time*(line - pixelCorrectionHeight/2);
		double delta_t_xyz = (is_infinity)? 0.0: delta_t;

		Vector3D d_vaz = new Vector3D(0.0,        delta_t,     0.0);       // is that correct?
		Vector3D d_vtl = new Vector3D(delta_t,    0.0,         0.0);       // is that correct?
		Vector3D d_vrl = new Vector3D(0.0,        0.0,         delta_t);
		
		Vector3D d_vx = new Vector3D(delta_t_xyz, 0.0,         0.0);
		Vector3D d_vy = new Vector3D(0.0,         delta_t_xyz, 0.0);
		Vector3D d_vz = new Vector3D(0.0,         0.0,         delta_t_xyz);
		
		/*
		Vector3D dw_dvaz0 = cam_orient_center.applyInverseTo(Vector3D.crossProduct(d_vaz, cam_center_now_local));
		Vector3D dw_dvtl0 = cam_orient_center.applyInverseTo(Vector3D.crossProduct(d_vtl, cam_center_now_local));
		Vector3D dw_dvrl0 = cam_orient_center.applyInverseTo(Vector3D.crossProduct(d_vrl, cam_center_now_local));
		
		Vector3D dw_dvx0 =  cam_orient_center.applyInverseTo(d_vx);
		Vector3D dw_dvy0 =  cam_orient_center.applyInverseTo(d_vy);
		Vector3D dw_dvz0 =  cam_orient_center.applyInverseTo(d_vz);
		*/
		Vector3D dw_dvaz = matrixTimesVector(rot_deriv_matrices_inverse[0], Vector3D.crossProduct(d_vaz, cam_center_now_local));
		Vector3D dw_dvtl = matrixTimesVector(rot_deriv_matrices_inverse[0], Vector3D.crossProduct(d_vtl, cam_center_now_local));
		Vector3D dw_dvrl = matrixTimesVector(rot_deriv_matrices_inverse[0], Vector3D.crossProduct(d_vrl, cam_center_now_local));
		
		Vector3D dw_dvx =  matrixTimesVector(rot_deriv_matrices_inverse[0], d_vx);
		Vector3D dw_dvy =  matrixTimesVector(rot_deriv_matrices_inverse[0], d_vy);
		Vector3D dw_dvz =  matrixTimesVector(rot_deriv_matrices_inverse[0], d_vz);
		
		// overall rotations and offsets		
		
		Vector3D dw_daz = matrixTimesVector(rot_deriv_matrices_inverse[1], cam_center_local);
		Vector3D dw_dtl = matrixTimesVector(rot_deriv_matrices_inverse[2], cam_center_local);
		Vector3D dw_drl = matrixTimesVector(rot_deriv_matrices_inverse[3], cam_center_local);
		
		Vector3D dw_dx = is_infinity ? Vector3D.ZERO : Vector3D.PLUS_I;
		Vector3D dw_dy = is_infinity ? Vector3D.ZERO : Vector3D.PLUS_J;
		Vector3D dw_dz = is_infinity ? Vector3D.ZERO : Vector3D.PLUS_K;
		
		Vector3D [] result = {
				world_xyz, dw_dpX,  dw_dpY, dw_dd,
				dw_dvaz,   dw_dvtl, dw_dvrl,
				dw_dvx,    dw_dvy,  dw_dvz,
				dw_daz,    dw_dtl,  dw_drl,
				dw_dx,     dw_dy,   dw_dz};
		return result;
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





