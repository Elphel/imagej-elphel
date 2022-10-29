/**
 **
 ** XyzAtr - store position+orientation (world XYZ, Azimuth, Tilt, Roll)
 ** of other scenes relative to the position of this camera.
 ** Positions/orientations are sampled during scanning of the center line 
 **
 **
 ** Copyright (C) 2022 Elphel, Inc.
 **
 ** -----------------------------------------------------------------------------**
 **
 **  XyzAtr.java is free software: you can redistribute it and/or modify
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

public class XyzAtr {
	double [] xyz;
	double [] atr;
	
	double [] ers_xyz_dt;  // world camera Vx, Vy, Vz (m/s)
	double [] ers_atr_dt;  // camera rotations (az, tilt, roll in radians/s, corresponding to the frame center)
	double [] ers_atr_d2t; // camera rotations (az, tilt, roll in radians/s, corresponding to the frame center)
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
		double [] d = ErsCorrection.parseDoublesCSV(s);
		xyz = new double [] {d[0], d[1], d[2]};
		atr = new double [] {d[3], d[4], d[5]};
	}

	public XyzAtr(String [] ss) {
		if (ss[0] != null) {
			double [] d = ErsCorrection.parseDoublesCSV(ss[0]);
			xyz = new double [] {d[0], d[1], d[2]};
			atr = new double [] {d[3], d[4], d[5]};
			this.ers_xyz_dt =  new double[3];
			this.ers_atr_dt =  new double[3];
			this.ers_xyz_d2t = new double[3];
			this.ers_atr_d2t = new double[3];
			if (ss.length > 1) {
				if (ss[1] != null) {
					d = ErsCorrection.parseDoublesCSV(ss[1]);
					ers_xyz_dt = new double [] {d[0], d[1], d[2]};
					ers_atr_dt = new double [] {d[3], d[4], d[5]};
				}
				if (ss.length > 2) {
					if (ss[2] != null) {
						d = ErsCorrection.parseDoublesCSV(ss[2]);
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