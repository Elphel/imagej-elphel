package com.elphel.imagej.tileprocessor;
/**
 **
 ** CorrVector - Geometry correction vector
 **
 ** Copyright (C) 2017 Elphel, Inc.
 **
 ** -----------------------------------------------------------------------------**
 **
 **  CorrVector.java is free software: you can redistribute it and/or modify
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
 */

import java.util.ArrayList;

import com.elphel.imagej.common.GenericJTabbedDialog;
//import com.elphel.imagej.tileprocessor.GeometryCorrection.CorrVector;

import Jama.Matrix;

public class CorrVector{ // TODO: Update to non-quad (extract to a file first)?
	
	private static final double ROT_AZ_SGN = -1.0; // sign of first sin for azimuth rotation
	private static final double ROT_TL_SGN =  1.0; // sign of first sin for tilt rotation
	private static final double ROT_RL_SGN =  1.0; // sign of first sin for roll rotation
	private double [] vector;
	private GeometryCorrection geometryCorrection;
	public static final String CORR_TILT =    "tilt";
	public static final String CORR_AZIMUTH = "azimuth";
	public static final String CORR_ROLL =    "roll";
	public static final String CORR_ZOOM =    "zoom";
	public static final String [] CORR_IMU = {"omega_tilt", "omega_azimuth", "omega_roll",
			                                  "velocity_x", "velocity_y", "velocity_z"};
	
	

	// replacing static constants (they will be different for diffrent numSensors)
	// 1 step - add methods, 2-nd update methods themselves
	public int getLength()       {return getLength      (getNumSensors());}
	public int getLengthAngles() {return getLengthAngles(getNumSensors());}
	public int getTiltIndex()    {return getTiltIndex   (getNumSensors());}
	public int getAzimuthIndex() {return getAzimuthIndex(getNumSensors());}
	public int getRollIndex()    {return getRollIndex   (getNumSensors());}
	public int getZoomIndex()    {return getZoomIndex   (getNumSensors());}
	public int getIMUIndex()     {return getIMUIndex    (getNumSensors());}
	public double getRotAzSgn()  {return ROT_AZ_SGN;}
	public double getRotTlSgn()  {return ROT_TL_SGN;}
	public double getRotRlzSgn() {return ROT_RL_SGN;}
	public double[] getVector()  {return vector;}
	
	// static that (will) use number of subcameras
	public static int getCamerasFromEV(int ev_length) {return (ev_length + 3 - CORR_IMU.length)/4;}
	public static int getLength       (int num_chn) {return 4 * num_chn - 3 + CORR_IMU.length;}
	public static int getLengthAngles (int num_chn) {return 3 * num_chn - 2;}
	public static int getTiltIndex    (int num_chn) {return 0;}
	public static int getAzimuthIndex (int num_chn) {return     num_chn - 1;}
	public static int getRollIndex    (int num_chn) {return 2 * num_chn - 2;}
	public static int getZoomIndex    (int num_chn) {return 3 * num_chn - 2;}
	public static int getIMUIndex     (int num_chn) {return 4 * num_chn - 3;}
//	public static double getRotAzSgn(int num_chn)  {return ROT_AZ_SGN;}
//	public static double getRotTlSgn(int num_chn)  {return ROT_TL_SGN;}
//	public static double getRotRlzSgn(int num_chn) {return ROT_RL_SGN;}
	
	public SymmVectorsSet symmVectorsSet;
	
	public int getNumSensors() {
		return geometryCorrection.getNumSensors();
	}
	
	public static int getNumSensors(int vector_length) {
		int num_sensors=4;
		for (; getLength(num_sensors) < vector_length; num_sensors++);
		if (getLength(num_sensors) != vector_length) {
			throw new IllegalArgumentException ("Invalid vector length = "+vector_length);
		}
		return num_sensors;
	}
	
//getLength()	
	public static String [] getCorrNames(int num_chn) {
		String [] corr_names = new String[getLength(num_chn)];
		for (int n = 0; n < num_chn; n++) {
			if (n < (num_chn - 1)) {
				corr_names[getTiltIndex(num_chn) + n]=    CORR_TILT    + n;
				corr_names[getAzimuthIndex(num_chn) + n]= CORR_AZIMUTH + n;
				corr_names[getZoomIndex(num_chn) + n]=    CORR_ZOOM    + n;
			}
			corr_names[getRollIndex(num_chn) + n]=        CORR_ROLL    + n;
		}
		for (int i = 0; i < CORR_IMU.length; i++) {
			corr_names[getIMUIndex(num_chn) + i]=         CORR_IMU[i];
		}
		return corr_names;
	}

	
	
	public String [] getCorrNames() {
		return getCorrNames(getNumSensors());
	}
	
	
	
	public Matrix [] getRotMatrices(Matrix rigMatrix)// USED in lwir
	{
		Matrix [] rots = getRotMatrices();
		if (rigMatrix != null) {
			for (int chn = 0; chn < rots.length; chn++) {
//				rots[chn] = rigMatrix.times(rots[chn]);
				rots[chn] = rots[chn].times(rigMatrix); // rots[chn] left matrix is rotation (closest to the camera), it should remain leftmost
			}
		}
		return rots;
	}

	// not yet used
	public Matrix [][] getRotDeriveMatrices(Matrix rigMatrix)// not used in lwir
	{
		Matrix [][] derivs = getRotDeriveMatrices();
		if (rigMatrix != null) {
			for (int chn = 0; chn < derivs.length; chn++) {
				for (int deriv_index=0; deriv_index < derivs[chn].length; deriv_index++) {
//					derivs[chn][deriv_index] = rigMatrix.times(derivs[chn][deriv_index]);
					derivs[chn][deriv_index] = derivs[chn][deriv_index].times(rigMatrix); // left matrix is rotation (closest to the camera), it should remain leftmost
				}
			}
		}
		return derivs;
	}
	public Matrix [] getRotMatricesDbg() {
		Matrix [] rots = new Matrix [getNumSensors()];
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
					{ ca,               0.0,  sa * ROT_AZ_SGN }, // -1.0
					{ 0.0,              1.0,  0.0},
					{ -sa* ROT_AZ_SGN,  0.0,  ca}};

			double [][] a_t =  { // inverted - OK
					{ 1.0,  0.0, 0.0},
					{ 0.0,  ct,               st * ROT_TL_SGN}, // +1.0
					{ 0.0, -st * ROT_TL_SGN,  ct}};

			double [][] a_r =  { // inverted OK
					{ cr,                sr * ROT_RL_SGN,  0.0}, // +1.0
					{ -sr * ROT_RL_SGN,  cr,               0.0},
					{ 0.0,               0.0,              1.0}};
			rots[chn] = (new Matrix(a_r).times(new Matrix(a_t).times(new Matrix(a_az))));
			System.out.println("Matrices for camera "+chn);
			System.out.println("Azimuth "+chn+" angle = "+ azimuths[chn]);
			(new Matrix(a_az)).print(10, 7);
			System.out.println("Tilt "+chn+" angle = "+ tilts[chn]);
			(new Matrix(a_t)).print(10, 7);
			System.out.println("Roll/Zoom "+chn+" angle = "+ rolls[chn]+" zoom="+zooms[chn]+
					" cr = "+ Math.cos(rolls[chn])+ " sr = "+ Math.sin(rolls[chn]));
			(new Matrix(a_r)).print(10, 7);
			System.out.println("Combined "+chn);
			rots[chn].print(10, 7);
		}
		return rots;
	}

	public Matrix [][] getRotDeriveMatricesDbg() // USED in lwir
	{
		Matrix [][] rot_derivs = new Matrix [getNumSensors()][4]; // channel, azimuth-tilt-roll-zoom
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

			System.out.println("\nMatrices for camera "+chn);
			System.out.println("Azimuth "+chn+" angle = "+ azimuths[chn]);
			(new Matrix(a_az)).print(10, 7);
			System.out.println("Tilt "+chn+" angle = "+ tilts[chn]);
			(new Matrix(a_t)).print(10, 7);

			System.out.println("Roll/Zoom "+chn+" angle = "+ rolls[chn]+" zoom="+zooms[chn]+
					" cr = "+ Math.cos(rolls[chn])+ " sr = "+ Math.sin(rolls[chn]));
			(new Matrix(a_r)).print(10, 7);

			System.out.println("Combined "+chn);
			(new Matrix(a_r).times(new Matrix(a_t).times(new Matrix(a_az)))).print(10, 7);

			System.out.println("a_t * a_az "+chn);
			(new Matrix(a_t).times(new Matrix(a_az))).print(10, 7);



			System.out.println("\n=== Matrices of Derivatives for camera "+chn);
			System.out.println("d/d_Azimuth "+chn+" angle = "+ azimuths[chn]+" ca="+ca+" sa="+sa);
			(new Matrix(a_daz)).print(10, 7);
			System.out.println("d/d_Tilt "+chn+" angle = "+ tilts[chn]+" ct="+ct+" st="+st);
			(new Matrix(a_dt)).print(10, 7);
			System.out.println("d/d_Roll "+chn+" angle = "+ rolls[chn]+" cr="+cr+" sr="+sr);
			(new Matrix(a_dr)).print(10, 7);
			System.out.println("d/d_Zoom "+chn+" zoom="+zooms[chn]);
			(new Matrix(a_dzoom)).print(10, 7);
			System.out.println("Combined d/daz"+chn);
			rot_derivs[chn][0].print(10, 7);
			System.out.println("Combined d/dt"+chn);
			rot_derivs[chn][1].print(10, 7);
			System.out.println("Combined d/dr"+chn);
			rot_derivs[chn][2].print(10, 7);
			System.out.println("Combined d/dzoom"+chn);
			rot_derivs[chn][3].print(10, 7);
		}
		return rot_derivs;
	}





	public Matrix [] getRotMatrices() // USED in lwir TODO: Update to non-quad!
	{
		Matrix [] rots = new Matrix [getNumSensors()];
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
	public Matrix [][] getRotDeriveMatrices() // USED in lwir
	{
		Matrix [][] rot_derivs = new Matrix [getNumSensors()][4]; // channel, azimuth-tilt-roll-zoom
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

	public CorrVector (GeometryCorrection geometryCorrection)// USED in lwir
	{
		int num_sensors = geometryCorrection.getNumSensors();
		this.symmVectorsSet = 	SymmVector.getSymmVectorsSet (num_sensors);
		this.vector = new double[getLength(num_sensors)];
		this.geometryCorrection = geometryCorrection;
	}

	public CorrVector (// USED in lwir
			GeometryCorrection geometryCorrection,
			double [] sym_vector,
			boolean [] par_mask)
	{
		int num_sensors = geometryCorrection.getNumSensors();
		this.symmVectorsSet = 	SymmVector.getSymmVectorsSet (num_sensors);
		this.vector = toTarArray(sym_vector, par_mask);
		this.geometryCorrection = geometryCorrection;
	}

	public CorrVector (
			GeometryCorrection geometryCorrection,
			double [] vector)
	{
		int num_sensors = geometryCorrection.getNumSensors();
		this.symmVectorsSet = 	SymmVector.getSymmVectorsSet (num_sensors);
		if (vector != null) {
			if (vector.length != getLength(num_sensors)) {
				throw new IllegalArgumentException("vector.length = "+vector.length+" != "+getLength(num_sensors));// not used in lwir
			}
			this.vector = vector;
		} else {
			this.vector = new double[getLength(num_sensors)];// not used in lwir
		}
		this.geometryCorrection = geometryCorrection;
	}
	
	public CorrVector (
			GeometryCorrection geometryCorrection,
			GeometryCorrection sourceGeometryCorrection)
	{
		int num_sensors = geometryCorrection.getNumSensors();
		this.symmVectorsSet = 	SymmVector.getSymmVectorsSet (num_sensors);
		this.geometryCorrection = geometryCorrection;
		this.vector = new double[getLength(num_sensors)];
		CorrVector scv = sourceGeometryCorrection.getCorrVector();
		int num_sens =  geometryCorrection.getNumSensors();
		int snum_sens = sourceGeometryCorrection.getNumSensors();
		int min_num_sens = (snum_sens < num_sens) ? snum_sens : num_sens; 
		System.arraycopy(scv.vector, scv.getTiltIndex(),    vector, getTiltIndex(),    min_num_sens - 1);
		System.arraycopy(scv.vector, scv.getAzimuthIndex(), vector, getAzimuthIndex(), min_num_sens - 1);
		System.arraycopy(scv.vector, scv.getRollIndex(),    vector, getRollIndex(),    min_num_sens);
		System.arraycopy(scv.vector, scv.getZoomIndex(),    vector, getZoomIndex(),    min_num_sens - 1);
		System.arraycopy(scv.vector, scv.getIMUIndex(),     vector, getIMUIndex(),     CORR_IMU.length);
	}
	
	
	public CorrVector getCorrVector(
			GeometryCorrection geometryCorrection,
			double [] vector){// not used in lwir
		return new CorrVector(geometryCorrection,vector);
	}

	public float [] toFloatArray() {
		if (vector == null) {
			return null;
		}
		float [] fvector = new float [vector.length];
		for (int i = 0; i < vector.length; i++) {
			fvector[i] = (float) vector[i];
		}
		return fvector;
	}

	public double [] toFullRollArray()
	{
		double [] v =       vector.clone();
		double [] rolls =   getFullRolls();
		for (int i = 0; i < (getZoomIndex() - getRollIndex()); i++) {
			v[getRollIndex() + i] = rolls[i];
		}
		return v;
	}

	public double [] toArray() // USED in lwir
	{
		return vector;
	}

	public double [] getTilts() // USED in lwir
	{
//		double [] tilts =     {vector[0], vector[1], vector[2], - (vector[0] + vector[1] +vector[2])};
		double [] tilts = new double [getNumSensors()];
		int indx = getTiltIndex();
		for (int i = 0; i < (tilts.length - 1); i++) {
			tilts[i] =              vector[indx + i];
			tilts[tilts.length - 1] -= tilts[i];
		}
		return tilts;
	}

	public double [] getAzimuths() // USED in lwir
	{
//		double [] azimuths =  {vector[3], vector[4], vector[5], -(vector[3] + vector[4] + vector[5])};
		double [] azimuths = new double [getNumSensors()];
		int indx = getAzimuthIndex();
		for (int i = 0; i < (azimuths.length - 1); i++) {
			azimuths[i] =              vector[indx + i];
			azimuths[azimuths.length - 1] -= azimuths[i];
		}
		return azimuths;
	}
	
	public double [] getRolls() // not used in lwir
	{
//		double [] rolls =     {vector[6],vector[7],vector[8], vector[9]};
		double [] rolls =     new double [getNumSensors()];
		int indx = getRollIndex();
		for (int i = 0; i < rolls.length; i++) {
			rolls[i] =        vector[indx + i];
		}
		return rolls;
	}


	public double [] getZooms() // USED in lwir
	{
//		double [] zooms =     {vector[10], vector[11], vector[12], - (vector[10] + vector[11] +vector[12])};
		double [] zooms = new double [getNumSensors()];
		int indx = getZoomIndex();
		for (int i = 0; i < (zooms.length - 1); i++) {
			zooms[i] =              vector[indx + i];
			zooms[zooms.length - 1] -= zooms[i];
		}
		return zooms;
	}
	
	public double [] getIMU() {
//		double [] imu = {vector[IMU_INDEX + 0], vector[IMU_INDEX + 1], vector[IMU_INDEX + 2],vector[IMU_INDEX + 3], vector[IMU_INDEX + 4], vector[IMU_INDEX + 5]};
		int indx = getIMUIndex(); // Was bug: getZoomIndex();
		double [] imu = new double[6];
		for (int i = 0; i < imu.length; i++) {
			imu[i] = vector[indx + i];
		}
		return imu;
	}

	public void setIMU(double [] imu) {
		int indx = getIMUIndex(); // Was bug: getZoomIndex();
		for (int i = 0; i < imu.length; i++) {
			vector [indx + i] = imu[i];
		}
		System.out.print("setIMU(");
		for (int i = 0; i < imu.length; i++) {
			System.out.print(imu[i]+",");
		}
		System.out.println(")");
	}

	public double setZoomsFromF(double [] f) { // USED in lwir
		int indx = getZoomIndex();
//		double f_avg = (f0+f1+f2+f3)/4;
		double f_avg = 0;
		for (int i = 0; i < f.length;i++) {
			f_avg += f[i];
		}
		f_avg /= f.length;
		for (int i = 0; i < (f.length-1); i++) {
			vector[indx+i] = (f[i] - f_avg)/f_avg; // here first time wrong number of sensors
		}
		return f_avg;
	}
	
	// Tilts in radians, theta in degrees
	public double setTiltsFromThetas(double [] t) { // USED in lwir
		int indx = getTiltIndex();
		double t_avg = 0;
		for (int i = 0; i < t.length; i++) {
			t_avg += t[i];
		}
		t_avg /= t.length;
		for (int i = 0; i < (t.length -1); i++) {
			vector[indx+i] = (t[i] - t_avg)*Math.PI/180.0;
		}
		return t_avg;
	}
	
	// Azimuths in radians, headings in degrees
	public double setAzimuthsFromHeadings(double [] h) { // USED in lwir
		int indx = getAzimuthIndex();
		double h_avg = 0;
		for (int i = 0; i < h.length; i++) {
			h_avg += h[i];
		}
		h_avg /= h.length;
		for (int i = 0; i < (h.length -1); i++) {
			vector[indx+i] = (h[i] - h_avg)*Math.PI/180.0;
		}
		return h_avg;
	}

	// Include factory calibration rolls
	public double [] getFullRolls() // USED in lwir
	{
		int indx =getRollIndex();
		double d2r= Math.PI/180.0;
		double [] rolls =  new double [geometryCorrection.roll.length];
		for (int i = 0; i < rolls.length; i++) {
			rolls[i] = vector[indx + i] + d2r * geometryCorrection.roll[i];
		}
		return rolls;
	}
	
	/**
	 * Return parameter value for reports
	 * @param indx parameter index (use CorrVector.XXX static integers)
	 * @param inPix show result in pixels , false - in radians (even for zooms)
	 * @return parameter value
	 */
	public double getExtrinsicParameterValue(int indx, boolean inPix) { // not used in lwir
		if (indx <0) return Double.NaN;
		if (indx <    getRollIndex()) return vector[indx]* (inPix? (1000.0*geometryCorrection.focalLength     /geometryCorrection.pixelSize): 1.0); // tilt and azimuth
		if (indx < getLengthAngles()) return vector[indx]* (inPix? (1000.0*geometryCorrection.distortionRadius/geometryCorrection.pixelSize): 1.0); //  rolls
		if (indx <     getIMUIndex()) return vector[indx]* (inPix? (1000.0*geometryCorrection.distortionRadius/geometryCorrection.pixelSize): 1.0); //  zooms
		if (indx < (getIMUIndex()+2)) return vector[indx]* (inPix? (1000.0*geometryCorrection.focalLength     /geometryCorrection.pixelSize): 1.0); // omega_tilt and omega_azimuth
		if (indx < (getIMUIndex()+3)) return vector[indx]* (inPix? (1000.0*geometryCorrection.distortionRadius/geometryCorrection.pixelSize): 1.0); // omega_roll
		if (indx <       getLength()) return vector[indx]* (inPix? (1000.0): 1.0); //  m/s or mm/s
		return Double.NaN;
	}
	public double getExtrinsicSymParameterValue(int indx, boolean inPix) { // not used in lwir
		double [] sym_vect = toSymArray(null);
		if (indx <0) return Double.NaN;
		if (indx <    getRollIndex()) return sym_vect[indx]* (inPix? (1000.0*geometryCorrection.focalLength     /geometryCorrection.pixelSize): 1.0); // tilt and azimuth
		if (indx < getLengthAngles()) return sym_vect[indx]* (inPix? (1000.0*geometryCorrection.distortionRadius/geometryCorrection.pixelSize): 1.0); //  rolls
		if (indx <     getIMUIndex()) return sym_vect[indx]* (inPix? (1000.0*geometryCorrection.distortionRadius/geometryCorrection.pixelSize): 1.0); //  zooms
		if (indx < (getIMUIndex()+2)) return sym_vect[indx]* (inPix? (1000.0*geometryCorrection.focalLength     /geometryCorrection.pixelSize): 1.0); // omega_tilt and omega_azimuth
		if (indx < (getIMUIndex()+3)) return sym_vect[indx]* (inPix? (1000.0*geometryCorrection.distortionRadius/geometryCorrection.pixelSize): 1.0); // omega_roll
		if (indx <        getLength()) return sym_vect[indx]* (inPix? (1000.0): 1.0);  // m/s mm/s
		return Double.NaN;
	}
	
	@Override
	public String toString() { // USED in lwir
		return toString(false);
	}
	
	public double getTiltAzPerPixel() {
		return 1.0/1000.0*geometryCorrection.focalLength/geometryCorrection.pixelSize;
	}
	
	public String toString(boolean short_out) {
		String s;
		double [] sym_vect = toSymArray(null);
		double [] v = new double [vector.length];
		double [] sv = new double [vector.length];
		for (int i = 0; i < getRollIndex(); i++){
			v[i] = vector[i]*1000.0*geometryCorrection.focalLength/geometryCorrection.pixelSize; // tilt and azimuth
			sv[i] = sym_vect[i]*1000.0*geometryCorrection.focalLength/geometryCorrection.pixelSize; // all 3 angles
		}
		for (int i = getRollIndex(); i < getLengthAngles(); i++){
			v[i] = vector[i]*1000.0*geometryCorrection.distortionRadius/geometryCorrection.pixelSize; // rolls
			sv[i] = sym_vect[i]*1000.0*geometryCorrection.distortionRadius/geometryCorrection.pixelSize; // combined rolls
		}
		for (int i = getLengthAngles(); i < getIMUIndex(); i++){
			v[i] = vector[i]*1000.0*geometryCorrection.distortionRadius/geometryCorrection.pixelSize; // zooms
			sv[i] = sym_vect[i]*1000.0*geometryCorrection.distortionRadius/geometryCorrection.pixelSize; // zooms
		}

		for (int i = getIMUIndex(); i < getIMUIndex()+2; i++){
			v[i] =  vector[i]*  1000.0*geometryCorrection.focalLength/geometryCorrection.pixelSize;    // omega_tilt and omega_azimuth
			sv[i] = sym_vect[i]*1000.0*geometryCorrection.focalLength/geometryCorrection.pixelSize;    // omega_tilt and omega_azimuth
		}
		for (int i = getIMUIndex()+2; i < getIMUIndex()+3; i++){
			v[i] = vector[i]*1000.0*geometryCorrection.distortionRadius/geometryCorrection.pixelSize;    // omega_roll
			sv[i] = sym_vect[i]*1000.0*geometryCorrection.distortionRadius/geometryCorrection.pixelSize; // omega_rolls
		}
		for (int i = getIMUIndex()+3; i < getLength(); i++){
			v[i] = vector[i] *1000.0; // *distortionRadius/pixelSize;    // movement mm/s
			sv[i] = sym_vect[i]*1000.0; // *distortionRadius/pixelSize; // movement mm/s
		}
		double [] tilts =    new double[getNumSensors()];
		double [] azimuths = new double[getNumSensors()];
		double [] rolls =    new double[getNumSensors()];
		double [] zooms =    new double[getNumSensors()];
		int n = getNumSensors();
		for (int i = 0; i < getNumSensors();i++) {
			rolls[i] = v[getRollIndex()+i];
			if (i < (n-1)) {
				tilts[i] = v[getTiltIndex()+i];
				tilts[n-1] -= tilts[i];
				azimuths[i] = v[getAzimuthIndex()+i];
				azimuths[n-1] -= azimuths[i];
				zooms[i] = v[getZoomIndex()+i];
				zooms[n-1] -= zooms[i];
			}
		}
		s = "";
		s+= "tilt    (up):   "; for (int i = 0; i < n; i++) s += String.format(" %8.4fpx", tilts[i]);    s+=" (shift of the image center)\n"; 
		s+= "azimuth (right):"; for (int i = 0; i < n; i++) s += String.format(" %8.4fpx", azimuths[i]); s+=" (shift of the image center)\n"; 
		s+= "roll    (CW):   "; for (int i = 0; i < n; i++) s += String.format(" %8.4fpx", rolls[i]);    s+=" (shift at the image half-width from the center)\n"; 
		s+= "diff zoom (in): "; for (int i = 0; i < n; i++) s += String.format(" %8.4fpx", zooms[i]);    s+=" (shift at the image half-width from the center)\n";
		if (short_out) {
			return s;
		}
		s += "Symmetrical vector:\n";
		if (getNumSensors() == 4) { // Use arrows for quad camera only (but update to match new
			// ← → ↑ ↓ ⇖ ⇗ ⇘ ⇙ ↔ ↺ ↻ 
			//		s += "    |⇘ ⇙|     |⇘ ⇗|     |⇗ ⇘|     |⇙ ⇘|      |⇙ ⇗|     |⇖  ⇘| 6: common roll 7:(r0-r3)/2,           |- +|     |- -|     |- +|\n";
			//		s += " 0: |⇗ ⇖|  1: |⇙ ⇖|  2: |⇖ ⇙|  3: |⇖ ⇗|  4:  |⇗ ⇙|  5: |⇘  ⇖| 8:(r1-r2)/2    9:(r0+r3-r1-r2)/4  10: |- +| 11: |+ +| 12: |+ -|\n";
			s += "  |⇘ ⇙|  |⇗ ⇘|  |⇘ ⇗|  |⇗ ⇖|   |⇘ ⇖|  |⇗ ⇙|  |↻ ↻|  |↺ ↻|  |↺ ↺|  |↺ ↻|   |- +|     |- -|     |- +|\n";
			s += "0:|⇗ ⇖|1:|⇖ ⇙|2:|⇙ ⇖|3:|⇘ ⇙|4: |⇖ ⇘|5:|⇙ ⇗|6:|↻ ↻|7:|↺ ↻|8:|↻ ↻|9:|↻ ↺|10:|- +| 11: |+ +| 12: |+ -|\n";

			s += String.format(" 0:%9.6fpx 1:%8.5fpx 2:%8.5fpx 3:%8.5fpx 4:%8.5fpx 5:%8.5fpx 6:%8.5fpx 7:%8.5fpx 8:%8.5fpx 9:%8.5fpx 10: %8.5fpx 11:%8.5fpx 12:%8.5fpx\n" ,
					sv[0], sv[1], sv[2], sv[3], sv[4], sv[5], sv[6], sv[7], sv[8], sv[9], sv[10], sv[11], sv[12] );
			s += String.format("omega_tilt =  %9.4f pix/s, omega_azimuth = %9.4f pix/s, omega_roll = %9.4f pix/s\n", sv[13],  sv[14],  sv[15]);
			s += String.format("velocity Vx = %9.4f mm/s,             Vy = %9.4f mm/s,          Vz = %9.4f mm/s\n",  sv[16],  sv[17],  sv[18]);
/*
proto_type1.length=4 proto_type1_extra.length=0 proto_type2.length=2 proto_type2_extra.length=0 proto_type3.length=0
Vector #  0 [0 | 0 | 0 | 0] 1.00<1.00  [1.0  | 1.0  | 1.0  | 1.0 ]  l= -1.0 n=-1(-1) type1:0
Vector #  1 [1 | 1 | 1 | 1] 1.00<1.00  [1.0  | 1.0  | 1.0  | 1.0 ]  l= 1.0 n=5(5) type1:1
Vector #  2 [0 | 2 | 0 | 2] 1.00<1.00  [1.0  | 1.0  | 1.0  | 1.0 ]  l= 1.0 n=4(4) type1:2
Vector #  3 [1 | 3 | 1 | 3] 1.00<1.00  [1.0  | 1.0  | 1.0  | 1.0 ]  l= 1.0 n=3(3) type1:3
Vector #  4 [0 | 3 | 2 | 1] 1.00<1.00  [1.0  | 1.0  | 1.0  | 1.0 ]  l= 1.0 n=2(2) type2:0
Vector #  5 [1 | 0 | 3 | 2] 1.00<1.00  [1.0  | 1.0  | 1.0  | 1.0 ]  l= 1.0 n=1(1) type2:1
No more orthogonal vectors, current number of independent vectors = 6(of 6 needed)
Vector #  0 [+|+|+|+] 1.00<1.00  [1.0  | 1.0  | 1.0  | 1.0 ]  l= -1.0 n=-1(-1)
Vector #  1 [-|+|+|-] 1.00<1.00  [1.0  | 1.0  | 1.0  | 1.0 ]  l= 1.0 n=3(3)
Vector #  2 [-|-|+|+] 1.00<1.00  [1.0  | 1.0  | 1.0  | 1.0 ]  l= 1.0 n=2(2)
Vector #  3 [-|+|-|+] 1.00<1.00  [1.0  | 1.0  | 1.0  | 1.0 ]  l= 1.0 n=1(1)
↻ 	↺	
 */
		} else {
			// azimuth/tilts
			String [] rt4 = {"↓","↷","↑", "↶"};
			int [][] rt_proto = symmVectorsSet.dir_rt;
			int indx = 0;
			for (int i = 0; i < rt_proto.length; i++) {
				s += String.format("%2d: [" , indx);
				for (int j = 0; j < n; j++) {
					if ((j>0) && (j%(n/4) == 0)) {
						 s += "|";
					}
					s += rt4[rt_proto[i][j]];
				}
				s+=String.format("] %9.6f px\n", sv[indx]);
				indx++;
			}
			// rolls
			double [][] rot_proto = symmVectorsSet.rots;
			for (int i = 0; i < rot_proto.length; i++) {
				s += String.format("%2d: [" , indx);
				for (int j = 0; j < n; j++) {
					if ((j>0) && (j%(n/4) == 0)) {
						 s += "|";
					}
					if (rot_proto[i][j] > 0) { // Index 30 out of bounds for length 16
						s+="↻";
					} else if (rot_proto[i][j] < 0) {
						s+="↺";
					} else {
						s+="0";
					}
				}
				s+=String.format("] %9.6f px\n", sv[indx]);
				indx++;
			}
			// zooms
			double [][] zoom_proto = symmVectorsSet.zooms;
			for (int i = 0; i < zoom_proto.length; i++) {
				s += String.format("%2d: [" , indx);
				for (int j = 0; j < n; j++) {
					if ((j>0) && (j%(n/4) == 0)) {
						 s += "|";
					}
					if (zoom_proto[i][j] > 0) {
						s+="+";
					} else if (zoom_proto[i][j] < 0) {
						s+="-";
					} else {
						s+="0";
					}
				}
				s+=String.format("] %9.6f px\n", sv[indx]);
				indx++;
			}
		}
		return s;
	}
	
	public boolean editVector() {
		int n = getNumSensors();
		boolean use_tabs = (n>4);
		System.out.println(toString());
		String lines[] = toString().split("\\n");
		double par_scales_tlaz = 1000.0*geometryCorrection.focalLength/geometryCorrection.pixelSize;
		double par_scales_roll = 1000.0*geometryCorrection.distortionRadius/geometryCorrection.pixelSize;
		double par_scales_linear = 1000.0;
			GenericJTabbedDialog gd = new GenericJTabbedDialog("Update extrinsic parameters",1500,1100);
			for (String s:lines) {
				gd.addMessage(s);
			}
		if (use_tabs) gd.addTab("Tilts","Vertical (around horizontal axis perpendicular to te camera view direction) rotations in pixels for current camera parameters");
		for (int i = getTiltIndex(); i < getAzimuthIndex(); i++){
			gd.addNumericField("Tilt "+(i -getTiltIndex())+" (up)",                               par_scales_tlaz * vector[i],                10, 13,"pix",
					"Vertical (around horizontal axis perpendicular to te camera view direction) rotations in pixels for current camera parameters");
		}
		if (use_tabs)gd.addTab("Azimuths","Horizontal (around vertical axis) rotations in pixels for current camera parameters");
		for (int i = getAzimuthIndex(); i < getRollIndex(); i++){
			gd.addNumericField("Azimuth " + (i - getAzimuthIndex()) + " (right)",   par_scales_tlaz * vector[i],                10, 13,"pix",
					"Horizontal (around vertical axis) rotations in pixels for current camera parameters");
		}
		if (use_tabs) gd.addTab("Rolls","Roll around view axis rotations in pixels for current camera parameters, 1 pixel corresponds to FoV edges");
		for (int i = getRollIndex(); i < getZoomIndex(); i++){
			gd.addNumericField("Roll " + (i - getRollIndex()) + " (clockwise)",     par_scales_roll * vector[i] ,               10, 13,"pix",
					"Roll around view axis rotations in pixels for current camera parameters, 1 pixel corresponds to FoV edges");
		}
		if (use_tabs) gd.addTab("Zooms","Relative (to the average) focal length modifications, 1 pixel corresponds to FoV edges");
		for (int i = getZoomIndex(); i < getIMUIndex(); i++){
			gd.addNumericField("Zoom " + (i - getZoomIndex()) + " (zoom in)",       par_scales_roll * vector[i],                10, 13,"pix",
					"Relative (to the average) focal length modifications, 1 pixel corresponds to FoV edges");
		}
		if (use_tabs) gd.addTab("ERS","ERS correction parameters");
		else gd.addMessage("--- ERS correction parameters ---");

		gd.addNumericField("Angular rotation azimuth  (right)",                 par_scales_tlaz   * vector[getIMUIndex() + 0],  10, 13,"pix/s",
				"Horizonatal camera rotation (reasonable value ~1000 pix/s)");
		gd.addNumericField("Angular rotation  tilt (up)",                       par_scales_tlaz   * vector[getIMUIndex() + 1],  10, 13,"pix/s",
				"Vertical camera rotation (reasonable value ~1000 pix/s)");
		gd.addNumericField("Angular rotation  tilt (clockwise)",                par_scales_roll   * vector[getIMUIndex() + 2],  10, 13,"pix/s",
				"Roll camera rotation (reasonable value ~1000 pix/s)");
		gd.addNumericField("Linear velocity (right)",                           par_scales_linear * vector[getIMUIndex() + 3],  10, 13,"mm/s",
				"Linear movement right");
		gd.addNumericField("Linear velocity (up)",                              par_scales_linear * vector[getIMUIndex() + 4],  10, 13,"mm/s",
				"Linear movement up");
		gd.addNumericField("Linear velocity (forward)",                         par_scales_linear * vector[getIMUIndex() + 5],  10, 13,"mm/s",
				"Linear movement forward");
			gd.showDialog();
		if (gd.wasCanceled()) return false;

		for (int i = 0; i < getRollIndex(); i++){
			vector[i] = gd.getNextNumber()/par_scales_tlaz;
		}

		for (int i = getRollIndex(); i < getIMUIndex(); i++){
			vector[i] = gd.getNextNumber()/par_scales_roll;
		}

		vector[getIMUIndex() + 0] = gd.getNextNumber()/par_scales_tlaz;
		vector[getIMUIndex() + 1] = gd.getNextNumber()/par_scales_tlaz;
		vector[getIMUIndex() + 2] = gd.getNextNumber()/par_scales_roll;
		vector[getIMUIndex() + 3] = gd.getNextNumber()/par_scales_linear;
		vector[getIMUIndex() + 4] = gd.getNextNumber()/par_scales_linear;
		vector[getIMUIndex() + 5] = gd.getNextNumber()/par_scales_linear;
		return true;
	}


	public boolean editIMU() {
		double par_scales_tlaz = 1000.0*geometryCorrection.focalLength/geometryCorrection.pixelSize;
		double par_scales_roll = 1000.0*geometryCorrection.distortionRadius/geometryCorrection.pixelSize;
		double par_scales_linear = 1000.0;
			GenericJTabbedDialog gd = new GenericJTabbedDialog("Set camera angular velocities (in pix/sec)",800,300);
		gd.addNumericField("Angular rotation azimuth  (positive - rotate camera  the right)", par_scales_tlaz   * vector[getIMUIndex() + 0],  3,6,"pix/s",
				"Horizonatal camera rotation (reasonable value ~1000 pix/s)");
		gd.addNumericField("Angular rotation  tilt (positive - rotate camera up)",            par_scales_tlaz   * vector[getIMUIndex() + 1],  3,6,"pix/s",
				"Vertical camera rotation (reasonable value ~1000 pix/s)");
		gd.addNumericField("Angular rotation  tilt (positive - clockwise)",                   par_scales_roll   * vector[getIMUIndex() + 2],  3,6,"pix/s",
				"Roll camera rotation (reasonable value ~1000 pix/s)");
		gd.addNumericField("Linear velocity (positive - right)",                              par_scales_linear * vector[getIMUIndex() + 3],  3,6,"mm/s",
				"Linear movement right");
		gd.addNumericField("Linear velocity (positive - up)",                                 par_scales_linear * vector[getIMUIndex() + 4],  3,6,"mm/s",
				"Linear movement up");
		gd.addNumericField("Linear velocity (positive - forward)",                            par_scales_linear * vector[getIMUIndex() + 5],  3,6,"mm/s",
				"Linear movement forward");
			gd.showDialog();
		if (gd.wasCanceled()) return false;
		vector[getIMUIndex() + 0] = gd.getNextNumber()/par_scales_tlaz;
		vector[getIMUIndex() + 1] = gd.getNextNumber()/par_scales_tlaz;
		vector[getIMUIndex() + 2] = gd.getNextNumber()/par_scales_roll;
		vector[getIMUIndex() + 3] = gd.getNextNumber()/par_scales_linear;
		vector[getIMUIndex() + 4] = gd.getNextNumber()/par_scales_linear;
		vector[getIMUIndex() + 5] = gd.getNextNumber()/par_scales_linear;
		return true;
	}


	// returns false if any component is NaN, in that case do not increment
	public boolean incrementVector(double [] incr, // USED in lwir
			double scale)
	{
		for (int i = 0; i < incr.length; i++){
			if (Double.isNaN(vector[i])) return false;
		}
		for (int i = 0; i < incr.length; i++){
			vector[i]+= incr[i] * scale;
		}
		return true;
	}

	public boolean incrementVector(CorrVector incr, double scale) // USED in lwir
	{
		return incrementVector(incr.toArray(), scale);
	}

	public CorrVector diffFromVector(CorrVector other_vector) {
		double [] aother = other_vector.toArray();
		double [] athis = toArray().clone();
		for (int i = 0; i < athis.length; i++) {
			athis[i] -= aother[i];
		}
		return new CorrVector(geometryCorrection, athis);
	}

	public double getNorm() {
		double s2 = 0;
		for (int i = 0; i < vector.length; i++) {
			double v = vector[i];
			// manually reducing weights of ERS parameters
			if (i >= (getIMUIndex()+3)) { 
				v *= 0.001; // imu mm/s
			} else if (i >= getIMUIndex()) {
				v *= 0.01; // imu mm/s
			}  
			s2 += v*v;
		}
		return Math.sqrt(s2); // add weights to compare apples and oranges?
	}


	@Override
	public CorrVector clone(){ // not used in lwir
		return new CorrVector(geometryCorrection, this.vector.clone());
	}

	/**
	 * Convert manual pixel shift between images to azimuth/tilt rotations (distortions ignored)
	 * and apply (add) them to the current vector (normally should be all 0.0)
	 * @param pXY_shift manula XY pixel corrections (shiftXY made of clt_parameters.fine_corr_[xy]_[0123])
	 */
	public void applyPixelShift(double [][] pXY_shift){ // not used in lwir
		double [] pXY_avg = {0.0,0.0};
		for (int i = 0; i < geometryCorrection.numSensors; i++){
			for (int j = 0; j < 2; j++) {
				pXY_avg[j] += pXY_shift[i][j]/geometryCorrection.numSensors;
			}
		}
		double pix_to_rads = (0.001*geometryCorrection.pixelSize)/geometryCorrection.focalLength;
		for (int i = 0; i < (geometryCorrection.numSensors -1); i++){
			vector[getTiltIndex()  + i]    -= pix_to_rads * (pXY_shift[i][1] - pXY_avg[1]);
			vector[getAzimuthIndex() + i] += pix_to_rads * (pXY_shift[i][0] - pXY_avg[0]);
		}
	}

	/**
	 * Get conversion matrix from symmetrical coordinates
	 * @return
	 *
	 *     |⇘ ⇙|     |⇘ ⇗|     |⇗ ⇘|     |⇙ ⇘|      |⇙ ⇗|     |⇖  ⇘|
	 *  0: |⇗ ⇖|  1: |⇙ ⇖|  2: |⇖ ⇙|  3: |⇖ ⇗|  4:  |⇗ ⇙|  5: |⇘  ⇖|
	 *
	 */
	/*
	public double [][] dSym_j_dTar_i() // USED in lwir
	{
		double [][] tar_to_sym = {
				{-2.0, -2.0,  2.0, -2.0,  0.0,  0.0, 0.0,   0.0,  0.0,  0.0,   0.0,  0.0,  0.0,   0.0,  0.0,  0.0,  0.0,  0.0,  0.0}, // t0
				{-2.0,  0.0,  0.0, -2.0,  2.0, -2.0, 0.0,   0.0,  0.0,  0.0,   0.0,  0.0,  0.0,   0.0,  0.0,  0.0,  0.0,  0.0,  0.0}, // t1
				{ 0.0, -2.0,  2.0,  0.0,  2.0, -2.0, 0.0,   0.0,  0.0,  0.0,   0.0,  0.0,  0.0,   0.0,  0.0,  0.0,  0.0,  0.0,  0.0}, // t2

				{ 2.0,  2.0,  2.0, -2.0,  0.0,  0.0, 0.0,   0.0,  0.0,  0.0,   0.0,  0.0,  0.0,   0.0,  0.0,  0.0,  0.0,  0.0,  0.0}, // a0
				{ 0.0,  2.0,  2.0,  0.0,  2.0,  2.0, 0.0,   0.0,  0.0,  0.0,   0.0,  0.0,  0.0,   0.0,  0.0,  0.0,  0.0,  0.0,  0.0}, // a1
				{ 2.0,  0.0,  0.0, -2.0,  2.0,  2.0, 0.0,   0.0,  0.0,  0.0,   0.0,  0.0,  0.0,   0.0,  0.0,  0.0,  0.0,  0.0,  0.0}, // a2

				{ 0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.25,  0.5,  0.0,  0.25,  0.0,  0.0,  0.0,   0.0,  0.0,  0.0,  0.0,  0.0,  0.0}, // roll 0
				{ 0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.25,  0.0,  0.5, -0.25,  0.0,  0.0,  0.0,   0.0,  0.0,  0.0,  0.0,  0.0,  0.0}, // roll 1
				{ 0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.25,  0.0, -0.5, -0.25,  0.0,  0.0,  0.0,   0.0,  0.0,  0.0,  0.0,  0.0,  0.0}, // roll 2
				{ 0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.25, -0.5,  0.0,  0.25,  0.0,  0.0,  0.0,   0.0,  0.0,  0.0,  0.0,  0.0,  0.0}, // roll 3

				{ 0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.0,   0.0,  0.0,  0.0,  -1.0, -1.0,  0.0,   0.0,  0.0,  0.0,  0.0,  0.0,  0.0}, // scale 0
				{ 0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.0,   0.0,  0.0,  0.0,   0.0, -1.0,  1.0,   0.0,  0.0,  0.0,  0.0,  0.0,  0.0}, // scale 1
				{ 0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.0,   0.0,  0.0,  0.0,  -1.0,  0.0,  1.0,   0.0,  0.0,  0.0,  0.0,  0.0,  0.0},  // scale 2

				{ 0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.0,   0.0,  0.0,  0.0,   0.0,  0.0,  0.0,   1.0,  0.0,  0.0,  0.0,  0.0,  0.0}, // omega_tilt
				{ 0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.0,   0.0,  0.0,  0.0,   0.0,  0.0,  0.0,   0.0,  1.0,  0.0,  0.0,  0.0,  0.0}, // omega_azimuth
				{ 0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.0,   0.0,  0.0,  0.0,   0.0,  0.0,  0.0,   0.0,  0.0,  1.0,  0.0,  0.0,  0.0},  // omega_roll

				{ 0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.0,   0.0,  0.0,  0.0,   0.0,  0.0,  0.0,   0.0,  0.0,  0.0,  1.0,  0.0,  0.0}, // velocity_x
				{ 0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.0,   0.0,  0.0,  0.0,   0.0,  0.0,  0.0,   0.0,  0.0,  0.0,  0.0,  1.0,  0.0}, // velocity_y
				{ 0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.0,   0.0,  0.0,  0.0,   0.0,  0.0,  0.0,   0.0,  0.0,  0.0,  0.0,  0.0,  1.0}  // velocity_z
		};

		return tar_to_sym;
	}
	public double [][] dTar_j_dSym_i(){
		//			Matrix sym_to_tar = new Matrix(symToTARArray());
		//			Matrix tar_to_sym = sym_to_tar.inverse();
		//			return tar_to_sym.getArray();
		 */
	
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
	/*
		double [][] sym_to_tar=	{ // USED in lwir
				// t0     t1     t2     a0     a1     a2     r0    r1    r2    r3    s0    s1    s2
				{-0.125,-0.125, 0.125, 0.125,-0.125, 0.125, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,  0.0,   0.0,  0.0,  0.0 },  // sym0
		        {-0.125, 0.125,-0.125, 0.125, 0.125,-0.125, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,  0.0,   0.0,  0.0,  0.0 },  // sym1
		        { 0.125,-0.125, 0.125, 0.125, 0.125,-0.125, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,  0.0,   0.0,  0.0,  0.0 },  // sym2
		        {-0.125,-0.125, 0.125,-0.125, 0.125,-0.125, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,  0.0,   0.0,  0.0,  0.0 },  // sym3
		        {-0.125, 0.125, 0.125,-0.125, 0.125, 0.125, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,  0.0,   0.0,  0.0,  0.0 },  // sym4
		        { 0.125,-0.125,-0.125,-0.125, 0.125, 0.125, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,  0.0,   0.0,  0.0,  0.0 },  // sym5
		        { 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 1.0,  1.0,  1.0,  1.0,  0.0,  0.0,  0.0,   0.0,  0.0,  0.0,   0.0,  0.0,  0.0 },  // sym6  = (r0+r1+r2+r3)/4
		        { 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 1.0,  0.0,  0.0, -1.0,  0.0,  0.0,  0.0,   0.0,  0.0,  0.0,   0.0,  0.0,  0.0 },  // sym7  = (r0-r3)/2
		        { 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0,  1.0, -1.0,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,  0.0,   0.0,  0.0,  0.0 },  // sym8  = (r1-r2)/2
		        { 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 1.0, -1.0, -1.0,  1.0,  0.0,  0.0,  0.0,   0.0,  0.0,  0.0,   0.0,  0.0,  0.0 },  // sym9  = (r0+r3-r1-r2)/4
		        { 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0,  0.0, -1.0,  1.0, -0.5,  0.5, -0.5,   0.0,  0.0,  0.0,   0.0,  0.0,  0.0 },  // sym10 = -s0 - s2
		        { 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0,  0.0, -1.0,  1.0, -0.5, -0.5,  0.5,   0.0,  0.0,  0.0,   0.0,  0.0,  0.0 },  // sym11 = -s0 - s1
		        { 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0,  0.0, -1.0,  1.0, -0.5,  0.5,  0.5,   0.0,  0.0,  0.0,   0.0,  0.0,  0.0 },  // sym12 =  s1 + s2
				{ 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,   1.0,  0.0,  0.0,   0.0,  0.0,  0.0 },  // sym13 =  omega_tilt
				{ 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,   0.0,  1.0,  0.0,   0.0,  0.0,  0.0 },  // sym14 =  omega_azimuth
				{ 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,  1.0,   0.0,  0.0,  0.0 },  // sym15 =  omega_roll
				{ 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,  0.0,   1.0,  0.0,  0.0 },  // sym16 =  velocity_x
				{ 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,  0.0,   0.0,  1.0,  0.0 },  // sym17 =  velocity_y
				{ 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,  0.0,   0.0,  0.0,  1.0 }   // sym18 =  velocity_z
		        } ;
		return sym_to_tar;
	}
*/
	public boolean [] getParMask( // USED in lwir
			boolean use_disparity,
			boolean use_aztilts,       // Adjust azimuths and tilts excluding disparity
			boolean use_diff_rolls,    // Adjust differential rolls (3 of 4 angles)
			boolean common_roll,
			boolean corr_focalLength,
			boolean ers_rot,           // Enable ERS correction of the camera rotation
			boolean ers_forw,      // Enable ERS correction of the camera linear movement in z direction
			boolean ers_side,      // Enable ERS correction of the camera linear movement in x direction
			boolean ers_vert,      // Enable ERS correction of the camera linear movement in y direction
	  		int     manual_par_sel)    // Manually select the parameter mask bit 0 - sym0, bit1 - sym1, ... (0 - use boolean flags, != 0 - ignore boolean flags)
	{
//		SymmVectorsSet
		int num_aztilts = symmVectorsSet.rt.length;
		int num_rolls =   symmVectorsSet.rots.length;
		int num_zooms =   symmVectorsSet.zooms.length;
		int num_imus =    6;
		boolean [] par_mask = new boolean[num_aztilts + num_rolls + num_zooms + num_imus];
		par_mask[0] =                                       use_disparity;
		for (int i = 1; i < num_aztilts; i++) par_mask[i] = use_aztilts;
		par_mask[num_aztilts] =                             common_roll;
		for (int i = 1; i < num_rolls; i++)   par_mask[i + num_aztilts] = use_diff_rolls;
		
		for (int i = 1; i < num_zooms; i++)   par_mask[i + num_aztilts + num_rolls] = corr_focalLength;
		int imu_offset = num_aztilts + num_rolls + num_zooms;
		par_mask[imu_offset + 0] = ers_rot;
		par_mask[imu_offset + 1] = ers_rot;
		par_mask[imu_offset + 2] = ers_rot;
		par_mask[imu_offset + 3] = ers_side;
		par_mask[imu_offset + 4] = ers_vert;
		par_mask[imu_offset + 5] = ers_forw;
		/*
		boolean [] par_mask = {
				use_disparity,    //sym0
				use_aztilts,      //sym1
				use_aztilts,      //sym2
				use_aztilts,      //sym3
				use_aztilts,      //sym4
				use_aztilts,      //sym5
				common_roll,      //sym6 // common roll
				use_diff_rolls,   //sym7
				use_diff_rolls,   //sym8
				use_diff_rolls,   //sym9
				corr_focalLength, //sym10
				corr_focalLength, //sym11
				corr_focalLength, //sym12
				ers_rot,          //sym13
				ers_rot,          //sym14
				ers_rot,          //sym15
				ers_side,         //sym16
				ers_vert,         //sym17
				ers_forw          //sym18
		};
		*/
		if (manual_par_sel != 0) { // not used in lwir
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
	 * are only responsible for the "lazy eye" (last 4  were roll, now differential zooms are added):
	 *
	 *     |⇘ ⇙|     |⇘ ⇗|     |⇗ ⇘|     |⇙ ⇘|      |⇙ ⇗|     |⇖  ⇘|
	 *  0: |⇗ ⇖|  1: |⇙ ⇖|  2: |⇖ ⇙|  3: |⇖ ⇗|  4:  |⇗ ⇙|  5: |⇘  ⇖|
	 *
	 * @param port_coord_deriv result of getPortsCoordinatesAndDerivatives: first index: port0x, port0y... port3y,
	 *                         second index - parameters [tilt0, tilt1, ..., roll3]
	 * @param par_mask array of 10->13->16 elements - which parameters to use (normally all true or all but first
	 *
	 *  UPDATE
	 * @return
	 */
	public double [][] getJtPartial( // USED in lwir
			double [][] port_coord_deriv,
			boolean []  par_mask)
	{
		Matrix from_sym = symmVectorsSet.from_sym;
		if (par_mask != null) {
			ArrayList<Integer> par_list = new ArrayList<Integer>();
			for (int i = 0; i < par_mask.length; i++) {
				if (par_mask[i]) {
					par_list.add(i);
				}
			}
			if (par_list.size() < par_mask.length) {
				int[] sub_pars = new int [par_list.size()];
				for (int i = 0; i < sub_pars.length; i++) {
					sub_pars[i] = par_list.get(i); 
				}
//				int[] sub_pars = par_list.stream().mapToInt(i -> i).toArray(); // Java8
				from_sym = from_sym.getMatrix(0, par_mask.length-1,sub_pars); // remove masked columns
			}
		}
		Matrix pcd = new Matrix (port_coord_deriv); // rows: px0,py0,... px[n-1], py[n-1], columns: tar
		Matrix mjt_part = pcd.times(from_sym).transpose(); // rows: sym0..., columns: px0,py0,... px[n-1], py[n-1]
		return mjt_part.getArray();
		/*
		int num_pars = 0;
		for (int npar = 0; npar < getLength(); npar++) if ((par_mask==null) || par_mask[npar]) {
			num_pars++;
		}
		double [][] sym_to_tar=	dTar_j_dSym_i();
		double [][] jt_part = new double[num_pars][port_coord_deriv.length];
		int opar = 0;
		for (int npar = 0; npar < getLength(); npar++) if ((par_mask==null) || par_mask[npar]) {
			for (int i = 0; i < port_coord_deriv.length; i++){ // pxy index (0..7)
				for (int k = 0; k < sym_to_tar[npar].length; k++){
					jt_part[opar][i] += sym_to_tar[npar][k]* port_coord_deriv[i][k]; // transposing port_coord_deriv
				}
			}
			opar++;
		}
		return jt_part;
		*/
	}

	// convert tilt0,... roll3 array to symmetrical coordinates [0] - to the center (disparity)
	// Use   public void setMatrix (int[] r, int[] c, Matrix X) ? 

	// symmVectorsSet.from_sym, symmVectorsSet.to_sym;
	public double [] toSymArray(boolean [] par_mask) // USED in lwir
	{
		return toSymArray(this.vector, par_mask);
	}
    /*
	public double [] toSymArray( // USED in lwir
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
		for (int npar = 0; npar < getLength(); npar++) if ((par_mask==null) || par_mask[npar]) {
			for (int i = 0; i < tar_array.length; i++){
				sym_array[opar] += tar_array[i] * tar_to_sym[i][npar];
			}
			opar++;
		}
		return sym_array;
	}

	public double [] toTarArray( // USED in lwir
			double [] sym_array,
			boolean [] par_mask)
	{
		double [][] sym_to_tar = dTar_j_dSym_i();
		double [] tar_array = new double  [sym_to_tar[0].length];
		for (int npar = 0; npar < getLength(); npar++)  {
			int spar = 0;
			for (int i = 0; i < getLength(); i++) if ((par_mask==null) || par_mask[i]){
				tar_array[npar] += sym_array[spar] * sym_to_tar[i][npar];
				spar++;
			}
		}
		return tar_array;
	}
	*/
	public double [] toSymArray( // USED in lwir
			double [] tar_array,
			boolean [] par_mask)
	{
		Matrix to_sym = symmVectorsSet.to_sym;
		if (par_mask != null) {
			ArrayList<Integer> par_list = new ArrayList<Integer>();
			for (int i = 0; i < par_mask.length; i++) {
				if (par_mask[i]) {
					par_list.add(i);
				}
			}
			if (par_list.size() < par_mask.length) {
				int[] sub_pars = par_list.stream().mapToInt(i -> i).toArray(); // Java8
				to_sym = to_sym.getMatrix(sub_pars, 0, par_mask.length-1); // remove masked rows
			}
		}
		Matrix tar = new Matrix(tar_array,tar_array.length);
		Matrix sym = to_sym.times(tar);
		return sym.getColumnPackedCopy();
	}

	public double [] toTarArray( // USED in lwir
			double [] sym_array,
			boolean [] par_mask)
	{
		Matrix from_sym = symmVectorsSet.from_sym;
		if (par_mask != null) {
			ArrayList<Integer> par_list = new ArrayList<Integer>();
			for (int i = 0; i < par_mask.length; i++) {
				if (par_mask[i]) {
					par_list.add(i);
				}
			}
			if (par_list.size() < par_mask.length) {
				int[] sub_pars = par_list.stream().mapToInt(i -> i).toArray(); // Java8
				from_sym = from_sym.getMatrix(0, par_mask.length-1,sub_pars); // remove masked columns
			}
		}
		Matrix sym = new Matrix(sym_array,sym_array.length);
		Matrix tar = from_sym.times(sym);
		return tar.getColumnPackedCopy();
	}
	
}
