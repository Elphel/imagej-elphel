package com.elphel.imagej.tileprocessor;

import java.util.ArrayList;
import java.util.HashMap;

/**
 **
 ** SymmVector - Generating symmetrical geometry correction vectors for multi-camera
 ** circular configurations.
 **
 ** Copyright (C) 2021 Elphel, Inc.
 **
 ** -----------------------------------------------------------------------------**
 **
 **  SymmVector.java is free software: you can redistribute it and/or modify
 **  it under the terms of the GNU General Public License as published by
 **  the Free Software Foundation, either version 3 of the License, or
 **  (at your option) any later version.
 **
 **  This program is distributed in the hope that it will be useful,
 **  but WITHOUT ANY WARRANTY; without even the implied warranty of
 **  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See theE
 **  GNU General Public License for more details.
 **
 **  You should have received a copy of the GNU General Public License
 **  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 ** -----------------------------------------------------------------------------**
 */

public class SymmVector {
	public static int R_MINUS = 0; // radial inwards
	public static int T_PLUS =  1; // tangential clockwise
	public static int R_PLUS =  2; // radial outwards
	public static int T_MINUS = 3; // tangential counter-clockwise
	public static double [][] AX_DIR = {{-1.0, 0.0}, {0.0, 1.0}, {1.0,0.0}, {0.0,-1.0}};
	public final  int N;  // For now assuming N % 4 == 0,
	public boolean full_type1 = false; 
	public boolean full_type2 = false; 
	public double fail_length = 1E-6; // minimal new vector height
	public double best_delta =  1E-6; // count almost same candidates that differ by less than this
	
	private int [][] proto_type1 =       new int [0][];        // second half is the same as first half, this balances sum(x) and sum(y)
	                                     // all are either R or T
	private int [][] proto_type1_extra = new int [0][];  // R and T are mixed, not included in proto_type1
	private int [][] proto_type2 =       new int [0][];        // first quarter - any, second quarter - T2=R1,R2=-T1. 
	                                     // Third quarter is R3 =-R1, T3=-T1. Fourth quarter T4=R3,R4=-T3 
	private int [][] proto_type2_extra = new int [0][];  // Similar to proto_type2, but quarter 3 is independent of quarter 1
	                                     // Does not include those of proto_type2.
	
	private int [][] proto_type3 =       new int [0][];        // non-quad 
	
	private int [][] proto_all;          //

	private double [][] dvectors;
	private int [] sym_indices;
	private int num_defined;
	private boolean [] used_indices;
	private double [] cumul_influences; // did not work - all candidates correlate the same
	
	private int [][]    proto_rz;       //+1 - zoom in, rot CW, -1 - zoom out, rot CCW, 0 - error (not initialized)
	private double [][] rz_vectors;     // normalized rot/zoom vectors
	private int []      rz_indices;     // indices of selected rz_vectors
	private int         num_rz_defined; // number of defined rz_indices
	private boolean []  used_rz_indices;// already used rz_vectors
	private double []   cumul_rz_influences; // will not work? - all candidates correlate the same
	public int debug_level =          -1;
	
	// caching results of vectors generation for used camera configuration not to re-generate each time
	// CorrVector is constructed
	public static HashMap <Integer,SymmVectorsSet> vectors_cache = new HashMap <Integer,SymmVectorsSet>();

	public static SymmVectorsSet getSymmVectorsSet (int num_cameras) {
		SymmVectorsSet rvs = vectors_cache.get(num_cameras);
		if (rvs == null) {
			rvs = newVectors (num_cameras);
			vectors_cache.put(num_cameras, rvs);
		}
		return rvs;
	}
	
	public static double [][] getSymmXY(int num_cameras){
		SymmVectorsSet rvs = vectors_cache.get(num_cameras);
		if (rvs == null) {
			rvs = newVectors (num_cameras);
			vectors_cache.put(num_cameras, rvs);
		}
		return rvs.xy;
	}
	
	public static SymmVectorsSet newVectors (int num_cameras) {
		boolean full_type1 = false;
		boolean full_type2 = false;
		int debug_level = -1;
		SymmVectorsSet rvs = new SymmVectorsSet();
		SymmVector sv = new SymmVector(
				num_cameras,
				full_type1,
				full_type2,
				debug_level);
		rvs.xy =     sv.exportXY();
		rvs.rt =     sv.exportRT();
		rvs.dir_rt = sv.exportDirRT();
		rvs.rots =   sv.exportRZ(false); // include common roll
		rvs.zooms =  sv.exportRZ(true);  // no common zoom
		return rvs;
	}
	
	public SymmVector (
			int num_cameras,
			boolean full_type1,  // false - all R or all T, true - mixed
			boolean full_type2,// false - quarter 3 is negated quarter 1, true - independent
			int debug_level) {
		this.N = num_cameras;
		this.full_type1 = full_type1;
		this.full_type2 = full_type2;
		this.debug_level = debug_level;
		generateType1();
		if (full_type1) generateType1Extra();
		generateType2();
		if (full_type2) generateType2Extra(full_type1);
		if (N >= 8) {
			generateType3(-1); // so far - only for 8, just generate length for proto_type3
		}
		if (debug_level > 0) System.out.println("=========== N = "+N+" ===========");

		if (debug_level > 0) System.out.println(
				"proto_type1.length="+        proto_type1.length+
				" proto_type1_extra.length="+ proto_type1_extra.length+
				" proto_type2.length="+       proto_type2.length+
				" proto_type2_extra.length="+ proto_type2_extra.length+
				" proto_type3.length="+       proto_type3.length);
		int num_all =               proto_type1.length +
				                    proto_type2.length +
				      (full_type1 ? proto_type1_extra.length:0) +
				      (full_type2 ? proto_type2_extra.length:0);
		int type3_offs = num_all; 
		num_all = num_all + proto_type3.length; 
		proto_all = new int [num_all][];
		
		int indx = 0;
		for (int i = 0; i < proto_type1.length; i++) {
			proto_all[indx++] = proto_type1[i];
		}
		for (int i = 0; i < proto_type2.length; i++) {
			proto_all[indx++] = proto_type2[i];
		}
		if (full_type1) for (int i = 0; i < proto_type1_extra.length; i++) {
			proto_all[indx++] = proto_type1_extra[i];
		} 
		if (full_type2) for (int i = 0; i < proto_type2_extra.length; i++) {
			proto_all[indx++] = proto_type2_extra[i];
		}
		for (int i = 0; i < proto_type3.length; i++) {
			proto_all[indx++] = proto_type3[i];
		}
		used_indices = new boolean [num_all];
		sym_indices =  new int[2*N - 2];
		num_defined = 0;
		dvectors = new double[num_all][2*N];
		double scale = 1.0/Math.sqrt(N); // all vectors have the same length
		for (int i = 0; i < num_all; i++) {
			for (int j = 0; j < N; j++) {
				if ((dvectors[i] == null) || (proto_all[i] == null)) {
					System.out.println("null pointer");
				}
				dvectors[i][2*j] =     scale * AX_DIR[proto_all[i][j]][0]; // null
				dvectors[i][2*j + 1] = scale * AX_DIR[proto_all[i][j]][1];
			}
		}
		if (proto_type3.length > 0) {
			generateType3(type3_offs);
		}
		// re-normalize all vectors (remove earlier normalization)
		for (int i = 0; i < num_all; i++) {
			double l2 = 0.0;
			for (int j = 0; j < 2 *N; j++) {
				l2 += dvectors[i][j] * dvectors[i][j]; 
			}
			scale = 1.0/Math.sqrt(l2);
			for (int j = 0; j < 2 *N; j++) {
				dvectors[i][j] *= scale; 
			}
		}
		boolean use_min_influence = true;
		cumul_influences = new double[N];
		for (int ivect = 0; ivect <= sym_indices.length; ivect++) { // <= to check for impossible
			int best_index = 0;
			int num_best0 = -1;
			int num_best  = -1;
			double best_height = -1.0;
//			double best_metrics = -1.0;
			if (ivect == 0) {
				best_index = 0;
			} else {
				double [] new_heights = getNewHeights(); // normailize?
				double [] new_mins =    use_min_influence ? getNewMins() : null;
//				double [] metrics = new_heights.clone();
//				best_index =        bestVector(metrics);
//				best_index =        bestVector(best_delta, new_heights, new_mins);
				best_index =        bestVector(best_delta, new_heights, best_delta, new_mins, cumul_influences, false);
				
				if (ivect == 1) {
					if (debug_level > -1) System.out.println("Vector # "+ivect+": overwriting best_index= "+best_index+" with 1");
					best_index = 1; // overwrite
				}
				best_height =    new_heights[best_index];
//				best_metrics =   metrics[best_index];
				num_best0 =      getNumBest(best_delta, new_heights);
				num_best =       getNumBest(best_delta, new_heights, best_delta, new_mins);
				
				if (best_height < fail_length) {
					if (ivect < sym_indices.length) {
						System.out.println("Failed to build orthogonal system, current number of independent vectors = "+
								num_defined+ "(of "+sym_indices.length+ " needed)");
						break; // continue; // return;
					} else {
						System.out.println("No more orthogonal vectors, current number of independent vectors = "+
								num_defined+ "(of "+sym_indices.length+ " needed)");
						break; // continue; // return;
						
					}
				}
			}
			if (debug_level > 0) {
				System.out.print(String.format("Vector # %2d",ivect)+" [");
				for (int j=0; j < N; j++) {
					System.out.print(proto_all[best_index][j]);
					if ((j+1)%(N/4) != 0) {
//						System.out.print(", ");
					} else if (j < (N-1)) {
						System.out.print(" | ");
					} else {
						System.out.print("] ");
					}
				}
				double [] sXY = getXY(dvectors[best_index]);
				double [] ni =  getNormInfluence(dvectors[best_index]);
				double [] mmi = minMaxInfluence(dvectors[best_index]);
				System.out.print(String.format("%4.2f<%4.2f ", mmi[0],mmi[1]));
				System.out.print(" [");
				for (int j=0; j < N; j++) {
					System.out.print(String.format("%3.1f ", ni[j]));
					if ((j+1)%(N/4) != 0) {
//						System.out.print(", ");
					} else if (j < (N-1)) {
						System.out.print(" | ");
					} else {
						System.out.print("] ");
					}
				}
				System.out.print(" l= " + best_height);
				System.out.print(String.format(" n=%d(%d)", num_best,num_best0));
				System.out.print(" "+getType(best_index));
				if (debug_level > 1) {
					System.out.print(" sXY=["+sXY[0]+","+sXY[1]+"]");
				}
				System.out.println();
			}
			add_vector(best_index);
		}
		
		// Generate vectors for rot and zoom (same)
		generateTypeRZ(); // generate candidates as both proto_rot_zoom and rz_vectors
		used_rz_indices = new boolean [num_all];
		rz_indices =      new int[N]; // for zoom the zero will be removed
		num_rz_defined = 0;
		// Not needed - re-normalize all vectors (remove earlier normalization)
		
		boolean use_min_influence_rz = true;
		cumul_rz_influences = new double[N];
		for (int ivect = 0; ivect <= rz_indices.length; ivect++) { // <= to check for impossible
			int best_index = 0;
			int num_best0 = -1;
			int num_best  = -1;
			double best_height = -1.0;
			if (ivect == 0) {
				best_index = 0;
			} else {
				double [] new_heights = getNewHeightsRZ();
				double [] new_mins =    use_min_influence_rz ? getNewMinsRZ() : null;
				best_index =        bestVector(best_delta, new_heights, best_delta, new_mins, cumul_influences, true);
				
				if (ivect < 4) {
					if (debug_level > -1) System.out.println("Vector # "+ivect+": overwriting best_index= "+best_index+" with "+ivect);
					best_index = ivect; // overwrite with pre-defined vectors
				}
				best_height =    new_heights[best_index];
				num_best0 =      getNumBest(best_delta, new_heights);
				num_best =       getNumBest(best_delta, new_heights, best_delta, new_mins);
				
				if (best_height < fail_length) {
					if (ivect < rz_indices.length) {
						System.out.println("Failed to build orthogonal system, current number of independent vectors = "+
								num_rz_defined+ "(of "+rz_indices.length+ " needed)");
						break; // continue; // return;
					} else {
						System.out.println("No more orthogonal vectors, current number of independent vectors = "+
								num_rz_defined+ "(of "+rz_indices.length+ " needed)");
						break; // continue; // return;
						
					}
				}
			}
			if (debug_level > 0) {
				System.out.print(String.format("Vector # %2d",ivect)+" [");
				for (int j=0; j < N; j++) {
					System.out.print((proto_rz[best_index][j] > 0)?("+"):((proto_rz[best_index][j] < 0) ? "-":"0"));
					if ((j+1)%(N/4) != 0) {
//						System.out.print(", ");
					} else if (j < (N-1)) {
						System.out.print("|");
					} else {
						System.out.print("] ");
					}
				}
				double s = 0.0;
				for (int i = 0; i < N; i++) {
					s+=rz_vectors[best_index][i];
				}
				
				double [] ni =  getNormInfluenceRZ(rz_vectors[best_index]);
				double [] mmi = minMaxInfluenceRZ (rz_vectors[best_index]);
				System.out.print(String.format("%4.2f<%4.2f ", mmi[0],mmi[1]));
				System.out.print(" [");
				for (int j=0; j < N; j++) {
					System.out.print(String.format("%3.1f ", ni[j]));
					if ((j+1)%(N/4) != 0) {
//						System.out.print(", ");
					} else if (j < (N-1)) {
						System.out.print(" | ");
					} else {
						System.out.print("] ");
					}
				}
				System.out.print(" l= " + best_height);
				System.out.print(String.format(" n=%d(%d)", num_best,num_best0));
//				System.out.print(" "+getType(best_index));
				if (debug_level > 1) {
					System.out.print(" s="+s);
				}
				System.out.println();
			}
			add_vector_rz(best_index);
		}
	}

	/**
	 * Return set of generated orthonormal symmetrical vectors X,Y components (X - right, Y - down).
	 * {x0,y0,x1,y1, ... x[N-1],y[N-1]}
	 * Cameras positioned on R=1 circle, #0 - up, others - clockwise looking where camera is pointed
	 * vector [0] - disparity, all po9inted to the center (similar to quad),
	 * vector [1] - tangential only, clockwise.  	
	 * @return
	 */
	public double [][] exportXY(){
		double [][] rslt = new double[sym_indices.length][2*N];
		for (int n = 0; n < sym_indices.length; n++){
			double [] vect = dvectors[sym_indices[n]]; // should be normalized
			for (int i = 0; i < N; i++) {
				double alpha = i*2*Math.PI/N;
				rslt[n][2*i + 0] =  vect[2*i + 0]*Math.sin(alpha) + vect[2*i + 1]*Math.cos(alpha);
				rslt[n][2*i + 1] = -vect[2*i + 0]*Math.cos(alpha) + vect[2*i + 1]*Math.sin(alpha);
			}
		}
		return rslt;
	}

	/**
	 * Return set of generated orthonormal symmetrical vectors radial/ tangential components (radial - outwards, X - right, Y - down) 	
	 * {r0,t0,r1,t1, ... r[N-1],t[N-1]}
	 * Cameras positioned on R=1 circle, #0 - up, others - clockwise looking where camera is pointed  	
	 * vector [0] - disparity, all po9inted to the center (similar to quad),
	 * vector [1] - tangential only, clockwise.  	
	 * @return
	 */
	public double [][] exportRT(){
		double [][] rslt = new double[sym_indices.length][];
		for (int n = 0; n < sym_indices.length; n++){
			rslt[n]= dvectors[sym_indices[n]]; // should be normalized
		}
		return rslt;
	}
	
	/**
	 * Export directions (0 - radial to center, 1 - tangential CW, 2 - radial out , 3 tangential CW
	 * @return
	 */
	public int [][] exportDirRT(){
		int [][] rslt = new int[sym_indices.length][];
		for (int n = 0; n < sym_indices.length; n++){
			rslt[n]= proto_all[sym_indices[n]];
		}
		return rslt;
	}
	
	/**
	 * 
	 * @param zoom_mode
	 * @return
	 */
	public double [][] exportRZ(boolean zoom_mode){
		int offs = zoom_mode ? 1 : 0;
		double [][] rslt = new double[rz_indices.length - offs][];
		for (int n = 0; n < rslt.length; n++) {
			rslt[n]= rz_vectors[rz_indices[n + offs]]; // should be normalized
		}
		return rslt;
	}
	
	
	
	
	/**
	 * To verify that center of mass does not move
	 * @param vect alternating radial/tangential components, length= 2*N 
	 * @return {Sx, Sy}
	 */
	public double [] getXY(double [] vect) {
		double [] SXY = {0.0,0.0};
		for (int i = 0; i < N; i++ ) {
			double alpha = i*2*Math.PI/N;
			SXY[0] +=  vect[2*i + 0]*Math.sin(alpha) + vect[2*i + 1]*Math.cos(alpha);
			SXY[1] += -vect[2*i + 0]*Math.cos(alpha) + vect[2*i + 1]*Math.sin(alpha);
		}
		return SXY;
	}
	
	public double [] getNormInfluence(double [] vect) {
		double l2 = 0.0;
		for (int i = 0; i < N; i++) {
			l2 += vect[2*i]*vect[2*i]+vect[2*i+1]*vect[2*i+1];
		}
		l2 /= N;
		double [] ni = new double [N];
		for (int i = 0; i < N; i++) {
			ni[i] = Math.sqrt((vect[2*i]*vect[2*i]+vect[2*i+1]*vect[2*i+1])/l2);
		}
		return ni;
	}
	public double [] minMaxInfluence(double [] vect) {
		double [] ni = getNormInfluence(vect);
		double mn = ni[0], mx = ni[0];
		for (int i = 1; i < N; i++) {
			if (ni[i] < mn) mn = ni[i];
			else if (ni[i] > mx) mx = ni[i];
		}
		return new double [] {mn,mx};
	}
	
	private String getType(int indx) {
		if (indx < proto_type1.length) {
			return "type1:"+indx;
		}
		indx -=  proto_type1.length;
		if (indx < proto_type2.length) {
			return "type2:"+indx;
		}
		indx -=  proto_type2.length;
		if (full_type1) {
			if (indx < proto_type1_extra.length) {
				return "type1_extra:"+indx;
			}
			indx -=  proto_type1_extra.length;
		}
		if (full_type2) {
			if (indx < proto_type2_extra.length) {
				return "type2_extra:"+indx;
			}
			indx -=  proto_type2_extra.length;
		}
		if (proto_type3.length > 0) {
			if (indx < proto_type3.length) {
				return "type3:"+(indx/N)+":"+(indx%N);
			}
			indx -=  proto_type3.length;
			
		}
		return "invalid:"+indx;
	}

	private void add_vector(int indx) {
		sym_indices[num_defined++] = indx;
		used_indices[indx] = true;
		// normalize selected vector (used when they are orthogonalized
		double l2 = 0.0;
		for (int j = 0; j < 2 *N; j++) {
			l2 += dvectors[indx][j] * dvectors[indx][j]; 
		}
		double scale = 1.0/Math.sqrt(l2);
		for (int j = 0; j < 2 *N; j++) {
			dvectors[indx][j] *= scale; 
		}
		double [] ni = getNormInfluence(dvectors[indx]);
		for (int j = 0; j < N; j++) {
			cumul_influences[j] += ni[j];
		}
	}
	
	private double[] remove_projection(double [] new_vect, double [] used_vect) { // |used_vect| === 1.0);
		double dotp = 0.0;
		for (int i = 0; i < new_vect.length; i++) {
			dotp += new_vect[i] * used_vect[i];
		}
		for (int i = 0; i < new_vect.length; i++) {
			new_vect[i] -= used_vect[i] * dotp;
		}
		return new_vect; // Math.sqrt(l2);
	}
	
	private double remove_projections(double [] new_vect) {
		for (int n = 0; n < num_defined; n++) {
			remove_projection(new_vect, dvectors[sym_indices[n]]);
		}
		return getL1(new_vect);
	}
	
	private double getL1(double [] vect) {
		double l2 = 0.0;
		for (int i = 0; i < vect.length; i++) {
			l2 += vect[i] * vect[i];
		}
		return Math.sqrt(l2);
	}

	private double [] getNewHeights() {
		double [] heights = new double [dvectors.length];
		for (int i = 0; i < heights.length; i++) {
			if (!used_indices[i]) {
//				heights[i] = remove_projections(dvectors[i].clone());
				heights[i] = remove_projections(dvectors[i]); // modify sources
			}
		}
		return heights;
	}

	private double [] getNewMins() {
		double [] mins = new double [dvectors.length];
		for (int i = 0; i < mins.length; i++) {
			if (!used_indices[i]) {
//				heights[i] = remove_projections(dvectors[i].clone());
				remove_projections(dvectors[i]); // modify sources
				mins[i] = minMaxInfluence(dvectors[i])[0];
			}
		}
		return mins;
	}
	
	
	private int bestVector(double [] primary) { // OK for RZ
		int indx = 0;
		for (int i = 1; i < primary.length; i++) {
			if (primary[i] > primary[indx]) {
				indx = i;
			}
		}
		return indx;
	}

	private int bestVector(double [] primary, boolean [] mask) { // OK for RZ
		if (mask == null) {
			return bestVector(primary);
		}
		int indx = -1;
		for (int i = 0; i < primary.length; i++) if (mask[i]){
			if ((indx < 0) || (primary[i] > primary[indx])) {
				indx = i;
			}
		}
		return indx;
	}

	private int bestVector(double delta, double [] primary, double [] secondary) { // OK for RZ
		if (secondary == null) {
			return bestVector(primary);
		}
		int ibest = bestVector(primary);
		double threshold = primary[ibest] * (1.0 - delta);
		boolean [] mask = new boolean [primary.length];
		for (int i = 0; i < primary.length; i++) {
			mask[i] = primary[i] >= threshold;
		}
		return bestVector(secondary, mask);
	}
	
	private int getNumBest(double delta_primary, double [] primary, double delta_secondary, double [] secondary) {  // OK for RZ
		if (secondary == null) {
			return getNumBest(delta_primary, primary);
		}
		int ibest = bestVector(primary);
		double threshold = primary[ibest] * (1.0 - delta_primary);
		boolean [] mask = new boolean [primary.length];
		for (int i = 0; i < primary.length; i++) {
			mask[i] = primary[i] >= threshold;
		}
		ibest = bestVector(secondary, mask);
		threshold = secondary[ibest] * (1.0 - delta_secondary);
		int num_best = 0;
		for (int i = 0; i < secondary.length; i++) {
			if (mask[i] && (secondary[i] >= threshold)) {
				num_best++;
			}
		}
		return num_best;
	}	

	private int bestVector(double delta_primary,
			               double [] primary,
			               double delta_secondary,
			               double [] secondary,
			               double[] cumul_influences,
			               boolean rz_mode) {
		if (cumul_influences == null) {
			return bestVector(delta_primary, primary, secondary);
		}
		int ibest = bestVector(primary);
		double threshold = primary[ibest] * (1.0 - delta_primary);
		boolean [] mask = new boolean [primary.length];
		for (int i = 0; i < primary.length; i++) {
			mask[i] = primary[i] >= threshold;
		}
		ibest = bestVector(secondary, mask);
		threshold = secondary[ibest] * (1.0 - delta_secondary);
		for (int i = 0; i < secondary.length; i++) {
			if (mask[i] ) {
				if (secondary[i] < threshold) {
					mask[i] = false;
				}
			}
		}
		boolean same_cumul = true;
		for (int j = 1; j < N; j++) {
			if (cumul_influences[j] != cumul_influences[0]) {
				same_cumul = false;
				break;
			}
		}
		if (!same_cumul) {
			System.out.print("");
		}
		//Balancing influences does not seem to work - all the remaining have them exactly the same
		double [] corr = new double[mask.length];
		for (int i = 0; i < mask.length; i++) if (mask[i]){
			double [] ni = rz_mode ? getNormInfluenceRZ(rz_vectors[i]): getNormInfluence(dvectors[i]);
			for (int j = 0; j < N; j++) {
				corr[i] += cumul_influences[j] * ni[j]; // best has minimal value (most negative)
			}
		}
		int best_index = -1;
		for (int i = 0; i < mask.length; i++) if (mask[i]){
			if ((best_index < 0) || (corr[i] < corr[best_index])) {
				best_index = i;
			}
		}		
		return best_index;
	}	
	
	
	
	private int getNumBest(double delta, double [] primary) {
		int ibest = bestVector(primary);
		double threshold = primary[ibest] * (1.0-delta);
		int num_best = 0;
		for (int i = 0; i < primary.length; i++) {
			if (primary[i] >= threshold) {
				num_best++;
			}
		}
		return num_best;
	}
	
	private void generateType1() {
		int numv = 1 << (N/2-1); // N/2 - first half, -1 - first choice is always 0 (inwards, clockwise)
		int [][] partial = new int [numv][N/2];
		for (int i = 0; i < numv; i++) {
			partial[i][0] = 0;
			for (int j = 1; j < N/2; j++) {
				partial[i][j] = (i >> (j - 1)) & 1;
			}
		}
		proto_type1 = new int [2*numv][N];
		for (int i = 0; i < numv; i++) {
			for (int j = 0; j < N/2; j++) {
				proto_type1[2 * i][j] =     2 * partial[i][j];  // radial in(0)/out(1)
				proto_type1[2 * i][N/2+j] = proto_type1[2 * i][j]; 
				proto_type1[2 * i + 1][j] =     2 * partial[i][j] + 1; // tangential CW(0)/CCW(1)
				proto_type1[2 * i + 1][N/2+j] = proto_type1[2 * i + 1][j]; 
			}
		}
	}

	private void generateType2() {
		int numv = 1 << (2 * (N/4 - 1) +1); // first 2 - choices, others - 4 choices
		int [][] partial = new int [numv][N/4];
		for (int i = 0; i < numv; i++) {
			partial[i][0] = (i & 1); // in or CW
			for (int j = 1; j < N/4; j++) {
				partial[i][j] = (i >> (2*(j - 1) + 1)) & 3;
			}
		}
		proto_type2 = new int [numv][N];
		for (int i = 0; i < numv; i++) {
			for (int j = 0; j < N/4; j++) {
				proto_type2[i][j + 0*N/4] =     partial[i][j];
				proto_type2[i][j + 1*N/4] =     (partial[i][j] + 3) % 4;
				proto_type2[i][j + 2*N/4] =     (partial[i][j] + 2) % 4;
				proto_type2[i][j + 3*N/4] =     (partial[i][j] + 1) % 4;
			}
		}
		if (debug_level > 2) System.out.println("proto_type2.length="+proto_type2.length);
	}
	
	private void generateType1Extra() {
		int numv = 1 << (N/2-1); // N/2 - first half, -1 - first choice is always 0 (inwards, clockwise)
		int numv_all = 1 << (2 * (N/2 - 1) +1); // first 2 - choices, others - 4 choices
		int [][] partial = new int [numv_all][N/2];
		for (int i = 0; i < numv_all; i++) {
			partial[i][0] = (i & 1); // in or CW
			for (int j = 1; j < N/2; j++) {
				partial[i][j] = (i >> (2*(j - 1) + 1)) & 3;
			}
		}
		proto_type1_extra = new int [numv_all - 2 * numv][N];
		int indx = 0;
		for (int i = 0; i < numv_all; i++) {
			int all_and = 3;
			int all_or = 0;
			for (int j = 0; j < N/2; j++) {
				all_and &= partial[i][j];
				all_or  |= partial[i][j];
			}
			if (((all_and & 1) == 1) || ((all_or & 1) == 0)) {
				continue; // all radial or tangential, should be included in proto_type1
			}
			for (int j = 0; j < N/2; j++) {
				proto_type1_extra[indx][j] =     partial[i][j];
				proto_type1_extra[indx][N/2+j] = proto_type1_extra[indx][j];
			}
			indx++;
		}
		if (debug_level > 2) System.out.println("proto_type1.length="+proto_type1.length+" proto_type1_extra.length="+proto_type1_extra.length);
	}
	
	private void generateType2Extra(boolean no_type1_extra) {
		int numv_all = 1 << (2 * (N/2 - 1) +1); // first 2 - choices, others - 4 choices (combined quarters 1 and 3,
		                                        // verify that Q3 is not -Q1
		int [][] partial = new int [numv_all][N/2]; // 0..N/4-1 -> Q1, N/4...N/2-1 ->Q3
		for (int i = 0; i < numv_all; i++) {
			partial[i][0] = (i & 1); // in or CW
			for (int j = 1; j < N/2; j++) {
				partial[i][j] = (i >> (2*(j - 1) + 1)) & 3;
			}
		}
		boolean [] new_vect = new boolean[numv_all];
		int num_new = 0;
		for (int i = 0; i < numv_all; i++) {
			boolean differs=false;
			for (int j = 0; j < N/4; j++) {
				if (((partial[i][j] - partial[i][j+N/4] + 2) %4) !=0) {
					differs = true;
					break;
				}
			}
			if (!differs) {
				continue; // Q3 = - Q1, included in proto_type2
			}
			if (no_type1_extra) {
				differs=false;
				for (int j = 0; j < N/4; j++) {
					if (((partial[i][j] - partial[i][j+N/4]) %4) !=0) {
						differs = true;
						break;
					}
				}
				if (!differs) {
					continue; // Q3 = - Q1, included in proto_type2
				}
			}
			new_vect[i]=true;
			num_new++;
		}		
		
		proto_type2_extra = new int [num_new][N];
		int indx = 0;
		for (int i = 0; i < numv_all; i++) if (new_vect[i]){
			for (int j = 0; j < N/4; j++) {
				proto_type2_extra[indx][j + 0*N/4] =     partial[i][j];
				proto_type2_extra[indx][j + 1*N/4] =     (partial[i][j] + 3) % 4;
				proto_type2_extra[indx][j + 2*N/4] =     partial[i][j+N/4];
				proto_type2_extra[indx][j + 3*N/4] =     (partial[i][j+N/4] + 3) % 4;
			}
			indx++;
		}
		if (debug_level > 2) System.out.println("proto_type2_extra.length="+proto_type2.length);
	}
	
	private void generateType3(int offset) { // in N >=8 (N%4=0)
		int ord = 0;
		for (ord = 0; (1 << ord) < (N/4); ord++);
		if (ord < 1) {
			proto_type3 = new int[0][];
			return;
		}
		proto_type3 = new int[N * ord][];
		int [] flip_tang = {0,3,2,1}; 
		for (int mode = ord; mode > 0; mode--) {
			int vect_offs = (ord - mode) * N; 
			int [] proto = new int [N];
			int [] up = new int [N]; // 1: camera moves up, 0: down  
			int period = 1 << mode;
			if (period > N/4) {
				period = N/4;
			}
			int period_high = (period+1) / 2;
			int up_index = period_high/2; 
			for (int i = 0; i < N/2; i++) {
				int i1 = (i + up_index) % (2 * period); // phase in 2 periods
				if (i1 < period) {
					up[i] = (i1 <       (period + 1) / 2) ? 1 : -1;
				} else {
					up[i] = ((i1 - period) < period / 2)?   1 :  -1;
				}

			}
			for (int i = 0; i < N/2; i++) {
				int i2 = 2*N - 2*i - (period_high / 2);
				if (up[i] < 0) {
					i2 += N; // 180 phase shift for down
				}
				i2 = (i2 + N/4) % (2*N);
				if (i2 < N/2) {
					proto[i] = R_PLUS;
				} else if (i2 < N) {
					proto[i] = T_PLUS;
				} else if (i2 < 3 * N/2) {
					proto[i] = R_MINUS;
				} else {
					proto[i] = T_MINUS;
				}
			}
			if ((period_high % 2) > 0) {
				proto[N/2] = R_MINUS;
				up[N/2] = 1;
				for (int i = 1; i < N/2; i++) {
					up[N - i] = up[i];
					proto[N - i] = flip_tang[proto[i]];
				}
			} else {
				for (int i = 0; i < N/2; i++) {
					up   [N - 1 - i] = up[i];
					proto[N -1 - i] =  flip_tang[proto[i]];
				}
			}
			proto_type3[vect_offs] = proto;
			for (int i = 1; i < N; i++) {
				proto_type3[vect_offs + i]= new int[N];
				for (int j=0; j< N; j++) {
					proto_type3[vect_offs + i][j] = proto_type3[vect_offs + 0][(N+j-i) % N]; 
				}
			}
			if (offset < 0) continue;
			// balance up/down 
			double [] scales = {0.0, 0.0 };
			double [] dbgx =   {0.0, 0.0 };
			if (debug_level > 1) {
				System.out.print(mode+"("+ord+"): ");
				for (int i = 0; i < N; i++) {
					System.out.print(proto[i]+" ");
				}			
				System.out.println();
			}
			for (int i = 0; i < N; i++) {
				double rad =  AX_DIR[proto[i]][0];
				double tang =  AX_DIR[proto[i]][1];
				if (debug_level > 1) System.out.print(String.format("(%4.1f, %4.1f) ",rad,tang));
				double alpha = (i*2 + (((period_high % 2)==0)? 1: 0))*Math.PI/N;
				double y = 	-rad*Math.cos(alpha) + tang*Math.sin(alpha);
				double x = 	rad*Math.sin(alpha) + tang*Math.cos(alpha);				
				if (debug_level > 1)System.out.print(String.format("[%5.2f, %5.2f] ",x,y));
				if (up[i] > 0) {
					scales[0] += y;
					dbgx[0] += x;
				} else {
					scales[1] += y;
					dbgx[1] += x;
				}
			}
			if (debug_level > 1) System.out.println();
			if (debug_level > 2) {
				System.out.println("dbgx[0], dbgx[1] = "+dbgx[0]+", "+dbgx[1]);
				System.out.println("dbgy[0], dbgy[1] = "+scales[0]+", "+scales[1]);
				System.out.println("sum_up/sum_down = "+(scales[0]/scales[1]));
			}
			
			scales[0] = Math.abs(Math.sqrt(-scales[0] * scales[1]) / scales[0]);  
			scales[1] = 1.0 / scales[0];
			if (debug_level > 2) {
				System.out.println("scale_down/scale_up = "+(scales[1]/scales[0]));
			}
			int mod_offs = offset + vect_offs; // N * (ord - mode);
			for (int i = 0; i <N; i++) {
				double l2 = 0.0;
				for (int j = 0; j < N; j++) {
					int j1 =(N+j-i) % N;
					double s = (up[j1] > 0) ? scales[0] : scales[1];
					dvectors[i + mod_offs][2*j] =     s * AX_DIR[proto_type3[i + vect_offs][j]][0];
					dvectors[i + mod_offs][2*j + 1] = s * AX_DIR[proto_type3[i + vect_offs][j]][1];
					l2 += dvectors[i + mod_offs][2*j] *     dvectors[i + mod_offs][2*j];
					l2 += dvectors[i + mod_offs][2*j + 1] * dvectors[i + mod_offs][2*j + 1];
					if (debug_level > 3) {
						System.out.print(String.format("%d:%4.2f(%4.1f, %4.1f) ",proto_type3[i + vect_offs][j], s,
								AX_DIR[proto_type3[i + vect_offs][j]][0], AX_DIR[proto_type3[i + vect_offs][j]][1]));
					}
				}
				if (debug_level > 3) {
					System.out.println();
				}
				double scale = 1.0/Math.sqrt(l2);
				for (int j = 0; j < 2* N; j++) {
					dvectors[i + mod_offs][j] *=     scale;
				}
				//debug
				if (debug_level > 2) {
					double [] sXY = getXY(dvectors[i + mod_offs]);
					System.out.println("mode="+mode+" ("+ord+"), i="+i+": sXY=["+sXY[0]+","+sXY[1]+"] scale ="+scale);
				}
			}
		}
		if (debug_level > 2) {
			System.out.println();
		}
	}
	class RZLine {
		int [] line;
		RZLine(int n){
			line = new int[n];
		}
	}
	
	private void generateTypeRZ() {
		ArrayList<RZLine> rz_list = new ArrayList<RZLine>();
		RZLine line = new RZLine(N);
		for (int i = 0; i < N; i++) line.line[i] = 1;
		rz_list.add(line);
		line = new RZLine(N);
		for (int i = 0; i < N; i++) line.line[i] = ((i + N/4) % N < N/2)? -1 : 1;
		rz_list.add(line);
		line = new RZLine(N);
		for (int i = 0; i < N; i++) line.line[i] = (i < N/2)?             -1 : 1;
		rz_list.add(line);
	    line = new RZLine(N);
		for (int i = 0; i < N; i++) line.line[i] = (((i + N/4) % N < N/2) ^ (i < N/2))? 1 : -1;
		rz_list.add(line);
		int ord = 0;
		for (ord = 0; (1 << ord) < (N/4); ord++);
		if (ord >= 1) {
			for (int mode = ord; mode > 0; mode--) {
				int period = 1 << mode;
				if (period > N/4) {
					period = N/4; 
				}
				boolean odd_period = (period & 1) != 0; 
				// rotate each pattern by 0..period-1
				for (int shft = 0; shft < period; shft++) {
					// for each period but first allow independent reversal
					int num_pers = (N/period);
					int num_reversals = 1 << (num_pers - 1);
					for (int rvrs_patt = 0; rvrs_patt <  num_reversals;  rvrs_patt++) {
						int [] proto = new int [N];
						for (int iper = 0; iper < num_pers; iper++) {
							boolean rvrs = ((rvrs_patt >> (num_pers - iper -1)) & 1) != 0; // for the first period - always false
							boolean odds =  odd_period && ((iper & 1) != 0); 
							for (int j = 0; j < period; j++) {
								int rj = rvrs ? (period - j -1): j;
								boolean plus = odds? (rj < (period - 1)/2) : (rj < (period + 1)/2);
								proto[(iper * period + j + shft) % N] = plus ? 1 : -1; 
							}
						}
						line = new RZLine(N);
						line.line = proto;
						rz_list.add(line);
/*
				int [] proto = new int [N/2];
				
				for (int i = 0; i < N/2; i++) {
					if (i < period) {
						proto[i] = (i <       (period + 1) / 2) ? 1 : -1;
					} else {
						proto[i] = ((i - period) < period / 2)?   1 :  -1;
					}
				}
				for (int i = 0; i < period; i++) {
					line = new RZLine(N);
					for (int j = 0; j< N/2; j++) {
						line.line[j] =     proto [(i + j ) % (N/2)];	
						line.line[j+N/2] = line.line[j];	
					}
					rz_list.add(line);
				}
*/				
					}
				}
			}
		}
		proto_rz = new int [rz_list.size()][];
		rz_vectors = new double [proto_rz.length][];
		used_rz_indices = new boolean [proto_rz.length];
		for (int i = 0; i < proto_rz.length; i++) {
			proto_rz[i] = rz_list.get(i).line;
			rz_vectors[i] = new double[proto_rz[i].length];
			// normalize and create rz_vectors
			double l2 = 0.0;
			for (int j = 0; j < proto_rz[i].length; j++) {
				l2 += proto_rz[i][j]*proto_rz[i][j];
			}
			double s = 1.0/Math.sqrt(l2);
			for (int j = 0; j < proto_rz[i].length; j++) {
				rz_vectors[i][j] = s * proto_rz[i][j];
			}
		}
	}

	private void add_vector_rz(int indx) {
		rz_indices[num_rz_defined++] = indx;
		used_rz_indices[indx] = true;
		// normalize selected vector (used when they are orthogonalized
		double l2 = 0.0;
		for (int j = 0; j < N; j++) {
			l2 += rz_vectors[indx][j] * rz_vectors[indx][j]; 
		}
		double scale = 1.0/Math.sqrt(l2);
		for (int j = 0; j < N; j++) {
			rz_vectors[indx][j] *= scale; 
		}
		double [] ni = getNormInfluenceRZ(rz_vectors[indx]);
		for (int j = 0; j < N; j++) {
			cumul_rz_influences[j] += ni[j];
		}
	}
	
	
	private double remove_projections_rz(double [] new_vect) {
		for (int n = 0; n < num_rz_defined; n++) {
			remove_projection(new_vect, rz_vectors[rz_indices[n]]);
		}
		return getL1(new_vect);
	}
	
	private double [] getNewHeightsRZ() {
		double [] heights = new double [rz_vectors.length];
		for (int i = 0; i < heights.length; i++) {
			if (!used_rz_indices[i]) {
				heights[i] = remove_projections_rz(rz_vectors[i]); // modify sources
			}
		}
		return heights;
	}

	private double [] getNewMinsRZ() {
		double [] mins = new double [rz_vectors.length];
		for (int i = 0; i < mins.length; i++) {
			if (!used_rz_indices[i]) {
//				heights[i] = remove_projections(dvectors[i].clone());
				remove_projections_rz(rz_vectors[i]); // modify sources
				mins[i] = minMaxInfluenceRZ(rz_vectors[i])[0];
			}
		}
		return mins;
	}
	
	public double [] getNormInfluenceRZ(double [] vect) {
		double l2 = 0.0;
		for (int i = 0; i < N; i++) {
			l2 += vect[i]*vect[i];
		}
		l2 /= N;
		double [] ni = new double [N];
		for (int i = 0; i < N; i++) {
			ni[i] = Math.sqrt((vect[i]*vect[i])/l2);
		}
		return ni;
	}
	public double [] minMaxInfluenceRZ(double [] vect) {
		double [] ni = getNormInfluenceRZ(vect);
		double mn = ni[0], mx = ni[0];
		for (int i = 1; i < N; i++) {
			if (ni[i] < mn) mn = ni[i];
			else if (ni[i] > mx) mx = ni[i];
		}
		return new double [] {mn,mx};
	}
	
	
	/*
	 * 
	public static int R_MINUS = 0; // radial inwards
	public static int T_PLUS =  1; // tangential clockwise
	public static int R_PLUS =  2; // radial outwards
	public static int T_MINUS = 3; // tangential counter-clockwise
	 * 
	 * Try 1) proto_type3[1] - rotated [0] to balance weights - yes, seems to work
	 * Verify Sx, Sy
	 * For 16: use same 2 proto_type3 pair, same - rotated by 1/16 and last pair - specific for 16
	 * 
	 * for 8/16 - try all rotations, not just 2 to find most ortho?
	 * 
	 * Build next decimation by alternating approximate directions. Use just 2 lengths
	 * 
	 * for 16: 2 one direction, 2 another, 2 - first, ...
	 * 
	 * Or use only radial? probably no
	 * 
	 * No need for _extra?
	 * 
	 */
	
}
class SymmVectorsSet {
//	int num_sensors;
	double [][] xy; 
	double [][] rt;
	int    [][] dir_rt;
	double [][] rots;
	double [][] zooms;
}

