package com.elphel.imagej.tileprocessor;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.ObjectOutputStream;
import java.io.Serializable;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
import java.util.ArrayList;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.atomic.DoubleAdder;

import javax.xml.bind.DatatypeConverter;

import Jama.Matrix;

public class IntersceneLma {
//	OpticalFlow opticalFlow = null;
	QuadCLT [] scenesCLT =    null; // now will use just 2 - 0 -reference scene, 1 - scene.  
	private double []         last_rms =        null; // {rms, rms_pure}, matching this.vector
	private double []         good_or_bad_rms = null; // just for diagnostics, to read last (failed) rms
	private double []         initial_rms =     null; // {rms, rms_pure}, first-calcualted rms
	private double []         y_vector =        null; // sum of fx(initial parameters) and correlation offsets
//	private double [][]       vector_XYS =      null;  // optical flow X,Y, confidence obtained from the correlate2DIterate()
	private double []         last_ymfx =       null;
	private double [][]       last_jt =         null;
	private double []         weights; // normalized so sum is 1.0 for all - samples and extra regularization 
	private double            pure_weight; // weight of samples only
	private boolean []        par_mask =        null;
	private int []            par_indices =     null;
	private double  []        parameters_full = null; // full parameters vector for the currenl LMA run
	private double  []        backup_parameters_full = null; // full parameters vector before first LMA run
	private double  []        parameters_vector = null; 
//	private double  []        parameters_initial = null;  // (will be used to pull for regularization)
	private double  [][]      macrotile_centers = null;  // (will be used to pull for regularization)
	private double            infinity_disparity = 0.1;  // treat lower as infinity
	private int               num_samples = 0;
	private boolean           thread_invariant = true; // Do not use DoubleAdder, provide results not dependent on threads
	public IntersceneLma(
//			OpticalFlow opticalFlow,
			boolean thread_invariant
			) {
		this.thread_invariant = thread_invariant;
//		this.opticalFlow = opticalFlow;
	}
	
	public double [][]       getLastJT(){
		return last_jt;
	}

	
	public double[] getLastRms() {
		return last_rms;
	}
	public double [] getSceneXYZ(boolean initial) {
		double [] full_vector = initial? backup_parameters_full: getFullVector(parameters_vector);
		return new double[] {
				full_vector[ErsCorrection.DP_DSX],full_vector[ErsCorrection.DP_DSY],full_vector[ErsCorrection.DP_DSZ]};
	}
	public double [] getSceneATR(boolean initial) {
		double [] full_vector = initial? backup_parameters_full: getFullVector(parameters_vector);
		return new double[] {
				full_vector[ErsCorrection.DP_DSAZ],full_vector[ErsCorrection.DP_DSTL],full_vector[ErsCorrection.DP_DSRL]};
	}
	public double [] getReferenceXYZ(boolean initial) {
		double [] full_vector = initial? backup_parameters_full: getFullVector(parameters_vector);
		return new double[] {
				full_vector[ErsCorrection.DP_DX],full_vector[ErsCorrection.DP_DY],full_vector[ErsCorrection.DP_DZ]};
	}
	public double [] getReferenceATR(boolean initial) {
		double [] full_vector = initial? backup_parameters_full: getFullVector(parameters_vector);
		return new double[] {
				full_vector[ErsCorrection.DP_DAZ],full_vector[ErsCorrection.DP_DTL],full_vector[ErsCorrection.DP_DRL]};
	}

	public double [] getSceneERSXYZ(boolean initial) { // never used
		double [] full_vector = initial? backup_parameters_full: getFullVector(parameters_vector);
		return new double[] {
				full_vector[ErsCorrection.DP_DSVX],full_vector[ErsCorrection.DP_DSVY],full_vector[ErsCorrection.DP_DSVZ]};
	}
	public double [] getSceneERSATR(boolean initial) { // never used
		double [] full_vector = initial? backup_parameters_full: getFullVector(parameters_vector);
		return new double[] {
				full_vector[ErsCorrection.DP_DSVAZ],full_vector[ErsCorrection.DP_DSVTL],full_vector[ErsCorrection.DP_DSVRL]};
	}
	public double [] getReferenceERSXYZ(boolean initial) { // never used
		double [] full_vector = initial? backup_parameters_full: getFullVector(parameters_vector);
		return new double[] {
				full_vector[ErsCorrection.DP_DSVX],full_vector[ErsCorrection.DP_DSVY],full_vector[ErsCorrection.DP_DSVZ]};
	}
	public double [] getReferenceERSATR(boolean initial) { // never used
		double [] full_vector = initial? backup_parameters_full: getFullVector(parameters_vector);
		return new double[] {
				full_vector[ErsCorrection.DP_DSVAZ],full_vector[ErsCorrection.DP_DSVTL],full_vector[ErsCorrection.DP_DSVRL]};
	}

	public String [] printOldNew(boolean allvectors) {
		return printOldNew(allvectors, 9, 6);
	}

	public String [] printOldNew(boolean allvectors, int w, int d) {
		String fmt1 = String.format("%%%d.%df", w+2,d+2); // more for the differences
		ArrayList<String> lines = new ArrayList<String>();
		for (int n = ErsCorrection.DP_DVAZ; n < ErsCorrection.DP_NUM_PARS; n+=3) {
			boolean adj = false;
			for (int i = 0; i <3; i++) adj |= par_mask[n+i];
			if (allvectors || adj) {
				String line = printNameV3(n, false, w,d)+"  (was "+printNameV3(n, true, w,d)+")";
				line += ", diff_last="+String.format(fmt1, getV3Diff(n)[0]);
				line += ", diff_first="+String.format(fmt1, getV3Diff(n)[1]);
				lines.add(line);
			}
		}
		return lines.toArray(new String[lines.size()]);
	}
	
	/**
	 * Calculate L2 for 3-vector difference between the adjusted values, this LMA start values
	 * and "backup" parameters - snapshot before the first adjustment (first calculation of 
	 * the motion vectors+LMA iteration  
	 * @param indx offset to the first element of the 3-vector in the full parameters array
	 * @return a pair of L2 differences {diff_with_this_LMA_start, diff_with_first_lma_start}
	 */
	public double [] getV3Diff(int indx) {
		double [] v_new = new double[3], v_backup = new double[3], v_start = new double[3];
		System.arraycopy(getFullVector(parameters_vector), indx, v_new, 0, 3);
		System.arraycopy(backup_parameters_full, indx, v_backup, 0, 3);
		System.arraycopy(parameters_full, indx, v_start, 0, 3);
		double l2_backup = 0;
		double l2_start = 0;
		for (int i = 0; i < 3; i++) {
			l2_backup += (v_new[i]-v_backup[i]) * (v_new[i]-v_backup[i]); 
			l2_start +=  (v_new[i]-v_start[i]) *  (v_new[i]-v_start[i]); 
		}
		return new double [] {Math.sqrt(l2_start), Math.sqrt(l2_backup)};
	}
	
	public String printNameV3(int indx, boolean initial, int w, int d) {
		double [] full_vector = initial? backup_parameters_full: getFullVector(parameters_vector);
		double [] vector = new double[3];
		for (int i = 0; i <3; i++) {
			vector[i] = full_vector[indx + i];
		}
		String name = ErsCorrection.DP_VECTORS_NAMES[indx];
		return printNameV3(name, vector, w, d);
	}
	
	public static String printNameV3(String name, double[] vector) {
		return printNameV3(name, vector, 10, 6);
	}
	
	public static String printNameV3(String name, double[] vector, int w, int d) {
		return String.format("%14s: %s", name, printV3(vector, w, d));
	}

	
	public static String printV3(double[] vector) {
		return printV3(vector, 8, 5);
	}
	public static String printV3(double[] vector, int w, int d) {
		String fmt = String.format("[%%%d.%df, %%%d.%df, %%%d.%df]", w,d,w,d,w,d);
		return String.format(fmt, vector[0], vector[1], vector[2]);
	}
	
	
	public void prepareLMA(
			final double []   scene_xyz0,     // camera center in world coordinates (or null to use instance)
			final double []   scene_atr0,     // camera orientation relative to world frame (or null to use instance)
			// reference atr, xyz are considered 0.0
			final QuadCLT     scene_QuadClt,
			final QuadCLT     reference_QuadClt,
			final boolean[]   param_select,
			final double []   param_regweights,
			final double [][] vector_XYS, // optical flow X,Y, confidence obtained from the correlate2DIterate()
			final double [][] centers,    // macrotile centers (in pixels and average disparities
			boolean           first_run,
			final int         debug_level)
	{
		scenesCLT = new QuadCLT [] {reference_QuadClt, scene_QuadClt};
		par_mask = param_select;
		macrotile_centers = centers;
		num_samples = 2 * centers.length;
		ErsCorrection ers_ref =   reference_QuadClt.getErsCorrection();
		ErsCorrection ers_scene = scene_QuadClt.getErsCorrection();
		final double []   scene_xyz = (scene_xyz0 != null) ? scene_xyz0 : ers_scene.camera_xyz;
		final double []   scene_atr = (scene_atr0 != null) ? scene_atr0 : ers_scene.camera_atr;
		final double []   reference_xyz = new double[3];
		final double []   reference_atr = new double[3];
		double [] full_parameters_vector = new double [] {
				0.0,                             0.0,                             0.0,
				ers_ref.ers_watr_center_dt[0],   ers_ref.ers_watr_center_dt[1],   ers_ref.ers_watr_center_dt[2],
				ers_ref.ers_wxyz_center_dt[0],   ers_ref.ers_wxyz_center_dt[1],   ers_ref.ers_wxyz_center_dt[2],
				reference_atr[0],                reference_atr[1],                reference_atr[2],
				reference_xyz[0],                reference_xyz[1],                reference_xyz[2],
				ers_scene.ers_watr_center_dt[0], ers_scene.ers_watr_center_dt[1], ers_scene.ers_watr_center_dt[2],
				ers_scene.ers_wxyz_center_dt[0], ers_scene.ers_wxyz_center_dt[1], ers_scene.ers_wxyz_center_dt[2],
				scene_atr[0],                    scene_atr[1],                    scene_atr[2],
				scene_xyz[0],                    scene_xyz[1],                    scene_xyz[2]};
		parameters_full = full_parameters_vector.clone();
		if ((vector_XYS != null) && (first_run || (backup_parameters_full == null))) {
			backup_parameters_full = full_parameters_vector.clone();
		}
		int num_pars = 0;
		for (int i = 0; i < par_mask.length; i++) if (par_mask[i]) num_pars++;
		par_indices = new int [num_pars];
		num_pars = 0;
		for (int i = 0; i < par_mask.length; i++) if (par_mask[i]) par_indices[num_pars++] = i;
		parameters_vector = new double [par_indices.length];
		for (int i = 0; i < par_indices.length; i++) parameters_vector[i] = full_parameters_vector[par_indices[i]];

		if (vector_XYS != null) {// skip when used for the motion blur vectors, not LMA
			setSamplesWeights(vector_XYS);  // not regularized yet !
		} else {
			weights = null; // new double[2 * centers.length];
		}

		last_jt = new double [parameters_vector.length][];
		if (debug_level > 1) {
			System.out.println("prepareLMA() 1");
		}
		double [] fx = getFxDerivs(
				parameters_vector, // double []         vector,
				last_jt,           // final double [][] jt, // should be null or initialized with [vector.length][]
				scenesCLT[1],      // final QuadCLT     scene_QuadClt,
				scenesCLT[0],      // final QuadCLT     reference_QuadClt,
				debug_level);      // final int         debug_level)
		if (vector_XYS == null) {
			return; // for MB vectors (noLMA)
		}
		
		double [][] wjtj = getWJtJlambda( // USED in lwir all NAN
				0.0,               // final double      lambda,
				last_jt);               // final double [][] jt) all 0???
		for (int i = 0; i < parameters_vector.length; i++) {
			int indx = num_samples + i;
			weights[indx] = param_regweights[par_indices[i]]/Math.sqrt(wjtj[i][i]);
		}
		normalizeWeights(); // make full weight == 1.0; pure_weight <= 1.0;
		// remeasure fx - now with regularization terms.

		if (debug_level > 1) {
			System.out.println("prepareLMA() 2");
		}
		fx = getFxDerivs(
				parameters_vector, // double []         vector,
				last_jt,                // final double [][] jt, // should be null or initialized with [vector.length][]
				scene_QuadClt,     // final QuadCLT     scene_QuadClt,
				reference_QuadClt, // final QuadCLT     reference_QuadClt,
				debug_level);      // final int         debug_level)
		y_vector = fx.clone();
		for (int i = 0; i < vector_XYS.length; i++) {
			if (vector_XYS[i] != null){
				y_vector[2 * i + 0] += vector_XYS[i][0];
				y_vector[2 * i + 1] += vector_XYS[i][1];
			}
		}
		last_rms = new double [2];
		last_ymfx = getYminusFxWeighted(
//				vector_XYS, // final double [][] vector_XYS,
				fx, // final double []   fx,
				last_rms); // final double []   rms_fp // null or [2]
		initial_rms = last_rms.clone();
		good_or_bad_rms = this.last_rms.clone();
	}
	
	public int runLma( // <0 - failed, >=0 iteration number (1 - immediately)
			double lambda,           // 0.1
			double lambda_scale_good,// 0.5
			double lambda_scale_bad, // 8.0
			double lambda_max,       // 100
			double rms_diff,         // 0.001
			int    num_iter,         // 20
			boolean last_run,
			int    debug_level)
	{
		boolean [] rslt = {false,false};
		this.last_rms = null; // remove?
		int iter = 0;
		for (iter = 0; iter < num_iter; iter++) {
			rslt =  lmaStep(
					lambda,
					rms_diff,
					debug_level);
			if (rslt == null) {
				return -1; // false; // need to check
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
				if (debug_level > 1) System.out.println("Step "+iter+": Improved, but number of steps exceeded maximal");
			} else {
				if (debug_level > 1) System.out.println("Step "+iter+": LMA: Success");
			}

		} else { // improved over initial ?
			if (last_rms[0] < initial_rms[0]) { // NaN
				rslt[0] = true;
				if (debug_level > 1) System.out.println("Step "+iter+": Failed to converge, but result improved over initial");
			} else {
				if (debug_level > 1) System.out.println("Step "+iter+": Failed to converge");
			}
		}
		boolean show_intermediate = true;
		if (show_intermediate && (debug_level > 0)) {
			System.out.println("LMA: full RMS="+last_rms[0]+" ("+initial_rms[0]+"), pure RMS="+last_rms[1]+" ("+initial_rms[1]+") + lambda="+lambda);
		}
		if (debug_level > 2){ 
			String [] lines1 = printOldNew(false); // boolean allvectors)
			System.out.println("iteration="+iter);
			for (String line : lines1) {
				System.out.println(line);
			}
		}
		if (debug_level > 0) {
			if ((debug_level > 1) ||  last_run) { // (iter == 1) || last_run) {
				if (!show_intermediate) {
					System.out.println("LMA: iter="+iter+",   full RMS="+last_rms[0]+" ("+initial_rms[0]+"), pure RMS="+last_rms[1]+" ("+initial_rms[1]+") + lambda="+lambda);
				}
				String [] lines = printOldNew(false); // boolean allvectors)
				for (String line : lines) {
					System.out.println(line);
				}
			}
		}
		if ((debug_level > -2) && !rslt[0]) { // failed
			if ((debug_level > 1) || (iter == 1) || last_run) {
				System.out.println("LMA failed on iteration = "+iter);
				String [] lines = printOldNew(true); // boolean allvectors)
				for (String line : lines) {
					System.out.println(line);
				}
			}
			System.out.println();
		}

		return rslt[0]? iter : -1;
	}
	
	private boolean [] lmaStep(
			double lambda,
			double rms_diff,
			int debug_level) {
		boolean [] rslt = {false,false};
		// maybe the following if() branch is not needed - already done in prepareLMA !
		if (this.last_rms == null) { //first time, need to calculate all (vector is valid)
			last_rms = new double[2];
			if (debug_level > 1) {
				System.out.println("lmaStep(): first step");
			}
			double [] fx = getFxDerivs(
					parameters_vector, // double []         vector,
					last_jt,           // final double [][] jt, // should be null or initialized with [vector.length][]
					scenesCLT[1],      // final QuadCLT     scene_QuadClt,
					scenesCLT[0],      // final QuadCLT     reference_QuadClt,
					debug_level);      // final int         debug_level)
			last_ymfx = getYminusFxWeighted(
//					vector_XYS, // final double [][] vector_XYS,
					fx, // final double []   fx,
					last_rms); // final double []   rms_fp // null or [2]
			this.initial_rms = this.last_rms.clone();
			this.good_or_bad_rms = this.last_rms.clone();

			if (debug_level > -1) { // temporary
				/*
				dbgYminusFxWeight(
						this.last_ymfx,
						this.weights,
						"Initial_y-fX_after_moving_objects");
                */
			}
			if (last_ymfx == null) {
				return null; // need to re-init/restart LMA
			}
			// TODO: Restore/implement
			if (debug_level > 3) {
				/*
				 dbgJacobians(
							corr_vector, // GeometryCorrection.CorrVector corr_vector,
							1E-5, // double delta,
							true); //boolean graphic)
				*/
			}
		}
		Matrix y_minus_fx_weighted = new Matrix(this.last_ymfx, this.last_ymfx.length);

		Matrix wjtjlambda = new Matrix(getWJtJlambda(
				lambda, // *10, // temporary
				this.last_jt)); // double [][] jt)
		if (debug_level > 3) {
			try {
				System.out.println("getFxDerivs(): getChecksum(this.y_vector)="+      getChecksum(this.y_vector));
				System.out.println("getFxDerivs(): getChecksum(this.weights)="+       getChecksum(this.weights));
				System.out.println("getFxDerivs(): getChecksum(this.last_ymfx)="+     getChecksum(this.last_ymfx));
				System.out.println("getFxDerivs(): getChecksum(y_minus_fx_weighted)="+getChecksum(y_minus_fx_weighted));
				System.out.println("getFxDerivs(): getChecksum(wjtjlambda)=         "+getChecksum(wjtjlambda));
			} catch (NoSuchAlgorithmException | IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		
		if (debug_level>2) {
			System.out.println("JtJ + lambda*diag(JtJ");
			wjtjlambda.print(18, 6);
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
		if (debug_level > 2) {
			try {
				System.out.println("getFxDerivs(): getChecksum(jtjl_inv)="+getChecksum(jtjl_inv));
				System.out.println("getFxDerivs(): getChecksum(jty)=     "+getChecksum(jty));
			} catch (NoSuchAlgorithmException | IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}

		
		
		Matrix mdelta = jtjl_inv.times(jty);
		if (debug_level>2) {
			System.out.println("mdelta");
			mdelta.print(18, 6);
		}

		double scale = 1.0;
		double []  delta =      mdelta.getColumnPackedCopy();
		double []  new_vector = parameters_vector.clone();
		for (int i = 0; i < parameters_vector.length; i++) {
			new_vector[i] += scale * delta[i];
		}
		if (debug_level > 2) {
			try {
				System.out.println("getFxDerivs(): getChecksum(mdelta)=            "+getChecksum(mdelta));
				System.out.println("getFxDerivs(): getChecksum(delta)=             "+getChecksum(delta));
				System.out.println("getFxDerivs(): getChecksum(parameters_vector)= "+getChecksum(parameters_vector));
				System.out.println("getFxDerivs(): getChecksum(new_vector)=        "+getChecksum(new_vector));
			} catch (NoSuchAlgorithmException | IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		
		
		double [] fx = getFxDerivs(
				new_vector, // double []         vector,
				last_jt,           // final double [][] jt, // should be null or initialized with [vector.length][]
				scenesCLT[1],      // final QuadCLT     scene_QuadClt,
				scenesCLT[0],      // final QuadCLT     reference_QuadClt,
				debug_level);      // final int         debug_level)
		double [] rms = new double[2];
		last_ymfx = getYminusFxWeighted(
//				vector_XYS, // final double [][] vector_XYS,
				fx, // final double []   fx,
				rms); // final double []   rms_fp // null or [2]
		if (debug_level > 2) {
			/*
			dbgYminusFx(this.last_ymfx, "next y-fX");
			dbgXY(new_vector, "XY-correction");
			*/
		}

		if (last_ymfx == null) {
			return null; // need to re-init/restart LMA
		}

		this.good_or_bad_rms = rms.clone();
		if (rms[0] < this.last_rms[0]) { // improved
			rslt[0] = true;
			rslt[1] = rms[0] >=(this.last_rms[0] * (1.0 - rms_diff));
			this.last_rms = rms.clone();

			this.parameters_vector = new_vector.clone();
			if (debug_level > 2) {
				// print vectors in some format
				/*
				System.out.print("delta: "+corr_delta.toString()+"\n");
				System.out.print("New vector: "+new_vector.toString()+"\n");
				System.out.println();
				*/
			}
		} else { // worsened
			rslt[0] = false;
			rslt[1] = false; // do not know, caller will decide
			// restore state
			fx = getFxDerivs(
					parameters_vector, // double []         vector,
					last_jt,           // final double [][] jt, // should be null or initialized with [vector.length][]
					scenesCLT[1],      // final QuadCLT     scene_QuadClt,
					scenesCLT[0],      // final QuadCLT     reference_QuadClt,
					debug_level);      // final int         debug_level)
			last_ymfx = getYminusFxWeighted(
//					vector_XYS, // final double [][] vector_XYS,
					fx, // final double []   fx,
					this.last_rms); // final double []   rms_fp // null or [2]
			if (last_ymfx == null) {
				return null; // need to re-init/restart LMA
			}
			if (debug_level > 2) {
				/*
				 dbgJacobians(
							corr_vector, // GeometryCorrection.CorrVector corr_vector,
							1E-5, // double delta,
							true); //boolean graphic)
							*/
			}
		}
		return rslt;
	}
	
	
	
	private void setSamplesWeights(
			final double [][] vector_XYS) // not regularizedn yet
	{
		this.weights = new double [num_samples + parameters_vector.length];
		
		final Thread[] threads = ImageDtt.newThreadArray(QuadCLT.THREADS_MAX);
		final AtomicInteger ai = new AtomicInteger(0);
		double sum_weights;
		if (thread_invariant) {
			final double [] sw_arr = new double [vector_XYS.length]; 
			for (int ithread = 0; ithread < threads.length; ithread++) {
				threads[ithread] = new Thread() {
					public void run() {
						for (int iMTile = ai.getAndIncrement(); iMTile < vector_XYS.length; iMTile = ai.getAndIncrement()) if (vector_XYS[iMTile] != null){
							double w = vector_XYS[iMTile][2];
							weights[2 * iMTile] = w;
//							asum_weight.add(w);
							sw_arr[iMTile] = w;
						}
					}
				};
			}		      
			ImageDtt.startAndJoin(threads);
			sum_weights = 0.0;
			for (double w:sw_arr) {
				sum_weights += w;
			}
		} else {
			final DoubleAdder asum_weight = new DoubleAdder(); 
			for (int ithread = 0; ithread < threads.length; ithread++) {
				threads[ithread] = new Thread() {
					public void run() {
						for (int iMTile = ai.getAndIncrement(); iMTile < vector_XYS.length; iMTile = ai.getAndIncrement()) if (vector_XYS[iMTile] != null){
							double w = vector_XYS[iMTile][2];
							weights[2 * iMTile] = w;
							asum_weight.add(w);
						}
					}
				};
			}		      
			ImageDtt.startAndJoin(threads);
			sum_weights = asum_weight.sum();
		}
		if (sum_weights <= 1E-8) {
			System.out.println("!!!!!! setSamplesWeights(): sum_weights=="+sum_weights+" <= 1E-8");
		}
		ai.set(0);
//		final double s = 0.5/asum_weight.sum();
		final double s = 0.5/sum_weights;
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				public void run() {
					for (int iMTile = ai.getAndIncrement(); iMTile < vector_XYS.length; iMTile = ai.getAndIncrement()) if (vector_XYS[iMTile] != null){
						weights[2 * iMTile] *= s;
						weights[2 * iMTile + 1] = weights[2 * iMTile];
					}
				}
			};
		}		      
		ImageDtt.startAndJoin(threads);
		pure_weight = 1.0;
	}

	/*
	@Deprecated
	private void normalizeWeights_old()
	{
//num_samples
		final Thread[] threads = ImageDtt.newThreadArray(opticalFlow.threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		final DoubleAdder asum_weight = new DoubleAdder(); 
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				public void run() {
					for (int i = ai.getAndIncrement(); i < num_samples; i = ai.getAndIncrement()){
						asum_weight.add(weights[i]);
					}
				}
			};
		}		      
		ImageDtt.startAndJoin(threads);
		final double s_pure = asum_weight.sum();
		for (int i = 0; i < par_indices.length; i++) {
			int indx =  num_samples + i;
			asum_weight.add(weights[indx]);
		}
		double full_weight = asum_weight.sum();
		pure_weight = s_pure/full_weight;
		final double s = 1.0/asum_weight.sum();
		if (Double.isNaN(s)) {
			System.out.println("normalizeWeights(): s == NaN : 1");
		}
		ai.set(0);
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				public void run() {
					for (int i = ai.getAndIncrement(); i < weights.length; i = ai.getAndIncrement()){
						weights[i] *= s;
					}
				}
			};
		}		      
		ImageDtt.startAndJoin(threads);
	}
	*/

	private void normalizeWeights()
	{
		final Thread[] threads = ImageDtt.newThreadArray(QuadCLT.THREADS_MAX);
		final AtomicInteger ai = new AtomicInteger(0);
		double full_weight, sum_weight_pure;
		if (thread_invariant) {
			sum_weight_pure = 00;
			for (int i = 0;  i < num_samples; i++) {
				sum_weight_pure += weights[i];
			}
			full_weight = sum_weight_pure;
			for (int i = 0; i < par_indices.length; i++) {
				int indx =  num_samples + i;
				full_weight += weights[indx];
			}
		} else {
			final DoubleAdder asum_weight = new DoubleAdder(); 
			for (int ithread = 0; ithread < threads.length; ithread++) {
				threads[ithread] = new Thread() {
					public void run() {
						for (int i = ai.getAndIncrement(); i < num_samples; i = ai.getAndIncrement()){
							asum_weight.add(weights[i]);
						}
					}
				};
			}		      
			ImageDtt.startAndJoin(threads);
			sum_weight_pure = asum_weight.sum();
			for (int i = 0; i < par_indices.length; i++) {
				int indx =  num_samples + i;
				asum_weight.add(weights[indx]);
			}
			full_weight = asum_weight.sum();
		}
		pure_weight = sum_weight_pure/full_weight;
		final double s = 1.0/full_weight;
		if (Double.isNaN(s)) {
			System.out.println("normalizeWeights(): s == NaN  : 2");
		}
		ai.set(0);
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				public void run() {
					for (int i = ai.getAndIncrement(); i < weights.length; i = ai.getAndIncrement()){
						weights[i] *= s;
					}
				}
			};
		}		      
		ImageDtt.startAndJoin(threads);
	}
	
	
	/**
	 * Combine unmodified parameters with these ones, using parameters mask
	 * @param vector variable-length adjustable vector, corresponding to a parameter mask
	 * @return full parameters vector combining adjusted and fixed parameters
	 */
	private double [] getFullVector(double [] vector) {
		double [] full_vector = parameters_full.clone();
		for (int i = 0; i < par_indices.length; i++) {
			full_vector[par_indices[i]] = vector[i];
		}
		return full_vector;
	}
	
	private double [] getFxDerivs(
			double []         vector,
			final double [][] jt, // should be null or initialized with [vector.length][]
			final QuadCLT     scene_QuadClt,
			final QuadCLT     reference_QuadClt,
			final int         debug_level)
	{
		final double [] full_vector = getFullVector(vector);
		final ErsCorrection ers_ref =   reference_QuadClt.getErsCorrection();
		final ErsCorrection ers_scene = scene_QuadClt.getErsCorrection();
		final double [] scene_xyz =     new double[3];;
		final double [] scene_atr  =    new double[3];
		final double [] reference_xyz = new double[3]; // will stay 0
		final double [] reference_atr = new double[3]; // will stay 0
		final boolean mb_mode = (weights == null);
		final int weights_length = mb_mode ?  (2 * macrotile_centers.length) : weights.length; 
		final double [] fx =         mb_mode ? null : (new double [weights_length]); // weights.length]; : weights.length :
		if (jt != null) {
			for (int i = 0; i < jt.length; i++) {
				jt[i] = new double [weights_length]; // weights.length];
			}
		}
		
		for (int i = 0; i <3; i++) {
			ers_ref.ers_watr_center_dt[i] =   full_vector[ErsCorrection.DP_DVAZ +  i];
			ers_ref.ers_wxyz_center_dt[i] =   full_vector[ErsCorrection.DP_DVX +   i];
			ers_scene.ers_watr_center_dt[i] = full_vector[ErsCorrection.DP_DSVAZ + i];
			ers_scene.ers_wxyz_center_dt[i] = full_vector[ErsCorrection.DP_DSVX +  i];
			reference_atr[i] = full_vector[ErsCorrection.DP_DAZ +  i];
			reference_xyz[i] = full_vector[ErsCorrection.DP_DX +  i];
			scene_atr[i] =     full_vector[ErsCorrection.DP_DSAZ +  i];
			scene_xyz[i] =     full_vector[ErsCorrection.DP_DSX +  i];
		}
		ers_scene.setupERS();
		ers_ref.setupERS();
		final Matrix [] reference_matrices_inverse = ErsCorrection.getInterRotDeriveMatrices(
				reference_atr, // double [] atr);
				true);         // boolean invert));

		final Matrix [] scene_matrices_inverse = ErsCorrection.getInterRotDeriveMatrices(
				scene_atr, // double [] atr);
				true);         // boolean invert));

		final Matrix  scene_rot_matrix = ErsCorrection.getInterRotDeriveMatrices(
				scene_atr, // double [] atr);
				false)[0]; // boolean invert));
		
		final Thread[] threads = ImageDtt.newThreadArray(QuadCLT.THREADS_MAX);
		final AtomicInteger ai = new AtomicInteger(0);
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				public void run() {
					for (int iMTile = ai.getAndIncrement(); iMTile < macrotile_centers.length; iMTile = ai.getAndIncrement()) {
						if ((macrotile_centers[iMTile]!=null) && (mb_mode || (weights[2*iMTile] > 0.0))){  // was: weights[iMTile]?
							//infinity_disparity
							boolean is_infinity = macrotile_centers[iMTile][2] < infinity_disparity;
							double [][] deriv_params = ers_ref.getDPxSceneDParameters(
									ers_scene,                   // ErsCorrection    ers_scene,
									true,                        // boolean correctDistortions,
									is_infinity,                 // boolean is_infinity,
									macrotile_centers[iMTile],   // double [] pXpYD_reference,
									reference_xyz,               // double [] reference_xyz, // reference (this) camera center in world coordinates
									scene_xyz,                   // double [] scene_xyz,    // (other, ers_scene) camera center in world coordinates
									reference_matrices_inverse,  // Matrix [] reference_matrices_inverse,
									scene_matrices_inverse,      // Matrix [] scene_matrices_inverse,
									scene_rot_matrix,            // Matrix    scene_rot_matrix,        // single rotation (direct) matrix for the scene
									debug_level);                // int debug_level);
							if (deriv_params!= null) {
								boolean bad_tile = false;
								if (!bad_tile) {
									if (!mb_mode) {
										fx[2 * iMTile + 0] = deriv_params[0][0]; // pX
										fx[2 * iMTile + 1] = deriv_params[0][1]; // pY
									}
									if (jt != null) {
										for (int i = 0; i < par_indices.length; i++) {
											int indx = par_indices[i] + 1;
											jt[i][2 * iMTile + 0] = deriv_params[indx][0]; // pX
											jt[i][2 * iMTile + 1] = deriv_params[indx][1]; // pY (disparity is not used)
										}
									}
								}
							} else if (mb_mode) {
								for (int i = 0; i < par_indices.length; i++) {
									jt[i][2 * iMTile + 0] = Double.NaN; // pX
									jt[i][2 * iMTile + 1] = Double.NaN; // ; // pY (disparity is not used)
								}
							}
						}
					}
				}
			};
		}		      
		ImageDtt.startAndJoin(threads);
		if (mb_mode) {
			return null;
		}
		// pull to the initial parameter values
		for (int i = 0; i < par_indices.length; i++) {
			fx [i + 2 * macrotile_centers.length] = vector[i]; // - parameters_initial[i]; // scale will be combined with weights
			jt[i][i + 2 * macrotile_centers.length] = 1.0; // scale will be combined with weights
		}
		if (debug_level > 3) {
			try {
				System.out.println    ("getFxDerivs(): getChecksum(fx)="+getChecksum(fx));
				if (jt != null) {
					System.out.println("getFxDerivs(): getChecksum(jt)="+getChecksum(jt));
				}
			} catch (NoSuchAlgorithmException | IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		return fx;
	}
	
	private double [][] getWJtJlambda( // USED in lwir
			final double      lambda,
			final double [][] jt)
	{
		final int num_pars = jt.length;
		final int num_pars2 = num_pars * num_pars;
		final int nup_points = jt[0].length;
		final double [][] wjtjl = new double [num_pars][num_pars];
		final Thread[] threads = ImageDtt.newThreadArray(QuadCLT.THREADS_MAX);
		final AtomicInteger ai = new AtomicInteger(0);
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				public void run() {
					for (int indx = ai.getAndIncrement(); indx < num_pars2; indx = ai.getAndIncrement()) {
						int i = indx / num_pars;
						int j = indx % num_pars;
						if (j >= i) {
							double d = 0.0;
							for (int k = 0; k < nup_points; k++) {
								if (jt[i][k] != 0) {
									d+=0;
								}
								d += weights[k]*jt[i][k]*jt[j][k];
							}
							wjtjl[i][j] = d;
							if (i == j) {
								wjtjl[i][j] += d * lambda;
							} else {
								wjtjl[j][i] = d;
							}
						}
					}
				}
			};
		}		      
		ImageDtt.startAndJoin(threads);
		return wjtjl;
	}
	
	private double [] getYminusFxWeighted(
//			final double [][] vector_XYS,
			final double []   fx,
			final double []   rms_fp // null or [2]
			) {
		final Thread[]      threads =     ImageDtt.newThreadArray(QuadCLT.THREADS_MAX);
		final AtomicInteger ai =          new AtomicInteger(0);
		final double []     wymfw =       new double [fx.length];
		double s_rms; 
		if (thread_invariant) {
			final double [] l2_arr = new double [num_samples]; 
			for (int ithread = 0; ithread < threads.length; ithread++) {
				threads[ithread] = new Thread() {
					public void run() {
						for (int i = ai.getAndIncrement(); i < num_samples; i = ai.getAndIncrement())  {
							double d = y_vector[i] - fx[i];
							double wd = d * weights[i];
							//double l2 = d * wd;
							l2_arr[i] = d * wd;
							wymfw[i] = wd;
						}
					}
				};
			}		      
			ImageDtt.startAndJoin(threads);
			s_rms = 0.0;
			for (double l2:l2_arr) {
				s_rms += l2;
			}
		} else {
			final DoubleAdder   asum_weight = new DoubleAdder();
			for (int ithread = 0; ithread < threads.length; ithread++) {
				threads[ithread] = new Thread() {
					public void run() {
						for (int i = ai.getAndIncrement(); i < num_samples; i = ai.getAndIncrement())  {
							double d = y_vector[i] - fx[i];
							double wd = d * weights[i];
							double l2 = d * wd;
							wymfw[i] = wd;
							asum_weight.add(l2);
						}
					}
				};
			}		      
			ImageDtt.startAndJoin(threads);
			s_rms = asum_weight.sum();
		}
		double rms_pure = Math.sqrt(s_rms/pure_weight);
		for (int i = 0; i < par_indices.length; i++) {
			int indx = i + num_samples;
			double d = y_vector[indx] - fx[indx]; // fx[indx] == vector[i]
			double wd = d * weights[indx];
			s_rms += d * wd;
			wymfw[indx] =   wd;
		}
		double rms = Math.sqrt(s_rms); // assuming sum_weights == 1.0; /pure_weight); shey should be re-normalized after adding regularization
		if (rms_fp != null) {
			rms_fp[0] = rms;
			rms_fp[1] = rms_pure;
		}
		return wymfw;
	}
	
	public static String getChecksum(Serializable object) throws IOException, NoSuchAlgorithmException {
	    ByteArrayOutputStream baos = null;
	    ObjectOutputStream oos = null;
	    try {
	        baos = new ByteArrayOutputStream();
	        oos = new ObjectOutputStream(baos);
	        oos.writeObject(object);
	        MessageDigest md = MessageDigest.getInstance("MD5");
	        byte[] thedigest = md.digest(baos.toByteArray());
	        return DatatypeConverter.printHexBinary(thedigest);
	    } finally {
	        oos.close();
	        baos.close();
	    }
	}
	
	
/*	
 * par_indices
		double [] rslt = {rms, rms_pure};

 * 
vector_XYS * 
	private double [] getFx(
			GeometryCorrection.CorrVector corr_vector)
	{
		int clusters = clustersX * clustersY;
		Matrix []   corr_rots =  corr_vector.getRotMatrices(); // get array of per-sensor rotation matrices
		Matrix [][] deriv_rots = corr_vector.getRotDeriveMatrices();
		double [] imu = corr_vector.getIMU(); // i)
		double [] y_minus_fx = new double  [clusters * POINTS_SAMPLE];
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
				//		    System.arraycopy(src_pixels, 0, dst_pixels, 0, src_pixels.length);  for the borders closer to 1/2 kernel size
				ddnd[0] = ddnd[0]; // ?
				for (int i = 0; i < NUM_SENSORS; i++) {
					ddnd[i + 1] = ddnd[i + 1];
					ddnd[i + 5] = ddnd[i + 5];
				}
				System.arraycopy(ddnd, 0, y_minus_fx, cluster*POINTS_SAMPLE, POINTS_SAMPLE);
			}
		}
		return y_minus_fx;
	}

	private double [][] getJacobianTransposed(
			GeometryCorrection.CorrVector corr_vector,
			double delta){
		int clusters = clustersX * clustersY;
		int num_pars = getNumPars();
		double [][] jt = new double  [num_pars][clusters * POINTS_SAMPLE ];
		double rdelta = 1.0/delta;
		for (int par = 0; par < num_pars; par++) {
			double [] pars = new double[num_pars];
			pars[par] =  delta;
			GeometryCorrection.CorrVector corr_delta = geometryCorrection.getCorrVector(pars, par_mask);
			GeometryCorrection.CorrVector corr_vectorp = corr_vector.clone();
			GeometryCorrection.CorrVector corr_vectorm = corr_vector.clone();
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
{}
*/	
	
}
