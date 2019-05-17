package com.elphel.imagej.tileprocessor;
import java.util.ArrayList;

import Jama.Matrix;

/**
 **
 ** Correlation2dLMA - Fit multi - baseline correaltion pairs to the model
 **
 ** Copyright (C) 2018 Elphel, Inc.
 **
 ** -----------------------------------------------------------------------------**
 **
 **  CorrelationRigLMA.java is free software: you can redistribute it and/or modify
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


/*
 * Fitting parabola (use value_to_weight to discount off-maximum samples) for multiple cameras, such as a pair of quad cams
 * It takes 2d correlations, uses common offsets dependence for x and y offsets on disparity
 *
 * Fits the following commonn (to all pairs) parameters:
 * A - overlall strength
 * Sx, Sy, Sxy - reverse widths
 * d - disparity
 *
 * Each pair may have 3 other parameters, with a commom mask and regularization coefficients (2) to punish for not being == 0:
 * M[i] - additive term (scaled with A
 * X[i], Y[i] - additional x, y ofsfet caused my mis-alignment of the cameras
 *
 * 2-d function to be fitted for each pair is (x,y - coordinates in 2d correlation, relative to the center):
 * f(x,y) = A * (M[i] + (1 - (Sx^2*(x-x0)^2 + Sy^2*(y-y0)^2+Sxy*(x-x0)*(y-y0)))), where
 * x0 = X[i] - d * Kx[i]
 * y0 = Y[i] - d * Ky[i]
 *
 *  Kx, Ky are calculated in advance, form geometry andrelative disparity (1.0 for the main camera, ~0.6 for aux, ~5.0 for inter)
 *  Kxy={
 *  {1, 0}, // top pair
 *  {1, 0}, // bottom pair
 *  {0, 1}, // left pair
 *  {0, 1}, // right pair
 *  {1, 1}, // main diagonal pair
 *  {1,-1}} // other diagonal pair
 */

public class CorrelationRigLMA {
	static final double [][] QUAD_KXY = {
			 {1.0, 0.0}, // top pair
			 {1.0, 0.0}, // bottom pair
			 {0.0, 1.0}, // left pair
			 {0.0, 1.0}, // right pair
			 {1.0, 1.0}, // main diagonal pair
			 {1.0,-1.0}}; // other diagonal pair
	static final int INDEX_DISP =  0;
	static final int INDEX_DA =    1;
	static final int INDEX_DSX =   2;
	static final int INDEX_DSY =   3;
	static final int INDEX_DSXY =  4;
	static final int INDEX_DM =    5;
	static final int INDEX_DX0 =   6;
	static final int INDEX_DY0 =   7;
	static final int INDEX_LEN =   8;

	static final String [] PAR_NAMES = {"Disparity","Amplitude","sharpness-X","sharpness-Y","Narrowness-XY", "shift", "dx", "dy"};
	static final String [] PAIR_NAMES = {"top","bottom","left","right","diagm","diago"};
	static final String [] CAM_NAMES = {"main","aux","inter"};

	double [][] Kxy = null;
	boolean     use_master =       true;
	boolean     use_aux =          true;
	boolean     use_inter =        true;
	double      default_width =    2.0; // pix
	boolean     adjust_M =         false;
	boolean     adjust_xy =        false;
	double      cost_width =       0.0;
	double      cost_max_value =   0.0;
	double      cost_xy =          0.0;
	double      weight_inter =     4.0; // relative weight of inter-camera correlation to a pair correlation.
	double [][] corr_data =       null;
	int         num_pairs;        // number of correlation pairs
	int         num_pair_pars;    // number of per-pair individual paramers (0,1,2,3 - m, xy)

	int [][]    pair_par_index =  null; // reverse table - from parameter number to pair (-1 - global) and derivative index ([0] - not used?)
	int [][]    mxy_pair_index = null;  // from pair to parameter number for M, and DX (-1 - none)
	ArrayList<ArrayList<Sample>> samples;
	double []   y_vector; // measured values and trailing zeros for regularization (0nly m, xy, no width)
	double []   w_vector; // weights maching y
	double      weights_pure; // sum of sampes weights (1.0 - regularization weights)
	double [][] xy_vector; // x and y coordinates in the 2d correlation arrays
	int []      pair_vector; // number of a pair for sample

	double []   vector;

	// next values are only updated after success
	double []   last_rms =        null; // {rms, rms_pure}, matching this.vector
	double []   good_or_bad_rms = null; // just for diagnostics, to read last (failed) rms
	double []   initial_rms =     null; // {rms, rms_pure}, first-calcualted rms
	double []   last_ymfx =       null;
	double [][] last_jt =         null;
//	boolean     input_diag =      false; // valid during adding samples, should be set before changing groups
//	double []   poly_coeff =      null;  // 6 elements - Xc, Yx, f(x,y), A, B, C (from A*x^2 + B*y^2 +C*x*y+...)
//	double []   poly_xyvwh =      null;  // result of 2-d polynomial approximation instead of the LMA - used for lazy eye correction



	class Sample{
		double x;      // x coordinate on the common scale (corresponding to the largest baseline), along the disparity axis
		double y;      // y coordinate (0 - disparity axis)
		double v;      // correlation value at that point
		double w;
		Sample (
				double x,      // x coordinate on the common scale (corresponding to the largest baseline), along the disparity axis
				double y,      // y coordinate (0 - disparity axis)
				double v,      // correlation value at that point
				double w)     // baseline scale index
			{
				this.x =  x;
				this.y =  y;
				this.v =  v;
				this.w =  w;
			}
	}


	public void initRigKx(
			double rel_disp_aux,
			double rel_disp_inter_x,  // for aux camera to the right of the main rel_disp_inter_x>0, rel_disp_inter_y = 0
			double rel_disp_inter_y) {
		int num_cam_pairs = QUAD_KXY.length;
		Kxy = new double [2 * num_cam_pairs+1][]; // 2 cameras with 6 pairs and one inter-camera pair
		for (int i = 0; i < num_cam_pairs; i++) {
			Kxy[i                ][0] = QUAD_KXY[i][0];
			Kxy[i                ][1] = QUAD_KXY[i][1];
			Kxy[i + num_cam_pairs][0] = QUAD_KXY[i][0] * rel_disp_aux;
			Kxy[i + num_cam_pairs][1] = QUAD_KXY[i][0] * rel_disp_aux;
		}
		Kxy[2 * num_cam_pairs][0] = rel_disp_inter_x; //
		Kxy[2 * num_cam_pairs][1] = rel_disp_inter_y; // 0.0 for horizontal pair
	}

	public void setup( // prepare once for all tiles
			boolean use_master,
			boolean use_aux,
			boolean use_inter,
			boolean adjust_max_value,
			boolean adjust_xy,
			double  default_width,
			double  cost_width,
			double  cost_max_value,
			double  cost_xy,
			double  weight_inter) {
		this.use_master =       use_master;
		this.use_aux =          use_aux;
		this.use_inter =        use_inter;

		this.adjust_M =         adjust_max_value;
		this.adjust_xy =        adjust_xy;
		this.default_width =    default_width;

		this.cost_width =       cost_width;
		this.cost_max_value =   cost_max_value;
		this.cost_xy =          cost_xy;

		this.weight_inter =     weight_inter;

		setupMask();
		corr_data = new double [pair_par_index.length][];
		samples = new ArrayList<ArrayList<Sample>>();
		for (int i = 0; i < Kxy.length; i++) samples.add((Kxy[i] != null)?(new ArrayList<Sample>()):null);
	}

	public void addSample(
			int         ncam,      // 0 - main, 1 - aux, 2 - inter
			int         npair,
			double      x,      // x coordinate on the common scale (corresponding to the largest baseline), along the disparity axis
			double      y,      // y coordinate (0 - disparity axis)
			double      v,      // correlation value at that point
			double      w) {
		if (ncam == 2) 	w *= weight_inter;
		samples.get(QUAD_KXY.length * ncam + npair).add(new Sample(x,y,v,w));
	}

	public void setYW() {
		int num_samples =0;
		for (ArrayList<Sample> als : samples) {
			if (als != null) num_samples += als.size();
		}
		int num_reg = 	num_pairs * num_pair_pars;

		w_vector = new double [num_samples + num_reg];
		y_vector = new double [w_vector.length];
		xy_vector = new double [w_vector.length][];
		int ns = 0;
		double sw = 0.0;
		double w_ind = (adjust_M? cost_max_value : 0.0) + 2 * (adjust_xy? cost_xy : 0.0); // for all tiles
		weights_pure = 1.0 - w_ind;
		double wm =  cost_max_value/num_pairs;
		double wxy = cost_xy/num_pairs;
		for (int pair = 0; pair < samples.size(); pair++) {
			ArrayList<Sample> als = samples.get(pair);
			if (als != null) {
				for (Sample sample : als) {
					w_vector[ns] = sample.w;
					sw +=          sample.w;
					y_vector[ns] = sample.v;
					xy_vector[ns] = new double[2];
					xy_vector[ns][0] = sample.x;
					xy_vector[ns][1] = sample.y;
					pair_vector[ns] = pair;
					ns++;
				}
			}
		}
		for (int np = 0; np < num_pairs; np++) {
			if (adjust_M) w_vector[ns++] = wm;
			if (adjust_xy) {
				w_vector[ns++] = wxy;
				w_vector[ns++] = wxy;
			}
		}
		if (sw > 0.0) {
			double kw =  weights_pure/sw;
			for (int i = 0; i < num_samples; i++) w_vector[i] *= kw;
		}
	}

	double [] getFxAndJacobian(
			double [] vector,
			double [][] jt) { // should be initialized as double [<number of parameters>][] or null
		int num_points = w_vector.length;
		double [] fx = new double [num_points];
		double [] derivs = null;
		if (jt != null) {
			for (int i= 0; i< jt.length; i++) jt[i] = new double [num_points];
			derivs = new double[INDEX_LEN];
		}
		int num_reg = 	num_pairs * num_pair_pars;
		int ns;
		for (ns = 0; ns < num_points - num_reg; ns ++) {
			int np = pair_vector[ns];
			double M = (mxy_pair_index[np][0] >= 0)? vector[mxy_pair_index[np][0]]:Double.NaN;
			double [] Dxy = null;
			int xy_index = mxy_pair_index[np][1]; // pair number to parameter index (0 - m, 1 - x in xy)
			if (xy_index >= 0) {
				Dxy = new double[2];
				Dxy[0] = vector[xy_index];
				Dxy[1] = vector[xy_index+1];
			}
			fx[ns] = getFxAndDerivatives(
					np,    // int         pair,
					xy_vector[ns][0],   // double         x,
					xy_vector[ns][1],   // double         y,
					vector[INDEX_DISP], // double      disparity,
					vector[INDEX_DA],   // double      A,
					vector[INDEX_DSX],  // double      Sx,
					vector[INDEX_DSY],  // double      Sy,
					vector[INDEX_DSXY], // double      Sxy,
					M,                  // double      M,  // NaN - assumed 0
					Dxy,                // double []   Dxy,// null - assumed {0,0}
					derivs);            // double []   deriv);
			// save jacobian results - just for the relevant parameters
			for (int nderiv = 0; nderiv <INDEX_DM; nderiv++) {
				jt[nderiv][ns] = derivs[nderiv]; // common parameters
			}
			if (mxy_pair_index[np][0] >= 0) {
				jt[mxy_pair_index[np][0]][ns] =  derivs[INDEX_DM];
			}
			if (xy_index >= 0) {
				jt[xy_index  ][ns] =  derivs[INDEX_DX0];
				jt[xy_index+1][ns] =  derivs[INDEX_DY0];
			}
		}
		// ns points to the first regularization "sample". Put diagonal of 1.0;
		for (int npar = INDEX_DM; npar < jt.length; npar++ ) {
			fx[ns] = vector[npar];
			jt[npar][ns++] = 1.0;
		}
		return fx;
	}

	public double [] getFxAndJacobian(
			double      delta, // for testing derivatives: calculates as delta-F/delta_x
			double []   vector,
			double [][] jt) { // should be either [vector.length][] or null - then only fx is calculated
		double [] fx0=getFxAndJacobian(vector,null);
		for (int np = 0; np < vector.length; np++) {
			double [] vector1 = vector.clone();
			vector1[np]+= delta;
			double [] fxp=getFxAndJacobian(vector1,null);
			vector1 = vector.clone();
			vector1[np]-= delta;
			double [] fxm=getFxAndJacobian(vector1,null);

			jt[np] = new double [fxp.length];
			for (int i = 0; i < fxp.length; i++) {
				jt[np][i] = (fxp[i] - fxm[i])/delta/2;
			}
		}
		return fx0;
	}

    public void debugJt(
    		double      delta,
    		double []   vector) {
//    	int num_points = this.values.length;
    	int num_pars = vector.length;
    	double [] max_diff = new double [num_pars];

//    	delta = 0.001;

    	double [][] jt =       new double [num_pars][];
    	double [][] jt_delta = new double [num_pars][];
    	double [] fx = getFxAndJacobian( vector,jt);
    	getFxAndJacobian(delta, vector,jt_delta);
    	System.out.println("Test of jt-jt_delta difference,  delta = "+delta);
    	System.out.print(String.format("%3s: %10s ", "#", "fx"));
    	for (int anp = 0; anp< pair_par_index.length; anp++){
    		int pair =     pair_par_index[anp][0];
    		int par_type = pair_par_index[anp][1];
    		int ncam = -1;
    		if (pair >= 0) ncam = pair/ QUAD_KXY.length;
    		String par_name = PAR_NAMES[par_type] + ((pair >= 0) ?("-"+PAIR_NAMES[pair]+"-"+CAM_NAMES[ncam]):"");
        	System.out.print(String.format("%17s ", par_name));
    	}
    	System.out.println();

    	for (int i = 0; i < fx.length; i++) {
        	System.out.print(String.format("%3d: %10.7f ", i, fx[i]));
        	for (int np = 0; np < num_pars; np++) {
            	System.out.print(String.format("%8.5f %8.5f ", jt_delta[np][i], 1000*(jt[np][i] - jt_delta[np][i])));
            	double adiff = Math.abs(jt[np][i] - jt_delta[np][i]);
            	if (adiff > max_diff[np]) {
            		max_diff[np] = adiff;
            	}
        	}
        	System.out.println();
    	}
    	System.out.print(String.format("%15s ", "Maximal diff:"));
    	for (int np = 0; np < num_pars; np++) {
        	System.out.print(String.format("%8s %8.5f ", "1/1000×",  1000*max_diff[np]));
    	}
    	System.out.println();
    }

//	int [][]    pair_par_index =  null; // reverse table - from parameter number to pair (-1 - global) and derivative index ([0] - not used?)


	double [] getYMinusFxWeighted(
			double [] fx,
			double [] rmses // {rms,rms_pure};
			) {
		if (rmses == null) {
			rmses = new double[2];
		}
		double [] y_minus_fx_w = new double [fx.length];
		double swd2=0.0;
		int num_pure = fx.length - num_pairs * num_pair_pars;
		for (int ns = 0; ns < num_pure; ns++) {
			double d =  y_vector[ns]-fx[ns];
			double dw = d * w_vector[ns];
			swd2 += d * dw;
			y_minus_fx_w[ns] = dw;
		}
		rmses[1] = Math.sqrt(swd2)/ weights_pure;
		for (int ns = num_pure; ns < fx.length; ns++) {
			double d = -fx[ns];
			double dw = d * w_vector[ns];
			swd2 += d * dw;
			y_minus_fx_w[ns] = dw;
		}
		rmses[0] = Math.sqrt(swd2);
		return y_minus_fx_w;
	}

// Calculate JtJ, JtYminusFx

	public double [][] getWJtJlambda(
			double      lambda,
			double [][] jt){
		int num_pars = jt.length;
		int nup_points = jt[0].length;
		double [][] wjtjl = new double [num_pars][num_pars];
		for (int i = 0; i < num_pars; i++) {
			for (int j = i; j < num_pars; j++) {
				double d = 0.0;
				for (int k = 0; k < nup_points; k++) {
					d += this.w_vector[k]*jt[i][k]*jt[j][k];
				}
				wjtjl[i][j] = d;
				if (i == j) {
					wjtjl[i][j] += d * lambda;
				} else {
					wjtjl[j][i] = d;
				}
			}
		}
		return wjtjl;
	}

	public double [] getJtWdiff(
			double []   wdiff,
			double [][] jt){
		int num_pars = jt.length;
		int nup_points = jt[0].length;
		double [] wjtymfx = new double [num_pars];
		for (int i = 0; i < num_pars; i++) {
			double d = 0;
			for (int j = 0; j < nup_points; j++) d += wdiff[j] + jt[i][j];
			wjtymfx[i] = d;
		}
		return wjtymfx;
	}

	public boolean runLma(
			double lambda,           // 0.1
			double lambda_scale_good,// 0.5
			double lambda_scale_bad, // 8.0
			double lambda_max,       // 100
			double rms_diff,         // 0.001
			int    num_iter,         // 20
			int    debug_level)
	{
		boolean [] rslt = {false,false};
		int iter = 0;
		for (iter = 0; iter < num_iter; iter++) {
			rslt =  lmaStep(
					lambda,
					rms_diff,
					debug_level);
			if (debug_level > 1) {
				System.out.println("LMA step "+iter+": {"+rslt[0]+","+rslt[1]+"} full RMS="+good_or_bad_rms[0]+
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
					break;
				}
			}
		}
		if (rslt[0]) { // better, but num tries exceeded
			if (iter >= num_iter) {
				if (debug_level > 0) System.out.println("Step "+iter+": Improved, but number of steps exceeded maximal");
			} else {
				if (debug_level > 0) System.out.println("Step "+iter+": LMA: Success");
			}

		} else { // improved over initial ?
			if (last_rms[0] < initial_rms[0]) {
				rslt[0] = true;
				if (debug_level > 0) System.out.println("Step "+iter+": Failed to converge, but result improved over initial");
			} else {
				if (debug_level > 0) System.out.println("Step "+iter+": Failed to converge");
			}
		}
		if (debug_level > 0) {
			System.out.println("LMA: full RMS="+last_rms[0]+" ("+initial_rms[0]+"), pure RMS="+last_rms[1]+" ("+initial_rms[1]+") + lambda="+lambda);
		}

		return rslt[0];
	}

	public boolean [] lmaStep(
			double lambda,
			double rms_diff,
			int debug_level) {
		int num_pars = vector.length;
		boolean [] rslt = {false,false};
		if (this.last_rms == null) { //first time, need to calculate all (vector is valid)
			this.last_jt = new double [num_pars][] ; // [num_points];
			double [] fx = getFxAndJacobian(
					this.vector, // double []   vector,
					this.last_jt); // double [][] jt) { // should be either [vector.length][samples.size()] or null - then only fx is calculated
			this.last_ymfx = getYMinusFxWeighted(
					fx, // double [] fx,
					this.last_rms); // double [] rmses) // {rms,rms_pure};

			this.initial_rms = this.last_rms.clone();
			this.good_or_bad_rms = this.last_rms.clone();
			if (debug_level > 3) {
				debugJt(
						0.000001, // double      delta,
						this.vector); // double []   vector);
			}
		}
		Matrix y_minus_fx_weighted = new Matrix(this.last_ymfx, this.last_ymfx.length);

		Matrix wjtjlambda = new Matrix(getWJtJlambda(
				lambda, // *10, // temporary
				this.last_jt)); // double [][] jt)
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
				System.out.println("Singular Matrix");
			}

			return rslt;
		}

		if (debug_level>2) {
			System.out.println("(JtJ + lambda*diag(JtJ).inv()");
			jtjl_inv.print(18, 6);
		}

		Matrix jty = (new Matrix(this.last_jt)).times(y_minus_fx_weighted);
		if (debug_level>2) {
			System.out.println("Jt * (y-fx)");
			jty.print(18, 6);
		}

		Matrix mdelta = jtjl_inv.times(jty);
		if (debug_level>2) {
			System.out.println("mdelta");
			mdelta.print(18, 6);
		}


		double []  delta = mdelta.getColumnPackedCopy();
		double [] new_vector = this.vector.clone();
		for (int i = 0; i < num_pars; i++) new_vector[i]+= delta[i];
		// being optimistic, modify jt and last_ymfx in place, restore if failed
		double[] fx = getFxAndJacobian(
				new_vector, // double []   vector,
				this.last_jt); // double [][] jt) { // should be either [vector.length][samples.size()] or null - then only fx is calculated
		double [] rms = new double[2];
		this.last_ymfx = getYMinusFxWeighted(
				fx, // double [] fx,
				this.last_rms); // double [] rmses) // {rms,rms_pure};

		this.good_or_bad_rms = rms.clone();
		if (rms[0] < this.last_rms[0]) { // improved
			rslt[0] = true;
			rslt[1] = rms[0] >=(this.last_rms[0] * (1.0 - rms_diff));
			this.last_rms = rms.clone();
			this.vector = new_vector.clone();
			if (debug_level > 2) {
				System.out.print("New vector: ");
				for (int np = 0; np < vector.length; np++) {
					System.out.print(this.vector[np]+" ");
				}
				System.out.println();
			}

		} else { // worsened
			rslt[0] = false;
			rslt[1] = false; // do not know, caller will decide
			// restore state
			fx = getFxAndJacobian( //recalculate fx
					new_vector, // double []   vector,
					this.last_jt); // double [][] jt) { // should be either [vector.length][samples.size()] or null - then only fx is calculated

			this.last_ymfx =  getYMinusFxWeighted(
					fx, // double [] fx,
					this.last_rms); // double [] rmses) // {rms,rms_pure};

			if (debug_level > 2) {
				debugJt(
						0.000001, // double      delta,
						this.vector); // double []   vector);
			}
		}
		return rslt;
	}




	private void setupMask() {
		if (!use_master) for (int i = 0; i < QUAD_KXY.length; i++) Kxy[i] =                   null;
		if (!use_aux)    for (int i = 0; i < QUAD_KXY.length; i++) Kxy[i + QUAD_KXY.length] = null;
		if (!use_inter)                                            Kxy[2 * QUAD_KXY.length] = null;
		num_pairs = (use_master ? QUAD_KXY.length:0) + (use_aux ? QUAD_KXY.length:0) + (use_inter ? 1:0);
		num_pair_pars = (adjust_M? 1 : 0) + (adjust_xy ? 2:0);
		int num_pars = INDEX_DM + num_pairs * num_pair_pars;
		//double [] vector = new double [num_pars];
		pair_par_index = new int [num_pars][2];
		mxy_pair_index = new int [num_pairs][2];  // from pair to parameter number for M, and DX (-1 - none)

		int npar = 0;
		for (int i = 0; i <= INDEX_DSXY; i++) { // common parameters
			pair_par_index[npar]  [0] = -1;
			pair_par_index[npar++][1] = i;
		}
		for (int pair = 0; pair < Kxy.length; pair++) if (Kxy[pair] != null) {
			if (adjust_M) {
				mxy_pair_index[pair][0] = npar;
				pair_par_index[npar  ][0] = pair;
				pair_par_index[npar++][1] = INDEX_DM;
			} else {
				mxy_pair_index[pair][0] = -1;
			}
			if (adjust_xy) {
				mxy_pair_index[pair  ][1] = npar;
				pair_par_index[npar  ][0] = pair;
				pair_par_index[npar++][1] = INDEX_DX0; //== 1
				pair_par_index[npar  ][0] = pair;
				pair_par_index[npar++][1] = INDEX_DY0; //== 2
			} else {
				mxy_pair_index[pair][1] = -1;
			}
		}
	}

	public void addCorrData(
			int         ncam,      // 0 - main, 1 - aux, 2 - inter
			double [][] corr_data) // 6 pairs for cameras, 1 pair - for inter
	{
		for (int i = 0; i < corr_data.length; i++) this.corr_data[QUAD_KXY.length * ncam + i] = corr_data[i];
	}
	public double [] init_vector() {
		double [] vector = new double [pair_par_index.length];
		vector[INDEX_DSX] = 1.0/default_width;
		vector[INDEX_DSY] = 1.0/default_width;
		// all other parameters start with 0.0;
		return vector;
	}


	public double getFxAndDerivatives(
			int         pair,
			double      x,
			double      y,
			double      disparity,
			double      A,
			double      Sx,
			double      Sy,
			double      Sxy,
			double      M,  // null - assumed 0
			double []   Dxy,// null - assumed {0,0}
			double []   deriv) // null - do not calculate. should be initialized to INDEX_LEN length
	                           // returns derivatives for d, A, Sx, Sy, Sxy, M[pair], Dxy[pair][0], Dxy[pair][1]
	{	double [] xy0 = {-disparity* Kxy[pair][0], -disparity* Kxy[pair][1]};
		if (Dxy != null) {
			xy0[0] += Dxy[0];
			xy0[1] += Dxy[1];
		}
		double dxn = x - xy0[0];
		double dyn = y - xy0[1];
		double dx = Sx * dxn;
		double dy = Sy * dyn;
		double nfxy = (1.0 - (dx*dx + dy*dy +2*Sxy*dx*dy));
		if (!Double.isNaN(M)) nfxy += M;
		double  fxy = A * nfxy;
		if (deriv != null) {
			// df/dA
			deriv[INDEX_DA] = nfxy;
			// d/dSxy[2]
			deriv[INDEX_DSXY] = -A * dx*dy;
			// d/dother
			double d_ddx = -2 * A* (dx + Sxy*dy);
			double d_ddy = -2 * A* (dy + Sxy*dx);
			// d/ddisparity = d_ddx * ddx_ddisparity + d_ddy * ddy_ddisparity
			deriv[INDEX_DISP] = (d_ddx * Sx * Kxy[pair][0]) + ( d_ddy * Sy * Kxy[pair][1]);
			// d/dSX =  d_ddx * ddx_dSX
			deriv[INDEX_DSX] = d_ddx * dxn;
			// d/dSY =  d_ddy * ddy_dSY
			deriv[INDEX_DSY] = d_ddy * dyn;
			// individual/maskable
			if (!Double.isNaN(M)) deriv[INDEX_DM] = A;
			if (Dxy != null) {
				// d/dDX =  d_ddx * ddx_dDX
				deriv[INDEX_DX0] = -d_ddx * Sx;
				deriv[INDEX_DY0] = -d_ddy * Sy;
			}
		}
		return fxy;
	}
/*
	static final int INDEX_DISP =  0;
	static final int INDEX_DA =    1;
	static final int INDEX_DSX =   2;
	static final int INDEX_DSY =   3;
	static final int INDEX_DSXY =  4;
	static final int INDEX_DM =    5;
	static final int INDEX_DX0 =   6;
	static final int INDEX_DY0 =   7;
	static final int INDEX_LEN =   8;

 */

}



