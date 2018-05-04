import java.util.ArrayList;
import java.util.HashMap;

import Jama.Matrix;

/**
 **
 ** Correlation2dLMA - Fit multi - baseline correaltion pairs to the model
 **
 ** Copyright (C) 2018 Elphel, Inc.
 **
 ** -----------------------------------------------------------------------------**
 **
 **  Correlation2dLMA.java is free software: you can redistribute it and/or modify
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
 * Fitting parabola for multiple grids
 * Difference between ortho and diagonals - just point coordinates and extra overall scale (weight account for number averaged)
 * Each group of compatible is already averaged, so each group has a single individual variable - scale.
 *
 * Parabolas (1 for each group) are Ag * (1 - ((x-x0)/Wx) ^ 2 - (y/Wy)^2), where As is a per-group scale
 * Wy = Wm * scale +Wyd
 * Wx = Wm * scale +Wyd + Wxy
 *
 * Wm is a correlation measurement parameter, it does not depend on x/y and on particular pair, it depends on the LPF, so the
 * total contribution is proportional to the baseline reduction (scale)
 *
 * Wyd is widening caused the image itself - probably noise and other factors of poor correlation contrast. When multiple
 * orthogonal directions are combined it influences equally all directions (x,y) so Wx includes that term also
 *
 * Wxy widens maximum in disparity direction, it is caused by multiple overlapping maximums for different disparities and for
 * strong enough matches can indicate miz of disparities in the same tile
 *
 * Fitting of a single scale groups (1 or 2) has to have Wm constant.
 *
 *

 *
 */

public class Correlations2dLMA {
	final static int X0_INDEX =  0;
	final static int WM_INDEX =  1; // may be frozen
	final static int WYD_INDEX = 2;
	final static int WXY_INDEX = 3;
	final static int AG_INDEX =  4; // and above
	final static String [] PAR_NAMES = {"X0","HALF-WIDTH","WIDTH-EXTRA","WIDTH_X/Y","AMPLITUDE"};
	double  [] all_pars;
	boolean [] par_mask;
	double  [] vector;
	double  [] scales = {1.0, 2.0, 4.0};
	ArrayList<Sample> samples = new ArrayList<Sample>();
	HashMap<Integer,NumDiag> groups = new HashMap<Integer,NumDiag>();
	double  [] group_weights =   null; // per-group weights of samples sum == 1.0
	double [] weights; // normalized so sum is 1.0 for all - samples and extra regularization terms
	double    pure_weight; // weight of samples only
	double [] values;
	// next values are only updated after success
	double []   last_rms =        null; // {rms, rms_pure}, matching this.vector
	double []   good_or_bad_rms = null; // just for diagnostics, to read last (failed) rms
	double []   initial_rms =     null; // {rms, rms_pure}, first-calcualted rms
	double []   last_ymfx =       null;
	double [][] last_jt =         null;
	boolean     input_diag =      false; // valid during adding samples, should be set before changing groups
	double []   poly_coeff =      null;  // 6 elements - Xc, Yx, f(x,y), A, B, C (from A*x^2 + B*y^2 +C*x*y+...)
	double []   poly_xyvwh =      null;  // result of 2-d polynomial approximation instead of the LMA - used for lazy eye correction

	public class NumDiag{
		int     num;
		boolean diag;
		public NumDiag(int num, boolean diag) {
			this.num =  num;
			this.diag = diag;
		}
	}

	public class Sample{
		double x;      // x coordinate on the common scale (corresponding to the largest baseline), along the disparity axis
		double y;      // y coordinate (0 - disparity axis)
		double v;      // correlation value at that point
		double w;
		int    si;     // baseline scale index
		int    gi;     // group index
		Sample (
				double x,      // x coordinate on the common scale (corresponding to the largest baseline), along the disparity axis
				double y,      // y coordinate (0 - disparity axis)
				double v,      // correlation value at that point
				double w,
				int    si,     // baseline scale index
				int    gi)     // baseline scale index
			{
				this.x =  x;
				this.y =  y;
				this.v =  v;
				this.w =  w;
				this.si = si;
				this.gi = gi;
			}
	}

	public void printParams() {
		for (int np = 0; np < all_pars.length; np++) {
			System.out.println(String.format("%2d%1s %22s %f",
					np,
					(par_mask[np]?"+":" "),
					((np < AG_INDEX)?PAR_NAMES[np]:(PAR_NAMES[AG_INDEX]+(np-AG_INDEX))),
					all_pars[np]));
		}
	}

	public void printInputDataFx(boolean show_fx){
		if 	(show_fx) {
			Sample s = null;
			double [] fx = getPolyFx();
			if (fx == null) fx = getFx();
			if (fx == null) return;
			for (int i = 0; i < fx.length; i++) {
				double fx_pos = (fx[i] >= 0)? fx[i]: Double.NaN;
				if (i < samples.size()) {
					s = samples.get(i);
					System.out.println(String.format("%3d: x=%8.4f y=%8.4f v=%9.6f fx=%9.6f w=%9.7f si=%1d gi=%1d", i, s.x, s.y, s.v, fx_pos, s.w, s.si, s.gi ));
				}
				else {
					System.out.println(String.format("%3d: %10s %10s v=%9.6f fx=%9.6f w=%9.7f", i, "---", "---", this.values[i], fx_pos, this.weights[i]));
				}
			}
		} else {
			int ns =0;
			for (Sample s:samples){
				System.out.println(String.format("%3d: x=%8.4f y=%8.4f v=%9.6f w=%9.7f si=%1d gi=%1d", ns++, s.x, s.y, s.v, s.w, s.si, s.gi ));
			}
		}
	}

	public double [] getRMS() {
		return last_rms;
	}
	public double [] getGoodOrBadRMS() {
		return good_or_bad_rms;
	}

	public double [] getAllPars() {
		return all_pars;
	}

	public double [] getDisparityStrength() {
		if (group_weights == null) return null;
		double disparity = -all_pars[X0_INDEX];
		double sum_amp = 0.0;
		for (int i = 0; i < (all_pars.length - AG_INDEX); i++) {
			sum_amp += group_weights[i]* all_pars[AG_INDEX + i]; // group_weights is normalized
		}
		// protect from weird fitting results
		double max_amp = 0.0;
		for (Sample s: samples) if (s.v > max_amp) max_amp = s.v;
		if (sum_amp > 1.25 * max_amp) sum_amp = max_amp;
		double [] ds = {disparity, sum_amp};
		return ds;
	}
	public double [] getDisparityStrengthWidth() {
		double [] ds = getDisparityStrength();
		if (ds == null) return null;
		double [] dsw = {ds[0], ds[1], all_pars[WM_INDEX], all_pars[WXY_INDEX]}; // asymmetry
		return dsw;
	}


	public Correlations2dLMA (
			double [] scales // null - use default table
			) {
		if (scales != null)	this.scales = scales.clone();
	}

	public void setDiag(boolean diag_in) {
		this.input_diag = diag_in;
	}
	public void addSample(
			double x,      // x coordinate on the common scale (corresponding to the largest baseline), along the disparity axis
			double y,      // y coordinate (0 - disparity axis)
			double v,      // correlation value at that point
			double w,
			int    si,     // baseline scale index
			int    gi){     // baseline scale index
		samples.add(new Sample(x,y,v,w,si,gi));
		if (!groups.containsKey(gi)) {
			groups.put(gi,new NumDiag(groups.size(),this.input_diag));
		}

	}
//NumDiag
// TODO: add auto x0, half-width?
// should be called ater all samples are entered (to list groups)
	public void initVector(
			boolean adjust_wm,
			boolean adjust_wy,
			boolean adjust_wxy,
			boolean adjust_Ag,
			double  x0,
			double  half_width,
			double  cost_wm,     // cost of non-zero this.all_pars[WM_INDEX]
			double  cost_wxy     // cost of non-zero this.all_pars[WXY_INDEX]
			) {
//		for (Sample s:samples) if (!groups.containsKey(s.gi)) groups.put(s.gi,groups.size());
		int num_groups = groups.size();
		this.all_pars =        new double[AG_INDEX + num_groups];
		this.all_pars[X0_INDEX] =  x0;
		this.all_pars[WM_INDEX] =  half_width;
		this.all_pars[WYD_INDEX] = 0.0;
		this.all_pars[WXY_INDEX] = 0.0;
		this.par_mask = new boolean[AG_INDEX + num_groups];
		this.par_mask[X0_INDEX] = true;
		this.par_mask[WM_INDEX] = adjust_wm;
		this.par_mask[WYD_INDEX] = adjust_wy;
		this.par_mask[WXY_INDEX] = adjust_wxy;
		for (int i = 0; i <num_groups; i++) {
			this.par_mask[AG_INDEX + i] = adjust_Ag;
		}
		setWeightsValues(half_width, cost_wm, cost_wxy);
		int imx = 0;
		int num_samples = samples.size();
		for (int i = 1; i<num_samples; i++) if (values[i] > values[imx]) imx = i;
		for (int i = 0; i < num_groups; i++) this.all_pars[AG_INDEX + i] = values[imx];
		toVector();
	}

	public void setWeightsValues(
			double half_width,   // expected width
			double  cost_wm,     // cost of non-zero this.all_pars[WYD_INDEX]
			double  cost_wxy) {  // cost of non-zero this.all_pars[WXY_INDEX]
		int np = samples.size();
		group_weights = new double[groups.size()];
		weights = new double [np+2];
		values =  new double [np+2];
		weights[np] =   cost_wm;
		weights[np+1] = cost_wxy;
		values[np] =    half_width;
		values[np+1] =  0.0;
		double sw = 0;
		for (int i = 0; i < np; i++) {
			Sample s = samples.get(i);
			weights[i] = s.w;
			values[i] =  s.v;
			group_weights[groups.get(s.gi).num] += s.w;
			if (Double.isNaN(values[i]) || Double.isNaN(weights[i])) {
				weights[i] = 0.0;
				values[i] = 0.0;
			}
			sw += weights[i];
		}
		pure_weight = sw;
		sw += weights[np] + weights[np+1];
		if (sw != 0.0) {
			sw = 1.0/sw;
			for (int i = 0; i < weights.length; i++) weights[i] *= sw;
		}

		if (pure_weight > 0.0) for (int i = 0; i < group_weights.length; i++) group_weights[i] /= pure_weight;
		pure_weight *= sw;
	}

	public void toVector() {
		int np = 0;
		for (int i = 0; i < par_mask.length; i++) if (par_mask[i]) np++;
		vector = new double[np];
		np = 0;
		for (int i = 0; i < par_mask.length; i++) if (par_mask[i]) vector[np++] = all_pars[i];
	}

	public void updateFromVector() {
		int np = 0;
		for (int i = 0; i < par_mask.length; i++) if (par_mask[i]) all_pars[i] = vector[np++];
	}

	public double [] fromVector(double [] vector) { // mix fixed and variable parameters
		if ( all_pars == null) return null;
		double [] ap = all_pars.clone();
		int np = 0;
		for (int i = 0; i < par_mask.length; i++) if (par_mask[i]) ap[i] = vector[np++];
		return ap;
	}

    public void debugJt(
    		double      delta,
    		double []   vector) {
    	int num_points = this.values.length;
    	int num_pars = vector.length;
    	double [] max_diff = new double [num_pars];

//    	delta = 0.001;

    	double [][] jt =       new double [num_pars][num_points];
    	double [][] jt_delta = new double [num_pars][num_points];
    	double [] fx = getFxJt( vector,jt);
    	getFxJt(delta, vector,jt_delta);
    	System.out.println("Test of jt-jt_delta difference,  delta = "+delta);
    	System.out.print(String.format("%3s: %10s ", "#", "fx"));
    	for (int anp = 0; anp< all_pars.length; anp++) if(par_mask[anp]){
        	System.out.print(String.format("%17s ", ((anp < AG_INDEX)?PAR_NAMES[anp]:(PAR_NAMES[AG_INDEX]+(anp-AG_INDEX)))));
    	}
    	System.out.println();

    	for (int i = 0; i < num_points; i++) {
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


	public double [] getFxJt(
			double      delta, // for testing derivatives: calculates as delta-F/delta_x
			double []   vector,
			double [][] jt) { // should be either [vector.length][samples.size()] or null - then only fx is calculated
		double [] fx0=getFxJt(vector,null);
		for (int np = 0; np < vector.length; np++) {
			double [] vector1 = vector.clone();
			vector1[np]+= delta;
			double [] fxp=getFxJt(vector1,null);
			vector1 = vector.clone();
			vector1[np]-= delta;
			double [] fxm=getFxJt(vector1,null);

			jt[np] = new double [fxp.length];
			for (int i = 0; i < fxp.length; i++) {
				jt[np][i] = (fxp[i] - fxm[i])/delta/2;
			}
		}
		return fx0;
	}



	public double [] getFx() {
		return getFxJt(this.vector, null);
	}

	public double [] getFxJt(
			double []   vector,
			double [][] jt) { // should be either [vector.length][samples.size()] or null - then only fx is calculated
		if (vector == null) return null;
		double [] av = fromVector(vector);
		int num_samples = samples.size();
		double [] fx= new double [num_samples + 2];
		int num_groups = groups.size();
		double sqrt2 = Math.sqrt(2.0);
		for (int ns = 0; ns < num_samples; ns++) {
			Sample s = samples.get(ns);
			int grp = groups.get(s.gi).num;
			boolean diag = groups.get(s.gi).diag;
			double wScale = diag?sqrt2:1.0;
			double Wy = av[WM_INDEX] * scales[s.si] + av[WYD_INDEX];
			double Wx = Wy + av[WXY_INDEX];
			double dx = s.x - av[X0_INDEX];
			double Ag = av[AG_INDEX+grp];
			double dxw = wScale*dx/Wx;
			double dyw = wScale*s.y/Wy;
			double d = (1.0 - dxw*dxw - dyw*dyw);
			fx[ns] = d * Ag;
			int np = 0;
			if (jt != null) {
//				if (par_mask[X0_INDEX])  jt[np++][ns] = Ag * 2 * dxw/Wx; // d/dx0
				if (par_mask[X0_INDEX])  jt[np++][ns] = wScale * Ag * 2 * dxw/Wx; // d/dx0
				double dfdWx = Ag * 2 * dxw * dxw / Wx;
				double dfdWy = Ag * 2 * dyw * dyw / Wy;
				if (par_mask[WM_INDEX])  jt[np++][ns] =  scales[s.si] * (dfdWx + dfdWy); // d/dWm
				if (par_mask[WYD_INDEX]) jt[np++][ns] = (dfdWx + dfdWy); // d/dWyd
				if (par_mask[WXY_INDEX]) jt[np++][ns] =  dfdWx;          // d/dWxy
				for (int i = 0; i < num_groups; i++) {
					if (par_mask[AG_INDEX + i]) jt[np++][ns] =  (i == grp) ? d : 0.0;    // d/dWAg
				}
			}
		}
		// add 2 extra samples - damping cost_wy, cost_wxy
		fx[num_samples] =     av[WM_INDEX];
		fx[num_samples + 1] = av[WXY_INDEX];
		// and derivatives
		if (jt != null) {
			int np = 0;
			for (int i = 0; i < par_mask.length; i++) if (par_mask[i]) {
				jt[np][num_samples] =     (i == WM_INDEX)? 1.0 : 0.0;
				jt[np++][num_samples+1] = (i == WXY_INDEX)? 1.0 : 0.0;
			}
		}
		return fx;
	}

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
					d += this.weights[k]*jt[i][k]*jt[j][k];
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

	// returns {rms, rms_pure}
	public double [] getWYmFxRms(
			double []   fx) { // will be replaced with y-fx
		int num_samples = samples.size();
		int num_points = fx.length; // includes 2 extra for regularization
		double rms = 0, rms_pure = 0;
		for (int i = 0; i < num_samples; i++) {
			double d = (values[i] - fx[i]);
			fx[i] = this.weights[i] * d;
			rms += fx[i]*d; // sum of weights
		}
		rms_pure = Math.sqrt(rms)/this.pure_weight;
		for (int i = num_samples; i < num_points; i++) {
			double d = (values[i] - fx[i]);
			fx[i] = this.weights[i] * d;
			rms += fx[i]*d; // sum of weights
		}
		rms = Math.sqrt(rms);
		double [] rslt = {rms, rms_pure};
		return rslt;
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
//			System.out.println("LMA: full RMS="+good_or_bad_rms[0]+" ("+initial_rms[0]+"), pure RMS="+good_or_bad_rms[1]+" ("+initial_rms[1]+") + lambda="+lambda);
		}

		return rslt[0];
	}

/*
	double []   last_rms =     null; // {rms, rms_pure}, matching this.vector
	double []   initial_rms =  null; // {rms, rms_pure}, first-calcualted rms

 */

	// returns {success, done}
	public boolean [] lmaStep(
			double lambda,
			double rms_diff,
			int debug_level) {
		int num_points = this.weights.length; // includes 2 extra for regularization
		int num_pars = vector.length;
		boolean [] rslt = {false,false};
		if (this.last_rms == null) { //first time, need to calculate all (vector is valid)
			this.last_jt = new double [num_pars][num_points];
			this.last_ymfx = getFxJt(
					this.vector, // double []   vector,
					this.last_jt); // double [][] jt) { // should be either [vector.length][samples.size()] or null - then only fx is calculated
			this.last_rms = getWYmFxRms(this.last_ymfx); // modifies this.last_ymfx
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
/*
				try {
					this.SYNC_COMMAND.wait();
				} catch (InterruptedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}

 */
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
		this.last_ymfx = getFxJt(
				new_vector, // double []   vector,
				this.last_jt); // double [][] jt) { // should be either [vector.length][samples.size()] or null - then only fx is calculated
		double [] rms = getWYmFxRms(this.last_ymfx); // modifies this.last_ymfx
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
			this.last_ymfx = getFxJt( // recalculate fx
					this.vector, // double []   vector,
					this.last_jt); // double [][] jt) { // should be either [vector.length][samples.size()] or null - then only fx is calculated
			this.last_rms = getWYmFxRms(this.last_ymfx); // modifies this.last_ymfx
			if (debug_level > 2) {
				debugJt(
						0.000001, // double      delta,
						this.vector); // double []   vector);
			}


		}
		return rslt;
	}

// modify to reuse Samples and apply polynomial approximation to resolve x0,y0 and strength?
	public double [] getMaxXYPoly( // get interpolated maximum coordinates using 2-nd degree polynomial
///			double  outside, // how much solution may be outside of the samples
			boolean debug
     ) {
		double [][][] mdata = new double[samples.size()][3][];
///		double min_x = Double.NaN, max_x = Double.NaN, min_y = Double.NaN, max_y =Double.NaN , min_v = Double.NaN;
		for (int i = 0; i < mdata.length; i++) {
			Sample s = samples.get(i);
			mdata[i][0] = new double [2];
			mdata[i][0][0] =  s.x;
			mdata[i][0][1] =  s.y;
			mdata[i][1] = new double [1];
			mdata[i][1][0] =  s.v;
			mdata[i][2] = new double [1];
			mdata[i][2][0] =  s.w;
			/*
			if (i == 0) {
				min_x = s.x;
				max_x = min_x;
				min_y = s.y;
				max_y = min_y;
				min_v = s.v;
			} else {
				if      (s.x > max_x) max_x = s.x;
				else if (s.x < min_x) min_x = s.x;
				if      (s.y > max_y) max_y = s.y;
				else if (s.y < min_y) min_y = s.y;
				if      (s.v < min_v) min_v = s.x;
			}
			*/
		}
		double [] rslt = (new PolynomialApproximation()).quadraticMaxV2dX2Y2XY( // 9 elements - Xc, Yx, f(x,y), A, B, C, D, E, F (from A*x^2 + B*y^2 +C*x*y+...)
				mdata,
				1.0E-30,//25, // 1.0E-15,
				debug? 4:0);
		this.poly_coeff = rslt;

		if ((rslt == null) || (rslt[2] < 0.0) || // negative strength
				(rslt[3] >= 0.0) || // x: min, not max
				(rslt[4] >= 0.0)) { // y: min, not max
			this.poly_coeff = null;
			this.poly_xyvwh = null;
			return null;
		}
//		if ()

		// calculate width_x and width_y
		double hwx = Double.NaN, hwy = Double.NaN;
		if ((rslt[2] > 0.0) && (rslt[3] <0.0) && (rslt[4] <0.0)) {
			hwx = Math.sqrt(-rslt[2]/rslt[3]);
			hwy = Math.sqrt(-rslt[2]/rslt[4]);
		}
		double [] xyvwh = {rslt[0], rslt[1], rslt[2], hwx, hwy};
		if (debug){
			System.out.println("lma.getMaxXYPoly()");
			for (int i = 0; i< mdata.length; i++){
				System.out.println(i+": "+mdata[i][0][0]+"/"+mdata[i][0][1]+" z="+mdata[i][1][0]+" w="+mdata[i][2][0]);
			}
			System.out.println("quadraticMax2d(mdata) --> "+((rslt==null)?"null":(rslt[0]+"/"+rslt[1])));
		}
		this.poly_xyvwh = xyvwh;
		return xyvwh; // rslt;
	}
	public double [] getPoly() {
		return poly_xyvwh;
	}
	public double [] getPolyFx() {return getPolyFx(this.poly_coeff);}
	public double [] getPolyFx(
			double [] coeff) { // 6 elements - Xc, Yx, f(x,y), A, B, C (from A*x^2 + B*y^2 +C*x*y+...)
		if (coeff == null) {
			return null;
		}
		int num_samples = samples.size();
		double [] fx= new double [num_samples];
		for (int ns = 0; ns < num_samples; ns++) {
			Sample s = samples.get(ns);
			fx[ns]= coeff[3]*s.x*s.x + coeff[4]*s.y*s.y +  coeff[5]*s.x*s.y + coeff[6]*s.x + coeff[7]*s.y + + coeff[8];
		}
		return fx;
	}


//	public double [] getValues(xyvwh) {
//		double [] values= new double [samples.size()];
//		for (int i = 0; i < values.length; i++) values[i] = samples.get(i).v;
//		return values;
//	}

	/*
	 *
		if (debugLevel>-1) {
			jtj.print(18, 6);
		}
		Matrix jtj_inv = jtj.inverse();
		Matrix jty = jt.times(y_minus_fx_weighted);
		Matrix mrslt = jtj_inv.times(jty);
		double []  drslt = mrslt.getColumnPackedCopy();
	 *
	 * Fitting parabola for multiple grids
	 * Difference between ortho and diagonals - just point coordinates and extra overall scale (weight account for number averaged)
	 * Each group of compatible is already averaged, so each group has a single individual variable - scale.
	 *
	 * Parabolas (1 for each group) are Ag * (1 - ((x-x0)/Wx) ^ 2 - (y/Wy)^2), where As is a per-group scale
	 * Wy = Wm * scale +Wyd
	 * Wx = Wm * scale +Wyd + Wxy
	 *
	 * Wm is a correlation measurement parameter, it does not depend on x/y and on particular pair, it depends on the LPF, so the
	 * total contribution is proportional to the baseline reduction (scale)
	 *
	 * Wyd is widening caused the image itself - probably noise and other factors of poor correlation contrast. When multiple
	 * orthogonal directions are combined it influences equally all directions (x,y) so Wx includes that term also
	 *
	 * Wxy widens maximum in disparity direction, it is caused by multiple overlapping maximums for different disparities and for
	 * strong enough matches can indicate miz of disparities in the same tile
	 *
	 * Fitting of a single scale groups (1 or 2) has to have Wm constant.
	 *
	 *

	 *
	 */




}
