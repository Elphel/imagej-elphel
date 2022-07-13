package com.elphel.imagej.tileprocessor;
/**
 **
 ** Corr2dLMA - Fit multi - baseline correlation pairs to the model
 **
 ** Copyright (C) 2018 Elphel, Inc.
 **
 ** -----------------------------------------------------------------------------**
 **
 **  Corr2dLMA.java is free software: you can redistribute it and/or modify
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
 * strong enough matches can indicate mix of disparities in the same tile
 *
 * Fitting of a single scale groups (1 or 2) has to have Wm constant.
 *
 *
Each camera has the same scale and aligned parallel to all others
Each camera has a 2x2 matrix where d_disp - increment along disparity, d_ndisp - along rotated disparity:
Mi=  | dx/d_disp , dx/d_ndisp|
     | dy/d_disp , dy_d_ndisp|
For each pair matrices are subtracted (6 for 4 cameras)
Each correlation pair has the same shape, just shifted (add mask WXY that fades?)
ndisp - array of 4 values or all 0
disp - array of 4 symmetrized values, where the first is common disparity, 3 other may be forced to 0

A0 - common gain, A(i>0) = A0*Ai
W(x,y) window function for 2D correlation
x0j, y0j - individual for each pair, calculated from disp_i, ndisp_i, and Mi1-Mi2 (i1 - first channel in pair, i2 - second)

Fx=A(i)*W(x,y)*(pa*(x-x0j)^2+2*pb*(x-x0j)*(y-y0j)+pc*(y-y0j)^2)


 *
 */
import java.util.ArrayList;
import java.util.Arrays;

import com.elphel.imagej.common.DoubleGaussianBlur;
import com.elphel.imagej.common.PolynomialApproximation;
import com.elphel.imagej.common.ShowDoubleFloatArrays;

import Jama.Matrix;



public class Corr2dLMA {

///	final static int          NUM_CAMS =     4; // not all have to be used, so it is maximal number of cameras
///	final static int          NUM_PAIRS =    NUM_CAMS* (NUM_CAMS -1)/2; // number of possible pairs
	final static int          NTILE0 = 0;
	final static int          DISP_INDEX =   0; // common/average disparity (negative!)
	final static int          A_INDEX =      1; // A*(x-x0)^2
	final static int          B_INDEX =      2; // 2*B*(x-x0)*(y-y0)
	final static int          CMA_INDEX =    3; // C*(y-y0)^2, encode C-A
	final static int          G0_INDEX =     4; // scale of correlation pair,
///	final static int          TILE_PARAMS =  G0_INDEX +    NUM_PAIRS; // number of tile-individual parameters
	final static String []    PAR_NAMES = {"DISP","A","B","C-A"};
	final static String       PAR_NAME_SCALE = "SCALE";
	final static String       PAR_NAME_CORRDISP = "CORR-DISP";
	final static String       PAR_NAME_CORRNDISP = "CORR-NDISP";

	final Correlation2d       correlation2d;
	final int                 num_cams;
	final int                 num_pairs;
	final int                 tile_params;
	double  []                all_pars;
	private boolean []        par_mask;
	private int     []        par_map;
	private double  []        vector;
	private ArrayList<Sample> samples = new ArrayList<Sample>();
	private double            total_weight; // total weight of all samples
	private int               total_tiles;
	private double []         weights; // normalized so sum is 1.0 for all - samples and extra regularization terms
	private double            pure_weight; // weight of samples only as a fraction of total weight
	private double []         values;

	// next values are only updated after success
	private int iter =        0; // to be able to read number of used iterations 

	private double []         last_rms =        null; // {rms, rms_pure}, matching this.vector
	private double []         good_or_bad_rms = null; // just for diagnostics, to read last (failed) rms
	private double []         initial_rms =     null; // {rms, rms_pure}, first-calculated rms
	private double []         last_ymfx =       null;
	private double [][]       last_jt =         null;

	private int               ddisp_index;   // =  G0_INDEX +    NUM_PAIRS; // disparity offset per camera (at least 1 should be disabled)
	private int               ndisp_index;   //  =  ddisp_index + NUM_CAMS;  // disparity offset per camera - none should be disabled
	private int               num_all_pars;  //  = ndisp_index+ NUM_CAMS; // maximal number of parameters
	private final int []      used_cams_map; //  =  new int[NUM_CAMS]; // for each camera index return used index ???

	private int []            used_cams_rmap; // variable-length list of used cameras numbers
	private int [][]          used_pairs_map; // [tile][pair] -1 for unused pairs, >=0 for used ones
	private boolean []        last_common_scale = null; // new boolean[1]; //When switching from common to individual the
	                                                    // scale[0] is cloned, reverse - averaged to [0]
	private boolean []        used_tiles;

	private final int         transform_size;
    private final double [][] corr_wnd;
    private boolean []        used_cameras;
    private boolean [][]      used_pairs;
    private int               ncam_used; // number of used cameras
    private int []            npairs; // number of used pairs per tile
    private int               last_cam; // index of the last camera (special treatment for disparity correction)
    private int               pre_last_cam; // index of the pre-last camera (special treatment for disparity correction)
    private Matrix [][]       m_disp;       
    private Matrix [][]       m_pairs;
    private Matrix [][]       m_pairs_inv;  // inverted m_pairs to calculate x,y -> dd,nd for initial disparity calculation
    
    private int               numTiles = 1;
    private final int         numMax; // number of maximums to process (normally 1, some now 2)
    private Matrix            mddnd; // Matrix to calculate 2 last corrections in disparity direction and 1 ortho from the first ones (normally 2+3=5)

    private int               gaussian_mode = 3; // 0 - parabola, 1 - gaussian, 2 - limited parabola, 4 - limited squared parabola
    private boolean           lazy_eye; // calculate parameters/derivatives for the "lazy eye" parameters
    private boolean []        lazy_eye_dir = null; // after first set will be boolean[2] {adjust_lazyeye_par, adjust_lazyeye_ortho}
    private double [][]       rXY;
    private int []            dd_indices; //normally 5-long (2 * ncam -3), absolute parameter indices for dd_pre_last, dd_last and nd_last

    public int                bad_tile=-1;   // bad tile - exponent got infinite, remove the tile and start over again
                                             // happens in gaussian mode with the convex area around wrong maximum
    public int getNumIter() {
    	return iter;
    }

	public class Sample{ // USED in lwir
		int    tile;   // tile in a cluster
		int    pair;   // pair index
		int    imax;
		int    ix;     // x coordinate in 2D correlation (0.. 2*transform_size-2, center: (transform_size-1)
		int    iy;     // y coordinate in 2D correlation (0.. 2*transform_size-2, center: (transform_size-1)
		double v;      // correlation value at that point
		double w;      // weight
		Sample (
				int    tile,
				int    pair,
				int    imax,
				int    x,      // x coordinate on the common scale (corresponding to the largest baseline), along the disparity axis
				int    y,      //  coordinate in 2D correlation (0.. 2*transform_size-2, center: (transform_size-1)
				double v,      // correlation value at that point
				double w)
			{
			    this.tile = tile;
			    this.pair = pair;
			    this.imax = imax;
				this.ix =  x;
				this.iy =  y;
				this.v =   v;
				this.w =   w;
			}
		@Override
		public String toString() {
			return String.format("tile=%d, p=%d, x=%d, y=%d, v=%f, w=%f", tile, pair, ix, iy, v, w);
		}
	}


	public Corr2dLMA (
				int         numTiles,
				int         numMax,
				Correlation2d correlation2d,
				int         ts, // null - use default table
				double [][] corr_wnd, // may be null
				double [][] rXY, // non-distorted X,Y offset per nominal pixel of disparity
				int         gaussian_mode // 0 - parabola, 1 - Gaussian, 2 - limited parabola, 3 - limited squared parabola
				) {
		// [(-Disparity, A, B, C-A, G0,... G119, {-Disparity2, A2, B2, C2-A2, G0_2,... G119_2}], ddisp_index..., ndisp_index
		this.correlation2d = correlation2d;
		this.num_cams =  correlation2d.getNumSensors();
		this.num_pairs = correlation2d.getNumPairs();
		this.rXY = rXY;
		this.gaussian_mode = gaussian_mode;
		this.used_cams_map  =  new int[num_cams];
		this.numMax =          numMax;
		this.tile_params = G0_INDEX + this.num_pairs;
		this.numTiles = numTiles;
		ddisp_index = this.numTiles * this.tile_params * this.numMax; // ;
		ndisp_index =  ddisp_index + num_cams;  // disparity offset per camera - none should be disable
		num_all_pars = ndisp_index+ num_cams; // maximal number of parameters
		this.transform_size = ts;
		this.corr_wnd = corr_wnd;
	}
	public static double [][] getCorrWnd(int transform_size, double pwr){ // sq = false
//		double pwr = 1.3;
		double [][] corr_wnd = new double[2 * transform_size - 1][2 * transform_size - 1];
		int tsm1 =  transform_size - 1; // 7
		corr_wnd[tsm1][tsm1] = 1.0;
		for (int i = 1; i < transform_size; i++) {
			if (pwr != 1.0) {
				corr_wnd[tsm1 + i][tsm1    ] = Math.pow(Math.cos(Math.PI*i/(2 * transform_size)),pwr);
				corr_wnd[tsm1 - i][tsm1    ] = Math.pow(Math.cos(Math.PI*i/(2 * transform_size)),pwr);
				corr_wnd[tsm1    ][tsm1 + i] = Math.pow(Math.cos(Math.PI*i/(2 * transform_size)),pwr);
				corr_wnd[tsm1    ][tsm1 - i] = Math.pow(Math.cos(Math.PI*i/(2 * transform_size)),pwr);
			} else {
				corr_wnd[tsm1 + i][tsm1    ] = Math.cos(Math.PI*i/(2 * transform_size));
				corr_wnd[tsm1 - i][tsm1    ] = Math.cos(Math.PI*i/(2 * transform_size));
				corr_wnd[tsm1    ][tsm1 + i] = Math.cos(Math.PI*i/(2 * transform_size));
				corr_wnd[tsm1    ][tsm1 - i] = Math.cos(Math.PI*i/(2 * transform_size));

			}
		}
		for (int i = 1; i < transform_size; i++) {
			for (int j = 1; j < transform_size; j++) {
				double d = corr_wnd[tsm1 + i][tsm1] * corr_wnd[tsm1 + j][tsm1];
				corr_wnd[tsm1 + i][tsm1 + j] = d;
				corr_wnd[tsm1 + i][tsm1 - j] = d;
				corr_wnd[tsm1 - i][tsm1 + j] = d;
				corr_wnd[tsm1 - i][tsm1 - j] = d;
			}
		}
		return corr_wnd;
	}
	public void addSample( // x = 0, y=0 - center
			int    tile,
			int    pair,
			int    imax,   // number of maximum (just a hint)
			int    x,      // x coordinate on the common scale (corresponding to the largest baseline), along the disparity axis
			int    y,      // y coordinate (0 - disparity axis)
			double v,      // correlation value at that point
			double w) {    // sample weight
		if ((w > 0) && !Double.isNaN(v)) samples.add(new Sample(tile, pair, imax, x, y, v, w));
	}
	
	public ArrayList<Sample>  filterSamples(
			double [][] disparity_strength) {
		ArrayList<Sample> filtered_samples = new ArrayList<Sample>();
		for (Sample s:samples) {
			if (disparity_strength[s.tile][1] > 0.0) {
				filtered_samples.add(s);
			}
		}
		samples = filtered_samples;
		return samples;
	}

	// maybe in the future - just remove a bad pair?
	private void removeBadTile(
			int ntile,
			int fcam,   // not yet used, maybe try removing just a pair, not all tile?
			int scam) {
		bad_tile = ntile;
		ArrayList<Sample> filtered_samples = new ArrayList<Sample>();
		for (Sample s:samples) {
			if (s.tile != ntile) {
//			if (s.tile != ntile || (s.fcam!= fcam) || (s.scam !=scam)) {
				filtered_samples.add(s);
			}
		}
		samples = filtered_samples;
	}

	public int getBadTile() { // after failed LMA check bad_tile. If >=0 - reinit/restart LMA, samples have already removed bad_tile
		return bad_tile;
	}



	boolean samplesExist() {
		return ((samples !=null) && !samples.isEmpty());
	}

	void setSamples(ArrayList<Sample> samples) {
		this.samples = samples;
	}

	public double [][][] dbgGetSamples(double [][] ds, int mode){
		int numPairs = getNumAllTilesUsedPairs();
		int size = 2* transform_size -1;
		int size2 = size*size;
		double [][][] rslt = new double [numTiles][numPairs][size2];
		for (int nTile = 0; nTile < numTiles; nTile++) {
			for (int np = 0; np < numPairs; np++) {
				for (int i = 0; i < size2; i++) {
					rslt[nTile][np][i] = Double.NaN;
				}
			}
		}
		double [] fx = null;
		if (mode == 2) {
			fx = getFx();
			if (fx == null) return null;
		}
		for (int ns = 0; ns < samples.size(); ns++) {
			Sample s = samples.get(ns);
			if ((ds != null) && Double.isNaN(ds[s.tile][0])) {
				continue;
			}
			double d = Double.NaN;
			if      (mode == 0) d = s.v;
			else if (mode == 1) d = s.w;
			else if (mode == 2) d = fx[ns];
			int inp = used_pairs_map[s.tile][s.pair];
			rslt[s.tile][inp][s.iy*size + s.ix] = d;
		}

		return rslt;
	}
	
	public String [] dbgGetSliceTitles() {
		int [] comb_map = getCombMap();
		int np = getNumAllTilesUsedPairs(); // comb_map[0][0];
		String [] srslt = new String [np];
		for (int npair = 0; npair < comb_map.length; npair++) {
			if (comb_map[npair] >= 0) {
				int [] fs = correlation2d.getPair(npair);
				srslt[comb_map[npair]] = ""+fs[0]+"->"+fs[1];
			}
		}
		return srslt;
	}
	
	public int [] getCombMap(){
		int []comb_map = new int [num_pairs];
		Arrays.fill(comb_map, -1);
		int np = 0;
		for (int npair = 0; npair < num_pairs; npair++) {
			for (int t = 0; t < numTiles; t++) {
				if (used_pairs_map[t][npair] >= 0) {
					comb_map[npair] = np++;
					break;
				}
			}
		}
		return comb_map;
	}

	public int getNumAllTilesUsedPairs() {
		int np = 0;
		for (int npair = 0; npair < num_pairs; npair++) {
			for (int t = 0; t < numTiles; t++) {
				if (used_pairs_map[t][npair] >= 0) {
					np++;
					break;
				}
			}
		}
		return np;
	}
	
	public static Matrix [] getMatrices(double [][] am_disp, int num_cams) {
		Matrix [] m_disp = new Matrix[num_cams];
		for (int n = 0; n < num_cams; n++) {
			double [][] am = {
					{am_disp[n][0],am_disp[n][1]},
					{am_disp[n][2],am_disp[n][3]}};
			m_disp[n] = new Matrix(am);
		}
		return m_disp;
	}

	
	public void setMatrices(double [][] am_disp) {
		m_disp = new Matrix[1][num_cams];
		for (int n = 0; n < num_cams; n++) {
			double [][] am = {
					{am_disp[n][0],am_disp[n][1]},
					{am_disp[n][2],am_disp[n][3]}};
			m_disp[0][n] = new Matrix(am);
		}
	}

	public void setMatrices(double [][][] am_disp) {
		m_disp = new Matrix[am_disp.length][num_cams];
		for (int nt = 0; nt < numTiles; nt++) if (used_tiles[nt]){
			for (int n = 0; n < num_cams; n++) {
				double [][] am = {
						{am_disp[nt][n][0], am_disp[nt][n][1]},
						{am_disp[nt][n][2], am_disp[nt][n][3]}};
				m_disp[nt][n] = new Matrix(am);
			}
		}
	}

	/**
	 * Calculate matrices to find two last offsets in the disparity direction (ddl1, ddl) and one last offset in ortho to disparity,so
	 * 1: sum all dd is zero (not changing actual disparity that is individual to each tile in a cluster)
	 * 2: dd, nd offsets keep average x,y ot zero
	 */
	private void setOffsetMatrices() {
		Matrix A33;
		Matrix A35;
		double [][] aA33 =
			{       {rXY[pre_last_cam][0], rXY[last_cam][0], -rXY[last_cam][1]},   // sum(x) = d
					{rXY[pre_last_cam][1], rXY[last_cam][1],  rXY[last_cam][0]},   // sum(y) = 0
					{1.0,                  1.0,               0.0             }};  // sum(dd) = 0
		A33 = new Matrix(aA33);
		double [][] aA35 = new double [3][2 * ncam_used-3];
		for (int i = 0; i < (ncam_used-2); i++) { // coefficients for dd[i]
			int n  = used_cams_rmap[i];
			aA35[0][i] = -rXY[n][0];
			aA35[1][i] = -rXY[n][1];
			aA35[2][i] = -1.0;
		}
		for (int i = 0; i < (ncam_used-1); i++) { // coefficients for nd[i]
			int n  = used_cams_rmap[i];
			int i1 = ncam_used - 2 + i;
			aA35[0][i1] =  rXY[n][1];
			aA35[1][i1] = -rXY[n][0];
			//				aA35[2][i] =  0.0; // already 0.0
		}
		A35 = new Matrix(aA35);
		mddnd = A33.inverse().times(A35);

	}
	
	public void initDisparity ( // USED in lwir
			double [][] disp_str,         // initial value of disparity
			int         nmax
			) {
		for (int nTile = 0; nTile < numTiles; nTile++) {
			int offs = (nTile * numMax + nmax) * tile_params;			
			if ((disp_str[nTile] != null) && (disp_str[nTile][1] > 0.0) && used_tiles[nTile]) {
				this.all_pars[DISP_INDEX + offs] = -disp_str[nTile][0]; // disp0;
			} else {
				this.all_pars[DISP_INDEX + offs] = 0.0; // disp0; // null pointer
			}
		}
	}

	public void initDisparity ( // USED in lwir
			double [][][] disp_str         // initial value of disparity [max][tile]{disp, strength}
			) {
		for (int nmax = 0; nmax < numMax; nmax++) {
			for (int nTile = 0; nTile < numTiles; nTile++) {
				int offs = (nTile * numMax + nmax) * tile_params;			
				if (    (disp_str[nmax] !=        null) &&
						(disp_str[nmax][nTile] != null) &&
						(disp_str[nmax][nTile][1] > 0.0) &&
						used_tiles[nTile]) {
					this.all_pars[DISP_INDEX + offs] = -disp_str[nmax][nTile][0]; // disp0;
				} else {
					this.all_pars[DISP_INDEX + offs] = 0.0; // disp0; // null pointer
				}
			}
		}
	}
	
	
	
	public boolean setInitialLYOffsets(
			double [][] pair_centers,
    		double      step_weight, // scale corrections
    		double      min_correction,  // exit when maximal XY correction is below
    		boolean debug) {
		if (pair_centers == null) {
			return false;
		}
		// verify that LMA is for a single maximum 
		if (numMax != 1) {
			System.out.println("Number of maximums is not 1, it is "+numMax);
			return false; // not all cameras present
		}
		// Verify that all cameras are used
		for (int ncam = 0; ncam <num_cams; ncam++) if (!used_cameras[ncam]){
			return false; // not all cameras present
		}
		
		boolean apply = true; // false;
		
		double min_corr2 = min_correction * min_correction;
		Matrix[] xy_cam_disp = new Matrix [num_cams]; // per camera - x,y without LY - just from disparity
//		Matrix ddnd_disp = new Matrix(new double [] {all_pars[DISP_INDEX + NTILE0*tile_params], 0.0}, 2);
		Matrix ddnd_disp = new Matrix(new double [] {-all_pars[DISP_INDEX + NTILE0*tile_params], 0.0}, 2);
		for (int ncam = 0; ncam < num_cams; ncam++) {
			xy_cam_disp[ncam] = m_disp[NTILE0][ncam].times(ddnd_disp);
		}
		Matrix[] m_pair_centers = new Matrix[pair_centers.length];
		for (int pair = 0; pair < num_pairs; pair++) if (pair_centers[pair] != null) {
			m_pair_centers[pair] = new Matrix (new double[] {pair_centers[pair][0], pair_centers[pair][1]},2);
			int [] se =correlation2d.getPair(pair);
			m_pair_centers[pair].plusEquals(xy_cam_disp[se[0]]).minusEquals(xy_cam_disp[se[1]]); // TODO: verify signs
		}
		// m_pair_centers - measured xy pair centers, corrected for disparity, so sum(dd) should be 0;
		Matrix [] xy_ly = new Matrix [num_cams];
		Matrix [] xy_ly_new = new Matrix [num_cams];
		for (int ncam = 0; ncam < num_cams; ncam++) {
			xy_ly[ncam] = new Matrix(2,1); // zero [x,Y] column
		}
		int max_try = 101;
		double [] sum_w = new double [num_cams];  
		for (int ntry = 0; ntry < max_try; ntry++) {
			// reset new vectors
			for (int ncam = 0; ncam < num_cams; ncam++) {
				xy_ly_new[ncam] = new Matrix(2,1); // zero [x,Y] column
			}
			for (int pair = 0; pair < num_pairs; pair++) if (m_pair_centers[pair] != null){
				int [] se =correlation2d.getPair(pair);
				xy_ly_new[se[1]].plusEquals(xy_ly[se[0]].plus (m_pair_centers[pair]).times(pair_centers[pair][2]));
				xy_ly_new[se[0]].plusEquals(xy_ly[se[1]].minus(m_pair_centers[pair]).times(pair_centers[pair][2]));
				if (ntry == 0) {
					sum_w[se[1]] += pair_centers[pair][2];
					sum_w[se[0]] += pair_centers[pair][2];
				}
			}
			double max_diff2 = 0;
			for (int ncam = 0; ncam < num_cams; ncam++) {
				xy_ly_new[ncam] =xy_ly[ncam].times(1.0 - step_weight).plus(xy_ly_new[ncam].times(step_weight/sum_w[ncam]));
				double dx =  xy_ly_new[ncam].get(0, 0) - xy_ly[ncam].get(0, 0);
				double dy =  xy_ly_new[ncam].get(1, 0) - xy_ly[ncam].get(1, 0);
				double d2 = dx*dx + dy*dy;
				if (d2 > max_diff2) {
					max_diff2 = d2;
				}
				xy_ly[ncam] = xy_ly_new[ncam];
			}
			if (max_diff2 < min_corr2) {
				break;
			}
		}
		// correct to zero dx, dy
		Matrix s_xy = new Matrix(2,1); // zero [x,Y] column
		for (int ncam = 0; ncam < num_cams; ncam++) {
			s_xy.plusEquals(xy_ly[ncam]);
		}		
		s_xy.timesEquals(1.0/num_cams);
		for (int ncam = 0; ncam < num_cams; ncam++) {
			xy_ly[ncam].minusEquals(s_xy);
		}
		if (debug) {
			System.out.println(String.format("%2s: %6s %6s","","   X", "   Y"));
			for (int ncam = 0; ncam <num_cams; ncam++) {
				System.out.println(String.format("%2d: %6.3f %6.3f", ncam, xy_ly[ncam].get(0, 0), xy_ly[ncam].get(1, 0)));
			}
		}
		// convert to ddnd
		Matrix [] ddnd = new Matrix [num_cams];
		Matrix s_dd = new Matrix(2,1); // zero [x,Y] column
		for (int ncam = 0; ncam < num_cams; ncam++) {
			ddnd[ncam] = m_disp[NTILE0][ncam].inverse().times(xy_ly[ncam]);
			s_dd.plusEquals(ddnd[ncam]); 
		}
		s_dd.set(1, 0, 0);
		s_dd.timesEquals(1.0/num_cams);
		for (int ncam = 0; ncam < num_cams; ncam++) {
			ddnd[ncam].minusEquals(s_dd) ;
		}
		
		if (debug) {
			System.out.println(String.format("%2s: %6s %6s","","  DD", "  ND"));
			for (int ncam = 0; ncam <num_cams; ncam++) {
				System.out.println(String.format("%2d: %6.3f %6.3f", ncam, ddnd[ncam].get(0, 0), ddnd[ncam].get(1, 0)));
			}
		}
		
		// apply corrections
		if (apply) {
			for (int ncam = 0; ncam <num_cams; ncam++) {
				this.all_pars[ddisp_index + ncam] = -ddnd[ncam].get(0, 0);
				this.all_pars[ndisp_index + ncam] = -ddnd[ncam].get(1, 0);
			}
			toVector();
		}
		return true;
	}
	
	public boolean initVector( // USED in lwir
			boolean [] adjust_disparities, // null - adjust all, otherwise - per maximum
			boolean adjust_width,          // adjust width of the maximum -                                           lma_adjust_wm
			boolean adjust_scales,         // adjust 2D correlation scales -                                          lma_adjust_ag
			boolean adjust_ellipse,        // allow non-circular correlation maximums                                 lma_adjust_wy
			boolean adjust_lazyeye_par,    // adjust disparity corrections parallel to disparities                    lma_adjust_wxy
			boolean adjust_lazyeye_ortho,  // obsolete - make == adjust_lazyeye_par adjust disparity corrections orthogonal to disparities                  lma_adjust_ly1
			double [][][] disp_str_all,    // initial value of disparity [max][tile]{disp, strength}
			double  half_width,            // A=1/(half_widh)^2                                                       lma_half_width
			double  cost_lazyeye_par,      // cost for each of the non-zero disparity corrections                     lma_cost_wy
			double  cost_lazyeye_odtho     // cost for each of the non-zero ortho disparity corrections               lma_cost_wxy
			) {
		if (adjust_disparities == null) {
			adjust_disparities = new boolean[numMax];
			Arrays.fill(adjust_disparities, true);
		}
		if (last_common_scale == null) { // first time - same as it was {false, ...,false}
			last_common_scale = new boolean[numMax];
			Arrays.fill(last_common_scale, false);
		}
		
///		double [][] disp_str = disp_str_all[0]; //FIXME: **************** 
		adjust_lazyeye_ortho = adjust_lazyeye_par; // simplify relations for the calculated/dependent parameters
		lazy_eye = adjust_lazyeye_par | adjust_lazyeye_ortho;
		bad_tile = -1;
		used_cameras = new boolean[num_cams];
		boolean [][] used_pairs = new boolean[numTiles][num_pairs];
		// 0-weight values and NaN-s should be filtered on input!
		used_pairs_map = new int [numTiles][num_pairs];
		for (int t = 0; t < numTiles; t++) {
			Arrays.fill(used_pairs_map[t], -1);
		}
		used_tiles = new boolean[numTiles];
		for (Sample s:samples) { // ignore zero-weight samples
			int pair = s.pair;
			int [] fscam = correlation2d.getPair(pair);
			used_cameras[fscam[0]]=true;
			used_cameras[fscam[1]]=true;
			used_tiles[s.tile] = true;
			used_pairs[s.tile][pair]=true; // throws < 0 - wrong pair, f==s
		}
		ncam_used = 0;
		npairs =new int [numTiles];  // pairs in each tile
		for (int i = 0; i < num_cams; i++) {
			used_cams_map[i] = ncam_used;
			if (used_cameras[i]) {
				ncam_used++;
			}
		}
		used_cams_rmap = new int [ncam_used];
		ncam_used = 0;
		for (int i = 0; i < num_cams; i++) {
			if (used_cameras[i]) {
				used_cams_rmap[ncam_used++] = i;
			}
		}
		if (ncam_used < 2) {
			return false;
		}
		last_cam =     (ncam_used >= 1)? used_cams_rmap[ncam_used - 1] :-1;
		pre_last_cam = (ncam_used >= 2)? used_cams_rmap[ncam_used - 2] :-1;

		setOffsetMatrices();

		for (int nTile = 0; nTile < numTiles; nTile++) {
			int [] upmam = new int[num_pairs];
			for (int i = 0; i < num_pairs; i++) {
				upmam[i] = npairs[nTile];
				if (used_pairs[nTile][i])  npairs[nTile]++;
			}
			for (int pair = 0; pair < num_pairs; pair++) {
				int npair = upmam[pair];
				if      (used_pairs[nTile][pair]) used_pairs_map[nTile][pair] = npair;
			}
		}

		this.all_pars = new double[num_all_pars];
		this.par_mask = new boolean[num_all_pars];
		total_tiles = 0; // now
		// per-tile parameters
		for (int nTile = 0; nTile < numTiles; nTile++) if (used_tiles[nTile]) {
			int num_good_max = 0;
			for (int nmax = 0; nmax < numMax; nmax++) {
				if ((disp_str_all[nmax][nTile] != null) && (disp_str_all[nmax][nTile][1] > 0.0)){
					num_good_max++;
				}
			}
			if (num_good_max == numMax) { // only use tiles with all (both) max
				for (int nmax = 0; nmax < numMax; nmax++) {
					int offs = (nTile * numMax + nmax) * tile_params;
					this.all_pars[DISP_INDEX + offs] = -disp_str_all[nmax][nTile][0]; // disp0;
					this.all_pars[A_INDEX    + offs] = 1.0/(half_width * half_width);
					this.all_pars[B_INDEX    + offs] = 0.0;
					this.all_pars[CMA_INDEX  + offs] = 0.0; // C-A
					this.par_mask[DISP_INDEX + offs] = adjust_disparities[nmax];//  true;
					this.par_mask[A_INDEX    + offs] = adjust_width;
					this.par_mask[B_INDEX    + offs] = adjust_ellipse;
					this.par_mask[CMA_INDEX  + offs] = adjust_ellipse;
					for (int i = 0; i <num_pairs; i++) {
						this.par_mask[G0_INDEX + i + offs] = used_pairs[nTile][i] & adjust_scales;
						this.all_pars[G0_INDEX + i + offs] = Double.NaN; // will be assigned later for used - should be for all !
					}
				}
				total_tiles++;
			}
		}
		// common for all tiles parameters
		for (int i = 0; i <num_cams; i++) {
			this.all_pars[ddisp_index + i] = 0.0; // C-A
			this.par_mask[ddisp_index + i] = used_cameras[i] & adjust_lazyeye_par & (i != last_cam) & (i != pre_last_cam);
			this.all_pars[ndisp_index + i] = 0.0; // C-A
			this.par_mask[ndisp_index + i] = used_cameras[i] & adjust_lazyeye_ortho & (i != last_cam);
		}
		int np = samples.size();

		weights = new double [np + 2 * num_cams]; // npairs];
		values =  new double [np + 2 * num_cams]; // npairs];
		for (int i = 0; i < num_cams; i++) {
			weights[np + i] =            (used_cameras[i] & adjust_lazyeye_par)?   (cost_lazyeye_par * numTiles) :   0.0; // ddisp - including last_cameras
			weights[np + num_cams + i] = (used_cameras[i] & adjust_lazyeye_ortho)? (cost_lazyeye_odtho * numTiles) : 0.0; // ndisp - including last_camera
			values [np + i] =            0.0;
			values [np + num_cams + i] = 0.0;
		}

		double sw = 0;
		total_weight = 0.0;

		// use per-tile, per maximum maximal measured value as initial amplitude
		for (int i = 0; i < np; i++) {
			Sample s = samples.get(i);
			weights[i] = s.w;
			total_weight += s.w;
			values[i] =  s.v;
			sw += weights[i];
			int indx = G0_INDEX + s.pair + (s.tile * numMax + s.imax) * tile_params;
			double d = s.v;
			if (this.corr_wnd !=null) {
				d /= this.corr_wnd[s.iy][s.ix];
			}
			if (!(d <= this.all_pars[indx])) this.all_pars[indx] = d; // to include Double.isNaN()
		}

		pure_weight = sw;
		for (int i = 0; i < 2 * num_cams; i++) { // weight of the regularization terms (twice number of cameras, some may be disabled by a mask)
			sw += weights[np + i];
		}
		if (sw != 0.0) {
			double kw = 1.0/sw;
			for (int i = 0; i < weights.length; i++) weights[i] *= kw;
			pure_weight *= kw; // it is now fraction (0..1.0), and weights are normalized
		}

		par_map = new int [par_mask.length];
		int par_indx = 0;
		for (int i = 0; i < par_mask.length; i++) {
			if (par_mask[i]) par_map[i] = par_indx++;
			else par_map[i] = -1;
		}
		dd_indices = lazy_eye? new int [2 * ncam_used -3] : null;
		if (dd_indices != null) {
			int pi = 0;
			for (int i = 0; i < num_cams; i++) {
				if (par_map[ddisp_index + i] >=0) dd_indices[pi++] = par_map[ddisp_index + i];
			}
			for (int i = 0; i < num_cams; i++) {
				if (par_map[ndisp_index + i] >=0) dd_indices[pi++] = par_map[ndisp_index + i];
			}
		}
		toVector();
		return true;
	}

	/**
	 * Initialize parameters values before calling setParMask() that will copy adjustable parameters to a vector
	 * @param disp_str_all initial value of disparity [max][tile]{disp, strength}
	 * @param half_width   initial half-width of a maximum A=1/(half_widh)^2 
	 * @return OK/failure
	 */
	public boolean preparePars(
			double [][][] disp_str_all,
			double  half_width
			) {
		bad_tile = -1;
		used_cameras = new boolean[num_cams];
		used_pairs = new boolean[numTiles][num_pairs];
		// 0-weight values and NaN-s should be filtered on input!
		used_pairs_map = new int [numTiles][num_pairs];
		for (int t = 0; t < numTiles; t++) {
			Arrays.fill(used_pairs_map[t], -1);
		}
		used_tiles = new boolean[numTiles];
		for (Sample s:samples) { // ignore zero-weight samples
			int pair = s.pair;
			int [] fscam = correlation2d.getPair(pair);
			used_cameras[fscam[0]]=true;
			used_cameras[fscam[1]]=true;
			used_tiles[s.tile] = true;
			used_pairs[s.tile][pair]=true; // throws < 0 - wrong pair, f==s
		}
		ncam_used = 0;
		npairs =new int [numTiles];  // pairs in each tile
		for (int i = 0; i < num_cams; i++) {
			used_cams_map[i] = ncam_used;
			if (used_cameras[i]) {
				ncam_used++;
			}
		}
		used_cams_rmap = new int [ncam_used];
		ncam_used = 0;
		for (int i = 0; i < num_cams; i++) {
			if (used_cameras[i]) {
				used_cams_rmap[ncam_used++] = i;
			}
		}
		if (ncam_used < 2) {
			return false;
		}
		last_cam =     (ncam_used >= 1)? used_cams_rmap[ncam_used - 1] :-1;
		pre_last_cam = (ncam_used >= 2)? used_cams_rmap[ncam_used - 2] :-1;

		setOffsetMatrices();

		for (int nTile = 0; nTile < numTiles; nTile++) {
			int [] upmam = new int[num_pairs];
			for (int i = 0; i < num_pairs; i++) {
				upmam[i] = npairs[nTile];
				if (used_pairs[nTile][i])  npairs[nTile]++;
			}
			for (int pair = 0; pair < num_pairs; pair++) {
				int npair = upmam[pair];
				if      (used_pairs[nTile][pair]) used_pairs_map[nTile][pair] = npair;
			}
		}

		this.all_pars = new double[num_all_pars];
		this.par_mask = new boolean[num_all_pars];
		total_tiles = 0; // now
		// per-tile parameters
		for (int nTile = 0; nTile < numTiles; nTile++) if (used_tiles[nTile]) {
			int num_good_max = 0;
			for (int nmax = 0; nmax < numMax; nmax++) {
				if ((disp_str_all[nmax][nTile] != null) && (disp_str_all[nmax][nTile][1] > 0.0)){
					num_good_max++;
				}
			}
			if (num_good_max == numMax) { // only use tiles with all (both) max
				for (int nmax = 0; nmax < numMax; nmax++) {
					int offs = (nTile * numMax + nmax) * tile_params;
					this.all_pars[DISP_INDEX + offs] = -disp_str_all[nmax][nTile][0]; // disp0;
					this.all_pars[A_INDEX    + offs] = 1.0/(half_width * half_width);
					this.all_pars[B_INDEX    + offs] = 0.0;
					this.all_pars[CMA_INDEX  + offs] = 0.0; // C-A
					for (int i = 0; i <num_pairs; i++) {
						this.all_pars[G0_INDEX + i + offs] = Double.NaN; // will be assigned later for used - should be for all !
					}
				}
				total_tiles++;
			} else {
				used_tiles[nTile] = false; /// Was not here before 02/13/2022 
			}
		}
		// common for all tiles parameters
		for (int i = 0; i <num_cams; i++) {
			this.all_pars[ddisp_index + i] = 0.0; // C-A
			this.all_pars[ndisp_index + i] = 0.0; // C-A
		}
		int np = samples.size();

		weights = new double [np + 2 * num_cams]; // npairs];
		values =  new double [np + 2 * num_cams]; // npairs];
		total_weight = 0.0;

		// use per-tile, per maximum maximal measured value as initial amplitude
		for (int i = 0; i < np; i++) {
			Sample s = samples.get(i);
			if (used_tiles[s.tile]) { // normally should always be true
				weights[i] = s.w;
				total_weight += s.w;
				values[i] =  s.v;
				int indx = G0_INDEX + s.pair + (s.tile * numMax + s.imax) * tile_params;
				double d = s.v;
				if (this.corr_wnd !=null) {
					d /= this.corr_wnd[s.iy][s.ix];
				}
				if (!(d <= this.all_pars[indx])) this.all_pars[indx] = d; // to include Double.isNaN()
			}
		}
		for (int i = 0; i < num_cams; i++) {
			values [np + i] =            0.0;
			values [np + num_cams + i] = 0.0;
		}
		// Normalize w/o LY
		if (total_weight != 0.0) {
			double kw = 1.0/total_weight;
			for (int i = 0; i < weights.length; i++) weights[i] *= kw;
			pure_weight = 1.0; // *= kw; // it is now fraction (0..1.0), and weights are normalized
		}
		return true;
	}
	
	/**
	 * Set/modify parameters mask. May be called after preparePars () or after updateFromVector() if LMA was ran
	 * @param adjust_disparities null to adjust all (1 or 2) disparities or a boolean array of per maximum
	 *                           individual disparity adjusts
	 * @param common_scales per-maximum, if true - common scale for all pairs. Null - all individual (old mode)                         
	 * @param adjust_width adjust correlation maximum width      
	 * @param adjust_scales adjust per-pair amplitude
	 * @param adjust_ellipse adjust per-pair maximum shape as an ellipse
	 * @param adjust_lazyeye_par adjust lazy eye parallel to disparity
	 * @param adjust_lazyeye_ortho adjust lazy eye orthogonal to disparity
	 * @param cost_lazyeye_par cost of lazy eye parallel to disparity
	 * @param cost_lazyeye_odtho cost of lazy eye orthogonal to disparity
	 * @return OK/failure
	 */
	public boolean setParMask(
			boolean [] adjust_disparities, // null - adjust all, otherwise - per maximum
			boolean [] common_scale,       // per-maximum, if true - common scale for all pairs
			boolean adjust_width,          // adjust width of the maximum -                                           lma_adjust_wm
			boolean adjust_scales,         // adjust 2D correlation scales -                                          lma_adjust_ag
			boolean adjust_ellipse,        // allow non-circular correlation maximums                                 lma_adjust_wy
			boolean adjust_lazyeye_par,    // adjust disparity corrections parallel to disparities                    lma_adjust_wxy
			boolean adjust_lazyeye_ortho,  // obsolete - make == adjust_lazyeye_par adjust disparity corrections orthogonal to disparities                  lma_adjust_ly1
			double  cost_lazyeye_par,      // cost for each of the non-zero disparity corrections                     lma_cost_wy
			double  cost_lazyeye_odtho     // cost for each of the non-zero ortho disparity corrections               lma_cost_wxy
			) {
		if (adjust_disparities == null) {
			adjust_disparities = new boolean[numMax];
			Arrays.fill(adjust_disparities, true);
		}
		if (common_scale == null) {
			common_scale = new boolean[numMax];
			Arrays.fill(common_scale, false);
		}
		if (last_common_scale == null) { // first time - same as it was {false, ...,false}
			last_common_scale = new boolean[numMax];
			Arrays.fill(last_common_scale, false);
		}
		
		total_tiles = 0; // now
		// per-tile parameters
		for (int nTile = 0; nTile < numTiles; nTile++) if (used_tiles[nTile]) {
			for (int nmax = 0; nmax < numMax; nmax++) {
				int offs = (nTile * numMax + nmax) * tile_params;
				// Do scales need to be changed??
				if (common_scale[nmax] != last_common_scale[nmax]) {
					double cs = 0.0;
					// Do scales need to be averaged into [0]?
					if (common_scale[nmax]) {
						int num_used_pairs = 0;
						for (int np = 0; np <num_pairs; np++) if (used_pairs[nTile][np]) {
							//this.all_pars[G0_INDEX + i + offs] = Double.NaN; // will be assigned later for used - should be for all !
							cs += this.all_pars[G0_INDEX + np + offs];
							num_used_pairs++;
						}
						if (num_used_pairs > 0) {
							cs /= num_used_pairs;
						}
					} else { // need to be copied
						cs = this.all_pars[G0_INDEX + 0 + offs];
					}
					for (int np = 0; np <num_pairs; np++) { //  if (used_pairs[nTile][np]) {
						this.all_pars[G0_INDEX + np + offs] = cs;
					}						
					last_common_scale[nmax] = common_scale[nmax];
				}
				
				this.par_mask[DISP_INDEX + offs] = adjust_disparities[nmax];//  true;
				this.par_mask[A_INDEX    + offs] = adjust_width;
				this.par_mask[B_INDEX    + offs] = adjust_ellipse;
				this.par_mask[CMA_INDEX  + offs] = adjust_ellipse;
				if (last_common_scale[nmax]) {
					this.par_mask[G0_INDEX + 0 + offs] = adjust_scales;
					for (int np = 1; np <num_pairs; np++) {
						this.par_mask[G0_INDEX + np + offs] = false;
					}
				} else  {
					for (int np = 0; np <num_pairs; np++) {
						this.par_mask[G0_INDEX + np + offs] = used_pairs[nTile][np] & adjust_scales;
					}
				}
			}
			total_tiles++;
		}
		// common for all tiles parameters
		for (int i = 0; i <num_cams; i++) {
			this.par_mask[ddisp_index + i] = used_cameras[i] & adjust_lazyeye_par & (i != last_cam) & (i != pre_last_cam);
			this.par_mask[ndisp_index + i] = used_cameras[i] & adjust_lazyeye_ortho & (i != last_cam);
		}
		
		int np = samples.size();
		// First time always, next only if adjust_lazyeye_par or adjust_lazyeye_ortho is changed
		if ((lazy_eye_dir == null) || (lazy_eye_dir[0] != adjust_lazyeye_par)|| (lazy_eye_dir[1] != adjust_lazyeye_ortho)) {
			for (int i = 0; i < num_cams; i++) {
				weights[np + i] =            (used_cameras[i] & adjust_lazyeye_par)?   (cost_lazyeye_par * numTiles) :   0.0; // ddisp - including last_cameras
				weights[np + num_cams + i] = (used_cameras[i] & adjust_lazyeye_ortho)? (cost_lazyeye_odtho * numTiles) : 0.0; // ndisp - including last_camera
			}
			double sw = total_weight;
			for (int i = 0; i < 2 * num_cams; i++) { // weight of the regularization terms (twice number of cameras, some may be disabled by a mask)
				sw += weights[np + i];
			}
			if (sw != 0.0) {
				double kw_pure = total_weight / pure_weight/ sw;
				double kw = 1.0/sw;
				for (int i = 0;  i < np;             i++) weights[i] *= kw_pure;
				for (int i = np; i < weights.length; i++) weights[i] *= kw;
				pure_weight = total_weight * kw; // *= kw; // it is now fraction (0..1.0), and weights are normalized
			}
		}
		par_map = new int [par_mask.length];
		int par_indx = 0;
		for (int i = 0; i < par_mask.length; i++) {
			if (par_mask[i]) par_map[i] = par_indx++;
			else par_map[i] = -1;
		}
		dd_indices = lazy_eye? new int [2 * ncam_used -3] : null;
		if (dd_indices != null) {
			int pi = 0;
			for (int i = 0; i < num_cams; i++) {
				if (par_map[ddisp_index + i] >=0) dd_indices[pi++] = par_map[ddisp_index + i];
			}
			for (int i = 0; i < num_cams; i++) {
				if (par_map[ndisp_index + i] >=0) dd_indices[pi++] = par_map[ndisp_index + i];
			}
		}
		toVector();
		return true;
	}

	
	
	

	public void initMatrices() { // should be called after initVector and after setMatrices
		m_pairs = new Matrix[used_pairs_map.length][num_pairs];
		m_pairs_inv = new Matrix[used_pairs_map.length][num_pairs];
		for (int nTile = 0; nTile < used_pairs_map.length; nTile++) if (used_tiles[nTile]){
			for (int npair = 0; npair < num_pairs; npair++) {
				int [] fs = correlation2d.getPair(npair); // TODO: change used_pairs_map?
				if (used_pairs_map[nTile][npair] >= 0) {
					m_pairs[nTile][npair] =     m_disp[nTile][fs[0]].minus(m_disp[nTile][fs[1]]);
					m_pairs_inv[nTile][npair] = m_pairs[nTile][npair].inverse();
				}
			}
		}
	}


	public void initInvertMatrices() { // should be called after initMatrices only if m_pairs_inv are needed
		m_pairs_inv = new Matrix[used_pairs_map.length][num_pairs];

		for (int nTile = 0; nTile < used_pairs_map.length; nTile++)  if (used_tiles[nTile]){
			for (int npair = 0; npair < num_pairs; npair++) {
				int [] fs = correlation2d.getPair(npair); // TODO: change used_pairs_map?
				if (used_pairs_map[nTile][npair] >= 0) {
					m_pairs_inv[nTile][npair] = m_pairs[nTile][npair].inverse();
				}
			}
		}
	}

	/**
	 * Calculate initial disparity by polynomial approximation. Only works with a single tile now and a single maximum
	 * Note: combo correlation scale is one pixel for 1 pixel disparity for diameter cameras, while standard 
	 * (pre-shift) disparity refers to a quad-camera (square) pixel offset. So the returned disparity should be divided by sqrt(2)
	 * @param corr_wnd_inv_limited window (same as for finding convex) to boost gain in peripheral areas, but not the very marginal ones
	 * @param max_offset maximal abs(disparity) value to trust
	 * @param dgg_title if !=null - generate debug image with this title
	 * @return estimated disparity, strength and half-window sizes
	 */
	public double [] polyDisparity(
			double [] corr_wnd_inv_limited,
			double max_offset, // 5?
			String dbg_title) {
		int center = transform_size -1;
		int corr_size = 2 * transform_size -1;
		initInvertMatrices();
		int nSamples = samples.size();
    	double [][][] mdata = new double [2 * nSamples][3][];
    	double bv;
    	for (int ns = 0; ns < nSamples; ns++) {
    		Sample s = samples.get(ns);
    		if (corr_wnd_inv_limited != null) {
        		bv = s.v * corr_wnd_inv_limited[corr_size * s.iy +s.ix];
    		} else {
        		bv = s.v / corr_wnd[s.iy][s.ix];
    		}
    		bv /=this.all_pars[G0_INDEX + s.pair + s.tile * tile_params]; // FIXME: numMax
    		//corr_wnd
    		int indx = 2 * ns;
    	    mdata[indx  ][0] = new double [2];
    	    mdata[indx+1][0] = new double [2];
    	    mdata[indx  ][1] = new double [1];
    	    mdata[indx+1][1] = new double [1];
    	    mdata[indx  ][2] = new double [1];
    	    mdata[indx+1][2] = new double [1];

    	    double [] aXY = {s.ix - center, s.iy - center};
    	    Matrix mXY = new Matrix(aXY,2);
    	    Matrix mDDND = m_pairs_inv[s.tile][s.pair].times(mXY);

    	    mdata[indx  ][0][0] =  mDDND.get(0, 0); // dd
    	    mdata[indx+1][0][0] =  mdata[indx  ][0][0];
    	    mdata[indx  ][0][1] =  mDDND.get(1, 0); // nd
    	    mdata[indx+1][0][1] = -mdata[indx  ][0][1];
    	    mdata[indx  ][1][0] =  bv;
    	    mdata[indx+1][1][0] =  bv;
    	    mdata[indx  ][2][0] =  s.w;
    	    mdata[indx+1][2][0] =  s.w;
    	}
        double[] approx2d =(new PolynomialApproximation()).quadraticMaxV2dX2Y2XY( // 9 elements - Xc, Yx, f(x,y), A, B, C, D, E, F (from A*x^2 + B*y^2 +C*x*y+...)
				mdata,
				1.0E-30,//25, // 1.0E-15,
				0);
		if ((approx2d == null) || (approx2d[2] < 0.0) || // negative strength
				(approx2d[3] >= 0.0) || // x: min, not max
				(approx2d[4] >= 0.0)) { // y: min, not max
			return null; // Double.NaN;
		}
		if (Math.abs(approx2d[0]) > max_offset) { // > 6 - too far for the maximum
			return null; // Double.NaN;
		}
		// calculate width_x and width_y
		double hwx = Math.sqrt(-approx2d[2]/approx2d[3]);
		double hwy = Math.sqrt(-approx2d[2]/approx2d[4]);
		double [] rslt = {
				-approx2d[0],
				approx2d[2] / Math.sqrt(hwx * hwy)};
//				hwx, hwy};

		if (dbg_title != null) {
			polyDisparityDebug(
					dbg_title,
					mdata,
					10.0, // double        scale,
					100, // int           irad,
					10.0);//sigma);
		}
		return rslt; // [0];
	}

	double [][] polyDisparityDebug(
			String        title,
			double [][][] mdata,
			double        scale,
			int           irad,
			double sigma){
		double [][] dbg_img =     new double [num_cams*num_cams+1][];
		double [][] dbg_weights = new double [num_cams*num_cams+1][];
		String [] titles = new String [num_cams*num_cams+1];
		int size = 2 * irad+1;
		int nSamples = samples.size();
		DoubleGaussianBlur gb = new DoubleGaussianBlur();
		titles[0] = "combined";
		dbg_img[0] =     new double [size*size];
		dbg_weights[0] = new double [size*size];

    	for (int ns = 0; ns < nSamples; ns++) {
    		int indx = 2 * ns;
    		Sample s = samples.get(ns);
    		int np = s.pair; // .fcam* num_cams + s.scam +1;
    		int [] fs = correlation2d.getPair(np);
    		if (dbg_img[np] == null) {
    			dbg_img[np] =     new double [size*size];
    			dbg_weights[np] = new double [size*size];
//    			titles[np]=String.format("%d->%d", s.fcam,s.scam);
    			titles[np]=String.format("%d->%d", fs[0], fs[1]);
    		}
    		double d = mdata[indx][1][0];
    		double w = mdata[indx][2][0];
    		if (w <= 0) {
    			System.out.println("w==0, should not happen");
    			continue;
    		}
    		int ix = (int) Math.round(mdata[indx  ][0][0] * scale + irad);
    		int iy = (int) Math.round(mdata[indx  ][0][1] * scale + irad);
    		if ((ix >= 0) && (iy >= 0) && (ix < size) && (iy <size)) {
    			int offs = iy * size + ix;
    			dbg_img[np][offs] +=     w * d;
    			dbg_weights[np][offs] += w;
    			dbg_img[0][offs] +=     w * d;
    			dbg_weights[0][offs] += w;
    		}
    	}
    	for (int np = 0; np < dbg_img.length; np++) if (dbg_img[np] != null) {
    		gb.blurDouble(dbg_img[np],     size, size, sigma, sigma, 0.01);
    		gb.blurDouble(dbg_weights[np], size, size, sigma, sigma, 0.01);
    		for (int i = 0; i < dbg_img[np].length; i++) {
    			dbg_img[np][i] /= dbg_weights[np][i];
    		}
    	}
    	if (title != null) {
    		(new ShowDoubleFloatArrays()).showArrays(
    				dbg_img,
    				size,
    				size,
    				true,
    				title,
    				titles);
/*
    		(new ShowDoubleFloatArrays()).showArrays(
    				dbg_weights,
    				size,
    				size,
    				true,
    				title+"-weights",
    				titles);
*/
    	}
		return dbg_img;
	}

	private double [] getFxJt(
			double []   vector,
			double [][] jt) { // should be either [vector.length][samples.size()] or null - then only fx is calculated
		if (this.gaussian_mode == 0) return getFxJt_parabola         (vector, jt);
		if (this.gaussian_mode == 1) return getFxJt_gaussian         (vector, jt);
		if (this.gaussian_mode == 2) return getFxJt_parabola_lim     (vector, jt);
		if (this.gaussian_mode == 3) return getFxJt_parabola_squared (vector, jt);
		else                         return getFxJt_parabola         (vector, jt);
	}

	private double [] getFxJt_parabola( // USED in lwir
			double []   vector,
			double [][] jt) { // should be either [vector.length][samples.size()] or null - then only fx is calculated
		if (vector == null) return null;
		double [] av = fromVector(vector);
		Matrix [][] xcam_ycam = new Matrix[numTiles][num_cams];
		double [][][][][] xp_yp = new double[numMax][numTiles][num_cams][num_cams][];
		double [] axc_yc = {transform_size - 1.0, transform_size-1.0};
		Matrix xc_yc = new Matrix(axc_yc, 2);
		double [][] AT = new double [numMax][numTiles]; // av[A_INDEX];
		double [][] BT = new double [numMax][numTiles]; // av[B_INDEX];
		double [][] CT = new double [numMax][numTiles]; // A + av[CMA_INDEX];
		for (int nmax = 0; nmax < numMax; nmax++) {
			for (int nTile = 0; nTile < numTiles; nTile++) if (used_tiles[nTile]){
				int offs = (nTile * numMax + nmax) * tile_params;
				for (int i = 0; i < num_cams; i++) if (used_cameras[i]) {
					double [] add_dnd = {av[DISP_INDEX+ offs]+ av[ddisp_index + i],  av[ndisp_index + i]};
					xcam_ycam[nTile][i] = m_disp[nTile][i].times(new Matrix(add_dnd,2));
				}
				for (int f = 0; f < num_cams; f++) if (used_cameras[f]) {
					for (int s = 0; s < num_cams; s++) if (used_cameras[s]) {
						xp_yp[nmax][nTile][f][s] =xcam_ycam[nTile][f].minus(xcam_ycam[nTile][s]).plus(xc_yc).getColumnPackedCopy();
					}
				}
				AT[nmax][nTile] = av[A_INDEX + offs];
				BT[nmax][nTile] = av[B_INDEX + offs];
				CT[nmax][nTile] = AT[nmax][nTile] + av[CMA_INDEX + offs];
			}
		}
		int num_samples = samples.size();
		double [] fx= new double [num_samples + 2 * num_cams];
//corr_wnd
		for (int ns = 0; ns < num_samples; ns++) {
			Sample s = samples.get(ns);
			//			int pair = pindx[s.fcam][s.scam]; // all pairs, not just used?
			int pair = s.pair; // all pairs, not just used?
			int [] fs = correlation2d.getPair(pair);
			for (int nmax = 0; nmax < numMax; nmax++) {
				int offs = (s.tile * numMax + nmax) * tile_params;
				double A = AT[nmax][s.tile];
				double B = BT[nmax][s.tile];
				double C = CT[nmax][s.tile];

				int cpair = last_common_scale[nmax]? 0: pair;
				double Gp = av[G0_INDEX + cpair + offs]; // either common or individual;

				double Wp = corr_wnd[s.ix][s.iy];
				double WGp = Wp * Gp;
				//			double xmxp = s.ix - xp_yp[s.tile][s.fcam][s.scam][0];
				//			double ymyp = s.iy - xp_yp[s.tile][s.fcam][s.scam][1];
				double xmxp = s.ix - xp_yp[nmax][s.tile][fs[0]][fs[1]][0]; // TODO - change format of xp_yp
				double ymyp = s.iy - xp_yp[nmax][s.tile][fs[0]][fs[1]][1];
				double xmxp2 = xmxp * xmxp;
				double ymyp2 = ymyp * ymyp;
				double xmxp_ymyp = xmxp * ymyp;
				double d = Wp*(1.0 - (A*xmxp2 + 2 * B * xmxp_ymyp + C * ymyp2));
				fx[ns] += d * Gp;
				if (Double.isNaN(fx[ns])) {
					System.out.println("fx["+ns+"]="+fx[ns]);
				}
				if (s.tile > 0) {
					System.out.print("");
				}
				if (jt != null) {
					if (par_map[DISP_INDEX + offs] >= 0) {
						jt[par_map[DISP_INDEX + offs]][ns] = 2 * WGp *
								((A * xmxp + B * ymyp) * m_pairs[s.tile][s.pair].get(0, 0)+
										(B * xmxp + C * ymyp) * m_pairs[s.tile][s.pair].get(1, 0));
					}
					if (par_map[A_INDEX + offs] >= 0) {
						jt[par_map[A_INDEX + offs]][ns] = -WGp*(xmxp2 + ymyp2);
					}
					if (par_map[B_INDEX + offs] >= 0) {
						jt[par_map[B_INDEX + offs]][ns] = -WGp* 2 * xmxp_ymyp;
					}
					if (par_map[CMA_INDEX + offs] >= 0) {
						jt[par_map[CMA_INDEX + offs]][ns] = -WGp* ymyp2;
					}
					
					if (last_common_scale[nmax]) {
						if (par_map[G0_INDEX + 0 + offs] >= 0) {
							jt[par_map[G0_INDEX + 0 + offs]][ns] = d; // (par_mask[G0_INDEX + pair])? d;
						}
					} else {
						for (int p = 0; p < num_pairs; p++) { // par_mask[G0_INDEX + p] as all pairs either used, or not - then npairs == 0
							if (par_map[G0_INDEX + p + offs] >= 0) {
								jt[par_map[G0_INDEX + p + offs]][ns] = (p== pair)? d : 0.0; // (par_mask[G0_INDEX + pair])? d;
							}
						}
					}
					
					if (lazy_eye) {
						if (nmax == 0) { // only zero during first pass, then accumulate only
							for (int f = 0; f < num_cams; f++)  { // -1 for the last_cam and pre_last_cam
								if (par_map[ddisp_index + f] >= 0) jt[par_map[ddisp_index + f]][ns] = 0.0;
								if (par_map[ndisp_index + f] >= 0) jt[par_map[ndisp_index + f]][ns] = 0.0;
							}
						}
						double [] dd_deriv = new double[3]; // derivatives by dependent dd_pre_lars, dd_last and nd_last (calculated on demand) with sign according to first/second in a pair
						//					if ((s.fcam == pre_last_cam)) {
						if ((fs[0] == pre_last_cam)) {
							dd_deriv[0] =   2 * WGp *
									(       (A * xmxp + B * ymyp) * m_disp[s.tile][pre_last_cam].get(0, 0)+
											(B * xmxp + C * ymyp) * m_disp[s.tile][pre_last_cam].get(1, 0));
							//					} else if ((s.scam == pre_last_cam)) {
						} else if ((fs[1] == pre_last_cam)) {
							dd_deriv[0] =  -2 * WGp *
									(       (A * xmxp + B * ymyp) * m_disp[s.tile][pre_last_cam].get(0, 0)+
											(B * xmxp + C * ymyp) * m_disp[s.tile][pre_last_cam].get(1, 0));
						}
						//					if ((s.fcam == last_cam)) {
						if ((fs[0] == last_cam)) {
							dd_deriv[1] =   2 * WGp *
									(       (A * xmxp + B * ymyp) * m_disp[s.tile][last_cam].get(0, 0)+
											(B * xmxp + C * ymyp) * m_disp[s.tile][last_cam].get(1, 0));
							dd_deriv[2] =   2 * WGp *
									(       (A * xmxp + B * ymyp) * m_disp[s.tile][last_cam].get(0, 1)+
											(B * xmxp + C * ymyp) * m_disp[s.tile][last_cam].get(1, 1));

							//					} else if ((s.scam == last_cam)) {
						} else if ((fs[1] == last_cam)) {
							dd_deriv[1] =  -2 * WGp *
									(       (A * xmxp + B * ymyp) * m_disp[s.tile][last_cam].get(0, 0)+
											(B * xmxp + C * ymyp) * m_disp[s.tile][last_cam].get(1, 0));
							dd_deriv[2] =  -2 * WGp *
									(       (A * xmxp + B * ymyp) * m_disp[s.tile][last_cam].get(0, 1)+
											(B * xmxp + C * ymyp) * m_disp[s.tile][last_cam].get(1, 1));
						}
						// now accumulate derivatives:
						// first calculate contributions of the dd, nd directly:
						if (par_map[ddisp_index + fs[0]] >= 0){ // par_map[ddisp_index + last_cam] always <0
							jt[par_map[ddisp_index + fs[0]]][ns] += 2 * WGp *
									((A * xmxp + B * ymyp) * m_disp[s.tile][fs[0]].get(0, 0)+
											(B * xmxp + C * ymyp) * m_disp[s.tile][fs[0]].get(1, 0));
						}
						if (par_map[ddisp_index + fs[1]]>= 0){ // par_map[ddisp_index + last_cam] always <0
							jt[par_map[ddisp_index + fs[1]]][ns] -= 2 * WGp *
									((A * xmxp + B * ymyp) * m_disp[s.tile][fs[1]].get(0, 0)+
											(B * xmxp + C * ymyp) * m_disp[s.tile][fs[1]].get(1, 0));
						}
						if (par_map[ndisp_index + fs[0]] >=0){
							jt[par_map[ndisp_index + fs[0]]][ns] += 2 * WGp *
									(       (A * xmxp + B * ymyp) * m_disp[s.tile][fs[0]].get(0, 1)+
											(B * xmxp + C * ymyp) * m_disp[s.tile][fs[0]].get(1, 1));
						}
						if (par_map[ndisp_index + fs[1]] >= 0) {

							jt[par_map[ndisp_index + fs[1]]][ns] -= 2 * WGp *
									(       (A * xmxp + B * ymyp) * m_disp[s.tile][fs[1]].get(0, 1)+
											(B * xmxp + C * ymyp) * m_disp[s.tile][fs[1]].get(1, 1));
						}
						// now calculate indirect ones through derivatives by dd_pre_last (dd_deriv[0]), dd_last (dd_deriv[1]) and nd_last (dd_deriv[2])
						////    private int []            dd_indices; //normally 5-long (2 * ncam -3), absolute parameter indices for dd_pre_last, dd_last and nd_last
						for (int ddn = 0; ddn < 3; ddn++) if (dd_deriv[ddn] != 0.0) {
							for (int i = 0; i < dd_indices.length; i++) {
								jt[dd_indices[i]][ns] += dd_deriv[ddn] * mddnd.get(ddn, i);
								//							if (Double.isNaN(jt[dd_indices[i]][ns])){
								//								System.out.println("getFxJt_parabola(): jt[dd_indices["+i+"]]["+ns+"] == NaN, dd_indices["+i+"]="+dd_indices[i]);
								//							}
							}
						}
					}
				}
			}
		}
		if (lazy_eye) {
			for (int n = 0; n < num_cams; n++) { // av[ddisp_index +last_cam] and other 2 are already populated
				fx[num_samples +            n] = av[ddisp_index + n];
				fx[num_samples + num_cams + n] = av[ndisp_index + n];
			}

			// and derivatives
			if (jt != null) {
				for (int i = 0; i < num_cams; i++) {
					if ((i != last_cam) && (i != pre_last_cam) && (par_map[ddisp_index + i] >= 0)) {
						jt[par_map[ddisp_index + i]][num_samples + i] = 1.0;
					}
					if ((i != last_cam) && (par_map[ndisp_index + i] >= 0)) {
						jt[par_map[ndisp_index + i] ][num_samples + num_cams + i] = 1.0;
					}
				}
				for (int i = 0; i < dd_indices.length; i++) {
					jt[dd_indices[i]][num_samples +        pre_last_cam] = mddnd.get(0, i);
					jt[dd_indices[i]][num_samples +            last_cam] = mddnd.get(1, i);
					jt[dd_indices[i]][num_samples + num_cams + last_cam] = mddnd.get(2, i);
				}
			}
		}
		return fx;
	}

	private double [] getFxJt_parabola_lim( // USED in lwir
			double []   vector,
			double [][] jt) { // should be either [vector.length][samples.size()] or null - then only fx is calculated
		if (vector == null) return null;
		double [] av = fromVector(vector);
		Matrix [][] xcam_ycam = new Matrix[numTiles][num_cams];
		double [][][][][] xp_yp = new double[numMax][numTiles][num_cams][num_cams][];
		double [] axc_yc = {transform_size - 1.0, transform_size-1.0};
		Matrix xc_yc = new Matrix(axc_yc, 2);
		double [][] AT = new double [numMax][numTiles]; // av[A_INDEX];
		double [][] BT = new double [numMax][numTiles]; // av[B_INDEX];
		double [][] CT = new double [numMax][numTiles]; // A + av[CMA_INDEX];
		for (int nmax = 0; nmax < numMax; nmax++) {
			for (int nTile = 0; nTile < numTiles; nTile++) if (used_tiles[nTile]){
				int offs = (nTile * numMax + nmax) * tile_params;
				for (int i = 0; i < num_cams; i++) if (used_cameras[i]) {
					double [] add_dnd = {av[DISP_INDEX+ offs]+ av[ddisp_index + i],  av[ndisp_index + i]};
					xcam_ycam[nTile][i] = m_disp[nTile][i].times(new Matrix(add_dnd,2));
				}
				for (int f = 0; f < num_cams; f++) if (used_cameras[f]) {
					for (int s = 0; s < num_cams; s++) if (used_cameras[s]) {
						xp_yp[nmax][nTile][f][s] =xcam_ycam[nTile][f].minus(xcam_ycam[nTile][s]).plus(xc_yc).getColumnPackedCopy();
					}
				}
				AT[nmax][nTile] = av[A_INDEX + offs];
				BT[nmax][nTile] = av[B_INDEX + offs];
				CT[nmax][nTile] = AT[nmax][nTile] + av[CMA_INDEX + offs];
			}
		}
		int num_samples = samples.size();
		double [] fx= new double [num_samples + 2 * num_cams];
//corr_wnd
		for (int ns = 0; ns < num_samples; ns++) {
			Sample s = samples.get(ns);
			//			int pair = pindx[s.fcam][s.scam]; // all pairs, not just used?
			int pair = s.pair; // all pairs, not just used?
			int [] fs = correlation2d.getPair(pair);
			for (int nmax = 0; nmax < numMax; nmax++) {
				int offs = (s.tile * numMax + nmax) * tile_params;
				double A = AT[nmax][s.tile];
				double B = BT[nmax][s.tile];
				double C = CT[nmax][s.tile];

				int cpair = last_common_scale[nmax]? 0: pair;
				double Gp = av[G0_INDEX + cpair + offs]; // either common or individual;

				double Wp = corr_wnd[s.ix][s.iy];
				double WGp = Wp * Gp;
				//			double xmxp = s.ix - xp_yp[s.tile][s.fcam][s.scam][0];
				//			double ymyp = s.iy - xp_yp[s.tile][s.fcam][s.scam][1];
				double xmxp = s.ix - xp_yp[nmax][s.tile][fs[0]][fs[1]][0]; // TODO - change format of xp_yp
				double ymyp = s.iy - xp_yp[nmax][s.tile][fs[0]][fs[1]][1];
				double xmxp2 = xmxp * xmxp;
				double ymyp2 = ymyp * ymyp;
				double xmxp_ymyp = xmxp * ymyp;
				double d = Wp*(1.0 - (A*xmxp2 + 2 * B * xmxp_ymyp + C * ymyp2));

				double lim_negative = (d < 0)? 0.0 : 1.0;
				d *= lim_negative;

				fx[ns] += d * Gp;
				if (Double.isNaN(fx[ns])) {
					System.out.println("fx["+ns+"]="+fx[ns]);
				}
				if (s.tile > 0) {
					System.out.print("");
				}
				if (jt != null) {
					if (par_map[DISP_INDEX + offs] >= 0) {
						jt[par_map[DISP_INDEX + offs]][ns] = lim_negative* 2 * WGp *
								//						((A * xmxp + B * ymyp) * m_pairs[s.tile][s.fcam][s.scam].get(0, 0)+
								//						 (B * xmxp + C * ymyp) * m_pairs[s.tile][s.fcam][s.scam].get(1, 0));
								((A * xmxp + B * ymyp) * m_pairs[s.tile][s.pair].get(0, 0)+
										(B * xmxp + C * ymyp) * m_pairs[s.tile][s.pair].get(1, 0));
					}
					if (par_map[A_INDEX + offs] >= 0) {
						jt[par_map[A_INDEX + offs]][ns] = -WGp * (xmxp2 + ymyp2) * lim_negative;
					}
					if (par_map[B_INDEX + offs] >= 0) {
						jt[par_map[B_INDEX + offs]][ns] = -WGp * 2 * xmxp_ymyp * lim_negative;
					}
					if (par_map[CMA_INDEX + offs] >= 0) {
						jt[par_map[CMA_INDEX + offs]][ns] = -WGp * ymyp2 * lim_negative;
					}
/*
					for (int p = 0; p < num_pairs; p++) { // par_mask[G0_INDEX + p] as all pairs either used, or not - then npairs == 0
						if (par_map[G0_INDEX + p + offs] >= 0) {
							jt[par_map[G0_INDEX + p + offs]][ns] = (p== pair)? d : 0.0; // (par_mask[G0_INDEX + pair])? d;
						}
					}
*/					
					if (last_common_scale[nmax]) {
						if (par_map[G0_INDEX + 0 + offs] >= 0) {
							jt[par_map[G0_INDEX + 0 + offs]][ns] = d; // (par_mask[G0_INDEX + pair])? d;
						}
					} else {
						for (int p = 0; p < num_pairs; p++) { // par_mask[G0_INDEX + p] as all pairs either used, or not - then npairs == 0
							if (par_map[G0_INDEX + p + offs] >= 0) {
								jt[par_map[G0_INDEX + p + offs]][ns] = (p== pair)? d : 0.0; // (par_mask[G0_INDEX + pair])? d;
							}
						}
					}
					
					
					
					if (lazy_eye) {
						if (nmax == 0) { // only zero during first pass, then accumulate only //   ****
							for (int f = 0; f < num_cams; f++)  { // -1 for the last_cam and pre_last_cam
								if (par_map[ddisp_index + f] >= 0) jt[par_map[ddisp_index + f]][ns] = 0.0;
								if (par_map[ndisp_index + f] >= 0) jt[par_map[ndisp_index + f]][ns] = 0.0;
							}
						}
						double [] dd_deriv = new double[3]; // derivatives by dependent dd_pre_last (0), dd_last (1) and nd_last (1)
						                                    // (calculated on demand) with sign according to first/second in a pair
						if ((fs[0] == pre_last_cam)) {
							dd_deriv[0] =   2 * WGp * lim_negative * 
									(       (A * xmxp + B * ymyp) * m_disp[s.tile][pre_last_cam].get(0, 0)+
											(B * xmxp + C * ymyp) * m_disp[s.tile][pre_last_cam].get(1, 0));
						} else if ((fs[1] == pre_last_cam)) {
							dd_deriv[0] =  -2 * WGp * lim_negative *
									(       (A * xmxp + B * ymyp) * m_disp[s.tile][pre_last_cam].get(0, 0)+
											(B * xmxp + C * ymyp) * m_disp[s.tile][pre_last_cam].get(1, 0));
						}
						if ((fs[0] == last_cam)) {
							dd_deriv[1] =   2 * WGp * lim_negative *
									(       (A * xmxp + B * ymyp) * m_disp[s.tile][last_cam].get(0, 0)+
											(B * xmxp + C * ymyp) * m_disp[s.tile][last_cam].get(1, 0));
							dd_deriv[2] =   2 * WGp * lim_negative *
									(       (A * xmxp + B * ymyp) * m_disp[s.tile][last_cam].get(0, 1)+
											(B * xmxp + C * ymyp) * m_disp[s.tile][last_cam].get(1, 1));

						} else if ((fs[1] == last_cam)) {
							dd_deriv[1] =  -2 * WGp * lim_negative *
									(       (A * xmxp + B * ymyp) * m_disp[s.tile][last_cam].get(0, 0)+
											(B * xmxp + C * ymyp) * m_disp[s.tile][last_cam].get(1, 0));
							dd_deriv[2] =  -2 * WGp * lim_negative *
									(       (A * xmxp + B * ymyp) * m_disp[s.tile][last_cam].get(0, 1)+
											(B * xmxp + C * ymyp) * m_disp[s.tile][last_cam].get(1, 1));
						}
						// now accumulate derivatives:
						// first calculate contributions of the dd, nd directly:
						if (par_map[ddisp_index + fs[0]] >= 0){ // par_map[ddisp_index + last_cam] always <0
							jt[par_map[ddisp_index + fs[0]]][ns] += 2 * WGp * lim_negative *
									((A * xmxp + B * ymyp) * m_disp[s.tile][fs[0]].get(0, 0)+
											(B * xmxp + C * ymyp) * m_disp[s.tile][fs[0]].get(1, 0));
						}
						if (par_map[ddisp_index + fs[1]]>= 0){ // par_map[ddisp_index + last_cam] always <0
							jt[par_map[ddisp_index + fs[1]]][ns] -= 2 * WGp * lim_negative *
									((A * xmxp + B * ymyp) * m_disp[s.tile][fs[1]].get(0, 0)+
											(B * xmxp + C * ymyp) * m_disp[s.tile][fs[1]].get(1, 0));
						}
						if (par_map[ndisp_index + fs[0]] >=0){
							jt[par_map[ndisp_index + fs[0]]][ns] += 2 * WGp * lim_negative *
									(       (A * xmxp + B * ymyp) * m_disp[s.tile][fs[0]].get(0, 1)+
											(B * xmxp + C * ymyp) * m_disp[s.tile][fs[0]].get(1, 1));
						}
						if (par_map[ndisp_index + fs[1]] >= 0) {
							jt[par_map[ndisp_index + fs[1]]][ns] -= 2 * WGp * lim_negative *
									(       (A * xmxp + B * ymyp) * m_disp[s.tile][fs[1]].get(0, 1)+
											(B * xmxp + C * ymyp) * m_disp[s.tile][fs[1]].get(1, 1));
						}
						// now calculate indirect ones through derivatives by dd_pre_last (dd_deriv[0]), dd_last (dd_deriv[1]) and nd_last (dd_deriv[2])
						////    private int []            dd_indices; //normally 5-long (2 * ncam -3), absolute parameter indices for dd_pre_last, dd_last and nd_last
						for (int ddn = 0; ddn < 3; ddn++) if (dd_deriv[ddn] != 0.0) {
							for (int i = 0; i < dd_indices.length; i++) {
								jt[dd_indices[i]][ns] += dd_deriv[ddn] * mddnd.get(ddn, i); // already includes  lim_negative *
								//							if (Double.isNaN(jt[dd_indices[i]][ns])){
								//								System.out.println("getFxJt_parabola(): jt[dd_indices["+i+"]]["+ns+"] == NaN, dd_indices["+i+"]="+dd_indices[i]);
								//							}
							}
						}
					}
				}
			}
		}
		if (lazy_eye) { // do not include  lim_negative *
			for (int n = 0; n < num_cams; n++) { // av[ddisp_index +last_cam] and other 2 are already populated
				fx[num_samples +            n] = av[ddisp_index + n];
				fx[num_samples + num_cams + n] = av[ndisp_index + n];
			}

			// and derivatives
			if (jt != null) {
				for (int i = 0; i < num_cams; i++) {
					if ((i != last_cam) && (i != pre_last_cam) && (par_map[ddisp_index + i] >= 0)) {
						jt[par_map[ddisp_index + i]][num_samples + i] = 1.0;
					}
					if ((i != last_cam) && (par_map[ndisp_index + i] >= 0)) {
						jt[par_map[ndisp_index + i] ][num_samples + num_cams + i] = 1.0;
					}
				}
				for (int i = 0; i < dd_indices.length; i++) {
					jt[dd_indices[i]][num_samples +        pre_last_cam] = mddnd.get(0, i);
					jt[dd_indices[i]][num_samples +            last_cam] = mddnd.get(1, i);
					jt[dd_indices[i]][num_samples + num_cams + last_cam] = mddnd.get(2, i);
				}
			}
		}
		return fx;
	}
	
	private double [] getFxJt_parabola_squared( // USED in lwir
			double []   vector,
			double [][] jt) { // should be either [vector.length][samples.size()] or null - then only fx is calculated
		if (vector == null) return null;
		double [] av = fromVector(vector);
		Matrix [][] xcam_ycam = new Matrix[numTiles][num_cams];
		double [][][][][] xp_yp = new double[numMax][numTiles][num_cams][num_cams][];
		double [] axc_yc = {transform_size - 1.0, transform_size-1.0};
		Matrix xc_yc = new Matrix(axc_yc, 2);
		double [][] AT = new double [numMax][numTiles]; // av[A_INDEX];
		double [][] BT = new double [numMax][numTiles]; // av[B_INDEX];
		double [][] CT = new double [numMax][numTiles]; // A + av[CMA_INDEX];
		for (int nmax = 0; nmax < numMax; nmax++) { 
			for (int nTile = 0; nTile < numTiles; nTile++) if (used_tiles[nTile]){
				int offs = (nTile * numMax + nmax) * tile_params;
				for (int i = 0; i < num_cams; i++) if (used_cameras[i]) {
					double [] add_dnd = {av[DISP_INDEX+ offs]+ av[ddisp_index + i],  av[ndisp_index + i]};
					xcam_ycam[nTile][i] = m_disp[nTile][i].times(new Matrix(add_dnd,2));
				}
				for (int f = 0; f < num_cams; f++) if (used_cameras[f]) {
					for (int s = 0; s < num_cams; s++) if (used_cameras[s]) {
						xp_yp[nmax][nTile][f][s] =xcam_ycam[nTile][f].minus(xcam_ycam[nTile][s]).plus(xc_yc).getColumnPackedCopy();
					}
				}
				AT[nmax][nTile] = av[A_INDEX + offs];
				BT[nmax][nTile] = av[B_INDEX + offs];
				CT[nmax][nTile] = AT[nmax][nTile] + av[CMA_INDEX + offs];
			}
		}

		int num_samples = samples.size();
		double [] fx= new double [num_samples + 2 * num_cams];
//corr_wnd
		for (int ns = 0; ns < num_samples; ns++) {
			Sample s = samples.get(ns);
			int pair = s.pair; // all pairs, not just used?
//			if (s.pair == 78) {
//				System.out.println("getFxJt_parabola_squared() s.pair="+s.pair);
//			}
			int [] fs = correlation2d.getPair(pair);
			for (int nmax = 0; nmax < numMax; nmax++) {
				int offs = (s.tile * numMax + nmax) * tile_params;
				double A = AT[nmax][s.tile];
				double B = BT[nmax][s.tile];
				double C = CT[nmax][s.tile];
				int cpair = last_common_scale[nmax]? 0: pair;
				double Gp = av[G0_INDEX + cpair + offs]; // either common or individual;
				double Wp = corr_wnd[s.ix][s.iy];
				double WGp = Wp * Gp;
				double xmxp = s.ix - xp_yp[nmax][s.tile][fs[0]][fs[1]][0]; // TODO - change format of xp_yp
				double ymyp = s.iy - xp_yp[nmax][s.tile][fs[0]][fs[1]][1];
				double xmxp2 = xmxp * xmxp;
				double ymyp2 = ymyp * ymyp;
				double xmxp_ymyp = xmxp * ymyp;
				double d = Wp*(1.0 - (A*xmxp2 + 2 * B * xmxp_ymyp + C * ymyp2));

				double lim_negative = (d < 0)? 0.0 : 1.0;
				d *= lim_negative;

				fx[ns] += d * d * Gp;

				// d(d^2)/dp = 2*d *dd/dp

				double lim_negative_2d = 2 * d * lim_negative;

				if (Double.isNaN(fx[ns])) {
					System.out.println("fx["+ns+"]="+fx[ns]);
				}
				if (s.tile > 0) {
					System.out.print("");
				}
				if (jt != null) {
					if (par_map[DISP_INDEX + offs] >= 0) {
						jt[par_map[DISP_INDEX + offs]][ns] = lim_negative_2d* 2 * WGp *
								((A * xmxp + B * ymyp) * m_pairs[s.tile][s.pair].get(0, 0)+
										(B * xmxp + C * ymyp) * m_pairs[s.tile][s.pair].get(1, 0));
					}
					if (par_map[A_INDEX + offs] >= 0) {
						jt[par_map[A_INDEX + offs]][ns] = -WGp * (xmxp2 + ymyp2) * lim_negative_2d;
					}
					if (par_map[B_INDEX + offs] >= 0) {
						jt[par_map[B_INDEX + offs]][ns] = -WGp * 2 * xmxp_ymyp * lim_negative_2d;
					}
					if (par_map[CMA_INDEX + offs] >= 0) {
						jt[par_map[CMA_INDEX + offs]][ns] = -WGp * ymyp2 * lim_negative_2d;
					}

					if (last_common_scale[nmax]) {
						if (par_map[G0_INDEX + 0 + offs] >= 0) {
							jt[par_map[G0_INDEX + 0 + offs]][ns] = d * d; // (par_mask[G0_INDEX + pair])? d;
						}
					} else {
						for (int p = 0; p < num_pairs; p++) { // par_mask[G0_INDEX + p] as all pairs either used, or not - then npairs == 0
							if (par_map[G0_INDEX + p + offs] >= 0) {
								jt[par_map[G0_INDEX + p + offs]][ns] = (p== pair)? (d * d) : 0.0; // (par_mask[G0_INDEX + pair])? d;
							}
						}
					}
					
					if (lazy_eye) {
						if (nmax == 0) { // only zero during first pass, then accumulate only //   ****
							for (int f = 0; f < num_cams; f++)  { // -1 for the last_cam and pre_last_cam
								if (par_map[ddisp_index + f] >= 0) jt[par_map[ddisp_index + f]][ns] = 0.0;
								if (par_map[ndisp_index + f] >= 0) jt[par_map[ndisp_index + f]][ns] = 0.0;
							}
						}
						double [] dd_deriv = new double[3]; // derivatives by dependent dd_pre_lars, dd_last and nd_last (calculated on demand) with sign according to first/second in a pair
						if ((fs[0] == pre_last_cam)) {
							dd_deriv[0] =   2 * WGp * lim_negative_2d * 
									(       (A * xmxp + B * ymyp) * m_disp[s.tile][pre_last_cam].get(0, 0)+
											(B * xmxp + C * ymyp) * m_disp[s.tile][pre_last_cam].get(1, 0));
						} else if ((fs[1] == pre_last_cam)) {
							dd_deriv[0] =  -2 * WGp * lim_negative_2d *
									(       (A * xmxp + B * ymyp) * m_disp[s.tile][pre_last_cam].get(0, 0)+
											(B * xmxp + C * ymyp) * m_disp[s.tile][pre_last_cam].get(1, 0));
						}
						if ((fs[0] == last_cam)) {
							dd_deriv[1] =   2 * WGp * lim_negative_2d *
									(       (A * xmxp + B * ymyp) * m_disp[s.tile][last_cam].get(0, 0)+
											(B * xmxp + C * ymyp) * m_disp[s.tile][last_cam].get(1, 0));
							dd_deriv[2] =   2 * WGp * lim_negative_2d *
									(       (A * xmxp + B * ymyp) * m_disp[s.tile][last_cam].get(0, 1)+
											(B * xmxp + C * ymyp) * m_disp[s.tile][last_cam].get(1, 1));

						} else if ((fs[1] == last_cam)) {
							dd_deriv[1] =  -2 * WGp * lim_negative_2d *
									(       (A * xmxp + B * ymyp) * m_disp[s.tile][last_cam].get(0, 0)+
											(B * xmxp + C * ymyp) * m_disp[s.tile][last_cam].get(1, 0));
							dd_deriv[2] =  -2 * WGp * lim_negative_2d *
									(       (A * xmxp + B * ymyp) * m_disp[s.tile][last_cam].get(0, 1)+
											(B * xmxp + C * ymyp) * m_disp[s.tile][last_cam].get(1, 1));
						}
						// now accumulate derivatives:
						// first calculate contributions of the dd, nd directly:
						if (par_map[ddisp_index + fs[0]] >= 0){ // par_map[ddisp_index + last_cam] always <0
							jt[par_map[ddisp_index + fs[0]]][ns] += 2 * WGp * lim_negative_2d *
									((A * xmxp + B * ymyp) * m_disp[s.tile][fs[0]].get(0, 0)+
											(B * xmxp + C * ymyp) * m_disp[s.tile][fs[0]].get(1, 0));
						}
						if (par_map[ddisp_index + fs[1]]>= 0){ // par_map[ddisp_index + last_cam] always <0
							jt[par_map[ddisp_index + fs[1]]][ns] -= 2 * WGp * lim_negative_2d *
									((A * xmxp + B * ymyp) * m_disp[s.tile][fs[1]].get(0, 0)+
											(B * xmxp + C * ymyp) * m_disp[s.tile][fs[1]].get(1, 0));
						}
						if (par_map[ndisp_index + fs[0]] >=0){
							jt[par_map[ndisp_index + fs[0]]][ns] += 2 * WGp * lim_negative_2d *
									(       (A * xmxp + B * ymyp) * m_disp[s.tile][fs[0]].get(0, 1)+
											(B * xmxp + C * ymyp) * m_disp[s.tile][fs[0]].get(1, 1));
						}
						if (par_map[ndisp_index + fs[1]] >= 0) {

							jt[par_map[ndisp_index + fs[1]]][ns] -= 2 * WGp * lim_negative_2d *
									(       (A * xmxp + B * ymyp) * m_disp[s.tile][fs[1]].get(0, 1)+
											(B * xmxp + C * ymyp) * m_disp[s.tile][fs[1]].get(1, 1));
						}
						// now calculate indirect ones through derivatives by dd_pre_last (dd_deriv[0]), dd_last (dd_deriv[1]) and nd_last (dd_deriv[2])
						////    private int []            dd_indices; //normally 5-long (2 * ncam -3), absolute parameter indices for dd_pre_last, dd_last and nd_last
						for (int ddn = 0; ddn < 3; ddn++) if (dd_deriv[ddn] != 0.0) {
							for (int i = 0; i < dd_indices.length; i++) {
								jt[dd_indices[i]][ns] += dd_deriv[ddn] * mddnd.get(ddn, i); // already includes  lim_negative *
								//							if (Double.isNaN(jt[dd_indices[i]][ns])){
								//								System.out.println("getFxJt_parabola(): jt[dd_indices["+i+"]]["+ns+"] == NaN, dd_indices["+i+"]="+dd_indices[i]);
								//							}
							}
						}
					}
				}
			} // for (int nmax = 0; nmax < numMax; nmax++) {//              ****
		} // for (int ns = 0; ns < num_samples; ns++) {
		if (lazy_eye) { // do not include  lim_negative *
			for (int n = 0; n < num_cams; n++) { // av[ddisp_index +last_cam] and other 2 are already populated
				fx[num_samples +            n] = av[ddisp_index + n];
				fx[num_samples + num_cams + n] = av[ndisp_index + n];
			}

			// and derivatives
			if (jt != null) {
				for (int i = 0; i < num_cams; i++) {
					if ((i != last_cam) && (i != pre_last_cam) && (par_map[ddisp_index + i] >= 0)) {
						jt[par_map[ddisp_index + i]][num_samples + i] = 1.0;
					}
					if ((i != last_cam) && (par_map[ndisp_index + i] >= 0)) {
						jt[par_map[ndisp_index + i] ][num_samples + num_cams + i] = 1.0;
					}
				}
				for (int i = 0; i < dd_indices.length; i++) {
					jt[dd_indices[i]][num_samples +        pre_last_cam] = mddnd.get(0, i);
					jt[dd_indices[i]][num_samples +            last_cam] = mddnd.get(1, i);
					jt[dd_indices[i]][num_samples + num_cams + last_cam] = mddnd.get(2, i);
				}
			}
		}
		return fx;
	}
	
	private double [] getFxJt_gaussian( // USED in lwir
			double []   vector,
			double [][] jt) { // should be either [vector.length][samples.size()] or null - then only fx is calculated
		if (vector == null) return null;
		double [] av = fromVector(vector);
		Matrix [][] xcam_ycam = new Matrix[numTiles][num_cams];
		double [][][][][] xp_yp = new double[numMax][numTiles][num_cams][num_cams][];
		double [] axc_yc = {transform_size - 1.0, transform_size-1.0};
		Matrix xc_yc = new Matrix(axc_yc, 2);
		double [][] AT = new double [numMax][numTiles]; // av[A_INDEX];
		double [][] BT = new double [numMax][numTiles]; // av[B_INDEX];
		double [][] CT = new double [numMax][numTiles]; // A + av[CMA_INDEX];
		for (int nmax = 0; nmax < numMax; nmax++) { 
			for (int nTile = 0; nTile < numTiles; nTile++)  if (used_tiles[nTile]){
				int offs = (nTile * numMax + nmax) * tile_params;
				for (int ncam = 0; ncam < num_cams; ncam++) if (used_cameras[ncam]) {
					double [] add_dnd = {av[DISP_INDEX+ offs]+ av[ddisp_index + ncam],  av[ndisp_index + ncam]};
					xcam_ycam[nTile][ncam] = m_disp[nTile][ncam].times(new Matrix(add_dnd,2));
				}
				for (int f = 0; f < num_cams; f++) if (used_cameras[f]) {
					for (int s = 0; s < num_cams; s++) if (used_cameras[s]) {
						xp_yp[nmax][nTile][f][s] =xcam_ycam[nTile][f].minus(xcam_ycam[nTile][s]).plus(xc_yc).getColumnPackedCopy();
					}
				}
				AT[nmax][nTile] = av[A_INDEX + offs];
				BT[nmax][nTile] = av[B_INDEX + offs];
				CT[nmax][nTile] = AT[nmax][nTile] + av[CMA_INDEX + offs];
			}
		}
		int num_samples = samples.size();
		double [] fx= new double [num_samples + 2 * num_cams];
//corr_wnd
		for (int ns = 0; ns < num_samples; ns++) {
			Sample s = samples.get(ns);
			//			int pair = pindx[s.fcam][s.scam]; // all pairs, not just used?
			int pair = s.pair; // all pairs, not just used?
			int [] fs = correlation2d.getPair(pair);
			for (int nmax = 0; nmax < numMax; nmax++) {
				int offs = (s.tile * numMax + nmax) * tile_params;
				double A = AT[nmax][s.tile];
				double B = BT[nmax][s.tile];
				double C = CT[nmax][s.tile];
				int cpair = last_common_scale[nmax]? 0: pair;
				double Gp = av[G0_INDEX + cpair + offs]; // either common or individual;
				double Wp = corr_wnd[s.ix][s.iy];
				double WGp = Wp * Gp;
				double xmxp = s.ix - xp_yp[nmax][s.tile][fs[0]][fs[1]][0];
				double ymyp = s.iy - xp_yp[nmax][s.tile][fs[0]][fs[1]][1];
				double xmxp2 = xmxp * xmxp;
				double ymyp2 = ymyp * ymyp;
				double xmxp_ymyp = xmxp * ymyp;
				////			double comm = Wp*(1.0 - (A*xmxp2 + 2 * B * xmxp_ymyp + C * ymyp2));
				double exp = Math.exp(-(A*xmxp2 + 2 * B * xmxp_ymyp + C * ymyp2));
				//			if ((exp > 1000.0) || (exp < 0.0000001)) {
				if (exp > 1000.0) {
					//				System.out.println("Unreasonable exp = "+exp);
					removeBadTile(s.tile, fs[0], fs[1]);
					return null; // should re-start LMA
				}

				if (Double.isInfinite(exp)) {
					removeBadTile(s.tile, fs[0], fs[1]);
					return null; // should re-start LMA
				}

				double comm = exp * Wp;
				double WGpexp = WGp*exp;

				fx[ns] += comm * Gp;
				if (Double.isNaN(fx[ns])) {
					System.out.println("fx["+ns+"]="+fx[ns]);
				}
				if (s.tile > 0) {
					System.out.print("");
				}


				if (jt != null) { // Need derivatives too, no just Fx
					// all those before lazy_eye for different nmax do not overlap, so OK to "=", not "+="
					if (par_map[DISP_INDEX + offs] >= 0)       jt[par_map[DISP_INDEX + offs]][ns] = 2 * WGpexp *
							((A * xmxp + B * ymyp) * m_pairs[s.tile][s.pair].get(0, 0)+
									(B * xmxp + C * ymyp) * m_pairs[s.tile][s.pair].get(1, 0));
					if (par_map[A_INDEX + offs] >= 0)          jt[par_map[A_INDEX + offs]][ns] = -WGpexp*(xmxp2 + ymyp2);
					if (par_map[B_INDEX + offs] >= 0)          jt[par_map[B_INDEX + offs]][ns] = -WGpexp* 2 * xmxp_ymyp;
					if (par_map[CMA_INDEX + offs] >= 0)        jt[par_map[CMA_INDEX + offs]][ns] = -WGpexp* ymyp2;
					/*
					for (int p = 0; p < npairs[s.tile]; p++) { // par_mask[G0_INDEX + p] as all pairs either used, or not - then npairs == 0
						if (par_map[G0_INDEX + p + offs] >= 0) jt[par_map[G0_INDEX + p + offs]][ns] = (p== pair)? comm : 0.0; // (par_mask[G0_INDEX + pair])? d;
					}
					*/
					if (last_common_scale[nmax]) {
						if (par_map[G0_INDEX + 0 + offs] >= 0) {
							jt[par_map[G0_INDEX + 0 + offs]][ns] = comm; // (par_mask[G0_INDEX + pair])? d;
						}
					} else {
						for (int p = 0; p < num_pairs; p++) { // par_mask[G0_INDEX + p] as all pairs either used, or not - then npairs == 0
							if (par_map[G0_INDEX + p + offs] >= 0) {
								jt[par_map[G0_INDEX + p + offs]][ns] = (p== pair)? comm : 0.0; // (par_mask[G0_INDEX + pair])? d;
							}
						}
					}
					
					// process ddisp (last camera not used, is equal to minus sum of others to make a sum == 0)
					if (lazy_eye) {
						if (nmax == 0) { // only zero during first pass, then accumulate only
							for (int f = 0; f < num_cams; f++) if (par_map[ddisp_index + f] >= 0) { // -1 for the last_cam
								jt[par_map[ddisp_index + f]][ns] = 0.0;
								jt[par_map[ndisp_index + f]][ns] = 0.0;
							}
						}
						double [] dd_deriv = new double[3]; // derivatives by dependent dd_pre_lars, dd_last and nd_last (calculated on demand) with sign according to first/second in a pair
						if ((fs[0] == pre_last_cam)) {
							dd_deriv[0] =   2 * WGpexp *
									(       (A * xmxp + B * ymyp) * m_disp[s.tile][pre_last_cam].get(0, 0)+
											(B * xmxp + C * ymyp) * m_disp[s.tile][pre_last_cam].get(1, 0));
						} else if ((fs[1] == pre_last_cam)) {
							dd_deriv[0] =  -2 * WGpexp *
									(       (A * xmxp + B * ymyp) * m_disp[s.tile][pre_last_cam].get(0, 0)+
											(B * xmxp + C * ymyp) * m_disp[s.tile][pre_last_cam].get(1, 0));
						}

						if ((fs[0] == last_cam)) {
							dd_deriv[1] =   2 * WGpexp *
									(       (A * xmxp + B * ymyp) * m_disp[s.tile][last_cam].get(0, 0)+
											(B * xmxp + C * ymyp) * m_disp[s.tile][last_cam].get(1, 0));
							dd_deriv[2] =   2 * WGpexp *
									(       (A * xmxp + B * ymyp) * m_disp[s.tile][last_cam].get(0, 1)+
											(B * xmxp + C * ymyp) * m_disp[s.tile][last_cam].get(1, 1));

						} else if ((fs[1] == last_cam)) {
							dd_deriv[1] =  -2 * WGpexp *
									(       (A * xmxp + B * ymyp) * m_disp[s.tile][last_cam].get(0, 0)+
											(B * xmxp + C * ymyp) * m_disp[s.tile][last_cam].get(1, 0));
							dd_deriv[2] =  -2 * WGpexp *
									(       (A * xmxp + B * ymyp) * m_disp[s.tile][last_cam].get(0, 1)+
											(B * xmxp + C * ymyp) * m_disp[s.tile][last_cam].get(1, 1));
						}
						// now accumulate derivatives:
						// first calculate contributions of the dd, nd directly:
						if (par_map[ddisp_index + fs[0]] >= 0){ // par_map[ddisp_index + last_cam] always <0
							jt[par_map[ddisp_index + fs[0]]][ns] += 2 * WGpexp *
									(       (A * xmxp + B * ymyp) * m_disp[s.tile][fs[0]].get(0, 0)+
											(B * xmxp + C * ymyp) * m_disp[s.tile][fs[0]].get(1, 0));
						}
						if (par_map[ddisp_index + fs[1]]>= 0){ // par_map[ddisp_index + last_cam] always <0
							jt[par_map[ddisp_index + fs[1]]][ns] -= 2 * WGpexp *
									(       (A * xmxp + B * ymyp) * m_disp[s.tile][fs[1]].get(0, 0)+
											(B * xmxp + C * ymyp) * m_disp[s.tile][fs[1]].get(1, 0));
						}
						if (par_map[ndisp_index + fs[0]] >=0){
							jt[par_map[ndisp_index + fs[0]]][ns] += 2 * WGpexp *
									(       (A * xmxp + B * ymyp) * m_disp[s.tile][fs[0]].get(0, 1)+
											(B * xmxp + C * ymyp) * m_disp[s.tile][fs[0]].get(1, 1));
						}
						if (par_map[ndisp_index + fs[1]] >= 0) {
							jt[par_map[ndisp_index + fs[1]]][ns] -= 2 * WGpexp *
									(       (A * xmxp + B * ymyp) * m_disp[s.tile][fs[1]].get(0, 1)+
											(B * xmxp + C * ymyp) * m_disp[s.tile][fs[1]].get(1, 1));
						}

						// now calculate indirect ones through derivatives by dd_pre_last (dd_deriv[0]), dd_last (dd_deriv[1]) and nd_last (dd_deriv[2])
						////    private int []            dd_indices; //normally 5-long (2 * ncam -3), absolute parameter indices for dd_pre_last, dd_last and nd_last
						for (int ddn = 0; ddn < 3; ddn++) if (dd_deriv[ddn] != 0.0) {
							for (int i = 0; i < dd_indices.length; i++) {
								jt[dd_indices[i]][ns] += dd_deriv[ddn] * mddnd.get(ddn, i);
								//							if (Double.isNaN(jt[dd_indices[i]][ns])){
								//								System.out.println("getFxJt_gaussian(): jt[dd_indices["+i+"]]["+ns+"] == NaN, dd_indices["+i+"]="+dd_indices[i]);
								//							}
							}
						}
					}
				}
			} // 		for (int nmax = 0; nmax < numMax; nmax++) { 

		} // for (int ns = 0; ns < num_samples; ns++) 
		
		
		if (lazy_eye) {
			for (int n = 0; n < num_cams; n++) { // av[ddisp_index +last_cam] and other 2 are already populated
				fx[num_samples +            n] = av[ddisp_index + n];
				fx[num_samples + num_cams + n] = av[ndisp_index + n];
			}

			// and derivatives
			if (jt != null) {
				for (int i = 0; i < num_cams; i++) {
					if ((i != last_cam) && (i != pre_last_cam) && (par_map[ddisp_index + i] >= 0)) {
						jt[par_map[ddisp_index + i]][num_samples + i] = 1.0;
					}
					if ((i != last_cam) && (par_map[ndisp_index + i] >= 0)) {
						jt[par_map[ndisp_index + i] ][num_samples + num_cams + i] = 1.0;
					}
				}
				for (int i = 0; i < dd_indices.length; i++) {
					jt[dd_indices[i]][num_samples +        pre_last_cam] = mddnd.get(0, i);
					jt[dd_indices[i]][num_samples +            last_cam] = mddnd.get(1, i);
					jt[dd_indices[i]][num_samples + num_cams + last_cam] = mddnd.get(2, i);
				}
			}
		}
		return fx;
	}

	@Deprecated
	private double [] getFxJt_parabola_lim_nodual( // USED in lwir
			double []   vector,
			double [][] jt) { // should be either [vector.length][samples.size()] or null - then only fx is calculated
		if (vector == null) return null;
		double [] av = fromVector(vector);
		Matrix [][] xcam_ycam = new Matrix[numTiles][num_cams];
		double [][][][] xp_yp = new double[numTiles][num_cams][num_cams][];
		double [] axc_yc = {transform_size - 1.0, transform_size-1.0};
		Matrix xc_yc = new Matrix(axc_yc, 2);
		double [] AT = new double [numTiles]; // av[A_INDEX];
		double [] BT = new double [numTiles]; // av[B_INDEX];
		double [] CT = new double [numTiles]; // A + av[CMA_INDEX];
		for (int nTile = 0; nTile < numTiles; nTile++) if (used_tiles[nTile]){
			for (int i = 0; i < num_cams; i++) if (used_cameras[i]) {
				double [] add_dnd = {av[DISP_INDEX+ nTile * tile_params]+ av[ddisp_index + i],  av[ndisp_index + i]};
				xcam_ycam[nTile][i] = m_disp[nTile][i].times(new Matrix(add_dnd,2));
			}
			for (int f = 0; f < num_cams; f++) if (used_cameras[f]) {
				for (int s = 0; s < num_cams; s++) if (used_cameras[s]) {
					xp_yp[nTile][f][s] =xcam_ycam[nTile][f].minus(xcam_ycam[nTile][s]).plus(xc_yc).getColumnPackedCopy();
				}
			}
			AT[nTile] = av[A_INDEX + nTile * tile_params];
			BT[nTile] = av[B_INDEX + nTile * tile_params];
			CT[nTile] = AT[nTile] + av[CMA_INDEX + nTile * tile_params];
		}

		int num_samples = samples.size();
		double [] fx= new double [num_samples + 2 * num_cams];
//corr_wnd
		for (int ns = 0; ns < num_samples; ns++) {
			Sample s = samples.get(ns);
//			int pair = pindx[s.fcam][s.scam]; // all pairs, not just used?
			int pair = s.pair; // all pairs, not just used?
			int [] fs = correlation2d.getPair(pair);
			double A = AT[s.tile];
			double B = BT[s.tile];
			double C = CT[s.tile];

			double Gp = av[G0_INDEX + pair + s.tile * tile_params];
			double Wp = corr_wnd[s.ix][s.iy];
			double WGp = Wp * Gp;
//			double xmxp = s.ix - xp_yp[s.tile][s.fcam][s.scam][0];
//			double ymyp = s.iy - xp_yp[s.tile][s.fcam][s.scam][1];
			double xmxp = s.ix - xp_yp[s.tile][fs[0]][fs[1]][0]; // TODO - change format of xp_yp
			double ymyp = s.iy - xp_yp[s.tile][fs[0]][fs[1]][1];
			double xmxp2 = xmxp * xmxp;
			double ymyp2 = ymyp * ymyp;
			double xmxp_ymyp = xmxp * ymyp;
			double d = Wp*(1.0 - (A*xmxp2 + 2 * B * xmxp_ymyp + C * ymyp2));
			
			double lim_negative = (d < 0)? 0.0 : 1.0;
			d *= lim_negative;
			
			fx[ns] = d * Gp;
			if (Double.isNaN(fx[ns])) {
				System.out.println("fx["+ns+"]="+fx[ns]);
			}
			if (s.tile > 0) {
				System.out.print("");
			}
			if (jt != null) {
				if (par_map[DISP_INDEX + s.tile*tile_params] >= 0) {
					jt[par_map[DISP_INDEX + s.tile*tile_params]][ns] = lim_negative* 2 * WGp *
//						((A * xmxp + B * ymyp) * m_pairs[s.tile][s.fcam][s.scam].get(0, 0)+
//						 (B * xmxp + C * ymyp) * m_pairs[s.tile][s.fcam][s.scam].get(1, 0));
						((A * xmxp + B * ymyp) * m_pairs[s.tile][s.pair].get(0, 0)+
						 (B * xmxp + C * ymyp) * m_pairs[s.tile][s.pair].get(1, 0));
				}
				if (par_map[A_INDEX + s.tile*tile_params] >= 0) {
					jt[par_map[A_INDEX + s.tile*tile_params]][ns] = -WGp * (xmxp2 + ymyp2) * lim_negative;
				}
				if (par_map[B_INDEX + s.tile*tile_params] >= 0) {
					jt[par_map[B_INDEX + s.tile*tile_params]][ns] = -WGp * 2 * xmxp_ymyp * lim_negative;
				}
				if (par_map[CMA_INDEX + s.tile*tile_params] >= 0) {
					jt[par_map[CMA_INDEX + s.tile*tile_params]][ns] = -WGp * ymyp2 * lim_negative;
				}
//				for (int p = 0; p < npairs[s.tile]; p++) { // par_mask[G0_INDEX + p] as all pairs either used, or not - then npairs == 0
				for (int p = 0; p < num_pairs; p++) { // par_mask[G0_INDEX + p] as all pairs either used, or not - then npairs == 0
					if (par_map[G0_INDEX + p + s.tile*tile_params] >= 0) {
						jt[par_map[G0_INDEX + p + s.tile*tile_params]][ns] = (p== pair)? d : 0.0; // (par_mask[G0_INDEX + pair])? d;
					}
				}
				if (lazy_eye) {
					for (int f = 0; f < num_cams; f++)  { // -1 for the last_cam and pre_last_cam
						if (par_map[ddisp_index + f] >= 0) jt[par_map[ddisp_index + f]][ns] = 0.0;
						if (par_map[ndisp_index + f] >= 0) jt[par_map[ndisp_index + f]][ns] = 0.0;
					}
					double [] dd_deriv = new double[3]; // derivatives by dependent dd_pre_lars, dd_last and nd_last (calculated on demand) with sign according to first/second in a pair
					if ((fs[0] == pre_last_cam)) {
						dd_deriv[0] =   2 * WGp * lim_negative * 
								(       (A * xmxp + B * ymyp) * m_disp[s.tile][pre_last_cam].get(0, 0)+
										(B * xmxp + C * ymyp) * m_disp[s.tile][pre_last_cam].get(1, 0));
					} else if ((fs[1] == pre_last_cam)) {
						dd_deriv[0] =  -2 * WGp * lim_negative *
								(       (A * xmxp + B * ymyp) * m_disp[s.tile][pre_last_cam].get(0, 0)+
										(B * xmxp + C * ymyp) * m_disp[s.tile][pre_last_cam].get(1, 0));
					}
					if ((fs[0] == last_cam)) {
						dd_deriv[1] =   2 * WGp * lim_negative *
								(       (A * xmxp + B * ymyp) * m_disp[s.tile][last_cam].get(0, 0)+
										(B * xmxp + C * ymyp) * m_disp[s.tile][last_cam].get(1, 0));
						dd_deriv[2] =   2 * WGp * lim_negative *
								(       (A * xmxp + B * ymyp) * m_disp[s.tile][last_cam].get(0, 1)+
										(B * xmxp + C * ymyp) * m_disp[s.tile][last_cam].get(1, 1));

					} else if ((fs[1] == last_cam)) {
						dd_deriv[1] =  -2 * WGp * lim_negative *
								(       (A * xmxp + B * ymyp) * m_disp[s.tile][last_cam].get(0, 0)+
										(B * xmxp + C * ymyp) * m_disp[s.tile][last_cam].get(1, 0));
						dd_deriv[2] =  -2 * WGp * lim_negative *
								(       (A * xmxp + B * ymyp) * m_disp[s.tile][last_cam].get(0, 1)+
										(B * xmxp + C * ymyp) * m_disp[s.tile][last_cam].get(1, 1));
					}
					// now accumulate derivatives:
					// first calculate contributions of the dd, nd directly:
					if (par_map[ddisp_index + fs[0]] >= 0){ // par_map[ddisp_index + last_cam] always <0
						jt[par_map[ddisp_index + fs[0]]][ns] += 2 * WGp * lim_negative *
								((A * xmxp + B * ymyp) * m_disp[s.tile][fs[0]].get(0, 0)+
										(B * xmxp + C * ymyp) * m_disp[s.tile][fs[0]].get(1, 0));
					}
					if (par_map[ddisp_index + fs[1]]>= 0){ // par_map[ddisp_index + last_cam] always <0
						jt[par_map[ddisp_index + fs[1]]][ns] -= 2 * WGp * lim_negative *
								((A * xmxp + B * ymyp) * m_disp[s.tile][fs[1]].get(0, 0)+
										(B * xmxp + C * ymyp) * m_disp[s.tile][fs[1]].get(1, 0));
					}
					if (par_map[ndisp_index + fs[0]] >=0){
						jt[par_map[ndisp_index + fs[0]]][ns] += 2 * WGp * lim_negative *
								(       (A * xmxp + B * ymyp) * m_disp[s.tile][fs[0]].get(0, 1)+
										(B * xmxp + C * ymyp) * m_disp[s.tile][fs[0]].get(1, 1));
					}
					if (par_map[ndisp_index + fs[1]] >= 0) {

						jt[par_map[ndisp_index + fs[1]]][ns] -= 2 * WGp * lim_negative *
								(       (A * xmxp + B * ymyp) * m_disp[s.tile][fs[1]].get(0, 1)+
										(B * xmxp + C * ymyp) * m_disp[s.tile][fs[1]].get(1, 1));
					}
					// now calculate indirect ones through derivatives by dd_pre_last (dd_deriv[0]), dd_last (dd_deriv[1]) and nd_last (dd_deriv[2])
					////    private int []            dd_indices; //normally 5-long (2 * ncam -3), absolute parameter indices for dd_pre_last, dd_last and nd_last
					for (int ddn = 0; ddn < 3; ddn++) if (dd_deriv[ddn] != 0.0) {
						for (int i = 0; i < dd_indices.length; i++) {
							jt[dd_indices[i]][ns] += dd_deriv[ddn] * mddnd.get(ddn, i); // already includes  lim_negative *
//							if (Double.isNaN(jt[dd_indices[i]][ns])){
//								System.out.println("getFxJt_parabola(): jt[dd_indices["+i+"]]["+ns+"] == NaN, dd_indices["+i+"]="+dd_indices[i]);
//							}
						}
					}
				}
			}
		}
		if (lazy_eye) { // do not include  lim_negative *
			for (int n = 0; n < num_cams; n++) { // av[ddisp_index +last_cam] and other 2 are already populated
				fx[num_samples +            n] = av[ddisp_index + n];
				fx[num_samples + num_cams + n] = av[ndisp_index + n];
			}

			// and derivatives
			if (jt != null) {
				for (int i = 0; i < num_cams; i++) {
					if ((i != last_cam) && (i != pre_last_cam) && (par_map[ddisp_index + i] >= 0)) {
						jt[par_map[ddisp_index + i]][num_samples + i] = 1.0;
					}
					if ((i != last_cam) && (par_map[ndisp_index + i] >= 0)) {
						jt[par_map[ndisp_index + i] ][num_samples + num_cams + i] = 1.0;
					}
				}
				for (int i = 0; i < dd_indices.length; i++) {
					jt[dd_indices[i]][num_samples +        pre_last_cam] = mddnd.get(0, i);
					jt[dd_indices[i]][num_samples +            last_cam] = mddnd.get(1, i);
					jt[dd_indices[i]][num_samples + num_cams + last_cam] = mddnd.get(2, i);
				}
			}
		}
		return fx;
	}
	
	/**
	 * Calculate each defined pair x,y expected offset, assuming only disparity, not lazy eye   
	 * @param corrs - pairs 2D correlations (each in scanline order) - just to determine null/non-null
	 * @param disparity - expected disparity (e.g. from CM)
	 * @param disp_dist - per camera disparity matrix as a 1d (linescan order))
	 * @return per pair x,y expected center offset in 2D correlations or nulls for undefined pairs (or null if already set)
	 */
	public double [][] getPairsOffsets(
			double [][]       corrs,
			boolean []        pair_mask,
			double            disparity,
			double [][]       disp_dist){ // 
		double [][] xy_offsets = new double [corrs.length][];
		if (disp_dist != null) {
			setMatrices(disp_dist);
		}
		for (int pair = 0; pair < xy_offsets.length; pair++) if ((pair < correlation2d.getNumPairs()) && (corrs[pair] != null) && ((pair_mask == null) || pair_mask[pair])){ // OK to calculate for each
 			int [] fscam = correlation2d.getPair(pair); // returns [first_cam, second_cam]
			Matrix mdd_dnd = new Matrix(new double[] {-disparity,  0.0},2);
			Matrix xcam_ycam_f = m_disp[0][fscam[0]].times(mdd_dnd);
			Matrix xcam_ycam_s = m_disp[0][fscam[1]].times(mdd_dnd);
			xy_offsets[pair] = xcam_ycam_f.minus(xcam_ycam_s).getColumnPackedCopy();	
		}
		return xy_offsets;
	}
	
	

	public void printParams() { // not used in lwir
		// to make sure it is updated
		System.out.println();
		for (int np = 0; np < all_pars.length; np++) {
			String parname;
			if      (np >= ndisp_index) parname = PAR_NAME_CORRNDISP + (np - ndisp_index);
			else if (np >= ddisp_index) parname = PAR_NAME_CORRDISP +  (np - ddisp_index);
			else {
				int ntile_max = np / tile_params;
				int ntile = ntile_max / numMax; 
				int nmax =  ntile_max % numMax; 
				int anpr =  np % tile_params;
				if (anpr < G0_INDEX) parname = PAR_NAMES[anpr]+nmax+"-"+ntile;
				else                 parname = PAR_NAME_SCALE +nmax+"-"+ntile + ":"+ (anpr - G0_INDEX);
			}

			System.out.println(String.format("%2d%1s %22s %f",
					np,
					(par_mask[np]?"+":" "),
					parname,
					all_pars[np]));
		}
	}

	public void printStats(
			double [][] ds,
			int clust_width) {
		String form_head, form_line;
		boolean scv = true;
		if (scv) {
			form_head = "%4s,%9s,%9s,%9s,%9s,%9s,%9s,%9s,%9s,%9s,%9s,%9s";
			form_line=  "%4d,%9.4f,%9.4f,%9.4f,%9.4f,%9.4f,%9.4f,%9.4f,%9.4f,%9.4f,%9.4f,%9.4f";
		}else {
			form_head = "%4s%9s%9s%9s%9s%9s%9s%9s%9s%9s%9s%9s";
			form_line=  "%4d%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f";
		}

		double []   rmstile =    getRmsTile();
		double [][] maxmin_amp = getMaxMinAmpTile();
		double [][] maxmin_val = getMaxMinValTile();
		double [][] abc =        getABCTile();
		System.out.println(String.format(form_head,
				"Tile","RMS","disparity","strength","max amp","min amp","diff_amp","max val","min val","diff_val",
				"max_ac","size"));
		for (int tile = 0; tile < ds.length; tile++) {
			double max_ac = Math.max(abc[tile][0], abc[tile][2]);
			double area = Double.NaN;
			if ((abc[tile][0] > 0.0) && (abc[tile][2] > 0.0)) {
				area = 1.0/abc[tile][0] + 1.0/abc[tile][2]; // area of a maximum
			}
			System.out.println(String.format(form_line,
					tile, rmstile[tile] / ((maxmin_amp[tile][0] + maxmin_amp[tile][1])/2), // normalize to average
					ds[tile][0], ds[tile][1],
					maxmin_amp[tile][0], maxmin_amp[tile][1],
					(maxmin_amp[tile][0] - maxmin_amp[tile][1]) / ((maxmin_amp[tile][0] + maxmin_amp[tile][1])/2), // normalize to average
					maxmin_val[tile][0], maxmin_val[tile][1],
					(maxmin_val[tile][0] - maxmin_val[tile][1])/((maxmin_val[tile][0] + maxmin_val[tile][1])/2),
					max_ac, Math.sqrt(area)));
			if ((clust_width > 1) && (((tile+1)%clust_width) == 0)) {
				System.out.println();
			}
		}
	}
	public double [][] getTileStats(){
		double []   rms =        getRmsTile();
		double [][] maxmin_amp = getMaxMinAmpTile();
//		double [][] maxmin_val = getMaxMinValTile();
		double [][] abc =        getABCTile();
		double [][] tileStats = new double [numTiles][];
		for (int tile = 0; tile < numTiles; tile++) if (!Double.isNaN(rms[tile]) && !Double.isNaN(maxmin_amp[tile][0])){
			tileStats[tile] = new double[4];
			double avg = 0.5*(maxmin_amp[tile][0]+maxmin_amp[tile][1]);
			double rrms = rms[tile]/avg;
			double strength = Math.sqrt(avg/rrms);
			double area = Double.NaN;
			if ((abc[tile][0] > 0.0) && (abc[tile][2] > 0.0)) {
				area = 1.0/abc[tile][0] + 1.0/abc[tile][2]; // area of a maximum
			}
			tileStats[tile][0] = rrms;
			tileStats[tile][1] = strength;
			tileStats[tile][2] = Math.max(abc[tile][0], abc[tile][2]);
			tileStats[tile][3] = area;
		}
		return tileStats;
	}



	public void printInputDataFx(boolean show_fx){ // not used in lwir
		if 	(show_fx) {
			Sample s = null;
			double [] fx = getFx();
			if (fx == null) return;
			for (int i = 0; i < fx.length; i++) {
				double fx_pos = fx[i];
				if (i < samples.size()) {
					s = samples.get(i);
					int [] fs = correlation2d.getPair(s.pair);
					System.out.println(String.format("%3d: x=%2d y=%2d v=%9.6f fx=%9.6f w=%9.7f pair=%3d fcam=%2d scam=%2d tile=%d mx=%d",
							i, s.ix, s.iy, s.v, fx_pos, s.w, s.pair, fs[0], fs[1], s.tile, s.imax));
				}
				else {
					System.out.println(String.format("%3d: %2s %2s v=%9.6f fx=%9.6f w=%9.7f", i, "-", "-", this.values[i], fx_pos, this.weights[i]));
				}
			}
		} else {
			int ns =0;
			for (Sample s:samples){
				int [] fs = correlation2d.getPair(s.pair);
				System.out.println(String.format("%3d: x=%2d y=%2d v=%9.6f w=%9.7f pair= %3d fcam=%2d scam=%2d tile=%d  mx=%d",
						ns++, s.ix, s.iy, s.v, s.w, s.pair, fs[0], fs[1], s.tile, s.imax));
			}
		}
	}

	public double [] getRMS() { // USED in lwir
		return last_rms;
	}
	public double [] getGoodOrBadRMS() { // not used in lwir
		return good_or_bad_rms;
	}

	public double [] getAllPars() { // not used in lwir
		return all_pars;
	}

	/*
	public double [] getDisparityStrength(int nTile) { // USED in lwir
		double disparity = -all_pars[DISP_INDEX + nTile * tile_params];
		double sum_amp = 0.0;
		for (int i = 0; i < num_pairs; i++) {
			sum_amp += all_pars[G0_INDEX + i + tile_params * nTile]; // group_weights is normalized
		}
		// protect from weird fitting results
		double max_amp = 0.0;
		for (Sample s: samples) if ((s.v > max_amp) && (s.tile == nTile)) max_amp = s.v;
		if (sum_amp > 1.25 * max_amp) sum_amp = max_amp;
		double [] ds = {disparity, sum_amp};
		return ds;
	}*/
	
	public double [] getDisparityStrength(int nTile) { // USED in lwir
		return getDisparityStrength(nTile, 0);
	}
	
	public double [] getDisparityStrength(int nTile, int nMax) { // USED in lwir
		int offs = (nTile * numMax + nMax) * tile_params;
		double disparity = -all_pars[DISP_INDEX + offs];
		double sum_amp = 0.0;
		for (int i = 0; i < num_pairs; i++) {
///			sum_amp += all_pars[G0_INDEX + i + tile_params * nTile]; // group_weights is normalized
			sum_amp += all_pars[G0_INDEX + i + offs]; // group_weights is normalized
		}
		// protect from weird fitting results
		double max_amp = 0.0;
		for (Sample s: samples) if ((s.v > max_amp) && (s.tile == nTile)) max_amp = s.v;
		if (sum_amp > 1.25 * max_amp) sum_amp = max_amp;
		double [] ds = {disparity, sum_amp};
		return ds;
	}

	public double [][] getDisparityStrength() { // USED in lwir
		double [][] ds = new double [numTiles][];
		for (int tile = 0; tile < numTiles; tile++) {
			ds[tile] = getDisparityStrength(tile);
		}
		return ds;
	}



	public double [] getDisparityStrengthABC(int nTile) {// width = 1/sqrt(all_pars[A_INDEX])
		return getDisparityStrengthABC(nTile, 0);
	}
	
	public double [] getDisparityStrengthABC(int nTile, int nMax) {// width = 1/sqrt(all_pars[A_INDEX])
		double [] ds = getDisparityStrength(nTile, nMax);
		if (ds == null) return null;
		int offs = (nTile * numMax + nMax) * tile_params;
		double [] dsw = {
				ds[0],
				ds[1],
				all_pars[A_INDEX +  offs],
				all_pars[B_INDEX +  offs],
				all_pars[CMA_INDEX+ offs]}; // asymmetry
		return dsw;
	}

	public double [] getLazyEye() {
		double [] rslt = new double [2 * num_cams];
		for (int i = 0; i < rslt.length; i++) {
			rslt[i] = all_pars[ddisp_index + i];
		}
		return rslt;
	}

	public void toVector() { // USED in lwir
		int np = 0;
		for (int i = 0; i < par_mask.length; i++) if (par_mask[i]) np++;
		vector = new double[np];
		np = 0;
		for (int i = 0; i < par_mask.length; i++) if (par_mask[i]) vector[np++] = all_pars[i];
	}

	public void updateFromVector() { // USED in lwir
		int np = 0;

		for (int i = 0; i < par_mask.length; i++) if (par_mask[i]) all_pars[i] = vector[np++];
		// just for reporting
		double [] a5 = new double [2 * ncam_used-3];
		for (int i = 0; i < (ncam_used - 2); i++) {
			a5[i] = all_pars[ddisp_index + used_cams_rmap[i]];
		}
		for (int i = 0; i < (ncam_used - 1); i++) {
			a5[ncam_used - 2 + i] = all_pars[ndisp_index + used_cams_rmap[i]];
		}
		Matrix m5 = new Matrix(a5,a5.length); // single column, normally 5 rows
		Matrix m3 = mddnd.times(m5);
		all_pars[ddisp_index + pre_last_cam] = m3.get(0, 0);
		all_pars[ddisp_index + last_cam] =     m3.get(1, 0);
		all_pars[ndisp_index + last_cam] =     m3.get(2, 0);

	}


	public double [] fromVector(double [] vector) { // mix fixed and variable parameters // USED in lwir
		if ( all_pars == null) return null;
		double [] ap = all_pars.clone();
		int np = 0;
		for (int i = 0; i < par_mask.length; i++) if (par_mask[i]) ap[i] = vector[np++];
		// Fill in missing values (2 last dd-s, 1 last nd)
		double [] a5 = new double [2 * ncam_used-3];
		for (int i = 0; i < (ncam_used - 2); i++) {
			a5[i] = ap[ddisp_index + used_cams_rmap[i]];
		}
		for (int i = 0; i < (ncam_used - 1); i++) {
			a5[ncam_used - 2 + i] = ap[ndisp_index + used_cams_rmap[i]];
		}
		Matrix m5 = new Matrix(a5,a5.length); // single column, normally 5 rows
		Matrix m3 = mddnd.times(m5);
		ap[ddisp_index + pre_last_cam] = m3.get(0, 0);
		ap[ddisp_index + last_cam] =     m3.get(1, 0);
		ap[ndisp_index + last_cam] =     m3.get(2, 0);
		return ap;
	}

////////////////////////////////////////////////////////////////////////

	private boolean debugJt(
    		double      delta,
    		double []   vector) {
    	int num_points = this.values.length;
    	int num_pars = vector.length;
    	double [] max_diff = new double [num_pars];

//    	delta = 0.001;
    	// FIXME numMax
    	double [][] jt =       new double [num_pars][num_points];
    	double [][] jt_delta = new double [num_pars][num_points];
    	double [] fx = getFxJt( vector,jt);
    	if (fx == null) return false;
    	if (getFxJt(delta, vector,jt_delta) == null) return false;
    	System.out.println("Test of jt-jt_delta difference,  delta = "+delta+ ":");
    	System.out.print(String.format("Til P %3s: %10s ", "#", "fx"));
    	for (int anp = 0; anp< all_pars.length; anp++) if(par_mask[anp]){
			String parname;
			if      (anp >= ndisp_index) parname = PAR_NAME_CORRNDISP + (anp - ndisp_index);
			else if (anp >= ddisp_index) parname = PAR_NAME_CORRDISP +  (anp - ddisp_index);
			else {
				int ntile_max = anp / tile_params;
				int ntile = ntile_max / numMax; 
				int nmax =  ntile_max % numMax; 
				int anpr =  anp % tile_params;
				if (anpr < G0_INDEX) parname = PAR_NAMES[anpr]+nmax+"-"+ntile;
				else                 parname = PAR_NAME_SCALE +nmax+"-"+ntile + ":"+ (anpr - G0_INDEX);
			}
        	System.out.print(String.format("| %16s ", parname));
    	}
    	System.out.println();
    	int npair0 = -1;
    	for (int i = 0; i < num_points; i++) {
    		if (i < samples.size()) {
//    			int [] fs = correlation2d.getPair(samples.get(i).pair);
//        		int npair = used_pairs_map[samples.get(i).tile][samples.get(i).fcam][samples.get(i).scam];
//        		int npair = used_pairs_map[samples.get(i).tile][fs[0]][fs[1]];
        		int npair = used_pairs_map[samples.get(i).tile][samples.get(i).pair];
        		if (npair !=npair0) {
        			if (npair0 >=0) System.out.println();
        			npair0 = npair;
        		}
    			System.out.print(String.format("%3d %1d %3d: %10.7f ",samples.get(i).tile, npair, i, fx[i]));
    		} else {
            	System.out.print(String.format(" -  - %3d: %10.7f ", i, fx[i]));
    		}
        	for (int np = 0; np < num_pars; np++) {
//            	System.out.print(String.format("|%8.5f %8.5f ", jt_delta[np][i], 1000*(jt[np][i] - jt_delta[np][i])));
            	System.out.print(String.format("|%8.5f %8.5f ", jt_delta[np][i], 1.0 * (jt[np][i] - jt_delta[np][i])));
            	double adiff = Math.abs(jt[np][i] - jt_delta[np][i]);
            	if (adiff > max_diff[np]) {
            		max_diff[np] = adiff;
            	}
        	}
        	System.out.println();
    	}
    	double tmd = 0.0;
    	for (int np = 0; np < num_pars; np++) {
    		if (max_diff[np] > tmd) tmd= max_diff[np];
    	}
//        System.out.print(String.format("      %15s ", "Maximal diff:"));
        System.out.print(String.format("Max diff.(%10.5f):", tmd));

    	for (int np = 0; np < num_pars; np++) {
        	System.out.print(String.format("|%8s %8.5f ", "1/1000×",  1000*max_diff[np]));
    	}
    	System.out.println();
    	return true;
    }


	public double [] getFxJt( // not used in lwir
			double      delta, // for testing derivatives: calculates as delta-F/delta_x
			double []   vector,
			double [][] jt) { // should be either [vector.length][samples.size()] or null - then only fx is calculated
		double [] fx0=getFxJt(vector,null);
		if (fx0 == null) return null;
		for (int np = 0; np < vector.length; np++) {
			double [] vector1 = vector.clone();
			vector1[np]+= delta;
			double [] fxp=getFxJt(vector1,null);
			if (fxp == null) return null;
			vector1 = vector.clone();
			vector1[np]-= delta;
			double [] fxm=getFxJt(vector1,null);
			if (fxm == null) return null;
			jt[np] = new double [fxp.length];
			for (int i = 0; i < fxp.length; i++) {
				jt[np][i] = (fxp[i] - fxm[i])/delta/2;
			}
		}
		return fx0;
	}



	public double [] getFx() { // not used in lwir
		return getFxJt(this.vector, null);
	}


	public double [][] getWJtJlambda( // USED in lwir
			double      lambda,
			double [][] jt){
		int num_pars = jt.length;
		int nup_points = (num_pars > 0)? jt[0].length : 0;
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
	public double [] getWYmFxRms( // USED in lwir
			double []   fx) { // will be replaced with y-fx
		int num_samples = samples.size();
		int num_points = fx.length; // includes some extra for regularization
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

	// has common threshold for scale ratios for foreground and background corr. maximums
	public double [][] lmaDisparityStrength(
			boolean bypass_tests,     // keep even weak for later analysis. Normally - only in test mode
			double  lmas_min_amp,     // minimal ratio of minimal pair correlation amplitude to maximal pair correlation amplitude
			double  lma_max_rel_rms,  // maximal relative (to average max/min amplitude LMA RMS) // May be up to 0.3)
			double  lma_min_strength, // minimal composite strength (sqrt(average amp squared over absolute RMS)
			double  lma_min_max_ac,   // minimal of A and C coefficients maximum (measures sharpest point/line)
			double  lma_min_min_ac,   // minimal of A and C coefficients minimum (measures sharpest point)
			double  lma_max_area,     // maximal half-area (if > 0.0)
			double  lma_str_scale,    // convert lma-generated strength to match previous ones - scale
			double  lma_str_offset,   // convert lma-generated strength to match previous ones - add to result
			double  lma_ac_offset     // add to a,c coefficients for near-lines where A,C could become negative because of window 
			){
		return lmaDisparityStrengths(
				bypass_tests, //boolean bypass_tests,     // keep even weak for later analysis. Normally - only in test mode
				lmas_min_amp,     // double  lmas_min_amp_fg,  // minimal ratio of minimal pair correlation amplitude to maximal pair correlation amplitude
				lmas_min_amp,     // double  lmas_min_amp_bg,  // Same for bg correlation max (only used for multi-max)
				lma_max_rel_rms,  // maximal relative (to average max/min amplitude LMA RMS) // May be up to 0.3)
				lma_min_strength, // minimal composite strength (sqrt(average amp squared over absolute RMS)
				lma_min_max_ac,   // minimal of A and C coefficients maximum (measures sharpest point/line)
				lma_min_min_ac,   // minimal of A and C coefficients minimum (measures sharpest point)
				lma_max_area,     // maximal half-area (if > 0.0)
				lma_str_scale,    // convert lma-generated strength to match previous ones - scale
				lma_str_offset,    // convert lma-generated strength to match previous ones - add to result
				false,             // boolean dbg_mode
				lma_ac_offset     // add to a,c coefficients for near-lines where A,C could become negative because of window 
				)[0];
	}

	// has separate thresholds for scale ratios for foreground and background corr. maximums
	public double [][][] lmaDisparityStrengths(
			boolean bypass_tests,     // keep even weak for later analysis. Normally - only in test mode
			double  lmas_min_amp_fg,  // minimal ratio of minimal pair correlation amplitude to maximal pair correlation amplitude
			double  lmas_min_amp_bg,  // Same for bg correlation max (only used for multi-max)
			double  lma_max_rel_rms,  // maximal relative (to average max/min amplitude LMA RMS) // May be up to 0.3)
			double  lma_min_strength, // minimal composite strength (sqrt(average amp squared over absolute RMS)
			double  lma_min_max_ac,   // maximal of A and C coefficients minimum (measures sharpest point/line)
			double  lma_min_min_ac,   // minimal of A and C coefficients minimum 
			double  lma_max_area,     // maximal half-area (if > 0.0) // 20
			double  lma_str_scale,    // convert lma-generated strength to match previous ones - scale
			double  lma_str_offset,   // convert lma-generated strength to match previous ones - add to result
			boolean dbg_mode,
			double  lma_ac_offset     // add to a,c coefficients for near-lines where A,C could become negative because of window 
			){
		
		double [][][] ds =         new double[numMax][numTiles][dbg_mode? 14 : 3];
		int reason_index = (dbg_mode && bypass_tests)? 13:-1;
		double []   rms =        getRmsTile();
		for (int nmax = 0; nmax < numMax; nmax++) {
			double [][] maxmin_amp = getMaxMinAmpTile(nmax); // nmax
			double [][] abc =        getABCTile(nmax);       // nmax
			for (int tile = 0; tile < numTiles; tile++) {
				abc[tile][0] += lma_ac_offset;
				abc[tile][2] += lma_ac_offset;
				int offs = (tile * numMax + nmax) * tile_params;
				ds[nmax][tile][0] = Double.NaN;
				if (Double.isNaN(maxmin_amp[tile][0])) {
					continue;
				}
				if (maxmin_amp[tile][1] < 0.0) {
					continue; // inverse maximum - discard tile
				}
				double lmas_min_amp = lmas_min_amp_fg;
				double disparity = -all_pars[DISP_INDEX + offs];
				if (numMax > 1) {
					for (int nmax1 = 0; nmax1 < numMax; nmax1++) if (nmax1 != nmax){
						int offs1 = (tile * numMax + nmax1) * tile_params;
						double disparity1 = -all_pars[DISP_INDEX + offs1];
						if (disparity1 > disparity) {
							lmas_min_amp = lmas_min_amp_bg;
							break;
						}
					}
				}
				if ((maxmin_amp[tile][1]/maxmin_amp[tile][0]) < lmas_min_amp) {
					continue; // inverse maximum - discard tile
				}
				double avg = 0.50*(maxmin_amp[tile][0]+maxmin_amp[tile][1]); // max_min[1] can be negative - filter it out?
				double rrms = rms[tile]/avg;
				if ((lma_max_rel_rms > 0.00) && (rrms > lma_max_rel_rms)) {
					if (bypass_tests) {
						if ((reason_index > 0) && (ds[nmax][tile][reason_index] == 0)) ds[nmax][tile][reason_index] = 1;
					} else {
						continue;
					}
				}
				if (Math.max(abc[tile][0], abc[tile][2]) < (lma_min_max_ac + lma_ac_offset)) { // so old lma_min_max_ac will stay
					if (bypass_tests) {
						if ((reason_index > 0) && (ds[nmax][tile][reason_index] == 0)) ds[nmax][tile][reason_index] = 2;
					} else {
						continue;
					}
				}
				if ((lma_min_min_ac > 0.0) && ((abc[tile][0] < lma_min_min_ac) || (abc[tile][2] < lma_min_min_ac))){
					if (bypass_tests) {
						if ((reason_index > 0) && (ds[nmax][tile][reason_index] == 0)) ds[nmax][tile][reason_index] = 3;
					} else {
						continue; // too large a or c (not sharp along at least one direction)
					}
				}
				double area = 1000.0;
				if ((abc[tile][0] > 0.0) && (abc[tile][2] > 0.0)) {
					//					double area_old = 1.0/abc[tile][0] + 1.0/abc[tile][2]; // area of a maximum
					area =  1.0/Math.sqrt(abc[tile][0] * abc[tile][2]);
					if ((lma_max_area > 0) && (area > lma_max_area)) {
						if (bypass_tests) {
						    if ((reason_index > 0) && (ds[nmax][tile][reason_index] == 0)) ds[nmax][tile][reason_index] = 4;
						} else {
						    continue; // too wide maximum
						}
					}
				} else {
					if (bypass_tests) {
						if ((reason_index > 0) && (ds[nmax][tile][reason_index] == 0)) ds[nmax][tile][reason_index] = 5;
					} else {
						continue;// not a maximum (keep area = 1000.0)
					}
				}
				double strength = Math.sqrt(avg/rrms);
//				double disparity = -all_pars[DISP_INDEX + offs];
				if ((strength < lma_min_strength) || Double.isNaN(disparity)) {
					if (bypass_tests) {
					    if ((reason_index > 0) && (ds[nmax][tile][reason_index] == 0)) ds[nmax][tile][reason_index] = 6;
					} else {
					    continue;
					}
				}
				double ac = abc[tile][0] * abc[tile][2];
				if (ac < 0) {
					if (bypass_tests) {
					    if ((reason_index > 0) && (ds[nmax][tile][reason_index] == 0)) ds[nmax][tile][reason_index] = 7;
					} else {
					    continue;
					}
				}
				double strength1 = Math.sqrt(strength * Math.sqrt(ac)); // / area ); // new strength
				ds[nmax][tile][0] = disparity;
				ds[nmax][tile][1] = (strength1 * lma_str_scale) + lma_str_offset;
				ds[nmax][tile][2] = strength1; // as is
				if (ds[nmax][tile].length > 3) {
					ds[nmax][tile][3] = area;
					ds[nmax][tile][4] = ac;
					ds[nmax][tile][5] = Math.min(abc[tile][0],abc[tile][2]);
					ds[nmax][tile][6] = Math.max(abc[tile][0],abc[tile][2]);
					ds[nmax][tile][7] = abc[tile][0];
					ds[nmax][tile][8] = abc[tile][2];
					ds[nmax][tile][9] = abc[tile][1];
					ds[nmax][tile][10] = strength;
					ds[nmax][tile][11] = rrms;
					ds[nmax][tile][12] = rms[tile];
				}
			}
		}
		return ds;
	}
	
	public double [][] lmaDisparityStrengthLY(
			double  lma_max_rel_rms,  // maximal relative (to average max/min amplitude LMA RMS) // May be up to 0.3)
			double  lma_min_strength, // minimal composite strength (sqrt(average amp squared over absolute RMS)
			double  lma_min_max_ac,   // minimal of A and C coefficients maximum (measures sharpest point/line)
			double  lma_min_min_ac,   // minimal of A and C coefficients minimum (measures sharpest point)
			double  lma_max_area,     // maximal half-area (if > 0.0)
			double  lma_str_scale,    // convert lma-generated strength to match previous ones - scale
			double  lma_str_offset    // convert lma-generated strength to match previous ones - add to result
			){
		double [][] ds =         new double[numTiles][2];
		double []   rms =        getRmsTile();
		double [][] maxmin_amp = getMaxMinAmpTile();
		double [][] abc =        getABCTile();
		if (numMax != 1) {
			System.out.println("lmaDisparityStrengthLY(): numMax = "+numMax+", while only numMax==1 is supported");
			return null;
		}
		for (int tile = 0; tile < numTiles; tile++) {
			ds[tile][0] = Double.NaN;
			if (Double.isNaN(maxmin_amp[tile][0])) {
				continue;
			}
			double avg = 0.5*(maxmin_amp[tile][0]+maxmin_amp[tile][1]);
			double rrms = rms[tile]/avg;
			if (((lma_max_rel_rms > 0.0) && (rrms > lma_max_rel_rms)) ||
					(Math.max(abc[tile][0], abc[tile][2]) < lma_min_max_ac) ||
					(Math.min(abc[tile][0], abc[tile][2]) < lma_min_min_ac)) {
				continue;
			}
			if (lma_max_area > 0) {
				if ((abc[tile][0] > 0.0) && (abc[tile][2] > 0.0)) {
//					double area_old = 1.0/abc[tile][0] + 1.0/abc[tile][2]; // area of a maximum
					double area =  1.0/Math.sqrt(abc[tile][0] * abc[tile][2]);
					if (area > lma_max_area) {
						continue; // too wide maximum
					}
				} else {
					continue; // not a maximum
				}


			}
			double strength = Math.sqrt(avg/rrms);
			double disparity = -all_pars[DISP_INDEX + tile*tile_params];
			if ((strength < lma_min_strength) || Double.isNaN(disparity)) {
				continue;
			}
			strength = Math.sqrt(strength * Math.sqrt(abc[tile][0] * abc[tile][2])); // / area ); // new strength
			ds[tile][0] = disparity;
			ds[tile][1] = (strength * lma_str_scale) + lma_str_offset;
		}
		return ds;
	}
	
	
	
	public double [][] lmaGetExtendedStats(
			double  lma_max_rel_rms,  // maximal relative (to average max/min amplitude LMA RMS) // May be up to 0.3)
			double  lma_min_strength, // minimal composite strength (sqrt(average amp squared over absolute RMS)
			double  lma_min_max_ac,   // minimal of A and C coefficients maximum (measures sharpest point/line)
			double  lma_min_min_ac,   // minimal of A and C coefficients minimum (measures sharpest point)
			double  lma_max_area,     // maximal half-area (if > 0.0)
			double  lma_str_scale,    // convert lma-generated strength to match previous ones - scale
			double  lma_str_offset    // convert lma-generated strength to match previous ones - add to result
			){
		if (numMax != 1) {
			System.out.println("lmaGetExtendedStats(): numMax = "+numMax+", while only numMax==1 is supported");
			return null;
		}
//		double [][] ds =         new double[numTiles][2];
		double [][] ext_stats =  new double[numTiles][11];
		double []   rms =        getRmsTile();
		double [][] maxmin_amp = getMaxMinAmpTile();
		double [][] abc =        getABCTile();
		for (int tile = 0; tile < numTiles; tile++) {
			if (Double.isNaN(maxmin_amp[tile][0])) {
				if (ext_stats[tile][10] == 0) ext_stats[tile][10] = 1;
			}
			ext_stats[tile][0] = maxmin_amp[tile][0];

			double avg = 0.5*(maxmin_amp[tile][0]+maxmin_amp[tile][1]);
			double rrms = rms[tile]/avg;
			ext_stats[tile][1] = avg;
			ext_stats[tile][2] = rrms; // -lma_max_rel_rms;
//			ext_stats[tile][3] = Math.max(abc[tile][0], abc[tile][2])-lma_min_ac;
			
			if (((lma_max_rel_rms > 0.0) && (rrms > lma_max_rel_rms)) ||
					(Math.max(abc[tile][0], abc[tile][2]) < lma_min_max_ac) ||
					(Math.min(abc[tile][0], abc[tile][2]) < lma_min_min_ac)) {
				if (ext_stats[tile][10] == 0) ext_stats[tile][10] = 2;
//				continue;
			}
			double area =  Double.NaN;
			double area_old = Double.NaN;
			if (lma_max_area > 0) {
				ext_stats[tile][5] = abc[tile][0];
				ext_stats[tile][6] = abc[tile][2];

				if ((abc[tile][0] > 0.0) && (abc[tile][2] > 0.0)) {
					area_old = 1.0/abc[tile][0] + 1.0/abc[tile][2]; // area of a maximum
					area =  1.0/Math.sqrt(abc[tile][0] * abc[tile][2]);
					if (area > lma_max_area) {
						if (ext_stats[tile][10] == 0) ext_stats[tile][10] = 3;
//						continue; // too wide maximum
					}
				} else {
					if (ext_stats[tile][10] == 0) ext_stats[tile][10] = 4;
//					continue; // not a maximum
				}
				ext_stats[tile][7] = 1.0/area; // area-lma_max_area;
				ext_stats[tile][4] = 1.0/area_old;
			}
			double strength = Math.sqrt(avg/rrms);
			ext_stats[tile][9] = (strength * lma_str_scale) + lma_str_offset;
			ext_stats[tile][3] = Math.sqrt(strength / area ); // new strength
			double disparity = -all_pars[DISP_INDEX + tile*tile_params];
			if ((strength < lma_min_strength) || Double.isNaN(disparity)) {
				if (ext_stats[tile][10] == 0) ext_stats[tile][10] = 5;
				continue;
			}
			
			ext_stats[tile][8] = -all_pars[DISP_INDEX + tile*tile_params];
			if (ext_stats[tile][10] == 0) ext_stats[tile][10] = 6;
		}
		return ext_stats;
	}
//			ext_stats[tile][3] = Math.max(abc[tile][0], abc[tile][2])-lma_min_ac;

	public double [][] getDdNd(){ // this.all_pars should be current
		double [][] ddnd = new double [num_cams][2];
		for (int nc = 0; nc < num_cams; nc++) {
			ddnd[nc][0] = all_pars[ddisp_index + nc];
			ddnd[nc][1] = all_pars[ndisp_index + nc];
		}
		return ddnd;
	}

	public double [] getStats(){
		double [] stats = {total_weight/numTiles, 1.0*total_tiles/numTiles, last_rms[0]};
		return stats;
	}

	public double [] getStats(int num_good_tiles){
		double [] stats = {total_weight/num_good_tiles, 1.0*total_tiles/num_good_tiles, last_rms[0]};
		return stats;
	}

//

	public void getDdNd(double [][] ddnd){ // this.all_pars should be current
		if (ddnd != null) {
			for (int nc = 0; nc < num_cams; nc++) {
				ddnd[nc] = new double [2];
				ddnd[nc][0] = all_pars[ddisp_index + nc];
				ddnd[nc][1] = all_pars[ndisp_index + nc];
			}
		}
	}

	public void getStats(double [] stats){
		if ((stats != null) && (stats.length > 0)) {
			stats[0] = total_weight;
			if (stats.length > 1) {
				stats[1] = total_tiles;
				if (stats.length > 2) {
					stats[2] = last_rms[0];
				}
			}
		}
	}
	public void getStats(double [][] stats){
		double [] stats0 = {total_weight, total_tiles, last_rms[0]};
		stats[0] = stats0;
	}



	public double [] getRmsTile() { //
		double [] tile_rms =     new double[numTiles];
		double [] tile_weights = new double[numTiles];
		int num_samples = samples.size();
		for (int i = 0; i < num_samples; i++) {
			int tile =	samples.get(i).tile;
			double wd = last_ymfx[i];
			tile_rms[tile] += wd*wd/weights[i];
			tile_weights[tile] += weights[i];
		}
		for (int tile = 0; tile < numTiles; tile++) {
			if (tile_weights[tile] > 0) {
				tile_rms[tile] = Math.sqrt(tile_rms[tile]/tile_weights[tile]);
			} else {
				tile_rms[tile] = Double.NaN;
			}
		}
		return tile_rms;
	}


	/**
	 * Get maximal/minimal (per pair) amplitudes and values for each tile
	 * @return array [tile][{amax, amin}]
	 */
	public double [][] getMaxMinAmpTile(){
		return getMaxMinAmpTile(0);
	}
	public double [][] getMaxMinAmpTile(int nmax){
		double [][] maxmins = new double[numTiles][2];
		for (int tile = 0; tile < numTiles; tile++) {
			int offs = (tile * numMax + nmax) * tile_params;
			maxmins[tile][0] = Double.NaN;
			maxmins[tile][1] = Double.NaN;
			for (int i = 0; i <num_pairs; i++) {
				if (par_mask[G0_INDEX + i + offs]) {
					if (!(all_pars[G0_INDEX + i + offs] <= maxmins[tile][0])) { // NaN will trigger true
						maxmins[tile][0] = all_pars[G0_INDEX + i + offs];
					}
					if (!(all_pars[G0_INDEX + i + offs] >= maxmins[tile][1])) {
						maxmins[tile][1] = all_pars[G0_INDEX + i + offs];
					}
				}
			}
		}
		return maxmins;
	}
	/**
	 * Get maximal/minimal (per pair) amplitudes and values for each tile
	 * @return array [tile][{vmax, vmin}]
	 */
	public double [][] getMaxMinValTile(){
		double [][] tile_pair_max = new double[numTiles][num_pairs];
		int num_samples = samples.size();
		for (int i = 0; i < num_samples; i++) {
			Sample s = samples.get(i);
			int tile =	s.tile;
//			int pair =  pindx[s.fcam][s.scam];
			int pair =  s.pair;
			if (s.v > tile_pair_max[tile][pair]) tile_pair_max[tile][pair] = s.v;
		}

		double [][] maxmins = new double[numTiles][2];
		for (int tile = 0; tile < numTiles; tile++) {
			maxmins[tile][0] = Double.NaN;
			maxmins[tile][1] = Double.NaN;
			for (int i = 0; i <num_pairs; i++) {
				if (tile_pair_max[tile][i] > 0) {
					if (!(tile_pair_max[tile][i] <= maxmins[tile][0])) { // NaN will trigger true
						maxmins[tile][0] = tile_pair_max[tile][i];
					}
					if (!(tile_pair_max[tile][i] >= maxmins[tile][1])) {
						maxmins[tile][1] = tile_pair_max[tile][i];
					}
				}
			}
		}
		return maxmins;
	}

	/*public double [][] getABCTile(){
		double [][] abc = new double[numTiles][3];
		for (int tile = 0; tile < numTiles; tile++) {
			abc[tile][0] =  all_pars[A_INDEX+ tile * tile_params];
			abc[tile][1] =  all_pars[B_INDEX+ tile * tile_params];
			abc[tile][2] =  abc[tile][0] + all_pars[CMA_INDEX+ tile * tile_params];
		}
		return abc;
	}*/

	public double [][] getABCTile(int nmax){
		double [][] abc = new double[numTiles][3];
		for (int tile = 0; tile < numTiles; tile++) {
			int offs = (tile * numMax + nmax) * tile_params;
			abc[tile][0] =  all_pars[A_INDEX+ offs];
			abc[tile][1] =  all_pars[B_INDEX+ offs];
			abc[tile][2] =  abc[tile][0] + all_pars[CMA_INDEX+ offs];
		}
		return abc;
	}
	
	public double [][] getABCTile(){
		return getABCTile(0);
	}	
	
	@Deprecated
	public static double [] getJtWdiff( // not used in lwir
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
		this.last_rms = null;
		for (iter = 0; iter < num_iter; iter++) {
			rslt =  lmaStep(
					lambda,
					rms_diff,
					debug_level);
			if (rslt == null) {
				return false; // need to check
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
			System.out.println("Corr2dLMA:LMA: full RMS="+last_rms[0]+" ("+initial_rms[0]+"), pure RMS="+last_rms[1]+" ("+initial_rms[1]+") + lambda="+lambda);
		}

		return rslt[0];
	}


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
			if (last_ymfx == null) {
				return null; // need to re-init/restart LMA
			}
			this.last_rms = getWYmFxRms(this.last_ymfx); // modifies this.last_ymfx
			this.initial_rms = this.last_rms.clone();
			this.good_or_bad_rms = this.last_rms.clone();
			if (debug_level > 3) {
				debugJt(
						0.000001, // double      delta, // 0.2, //
						this.vector); // double []   vector);
			}
		}
		Matrix y_minus_fx_weighted = new Matrix(this.last_ymfx, this.last_ymfx.length);
		double [][] awjtjlambda = getWJtJlambda( // null
				lambda, // *10, // temporary
				this.last_jt); // double [][] jt)
		if ((awjtjlambda ==  null) || (awjtjlambda.length==0)) {
			return null;
		}
		Matrix wjtjlambda = new Matrix(awjtjlambda); // double [][] jt)
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
			System.out.println("(JtJ + lambda*diag(JtJ)).inv()");
			jtjl_inv.print(18, 6);
		}
//last_jt has NaNs
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
		if (last_ymfx == null) {
			return null; // need to re-init/restart LMA
		}
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
			if (last_ymfx == null) {
				return null; // need to re-init/restart LMA
			}
			this.last_rms = getWYmFxRms(this.last_ymfx); // modifies this.last_ymfx
			if (debug_level > 2) {
				debugJt(
						0.000001, // double      delta,
						this.vector); // double []   vector);
			}


		}
		return rslt;
	}
}
