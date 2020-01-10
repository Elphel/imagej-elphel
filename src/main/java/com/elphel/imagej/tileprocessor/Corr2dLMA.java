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

import Jama.Matrix;



public class Corr2dLMA {

	final static int NUM_CAMS =     4; // not all have to be used, so it is maximal number of cameras
	final static int NUM_PAIRS =    NUM_CAMS* (NUM_CAMS -1)/2; // number of possible pairs

	final static int DISP_INDEX =   0; // common/average disparity
	final static int A_INDEX =      1; // A*(x-x0)^2
	final static int B_INDEX =      2; // 2*B*(x-x0)*(y-y0)
	final static int CMA_INDEX =    3; // C*(y-y0)^2, encode C-A
	final static int G0_INDEX =     4; // scale of correlation pair,
	final static int TILE_PARAMS =  G0_INDEX +    NUM_PAIRS; // number of tile-individual parameters
//	final static int DDISP_INDEX =  G0_INDEX +    NUM_PAIRS; // disparity offset per camera (at least 1 should be disabled)
//	final static int NDISP_INDEX =  DDISP_INDEX + NUM_CAMS;  // disparity offset per camera - none should be disable
//	final static int NUM_ALL_PARS = NDISP_INDEX+ NUM_CAMS; // maximal number of parameters

	private int DDISP_INDEX;   // =  G0_INDEX +    NUM_PAIRS; // disparity offset per camera (at least 1 should be disabled)
	private int NDISP_INDEX;   //  =  DDISP_INDEX + NUM_CAMS;  // disparity offset per camera - none should be disable
	private int NUM_ALL_PARS;  //  = NDISP_INDEX+ NUM_CAMS; // maximal number of parameters


	final int []     USED_CAMS_MAP =  new int[NUM_CAMS]; // for each camera index return used index ???
//	final int [][]   USED_PAIRS_MAP = new int[NUM_CAMS][NUM_CAMS]; // for each camera index return used index ??
	private int [][][] USED_PAIRS_MAP; // for each camera index return used index ??

	final static String [] PAR_NAMES = {"DISP","A","B","C-A"};
	final static String PAR_NAME_SCALE = "SCALE";
	final static String PAR_NAME_CORRDISP = "CORR-DISP";
	final static String PAR_NAME_CORRNDISP = "CORR-NDISP";
	double  [] all_pars;
	boolean [] par_mask;
	int     [] par_map;
	double  [] vector;
	double  [] scales = {1.0, 2.0, 4.0};
	ArrayList<Sample> samples = new ArrayList<Sample>();
////	double [] pair_weights = null; // per pair weights (sum == 1.0) Not really needed?
	double [] weights; // normalized so sum is 1.0 for all - samples and extra regularization terms
	double    pure_weight; // weight of samples only
	double [] values;
	// next values are only updated after success
	double []   last_rms =        null; // {rms, rms_pure}, matching this.vector
	double []   good_or_bad_rms = null; // just for diagnostics, to read last (failed) rms
	double []   initial_rms =     null; // {rms, rms_pure}, first-calcualted rms
	double []   last_ymfx =       null;
	double [][] last_jt =         null;
	double []   poly_coeff =      null;  // 6 elements - Xc, Yx, f(x,y), A, B, C (from A*x^2 + B*y^2 +C*x*y+...)
	double []   poly_xyvwh =      null;  // result of 2-d polynomial approximation instead of the LMA - used for lazy eye correction
	private final int transform_size;
    private final double [][] corr_wnd;
    private boolean [] used_cameras;
//    private final Matrix [] m_disp = new Matrix[NUM_CAMS];
    private Matrix [][] m_disp;
    private int ncam; // number of used cameras
    private int [] npairs; // number of used pairs per tile
    private int     last_cam; // index of the last camera (special treatment for disparity correction)
//    private boolean second_last; // there is a pair where the second camera is the last one (false: first in a pair is the last one)
//    private final Matrix [][] m_pairs =      new Matrix[NUM_CAMS][NUM_CAMS];
    private Matrix [][][] m_pairs;
//    private final Matrix [][] m_pairs_last = new Matrix[NUM_CAMS][NUM_CAMS];
    private final int [][] pindx =           new int [NUM_CAMS][NUM_CAMS];
    private int            numTiles = 1;

    public boolean gaussian_mode = true;

	public class Sample{ // USED in lwir
		int    tile;   // tile in a cluster
		int    fcam;   // first  camera index
		int    scam;   // second camera index
		int    ix;     // x coordinate in 2D correlation (0.. 2*transform_size-2, center: (transform_size-1)
		int    iy;     // y coordinate in 2D correlation (0.. 2*transform_size-2, center: (transform_size-1)
		double v;      // correlation value at that point
		double w;      // weight
		Sample (
				int    tile,
				int    fcam,   // first  camera index
				int    scam,   // second camera index
				int    x,      // x coordinate on the common scale (corresponding to the largest baseline), along the disparity axis
				int    y,      //  coordinate in 2D correlation (0.. 2*transform_size-2, center: (transform_size-1)
				double v,      // correlation value at that point
				double w)
			{
			    this.tile = tile;
				this.fcam = fcam;
				this.scam = scam;
				this.ix =  x;
				this.iy =  y;
				this.v =   v;
				this.w =   w;
			}
	}


	public Corr2dLMA (
				int numTiles,
				int ts, // null - use default table
				double [][] corr_wnd, // may be null
				boolean gaussian_mode
				) {
		this.gaussian_mode = gaussian_mode;
		for (int f = 0; f < NUM_CAMS; f++) {
			pindx[f][f]=-1;
			for (int s = f+1; s < NUM_CAMS; s++) {
				pindx[f][s] = getPairIndex(f,s);
				pindx[s][f] = pindx[f][s];
			}
		}
		this.numTiles = numTiles;
		DDISP_INDEX = this.numTiles * TILE_PARAMS;
		NDISP_INDEX =  DDISP_INDEX + NUM_CAMS;  // disparity offset per camera - none should be disable
		NUM_ALL_PARS = NDISP_INDEX+ NUM_CAMS; // maximal number of parameters

		boolean sq = false;
		this.transform_size = ts;
		if (corr_wnd != null) {
		    this.corr_wnd = corr_wnd;
		} else {
			this.corr_wnd = new double[2 * transform_size - 1][2 * transform_size - 1];
			int tsm1 =  transform_size - 1; // 7
			int dtsm1 = 2 * transform_size - 1; // 15
			this.corr_wnd[tsm1][tsm1] = 1.0;
			for (int i = 1; i < transform_size; i++) {
				this.corr_wnd[tsm1 + i][tsm1    ] = Math.cos(Math.PI*i/(2 * transform_size));
				this.corr_wnd[tsm1 - i][tsm1    ] = Math.cos(Math.PI*i/(2 * transform_size));
				this.corr_wnd[tsm1    ][tsm1 + i] = Math.cos(Math.PI*i/(2 * transform_size));
				this.corr_wnd[tsm1    ][tsm1 - i] = Math.cos(Math.PI*i/(2 * transform_size));
			}
			for (int i = 1; i < transform_size; i++) {
				for (int j = 1; j < transform_size; j++) {
					double d = this.corr_wnd[tsm1 + i][tsm1] * this.corr_wnd[tsm1 + j][tsm1];
					this.corr_wnd[tsm1 + i][tsm1 + j] = d;
					this.corr_wnd[tsm1 + i][tsm1 - j] = d;
					this.corr_wnd[tsm1 - i][tsm1 + j] = d;
					this.corr_wnd[tsm1 - i][tsm1 - j] = d;
				}
			}
			if (sq) {
				for (int i = 0; i < dtsm1; i++) {
					for (int j = 0; j < dtsm1; j++) {
						this.corr_wnd[i][j] *=this.corr_wnd[i][j];
					}
				}
			}
		}
	}

	public double[][] getCorrWnd() {
		return this.corr_wnd;
	}

	public void addSample( // x = 0, y=0 - center
			int    tile,
			int    fcam,   // first  camera index
			int    scam,   // second camera index
			int    x,      // x coordinate on the common scale (corresponding to the largest baseline), along the disparity axis
			int    y,      // y coordinate (0 - disparity axis)
			double v,      // correlation value at that point
			double w){     // sample weight
		if ((w > 0) && !Double.isNaN(v)) samples.add(new Sample(tile,fcam,scam,x,y,v,w));
	}

	public double [][][] dbgGetSamples(int mode){
		int [][] comb_map = getCombMap();
		int numPairs = comb_map[0][0];
		comb_map[0][0] = -1;

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
			double d = Double.NaN;
			if      (mode == 0) d = s.v;
			else if (mode == 1) d = s.w;
			else if (mode == 2) d = fx[ns];
//			int np = USED_PAIRS_MAP[0][s.fcam][s.scam]; ////////////////////
			int np = comb_map[s.fcam][s.scam]; ////////////////////
			rslt[s.tile][np][s.iy*size + s.ix] = d;
		}

		return rslt;
	}
	/*
	public String [] dbgGetSliceTiles(int ntile) {
		String [] srslt = new String [npairs[ntile]];
		for (int f = 0; f < NUM_CAMS; f++) for (int s = 0; s < NUM_CAMS; s++) {
			if (USED_PAIRS_MAP[ntile][f][s] >= 0) {
				srslt[USED_PAIRS_MAP[ntile][f][s]] = ""+f+"->"+s;
			}
		}
		return srslt;
	}
	*/
	public String [] dbgGetSliceTiles() {
		int [][] comb_map = getCombMap();
		int np = comb_map[0][0];
		comb_map[0][0] = -1;

		String [] srslt = new String [np];
		for (int f = 0; f < NUM_CAMS; f++) for (int s = 0; s < NUM_CAMS; s++) {
			if (comb_map[f][s] >= 0) {
				srslt[comb_map[f][s]] = ""+f+"->"+s;
			}
		}
		return srslt;
	}

	public int [][] getCombMap(){
		boolean [][]  comb_pairs = new boolean[NUM_CAMS][NUM_CAMS];

		for (int t = 0; t < numTiles; t++) {
			for (int f = 0; f < NUM_CAMS; f++) for (int s = 0; s < NUM_CAMS; s++) {
				comb_pairs[f][s] |= USED_PAIRS_MAP[t][f][s] >= 0;
			}
		}
		int np = 0;
		int [][] comb_map = new int [NUM_CAMS][NUM_CAMS];
		for (int f = 0; f < NUM_CAMS; f++) for (int s = 0; s < NUM_CAMS; s++) {
			if  (comb_pairs[f][s]) comb_map[f][s] = np++;
			else comb_map[f][s] = -1;
		}
		comb_map[0][0] = np;
		return comb_map;

	}



	public int getPairIndex(int f, int s) {
		if (f > s) {
			int t = f;
			f = s;
			s = t;
		}
		return (NUM_CAMS * f) - (f + 1)*f/2 - f - 1 + s ; // return n*i - i*(i+1)//2 - i + j -1
	}

	public void setMatrices(double [][] am_disp) {
		m_disp = new Matrix[1][NUM_CAMS];
		for (int n = 0; n < NUM_CAMS; n++) {
			double [][] am = {
					{am_disp[n][0],am_disp[n][1]},
					{am_disp[n][2],am_disp[n][3]}};
			m_disp[0][n] = new Matrix(am);
		}
	}

	public void setMatrices(double [][][] am_disp) {
		m_disp = new Matrix[am_disp.length][NUM_CAMS];
		for (int nt = 0; nt < numTiles; nt++) {
			for (int n = 0; n < NUM_CAMS; n++) {
				double [][] am = {
						{am_disp[nt][n][0], am_disp[nt][n][1]},
						{am_disp[nt][n][2], am_disp[nt][n][3]}};
				m_disp[nt][n] = new Matrix(am);
			}
		}
	}

	public void initVector( // USED in lwir
			boolean adjust_width,         // adjust width of the maximum -                                           lma_adjust_wm
			boolean adjust_scales,        // adjust 2D correlation scales -                                          lma_adjust_ag
			boolean adjust_ellipse,       // allow non-circular correlation maximums                                 lma_adjust_wy
			boolean adjust_lazyeye_par,   // adjust disparity corrections parallel to disparities                    lma_adjust_wxy
			boolean adjust_lazyeye_ortho, // adjust disparity corrections orthogonal to disparities                  lma_adjust_ly1
			double [][] disp_str,         // initial value of disparity
//			double  disp0,                // initial value of disparity
			double  half_width,           // A=1/(half_widh)^2                                                       lma_half_width
			double  cost_lazyeye_par,     // cost for each of the non-zero disparity corrections                     lma_cost_wy
			double  cost_lazyeye_odtho    // cost for each of the non-zero ortho disparity corrections               lma_cost_wxy
			) {
//		int [][] pindx = new int [NUM_CAMS][NUM_CAMS];
		USED_PAIRS_MAP = new int [numTiles][NUM_CAMS][NUM_CAMS];
		used_cameras = new boolean[NUM_CAMS];
		boolean [][] used_pairs = new boolean[numTiles][NUM_PAIRS];
		// 0-weight values and NaN-s should be filtered on input!
		last_cam = -1;
		for (int t = 0; t < numTiles; t++) for (int f = 0; f < NUM_CAMS; f++) for (int s = 0; s < NUM_CAMS; s++) {
			USED_PAIRS_MAP[t][f][s] = -1;
		}
		boolean [][][] used_pairs_dir = new boolean [numTiles][NUM_CAMS][NUM_CAMS];
		for (Sample s:samples) { // ignore zero-weight samples
			used_cameras[s.fcam]=true;
			used_cameras[s.scam]=true;
			used_pairs[s.tile][pindx[s.fcam][s.scam]]=true; // throws < 0 - wrong pair, f==s
			used_pairs_dir[s.tile][s.fcam][s.scam] = true;
		}
		ncam = 0;
		npairs =new int [numTiles];
		for (int i = 0; i < NUM_CAMS; i++) {
			USED_CAMS_MAP[i] = ncam;
			if (used_cameras[i]) {
				last_cam = ncam;
				ncam++;
			}
		}
		for (int nTile = 0; nTile < numTiles; nTile++) {
			int [] upmam = new int[NUM_PAIRS];
			for (int i = 0; i < NUM_PAIRS; i++) {
				upmam[i] = npairs[nTile];
				if (used_pairs[nTile][i])  npairs[nTile]++;
			}
			for (int f = 0; f < NUM_CAMS; f++) {
				for (int s = f+1; s < NUM_CAMS; s++) {
					int npair = upmam[pindx[f][s]];
					if      (used_pairs_dir[nTile][f][s]) USED_PAIRS_MAP[nTile][f][s] = npair;  // either or, can not be f,s and s,f pairs
					else if (used_pairs_dir[nTile][s][f]) USED_PAIRS_MAP[nTile][s][f] = npair;
				}
			}
		}

		this.all_pars = new double[NUM_ALL_PARS];
		this.par_mask = new boolean[NUM_ALL_PARS];
		// per-tile parameters
		for (int nTile = 0; nTile < numTiles; nTile++) {
			this.all_pars[DISP_INDEX + nTile*TILE_PARAMS] = disp_str[nTile][0]; // disp0;
			this.all_pars[A_INDEX    + nTile*TILE_PARAMS] = 1.0/(half_width * half_width);
			this.all_pars[B_INDEX    + nTile*TILE_PARAMS] = 0.0;
			this.all_pars[CMA_INDEX  + nTile*TILE_PARAMS] = 0.0; // C-A
			this.par_mask[DISP_INDEX + nTile*TILE_PARAMS] = true;
			this.par_mask[A_INDEX    + nTile*TILE_PARAMS] = adjust_width;
			this.par_mask[B_INDEX    + nTile*TILE_PARAMS] = adjust_ellipse;
			this.par_mask[CMA_INDEX  + nTile*TILE_PARAMS] = adjust_ellipse;
			for (int i = 0; i <NUM_PAIRS; i++) {
				this.par_mask[G0_INDEX + i + nTile*TILE_PARAMS] = used_pairs[nTile][i] & adjust_scales;
				this.all_pars[G0_INDEX + i + nTile*TILE_PARAMS] = Double.NaN; // will be assigned later for used - should be for all !
			}
		}
		// common for all tiles parameters
		for (int i = 0; i <NUM_CAMS; i++) {
			this.all_pars[DDISP_INDEX + i] = 0.0; // C-A
			this.par_mask[DDISP_INDEX + i] = used_cameras[i] & adjust_lazyeye_par & (i != last_cam);
			this.all_pars[NDISP_INDEX + i] = 0.0; // C-A
			this.par_mask[NDISP_INDEX + i] = used_cameras[i] & adjust_lazyeye_ortho;
		}
		int np = samples.size();

		weights = new double [np + 2 * NUM_CAMS]; // npairs];
		values =  new double [np + 2 * NUM_CAMS]; // npairs];
		for (int i = 0; i < NUM_CAMS; i++) {
			weights[np + i] =            (used_cameras[i] & adjust_lazyeye_par)?   (cost_lazyeye_par * numTiles) :   0.0; // ddisp - including last_camera
			weights[np + NUM_CAMS + i] = (used_cameras[i] & adjust_lazyeye_ortho)? (cost_lazyeye_odtho * numTiles) : 0.0; // ndisp
			values [np + i] =            0.0;
			values [np + NUM_CAMS + i] = 0.0;
		}

		double sw = 0;
////		this.pair_weights = new double[NUM_PAIRS];
		for (int i = 0; i < np; i++) {
			Sample s = samples.get(i);
			weights[i] = s.w;
			values[i] =  s.v;
			sw += weights[i];
			int indx = G0_INDEX + pindx[s.fcam][s.scam] + s.tile * TILE_PARAMS;
			double d = s.v;
			if (this.corr_wnd !=null) {
				d /= this.corr_wnd[s.iy][s.ix];
			}
			if (!(d <= this.all_pars[indx])) this.all_pars[indx] = d; // to include Double.isNaN()
		}
		pure_weight = sw;
		for (int i = 0; i < 2 * NUM_CAMS; i++) { // weight of the regularization terms (twice number of cameras, some may be disabled by a mask)
			sw += weights[np + i];
		}
		if (sw != 0.0) {
			double kw = 1.0/sw;
			for (int i = 0; i < weights.length; i++) weights[i] *= kw;
			pure_weight *= kw; // it is now fraction (0..1.0), and weights are normalized
		}
		/*****
		double spw = 0;
		for (int i = 0; i < NUM_PAIRS; i++) {
			spw += pair_weights[i];
		}
		if (spw > 0) {
			double rspw = 1.0/spw;
			for (int i = 0; i < NUM_PAIRS; i++) {
				pair_weights[i]*=rspw;
			}
		}
		*/
		par_map = new int [par_mask.length];
		int par_indx = 0;
		for (int i = 0; i < par_mask.length; i++) {
			if (par_mask[i]) par_map[i] = par_indx++;
			else par_map[i] = -1;
		}

		toVector();
	}



	public void initMatrices() { // should be called after initVector and after setMatrices
		m_pairs = new Matrix[USED_PAIRS_MAP.length][NUM_CAMS][NUM_CAMS];
		for (int nTile = 0; nTile < USED_PAIRS_MAP.length; nTile++) {
			for (int f = 0; f < NUM_CAMS; f++) for (int s = 0; s < NUM_CAMS; s++) {
				m_pairs[nTile][f][s] =      null;
				//			m_pairs_last[f][s] = null;
				if (USED_PAIRS_MAP[nTile][f][s] >= 0) {
					m_pairs[nTile][f][s] = m_disp[nTile][f].minus(m_disp[nTile][s]);
				}
				/*
			if (f == last_cam) {
				m_pairs_last[f][s] = m_disp[s].uminus();
				for (int i = 0; i < NUM_CAMS; i++) if (used_cameras[i] && (i != last_cam) ){
					m_pairs_last[f][s].minusEquals(m_disp[i]);
				}
			} else if (s == last_cam) {
				m_pairs_last[f][s] = m_disp[f].copy();
				for (int i = 0; i < NUM_CAMS; i++) if (used_cameras[i] && (i != last_cam) ){
					m_pairs_last[f][s].plusEquals(m_disp[i]);
				}
			}
				 */
			}
		}
	}
	public double [] getFxJt( // USED in lwir
			double []   vector,
			double [][] jt) { // should be either [vector.length][samples.size()] or null - then only fx is calculated
		if (this.gaussian_mode) return getFxJt_gaussian(vector, jt);
		else                    return getFxJt_parabola(vector, jt);
	}

	public double [] getFxJt_parabola( // USED in lwir
			double []   vector,
			double [][] jt) { // should be either [vector.length][samples.size()] or null - then only fx is calculated
		if (vector == null) return null;
		double [] av = fromVector(vector);
		Matrix [][] xcam_ycam = new Matrix[numTiles][NUM_CAMS];
		double [][][][] xp_yp = new double[numTiles][NUM_CAMS][NUM_CAMS][];
		double [] axc_yc = {transform_size - 1.0, transform_size-1.0};
		Matrix xc_yc = new Matrix(axc_yc, 2);
		double [] AT = new double [numTiles]; // av[A_INDEX];
		double [] BT = new double [numTiles]; // av[B_INDEX];
		double [] CT = new double [numTiles]; // A + av[CMA_INDEX];
		for (int nTile = 0; nTile < numTiles; nTile++) {
			for (int i = 0; i < NUM_CAMS; i++) if (used_cameras[i]) {
				double [] add_dnd = {av[DISP_INDEX+ nTile * TILE_PARAMS]+ av[DDISP_INDEX + i],  av[NDISP_INDEX + i]};
				xcam_ycam[nTile][i] = m_disp[nTile][i].times(new Matrix(add_dnd,2));
			}
			for (int f = 0; f < NUM_CAMS; f++) if (used_cameras[f]) {
				for (int s = 0; s < NUM_CAMS; s++) if (used_cameras[s]) {
					xp_yp[nTile][f][s] =xcam_ycam[nTile][f].minus(xcam_ycam[nTile][s]).plus(xc_yc).getColumnPackedCopy();
				}
			}
			AT[nTile] = av[A_INDEX + nTile * TILE_PARAMS];
			BT[nTile] = av[B_INDEX + nTile * TILE_PARAMS];
			CT[nTile] = AT[nTile] + av[CMA_INDEX + nTile * TILE_PARAMS];
		}

		int num_samples = samples.size();
		double [] fx= new double [num_samples + 2 * NUM_CAMS];
//		double A = av[A_INDEX];
//		double B = av[B_INDEX];
//		double C = A + av[CMA_INDEX];
//corr_wnd
		for (int ns = 0; ns < num_samples; ns++) {
			Sample s = samples.get(ns);
			int pair = pindx[s.fcam][s.scam]; // all pairs, noit just used?
			double A = AT[s.tile];
			double B = BT[s.tile];
			double C = CT[s.tile];

			double Gp = av[G0_INDEX + pair + s.tile * TILE_PARAMS];
			double Wp = corr_wnd[s.ix][s.iy];
			double WGp = Wp * Gp;
			double xmxp = s.ix - xp_yp[s.tile][s.fcam][s.scam][0];
			double ymyp = s.iy - xp_yp[s.tile][s.fcam][s.scam][1];
			double xmxp2 = xmxp * xmxp;
			double ymyp2 = ymyp * ymyp;
			double xmxp_ymyp = xmxp * ymyp;
			double d = Wp*(1.0 - (A*xmxp2 + 2 * B * xmxp_ymyp + C * ymyp2));
			fx[ns] = d * Gp;
			if (Double.isNaN(fx[ns])) {
				System.out.println("fx["+ns+"]="+fx[ns]);
			}
			if (s.tile > 0) {
				System.out.print("");
			}
			if (jt != null) {
				if (par_map[DISP_INDEX + s.tile*TILE_PARAMS] >= 0)  jt[par_map[DISP_INDEX + s.tile*TILE_PARAMS]][ns] = 2 * WGp *
						((A * xmxp + B * ymyp) * m_pairs[s.tile][s.fcam][s.scam].get(0, 0)+
						 (B * xmxp + C * ymyp) * m_pairs[s.tile][s.fcam][s.scam].get(1, 0));
				if (par_map[A_INDEX + s.tile*TILE_PARAMS] >= 0)     jt[par_map[A_INDEX + s.tile*TILE_PARAMS]][ns] = -WGp*(xmxp2 + ymyp2);
				if (par_map[B_INDEX + s.tile*TILE_PARAMS] >= 0)     jt[par_map[B_INDEX + s.tile*TILE_PARAMS]][ns] = -WGp* 2 * xmxp_ymyp;
				if (par_map[CMA_INDEX + s.tile*TILE_PARAMS] >= 0)   jt[par_map[CMA_INDEX + s.tile*TILE_PARAMS]][ns] = -WGp* ymyp2;
				for (int p = 0; p < npairs[s.tile]; p++) { // par_mask[G0_INDEX + p] as all pairs either used, or not - then npairs == 0
					if (par_map[G0_INDEX + p + s.tile*TILE_PARAMS] >= 0) jt[par_map[G0_INDEX + p + s.tile*TILE_PARAMS]][ns] = (p== pair)? d : 0.0; // (par_mask[G0_INDEX + pair])? d;
				}
				// process ddisp (last camera not used, is equal to minus sum of others to make a sum == 0)

				for (int f = 0; f < NUM_CAMS; f++) if (par_map[DDISP_INDEX + f] >= 0) { // -1 for the last_cam
					jt[par_map[DDISP_INDEX + f]][ns] = 0.0;
				}
				if (par_map[DDISP_INDEX + s.fcam] >= 0){ // par_map[DDISP_INDEX + last_cam] always <0
						jt[par_map[DDISP_INDEX + s.fcam]][ns] += 2 * WGp *
									((A * xmxp + B * ymyp) * m_disp[s.tile][s.fcam].get(0, 0)+
									 (B * xmxp + C * ymyp) * m_disp[s.tile][s.fcam].get(1, 0));
				} else if (s.fcam == last_cam) {
					for (int c = 0; c < NUM_CAMS; c++) if ((c != last_cam) && (par_map[DDISP_INDEX + c] >=0)) {
						jt[par_map[DDISP_INDEX + c]][ns] -= 2 * WGp *
								(       (A * xmxp + B * ymyp) * m_disp[s.tile][s.fcam].get(0, 0)+
										(B * xmxp + C * ymyp) * m_disp[s.tile][s.fcam].get(1, 0));
					}
				}
				if (par_map[DDISP_INDEX + s.scam]>= 0){ // par_map[DDISP_INDEX + last_cam] always <0
					 jt[par_map[DDISP_INDEX + s.scam]][ns] -= 2 * WGp *
								((A * xmxp + B * ymyp) * m_disp[s.tile][s.scam].get(0, 0)+
								 (B * xmxp + C * ymyp) * m_disp[s.tile][s.scam].get(1, 0));

				} else if (s.scam == last_cam) {
					for (int c = 0; c < NUM_CAMS; c++) if ((c != last_cam) && (par_map[DDISP_INDEX + c] >= 0)) {
						 jt[par_map[DDISP_INDEX + c]][ns] += 2 * WGp *
							        ((A * xmxp + B * ymyp) * m_disp[s.tile][s.scam].get(0, 0)+
									 (B * xmxp + C * ymyp) * m_disp[s.tile][s.scam].get(1, 0));
					}
				}

				// process ndisp
				for (int f = 0; f < ncam; f++) if (par_map[NDISP_INDEX + f] >= 0) {
					jt[par_map[NDISP_INDEX + f]][ns] = 0.0;
				}
				if (par_map[NDISP_INDEX + s.fcam] >=0){
					jt[par_map[NDISP_INDEX + s.fcam]][ns] += 2 * WGp *
							(       (A * xmxp + B * ymyp) * m_disp[s.tile][s.fcam].get(0, 1)+
									(B * xmxp + C * ymyp) * m_disp[s.tile][s.fcam].get(1, 1));
				}

				if (par_map[NDISP_INDEX + s.scam] >= 0) {

					jt[par_map[NDISP_INDEX + s.scam]][ns] -= 2 * WGp *
							(       (A * xmxp + B * ymyp) * m_disp[s.tile][s.scam].get(0, 1)+
									(B * xmxp + C * ymyp) * m_disp[s.tile][s.scam].get(1, 1));
				}
			}
		}
		for (int n = 0; n < NUM_CAMS; n++) { // av[DDISP_INDEX +last_cam] is already populated
			fx[num_samples +            n] = av[DDISP_INDEX + n];
			fx[num_samples + NUM_CAMS + n] = av[NDISP_INDEX + n];
		}

// and derivatives
		if (jt != null) {
			for (int i = 0; i < NUM_CAMS; i++) {
				if ((i != last_cam) && (par_map[DDISP_INDEX + i] >= 0)) {
					for (int j = 0; j < NUM_CAMS; j++) { // j - column
						jt[par_map[DDISP_INDEX + i]][num_samples + j] = (i==j)? 1.0 : 0.0;
					}
					jt[par_map[DDISP_INDEX + i]][num_samples + last_cam] = -1.0;
				}
			}
			for (int i = 0; i < NUM_CAMS; i++) {
				if (par_map[NDISP_INDEX + i] >= 0) {
					for (int j = 0; j < NUM_CAMS; j++) { // j - column
						jt[par_map[NDISP_INDEX + i] ][num_samples + NUM_CAMS + j] = (i==j)? 1.0 : 0.0;
					}
				}
			}
		}
		return fx;
	}

	public double [] getFxJt_gaussian( // USED in lwir
			double []   vector,
			double [][] jt) { // should be either [vector.length][samples.size()] or null - then only fx is calculated
		if (vector == null) return null;
		double [] av = fromVector(vector);
		Matrix [][] xcam_ycam = new Matrix[numTiles][NUM_CAMS];
		double [][][][] xp_yp = new double[numTiles][NUM_CAMS][NUM_CAMS][];
		double [] axc_yc = {transform_size - 1.0, transform_size-1.0};
		Matrix xc_yc = new Matrix(axc_yc, 2);
		double [] AT = new double [numTiles]; // av[A_INDEX];
		double [] BT = new double [numTiles]; // av[B_INDEX];
		double [] CT = new double [numTiles]; // A + av[CMA_INDEX];
		for (int nTile = 0; nTile < numTiles; nTile++) {
			for (int i = 0; i < NUM_CAMS; i++) if (used_cameras[i]) {
				double [] add_dnd = {av[DISP_INDEX+ nTile * TILE_PARAMS]+ av[DDISP_INDEX + i],  av[NDISP_INDEX + i]};
				xcam_ycam[nTile][i] = m_disp[nTile][i].times(new Matrix(add_dnd,2));
			}
			for (int f = 0; f < NUM_CAMS; f++) if (used_cameras[f]) {
				for (int s = 0; s < NUM_CAMS; s++) if (used_cameras[s]) {
					xp_yp[nTile][f][s] =xcam_ycam[nTile][f].minus(xcam_ycam[nTile][s]).plus(xc_yc).getColumnPackedCopy();
				}
			}
			AT[nTile] = av[A_INDEX + nTile * TILE_PARAMS];
			BT[nTile] = av[B_INDEX + nTile * TILE_PARAMS];
			CT[nTile] = AT[nTile] + av[CMA_INDEX + nTile * TILE_PARAMS];
		}

		int num_samples = samples.size();
		double [] fx= new double [num_samples + 2 * NUM_CAMS];
//corr_wnd
		for (int ns = 0; ns < num_samples; ns++) {
			Sample s = samples.get(ns);
			int pair = pindx[s.fcam][s.scam]; // all pairs, noit just used?
			double A = AT[s.tile];
			double B = BT[s.tile];
			double C = CT[s.tile];

			double Gp = av[G0_INDEX + pair + s.tile * TILE_PARAMS];
			double Wp = corr_wnd[s.ix][s.iy];
			double WGp = Wp * Gp;
			double xmxp = s.ix - xp_yp[s.tile][s.fcam][s.scam][0];
			double ymyp = s.iy - xp_yp[s.tile][s.fcam][s.scam][1];
			double xmxp2 = xmxp * xmxp;
			double ymyp2 = ymyp * ymyp;
			double xmxp_ymyp = xmxp * ymyp;
////			double comm = Wp*(1.0 - (A*xmxp2 + 2 * B * xmxp_ymyp + C * ymyp2));
			double exp = Math.exp(-(A*xmxp2 + 2 * B * xmxp_ymyp + C * ymyp2));
			double comm = exp * Wp;
			double WGpexp = WGp*exp;

			fx[ns] = comm * Gp;
			if (Double.isNaN(fx[ns])) {
				System.out.println("fx["+ns+"]="+fx[ns]);
			}
			if (s.tile > 0) {
				System.out.print("");
			}
			if (jt != null) {
				if (par_map[DISP_INDEX + s.tile*TILE_PARAMS] >= 0)  jt[par_map[DISP_INDEX + s.tile*TILE_PARAMS]][ns] = 2 * WGpexp *
						((A * xmxp + B * ymyp) * m_pairs[s.tile][s.fcam][s.scam].get(0, 0)+
						 (B * xmxp + C * ymyp) * m_pairs[s.tile][s.fcam][s.scam].get(1, 0));


				if (par_map[A_INDEX + s.tile*TILE_PARAMS] >= 0)     jt[par_map[A_INDEX + s.tile*TILE_PARAMS]][ns] = -WGpexp*(xmxp2 + ymyp2);
				if (par_map[B_INDEX + s.tile*TILE_PARAMS] >= 0)     jt[par_map[B_INDEX + s.tile*TILE_PARAMS]][ns] = -WGpexp* 2 * xmxp_ymyp;
				if (par_map[CMA_INDEX + s.tile*TILE_PARAMS] >= 0)   jt[par_map[CMA_INDEX + s.tile*TILE_PARAMS]][ns] = -WGp* ymyp2 * exp;
				for (int p = 0; p < npairs[s.tile]; p++) { // par_mask[G0_INDEX + p] as all pairs either used, or not - then npairs == 0
					if (par_map[G0_INDEX + p + s.tile*TILE_PARAMS] >= 0) jt[par_map[G0_INDEX + p + s.tile*TILE_PARAMS]][ns] = (p== pair)? comm : 0.0; // (par_mask[G0_INDEX + pair])? d;
				}
				// process ddisp (last camera not used, is equal to minus sum of others to make a sum == 0)

				for (int f = 0; f < NUM_CAMS; f++) if (par_map[DDISP_INDEX + f] >= 0) { // -1 for the last_cam
					jt[par_map[DDISP_INDEX + f]][ns] = 0.0;
				}
				if (par_map[DDISP_INDEX + s.fcam] >= 0){ // par_map[DDISP_INDEX + last_cam] always <0
						jt[par_map[DDISP_INDEX + s.fcam]][ns] += 2 * WGpexp *
									((A * xmxp + B * ymyp) * m_disp[s.tile][s.fcam].get(0, 0)+
									 (B * xmxp + C * ymyp) * m_disp[s.tile][s.fcam].get(1, 0));
				} else if (s.fcam == last_cam) {
					for (int c = 0; c < NUM_CAMS; c++) if ((c != last_cam) && (par_map[DDISP_INDEX + c] >=0)) {
						jt[par_map[DDISP_INDEX + c]][ns] -= 2 * WGpexp *
								(       (A * xmxp + B * ymyp) * m_disp[s.tile][s.fcam].get(0, 0)+
										(B * xmxp + C * ymyp) * m_disp[s.tile][s.fcam].get(1, 0));
					}
				}
				if (par_map[DDISP_INDEX + s.scam]>= 0){ // par_map[DDISP_INDEX + last_cam] always <0
					 jt[par_map[DDISP_INDEX + s.scam]][ns] -= 2 * WGpexp *
								((A * xmxp + B * ymyp) * m_disp[s.tile][s.scam].get(0, 0)+
								 (B * xmxp + C * ymyp) * m_disp[s.tile][s.scam].get(1, 0));

				} else if (s.scam == last_cam) {
					for (int c = 0; c < NUM_CAMS; c++) if ((c != last_cam) && (par_map[DDISP_INDEX + c] >= 0)) {
						 jt[par_map[DDISP_INDEX + c]][ns] += 2 * WGpexp *
							        ((A * xmxp + B * ymyp) * m_disp[s.tile][s.scam].get(0, 0)+
									 (B * xmxp + C * ymyp) * m_disp[s.tile][s.scam].get(1, 0));
					}
				}

				// process ndisp
				for (int f = 0; f < ncam; f++) if (par_map[NDISP_INDEX + f] >= 0) {
					jt[par_map[NDISP_INDEX + f]][ns] = 0.0;
				}
				if (par_map[NDISP_INDEX + s.fcam] >=0){
					jt[par_map[NDISP_INDEX + s.fcam]][ns] += 2 * WGpexp *
							(       (A * xmxp + B * ymyp) * m_disp[s.tile][s.fcam].get(0, 1)+
									(B * xmxp + C * ymyp) * m_disp[s.tile][s.fcam].get(1, 1));
				}

				if (par_map[NDISP_INDEX + s.scam] >= 0) {

					jt[par_map[NDISP_INDEX + s.scam]][ns] -= 2 * WGpexp *
							(       (A * xmxp + B * ymyp) * m_disp[s.tile][s.scam].get(0, 1)+
									(B * xmxp + C * ymyp) * m_disp[s.tile][s.scam].get(1, 1));
				}
			}
		}
		for (int n = 0; n < NUM_CAMS; n++) { // av[DDISP_INDEX +last_cam] is already populated
			fx[num_samples +            n] = av[DDISP_INDEX + n];
			fx[num_samples + NUM_CAMS + n] = av[NDISP_INDEX + n];
		}

// and derivatives
		if (jt != null) {
			for (int i = 0; i < NUM_CAMS; i++) {
				if ((i != last_cam) && (par_map[DDISP_INDEX + i] >= 0)) {
					for (int j = 0; j < NUM_CAMS; j++) { // j - column
						jt[par_map[DDISP_INDEX + i]][num_samples + j] = (i==j)? 1.0 : 0.0;
					}
					jt[par_map[DDISP_INDEX + i]][num_samples + last_cam] = -1.0;
				}
			}
			for (int i = 0; i < NUM_CAMS; i++) {
				if (par_map[NDISP_INDEX + i] >= 0) {
					for (int j = 0; j < NUM_CAMS; j++) { // j - column
						jt[par_map[NDISP_INDEX + i] ][num_samples + NUM_CAMS + j] = (i==j)? 1.0 : 0.0;
					}
				}
			}
		}
		return fx;
	}

	public void printParams() { // not used in lwir
		for (int np = 0; np < all_pars.length; np++) {
			String parname;
//			if      (np < G0_INDEX)    parname = PAR_NAMES[np];
//			else if (np < DDISP_INDEX) parname = PAR_NAME_SCALE;
//			else if (np < NDISP_INDEX) parname = PAR_NAME_CORRDISP;
//			else                       parname = PAR_NAME_CORRNDISP;
			if      (np >= NDISP_INDEX) parname = PAR_NAME_CORRNDISP + (np - NDISP_INDEX);
			else if (np >= DDISP_INDEX) parname = PAR_NAME_CORRDISP +  (np - DDISP_INDEX);
			else {
				int ntile = np / TILE_PARAMS;
				int anpr =  np % TILE_PARAMS;
				if (anpr < G0_INDEX) parname = PAR_NAMES[anpr]+"-"+ntile;
				else                 parname = PAR_NAME_SCALE +"-"+ntile + ":"+ (anpr - G0_INDEX);
			}

			System.out.println(String.format("%2d%1s %22s %f",
					np,
					(par_mask[np]?"+":" "),
					parname,
					all_pars[np]));
		}
	}

	public void printInputDataFx(boolean show_fx){ // not used in lwir
		if 	(show_fx) {
			Sample s = null;
//			double [] fx = getPolyFx();
//			if (fx == null) fx = getFx();
			double [] fx = getFx();
			if (fx == null) return;
			for (int i = 0; i < fx.length; i++) {
//				double fx_pos = (fx[i] >= 0)? fx[i]: Double.NaN;
				double fx_pos = fx[i];
				if (i < samples.size()) {
					s = samples.get(i);
					System.out.println(String.format("%3d: x=%2d y=%2d v=%9.6f fx=%9.6f w=%9.7f fcam=%1d scam=%1d tile=%d", i, s.ix, s.iy, s.v, fx_pos, s.w, s.fcam, s.scam, s.tile));
				}
				else {
					System.out.println(String.format("%3d: %2s %2s v=%9.6f fx=%9.6f w=%9.7f", i, "-", "-", this.values[i], fx_pos, this.weights[i]));
				}
			}
		} else {
			int ns =0;
			for (Sample s:samples){
				System.out.println(String.format("%3d: x=%2d y=%2d v=%9.6f w=%9.7f fcam=%1d scam=%1d tile=%d", ns++, s.ix, s.iy, s.v, s.w, s.fcam, s.scam, s.tile));
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

	public double [] getDisparityStrength(int nTile) { // USED in lwir
////		if (pair_weights == null) return null;
		double disparity = -all_pars[DISP_INDEX];
		double sum_amp = 0.0;
		for (int i = 0; i < NUM_PAIRS; i++) {
////			sum_amp += pair_weights[i] * all_pars[G0_INDEX + i]; // group_weights is normalized
			sum_amp += all_pars[G0_INDEX + i + TILE_PARAMS * nTile]; // group_weights is normalized
		}
		// protect from weird fitting results
		double max_amp = 0.0;
		for (Sample s: samples) if ((s.v > max_amp) && (s.tile == nTile)) max_amp = s.v;
		if (sum_amp > 1.25 * max_amp) sum_amp = max_amp;
		double [] ds = {disparity, sum_amp};
		return ds;
	}

	public double [] getDisparityStrengthABC(int nTile) {// width = 1/sqrt(all_pars[A_INDEX])
		double [] ds = getDisparityStrength(nTile);
		if (ds == null) return null;
		double [] dsw = {ds[0], ds[1], all_pars[A_INDEX + TILE_PARAMS * nTile], all_pars[B_INDEX+ TILE_PARAMS * nTile],all_pars[CMA_INDEX+ TILE_PARAMS * nTile]}; // asymmetry
		return dsw;
	}

	public double [] getLazyEye() {
		double [] rslt = new double [2 * NUM_CAMS];
		for (int i = 0; i < rslt.length; i++) {
			rslt[i] = all_pars[DDISP_INDEX + i];
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
		all_pars[DDISP_INDEX + last_cam] = 0.0;
		for (int i = 0; i < NUM_CAMS; i++) {
			if (used_cameras[i] & (i != last_cam)) {
				all_pars[DDISP_INDEX + last_cam] -= all_pars[DDISP_INDEX + i];
			}
		}

	}



	public double [] fromVector(double [] vector) { // mix fixed and variable parameters // USED in lwir
		if ( all_pars == null) return null;
		double [] ap = all_pars.clone();
		int np = 0;
		for (int i = 0; i < par_mask.length; i++) if (par_mask[i]) ap[i] = vector[np++];

		ap[DDISP_INDEX + last_cam] = 0.0;
		for (int i = 0; i < NUM_CAMS; i++) {
			if (used_cameras[i] & (i != last_cam)) {
				ap[DDISP_INDEX + last_cam] -= ap[DDISP_INDEX + i];
			}
		}
		return ap;
	}

////////////////////////////////////////////////////////////////////////

    public void debugJt( // not used in lwir
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
    	System.out.println("Test of jt-jt_delta difference,  delta = "+delta+ ":");
    	System.out.print(String.format("Til P %3s: %10s ", "#", "fx"));
    	for (int anp = 0; anp< all_pars.length; anp++) if(par_mask[anp]){
			String parname;
			if      (anp >= NDISP_INDEX) parname = PAR_NAME_CORRNDISP + (anp - NDISP_INDEX);
			else if (anp >= DDISP_INDEX) parname = PAR_NAME_CORRDISP +  (anp - DDISP_INDEX);
			else {
				int ntile = anp / TILE_PARAMS;
				int anpr =  anp % TILE_PARAMS;
				if (anpr < G0_INDEX) parname = PAR_NAMES[anpr]+"-"+ntile;
				else                 parname = PAR_NAME_SCALE +"-"+ntile + ":"+ (anpr - G0_INDEX);
			}
        	System.out.print(String.format("| %16s ", parname));
    	}
    	System.out.println();
    	int npair0 = -1;
    	for (int i = 0; i < num_points; i++) {
    		if (i < samples.size()) {
        		int npair = USED_PAIRS_MAP[samples.get(i).tile][samples.get(i).fcam][samples.get(i).scam];
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
    	System.out.print(String.format("      %15s ", "Maximal diff:"));
    	for (int np = 0; np < num_pars; np++) {
        	System.out.print(String.format("|%8s %8.5f ", "1/1000×",  1000*max_diff[np]));
    	}
    	System.out.println();
    }


	public double [] getFxJt( // not used in lwir
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



	public double [] getFx() { // not used in lwir
		return getFxJt(this.vector, null);
	}


	public double [][] getWJtJlambda( // USED in lwir
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

	public double [] getJtWdiff( // not used in lwir
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

	public boolean runLma( // USED in lwir
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
			System.out.println("LMA: full RMS="+last_rms[0]+" ("+initial_rms[0]+"), pure RMS="+last_rms[1]+" ("+initial_rms[1]+") + lambda="+lambda);
		}

		return rslt[0];
	}


	// returns {success, done}
	public boolean [] lmaStep( // USED in lwir
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

}
