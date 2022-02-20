package com.elphel.imagej.tileprocessor;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Iterator;
import java.util.concurrent.atomic.AtomicInteger;

import com.elphel.imagej.common.DoubleGaussianBlur;
import com.elphel.imagej.common.PolynomialApproximation;
import com.elphel.imagej.common.ShowDoubleFloatArrays;
import com.elphel.imagej.tileprocessor.Corr2dLMA.Sample;

import Jama.Matrix;

/**
 **
 ** Correlation2d - Handle 2-d (phase) correlations, combining multiple-pair data
 **
 ** Copyright (C) 2018 Elphel, Inc.
 **
 ** -----------------------------------------------------------------------------**
 **
 **  Correlation2d.java is free software: you can redistribute it and/or modify
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

public class Correlation2d {
	public final static int PAIR_HORIZONTAL =      0;
	public final static int PAIR_VERTICAL =        1;
	public final static int PAIR_DIAGONAL_MAIN =   2;
	public final static int PAIR_DIAGONAL_OTHER =  3;
	
	public final static int CAMEL_BOTH =      0;
	public final static int CAMEL_STRONGEST = 1;
	public final static int CAMEL_NEAREST =   2;
	public final static int CAMEL_FG =        3;
	public final static int CAMEL_BG =        4;
	
	
	private final DttRad2 dtt;
	private final int transform_size;
	private final int transform_len;
	private final int corr_size;
	private final int [] transpose_all_ortho;
	private final int [] transpose_all_diagonal;
	private final double [] ortho_notch_filter;
	private final double [] corr_wndx;
	private final double [] corr_wndy;
	private final double [] corr_wndy_notch;
	private final boolean monochrome;

	public static int CORR_SEL_BIT_ALL =   0; 
	public static int CORR_SEL_BIT_DIA =   1; 
	public static int CORR_SEL_BIT_SQ =    2;
	public static int CORR_SEL_BIT_NEIB =  3;
	public static int CORR_SEL_BIT_HOR =   4;
	public static int CORR_SEL_BIT_VERT =  5;
	
	public static int CORR_SEL_BIT_LIMIT =  8;
	public static int CORR_SEL_BITS_LIMIT = 3;
// Debug feature to limit selections to:
// 1: 2 sensors - horizontal diameter, sensors  4,12
// 2: 4 sensors - original quad selection, sensors 2,6,10,14
// 3: 8 sensors - 0,2,4,6,8,10,12,14	

	public static int corrSelEncode(
			boolean sel_all,
			boolean sel_dia,
			boolean sel_sq,
			boolean sel_neib,
			boolean sel_hor,
			boolean sel_vert,
			int     limit_sensors) {
		return  (sel_all ?  (1 << CORR_SEL_BIT_ALL):  0) |
				(sel_dia ?  (1 << CORR_SEL_BIT_DIA):  0) |
				(sel_sq ?   (1 << CORR_SEL_BIT_SQ):   0) |
				(sel_neib ? (1 << CORR_SEL_BIT_NEIB): 0) |
				(sel_hor ?  (1 << CORR_SEL_BIT_HOR):  0) |
				(sel_vert ? (1 << CORR_SEL_BIT_VERT): 0) |
				((limit_sensors & ((1 << CORR_SEL_BITS_LIMIT) - 1)) << CORR_SEL_BIT_LIMIT);
	}
	public static boolean isCorrAll (int sel) { return ((sel >> CORR_SEL_BIT_ALL)  & 1) != 0;}
	public static boolean isCorrDia (int sel) { return ((sel >> CORR_SEL_BIT_DIA)  & 1) != 0;}
	public static boolean isCorrSq  (int sel) { return ((sel >> CORR_SEL_BIT_SQ)   & 1) != 0;}
	public static boolean isCorrNeib(int sel) { return ((sel >> CORR_SEL_BIT_NEIB) & 1) != 0;}
	public static boolean isCorrHor (int sel) { return ((sel >> CORR_SEL_BIT_HOR)  & 1) != 0;}
	public static boolean isCorrVert(int sel) { return ((sel >> CORR_SEL_BIT_VERT) & 1) != 0;}
	public static int getSensorsLimit(int sel){ return ((sel >> CORR_SEL_BIT_LIMIT) & ((1 << CORR_SEL_BITS_LIMIT) - 1));}

	public static int corrSelEncode(ImageDttParameters  img_dtt, int num_sensors) {
		return corrSelEncode(
				img_dtt.getMcorrAll  (num_sensors),  // boolean sel_all,
				img_dtt.getMcorrDia  (num_sensors),  // boolean sel_dia,
				img_dtt.getMcorrSq   (num_sensors),  // boolean sel_sq,
				img_dtt.getMcorrNeib (num_sensors),  // boolean sel_neib,
				img_dtt.getMcorrHor  (num_sensors),  // boolean sel_hor,
				img_dtt.getMcorrVert (num_sensors),  // boolean sel_vert);
				img_dtt.mcorr_limit_sensors); // 0 - no limit, 1 - (4,12) 2 - (2, 6, 10, 14), 3 - (0,2,4,6,8,10,12,14)
	}

	public static int corrSelEncodeAll(int limit_sensors) { // 0 - no limit, 1 - (4,12) 2 - (2, 6, 10, 14), 3 - (0,2,4,6,8,10,12,14)
		return corrSelEncode(
				true,  // boolean sel_all,
				false,  // boolean sel_dia,
				false,  // boolean sel_sq,
				false,  // boolean sel_neib,
				false,  // boolean sel_hor,
				false,  // boolean sel_vert);
				limit_sensors); // 0 - no limit, 1 - (4,12) 2 - (2, 6, 10, 14), 3 - (0,2,4,6,8,10,12,14)
	}
	
	  
	// configuration for 8-lens and 4-lens cameras. 8-lens has baseline = 1 for 1..4 and 1/2 for 4..7
/*0        1
     4  5
     6  7
  2        3 */
	final static int [][] PAIRS ={ // {first, second, orientation, scale}
			{0, 1, PAIR_HORIZONTAL,     1},
			{2, 3, PAIR_HORIZONTAL,     1},
			{0, 2, PAIR_VERTICAL,       1},
			{1, 3, PAIR_VERTICAL,       1},
			{0, 3, PAIR_DIAGONAL_MAIN,  1},
			{2, 1, PAIR_DIAGONAL_OTHER, 1},

			{4, 5, PAIR_HORIZONTAL,     2},
			{6, 7, PAIR_HORIZONTAL,     2},
			{4, 6, PAIR_VERTICAL,       2},
			{5, 7, PAIR_VERTICAL,       2},
			{4, 7, PAIR_DIAGONAL_MAIN,  2},
			{6, 5, PAIR_DIAGONAL_OTHER, 2},

			{0, 4, PAIR_DIAGONAL_MAIN,  4},
			{7, 3, PAIR_DIAGONAL_MAIN,  4},
			{2, 6, PAIR_DIAGONAL_OTHER, 4},
			{5, 1, PAIR_DIAGONAL_OTHER, 4},
			};
	final static int [][] GROUPS = { // {diagonal, scale}
			{0, 1},
			{1, 1},
			{0, 2},
			{1, 2},
			{0, 4},
			{1, 4}};

	final double[][] port_offsets = {
			{-0.5, -0.5},
			{ 0.5, -0.5},
			{-0.5,  0.5},
			{ 0.5,  0.5}};
	// This table is used for CLT-based transform of the correlation results to the pixel domain
	final static int [][] ZI =
		{{ 0,  1,  2,  3},
		 {-1,  0, -3,  2},
		 {-2, -3,  0,  1},
		 { 3, -2, -1,  0}};
	
    public final int numSensors;
    public final boolean top_is_0;


	public boolean isMonochrome() {return monochrome;} // not used in lwir
// for 8 cameras and 16 pairs. Following data moved from ImageDtt
    // which images to use (0..3 - external, 4..7 - internal)
    public static int getImgMask  (int data){ return (data & 0xff);} // not used in lwir
    // which pairs to combine in the combo see PAIRS data
    public static int getPairMask (int data){ return ((data >> 8) & 0xffff);} // not used in lwir
    public static int setImgMask  (int data, int mask) {return (data & ~0xff) | (mask & 0xff);} // not used in lwir
    public static int setPairMask (int data, int mask) {return (data & ~0xffff00) | ((mask & 0xffff) << 8);} // not used in lwir
    public static boolean getForcedDisparity (int data){return (data & 0x1000000) != 0;} // not used in lwir
    public static int     setForcedDisparity (int data, boolean force) {return (data & ~0x1000000) | (force?0x1000000:0);} // not used in lwir
    public static boolean getOrthoLines (int data){return (data & 0x2000000) != 0;} // not used in lwir
    public static int     setOrthoLines (int data, boolean ortho_lines) {return (data & ~0x2000000) | (ortho_lines?0x2000000:0);} // not used in lwir
    public static boolean isOrthoPair         (int npair) {return (PAIRS[npair][2]== PAIR_HORIZONTAL) || (PAIRS[npair][2]== PAIR_VERTICAL);} // not used in lwir
    public static boolean isDiagonalPair      (int npair) {return (PAIRS[npair][2]== PAIR_DIAGONAL_MAIN) || (PAIRS[npair][2]== PAIR_DIAGONAL_OTHER);} // USED in lwir
    public static boolean isHorizontalPair    (int npair) {return  PAIRS[npair][2]== PAIR_HORIZONTAL;} // USED in lwir
    public static boolean isVerticalPair      (int npair) {return  PAIRS[npair][2]== PAIR_VERTICAL;} // USED in lwir
    public static boolean isDiagonalMainPair  (int npair) {return  PAIRS[npair][2]== PAIR_DIAGONAL_MAIN;} // USED in lwir
    public static boolean isDiagonalOtherPair (int npair) {return  PAIRS[npair][2]== PAIR_DIAGONAL_OTHER;} // USED in lwir
    public static int     getScaleOfPair      (int npair) {return  PAIRS[npair][3];} // not used in lwir

    public static int     getMaskHorizontal(int scale)    {return getMaskType(PAIR_HORIZONTAL, scale);} // USED in lwir
    public static int     getMaskVertical(int scale)      {return getMaskType(PAIR_VERTICAL, scale);} // USED in lwir
    public static int     getMaskDiagonalMain(int scale)  {return getMaskType(PAIR_DIAGONAL_MAIN, scale);} // not used in lwir
    public static int     getMaskDiagonalOther(int scale) {return getMaskType(PAIR_DIAGONAL_OTHER, scale);} // not used in lwir
    public static int     getMaskOrtho(int scale)         {return getMaskHorizontal(scale) | getMaskVertical(scale);} // not used in lwir
    public static int     getMaskDiagonal(int scale)      {return getMaskDiagonalMain(scale) | getMaskDiagonalOther(scale);} // not used in lwir

    private static int     getMaskType(int type, int scale) { // scale <0 = any  // USED in lwir
    	int bm = 0;
    	for (int i = 0; i <PAIRS.length; i++) if ((PAIRS[i][2]==type) && ((scale < 0) || (scale == PAIRS[i][3]))) bm |= 1 << i;
    	return bm;
    }
    
    private final int [][]  pair_start_end; // start_index, end_index
    private final int [][]  pair_orient;    // second pair parallel to first pair - 0, CW90 - 1, 180 - 2, CW270 - 3, other - -1
    private final int []    pair_length;    // number of sensors from first to second
    private boolean []      corr_pairs;     // Which pairs to calculate
    private boolean []      corr_pairs_filter; // Remove pairs that do not have measurement
    public static enum      MCORR_COMB {ALL, DIA, SQ, NEIB, HOR, VERT}; // MCORR_COMB.SQ.ordinal() == 2
    private static String [] CORR_TITLES_EXTRA={"*all","*dia","*sq","*neib","*hor","*vert"};
    
    public final String []  cor_titles;       // per-pair correlation slices titles
    public final String []  cor_titles_combo; // combo correlation slices titles
    
    
	private int             mcorr_comb_width;  // combined correlation tile width
	private int             mcorr_comb_height; // combined correlation tile full height
	private int             mcorr_comb_offset; // combined correlation tile height offset: 0 - centered (-height/2 to height/2), height/2 - only positive (0 to height)
	private double          mcorr_comb_disp;   // Combined tile per-pixel disparity for baseline == side of a square
//	private double [][][][] resample;    // [num_pair][out_index]<null or variable number of {source_index, weight} pairs>
	
	private int    [][][]   resample_indices;  // should have the same length as number of pairs
	private double [][][]   resample_weights;  // should have the same length as number of pairs

	
	private static int SUB_SAMPLE = 16;  // subsample source pixel in each direction when generating
	public static int THREADS_MAX = 100; 
	public int getPairLength(int pair) {
		return pair_length[pair];
	}
	public int [] getPairLengths() {
		return pair_length;
	}
	
	// All used pairs (but the diameters) are clockwise (end is clockwise of start)
	// Orientation calculations are valid for clockwise only
    private void setupPairs() {
    	int indx = 0;
    	for (int i = 1; i <= numSensors/2; i++) { // CW length
    		for (int j = 0; j  < ((i < (numSensors/2))? numSensors : (numSensors/2) ); j++) {
    			pair_length[indx] = i;
    			pair_start_end[indx][0] = j;
    			pair_start_end[indx][1] = (i + j) % numSensors;
    			indx++;
 //   			if (indx >= pair_length.length) {
 //   				break; // BUG
 //   			}
    		}
    	}
    	for (int i = 0; i < pair_start_end.length; i++ ) {
    		int n1 = (2 * pair_start_end[i][0] + pair_length[i]) % (2 * numSensors); 
        	for (int j = 0; j < pair_start_end.length; j++ ) {
        		int n2 = (2 * pair_start_end[j][0] + pair_length[j]) % (2 * numSensors);
//        		if ((i >= 116) && (j >= 112)) {
//        			System.out.println ("i="+i+", j="+j);
//        			System.out.println ("(n2 - n1) % (numSensors / 2)="+((n2 - n1) % (numSensors / 2)));
//        			System.out.println ("((n2 - n1) / (numSensors / 2)) % 4 = "+(((n2 - n1) / (numSensors / 2)) % 4));
//        		}
        		if ((n2 - n1) % (numSensors / 2) != 0) {
        			pair_orient[i][j] = -1; // non-ortho
        		} else {
        			pair_orient[i][j] = Math.floorMod(((n2 - n1) / (numSensors / 2)) , 4);
        		}
        	}
    	}
//    	cor_titles = new String[num_pairs + CORR_TITLES_EXTRA.length];
    	for (int i = 0; i < pair_start_end.length; i++ ) {
    		cor_titles[i] = getPair(i)[0]+"-"+getPair(i)[1]; //getPair(i) to correct sensor numbers for quad (Z-order)
    	}
    	for (int i = 0; i < CORR_TITLES_EXTRA.length; i++ ) {
    		cor_titles_combo[i] = CORR_TITLES_EXTRA[i];
    	}
    	return;
    }
    
    public boolean [] selectLength(boolean [] cur_sel, int l) {
    	boolean [] sel = (cur_sel == null) ? (new boolean[pair_length.length]) : cur_sel;
    	for (int i = 0; i < pair_length.length; i++) {
    		sel[i] |= pair_length[i] == l;
    	}
    	return sel; 
    }

    public boolean [] selectDiameters(boolean [] cur_sel) {
    	return selectLength(cur_sel, numSensors / 2); 
    }

    public boolean [] selectSquares(boolean [] cur_sel) {
    	return selectLength(cur_sel, numSensors / 4); 
    }

    public boolean [] selectNeibs(boolean [] cur_sel) {
    	return selectLength(cur_sel, 1); 
    }
    
    public boolean [] selectParallel(boolean [] cur_sel, int start, int end) {
    	boolean [] sel = (cur_sel == null) ? (new boolean[pair_length.length]) : cur_sel;
    	int l = (end + numSensors - start) % numSensors; // CCW pairfs will have l > numSensors/2
    	
		int n1 = (2 * start + l ) % (2 * numSensors); 
    	for (int j = 0; j < pair_start_end.length; j++ ) { // lengths are always positive
    		int n2 = (2 * pair_start_end[j][0] + pair_length[j]) % (2 * numSensors);
    		if ((n2 - n1) % numSensors == 0) { // only parallel or anti-parallel
    			sel[j] = true;
    		}
    	}
    	return sel; 
    }

    public boolean [] selectHorizontal(boolean [] cur_sel) {
    	if (top_is_0) {
    		return selectParallel(cur_sel, numSensors - 1, 1);
    	} else {
    		return selectParallel(cur_sel, numSensors - 1, 0);
    	}
    }
    
    public boolean [] selectVertical(boolean [] cur_sel) {
    	if (top_is_0) {
    		return selectParallel(cur_sel, 0 , numSensors/2);
    	} else {
    		return selectParallel(cur_sel, 0, numSensors/2 - 1);
    	}
    }
    
    // combine with available correlation results
    
    
    public boolean [] selectAll() {
    	boolean [] sel = new boolean[pair_length.length];
    	Arrays.fill(sel, true);
    	return sel; 
    }

    public int [] getPair(int indx) {
    	int [] quad_indx = {1,3,2,0}; // CW to EO numbers
    	if (numSensors > 4) {
    		return pair_start_end[indx];
    	} else {
    		return new int [] {quad_indx[pair_start_end[indx][0]], quad_indx[pair_start_end[indx][1]]};
    	}
    }
    public int getNumPairs() {
    	return pair_start_end.length;
    }
    
    public static int getNumPairs(int numSensors) {
     return  numSensors * (numSensors-1) /2;
    }
    
    public static boolean [] boolCorrPairs(int mcorr_sel, int numSensors) {
    	Correlation2d c2d = new Correlation2d(
    			numSensors, // other parameters do not matter
    			8, 
    			true,
    			false);
		return c2d.setCorrPairs(mcorr_sel);
    }
    
    public static int [] intCorrPairs(int mcorr_sel, int numSensors, int num_out) {
    	boolean [] boolCorrPairs = boolCorrPairs(mcorr_sel, numSensors);
    	int [] int_pairs = new int [num_out]; // normally 4
    	for (int i = 0; i < boolCorrPairs.length; i++) if (boolCorrPairs[i]){
    		int_pairs[i >> 5] |= 1 << (i & 31);
    	}
    	return int_pairs;
    	
    }
    
    public boolean[] setCorrPairs(int mcorr_sel) {
		boolean [] corr_calculate = null;
		if (isCorrAll  (mcorr_sel)) corr_calculate = selectAll();
		if (isCorrDia  (mcorr_sel)) corr_calculate = selectDiameters  (corr_calculate);
		if (isCorrSq   (mcorr_sel)) corr_calculate = selectSquares    (corr_calculate);
		if (isCorrNeib (mcorr_sel)) corr_calculate = selectNeibs      (corr_calculate);
		if (isCorrHor  (mcorr_sel)) corr_calculate = selectHorizontal (corr_calculate);
		if (isCorrVert (mcorr_sel)) corr_calculate = selectVertical   (corr_calculate);
		
		int sensors_limit = getSensorsLimit(mcorr_sel);
		if ((sensors_limit != 0) && (numSensors == 16)) { // so far only for 16 sensors
			boolean [] sel_sensors = new boolean [numSensors];
			int [][] indices = {
					{4,12},                // binocular,
					{2,6,10,14},           // quad
					{0,2,4,6,8,10,12,14}}; // octal
			for (int i = 0; i < indices[sensors_limit-1].length; i++) {
				sel_sensors[indices[sensors_limit-1][i]] = true;
			}
//			boolean [] pair_sel = getCorrPairs();
			for (int i = 0; i < corr_calculate.length; i++) if (corr_calculate[i]){
				int [] se = getPair(i);
				if (!sel_sensors[se[0]] || !sel_sensors[se[1]]) {
					corr_calculate[i] = false;
				}
			}
//			setCorrPairs(pair_sel);
			//*********************************** limit pairs
		}
		setCorrPairs(corr_calculate);
		return corr_calculate;
    }
    
    
    public void setCorrPairs(boolean [] sel) { // these pairs will be correlated
    	corr_pairs = sel.clone();
    }
    
    public int getNumSensors() {
    	return this.numSensors;
    }
    
    public void setCorrPairsFilter(double [][][][][][] clt_data, int tileY, int tileX ) { // these pairs will be correlated
    	corr_pairs_filter = new boolean [pair_start_end.length];
    	boolean [] en = new boolean [numSensors];
    	if (clt_data != null) {
    		for (int i = 0; i < numSensors; i++) if (clt_data[i] != null){
    			for (int ncol = 0; ncol < clt_data[i].length; ncol++) {
    				if ((clt_data[i][ncol] != null) && (clt_data[i][ncol][tileY] != null) && (clt_data[i][ncol][tileY][tileX] != null)) {
    					en[i] = true;
    					break;
    				}
    			}
    		}
    		for (int i = 0; i < pair_start_end.length; i++) {
    			corr_pairs_filter[i] = en[pair_start_end[i][0]] && en[pair_start_end[i][1]];
    		}
    	}
    }
    
    /**
     * Filter sel array by removing unavailable sensors and pairs that are not calculated 
     * @param sel selection to be filtered
     * @return number of remaining selected pairs 
     */
    public int numRemainingPairs(boolean [] sel ) {
    	int num_sel = 0;
    	for (int i = 0; i < sel.length; i++) {
    		sel[i] &= corr_pairs[i] &&  corr_pairs_filter[i];
    		if (sel[i]) num_sel++;
    	}
    	return num_sel;
    }
    

    public boolean [] getCorrPairs() { // these pairs will be correlated
    	return corr_pairs;
    }
    
    public String[] getCorrTitles() {
    	return cor_titles;
    }
    public String[] getComboTitles() {
    	return cor_titles_combo;
    }
    
    public static boolean [] longToArray(long sel_bits, int num_pairs) {
    	boolean [] sel = new boolean [num_pairs];
    	for (int i = 0; i < sel.length;i++) {
    		sel[i] = (sel_bits & 1) > 0;
    		sel_bits >>= 1;
    	}
    	return sel;
    }

    public boolean [] longToArray(long sel_bits) {
    	return longToArray (sel_bits, pair_start_end.length);
    }
    
    
    /**
     * Add 2D correlation for a pair to the combined correlation tile, applying rotation/scaling
     * @param accum_tile tile for accumulation in line-scan order, same dimension as during generateResample()
     * @param corr_tile correlation tile (currently 15x15)
     * @param num_pair number of correlation pair for which resampling data exists
     * @param weight multiply added tile data by this coefficient before accumulation
     */
    public void accummulatePair(
    		double [] accum_tile,
    		double [] corr_tile,
    		int       num_pair,
    		double    weight) {
    	if ((resample_indices == null) || (resample_indices[num_pair] == null)) {
    		throw new IllegalArgumentException ("No resample data for num_pair = "+num_pair);
    	}
    	for (int i = 0; i < accum_tile.length; i++) if (resample_indices[num_pair][i] != null) {
    		for (int j = 0; j < resample_indices[num_pair][i].length; j++) {
    			accum_tile[i] += weight * corr_tile[resample_indices[num_pair][i][j]] * resample_weights[num_pair][i][j];
    		}
    	}
    }
    
    /**
     * Add multiple 2D correlation tiles to a combined correlation tile, applying rotation/scaling
     * @param corr_tiles (sparse) array of correlation tile, some tiles may be null
     * @param selection boolean selection array that specifies which (of existing) correlation tiles to add
     * @param weight  multiply added tile data by this coefficient before accumulation, common for all added tiles
     */
    public double accummulatePairs(
    		double []   accum_tile,
    		double [][] corr_tiles,
    		boolean []  selection,
    		double      weight) { // same weights
    	double sumw = 0.0;
    	for (int num_pair = 0; num_pair < selection.length; num_pair++) {
    		if (selection[num_pair] && (corr_tiles[num_pair] != null)) {
    			accummulatePair(
    					accum_tile,
    		    		corr_tiles[num_pair],
    		    		num_pair,
    		    		weight);
    			sumw+=weight;
    		}
    	}
    	return sumw;
    }
    
    public double accummulatePairs(
    		double []   accum_tile,
    		double [][] corr_tiles,
    		double  []   weights) {
    	double sumw = 0.0;
    	for (int num_pair = 0; num_pair < weights.length; num_pair++) {
    		if ((weights[num_pair]>0.0) && (corr_tiles[num_pair] != null)) {
    			accummulatePair(
    					accum_tile,
    		    		corr_tiles[num_pair],
    		    		num_pair,
    		    		weights[num_pair]);
    			sumw+=weights[num_pair];
    		}
    	}
    	return sumw;
    }
    
    
    
    public void normalizeAccumulatedPairs(
    		double []   accum_tile,
    		double      sumw) {
    	if ((accum_tile != null) && (sumw > 0)) {
    		double k = 1.0/sumw;
    		for (int i= 0; i < accum_tile.length; i++) {
    			accum_tile[i] *= k;
    		}
    	}
    }
    
    /**
     * Calculate correlation pair width in the disparity direction to boost weights of the most
     * informative (for disparity) pairs 
     * @param corr_tile for selected pair 
     * @param num_pair selected pair number
     * @return pair width (in pixels) in the disparity direction after applying rotation and scaling
     */
    public double getPairWidth(
    		double [] corr_tile,
    		int       num_pair){ // Center line only
     	if ((resample_indices == null) || (resample_indices[num_pair] == null)) {
    		throw new IllegalArgumentException ("getPairWidth(): No resample data for num_pair = "+num_pair);
    	}
    	double s0 = 0.0, s1 = 0.0, s2 = 0.0;
		int indx = mcorr_comb_width * (mcorr_comb_offset + ((mcorr_comb_height - 1) / 2));
    	for (int ix = 0; ix < mcorr_comb_width; ix ++ ) if (resample_indices[num_pair][indx] != null){
    		for (int j = 0; j < resample_indices[num_pair][indx].length; j++) {
    			double d = corr_tile[resample_indices[num_pair][indx][j]] * resample_weights[num_pair][indx][j];
    			if (d < 0) {
    				d = 0.0; // was get5ting negative (s0*s2- s1*s1) !
    			}
    			s0 += d;
    			s1 += d * ix; 
    			s2 += d * ix * ix; 
    		}
    		indx++;
    	}
    	return Math.sqrt(s0*s2- s1*s1)/s0;
    }
    
    
    public double [] accumulateInit() {return new double [mcorr_comb_width * mcorr_comb_height]; }
    public int getCombWidth() {return mcorr_comb_width;}
    public int getCombHeight() {return mcorr_comb_height;}
    public int getCombOffset() {return mcorr_comb_offset;}
    public double getCombDisp() {return mcorr_comb_disp;}
    
    @Deprecated
    public void generateResampleOld( // should be called before
			final int           mcorr_comb_width,  // combined correlation tile width
			final int           mcorr_comb_height, // combined correlation tile full height
			final int           mcorr_comb_offset, // combined correlation tile height offset: 0 - centered (-height/2 to height/2), height/2 - only positive (0 to height)
			final double        mcorr_comb_disp){  // Combined tile per-pixel disparity for baseline == side of a square
    	
		this.mcorr_comb_width =  mcorr_comb_width;  // combined correlation tile width
		this.mcorr_comb_height = mcorr_comb_height; // combined correlation tile full height
		this.mcorr_comb_offset = mcorr_comb_offset; // combined correlation tile height offset: 0 - centered (-height/2 to height/2), height/2 - only positive (0 to height)
		this.mcorr_comb_disp =   mcorr_comb_disp;  // Combined tile per-pixel disparity for baseline == side of a square
		
		resample_indices = new int     [corr_pairs.length][][];
		resample_weights = new double  [corr_pairs.length][][];
		
		
		final int weights_size = 2*SUB_SAMPLE-1;
		final double [] weights1 = new double [weights_size];
		final double [] weights = new double [weights_size * weights_size];
		for (int i = 0; i < weights1.length; i++) {
			double w = Math.sin(Math.PI*(i+1)/(2*SUB_SAMPLE));
			weights1[i] = w * w;
		}
		double s = 0.0;
		int indx = 0;
		for (int i = 0; i < weights1.length; i++) {
			for (int j = 0; j < weights1.length; j++) {
				double w = weights1[i] * weights1[j];
				s+=w;
				weights[indx++] = w;
			}
		}
		// normalize sum == 1.0;
		double k = 1.0/s;
		for (int i = 0; i < weights.length; i++) {
			weights[i] *= k;
		}
		
		
//    	final double [][][][] resample = new double [pair_start_end.length][mcorr_comb_width * mcorr_comb_height][][];
    	// use multithreading?
		final Thread[] threads = ImageDtt.newThreadArray(THREADS_MAX);
		final AtomicInteger ai = new AtomicInteger(0);
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				public void run() {
					double [] contrib = new double [corr_size*corr_size];
					HashSet<Integer> contrib_set = new HashSet<Integer>();
					Iterator<Integer> contrib_itr;
					for (int num_pair = ai.getAndIncrement(); num_pair < corr_pairs.length; num_pair = ai.getAndIncrement()) {
//						if (num_pair == 62) {
//							System.out.println("num_pair="+num_pair);
//							System.out.println("num_pair="+num_pair);
//						}
						if (corr_pairs[num_pair]) {
							resample_indices[num_pair] = new int    [mcorr_comb_width * mcorr_comb_height][];
							resample_weights[num_pair] = new double [mcorr_comb_width * mcorr_comb_height][];
							// scale and angle
							int istart = pair_start_end[num_pair][0];
							int iend =   pair_start_end[num_pair][1];
							double cwrot = (((istart+iend) % numSensors) + (top_is_0 ? 0.0 : 0.5)) * Math.PI/numSensors; //  +(top_is_0 ? 0.0 : 0.5)
							if ((iend > istart) ^ ((iend + istart) <16)) {
								cwrot += Math.PI;
							}
							double pl = 2 * Math.sin(Math.PI* pair_length[num_pair] / numSensors); // length for R=1
//							double scale = pl/Math.sqrt(2.0)* mcorr_comb_disp; // Math.sqrt(2.0) - relative to side of a square - may be change later?
							double scale = pl/2.0* mcorr_comb_disp; // Math.sqrt(2.0) - relative to diameter
							Matrix toPair = new Matrix(new double[][] {
								{scale * Math.cos(cwrot), -scale*Math.sin(cwrot)},
								{scale * Math.sin(cwrot),  scale*Math.cos(cwrot)}});
							
							// scale - to get pair (source) radius from combo (destination) radius
							double ksub = 1.0/SUB_SAMPLE;
							Matrix mxy = new Matrix(2,1);
							for (int i = 0; i < mcorr_comb_height; i++) {
								int iy = i + mcorr_comb_offset - mcorr_comb_height/2;
								for (int j = 0; j < mcorr_comb_width; j++) {
									int ix = j - mcorr_comb_width/2;
									Arrays.fill(contrib, 0.0);
									contrib_set.clear();
//									if ((num_pair == 62) && (i==7) && (j==7)) {
//										System.out.println("num_pair="+num_pair+", i="+i+", j="+j);
//										System.out.println("num_pair="+num_pair+", i="+i+", j="+j);
//									}
									
									for (int idy = 0; idy < weights_size; idy++) {
										mxy.set(1, 0,  iy+ (idy - SUB_SAMPLE + 1) * ksub); //idy == (SUB_SAMPLE -1) - no fractional pixel
										for (int idx = 0; idx < weights_size; idx++) {
											mxy.set(0, 0,  ix+ (idx - SUB_SAMPLE + 1) * ksub); // idy == (SUB_SAMPLE -1) - no fractional pixel
											double [] pxy = toPair.times(mxy).getColumnPackedCopy();
//											int ipx = (int) Math.round(pxy[0]+transform_size -1);
//											int ipy = (int) Math.round(pxy[1]+transform_size -1);
											// round symmetrically (away from zero)
											double dpx = pxy[0]+transform_size -1;
											double dpy = pxy[1]+transform_size -1;
											int ipx = (int) Math.round(Math.abs(dpx));
											int ipy = (int) Math.round(Math.abs(dpy));
											if (dpx < 0) ipx = -ipx;
											if (dpy < 0) ipy = -ipy;
											
											if ((ipx >= 0) && (ipy >= 0) && (ipx < corr_size) && (ipy < corr_size)) {
												int indx_src = ipy * corr_size + ipx;
												contrib_set.add(indx_src);
												contrib[indx_src] += weights[idy * weights_size + idx];
											}
										}
									}
									if (!contrib_set.isEmpty()) {
										int indx = i * mcorr_comb_width + j;
										resample_indices[num_pair][indx] = new int    [contrib_set.size()];
										resample_weights[num_pair][indx] = new double [contrib_set.size()];
										contrib_itr = contrib_set.iterator();
										int contrib_num = 0;
										while(contrib_itr.hasNext()) {
											int indx_src = contrib_itr.next();
											resample_indices[num_pair][indx][contrib_num] = indx_src;
											resample_weights[num_pair][indx][contrib_num] = contrib[indx_src];
											contrib_num++;
										}
									}
								}
							}
						} else {
							resample_indices[num_pair] = null;
							resample_weights[num_pair] = null;
						}
					}
				}
			};
		}		      
		ImageDtt.startAndJoin(threads);
    }

    public void generateResample( // should be called before
 			final int           mcorr_comb_width,  // combined correlation tile width
 			final int           mcorr_comb_height, // combined correlation tile full height
 			final int           mcorr_comb_offset, // combined correlation tile height offset: 0 - centered (-height/2 to height/2), height/2 - only positive (0 to height)
 			final double        mcorr_comb_disp){  // Combined tile per-pixel disparity for baseline == side of a square -> diameter !
     	final double ignore_contrib = 0.001; // ignore contributors with weight below 
 		this.mcorr_comb_width =  mcorr_comb_width;  // combined correlation tile width
 		this.mcorr_comb_height = mcorr_comb_height; // combined correlation tile full height
 		this.mcorr_comb_offset = mcorr_comb_offset; // combined correlation tile height offset: 0 - centered (-height/2 to height/2), height/2 - only positive (0 to height)
 		this.mcorr_comb_disp =   mcorr_comb_disp;  // Combined tile per-pixel disparity for baseline == side of a square
 		
 		resample_indices = new int     [corr_pairs.length][][];
 		resample_weights = new double  [corr_pairs.length][][];
 		final double [][] four_corners = {{-1,-1},{1,-1},{-1, 1},{1, 1}};
// 		final int corr_size = 2 * transform_size - 1;
 		final int corr_center_offs =  2 * transform_size * (transform_size -1);
//     	final double [][][][] resample = new double [pair_start_end.length][mcorr_comb_width * mcorr_comb_height][][];
     	// use multithreading?
 		final Thread[] threads = ImageDtt.newThreadArray(THREADS_MAX);
 		final AtomicInteger ai = new AtomicInteger(0);
 		for (int ithread = 0; ithread < threads.length; ithread++) {
 			threads[ithread] = new Thread() {
 				public void run() {
 					ArrayList<Integer> contrib_list = new ArrayList<Integer>();
 					ArrayList<Double>  contrib_weights_list = new ArrayList<Double>();
 					for (int num_pair = ai.getAndIncrement(); num_pair < corr_pairs.length; num_pair = ai.getAndIncrement()) {
// 						if (num_pair == 62) {
// 							System.out.println("num_pair="+num_pair);
// 							System.out.println("num_pair="+num_pair);
// 						}
 						if (corr_pairs[num_pair]) {
 							resample_indices[num_pair] = new int    [mcorr_comb_width * mcorr_comb_height][];
 							resample_weights[num_pair] = new double [mcorr_comb_width * mcorr_comb_height][];
 							// scale and angle
 							int istart = pair_start_end[num_pair][0];
 							int iend =   pair_start_end[num_pair][1];
 							double cwrot = (((istart+iend) % numSensors) + (top_is_0 ? 0.0 : 0.5)) * Math.PI/numSensors; //  +(top_is_0 ? 0.0 : 0.5)
 							if ((iend > istart) ^ ((iend + istart) <16)) {
 								cwrot += Math.PI;
 							}
 							double pl = 2 * Math.sin(Math.PI* pair_length[num_pair] / numSensors); // length for R=1
 							// scale - to get pair (source) radius from combo (destination) radius
//							double scale = pl/Math.sqrt(2.0)* mcorr_comb_disp; // Math.sqrt(2.0) - relative to side of a square - may be change later?
							double scale = pl/2.0* mcorr_comb_disp; // Math.sqrt(2.0) - relative to diameter
 							Matrix toPair = new Matrix(new double[][] {
 								{scale * Math.cos(cwrot), -scale*Math.sin(cwrot)},
 								{scale * Math.sin(cwrot),  scale*Math.cos(cwrot)}});
 							Matrix fromPair = toPair.inverse();
 							Matrix mxy = new Matrix(2,1);
 							Matrix mpxy = new Matrix(2,1);
 							// 2 methods depending on scale
 							// try with >= 0.9999 to switch branches and compare results
/// 							if (scale >= 0.9999) { // 1.0) { // pair grid is finer, than accumulated grid
  							if (scale >= 1.0) { // pair grid is finer, than accumulated grid
 								for (int i = 0; i < mcorr_comb_height; i++) {
 									int iy = i + mcorr_comb_offset - mcorr_comb_height/2;
 									for (int j = 0; j < mcorr_comb_width; j++) {
 										int ix = j - mcorr_comb_width/2;
// 										if ((num_pair == 62) && (i==7) && (j==7)) {
// 											System.out.println("num_pair="+num_pair+", i="+i+", j="+j);
// 											System.out.println("num_pair="+num_pair+", i="+i+", j="+j);
// 										}
 		 								// convert +/- 1 pixel to pair
 										double minXPair = Double.NaN,minYPair = Double.NaN, maxXPair = Double.NaN, maxYPair = Double.NaN;
 										for (int d = 0; d < four_corners.length; d++) {
 											mxy.set(0, 0,  ix + four_corners[d][0]);
 											mxy.set(1, 0,  iy + four_corners[d][1]);
 											double [] xy_pair = toPair.times(mxy).getColumnPackedCopy();
 											if (d == 0) {
 												minXPair = xy_pair[0];
 												maxXPair = minXPair;
 												minYPair = xy_pair[1];
 												maxYPair = minYPair;
 											} else {
 												minXPair = Math.min(minXPair,xy_pair[0]);
 												minYPair = Math.min(minYPair,xy_pair[1]);
 												maxXPair = Math.max(maxXPair,xy_pair[0]);
 												maxYPair = Math.max(maxYPair,xy_pair[1]);
 											}
 										}
 										int iMinXPair = (int) Math.floor(minXPair);
 										int iMinYPair = (int) Math.floor(minYPair);
 										int iMaxXPair = (int) Math.ceil (maxXPair);
 										int iMaxYPair = (int) Math.ceil (maxYPair);
 										// limit by available data
 										if (iMinXPair < (1 - transform_size)) iMinXPair = 1- transform_size; 
 										if (iMaxXPair > (transform_size - 1)) iMaxXPair = transform_size -1; 
 										if (iMinYPair < (1 - transform_size)) iMinYPair = 1- transform_size; 
 										if (iMaxYPair > (transform_size - 1)) iMaxYPair = transform_size -1; 
 										// corr_center_offs , corr_len
 										if ((iMaxXPair >= iMinXPair) && (iMaxYPair >= iMinYPair)) {
 	 			 							contrib_list.clear();
 	 			 							contrib_weights_list.clear();
 	 			 							double sumw = 0.0;
 											for (int ipy= iMinYPair; ipy <= iMaxYPair; ipy ++) {
 												mpxy.set(1, 0,  ipy);
 												for (int ipx= iMinXPair; ipx <= iMaxXPair; ipx ++) {
 													mpxy.set(0, 0,  ipx);
 													double [] xy_acc = fromPair.times(mpxy).getColumnPackedCopy();
 													// did it get to square of influence?
 													if (    (xy_acc[0] > (ix - 1)) &&
 															(xy_acc[0] < (ix + 1)) &&
 															(xy_acc[1] > (iy - 1)) &&
 															(xy_acc[1] < (iy + 1))) {
 														double w = Math.cos(0.5 * Math.PI * (xy_acc[0] - ix)) *
 																Math.cos(0.5 * Math.PI * (xy_acc[1] - iy));
 														w *= w;
 														if (w >= ignore_contrib) {
 															int pair_indx = ipy * corr_size + ipx + corr_center_offs;
 															contrib_list.add(pair_indx);
 															sumw += w;
 															contrib_weights_list.add(w);
 														}
 													}
 												}
 											}
 											if (sumw > 0) {
 	 											// normalize and store contributions
 		 										int indx = i * mcorr_comb_width + j;
 		 										int ncontrib = contrib_list.size();
 		 										resample_indices[num_pair][indx] = new int    [ncontrib];
 		 										resample_weights[num_pair][indx] = new double [ncontrib];
 		 										for (int icontrib = 0; icontrib < ncontrib; icontrib++) {
 		 											resample_indices[num_pair][indx][icontrib] = contrib_list.get(icontrib);
 		 											resample_weights[num_pair][indx][icontrib] = contrib_weights_list.get(icontrib) / sumw;
 		 										}
 											}
 										}
 									}
 								}
 							} else { // if (scale >= 1.0) { // pair grid is coarser, than accumulated grid
 								for (int i = 0; i < mcorr_comb_height; i++) {
 									int iy = i + mcorr_comb_offset - mcorr_comb_height/2;
 									for (int j = 0; j < mcorr_comb_width; j++) {
 										int ix = j - mcorr_comb_width/2;
 										mxy.set(0, 0,  ix);
 										mxy.set(1, 0,  iy);
 										double [] xy_pair = toPair.times(mxy).getColumnPackedCopy();
// 										if ((num_pair == 62) && (i==7) && (j==7)) {
// 											System.out.println("num_pair="+num_pair+", i="+i+", j="+j);
// 											System.out.println("num_pair="+num_pair+", i="+i+", j="+j);
// 										}
 										// find 4 corners in the pair array
 										if (    (xy_pair[0] >= (1 - transform_size)) &&
 												(xy_pair[0] <= (transform_size - 1)) &&
 												(xy_pair[1] >= (1 - transform_size)) &&
 												(xy_pair[1] <= (transform_size - 1))) {
 	 			 							contrib_list.clear();
 	 			 							contrib_weights_list.clear();
 	 			 							double [] floor_xy = {Math.floor(xy_pair[0]),Math.floor(xy_pair[1])};
 	 			 							double [] ceil_xy =  {Math.ceil(xy_pair[0]), Math.ceil(xy_pair[1])};
 	 			 							int px0 = (int) floor_xy[0] + (transform_size - 1);  
 	 			 							int py0 = (int) floor_xy[1] + (transform_size - 1);  
 	 			 							double sumw = 0.0;
 	 			 							/*
 	 			 							if ((num_pair == 62) && (i==7) && (j==7)) {
 	 			 								System.out.println("Math.floor(xy_pair[0]="+Math.floor(xy_pair[0]));
 	 			 								System.out.println("Math.floor(xy_pair[1]="+Math.floor(xy_pair[1]));
 	 			 								System.out.println("Math.ceil(xy_pair[0]="+Math.ceil(xy_pair[0]));
 	 			 								System.out.println("Math.ceil(xy_pair[1]="+Math.ceil(xy_pair[1]));
 	 			 							}
 	 			 							*/
 	 			 							double wx0 = Math.cos(0.5*Math.PI*(xy_pair[0] - floor_xy[0]));
 	 			 							double wx1 = Math.cos(0.5*Math.PI*(ceil_xy[0] - xy_pair[0]));
 	 			 							if (ceil_xy[0] == floor_xy[0]) wx1 = 0.0;
 	 			 							double wy0 = Math.cos(0.5*Math.PI*(xy_pair[1] - floor_xy[1]));
 	 			 							double wy1 = Math.cos(0.5*Math.PI*(ceil_xy[1] - xy_pair[1]));
 	 			 							if (ceil_xy[1] == floor_xy[1]) wy1 = 0.0;
 	 			 							double [] wxy = {wx0*wy0, wx1*wy0, wx0*wy1, wx1*wy1};
 	 			 							int [] pair_ind = {
 	 			 									py0 * corr_size + px0,
 	 			 									py0 * corr_size + px0 + 1, 
 	 			 									(py0 + 1) * corr_size + px0,
 	 			 									(py0 + 1) * corr_size + px0 + 1};
 	 			 							for (int d = 0; d < wxy.length; d++) {
 	 			 								double w = wxy[d]*wxy[d];
 	 			 								if (w >= ignore_contrib) {
 	 			 									contrib_list.add(pair_ind[d]);
 	 			 									contrib_weights_list.add(w);
 	 			 									sumw += w;
 	 			 								}
 	 			 							}
 											if (sumw > 0) {
 	 											// normalize and store contributions
 		 										int indx = i * mcorr_comb_width + j;
 		 										int ncontrib = contrib_list.size();
 		 										resample_indices[num_pair][indx] = new int    [ncontrib];
 		 										resample_weights[num_pair][indx] = new double [ncontrib];
 		 										for (int icontrib = 0; icontrib < ncontrib; icontrib++) {
 		 											resample_indices[num_pair][indx][icontrib] = contrib_list.get(icontrib);
 		 											resample_weights[num_pair][indx][icontrib] = contrib_weights_list.get(icontrib) / sumw;
 		 										}
 											}
 										}
 									}
 								}
 							}
 						}
 					}
 				}
 			};
 		}		      
 		ImageDtt.startAndJoin(threads);
     }
  
    public Correlation2d ( // USED in lwir
    		int                numSensors,
    		ImageDttParameters imgdtt_params,
    		int                transform_size,
    		double             wndx_scale, // (wndy scale is always 1.0)
    		boolean            monochrome,
    		boolean            debug) {
    	this.numSensors = numSensors;
    	this.monochrome = monochrome;
    	this.dtt = new DttRad2(transform_size);
    	this.dtt.setup_CSIIe(transform_size); // common for multiple threads, initializing on demand fails
    	this.transform_size = transform_size;
    	this.transform_len = transform_size * transform_size;
    	this.corr_size = transform_size * 2 -1;
    	// not initialized until needed
    	this.transpose_all_ortho =     new int [corr_size*corr_size];
    	this.transpose_all_diagonal =  new int [corr_size*corr_size];
    	this.ortho_notch_filter = new double [corr_size];
    	this.corr_wndy = halfFlatTopWindow(
				imgdtt_params.corr_wndy_size,   // int     ihwidth,
				imgdtt_params.corr_wndy_hwidth, // double  hwidth,
				imgdtt_params.corr_wndy_blur,   // double  blur,
				true,                           // boolean normalize,
				false,                          // boolean notch,
				1.0); // double  scale);
    	this.corr_wndx = halfFlatTopWindow(
				imgdtt_params.corr_wndx_size,   // int     ihwidth,
				imgdtt_params.corr_wndx_hwidth, //  double  hwidth,
				imgdtt_params.corr_wndx_blur,   // double  blur,
				true,                           // boolean normalize,
				false,                          // boolean notch,
				wndx_scale);                    // double  scale);
    	this.corr_wndy_notch = halfFlatTopWindow(
				imgdtt_params.corr_strip_notch, // int     ihwidth,
				imgdtt_params.corr_notch_hwidth,//  double  hwidth,
				imgdtt_params.corr_notch_blur,  // double  blur,
				true,                            // boolean normalize,
				true,                            // boolean notch,
				wndx_scale);                     // double  scale);
    	int num_pairs = getNumPairs(numSensors);//  numSensors * (numSensors-1) /2;
    	pair_start_end = new int [num_pairs][2];
    	pair_orient =    new int [num_pairs][num_pairs];
    	pair_length =    new int [num_pairs];
    	top_is_0 = imgdtt_params.getTopIs0(numSensors); // sensor 0 is straight up, false half-rotated 
    	cor_titles =       new String[num_pairs];
    	cor_titles_combo = new String[CORR_TITLES_EXTRA.length];
    	setupPairs();
      }

    public Correlation2d ( // USED in lwir
    		int     numSensors,
    		int     transform_size,
    		boolean monochrome,
    		boolean debug) {
    	this.numSensors = numSensors;
    	this.monochrome = monochrome;
    	this.dtt = new DttRad2(transform_size);
    	this.dtt.setup_CSIIe(transform_size); // common for multiple threads, initializing on demand fails
    	this.transform_size = transform_size;
    	this.transform_len = transform_size * transform_size;
    	this.corr_size = transform_size * 2 -1;
    	// not initialized until needed
    	this.transpose_all_ortho =     new int [corr_size*corr_size];
    	this.transpose_all_diagonal =  new int [corr_size*corr_size];
    	this.ortho_notch_filter =      new double [corr_size];
    	this.corr_wndy =               null; // will not be used
    	this.corr_wndx =               null; // will not be used
    	this.corr_wndy_notch =         null; // will not be used
    	int num_pairs = numSensors * (numSensors-1) /2;
    	pair_start_end = new int [num_pairs][2];
    	pair_orient =    new int [num_pairs][num_pairs];
    	pair_length =    new int [num_pairs];
    	top_is_0 = numSensors > 4; // sensor 0 is straight up, false half-rotated
//     	top_is_0 = imgdtt_params.getTopIs0(numSensors); // sensor 0 is straight up, false half-rotated 
    	cor_titles =       new String[num_pairs];
    	cor_titles_combo = new String[CORR_TITLES_EXTRA.length];
    	setupPairs();
    }
    
      public int [] getTransposeAll(boolean diagonal){ // USED in lwir
    	  if (diagonal) return getTransposeAllDiagonal();
    	  else          return getTransposeAllOrtho();
      }

      public int [] getTransposeAllOrtho(){ // USED in lwir
    	  if (this.transpose_all_ortho[0] == this.transpose_all_ortho[1]) {
    		  for (int i =0; i < corr_size; i++){
    			  for (int j =0; j < corr_size; j++){
    				  this.transpose_all_ortho[i * corr_size + j] = j * corr_size + i;
    			  }
    		  }
    	  }
    	  return this.transpose_all_ortho;
      }

      public int [] getTransposeAllDiagonal(){ // USED in lwir
    	  if (this.transpose_all_diagonal[0] == this.transpose_all_diagonal[1]) {
    		  for (int i =0; i < corr_size; i++){
    			  for (int j = 0; j < corr_size; j++){
    				  this.transpose_all_diagonal[i * corr_size + j] = (corr_size - i -1) * corr_size + j;
    			  }
    		  }
    	  }
    	  return this.transpose_all_diagonal;
      }

      /**
       * Multiply CLT data of two channels, OK with null inputs (missing colors for monochrome images)
       * @param clt_data1 first operand FD CLT data[4][transform_len]
       * @param clt_data2 second operand FD CLT data[4][transform_len]
       * @return [4][transform_len] FD CLT data
       */
      public double[][] correlateSingleColorFD(
      		double [][] clt_data1,
      		double [][] clt_data2,
      		double [][] tcorr){ // null or initialized to [4][transform_len]

      	if (tcorr == null) tcorr = new double [4][transform_len];
    	if ((clt_data1 == null) || (clt_data1 == null)) return null; // to work with missing colors for monochrome
  		for (int i = 0; i < transform_len; i++) {
  			for (int n = 0; n<4; n++){
  				tcorr[n][i] = 0;
  				for (int k=0; k<4; k++){
  					if (ZI[n][k] < 0)
  						tcorr[n][i] -=
  								clt_data1[-ZI[n][k]][i] * clt_data2[k][i];
  					else
  						tcorr[n][i] +=
  								clt_data1[ZI[n][k]][i] * clt_data2[k][i];
  				}
  			}
  		}
  		return tcorr;
      }

      /**
       * Normalize 2D correlation in FD, LPF (if not null) and convert to pixel domain and trim
       * @param tcorr FD representation of the correlation[4][64]
       * @param lpf LPF [64] or null
       * @param afat_zero2 fat zero to add during normalization, units of squared values
       * @param corr_radius if >=0 and < 7 - extract only the central part of the 15x15 square
       * @return 2D phase correlation in linescan order
       */
      public double[] normalizeConvertCorr(
    		  double [][] tcorr, // null or initialized to [4][transform_len]
    		  double []   lpf,
    		  double      afat_zero2, // absolute fat zero, same units as components squared values
    		  int corr_radius,
    		  boolean debug_gpu){
    	  if (tcorr == null) return null;
    	  double afat_zero4 = afat_zero2*afat_zero2;

    	  for (int i = 0; i < transform_len; i++) {
    		  double s = afat_zero4;
    		  for (int n = 0; n< 4; n++){
    			  s += tcorr[n][i]*tcorr[n][i];
    		  }
    		  double k = 1.0/ Math.sqrt(s);
    		  for (int n = 0; n< 4; n++){
    			  tcorr[n][i]*= k;
    		  }
    	  }
    	  if (debug_gpu) {
    		  System.out.println("=== NORMALIZED CORRELATION , afat_zero2="+afat_zero2+", afat_zero4="+afat_zero4+" ===");
    		  for (int dct_mode = 0; dct_mode < 4; dct_mode++) {
    			  System.out.println("------dct_mode="+dct_mode);
    			  for (int i = 0; i < transform_size; i++) {
    				  for (int j = 0; j < transform_size; j++) {
    					  System.out.print(String.format("%10.5f ", tcorr[dct_mode][transform_size * i + j]));
    				  }
    				  System.out.println();
    			  }
    		  }
    	  }
    	  if (lpf != null) {
        	  if (debug_gpu) {
        		  System.out.println("=== LPF for CORRELATION ===");
/*
        		  for (int i = 0; i < transform_size; i++) {
        			  System.out.print("\t\t");
        			  for (int j = 0; j < transform_size; j++) {
        				  System.out.print(String.format("%10.8ff", lpf[transform_size * i + j]));
        				  if ((j < (transform_size-1)) || (j < (transform_size-1))){
        					  System.out.print(", ");
        				  }
        			  }
        			  System.out.println();
        		  }
*/
      			System.out.print("__constant__ float lpf_corr[64]={");
    			for (int i=0; i<lpf.length;i++){
    				System.out.print(String.format("%10.8ff", lpf[i]));
    				if (i == 63) {
    					System.out.println("};");
    				} else {
    					System.out.print(", ");
    					if ((i % 8) == 7) {
    						System.out.print("\n                                 ");
    					}
    				}
    			}

        	  }
    		  for (int n = 0; n<4; n++) {
    			  for (int i = 0; i < transform_len; i++) {
    				  tcorr[n][i] *= lpf[i];
    			  }
    		  }
    	  }
    	  if (debug_gpu) {
    		  System.out.println("=== LPF-ed CORRELATION ===");
    		  for (int dct_mode = 0; dct_mode < 4; dct_mode++) {
    			  System.out.println("------dct_mode="+dct_mode);
    			  for (int i = 0; i < transform_size; i++) {
    				  for (int j = 0; j < transform_size; j++) {
    					  System.out.print(String.format("%10.5f ", tcorr[dct_mode][transform_size * i + j]));
    				  }
    				  System.out.println();
    			  }
    		  }
    	  }
		  for (int quadrant = 0; quadrant < 4; quadrant++){
			  int mode = ((quadrant << 1) & 2) | ((quadrant >> 1) & 1); // transpose
			  tcorr[quadrant] = dtt.dttt_iie(tcorr[quadrant], mode, transform_size, debug_gpu); // not orthogonal, term[0] is NOT *= 1/sqrt(2)
		  }
    	  if (debug_gpu) {
    		  System.out.println("=== CONVERTED CORRELATION ===");
    		  for (int dct_mode = 0; dct_mode < 4; dct_mode++) {
    			  System.out.println("------dct_mode="+dct_mode);
    			  for (int i = 0; i < transform_size; i++) {
    				  for (int j = 0; j < transform_size; j++) {
    					  System.out.print(String.format("%10.5f ", tcorr[dct_mode][transform_size * i + j]));
    				  }
    				  System.out.println();
    			  }
    		  }
    	  }

    	  // convert from 4 quadrants to 15x15 centered tiles (only composite)
    	  double [] corr_pd =  dtt.corr_unfold_tile(tcorr,	transform_size);
    	  if (debug_gpu) {
    		  int corr_size = 2* transform_size -1;
    		  System.out.println("=== UNFOLDED CORRELATION ===");
    		  for (int i = 0; i < corr_size; i++) {
    			  for (int j = 0; j < corr_size; j++) {
    				  System.out.print(String.format("%10.5f ", corr_pd[corr_size * i + j]));
    			  }
    			  System.out.println();
    		  }
    	  }

    	  if ((corr_radius <= 0) || (corr_radius >= (transform_size - 1))) {
    		  return corr_pd;
    	  }
    	  int full_size = 2 * transform_size - 1;
    	  int trimmed_size = 2 * corr_radius + 1;
    	  int trim =  transform_size - 1 - corr_radius;
    	  double [] trimmed_pd = new double [trimmed_size * trimmed_size];
    	  int ioffs = (full_size + 1)*trim;
    	  for (int orow = 0; orow < trimmed_size; orow++) {
    		  System.arraycopy(corr_pd, orow*full_size + ioffs, trimmed_pd, orow*trimmed_size, trimmed_size);
    	  }
    	  return trimmed_pd;
      }


    /**
     * Multiply CLT data of two channels, normalize amplitude, OK with null inputs (missing colors for monochrome images)
     * @param clt_data1 first operand FD CLT data[4][transform_len]
     * @param clt_data2 second operand FD CLT data[4][transform_len]
     * @param fat_zero add to normalization amplitude
     * @return [4][transform_len] FD CLT data
     */
	public double[][] correlateSingleColorFD( // USED in lwir
    		double [][] clt_data1,
    		double [][] clt_data2,
    		double [][] tcorr, // null or initialized to [4][transform_len]
    		double      fat_zero) {

    	if (tcorr == null) tcorr = new double [4][transform_len];
    	if ((clt_data1 == null) || (clt_data1 == null)) return null; // to work with missing colors for monochrome
    	double [] a2 = new double[transform_len];
    	double sa2 = 0.0;
		for (int i = 0; i < transform_len; i++) {
			double s1 = 0.0, s2=0.0;
			for (int n = 0; n< 4; n++){
				s1+=clt_data1[n][i] * clt_data1[n][i];
				s2+=clt_data2[n][i] * clt_data2[n][i];
			}
			a2[i] = Math.sqrt(s1*s2);
			sa2 += a2[i];
		}
		double fz2 = sa2/transform_len * fat_zero * fat_zero; // fat_zero squared to match units
		for (int i = 0; i < transform_len; i++) {
			double scale = 1.0 / (a2[i] + fz2);
			for (int n = 0; n<4; n++){
				tcorr[n][i] = 0;
				for (int k=0; k<4; k++){
					if (ZI[n][k] < 0)
						tcorr[n][i] -=
								clt_data1[-ZI[n][k]][i] * clt_data2[k][i];
					else
						tcorr[n][i] +=
								clt_data1[ZI[n][k]][i] * clt_data2[k][i];
				}
				tcorr[n][i] *= scale;
			}
		}
		return tcorr;
    }

    
    /**
     * Calculate color channels FD phase correlations, mix results with weights, apply optional low-pass filter
     * and convert to the pixel domain  as [(2*transform_size-1) * (2*transform_size-1)] tiles (15x15)
     * No transposing or rotation
     * @param clt_data1 [3][4][transform_len] first operand data. First index - RBG color
     * @param clt_data2 [3][4][transform_len] first operand data. First index - RBG color
     * @param lpf   optional [transform_len] LPF filter data
     * @param scale_value   scale correlation results to compensate for lpf changes and other factors
     * @param col_weights [3] - color weights {R, B, G} - green is last, normalized to sum =1.0
     * @param fat_zero fat zero for phase correlation (0 seems to be OK)
     * @return correlation result [(2*transform_size-1) * (2*transform_size-1)]
     */
    public double[]       correlateCompositeFD( // USED in lwir
    		double [][][] clt_data1,
    		double [][][] clt_data2,
    		double []     lpf,
    		double        scale_value, // scale correlation value
    		double []     col_weights_in, // should have the same dimension as clt_data1 and clt_data2
    		double        fat_zero) {

//    	if ((clt_data1 == null) || (clt_data1 == null)) return null;
    	// work with sparse clt
    	double [] col_weights = col_weights_in.clone();
    	double s = 0.0;
    	for (int i = 0; i < col_weights.length; i++) {
    		if ((clt_data1[i] == null) || (clt_data2[i] == null)) {
    			col_weights[i]= 0.0;
    		}
    		s+=col_weights[i];
    	}
    	for (int i = 0; i < col_weights.length; i++) {
    		if (col_weights[i] != 0.0) col_weights[i] *= scale_value;
    	}

    	if (clt_data1.length == 1) { // monochrome  // not used in lwir
    		col_weights = new double[1];
    		col_weights[0] = 1.0;
    	}
    	for (int i = 0; i < col_weights.length; i++) {
    		if (col_weights[i] != 0.0) col_weights[i]/=s; // will have 1.0 for the single color
    	}


    	double [][][]tcorr = new double [clt_data1.length][4][transform_len];
    	int first_col = -1;
    	for (int col = 0; col < tcorr.length; col++) if (col_weights[col] > 0.0 ) {
    		 correlateSingleColorFD(
    		    		clt_data1[col],
    		    		clt_data2[col],
    		    		tcorr[col],
    		    		fat_zero/scale_value);

    		 if (first_col < 0) {// accumulate all channels in first non-null color ( 0 for color, 2 for mono?)
    			 first_col = col; // first non-empty color (2, green) or 0 for color images
    			 for (int n = 0; n < 4; n++) {
    				 for (int i = 0; i < transform_len; i++) {
    					 tcorr[first_col][n][i] *= col_weights[col];
    				 }
    			 }
    		 } else {
    			 for (int n = 0; n < 4; n++) {
    				 for (int i = 0; i < transform_len; i++) {
    					 tcorr[first_col][n][i] += tcorr[col][n][i] * col_weights[col];
    				 }
    			 }
    		 }
    	}
    	if (lpf != null) {
    		for (int n = 0; n<4; n++) {
    			for (int i = 0; i < transform_len; i++) {
    				tcorr[first_col][n][i] *= lpf[i];
    			}
    		}
    	}

    	for (int quadrant = 0; quadrant < 4; quadrant++){
    		int mode = ((quadrant << 1) & 2) | ((quadrant >> 1) & 1); // transpose
    		tcorr[first_col][quadrant] = dtt.dttt_iie(tcorr[first_col][quadrant], mode, transform_size);
    	}
		// convert from 4 quadrants to 15x15 centered tiles (only composite)
    	double [] corr_pd =  dtt.corr_unfold_tile(tcorr[first_col],	transform_size);
    	return corr_pd;
    }
    
    
    /**
     * Combine compatible correlations in TD (for 6 pairs requires 2 calls): vertical will be transposed,
     * (all quadrants transposed, quadrants 1 and 2 - swapped), "diagonal other" (#5) will be flipped
     * vertically (quadrants 2 and 3 - negated) 
     * @param corr_combo_td per-pair array of the TD representation of correlation pairs, mixed color components
     * @param pairs_mask bitmask of the pairs to combinde
     * @return TD representation of the combined 2D correlation
     */
    public double [][] combineCorrelationsTD(
    		double [][][] corr_combo_td, // vertical will be transposed, other diagonal flipped vertically (should be separate)
    		int           pairs_mask  // vertical
    		) {
    	return combineCorrelationsTD(corr_combo_td, pairs_mask, false);
    }

    /**
     * Combine compatible correlations in TD (for 6 pairs requires 2 calls): vertical will be transposed,
     * (all quadrants transposed, quadrants 1 and 2 - swapped), "diagonal other" (#5) will be flipped
     * vertically (quadrants 2 and 3 - negated) 
     * @param corr_combo_td per-pair array of the TD representation of correlation pairs, mixed color components
     * @param pairs_mask bitmask of the pairs to combinde
     * @param no_transpose_vertical Do not transpose vertical pairs (used when combining just vertical together)
     * @return TD representation of the combined 2D correlation
     */
    public double [][] combineCorrelationsTD(
    		double [][][] corr_combo_td, // vertical will be transposed, other diagonal flipped vertically (should be separate)
    		int           pairs_mask, // vertical
    		boolean       no_transpose_vertical
    		) {
    	if (corr_combo_td != null) {
    		double [][]tcorr = new double [4][transform_len];
    		int num_combined = 0;
    		for (int npair = 0; npair < corr_combo_td.length; npair++) if ((((pairs_mask >> npair) & 1) != 0 ) && (corr_combo_td[npair] !=null)) {
    			if (isHorizontalPair(npair) || isDiagonalMainPair(npair) || (isVerticalPair(npair) && no_transpose_vertical)) {
    				for (int q = 0; q < 4; q++) {
    					for (int i = 0; i < transform_len; i++) {
    						tcorr[q][i] += corr_combo_td[npair][q][i];
    					}
    				}
    				num_combined++;
    			} else if (isVerticalPair(npair)) {
    				for (int q = 0; q < 4; q++) {
    					int qr = ((q & 1) << 1) | ((q >> 1) & 1);
    					int i = 0;
    					for (int iy = 0; iy < transform_size; iy++ ) {
        					for (int ix = 0; ix < transform_size; ix++ ) {
        						int ir = ix*transform_size + iy;
        						tcorr[qr][ir] += corr_combo_td[npair][q][i];
        						i++;
        					}
    					}
    				}
    				num_combined++;
    			} else if (isDiagonalOtherPair(npair)) {
    				for (int q = 0; q < 2; q++) {
    					for (int i = 0; i < transform_len; i++) {
    						tcorr[q][i] += corr_combo_td[npair][q][i];
    					}
    				}
    				for (int q = 2; q < 4; q++) {
    					for (int i = 0; i < transform_len; i++) {
    						tcorr[q][i] -= corr_combo_td[npair][q][i];
    					}
    				}
    				num_combined++;
    			}
    		}
    		if (num_combined > 0) {
				for (int q = 0; q < 4; q++) {
					for (int i = 0; i < transform_len; i++) {
						tcorr[q][i] /= num_combined;
					}
				}
    			return tcorr;
    		}
    	}
    	return null;
    }
    
/**
 * Normalize TD tile (before converting to pixel domain as a phase correlation result and optionally LOF
 * @param td 4-quadrant TD tile, will be normalized in-place
 * @param lpf optional low-pass filter (or null)
 * @param fat_zero relative fat zero to apply during normalization
 */
    public void normalize_TD(
    		double [][] td,
    		double []   lpf, // or null
    		double      fat_zero) {
    	if (td == null)    		return;
    	double [] a2 = new double[transform_len];
    	double sa2 = 0.0;
		for (int i = 0; i < transform_len; i++) {
			double s1 = 0.0;
			for (int n = 0; n< 4; n++){
				s1+=td[n][i] * td[n][i];
			}
			a2[i] = Math.sqrt(s1);
			sa2 += a2[i];
		}
		double fz2 = sa2/transform_len * fat_zero * fat_zero; // fat_zero squared to match units
		for (int i = 0; i < transform_len; i++) {
			double scale = 1.0 / (a2[i] + fz2);
			for (int n = 0; n<4; n++){
				td[n][i] *= scale;
			}
		}
    	if (lpf != null) {
    		for (int n = 0; n<4; n++) {
    			for (int i = 0; i < transform_len; i++) {
    				td[n][i] *= lpf[i];
    			}
    		}
    	}
    }
    
    /**
     * Convert 2D correlation tile from the TD to pixel domain
     * @param td 4-quadrant TD representation of the tile
     * @return normally a 15x15 pixel 2D correlation
     */
    public double [] convertCorrToPD(
    		double [][] td) {
    	if (td == null)	return null;
    	for (int quadrant = 0; quadrant < 4; quadrant++){
    		int mode = ((quadrant << 1) & 2) | ((quadrant >> 1) & 1); // transpose
    		td[quadrant] = dtt.dttt_iie(td[quadrant], mode, transform_size);
    	}
		// convert from 4 quadrants to 15x15 centered tiles (only composite)
    	double [] corr_pd =  dtt.corr_unfold_tile(td,	transform_size);
    	return corr_pd;
    }
    
    
    /**
     * Calculate all required image pairs phase correlation, stay in Transform Domain
     * @param clt_data aberration-corrected FD CLT data [camera][color][tileY][tileX][quadrant][index]
     * @param tileX tile to extract X index
     * @param tileY tile to extract Y index
     * @param pairs_mask bimask of required pairs
     * @param lpf_rb optional low-pass filter - extra LPF for red and blue
     * @param scale_value   scale correlation results to compensate for lpf changes and other factors
     * @param col_weights RBG color weights
     * @return [pair][quadrant][index]
     */
    @Deprecated
    public double [][][]  correlateCompositeTD(
    		double [][][][][][] clt_data,
    		int                 tileX,
    		int                 tileY,
    		int                 pairs_mask,
    		double []           lpf_rb, // extra lpf for red and blue (unused for mono) or null
    		double              scale_value, // scale correlation value
    		double []           col_weights) {
    	double [][][][]     clt_data_tile = new double[clt_data.length][][][]; // [camera][color][quadrant][index]
    	for (int ncam = 0; ncam < clt_data.length; ncam++) if (clt_data[ncam] != null){
    		clt_data_tile[ncam] = new double[clt_data[ncam].length][][];
        	for (int ncol = 0; ncol < clt_data[ncam].length; ncol++) if ((clt_data[ncam][ncol] != null) && (clt_data[ncam][ncol][tileY] != null)){
        		clt_data_tile[ncam][ncol] = clt_data[ncam][ncol][tileY][tileX];
        	}
    	}
    	return correlateCompositeTD(
    		    		clt_data_tile,
    		    		pairs_mask, // already decoded so bit 0 - pair 0
    		    		lpf_rb,
    		    		scale_value,
    		    		col_weights);
    }

    /**
     * Calculate all required image pairs phase correlation, stay in Transform Domain
     * @param clt_data aberration-corrected FD CLT data [camera][color][tileY][tileX][quadrant][index]
     * @param tileX tile to extract X index
     * @param tileY tile to extract Y index
     * @param pairs_mask bimask of required pairs
     * @param lpf_rb optional low-pass filter - extra LPF for red and blue
     * @param scale_value   scale correlation results to compensate for lpf changes and other factors
     * @param col_weights RBG color weights
     * @return [pair][quadrant][index]
     */
    public double [][][]  correlateCompositeTD(
    		double [][][][][][] clt_data,
    		int                 tileX,
    		int                 tileY,
    		boolean []          pairs_mask,
    		double []           lpf_rb, // extra lpf for red and blue (unused for mono) or null
    		double              scale_value, // scale correlation value
    		double []           col_weights) {
    	double [][][][]     clt_data_tile = new double[clt_data.length][][][]; // [camera][color][quadrant][index]
    	for (int ncam = 0; ncam < clt_data.length; ncam++) if (clt_data[ncam] != null){
    		clt_data_tile[ncam] = new double[clt_data[ncam].length][][];
        	for (int ncol = 0; ncol < clt_data[ncam].length; ncol++) if ((clt_data[ncam][ncol] != null) && (clt_data[ncam][ncol][tileY] != null)){
        		clt_data_tile[ncam][ncol] = clt_data[ncam][ncol][tileY][tileX];
        	}
    	}
    	return correlateCompositeTD(
    		    		clt_data_tile,
    		    		pairs_mask, // already decoded so bit 0 - pair 0
    		    		lpf_rb,
    		    		scale_value,
    		    		col_weights);
    }
    
    
    /**
     * Calculate all required image pairs phase correlation, stay in Transform Domain
     * @param clt_data_tile aberration-corrected FD CLT data for one tile [camera][color][quadrant][index]
     * @param pairs_mask bimask of required pairs NOW USE boolean array and new pairs
     * @param lpf optional final low-pass filter 
     * @param lpf_rb optional low-pass filter (extra) for R,B components
     * @param scale_value   scale correlation results to compensate for lpf changes and other factors
     * @param col_weights RBG color weights
     * @return [pair][quadrant][index]
     */
    @Deprecated
    public double [][][]  correlateCompositeTD(
    		double [][][][]     clt_data_tile,
    		int                 pairs_mask, // already decoded so bit 0 - pair 0
    		double []           lpf_rb, // extra lpf for red and blue (unused for mono) or null
    		double              scale_value, // scale correlation value
    		double []           col_weights) {
    	if (clt_data_tile == null) return null;
    	double [][][] pairs_corr = new double [PAIRS.length][][];
    	for (int npair = 0; npair < pairs_corr.length; npair++) if (((pairs_mask >> npair) & 1) != 0 ) {
    		int ncam1 = PAIRS[npair][0];
    		int ncam2 = PAIRS[npair][1];
    		if ((ncam1 < clt_data_tile.length) && (clt_data_tile[ncam1] != null) && (ncam2 < clt_data_tile.length) && (clt_data_tile[ncam2] != null)) {
    			pairs_corr[npair] =  correlateCompositeTD(
    					clt_data_tile[ncam1], // double [][][] clt_data1,
    					clt_data_tile[ncam2], // double [][][] clt_data2,
    		    		lpf_rb,               // double []     lpf_rb,
    		    		scale_value,
    		    		col_weights);         // double []     col_weights,
    		}
    	}
    	return pairs_corr;
    }

    
    /**
     * Calculate all required image pairs phase correlation, stay in Transform Domain
     * @param clt_data_tile aberration-corrected FD CLT data for one tile [camera][color][quadrant][index]
     * @param pairs_mask bimask of required pairs
     * @param lpf optional final low-pass filter 
     * @param lpf_rb optional low-pass filter (extra) for R,B components
     * @param scale_value   scale correlation results to compensate for lpf changes and other factors
     * @param col_weights RBG color weights
     * @return [pair][quadrant][index]
     */
    public double [][][]  correlateCompositeTD(
    		double [][][][]     clt_data_tile,
    		boolean[]           pairs_mask, // already decoded so bit 0 - pair 0
    		double []           lpf_rb, // extra lpf for red and blue (unused for mono) or null
    		double              scale_value, // scale correlation value
    		double []           col_weights) {
    	if (clt_data_tile == null) return null;
    	double [][][] pairs_corr = new double [getNumPairs()][][];
    	for (int npair = 0; npair < pairs_corr.length; npair++) if (pairs_mask[npair]) {
    		int [] pair = getPair(npair);
    		int ncam1 = pair[0]; // start 
    		int ncam2 = pair[1]; // end
    		if ((ncam1 < clt_data_tile.length) && (clt_data_tile[ncam1] != null) && (ncam2 < clt_data_tile.length) && (clt_data_tile[ncam2] != null)) {
    			pairs_corr[npair] =  correlateCompositeTD(
    					clt_data_tile[ncam1], // double [][][] clt_data1,
    					clt_data_tile[ncam2], // double [][][] clt_data2,
    		    		lpf_rb,               // double []     lpf_rb,
    		    		scale_value,
    		    		col_weights);         // double []     col_weights,
    		}
    	}
    	return pairs_corr;
    }
    
    
    /**
     * Calculate color channels FD phase correlations, mix results with weights, apply optional low-pass filter
     * No transposing or rotation
     * @param clt_data1 [3][4][transform_len] first operand data. First index - RBG color
     * @param clt_data2 [3][4][transform_len] first operand data. First index - RBG color
     * @param lpf optional final low-pass filter applied to mixed colors
     * @param lpf_rb optional low-pass filter (extra) for R,B components only
     * @param scale_value scale correlation results to compensate for lpf changes and other factors
     * @param col_weights [3] - color weights {R, B, G} - green is last, normalized to sum =1.0 
     * @return correlation result (Transform Domain) [quadrant][index]
     */
    public double[][]   correlateCompositeTD(
    		double [][][] clt_data1,
    		double [][][] clt_data2,
    		double []     lpf_rb, // extra lpf for red and blue (unused for mono) or null
    		double        scale_value, // scale correlation value
    		double []     col_weights_in){ // should have the same dimension as clt_data1 and clt_data2
    	double [] col_weights = col_weights_in.clone();
    	double s = 0.0;
    	for (int i = 0; i < col_weights.length; i++) {
    		if ((clt_data1[i] == null) || (clt_data2[i] == null)) {
    			col_weights[i]= 0.0;
    		}
    		s+=col_weights[i];
    	}
    	for (int i = 0; i < col_weights.length; i++) {
    		if (col_weights[i] != 0.0) col_weights[i] *= scale_value;
    	}

    	if (clt_data1.length == 1) { // monochrome  // not used in lwir
    		col_weights = new double[1];
    		col_weights[0] = 1.0;
    	}
    	for (int i = 0; i < col_weights.length; i++) {
    		if (col_weights[i] != 0.0) col_weights[i]/=s; // will have 1.0 for the single color
    	}


    	double [][][]tcorr = new double [clt_data1.length][4][transform_len];
    	int first_col = -1;
    	for (int col = 0; col < tcorr.length; col++) if (col_weights[col] > 0.0 ) {
    		 correlateSingleColorFD( // no amplitude normalization
    		    		clt_data1[col],
    		    		clt_data2[col],
    		    		tcorr[col]);

    		 if (first_col < 0) {// accumulate all channels in first non-null color ( 0 for color, 2 for mono?)
    			 first_col = col; // first non-empty color (2, green) or 0 for color images
    			 for (int n = 0; n < 4; n++) {
    				 for (int i = 0; i < transform_len; i++) {
    					 tcorr[first_col][n][i] *= col_weights[col];
    				 }
    			 }
    		 } else {
    			 for (int n = 0; n < 4; n++) {
    				 for (int i = 0; i < transform_len; i++) {
    					 tcorr[first_col][n][i] += tcorr[col][n][i] * col_weights[col];
    				 }
    			 }
    			 if (col == 1) { // after blue, before red
    				 if (lpf_rb != null) {
    					 for (int n = 0; n<4; n++) {
    						 for (int i = 0; i < transform_len; i++) {
    							 tcorr[first_col][n][i] *= lpf_rb[i];
    						 }
    					 }
    				 }
    			 }
    		 }
    	}
    	return tcorr[first_col];
    }
    
    /**
     * Calculate all required image pairs phase correlation
     * @param clt_data aberration-corrected FD CLT data [camera][color][tileY][tileX][quadrant][index]
     * @param tileX tile to extract X index
     * @param tileY tile to extract Y index
     * @param pairs_mask bimask of required pairs
     * @param lpf optional low-pass filter
     * @param scale_value   scale correlation results to compensate for lpf changes and other factors
     * @param col_weights RBG color weights
     * @param fat_zero fat zero for phase correlations
     * @return [pair][corr_index]
     */
    public double [][]  correlateCompositeFD( // USED in lwir
    		double [][][][][][] clt_data,
    		int                 tileX,
    		int                 tileY,
    		int                 pairs_mask,
    		double []           lpf,
    		double              scale_value, // scale correlation value
    		double []           col_weights,
    		double              fat_zero) {
    	double [][][][]     clt_data_tile = new double[clt_data.length][][][]; // [camera][color][quadrant][index]
    	for (int ncam = 0; ncam < clt_data.length; ncam++) if (clt_data[ncam] != null){
    		clt_data_tile[ncam] = new double[clt_data[ncam].length][][];
        	for (int ncol = 0; ncol < clt_data[ncam].length; ncol++) if ((clt_data[ncam][ncol] != null) && (clt_data[ncam][ncol][tileY] != null)){
        		clt_data_tile[ncam][ncol] = clt_data[ncam][ncol][tileY][tileX];
        	}
    	}
    	return correlateCompositeFD(
    		    		clt_data_tile,
    		    		pairs_mask, // already decoded so bit 0 - pair 0
    		    		lpf,
    		    		scale_value,
    		    		col_weights,
    		    		fat_zero);
    }

    /**
     * Calculate all required image pairs phase correlation
     * @param clt_data aberration-corrected FD CLT data [camera][color][tileY][tileX][quadrant][index]
     * @param tileX tile to extract X index
     * @param tileY tile to extract Y index
     * @param pairs_mask boolean array of required pairs
     * @param lpf optional low-pass filter
     * @param scale_value   scale correlation results to compensate for lpf changes and other factors
     * @param col_weights RBG color weights
     * @param fat_zero fat zero for phase correlations
     * @return [pair][corr_index]
     */
    public double [][]  correlateCompositeFD( // USED in lwir
    		double [][][][][][] clt_data,
    		int                 tileX,
    		int                 tileY,
    		boolean[]           pairs_mask,
    		double []           lpf,
    		double              scale_value, // scale correlation value
    		double []           col_weights,
    		double              fat_zero) {
    	double [][][][]     clt_data_tile = new double[clt_data.length][][][]; // [camera][color][quadrant][index]
    	for (int ncam = 0; ncam < clt_data.length; ncam++) if (clt_data[ncam] != null){
    		clt_data_tile[ncam] = new double[clt_data[ncam].length][][];
        	for (int ncol = 0; ncol < clt_data[ncam].length; ncol++) if ((clt_data[ncam][ncol] != null) && (clt_data[ncam][ncol][tileY] != null)){
        		clt_data_tile[ncam][ncol] = clt_data[ncam][ncol][tileY][tileX];
        	}
    	}
    	return correlateCompositeFD(
    		    		clt_data_tile,
    		    		pairs_mask,
    		    		lpf,
    		    		scale_value,
    		    		col_weights,
    		    		fat_zero);
    }
    
    
    /**
     * Calculate all required image pairs phase correlation
     * @param clt_data_tile aberration-corrected FD CLT data for one tile [camera][color][quadrant][index]
     * @param pairs_mask bimask of required pairs
     * @param lpf optional low-pass filter
     * @param scale_value   scale correlation results to compensate for lpf changes and other factors
     * @param col_weights RBG color weights
     * @param fat_zero fat zero for phase correlations
     * @return [pair][corr_index]
     */
    public double [][]  correlateCompositeFD( // USED in lwir
    		double [][][][]     clt_data_tile,
    		int                 pairs_mask, // already decoded so bit 0 - pair 0
    		double []           lpf,
    		double              scale_value, // scale correlation value
    		double []           col_weights,
    		double              fat_zero) {
    	if (clt_data_tile == null) return null;
    	double [][] pairs_corr = new double [PAIRS.length][];
    	for (int npair = 0; npair < pairs_corr.length; npair++) if (((pairs_mask >> npair) & 1) != 0 ) {
    		int ncam1 = PAIRS[npair][0];
    		int ncam2 = PAIRS[npair][1];
    		if ((ncam1 < clt_data_tile.length) && (clt_data_tile[ncam1] != null) && (ncam2 < clt_data_tile.length) && (clt_data_tile[ncam2] != null)) {
    			pairs_corr[npair] =  correlateCompositeFD(
    					clt_data_tile[ncam1], // double [][][] clt_data1,
    					clt_data_tile[ncam2], // double [][][] clt_data2,
    		    		lpf,                  // double []     lpf,
    		    		scale_value,
    		    		col_weights,          // double []     col_weights,
    		    		fat_zero);            // double        fat_zero)
    		}
    	}
    	return pairs_corr;
    }
    
    /**
     * Calculate all required image pairs phase correlation
     * @param clt_data_tile aberration-corrected FD CLT data for one tile [camera][color][quadrant][index]
     * @param pairs_mask boolean array of required pairs
     * @param lpf optional low-pass filter
     * @param scale_value   scale correlation results to compensate for lpf changes and other factors
     * @param col_weights RBG color weights
     * @param fat_zero fat zero for phase correlations
     * @return [pair][corr_index]
     */

    public double [][]  correlateCompositeFD( // USED in lwir
    		double [][][][]     clt_data_tile,
    		boolean []          pairs_mask, 
    		double []           lpf,
    		double              scale_value, // scale correlation value
    		double []           col_weights,
    		double              fat_zero) {
    	if (clt_data_tile == null) return null;
    	double [][] pairs_corr = new double [getNumPairs()][];
    	for (int npair = 0; npair < pairs_corr.length; npair++) if (pairs_mask[npair]) {
    		int [] pair = getPair(npair);
    		int ncam1 = pair[0]; // start 
    		int ncam2 = pair[1]; // end
    		if ((ncam1 < clt_data_tile.length) && (clt_data_tile[ncam1] != null) && (ncam2 < clt_data_tile.length) && (clt_data_tile[ncam2] != null)) {
    			pairs_corr[npair] =  correlateCompositeFD(
    					clt_data_tile[ncam1], // double [][][] clt_data1,
    					clt_data_tile[ncam2], // double [][][] clt_data2,
    		    		lpf,                  // double []     lpf,
    		    		scale_value,
    		    		col_weights,          // double []     col_weights,
    		    		fat_zero);            // double        fat_zero)
    		}
    	}
    	return pairs_corr;
    }
    
    
    /**
     * Calculate FD phase correlation between averaged FD data from two quad (or octal/mixed)
     * cameras, each should be pre-shifted the same disparity
     * @param clt_data_tile_main aberration-corrected FD CLT data for one tile of the main quad camera  [sub-camera][color][quadrant][index]
     * @param clt_data_tile_aux aberration-corrected FD CLT data for one tile of the auxiliary quad camera  [sub-camera][color][quadrant][index]
     * @param lpf optional low-pass filter
     * @param scale_value   scale correlation results to compensate for lpf changes and other factors
     * @param col_weights RBG color weights
     * @param fat_zero fat zero for phase correlations
     * @return 2-d correlation array in line scan order
     */
    public double []  correlateInterCamerasFD( // not used in lwir
    		double [][][][]     clt_data_tile_main,
    		double [][][][]     clt_data_tile_aux,
    		double []           lpf,
    		double              scale_value, // scale correlation value
    		double []           col_weights,
    		double              fat_zero) {
    	if ((clt_data_tile_main == null) || (clt_data_tile_aux == null)) return null;
    	double [][][] clt_mix_main = cltMixCameras(clt_data_tile_main);
    	double [][][] clt_mix_aux =  cltMixCameras(clt_data_tile_aux);
    	double [] inter_cam_corr = correlateCompositeFD(
    			clt_mix_main,         // double [][][] clt_data1,
    			clt_mix_aux,          // double [][][] clt_data2,
	    		lpf,                  // double []     lpf,
	    		scale_value,          // scale correlation value
	    		col_weights,          // double []     col_weights,
	    		fat_zero);            // double        fat_zero)
    	return inter_cam_corr;
    }



    /**
     * Average FD data from 4 sub-cameras (rendered for the same specific disparity/distance),
     * each color component separately. Used to correlate a pair of quad-camera composite images
     * @param clt_data_tile aberration-corrected FD CLT data for one tile [camera][color][quadrant][index]
     * @return averaged for all cameras FD data [color][quadrant][index]
     */
    public double [][][] cltMixCameras( // not used in lwir
    		double [][][][]     clt_data_tile){
    	int tlen = transform_size * transform_size;
    	double [][][] clt_mix = new double [clt_data_tile[0].length][4][tlen];
    	for (int color = 0; color < clt_mix.length; color++) {

    		for (int cltq = 0; cltq <4; cltq++) {
    			for (int i = 0; i < tlen; i++) {
    				for (int cam = 0; cam < clt_data_tile.length; cam++) {
    					if (clt_data_tile[cam][color] != null) {
    						clt_mix[color][cltq][i] += clt_data_tile[cam][color][cltq][i];
    					} else {
    						clt_mix[color] = null; // delete it
    					}
    				}
    			}
    		}
    	}
    	double k = 1.0/clt_data_tile.length;
    	for (int color = 0; color < clt_mix.length; color++) if (clt_mix[color] != null) {
    		for (int cltq = 0; cltq <4; cltq++) {
    			for (int i = 0; i < tlen; i++) {
    				clt_mix[color][cltq][i] *= k;
    			}
    		}
    	}
    	return clt_mix;
    }


   /**
    * Combine (average) several specified correlation pairs that have the same grid (ortho/diagonal, different baselines)
    * Ortho pairs will be transposed as needed to match horizontal pairs, diagonal ones - to match main diagonal (0->3)
    * @param correlations per-pair correlations (according to PAIRS), some may be nulls
    * @param pairs_mask bitmask of the pairs to combine
    * @param diagonal use only pairs that are ortho (false) or diagonal (true)
    * @param baseline_scale use only pairs with this ortho (diagonals have the same scale as ortho) baseline scale
    *        (1 - largest, 2 - half, 4 - quarter)
    * @return single square correlation array, same dimension as the input (now 15x15)
    */
    @Deprecated
    public double [] combineCompatiblePairs(// USED in lwir
    		double [][] correlations,
    		boolean []  pairs_mask,
        	boolean     diagonal,
        	int         baseline_scale
    		) {
    	int width = 2 * transform_size - 1;
    	double [] combo = new double [width * width];
    	int number_combined = 0;
    	// find diagonal/ortho and scale that determine compatible correlations
    	//       			if (pairs_mask[npair] && (correlations[npair]!=null)){
//    	for (int npair = 0; npair < PAIRS.length; npair++) if ((((pairs_mask >> npair) & 1) != 0 ) && (correlations[npair]!=null) &&
    	for (int npair = 0; npair < PAIRS.length; npair++) if (pairs_mask[npair] && (correlations[npair]!=null) &&
    		(isDiagonalPair(npair) == diagonal) && (PAIRS[npair][3] == baseline_scale)){
    		if (isHorizontalPair(npair) || isDiagonalMainPair(npair)) {
    			for (int i = 0; i < combo.length; i++) combo[i]+= correlations[npair][i];
    		} else {
    			int [] transpose_indices = getTransposeAll(isDiagonalOtherPair(npair));
				for (int i = 0; i < transpose_indices.length; i++) {
		  			combo[i]+= correlations[npair][transpose_indices[i]];
				}
    		}
    		number_combined++;
    	}
    	if (number_combined == 0) return null;
    	else if (number_combined > 1) {
    		for (int i = 0; i < combo.length; i++) combo[i] /= number_combined;

    	}
    	return combo;
    }


    /**
     * Get number of compatible pairs among the selection (to be used as weights)
     * @param correlations per-pair correlations (according to PAIRS), some may be nulls
     * @param pairs_mask bitmask of the pairs to combine
     * @param diagonal use only pairs that are ortho (false) or diagonal (true)
     * @param baseline_scale use only pairs with this ortho (diagonals have the same scale as ortho) baseline scale
     *        (1 - largest, 2 - half, 4 - quarter)
     * @return {number of compatible pairs among the selection, index of the base pair}
     */
    @Deprecated
    public int [] getNumberBaseOfCompatiblePairs(// USED in lwir
    		double [][] correlations,
    		boolean []  pairs_mask,
        	boolean     diagonal,
        	int         baseline_scale
    		) {
    	int number_combined = 0;
    	// find diagonal/ortho and scale that determine compatible correlations
    	int base_pair = -1;
    	for (int npair = 0; npair < PAIRS.length; npair++) {
    		if ((isDiagonalPair(npair) == diagonal) && (PAIRS[npair][3] == baseline_scale)){
    			if (base_pair < 0) base_pair = npair;
//    			if ((((pairs_mask >> npair) & 1) != 0 ) && (correlations[npair]!=null)){
       			if (pairs_mask[npair] && (correlations[npair]!=null)){
    				number_combined++;
    			}
    		}
    	}
    	int [] rslt = {number_combined, base_pair};
    	return rslt;
    }




    /**
     * Find minimal/maximal subsampling for selected correlation pairs
     * @param pairs_mask bitmask of selected pairs
     * @return subsampling {min, max}: 1 - no subsampling, also possible 2 and 4 (8-camera)
     */
    public int [] getMinMaxSubSample(int pairs_mask) {// USED in lwir
    	int ss_max = 0;
    	for (int npair = 0; npair < PAIRS.length; npair++) if (((pairs_mask >> npair) & 1) != 0 ) {
    		if (PAIRS[npair][3] > ss_max) {
    			ss_max = PAIRS[npair][3];
    		}
    	}
    	int ss_min = ss_max;
    	for (int npair = 0; npair < PAIRS.length; npair++) if (((pairs_mask >> npair) & 1) != 0 ) {
    		if (PAIRS[npair][3] < ss_min) { // not used in lwir
    			ss_min = PAIRS[npair][3];
    		}
    	}
    	int [] mm = {ss_min,ss_max};
    	return mm;
    }

    /**
     * Process multiple correlation pairs according to pairs_mask (see scaleRotateInterpoateSingleCorrelation() method)
     * This method has limited sub-pixel resolution, it is used to prevent false positives on periodic structures
     * @param correlations array of per-pair correlations (some elements may be nulls)
     * @param pairs_mask bitmask of selected pairs
     * @param hwidth number of the result rows (1 - only main diagonal, 2 - main diagonal and 2 other color ones
     * @return transformed array of correlation arrays [hwidth][2*transform_size-1] (some may be nulls)
     */
    public double [][] scaleRotateInterpoateCorrelations(// USED in lwir
    		double [][] correlations,
    		int         pairs_mask,
    		int         hwidth,
    		int         debug_mask
    		) {
    	double [][] strips = new double [correlations.length][];
    	int [] ss_mm = getMinMaxSubSample(pairs_mask);
    	for (int npair = 0; npair < correlations.length; npair++) if (((pairs_mask & (1 << npair)) != 0) && (correlations[npair] != null)){
    		strips[npair] = scaleRotateInterpoateSingleCorrelation(
    	    		correlations,
    	    		npair,
    	    		ss_mm[0], // sub_sampling,
    	    		hwidth,
    	    		((debug_mask & (1 << npair)) != 0) );
    	}

    	return strips;
    }

    /**
     * Assuming bi-quad camera configuration orthogonal pairs are mapped to a checkerboard cells,
     * while diagonal - to all cells. Correlation is mirrored and averaged around disparity axis
     * orthogonal (hor/vert) pairs are rotated 45 degrees so disparity axis  is along main diagonal
     * that corresponds to first (zero) row of the result.
     * Second (#1) and all odd rows are shifted by 0.5 pix and correspond to other checkerboard color.
     * As this operation is performed only to  locate intersection of the features in all pairs
     * (by offset multiplication) it is only defined for the disparity range present for all selected
     * pairs, data for the lower baseline pairs is truncated.
     * This program uses simple bi-linear interpolation, it is possible to use slower and more precise
     * polynomial interpolation too.
     *
     * @param correlations array of per-pair correlations (some elements may be nulls)
     * @param npair number of this pair to extract
     * @param sub_sampling minimal subsampling for the selected pairs (divide by it)
     * @param hwidth number of the result rows (1 - only main diagonal, 2 - main diagonal and 2 other color ones
     * @return transformed correlation array [hwidth][2*transform_size-1]
     */

    public double [] scaleRotateInterpoateSingleCorrelation(// USED in lwir
    		double [][] correlations,
    		int         npair,
    		int         sub_sampling,
    		int         hwidth,
    		boolean     debug
    		) {
   return scaleRotateInterpoateSingleCorrelation(
        		correlations[npair],
        		hwidth,
            	PAIRS[npair][2], // 0 - hor, 1 - vert, 2 - parallel to row = col (main) diagonal (0->3), 3 -2->1
            	PAIRS[npair][3]/sub_sampling,
        		debug);
    }
    public double [] scaleRotateInterpoateSingleCorrelation(// USED in lwir
    		double []   corr,
    		int         hwidth,
    		int         dir, // 0 - hor, 1 - vert, 2 - parallel to row = col (main) diagonal (0->3), 3 -2->1
    		int         ss,
    		boolean     debug
    		) {
    	int center = transform_size - 1;
    	int width = 2 * center + 1;
    	double [] strip = new double [hwidth * width];
    	int xnum=0,ynum=0;
    	int denom =  ss * ((dir > 1)?1:2);
    	double rdenom = denom;
    	int ilimit = center * denom;
		if (debug) {
			System.out.println("\n============== scaleRotateInterpoateSingleCorrelation() ===============");
		}

    	for (int row = 0; row < hwidth; row++) {
    		for (int scol = -center; scol <= center; scol++) {
    			int indx = row * width + scol + center;
    			for (int down = -1; down < ((row == 0)?0:2); down += 2) {
    				switch (dir) {
    				case 0:
    					xnum = 2 * scol + (row & 1);
    					ynum = down * row;
    					break;
    				case 1:
    					xnum = down * row;
    					ynum = 2 * scol + (row & 1);
    					break;
    				case 2:
    					xnum =  scol - (( down * row    ) >> 1);
    					ynum =  scol + (( down * row + 1) >> 1);
    					break;
    				case 3:
    					xnum =  scol - (( down * row    ) >> 1);
    					ynum = -scol - (( down * row + 1) >> 1);
    					break;
    				}
    				// see if it is within limits
    				if (debug) {
    					System.out.print(String.format("%2d/%2d/%2d: %3d|%3d|%1d   ", row,scol,down,xnum,ynum,denom));
    				}
    				if ((xnum > ilimit) || (xnum < -ilimit) || (ynum > ilimit) || (ynum < -ilimit)) {
    					strip[indx] = Double.NaN;
    				} else {
    					int ix0 = (xnum + ilimit) / denom;
    					int iy0 = (ynum + ilimit) / denom;

    					boolean int_x = (xnum % denom) == 0;
    					boolean int_y = (ynum % denom) == 0;
    					double d;
    					if (int_x && int_y) {
    						d = corr[iy0*width + ix0];
    					} else { //interpolate x only
    						double dx = (xnum + ilimit) / rdenom - ix0;
    						double dy = (ynum + ilimit) / rdenom - iy0;
    						if (int_y) { // not used in lwir
    							d = (1.0-dx)*corr[iy0*width + ix0] + dx*corr[iy0 * width + ix0 + 1];
    						} else if (int_x) { // not used in lwir
    							d = (1.0-dy)*corr[iy0*width + ix0] + dy*corr[iy0 * width + ix0 + width];
    						} else { // bilinear - // USED in LWIR
    							d = (   (1.0 - dx) * (1.0 - dy) * corr[iy0 * width + ix0]) +
    									(       dx  * (1.0 - dy) * corr[iy0 * width + ix0 + 1]) +
    									((1.0 - dx) *        dy  * corr[iy0 * width + ix0 + width]) +
    									(       dx  *        dy  * corr[iy0 * width + ix0 + width + 1]);
    						}
    					}
    					if (row == 0) {
    						strip[indx] = d;
    	    				if (debug) {
    	    					System.out.print(String.format("[%3d:%7.4f] ", indx, d));
    	    				}

    					} else {
    						strip[indx] += 0.5 * d; // average between symmetrical around disparity
    	    				if (debug) {
    	    					System.out.print(String.format("{%3d:%7.4f} ", indx, d));
    	    				}
    					}
    				}
    			}
    		}
			if (debug) {
				System.out.println();
			}
    	}
		if (debug) {
			System.out.println("\n--- Strip ---");
	    	for (int row = 0; row < hwidth; row++) {
	    		if ((row & 1) != 0) {
	    			System.out.print("    ");
	    		}
	    		for (int col = 0; col < width; col++) {
	    			int indx = row * width + col;
	    			System.out.print(String.format("%7.4f ", strip[indx]));
	    		}
    			System.out.println();

	    	}
		}
    	return strip;
// TODO: if there are no diagonals - why interpolate?
    }

    /**
     * Combine calculated (rotated, scaled, interpolated) correlation half-strips
     * @param strips array of per-correlation
     * @param pairs_mask which pairs to combine
     * @param offset add before multiplication, subtract in the end. If negative - use averaging
     *  instead of the shifted multiplication
     * @param twice_diagonal diagonal pairs provide twice denser samples, when true it doubles
     *  the weight of diagonal pairs
     * @return combined correlation data
     */
    public double [] combineInterpolatedCorrelations(  // USED in lwir
    		double [][] strips,
    		int         pairs_mask,
    		double      offset,
    		boolean     twice_diagonal){
    	double [] combo = null;
    	int ncombined = 0;
    	if (offset >= 0) { // use shifted multiplication
    		for (int npair = 0; npair < strips.length; npair++) if (((pairs_mask & (1 << npair)) != 0) && (strips[npair] != null)){
    			if (combo == null) {
    				combo = new double [strips[npair].length];
    				for (int i = 0; i < combo.length; i++) combo[i] = 1.0;
    			}
    			if (twice_diagonal) { // USED in lwir
    				for (int i = 0; i < combo.length; i++) {
    					double d = strips[npair][i] + offset;
    					if (d < 0.0) d = 0.0;
    					combo[i] *= d * d;
    				}
    				ncombined++;
    			} else { // not used in lwir
    				for (int i = 0; i < combo.length; i++) {
    					double d = strips[npair][i] + offset;
    					if (d < 0.0) d = 0.0;
    					combo[i] *= d;
    				}
    			}
    			ncombined++;
    		}
    		double pwr = 1.0/ncombined;
    		if (combo != null) {
    			for (int i = 0; i < combo.length; i++) {
    				if (combo[i] > 0.0) combo[i] = Math.pow(combo[i], pwr) - offset;
    			}
    		}
    	} else { // use addition // not used in lwir
    		for (int npair = 0; npair < strips.length; npair++) if (((pairs_mask & (1 << npair)) != 0) && (strips[npair] != null)){
    			if (combo == null) {
    				combo = new double [strips[npair].length];
    			}
    			if (twice_diagonal) {
    				for (int i = 0; i < combo.length; i++) {
    					combo[i] += 2 * strips[npair][i];
    				}
    				ncombined++;
    			} else {
    				for (int i = 0; i < combo.length; i++) {
    					combo[i] += strips[npair][i];
    				}
    			}
    			ncombined++;
    		}
    		double scale = 1.0/ncombined;
    		for (int i = 0; i < combo.length; i++) {
    			combo[i] *= scale;
    		}
    	}
    	return combo;
    }
//isDiagonalPair
/**
 * Find maximum correlation on the grid (before 08/28/2021)
 * @param data correlation data - full square or just a combined strip with axis at row 0
 * @param axis_only look for the maximum on the disparity axis only. For the rectangular strips and symmetrical forced to be true.
 * @param minMax minimal value of the maximum to be considered valid
 * @param debug print debug data
 * @return a pair of {x,y} or null. x, y are 0 in the center, disparity is -x
 */
    @Deprecated
	public int [] getMaxXYInt( // find integer pair or null if below threshold // USED in lwir
			double [] data,      // [data_size * data_size]
			boolean   axis_only,
			double    minMax,    // minimal value to consider (at integer location, not interpolated)
			boolean   debug)
	{
		int data_width = 2 * transform_size - 1;
		int data_height = data.length / data_width;
		int center_row = 0;
		int center = transform_size - 1;
		if (data_height == data_width) {
			center_row = center; // not used in lwir
		} else {
			axis_only = true;
		}
		return getMaxXYInt( // find integer pair or null if below threshold // USED in lwir
				data,      // [data_size * data_size]
				null,
				data_width,
				center_row,
				axis_only,
				minMax,    // minimal value to consider (at integer location, not interpolated)
				debug);
		
	}
    
/**
 * Find maximum correlation on the grid (before 08/28/2021)
 * @param data correlation data - full square or just a combined strip with axis at row 0
 * @param data_width width of the data array
 * @param center_row for axis-only row0 is center, for square arrays - center row
 * @param axis_only look for the maximum on the disparity axis only. For the rectangular strips and symmetrical forced to be true.
 * @param minMax minimal value of the maximum to be considered valid
 * @param debug print debug data
 * @return a pair of {x,y} or null. x, y are 0 in the center, disparity is -x
 */
	
	public int [] getMaxXYInt( // find integer pair or null if below threshold // USED in lwir
			double [] data,      // [data_size * data_size]
			double [] disp_str,  // if not null, will return {disparity, strength}
			int       data_width,
			int       center_row,
			boolean   axis_only,
			double    minMax,    // minimal value to consider (at integer location, not interpolated)
			boolean   debug)
	{
		//mcorr_comb_disp
//		int data_width = 2 * transform_size - 1;
//		int data_height = data.length / data_width;
		int center = data_width / 2; // transform_size - 1;
		/*
		int center_row = 0;
		if (data_height == data_width) {
			center_row = center; // not used in lwir
		} else {
			axis_only = true;
		}
		*/
		int    imx = 0;
		if (debug){
			System.out.println("getMaxXYInt(): axis_only="+axis_only+", minMax="+minMax);
		}
		if (axis_only) {
			int sol = data_width * center_row;
			imx = sol;
			for (int i = 1; i < data_width; i++) {
				if (Double.isNaN(data[imx]) || (data[sol+i] > data[imx])) imx = sol+i;
			}
		} else { // only for the the square tile // not used in lwir
			for (int i = 0; i < data.length; i++) {
				if (Double.isNaN(data[imx]) || (data[i] > data[imx])) imx = i;
			}
		}
		if (!(data[imx] >= minMax)) {
			if (debug){
				System.out.println("getMaxXYInt() -> null (data["+imx+"] = "+data[imx]+" < "+minMax);
			}
			return null;
		}
		int [] rslt = {imx %  data_width - center, axis_only ? 0 : (imx / data_width - center)};
		if (debug){
			System.out.println("getMaxXYInt() -> "+rslt[0]+"/"+rslt[1]);
		}
		if (disp_str != null) {
			disp_str[0] = -rslt[0] * mcorr_comb_disp;
			disp_str[1] = data[imx];
		}
		return rslt;
	}
	
	
	
	/**
	 * Get fractional center as a "center of mass" inside circle/square from the integer max. Works on the square 2d phase
	 * correlation results to provide data for channel x/y offset
	 * @param data correlation data [(2 * transform_size - 1) * (2 * transform_size - 1)]
	 * @param icenter integer coordinates of the maximal point, relative to the correlation tile center
	 * @param radius positive - within that distance, negative - within 2*(-radius)+1 square
	 * @param debug
	 * @return a pair of {x,y} offsets of the center of masses to the tile center.
	 */
	public double [] getMaxXYCm( // not used in lwir
			double [] data,
			int []    icenter,
			double    radius,    // positive - within that distance, negative - within 2*(-radius)+1 square
			boolean   debug)
	{
		if (icenter == null) {
			double [] rslt = {Double.NaN,Double.NaN};
			return rslt; //gigo
		}
		int center = transform_size - 1;
		int data_size = 2 * transform_size - 1;
		//calculate as "center of mass"
		int iradius = (int) Math.abs(radius);
		int ir2 = (int) (radius*radius);
		boolean square = radius <0;
		double s0 = 0, sx=0,sy = 0;
		for (int y = - iradius ; y <= iradius; y++){
			int dataY = icenter[1] +y;
			int iy = dataY + center;
			if ((iy >= 0) && (iy < data_size)){
				int y2 = y*y;
				for (int x = - iradius ; x <= iradius; x++){
					int dataX = icenter[0] +x;
					int ix = dataX + center;
					double r2 = y2 + x * x;
					if ((ix >= 0) && (ix < data_size) && (square || (r2 <= ir2))){
						double d =  data[iy * data_size + ix];
						s0 += d;
						sx += d * dataX;
						sy += d * dataY;
					}
				}
			}
		}
		double [] rslt = {sx / s0, sy / s0};
		if (debug){
			System.out.println("getMaxXYCm() -> "+rslt[0]+"/"+rslt[1]);
		}
		return rslt;
	}
	
	/**
	 * Get fractional center as a "center of mass" for all pixels
	 * Data should be masked by the caller 
	 * @param data correlation data [(2 * transform_size - 1) * (2 * transform_size - 1)]
	 * @param data_width
	 * @param center_row
	 * @param debug
	 * @return argmax() relative to the tile center
	 */
	
	public double [] getMaxXYCm( // not used in lwir
			double [] data,
			int       data_width,      //  = 2 * transform_size - 1;
			int       center_row,
			boolean   debug)
	{
		int data_height = data.length/data_width;
		int center_x = (data_width - 1)/2; //  = transform_size - 1;
		//calculate as "center of mass"
		double s0 = 0, sx=0,sy = 0;
		for (int iy = 0; iy < data_height; iy++) {
			double y = iy - center_row;
			for (int ix = 0; ix < data_width; ix++) {
				double x = ix - center_x;
				double d =  data[iy * data_width + ix];
				s0 += d;
				sx += d * x;
				sy += d * y;
			}
		}
		double [] rslt = {sx / s0, sy / s0};
		if (debug){
			System.out.println("getMaxXYCm() -> "+rslt[0]+":"+rslt[1]);
		}
		return rslt;
	}
	
	/**
	 * Analyze 1d correlation (single centerline of the 3D phase correlation combined output
	 * from all pairs or only diameters to detect double maximum (simultaneous FG+BG)
	 * @param disparity_scale multiply pixel coordinates in combo_corrs to match bominal disparity
	 *        pixels (measured for quad square camera). It is now 1/sqrt(2) 
	 * @param combo_corrs 2D phase correlation, now 15x15 = 255 pixels long 
	 * @param min_fraction minimal ratio of weaker maximum to the strongest
	 * @return array of 1 or 2 {disparity, strength} pairs (zero pairs if no local max)
	 */
	public double [][] getDoublePoly(
			double    disparity_scale,
			double [] combo_corrs,
			double    min_fraction
			){
		int center_x = 2* (transform_size - 1) * transform_size; // 112
		int [] imx = new int[2];
		for (int i = center_x - transform_size + 2; i < center_x + transform_size - 1; i++) {
			double c = combo_corrs[i];
			if ((c > combo_corrs[i - 1]) && (c > combo_corrs[i + 1])) {
				if ((imx[0] == 0) || (c > combo_corrs[imx[0]]))  {
					imx[1] = imx[0];
					imx[0] = i;
				} else if ((imx[1] == 0) || (c > combo_corrs[imx[1]]))  {
					imx[1] = i;
				}
				i++; // skip next after max
			}
		}
		if (imx[0] == 0) return new double[0][];
		int nm = 1;
		if ((imx[1] > 0) && (combo_corrs[imx[1]]/combo_corrs[imx[0]] > min_fraction)) {
			nm++; 
		}
		double [][] maxes = new double [nm][2];
		
		for (int i = 0; i < nm; i++) {
			double c = combo_corrs[imx[i]];
			double a = (combo_corrs[imx[i] + 1] + combo_corrs[imx[i] - 1]) / 2 - c;
			double b = (combo_corrs[imx[i] + 1] - combo_corrs[imx[i] - 1]);
			maxes[i][0] = -disparity_scale * (imx[i]- center_x - 0.5 * b / a); // disparity, not x!  
			maxes[i][1] = c - 0.25 * b * b / a;
		}
		return maxes;
	}
	
//corr_size	transform_size
	/**
	 * Calculate 1-d maximum location, strength and half-width for the special strip (odd rows shifted by 0.5
	 * Negative values are ignored!
	 * Both x and y half-windows can be variable length (to reduce calculations with 0.0 elements), normalized
	 * so sums of zero element and twice all others are 1.0
	 * Window in Y direction corresponds to correlation strip rows, corresponding to sqrt(2)/2 sensor pixels
	 * for the largest baseline pairs,
	 * Window in X direction has the same sqrt(2)/2 step, but it is half of the horizontal steps of the correlation
	 * results strip
	 * @param data  special diagonal checkerboard array (step in y is 0.5 step in x, odd rows are shifted by 0.5)
	 * of (2 * transform_size - 1)*numrows length, row0 corresponds to the centerline (disparity axis)
	 * @param ixcenter integer argmax on x-axis, relative to the center
	 * @param window_y
	 * @param window_x
	 * @param debug not yet used
	 * @return {argmax from center, weight, half_width} or null
	 */
	public double [] getMaxXCm( // get fractional center as a "center of mass" inside circle/square from the integer max // USED in lwir
			double [] data,      // [data_size * data_size]
			int       data_width,      //  = 2 * transform_size - 1;
			int       center_row,
			int       ixcenter,  // integer center x
			boolean   debug) {
		return getMaxXCm( // get fractional center as a "center of mass" inside circle/square from the integer max
				data,      // double [] data,      // [data_size * data_size]
				data_width,           // int       data_width,      //  = 2 * transform_size - 1;
				center_row,           // int       center_row,
				ixcenter,  // int       ixcenter,  // integer center x
				this.corr_wndy, // double [] window_y,  // (half) window function in y-direction(perpendicular to disparity: for row0  ==1
				this.corr_wndx, // double [] window_x,  // half of a window function in x (disparity) direction
				debug);// boolean   debug);
	}
	
	public double [] getMaxXCmNotch( // get fractional center as a "center of mass" inside circle/square from the integer max // not used in lwir
			double [] data,      // [data_size * data_size]
			int       data_width,      //  = 2 * transform_size - 1;
			int       center_row,
			int       ixcenter,  // integer center x
			boolean   debug) {
		return getMaxXCm(             // get fractional center as a "center of mass" inside circle/square from the integer max
				data,                 // double [] data,      // [data_size * data_size]
				data_width,           // int       data_width,      //  = 2 * transform_size - 1;
				center_row,           // int       center_row,
				ixcenter,             // int       ixcenter,  // integer center x
				this.corr_wndy_notch, // double [] window_y,  // (half) window function in y-direction(perpendicular to disparity: for row0  ==1
				this.corr_wndx,       // double [] window_x,  // half of a window function in x (disparity) direction
				debug);               // boolean   debug);
	}
	
	@Deprecated
	public double [] getMaxXCm( // get fractional center as a "center of mass" inside circle/square from the integer max // not used in lwir
			double [] data,      // [data_size * data_size]
			int       ixcenter,  // integer center x
			boolean   debug) {
		return getMaxXCm( // get fractional center as a "center of mass" inside circle/square from the integer max // not used in lwir
				data,      // [data_size * data_size]
				2 * transform_size - 1, // int       data_width, 
				transform_size - 1,     // int       center_row,
				ixcenter,  // integer center x
				debug);
	}
	
	@Deprecated
	public double [] getMaxXCmNotch( // get fractional center as a "center of mass" inside circle/square from the integer max // not used in lwir
			double [] data,      // [data_size * data_size]
			int       ixcenter,  // integer center x
			boolean   debug) {
		return getMaxXCmNotch( // get fractional center as a "center of mass" inside circle/square from the integer max // not used in lwir
				data,      // [data_size * data_size]
				2 * transform_size - 1, // int       data_width, 
				transform_size - 1,     // int       center_row,
				ixcenter,  // integer center x
				debug);
	}
	
	// No shift by 0.5 for 2021
	public double [] getMaxXCm( // get fractional center as a "center of mass" inside circle/square from the integer max // USED in lwir
			double [] data,      // rectangular strip of 1/2 of the correlation are with odd rows shifted by 1/2 pixels
			int       data_width,      //  = 2 * transform_size - 1;
			int       center_row,
			int       ixcenter,  // integer center x
			double [] window_y,  // (half) window function in y-direction(perpendicular to disparity: for row0  ==1
			double [] window_x,  // half of a window function in x (disparity) direction
			boolean   debug) {
		int data_height = data.length/data_width;
		int center_x = (data_width - 1)/2; //  = transform_size - 1;
		int x0 = center_x + ixcenter; // index of the argmax, starting with 0

//		int data_width = 2 * transform_size - 1;
//		int data_height = data.length/data_width;
//		double wy_scale = 1.0;
		/*
		if (data_height > window_y.length) {
			data_height = window_y.length;
		} else if (data_height < window_y.length) { // re-  // not used in lwir
			double swy = window_y[0];
			for (int i = 1; i < data_height; i++) swy += window_y[i];
			wy_scale = 1.0/swy;
		}
	    */
		double w_scale = 1.0;
		if ((center_row + window_y.length > data_height) || (center_row - window_y.length < 0)) {
			double sw = 0.0;
			for (int i = 0; i < data_height; i++) {
				int dy = i - center_row;
				int ady = (dy > 0) ? dy : -dy;
				if (ady < window_y.length) {
					sw += window_y[ady];
				}
			}
			w_scale /= sw; 
		}

		if ((x0 + window_x.length > data_width) || (x0 - window_x.length < 0)) {
			double sw = 0.0;
			for (int i = 0; i < data_width; i++) {
				int dx = i - x0;
				int adx = (dx > 0) ? dx : -dx;
				if (adx < window_x.length) {
					sw += window_x[adx];
				}
			}
			w_scale /= sw; 
		}
		
		
		//		double [][]dbg_data = null;
		double s0=0.0, sx=0.0, sx2 = 0.0;
		for (int iy = 0; iy < data_height; iy++) {
			int dy = iy - center_row;
			int ady = (dy > 0) ? dy : -dy;
			if (ady < window_y.length) {
				double wy = w_scale * window_y[ady];
				int indx0 = data_width * iy;
				for (int ix = 0; ix < data_width; ix++) {
					int dx = ix - x0; // 0 at argmax
					int adx = (dx > 0) ? dx : -dx;
					if (adx < window_x.length) {
						double d = data[indx0 + ix];
						if (debug)	System.out.print(String.format(" %2d:%2d:%8.5f ", dy, dx, d));
						if (!Double.isNaN(d) && (d > 0.0)) { // with negative d s0 can get very low value (or even negative)
							double w = wy * window_x[adx];
							d*=    w;
							s0+=   d;
							sx +=  d * dx;
							sx2 += d * dx * dx;
						}						
					}
				}
			}
			if (debug)	System.out.println();
		}
		if (debug){
			System.out.println("getMaxXCm() -> s0="+s0+", sx="+sx+", sx2="+sx2+", ixcenter="+ixcenter);
		}

		if (s0 == 0.0) return null;

		double [] rslt = {
				(ixcenter + sx/s0)* mcorr_comb_disp, // /2,              // new center in disparity units, relative to the correlation center
				s0,                                  // total "weight"
				(Math.sqrt(s0*sx2 - sx*sx)/s0)* mcorr_comb_disp}; // /2}; // standard deviation in disparity units (divide weight by the standard deviation for quality?)
		if (debug){
			System.out.println("getMaxXCm() -> "+rslt[0]+"/"+rslt[1]+"/"+rslt[2]);
		}
		return rslt;
	}

	/*
	public double [] getMaxXCm( // get fractional center as a "center of mass" inside circle/square from the integer max // USED in lwir
			double [] data,      // rectangular strip of 1/2 of the correlation are with odd rows shifted by 1/2 pixels
			int       center, //  = transform_size - 1;
			int       data_width, //  = 2 * transform_size - 1;
			int       ixcenter,  // integer center x
			double [] window_y,  // (half) window function in y-direction(perpendicular to disparity: for row0  ==1
			double [] window_x,  // half of a window function in x (disparity) direction
			boolean   debug) {
//		int center = transform_size - 1;
//		int data_width = 2 * transform_size - 1;
		int data_height = data.length/data_width;
		double wy_scale = 1.0;
		if (data_height > window_y.length) {
			data_height = window_y.length;
		} else if (data_height < window_y.length) { // re-  // not used in lwir
			double swy = window_y[0];
			for (int i = 1; i < data_height; i++) swy += window_y[i];
			wy_scale = 1.0/swy;
		}
		double [][]dbg_data = null;
		if (debug) {
			String [] dbg_titles = {"strip","*wnd_y"};
			dbg_data = new double [2][];
			dbg_data[0] =  debugStrip3(data);
			double [] data_0 = data.clone();
			for (int i = 0; i < data_height; i++) {
				for (int j = 0; j <  data_width; j++) {
					data_0[i * data_width + j] *= (i < window_y.length) ? (wy_scale * window_y[i]): 0.0;
				}
			}
			dbg_data[1] =  debugStrip3(data_0);
			int long_width = 2 * data_width; // (2 * transform_size-1);
			if (dbg_data[0] != null) {
				(new ShowDoubleFloatArrays()).showArrays(
						dbg_data,
						long_width,
						dbg_data[0].length/long_width,
						true,
						"Strip",
						dbg_titles);
			}

			System.out.println("getMaxXCm(), ixcenter = "+ixcenter);
			for (int dy = 0; dy < data_height; dy++) {
				if ((dy & 1) != 0) System.out.print("    ");
				for (int dx = 0; dx < data_width; dx++) {
					System.out.print(String.format(" %8.5f", data[dy * data_width + dx]));
				}
				System.out.println();
			}
			System.out.println();
		}
		double s0=0.0, sx=0.0, sx2 = 0.0;
		int x0 = center + ixcenter; // index of the argmax, starting with 0
		for (int dy = 0; dy < data_height; dy++) {
			int odd = dy & 1;
			double wy = ((dy == 0)? wy_scale: (2.0 * wy_scale))*window_y[dy];
			int indx0 = data_width * dy;
			for (int adx = odd; adx < window_x.length; adx+=2) { // index in window_x
				for (int dir = (adx == 0)?1:-1; dir <= 1; dir+=2) {
					// calculate data index
					int idx = (adx * dir) >> 1;
					int x = 2 * idx + odd;
					int x1 = x0 + idx; // correct
					if (debug)	System.out.print(String.format(" %2d:%2d:%d %3d", dy,adx,dir,x));
					if ((x1 >= 0 ) && (x1 < data_width)) {
						double d = data[indx0+x1];
///						if (!Double.isNaN(d)) {
						if (!Double.isNaN(d) && (d > 0.0)) { // with negative d s0 can get very low value (or even negative)
							d*= wy*window_x[adx];
							s0+= d;
							sx +=  d * x; // result x is twice larger (corresponds to window_x)
							sx2 += d * x * x;
							if (debug)	System.out.print(String.format("%8.5f", data[indx0+x1])); //d));
						}
					} else {
						if (debug)	System.out.print("********");
					}
				}
			}
			if (debug)	System.out.println();

		}
		if (debug){
			System.out.println("getMaxXCm() -> s0="+s0+", sx="+sx+", sx2="+sx2+", ixcenter="+ixcenter);
		}

		if (s0 == 0.0) return null;

		double [] rslt = {
				ixcenter + sx/s0/2,              // new center in disparity units, relative to the correlation center
				s0,                              // total "weight"
				Math.sqrt(s0*sx2 - sx*sx)/s0/2}; // standard deviation in disparity units (divide weight by the standard deviation for quality?)
		if (debug){
			System.out.println("getMaxXCm() -> "+rslt[0]+"/"+rslt[1]+"/"+rslt[2]);
		}
		return rslt;
	}
	*/
	
	
	
	
	/**
	 * Generate a half-window for correlation center of mass calculation. Window has a flat top,
	 * then half-cosine fade to zero
	 * @param ihwidth number of samples to keep
	 * @param hwidth half-width corresponding to 50% of the top value
	 * @param blur full width of transition from top value to zero
	 * @param normalize: false - no normalization, [0] = 1.0 (blur permitting),
	 *                   true  - sum of zero element and twice each other == 1.0
	 * @param notch: true - make it a notch filter
	 * @param scale  multiply each value
	 * @return half-window array [ihwidth]
	 */

	public double [] halfFlatTopWindow( // USED in lwir
			int     ihwidth,
			double  hwidth,
			double  blur,
			boolean normalize,
			boolean notch,
			double  scale) {
		double [] wnd = new double [ihwidth];
		for (int i = 0; i < ihwidth; i++) {
			if (blur <= 0.0) { // not used in lwir
				wnd[i] = (i < hwidth)? 1.0:0.0;
			} else {
				if (i < hwidth - blur/2)       wnd[i] = 1.0;
				 else if (i > hwidth + blur/2) wnd[i] = 0.0;
				 else                          wnd[i] = 0.5*(1.0 - Math.sin(Math.PI * (i-hwidth)/blur));
			}
			if (notch) {
				wnd[i] = 1.0 - wnd[i];
			}
		}
		if (normalize) {
			double s = 0.0;
			for (int i = 0; i < ihwidth; i++) {
				s += ((i==0)?1.0:2.0) * wnd[i];
			}
			s = 1/s;
			for (int i = 0; i < ihwidth; i++) {
				wnd[i] *= s;
			}
		}
		for (int i = 0; i < ihwidth; i++) {
			wnd[i] *= scale;
		}
		return wnd;
	}


	/**
	 * Extract center 2-d correlation around zero from the full (now 15x15)
	 * @param hwidth half width of the output tile (0 -> 1x1, 1-> 3x3, 2->5x5)
	 * @param full_corr  full pixel-domain correlation (now 15x15=225 long)
	 * @param center_corr - output array [(2*hwidth+1)*(2*hwidth+1)] or null
	 * @return center_corr - center part of the correlation in linescan order
	 */
    public double [] corrCenterValues( // USED in lwir
    		int       hwidth,
    		double [] full_corr,
    		double [] center_corr) {
    	if (full_corr == null) return null;
    	int center = transform_size - 1;
    	int width = 2 * center + 1;
    	int owidth = 2*hwidth+1;
    	if (center_corr == null) center_corr = new double [owidth*owidth];
    	int indx = 0;
    	int findx = (center - hwidth) * (width + 1); // top left corner
    	for (int row = 0; row < owidth; row++) {
        	for (int col = 0; col < owidth; col++) {
        		center_corr[indx++] = full_corr[findx++];
        	}
        	findx += width-owidth;
    	}
    	return center_corr;
    }

    /**
     * Extract center 2-d correlations around zero from the full (now 15x15) correlations
     * for each of the 6 pairs and 2 combined directions (horizontal, vertical)
     * @param hwidth half width of the output tile (0 -> 1x1, 1-> 3x3, 2->5x5)
     * @param offset add before multiplication, subtract in the end. If negative - use averaging
     *  instead of the shifted multiplication
     * @param full_corr  full pixel-domain correlation (now 15x15=225 long)for each of 6 pairs
     * @param center_corr - output array [(2*hwidth+1)*(2*hwidth+1)]. should be [8][]
     */
    public void corrCenterValues( // USED in lwir
    		int         hwidth,
    		double      offset,
    		double [][] full_corr,
    		double [][] center_corr) {
    	// first 6 layers - directly correspond to pairs (top, bottom, left, right, diagonal main, diagonal other)
    	for (int i = 0; i < 6; i++) {
    		center_corr[i] = corrCenterValues(
    				hwidth,
    				full_corr[i],
    				center_corr[i]);
    	}

    	// combine vertical and horizontal pairs
    	int center = transform_size - 1;
    	int width =  2 * center + 1;
    	int owidth = 2 * hwidth + 1;
    	for (int ndir = 0; ndir < 2; ndir++) { // 0- hor, 1- vert
    		if (center_corr[ndir] == null) center_corr[ndir] = new double [owidth*owidth];
    		int indx = 0;
    		int findx = (center - hwidth) * (width + 1); // top left corner
    		for (int row = 0; row < owidth; row++) {
    			for (int col = 0; col < owidth; col++) {
    				double fc0 = full_corr[2 * ndir + 0][findx]; // 0 (top), 2 (left)
    				double fc1 = full_corr[2 * ndir + 1][findx++]; // 1 - bottom, 3 (right)
    				double cc = 0.0;
    				if (offset >= 0.0) {
    					if ((fc0 > -offset) && (fc1 > -offset)) {
    						cc = Math.sqrt((fc0 + offset) * (fc1 + offset)) - offset;
    					} else {
    						cc = -offset; // smallest value
    					}
    				} else { // not used in lwir
    					cc =  0.5*(fc0+fc1);
    				}

    				center_corr[ndir + 6][indx++] = cc; // save to 6-th and 7-th layer
    			}
    			findx += width-owidth;
    		}
    	}
    }

    /**
     * Save 2d correlation data for one layer, one tile into the combined multi-layer ML array, viewable as an image
     * @param tileX horizontal tile index
     * @param tileY vertical tile index
     * @param ml_hwidth half-width of the preserved 2d correlation (0 - single point, 1 -> 3x3, 2 -> 5x5, 7 - all data)
     * @param ml_data multi-layer array, each layer matches an image of ((2 * ml_hwidth + 1) * tilesX) by ((2 * ml_hwidth + 1) * tilesY) in scanline order
     * Each tile corresponds to  (2 * ml_hwidth + 1) * (2 * ml_hwidth + 1) square in the image. Only selected tiles will be updated, so it is good to initialize array
     * with all Double.NaN values
     * @param ml_layer layer to save tile data
     * @param ml_tile (2 * ml_hwidth + 1) * (2 * ml_hwidth + 1) tile data to be saved
     * @param tilesX image width in tiles
     */
    public void saveMlTile( // USED in lwir
    		int         tileX,
    		int         tileY,
    		int         ml_hwidth,
    		double [][] ml_data,
    		int         ml_layer,
    		double []   ml_tile,
    		int         tilesX) {
    	int tile_width = 2 * ml_hwidth + 1;
    	int full_width = tile_width * tilesX;
    	int oindex = tileY *tile_width * full_width  + tileX * tile_width;
    	for (int row = 0; row < tile_width; row++) {
			System.arraycopy(ml_tile, row * tile_width, ml_data[ml_layer], oindex, tile_width);
			oindex += full_width;
    	}
    }
    /**
     * Save a single value to the combined multi-layer ML array, viewable as an image
     * @param tileX horizontal tile index
     * @param tileY vertical tile index
     * @param ml_hwidth half-width of the preserved 2d correlation (0 - single point, 1 -> 3x3, 2 -> 5x5, 7 - all data)
     * @param ml_data multi-layer array, each layer matches an image of ((2 * ml_hwidth + 1) * tilesX) by ((2 * ml_hwidth + 1) * tilesY) in scanline order
     * Each tile corresponds to  (2 * ml_hwidth + 1) * (2 * ml_hwidth + 1) square in the image. Only selected tiles will be updated, so it is good to initialize array
     * with all Double.NaN values
     * @param ml_layer layer to save tile data
     * @param ml_index data index within tile
     * @param ml_value value to set
     * @param tilesX image width in tiles
     */
    public void saveMlTilePixel(// USED in lwir
    		int         tileX,
    		int         tileY,
    		int         ml_hwidth,
    		double [][] ml_data,
    		int         ml_layer,
    		int         ml_index,
    		double      ml_value,
    		int         tilesX) {
    	int tile_width = 2 * ml_hwidth + 1;
    	int full_width = tile_width * tilesX;
    	int oindex = tileY *tile_width * full_width  + tileX * tile_width + (ml_index/tile_width)*full_width + (ml_index%tile_width) ;
    	ml_data[ml_layer][oindex] = ml_value;
    }
    /**
     * Get a single value from the combined multi-layer ML array, viewable as an image
     * @param tileX horizontal tile index
     * @param tileY vertical tile index
     * @param ml_hwidth half-width of the preserved 2d correlation (0 - single point, 1 -> 3x3, 2 -> 5x5, 7 - all data)
     * @param ml_data multi-layer array, each layer matches an image of ((2 * ml_hwidth + 1) * tilesX) by ((2 * ml_hwidth + 1) * tilesY) in scanline order
     * Each tile corresponds to  (2 * ml_hwidth + 1) * (2 * ml_hwidth + 1) square in the image. Only selected tiles will be updated, so it is good to initialize array
     * with all Double.NaN values
     * @param ml_layer layer to save tile data
     * @param ml_index data index within tile
     * @param tilesX image width in tiles
     * @return value indexed by tileX, tileY, ml_layer and ml_index
     */
    public double restoreMlTilePixel( // not used in lwir
    		int         tileX,
    		int         tileY,
    		int         ml_hwidth,
    		double [][] ml_data,
    		int         ml_layer,
    		int         ml_index,
    		int         tilesX) {
    	int tile_width = 2 * ml_hwidth + 1;
    	int full_width = tile_width * tilesX;
    	int oindex = tileY *tile_width * full_width  + tileX * tile_width + (ml_index/tile_width)*full_width + (ml_index%tile_width) ;
    	return ml_data[ml_layer][oindex];
    }

    /**
     * Get an index of the selected tile+index in a ML array layer
     * @param tileX horizontal tile index
     * @param tileY vertical tile index
     * @param ml_hwidth half-width of the preserved 2d correlation (0 - single point, 1 -> 3x3, 2 -> 5x5, 7 - all data)
     * @param ml_index data index within tile
     * @param tilesX image width in tiles
     * @return index of teh selected pixel in thye whole image (specified by  tileX, tileY, and ml_index)
     */

    public  int getMlTilePixelIndex( // not used in lwir
    		int         tileX,
    		int         tileY,
    		int         ml_hwidth,
    		int         ml_index,
    		int         tilesX) {
    	int tile_width = 2 * ml_hwidth + 1;
    	int full_width = tile_width * tilesX;
    	int oindex = tileY *tile_width * full_width  + tileX * tile_width + (ml_index/tile_width)*full_width + (ml_index%tile_width) ;
    	return oindex;
    }



    public double [] debugStrip( // USED in lwir
    		double [] strip) {
    	if (strip == null) return null;
    	int center = transform_size - 1;
    	int width = 2 * center + 1;
    	int height =  strip.length/width;
    	double [] padded_strip = new double [width*width];
    	for (int row = 0; row < width; row++) {
    		if ((row <= center - height) || (row >= center+ height)) { // not used in lwir
    			for (int j = 0; j<width; j++) padded_strip[row*width+j] = Double.NaN;
    		} else {
    			int srow = (row >= center)? (row - center) : (center - row);
    			for (int j = 0; j<width; j++) padded_strip[row*width+j] = strip[srow*width+j];
    		}
    	}

    	return padded_strip;
    }

    // only show center part, but with correct shift
    public double [] debugStrip2( // USED in lwir
    		double [] strip) {
    	if (strip == null) return null;
    	int center = transform_size - 1;
    	int width = 2 * center + 1;
    	int height =  strip.length/width;
    	double [] padded_strip = new double [width*width];
    	for (int row = 0; row < width; row++) {
    		if ((row <= center - height) || (row >= center+ height)) { // not used in lwir
    			for (int j = 0; j<width; j++) padded_strip[row*width+j] = Double.NaN;
    		} else {
    			int srow = (row >= center)? (row - center) : (center - row);
        		int odd = srow & 1;
    			for (int j = 0; j<width; j++) {
    				int j1 = transform_size / 2 + ((j - odd) >> 1);
    				padded_strip[row * width + j] = strip[srow * width + j1];
    			}
    		}
    	}

    	return padded_strip;
    }

    // Full size/resolution.but on a larger rectangle
    public double [] debugStrip3( // not used in lwir
    		double [] strip) {
    	if (strip == null) return null;
    	int center = transform_size - 1;
    	int width = 2 * center + 1;// 2* transform_siza-1
    	int long_width = 2 * width;
    	int height =  strip.length/width;
    	int long_height = 2 * height - 1;
    	double [] padded_strip = new double [long_height * long_width];
    	for (int row = 0; row < long_height; row++) {
			int srow = (row >= height)? (row - height) : (height - row);
       		int odd = srow & 1;
       		for (int j = 0; j<width; j++) {
       			int j1 = transform_size / 2 + ((j - odd) >> 1);
       			if ((j1 < 0) || (j1 >- width)) {
       				padded_strip[row * width + j] = Double.NaN;
       			} else {
       				padded_strip[row * width + j] = strip[srow * width + j1];
       			}
       		}
    	}

    	return padded_strip;
    }

    public double [] mismatchPairsCM( // returns x-xcenter, y, strength (sign same as disparity) // USED in lwir
    		ImageDttParameters  imgdtt_params,
    		double [][]         corrs,
    		int                 pair_mask, // which pairs to process
    		double              xcenter,   // -disparity to compare
			double              radius,    // positive - within that distance, negative - within 2*(-radius)+1 square
    		int                 debug_level,
    		int                 tileX, // just for debug output
    		int                 tileY)
    {
    	boolean debug = debug_level > 0;
    	int width = 2 * transform_size - 1;
    	int center = transform_size - 1;
    	int center_index = (width + 1) * center; //
    	int num_pairs = 0;
    	for (int i = 0; i <  corrs.length; i++)  if ((corrs[i] != null) && (((1 << i) & pair_mask) != 0)) {
    		num_pairs++;
    	}
    	double [] rslt = new double[3 * num_pairs];

    	int np= 0;
    	int [] icenter = new int[2];
    	int ixcenter = (int) Math.round(xcenter);
    	for (int pair = 0; pair <  corrs.length; pair++)  if ((corrs[pair] != null) && (((1 << pair) & pair_mask) != 0)) {
//    		int this_mask = 1 << pair;
			if (       isHorizontalPair(pair)) {
				icenter[0] =  ixcenter;
				icenter[1] = 0;
			} else if (isVerticalPair(pair)) {
				icenter[0] = 0;
				icenter[1] =  ixcenter;
			} else if (isDiagonalMainPair(pair)) {
				icenter[0] =  ixcenter;
				icenter[1] =  ixcenter;
			} else if (isDiagonalOtherPair(pair)) {
				icenter[0] =  ixcenter;
				icenter[1] = -ixcenter;
			} else { // not used in lwir
				System.out.println("************ BUG: illegal pair type for pair1"+pair);
				return null;
			}
			//calculate as "center of mass"
			int iradius = (int) Math.abs(radius);
			int ir2 = (int) (radius*radius);
			boolean square = radius <0;
			double s0 = 0, sx=0, sy = 0;
			for (int y = - iradius ; y <= iradius; y++){
				int dataY = icenter[1] +y;
				if ((dataY >= -center) && (dataY <= center)){
					int y2 = y*y;
					for (int x = - iradius ; x <= iradius; x++){
						int dataX = icenter[0] +x;
						double r2 = y2 + x * x;
						if ((dataX >= -center) && (dataX <= center) && (square || (r2 <= ir2))){
							double d =  corrs[pair][dataY * width + dataX + center_index];
							if (d > 0.0) {
								s0 += d;
								sx += d * dataX;
								sy += d * dataY;
							}
						}
					}
				}
			}
			double xm = sx / s0;
			double ym = sy / s0;
			int ixm = (int) Math.round(xm);
			int iym = (int) Math.round(ym);
			double s = corrs[pair][iym * width + ixm + center_index];
			if (       isHorizontalPair(pair)) {
				rslt[3 * np + 0] =  xcenter - xm;
				rslt[3 * np + 1] = -ym;
			} else if (isVerticalPair(pair)) {
				rslt[3 * np + 0] = -xm;
				rslt[3 * np + 1] =  xcenter - ym;
			} else if (isDiagonalMainPair(pair)) {
				rslt[3 * np + 0] =  xcenter - xm;
				rslt[3 * np + 1] =  xcenter - ym;
			} else if (isDiagonalOtherPair(pair)) {
				rslt[3 * np + 0] =  xcenter - xm;
				rslt[3 * np + 1] = -xcenter - ym;
			} else { // not used in lwir
				System.out.println("************ BUG: illegal pair type for pair "+pair);
				return null;
			}

			rslt[3 * np + 2] = s;
			if (debug){
				System.out.println("getMaxXYInt() -> "+rslt[0]+"/"+rslt[1]);
			}

    		np++;
    	}
    	return rslt;
    }

    // returns array 3*num_pairs long
    // TODO: now works for small offsets. Maybe add re-calculate int argmax for each pair? xcenter is still needed to subtract Add switch? (small/large correction)
    public double [] mismatchPairs( // returns x-xcenter, y, strength (sign same as disparity) // not used in lwir
    		ImageDttParameters  imgdtt_params,
    		double [][]         corrs,
    		int                 pair_mask, // which pairs to process
    		double              xcenter,   // -disparity to compare
    		double              vasw_pwr,  // value as weight to this power,
    		int                 debug_level,
    		int                 tileX, // just for debug output
    		int                 tileY
) {
    	int num_pairs = 0;
    	for (int i = 0; i <  corrs.length; i++)  if ((corrs[i] != null) && (((1 << i) & pair_mask) != 0)) {
    		num_pairs++;
    	}
    	double [] rslt = new double[3 * num_pairs];

    	int np= 0;
    	for (int pair = 0; pair <  corrs.length; pair++)  if ((corrs[pair] != null) && (((1 << pair) & pair_mask) != 0)) {
    		int this_mask = 1 << pair;
    		if (debug_level > -1) {
    			System.out.println(String.format("mismatchPairs(), np = %d pairs mask = 0x%x", np, this_mask));
    		}
    		Correlations2dLMA lma=corrLMA(
    				imgdtt_params, // ImageDttParameters  imgdtt_params,
    				corrs,         // double [][]         corrs,
    				longToArray(this_mask), //this_mask,     // int                 pair_mask, // which pairs to process
    				true,          // boolean             run_poly_instead, // true - run LMA, false - run 2d polynomial approximation
    				xcenter,       // double              xcenter,   // preliminary center x in pixels for largest baseline
    				vasw_pwr,      // double              vasw_pwr,  // value as weight to this power,
    				debug_level, // -1, // int                 debug_level,
    				tileX,         // int                 tileX, // just for debug output
    				tileY);        //int                 tileY
    		if ((lma == null) || (lma.getPoly() == null)) {
    			rslt[3 * np + 0] = Double.NaN;
    			rslt[3 * np + 1] = Double.NaN;
    			rslt[3 * np + 2] = 0.0;
    		} else {
    			double [] poly_xyvwh = lma.getPoly();
    			if (       isHorizontalPair(pair)) {
    				rslt[3 * np + 0] =  xcenter - poly_xyvwh[0];
    				rslt[3 * np + 1] = -poly_xyvwh[1];
    			} else if (isVerticalPair(pair)) {
    				rslt[3 * np + 0] = -poly_xyvwh[1];
    				rslt[3 * np + 1] =  xcenter - poly_xyvwh[0];
    			} else if (isDiagonalMainPair(pair)) {
    				rslt[3 * np + 0] =  xcenter - poly_xyvwh[0] + poly_xyvwh[1]; // x - y
    				rslt[3 * np + 1] =  xcenter - poly_xyvwh[0] - poly_xyvwh[1]; // x + y
    			} else if (isDiagonalOtherPair(pair)) {
    				rslt[3 * np + 0] =  xcenter - poly_xyvwh[0] + poly_xyvwh[1]; // x - y
    				rslt[3 * np + 1] = -xcenter + poly_xyvwh[0] + poly_xyvwh[1]; // x + y
    			} else {
    				System.out.println("************ BUG: illegal pair type for pair "+pair);
    				return null;
    			}
    			rslt[3 * np + 2] = Double.isNaN(poly_xyvwh[2])?0.0: poly_xyvwh[2];
    		}
    		if ((Double.isNaN(rslt[3 * np + 0]) || Double.isNaN(rslt[3 * np + 1])) && (rslt[3 * np + 2] > 0.0)) {
    			System.out.println("************ mismatchPairs() Fix NaN!!!!! **************");
    			System.out.println(String.format("(), np = %d pairs mask = 0x%x, dx=%f, dy=%f, strength=%f", np, this_mask, rslt[3 * np + 0], rslt[3 * np + 1], rslt[3 * np + 2]));
    		}
    		if (debug_level > -1) {
    			System.out.println(String.format("(), np = %d pairs mask = 0x%x, dx=%f, dy=%f, strength=%f", np, this_mask, rslt[3 * np + 0], rslt[3 * np + 1], rslt[3 * np + 2]));
    		}
    		np++;
    	}
    	return rslt;
    }

// run a single correlation poly
    public double [] single2dPoly( // returns x-xcenter, y, strength (sign same as disparity) // not used in lwir
    		ImageDttParameters  imgdtt_params,
    		double []           corr,
    		double              xcenter,   // -disparity to compare. use 0?
    		double              vasw_pwr,  // value as weight to this power,
    		int                 debug_level,
    		int                 tileX, // just for debug output
    		int                 tileY
    		) {
    	double [] rslt = {Double.NaN,Double.NaN, 0.0};
    	Correlations2dLMA lma=corrLMA(
    			imgdtt_params, // ImageDttParameters  imgdtt_params,
    			corr,          // double [][]         corrs,
    			true,          // boolean             run_poly_instead, // true - run LMA, false - run 2d polynomial approximation
    			xcenter,       // double              xcenter,   // preliminary center x in pixels for largest baseline
    			vasw_pwr,      // double              vasw_pwr,  // value as weight to this power,
    			debug_level,   // -1, // int                 debug_level,
    			tileX,         // int                 tileX, // just for debug output
    			tileY);        //int                 tileY
    	if ((lma != null) && (lma.getPoly() != null)) {
    		double [] poly_xyvwh = lma.getPoly();
    		rslt[0] =  xcenter - poly_xyvwh[0];
    		rslt[1] = -poly_xyvwh[1];
    		rslt[2] = Double.isNaN(poly_xyvwh[2])?0.0: poly_xyvwh[2];
    	}
    	return rslt;
    }

    // ignores negative values
    public double [] single2dCM( // returns x-xcenter, y, strength (sign same as disparity) // not used in lwir
    		ImageDttParameters  imgdtt_params,
    		double []           corr,
    		double              xcenter,   // -disparity to compare. use 0?
			double              radius,    // positive - within that distance, negative - within 2*(-radius)+1 square
    		int                 debug_level,
    		int                 tileX, // just for debug output
    		int                 tileY
    		) {
    	int width = 2 * transform_size - 1;
    	int center = transform_size - 1;
    	int center_index = (width + 1) * center; //
    	double [] rslt = {Double.NaN,Double.NaN, 0.0};
    	int ixcenter = (int) Math.round(xcenter);
    	int [] icenter = {ixcenter, 0};

		//calculate as "center of mass"
		int iradius = (int) Math.abs(radius);
		int ir2 = (int) (radius*radius);
		boolean square = radius <0;
		double s0 = 0, sx=0, sy = 0;
		for (int y = - iradius ; y <= iradius; y++){
			int dataY = icenter[1] +y;
			if ((dataY >= -center) && (dataY <= center)){
				int y2 = y*y;
				for (int x = - iradius ; x <= iradius; x++){
					int dataX = icenter[0] +x;
					double r2 = y2 + x * x;
					if ((dataX >= -center) && (dataX <= center) && (square || (r2 <= ir2))){
						double d =  corr[dataY * width + dataX + center_index];
						if (d > 0.0) {
							s0 += d;
							sx += d * dataX;
							sy += d * dataY;
						}
					}
				}
			}
		}
		double xm = sx / s0;
		double ym = sy / s0;
		int ixm = (int) Math.round(xm);
		int iym = (int) Math.round(ym);
		double s = corr[iym * width + ixm + center_index];

		rslt[0] =  xcenter - xm;
		rslt[1] = -ym;
		rslt[2] = s;
    	return rslt;
    }



    public double [][] corr4dirsLMA( // USED in lwir
    		ImageDttParameters  imgdtt_params,
    		double [][]         corrs,
    		int                 pair_mask, // which pairs to process
    		double              xcenter,   // preliminary center x in pixels for largest baseline
    		double              vasw_pwr,  // value as weight to this power,
    		int                 debug_level,
    		int                 tileX, // just for debug output
    		int                 tileY) {
		if (debug_level > 0) {
			System.out.println("corr4dirsLMA()");
		}
    	double [][] rslt = new double[4][];
    	int [] types = {PAIR_HORIZONTAL, PAIR_VERTICAL, PAIR_DIAGONAL_MAIN, PAIR_DIAGONAL_OTHER};
    	for (int dir = 0; dir < types.length; dir++) {
    		int this_mask = 0;
    		for (int i = 0; i < PAIRS.length; i++) if (((pair_mask & (1 << i)) != 0) && (PAIRS[i][2] == types[dir])) {
    			this_mask |= (1 << i);
    		}
    		if (this_mask == 0) { // not used in lwir
    			rslt[dir] = null;
    		} else {
    			Correlations2dLMA lma=corrLMA(
    					imgdtt_params, // ImageDttParameters  imgdtt_params,
    					corrs,         // double [][]         corrs,
    					longToArray(this_mask), // this_mask,     // int                 pair_mask, // which pairs to process
    		    		false,         // boolean             run_poly_instead, // true - run LMA, false - run 2d polynomial approximation
    		    		xcenter,       // double              xcenter,   // preliminary center x in pixels for largest baseline
    		    		vasw_pwr,      // double              vasw_pwr,  // value as weight to this power,
    		    		debug_level-3, // int                 debug_level,
    		    		tileX,         // int                 tileX, // just for debug output
    		    		tileY);        //int                 tileY
    			if (lma == null) {
    				rslt[dir] = null;
    			} else {
    				rslt[dir] = lma.getDisparityStrengthWidth();
    			}
    		}
    		if (debug_level > 0) {
    			if (rslt[dir] != null) {

    			System.out.println(String.format("corr4dirsLMA() dir = %d, pairs mask = 0x%2x, disparity = %8.5f, strength = %8.5f, width_d = %8.5f, width_y = %8.5f, eff_strength= %8.5f",
    					dir,
    					this_mask,
    					rslt[dir][0],
    					rslt[dir][1],
    					rslt[dir][2] + rslt[dir][3],
    					rslt[dir][2]- rslt[dir][3],
    					rslt[dir][1] * (rslt[dir][2]- rslt[dir][3])/(rslt[dir][2]+ rslt[dir][3]) ));
    			} else {
        			System.out.println(String.format("corr4dirsLMA() rslt[%d] is null", dir));
    			}
    		}
    	}
    	return rslt;
    }

 // TODO: if max/min is strong enough, but the other is not (or not ortho) - just use with  overcorrection = 0.0
 // should always return max-min (NaN OK)
/*
    		if (debug_level > 0) {
    			System.out.println(String.format("corr4dirsLMA() dir = %d, pairs mask = 0x%2x, disparity = %8.5f, strength = %8.5f, width_d = %8.5f, width_y = %8.5f, eff_strength= %8.5f",
    					dir,
    					this_mask,
    					rslt[dir][0],
    					rslt[dir][1],
    					rslt[dir][2] + rslt[dir][3],
    					rslt[dir][2]- rslt[dir][3],
    					rslt[dir][1] * (rslt[dir][2]- rslt[dir][3])/(rslt[dir][2]+ rslt[dir][3]) ));
    		}
// each element may be null, data may contain NaN
 */
    public double [] foregroundCorrect( // USED in lwir
    		boolean     bg,
    		boolean     ortho,
    		double [][] dir_disp_strength, //
    		double      full_disp,
    		double      min_strength,
    		double      min_eff,
    		double      min_eff_ratio,
    		double      max_hwidth, //  =          3.0;  // maximal half-width disparity  direction to try to correct far objects
    		double      min_diff,
    		double      overcorrection,
    		double      lim_overcorr,
    		boolean debug) {
    	double mn = Double.NaN, mx = Double.NaN;
    	int imn = -1, imx = -1;
    	double [] eff_strength = new double [4];
    	double [] width_d =      new double [4];

    	boolean [] strong = {false,false,false,false};
    	for (int i = 0; i < dir_disp_strength.length; i++) {
    		if (dir_disp_strength[i] != null){
    			if (Double.isNaN(mn) || (dir_disp_strength[i][0] < mn)) {mn = dir_disp_strength[i][0]; imn = i;}
    			if (Double.isNaN(mx) || (dir_disp_strength[i][0] > mx)) {mx = dir_disp_strength[i][0]; imx = i;}
    			eff_strength[i] = dir_disp_strength[i][1] * (dir_disp_strength[i][2]- dir_disp_strength[i][3])/(dir_disp_strength[i][2]+ dir_disp_strength[i][3]);
    			width_d[i] = dir_disp_strength[i][2] + dir_disp_strength[i][3]; // width in disparity direction
    			strong[i] = (dir_disp_strength[i][1] >= min_strength) && (eff_strength[i] >= min_eff);
    			// TODO what about strong ortho? yes, use - if not strong - overcorr = 0
    		} else {
    			eff_strength[i] = Double.NaN;
    		}
    	}
    	boolean are_ortho = (imn ^ imx) == 1;
		double corr = overcorrection * (mx - mn);

		int isel =   bg ? imn : imx; // isel may be -1 !!!
		int iother = bg ? imx : imn;
		int iortho = isel ^ 1;

    	if (debug) {
    		System.out.println("foregroundCorrect("+bg+", "+ortho+",...): full: "+full_disp+", min = "+mn+", mx = "+mx+", diff="+(mx-mn)+", are_ortho="+are_ortho);
    		System.out.println("foregroundCorrect() min_strength="+min_strength+", min_diff="+min_diff+", overcorrection="+overcorrection+",lim_overcorr="+lim_overcorr);
    		System.out.println("strong = ["+strong[0]+", "+strong[1]+", "+strong[2]+", "+strong[3]+"]");
    		System.out.println(String.format("eff_strength = [%8.5f, %8.5f, %8.5f, %8.5f]", eff_strength[0], eff_strength[1], eff_strength[2], eff_strength[3]));
    		System.out.println(String.format("width_d =      [%8.5f, %8.5f, %8.5f, %8.5f]", width_d[0],      width_d[1],      width_d[2],      width_d[3]));
    	}
		if ((isel <0) || !strong[isel]) {
			corr = Double.NaN;
			if (debug) System.out.println("Direction with "+(bg?"min":"max")+" disparity is not strong enough -> no correction");
		} else 		if (width_d[isel] > max_hwidth) {
			corr = Double.NaN;
			if (debug) System.out.println("Direction with "+(bg?"min":"max")+" has too wide correlation maximum in disparity direction -> no correction ("+
					"width_d["+isel+"] = "+width_d[isel]+" > "+max_hwidth+")");
		} else 		if ((eff_strength[isel]/eff_strength[iortho] < min_eff_ratio)) {
			corr = Double.NaN;
			if (debug) System.out.println("Direction with "+(bg?"min":"max")+" has effective strength ratio too small -> no correction ("+
					"(eff_strength["+isel+"] = "+eff_strength[isel]+")/(eff_strength["+iortho+"] = "+eff_strength[iortho]+") < "+min_eff_ratio+")");
		} else {
			if (ortho && !are_ortho) {
				if (debug) System.out.println("Min/max orthogonal required, but they are not ("+imn+","+imx+") -> overcorrection = 0 (just using "+(bg?"min":"max") );
				corr = 0.0;
			} else { // see if min/max is strong enough
				if (!strong[iother]) {
					corr = 0.0;
					if (debug) System.out.println("Orthogonal direction is not strong enough -> overcorrection = 0 (just using "+(bg?"min":"max") );
				}
			}
			if (!(mx >= (mn + min_diff))) { // so NaN are OK // not used in lwir
				corr = Double.NaN;
	    		if (debug) System.out.println("difference max - min="+(mx - mn)+" < "+min_diff+" -> no fo correction");
			}
		}
    	double disp = full_disp;
    	if (!Double.isNaN(corr)) {
    		double lim;
    		if (bg) { // not used in lwir
    			lim = full_disp - (full_disp - mn) * lim_overcorr;
    			disp =  Math.max(mn - corr, lim);
    		} else {
    			lim = full_disp + (mx - full_disp) * lim_overcorr;
    			disp =  Math.min(mx + corr, lim);
    		}
    		if (debug) System.out.println("lim =  "+lim+", disp = "+disp);
    	}
    	double [] rslt = {disp, (isel >=0) ? eff_strength[isel]:0.0, mx - mn, are_ortho ?1.0 : 0.0};
		if (debug) System.out.println(String.format("foregroundCorrect() -> [%8.5f, %8.5f, %8.5f, %3.1f]", rslt[0], rslt[1], rslt[2], rslt[3]));

    	return    rslt;
    }

    double [][] repackCluster(
    		double [][][] data,
    		int clust_width){
		int  clust_height = data.length/clust_width;
    	int corr_size = 2 * transform_size - 1;
    	int nslices = 0;
    	for (int i = 0; i < data.length; i++) {
    		if (data[i] != null) {
    			nslices = data[i].length;
    			break;
    		}
    	}
    	int out_width =  (corr_size+1)*clust_width;
    	int out_height = (corr_size+1)*clust_height;
    	double [][] mosaic = new double [nslices][];
    	for (int ns = 0; ns <  nslices; ns++) {
    		boolean slice_exists = false;
    		for (int i = 0; i < data.length; i++){
    			if ((data[i] != null) && (data[i][ns] != null)) {
    				slice_exists = true;
    				break;
    			}
    		}
    		if (slice_exists) {
    			mosaic[ns] = new double [out_height*out_width];
    			for (int cy = 0; cy < clust_height; cy++) {
    				for (int cx = 0; cx < clust_width; cx++) {
    					int nclust = cy * clust_width + cx;
    					int indx = cy * (corr_size+1)*out_width + cx * (corr_size+1);
    					if ((data[nclust] == null) || (data[nclust][ns] == null)) {
    						for (int i = 0; i < corr_size; i ++) {
    							for (int j = 0; j < corr_size; j ++) {
    								mosaic[ns][indx + j] = Double.NaN;
    							}
    							indx += out_width;
    						}
    					} else {
    						for (int i = 0; i < corr_size; i ++) {
    							for (int j = 0; j < corr_size; j ++) {
    								mosaic[ns][indx + j] = data[nclust][ns][i * corr_size + j];
    							}
    							indx += out_width;
    						}
    					}
    				}
    			}
    		}
    	}
    	return mosaic;
    }

    public Corr2dLMA corrLMA2Multi( // multi-tile
    		ImageDttParameters  imgdtt_params,
    		int                 clust_width,
    		double [][]         corr_wnd, // correlation window to save on re-calculation of the window
    		double []           corr_wnd_inv_limited, // correlation window, limited not to be smaller than threshold - used for finding max/convex areas (or null)
    		double [][][]       corrs, // per tile, per pair, 2 correlation in line-scan order
    		double [][][]       disp_dist, // per tile, per camera disparity matrix as a 1d (linescan order)
    		double [][]         rXY, // non-distorted X,Y offset per nominal pixel of disparity
    		boolean []          pair_mask, // which pairs to process
    		double[][]          disp_str,   // -preliminary center x in pixels for largest baseline
    		double              vasw_pwr,  // value as weight to this power,
    		int                 debug_level,
    		int                 tileX, // just for debug output
    		int                 tileY)
    {
    	// corrs are organized as PAIRS, some are null if not used
    	// for each enabled and available pair find a maximum, filter convex and create sample list
    	boolean debug_graphic = imgdtt_params.lma_debug_graphic && (debug_level > -1) ;
//    	boolean debug_graphic = true;
//    	boolean debug_second_all = false; // true; // alse; // true;
		int  clust_height = corrs.length/clust_width;
    	int ntiles = corrs.length;
    	DoubleGaussianBlur gb = null;
    	if (imgdtt_params.lma_sigma > 0) gb = new DoubleGaussianBlur();
    	int center =       transform_size - 1;
    	int corr_size = 2 * transform_size - 1;
//    	double [][] blur_max = new double [corrs.length][];
    	if (debug_level > -2) {
    		System.out.println("Debugging corrLMA2()// multi-tile");
    	}
    		Corr2dLMA lma = new Corr2dLMA(
    			corrs.length,
				1, // int         numMax,
				this, // Correlation2d correlation2d,
    			transform_size,
    			corr_wnd,
    			rXY, //double [][] rXY, // non-distorted X,Y offset per nominal pixel of disparity
    			imgdtt_params.lma_gaussian//boolean gaussian_mode
    			);

    	double [][][] dbg_corr =    debug_graphic ? new double [corrs.length][][] : null;
    	double [][][] dbg_weights = debug_graphic ? new double [corrs.length][][] : null;
    	int dbg_out_width =  (corr_size + 1) * clust_width;
    	int dbg_out_height = (corr_size + 1) * clust_height;

    	if (debug_graphic) {
    	    double [][] dbg_corrs = repackCluster(
    	    		corrs,
    	    		clust_width);

    		(new ShowDoubleFloatArrays()).showArrays(
    				dbg_corrs,
    				dbg_out_width,
    				dbg_out_height,
    				true,
    				"corr_pairs-"+"_x"+tileX+"_y"+tileY);

    	}
    	int numpairs = 0;
    	for (int ntile = 0; ntile < ntiles; ntile++) if ((corrs[ntile] != null) &&
    			((disp_str == null) || ((disp_str[ntile] != null) && (disp_str[ntile][1] > 0.0)))){
    		double[][] corr = new double[corrs[ntile].length][];
    		double [][] filtWeight = new double [corrs[ntile].length][];
//    		for (int npair = 0; npair < corrs[ntile].length; npair++) if ((corrs[ntile][npair] != null) && (((pair_mask >> npair) & 1) !=0)){
       		for (int npair = 0; npair < corrs[ntile].length; npair++) if ((corrs[ntile][npair] != null) && (pair_mask[npair])){
    			corr[npair] = corrs[ntile][npair].clone();
    			if (corr_wnd_inv_limited != null) {
    				for (int i = 0; i < corr.length; i++) {
    					corr[npair][i] *= corr_wnd_inv_limited[i];
    				}
    			}
    			if (imgdtt_params.lma_sigma > 0) {
    				gb.blurDouble(corr[npair], corr_size, corr_size, imgdtt_params.lma_sigma, imgdtt_params.lma_sigma, 0.01);
    			}
    			int imx = imgdtt_params.lma_soft_marg * (corr_size + 1);
    			for (int iy = imgdtt_params.lma_soft_marg; iy < (corr_size - imgdtt_params.lma_soft_marg); iy++) {
        			for (int ix = imgdtt_params.lma_soft_marg; ix < (corr_size - imgdtt_params.lma_soft_marg); ix++) {
        				int indx = iy * corr_size + ix;
        				if (corr[npair][indx] > corr[npair][imx]) imx = indx;

        			}
    			}
    			// filter convex
    			int ix0 = (imx % corr_size) - center; // signed, around center to match filterConvex
    			int iy0 = (imx / corr_size) - center; // signed, around center to match filterConvex
    			filtWeight[npair] =  filterConvex(
    					corr[npair],                  // double [] corr_data,
    					imgdtt_params.cnvx_hwnd_size, // int       hwin,
    					ix0,                          // int       x0,
    					iy0,                          // int       y0,
    					imgdtt_params.cnvx_add3x3,    // boolean   add3x3,
    					imgdtt_params.cnvx_weight,    // double    nc_cost,
    					(debug_level > 2));           // boolean   debug);
    			// remove too small clusters, or all collinear
    		    if (!checkCluster(
    		    		filtWeight[npair], // double [] weights,
    		    		corr_size, // int       size,
    		    		imgdtt_params.cnvx_min_samples, // int min_samples,
    		    		imgdtt_params.cnvx_non_coll // boolean non_coll
    		    		)) {
    		    	filtWeight[npair] = null;
    		    } else {
    		    	numpairs++;
    		    	if (dbg_corr    != null) 	{
    		    		if (dbg_corr[ntile] == null) {
    		    			dbg_corr[ntile] = new double [corrs[ntile].length][];
    		    		}
    		    		dbg_corr   [ntile][npair] = corr[npair];
    		    	}
    		    	if (dbg_weights != null) {
    		    		if (dbg_weights[ntile] == null) {
    		    			dbg_weights[ntile] = new double [corrs[ntile].length][];
    		    		}
    		    		dbg_weights[ntile][npair] = filtWeight[npair];
    		    	}
    		    }
    			// Normalize weight for each pair to compensate for different number of convex samples?
    		}

    		// numpairs
    		if (numpairs >= imgdtt_params.cnvx_min_pairs) {
    			for (int npair = 0; npair < corrs[ntile].length; npair++) if (filtWeight[npair] != null){
//    				int fcam = PAIRS[npair][0];
//    				int scam = PAIRS[npair][1];
    				for (int i = 1; i < filtWeight[npair].length; i++) if (filtWeight[npair][i] > 0.0) {
    					int ix = i % corr_size; // >=0
    					int iy = i / corr_size; // >=0
    					double v = corrs[ntile][npair][i]; // not blurred
    					double w = filtWeight[npair][i];
    					if (vasw_pwr != 0) {
    						w *= Math.pow(Math.abs(v), vasw_pwr);
    					}
    					lma.addSample( // x = 0, y=0 - center
    							ntile, // tile
    							npair,
    							0,     // int    imax,   // number of maximum (just a hint)
//    							fcam,  // int    fcam, // first  camera index
//    							scam,  // int    scam, // second camera index
    							ix,    // int    x,      // x coordinate on the common scale (corresponding to the largest baseline), along the disparity axis
    							iy,    // int    y,      // y coordinate (0 - disparity axis)
    							v,     // double v,       // correlation value at that point
    							w);     //double w){      // sample weight
    				}
    			}
    		}
    	}


    	if (debug_graphic && (dbg_corr != null)) {

    	    double [][] dbg_corrs = repackCluster(
    	    		dbg_corr,
    	    		clust_width);

    		(new ShowDoubleFloatArrays()).showArrays( // empty array
    				dbg_corrs,
    				dbg_out_width,
    				dbg_out_height,
    				true,
    				"corr_blurred"+"_x"+tileX+"_y"+tileY);

    	    double [][] dbg_w = repackCluster(
    	    		dbg_weights,
    	    		clust_width);

    		(new ShowDoubleFloatArrays()).showArrays(
    				dbg_w,
    				dbg_out_width,
    				dbg_out_height,
    				true,
    				"corr_weights"+"_x"+tileX+"_y"+tileY);
    	}

    	boolean lmaSuccess = false;
    	while (!lmaSuccess) {
    		boolean OK = lma.initVector( // USED in lwir
    				null, // 			boolean [] adjust_disparities, // null - adjust all, otherwise - per maximum
    				imgdtt_params.lma_adjust_wm,  // boolean adjust_width,     // adjust width of the maximum - lma_adjust_wm
    				imgdtt_params.lma_adjust_ag,  // boolean adjust_scales,    // adjust 2D correlation scales - lma_adjust_ag
    				imgdtt_params.lma_adjust_wy,  // boolean adjust_ellipse,   // allow non-circular correlation maximums lma_adjust_wy
    				imgdtt_params.lma_adjust_wxy, // boolean adjust_lazyeye_par,   // adjust disparity corrections parallel to disparities  lma_adjust_wxy
    				imgdtt_params.lma_adjust_ly1, // boolean adjust_lazyeye_ortho, // adjust disparity corrections orthogonal to disparities lma_adjust_ly1
    				new double[][][] {disp_str}, // xcenter, 
    				imgdtt_params.lma_half_width, // double  half_width,       // A=1/(half_widh)^2   lma_half_width
    				imgdtt_params.lma_cost_wy,     // double  cost_lazyeye_par,     // cost for each of the non-zero disparity corrections        lma_cost_wy
    				imgdtt_params.lma_cost_wxy     // double  cost_lazyeye_odtho    // cost for each of the non-zero ortho disparity corrections  lma_cost_wxy
    				);
    		if (!OK) {
    			return null; // empty cluster
    		}
    		lma.setMatrices(disp_dist);
    		lma.initMatrices(); // should be called after initVector and after setMatrices
    		//    	boolean lmaSuccess = false;
    		if (debug_level > 1) {
    			System.out.println("Input data:");
    			lma.printInputDataFx(false);
    			lma.printParams();

    		}
    		lmaSuccess = 	lma.runLma(
    				imgdtt_params.lma_lambda_initial,      // double lambda,            // 0.1
    				imgdtt_params.lma_lambda_scale_good,   // double lambda_scale_good, // 0.5
    				imgdtt_params.lma_lambda_scale_bad,    // double lambda_scale_bad,  // 8.0
    				imgdtt_params.lma_lambda_max,          // double lambda_max,        // 100
    				imgdtt_params.lma_rms_diff,            // double rms_diff,          // 0.001
    				imgdtt_params.lma_num_iter,            // int    num_iter,          // 20
    				imgdtt_params.lma_debug_level1); //4); // debug_level);             // int    debug_level) // > 3
    		if (!lmaSuccess && (lma.getBadTile() >= 0)) {
        		if (debug_level > -1) { // -20) {
        			System.out.println("Found bad tile/pair (probably wrong initial maximum - try around preliminary? "+lma.getBadTile());
        		}
    		} else {
    			break;
    		}
    	}


    	double [][] ds = null;

    	if (lmaSuccess) {
    		lma.updateFromVector();
    		double [] rms = lma.getRMS();
    		if (debug_level > 0) {
    			System.out.println("LMA -> "+lmaSuccess+" RMS= "+rms[0]+", pure RMS= "+rms[1]);
    			lma.printParams();
    		}

    		//    		double [][] ds =         lma.getDisparityStrength();
    		ds = lma.lmaDisparityStrength(
    				imgdtt_params.lmas_min_amp,      //  minimal ratio of minimal pair correlation amplitude to maximal pair correlation amplitude
    				imgdtt_params.lma_max_rel_rms,  // maximal relative (to average max/min amplitude LMA RMS) // May be up to 0.3)
    				imgdtt_params.lma_min_strength, // minimal composite strength (sqrt(average amp squared over absolute RMS)
    				imgdtt_params.lma_min_ac,       // minimal of A and C coefficients maximum (measures sharpest point/line)
    				imgdtt_params.lma_min_min_ac,   // minimal of A and C coefficients minimum (measures sharpest point)
    				imgdtt_params.lma_max_area,     //double  lma_max_area,     // maximal half-area (if > 0.0)
    				imgdtt_params.lma_str_scale,    // convert lma-generated strength to match previous ones - scale
    				imgdtt_params.lma_str_offset    // convert lma-generated strength to match previous ones - add to result
    				);
    		if (debug_level > 0) lma.printStats(ds,clust_width);


    		if (debug_level > 1) {
    			System.out.println("Input data and approximation:");
    			lma.printInputDataFx(true);
    		}
    		//    	boolean debug_second_all = true;

    		if (lmaSuccess && imgdtt_params.lma_second) {
    			double [] all_pars_save = null;
    			ArrayList<Sample> samples = lma.filterSamples(ds);
    			if (!lma.samplesExist()) {
    				if (debug_level > 1) {
    					System.out.println("No tiles left in a cluster!");
    				}
    				return null;
    			}
    			lmaSuccess = false;
    			while (!lmaSuccess ) {
    				if (imgdtt_params.lma_gaussian == imgdtt_params.lma_second_gaussian) {
    					all_pars_save = lma.all_pars.clone();
    				} else { // have to restart LMA initialization
    					lma = new Corr2dLMA(
    							corrs.length,
    							1, // int         numMax,
    							this, // 				Correlation2d correlation2d,
    							transform_size,
    							corr_wnd,
    							rXY, //double [][] rXY, // non-distorted X,Y offset per nominal pixel of disparity
    							imgdtt_params.lma_second_gaussian//boolean gaussian_mode
    							);
    					lma.setSamples(samples); // restore samples
    				}
    				lma.initVector( // USED in lwir
    	    				null, // 			boolean [] adjust_disparities, // null - adjust all, otherwise - per maximum
    						imgdtt_params.lma_adjust_wm,  // boolean adjust_width,     // adjust width of the maximum - lma_adjust_wm
    						imgdtt_params.lma_adjust_ag,  // boolean adjust_scales,    // adjust 2D correlation scales - lma_adjust_ag
    						imgdtt_params.lma_adjust_wy,  // boolean adjust_ellipse,   // allow non-circular correlation maximums lma_adjust_wy
    						true, // imgdtt_params.lma_adjust_wxy, // boolean adjust_lazyeye_par,   // adjust disparity corrections parallel to disparities  lma_adjust_wxy
    						true, // imgdtt_params.lma_adjust_ly1, // boolean adjust_lazyeye_ortho, // adjust disparity corrections orthogonal to disparities lma_adjust_ly1
//    						ds, // disp_str, // xcenter_str,                  // double [][] disp_str,         // initial value of disparity/strength/?
    	    				new double[][][] {ds}, // xcenter, 
    						imgdtt_params.lma_half_width, // double  half_width,       // A=1/(half_widh)^2   lma_half_width
    						imgdtt_params.lma_cost_wy,     // double  cost_lazyeye_par,     // cost for each of the non-zero disparity corrections        lma_cost_wy
    						imgdtt_params.lma_cost_wxy     // double  cost_lazyeye_odtho    // cost for each of the non-zero ortho disparity corrections  lma_cost_wxy
    						);
    				if (all_pars_save != null) {
    					lma.all_pars = all_pars_save;
    					lma.toVector();
    				}

    				lma.setMatrices(disp_dist);
    				lma.initMatrices(); // should be called after initVector and after setMatrices
    				//boolean
    				if (debug_level > 1) {
    					System.out.println("Input data:");
    					lma.printInputDataFx(false);
    					lma.printParams();

    				}

    				lmaSuccess = 	lma.runLma(
    						imgdtt_params.lma_lambda_initial,     // double lambda,           // 0.1
    						imgdtt_params.lma_lambda_scale_good,  // double lambda_scale_good,// 0.5
    						imgdtt_params.lma_lambda_scale_bad,   // double lambda_scale_bad, // 8.0
    						imgdtt_params.lma_lambda_max,         // double lambda_max,       // 100
    						imgdtt_params.lma_rms_diff,           // double rms_diff,         // 0.001
    						imgdtt_params.lma_num_iter,           // int    num_iter,         // 20
    						debug_level); // 2); //4); // debug_level);       // int    debug_level) // > 3
    	    		if (!lmaSuccess && (lma.getBadTile() >= 0)) {
    					if (debug_level > -2) {
    						System.out.println("Found bad tile/pair in a second pass (probably wrong initial maximum - try around preliminary? "+lma.getBadTile());
    					}
    	    		} else {
    	    			break;
    	    		}
    			}

    			//double []
    			rms = lma.getRMS();
    			if (debug_level > 0) {
    				System.out.println("LMA -> "+lmaSuccess+" RMS="+rms[0]+", pure RMS="+rms[1]);
    				lma.printParams();
    			}

    			if (lmaSuccess) {
        			lma.updateFromVector();
        			if (debug_level > 1) {
        				System.out.println("Input data and approximation:");
        				lma.printInputDataFx(true);
        			}
    				//        		double [][] ds =         lma.getDisparityStrength();
    				ds = lma.lmaDisparityStrength(
    	    				imgdtt_params.lmas_min_amp,      //  minimal ratio of minimal pair correlation amplitude to maximal pair correlation amplitude
    						imgdtt_params.lma_max_rel_rms,  // maximal relative (to average max/min amplitude LMA RMS) // May be up to 0.3)
    						imgdtt_params.lma_min_strength, // minimal composite strength (sqrt(average amp squared over absolute RMS)
    						imgdtt_params.lma_min_ac,       // minimal of A and C coefficients maximum (measures sharpest point/line)
    	    				imgdtt_params.lma_min_min_ac,   // minimal of A and C coefficients minimum (measures sharpest point)
    	    				imgdtt_params.lma_max_area,      //double  lma_max_area,     // maximal half-area (if > 0.0)
    						imgdtt_params.lma_str_scale,    // convert lma-generated strength to match previous ones - scale
    						imgdtt_params.lma_str_offset    // convert lma-generated strength to match previous ones - add to result
    						);
    				if (debug_level > 0) { // -2) {
    					lma.printStats(ds,clust_width);
    				}
    			}
    		}
    	}

    	// ds[tile]{disparity, strength}

    	if (debug_graphic && lmaSuccess) {
    		String [] sliceTitles = lma.dbgGetSliceTitles();
    		if (corrs.length == 1) {
    			(new ShowDoubleFloatArrays()).showArrays(
    					lma.dbgGetSamples(ds, 0)[0],
    					corr_size,
    					corr_size,
    					true,
    					"corr_values"+"_x"+tileX+"_y"+tileY, sliceTitles);
    			(new ShowDoubleFloatArrays()).showArrays(
    					lma.dbgGetSamples(ds, 1)[0],
    					corr_size,
    					corr_size,
    					true,
    					"corr_weights"+"_x"+tileX+"_y"+tileY, sliceTitles);
    			(new ShowDoubleFloatArrays()).showArrays(
    					lma.dbgGetSamples(ds, 2)[0],
    					corr_size,
    					corr_size,
    					true,
    					"corr_fx"+"_x"+tileX+"_y"+tileY, sliceTitles);
    		} else {
    			double [][] repacked_y = repackCluster(lma.dbgGetSamples(ds, 0),clust_width);
    			double [][] repacked_fx = repackCluster(lma.dbgGetSamples(ds, 2),clust_width);
    			double [][] y_minux_fx = new double [repacked_y.length][];
    			for (int i = 0; i < repacked_y.length; i++) if ((repacked_y[i] != null) && (repacked_fx[i] != null)){
    				y_minux_fx[i] = new double [repacked_y[i].length];
    				for (int j = 0; j < y_minux_fx[i].length; j++) y_minux_fx[i][j] = repacked_y[i][j] - repacked_fx[i][j];
    			}
    			(new ShowDoubleFloatArrays()).showArrays(
    					repacked_y,
    					dbg_out_width,
    					dbg_out_height,
    					true,
    					"y"+"_x"+tileX+"_y"+tileY, sliceTitles);
    			(new ShowDoubleFloatArrays()).showArrays(
    					repacked_fx,
    					dbg_out_width,
    					dbg_out_height,
    					true,
    					"fx"+"_x"+tileX+"_y"+tileY, sliceTitles);
    			(new ShowDoubleFloatArrays()).showArrays(
    					y_minux_fx,
    					dbg_out_width,
    					dbg_out_height,
    					true,
    					"y-fx"+"_x"+tileX+"_y"+tileY, sliceTitles);
    			(new ShowDoubleFloatArrays()).showArrays(
    					repackCluster(lma.dbgGetSamples(ds, 1),clust_width),
    					dbg_out_width,
    					dbg_out_height,
    					true,
    					"weights_late"+"_x"+tileX+"_y"+tileY, sliceTitles);
    		}
    	}

    	return lmaSuccess? lma: null;
    }

    public Corr2dLMA corrLMA2Single( // single tile
    		ImageDttParameters  imgdtt_params,
    		boolean             adjust_ly, // adjust Lazy Eye
    		double [][]         corr_wnd, // correlation window to save on re-calculation of the window
    		double []           corr_wnd_inv_limited, // correlation window, limited not to be smaller than threshold - used for finding max/convex areas (or null)
    		double [][]         corrs, // may have more elements than pair_mask (corrs may have combo as last elements)
    		double [][]         disp_dist, // per camera disparity matrix as a 1d (linescan order)
    		double [][]         rXY, // non-distorted X,Y offset per nominal pixel of disparity
    		boolean []          pair_mask, // which pairs to process
    		double[]            disp_str,   // -preliminary center x in pixels for largest baseline
    		double[]            poly_ds,    // null or pair of disparity/strength
    		double              vasw_pwr,  // value as weight to this power,
    		int                 debug_level,
    		int                 tileX, // just for debug output
    		int                 tileY
    		)
    {
    	return  corrLMA2Single(       // single tile
        		imgdtt_params,
        		adjust_ly,            // adjust Lazy Eye
        		corr_wnd,             // correlation window to save on re-calculation of the window
        		corr_wnd_inv_limited, // correlation window, limited not to be smaller than threshold - used for finding max/convex areas (or null)
        		corrs,                // may have more elements than pair_mask (corrs may have combo as last elements)
        		disp_dist,            // per camera disparity matrix as a 1d (linescan order)
        		rXY,                  // non-distorted X,Y offset per nominal pixel of disparity
        		pair_mask,            // which pairs to process
        		disp_str,             // -preliminary center x in pixels for largest baseline
        		poly_ds,              // null or pair of disparity/strength
        		vasw_pwr,             // value as weight to this power,
        		null,                 //double []           debug_lma_tile,
        		debug_level,
        		tileX,                // just for debug output
        		tileY);
    }
    
    
    public Corr2dLMA corrLMA2Single( // single tile
    		ImageDttParameters  imgdtt_params,
    		boolean             adjust_ly, // adjust Lazy Eye
    		double [][]         corr_wnd, // correlation window to save on re-calculation of the window
    		double []           corr_wnd_inv_limited, // correlation window, limited not to be smaller than threshold - used for finding max/convex areas (or null)
    		double [][]         corrs, // may have more elements than pair_mask (corrs may have combo as last elements)
    		double [][]         disp_dist, // per camera disparity matrix as a 1d (linescan order)
    		double [][]         rXY, // non-distorted X,Y offset per nominal pixel of disparity
    		boolean []          pair_mask, // which pairs to process
    		double[]            disp_str,   // -preliminary center x in pixels for largest baseline
    		// Not anymore, it is expressed in quad camera disparity (-x offset divided by sqrt(2)!
    		double[]            poly_ds,    // null or pair of disparity/strength
    		double              vasw_pwr,  // value as weight to this power,
    		double []           debug_lma_tile,
    		int                 debug_level,
    		int                 tileX, // just for debug output
    		int                 tileY
    		)
    {
    	// corrs are organized as PAIRS, some are null if not used
    	// for each enabled and available pair find a maximum, filter convex and create sample list
    	boolean need_poly = (disp_str == null); // true; // find initial disparity by polynomial approximation
    	boolean debug_graphic = imgdtt_params.lma_debug_graphic && (imgdtt_params.lma_debug_level1 > 3) && (debug_level > 0) ;
    	debug_graphic |= imgdtt_params.lmamask_dbg && (debug_level > 0);
		String dbg_title = null;
    	if (debug_graphic) {
			dbg_title = String.format("tX%d_tY%d",tileX,tileY);
		}
    	DoubleGaussianBlur gb = null;
    	if (imgdtt_params.lma_sigma > 0) gb = new DoubleGaussianBlur();
    	int center =       transform_size - 1;
    	int corr_size = 2 * transform_size - 1;
    	Corr2dLMA lma = new Corr2dLMA(
    			1,
				1, // int         numMax,
				this, // Correlation2d correlation2d,
    			transform_size,
    			corr_wnd,
    			rXY,                       //double [][] rXY, // non-distorted X,Y offset per nominal pixel of disparity
    			imgdtt_params.lmas_gaussian //boolean gaussian_mode
    			);
    	double []   corr_shape = null;
    	double []   corr_shape_dia = null;
    	double []   norm_shape = null;
    	double [][] pair_shape_masks = null;
    	double [][] pair_offsets = null;
///    	if (imgdtt_params.lmamask_en && (disp_str != null)) {
       	if (disp_str != null) {
    		pair_offsets = lma.getPairsOffsets(
    				corrs,                                   // double [][]       corrs,
    				pair_mask,                               // boolean []        pair_mask,
//    				disp_str[0]/imgdtt_params.lmamask_magic, // double            disparity,
    				// Moved to the caller
///    				disp_str[0]/Math.sqrt(2), // *Math.sqrt(2), // double            disparity,
    				disp_str[0],                             // *Math.sqrt(2), // double            disparity,
    				disp_dist);                              // double [][]       disp_dist);
    		corr_shape = getCorrShape(
    				corrs,  // double [][] corrs,
    				pair_offsets); // double [][] xy_offsets)
    		if (debug_graphic) {
    			int min_dia = 96;
    			double [][] corrs_dia = new double[corrs.length][];
    			for (int i = min_dia; i < corrs.length; i++) {
    				corrs_dia[i] = corrs[i];
    			}
    			corr_shape_dia = getCorrShape(
    					corrs_dia,     // double [][] corrs,
    					pair_offsets); // double [][] xy_offsets)
    		}
    		norm_shape = conditionCorrShape(
    				corr_shape,        // double []   corrs_shape,
    				imgdtt_params.lmamask_min_main,          // double      min_main,
    				imgdtt_params.lmamask_min_neib,          // double      min_neib,
    				imgdtt_params.lmamask_weight_neib,       // double      weight_neib);
    				imgdtt_params.lmamask_weight_neib_neib); // double weight_neib_neib
    		pair_shape_masks = applyCorrShape(
    				norm_shape,    // double []   corrs_shape,
    				pair_offsets); // double [][] xy_offsets)
    	}
    	
    	double [][] dbg_corr =    debug_graphic ? new double [corrs.length][] : null;
    	if (debug_graphic) {
    		(new ShowDoubleFloatArrays()).showArrays(
    				corrs,
    				corr_size,
    				corr_size,
    				true,
    				"corr_pairs"+"_x"+tileX+"_y"+tileY,
    				getCorrTitles());
    		if (corr_shape != null) {
    		(new ShowDoubleFloatArrays()).showArrays(
    				new double [][] {corr_shape,norm_shape, corr_shape_dia},
    				corr_size,
    				corr_size,
    				true,
    				"corr_shape"+"_x"+tileX+"_y"+tileY,
    				new String [] {"corr_shape","norm_shape","corr_shape_dia"});
    		}
    		if (pair_shape_masks != null) {
    		(new ShowDoubleFloatArrays()).showArrays(
    				pair_shape_masks,
    				corr_size,
    				corr_size,
    				true,
    				"corr_shape_masks"+"_x"+tileX+"_y"+tileY,
    				getCorrTitles());
    		}
    	}
    	// try alternative mask generation by accumulation of the pre-shifted (from CM estimation with magic 0.85) correlations

    	double [][] filtWeight =    new double [corrs.length][];
    	double [][] samplesWeight = new double [corrs.length][];
    	int num_disp_samples = 0;
    	int num_cnvx_samples = 0;
    	int num_comb_samples = 0;
    	
   		for (int npair = 0; npair < pair_mask.length; npair++) if ((corrs[npair] != null) && (pair_mask[npair])){
				double [] corr_blur = null;
   			if (imgdtt_params.cnvx_en || (pair_shape_masks == null)) {
   				corr_blur = corrs[npair].clone();
   				if (corr_wnd_inv_limited != null) {
   					for (int i = 0; i < corr_blur.length; i++) {
   						corr_blur[i] *= corr_wnd_inv_limited[i];
   					}
   				}
   				if (imgdtt_params.lma_sigma > 0) {
   					gb.blurDouble(corr_blur, corr_size, corr_size, imgdtt_params.lma_sigma, imgdtt_params.lma_sigma, 0.01);
   				}
   				int imx = imgdtt_params.lma_soft_marg * (corr_size + 1);
   				for (int iy = imgdtt_params.lma_soft_marg; iy < (corr_size - imgdtt_params.lma_soft_marg); iy++) {
   					for (int ix = imgdtt_params.lma_soft_marg; ix < (corr_size - imgdtt_params.lma_soft_marg); ix++) {
   						int indx = iy * corr_size + ix;
   						if (corr_blur[indx] > corr_blur[imx]) imx = indx;
   					}
   				}


   				// filter convex
   				int ix0 = (imx % corr_size) - center; // signed, around center to match filterConvex
   				int iy0 = (imx / corr_size) - center; // signed, around center to match filterConvex
   				filtWeight[npair] =  filterConvex(
   						corr_blur,             // double [] corr_data,
   						imgdtt_params.cnvx_hwnd_size, // int       hwin,
   						ix0,                          // int       x0,
   						iy0,                          // int       y0,
   						imgdtt_params.cnvx_add3x3,    // boolean   add3x3,
   						imgdtt_params.cnvx_weight,    // double    nc_cost,
   						(debug_level > 2));           // boolean   debug);
   			}
   			if (dbg_corr    != null) 	dbg_corr   [npair] = corr_blur;

    	    // Normalize weight for each pair to compensate for different number of convex samples?
   			// Combine/use window masks
   			if (filtWeight[npair] == null) {
   				samplesWeight[npair] = (pair_shape_masks != null)? pair_shape_masks[npair] : null;
   			} else if ((pair_shape_masks == null) || (pair_shape_masks[npair] == null) || !imgdtt_params.lmamask_en) {
   				samplesWeight[npair] = filtWeight[npair]; 
   			} else {
   				samplesWeight[npair] = filtWeight[npair].clone();
   				if (imgdtt_params.cnvx_or) {
   					for (int i = 0; i < samplesWeight[npair].length; i++) {
   						samplesWeight[npair][i] = Math.max(samplesWeight[npair][i], pair_shape_masks[npair][i]); 
   					}
   				} else {
   					for (int i = 0; i < samplesWeight[npair].length; i++) {
   						samplesWeight[npair][i] *= pair_shape_masks[npair][i]; 
   					}
   				}
   			}
   			if (debug_lma_tile != null) { // calculate and return number of non-zero tiles
   				if (pair_shape_masks[npair] != null) {
   					for (int i = 0; i < samplesWeight[npair].length; i++) if (samplesWeight[npair][i] > 0.0) num_disp_samples++;
   				}
   				if (filtWeight[npair] != null) {
   					for (int i = 0; i < filtWeight[npair].length; i++) if (filtWeight[npair][i] > 0.0) num_cnvx_samples++;
   				}
   				if (samplesWeight[npair] != null) {
   					for (int i = 0; i < samplesWeight[npair].length; i++) if (samplesWeight[npair][i] > 0.0) num_comb_samples++;
   				}
   			}
   			
   			
   			
    	    for (int i = 1; i < samplesWeight[npair].length; i++) if (samplesWeight[npair][i] > 0.0) {
    	    	int ix = i % corr_size; // >=0
    	    	int iy = i / corr_size; // >=0
    	    	double v = corrs[npair][i]; // not blurred
    	    	double w = samplesWeight[npair][i];
    	        if (vasw_pwr != 0) {
    	            w *= Math.pow(Math.abs(v), vasw_pwr);
    	        }
    	    	lma.addSample( // x = 0, y=0 - center
    	    			0,    // tile
    	    			npair,
    	    			0,    // int    imax,   // number of maximum (just a hint)
    	    			ix,   // int    x,     // x coordinate on the common scale (corresponding to the largest baseline), along the disparity axis
    	    			iy,   // int    y,     // y coordinate (0 - disparity axis)
    	    			v,    // double v,     // correlation value at that point
    	    			w);   //double w)      // sample weight
    	    }

    	}
   		if (debug_lma_tile != null) { // calculate and return number of non-zero tiles
   			debug_lma_tile[0] = num_disp_samples; 
   			debug_lma_tile[1] = num_cnvx_samples; 
   			debug_lma_tile[2] = num_comb_samples; 
   			debug_lma_tile[3] = -1; // number of LMA iterations 
   			debug_lma_tile[4] = -1; // last number of LMA iterations 
   			debug_lma_tile[5] = -1; // LMA RMA  
   		}
   		if (debug_graphic) {
   			if (dbg_corr != null) {
   				(new ShowDoubleFloatArrays()).showArrays(
   						dbg_corr,
   						corr_size,
   						corr_size,
   						true,
   						"corr_blurred"+"_x"+tileX+"_y"+tileY,
   						getCorrTitles());
   			}
   			if (filtWeight != null) {
   				(new ShowDoubleFloatArrays()).showArrays(
   						filtWeight,
   						corr_size,
   						corr_size,
   						true,
   						"filt_weight"+"_x"+tileX+"_y"+tileY,
   						getCorrTitles());
   			}
   			if (samplesWeight != null) {
   				(new ShowDoubleFloatArrays()).showArrays(
   						samplesWeight,
   						corr_size,
   						corr_size,
   						true,
   						"samples_weights"+"_x"+tileX+"_y"+tileY,
   						getCorrTitles());
   			}
   		}
//    	double [][] disp_str = {{xcenter, 1.0}}; // temporary
    	double [][] disp_str2 = {{0.0, 1.0}}; // temporary // will be calculated/set later
    	if (disp_str != null) {
    		disp_str2[0] = disp_str;
    	}
    	
    	boolean lmaSuccess = false;
    	int num_lma_retries = 0;
		double [] disp = null;
		// adjust_ly
		double [][] ly_offsets_pairs = null;
		if (adjust_ly) {
			ly_offsets_pairs = getPairsCenters(
					corrs, // 		double [][] corrs,
					samplesWeight); // double [][] weights)
		}
		double      step_weight = 0.5; // scale corrections
		double      min_correction = 0.1; //  exit when maximal XY correction is below
		
		while (!lmaSuccess) {
			num_lma_retries ++; // debug
			// FIXME: ugly fix
			double [][] disp_str2_scaled = disp_str2.clone();
			for (int i = 0; i < disp_str2_scaled.length; i++) {
				if (disp_str2_scaled[i] != null) {
					disp_str2_scaled[i] = disp_str2_scaled[i].clone();
					disp_str2_scaled[i][0] /= imgdtt_params.lmamask_magic;
				}
			}
			
    		lma.initVector(
    				null, // 			boolean [] adjust_disparities, // null - adjust all, otherwise - per maximum
    				imgdtt_params.lmas_adjust_wm,  // boolean adjust_width,     // adjust width of the maximum - lma_adjust_wm
    				imgdtt_params.lmas_adjust_ag,  // boolean adjust_scales,    // adjust 2D correlation scales - lma_adjust_ag
    				imgdtt_params.lmas_adjust_wy,  // boolean adjust_ellipse,   // allow non-circular correlation maximums lma_adjust_wy
    				(adjust_ly ? imgdtt_params.lma_adjust_wxy : false), //imgdtt_params.lma_adjust_wxy, // boolean adjust_lazyeye_par,   // adjust disparity corrections parallel to disparities  lma_adjust_wxy
    				(adjust_ly ? imgdtt_params.lma_adjust_ly1: false), // imgdtt_params.lma_adjust_ly1, // boolean adjust_lazyeye_ortho, // adjust disparity corrections orthogonal to disparities lma_adjust_ly1
    				new double[][][] {disp_str2_scaled}, // xcenter, 
    				imgdtt_params.lma_half_width, // double  half_width,       // A=1/(half_widh)^2   lma_half_width
    				(adjust_ly ? imgdtt_params.lma_cost_wy : 0.0), // imgdtt_params.lma_cost_wy,     // double  cost_lazyeye_par,     // cost for each of the non-zero disparity corrections        lma_cost_wy
    				(adjust_ly ? imgdtt_params.lma_cost_wxy : 0.0) //imgdtt_params.lma_cost_wxy     // double  cost_lazyeye_odtho    // cost for each of the non-zero ortho disparity corrections  lma_cost_wxy
    				);
    		
    		lma.setMatrices(disp_dist);
    		lma.initMatrices(); // should be called after initVector and after setMatrices
    		boolean all_sensors_used = lma.setInitialLYOffsets(
    				ly_offsets_pairs,  // double [][] pair_centers,
    				step_weight,       // double      step_weight, // scale corrections
    				min_correction,    // double      min_correction ){  // exit when maximal XY correction is below
    				(debug_level > 0)); // 
    		if (adjust_ly && !all_sensors_used) {
    			return null; //LY requested, but not all sensors present
    		}

    		//center
    		disp = null;
    		if (need_poly) {
    			disp = lma.polyDisparity(
    					corr_wnd_inv_limited,
    					transform_size-1-imgdtt_params.lma_soft_marg,//double max_offset, // 5?
    					debug_graphic?dbg_title:null); // 		double [] rslt = {-approx2d[0], approx2d[2], hwx, hwy};

    			if (disp == null) {
    				if (imgdtt_params.lmas_poly_continue && (disp == null)) {
    					disp = disp_str2[0];
    					if (debug_level > 0) {
    						System.out.println("Poly disparity=NULL, using tile center for initial LMA");
    					}
    				} else {
    					if (debug_level > 0) {
    						System.out.println("Poly disparity=NULL, set lmas_poly_continue to true to use tile center instead");
    					}
    				}
    			} else {
    				disp[1] *= imgdtt_params.lmas_poly_str_scale;
    				disp[0] /= Math.sqrt(2); // disparity is expressed in pixels of a quad camera, combo correlation is for diameter cameras
    				if (debug_level > 0) {
    					System.out.println(String.format("Poly disparity (quad camera scale) =%8.5f , str=%8.5f, disp_str2[0][0]=%8.5f, disp_str2[0][1]=%8.5f",
    							disp[0],disp[1],disp_str2[0][0],disp_str2[0][1]));
    				}
    				if (disp[1] < imgdtt_params.lmas_poly_str_min) {
    					if (debug_level > 0) {
    						System.out.println("Poly strength too low ("+disp[1]+" < "+imgdtt_params.lmas_poly_str_min+")");
    					}
    					disp = null;
    				}
    			}
    			//    		double[]            poly_ds,    // null or pair of disparity/strength
    			if (poly_ds != null) {
    				poly_ds[0] = (disp==null) ? Double.NaN: disp[0];
    				poly_ds[1] = (disp==null) ? 0.0:        disp[1];
    			}
    		} else {
    			disp = disp_str;
    		}
    		
    		
    		
    		if (disp != null) {
    			disp_str2[0] = disp;
    			lma.initDisparity( // USED in lwir null pointer
    					disp_str2, // double [][] disp_str         // initial value of disparity
        				0); // int         nmax); //

    			if (debug_level > 1) {
    				System.out.println("Input data:");
    				lma.printInputDataFx(false);
    				lma.printParams();

    			}
    			lmaSuccess = 	lma.runLma(
    					imgdtt_params.lmas_lambda_initial,     // double lambda,           // 0.1
    					imgdtt_params.lma_lambda_scale_good,   // double lambda_scale_good,// 0.5
    					imgdtt_params.lma_lambda_scale_bad,    // double lambda_scale_bad, // 8.0
    					imgdtt_params.lma_lambda_max,          // double lambda_max,       // 100
    					imgdtt_params.lmas_rms_diff,           // double rms_diff,         // 0.001
    					imgdtt_params.lmas_num_iter,           // int    num_iter,         // 20
    					debug_level);  // imgdtt_params.lma_debug_level1);      // 4); // int    debug_level) // > 3
	    		if (!lmaSuccess && (lma.getBadTile() >= 0)) {
    				if (debug_level > -2) {
    					System.out.println("Found bad tile/pair during single (probably wrong initial maximum - try around preliminary? "+lma.getBadTile());
    				}
	    		} else {
	    			break;
	    		}

    		} else {
    			break;
    		}
    	}

    	if (lmaSuccess) {

    		lma.updateFromVector();


    		double [][] dispStr = lma.lmaDisparityStrength( //TODO: add parameter to filter out negative minimums ?
    				imgdtt_params.lmas_min_amp,      //  minimal ratio of minimal pair correlation amplitude to maximal pair correlation amplitude
    				imgdtt_params.lmas_max_rel_rms,  // maximal relative (to average max/min amplitude LMA RMS) // May be up to 0.3)
    				imgdtt_params.lmas_min_strength, // minimal composite strength (sqrt(average amp squared over absolute RMS)
    				imgdtt_params.lmas_min_ac,       // minimal of A and C coefficients maximum (measures sharpest point/line)
    				imgdtt_params.lmas_min_min_ac,   // minimal of A and C coefficients minimum (measures sharpest point)
    				imgdtt_params.lmas_max_area,     //double  lma_max_area,     // maximal half-area (if > 0.0)
    				imgdtt_params.lma_str_scale,     // convert lma-generated strength to match previous ones - scale
    				imgdtt_params.lma_str_offset     // convert lma-generated strength to match previous ones - add to result
    				);
    		if (dispStr[0][1] <= 0) {
    			lmaSuccess = false;
    			if (debug_level > -2) { // 0
    				System.out.println(String.format("Poly disparity=%8.5f , str=%8.5f", disp[0],disp[1]));
    			}
    			if (debug_lma_tile != null) {
    	    		debug_lma_tile[3] = num_lma_retries; // number of wasted attempts
    	    		debug_lma_tile[4] = lma.getNumIter(); 
    	    	}    			
    		} else {
    			if (debug_level > -2) {
    				System.out.println(String.format("Poly disparity=%8.5f , str=%8.5f, LMA disparity=%8.5f, str=%8.5f",
    						disp[0],disp[1],dispStr[0][0],dispStr[0][1]));
    			}
    			//    		System.out.println("dispStr[0][0]="+dispStr[0][0]+" dispStr[0][1]="+dispStr[0][1]);

    			double [] rms = lma.getRMS();
    			if (debug_lma_tile != null) {
    	    		debug_lma_tile[3] = num_lma_retries; // number of wasted attempts
    	    		debug_lma_tile[4] = lma.getNumIter(); 
    	    		debug_lma_tile[5] = rms[1]; // pure rms 
    	    	}    			
    			
    			if (debug_level > 0) {
    				System.out.println("LMA -> "+lmaSuccess+" RMS="+rms[0]+", pure RMS="+rms[1]);
    				lma.printParams();
    			}

    			if (debug_level > 1) {
    				System.out.println("Input data and approximation:");
    				lma.printInputDataFx(true);
    			}
    			//    				double [][] ds = null;
    			if (debug_graphic && lmaSuccess) {
    				String [] sliceTitles = lma.dbgGetSliceTitles();
 //   				if (corrs.length == 1) { // only for single-tile cluster (here it is always single, and corrs is double [][], not double [][][]
    					(new ShowDoubleFloatArrays()).showArrays(
    							lma.dbgGetSamples(null,0)[0],
    							corr_size,
    							corr_size,
    							true,
    							"corr_values"+"_x"+tileX+"_y"+tileY, sliceTitles);
    					(new ShowDoubleFloatArrays()).showArrays(
    							lma.dbgGetSamples(null,2)[0],
    							corr_size,
    							corr_size,
    							true,
    							"corr_fx"+"_x"+tileX+"_y"+tileY, sliceTitles);
    					(new ShowDoubleFloatArrays()).showArrays(
    							lma.dbgGetSamples(null,1)[0],
    							corr_size,
    							corr_size,
    							true,
    							"corr_weights_late"+"_x"+tileX+"_y"+tileY, sliceTitles);
//    				}
    			}
    		}

    	} else if (debug_lma_tile != null) {
    		debug_lma_tile[3] = num_lma_retries; // number of wasted attempts 
    	}
    	return lmaSuccess? lma: null;
    }
    
    /**
     * Sort array of {disparity,strength} pairs in descending disparity order. Maintain is_lma if not null  
     * @param disp_str_in array of {disparity, strength} pairs (will not be modified)
     * @param is_lma optional array of boolean is_lma elements or null
     * @return reordered array of {disparity, strength} pairs. If not null, is_lma will be updated accordingly
     */
    public static double [][] sortFgFirst(
    		double[][] disp_str_in0,
    		boolean [] is_lma
    		){
    	double[][] disp_str_in = disp_str_in0.clone();
    	boolean [] is_lma1 = (is_lma == null)? null : is_lma.clone();
		int nmax = 0;
		for (int i = 0; i < disp_str_in.length; i++) if (disp_str_in[i] != null) nmax++;
		double [][] disp_str = new double [nmax][];
		for (int n = 0; n < nmax; n++) {
			int idmx = -1;
			for (int i = 0; i < disp_str_in.length; i++) if ((disp_str_in[i] != null) && !Double.isNaN(disp_str_in[i][0])){
				if ((idmx < 0) || ( disp_str_in[i][0] > disp_str_in[idmx][0])) {
					idmx = i;
				}
			}
			if (idmx < 0) {
				// truncate
				double [][] disp_str_trunc = new double [n][];
				for (int i = 0; i <n; i++) {
					disp_str_trunc[i]=disp_str[i];
				}
				return disp_str_trunc;
			}
			disp_str[n] = disp_str_in[idmx];
			if (is_lma != null) {
				is_lma[n] = is_lma1[idmx];
			}
			disp_str_in[idmx] = null;
		}
		return  disp_str;
    }
    
    
    /**
     * Select {disparity, strength} pair from potentially multiple 
     * @param combine_mode processing mode:  0 - both,  1 - strongest, 2 - nearest to zero, 3 - FG, 4 - BG
     * @param disp_str_dual [nmax]{disparity, strength}. One or two pairs of {disparity, strength}, 
     *                      in descending strength order
     * @return selected [nmax]{disparity, strength}. If two are returned, they have the same order as the input one
     */
    public static double [][] selDispStr( // single tile
        	int                 combine_mode,   // 0 - both,  1 - strongest, 2 - nearest to zero, 3 - FG, 4 - BG
    		double[][]          disp_str_dual)  // -preliminary center x in pixels for largest baseline
    {
    	double [][]   disp_str_all; 
    	if (disp_str_dual.length < 2) {
    		disp_str_all =     disp_str_dual;
    	} else { // only 2 max are supported
    		if (disp_str_dual.length > 2) {
    			System.out.println("selDispStr(): Only 2 correlation maximums are currently supported, all but 2 strongest are discarded");
    		}
    		int nearest_max = (Math.abs(disp_str_dual[0][0]) < Math.abs(disp_str_dual[1][0]))? 0 : 1;
    		int fg_max =  (disp_str_dual[0][0] > disp_str_dual[1][0]) ? 0 : 1;
    		switch (combine_mode) {
    		case CAMEL_BOTH: // keep both
    			disp_str_all =     disp_str_dual;
    			break; 
    		case CAMEL_STRONGEST: // keep strongest
    			disp_str_all =     new double [][]   {disp_str_dual    [0]};
    			break;
    		case CAMEL_NEAREST: // keep nearest
    			disp_str_all =     new double [][]   {disp_str_dual    [nearest_max]};
    			break;
    		case CAMEL_FG: // keep foreground
    			disp_str_all =     new double [][]   {disp_str_dual    [fg_max]};
    			break;
    		case CAMEL_BG: // keep background
    			disp_str_all =     new double [][]   {disp_str_dual    [1-fg_max]};
    			break;
    		default: // keep both
    			disp_str_all =     disp_str_dual;
    		}
    	}
        return disp_str_all;
    }
  
    /**
     * Select index of {disparity, strength} pair from potentially multiple 
     * @param combine_mode processing mode:  0 - both,  1 - strongest, 2 - nearest to zero, 3 - FG, 4 - BG
     * @param disp_str_dual [nmax]{disparity, strength}. One or two pairs of {disparity, strength}, 
     *                      in descending strength order
     * @return index of selected pair or -1 if all are selected
     */
    public static int selDispStrIndex( // single tile
    		int                 combine_mode,   // 0 - both,  1 - strongest, 2 - nearest to zero, 3 - FG, 4 - BG
    		double[][]          disp_str_dual) {  // -preliminary center x in pixels for largest baseline
    	if (disp_str_dual.length < 2) {
    		return 0;
    	} else { // only 2 max are supported
    		if (disp_str_dual.length > 2) {
    			System.out.println("selDispStr(): Only 2 correlation maximums are currently supported, all but 2 strongest are discarded");
    		}
    		int nearest_max = (Math.abs(disp_str_dual[0][0]) < Math.abs(disp_str_dual[1][0]))? 0 : 1;
    		int fg_max =  (disp_str_dual[0][0] > disp_str_dual[1][0]) ? 0 : 1;
    		switch (combine_mode) {
    		case CAMEL_BOTH: // keep both
    			return -1;
    		case CAMEL_STRONGEST:// keep strongest
    			return 0;
    		case CAMEL_NEAREST:  // keep nearest
    			return nearest_max;
    		case CAMEL_FG:       // keep foreground
    			return fg_max;
    		case CAMEL_BG:       // keep background
    			return 1-fg_max;
    		default: // keep both
    			return -1;
    		}
    	}
    }
    
    
    /**
     * Process multiple 2D correlation pairs with LMA and extract {disparity,strength} data with
     * Levenberg-Marquardt Algorithm (LMA). The data should be pre-processed (e.g. with
     * rotation+scaling+accumulation of the correlation pairs followed by the polynomial approximation
     * and expected disparity value (or up to 2 values) should be known. In the case of multiple
     * (now just 2) disparities that correspond to foreground (FG) and background (BG) objects in the
     * same tile, this method may process only one of them (strongest, nearest to 0, FG or BG) or both
     * depending on the combine_mode.  
     * @param imgdtt_params multiple processing parameters
     * @param combine_mode processing mode:  0 - both,  1 - strongest, 2 - nearest to zero, 3 - FG, 4 - BG
     * @param corr_wnd correlation window to save on re-calculation of the window
     * @param corr_wnd_inv_limited correlation window, limited not to be smaller than threshold - used
     *                             to find max/convex areas (or null)
     * @param corrs     [pair][pixel] correlation data (per pair, per pixel). May have nulls for
     *                  unused pairs as well as extra elements in the end (such as composite "combo" image
     *                  containing result of accumulation of rotated+scaled individual correlation pairs).
     * @param disp_dist per camera disparity matrix as a 1d (linescan order))
     * @param rXY       non-distorted X,Y offset per nominal pixel of disparity
     * @param pair_mask per pair boolean array, false elements disable corresponding correlation pairs in corrs  
     * @param disp_str_dual [nmax]{disparity, strength}. One or two pairs of {disparity, strength}, 
     *                      in descending strength order
     * @param debug_lma_tile array for per-tile debug data or null if not needed
     * @param debug_level debug level
     * @param tileX     debug tile X (just for the debug images titles
     * @param tileY     debug tile X (just for the debug images titles
     * @return          Corr2dLMA instance if success, or null in case of failure. Corr2dLMA instance may be
     *                  used to extract result (refined disparity/strength pair(s) and other values to judge
     *                  the quality of LMA fitting
     */
    public Corr2dLMA corrLMA2DualMax( // single tile
    		ImageDttParameters  imgdtt_params,
        	int                 combine_mode,   // 0 - both,  1 - strongest, 2 - nearest to zero, 3 - FG, 4 - BG
    		double [][]         corr_wnd, // correlation window to save on re-calculation of the window
    		double []           corr_wnd_inv_limited, // correlation window, limited not to be smaller than threshold - used for finding max/convex areas (or null)
    		double [][]         corrs, // may have more elements than pair_mask (corrs may have combo as last elements)
    		double [][]         disp_dist, // per camera disparity matrix as a 1d (linescan order)
    		double [][]         rXY, // non-distorted X,Y offset per nominal pixel of disparity
    		boolean []          pair_mask, // which pairs to process
    		// should never be null
    		double[][]          disp_str_dual,          // -preliminary center x in pixels for largest baseline
    		double []           debug_lma_tile,
    		int                 debug_level,
    		int                 tileX, // just for debug output
    		int                 tileY)
    {
		int sel_max = selDispStrIndex( // single tile
				imgdtt_params.bimax_combine_mode, // int                 combine_mode,   // 0 - both,  1 - strongest, 2 - nearest to zero, 3 - FG, 4 - BG
				disp_str_dual);                   // double[][]          disp_str_dual)  // -preliminary center x in pixels for largest baseline
    	
    	boolean debug_graphic = imgdtt_params.lma_debug_graphic && (imgdtt_params.lma_debug_level1 > 3) && (debug_level > 0) ;
    	debug_graphic |= imgdtt_params.lmamask_dbg && (debug_level > 0) ;

    	double [][] dbg_corr =    debug_graphic ? new double [corrs.length][] : null;
    	DoubleGaussianBlur gb = null;
    	if (imgdtt_params.lma_sigma > 0) gb = new DoubleGaussianBlur();
    	int center =       transform_size - 1;
    	int corr_size = 2 * transform_size - 1;
    	double [][][]       pair_offsets0 = getPairsOffsets(
    			corrs,         // double [][]       corrs,
    			pair_mask,     // boolean []        pair_mask,
    			disp_str_dual, // double  [][]        disparity_strengths,
    			disp_dist);    // double [][]       disp_dist);
    	double [][][]       own_masks0 = new double [pair_offsets0.length][][];
    	double [][][] lma_corr_weights0 = getLmaWeights(
    			imgdtt_params,                 // ImageDttParameters  imgdtt_params,
    			imgdtt_params.bimax_lpf_neib,  // double              lpf_neib,  // if >0, add ortho neibs (corners - squared)
    			imgdtt_params.bimax_notch_pwr, // double              notch_pwr, //  = 4.00;
    			imgdtt_params.bimax_adv_power, // double              adv_power,    // reduce weight from overlap with adversarial maximum 
    			corrs,                         // double [][]         corrs, // may have more elements than pair_mask (corrs may have combo as last elements)
    			disp_dist,                     // double [][]         disp_dist, // per camera disparity matrix as a 1d (linescan order)
    			disp_dist,                     // double [][]         rXY, // non-distorted X,Y offset per nominal pixel of disparity
    			pair_mask,                     // boolean []          pair_mask, // which pairs to process
    			pair_offsets0,                 // double [][][]       pair_offsets,
    			own_masks0,                    // double [][][]       own_masks, // individual per-maximum, per-pair masks regardless of adversaries or null
    			disp_str_dual,                 // double[][]          disp_str_dual, // single or a pair of {disparity, strength} pairs. First is the strongest of two.
    			debug_level,                   // int                 debug_level,
    			tileX,                         // int                 tileX, // just for debug output
    			tileY);                        // int                 tileY
    	double [][][] lma_corr_weights;
    	double [][]   disp_str_all;
    	double [][][] own_masks;
    	double [][][] pair_offsets;
		boolean [] common_scale; //  = {false,false}; // true}; {true,true}; // // TODO: implement
//		int combine_mode_pre = imgdtt_params.bimax_post_LMA ? 0 : imgdtt_params.bimax_combine_mode;
    	
    	if (lma_corr_weights0.length < 2) {
    		pair_offsets =     pair_offsets0;
    		lma_corr_weights = lma_corr_weights0;
    		disp_str_all =     disp_str_dual;
    		own_masks =        own_masks0;
    		common_scale = new boolean[] {imgdtt_params.bimax_common_fg};
    	} else { // only 2 max are supported
    		if (lma_corr_weights0.length > 2) {
    			System.out.println("corrLMA2DualMax(): Only 2 correlation maximums are currently supported, all but 2 strongest are discarded");
    		}
    		int nearest_max = (Math.abs(disp_str_dual[0][0]) < Math.abs(disp_str_dual[1][0]))? 0 : 1;
    		int fg_max =  (disp_str_dual[0][0] > disp_str_dual[1][0]) ? 0 : 1;
    		if (imgdtt_params.bimax_post_LMA || (sel_max < 0)) { // keep multi
        		pair_offsets =     pair_offsets0;
    			lma_corr_weights = lma_corr_weights0; 
    			disp_str_all =     disp_str_dual;
        		own_masks =        own_masks0;
        		common_scale = new boolean[disp_str_dual.length];
        		for (int i = 0; i < common_scale.length; i++) {
        			common_scale[i] = (i == fg_max) ? imgdtt_params.bimax_common_fg : imgdtt_params.bimax_common_bg; 
        		}
    		} else { // select one max right now
    			//(sel_max == fg_max)
        		pair_offsets =     new double [][][] {pair_offsets0    [sel_max]};
    			lma_corr_weights = new double [][][] {lma_corr_weights0[sel_max]};
    			disp_str_all =     new double [][]   {disp_str_dual    [sel_max]};
    			own_masks =        new double [][][] {own_masks0       [sel_max]};
        		common_scale =     new boolean[]     {(sel_max == fg_max) ? imgdtt_params.bimax_common_fg : imgdtt_params.bimax_common_bg}; 
    		}
    		
    		/*
    		int nearest_max = (Math.abs(disp_str_dual[0][0]) < Math.abs(disp_str_dual[1][0]))? 0 : 1;
    		int fg_max =  (disp_str_dual[0][0] > disp_str_dual[1][0]) ? 0 : 1;
    		
    		switch (combine_mode) {
    		case 0: // keep both
        		pair_offsets =     pair_offsets0;
    			lma_corr_weights = lma_corr_weights0; 
    			disp_str_all =     disp_str_dual;
        		own_masks =        own_masks0;
        		common_scale = new boolean[disp_str_dual.length];
        		for (int i = 0; i < common_scale.length; i++) {
        			common_scale[i] = (i == fg_max) ? imgdtt_params.bimax_common_fg : imgdtt_params.bimax_common_bg; 
        		}
    			break; 
    		case 1: // keep strongest
        		pair_offsets =     new double [][][] {pair_offsets0    [0]};
    			lma_corr_weights = new double [][][] {lma_corr_weights0[0]};
    			disp_str_all =     new double [][]   {disp_str_dual    [0]};
    			own_masks =        new double [][][] {own_masks0       [0]};
        		common_scale =     new boolean[]     {(0 == fg_max) ? imgdtt_params.bimax_common_fg : imgdtt_params.bimax_common_bg}; 
    			break;
    		case 2: // keep nearest
        		pair_offsets =     new double [][][] {pair_offsets0    [nearest_max]};
    			lma_corr_weights = new double [][][] {lma_corr_weights0[nearest_max]};
    			disp_str_all =     new double [][]   {disp_str_dual    [nearest_max]};
    			own_masks =        new double [][][] {own_masks0       [nearest_max]};
        		common_scale =     new boolean[]     {(nearest_max == fg_max) ? imgdtt_params.bimax_common_fg : imgdtt_params.bimax_common_bg}; 
    			break;
    		case 3: // keep foreground
        		pair_offsets =     new double [][][] {pair_offsets0    [fg_max]};
    			lma_corr_weights = new double [][][] {lma_corr_weights0[fg_max]};
    			disp_str_all =     new double [][]   {disp_str_dual    [fg_max]};
    			own_masks =        new double [][][] {own_masks0       [fg_max]};
        		common_scale =     new boolean[]     {imgdtt_params.bimax_common_fg}; 
    			break;
    		case 4: // keep background
        		pair_offsets =     new double [][][] {pair_offsets0    [1-fg_max]};
    			lma_corr_weights = new double [][][] {lma_corr_weights0[1-fg_max]};
    			disp_str_all =     new double [][]   {disp_str_dual    [1-fg_max]};
    			own_masks =        new double [][][] {own_masks0       [1-fg_max]};
        		common_scale =     new boolean[]     {imgdtt_params.bimax_common_bg}; 
    			break;
    		default: // keep both
        		pair_offsets =     pair_offsets0;
    			lma_corr_weights = lma_corr_weights0;
    			disp_str_all =     disp_str_dual;
    			own_masks =        own_masks0;
        		common_scale = new boolean[disp_str_dual.length];
        		for (int i = 0; i < common_scale.length; i++) {
        			common_scale[i] = (i == fg_max) ? imgdtt_params.bimax_common_fg : imgdtt_params.bimax_common_bg; 
        		}
    		}
    		*/
    	}
    	
     	double [][][] filtWeight =    new double [lma_corr_weights.length][corrs.length][];
    	int num_disp_samples = 0;
    	int num_cnvx_samples = 0;
    	int num_comb_samples = 0;
    	int high_marg = corr_size - imgdtt_params.lma_soft_marg -1;
    	boolean [] used_pairs = pair_mask.clone();
    	int num_used_pairs = 0;
   		for (int npair = 0; npair < pair_mask.length; npair++) if ((corrs[npair] != null) && (pair_mask[npair])){
				double [] corr_blur = null;
   			if (imgdtt_params.cnvx_en) { //  || (pair_shape_masks == null)) {
   				corr_blur = corrs[npair].clone();
   				if (corr_wnd_inv_limited != null) {
   					for (int i = 0; i < corr_blur.length; i++) {
   						corr_blur[i] *= corr_wnd_inv_limited[i];
   					}
   				}
   				if (imgdtt_params.lma_sigma > 0) {
   					gb.blurDouble(corr_blur, corr_size, corr_size, imgdtt_params.lma_sigma, imgdtt_params.lma_sigma, 0.01);
   				}
   				if (dbg_corr != null) {
   					dbg_corr[npair] = corr_blur;
   				}
   				for (int nmax = 0; nmax < lma_corr_weights.length; nmax++) { // use same blurred version for all max-es
   					int x00 = (int) Math.round(pair_offsets[nmax][npair][0])+center;
   					int y00 = (int) Math.round(pair_offsets[nmax][npair][1])+center;
   					
   					int x_min = Math.max (x00 - imgdtt_params.bimax_rad_convex_search,imgdtt_params.lma_soft_marg);
   					int x_max = Math.min(x00 + imgdtt_params.bimax_rad_convex_search , high_marg);
   					int y_min = Math.max (y00 - imgdtt_params.bimax_rad_convex_search,imgdtt_params.lma_soft_marg);
   					int y_max = Math.min(y00 + imgdtt_params.bimax_rad_convex_search, high_marg);
   					int imx = x_min + y_min * corr_size;
   					double [] own_mask = own_masks[nmax][npair];
   	   				for (int iy = y_min; iy <= y_max; iy++) {
   	   					for (int ix = x_min; ix <= x_max; ix++) {
   	   						int indx = iy * corr_size + ix;
   	   						if (corr_blur[indx] * own_mask[indx] > corr_blur[imx] * own_mask[imx]) {
   	   							imx = indx;
   	   						}
   	   					}
   	   				}
   	   				// filter convex
   	   				int ix0 = (imx % corr_size) - center; // signed, around center to match filterConvex
   	   				int iy0 = (imx / corr_size) - center; // signed, around center to match filterConvex
   	   				filtWeight[nmax][npair] =  filterConvex(
   	   						corr_blur,             // double [] corr_data,
   	   						imgdtt_params.cnvx_hwnd_size, // int       hwin,
   	   						ix0,                          // int       x0,
   	   						iy0,                          // int       y0,
   	   						imgdtt_params.cnvx_add3x3,    // boolean   add3x3,
   	   						imgdtt_params.cnvx_weight,    // double    nc_cost,
   	   						(debug_level > 2));           // boolean   debug);

   	   				for (int i = 0; i < filtWeight[nmax][npair].length; i++){
   	   					lma_corr_weights[nmax][npair][i] *= filtWeight[nmax][npair][i];
   	   				}
   				}
   			}
   			used_pairs[npair] = true;
   			for (int nmax = 0; nmax < lma_corr_weights.length; nmax++) { // use same blurred version for all max-es
   				int num_pair_samples = 0;
   				for (int i = 0; i < lma_corr_weights[nmax][npair].length; i++) if (lma_corr_weights[nmax][npair][i] > 0.0) {
   					num_pair_samples++;
   				}
   				if (num_pair_samples < imgdtt_params.bimax_min_num_samples) {
   					lma_corr_weights[nmax][npair] = null;
   					used_pairs[npair] = false;
   				}
   			}
   			if (used_pairs[npair]) {
   				num_used_pairs++;
   			}
    	}
   		if (num_used_pairs < imgdtt_params.bimax_min_num_pairs) {
   			return null;
   		}
    	Corr2dLMA lma = new Corr2dLMA(
    			1,
				lma_corr_weights.length, // int         numMax,
				this, // Correlation2d correlation2d,
    			transform_size,
    			corr_wnd,
    			rXY,                       //double [][] rXY, // non-distorted X,Y offset per nominal pixel of disparity
    			imgdtt_params.lmas_gaussian //boolean gaussian_mode
    			);
   		//imgdtt_params.ortho_vasw_pwr
    	double vasw_pwr = imgdtt_params.ortho_vasw_pwr;  // value as weight to this power,
   		for (int npair = 0; npair < pair_mask.length; npair++) if ((corrs[npair] != null) && (used_pairs[npair])){
   			for (int nmax = 0; nmax < lma_corr_weights.length; nmax++) { // use same blurred version for all max-es
   				for (int i = 1; i < lma_corr_weights[nmax][npair].length; i++) if (lma_corr_weights[nmax][npair][i] > 0.0) {
   					int ix = i % corr_size; // >=0
   					int iy = i / corr_size; // >=0
   					double v = corrs[npair][i]; // not blurred
   					double w = lma_corr_weights[nmax][npair][i];
   					if (vasw_pwr != 0) {
   						w *= Math.pow(Math.abs(v), vasw_pwr);
   					}
					lma.addSample( // x = 0, y=0 - center
   							0,     // tile
   							npair,
   							nmax,  // int    imax,   // number of maximum (just a hint)   							
   							ix,    // int    x,     // x coordinate on the common scale (corresponding to the largest baseline), along the disparity axis
   							iy,    // int    y,     // y coordinate (0 - disparity axis)
   							v,     // double v,     // correlation value at that point
   							w);    //double w)      // sample weight
   				}
   			}
   		}
   		
   	// Just debug - modify/restore	
   		if (debug_lma_tile != null) { // calculate and return number of non-zero tiles
   			debug_lma_tile[0] = num_disp_samples; 
   			debug_lma_tile[1] = num_cnvx_samples; 
   			debug_lma_tile[2] = num_comb_samples; 
   			debug_lma_tile[3] = -1; // number of LMA iterations 
   			debug_lma_tile[4] = -1; // last number of LMA iterations 
   			debug_lma_tile[5] = -1; // LMA RMA  
   		}
   		if (debug_graphic) {
   			if (dbg_corr != null) {
   				(new ShowDoubleFloatArrays()).showArrays(
   						dbg_corr,
   						corr_size,
   						corr_size,
   						true,
   						"corr_blurred"+"_x"+tileX+"_y"+tileY,
   						getCorrTitles());
   			}
   			for (int nmax = 0; nmax < filtWeight.length; nmax++) {
   				if (filtWeight[nmax] != null) {
   					(new ShowDoubleFloatArrays()).showArrays(
   							filtWeight[nmax],
   							corr_size,
   							corr_size,
   							true,
   							"filt_weight"+"_x"+tileX+"_y"+tileY+"_n"+nmax,
   							getCorrTitles());
   				}
   			}
   			for (int nmax = 0; nmax < lma_corr_weights.length; nmax++) {
   				if (lma_corr_weights[nmax] != null) {
   					(new ShowDoubleFloatArrays()).showArrays(
   							lma_corr_weights[nmax],
   							corr_size,
   							corr_size,
   							true,
   							"samples_weights"+"_x"+tileX+"_y"+tileY+"_nmx"+nmax,
   							getCorrTitles());
   				}
   			}
   		}

   		// initially will use just first (selected) maximum, both - whenLMA will be modified to support 3 maximums.
    	double [][][] disp_str2 = new double[disp_str_all.length][1][];// = {disp_str_all[0]}; // temporary // will be calculated/set later
    	for (int i = 0; i <disp_str_all.length; i++) if (disp_str_all[i] != null){
    		disp_str2[i][0] = disp_str_all[i];
    	}
    	boolean lmaSuccess = false;
		
    	// When running LMA - first do not touch disparity?
    	boolean [] adjust_disparities = new boolean [disp_str_all.length]; // all false;
    	boolean needprep = true; //t npass = 0;
    	for (int npass = (imgdtt_params.bimax_dual_pass? 0 : 1); npass < 2; npass++) { // may break while (!lmaSuccess) {
    		if (needprep) {
    			lma.preparePars(
    					disp_str2, // double [][][] disp_str_all,  initial value of disparity [max][tile]{disp, strength}
    					imgdtt_params.lma_half_width); // double  half_width,       // A=1/(half_widh)^2   lma_half_width
    			needprep = false;
        		lma.setMatrices (disp_dist);
        		lma.initMatrices (); // should be called after initVector and after setMatrices
        		lma.initDisparity(disp_str2); // initial value of disparity [max][tile]{disp, strength}
    		} else {
    			lma.updateFromVector();
    		}
    		if (npass > 0) {
    			adjust_disparities = null;
    		}
    		lma.setParMask(
    				adjust_disparities, // null, // 			boolean [] adjust_disparities, // null - adjust all, otherwise - per maximum
    				common_scale,       //boolean [] common_scale,       // per-maximum, if true - common scale for all pairs
    				imgdtt_params.lmas_adjust_wm,  // boolean adjust_width,     // adjust width of the maximum - lma_adjust_wm
    				imgdtt_params.lmas_adjust_ag,  // boolean adjust_scales,    // adjust 2D correlation scales - lma_adjust_ag
    				imgdtt_params.lmas_adjust_wy,  // boolean adjust_ellipse,   // allow non-circular correlation maximums lma_adjust_wy
    				false, // (adjust_ly ? imgdtt_params.lma_adjust_wxy : false), //imgdtt_params.lma_adjust_wxy, // boolean adjust_lazyeye_par,   // adjust disparity corrections parallel to disparities  lma_adjust_wxy
    				false, // (adjust_ly ? imgdtt_params.lma_adjust_ly1: false), // imgdtt_params.lma_adjust_ly1, // boolean adjust_lazyeye_ortho, // adjust disparity corrections orthogonal to disparities lma_adjust_ly1
    				0.0, // (adjust_ly ? imgdtt_params.lma_cost_wy : 0.0), // imgdtt_params.lma_cost_wy,     // double  cost_lazyeye_par,     // cost for each of the non-zero disparity corrections        lma_cost_wy
    				0.0);  // (adjust_ly ? imgdtt_params.lma_cost_wxy : 0.0) //imgdtt_params.lma_cost_wxy     // double  cost_lazyeye_odtho    // cost for each of the non-zero ortho disparity corrections  lma_cost_wxy
    		if (debug_level > 1) { // 1) {
    			System.out.println("Input data:");
    			lma.printInputDataFx(false);
    			lma.printParams();
    		}
    		lmaSuccess = lma.runLma(
    				imgdtt_params.lmas_lambda_initial,     // double lambda,           // 0.1
    				imgdtt_params.lma_lambda_scale_good,   // double lambda_scale_good,// 0.5
    				imgdtt_params.lma_lambda_scale_bad,    // double lambda_scale_bad, // 8.0
    				imgdtt_params.lma_lambda_max,          // double lambda_max,       // 100
    				imgdtt_params.lmas_rms_diff,           // double rms_diff,         // 0.001
    				imgdtt_params.lmas_num_iter,           // int    num_iter,         // 20
    				debug_level);  // imgdtt_params.lma_debug_level1);      // 4); // int    debug_level) // > 3
    		if (!lmaSuccess && (lma.getBadTile() >= 0)) {
    			if (debug_level > -2) {
    				System.out.println("Found bad tile/pair during single (probably wrong initial maximum - try around preliminary? "+lma.getBadTile());
    			}
    		}
    	}

		if (lmaSuccess) {
			lma.updateFromVector();
			double [][][] dispStrs = lma.lmaDisparityStrengths( //TODO: add parameter to filter out negative minimums ?
					imgdtt_params.lmas_min_amp,      //  minimal ratio of minimal pair correlation amplitude to maximal pair correlation amplitude
					imgdtt_params.lmas_min_amp_bg,   //  minimal ratio of minimal pair correlation amplitude to maximal pair correlation amplitude
					imgdtt_params.lmas_max_rel_rms,  // maximal relative (to average max/min amplitude LMA RMS) // May be up to 0.3)
					imgdtt_params.lmas_min_strength, // minimal composite strength (sqrt(average amp squared over absolute RMS)
					imgdtt_params.lmas_min_ac,       // minimal of A and C coefficients maximum (measures sharpest point/line)
					imgdtt_params.lmas_min_min_ac,   // minimal of A and C coefficients minimum (measures sharpest point)
					imgdtt_params.lmas_max_area,     //double  lma_max_area,     // maximal half-area (if > 0.0)
					imgdtt_params.lma_str_scale,     // convert lma-generated strength to match previous ones - scale
					imgdtt_params.lma_str_offset     // convert lma-generated strength to match previous ones - add to result
					);
			for (int nmax = 00; nmax < dispStrs.length; nmax++) if (dispStrs[nmax][0][1] <= 0) {
				lmaSuccess = false;
				break;
			}
			if (!lmaSuccess) {
				if (debug_lma_tile != null) {
//					debug_lma_tile[3] = num_lma_retries; // number of wasted attempts
					debug_lma_tile[4] = lma.getNumIter(); 
				}    			
			} else {
				if (debug_level > -2) {
					for (int nmax = 0; nmax < dispStrs.length; nmax++) {
					System.out.println(String.format("LMA max-%d: disparity=%8.5f, str=%8.5f",
							nmax,dispStrs[nmax][0][0],dispStrs[nmax][0][1]));
					}
				}
				double [] rms = lma.getRMS();
				if (debug_lma_tile != null) {
//					debug_lma_tile[3] = num_lma_retries; // number of wasted attempts
					debug_lma_tile[4] = lma.getNumIter(); 
					debug_lma_tile[5] = rms[1]; // pure rms 
				}    			

				if (debug_level > 0) {
					System.out.println("LMA -> "+lmaSuccess+" RMS="+rms[0]+", pure RMS="+rms[1]);
					lma.printParams();
				}

				if (debug_level > 1) {
					System.out.println("Input data and approximation:");
					lma.printInputDataFx(true);
				}
				if (debug_graphic && lmaSuccess) {
					String [] sliceTitles = lma.dbgGetSliceTitles();
					(new ShowDoubleFloatArrays()).showArrays(
							lma.dbgGetSamples(null,0)[0],
							corr_size,
							corr_size,
							true,
							"corr_values"+"_x"+tileX+"_y"+tileY, sliceTitles);
					(new ShowDoubleFloatArrays()).showArrays(
							lma.dbgGetSamples(null,2)[0],
							corr_size,
							corr_size,
							true,
							"corr_fx"+"_x"+tileX+"_y"+tileY, sliceTitles);
					(new ShowDoubleFloatArrays()).showArrays(
							lma.dbgGetSamples(null,1)[0],
							corr_size,
							corr_size,
							true,
							"corr_weights_late"+"_x"+tileX+"_y"+tileY, sliceTitles);
				}
			}

		} else if (debug_lma_tile != null) {
 //   		debug_lma_tile[3] = num_lma_retries; // number of wasted attempts 
    	}
    	return lmaSuccess? lma: null;
    }
   /**
    * Low-pass filtering of correlation data with a 3x3 kernel
    * @param data correlation data in scan-line order (will be modified), currently 15x15=225 long
    * @param orth_val relative (to the center pixel) weight of the 4 ortho pixels (corner weights
    *                 are orth_val * orth_val 
    */
    void lpf_neib(
    		double [] data,
    		double orth_val) {
    	double scale = 1.0/(1.0 + 4 * orth_val * (1.0 + orth_val));
    	int [] ortho = {-1, - corr_size, 1, corr_size};
    	int [] corners = {- corr_size - 1, - corr_size + 1,  corr_size -1, corr_size + 1};
    	double [] scales = {scale, scale*orth_val, scale * orth_val * orth_val}; // center, ortho,corner
    	double [] d = data.clone();
    	Arrays.fill(data, 0.0);
    	int i_last = data.length - corr_size - 1;
    	for (int i = corr_size+1; i < i_last; i++) { // assuming nothing in a 1-pixel frame
    		if (d[i] != 0.0) {
    			data[i] += d[i] * scales[0];
				double a = d[i]* scales[1];
    			for (int offs:ortho) {
    				data[i+offs] += a;	
    			}
				a = d[i]* scales[2];
    			for (int offs:corners) {
    				data[i+offs] += a;	
    			}
    		}
    	}
    }
    
	/**
	 * Calculate offsets for multiple disparities, for each defined pair x,y expected offset, assuming only disparity, not lazy eye   
	 * @param corrs - pairs 2D correlations (each in scanline order) - just to determine null/non-null
	 * @param disparity_strengths - expected {disparity, strength} pairs (e.g. from CM)
	 * @param disp_dist - per camera disparity matrix as a 1d (linescan order))
	 * @return array of per-disparity, per pair x,y expected center offset in 2D correlations or nulls for undefined pairs (or null if already set)
	 */
	public double [][][] getPairsOffsets(
			double [][]       corrs,
			boolean []        pair_mask,
			double [][]       disparity_strengths,
			double [][]       disp_dist){ //
		double [][][] pairs_offsets = new double [disparity_strengths.length][][];
		for (int i = 0; i < pairs_offsets.length; i++) {
			pairs_offsets[i] = getPairsOffsets(
					corrs,                     // double [][]       corrs,
					pair_mask,                 // boolean []        pair_mask,
					disparity_strengths[i][0], // double            disparity,
					disp_dist);                // double [][]       disp_dist)
		}
		return pairs_offsets;
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
		Matrix [] m_disp = Corr2dLMA.getMatrices(
				disp_dist, // double [][] am_disp,
				getNumSensors()); // int num_cams)
		for (int pair = 0; pair < xy_offsets.length; pair++) if ((pair < getNumPairs()) && (corrs[pair] != null) && ((pair_mask == null) || pair_mask[pair])){ // OK to calculate for each
 			int [] fscam = getPair(pair); // returns [first_cam, second_cam]
			Matrix mdd_dnd = new Matrix(new double[] {-disparity,  0.0},2);
			Matrix xcam_ycam_f = m_disp[fscam[0]].times(mdd_dnd);
			Matrix xcam_ycam_s = m_disp[fscam[1]].times(mdd_dnd);
			xy_offsets[pair] = xcam_ycam_f.minus(xcam_ycam_s).getColumnPackedCopy();	
		}
		return xy_offsets;
	}

    /**
     * Generate masks (analog weights) per correlation maximum, per pair taking into account
     * other (adversarial) correlation maximum(s) and reducing weigh.
     * @param imgdtt_params multiple processing parameters
     * @param lpf_neib low pass filtering parameter
     * @param notch_pwr power to raise notch filter to generate an accumulated over all pairs
     *                  notch filter that later is shifted according to the expected disparity
     *                  for each individual pair. The higher the power, the wider and sharper
     *                  is the filter.
     * @param adv_power power to raise relative (to adversarial maximum) this maximum strength
     *                  for each pixel. Higher power make filtering sharper
     * @param corrs     [pair][pixel] correlation data (per pair, per pixel). May have nulls for
     *                  unused pairs.  
     * @param disp_dist per camera disparity matrix as a 1d (linescan order))
     * @param rXY       non-distorted X,Y offset per nominal pixel of disparity
     * @param pair_mask per pair boolean array, false elements disable corresponding correlation pairs in corrs  
     * @param pair_offsets [nmax][pair]{x,y} per-maximum, per pair X,Y offset of the correlation maximum 
     * @param own_masks [nmax][pair][pix] per-maximum, per-pair correlation pixel array calculated regardless of
     *                  the adversaries. Should be initialized to double[nmax][][] or null (to skip calculation)
     * @param disp_str_dual [nmax]{disparity, strength}. One or two pairs of {disparity, strength}, 
     *                  in descending strength order
     * @param debug_level debug level
     * @param tileX     debug tile X (just for the debug images titles
     * @param tileY     debug tile X (just for the debug images titles
     * @return [nmax][pair][pixel] corresponding by the first dimension to the input disp_str_dual array
     */
    
	public double[][][] getLmaWeights(
    		ImageDttParameters  imgdtt_params,
    		double              lpf_neib,  // if >0, add ortho neibs (corners - squared)
    		double              notch_pwr, //  = 4.00;
    		double              adv_power,    // reduce weight from overlap with adversarial maximum 
    		double [][]         corrs, // may have more elements than pair_mask (corrs may have combo as last elements)
    		double [][]         disp_dist, // per camera disparity matrix as a 1d (linescan order)
    		double [][]         rXY, // non-distorted X,Y offset per nominal pixel of disparity
    		boolean []          pair_mask, // which pairs to process
    		double [][][]       pair_offsets,
    		double [][][]       own_masks, // individual per-maximum, per-pair masks regardless of adversaries or null
    		double[][]          disp_str_dual, // single or a pair of {disparity, strength} pairs. First is the strongest of two.
    		int                 debug_level,
    		int                 tileX, // just for debug output
    		int                 tileY){
    	boolean debug_graphic = imgdtt_params.lma_debug_graphic && (imgdtt_params.lma_debug_level1 > 3) && (debug_level > 0) ;
    	debug_graphic |= imgdtt_params.lmamask_dbg && (debug_level > 0); //  || true;
    	String dbg_title = null;
    	if (debug_graphic) {
    		dbg_title = String.format("tX%d_tY%d",tileX,tileY);
    	}
    	double [][]   corr_shapes =       new double [disp_str_dual.length][];
    	double [][]   corr_shapes_dia =   new double [disp_str_dual.length][];
    	double [][]   norm_shapes =       new double [disp_str_dual.length][];
    	double [][]   corr_sharps =       new double [disp_str_dual.length][]; // combination of corr_shapes and norm_shapes
    	double [][][] pair_shapes_masks = new double [disp_str_dual.length][][];
    	double [][][] pair_sharp_masks =  new double [disp_str_dual.length][][];
    	boolean is_multi = disp_str_dual.length > 1;
    	double [][]   notch = is_multi ? (new double [corrs.length][]) : null;
    	double [][]   corrs_notched;
    	double []     notch_shape =          null;
    	for (int nmax = 0; nmax < disp_str_dual.length; nmax++) {
    		if (is_multi && (nmax > 0)) {
    			corrs_notched = new double[notch.length][]; 
    			for (int nt = 0; nt < notch.length; nt++) {
    				if (corrs[nt] != null) {
    					if (notch[nt] != null) {
    						corrs_notched[nt] = corrs[nt].clone();
    						for (int i = 0; i < notch[nt].length; i++) {
    							corrs_notched[nt][i] *= notch[nt][i];
    						}
    					} else {
    						corrs_notched[nt] = corrs[nt];
    					}
    				}
    			}
    		} else {
    			corrs_notched = corrs;
    		}
    		corr_shapes[nmax] = getCorrShape(
    				corrs_notched,        // double [][] corrs,
    				pair_offsets[nmax]);  // double [][] xy_offsets)
    		if (debug_graphic) {
    			int min_dia = 96;
    			double [][] corrs_dia = new double[corrs.length][];
    			for (int i = min_dia; i < corrs.length; i++) {
    				corrs_dia[i] = corrs[i];
    			}
    			corr_shapes_dia[nmax] = getCorrShape(
    					corrs_dia,     // double [][] corrs,
    					pair_offsets[nmax]); // double [][] xy_offsets)
    		}
    		norm_shapes[nmax] = conditionCorrShape(
    				corr_shapes[nmax],        // double []   corrs_shape,
    				imgdtt_params.lmamask_min_main,          // double      min_main,
    				imgdtt_params.lmamask_min_neib,          // double      min_neib,
    				imgdtt_params.lmamask_weight_neib,       // double      weight_neib);
    				imgdtt_params.lmamask_weight_neib_neib); // double weight_neib_neib
    		pair_shapes_masks[nmax] = applyCorrShape(
    				norm_shapes[nmax],    // double []   corrs_shape,
    				pair_offsets[nmax]); // double [][] xy_offsets)
    		
    		// multiply shape and norm_shape and normalize
    		corr_sharps[nmax] = new double [norm_shapes[nmax].length];
       		for (int i = 0; i < corr_sharps[nmax].length; i++) if (norm_shapes[nmax][i] > 0.0){
    			corr_sharps[nmax][i] = norm_shapes[nmax][i] * Math.abs(corr_shapes[nmax][i]); // abs() as we are interested in adversarial influence
    		}
    		if (lpf_neib > 0) { // lpf filter
    			lpf_neib(
    					corr_sharps[nmax], // double [] data,
    					lpf_neib);         // double orth_val)     			
        		for (int i = 0; i < corr_sharps[nmax].length; i++) if ( norm_shapes[nmax][i] <= 0.0){
        			corr_sharps[nmax][i] = 0.0; // remove possible roll-overs in a first/last line
        		}
    		}
    		double dmax = 0.0;
       		for (int i = 0; i < corr_sharps[nmax].length; i++) if (norm_shapes[nmax][i] > 0.0){
    			if (corr_sharps[nmax][i] > dmax) {
    				dmax = corr_sharps[nmax][i]; 
    			}
    		}
    		for (int i = 0; i < corr_sharps[nmax].length; i++){
    			corr_sharps[nmax][i]/= dmax;
    		}
    		
    		pair_sharp_masks[nmax]= applyCorrShape(
    				corr_sharps[nmax],    // double []   corrs_shape,
    				pair_offsets[nmax]); // double [][] xy_offsets)
    		// FIXME: Some duplicate below - use corr_sharps[nmax]
    		if (is_multi && (nmax == 0)) {
//    			notch_shape = corr_shapes[0].clone();
    			notch_shape = new double [corr_shapes[0].length];
    			for (int i = 0; i < notch_shape.length; i++) {
    				if (corr_sharps[nmax][i] <= 0.0) {
    					notch_shape[i] = 1.0;
    				} else {
//    					notch_shape[i] = Math.pow(( 1.0 - notch_shape[i] / dmax * norm_shapes[0][i]), notch_pwr);
    					notch_shape[i] = Math.pow( 1.0 - corr_sharps[nmax][i],  notch_pwr);
    				}
    			}
    			notch = applyCorrShape(
    					notch_shape,         // double []   corrs_shape,
    					pair_offsets[0]); // double [][] xy_offsets)
    		}
    	}
    	// Apply adversarial selections. Only two maximums are considered here
    	if (disp_str_dual.length == 2) {
    		for (int nmax = 0; nmax < 2; nmax ++) {
        		double over_avd_scale = disp_str_dual[nmax][1] / disp_str_dual[1 - nmax][1];
    			for (int npair = 0; npair < pair_shapes_masks[nmax].length; npair++) {
    				if ((pair_shapes_masks[nmax][npair] != null) && (pair_sharp_masks[nmax][npair] != null) && (pair_sharp_masks[1-nmax][npair] != null)) {
    					for (int i = 0; i < pair_shapes_masks[nmax][npair].length; i++) if (pair_shapes_masks[nmax][npair][i] > 0.0){
    						if (pair_sharp_masks[1 - nmax][npair][i] > 0.0) { // adversary exists here
        						double relative_strength = over_avd_scale * 
        								pair_sharp_masks[nmax][npair][i] /
        								pair_sharp_masks[1 - nmax][npair][i];
        						double scale = relative_strength / (1.0 + relative_strength);
        						if (adv_power != 1.0) {
        							scale = Math.pow(scale, adv_power);
        							pair_shapes_masks[nmax][npair][i] *= scale;
        						}
    						}
    					}
    				}
    			}
    		}
    	}
    	
    	if (debug_graphic) {
    		(new ShowDoubleFloatArrays()).showArrays(
    				corrs,
    				corr_size,
    				corr_size,
    				true,
    				"corr_pairs"+"_x"+tileX+"_y"+tileY,
    				getCorrTitles());
    		if (corr_shapes != null) {
    			for (int nmax = 0; nmax < disp_str_dual.length; nmax++) {
    				(new ShowDoubleFloatArrays()).showArrays(
    						new double [][] {corr_shapes[nmax],corr_sharps[nmax], norm_shapes[nmax], corr_shapes_dia[nmax],notch_shape},
    						corr_size,
    						corr_size,
    						true,
    						"corr_shape"+"_x"+tileX+"_y"+tileY+"_M"+nmax,
    						new String [] {"corr_shape","corr_sharp","norm_shape","corr_shape_dia","notch_shape"});
    			}
    		}
    		//    		if (pair_shapes_masks != null) {
    		for (int nmax = 0; nmax < disp_str_dual.length; nmax++) {
    			(new ShowDoubleFloatArrays()).showArrays(
    					pair_shapes_masks[nmax],
    					corr_size,
    					corr_size,
    					true,
    					"corr_shape_masks"+"_x"+tileX+"_y"+tileY+"_M"+nmax,
    					getCorrTitles());
    		}

    		for (int nmax = 0; nmax < disp_str_dual.length; nmax++) {
    			(new ShowDoubleFloatArrays()).showArrays(
    					pair_sharp_masks[nmax],
    					corr_size,
    					corr_size,
    					true,
    					"pair_sharp_masks"+"_x"+tileX+"_y"+tileY+"_M"+nmax,
    					getCorrTitles());
    		}
    		if (notch != null) {
    				(new ShowDoubleFloatArrays()).showArrays(
    						notch,
    						corr_size,
    						corr_size,
    						true,
    						"notch"+"_x"+tileX+"_y"+tileY,
    						getCorrTitles());
    		}
    	}
    	if (own_masks != null) {
    		for (int nmax = 0; nmax < own_masks.length; nmax++) {
    			own_masks[nmax] = pair_sharp_masks[nmax];
    		}
    	}
    	return pair_shapes_masks;

    }
    
    
    
    /**
     * Condition correlation shape to use as a mask
     * @param corrs_shape raw centered correlation shape calculated with getCorrShape()
     * @param min_main threshold to keep points regardless of neighbors (fraction of maximum)
     * @param min_neib threshold to keep points if they have neighbors above min_main (fraction of maximum)
     * @param weight_neib weight of the points that satisfy min_neib, but not min_main condition, < 1.0
     * @param weight_neib_neib weight of neighbors of neighbors, regardless of value, < weight_neib
     * @return weights array with 1.0 for points above min_main, weight_neib for satisfying min_neib and 0 - for others
     */
    double [] conditionCorrShape(
    		double []   corrs_shape,
    		double      min_main,
    		double      min_neib,
    		double      weight_neib,
    		double      weight_neib_neib) {
    	
    	double []   cond_shape = new double [corrs_shape.length];
    	double weight_main = 1.0;
    	int corr_size = 2 * transform_size - 1;
//    	TileNeibs tn = new TileNeibs(corr_size,corr_size);
    	int [][] offs_xy= {{0,-1},{1,-1},{1,0},{1,1},{0,1},{-1,1},{-1,0},{-1,-1}}; 
    	double max =  corrs_shape[0];
    	for (int i = 1; i < cond_shape.length; i++) {
    		if (corrs_shape[i] > max) max = corrs_shape[i];  
    	}
    	min_main *= max;
    	min_neib *= max;
    	for (int i = 0; i < cond_shape.length; i++) if (corrs_shape[i] >= min_main) {
    			cond_shape[i] = weight_main;
   		}
    	if (weight_neib > 0.0) {
    		for (int y = 0; y < corr_size; y++) {
    			for (int x = 0; x < corr_size; x++) {
    				int indx = y * corr_size + x;
    				if ((cond_shape[indx] == 0.0) &&(corrs_shape[indx] >= min_neib)) {
    					for (int [] dxy : offs_xy) {
    						int y1 = y + dxy[1];
    						if ((y1 >= 0) && (y1 < corr_size)) { 
    							int x1 = x + dxy[0];
    							if ((x1 >= 0) && (x1 < corr_size)) {
    								if (cond_shape[y1 * corr_size + x1] == weight_main) {
    									cond_shape[y * corr_size + x] = weight_neib;
    									break;
    								}
    							}
    						}
    					}
    				}
    			}
    		}
    	}
    	if (weight_neib_neib > 0) {
    		for (int y = 0; y < corr_size; y++) {
    			for (int x = 0; x < corr_size; x++) {
    				int indx = y * corr_size + x;
    				if (cond_shape[indx] == 0.0) {
    					for (int [] dxy : offs_xy) {
    						int y1 = y + dxy[1];
    						if ((y1 >= 0) && (y1 < corr_size)) { 
    							int x1 = x + dxy[0];
    							if ((x1 >= 0) && (x1 < corr_size)) {
    								if (cond_shape[y1 * corr_size + x1] >= weight_neib) {
    									cond_shape[y * corr_size + x] = weight_neib_neib;
    									break;
    								}
    							}
    						}
    					}
    				}
    			}
    		}
    	}
    	return cond_shape;
    }
    
    /**
     * Calculate averaged (for defined pairs) 2D correlation shape assuming only
     * disparity, no lazy eye. Disparity may be estimated with averaging + CM 
     * @param corrs pairs 2D correlations (each in scanline order)  
     * @param xy_offsets per pair x,y offsets calculated with Corr2dLMA.getPairsCenters()
     * @return averaged and centered (by applying disparity-defined shift) 2D correlations
     * that can be used as mask for LMA input samples
     */
    double [] getCorrShape(
    		double [][] corrs,
			double [][] xy_offsets){ // 
    	int corr_size = 2 * transform_size - 1;
    	double [] corr_shape =   new double [corr_size*corr_size];
    	double [] corr_weights = new double [corr_size*corr_size];
    	for (int np = 0; np < corrs.length; np++) if ((xy_offsets[np] != null) && (corrs[np] != null)) {
    		int ix0 = (int) Math.floor(xy_offsets[np][0]);
    		int iy0 = (int) Math.floor(xy_offsets[np][1]);
    		for (int dy = 0; dy < 2; dy++) {
    			int iy = iy0 + dy;
    			int y0 = (iy > 0) ? 0 : -iy;
    			int y1 = (iy > 0) ? (corr_size - iy) : corr_size;
    			double ky = (dy > 0)? (xy_offsets[np][1] - iy0) : (iy0 + 1 - xy_offsets[np][1]);
        		for (int dx = 0; dx < 2; dx++) {
        			int ix = ix0 + dx;
        			int x0 = (ix > 0) ? 0 : -ix;
        			int x1 = (ix > 0) ? (corr_size - ix) : corr_size;
        			double kx = (dx > 0)? (xy_offsets[np][0] - ix0) : (ix0 + 1 - xy_offsets[np][0]);
        			double k = ky*kx;
        			int dsrc = iy * corr_size + ix;
        			for (int y = y0; y < y1; y++) {
            			for (int x = x0; x < x1; x++) {
            				int idst = y * corr_size + x;
            				int isrc = idst + dsrc;
            				corr_weights[idst] += k;
            				corr_shape[idst] += corrs[np][isrc] * k;
            			}
        			}
        		}
    		}
    	}
    	for (int i = 0; i < corr_shape.length; i++) if (corr_weights[i] > 0.0) {
    		corr_shape[i] /= corr_weights[i];
    	}
    	return corr_shape;
    }
    
    /**
     * Shift common corrs_shape in reverse direction of xy_offsets to create per-pair selection window for LMA fitting
     * @param corrs_shape centered correlation shape calculated by averaging shifted (according to disparity) 2D correlations
     * @param xy_offsets pairs of per-pair shift of 2D correlations according to common disparity
     * @return per-pair selection windows for LMA
     */
    double [][] applyCorrShape(
    		double []   corrs_shape,
			double [][] xy_offsets){ //
    	double [][] shifted_shapes = new double [xy_offsets.length][];
    	int corr_size = 2 * transform_size - 1;
    	int corr_len = corr_size*corr_size; 

    	for (int np = 0; np < xy_offsets.length; np++) if (xy_offsets[np]!= null){
    		double [] corr_weights = new double [corr_len];
    		shifted_shapes[np] = new double [corr_len];
    		int ix0 = (int) Math.floor(-xy_offsets[np][0]);
    		int iy0 = (int) Math.floor(-xy_offsets[np][1]);
    		for (int dy = 0; dy < 2; dy++) {
    			int iy = iy0 + dy;
    			int y0 = (iy > 0) ? 0 : -iy;
    			int y1 = (iy > 0) ? (corr_size - iy) : corr_size;
    			double ky = (dy > 0)? (-xy_offsets[np][1] - iy0) : (iy0 + 1 + xy_offsets[np][1]);
        		for (int dx = 0; dx < 2; dx++) {
        			int ix = ix0 + dx;
        			int x0 = (ix > 0) ? 0 : -ix;
        			int x1 = (ix > 0) ? (corr_size - ix) : corr_size;
        			double kx = (dx > 0)? (-xy_offsets[np][0] - ix0) : (ix0 + 1 + xy_offsets[np][0]);
        			double k = ky*kx;
        			int dsrc = iy * corr_size + ix;
        			for (int y = y0; y < y1; y++) {
            			for (int x = x0; x < x1; x++) {
            				int idst = y * corr_size + x;
            				int isrc = idst + dsrc;
            				corr_weights[idst] += k;
            				shifted_shapes[np][idst] += corrs_shape[isrc]  *k; 
            			}
        			}
        		}
    		}
    		for (int i = 0; i < corr_len; i++) if (corr_weights[i] > 0) {
    			shifted_shapes[np][i] /= corr_weights[i];
    		}
    	}
    	return shifted_shapes;
    }
    
    
    /**
     * Find individual pair centers for Lazy Eye initialization
     * @param corrs per-pair 2d correlations
     * @param weights per-pair, per-point sample weights 
     * @return per-pair x,y offsets
     */
    public double [][] getPairsCenters(
    		double [][] corrs,
    		double [][] weights){
    	double [][] xy_offsets_pairs = new double[corrs.length][]; // 
    	for (int np = 0; np < corrs.length; np++) if (corrs[np]!= null) {
            double [] wcorr = new double [corrs[np].length];
        	double  pair_weights = 0.0;
            for (int i = 0; i < wcorr.length; i++) if (weights[np][i] > 0.0){
            	wcorr[i] = corrs[np][i]*weights[np][i];
            	pair_weights += wcorr [i];
            }
            double [] dxy = getMaxXYCm(
            		wcorr,   // double [] data,
                    getCombWidth(), //       data_width,
                    getCombHeight()/2 - getCombOffset(), // int       center_row
        			false); // boolean   debug)
            xy_offsets_pairs[np] = new double[3];
            xy_offsets_pairs[np][0] = dxy[0];
            xy_offsets_pairs[np][1] = dxy[1];
            xy_offsets_pairs[np][2] = pair_weights;
    	}
    	return xy_offsets_pairs;
    }
   
    
    
    public Correlations2dLMA corrLMA( // USED in lwir
    		ImageDttParameters  imgdtt_params,
    		double [][]         corrs,
//    		int                 pair_mask, // which pairs to process
    		boolean []          pair_mask, // which pairs to process
    		boolean             run_poly_instead, // true - run LMA, false - run 2d polynomial approximation
    		double              xcenter,   // preliminary center x in pixels for largest baseline
    		double              vasw_pwr,  // value as weight to this power,
    		int                 debug_level,
    		int                 tileX, // just for debug output
    		int                 tileY
    		)
    {
    	// for quad camera
    	int [][] groups = new int [GROUPS.length][];
    	int   [] scale_ind = new int  [GROUPS.length];
    	// See which groups exist for current pairs mask
    	int ng = 0;
    	ArrayList<Integer> sl = new ArrayList<Integer>();
    	for (int i = 0; i < GROUPS.length; i++) {
    		groups[i] = getNumberBaseOfCompatiblePairs(
    				corrs,                       // double [][] correlations,
    				pair_mask,                   // int         pairs_mask,
    				(GROUPS[i][0] > 0),          // boolean     diagonal,
    				GROUPS[i][1]);               // int         baseline_scale
    		if (groups[i][0] > 0) {
    			ng++;
    			if (!sl.contains(GROUPS[i][1])) {
    				sl.add(GROUPS[i][1]);
    			}
    			scale_ind[i] = sl.indexOf(GROUPS[i][1]);
    		}
    	}
    	if (debug_level > 1) {
    		System.out.println("corrLMA(): found "+ng+" groups, "+sl.size()+" scales");
    	}
    	double [][] groups_LMA =      new double [ng][];
    	int    [][] groups_pairs =    new int [ng][]; // number of combined pairs, index of base pair
    	int    []   group_scale_ind = new int [ng]; // number of combined pairs, index of base pair
    	{
    		int ig = 0;
    		for (int i = 0; i < groups.length; i++) if (groups[i][0] >0){
    			groups_LMA[ig] = combineCompatiblePairs(
    					corrs,                       // double [][] correlations,
    					pair_mask,                   // int         pairs_mask,
    					(GROUPS[i][0] > 0),          // boolean     diagonal,
    					GROUPS[i][1]);               // int         baseline_scale
    			groups_pairs[ig] = groups[i]; // {number, base_pair_index}
    			group_scale_ind[ig] = scale_ind[i];
    			ig++;
    		}
    	}
    	double [] scales = new double [sl.size()];
    	for (int i = 0; i < scales.length; i++) scales[i] = sl.get(i); // from int to double
    	Correlations2dLMA lma = new Correlations2dLMA(scales);
    	if (debug_level > 1) {
    		for (int i = 0; i < groups_pairs.length; i++) {
    			System.out.println("Group #"+i+" - "+groups_pairs[i][0]+", type:"+groups_pairs[i][1]);
    		}
    	}
        addSamples(
        		xcenter,         // double            xcenter,   // preliminary center x in pixels for largest baseline
        		imgdtt_params.cnvx_hwnd_size, // int  hwindow_y, //  = window_y.length; // should actually be the same?
        		imgdtt_params.cnvx_hwnd_size, //int   hwindow_x, // = window_x.length;
        		vasw_pwr,                     // double            vasw_pwr,  // value as weight to this power,
        		groups_LMA,                   // double [][]       groups_LMA,
        		groups_pairs,                 // int    [][]       groups_pairs,
        		scales,                       // double []         scales,
        		group_scale_ind,              // int    []         group_scale_ind,
        		lma,                          // Correlations2dLMA lma,
        		imgdtt_params.cnvx_add3x3,    // boolean   add3x3,
        		imgdtt_params.cnvx_weight,    // double    nc_cost,
        		debug_level);    // int               debug_level
        boolean lmaSuccess;
        if (run_poly_instead) { // not used in lwir
        	lma.getMaxXYPoly( // get interpolated maximum coordinates using 2-nd degree polynomial
        			debug_level>3); // boolean debug
        	lmaSuccess = lma.getPolyFx() != null;
        } else {
        	lma.initVector(
        			imgdtt_params.lma_adjust_wm,   // boolean adjust_wm,
        			imgdtt_params.lma_adjust_wy,   // boolean adjust_wy,
        			imgdtt_params.lma_adjust_wxy,  // boolean adjust_wxy,
        			imgdtt_params.lma_adjust_ag,   // boolean adjust_Ag,
        			xcenter,                       // double  x0,
        			imgdtt_params.lma_half_width,  // double  half_width,
        			imgdtt_params.lma_cost_wy,     // double  cost_wy,     // cost of non-zero this.all_pars[WYD_INDEX]
        			imgdtt_params.lma_cost_wxy     //double  cost_wxy     // cost of non-zero this.all_pars[WXY_INDEX]
        			);
        	if (debug_level > 1) {
        		System.out.println("Input data:");
        		lma.printInputDataFx(false);
        	}

        	lmaSuccess = 	lma.runLma(
        			imgdtt_params.lma_lambda_initial,     // double lambda,           // 0.1
        			imgdtt_params.lma_lambda_scale_good,  // double lambda_scale_good,// 0.5
        			imgdtt_params.lma_lambda_scale_bad,   // double lambda_scale_bad, // 8.0
        			imgdtt_params.lma_lambda_max,         // double lambda_max,       // 100
        			imgdtt_params.lma_rms_diff,           // double rms_diff,         // 0.001
        			imgdtt_params.lma_num_iter,           // int    num_iter,         // 20
        			debug_level);       // int    debug_level)

        	lma.updateFromVector();
        	double [] rms = lma.getRMS();
        	if (debug_level > 0) {
        		System.out.println("LMA ->"+lmaSuccess+" RMS="+rms[0]+", pure RMS="+rms[1]);
        		lma.printParams();
        	}
        }
    	if ((debug_level > 1) && (groups_LMA !=null) && (groups_LMA.length > 0)) {
    		double [][] y_and_fx = new double [groups_LMA.length * 2][];
    		double [][] groups_fx =  getFitSamples( // just for debug to compare LMA-fitted fx with original data
    				xcenter,                      // double            xcenter,   // preliminary center x in pixels for largest baseline
            		imgdtt_params.cnvx_hwnd_size, // int  hwindow_y, //  = window_y.length; // should actually be the same?
            		imgdtt_params.cnvx_hwnd_size, //int   hwindow_x, // = window_x.length;
    				groups_pairs,                 // int    [][]       groups_pairs,
    				scales,                       // double []         scales,
    				group_scale_ind,              // int    []         group_scale_ind,
    				lma,                          // Correlations2dLMA lma,
            		groups_LMA,                   // double [][]       groups_LMA,
            		imgdtt_params.cnvx_add3x3,    // boolean   add3x3,
            		imgdtt_params.cnvx_weight,    // double    nc_cost,
            		debug_level);    // int               debug_level

    		String [] titles = new String [groups_LMA.length * 2];
    		for (int i = 0; i < groups_LMA.length; i++) if (groups_pairs[i][0] > 0){
    			int base_pair = groups_pairs[i][1];
    			titles[2 * i] =     (isDiagonalPair(base_pair)?"diag":"ortho")+getScaleOfPair(base_pair);
    			titles[2 * i + 1] = (isDiagonalPair(base_pair)?"diag":"ortho")+getScaleOfPair(base_pair)+"-fit";
    			y_and_fx[2 * i] =     groups_LMA[i];
    			y_and_fx[2 * i + 1] = groups_fx[i];
    		}

    		//    		String [] titles = {"ortho","diagonal"};
    		(new ShowDoubleFloatArrays()).showArrays(
    				y_and_fx,
    				2 * transform_size-1,
    				2 * transform_size-1,
    				true, (run_poly_instead?("mismatch"+pair_mask):"groups")+"_x"+tileX+"_y"+tileY, titles);
    	}
    	if (debug_level > 1) {
    		System.out.println("Input data and approximation:");
    		lma.printInputDataFx(true);
    	}
    	return lmaSuccess? lma: null;
    }




    public Correlations2dLMA corrLMA( // New for non-square
    		ImageDttParameters  imgdtt_params,
    		double [][]         disparity_distortions, // {d_disp/dx, d_ndisp/dx, d_disp/dy, d_ndisp/dy} for each camera
    		double [][]         corrs,
    		boolean[]           pair_mask, // which pairs to process
    		boolean             run_poly_instead, // true - run LMA, false - run 2d polynomial approximation
    		double              xcenter,   // preliminary center x in pixels for largest baseline
    		double              vasw_pwr,  // value as weight to this power,
    		int                 debug_level,
    		int                 tileX, // just for debug output
    		int                 tileY
    		)
    {
    	// for quad camera
    	int [][] groups = new int [GROUPS.length][];
    	int   [] scale_ind = new int  [GROUPS.length];
    	// See which groups exist for current pairs mask
    	int ng = 0;
    	ArrayList<Integer> sl = new ArrayList<Integer>();
    	for (int i = 0; i < GROUPS.length; i++) {
    		groups[i] = getNumberBaseOfCompatiblePairs(
    				corrs,                       // double [][] correlations,
    				pair_mask,                   // int         pairs_mask,
    				(GROUPS[i][0] > 0),          // boolean     diagonal,
    				GROUPS[i][1]);               // int         baseline_scale
    		if (groups[i][0] > 0) {
    			ng++;
    			if (!sl.contains(GROUPS[i][1])) {
    				sl.add(GROUPS[i][1]);
    			}
    			scale_ind[i] = sl.indexOf(GROUPS[i][1]);
    		}
    	}
    	if (debug_level > 1) {
    		System.out.println("corrLMA(): found "+ng+" groups, "+sl.size()+" scales");
    	}
    	double [][] groups_LMA =      new double [ng][];
    	int    [][] groups_pairs =    new int [ng][]; // [ng]{number of combined pairs, index of base pair}
    	int    []   group_scale_ind = new int [ng]; // number of combined pairs, index of base pair
    	{
    		int ig = 0;
    		for (int i = 0; i < groups.length; i++) if (groups[i][0] >0){
    			groups_LMA[ig] = combineCompatiblePairs(
    					corrs,                       // double [][] correlations,
    					pair_mask,                   // int         pairs_mask,
    					(GROUPS[i][0] > 0),          // boolean     diagonal,
    					GROUPS[i][1]);               // int         baseline_scale
    			groups_pairs[ig] = groups[i]; // {number, base_pair_index}
    			group_scale_ind[ig] = scale_ind[i];
    			ig++;
    		}
    	}
    	double [] scales = new double [sl.size()];
    	for (int i = 0; i < scales.length; i++) scales[i] = sl.get(i); // from int to double


    	Correlations2dLMA lma = new Correlations2dLMA(scales);
    	if (debug_level > 1) {
    		for (int i = 0; i < groups_pairs.length; i++) {
    			System.out.println("Group #"+i+" - "+groups_pairs[i][0]+", type:"+groups_pairs[i][1]);
    		}
    	}
        addSamples(
        		xcenter,         // double            xcenter,   // preliminary center x in pixels for largest baseline
        		imgdtt_params.cnvx_hwnd_size, // int  hwindow_y, //  = window_y.length; // should actually be the same?
        		imgdtt_params.cnvx_hwnd_size, //int   hwindow_x, // = window_x.length;
        		vasw_pwr,                     // double            vasw_pwr,  // value as weight to this power,
        		groups_LMA,                   // double [][]       groups_LMA,
        		groups_pairs,                 // int    [][]       groups_pairs,
        		scales,                       // double []         scales,
        		group_scale_ind,              // int    []         group_scale_ind,
        		lma,                          // Correlations2dLMA lma,
        		imgdtt_params.cnvx_add3x3,    // boolean   add3x3,
        		imgdtt_params.cnvx_weight,    // double    nc_cost,

        		debug_level);    // int               debug_level
        boolean lmaSuccess;
        if (run_poly_instead) { // not used in lwir
        	lma.getMaxXYPoly( // get interpolated maximum coordinates using 2-nd degree polynomial
        			debug_level>3); // boolean debug
        	lmaSuccess = lma.getPolyFx() != null;
        } else {
        	lma.initVector(
        			imgdtt_params.lma_adjust_wm,   //  boolean adjust_wm,
        			imgdtt_params.lma_adjust_wy,   // boolean adjust_wy,
        			imgdtt_params.lma_adjust_wxy,  // boolean adjust_wxy,
        			imgdtt_params.lma_adjust_ag,   // boolean adjust_Ag,
        			xcenter,                       // double  x0,
        			imgdtt_params.lma_half_width,  // double  half_width,
        			imgdtt_params.lma_cost_wy,     // double  cost_wy,     // cost of non-zero this.all_pars[WYD_INDEX]
        			imgdtt_params.lma_cost_wxy     //double  cost_wxy     // cost of non-zero this.all_pars[WXY_INDEX]
        			);
        	if (debug_level > 1) {
        		System.out.println("Input data:");
        		lma.printInputDataFx(false);
        	}

        	lmaSuccess = 	lma.runLma(
        			imgdtt_params.lma_lambda_initial,     // double lambda,           // 0.1
        			imgdtt_params.lma_lambda_scale_good,  // double lambda_scale_good,// 0.5
        			imgdtt_params.lma_lambda_scale_bad,   // double lambda_scale_bad, // 8.0
        			imgdtt_params.lma_lambda_max,         // double lambda_max,       // 100
        			imgdtt_params.lma_rms_diff,           // double rms_diff,         // 0.001
        			imgdtt_params.lma_num_iter,           // int    num_iter,         // 20
        			debug_level);       // int    debug_level)

        	lma.updateFromVector();
        	double [] rms = lma.getRMS();
        	if (debug_level > 0) {
        		System.out.println("LMA ->"+lmaSuccess+" RMS="+rms[0]+", pure RMS="+rms[1]);
        		lma.printParams();
        	}
        }
    	if ((debug_level > 1) && (groups_LMA !=null) && (groups_LMA.length > 0)) {
    		double [][] y_and_fx = new double [groups_LMA.length * 2][];
    		double [][] groups_fx =  getFitSamples( // just for debug to compare LMA-fitted fx with original data
    				xcenter,                      // double            xcenter,   // preliminary center x in pixels for largest baseline
            		imgdtt_params.cnvx_hwnd_size, // int  hwindow_y, //  = window_y.length; // should actually be the same?
            		imgdtt_params.cnvx_hwnd_size, //int   hwindow_x, // = window_x.length;
    				groups_pairs,                 // int    [][]       groups_pairs,
    				scales,                       // double []         scales,
    				group_scale_ind,              // int    []         group_scale_ind,
    				lma,                          // Correlations2dLMA lma,
            		groups_LMA,                   // double [][]       groups_LMA,
            		imgdtt_params.cnvx_add3x3,    // boolean   add3x3,
            		imgdtt_params.cnvx_weight,    // double    nc_cost,
            		debug_level);    // int               debug_level

    		String [] titles = new String [groups_LMA.length * 2];
    		for (int i = 0; i < groups_LMA.length; i++) if (groups_pairs[i][0] > 0){
    			int base_pair = groups_pairs[i][1];
    			titles[2 * i] =     (isDiagonalPair(base_pair)?"diag":"ortho")+getScaleOfPair(base_pair);
    			titles[2 * i + 1] = (isDiagonalPair(base_pair)?"diag":"ortho")+getScaleOfPair(base_pair)+"-fit";
    			y_and_fx[2 * i] =     groups_LMA[i];
    			y_and_fx[2 * i + 1] = groups_fx[i];
    		}

    		//    		String [] titles = {"ortho","diagonal"};
    		(new ShowDoubleFloatArrays()).showArrays(
    				y_and_fx,
    				2 * transform_size-1,
    				2 * transform_size-1,
    				true, (run_poly_instead?("mismatch"+pair_mask):"groups")+"_x"+tileX+"_y"+tileY, titles);
    	}
    	if (debug_level > 1) {
    		System.out.println("Input data and approximation:");
    		lma.printInputDataFx(true);
    	}
    	return lmaSuccess? lma: null;
    }

// Run for a single horizontal 2d correlation array
    public Correlations2dLMA corrLMA( // not used in lwir
    		ImageDttParameters  imgdtt_params,
    		double []           corr,
    		boolean             run_poly_instead, // true - run LMA, false - run 2d polynomial approximation
    		double              xcenter,   // preliminary center x in pixels for largest baseline
    		double              vasw_pwr,  // value as weight to this power,
    		int                 debug_level,
    		int                 tileX, // just for debug output
    		int                 tileY)
    {
    	double [][] groups_LMA =       {corr};     // new double [ng][];
    	int    [][] groups_pairs =     {{1,0}};     // new int[ng][]; // number of combined pairs, index of base pair
    	int    []   group_scale_ind =  {0}; // new int [ng]; // number of combined pairs, index of base pair

    	double [] scales = {1.0};

    	Correlations2dLMA lma = new Correlations2dLMA(scales);
    	if (debug_level > 1) {
    		for (int i = 0; i < groups_pairs.length; i++) {
    			System.out.println("Group #"+i+" - "+groups_pairs[i][0]+", type:"+groups_pairs[i][1]);
    		}
    	}
        addSamples(
        		xcenter,         // double            xcenter,   // preliminary center x in pixels for largest baseline
        		imgdtt_params.cnvx_hwnd_size, // int  hwindow_y, //  = window_y.length; // should actually be the same?
        		imgdtt_params.cnvx_hwnd_size, //int   hwindow_x, // = window_x.length;
        		vasw_pwr,                     // double            vasw_pwr,  // value as weight to this power,
        		groups_LMA,                   // double [][]       groups_LMA,
        		groups_pairs,                 // int    [][]       groups_pairs,
        		scales,                       // double []         scales,
        		group_scale_ind,              // int    []         group_scale_ind,
        		lma,                          // Correlations2dLMA lma,
        		imgdtt_params.cnvx_add3x3,    // boolean   add3x3,
        		imgdtt_params.cnvx_weight,    // double    nc_cost,

        		debug_level);    // int               debug_level
        boolean lmaSuccess;
        if (run_poly_instead) {
        	lma.getMaxXYPoly( // get interpolated maximum coordinates using 2-nd degree polynomial
        			debug_level>3); // boolean debug
        	lmaSuccess = lma.getPolyFx() != null;
        } else {
        	lma.initVector(
        			imgdtt_params.lma_adjust_wm,   //  boolean adjust_wm,
        			imgdtt_params.lma_adjust_wy,   // boolean adjust_wy,
        			imgdtt_params.lma_adjust_wxy,  // boolean adjust_wxy,
        			imgdtt_params.lma_adjust_ag,   // boolean adjust_Ag,
        			xcenter,                       // double  x0,
        			imgdtt_params.lma_half_width,  // double  half_width,
        			imgdtt_params.lma_cost_wy,     // double  cost_wy,     // cost of non-zero this.all_pars[WYD_INDEX]
        			imgdtt_params.lma_cost_wxy     //double  cost_wxy     // cost of non-zero this.all_pars[WXY_INDEX]
        			);
        	if (debug_level > 1) {
        		System.out.println("Input data:");
        		lma.printInputDataFx(false);
        	}

        	lmaSuccess = 	lma.runLma(
        			imgdtt_params.lma_lambda_initial,     // double lambda,           // 0.1
        			imgdtt_params.lma_lambda_scale_good,  // double lambda_scale_good,// 0.5
        			imgdtt_params.lma_lambda_scale_bad,   // double lambda_scale_bad, // 8.0
        			imgdtt_params.lma_lambda_max,         // double lambda_max,       // 100
        			imgdtt_params.lma_rms_diff,           // double rms_diff,         // 0.001
        			imgdtt_params.lma_num_iter,           // int    num_iter,         // 20
        			debug_level);       // int    debug_level)

        	lma.updateFromVector();
        	double [] rms = lma.getRMS();
        	if (debug_level > 0) {
        		System.out.println("LMA ->"+lmaSuccess+" RMS="+rms[0]+", pure RMS="+rms[1]);
        		lma.printParams();
        	}
        }
    	if ((debug_level > 1) && (groups_LMA !=null) && (groups_LMA.length > 0)) {
    		double [][] y_and_fx = new double [groups_LMA.length * 2][];
    		double [][] groups_fx =  getFitSamples( // just for debug to compare LMA-fitted fx with original data
    				xcenter,                      // double            xcenter,   // preliminary center x in pixels for largest baseline
            		imgdtt_params.cnvx_hwnd_size, // int  hwindow_y, //  = window_y.length; // should actually be the same?
            		imgdtt_params.cnvx_hwnd_size, //int   hwindow_x, // = window_x.length;
    				groups_pairs,                 // int    [][]       groups_pairs,
    				scales,                       // double []         scales,
    				group_scale_ind,              // int    []         group_scale_ind,
    				lma,                          // Correlations2dLMA lma,
            		groups_LMA,                   // double [][]       groups_LMA,
            		imgdtt_params.cnvx_add3x3,    // boolean   add3x3,
            		imgdtt_params.cnvx_weight,    // double    nc_cost,
            		debug_level);    // int               debug_level

    		String [] titles = new String [groups_LMA.length * 2];
    		for (int i = 0; i < groups_LMA.length; i++) if (groups_pairs[i][0] > 0){
    			int base_pair = groups_pairs[i][1];
    			titles[2 * i] =     (isDiagonalPair(base_pair)?"diag":"ortho")+getScaleOfPair(base_pair);
    			titles[2 * i + 1] = (isDiagonalPair(base_pair)?"diag":"ortho")+getScaleOfPair(base_pair)+"-fit";
    			y_and_fx[2 * i] =     groups_LMA[i];
    			y_and_fx[2 * i + 1] = groups_fx[i];
    		}

    		//    		String [] titles = {"ortho","diagonal"};
    		(new ShowDoubleFloatArrays()).showArrays(
    				y_and_fx,
    				2 * transform_size-1,
    				2 * transform_size-1,
    				true, (run_poly_instead?"poly":"lma")+"_x"+tileX+"_y"+tileY, titles);
    	}
    	if (debug_level > 1) {
    		System.out.println("Input data and approximation:");
    		lma.printInputDataFx(true);
    	}
    	return lmaSuccess? lma: null;
    }

    /**
     * Verify that 2D correlation pair has enough enabled samples and they are not all collinear (if requested)
     * @param weights 2D correlation weights (selected if > 0) in line-scan order
     * @param size of the square (15 for transform size 8)
     * @param min_samples minimal required number of non-zero samples
     * @param non_coll requre non-collinear samples
     * @return true if the data matches requirements
     */
    public boolean checkCluster(
    		double [] weights,
    		int       size,
    		int min_samples,
    		boolean non_coll
    		) {
    	int num_samples = 0;
    	int [][] first2 = new int [2][2];
    	boolean noncolinOK = !non_coll;
    	for (int iy = 0; iy < size; iy++) {
    		for (int ix = 0; ix < size; ix++) {
    			if (weights[iy * size + ix] > 0) {
    				if (num_samples == 0) {
    					first2[0][0] = ix; // x0
    					first2[0][1] = iy; // y0
    				} else if (num_samples == 1) {
    					first2[1][0] = ix - first2[0][0]; // dx1
    					first2[1][1] = iy - first2[0][1]; // dy1
    				} else {
    					noncolinOK |= ((ix-first2[0][0])*first2[1][1] != (iy-first2[0][1])*first2[1][0]);
    				}
    				num_samples++;
    				if (noncolinOK && (num_samples >= min_samples)) {
    					return true;
    				}
    			}
    		}
    	}
    	return false;
    }

    /**
     * Create mask of usable points, allowing only first non bi-convex away from the center
     * @param corr_data correlation data, packed in linescan order
     * @param hwin half-window width, in pixels
     * @param x0c start point x (-disparity)
     * @param y0c start point y (for ortho == 0, for diagonal ==x0
     * @param add3x3 always enable 3x3 points around start point
     * @param nc_cost - cost of non-convex points
     * @param debug
     * @return cost packed array, corresponding to the input. selected convex points have weight
     * 1.0, other selected - nc_cost
     */

    public double [] filterConvex(// USED in lwir
    		double [] corr_data,
    		int       hwin,
    		int       x0c,
    		int       y0c,
    		boolean   add3x3,
    		double    nc_cost,
    		boolean   debug) {

    	int center =       transform_size - 1;
    	int x0 =           x0c + center;
    	int y0=            y0c +center;
    	int width =        2 * center + 1;
    	int dlen =         width * width;
    	double [] weights = new double [dlen];
    	boolean [] convex = new boolean[dlen];
    	int min_row = y0 - hwin;
    	if (min_row < 0) min_row = 0;
    	int max_row = y0 + hwin;
    	if (max_row >= width) max_row = width -1 ;
    	int min_col = x0 - hwin;
    	if (min_col < 0) min_col = 0;
    	int max_col = x0 + hwin;
    	if (max_col >= width) max_col = width -1 ;

    	for (int row = min_row + 1; row < max_row; row ++) {
        	for (int col = min_col + 1; col < max_col; col ++) {
        		int indx = row * width + col;
        		convex[indx] = (2 * corr_data[indx] > (corr_data[indx - 1] + corr_data[indx + 1])) &&
        				( 2 * corr_data[indx] > (corr_data[indx - width] + corr_data[indx + width]));
        	}
    	}
    	boolean [] sel = new boolean [dlen];
    	if ((y0 >= min_row) && (y0 <= max_row) && (x0 >= min_col) && (x0 <= max_col)) {
    		sel [y0 * width + x0] = true; // start point
    	}

        if (debug) debug_convex(convex, sel,"initial");


    	for (int istep = 0; istep < hwin; istep ++) {
    		boolean [] sel_new = sel.clone();
    		boolean first = (istep == 0); // consider initial cell as if convex
    		int num_new = 0;
        	int step_min_row = y0 - istep;
        	if (step_min_row < 0) step_min_row = 0;
        	int step_max_row = y0 + istep;
        	if (step_max_row >= width) step_max_row = width -1 ;
        	int step_min_col = x0 - istep;
        	if (step_min_col < 0) step_min_col = 0;
        	int step_max_col = x0 + istep;
        	if (step_max_col >= width) step_max_col = width -1 ;
    		// expand up
    		if (istep < y0) {
    			int row = y0 - istep;
    			for (int col = step_min_col; col <= step_max_col; col++){
    				int indx = row * width + col;
    				if (sel[indx] && (convex[indx] || first)) {
    					if (!sel_new[indx - width]) {
    						sel_new[indx - width] = true;
    						num_new++;
    					}
    					// try next/prev
    					indx -= width;
    					if (convex[indx]) {
    						if ((col <= x0) && (col > min_col) && !sel_new[indx - 1]) {
    	    					sel_new[indx - 1] = true;
    	    					num_new++;
    						}
    						if ((col >= x0) && (col < max_col) && !sel_new[indx + 1]) {
    	    					sel_new[indx + 1] = true;
    	    					num_new++;
    						}
    					}
    				}
    			}
    		}
            if (debug) debug_convex(convex, sel_new,"expanded up, step="+istep+",  num_new="+num_new);
    		// expand down
    		if ((y0 + istep) < max_row) {
    			int row = y0 + istep;
    			for (int col = step_min_col; col <= step_max_col; col++){
    				int indx = row * width + col;
    				if (sel[indx] && (convex[indx] || first)) {
    					if (!sel_new[indx + width]) {
    						sel_new[indx + width] = true;
    						num_new++;
    					}
    					// try next/prev
    					indx += width;
    					if (convex[indx]) {
    						if ((col <= x0) && (col > min_col) && !sel_new[indx - 1]) {
    	    					sel_new[indx - 1] = true;
    	    					num_new++;
    						}
    						if ((col >= x0) && (col < max_col) && !sel_new[indx + 1]) {
    	    					sel_new[indx + 1] = true;
    	    					num_new++;
    						}
    					}
    				}
    			}
    		}
            if (debug) debug_convex(convex, sel_new,"expanded down, step="+istep+", num_new="+num_new);
    		// expand left
    		if (istep < x0) {
    			int col = x0 - istep;
    			for (int row = step_min_row; row <= step_max_row; row++){
    				int indx = row * width + col;
    				if (sel[indx] && (convex[indx] || first)) {
    					if (!sel_new[indx - 1]) {
    						sel_new[indx - 1] = true;
    						num_new++;
    					}
    					// try next/prev
    					indx -= 1;
    					if (convex[indx]) {
    						if ((row <= y0) && (row > min_row) && !sel_new[indx - width]) {
    	    					sel_new[indx - width] = true;
    	    					num_new++;
    						}
    						if ((row >= y0) && (row < max_row) && !sel_new[indx + width]) {
    	    					sel_new[indx + width] = true;
    	    					num_new++;
    						}
    					}
    				}
    			}
    		}
            if (debug) debug_convex(convex, sel_new,"expanded left, step="+istep+", num_new="+num_new);
    		// expand right
    		if ((x0 + istep) < max_col) {
    			int col = x0 + istep;
    			for (int row = step_min_row; row <= step_max_row; row++){
    				int indx = row * width + col;
    				if (sel[indx] && (convex[indx] || first)) {
    					if (!sel_new[indx + 1]) {
    						sel_new[indx + 1] = true;
    						num_new++;
    					}
    					// try next/prev
    					indx += 1;
    					if (convex[indx]) {
    						if ((row <= y0) && (row > min_row) && !sel_new[indx - width]) {
    	    					sel_new[indx - width] = true;
    	    					num_new++;
    						}
    						if ((row >= y0) && (row < max_row) && !sel_new[indx + width]) {
    	    					sel_new[indx + width] = true;
    	    					num_new++;
    						}
    					}
    				}
    			}
    		}
            if (debug) debug_convex(convex, sel_new,"expanded right, step="+istep+", num_new="+num_new);
    		if (num_new == 0) break;
    		sel = sel_new; // no need to clone
    	}
    	for (int i = 0; i < sel.length;i++) if (sel[i]){
    		weights[i] = convex[i]? 1.0: nc_cost;
    	}

    	// add central square
    	if (add3x3) {
    		for (int row = y0 - 1; row <= y0+1; row++)     if ((row >= 0) && (row < width)) {
        		for (int col = x0 - 1; col <= x0+1; col++) if ((col >= 0) && (col < width)) {
        			int indx = row * width + col;
        			if (weights[indx] == 0.0) weights[indx] = nc_cost; // non-connected convex should not matter
//        			sel[row * width + col] = true;
        		}
    		}
    	}
    	return weights;
    }

    public void debug_convex( // not used in lwir
    		boolean [] convex,
    		boolean [] sel,
    		String title) {
    	int center =       transform_size - 1;
    	int width =        2 * center + 1;
    	System.out.println(title);
		for (int row = 0; row < width; row++) {
			System.out.print(String.format("%3d: ", row));
    		for (int col = 0; col < width; col++){
    			int indx = row * width + col;
    			String s = sel[indx]? "*": ".";
    			if (convex[indx]) System.out.print(String.format("(%1s)", s));
    			else              System.out.print(String.format(" %1s ", s));
    		}
    		System.out.println();
		}

    }

    public void addSamples(// USED in lwir
    		double            xcenter,   // preliminary center x in pixels for largest baseline
        	int               hwindow_y, //  = window_y.length; // should actually be the same?
        	int               hwindow_x, // = window_x.length;
    		double            vasw_pwr,  // value as weight to this power,
    		double [][]       groups_LMA,
    		int    [][]       groups_pairs,
    		double []         scales,
    		int    []         group_scale_ind,
    		Correlations2dLMA lma,
    		boolean           add3x3,
    		double            nc_cost,
    		int               debug_level
    		) {
    	int center =       transform_size - 1;
    	int width =        2 * center + 1;
    	int center_index = (width + 1) * center; // in
    	// convex filter expects half window in pixels, arow/acol - half-pixel grid
    	int hwindow_y2 = 2 * hwindow_y + 1;
    	int hwindow_x2 = 2 * hwindow_x + 1;
    	int [][] quad_signs = {{-1,-1},{1,-1},{-1,1},{1,1}}; // {sign_x, sign_y} per quadrant
       	for (int ig = 0; ig < groups_pairs.length; ig++) if (groups_pairs[ig][0] > 0) {
    		double scale = scales[group_scale_ind[ig]];
    		double scale05 = scale/2.0;
    		boolean diagonal = isDiagonalPair(groups_pairs[ig][1]);
    		int ixcenter = (int) Math.round(xcenter / scale);
    		double xcs = ixcenter*scale;

    		if (debug_level > 0) {
    			System.out.println("\nCombined correlation data, diagonal = "+diagonal);
    			for (int row = 0; row < width; row++) {
        			System.out.print(String.format("%3d: ", row));
        			for (int col = 0; col < width; col++) {
        				if ((row == center) && (col == center)) {
                			System.out.print(String.format("[%8.5f]", groups_LMA[ig][width * row + col]));
        				} else {
        					System.out.print(String.format(" %8.5f ", groups_LMA[ig][width * row + col]));
        				}
        			}
        			System.out.println();
    			}
    			System.out.println("\nAll convex");
    			for (int row = 0; row < width; row++) {
        			System.out.print(String.format("%3d: ", row));
        			for (int col = 0; col < width; col++) {
        				if ((col==0) || (row==0) || (row == (width -1)) || (col == (width -1))) {
    						System.out.print(" - ");
        				} else {
        					boolean convex_x =  2 * groups_LMA[ig][width * row + col] > (groups_LMA[ig][width * row + col -1] +      groups_LMA[ig][width * row + col+1]);
        					boolean convex_y =  2 * groups_LMA[ig][width * row + col] > (groups_LMA[ig][width * row + col - width] + groups_LMA[ig][width * row + col + width]);
        					if ((row == center) && (col == center)) {
        						System.out.print(String.format("[%1s]", (convex_x && convex_y)?"*":"."));
        					} else {
        						System.out.print(String.format(" %1s ", (convex_x && convex_y)?"*":"."));
        					}
        				}
        			}
        			System.out.println();
    			}
    			// test filter convex
    		}
			double [] filtWeight =  filterConvex(
					groups_LMA[ig],         // double [] corr_data,
					hwindow_x,              // int       hwin,
					ixcenter,               // int       x0,
					(diagonal? ixcenter:0), // int       y0,
					add3x3,                 // boolean   add3x3,
					nc_cost,                // double    nc_cost,
		    		(debug_level > 2));     // boolean   debug);
			if (debug_level > 1) {
				System.out.println("\nConvex weights");
				for (int row = 0; row < width; row++) {
					System.out.print(String.format("%3d: ", row));
					for (int col = 0; col < width; col++) {
						double d = filtWeight[width * row + col];
						if ((row == center) && (col == center)) {
							if (d == 0.0) {
								System.out.print(String.format("[%5s]","."));
							} else {
								System.out.print(String.format("[%5.2f]",d));
							}
						} else {
							if (d == 0.0) {
								System.out.print(String.format(" %5s ","."));
							} else {
								System.out.print(String.format(" %5.2f ", d));
							}
						}
					}
					System.out.println();
				}
			}
			lma.setDiag(diagonal);
    		if (diagonal) {
    			for (int arow =  0; arow < hwindow_y2; arow ++) {
    				int odd = arow & 1;
    				for (int acol =  odd; acol < hwindow_x2; acol +=2) {
    					for (int quad = 0; quad < 4; quad ++) if (((arow > 0) || ((quad & 2) !=0 )) && ((acol > 0) || ((quad & 1) !=0 ))){
    						int cx = (quad_signs[quad][0] * acol - quad_signs[quad][1] * arow)/2 + ixcenter; // ix0;
    						int cy = (quad_signs[quad][0] * acol + quad_signs[quad][1] * arow)/2 + ixcenter; // ix0;
    						// calculate coordinates in the correlation array
    						if ((cx >= -center) && (cx <= center) && (cy >= -center) && (cy <= center)) {
    							double w = filtWeight[center_index + width * cy + cx] * groups_pairs[ig][0]; // number of pair averaged
    							if (w > 0.0) {
    								double v = groups_LMA[ig][center_index + width * cy + cx];
    								if (vasw_pwr != 0) {
    									w *= Math.pow(Math.abs(v), vasw_pwr);
    								}
    								lma.addSample(
    										quad_signs[quad][0] * acol * scale05 + xcs, // x * scale, // double x,      // x coordinate on the common scale (corresponding to the largest baseline), along the disparity axis
    										quad_signs[quad][1] * arow * scale05,       // y * scale, // double y,      // y coordinate (0 - disparity axis)
    										v,                   // double v,      // correlation value at that point
    										w,                 // double w,
    										group_scale_ind[ig], // int    si,     // baseline scale index
    										ig);                 // int    gi);
    							}
    						}
    					}
    				}
    			}
    		} else { // ortho
    			for (int arow =  0; arow < hwindow_y2; arow += 2) {
    				for (int acol =  0; acol < hwindow_x2; acol +=2) {
    					for (int quad = 0; quad < 4; quad ++) if (((arow > 0) || ((quad & 2) !=0 )) && ((acol > 0) || ((quad & 1) !=0 ))){
    						int cx = (quad_signs[quad][0] * acol)/2 + ixcenter; // ix0;
    						int cy = (quad_signs[quad][1] * arow)/2;
    						// calculate coordinates in the correlation array
    						if ((cx >= -center) && (cx <= center) && (cy >= -center) && (cy <= center)) {
    							double w = filtWeight[center_index + width * cy + cx] * groups_pairs[ig][0]; // number of pair averaged;
    							if (w > 0.0) {
    								double v = groups_LMA[ig][center_index + width * cy + cx];
    								if (vasw_pwr != 0) {
    									w *= Math.pow(Math.abs(v), vasw_pwr);
    								}
    								lma.addSample(
    										quad_signs[quad][0] * acol * scale05 + xcs, // x * scale, // double x,      // x coordinate on the common scale (corresponding to the largest baseline), along the disparity axis
    										quad_signs[quad][1] * arow * scale05,       // y * scale, // double y,      // y coordinate (0 - disparity axis)
    										v,                   // double v,      // correlation value at that point
    										w,                 // double w,
    										group_scale_ind[ig], // int    si,     // baseline scale index
    										ig);                 // int    gi);
    							}
    						}
    					}
    				}
    			}

    		}
    	}
    }

    // Mimics addSamples, but reads f(x) values instead of setting them
    public double [][] getFitSamples( // just for debug to compare LMA-fitted fx with original data // not used in lwir
    		double            xcenter,   // preliminary center x in pixels for largest baseline
        	int hwindow_y, // should actually be the same?
        	int hwindow_x,
    		int    [][]       groups_pairs,
    		double []         scales,
    		int    []         group_scale_ind,
    		Correlations2dLMA lma,
    		double [][]       groups_LMA,
    		boolean           add3x3,
    		double            nc_cost,
    		int               debug_level) {
    	double [] fx = lma.getPolyFx();  // try if poly results are available, if not use LMA results
    	if (fx == null) fx=lma.getFx(); // fx will be 2 samples longer because of the regularization terms.
    	int center =       transform_size - 1;
    	int width =        2 * center + 1;
    	double [][] groups_fitted = new double [groups_pairs.length][width * width];

    	int center_index = (width + 1) * center; // in
    	// convex filter expects half window in pixels, arow/acol - half-pixel grid
    	int hwindow_y2 = 2 * hwindow_y + 1;
    	int hwindow_x2 = 2 * hwindow_x + 1;
    	int [][] quad_signs = {{-1,-1},{1,-1},{-1,1},{1,1}}; // {sign_x, sign_y} per quadrant
    	int numSample = 0;
       	for (int ig = 0; ig < groups_pairs.length; ig++) if (groups_pairs[ig][0] > 0) {
       		for (int i = 0; i < groups_fitted[ig].length; i++) {
       			groups_fitted[ig][i] = Double.NaN;
       		}
    		double scale = scales[group_scale_ind[ig]];
    		boolean diagonal = isDiagonalPair(groups_pairs[ig][1]);
    		int ixcenter = (int) Math.round(xcenter / scale);
			double [] filtWeight =  filterConvex(
					groups_LMA[ig],       // double [] corr_data,
					hwindow_x,            // int       hwin,
					ixcenter,             // int       x0,
					(diagonal? ixcenter:0), // int       y0,
					add3x3, // false, // true,                 // boolean   add3x3,
					nc_cost,                  // double    nc_cost,
		    		false);// boolean   debug);

    		if (diagonal) {
    			for (int arow =  0; arow < hwindow_y2; arow ++) {
    				int odd = arow & 1;
    				for (int acol =  odd; acol < hwindow_x2; acol +=2) {
    					for (int quad = 0; quad < 4; quad ++) if (((arow > 0) || ((quad & 2) !=0 )) && ((acol > 0) || ((quad & 1) !=0 ))){
    						int cx = (quad_signs[quad][0] * acol - quad_signs[quad][1] * arow)/2 + ixcenter; // ix0;
    						int cy = (quad_signs[quad][0] * acol + quad_signs[quad][1] * arow)/2 + ixcenter; // ix0;
    						// calculate coordinates in the correlation array
    						if ((cx >= -center) && (cx <= center) && (cy >= -center) && (cy <= center)) {
    							double w = filtWeight[center_index + width * cy + cx];
    							if (w > 0.0) {
    								double v = fx[numSample++];
    								if (v > 0.0) groups_fitted[ig][center_index + width * cy + cx] = v;
    								//    							double v = groups_LMA[ig][center_index + width * cy + cx];
    							}

    						}
    					}
    				}
    			}
    		} else { // ortho
    			for (int arow =  0; arow < hwindow_y2; arow += 2) {
    				for (int acol =  0; acol < hwindow_x2; acol +=2) {
    					for (int quad = 0; quad < 4; quad ++) if (((arow > 0) || ((quad & 2) !=0 )) && ((acol > 0) || ((quad & 1) !=0 ))){
    						int cx = (quad_signs[quad][0] * acol)/2 + ixcenter; // ix0;
    						int cy = (quad_signs[quad][1] * arow)/2;
    						// calculate coordinates in the correlation array
    						if ((cx >= -center) && (cx <= center) && (cy >= -center) && (cy <= center)) {
    							double w = filtWeight[center_index + width * cy + cx];
    							if (w > 0.0) {
    								double v = fx[numSample++];
    								if (v > 0.0) groups_fitted[ig][center_index + width * cy + cx] = v;
    								//    							double v = groups_LMA[ig][center_index + width * cy + cx];
    							}

    						}
    					}
    				}
    			}

    		}
    	}
       	return groups_fitted;
    }

    public int [] listPairs( // USED in lwir
    		double [][] correlations,
    		int         pairs_mask) {
		ArrayList<Integer> pairs = new ArrayList<Integer>();
		for (int np = 0; np < correlations.length; np++) if ((correlations[np] != null) && ((( 1 << np) & pairs_mask) != 0)){
			pairs.add(np);
		}
    	int [] rslt = new int[pairs.size()];
    	for (int i = 0; i < rslt.length; i++) {
    		rslt[i] = pairs.get(i);
    	}
    	return rslt;
    }

    /**
     * Reproducing original hor/vert feature detection using a notch filter. Will need to extend
     * to use diagonals too, as well as variable baselines similar to LMA/poly
     * Combines several (same-oriented) correlation arrays
     * Uses 3-point polynomial interpolation
     * @param correlations correlation data packed in linescan order, first index is for pairs
     * @param pairs_mask bitmask of pairs to use (should only include same-scale, same orientation)
     * enhortho_scales notch filter array is taken from this class global array
     * @param radius positive - within that distance, negative - within 2*(-radius)+1 square
     * @param debug
     * @return {center, strength} pair (center is 0 for the correlation center)
     */
	public double [] getMaxXSOrtho( // // get fractional center using a quadratic polynomial // USED in lwir
    		double [][] correlations,
    		int         pairs_mask,
			double      offset,      // double      offset);
			boolean     symmetric,   // for comparing with old implementation average with symmetrical before multiplication
			boolean     is_vert,      // transpose X/Y
			boolean     debug)
	{
		if (debug) {
			System.out.println("getMaxXSOrtho()");
		}
		int corr_size = transform_size * 2 -1;
		double [] corr_1d = new double [corr_size];
		int [] pairs = listPairs(correlations,pairs_mask);
		if (pairs.length==0) return null;
		double w = 1.0/pairs.length;

		if (debug){
			System.out.println ("Before combining tiles");
			for (int ip = 0; ip < pairs.length; ip++) {
				int np = pairs[ip];
				System.out.println("pair # "+np);
				for (int i = 0; i < corr_size; i++) {
					System.out.print(String.format("%2d:", i));
					for (int j = 0; j < corr_size; j++) {
						System.out.print(String.format(" %8.5f", correlations[np][i * corr_size + j]));
					}
					System.out.println();
				}
				System.out.println();
			}
		}



		if (debug){
			System.out.println ("combined tiles (transposed):");
		}

		if (is_vert) {
			if (offset >= 0) { // use shifted multiplication (similar to combineInterpolatedCorrelations(), but no need to use twice_diagonal
				for (int j = 0; j < corr_size; j++){
					if (debug) System.out.print(String.format("transposed %2d:", j));
					corr_1d[j] = 0;
					for (int i = 0; i < corr_size; i++){
						double combo = 1.0;
						for (int ip = 0; ip < pairs.length; ip++) {
							int np = pairs[ip];
							double d = (symmetric ? (0.5*(correlations[np][j * corr_size + i] +  correlations[np][j * corr_size + (corr_size - i - 1)])): (correlations[np][j * corr_size + i])) + offset;
							/// if (d < 0.0) d = 0.0; // It is how it was before, and for combining exactly 2 pairs it may be OK
							combo *= d;
							//combo += w * correlations[np][i * corr_size + j] * this.ortho_notch_filter[i];
						}
						if (combo > 0.0) {
							combo = Math.pow(combo, w) - offset;
						} else {
							combo =  - offset; // how it was before, just to compare
						}
						corr_1d[j] += combo * this.ortho_notch_filter[i];
						if (debug) System.out.print(String.format(" %8.5f", combo));
					}
					if (debug) System.out.println();
				}
			} else { // use averaging (linear) // not used in lwir
				for (int j = 0; j < corr_size; j++){
					corr_1d[j] = 0;
					for (int ip = 0; ip < pairs.length; ip++) {
						int np = pairs[ip];
						for (int i = 0; i < corr_size; i++){
							corr_1d[j] += w * correlations[np][j * corr_size + i] * this.ortho_notch_filter[i];
						}
					}
				}
			}

		} else {
			if (offset >= 0) { // use shifted multiplication (similar to combineInterpolatedCorrelations(), but no need to use twice_diagonal
				for (int j = 0; j < corr_size; j++){
					if (debug) System.out.print(String.format("transposed %2d:", j));
					corr_1d[j] = 0;
					for (int i = 0; i < corr_size; i++){
						double combo = 1.0;
						for (int ip = 0; ip < pairs.length; ip++) {
							int np = pairs[ip];
							double d = (symmetric ? (0.5*(correlations[np][i * corr_size + j] +  correlations[np][(corr_size - i - 1) * corr_size + j])): (correlations[np][i * corr_size + j])) + offset;
							/// if (d < 0.0) d = 0.0; // It is how it was before, and for combining exactly 2 pairs it may be OK
							combo *= d;
							//combo += w * correlations[np][i * corr_size + j] * this.ortho_notch_filter[i];
						}
						if (combo > 0.0) {
							combo = Math.pow(combo, w) - offset;
						} else {
							combo =  - offset; // how it was before, just to compare
						}
						corr_1d[j] += combo * this.ortho_notch_filter[i];
						if (debug) System.out.print(String.format(" %8.5f", combo));
					}
					if (debug) System.out.println();
				}
			} else { // use averaging (linear) // not used in lwir
				for (int j = 0; j < corr_size; j++){
					corr_1d[j] = 0;
					for (int ip = 0; ip < pairs.length; ip++) {
						int np = pairs[ip];
						for (int i = 0; i < corr_size; i++){
							corr_1d[j] += w * correlations[np][i * corr_size + j] * this.ortho_notch_filter[i];
						}
					}
				}
			}
		}




		if (debug) {
			System.out.println();
			System.out.print ("corr_1d = ");
			for (int j = 0; j < corr_size; j++){
				if (debug) System.out.print(String.format(" %8.5f", corr_1d[j]));
			}
			System.out.println();
		}

		int icenter = 0;
		for (int i = 1; i < corr_size; i++){
			if (corr_1d[i] > corr_1d[icenter]) icenter = i;
		}

		double [] coeff = null;
		double xcenter = icenter;
		double [][] pa_data=null;
		// try 3-point parabola
		if ((icenter >0) && (icenter < (corr_size - 1))) {
			PolynomialApproximation pa = new PolynomialApproximation(debug?5:0); // debugLevel
			double [][] pa_data0 = {
					{icenter - 1,  corr_1d[icenter - 1]},
					{icenter,      corr_1d[icenter    ]},
					{icenter + 1,  corr_1d[icenter + 1]}};
			pa_data = pa_data0;
			coeff = pa.polynomialApproximation1d(pa_data, 2);
			if (coeff != null){
				xcenter = - coeff[1]/(2* coeff[2]);
			}
		}
		icenter = (int) Math.round(xcenter);
		double strength = corr_1d[icenter] / ((corr_size+1) / 2);// scale to ~match regular strength
		double [] rslt1 = {xcenter - (transform_size - 1), strength};
		return rslt1;
	}


	public void createOrtoNotch( // USED in lwir
			double enhortho_width,
			double enhortho_scale,
			boolean debug) {
		for (int i = 0; i < corr_size; i++){
			if ((i < (transform_size - enhortho_width)) || (i > (transform_size - 2 + enhortho_width))) {
				this.ortho_notch_filter[i] = 1.0;
			} else {
				this.ortho_notch_filter[i] = enhortho_scale;
			}
			if (i == (transform_size-1)) this.ortho_notch_filter[i] = 0.0 ; // hardwired 0 in the center
			this.ortho_notch_filter[i] *= Math.sin(Math.PI*(i+1.0)/(2*transform_size));
		}
		if (debug){
			System.out.println("enhortho_width="+ enhortho_width+" enhortho_scale="+ enhortho_scale);
			for (int i = 0; i < corr_size; i++){
				System.out.println(" enh_ortho_scale["+i+"]="+ this.ortho_notch_filter[i]);
			}
		}
	}



}
