package com.elphel.imagej.tileprocessor;
import java.util.ArrayList;

import com.elphel.imagej.common.PolynomialApproximation;
import com.elphel.imagej.common.ShowDoubleFloatArrays;

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
	// This table is used for CLT-based transform of teh correlation results to the pixel domain
	final static int [][] ZI =
		{{ 0,  1,  2,  3},
		 {-1,  0, -3,  2},
		 {-2, -3,  0,  1},
		 { 3, -2, -1,  0}};

// for 8 cameras and 16 pairs. Following data moved from ImageDtt
    // which images to use (0..3 - external, 4..7 - internal)
    public static int getImgMask  (int data){ return (data & 0xff);}
    // which pairs to combine in the combo see PAIRS data
    public static int getPairMask (int data){ return ((data >> 8) & 0xffff);}
    public static int setImgMask  (int data, int mask) {return (data & ~0xff) | (mask & 0xff);}
    public static int setPairMask (int data, int mask) {return (data & ~0xffff00) | ((mask & 0xffff) << 8);}
    public static boolean getForcedDisparity (int data){return (data & 0x1000000) != 0;}
    public static int     setForcedDisparity (int data, boolean force) {return (data & ~0x1000000) | (force?0x1000000:0);}
    public static boolean getOrthoLines (int data){return (data & 0x2000000) != 0;}
    public static int     setOrthoLines (int data, boolean ortho_lines) {return (data & ~0x2000000) | (ortho_lines?0x2000000:0);}
    public static boolean isOrthoPair         (int npair) {return (PAIRS[npair][2]== PAIR_HORIZONTAL) || (PAIRS[npair][2]== PAIR_VERTICAL);}
    public static boolean isDiagonalPair      (int npair) {return (PAIRS[npair][2]== PAIR_DIAGONAL_MAIN) || (PAIRS[npair][2]== PAIR_DIAGONAL_OTHER);}
    public static boolean isHorizontalPair    (int npair) {return  PAIRS[npair][2]== PAIR_HORIZONTAL;}
    public static boolean isVerticalPair      (int npair) {return  PAIRS[npair][2]== PAIR_VERTICAL;}
    public static boolean isDiagonalMainPair  (int npair) {return  PAIRS[npair][2]== PAIR_DIAGONAL_MAIN;}
    public static boolean isDiagonalOtherPair (int npair) {return  PAIRS[npair][2]== PAIR_DIAGONAL_OTHER;}
    public static int     getScaleOfPair      (int npair) {return  PAIRS[npair][3];}

    public static int     getMaskHorizontal(int scale)    {return getMaskType(PAIR_HORIZONTAL, scale);}
    public static int     getMaskVertical(int scale)      {return getMaskType(PAIR_VERTICAL, scale);}
    public static int     getMaskDiagonalMain(int scale)  {return getMaskType(PAIR_DIAGONAL_MAIN, scale);}
    public static int     getMaskDiagonalOther(int scale) {return getMaskType(PAIR_DIAGONAL_OTHER, scale);}
    public static int     getMaskOrtho(int scale)         {return getMaskHorizontal(scale) | getMaskVertical(scale);}
    public static int     getMaskDiagonal(int scale)      {return getMaskDiagonalMain(scale) | getMaskDiagonalOther(scale);}

    private static int     getMaskType(int type, int scale) { // scale <0 = any
    	int bm = 0;
    	for (int i = 0; i <PAIRS.length; i++) if ((PAIRS[i][2]==type) && ((scale < 0) || (scale == PAIRS[i][3]))) bm |= 1 << i;
    	return bm;
    }



    public Correlation2d (
    		  ImageDttParameters  imgdtt_params,
    		  int transform_size,
    		  double wndx_scale, // (wndy scale is always 1.0)
    		  boolean debug) {
    	this.dtt = new DttRad2(transform_size);
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
      }


      public int [] getTransposeAll(boolean diagonal){
    	  if (diagonal) return getTransposeAllDiagonal();
    	  else          return getTransposeAllOrtho();
      }

      public int [] getTransposeAllOrtho(){
    	  if (this.transpose_all_ortho[0] == this.transpose_all_ortho[1]) {
    		  for (int i =0; i < corr_size; i++){
    			  for (int j =0; j < corr_size; j++){
    				  this.transpose_all_ortho[i * corr_size + j] = j * corr_size + i;
    			  }
    		  }
    	  }
    	  return this.transpose_all_ortho;
      }

      public int [] getTransposeAllDiagonal(){
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
     * Multiply CLT data of two channels, normalize amplitude
     * @param clt_data1 first operand FD CLT data[4][transform_len]
     * @param clt_data2 second operand FD CLT data[4][transform_len]
     * @param fat_zero add to normalization amplitude
     * @return [4][transform_len] FD CLT data
     */
    public double[][] correlateSingleColorFD(
    		double [][] clt_data1,
    		double [][] clt_data2,
    		double [][] tcorr, // null or initialized to [4][transform_len]
    		double      fat_zero) {
    	if (tcorr == null) tcorr = new double [4][transform_len];
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

    public double[][] correlateSingleColorFD_old(
    		double [][] clt_data1,
    		double [][] clt_data2,
    		double [][] tcorr, // null or initialized to [4][transform_len]
    		double      fat_zero) {
    	if (tcorr == null) tcorr = new double [4][transform_len];
		for (int i = 0; i < transform_len; i++) {
			double s1 = 0.0, s2=0.0;
			for (int n = 0; n< 4; n++){
				s1+=clt_data1[n][i] * clt_data1[n][i];
				s2+=clt_data2[n][i] * clt_data2[n][i];
			}
			double scale = 1.0 / (Math.sqrt(s1*s2) + fat_zero*fat_zero); // squared to match units
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
     * @param col_weights [3] - color weights {R, B, G} - green is last, normalized to sum =1.0
     * @param fat_zero fat zero for phase correlation (0 seems to be OK)
     * @return correlation result [(2*transform_size-1) * (2*transform_size-1)]
     */
    public double[]       correlateCompositeFD(
    		double [][][] clt_data1,
    		double [][][] clt_data2,
    		double []     lpf,
    		double []     col_weights,
    		double        fat_zero) {

//    	if ((clt_data1 == null) || (clt_data1 == null)) return null;
    	if (clt_data1.length == 1) { // monochrome
    		col_weights = new double[1];
    		col_weights[0] = 1.0;
    	}
    	double [][][]tcorr = new double [clt_data1.length][4][transform_len];
    	for (int col = 0; col < tcorr.length; col++) {
    		 correlateSingleColorFD(
    		    		clt_data1[col],
    		    		clt_data2[col],
    		    		tcorr[col],
    		    		fat_zero);
    		 if (col == 0) { // accummulate all channels in color 0
    			 for (int n = 0; n<4; n++) {
    				 for (int i = 0; i < transform_len; i++) {
    					 tcorr[0][n][i] *= col_weights[col];
    				 }
    			 }
    		 } else {
    			 for (int n = 0; n<4; n++) {
    				 for (int i = 0; i < transform_len; i++) {
    					 tcorr[0][n][i] += tcorr[col][n][i] * col_weights[col];
    				 }
    			 }
    		 }
    	}
    	if (lpf != null) {
    		for (int n = 0; n<4; n++) {
    			for (int i = 0; i < transform_len; i++) {
    				tcorr[0][n][i] *= lpf[i];
    			}
    		}
    	}

    	for (int quadrant = 0; quadrant < 4; quadrant++){
    		int mode = ((quadrant << 1) & 2) | ((quadrant >> 1) & 1); // transpose
    		tcorr[0][quadrant] = dtt.dttt_iie(tcorr[0][quadrant], mode, transform_size);
    	}
		// convert from 4 quadrants to 15x15 centered tiles (only composite)
    	double [] corr_pd =  dtt.corr_unfold_tile(tcorr[0],	transform_size);

    	return corr_pd;
    }
    /**
     * Calculate all required image pairs phase correlation
     * @param clt_data aberration-corrected FD CLT data [camera][color][tileY][tileX][quadrant][index]
     * @param tileX tile to extract X index
     * @param tileY tile to extract Y index
     * @param pairs_mask bimask of required pairs
     * @param lpf optional low-pass filter
     * @param col_weights RBG color weights
     * @param fat_zero fat zero for phase correlations
     * @return [pair][corr_index]
     */
    public double [][]  correlateCompositeFD(
    		double [][][][][][] clt_data,
    		int                 tileX,
    		int                 tileY,
    		int                 pairs_mask,
    		double []           lpf,
    		double []           col_weights,
    		double              fat_zero) {
    	double [][][][]     clt_data_tile = new double[clt_data.length][][][];
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
    		    		col_weights,
    		    		fat_zero);
    }

    /**
     * Calculate all required image pairs phase correlation
     * @param clt_data aberration-corrected FD CLT data for one tile [camera][color][quadrant][index]
     * @param pairs_mask bimask of required pairs
     * @param lpf optional low-pass filter
     * @param col_weights RBG color weights
     * @param fat_zero fat zero for phase correlations
     * @return [pair][corr_index]
     */
    public double [][]  correlateCompositeFD(
    		double [][][][]     clt_data_tile,
    		int                 pairs_mask, // already decoded so bit 0 - pair 0
    		double []           lpf,
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
     * @param col_weights RBG color weights
     * @param fat_zero fat zero for phase correlations
     * @return 2-d correlation array in line scan order
     */
    public double []  correlateInterCamerasFD(
    		double [][][][]     clt_data_tile_main,
    		double [][][][]     clt_data_tile_aux,
    		double []           lpf,
    		double []           col_weights,
    		double              fat_zero) {
    	if ((clt_data_tile_main == null) || (clt_data_tile_aux == null)) return null;
    	double [][][] clt_mix_main = cltMixCameras(clt_data_tile_main);
    	double [][][] clt_mix_aux =  cltMixCameras(clt_data_tile_aux);
    	double [] inter_cam_corr = correlateCompositeFD(
    			clt_mix_main,         // double [][][] clt_data1,
    			clt_mix_aux,          // double [][][] clt_data2,
	    		lpf,                  // double []     lpf,
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
    public double [][][] cltMixCameras(
    		double [][][][]     clt_data_tile){
    	int tlen = transform_size * transform_size;
    	double [][][] clt_mix = new double [clt_data_tile[0].length][4][tlen];
    	for (int color = 0; color < clt_mix.length; color++) {
    		for (int cltq = 0; cltq <4; cltq++) {
    			for (int i = 0; i < tlen; i++) {
    				for (int cam = 0; cam < clt_data_tile.length; cam++)
    				clt_mix[color][cltq][i] += clt_data_tile[cam][color][cltq][i];
    			}
    		}
    	}
    	double k = 1.0/clt_data_tile.length;
    	for (int color = 0; color < clt_mix.length; color++) {
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

    public double [] combineCompatiblePairs(
    		double [][] correlations,
    		int         pairs_mask,
        	boolean     diagonal,
        	int         baseline_scale
    		) {
    	int width = 2 * transform_size - 1;
    	double [] combo = new double [width * width];
    	int number_combined = 0;
    	// find diagonal/ortho and scale that determine compatible correlations
    	for (int npair = 0; npair < PAIRS.length; npair++) if ((((pairs_mask >> npair) & 1) != 0 ) && (correlations[npair]!=null) &&
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
    public int [] getNumberBaseOfCompatiblePairs(
    		double [][] correlations,
    		int         pairs_mask,
        	boolean     diagonal,
        	int         baseline_scale
    		) {
    	int number_combined = 0;
    	// find diagonal/ortho and scale that determine compatible correlations
    	int base_pair = -1;
    	for (int npair = 0; npair < PAIRS.length; npair++) {
    		if ((isDiagonalPair(npair) == diagonal) && (PAIRS[npair][3] == baseline_scale)){
    			if (base_pair < 0) base_pair = npair;
    			if ((((pairs_mask >> npair) & 1) != 0 ) && (correlations[npair]!=null)){
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
    public int [] getMinMaxSubSample(int pairs_mask) {
    	int ss_max = 0;
    	for (int npair = 0; npair < PAIRS.length; npair++) if (((pairs_mask >> npair) & 1) != 0 ) {
    		if (PAIRS[npair][3] > ss_max) {
    			ss_max = PAIRS[npair][3];
    		}
    	}
    	int ss_min = ss_max;
    	for (int npair = 0; npair < PAIRS.length; npair++) if (((pairs_mask >> npair) & 1) != 0 ) {
    		if (PAIRS[npair][3] < ss_min) {
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
     * @param hwidth number of the result rows (1 - only main diagonal, 2 - main diagonal end 2 other color ones
     * @return transformed array of correlation arrays [hwidth][2*transform_size-1] (some may be nulls)
     */
    public double [][] scaleRotateInterpoateCorrelations(
    		double [][] correlations,
    		int         pairs_mask,
//    		int         sub_sampling,
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
     * @param hwidth number of the result rows (1 - only main diagonal, 2 - main diagonal end 2 other color ones
     * @return transformed correlation array [hwidth][2*transform_size-1]
     */

    public double [] scaleRotateInterpoateSingleCorrelation(
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
    public double [] scaleRotateInterpoateSingleCorrelation(
    		double []   corr,
    		int         hwidth,
    		int         dir, // 0 - hor, 1 - vert, 2 - parallel to row = col (main) diagonal (0->3), 3 -2->1
    		int         ss,
    		boolean     debug
    		) {
//    	int dir = PAIRS[npair][2]; // 0 - hor, 1 - vert, 2 - parallel to row = col (main) diagonal (0->3), 3 -2->1
//   	int ss =  PAIRS[npair][3]/sub_sampling;
    	int center = transform_size - 1;
    	int width = 2 * center + 1;
    	double [] strip = new double [hwidth * width];
    	int xnum=0,ynum=0;
    	int denom =  ss * ((dir > 1)?1:2);
    	double rdenom = denom;
    	int ilimit = center * denom;
//    	double [] corr = correlations[npair];
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
//    					ynum = scol + ((-down * row - 1) >> 1);
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
    						if (int_y) {
    							d = (1.0-dx)*corr[iy0*width + ix0] + dx*corr[iy0 * width + ix0 + 1];
    						} else if (int_x) {
    							d = (1.0-dy)*corr[iy0*width + ix0] + dy*corr[iy0 * width + ix0 + width];
    						} else { // bilinear
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
// todo: if there are no diagonals - why interpolate?
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
    public double [] combineInterpolatedCorrelations(
    		double [][] strips,
    		int         pairs_mask,
    		double      offset,
    		boolean     twice_diagonal){
//    	int center = transform_size - 1;
 //   	int width = 2 * center + 1;
    	double [] combo = null;
    	int ncombined = 0;
    	if (offset >= 0) { // use shifted multiplication
    		for (int npair = 0; npair < strips.length; npair++) if (((pairs_mask & (1 << npair)) != 0) && (strips[npair] != null)){
    			if (combo == null) {
    				combo = new double [strips[npair].length];
    				for (int i = 0; i < combo.length; i++) combo[i] = 1.0;
    			}
    			if (twice_diagonal) {
    				for (int i = 0; i < combo.length; i++) {
    					double d = strips[npair][i] + offset;
    					if (d < 0.0) d = 0.0;
    					combo[i] *= d * d;
    				}
    				ncombined++;
    			} else {
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
    	} else { // use addition
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
 * Find maximum correlation on the grid
 * @param data correlation data - full square or just a combined strip with axis at row 0
 * @param axis_only look for the maximum on the disparity axis only. For the rectangular strips and symmetrical forced to be true.
 * @param minMax minimal value of the maximum to be considered valid
 * @param debug print debug data
 * @return a pair of {x,y} or null. x, y are 0 in the center, disparity is -x
 */


	public int [] getMaxXYInt( // find integer pair or null if below threshold
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
			center_row = center;
		} else {
			axis_only = true;
		}
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
		} else { // only for the the square tile
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
	public double [] getMaxXYCm(
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
//		int [] icenter0 = {icenter[0]+center,icenter[1]+center};
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
	 * Calculate 1-d maximum location, strength and half-width for the special strip (odd rows shifted by 0.5
	 * Negative values are ignored!
	 * Both x and y half-windows can be variable length (to reduce calculations with 0.0 elements), normalized
	 * so sums of zero element and twice all others are 1.0
	 * Window in Y direction corresponds to correlation stripe rows, corresponding to sqrt(2)/2 sensor pixels
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
	public double [] getMaxXCm( // get fractional center as a "center of mass" inside circle/square from the integer max
			double [] data,      // [data_size * data_size]
			int       ixcenter,  // integer center x
			boolean   debug) {
		return getMaxXCm( // get fractional center as a "center of mass" inside circle/square from the integer max
				data,      // double [] data,      // [data_size * data_size]
				ixcenter,  // int       ixcenter,  // integer center x
				this.corr_wndy, // double [] window_y,  // (half) window function in y-direction(perpendicular to disparity: for row0  ==1
				this.corr_wndx, // double [] window_x,  // half of a window function in x (disparity) direction
				debug);// boolean   debug);
	}
	public double [] getMaxXCmNotch( // get fractional center as a "center of mass" inside circle/square from the integer max
			double [] data,      // [data_size * data_size]
			int       ixcenter,  // integer center x
			boolean   debug) {
		return getMaxXCm(             // get fractional center as a "center of mass" inside circle/square from the integer max
				data,                 // double [] data,      // [data_size * data_size]
				ixcenter,             // int       ixcenter,  // integer center x
				this.corr_wndy_notch, // double [] window_y,  // (half) window function in y-direction(perpendicular to disparity: for row0  ==1
				this.corr_wndx,       // double [] window_x,  // half of a window function in x (disparity) direction
				debug);               // boolean   debug);
	}
	public double [] getMaxXCm( // get fractional center as a "center of mass" inside circle/square from the integer max
			double [] data,      // rectangular strip of 1/2 of the correlation are with odd rows shifted by 1/2 pixels
			int       ixcenter,  // integer center x
			double [] window_y,  // (half) window function in y-direction(perpendicular to disparity: for row0  ==1
			double [] window_x,  // half of a window function in x (disparity) direction
			boolean   debug) {
		int center = transform_size - 1;
		int data_width = 2 * transform_size - 1;
		int data_height = data.length/data_width;
		double wy_scale = 1.0;
		if (data_height > window_y.length) {
			data_height = window_y.length;
		} else if (data_height < window_y.length) { // re-
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
			int long_width = 2 * (2 * transform_size-1);
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

	public double [] halfFlatTopWindow(
			int     ihwidth,
			double  hwidth,
			double  blur,
			boolean normalize,
			boolean notch,
			double  scale) {
		double [] wnd = new double [ihwidth];
		for (int i = 0; i < ihwidth; i++) {
			if (blur <= 0.0) {
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
    public double [] corrCenterValues(
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
    public void corrCenterValues(
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
    					if ((fc0 > 0.0) && (fc1 > 0.0)) {
    						cc = Math.sqrt((fc0+offset)*(fc1+offset)) - offset;
    					}
    				} else {
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
    public void saveMlTile(
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
    public void saveMlTilePixel(
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
    public double restoreMlTilePixel(
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

    public  int getMlTilePixelIndex(
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



    public double [] debugStrip(
    		double [] strip) {
    	if (strip == null) return null;
    	int center = transform_size - 1;
    	int width = 2 * center + 1;
    	int height =  strip.length/width;
    	double [] padded_strip = new double [width*width];
    	for (int row = 0; row < width; row++) {
    		if ((row <= center - height) || (row >= center+ height)) {
    			for (int j = 0; j<width; j++) padded_strip[row*width+j] = Double.NaN;
    		} else {
    			int srow = (row >= center)? (row - center) : (center - row);
    			for (int j = 0; j<width; j++) padded_strip[row*width+j] = strip[srow*width+j];
    		}
    	}

    	return padded_strip;
    }

    // only show center part, but with correct shift
    public double [] debugStrip2(
    		double [] strip) {
    	if (strip == null) return null;
    	int center = transform_size - 1;
    	int width = 2 * center + 1;
    	int height =  strip.length/width;
    	double [] padded_strip = new double [width*width];
    	for (int row = 0; row < width; row++) {
    		if ((row <= center - height) || (row >= center+ height)) {
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
    public double [] debugStrip3(
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

    public double [] mismatchPairsCM( // returns x-xcenter, y, strength (sign same as disparity)
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
			} else {
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
			} else {
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
    public double [] mismatchPairs( // returns x-xcenter, y, strength (sign same as disparity)
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
//    	for (int np = 0; np < num_pairs; np++) {
    		int this_mask = 1 << pair;
    		if (debug_level > -1) {
    			System.out.println(String.format("mismatchPairs(), np = %d pairs mask = 0x%x", np, this_mask));
    		}
    		Correlations2dLMA lma=corrLMA(
    				imgdtt_params, // ImageDttParameters  imgdtt_params,
    				corrs,         // double [][]         corrs,
    				this_mask,     // int                 pair_mask, // which pairs to process
    				true,          // boolean             run_poly_instead, // true - run LMA, false - run 2d polynomial approximation
    				xcenter,       // double              xcenter,   // preliminary center x in pixels for largest baseline
    				vasw_pwr,      // double              vasw_pwr,  // value as weight to this power,
//    				debug_level-3, // int                 debug_level,
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
//    				rslt[3 * np + 0] = xcenter - poly_xyvwh[0];
//    				rslt[3 * np + 1] = -poly_xyvwh[1];
    				rslt[3 * np + 0] =  xcenter - poly_xyvwh[0] + poly_xyvwh[1]; // x - y
    				rslt[3 * np + 1] =  xcenter - poly_xyvwh[0] - poly_xyvwh[1]; // x + y
    			} else if (isDiagonalOtherPair(pair)) {
//    				rslt[3 * np + 0] = xcenter - poly_xyvwh[0];
//    				rslt[3 * np + 1] =  poly_xyvwh[1];
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
    public double [] single2dPoly( // returns x-xcenter, y, strength (sign same as disparity)
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
    public double [] single2dCM( // returns x-xcenter, y, strength (sign same as disparity)
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



    public double [][] corr4dirsLMA(
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
    		if (this_mask == 0) {
    			rslt[dir] = null;
    		} else {
    			Correlations2dLMA lma=corrLMA(
    					imgdtt_params, // ImageDttParameters  imgdtt_params,
    					corrs,         // double [][]         corrs,
    		    		this_mask,     // int                 pair_mask, // which pairs to process
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
    public double [] foregroundCorrect(
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
			if (debug) System.out.println("Direction with "+(bg?"min":"max")+" has effective strenbgth ratio too small -> no correction ("+
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
			if (!(mx >= (mn + min_diff))) { // so NaN are OK
				corr = Double.NaN;
	    		if (debug) System.out.println("difference max - min="+(mx - mn)+" < "+min_diff+" -> no fo correction");
			}
		}
    	double disp = full_disp;
    	if (!Double.isNaN(corr)) {
    		double lim;
    		if (bg) {
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



    public Correlations2dLMA corrLMA(
    		ImageDttParameters  imgdtt_params,
    		double [][]         corrs,
    		int                 pair_mask, // which pairs to process
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
    				true, (run_poly_instead?("mismatch"+pair_mask):"groups")+"_x"+tileX+"_y"+tileY, titles);
    	}
    	if (debug_level > 1) {
    		System.out.println("Input data and approximation:");
    		lma.printInputDataFx(true);
    	}
    	return lmaSuccess? lma: null;
    }

// Run for a single horizontal 2d correlation array
    public Correlations2dLMA corrLMA(
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
    public double [] filterConvex(
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
//    	int width_m1 =     width - 1;
//    	int win = 2 * hwin -1;
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

    	// modify

//    	int [][] dirs8 = {
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

    public void debug_convex(
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

    public void addSamples(
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
//    				double wy = window_y[arow] * groups_pairs[0][0]; // number of pair averaged
    				for (int acol =  odd; acol < hwindow_x2; acol +=2) {
//    					double wxy = window_x[acol] * wy;  // full weight before value as weight
    					for (int quad = 0; quad < 4; quad ++) if (((arow > 0) || ((quad & 2) !=0 )) && ((acol > 0) || ((quad & 1) !=0 ))){
    						int cx = (quad_signs[quad][0] * acol - quad_signs[quad][1] * arow)/2 + ixcenter; // ix0;
    						int cy = (quad_signs[quad][0] * acol + quad_signs[quad][1] * arow)/2 + ixcenter; // ix0;
    						// calculate coordinates in the correlation array
    						if ((cx >= -center) && (cx <= center) && (cy >= -center) && (cy <= center)) {
//    							double w = wxy;  // full weight before value as weight
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
//    				double wy = window_y[arow] * groups_pairs[0][0]; // number of pair averaged
    				for (int acol =  0; acol < hwindow_x2; acol +=2) {
//    					double wxy = window_x[acol] * wy;  // full weight before value as weight
    					for (int quad = 0; quad < 4; quad ++) if (((arow > 0) || ((quad & 2) !=0 )) && ((acol > 0) || ((quad & 1) !=0 ))){
    						int cx = (quad_signs[quad][0] * acol)/2 + ixcenter; // ix0;
    						int cy = (quad_signs[quad][1] * arow)/2;
    						// calculate coordinates in the correlation array
    						if ((cx >= -center) && (cx <= center) && (cy >= -center) && (cy <= center)) {
//    							double w = wxy;  // full weight before value as weight
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
    public double [][] getFitSamples( // just for debug to compare LMA-fitted fx with original data
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
//    	int hwindow_y = window_y.length; // should actually be the same?
//    	int hwindow_x = window_x.length;
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

    public int [] listPairs(
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
//  		public double     max_corr_radius =   3.9;  // maximal distance from int max to consider
	public double [] getMaxXSOrtho( // // get fractional center using a quadratic polynomial
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
			} else { // use averaging (linear)
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
			} else { // use averaging (linear)
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


	public void createOrtoNotch(
			double enhortho_width,
			double enhortho_scale,
			boolean debug) {
//		int corr_size = transform_size * 2 -1;
//		double [] ortho_notch = new double [corr_size];
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
