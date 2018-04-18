/**
 **
 ** Correlation2d - Handle 2-d (phase) correlations, combining multiple-pair data
 **
 ** Copyright (C) 2018 Elphel, Inc.
 **
 ** -----------------------------------------------------------------------------**
 **
 **  ImageDtt.java is free software: you can redistribute it and/or modify
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
	private final int [][] transpose_indices_ortho;
	private final int [][] transpose_indices_diagonal;

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

//    public final ImageDtt imagedtt;

    // imagedtt.transform_size
//    public Correlation2d (ImageDtt imagedtt, int transform_size) {
      public Correlation2d (int transform_size) {
//    	this.imagedtt = imagedtt;
    	this.dtt = new DttRad2(transform_size);
    	this.transform_size = transform_size;
    	this.transform_len = transform_size * transform_size;
    	this.corr_size = transform_size * 2 -1;
    	// not initialized until needed
    	this.transpose_indices_ortho =  new int [corr_size*(corr_size-1)/2][];
    	this.transpose_indices_diagonal =  new int [corr_size*(corr_size-1)/2][];

      }

      public int [][] getTransposeIndices(boolean diagonal){
    	  if (diagonal) return getTransposeIndicesDiagonal();
    	  else          return getTransposeIndicesOrtho();
      }

      public int [][] getTransposeIndicesOrtho(){
    	  if (this.transpose_indices_ortho[0] == null) {
    		  int indx = 0;
    		  for (int i =0; i < corr_size-1; i++){
    			  for (int j = i+1; j < corr_size; j++){
    				  transpose_indices_ortho[indx  ] = new int[2];
    				  transpose_indices_ortho[indx  ][0] = i * corr_size + j;
    				  transpose_indices_ortho[indx++][1] = j * corr_size + i;
    			  }
    		  }
    	  }
    	  return this.transpose_indices_ortho;
      }

      public int [][] getTransposeIndicesDiagonal(){
    	  if (this.transpose_indices_diagonal[0] == null) {

    		  int indx = 0;
    		  for (int i =0; i < transform_size-1; i++){
    			  for (int j = 0; j < corr_size; j++){
    				  transpose_indices_diagonal[indx  ] = new int[2];
    				  transpose_indices_diagonal[indx  ][0] = i                  * corr_size + j;
    				  transpose_indices_diagonal[indx++][1] = (corr_size - i -1) * corr_size + j;
    			  }
    		  }
    	  }
    	  return this.transpose_indices_diagonal;
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
     * @param clt_data aberration-corrected FD CLT data [camera][color][quadrant][index]
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
    			int [][] transpose_indices = getTransposeIndices(isDiagonalOtherPair(npair));
				for (int i = 0; i < transpose_indices.length; i++) {
					double d = correlations[npair][transpose_indices[i][0]];
					correlations[npair][transpose_indices[i][0]] = correlations[npair][transpose_indices[i][1]];
					correlations[npair][transpose_indices[i][1]] = d;
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
     * @return number of compatible pairs among the selection
     */
    public int getNumberOfCompatiblePairs(
    		double [][] correlations,
    		int         pairs_mask,
        	boolean     diagonal,
        	int         baseline_scale
    		) {
    	int number_combined = 0;
    	// find diagonal/ortho and scale that determine compatible correlations
    	for (int npair = 0; npair < PAIRS.length; npair++) if ((((pairs_mask >> npair) & 1) != 0 ) && (correlations[npair]!=null) &&
    		(isDiagonalPair(npair) == diagonal) && (PAIRS[npair][3] == baseline_scale)){
    		number_combined++;
    	}
    	return number_combined;
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
    	int dir = PAIRS[npair][2]; // 0 - hor, 1 - vert, 2 - parallel to row = col (main) diagonal (0->3), 3 -2->1
    	int ss =  PAIRS[npair][3]/sub_sampling;
    	int center = transform_size - 1;
    	int width = 2 * center + 1;
    	double [] strip = new double [hwidth * width];
    	int xnum=0,ynum=0;
    	int denom =  ss * ((dir > 1)?1:2);
    	double rdenom = denom;
    	int ilimit = center * denom;
    	double [] corr = correlations[npair];
		if (debug) {
			System.out.println("\n============== scaleRotateInterpoateSingleCorrelation(),npair ="+ npair+", sub_sampling="+sub_sampling+" ===============");
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
     * @param pairs_mask
     * @param offset
     * @return
     */
    public double [] combineInterpolatedCorrelations(
    		double [][] strips,
    		int         pairs_mask,
    		double      offset){
//    	int center = transform_size - 1;
 //   	int width = 2 * center + 1;
    	double [] combo = null;
    	int ncombined = 0;
    	for (int npair = 0; npair < strips.length; npair++) if (((pairs_mask & (1 << npair)) != 0) && (strips[npair] != null)){
    		if (combo == null) {
    			combo = new double [strips[npair].length];
    			for (int i = 0; i < combo.length; i++) combo[i] = 1.0;
    		}
			for (int i = 0; i < combo.length; i++) {
				combo[i] *= (strips[npair][i] + offset);
			}
    		ncombined++;
    	}
    	double pwr = 1.0/ncombined;
		for (int i = 0; i < combo.length; i++) {
			combo[i] = Math.pow(combo[i], pwr) - offset;
		}
    	return combo;
    }

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
		if (data[imx] < minMax) {
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

		if (debug) {
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
						if (!Double.isNaN(d)) {
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
	 * @param scale  multiply each value
	 * @return half-window array [ihwidth]
	 */

	public double [] halfFlatTopWindow(
			int     ihwidth,
			double  hwidth,
			double  blur,
			boolean normalize,
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




    public double [] debugStrip(
    		double [] strip) {
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

//    				padded_strip[row*width+j] = strip[srow*width+j];
    			}
    		}
    	}

    	return padded_strip;
    }

/*
    public static int getPairMask (int data){ return ((data >> 8) & 0xffff);}
      * Calculate color channels FD phase correlations, mix results with weights, apply optional low-pass filter
     * and convert to the pixel domain  as [(2*transform_size-1) * (2*transform_size-1)] tiles (15x15)
     * No transposing or rotation
     * @param clt_data1 [3][4][transform_len] first operand data. First index - RBG color
     * @param clt_data2 [3][4][transform_len] first operand data. First index - RBG color
     * @param lpf   optional [transform_len] LPF filter data
     * @param col_weights [3] - color weights {R, B, G} - green is last, normalized to sum =1.0
     * @param fat_zero fat zero for phase correlation (0 seems to be OK)
     * @return correlation result [(2*transform_size-1) * (2*transform_size-1)]

		final double [][][][][][] clt_data = new double[quad][nChn][tilesY][tilesX][][];

								// transpose vertical pairs
								if (corr_pairs[pair][2] != 0) {
									for (int chn = firstColor; chn <= numcol; chn++){
										for (int i = 0; i < transpose_indices.length; i++) {
											double d = tcorr_partial[pair][chn][transpose_indices[i][0]];
											tcorr_partial[pair][chn][transpose_indices[i][0]] = tcorr_partial[pair][chn][transpose_indices[i][1]];
											tcorr_partial[pair][chn][transpose_indices[i][1]] = d;
											//transpose_indices
										}
									}
								}
								// make symmetrical around the disparity direction (horizontal) (here using just average, not mul/sum mixture)
								// symmetry can be added to result, not individual (if sum - yes, but with multiplication - not)
								if (corr_sym && (clt_mismatch == null)){ // when measuring clt_mismatch symmetry should be off !
									for (int chn = firstColor; chn <= numcol; chn++){
										for (int i = 1 ; i < transform_size; i++){
											int indx1 = (transform_size - 1 - i) * corr_size;
											int indx2 = (transform_size - 1 + i) * corr_size;
											for (int j = 0; j< corr_size; j++){
												int indx1j = indx1 + j;
												int indx2j = indx2 + j;
												tcorr_partial[pair][chn][indx1j] =
														0.5* (tcorr_partial[pair][chn][indx1j] + tcorr_partial[pair][chn][indx2j]);
												tcorr_partial[pair][chn][indx2j] = tcorr_partial[pair][chn][indx1j];
											}
										}
									}
								}

------
								for (int chn = 0; chn <numcol; chn++){
									double [][] data1 = clt_data[corr_pairs[pair][0]][chn][tileY][tileX];
									double [][] data2 = clt_data[corr_pairs[pair][1]][chn][tileY][tileX];
									for (int i = 0; i < transform_len; i++) {
										double s1 = 0.0, s2=0.0;
										for (int n = 0; n< 4; n++){
											s1+=data1[n][i] * data1[n][i];
											s2+=data2[n][i] * data2[n][i];
										}
										double scale = 1.0 / (Math.sqrt(s1*s2) + corr_fat_zero*corr_fat_zero); // squared to match units
										for (int n = 0; n<4; n++){
											tcorr_tpartial[pair][chn][n][i] = 0;
											for (int k=0; k<4; k++){
												if (zi[n][k] < 0)
													tcorr_tpartial[pair][chn][n][i] -=
															data1[-zi[n][k]][i] * data2[k][i];
												else
													tcorr_tpartial[pair][chn][n][i] +=
													data1[zi[n][k]][i] * data2[k][i];
											}
											tcorr_tpartial[pair][chn][n][i] *= scale;
										}
									}
									// got transform-domain correlation for the pair, 1 color
								}

								// calculate composite color
								for (int i = 0; i < transform_len; i++) {
									for (int n = 0; n<4; n++) {
										tcorr_tpartial[pair][numcol][n][i] =
												col_weights[0]* tcorr_tpartial[pair][0][n][i] +
												col_weights[1]* tcorr_tpartial[pair][1][n][i] +
												col_weights[2]* tcorr_tpartial[pair][2][n][i];
									}
								}
								// now lpf (only last/composite color if do not preserve intermediate
								int firstColor = (clt_corr_partial == null)? numcol : 0;
								if (corr_sigma >0) {
									for (int chn = firstColor; chn <= numcol; chn++){
										for (int i = 0; i < transform_len; i++) {
											for (int n = 0; n<4; n++) {
												tcorr_tpartial[pair][chn][n][i] *= filter[i];
											}
										}
									}
								}
								// convert to pixel domain - all or just composite color
								for (int chn = firstColor; chn <= numcol; chn++){
									for (int quadrant = 0; quadrant < 4; quadrant++){
										int mode = ((quadrant << 1) & 2) | ((quadrant >> 1) & 1); // transpose
										tcorr_tpartial[pair][chn][quadrant] =
												dtt.dttt_iie(tcorr_tpartial[pair][chn][quadrant], mode, transform_size);
									}
								}
								// convert from 4 quadrants to 15x15 centered tiles (each color or only composite)
								for (int chn = firstColor; chn <= numcol; chn++){
									tcorr_partial[pair][chn] = corr_unfold_tile(
											tcorr_tpartial[pair][chn],
											transform_size);
								}


 */


}
