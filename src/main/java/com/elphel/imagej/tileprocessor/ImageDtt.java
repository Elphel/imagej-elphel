package com.elphel.imagej.tileprocessor;
/**
 **
 ** ImageDtt - Process images with DTT-based methods
 **
 ** Copyright (C) 2016 Elphel, Inc.
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
import java.util.concurrent.atomic.AtomicInteger;

import com.elphel.imagej.cameras.CLTParameters;
import com.elphel.imagej.cameras.EyesisCorrectionParameters;
import com.elphel.imagej.common.PolynomialApproximation;
import com.elphel.imagej.common.ShowDoubleFloatArrays;
import com.elphel.imagej.correction.EyesisDCT;

import Jama.Matrix;
import ij.ImageStack;

public class ImageDtt {
	  static boolean FPGA_COMPARE_DATA= false; // true; // false; //
	  static int     FPGA_SHIFT_BITS =  7; // number of bits for fractional pixel shift
	  static int     FPGA_PIXEL_BITS = 15; // bits to represent pixel data (positive)
	  static int     FPGA_WND_BITS =   17; // bits to represent mclt window (positive for 18-bit signed mpy input)
	  static int     FPGA_DTT_IN =     22; // bits to represent maximal value after folding (input to DTT)
	  static int     FPGA_TILE_SIZE =  22; // size of square side for the composite colors tile (16..22)
	  static double [] kern_g={
			  0.0,   0.125,  0.0  ,
			  0.125, 0.5,    0.125,
			  0.0,   0.125,  0.0  };
	  static double [] kern_rb={
			  0.0625,  0.125, 0.0625,
			  0.125,   0.25,  0.125,
			  0.0625,  0.125, 0.0625};
//	  static double [][] kerns = {kern_rb,kern_rb,kern_g};
	  static int [][] corn_side_indices = { // of 012/345/678 3x3 square kernel
			  {4,5,7,8},           // top left corner
			  {3,4,5,6,7,8},       // top middle
			  {3,4,6,7},           // top right
			  {1,2,4,5,7,8},       // middle left
			  {0,1,2,3,4,5,6,7,8}, // middle
			  {0,1,3,4,6,7},       // middle right
			  {1,2,4,5},           // bottom left
			  {0,1,2,3,4,5},       // bottom middle
			  {0,1,3,4}};          // bottom right
//	 public static int FORCE_DISPARITY_BIT = 8; // move to parameters?

	  static int  QUAD =                           4; // number of cameras in camera
	  static int  GREEN_CHN =                      2; // index of green channel
	  static int  MONO_CHN =                       2; // index of channel used in monochrome mode

	  static int  DISPARITY_INDEX_INT =            0; // 0 - disparity from correlation integer pixels, 1 - ortho
	  // Now DISPARITY_INDEX_CM may be POLY with backup from CM (for bad correlation)
	  static int  DISPARITY_INDEX_CM =             2; // 2 - disparity from correlation "center mass", 3 - ortho (only used for fine correction)
	  static int  DISPARITY_INDEX_HOR =            4; // disparity from correlation of the horizontal pairs with center suppressed
	  static int  DISPARITY_INDEX_HOR_STRENGTH =   5; // strength for hor mode (emphasis on vertical lines)
	  static int  DISPARITY_INDEX_VERT =           6; // disparity from correlation of the vertical pairs with center suppressed
	  static int  DISPARITY_INDEX_VERT_STRENGTH =  7; // strength in vert mode (horizontal lines detection)
	  static int  DISPARITY_INDEX_POLY =           8; // index of disparity value in disparity_map == 2 (0,2 or 4)
	  static int  DISPARITY_STRENGTH_INDEX =      10; // index of strength data in disparity map ==6
	  static int  DISPARITY_VARIATIONS_INDEX =    11; // index of strength data in disparity map ==6
	  static int  IMG_DIFF0_INDEX =               12; // index of noise- normalized image difference for port 0 in disparity map
	  static int  OVEREXPOSED =                   16; // index of overexposed fraction of all pixels
	  static int  IMG_TONE_RGB =                  17; // 12 entries of r0,r1,r2,r3,g0,g1,g2,g3,b0,b1,b2,b3
// remove when not needed
	  static String [] DISPARITY_TITLES = {
			  "int_disp","int_y_disp","cm_disp","cm_y_disp","hor_disp","hor_strength","vert_disp","vert_strength",
			  "poly_disp", "poly_y_disp", "strength_disp", "vary_disp","diff0","diff1","diff2","diff3","overexp",
			  "r0","r1","r2","r3",
			  "g0","g1","g2","g3",
			  "b0","b1","b2","b3",
			  };
//			  "dbg0","dbg1","dbg2","dbg3","dbg4","dbg5","dbg6","dbg7","dbg8","dbg9","dbg10","dbg11","dbg12","dbg13","dbg14","dbg15","dbg16","dbg17","dbg18"};

	  static int  BI_DISP_FULL_INDEX =            0;  // 0 - disparity for all directions of the main camera
	  static int  BI_DISP_HOR_INDEX =             1;  // 1 - disparity for 2 horizontal pairs of the main camera
	  static int  BI_DISP_VERT_INDEX =            2;  // 2 - disparity for 2 vertical pairs of the main camera
	  static int  BI_DISP_DIAGM_INDEX =           3;  // 3 - disparity for main diagonal pair of the main camera
	  static int  BI_DISP_DIAGO_INDEX =           4;  // 4 - disparity for main diagonal pair of the main camera
	  static int  BI_ADISP_FULL_INDEX =           5;  // 5 - disparity for all directions of the aux camera
	  static int  BI_ADISP_HOR_INDEX =            6;  // 6 - disparity for 2 horizontal pairs of the aux camera
	  static int  BI_ADISP_VERT_INDEX =           7;  // 7 - disparity for 2 vertical pairs of the aux camera
	  static int  BI_ADISP_DIAGM_INDEX =          8;  // 8 - disparity for main diagonal pair of the aux camera
	  static int  BI_ADISP_DIAGO_INDEX =          9;  // 9 - disparity for main diagonal pair of the aux camera
	  static int  BI_DISP_CROSS_INDEX =          10;  //10 - disparity between the main the aux camera
	  static int  BI_DISP_CROSS_DX_INDEX =       11;  //11 - delta disparity between the main the aux camera (horizontal)
	  static int  BI_DISP_CROSS_DY_INDEX =       12;  //12 - delta disparity between the main the aux camera (vertical)
	  static int  BI_STR_FULL_INDEX =            13;  //13 - strength for all directions of the main camera
	  static int  BI_STR_HOR_INDEX =             14;  //14 - strength for 2 horizontal pairs of the main camera
	  static int  BI_STR_VERT_INDEX =            15;  //15 - strength for 2 vertical pairs of the main camera
	  static int  BI_STR_DIAGM_INDEX =           16;  //16 - strength for main diagonal pair of the main camera
	  static int  BI_STR_DIAGO_INDEX =           17;  //17 - strength for main diagonal pair of the main camera
	  static int  BI_ASTR_FULL_INDEX =           18;  //18 - strength for all directions of the aux camera
	  static int  BI_ASTR_HOR_INDEX =            19;  //19 - strength for 2 horizontal pairs of the aux camera
	  static int  BI_ASTR_VERT_INDEX =           20;  //20 - strength for 2 vertical pairs of the aux camera
	  static int  BI_ASTR_DIAGM_INDEX =          21;  //21 - strength for main diagonal pair of the aux camera
	  static int  BI_ASTR_DIAGO_INDEX =          22;  //22 - strength for main diagonal pair of the aux camera
	  static int  BI_STR_CROSS_INDEX =           23;  //23 - strength between the main the aux camera
	  static int  BI_STR_ALL_INDEX =             24;  //24 - average strength (product of strengths to 1/3 power), TODO: strength at cross disparity
	  static int  BI_TARGET_INDEX =              25;  //25 - target disparity
	  static int  BI_DBG1_INDEX =                26;  //26 - debug layer 1
	  static int  BI_DBG2_INDEX =                27;  //27 - debug layer 2
	  static int  BI_DBG3_INDEX =                28;  //28 - debug layer 2
	  static int  BI_DBG4_INDEX =                29;  //29 - debug layer 2

	  static String [] BIDISPARITY_TITLES = {
			  "disparity","disp_hor","disp_vert","disp_diagm","disp_diago",
			  "adisparity","adisp_hor","adisp_vert","adisp_diagm","adisp_diago",
			  "bi-disparity","bi-disparity-dx","bi-disparity-dy",
			  "strength", "str_hor", "str_vert", "str_diagm", "str_diago",
			  "astrength", "astr_hor", "astr_vert", "astr_diagm", "astr_diago",
			  "bi-strength", "all-strength", "target", "dbg1", "dbg2", "dbg3", "dbg4"};
	  static int [] BIDISPARITY_STRENGTHS= {
			  BI_STR_FULL_INDEX,   BI_STR_VERT_INDEX,  BI_STR_DIAGM_INDEX,  BI_STR_DIAGO_INDEX,
			  BI_ASTR_FULL_INDEX,  BI_ASTR_HOR_INDEX,  BI_ASTR_VERT_INDEX,  BI_ASTR_DIAGM_INDEX,
			  BI_ASTR_DIAGO_INDEX, BI_STR_CROSS_INDEX, BI_STR_ALL_INDEX};

	  static int  DISP_FULL_INDEX =            0;  // 0 - disparity for all directions of the main camera
	  static int  DISP_HOR_INDEX =             1;  // 1 - disparity for 2 horizontal pairs of the main camera
	  static int  DISP_VERT_INDEX =            2;  // 2 - disparity for 2 vertical pairs of the main camera
	  static int  DISP_DIAGM_INDEX =           3;  // 3 - disparity for main diagonal pair of the main camera
	  static int  DISP_DIAGO_INDEX =           4;  // 4 - disparity for main diagonal pair of the main camera
	  static int  STR_FULL_INDEX =             5;  //11 - strength for all directions of the main camera
	  static int  STR_HOR_INDEX =              6;  //12 - strength for 2 horizontal pairs of the main camera
	  static int  STR_VERT_INDEX =             7;  //13 - strength for 2 vertical pairs of the main camera
	  static int  STR_DIAGM_INDEX =            8;  //14 - strength for main diagonal pair of the main camera
	  static int  STR_DIAGO_INDEX =            9;  //15 - strength for main diagonal pair of the main camera

	  // ML data

	  static int  ML_TOP_INDEX =               0;  // 0 - top pair 2d correlation center area
	  static int  ML_BOTTOM_INDEX =            1;  // 1 - bottom pair 2d correlation center area
	  static int  ML_LEFT_INDEX =              2;  // 2 - left pair 2d correlation center area
	  static int  ML_RIGHT_INDEX =             3;  // 3 - right pair 2d correlation center area
	  static int  ML_DIAGM_INDEX =             4;  // 4 - main diagonal (top-left to bottom-right) pair 2d correlation center area
	  static int  ML_DIAGO_INDEX =             5;  // 5 - other diagonal (bottom-left to top-right) pair 2d correlation center area
	  static int  ML_HOR_INDEX =               6;  // 6 - horizontal pairs combined 2d correlation center area
	  static int  ML_VERT_INDEX =              7;  // 7 - vertical pairs combined 2d correlation center area

	  static int  ML_TOP_AUX_INDEX =           8;  // 8 - top pair 2d correlation center area (auxiliary camera)
	  static int  ML_BOTTOM_AUX_INDEX =        9;  // 9 - bottom pair 2d correlation center area (auxiliary camera)
	  static int  ML_LEFT_AUX_INDEX =         10;  //10 - left pair 2d correlation center area (auxiliary camera)
	  static int  ML_RIGHT_AUX_INDEX =        11;  //11 - right pair 2d correlation center area (auxiliary camera)
	  static int  ML_DIAGM_AUX_INDEX =        12;  //12 - main diagonal (top-left to bottom-right) pair 2d correlation center area (auxiliary camera)
	  static int  ML_DIAGO_AUX_INDEX =        13;  //13 - other diagonal (bottom-left to top-right) pair 2d correlation center area (auxiliary camera)
	  static int  ML_HOR_AUX_INDEX =          14;  //14 - horizontal pairs combined 2d correlation center area (auxiliary camera)
	  static int  ML_VERT_AUX_INDEX =         15;  //15 - vertical pairs combined 2d correlation center area (auxiliary camera)

	  static int  ML_INTER_INDEX =            16;  //16 - inter-camera (between two quad ones) correlation center area
	  static int  ML_OTHER_INDEX =            17;  //17 - other data: 0 (top left tile corner) - preset disparity of the tile, 1: (next element) - ground trouth data, 2:
	                                               //      ground truth confidence
	  static int  ML_DBG1_INDEX =             18;  //18 - just debug data (first - auto phase correlation)

	  static String [] ML_TITLES = {
			  "top-pair", "bottom-pair", "left_pair", "right-pair", "diagm-pair", "diago-pair","hor-pairs","vert-pairs",
			  "top-aux",  "bottom-aux",  "left_aux",  "right-aux",  "diagm-aux",  "diago-aux", "hor-aux",  "vert-aux",
			  "inter", "other", "dbg1"};

	  public static int  ML_OTHER_TARGET =            0;  // Offset to target disparity data in  ML_OTHER_INDEX layer tile
	  public static int  ML_OTHER_GTRUTH =            2;  // Offset to ground truth disparity data in ML_OTHER_INDEX layer tile
	  public static int  ML_OTHER_GTRUTH_STRENGTH =   4;  // Offset to ground truth confidence data in  ML_OTHER_INDEX layer tile


	  // indices in cross-camera correlation results
	  static int INDEX_DISP =     0;
	  static int INDEX_STRENGTH = 1;
	  static int INDEX_DX =       2;
	  static int INDEX_DY =       3;


	  static String [] SNGL_DISPARITY_TITLES = {
			  "disparity","disp_hor","disp_vert","disp_diagm","disp_diago",
			  "strength", "str_hor", "str_vert", "str_diagm", "str_diago"};
      static int [] SNGL_DISPARITY_NAN = {DISP_FULL_INDEX, DISP_HOR_INDEX, DISP_VERT_INDEX,  DISP_DIAGM_INDEX, DISP_DIAGO_INDEX};

      static int [][] SNGL_TO_BI = {
    		  {BI_DISP_FULL_INDEX, BI_DISP_HOR_INDEX, BI_DISP_VERT_INDEX,BI_DISP_DIAGM_INDEX,BI_DISP_DIAGO_INDEX,
    			  BI_STR_FULL_INDEX, BI_STR_HOR_INDEX, BI_STR_VERT_INDEX,BI_STR_DIAGM_INDEX,BI_STR_DIAGO_INDEX},
    		  {BI_ADISP_FULL_INDEX, BI_ADISP_HOR_INDEX, BI_ADISP_VERT_INDEX,BI_ADISP_DIAGM_INDEX,BI_ADISP_DIAGO_INDEX,
    				  BI_ASTR_FULL_INDEX, BI_ASTR_HOR_INDEX, BI_ASTR_VERT_INDEX,BI_ASTR_DIAGM_INDEX,BI_ASTR_DIAGO_INDEX}};

	  static int  TCORR_COMBO_RSLT =  0; // normal combined correlation from all   selected pairs (mult/sum)
	  static int  TCORR_COMBO_SUM =   1; // sum of channel correlations from all   selected pairs
	  static int  TCORR_COMBO_HOR =   2; // combined correlation from 2 horizontal pairs (0,1). Used to detect vertical features
	  static int  TCORR_COMBO_VERT =  3; // combined correlation from 2 vertical   pairs (0,1). Used to detect horizontal features
	  static String [] TCORR_TITLES = {"combo","sum","hor","vert"};

	  private final boolean monochrome;
	  private final double scale_strengths; // scale all correlation strengths (to compensate for LPF sigma changes)


     public static int getImgMask  (int data){ return (data & 0xf);}      // which images to use
     public static int getPairMask (int data){ return ((data >> 4) & 0xf);} // which pairs to combine in the combo:  1 - top, 2 bottom, 4 - left, 8 - right
     public static int setImgMask  (int data, int mask) {return (data & ~0xf) | (mask & 0xf);}
     public static int setPairMask (int data, int mask) {return (data & ~0xf0) | ((mask & 0xf) << 4);}
     public static boolean getForcedDisparity (int data){return (data & 0x100) != 0;}
     public static int     setForcedDisparity (int data, boolean force) {return (data & ~0x100) | (force?0x100:0);}
     public static boolean getOrthoLines (int data){return (data & 0x200) != 0;}
     public static int     setOrthoLines (int data, boolean force) {return (data & ~0x200) | (force?0x200:0);}

	public ImageDtt(
			boolean mono,
			double scale_strengths){
		this.monochrome = mono;
		this.scale_strengths = scale_strengths;
	}

	public boolean isMonochrome() {
		return monochrome;
	}
	// maybe change in the future
	public boolean isAux() {
		return monochrome;
	}


	public double [][][][] mdctStack(
			final ImageStack                                 imageStack,
			final int                                        subcamera, //
			final EyesisCorrectionParameters.DCTParameters   dctParameters, //
			final EyesisDCT                                  eyesisDCT,
			final int                                        threadsMax, // maximal step in pixels on the maxRadius for 1 angular step (i.e. 0.5)
			final int                                        debugLevel,
			final boolean                                    updateStatus) // update status info

	{
	  	  if (imageStack==null) return null;
		  final int imgWidth=imageStack.getWidth();
		  final int nChn=imageStack.getSize();
		  double [][][][] dct_data = new double [nChn][][][];
		  float [] fpixels;
		  int i,chn; //tileX,tileY;
		  /* find number of the green channel - should be called "green", if none - use last */
		  // Extract float pixels from inage stack, convert each to double

		  EyesisDCT.DCTKernels dct_kernels = null;
		  dct_kernels = ((eyesisDCT==null) || (eyesisDCT.kernels==null))?null:eyesisDCT.kernels[subcamera];
		  if (dct_kernels == null){
			  System.out.println("No DCT kernels available for subcamera # "+subcamera);
		  } else if (debugLevel>0){
			  System.out.println("Using DCT kernels for subcamera # "+subcamera);
		  }
//		  if (dctParameters.kernel_chn >=0 ){
//			  dct_kernels = eyesisDCT.kernels[dctParameters.kernel_chn];
//		  }

		  for (chn=0;chn<nChn;chn++) {
			  fpixels= (float[]) imageStack.getPixels(chn+1);
			  double[] dpixels = new double[fpixels.length];
			  for (i = 0; i <fpixels.length;i++) dpixels[i] = fpixels[i];
			  // convert each to DCT tiles
			  dct_data[chn] =lapped_dct(
						dpixels,
						imgWidth,
						dctParameters.dct_size,
						0, //     dct_mode,    // 0: dct/dct, 1: dct/dst, 2: dst/dct, 3: dst/dst
						dctParameters.dct_window, // final int       window_type,
						chn,
						dct_kernels,
						dctParameters.skip_sym,
						dctParameters.convolve_direct,
						dctParameters.tileX,
						dctParameters.tileY,
						dctParameters.dbg_mode,
						threadsMax,  // maximal number of threads to launch
						debugLevel);
		  }
		return dct_data;
	}

	public double [][][] lapped_dct(
			final double [] dpixels,
			final int       width,
			final int       dct_size,
			final int       dct_mode,    // 0: dct/dct, 1: dct/dst, 2: dst/dct, 3: dst/dst
			final int       window_type,
			final int       color,
			final EyesisDCT.DCTKernels dct_kernels,
			final boolean   skip_sym,
			final boolean   convolve_direct, // test feature - convolve directly with the symmetrical kernel
			final int       debug_tileX,
			final int       debug_tileY,
			final int       debug_mode,
			final int       threadsMax,  // maximal number of threads to launch
			final int       globalDebugLevel)
	{
		final int kernel_margin = 1; //move to parameters?
		final int height=dpixels.length/width;
		final int tilesX=width/dct_size-1;
		final int tilesY=height/dct_size-1;
		final int nTiles=tilesX*tilesY;
		final double [][][] dct_data = new double[tilesY][tilesX][dct_size*dct_size];
		final Thread[] threads = newThreadArray(threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		for (int tileY = 0; tileY < tilesY; tileY++){
			for (int tileX = 0; tileX < tilesX; tileX++){
				for (int i=0; i<dct_data[tileY][tileX].length;i++) dct_data[tileY][tileX][i]= 0.0; // actually not needed, Java initializes arrays
			}
		}
		double [] dc = new double [dct_size*dct_size];
		for (int i = 0; i<dc.length; i++) dc[i] = 1.0;
		DttRad2 dtt0 = new DttRad2(dct_size);
		dtt0.set_window(window_type);
		final double [] dciii = dtt0.dttt_iii  (dc, dct_size);
		final double [] dciiie = dtt0.dttt_iiie  (dc, 0, dct_size);
		if ((globalDebugLevel > 0) && (color ==2)) {
			double [][]dcx = {dc,dciii,dciiie, dtt0.dttt_ii(dc, dct_size),dtt0.dttt_iie(dc, 0, dct_size)};
			ShowDoubleFloatArrays sdfa_instance0 = new ShowDoubleFloatArrays(); // just for debugging?
			sdfa_instance0.showArrays(dcx,  dct_size, dct_size, true, "dcx");
		}


		if (globalDebugLevel > 0) {
			System.out.println("lapped_dctdc(): width="+width+" height="+height);
		}

		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				@Override
				public void run() {
					DttRad2 dtt = new DttRad2(dct_size);
					dtt.set_window(window_type);
					double [] tile_in = new double[4*dct_size * dct_size];
					double [] sym_conv= null;
					if ((dct_kernels != null) && convolve_direct){ // debug feature - directly convolve with symmetrical kernel
						sym_conv = new double[4*dct_size * dct_size];
					}
					double [] tile_folded;
					double [] tile_out; // = new double[dct_size * dct_size];
					int tileY,tileX;
					int n2 = dct_size * 2;
					double [] tile_out_copy = null;
					ShowDoubleFloatArrays sdfa_instance = new ShowDoubleFloatArrays(); // just for debugging?
					for (int nTile = ai.getAndIncrement(); nTile < nTiles; nTile = ai.getAndIncrement()) {
						tileY = nTile/tilesX;
						tileX = nTile - tileY * tilesX;
						int kernelTileY=0;
						int kernelTileX=0;
						//readDCTKernels() debugLevel = 1 kernels[0].size = 8 kernels[0].img_step = 16 kernels[0].asym_nonzero = 4 nColors = 3 numVert = 123 numHor =  164
						if (dct_kernels != null){ // convolve directly with asym_kernel
							int asym_center = dct_kernels.asym_size/2; // 7 for 15
							kernelTileY = kernel_margin + (tileY * dct_size) / dct_kernels.img_step;
							kernelTileX = kernel_margin + (tileX * dct_size) / dct_kernels.img_step;
							if ((globalDebugLevel > 0) && (tileY == debug_tileY) && (tileX == debug_tileX) && (color == 2)) {
								System.out.println("kernelTileY="+kernelTileY+" kernelTileX="+kernelTileX+" width="+width);
							}
							for (int i = 0; i < n2; i++){
								for (int j = 0; j < n2; j++){
									tile_in[i*n2 + j] = 0.0;
									// convolve list
									int [] asym_indx =   dct_kernels.asym_indx[color][kernelTileY][kernelTileX];
									double [] asym_val = dct_kernels.asym_val[color][kernelTileY][kernelTileX];
									for (int indx = 0; indx < asym_indx.length; indx++){
										int xy = asym_indx[indx];
										if ((globalDebugLevel > 0) && (tileY == debug_tileY) && (tileX == debug_tileX) && (color == 2)) {
											System.out.println("i="+i+" j="+j+" indx="+indx+" xy="+xy);
										}
										if (xy >= 0) {
											int dy = (xy / dct_kernels.asym_size) - asym_center;
											int dx = (xy % dct_kernels.asym_size) - asym_center;
											int y = tileY*dct_size - dy + i;
											int x = tileX*dct_size - dx + j;
											if (y < 0) y &= 1;
											if (x < 0) x &= 1;
											if (y >= height) y = (height - 2) + (y & 1);
											if (x >= width)  x = (width - 2) +  (x & 1);
											tile_in[i*n2 + j] += asym_val[indx] * dpixels[ y * width + x];
											if ((globalDebugLevel > 0) && (tileY == debug_tileY) && (tileX == debug_tileX) && (color == 2)) {
												System.out.println("dy= "+dy+" dx="+dx+" x = "+x+" y="+y+" y*width + x="+(y*width + x));
												System.out.println("asym_val["+indx+"]="+asym_val[indx]+
														"  dpixels["+(y * width + x)+"]="+ dpixels[ y * width + x]+
														"tile_in["+(i*n2 + j)+"]="+tile_in[i*n2 + j]);
											}
										}
									}
								}
							}
							// directly convolve with symmetrical kernel (debug feature
							if ((dct_kernels != null) && convolve_direct){
								double [] dir_sym = dct_kernels.st_direct[color][kernelTileY][kernelTileX];
								double s0 = 0;
								for (int i = 0; i < n2; i++){
									for (int j = 0; j < n2; j++){
										int indx = i*n2+j;
										sym_conv[indx] = 0.0; // dir_sym[0]* tile_in[indx];
										for (int dy = -dct_size +1; dy < dct_size; dy++){
											int ady = (dy>=0)?dy:-dy;
											int sgny = 1;
											int y = i - dy;
											if (y < 0){
												y = -1 -y;
												sgny = -sgny;
											}
											if (y >= n2){
												y = 2*n2 - y -1;
												sgny = -sgny;
											}
											for (int dx = -dct_size +1; dx < dct_size; dx++){
												int adx = (dx >= 0)? dx:-dx;
												int sgn = sgny;
												int x = j - dx;
												if (x < 0){
													x = -1 -x;
													sgn = -sgn;
												}
												if (x >= n2){
													x = 2*n2 - x -1;
													sgn = -sgn;
												}
												sym_conv[indx] += sgn*dir_sym[ady * dct_size + adx] * tile_in[y * n2 + x];
												s0+=dir_sym[ady * dct_size + adx];
												if ((globalDebugLevel > 0) && (tileY == debug_tileY) && (tileX == debug_tileX) && (color == 2) &&
														(i == dct_size) && (j== dct_size)) {
													System.out.println("i="+i+" j="+j+" dy="+dy+" dx="+dx+" ady="+ady+" adx="+adx+
															" y="+y+" x="+x+" sgny="+sgny+" sgn="+sgn+
															"sym_conv["+indx+"] += "+sgn+"* dir_sym["+(ady * dct_size + adx)+"] * tile_in["+(y * n2 + x)+"] +="+
															sgn+"* "+ dir_sym[ady * dct_size + adx]+" * "+tile_in[y * n2 + x]+" +="+
															(sgn*dir_sym[ady * dct_size + adx] * tile_in[y * n2 + x])+" ="+sym_conv[indx]);
												}

											}
										}
									}
								}


								if ((globalDebugLevel > 0) && (tileY == debug_tileY) && (tileX == debug_tileX) && (color == 2)) {
									//								if ((tileY == debug_tileY) && (tileX == debug_tileX)) {
									double [][] pair = {tile_in, sym_conv};
									sdfa_instance.showArrays(pair,  n2, n2, true, "dconv-X"+tileX+"Y"+tileY+"C"+color);
									sdfa_instance.showArrays(dir_sym,  dct_size, dct_size, "dk-X"+tileX+"Y"+tileY+"C"+color);
									double s1=0,s2=0;
									for (int i = 0; i<tile_in.length; i++){
										s1 +=tile_in[i];
										s2 +=sym_conv[i];
									}
									double s3 = 0.0;
									for (int i=0; i<dct_size;i++){
										for (int j=0; j<dct_size;j++){
											double d = dir_sym[i*dct_size+j];
											if (i > 0) d*=2;
											if (j > 0) d*=2;
											s3+=d;
										}
									}
									System.out.println("s1="+s1+" s2="+s2+" s1/s2="+(s1/s2)+" s0="+s0+" s3="+s3);
								}
//								tile_in = sym_conv.clone();
								System.arraycopy(sym_conv, 0, tile_in, 0, n2*n2);
							}
						} else { // no aberration correction, just copy data
							for (int i = 0; i < n2;i++){
								System.arraycopy(dpixels, (tileY*width+tileX)*dct_size + i*width, tile_in, i*n2, n2);
							}
						}
						tile_folded=dtt.fold_tile(tile_in, dct_size, 0); // DCCT
						tile_out=dtt.dttt_iv  (tile_folded, dct_mode, dct_size);
						if ((dct_kernels != null) && !skip_sym){ // convolve in frequency domain with sym_kernel
							double s0 =0;

							if (debug_mode == 2){
								for (int i=0;i<dct_kernels.st_kernels[color][kernelTileY][kernelTileX].length; i++){
									s0+=dct_kernels.st_kernels[color][kernelTileY][kernelTileX][i];
								}
								s0 = dct_size*dct_size/s0;
							} else if (debug_mode == 3){
								for (int i=0;i<dct_size;i++){
									double scale0 = (i>0)?2.0:1.0;
									for (int j=0;j<dct_size;j++){
										double scale = scale0*((j>0)?2.0:1.0);
										int indx = i*dct_size+j;
										s0+=scale*dct_kernels.st_kernels[color][kernelTileY][kernelTileX][indx];
									}
								}
								s0 = (2*dct_size-1)*(2*dct_size-1)/s0;
							}else if (debug_mode == 4){
								//dciii
								for (int i=0;i<dct_kernels.st_kernels[color][kernelTileY][kernelTileX].length; i++){
									s0+=dciii[i]* dct_kernels.st_kernels[color][kernelTileY][kernelTileX][i];
								}
								s0 = dct_size*dct_size/s0;
							} else s0 = 1.0;

							for (int i = 0; i < tile_out.length; i++){
								tile_out[i] *= s0;
							}
						}

						if ((tileY == debug_tileY) && (tileX == debug_tileX) && (color == 2)) {
							tile_out_copy = tile_out.clone();
						}


						if ((dct_kernels != null) && !skip_sym){ // convolve in frequency domain with sym_kernel
							for (int i = 0; i < tile_out.length; i++){
								tile_out[i] *=dct_kernels.st_kernels[color][kernelTileY][kernelTileX][i];
							}
						}


						if ((dct_kernels!=null) && (tileY == debug_tileY) && (tileX == debug_tileX) && (color == 2)) {
							double [][] dbg_tile = {
									dct_kernels.st_direct[color][kernelTileY][kernelTileX],
									dct_kernels.st_kernels[color][kernelTileY][kernelTileX],
									tile_out_copy,
									tile_out};
							if (globalDebugLevel > 0){
								sdfa_instance.showArrays(tile_in,  n2, n2, "tile_in-X"+tileX+"Y"+tileY+"C"+color);
								sdfa_instance.showArrays(dbg_tile,  dct_size, dct_size, true, "dbg-X"+tileX+"Y"+tileY+"C"+color);
								System.out.println("tileY="+tileY+" tileX="+tileX+" kernelTileY="+kernelTileY+" kernelTileX="+kernelTileX);
								double s0=0.0, s1=0.0, s2=0.0, s3=0.0;
								for (int i=0;i<dct_size;i++){
									double scale0 = (i>0)?2.0:1.0;
									for (int j=0;j<dct_size;j++){
										double scale = scale0*((j>0)?2.0:1.0);
										int indx = i*dct_size+j;
										s0+=scale*dct_kernels.st_direct[color][kernelTileY][kernelTileX][indx];
										s1+=scale*dct_kernels.st_kernels[color][kernelTileY][kernelTileX][indx];
										s2+=      dct_kernels.st_kernels[color][kernelTileY][kernelTileX][indx];
										s3+=dciii[indx]*dct_kernels.st_kernels[color][kernelTileY][kernelTileX][indx];
									}
								}
								System.out.println("s0="+s0+" s1="+s1+" s2="+s2+" s3="+s3);
							}
						}
						System.arraycopy(tile_out, 0, dct_data[tileY][tileX], 0, tile_out.length);
					}
				}
			};
		}
		startAndJoin(threads);
		return dct_data;
	}

	// extract DCT transformed parameters in linescan order (for visualization)
	public double [] lapped_dct_dbg(
			final double [][][] dct_data,
			final int           threadsMax,     // maximal number of threads to launch
			final int           globalDebugLevel)
	{
		final int tilesY=dct_data.length;
		final int tilesX=dct_data[0].length;
		final int nTiles=tilesX*tilesY;
		final int dct_size = (int) Math.round(Math.sqrt(dct_data[0][0].length));
		final int dct_len = dct_size*dct_size;
		final double [] dct_data_out = new double[tilesY*tilesX*dct_len];
		final Thread[] threads = newThreadArray(threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		for (int i=0; i<dct_data_out.length;i++) dct_data_out[i]= 0;

		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				@Override
				public void run() {
					int tileY,tileX;
					for (int nTile = ai.getAndIncrement(); nTile < nTiles; nTile = ai.getAndIncrement()) {
						tileY = nTile/tilesX;
						tileX = nTile - tileY * tilesX;
						for (int i = 0; i < dct_size;i++){
							System.arraycopy(dct_data[tileY][tileX], dct_size* i, dct_data_out, ((tileY*dct_size + i) *tilesX + tileX)*dct_size , dct_size);
						}
					}
				}
			};
		}
		startAndJoin(threads);
		return dct_data_out;
	}

	public void dct_lpf(
			final double sigma,
			final double [][][] dct_data,
			final int       threadsMax,     // maximal number of threads to launch
			final int       globalDebugLevel)
	{
		final int tilesY=dct_data.length;
		final int tilesX=dct_data[0].length;
		final int nTiles=tilesX*tilesY;
		final int dct_size = (int) Math.round(Math.sqrt(dct_data[0][0].length));
		final int dct_len = dct_size*dct_size;
		final double [] filter_direct= new double[dct_len];
		if (sigma == 0) {
			filter_direct[0] = 1.0;
			for (int i= 1; i<filter_direct.length;i++) filter_direct[i] =0;
		} else {
			for (int i = 0; i < dct_size; i++){
				for (int j = 0; j < dct_size; j++){
					filter_direct[i*dct_size+j] = Math.exp(-(i*i+j*j)/(2*sigma));
				}
			}
		}
		// normalize
		double sum = 0;
		for (int i = 0; i < dct_size; i++){
			for (int j = 0; j < dct_size; j++){
				double d = 	filter_direct[i*dct_size+j];
				d*=Math.cos(Math.PI*i/(2*dct_size))*Math.cos(Math.PI*j/(2*dct_size));
				if (i > 0) d*= 2.0;
				if (j > 0) d*= 2.0;
				sum +=d;
			}
		}
		for (int i = 0; i<filter_direct.length; i++){
			filter_direct[i] /= sum;
		}

		if (globalDebugLevel > 0) {
			for (int i=0; i<filter_direct.length;i++){
				System.out.println("dct_lpf_psf() "+i+": "+filter_direct[i]);
			}
		}
		DttRad2 dtt = new DttRad2(dct_size);
		final double [] filter= dtt.dttt_iiie(filter_direct);
		final double [] dbg_filter= dtt.dttt_ii(filter);

//		for (int i=0; i < filter.length;i++) filter[i] *= dct_size;
		for (int i=0; i < filter.length;i++) filter[i] *= 2*dct_size;

		if (globalDebugLevel > 0) {
			for (int i=0; i<filter.length;i++){
				System.out.println("dct_lpf_psf() "+i+": "+filter[i]);
			}
			ShowDoubleFloatArrays sdfa_instance = new ShowDoubleFloatArrays(); // just for debugging?
			double [][] ff = {filter_direct,filter,dbg_filter};
			sdfa_instance.showArrays(ff,  dct_size,dct_size, true, "filter_lpf");
		}

		final Thread[] threads = newThreadArray(threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);

		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				@Override
				public void run() {
					int tileY,tileX;
					for (int nTile = ai.getAndIncrement(); nTile < nTiles; nTile = ai.getAndIncrement()) {
						tileY = nTile/tilesX;
						tileX = nTile - tileY * tilesX;
						for (int i = 0; i < filter.length; i++){
							dct_data[tileY][tileX][i] *= filter[i];
						}
					}
				}
			};
		}
		startAndJoin(threads);
	}

	public double [][][][] dct_color_convert(
			final double [][][][] dct_data,
			final double kr,
			final double kb,
			final double sigma_rb,        // blur of channels 0,1 (r,b) in addition to 2 (g)
			final double sigma_y,         // blur of Y from G
			final double sigma_color,     // blur of Pr, Pb in addition to Y
			final int       threadsMax,     // maximal number of threads to launch
			final int       globalDebugLevel)
	{
		final int tilesY=dct_data[0].length;
		final int tilesX=dct_data[0][0].length;
		final int nTiles=tilesX*tilesY;
		final int dct_size = (int) Math.round(Math.sqrt(dct_data[0][0][0].length));
		final int dct_len = dct_size*dct_size;
		final double [][][][] yPrPb = new double [3][tilesY][tilesX][dct_len];
		final double [][][] filters = new double [3][3][dct_len];
		final double kg = 1.0 - kr - kb;
		final double [][] filters_proto_direct = new double[3][dct_len];
		final double [][] filters_proto = new double[3][];
		System.out.println("dct_color_convert(): kr="+kr+" kg="+kg+" kb="+kb);
		final double [] sigmas = {sigma_rb,sigma_y,sigma_color};
		double [] norm_sym_weights = new double [dct_size*dct_size];
		for (int i = 0; i < dct_size; i++){
			for (int j = 0; j < dct_size; j++){
				double d = 	Math.cos(Math.PI*i/(2*dct_size))*Math.cos(Math.PI*j/(2*dct_size));
				if (i > 0) d*= 2.0;
				if (j > 0) d*= 2.0;
				norm_sym_weights[i*dct_size+j] = d;
			}
		}

		for (int n = 0; n<3; n++) {

			double s = 0.0;
			for (int i = 0; i < dct_size; i++){
				for (int j = 0; j < dct_size; j++){
					double d;
					if (sigmas[n] == 0.0)   d = ((i == 0) && (j==0))? 1.0:0.0;
					else                    d = Math.exp(-(i*i+j*j)/(2*sigmas[n]));
					filters_proto_direct[n][i*dct_size+j] = d;
				}

			}
			for (int i = 0; i< dct_len; i++){
				s += norm_sym_weights[i]*filters_proto_direct[n][i];
			}

			if (globalDebugLevel>0) System.out.println("dct_color_convert(): sigmas["+n+"]="+sigmas[n]+", sum="+s);
			for (int i = 0; i < dct_len; i++){
				filters_proto_direct[n][i] /=s;
			}
		}

		DttRad2 dtt = new DttRad2(dct_size);
		for (int i = 0; i < filters_proto.length; i++){
			filters_proto[i] = dtt.dttt_iiie(filters_proto_direct[i]);
			if (globalDebugLevel > 0)  System.out.println("filters_proto.length="+filters_proto.length+" filters_proto["+i+"].length="+filters_proto[i].length+" dct_len="+dct_len+" dct_size="+dct_size);
			for (int j=0; j < dct_len; j++) filters_proto[i][j] *= 2*dct_size;

		}
		if (globalDebugLevel > 0) {
			ShowDoubleFloatArrays sdfa_instance = new ShowDoubleFloatArrays(); // just for debugging?
			double [][] ff = {filters_proto_direct[0],filters_proto_direct[1],filters_proto_direct[2],filters_proto[0],filters_proto[1],filters_proto[2]};
			sdfa_instance.showArrays(ff,  dct_size,dct_size, true, "filters_proto");
		}

		double [][] coeff_arr ={
				{ kr,              kb,               kg             },  // Y = R*Kr+G*Hg+B*Kb
				{ 0.5,            -kb/(2.0*(1-kr)), -kg/(2.0*(1-kr))},  // Pr =  R* 0.5  - G* Kg/(2.0*(1-Kr)) - B *Kb/(2.0*(1-Kr))
				{-kr/(2.0*(1-kb)), 0.5,             -kg/(2.0*(1-kb))}}; // Pb =  B* 0.5  - G* Kg/(2.0*(1-Kb)) - R *Kr/(2.0*(1-Kb))
		for (int k = 0; k < dct_len; k++){
			for (int i = 0; i < coeff_arr.length; i++){
				for (int j = 0; j < coeff_arr.length; j++){
					filters[i][j][k] = coeff_arr[i][j]* filters_proto[1][k];      // minimal blur - for all sigma_y
					if (i > 0){
						filters[i][j][k] *= filters_proto[2][k]; // for Pr, Pb sigma_color
					}
					if (j <2){ // all but green
						filters[i][j][k] *= filters_proto[0][k]; // for R,B sigma_rb
					}

				}
			}
		}
		if (globalDebugLevel > 0) {
			ShowDoubleFloatArrays sdfa_instance = new ShowDoubleFloatArrays(); // just for debugging?
			double [][] ff = {
					filters[0][0], filters[0][1], filters[0][2],
					filters[1][0], filters[1][1], filters[1][2],
					filters[2][0], filters[2][1], filters[2][2]};
			sdfa_instance.showArrays(ff,  dct_size,dct_size, true, "filters");
		}

		final Thread[] threads = newThreadArray(threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);

		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				@Override
				public void run() {
					int tileY,tileX;
					for (int nTile = ai.getAndIncrement(); nTile < nTiles; nTile = ai.getAndIncrement()) {
						tileY = nTile/tilesX;
						tileX = nTile - tileY * tilesX;
						for (int i = 0; i < filters.length; i++){
							for (int k = 0; k <dct_len; k++){
								yPrPb[i][tileY][tileX][k]=0.0;
								for (int j = 0; j < filters[i].length; j++){
									yPrPb[i][tileY][tileX][k] += filters[i][j][k] * dct_data[j][tileY][tileX][k];
								}
							}

						}
					}
				}
			};
		}
		startAndJoin(threads);
		return yPrPb;
	}






	public double [] lapped_idct(
//			final double [][][] dctdc_data,  // array [tilesY][tilesX][dct_size*dct_size+1] - last element is DC value
			final double [][][] dct_data,  // array [tilesY][tilesX][dct_size*dct_size]
			final int       dct_size,
			final int       window_type,
			final int       threadsMax,  // maximal number of threads to launch
			final int       globalDebugLevel)
	{
//		final int tilesX=dct_width/dct_size;
//		final int tilesY=dct_data.length/(dct_width*dct_size);
		final int tilesY=dct_data.length;
		final int tilesX=dct_data[0].length;

		final int width=  (tilesX+1)*dct_size;
		final int height= (tilesY+1)*dct_size;
		if (globalDebugLevel > 0) {
			System.out.println("lapped_idct():tilesX=   "+tilesX);
			System.out.println("lapped_idct():tilesY=   "+tilesY);
			System.out.println("lapped_idct():width=    "+width);
			System.out.println("lapped_idct():height=   "+height);
		}
		final double [] dpixels = new double[width*height];
		final Thread[] threads = newThreadArray(threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		final AtomicInteger nser = new AtomicInteger(0);
		final int [][][] tiles_list = new int[4][][];
		for (int n=0; n<4; n++){
			int nx = (tilesX + 1 - (n &1)) / 2;
			int ny = (tilesY + 1 - ((n>>1) & 1)) / 2;
			tiles_list[n] = new int [nx*ny][2];
			int indx = 0;
			for (int i = 0;i < ny; i++) for (int j = 0; j < nx; j++){
				tiles_list[n][indx][0]=2*j+(n &1);
				tiles_list[n][indx++][1]=2*i+((n>>1) & 1);
			}
		}
		for (int i=0; i<dpixels.length;i++) dpixels[i]= 0;
		for (int n=0; n<4; n++){
			nser.set(n);
			ai.set(0);
			for (int ithread = 0; ithread < threads.length; ithread++) {
				threads[ithread] = new Thread() {
					@Override
					public void run() {
						DttRad2 dtt = new DttRad2(dct_size);
						dtt.set_window(window_type);
						double [] tile_in = new double[dct_size * dct_size];
						double [] tile_dct; // = new double[dct_size * dct_size];
						double [] tile_out; //  = new double[4*dct_size * dct_size];
						int tileY,tileX;
						int n2 = dct_size * 2;
						for (int nTile = ai.getAndIncrement(); nTile < tiles_list[nser.get()].length; nTile = ai.getAndIncrement()) {
							tileX = tiles_list[nser.get()][nTile][0];
							tileY = tiles_list[nser.get()][nTile][1];
							System.arraycopy(dct_data[tileY][tileX], 0, tile_in, 0, tile_in.length);
							tile_dct=dtt.dttt_iv  (tile_in, 0, dct_size);
							tile_out=dtt.unfold_tile(tile_dct, dct_size, 0); // mpode=0 - DCCT
							for (int i = 0; i < n2;i++){
								int start_line = ((tileY*dct_size + i) *(tilesX+1) + tileX)*dct_size;
								for (int j = 0; j<n2;j++) {
									dpixels[start_line + j] += tile_out[n2 * i + j]; //  +1.0;
								}
							}
						}
					}
				};
			}
			startAndJoin(threads);
		}
		return dpixels;
	}

	// perform 2d clt and apply aberration corrections, all colors
	public double [][][][][] clt_aberrations(
			final double [][]       image_data,
			final int               width,
			final double [][][][][] clt_kernels, // [color][tileY][tileX][band][pixel] , size should match image (have 1 tile around)
			final int               kernel_step,
			final int               transform_size,
			final int               window_type,
			final double            shiftX, // shift image horizontally (positive - right) - just for testing
			final double            shiftY, // shift image vertically (positive - down)
			final int               debug_tileX,
			final int               debug_tileY,
			final boolean           no_fract_shift,
			final boolean           no_deconvolution,
			final boolean           transpose,
			final int               threadsMax,  // maximal number of threads to launch
			final int               globalDebugLevel)
	{
		final int nChn = image_data.length;
		final int height=image_data[0].length/width;
		final int tilesX=width/transform_size;
		final int tilesY=height/transform_size;
		final int nTilesInChn=tilesX*tilesY;
		final int nTiles=tilesX*tilesY*nChn;
		final double [][][][][] clt_data = new double[nChn][tilesY][tilesX][4][];
		final Thread[] threads = newThreadArray(threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		if (globalDebugLevel > 0) {
			System.out.println("clt_aberrations(): width="+width+" height="+height+" transform_size="+transform_size+
					" debug_tileX="+debug_tileX+" debug_tileY="+debug_tileY+" globalDebugLevel="+globalDebugLevel);
		}
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				@Override
				public void run() {
					DttRad2 dtt = new DttRad2(transform_size);
					dtt.set_window(window_type);
					int tileY,tileX, chn;
					//						showDoubleFloatArrays sdfa_instance = new showDoubleFloatArrays(); // just for debugging?
					double centerX; // center of aberration-corrected (common model) tile, X
					double centerY; //
					for (int nTile = ai.getAndIncrement(); nTile < nTiles; nTile = ai.getAndIncrement()) {
						chn=nTile/nTilesInChn;
						tileY =(nTile % nTilesInChn)/tilesX;
						tileX = nTile % tilesX;
//						centerX = tileX * transform_size - transform_size/2 - shiftX;
//						centerY = tileY * transform_size - transform_size/2 - shiftY;
						centerX = tileX * transform_size + transform_size/2 - shiftX;
						centerY = tileY * transform_size + transform_size/2 - shiftY;

						double [] fract_shiftXY = extract_correct_tile( // return a pair of resudual offsets
								image_data,
								width,       // image width
								clt_kernels, // [color][tileY][tileX][band][pixel]
								clt_data[chn][tileY][tileX], //double  [][]        clt_tile,    // should be double [4][];
								kernel_step,
								transform_size,
								dtt,
								chn,
								centerX, // center of aberration-corrected (common model) tile, X
								centerY, //
								((globalDebugLevel > 0) && (tileX == debug_tileX) && (tileY == debug_tileY) && (chn == 2)) ? 1 : 0, // external tile compare
								no_deconvolution,
								transpose,
								// no saturation processing
								null, // boolean []          saturation_imp, // (near) saturated pixels or null
								null); // int []              overexp_all ) // {number of overexposed,  number of all tiles} or null

						if ((globalDebugLevel > 0) && (debug_tileX == tileX) && (debug_tileY == tileY)  && (chn == 2)) {
							ShowDoubleFloatArrays sdfa_instance = new ShowDoubleFloatArrays(); // just for debugging?
							String [] titles = {"CC","SC","CS","SS"};
							sdfa_instance.showArrays(clt_data[chn][tileY][tileX],  transform_size, transform_size, true, "pre-shifted_x"+tileX+"_y"+tileY, titles);
						}

						if ((globalDebugLevel > -1) && (tileX >= debug_tileX - 2) && (tileX <= debug_tileX + 2) &&
								(tileY >= debug_tileY - 2) && (tileY <= debug_tileY+2)) {
							System.out.println("clt_aberrations(): color="+chn+", tileX="+tileX+", tileY="+tileY+
									" fract_shiftXY[0]="+fract_shiftXY[0]+" fract_shiftXY[1]="+fract_shiftXY[1]);
						}

						if (!no_fract_shift) {
							// apply residual shift
							fract_shift(    // fractional shift in transform domain. Currently uses sin/cos - change to tables with 2? rotations
									clt_data[chn][tileY][tileX], // double  [][]  clt_tile,
									transform_size,
									fract_shiftXY[0],            // double        shiftX,
									fract_shiftXY[1],            // double        shiftY,
//									(globalDebugLevel > 0) && (tileX == debug_tileX) && (tileY == debug_tileY)); // external tile compare
									((globalDebugLevel > 0) && (chn==0) && (tileX >= debug_tileX - 2) && (tileX <= debug_tileX + 2) &&
											(tileY >= debug_tileY - 2) && (tileY <= debug_tileY+2)));
							if ((globalDebugLevel > 0) && (debug_tileX == tileX) && (debug_tileY == tileY)) {
								ShowDoubleFloatArrays sdfa_instance = new ShowDoubleFloatArrays(); // just for debugging?
								String [] titles = {"CC","SC","CS","SS"};
								sdfa_instance.showArrays(clt_data[chn][tileY][tileX],  transform_size, transform_size, true, "shifted_x"+tileX+"_y"+tileY, titles);
							}
						}
					}
				}
			};
		}
		startAndJoin(threads);
		return clt_data;
	}


	public double [][][][][][] clt_aberrations_quad(
			final double              disparity,
			final double [][][]       image_data, // first index - number of image in a quad
			final int                 width,
			final GeometryCorrection  geometryCorrection,
			final double [][][][][][] clt_kernels, // [channel_in_quad][color][tileY][tileX][band][pixel] , size should match image (have 1 tile around)
			final int                 kernel_step,
			final int                 transform_size,
			final int                 window_type,
			final double              shiftX, // shift image horizontally (positive - right) - just for testing
			final double              shiftY, // shift image vertically (positive - down)
			final int                 debug_tileX,
			final int                 debug_tileY,
			final boolean             no_fract_shift,
			final boolean             no_deconvolution,
			final boolean             transpose,
			final int                 threadsMax,  // maximal number of threads to launch
			final int                 globalDebugLevel)
	{
		final int quad = 4;   // number of subcameras
		final int numcol = 3; // number of colors
		final int nChn = image_data[0].length;
		final int height=image_data[0][0].length/width;
		final int tilesX=width/transform_size;
		final int tilesY=height/transform_size;
		final int nTilesInChn=tilesX*tilesY;
//		final int nTiles=tilesX*tilesY*nChn;
		final double [][][][][][] clt_data = new double[quad][nChn][tilesY][tilesX][4][];
		final Thread[] threads = newThreadArray(threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		final Matrix [] corr_rots = geometryCorrection.getCorrVector().getRotMatrices(); // get array of per-sensor rotation matrices


		if (globalDebugLevel > 0) {
			System.out.println("clt_aberrations(): width="+width+" height="+height+" transform_size="+transform_size+
					" debug_tileX="+debug_tileX+" debug_tileY="+debug_tileY+" globalDebugLevel="+globalDebugLevel);
		}
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				@Override
				public void run() {
					DttRad2 dtt = new DttRad2(transform_size);
					dtt.set_window(window_type);
					int tileY,tileX; // , chn;
					//						showDoubleFloatArrays sdfa_instance = new showDoubleFloatArrays(); // just for debugging?
					double centerX; // center of aberration-corrected (common model) tile, X
					double centerY; //
					double [][] fract_shiftsXY = new double[quad][];

//					for (int nTile = ai.getAndIncrement(); nTile < nTiles; nTile = ai.getAndIncrement()) {
					for (int nTile = ai.getAndIncrement(); nTile < nTilesInChn; nTile = ai.getAndIncrement()) {
						// TODO: make all color channels to be processed here (atomically)
						//						chn=nTile/nTilesInChn;
						//						tileY =(nTile % nTilesInChn)/tilesX;
						//						tileX = nTile % tilesX;
						tileY = nTile /tilesX;
						tileX = nTile % tilesX;



						centerX = tileX * transform_size + transform_size/2 - shiftX;
						centerY = tileY * transform_size + transform_size/2 - shiftY;
//						double [][] centersXY = geometryCorrection.getPortsCoordinates(
//								centerX,
//								centerY,
//								disparity);
						double [][] centersXY = geometryCorrection.getPortsCoordinatesAndDerivatives(
								geometryCorrection, //			GeometryCorrection gc_main,
								false,          // boolean use_rig_offsets,
								corr_rots, // Matrix []   rots,
								null,      //  Matrix [][] deriv_rots,
								null,      // double [][] pXYderiv, // if not null, should be double[8][]
								centerX,
								centerY,
								disparity);

						if ((globalDebugLevel > 0) && (tileX >= debug_tileX - 2) && (tileX <= debug_tileX + 2) &&
								(tileY >= debug_tileY - 2) && (tileY <= debug_tileY+2)) {
							for (int i = 0; i < quad; i++) {
								System.out.println("clt_aberrations_quad():  tileX="+tileX+", tileY="+tileY+
										" centerX="+centerX+" centerY="+centerY+" disparity="+disparity+
										" centersXY["+i+"][0]="+centersXY[i][0]+" centersXY["+i+"][1]="+centersXY[i][1]);
							}
						}


						for (int chn = 0; chn <numcol; chn++) {

							for (int i = 0; i < quad; i++) {
								fract_shiftsXY[i] = extract_correct_tile( // return a pair of resudual offsets
										image_data[i],
										width,       // image width
										clt_kernels[i], // [color][tileY][tileX][band][pixel]
										clt_data[i][chn][tileY][tileX], //double  [][]        clt_tile,    // should be double [4][];
										kernel_step,
										transform_size,
										dtt,
										chn,
										centersXY[i][0], // centerX, // center of aberration-corrected (common model) tile, X
										centersXY[i][1], // centerY, //
										((globalDebugLevel > 0) && (tileX == debug_tileX) && (tileY == debug_tileY) && (chn == 2)) ? 1 : 0, // external tile compare
										no_deconvolution,
										transpose,
										// no saturation processing
										null, // boolean []          saturation_imp, // (near) saturated pixels or null
										null); // int []              overexp_all ) // {number of overexposed,  number of all tiles} or null
							}
							if ((globalDebugLevel > 0) && (debug_tileX == tileX) && (debug_tileY == tileY)  && (chn == 2)) {
								ShowDoubleFloatArrays sdfa_instance = new ShowDoubleFloatArrays(); // just for debugging?
								String [] titles = {"CC0","SC0","CS0","SS0","CC1","SC1","CS1","SS1","CC2","SC2","CS2","SS2","CC3","SC3","CS3","SS3"};
								double [][] dbg_tile = new double [16][];
								for (int i = 0; i < 16; i++) dbg_tile[i]=clt_data[i>>2][chn][tileY][tileX][i & 3];
								sdfa_instance.showArrays(dbg_tile,  transform_size, transform_size, true, "pre-shifted_x"+tileX+"_y"+tileY, titles);
							}

							if ((globalDebugLevel > 0) && (tileX >= debug_tileX - 2) && (tileX <= debug_tileX + 2) &&
									(tileY >= debug_tileY - 2) && (tileY <= debug_tileY+2)) {
								for (int i = 0; i < quad; i++) {
									System.out.println("clt_aberrations_quad(): color="+chn+", tileX="+tileX+", tileY="+tileY+
											" fract_shiftsXY["+i+"][0]="+fract_shiftsXY[i][0]+" fract_shiftsXY["+i+"][1]="+fract_shiftsXY[i][1]);
								}
							}

							if (!no_fract_shift) {
								// apply residual shift
								for (int i = 0; i < quad; i++) {
									fract_shift(    // fractional shift in transform domain. Currently uses sin/cos - change to tables with 2? rotations
											clt_data[i][chn][tileY][tileX], // double  [][]  clt_tile,
											transform_size,
											fract_shiftsXY[i][0],            // double        shiftX,
											fract_shiftsXY[i][1],            // double        shiftY,
											//									(globalDebugLevel > 0) && (tileX == debug_tileX) && (tileY == debug_tileY)); // external tile compare
											((globalDebugLevel > 0) && (chn==0) && (tileX >= debug_tileX - 2) && (tileX <= debug_tileX + 2) &&
													(tileY >= debug_tileY - 2) && (tileY <= debug_tileY+2)));
								}
								if ((globalDebugLevel > 0) && (debug_tileX == tileX) && (debug_tileY == tileY)) {
									ShowDoubleFloatArrays sdfa_instance = new ShowDoubleFloatArrays(); // just for debugging?
									String [] titles = {"CC0","SC0","CS0","SS0","CC1","SC1","CS1","SS1","CC2","SC2","CS2","SS2","CC3","SC3","CS3","SS3"};
									double [][] dbg_tile = new double [16][];
									for (int i = 0; i < 16; i++) dbg_tile[i]=clt_data[i>>2][chn][tileY][tileX][i & 3];
									sdfa_instance.showArrays(dbg_tile,  transform_size, transform_size, true, "shifted_x"+tileX+"_y"+tileY, titles);
								}
							}
						}
						// all color channels are done here
					}
				}
			};
		}
		startAndJoin(threads);
		return clt_data;
	}

	public void printSignsFPGA (
			DttRad2               dtt
			){
		double [][][] fold_coeff = dtt.getFoldK();
		int [][] fpga_fi = dtt.getFoldIndex();
// For R/B color channels (1 - 4 non-zero) show signs during folding of a single pixel contributor (4 Bayer variants) per each mode
		String [] mode_names={"CC","SC", "CS","SS"};
		int [] bayer_patterns = {0x1, 0x2, 0x4, 0x8}; // , 0x9, 0x6};
		boolean [][][] signs = new boolean [bayer_patterns.length][4][64];

//		for (int bp:bayer_patterns){
	    for (int ibp = 0; ibp < bayer_patterns.length; ibp++){
	    	int bp = bayer_patterns[ibp];
			System.out.println("\nPattern (row/col) "+bp+":");
			System.out.println("| "+(((bp & 1) !=0) ? "X ":"  ")+(((bp & 2) !=0) ? "X ":"  ")+"|");
			System.out.println("| "+(((bp & 4) !=0) ? "X ":"  ")+(((bp & 8) !=0) ? "X ":"  ")+"|");
			for (int mode = 0; mode < 4; mode++){
				if (mode == 0) 	System.out.println("DTT mode = "+mode+" ("+ mode_names[mode]+"): term sign");
				else 	        System.out.println("DTT mode = "+mode+" ("+ mode_names[mode]+"): term inverse relative to CC ");
				for (int i = 0; i < 64; i++){
					for (int k = 0; k < 4; k++){
						int row = (fpga_fi[i][k] >> 4);
						int col = (fpga_fi[i][k] & 0xf);
						int indx = (row & 1) + 2 * (col & 1);
						if (((1 << indx) & bp) != 0) { // only use non-zero pixels, for 1 in 4 - only one k would match
							signs[ibp][mode][i] = fold_coeff[mode][i][k] < 0;
							if (mode == 0) 	{
								if (fold_coeff[mode][i][k] < 0) System.out.print("- ");
								else                            System.out.print("+ ");
							} else {
								boolean sgn = signs[ibp][mode][i] ^ signs[ibp][0][i];
								if (sgn) System.out.print("* ");
								else     System.out.print(". ");
							}
//							continue;
						}
					}
					if ((i+1)%8 == 0) System.out.println();
				}
			}
		}
	}

	public void generateFPGACompareData(
			final double [][]         image_data, // for selected subcamera
			final double [][]         colorCentersXY, // pixel centers per color (2 - green)
			final int                 transform_size,
			final int                 width,
			DttRad2                   dtt
			){

		printSignsFPGA(dtt);

		int height = image_data[0].length/width;
		double [][][] fpga_clt_data_in = new double [3][4][];
		double [][][] fpga_clt_data_out = new double [3][4][];
		double [][][] fpga_clt_data_rot = new double [3][4][];
//		double [][] fpga_fract_shiftsXY = new double[3][];
		double [][] fpga_centersXY = new double [3][2];
//		int    [][] color_int_shifts = new int [3][2];

		double [][][][] fold_coeff = new double[3][][][];
		int [] ctile_left = new int [3];
		int [] ctile_top =  new int [3];
		double [][] residual_shift = new double[3][2];
		int    []   ishx = new int[3];
		int    []   ishy = new int[3];

		double [][] fpga_full_tile = new double [3][FPGA_TILE_SIZE * FPGA_TILE_SIZE];
		double [][] fpga_tile =      new double [3][4*transform_size*transform_size];
		for (int chn = 0; chn<3; chn++) for (int j = 0; j < 2; j++) {
			fpga_centersXY[chn][j] = colorCentersXY[chn][j];
			// Round to FPGA precision
			fpga_centersXY[chn][j] = Math.round(128*fpga_centersXY[chn][j])/128.0;
		}
		for (int chn = 0; chn<3; chn++) {
			double px = fpga_centersXY[chn][0] - transform_size;
			double py = fpga_centersXY[chn][1] - transform_size;
			// Was wrong rounding, fractional part gets to +0.5
			ctile_left[chn] = (int) -Math.round(-px);
			ctile_top[chn] =  (int) -Math.round(-py);
			residual_shift[chn][0] = -(px - ctile_left[chn]);
			residual_shift[chn][1] = -(py - ctile_top[chn]);
		}
		int lt = (FPGA_TILE_SIZE - 2 * transform_size)/2;
		for (int chn = 0; chn < 3; chn++){
			for (int i = 0; i < FPGA_TILE_SIZE; i++){
				System.arraycopy(
						image_data[chn],
						((ctile_top[GREEN_CHN] - lt) + i) * width + (ctile_left[GREEN_CHN] - lt),
						fpga_full_tile[chn], FPGA_TILE_SIZE * i,
						FPGA_TILE_SIZE);
			}
		}
		for (int chn = 0; chn < 3; chn++){

	 		if ((ctile_left[chn] >= 0) && (ctile_left[chn] < (width - transform_size * 2)) &&
					(ctile_top[chn] >= 0) && (ctile_top[chn] < (height - transform_size * 2))) {
				for (int i = 0; i < transform_size * 2; i++){
					System.arraycopy(image_data[chn], (ctile_top[chn] + i) * width + ctile_left[chn], fpga_tile[chn], transform_size * 2 * i, transform_size* 2);
				}
			} else { // copy by 1
				for (int i = 0; i < transform_size* 2; i++){
					int pi = ctile_top[chn] + i;
					if      (pi < 0)       pi &= 1;
					else if (pi >= height) pi = height - 2 + (pi & 1);
					for (int j = 0; j < transform_size* 2; j++){
						int pj = ctile_left[chn] + j;
						if      (pj < 0)      pj &= 1;
						else if (pj >= width) pj = width - 2 + (pj & 1);
						fpga_tile[chn][transform_size * 2 * i + j] = image_data[chn][pi * width + pj];
					}
				}
			}
		}
		// Fold and transform
		for (int chn = 0; chn < 3; chn++){
			fold_coeff[chn] = dtt.get_shifted_fold_2d ( // get_shifted_fold_2d(
					transform_size,
					residual_shift[chn][0],
					residual_shift[chn][1],
					0); // debug level
		}
		for (int chn = 0; chn < 3; chn++){
			for (int dct_mode = 0; dct_mode <4; dct_mode++) {
//				fpga_clt_data_in[chn][dct_mode] = dtt.fold_tile_debug (fpga_tile[chn], transform_size, dct_mode, fold_coeff[chn]); // DCCT, DSCT, DCST, DSST
				fpga_clt_data_in[chn][dct_mode] = dtt.fold_tile (fpga_tile[chn], transform_size, dct_mode, fold_coeff[chn]); // DCCT, DSCT, DCST, DSST
			}
		}

		for (int chn = 0; chn < 3; chn++) if (chn != GREEN_CHN) {
			for (int dct_mode = 0; dct_mode <4; dct_mode++) {
				for (int i = 0; i < 64; i++){
					fpga_clt_data_in[chn][dct_mode][i] *= 2.0; //adding twice each number in FPGA for R and B
				}
			}
		}


		for (int chn = 0; chn < 3; chn++) {
			for (int dct_mode = 0; dct_mode <4; dct_mode++) {
				fpga_clt_data_out[chn][dct_mode] = fpga_clt_data_in[chn][dct_mode].clone();
			}
		}

		double scale1 = (1 << (FPGA_DTT_IN - 9)); //  -1;
		scale1 *= 1.0 *((1 << FPGA_WND_BITS) -1) / (1 << FPGA_WND_BITS);
		scale1 *= 1.0 *((1 << FPGA_WND_BITS) -1) / (1 << FPGA_WND_BITS);
//		scale1 *= 2.0;
		System.out.println("scale1="+scale1);
		for (int chn = 0; chn < 3; chn++) {
			for (int dct_mode = 0; dct_mode <4; dct_mode++) {
				System.out.println("// Color="+chn+" fpga_clt_data_out[chn][dct_mode] = dtt.dttt_iv(..., scale1="+scale1);
				fpga_clt_data_out[chn][dct_mode] = dtt.dttt_iv   (fpga_clt_data_out[chn][dct_mode], dct_mode, transform_size, scale1, ((1 << 25) -1)); // debug level
//				fpga_clt_data_out[chn][dct_mode] = dtt.dttt_iv   (fpga_clt_data_out[chn][dct_mode], dct_mode, transform_size);
			}
		}

		for (int chn = 0; chn < 3; chn++) {
			for (int dct_mode = 0; dct_mode <4; dct_mode++) {
				fpga_clt_data_rot[chn][dct_mode] = fpga_clt_data_out[chn][dct_mode].clone();
			}
		}

// Rotate for fractional shift:

		for (int chn = 0; chn < 3; chn++) {
			fract_shift(    // fractional shift in transform domain. Currently uses sin/cos - change to tables with 2? rotations
					fpga_clt_data_rot[chn], // double  [][]  clt_tile,
					transform_size,
					residual_shift[chn][0],            // double        shiftX,
					residual_shift[chn][1],            // double        shiftY,
					true); // debug
		}

//		int byr_shift = ((ctile_top[GREEN_CHN] & 1) <<1) | (ctile_left[GREEN_CHN] & 1);
//GREEN_CHN

		// Printout
		System.out.println("// Debugging FPGA implementation");
		for (int chn = 0; chn<3; chn++) {
			System.out.println("// residual_shift["+chn+"][0]="+residual_shift[chn][0]+", residual_shift["+chn+"][1]="+residual_shift[chn][1]);
			ishx[chn] = (int) Math.round((1 << (FPGA_SHIFT_BITS)) * residual_shift[chn][0]);
			ishy[chn] = (int) Math.round((1 << (FPGA_SHIFT_BITS)) * residual_shift[chn][1]);
			if (ishx[chn] >= (1 << (FPGA_SHIFT_BITS-1))) ishx[chn] = (1 << (FPGA_SHIFT_BITS-1)) - 1;
			if (ishy[chn] >= (1 << (FPGA_SHIFT_BITS-1))) ishy[chn] = (1 << (FPGA_SHIFT_BITS-1)) - 1;
			if (ishx[chn] < -(1 << (FPGA_SHIFT_BITS-1))) ishx[chn] = -(1 << (FPGA_SHIFT_BITS-1));
			if (ishy[chn] < -(1 << (FPGA_SHIFT_BITS-1))) ishy[chn] = -(1 << (FPGA_SHIFT_BITS-1));
			residual_shift[chn][0] = ishx[chn] * (1.0/(1 << (FPGA_SHIFT_BITS)));
			residual_shift[chn][1] = ishy[chn] * (1.0/(1 << (FPGA_SHIFT_BITS)));
			System.out.println(String.format("%4x // color %d shift_x, %d bits", ishx[chn] & ((1 << (FPGA_SHIFT_BITS)) - 1),chn,FPGA_SHIFT_BITS));
			System.out.println(String.format("%4x // color %d shift_y, %d bits", ishy[chn] & ((1 << (FPGA_SHIFT_BITS)) - 1),chn,FPGA_SHIFT_BITS));
			System.out.println(String.format("%4x // color %d ctile_left", ctile_left[chn],chn));
			System.out.println(String.format("%4x // color %d ctile_top",  ctile_top[chn], chn));
		}






		System.out.println("\n// Full Bayer fpga tile data");
		int id = (1 << (FPGA_PIXEL_BITS - 9)); // 8
		for (int i = 0; i < FPGA_TILE_SIZE*FPGA_TILE_SIZE; i++) {
			double d = 0.0;
			for (int fpga_chn = 0; fpga_chn < 3; fpga_chn++){
				d +=  fpga_full_tile[fpga_chn][i];
			}
			System.out.print(String.format("%4x ",(int) Math.round(id * d)));
			if (((i+1) %FPGA_TILE_SIZE) == 0) {
				System.out.println();
			}
		}
		System.out.println();

		for (int chn = 0; chn<3; chn++) {
			double [] fpga_pix_lim = {0.0,0.0};
			for (int i = 0; i < 256; i++) if (fpga_tile[chn][i] != 0){
				if (fpga_tile[chn][i] > fpga_pix_lim[0]) fpga_pix_lim[0] = fpga_tile[chn][i];
				if (fpga_tile[chn][i] < fpga_pix_lim[1]) fpga_pix_lim[1] = fpga_tile[chn][i];
			}
			System.out.println(String.format("\n// Color # %d: Pixels input range: %f ... %f", chn, fpga_pix_lim[1], fpga_pix_lim[0]));
			System.out.println(String.format("//%x // shift_x, %d bits",ishx[chn] & ((1 << (FPGA_SHIFT_BITS)) - 1),FPGA_SHIFT_BITS));
			System.out.println(String.format("//%x // shift_y, %d bits",ishy[chn] & ((1 << (FPGA_SHIFT_BITS)) - 1),FPGA_SHIFT_BITS));
			for (int row = 0; row <16; row++){
				for (int col = 0; col <16; col++){
					System.out.print(String.format("%4x ",(int) Math.round(id * fpga_tile[chn][row*16 + col])));
				}
				System.out.println();
			}
			System.out.println();
		}
		System.out.println();
		for (int chn = 0; chn<3; chn++) {
			System.out.println("// Color="+chn+", signs table (per mode, per index - bitstring of variants, 0 - positive, 1 - negative)");
			for (int dct_mode = 0; dct_mode <4; dct_mode++) {
				for (int i = 0; i < 64; i++){
					int d = 0;
					for (int b = 0; b < 4; b++){
						if (fold_coeff[chn][dct_mode][i][b] < 0){
							d |= (1 << b);
						}
					}
					System.out.print(String.format("%x ",d));
					if ((i % 16) == 15){
						System.out.println();
					}
				}
			}
			System.out.println();
		}
		System.out.println();
		for (int chn = 0; chn<3; chn++) {
			System.out.println("// Color = "+chn+", absolute values, mode0 (CC), others are the same");
			//		for (int dct_mode = 0; dct_mode <4; dct_mode++) {
			int dct_mode = 0;
			for (int i = 0; i < 64; i++){
				for (int b = 0; b < 4; b++){
					int d = (int) Math.round(((1 << FPGA_WND_BITS) -1)* Math.abs(fold_coeff[chn][dct_mode][i][b]));
					System.out.print(String.format("%5x ",d & ((1 << (FPGA_WND_BITS)) - 1)));
				}
				if ((i % 4) == 3){
					System.out.println();
				}
			}
			System.out.println();
			//		}
		}
		System.out.println();



		for (int chn = 0; chn<3; chn++) {
			double [] fpga_dtt_lim = {0.0,0.0};
			for (int dct_mode = 0; dct_mode <4; dct_mode++) {
				for (int i = 0; i < 64; i++){
					if (fpga_clt_data_in[chn][dct_mode][i] > fpga_dtt_lim[0]) fpga_dtt_lim[0] = fpga_clt_data_in[chn][dct_mode][i];
					if (fpga_clt_data_in[chn][dct_mode][i] < fpga_dtt_lim[1]) fpga_dtt_lim[1] = fpga_clt_data_in[chn][dct_mode][i];
				}
			}
			System.out.println(String.format("// Color= %d, DTT input range: %f ... %f", chn, fpga_dtt_lim[1], fpga_dtt_lim[0]));
			double scale = (1 << (FPGA_DTT_IN - 9)); //  -1;
			scale *= 1.0 *((1 << FPGA_WND_BITS) -1) / (1 << FPGA_WND_BITS);
			scale *= 1.0 *((1 << FPGA_WND_BITS) -1) / (1 << FPGA_WND_BITS);
			scale *= 2; // Increased twice in FPGA adding twice each number in FPGA
			System.out.println("// Color="+chn+" fpga_clt_data_out[chn][dct_mode] = dtt.dttt_iv(..., scale="+scale);

//			if (chn != GREEN_CHN) scale *= 2; // adding twice each number in FPGA for R and B - done before
			for (int dct_mode = 0; dct_mode <4; dct_mode++) {
				for (int i = 0; i < 64; i++){
					int idd = (int) Math.round(scale * fpga_clt_data_in[chn][dct_mode][i]);
					System.out.print(String.format("%7x ", idd & ((1 << 25) -1)));
					if ((i % 8) == 7) System.out.println();
				}
				System.out.println();
			}
			System.out.println();
		}
		System.out.println();

		for (int chn = 0; chn<3; chn++) {
			double [] fpga_dtt_lim = {0.0,0.0};
			for (int dct_mode = 0; dct_mode <4; dct_mode++) {
				for (int i = 0; i < 64; i++){
					if (fpga_clt_data_out[chn][dct_mode][i] > fpga_dtt_lim[0]) fpga_dtt_lim[0] = fpga_clt_data_out[chn][dct_mode][i];
					if (fpga_clt_data_out[chn][dct_mode][i] < fpga_dtt_lim[1]) fpga_dtt_lim[1] = fpga_clt_data_out[chn][dct_mode][i];
				}
			}
			System.out.println(String.format("// Color = %d: DTT output range: %f ... %f", chn, fpga_dtt_lim[1], fpga_dtt_lim[0]));
			// scale = (1 << (FPGA_DTT_IN - 9)); //  -1;
//			double scale = (1 << (FPGA_DTT_IN - 8)); //  increased twice
			double scale = (1 << (FPGA_DTT_IN - 9)); //  Do not increase - lead to overflow !
			scale *= 1.0 *((1 << FPGA_WND_BITS) -1) / (1 << FPGA_WND_BITS);
			scale *= 1.0 *((1 << FPGA_WND_BITS) -1) / (1 << FPGA_WND_BITS);

			for (int dct_mode = 0; dct_mode <4; dct_mode++) {
				for (int i = 0; i < 64; i++){
					int idd = (int) Math.round(scale * fpga_clt_data_out[chn][dct_mode][i]);
					System.out.print(String.format("%7x ", idd & ((1 << 25) -1)));
					if ((i % 8) == 7) System.out.println();
				}
				System.out.println();
			}
			System.out.println();
			System.out.println("// Color = "+chn+" Testing symmetry of checkerboard patterns");
			for (int dct_mode = 0; dct_mode < 2; dct_mode++) {
				for (int i = 0; i < 64; i++){
					if ((i % 8) == 0) System.out.print("// ");
					int idd = (int) Math.round(scale * fpga_clt_data_out[chn][dct_mode][i]);
					int idd1 = (int) Math.round(scale * fpga_clt_data_out[chn][3-dct_mode][63-i]);
					System.out.print(String.format("%7x ", (idd-idd1) & ((1 << 25) -1)));
					if ((i % 8) == 7) System.out.println();
				}
				System.out.println();
			}
			System.out.println();
			System.out.println("// Color = "+chn+" Testing antisymmetry of checkerboard patterns");
			for (int dct_mode = 0; dct_mode < 2; dct_mode++) {
				for (int i = 0; i < 64; i++){
					if ((i % 8) == 0) System.out.print("// ");
					int idd = (int) Math.round(scale * fpga_clt_data_out[chn][dct_mode][i]);
					int idd1 = (int) Math.round(scale * fpga_clt_data_out[chn][3-dct_mode][63-i]);
					System.out.print(String.format("%7x ", (idd+idd1) & ((1 << 25) -1)));
					if ((i % 8) == 7) System.out.println();
				}
				System.out.println();
			}
			System.out.println();
		}
		System.out.println();

		for (int chn = 0; chn<3; chn++) {
//			double scale = (1 << (FPGA_DTT_IN - 9)); //  -1;
//			double scale = (1 << (FPGA_DTT_IN - 8)); //
			double scale = (1 << (FPGA_DTT_IN - 9)); //  Do not increase - lead to overflow !

			// compensate for DTT scale
			scale *= 1.0 *((1 << FPGA_WND_BITS) -1) / (1 << FPGA_WND_BITS);
			scale *= 1.0 *((1 << FPGA_WND_BITS) -1) / (1 << FPGA_WND_BITS);
			// compensate for rotator scale:
			scale *= 1.0 *((1 << FPGA_WND_BITS) -1) / (1 << FPGA_WND_BITS);
			scale *= 1.0 *((1 << FPGA_WND_BITS) -1) / (1 << FPGA_WND_BITS);
			double [] fpga_dtt_lim = {0.0,0.0};
			for (int dct_mode = 0; dct_mode <4; dct_mode++) {
				for (int j = 0; j < 64; j++){
					if (fpga_clt_data_rot[chn][dct_mode][j] > fpga_dtt_lim[0]) fpga_dtt_lim[0] = fpga_clt_data_rot[chn][dct_mode][j];
					if (fpga_clt_data_rot[chn][dct_mode][j] < fpga_dtt_lim[1]) fpga_dtt_lim[1] = fpga_clt_data_rot[chn][dct_mode][j];
				}
			}

			System.out.println(String.format("// Color = %d: DTT rotated, shift_x=%f. shift_y = %f", chn, residual_shift[chn][0],residual_shift[chn][1]));
			System.out.println(String.format("// DTT rotated  range: %f ... %f", fpga_dtt_lim[1], fpga_dtt_lim[0]));
//			scale = (1 << (FPGA_DTT_IN - 9)); //  -1;
			for (int dct_mode = 0; dct_mode <4; dct_mode++) {
				for (int j = 0; j < 64; j++){
					int idd = (int) Math.round(scale * fpga_clt_data_rot[chn][dct_mode][j]);
					System.out.print(String.format("%7x ", idd & ((1 << 25) -1)));
					if ((j % 8) == 7) System.out.println();
				}
				System.out.println();
			}
		}
	}


	public double [][][][][][] clt_aberrations_quad_corr_new(
			final ImageDttParameters  imgdtt_params,   // Now just extra correlation parameters, later will include, most others
			final int                 macro_scale,     // to correlate tile data instead of the pixel data: 1 - pixels, 8 - tiles
			final int [][]            tile_op,         // [tilesY][tilesX] - what to do - 0 - nothing for this tile
			final double [][]         disparity_array, // [tilesY][tilesX] - individual per-tile expected disparity
			final double [][][]       image_data, // first index - number of image in a quad
		    final boolean [][]        saturation_imp, // (near) saturated pixels or null
			 // correlation results - final and partial
			final double [][][][]     clt_corr_combo,  // [type][tilesY][tilesX][(2*transform_size-1)*(2*transform_size-1)] // if null - will not calculate
			                                           // [type][tilesY][tilesX] should be set by caller
													   // types: 0 - selected correlation (product+offset), 1 - sum

			final double [][][][][]   clt_corr_partial,// [tilesY][tilesX][quad]color][(2*transform_size-1)*(2*transform_size-1)] // if null - will not calculate
                                                       // [tilesY][tilesX] should be set by caller
			// When clt_mismatch is non-zero, no far objects extraction will be attempted
			final double [][]         clt_mismatch,    // [12][tilesY * tilesX] // ***** transpose unapplied ***** ?. null - do not calculate
			                                           // values in the "main" directions have disparity (*_CM) subtracted, in the perpendicular - as is

			final double [][]         disparity_map,   // [8][tilesY][tilesX], only [6][] is needed on input or null - do not calculate
			                                           // last 2 - contrast, avg/ "geometric average)
			final double [][][][]     texture_tiles,   // [tilesY][tilesX]["RGBA".length()][];  null - will skip images combining

			final int                 width,
			final double              corr_fat_zero,    // add to denominator to modify phase correlation (same units as data1, data2). <0 - pure sum
			final boolean             corr_sym,
			final double              corr_offset,
			final double              corr_red,
			final double              corr_blue,
			final double              corr_sigma,
			final boolean             corr_normalize,  // normalize correlation results by rms
	  		final double              min_corr,        // 0.02; // minimal correlation value to consider valid
			final double              max_corr_sigma,  // 1.2;  // weights of points around global max to find fractional
			final double              max_corr_radius, // 3.9;
			final boolean 			  max_corr_double, //"Double pass when masking center of mass to reduce preference for integer values
			final int                 corr_mode, // Correlation mode: 0 - integer max, 1 - center of mass, 2 - polynomial
			final double              min_shot,        // 10.0;  // Do not adjust for shot noise if lower than
			final double              scale_shot,      // 3.0;   // scale when dividing by sqrt ( <0 - disable correction)
			final double              diff_sigma,      // 5.0;//RMS difference from average to reduce weights (~ 1.0 - 1/255 full scale image)
			final double              diff_threshold,  // 5.0;   // RMS difference from average to discard channel (~ 1.0 - 1/255 full scale image)
			final boolean             diff_gauss,      // true;  // when averaging images, use Gaussian around average as weight (false - sharp all/nothing)
			final double              min_agree,       // 3.0;   // minimal number of channels to agree on a point (real number to work with fuzzy averages)
			final boolean             dust_remove,     // Do not reduce average weight when only one image differs much from the average
			final boolean             keep_weights,    // Add port weights to RGBA stack (debug feature)
			final GeometryCorrection  geometryCorrection,
			final GeometryCorrection  geometryCorrection_main, // if not null correct this camera (aux) to the coordinates of the main
			final double [][][][][][] clt_kernels, // [channel_in_quad][color][tileY][tileX][band][pixel] , size should match image (have 1 tile around)
			final int                 kernel_step,
			final int                 transform_size,
			final int                 window_type,
			final double [][]         shiftXY, // [port]{shiftX,shiftY}
			final double              disparity_corr, // disparity at infinity
			final double [][][]       fine_corr, // quadratic coefficients for fine correction (or null)
			final double              corr_magic_scale, // still not understood coefficient that reduces reported disparity value.  Seems to be around 0.85
			final double              shiftX, // shift image horizontally (positive - right) - just for testing
			final double              shiftY, // shift image vertically (positive - down)
			final int                 debug_tileX,
			final int                 debug_tileY,
			final boolean             no_fract_shift,
			final boolean             no_deconvolution,
			final int                 threadsMax,  // maximal number of threads to launch
			final int                 globalDebugLevel)
	{
		final boolean macro_mode = macro_scale != 1;      // correlate tile data instead of the pixel data
		final int quad = 4;   // number of subcameras

		final int numcol = 3; // number of colors // keep the same, just do not use [0] and [1], [2] - green

//		final int numColors = image_data[0].length;
		final int height=image_data[0][0].length/width;
		final int tilesX=width/transform_size;
		final int tilesY=height/transform_size;
		final int nTilesInChn=tilesX*tilesY;
		final double [][][][][][] clt_data = new double[quad][numcol][tilesY][tilesX][][];
		final Thread[] threads = newThreadArray(threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		final double [] col_weights= new double [numcol]; // colors are RBG

		// keep for now for mono, find out  what do they mean for macro mode
		if (macro_mode) { // all the same as they now mean different
			//compensating Bayer correction
			col_weights[0] = 0.25; //  1.0/3;
			col_weights[1] = 0.25; //  1.0/3;
			col_weights[2] = 0.5; // 1.0/3;
		} else {
			if (isMonochrome()) {
				col_weights[2] = 1.0;// green color/mono
				col_weights[0] = 0;
				col_weights[1] = 0;
			} else {
				col_weights[2] = 1.0/(1.0 + corr_red + corr_blue);    // green color
				col_weights[0] = corr_red *  col_weights[2];
				col_weights[1] = corr_blue * col_weights[2];
			}
		}

		final int corr_size = transform_size * 2 -1;
		final int [][] transpose_indices = new int [corr_size*(corr_size-1)/2][2];
		int indx = 0;
		if ((globalDebugLevel > -10) && (disparity_corr != 0.0)){
			System.out.println(String.format("Using manual infinity disparity correction of %8.5f pixels",disparity_corr));
		}
		for (int i =0; i < corr_size-1; i++){
			for (int j = i+1; j < corr_size; j++){
				transpose_indices[indx  ][0] = i * corr_size + j;
				transpose_indices[indx++][1] = j * corr_size + i;
			}
		}


		final int first_color = isMonochrome()? MONO_CHN : 0; // color that is non-zero

		// reducing weight of on-axis correlation values to enhance detection of vertical/horizontal lines
		// multiply correlation results inside the horizontal center strip  2*enhortho_width - 1 wide by enhortho_scale

		final double [] enh_ortho_scale = new double [corr_size];
		for (int i = 0; i < corr_size; i++){
			if ((i < (transform_size - imgdtt_params.getEnhOrthoWidth(isAux()))) || (i > (transform_size - 2 + imgdtt_params.getEnhOrthoWidth(isAux())))) {
				enh_ortho_scale[i] = 1.0;
			} else {
				enh_ortho_scale[i] = imgdtt_params.getEnhOrthoScale(isAux());
			}
			if (i == (transform_size-1)) enh_ortho_scale[i] = 0.0 ; // hardwired 0 in the center
			enh_ortho_scale[i] *= Math.sin(Math.PI*(i+1.0)/(2*transform_size));
		}
		if (globalDebugLevel > 1){
			System.out.println("getEnhOrthoWidth(isAux())="+ imgdtt_params.getEnhOrthoWidth(isAux())+" getEnhOrthoScale(isAux())="+ imgdtt_params.getEnhOrthoScale(isAux()));
			for (int i = 0; i < corr_size; i++){
				System.out.println(" enh_ortho_scale["+i+"]="+ enh_ortho_scale[i]);

			}
		}

		// Create window  to select center correlation strip using
		// ortho_height - full width of non-zero elements
		// ortho_eff_height - effective height (ration of the weighted column sum to the center value)
		int wcenter = transform_size - 1;
		final double [] ortho_weights = new double [corr_size]; // [15]
		for (int i = 0; i < corr_size; i++){
			if ((i >= wcenter - imgdtt_params.ortho_height/2) && (i <= wcenter + imgdtt_params.ortho_height/2)) {
				double dx = 1.0*(i-wcenter)/(imgdtt_params.ortho_height/2 + 1);
				ortho_weights[i] = 0.5*(1.0+Math.cos(Math.PI*dx))/imgdtt_params.ortho_eff_height;
			}
		}
		if (globalDebugLevel > 0){
			System.out.println("ortho_height="+ imgdtt_params.ortho_height+" ortho_eff_height="+ imgdtt_params.ortho_eff_height);
			for (int i = 0; i < corr_size; i++){
				System.out.println(" ortho_weights["+i+"]="+ ortho_weights[i]);
			}
		}




		if (globalDebugLevel > 0) {
			System.out.println("clt_aberrations_quad_corr(): width="+width+" height="+height+" transform_size="+transform_size+
					" debug_tileX="+debug_tileX+" debug_tileY="+debug_tileY+" globalDebugLevel="+globalDebugLevel);
		}
		final int [][] zi =
			{{ 0,  1,  2,  3},
			 {-1,  0, -3,  2},
			 {-2, -3,  0,  1},
			 { 3, -2, -1,  0}};
		final int [][] corr_pairs ={ // {first, second, rot} rot: 0 - as is, 1 - swap y,x
				{0,1,0},
				{2,3,0},
				{0,2,1},
				{1,3,1}};

		final double[][] port_offsets = {
				{-0.5, -0.5},
				{ 0.5, -0.5},
				{-0.5,  0.5},
				{ 0.5,  0.5}};
		final int transform_len = transform_size * transform_size;



		final double [] filter_direct= new double[transform_len];
		if (corr_sigma == 0) {
			filter_direct[0] = 1.0;
			for (int i= 1; i<filter_direct.length;i++) filter_direct[i] =0;
		} else {
			for (int i = 0; i < transform_size; i++){
				for (int j = 0; j < transform_size; j++){
					filter_direct[i*transform_size+j] = Math.exp(-(i*i+j*j)/(2*corr_sigma)); // FIXME: should be sigma*sigma !
				}
			}
		}
		// normalize
		double sum = 0;
		for (int i = 0; i < transform_size; i++){
			for (int j = 0; j < transform_size; j++){
				double d = 	filter_direct[i*transform_size+j];
				d*=Math.cos(Math.PI*i/(2*transform_size))*Math.cos(Math.PI*j/(2*transform_size));
				if (i > 0) d*= 2.0;
				if (j > 0) d*= 2.0;
				sum +=d;
			}
		}
		for (int i = 0; i<filter_direct.length; i++){
			filter_direct[i] /= sum;
		}

		DttRad2 dtt = new DttRad2(transform_size);
		final double [] filter= dtt.dttt_iiie(filter_direct);
		for (int i=0; i < filter.length;i++) filter[i] *= 2*transform_size;

		// prepare disparity maps and weights
		final int max_search_radius = (int) Math.abs(max_corr_radius); // use negative max_corr_radius for squares instead of circles?
		final int max_search_radius_poly = 1;
		if (globalDebugLevel > 0){
			System.out.println("max_corr_radius=       "+max_corr_radius);
			System.out.println("max_search_radius=     "+max_search_radius);
			System.out.println("max_search_radius_poly="+max_search_radius_poly);
			System.out.println("corr_fat_zero=         "+corr_fat_zero);
			System.out.println("disparity_array[0][0]= "+disparity_array[0][0]);


		}
		// add optional initialization of debug layers here
		if (disparity_map != null){
			for (int i = 0; i<disparity_map.length;i++){
				if (i < OVEREXPOSED) {
					disparity_map[i] = new double [tilesY*tilesX];
				} else if (i == OVEREXPOSED) {
					if (saturation_imp!= null) {
					disparity_map[i] = new double [tilesY*tilesX];
					}
				} else if (i >= IMG_TONE_RGB) {
					if (texture_tiles != null) { // for now - enable 12 tone layers only together with texture tiles
						disparity_map[i] = new double [tilesY*tilesX];
					}
				}
			}
		}
		if (clt_mismatch != null){
			for (int i = 0; i<clt_mismatch.length;i++){
				clt_mismatch[i] = new double [tilesY*tilesX]; // will use only "center of mass" centers
			}
		}

		dtt.set_window(window_type);
		final double [] lt_window = dtt.getWin2d();	// [256]
		final double [] lt_window2 = new double [lt_window.length]; // squared
		for (int i = 0; i < lt_window.length; i++) lt_window2[i] = lt_window[i] * lt_window[i];

		if (globalDebugLevel > 1) {
			ShowDoubleFloatArrays sdfa_instance = new ShowDoubleFloatArrays(); // just for debugging?
			sdfa_instance.showArrays(lt_window,  2*transform_size, 2*transform_size, "lt_window");
		}
		if (globalDebugLevel > 0) {
			System.out.println("macro_mode="+macro_mode);
		}

		Matrix [] corr_rots_aux = null;
		if (geometryCorrection_main != null) {
			corr_rots_aux = geometryCorrection.getCorrVector().getRotMatrices(geometryCorrection.getRotMatrix(true));
		}

		final boolean use_main = corr_rots_aux != null;
		final Matrix [] corr_rots = use_main ? corr_rots_aux : geometryCorrection.getCorrVector().getRotMatrices(); // get array of per-sensor rotation matrices
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				@Override
				public void run() {
					DttRad2 dtt = new DttRad2(transform_size);
					dtt.set_window(window_type);
					int tileY,tileX,tIndex; // , chn;
					//						showDoubleFloatArrays sdfa_instance = new showDoubleFloatArrays(); // just for debugging?
					double centerX; // center of aberration-corrected (common model) tile, X
					double centerY; //
					double [][] fract_shiftsXY = new double[quad][];
					double [][]     tcorr_combo =    null; // [15*15] pixel space
					double [][][]   tcorr_partial =  null; // [quad][numcol+1][15*15]
					double [][][][] tcorr_tpartial = null; // [quad][numcol+1][4][8*8]
					double [] ports_rgb = null;
					Correlation2d corr2d = new Correlation2d(
							imgdtt_params,              // ImageDttParameters  imgdtt_params,
							transform_size,             // int transform_size,
							2.0,                        //  double wndx_scale, // (wndy scale is always 1.0)
							isMonochrome(), // boolean monochrome,
							(globalDebugLevel > -1));   //   boolean debug)
					corr2d.createOrtoNotch(
							imgdtt_params.getEnhOrthoWidth(isAux()), // double getEnhOrthoWidth(isAux()),
							imgdtt_params.getEnhOrthoScale(isAux()), //double getEnhOrthoScale(isAux()),
							(imgdtt_params.lma_debug_level > 1)); // boolean debug);

					for (int nTile = ai.getAndIncrement(); nTile < nTilesInChn; nTile = ai.getAndIncrement()) {

						tileY = nTile /tilesX;
						tileX = nTile % tilesX;
						tIndex = tileY * tilesX + tileX;
						if (tile_op[tileY][tileX] == 0) continue; // nothing to do for this tile
						int                 img_mask =  getImgMask(tile_op[tileY][tileX]);         // which images to use
						int                 corr_mask = getPairMask(tile_op[tileY][tileX]);       // which pairs to combine in the combo:  1 - top, 2 bottom, 4 - left, 8 - right
						// mask out pairs that use missing channels
						for (int i = 0; i< corr_pairs.length; i++){
							if ((((1 << corr_pairs[i][0]) & img_mask) == 0) || (((1 << corr_pairs[i][1]) & img_mask) == 0)) {
								corr_mask &= ~ (1 << i);
							}
						}
						boolean debugTile =(tileX == debug_tileX) && (tileY == debug_tileY);

						final int [] overexp_all = (saturation_imp != null) ? ( new int [2]): null;

						// Moved from inside chn loop
						centerX = tileX * transform_size + transform_size/2 - shiftX;
						centerY = tileY * transform_size + transform_size/2 - shiftY;
						// TODO: move port coordinates out of color channel loop
						double [][] centersXY;


						if (macro_mode){
							if ((globalDebugLevel > -1) && (tileX == debug_tileX) && (tileY == debug_tileY)) { // before correction
								System.out.println("\nUsing MACRO mode, centerX="+centerX+", centerY="+centerY);
							}
							centersXY = geometryCorrection.getPortsCoordinatesIdeal(
									macro_scale,
									centerX,
									centerY,
									macro_scale* disparity_array[tileY][tileX] + disparity_corr);

						} else {
							if (use_main) { // this is AUX camera that uses main coordinates
								centersXY =  geometryCorrection.getPortsCoordinatesAndDerivatives(
										geometryCorrection_main, //			GeometryCorrection gc_main,
										true,            // boolean use_rig_offsets,
										corr_rots,       // Matrix []   rots,
										null,            //  Matrix [][] deriv_rots,
										null,            // double [][] pXYderiv, // if not null, should be double[8][]
										centerX,
										centerY,
										disparity_array[tileY][tileX] + disparity_corr); // _aux); //  + disparity_corr);


							} else {
								centersXY = geometryCorrection.getPortsCoordinatesAndDerivatives(
										geometryCorrection, //			GeometryCorrection gc_main,
										false,          // boolean use_rig_offsets,
										corr_rots, // Matrix []   rots,
										null,      //  Matrix [][] deriv_rots,
										null,      // double [][] pXYderiv, // if not null, should be double[8][]
										centerX,
										centerY,
										disparity_array[tileY][tileX] + disparity_corr);
							}

							if ((globalDebugLevel > 0) && (tileX == debug_tileX) && (tileY == debug_tileY)) {
								for (int i = 0; i < quad; i++) {
									System.out.println("clt_aberrations_quad_corr():  tileX="+tileX+", tileY="+tileY+
											" centerX="+centerX+" centerY="+centerY+" disparity="+disparity_array[tileY][tileX]+
											" centersXY["+i+"][0]="+centersXY[i][0]+" centersXY["+i+"][1]="+centersXY[i][1]);
								}
							}

							if ((globalDebugLevel > -1) && (tileX == debug_tileX) && (tileY == debug_tileY)) { // before correction
								System.out.print(disparity_array[tileY][tileX]+"\t"+
										centersXY[0][0]+"\t"+centersXY[0][1]+"\t"+
										centersXY[1][0]+"\t"+centersXY[1][1]+"\t"+
										centersXY[2][0]+"\t"+centersXY[2][1]+"\t"+
										centersXY[3][0]+"\t"+centersXY[3][1]+"\t");
							}

							for (int ip = 0; ip < centersXY.length; ip++){
								centersXY[ip][0] -= shiftXY[ip][0];
								centersXY[ip][1] -= shiftXY[ip][1];
							}
							// TODO: use correction after disparity applied (to work for large disparity values)
							if (fine_corr != null){

								for (int ip = 0; ip < centersXY.length; ip++){
									double [] tXY = geometryCorrection.getRelativeCoords(centersXY[ip]);
									for (int d = 0; d <2; d++) {
										centersXY[ip][d] -= (
												fine_corr[ip][d][0]*tXY[0]*tXY[0]+
												fine_corr[ip][d][1]*tXY[1]*tXY[1]+
												fine_corr[ip][d][2]*tXY[0]*tXY[1]+
												fine_corr[ip][d][3]*tXY[0]+
												fine_corr[ip][d][4]*tXY[1]+
												fine_corr[ip][d][5]);
									}
								}
							}
						} // if (macro_mode) ... else
						if (FPGA_COMPARE_DATA && (globalDebugLevel > 0) && (tileX == debug_tileX) && (tileY == debug_tileY)) {
							final int fpga_cam = 0;
							double [][] manual_offsets={
//									{ 1.3, -2.7},
//									{-1.3,  2.7},
//									{ 0.0,  0.0}};

//							{ 2.3, -2.7},
//							{-0.3,  2.7},
//							{ 0.0,  0.0}};

							{ 2.3, -2.7},
							{-0.3,  2.7},
							{ 1.0,  0.0}};

							double [][] colorCentersXY = {
									{centersXY[fpga_cam][0] + manual_offsets[0][0], centersXY[fpga_cam][1] + manual_offsets[0][1]}, // add manual offsets here
									{centersXY[fpga_cam][0] + manual_offsets[1][0], centersXY[fpga_cam][1] + manual_offsets[1][1]},
									{centersXY[fpga_cam][0] + manual_offsets[2][0], centersXY[fpga_cam][1] + manual_offsets[2][1]}
							};
							generateFPGACompareData(
									image_data[fpga_cam], // final double [][]                  image_data, // for selected subcamera
									colorCentersXY,       // final double [][]         colorCentersXY, // pixel centers per color (2 - green)
									transform_size,       // final int                 transform_size,
									width,                 // final int                 width
									dtt
									);
						}
// See if macro_mode uses color channels for non-color?
						for (int ncol = 0; ncol <numcol; ncol++) {
							if (!isMonochrome() || (ncol == MONO_CHN) || macro_mode) { // in monochrome mode skip all non-mono (green) channels
								boolean debug_for_fpga = FPGA_COMPARE_DATA && (globalDebugLevel > 0) && (tileX == debug_tileX) && (tileY == debug_tileY) && (ncol == 2);
								if ((globalDebugLevel > -1) && (tileX == debug_tileX) && (tileY == debug_tileY) && (ncol == 2)) {
									System.out.println("\nUsing "+(macro_mode?"MACRO":"PIXEL")+" mode, centerX="+centerX+", centerY="+centerY);
									System.out.println(disparity_array[tileY][tileX]+"\t"+
											centersXY[0][0]+"\t"+centersXY[0][1]+"\t"+
											centersXY[1][0]+"\t"+centersXY[1][1]+"\t"+
											centersXY[2][0]+"\t"+centersXY[2][1]+"\t"+
											centersXY[3][0]+"\t"+centersXY[3][1]+"\t");
								}

								for (int i = 0; i < quad; i++) {
									if (debug_for_fpga && (i==0)){
										double [][] fpga_clt_data = new double [4][];
										double [] fpga_fract_shiftsXY;
										double [] fpga_centersXY = {centersXY[i][0],centersXY[i][1]};
										int fpga_chn = ncol; // ==2, green
										// round to nearest 1/128 pix (supported by FPGA code)
										System.out.println(String.format("Center X= %f, center Y = %f", fpga_centersXY[0],fpga_centersXY[1]));
										for (int j=0; j<2;j++){
											fpga_centersXY[j] = Math.round(128*fpga_centersXY[j])/128.0;
										}


										for (int j=0; j<2;j++){
											fpga_centersXY[j] = Math.round(fpga_centersXY[j]);
										}


										//									fpga_centersXY[0]+=0.5; // half pixel shift horizontal zero pixel shift vertical
										//									fpga_centersXY[1]+=0.5; // half pixel shift vertical, zero pixel shift horizontal

										//									fpga_centersXY[0]+=1.0; //

										fpga_chn =           2;


										System.out.println(String.format("Manually changing offset: center X= %f, center Y = %f", fpga_centersXY[0],fpga_centersXY[1]));
										System.out.println(String.format("Manually changing color to %d (was %d)", fpga_chn, ncol));



										fpga_fract_shiftsXY = extract_correct_tile( // return a pair of residual offsets
												image_data[i],
												width,       // image width
												null,
												fpga_clt_data, //double  [][]        clt_tile,    // should be double [4][];
												kernel_step,
												transform_size,
												dtt,
												fpga_chn, // chn,
												fpga_centersXY[0], // centersXY[i][0], // centerX, // center of aberration-corrected (common model) tile, X
												fpga_centersXY[1], // centersXY[i][1], // centerY, //
												-10, // globalDebugLevel,
												true, // no_deconvolution,
												false, // ); // transpose);
												null,
												null);
										ShowDoubleFloatArrays sdfa_instance = new ShowDoubleFloatArrays(); // just for debugging?
										String [] titles = {"CC","SC","CS","SS"};
										double [][] dbg_tile = new double [4][];
										for (int im = 0; im < 4; im++) dbg_tile[im]=fpga_clt_data[im];
										sdfa_instance.showArrays(dbg_tile,  transform_size, transform_size, true, "fpre-shifted_x"+tileX+"_y"+tileY+"-z", titles);
										fract_shift(    // fractional shift in transform domain. Currently uses sin/cos - change to tables with 2? rotations
												fpga_clt_data, // double  [][]  clt_tile,
												transform_size,
												fpga_fract_shiftsXY[0],            // double        shiftX,
												fpga_fract_shiftsXY[1],            // double        shiftY,
												true); // debug
										for (int im = 0; im < 4; im++) dbg_tile[im]=fpga_clt_data[im];
										sdfa_instance.showArrays(dbg_tile,  transform_size, transform_size, true, "f-shifted_x"+tileX+"_y"+tileY+"-z", titles);
										System.out.println("Debugging for FPGA data, globalDebugLevel = "+globalDebugLevel+", tileX="+tileX+", tileY="+tileY+", sesnlor="+i+", color="+ncol);
										System.out.println("Debugging for FPGA data, fpga_fract_shiftsXY[0] = "+fpga_fract_shiftsXY[0]+", fpga_fract_shiftsXY[1]="+fpga_fract_shiftsXY[1]);
										System.out.println();

										double scale = (1 << (FPGA_DTT_IN - 9)); //  -1;
										// compensate for DTT scale
										scale *= 1.0 *((1 << FPGA_WND_BITS) -1) / (1 << FPGA_WND_BITS);
										scale *= 1.0 *((1 << FPGA_WND_BITS) -1) / (1 << FPGA_WND_BITS);
										// compensate for rotator scale:
										scale *= 1.0 *((1 << FPGA_WND_BITS) -1) / (1 << FPGA_WND_BITS);
										scale *= 1.0 *((1 << FPGA_WND_BITS) -1) / (1 << FPGA_WND_BITS);
										double [] fpga_dtt_lim = {0.0,0.0};
										for (int dct_mode = 0; dct_mode <4; dct_mode++) {
											for (int j = 0; j < 64; j++){
												if (fpga_clt_data[dct_mode][j] > fpga_dtt_lim[0]) fpga_dtt_lim[0] = fpga_clt_data[dct_mode][j];
												if (fpga_clt_data[dct_mode][j] < fpga_dtt_lim[1]) fpga_dtt_lim[1] = fpga_clt_data[dct_mode][j];
											}
										}

										System.out.println(String.format("// DTT rotated, shift_x=%f. shift_y = %f", fpga_fract_shiftsXY[0],fpga_fract_shiftsXY[1]));
										System.out.println(String.format("// DTT rotated  range: %f ... %f", fpga_dtt_lim[1], fpga_dtt_lim[0]));
										//									scale = (1 << (FPGA_DTT_IN - 9)); //  -1;
										for (int dct_mode = 0; dct_mode <4; dct_mode++) {
											for (int j = 0; j < 64; j++){
												int id = (int) Math.round(scale * fpga_clt_data[dct_mode][j]);
												System.out.print(String.format("%7x ", id & ((1 << 25) -1)));
												if ((j % 8) == 7) System.out.println();
											}
											System.out.println();
										}
									} // end of debug_for_fpga

									clt_data[i][ncol][tileY][tileX] = new double [4][];
									fract_shiftsXY[i] = extract_correct_tile( // return a pair of residual offsets
											image_data[i],
											width,       // image width
											(clt_kernels == null) ? null : clt_kernels[i], // [color][tileY][tileX][band][pixel]
													clt_data[i][ncol][tileY][tileX], //double  [][]        clt_tile,    // should be double [4][];
													kernel_step,
													transform_size,
													dtt,
													ncol,
													centersXY[i][0], // centerX, // center of aberration-corrected (common model) tile, X
													centersXY[i][1], // centerY, //
													(!FPGA_COMPARE_DATA && (globalDebugLevel > -1) && (tileX == debug_tileX) && (tileY == debug_tileY) && (ncol == 2) && (i==0)) ? (globalDebugLevel + 0) : 0, // external tile compare
															no_deconvolution,
															false, // ); // transpose);
															((saturation_imp != null) ? saturation_imp[i] : null), //final boolean [][]        saturation_imp, // (near) saturated pixels or null
															((saturation_imp != null) ? overexp_all: null)); // final double [] overexposed)
								} // for (int i = 0; i < quad; i++)
								if ((globalDebugLevel > -1) && (tileX == debug_tileX) && (tileY == debug_tileY) && (ncol == 2)) {
									System.out.println();
								}
								if ((globalDebugLevel > 0) && (debug_tileX == tileX) && (debug_tileY == tileY)  && (ncol == 2) && !FPGA_COMPARE_DATA) {
									ShowDoubleFloatArrays sdfa_instance = new ShowDoubleFloatArrays(); // just for debugging?
									String [] titles = {"CC0","SC0","CS0","SS0","CC1","SC1","CS1","SS1","CC2","SC2","CS2","SS2","CC3","SC3","CS3","SS3"};
									double [][] dbg_tile = new double [16][];
									for (int i = 0; i < 16; i++) dbg_tile[i]=clt_data[i>>2][ncol][tileY][tileX][i & 3];
									sdfa_instance.showArrays(dbg_tile,  transform_size, transform_size, true, "pre-shifted_x"+tileX+"_y"+tileY, titles);
								}

								if ((globalDebugLevel > 0) && (tileX >= debug_tileX - 2) && (tileX <= debug_tileX + 2) &&
										(tileY >= debug_tileY - 2) && (tileY <= debug_tileY+2)) {
									for (int i = 0; i < quad; i++) {
										System.out.println("clt_aberrations_quad(): color="+ncol+", tileX="+tileX+", tileY="+tileY+
												" fract_shiftsXY["+i+"][0]="+fract_shiftsXY[i][0]+" fract_shiftsXY["+i+"][1]="+fract_shiftsXY[i][1]);
									}
								}

								if (!no_fract_shift) {
									// apply residual shift
									for (int i = 0; i < quad; i++) {
										fract_shift(    // fractional shift in transform domain. Currently uses sin/cos - change to tables with 2? rotations
												clt_data[i][ncol][tileY][tileX], // double  [][]  clt_tile,
												transform_size,
												fract_shiftsXY[i][0],            // double        shiftX,
												fract_shiftsXY[i][1],            // double        shiftY,
												//									(globalDebugLevel > 0) && (tileX == debug_tileX) && (tileY == debug_tileY)); // external tile compare
												((globalDebugLevel > 1) &&
														((ncol==0) || isMonochrome()) &&
														(tileX >= debug_tileX - 2) && (tileX <= debug_tileX + 2) &&
														(tileY >= debug_tileY - 2) && (tileY <= debug_tileY+2)));
									}
									if ((globalDebugLevel > 0) && (debug_tileX == tileX) && (debug_tileY == tileY)) {
										ShowDoubleFloatArrays sdfa_instance = new ShowDoubleFloatArrays(); // just for debugging?
										String [] titles = {"CC0","SC0","CS0","SS0","CC1","SC1","CS1","SS1","CC2","SC2","CS2","SS2","CC3","SC3","CS3","SS3"};
										double [][] dbg_tile = new double [16][];
										for (int i = 0; i < 16; i++) dbg_tile[i]=clt_data[i>>2][ncol][tileY][tileX][i & 3];
										sdfa_instance.showArrays(dbg_tile,  transform_size, transform_size, true, "shifted_x"+tileX+"_y"+tileY+"-z", titles);
									}



								}
							} else { // if (!isMonochrome() || (chn == MONO_CHN) || macro_mode) { // in monochrome mode skip all non-mono (green) channels
								for (int i = 0; i < quad; i++) {
									clt_data[i][ncol] = null; // erase unused clt_data
								}
							}
						}// end of for (int chn = 0; chn <numcol; chn++)

						int tile_lma_debug_level =  ((tileX == debug_tileX) && (tileY == debug_tileY))? imgdtt_params.lma_debug_level : -1;

						// all color channels are done here
						double extra_disparity = 0.0; // used for textures:  if allowed, shift images extra before trying to combine

						// fill clt_corr_combo if it exists
						if (disparity_map != null){ // not null - calculate correlations
							for (int i = 0; i < disparity_map.length; i++) {
								if (disparity_map[i] != null) disparity_map[i][nTile] = (
										(i == DISPARITY_STRENGTH_INDEX) ||
										(i == DISPARITY_INDEX_HOR_STRENGTH) ||
										(i == DISPARITY_INDEX_VERT_STRENGTH)) ? 0.0 : Double.NaN; // once and for all
							}
							//clt_mismatch should only be used with disparity_map != null;
							if (clt_mismatch != null) {
								for (int np = 0; np < clt_mismatch.length/3; np++) {
									clt_mismatch[3 * np + 0 ][tIndex] = Double.NaN;
									clt_mismatch[3 * np + 1 ][tIndex] = Double.NaN;
									clt_mismatch[3 * np + 2 ][tIndex] = 0;
								}
							}

							// calculate overexposed fraction
							if (saturation_imp != null){
								disparity_map[OVEREXPOSED][nTile] = (1.0 * overexp_all[0]) / overexp_all[1];
							}

							// calculate all selected pairs correlations
							int all_pairs = imgdtt_params.dbg_pair_mask; //TODO: use tile tasks

							double [][]  corrs = corr2d.correlateCompositeFD( // now works with nulls for some clt_data colors
						    		clt_data,        // double [][][][][][] clt_data,
						    		tileX,           // int                 tileX,
						    		tileY,           // int                 tileY,
						    		all_pairs,       // int                 pairs_mask,
						    		filter,          // double []           lpf,
						    		scale_strengths, // double              scale_value, // scale correlation value
						    		col_weights,     // double []           col_weights,
						    		corr_fat_zero);  // double              fat_zero)

						    // calculate interpolated "strips" to match different scales and orientations (ortho/diagonal) on the
						    // fine (0.5 pix) grid. ortho for scale == 1 provide even/even samples (1/4 of all), diagonal ones -
						    // checkerboard pattern

						    double [][] strips = corr2d.scaleRotateInterpoateCorrelations(
						    		corrs,                          // double [][] correlations,
						    		all_pairs,                      // int         pairs_mask,
						    		imgdtt_params.corr_strip_hight, //);    // int         hwidth);
						    		(tile_lma_debug_level > 0) ? all_pairs:0); // debugMax);

						    // Combine strips for selected pairs. Now using only for all available pairs.
						    // Other combinations are used only if requested (clt_corr_partial != null)

						    double [] strip_combo = corr2d.combineInterpolatedCorrelations(
						    		strips,                        // double [][] strips,
						    		all_pairs,                     // int         pairs_mask,
						    		imgdtt_params.corr_offset,     // double      offset);
						    		imgdtt_params.twice_diagonal); //    		boolean     twice_diagonal)

						    // Debug feature - only calculated if requested
						    if ((clt_corr_partial != null) && imgdtt_params.corr_mode_debug) {
							    double [] strip_ortho = corr2d.combineInterpolatedCorrelations(
							    		strips,                         // double [][] strips,
							    		0x0f,                           // int         pairs_mask,
							    		imgdtt_params.corr_offset,      // double      offset);
							    		imgdtt_params.twice_diagonal);  //    		boolean     twice_diagonal)

							    double [] strip_diag = corr2d.combineInterpolatedCorrelations(
							    		strips,                         // double [][] strips,
							    		0x30,                           // int         pairs_mask,
							    		imgdtt_params.corr_offset,      // double      offset);
							    		imgdtt_params.twice_diagonal);  //    		boolean     twice_diagonal)

							    double [] strip_all = corr2d.combineInterpolatedCorrelations(
							    		strips,                         // double [][] strips,
							    		0x3f,                           // int         pairs_mask,
							    		imgdtt_params.corr_offset,      // double      offset);
							    		imgdtt_params.twice_diagonal);  //    		boolean     twice_diagonal)

							    double [] strip_hor = corr2d.combineInterpolatedCorrelations(
							    		strips,                         // double [][] strips,
							    		0x03,                           // int         pairs_mask,
							    		imgdtt_params.corr_offset,      // double      offset);
							    		imgdtt_params.twice_diagonal);  //    		boolean     twice_diagonal)

							    double [] strip_vert = corr2d.combineInterpolatedCorrelations(
							    		strips,                         // double [][] strips,
							    		0x0c,                           // int         pairs_mask,
							    		imgdtt_params.corr_offset,      // double      offset);
							    		imgdtt_params.twice_diagonal);  //    		boolean     twice_diagonal)

							    // re-using arrays that were made for color channels
							    clt_corr_partial[tileY][tileX] = new double[quad][numcol+1][];
						    	clt_corr_partial[tileY][tileX][0][0] = corrs[0];                        // 1
						    	clt_corr_partial[tileY][tileX][0][1] = corrs[1];                        // 2
						    	clt_corr_partial[tileY][tileX][0][2] = corrs[2];                        // 3
						    	clt_corr_partial[tileY][tileX][0][3] = corrs[3];                        // 4
						    	clt_corr_partial[tileY][tileX][1][0] = corrs[4];                        // 5
						    	clt_corr_partial[tileY][tileX][1][1] = corrs[5];                        // 6
						    	clt_corr_partial[tileY][tileX][1][2] = corr2d.debugStrip(strip_hor);    // 7
						    	clt_corr_partial[tileY][tileX][1][3] = corr2d.debugStrip(strip_vert);   // 8
						    	clt_corr_partial[tileY][tileX][2][0] = corr2d.debugStrip(strips[4]);    // 9
						    	clt_corr_partial[tileY][tileX][2][1] = corr2d.debugStrip(strips[5]);    // 10
						    	clt_corr_partial[tileY][tileX][2][2] = corr2d.debugStrip2(strip_hor);   // 11
						    	clt_corr_partial[tileY][tileX][2][3] = corr2d.debugStrip2(strip_vert);  // 12
						    	clt_corr_partial[tileY][tileX][3][0] = corr2d.debugStrip2(strip_ortho); // 13
						    	clt_corr_partial[tileY][tileX][3][1] = corr2d.debugStrip2(strip_diag);  // 14
						    	clt_corr_partial[tileY][tileX][3][2] = corr2d.debugStrip(strip_all);    // 15
						    	clt_corr_partial[tileY][tileX][3][3] = corr2d.debugStrip2(strip_combo); // 16
						    }
						    if ((clt_corr_combo != null) && imgdtt_params.corr_mode_debug) {
						    	// reuse it too?

						    }
						    // calculate CM maximums for all mixed channels
						    // First get integer correlation center, relative to the center
							int [] ixy =  corr2d.getMaxXYInt( // find integer pair or null if below threshold
									strip_combo,              // double [] data,
									true,                     // boolean   axis_only,
									imgdtt_params.min_corr,   //  double    minMax,    // minimal value to consider (at integer location, not interpolated)
									tile_lma_debug_level > 0); // boolean   debug);

							double [] corr_stat = null;

							// if integer argmax was strong enough, calculate CM argmax
							// will not fill out DISPARITY_INDEX_INT+1, DISPARITY_INDEX_CM+1, DISPARITY_INDEX_POLY+1
							// use clt_mismatch for that
							double strength = 0.0;
							double disparity = 0.0;
							if (ixy != null) {
								strength = strip_combo[ixy[0]+transform_size-1]; // strength at integer max on axis
								disparity_map[DISPARITY_INDEX_INT][tIndex] =      -ixy[0];
//								disparity_map[DISPARITY_INDEX_INT + 1][tIndex] =
								disparity_map[DISPARITY_STRENGTH_INDEX][tIndex] = strength;
								if (Double.isNaN(disparity_map[DISPARITY_STRENGTH_INDEX][tIndex])) {
									System.out.println("BUG: 1. disparity_map[DISPARITY_STRENGTH_INDEX]["+tIndex+"] should not be NaN");
								}
								corr_stat = corr2d.getMaxXCm(   // get fractional center as a "center of mass" inside circle/square from the integer max
										strip_combo,                      // double [] data,      // [data_size * data_size]
										ixy[0],                           // int       ixcenter,  // integer center x
										// corr_wndy,                        // double [] window_y,  // (half) window function in y-direction(perpendicular to disparity: for row0  ==1
										// corr_wndx,                        // double [] window_x,  // half of a window function in x (disparity) direction
										(tile_lma_debug_level > 0)); // boolean   debug);
							}

							// for compatibility with old code executed unconditionally. TODO: Move to if (corr_stat != null) ... condition below
							double [] hor_pair1 = corr2d.getMaxXSOrtho(
									corrs,                              // double [][] correlations,
									Correlation2d.getMaskHorizontal(1), // int         pairs_mask,
									imgdtt_params.corr_offset,          // double      corr_offset,
									true,                               // boolean     symmetric,   // for comparing with old implementation average with symmetrical before multiplication
									false,                              // boolean     is_vert,      // transpose X/Y
									tile_lma_debug_level > 0);          // boolean   debug);
							if (hor_pair1 != null) {
								disparity_map[DISPARITY_INDEX_HOR][tIndex] =          -hor_pair1[0];
								disparity_map[DISPARITY_INDEX_HOR_STRENGTH][tIndex] =  hor_pair1[1];
							}

							double [] vert_pair1 = corr2d.getMaxXSOrtho(
									corrs,                              // double [][] correlations,
									Correlation2d.getMaskVertical(1), // int         pairs_mask,
									imgdtt_params.corr_offset,        // double      corr_offset,
									true,                             // boolean     symmetric,   // for comparing with old implementation average with symmetrical before multiplication
									true,                             // boolean     is_vert,      // transpose X/Y
									tile_lma_debug_level > 0); // boolean   debug);
							if (vert_pair1 != null) {
								disparity_map[DISPARITY_INDEX_VERT][tIndex] =         -vert_pair1[0];
								disparity_map[DISPARITY_INDEX_VERT_STRENGTH][tIndex] = vert_pair1[1];
							}


							// proceed only if CM correlation result is non-null // for compatibility with old code we need it to run regardless of the strength of the normal correlation
							if (corr_stat != null) {
// skipping DISPARITY_VARIATIONS_INDEX - it was not used
								disparity = -corr_stat[0];
								disparity_map[DISPARITY_INDEX_CM][tIndex] = disparity; // disparity is negative X
								if (tile_lma_debug_level > 0) {
									System.out.println("Will run getMaxXSOrtho( ) for tileX="+tileX+", tileY="+tileY);
								}
//								disparity_map[DISPARITY_INDEX_CM + 1][tIndex] = // y not available here
								// calculate/fill out hor and vert
								// convert to multi-baseline combining results from several integer scales

								// see if strength is enough to proceed with LMA/poly (otherwise keep disp/strength
								if (strength > imgdtt_params.min_poly_strength) {
									// create LMA instance, calculate LMA composite argmax
							    	// Create 2 groups: ortho & diag
							    	Correlations2dLMA lma = corr2d.corrLMA(
							    			imgdtt_params,                // ImageDttParameters  imgdtt_params,
							    			corrs,                        // double [][]         corrs,
							    			imgdtt_params.dbg_pair_mask,  // int                 pair_mask, // which pairs to process
							        		false,                        // boolean             run_poly_instead, // true - run LMA, false - run 2d polynomial approximation
							    			corr_stat[0],                 // double    xcenter,   // preliminary center x in pixels for largest baseline
							    			imgdtt_params.ortho_vasw_pwr, // double    vasw_pwr,  // value as weight to this power,
							    			tile_lma_debug_level,         // int                 debug_level,
							        		tileX,                        // int                 tileX, // just for debug output
							        		tileY );                      // int                 tileY
							    	double [] lma_disparity_strength = null;
							    	if (lma != null) {
							    		double []   mod_disparity_diff = null;
							    		double [][] dir_corr_strength =  null;
								    	lma_disparity_strength = lma.getDisparityStrength();
							    		if (tile_lma_debug_level > 0){
							    				System.out.println(String.format("Tile X/Y = %d/%d LMA disparity = %7.4f, strength = %7.4f",
							    						tileX, tileY,
							    						lma_disparity_strength[0],lma_disparity_strength[1]));
							    		}
										disparity_map[DISPARITY_INDEX_POLY]         [tIndex] = lma_disparity_strength[0];

										// if enabled overwrite - replace  DISPARITY_INDEX_CM and DISPARITY_STRENGTH_INDEX
										if (imgdtt_params.mix_corr_poly) {
											disparity = lma_disparity_strength[0];
											strength =  lma_disparity_strength[1];
											disparity_map[DISPARITY_INDEX_CM]       [tIndex] = disparity;
											disparity_map[DISPARITY_STRENGTH_INDEX] [tIndex] = strength;
											if (Double.isNaN(disparity_map[DISPARITY_STRENGTH_INDEX][tIndex])) {
												System.out.println("BUG: 2. disparity_map[DISPARITY_STRENGTH_INDEX][tIndex] should not be NaN");
											}
										}
										// store debug data
										// if strong enough and enabled - try to improve far objects
										// temporarily removed strength requirement for debugging, restore later to make faster
///										if ((imgdtt_params.fo_correct && (strength > imgdtt_params.fo_min_strength)) || (clt_mismatch != null)) {

										if ((imgdtt_params.fo_correct && (strength > 0 * imgdtt_params.fo_min_strength)) || (clt_mismatch != null)) {
								    		// try all dirs:
											dir_corr_strength = corr2d.corr4dirsLMA(
								    				imgdtt_params,                // ImageDttParameters  imgdtt_params,
								    				corrs,                        // double [][]         corrs,
								    				imgdtt_params.dbg_pair_mask,  // int                 pair_mask, // which pairs to process
								    				-disparity,                   // double    xcenter,   // preliminary center x in pixels for largest baseline
								    				imgdtt_params.ortho_vasw_pwr, // double    vasw_pwr,  // value as weight to this power,
								    				tile_lma_debug_level,         // int                 debug_level,
								    				tileX,                        // int                 tileX, // just for debug output
								    				tileY );                      // int                 tileY
								    		if ((tile_lma_debug_level > 0) && (dir_corr_strength != null)) {
								    			double [] nan2 = {Double.NaN, Double.NaN, Double.NaN, Double.NaN};
								    			for (int ii = 0; ii < dir_corr_strength.length; ii++) {
								    				if (dir_corr_strength[ii] == null) dir_corr_strength[ii] = nan2;
								    			}
								    			System.out.println(String.format("corr4dirsLMA -> ↔: %7.4f (%7.4f) ↕:%7.4f (%7.4f) ⤡:%7.4f (%7.4f) ⤢:%7.4f (%7.4f)",
								    					dir_corr_strength[0][0],dir_corr_strength[0][1],dir_corr_strength[1][0],dir_corr_strength[1][1],
								    					dir_corr_strength[2][0],dir_corr_strength[2][1],dir_corr_strength[3][0],dir_corr_strength[3][1]));
								    		}

								    		mod_disparity_diff =     corr2d.foregroundCorrect(
								    				imgdtt_params.fo_far,            // boolean   bg,
								    				imgdtt_params.fo_ortho,          // boolean   ortho,
								    				dir_corr_strength,               // double [] dir_disp,
								    	    		disparity,                       // double    full_disp,
								    	    		imgdtt_params.fo_min_strength,   // double      min_strength,
								    	    		imgdtt_params.fo_min_eff,        // double      min_eff,
								    	    		imgdtt_params.fo_min_eff_ratio,  // double      min_eff_ratio,
								    	    		imgdtt_params.fo_max_hwidth,    // double      max_hwidth, //  =          3.0;  // maximal half-width disparity  direction to try to correct far objects
								    	    		imgdtt_params.fo_min_diff,       // double    fo_min_diff,
								    	    		imgdtt_params.fo_overcorrection, // double    fo_overcorrection,
								    	    		imgdtt_params.fo_lim_overcorr,   // double    fo_lim_overcorr,
								    	    		(tile_lma_debug_level > 0) );    // boolean debug);

// Do not use modified far object distance when mismatch is measured
											if ((mod_disparity_diff[0] != disparity) && (clt_mismatch == null)){ // if it changed
												if (imgdtt_params.fo_correct && (strength > imgdtt_params.fo_min_strength)) { // always
													disparity = mod_disparity_diff[0];
													disparity_map[DISPARITY_INDEX_CM]       [tIndex] = disparity;
												}
											}
										}
										if (tile_lma_debug_level > -1) {
											System.out.println("debug12348973591");
										}
										if (clt_mismatch != null) { // mod_disparity_diff should be calculated
											// bypass difference or zero strength if disparity difference is too high (will influence mismatch correction)
											// but setting it too low will make it impossible to correct larger mismatches. Maybe multi-pass?
											if (mod_disparity_diff[2] <= imgdtt_params.mismatch_max_diff) { // may be NaN, will fail test as intended
												if (tile_lma_debug_level > -1) {
													System.out.println("debug12348973590");
												}
												double [] mismatch_result = null;
												boolean need_CM = true;
												if (imgdtt_params.ly_poly) {
													mismatch_result = corr2d.mismatchPairs( // returns x-xcenter, y, strength (sign same as disparity)
															imgdtt_params,                // ImageDttParameters  imgdtt_params,
															corrs,                        // double [][]         corrs,
															all_pairs,                    // int                 pair_mask, // which pairs to process
															-disparity,                   // double    xcenter,   // preliminary center x in pixels for largest baseline
															max_corr_radius,              // double    vasw_pwr,  // value as weight to this power,
															tile_lma_debug_level,// int                 debug_level,
															tileX,         // int                 tileX, // just for debug output
															tileY );       // int                 tileY
													// check if got crazy poly, then retry with CM
													boolean has_NaN = false;
													boolean need_dbg = false;
													for (int dir = 0; dir < (mismatch_result.length/3); dir ++) {
														if (Double.isNaN(mismatch_result[3*dir + 0]) || Double.isNaN(mismatch_result[3*dir + 1])) {
															has_NaN = true;
														} else if ((mismatch_result[3*dir + 2] != 0.0) &&
																((Math.abs(mismatch_result[3*dir + 0]) > imgdtt_params.ly_crazy_poly) ||
																(Math.abs(mismatch_result[3*dir + 1]) > imgdtt_params.ly_crazy_poly))) {
															mismatch_result[3*dir + 2] = 0;
															has_NaN = true;
															need_dbg = true;
														}
													}
													need_CM = imgdtt_params.ly_poly_backup && has_NaN;
													if (need_dbg && (imgdtt_params.lma_debug_level > 0)) {
														System.out.println("Failed polynomial mismatch for tileX="+tileX+", tileY="+tileY);
														for (int dir = 0; dir < (mismatch_result.length/3); dir ++) {
															System.out.println(String.format("%d: dxy[%d]=%f8.5, dxy[%d]=%f8.5 strength=%7.5f",
																	dir, (dir*2+0), mismatch_result[dir*3 + 0], (dir*2 + 1), mismatch_result[dir*3 + 1], mismatch_result[dir*3 + 2]));
														}
													}
												}
												// TODO: use magic_scale for CM?
												if (need_CM) { // if poly was off or gave crazy poly
													mismatch_result = corr2d.mismatchPairsCM( // returns x-xcenter, y, strength (sign same as disparity)
															imgdtt_params,                // ImageDttParameters  imgdtt_params,
															corrs,                        // double [][]         corrs,
															all_pairs,                    // int                 pair_mask, // which pairs to process
															-disparity,                   // double    xcenter,   // preliminary center x in pixels for largest baseline
															max_corr_radius, // imgdtt_params.ortho_vasw_pwr, // radius,    // positive - within that distance, negative - within 2*(-radius)+1 square
															tile_lma_debug_level,// int                 debug_level,
															tileX,         // int                 tileX, // just for debug output
															tileY );       // int                 tileY

													if (imgdtt_params.ly_poly && (imgdtt_params.lma_debug_level > 0)) {
														System.out.println("Corrected by CM failed polynomial mismatch for tileX="+tileX+", tileY="+tileY);
														for (int dir = 0; dir < (mismatch_result.length/3); dir ++) {
															System.out.println(String.format("%d: dxy[%d]=%f8.5, dxy[%d]=%f8.5 strength=%7.5f",
																	dir, (dir*2+0), mismatch_result[dir*3 + 0], (dir*2 + 1), mismatch_result[dir*3 + 1], mismatch_result[dir*3 + 2]));
														}
													}
												}
									    		if (tile_lma_debug_level > 0) {
									    			System.out.println("Lazy eye mismatch:");
									    			for (int np = 0; np < mismatch_result.length/3; np++) {
									    				System.out.println(String.format("%2d: dx = %7.4f, dy = %7.4f, strength = %7.4f,",
									    						np, mismatch_result[3 * np + 0], mismatch_result[3 * np + 1], mismatch_result[3 * np + 2]));
									    			}
									    		}

												for (int np = 0; np < clt_mismatch.length/3; np++) if (np < mismatch_result.length/3){
													clt_mismatch[3 * np + 0 ][tIndex] = mismatch_result[3 * np + 0 ];
													clt_mismatch[3 * np + 1 ][tIndex] = mismatch_result[3 * np + 1 ];
													clt_mismatch[3 * np + 2 ][tIndex] = mismatch_result[3 * np + 2 ];
												}
											}
										}
							    	}
								}
							} // end of if (corr_stat != null)
							if      (corr_mode == 0) extra_disparity = disparity_map[DISPARITY_INDEX_INT][tIndex];
							else if (corr_mode == 1) extra_disparity = disparity_map[DISPARITY_INDEX_CM][tIndex];
							else if (corr_mode == 2) extra_disparity = disparity_map[DISPARITY_INDEX_POLY][tIndex];
							else if (corr_mode == 3) extra_disparity = disparity_map[DISPARITY_INDEX_HOR][tIndex];
							else if (corr_mode == 4) extra_disparity = disparity_map[DISPARITY_INDEX_VERT][tIndex];
							if (Double.isNaN(extra_disparity)) extra_disparity = 0;

							if (Double.isNaN(disparity_map[DISPARITY_STRENGTH_INDEX][tIndex])) {
								System.out.println("BUG: 3. disparity_map[DISPARITY_STRENGTH_INDEX][tIndex] should not be NaN");
							}
						} // if (disparity_map != null){ // not null - calculate correlations
						// only debug is left
						// old (per-color correlation)
						if ((clt_corr_combo != null)  && !imgdtt_params.corr_mode_debug){ // not null - calculate correlations
							tcorr_tpartial=  new double[corr_pairs.length][numcol+1][4][transform_len];
							tcorr_partial =  new double[quad][numcol+1][];

							for (int pair = 0; pair < corr_pairs.length; pair++){
								for (int ncol = 0; ncol <numcol; ncol++) if (clt_data[ncol] != null){
									double [][] data1 = clt_data[corr_pairs[pair][0]][ncol][tileY][tileX];
									double [][] data2 = clt_data[corr_pairs[pair][1]][ncol][tileY][tileX];
									if ((data1 != null) && (data2 != null)) {

										double [] a2 = new double[transform_len];
										double sa2 = 0.0;
										for (int i = 0; i < transform_len; i++) {
											double s1 = 0.0, s2=0.0;
											for (int n = 0; n< 4; n++){
												s1+=data1[n][i] * data1[n][i];
												s2+=data2[n][i] * data2[n][i];
											}
											a2[i] = Math.sqrt(s1*s2);
											sa2 += a2[i];
										}
										double fz2 = sa2/transform_len * corr_fat_zero * corr_fat_zero; // fat_zero squared to match units
										for (int i = 0; i < transform_len; i++) {
											double scale = 1.0 / (a2[i] + fz2);
											for (int n = 0; n<4; n++){
												tcorr_tpartial[pair][ncol][n][i] = 0;
												for (int k=0; k<4; k++){
													if (zi[n][k] < 0)
														tcorr_tpartial[pair][ncol][n][i] -=
														data1[-zi[n][k]][i] * data2[k][i];
													else
														tcorr_tpartial[pair][ncol][n][i] +=
														data1[zi[n][k]][i] * data2[k][i];
												}
												tcorr_tpartial[pair][ncol][n][i] *= scale;
											}
										}
									} else {
										tcorr_tpartial[pair][ncol] = null;
									}
									// got transform-domain correlation for the pair, 1 color
								}
								// calculate composite color
								for (int i = 0; i < transform_len; i++) {
									for (int n = 0; n<4; n++) {
										tcorr_tpartial[pair][numcol][n][i] = 0.0;
										for (int ncol= 0; ncol < tcorr_tpartial[pair].length; ncol++) {
											if (tcorr_tpartial[pair][ncol] != null) {
												tcorr_tpartial[pair][numcol][n][i] += col_weights[ncol] * tcorr_tpartial[pair][0][n][i];
											}
										}
									}
								}
								// now lpf (only last/composite color if do not preserve intermediate
								int firstColor = (clt_corr_partial == null)? numcol : 0;
								if (corr_sigma >0) {
									for (int ncol = firstColor; ncol <= numcol; ncol++) if (tcorr_tpartial[pair][ncol] != null){
										for (int i = 0; i < transform_len; i++) {
											for (int n = 0; n<4; n++) {
												tcorr_tpartial[pair][ncol][n][i] *= filter[i];
											}
										}
									}
								}
								// convert to pixel domain - all or just composite color
								for (int ncol = firstColor; ncol <= numcol; ncol++) if (tcorr_tpartial[pair][ncol] != null) {
									for (int quadrant = 0; quadrant < 4; quadrant++){
										int mode = ((quadrant << 1) & 2) | ((quadrant >> 1) & 1); // transpose
										tcorr_tpartial[pair][ncol][quadrant] =
												dtt.dttt_iie(tcorr_tpartial[pair][ncol][quadrant], mode, transform_size);
									}
								}
								// convert from 4 quadrants to 15x15 centered tiles (each color or only composite)
								for (int ncol = firstColor; ncol <= numcol; ncol++) if (tcorr_tpartial[pair][ncol] != null) {
									tcorr_partial[pair][ncol] = dtt.corr_unfold_tile(
											tcorr_tpartial[pair][ncol],
											transform_size);
								}
								// transpose vertical pairs
								if (corr_pairs[pair][2] != 0) {
									for (int ncol = firstColor; ncol <= numcol; ncol++) if (tcorr_tpartial[pair][ncol] != null) {
										for (int i = 0; i < transpose_indices.length; i++) {
											double d = tcorr_partial[pair][ncol][transpose_indices[i][0]];
											tcorr_partial[pair][ncol][transpose_indices[i][0]] = tcorr_partial[pair][ncol][transpose_indices[i][1]];
											tcorr_partial[pair][ncol][transpose_indices[i][1]] = d;
											//transpose_indices
										}
									}
								}
								// make symmetrical around the disparity direction (horizontal) (here using just average, not mul/sum mixture)
								// symmetry can be added to result, not individual (if sum - yes, but with multiplication - not)
								if (corr_sym && (clt_mismatch == null)){ // when measuring clt_mismatch symmetry should be off !
									for (int ncol = firstColor; ncol <= numcol; ncol++) if (tcorr_tpartial[pair][ncol] != null) {
										for (int i = 1 ; i < transform_size; i++){
											int indx1 = (transform_size - 1 - i) * corr_size;
											int indx2 = (transform_size - 1 + i) * corr_size;
											for (int j = 0; j< corr_size; j++){
												int indx1j = indx1 + j;
												int indx2j = indx2 + j;
												tcorr_partial[pair][ncol][indx1j] =
														0.5* (tcorr_partial[pair][ncol][indx1j] + tcorr_partial[pair][ncol][indx2j]);
												tcorr_partial[pair][ncol][indx2j] = tcorr_partial[pair][ncol][indx1j];
											}
										}
									}
								}
							} // all pairs calculated


							tcorr_combo = new double [TCORR_TITLES.length][corr_size * corr_size];

							int numPairs = 	0, numPairsHor = 0, numPairsVert = 0;
							for (int pair = 0; pair < corr_pairs.length; pair++) if (((corr_mask >> pair) & 1) != 0){
								numPairs++;
								if (corr_pairs[pair][2] == 0) { // horizontal pair)
									numPairsHor++;
								} else {
									numPairsVert++;
								}
							}
							double avScale = 0.0, avScaleHor = 0.0, avScaleVert = 0.0;
							if (numPairs > 0) {
								boolean debugMax =  tile_lma_debug_level > 1;
								avScale = 1.0/numPairs;
								if (numPairsHor > 0)  avScaleHor = 1.0/numPairsHor;
								if (numPairsVert > 0) avScaleVert = 1.0/numPairsVert;
								if (debugMax) {
									System.out.println("avScale = "+avScale+", avScaleHor = "+avScaleHor+", avScaleVert = "+avScaleVert+", corr_offset = "+corr_offset);
								}
								if (corr_offset < 0) { // just add all partial correlations for composite color
									for (int i = 0; i < tcorr_combo[TCORR_COMBO_RSLT].length; i++){
										tcorr_combo[TCORR_COMBO_RSLT][i] = 0.0;
										tcorr_combo[TCORR_COMBO_HOR][i] = 0.0;
										tcorr_combo[TCORR_COMBO_VERT][i] = 0.0;
										for (int pair = 0; pair < corr_pairs.length; pair++) if (((corr_mask >> pair) & 1) != 0){
											tcorr_combo[TCORR_COMBO_RSLT][i] += avScale*tcorr_partial[pair][numcol][i]; // only composite color channel
											if (corr_pairs[pair][2] == 0) { // horizontal pair
												tcorr_combo[TCORR_COMBO_HOR][i] +=  avScaleHor*tcorr_partial[pair][numcol][i]; // only composite color channel
											} else { //vertical pair
												tcorr_combo[TCORR_COMBO_VERT][i] += avScaleVert*tcorr_partial[pair][numcol][i]; // only composite color channel
											}
											if (debugMax) {
												System.out.println("tcorr_combo[TCORR_COMBO_RSLT]["+i+"]="+tcorr_combo[TCORR_COMBO_RSLT][i]+" tcorr_partial["+pair+"]["+numcol+"]["+i+"]="+tcorr_partial[pair][numcol][i]);
											}
										}
									}
								} else {
									for (int i = 0; i < tcorr_combo[TCORR_COMBO_RSLT].length; i++){
										tcorr_combo[TCORR_COMBO_RSLT][i] = 1.0;
										tcorr_combo[TCORR_COMBO_HOR][i] =  1.0;
										tcorr_combo[TCORR_COMBO_VERT][i] = 1.0;
										for (int pair = 0; pair < corr_pairs.length; pair++) if (((corr_mask >> pair) & 1) != 0){
											tcorr_combo[TCORR_COMBO_RSLT][i] *= (tcorr_partial[pair][numcol][i] + corr_offset); // only composite color channel
											if (corr_pairs[pair][2] == 0) { // horizontal pair
												tcorr_combo[TCORR_COMBO_HOR][i] *= (tcorr_partial[pair][numcol][i] + corr_offset); // only composite color channel
											} else { //vertical pair
												tcorr_combo[TCORR_COMBO_VERT][i] *= (tcorr_partial[pair][numcol][i] + corr_offset); // only composite color channel
											}
											if (debugMax) {
												System.out.println("tcorr_combo[TCORR_COMBO_RSLT]["+i+"]="+tcorr_combo[TCORR_COMBO_RSLT][i]+" tcorr_partial["+pair+"]["+numcol+"]["+i+"]="+tcorr_partial[pair][numcol][i]);
											}
										}
										if (corr_normalize) {
											if (tcorr_combo[TCORR_COMBO_RSLT][i] > 0.0){
												tcorr_combo[TCORR_COMBO_RSLT][i] = Math.pow(tcorr_combo[TCORR_COMBO_RSLT][i],avScale) - corr_offset;
											} else {
												tcorr_combo[TCORR_COMBO_RSLT][i] =  -corr_offset;
											}

											if (tcorr_combo[TCORR_COMBO_HOR][i] > 0.0){
												tcorr_combo[TCORR_COMBO_HOR][i] = Math.pow(tcorr_combo[TCORR_COMBO_HOR][i],avScaleHor) - corr_offset;
											} else {
												tcorr_combo[TCORR_COMBO_HOR][i] =  -corr_offset;
											}

											if (tcorr_combo[TCORR_COMBO_VERT][i] > 0.0){
												tcorr_combo[TCORR_COMBO_VERT][i] = Math.pow(tcorr_combo[TCORR_COMBO_VERT][i],avScaleVert) - corr_offset;
											} else {
												tcorr_combo[TCORR_COMBO_VERT][i] =  -corr_offset;
											}
										}
									}
								}
								// calculate sum also
								for (int i = 0; i < tcorr_combo[TCORR_COMBO_SUM].length; i++){
									tcorr_combo[TCORR_COMBO_SUM][i] = 0.0;
									for (int pair = 0; pair < corr_pairs.length; pair++) if (((corr_mask >> pair) & 1) != 0){
										tcorr_combo[TCORR_COMBO_SUM][i] += avScale*tcorr_partial[pair][numcol][i]; // only composite color channel
										if (debugMax) {
											System.out.println("tcorr_combo[TCORR_COMBO_SUM]["+i+"]="+tcorr_combo[TCORR_COMBO_SUM][i]+" tcorr_partial["+pair+"]["+numcol+"]["+i+"]="+tcorr_partial[pair][numcol][i]);
										}
									}
								}
								// return results
								for (int n = 0; n < clt_corr_combo.length; n++){ // tcorr_combo now may be longer than clt_corr_combo
									clt_corr_combo[n][tileY][tileX] = tcorr_combo[n];
								}
								if (clt_corr_partial != null){
									clt_corr_partial[tileY][tileX] = tcorr_partial;
								}

							} // if (numPairs > 0) {

						} // end of if (clt_corr_combo != null)




						if (texture_tiles !=null) {
							if ((extra_disparity != 0) && !getForcedDisparity(tile_op[tileY][tileX])){ // 0 - adjust disparity, 1 - use provided
								// shift images by 0.5 * extra disparity in the diagonal direction
								for (int ncol = 0; ncol <numcol; ncol++) { // color
									for (int i = 0; i < quad; i++) {
										if (clt_data[i][ncol] != null) {
											fract_shift(    // fractional shift in transform domain. Currently uses sin/cos - change to tables with 2? rotations
													clt_data[i][ncol][tileY][tileX], // double  [][]  clt_tile,
													transform_size,
													extra_disparity * port_offsets[i][0] / corr_magic_scale,     // double        shiftX,
													extra_disparity * port_offsets[i][1] / corr_magic_scale,     // double        shiftY,
													//									(globalDebugLevel > 0) && (tileX == debug_tileX) && (tileY == debug_tileY)); // external tile compare
													((globalDebugLevel > 0) && (ncol==0) && (tileX >= debug_tileX - 2) && (tileX <= debug_tileX + 2) &&
															(tileY >= debug_tileY - 2) && (tileY <= debug_tileY+2)));
										}
									}
								}
							}
							// lpf tiles (same as images before)
							// iclt tiles
							double [][][] iclt_tile = new double [quad][numcol][]; // in mono some may remain null
							double [] clt_tile;
							double scale = 0.25;  // matching iclt_2d
							for (int i = 0; i < quad; i++) {
								for (int ncol = 0; ncol <numcol; ncol++) if (clt_data[i][ncol] != null) { // color
									// double [] clt_tile = new double [transform_size*transform_size];
									for (int dct_mode = 0; dct_mode < 4; dct_mode++){
										clt_tile = clt_data[i][ncol][tileY][tileX][dct_mode].clone();
										// lpf each of the 4 quadrants before idct
										for (int j = 0; j < filter.length; j++){
											clt_tile[j] *= scale*filter[j];
										}
										// IDCT-IV should be in reversed order: CC->CC, SC->CS, CS->SC, SS->SS
										int idct_mode = ((dct_mode << 1) & 2) | ((dct_mode >> 1) & 1);

										if ((globalDebugLevel > 0) && debugTile) {
											double [] clt_tile_dbg =  clt_tile.clone();
											double [] clt_tile_dbg1 = clt_tile.clone();
											clt_tile_dbg = dtt.dttt_ivn   (clt_tile_dbg,  idct_mode, transform_size, false);
											clt_tile_dbg1 = dtt.dttt_ivn  (clt_tile_dbg1, idct_mode, transform_size, true);
											clt_tile = dtt.dttt_iv  (clt_tile, idct_mode, transform_size);
											System.out.println("\n-------- tileX="+tileX+", tileY="+tileY+", idct_mode="+idct_mode+" ---#, standard, good, bad, diff---");
											for (int nt = 0; nt < clt_tile_dbg.length; nt++) {
												System.out.println(String.format("%2d: %8f %8f %8f %8f", //  %8f",
														clt_tile[nt],
														clt_tile_dbg[nt],
														clt_tile_dbg1[nt],
														clt_tile_dbg[nt] - clt_tile_dbg[nt],
														clt_tile_dbg1[nt] - clt_tile_dbg[nt]));
											}
										} else {
											clt_tile = dtt.dttt_iv  (clt_tile, idct_mode, transform_size);
										}

										// iclt_tile[i][chn] = dtt.dttt_iv  (clt_data[i][chn][tileY][tileX][dct_mode], idct_mode, transform_size);
										double [] tile_mdct = dtt.unfold_tile(clt_tile, transform_size, dct_mode); // mode=0 - DCCT 16x16
										// accumulate partial mdct results
										if (dct_mode == 0){
											iclt_tile[i][ncol] = tile_mdct;
										} else{
											for (int j = 0; j<tile_mdct.length; j++){
												iclt_tile[i][ncol][j] += tile_mdct[j]; // matching iclt_2d
											}
										}
									}
								}
							}
							// iclt here: [quad][color][256]
							if ((globalDebugLevel > 0) && debugTile) {
								ShowDoubleFloatArrays sdfa_instance = new ShowDoubleFloatArrays(); // just for debugging?
								String [] titles = {"red0","blue0","green0","red1","blue1","green1","red2","blue2","green2","red3","blue3","green3"};
								double [][] dbg_tile = new double [quad*numcol][];
								for (int i = 0; i < quad; i++) {
									for (int ncol = 0; ncol <numcol; ncol++) if (iclt_tile[i][ncol] != null) { // color
										dbg_tile[i * numcol + ncol] = iclt_tile[i][ncol];
									}
								}
								sdfa_instance.showArrays(dbg_tile, 2* transform_size, 2* transform_size, true, "iclt_x"+tileX+"_y"+tileY, titles);
							}


							// "de-bayer" tiles for matching, use original data for output
							double [][][] tiles_debayered = new double [quad][numcol][];
							for (int i =0; i<quad; i++){
								for (int ncol = 0; ncol < numcol; ncol++) if (iclt_tile[i][ncol] != null) {
									if (isMonochrome()) {
										tiles_debayered[i][ncol] =  iclt_tile[i][ncol];
									} else {
										tiles_debayered[i][ncol] =  tile_debayer_shot_corr(
												(ncol != 2), // red or blue (false - green)
												iclt_tile[i][ncol],
												2 * transform_size,
												lt_window2, // squared lapping window
												min_shot,   // 10.0;  // Do not adjust for shot noise if lower than
												scale_shot,  //3.0;   // scale when dividing by sqrt
												lt_window2); // re-apply window to the result
									}
								}
							}
							if ((globalDebugLevel > 0) && debugTile) {
								ShowDoubleFloatArrays sdfa_instance = new ShowDoubleFloatArrays(); // just for debugging?
								String [] titles = {"red0","blue0","green0","red1","blue1","green1","red2","blue2","green2","red3","blue3","green3"};
								double [][] dbg_tile = new double [quad*numcol][];
								for (int i = 0; i < quad; i++) {
									for (int chn = 0; chn <numcol; chn++) { // color
										dbg_tile[i * numcol + chn] = tiles_debayered[i][chn];
									}
								}
								sdfa_instance.showArrays(dbg_tile, 2* transform_size, 2* transform_size, true, "tiles_debayered_x"+tileX+"_y"+tileY, titles);
							}

							double []     max_diff = null;
							if ((disparity_map != null) && (disparity_map.length >= (IMG_DIFF0_INDEX + quad))){
								max_diff = new double[quad];
							}
							int ports_rgb_len = quad*numcol;  // 12
							if ((disparity_map != null) && (disparity_map.length >= (IMG_TONE_RGB + ports_rgb_len))) {
								ports_rgb = new double[ports_rgb_len];
							}
							texture_tiles[tileY][tileX] =  tile_combine_rgba(
									tiles_debayered, // iclt_tile,      // [port][numcol][256]
									ports_rgb, // double []     ports_rgb,      // average values of R,G,B for each camera (R0,R1,...,B2,B3)
									max_diff,        // maximal (weighted) deviation of each channel from the average
									lt_window2,      // [256]
									port_offsets,    // [port]{x_off, y_off}
									img_mask,        // which port to use, 0xf - all 4 (will modify as local variable)
									diff_sigma,      // pixel value/pixel change
									diff_threshold,  // pixel value/pixel change
									diff_gauss,      // when averaging images, use gaussian around average as weight (false - sharp all/nothing)
									min_agree,       // minimal number of channels to agree on a point (real number to work with fuzzy averages)
									col_weights,     // color channel weights, sum == 1.0 (already set to 0/0/1 for monochrome
									dust_remove,     // boolean dust_remove,    // Do not reduce average weight when only one image differes much from the average
									keep_weights,    // keep_weights);   // return channel weights after A in RGBA
									(globalDebugLevel > 0) && debugTile);

							// mix RGB from iclt_tile, mix alpha with - what? correlation strength or 'don't care'? good correlation or all > min?
							for (int i = 0; i < iclt_tile[0][first_color].length; i++ ) {
								double sw = 0.0;
								for (int ip = 0; ip < quad; ip++) {
									sw += texture_tiles[tileY][tileX][numcol+1+ip][i];
								}
								if (sw != 0 ) sw = 1.0/sw;
								for (int ncol = 0; ncol < numcol; ncol++) if (iclt_tile[0][ncol] !=null){ // color
									texture_tiles[tileY][tileX][ncol][i] = 0.0; //iclt[tileY][tileX][chn]
									for (int ip = 0; ip < quad; ip++) {
										texture_tiles[tileY][tileX][ncol][i] += sw * texture_tiles[tileY][tileX][numcol+1+ip][i] * iclt_tile[ip][ncol][i];
									}
								}
							}
							if ((disparity_map != null) && (disparity_map.length >= (IMG_DIFF0_INDEX + quad))){
								for (int i = 0; i < max_diff.length; i++){
									disparity_map[IMG_DIFF0_INDEX + i][tIndex] = max_diff[i];
								}
							}
							if (ports_rgb != null) {
								for (int i = 0; i < ports_rgb.length; i++){
									disparity_map[IMG_TONE_RGB + i][tIndex] = ports_rgb[i];
								}
							}
						}
					}
				}
			};
		}
		startAndJoin(threads);
/*
		if (dbg_ports_coords != null) {
			(new showDoubleFloatArrays()).showArrays(dbg_ports_coords,  tilesX, tilesY, true, "ports_coordinates", dbg_titles);
		}
*/
		return clt_data;
	}

	public boolean dmExists(double [][] dm, int indx) {
		return (dm != null) && (dm.length > indx) && (dm[indx]!= null);
	}

	public double [][] tile_combine_rgba(
			double [][][] iclt_tile,     // [port][numcol][256] // in mono some are null
			double []     ports_rgb,     // average values of R,G,B for each camera (R0,R1,...,B2,B3)
			double []     max_diff,      // maximal (weighted) deviation of each channel from the average
			double []     lt_window,     // [256]
			double [][]   port_offsets,  // [port]{x_off, y_off} - just to scale pixel value differences
			int           port_mask,      // which port to use, 0xf - all 4 (will modify as local variable)
			double        diff_sigma,     // pixel value/pixel change
			double        diff_threshold, // pixel value/pixel change
			boolean       diff_gauss,     // when averaging images, use gaussian around average as weight (false - sharp all/nothing)
			double        min_agree,      // minimal number of channels to agree on a point (real number to work with fuzzy averages)
			double []     chn_weights,     // color channel weights, sum == 1.0
			boolean       dust_remove,    // Do not reduce average weight when only one image differes much from the average
			boolean       keep_weights,   // return channel weights after A in RGBA
			boolean       debug)
	{
		int ports =  iclt_tile.length;

		int numcol =  iclt_tile[0].length;
		int tile_len = 0;
		for (int ncol = 0; ncol < numcol; ncol++) {
			if (iclt_tile[0][ncol] != null) {
				tile_len = iclt_tile[0][ncol].length;
				break;
			}
		}

		int usedPorts = ((port_mask >> 0) & 1) + ((port_mask >> 1) & 1) + ((port_mask >> 2) & 1) + ((port_mask >> 3) & 1);


		double [][] port_weights = new double[ports][tile_len];
		double [][] color_avg =    new double[numcol][tile_len];
		double [][] rgba = new double[numcol + 1 + (keep_weights?(ports + 4):0)][];
		int rms_start = numcol + 1 + ports;
		if (keep_weights){
			for (int ncol = 0; ncol <= numcol ; ncol++) if ((ncol == numcol) || (iclt_tile[0][ncol] != null)) {
				rgba[rms_start + ncol] = new double [tile_len]; // rms for each color, then - weighted
			}
			for (int i = 0; i <tile_len; i++){
				double sw = 0.0;
				for (int ncol = 0; ncol < numcol; ncol ++ )  if (iclt_tile[0][ncol] != null){
					double s0 = 0, s1 = 0, s2 = 0;
					for (int ip = 0; ip < ports; ip++)if ((port_mask & ( 1 << ip)) != 0){
						s0 += 1;
						s1 += iclt_tile[ip][ncol][i];
						s2 += iclt_tile[ip][ncol][i]*iclt_tile[ip][ncol][i];
					}
					rgba[rms_start+ncol][i] = Math.sqrt(s0*s2 - s1*s1) / s0;
					sw += chn_weights[ncol]*rgba[rms_start+ncol][i]*rgba[rms_start+ncol][i];
				}
				rgba[rms_start+numcol][i] = Math.sqrt(sw); // will fade as window
			}
		}


		double []  alpha = new double[tile_len];
		double threshold2 = diff_sigma * diff_threshold;
		threshold2 *= threshold2; // squared to compare with diff^2
		if (usedPorts > 1) {
			double [] pair_dist2r =  new double [ports*(ports-1)/2]; // reversed squared distance between images - to be used with gaussian
			int [][]  pair_ports = new int [ports*(ports-1)/2][2];
			int indx = 0;
			double ksigma = 1.0/(2.0*diff_sigma*diff_sigma); // multiply by a weighted sum of squares of the differences
			for (int i = 0; i < ports; i++) if ((port_mask & ( 1 << i)) != 0){
				for (int j = i+1; j < ports; j++)  if ((port_mask & ( 1 << j)) != 0){
					double dx = port_offsets[j][0] - port_offsets[i][0];
					double dy = port_offsets[j][1] - port_offsets[i][1];
					pair_ports[indx][0] = i;
					pair_ports[indx][1] = j;
					pair_dist2r[indx++] = ksigma/(dx*dx+dy*dy); // 2*sigma^2 * r^2
				}
			}
			// there will be no pairs for a single used port
			for (int i = 0; i < tile_len; i++){
				for (int ip = 0; ip < ports; ip++) port_weights[ip][i] = 0.0;
				for (int ip = 0; ip<pair_ports.length; ip++){
					double d = 0;
					for (int ncol = 0; ncol < numcol; ncol++) if (iclt_tile[0][ncol] != null){
						double dc = iclt_tile[pair_ports[ip][0]][ncol][i] - iclt_tile[pair_ports[ip][1]][ncol][i];
						dc /= lt_window[i]; // to compensate fading near the edges
						d+= chn_weights[ncol]*dc*dc;

					}
					d = Math.exp(-pair_dist2r[ip]*d); // 0.5 for exact match, lower for mismatch. Add this weight to both ports involved
					// Add weight to both channels in a pair
					port_weights[pair_ports[ip][0]][i] +=d;
					port_weights[pair_ports[ip][1]][i] +=d;
				}
				// find 2 best ports (resolving 2 pairs of close values)
				int bestPort1=0;
				for (int ip = bestPort1+1; ip < ports; ip++) if (port_weights[ip][i] > port_weights[bestPort1][i]) bestPort1 = ip;
				int bestPort2 = (bestPort1 == 0)?1:0;
				for (int ip = bestPort2+1; ip < ports; ip++) if ((ip != bestPort1) && (port_weights[ip][i] > port_weights[bestPort2][i])) bestPort2 = ip;
				// find weighted average between these 2 ports
				double w1 = port_weights[bestPort1][i]/(port_weights[bestPort1][i]+port_weights[bestPort2][i]);
				double w2 = 1.0 - w1;
				for (int ncol = 0; ncol < numcol; ncol++) if (iclt_tile[0][ncol] != null) {
					color_avg[ncol][i] = w1 * iclt_tile[bestPort1][ncol][i] + w2 * iclt_tile[bestPort2][ncol][i];
				}
				// recalculate all weights using difference from this average of the best pair
				double [] d2 = new double [ports]; //weighted squared differences
				for (int ip = 0; ip < ports; ip++) if ((port_mask & ( 1 << ip)) != 0){
					d2[ip] = 0;
					for (int ncol = 0; ncol < numcol; ncol++)  if (iclt_tile[0][ncol] != null){
						double dc = iclt_tile[ip][ncol][i] - color_avg[ncol][i];
						dc /= lt_window[i]; // to compensate fading near the edges
						d2[ip]+= chn_weights[ncol]*dc*dc;
					}
					port_weights[ip][i] = Math.exp(-ksigma * d2[ip]);
				}
				// and now make a new average with those weights
				double k = 0.0;
				for (int ip = 0; ip < ports; ip++) k+=port_weights[ip][i];
				k = 1.0/k;
				for (int ncol = 0; ncol < numcol; ncol++) if (iclt_tile[0][ncol] != null) {
					color_avg[ncol][i] = 0;
					for (int ip = 0; ip < ports; ip++) {
						color_avg[ncol][i] += k* port_weights[ip][i] * iclt_tile[ip][ncol][i];
					}
				}

			} // or (int i = 0; i < tile_len; i++){

		} else if (usedPorts > 0){ // just copy from a single channel
			for (int ip = 0; ip < ports; ip++) if ((port_mask & ( 1 << ip)) != 0){
				for (int i = 0; i < tile_len; i++){
					for (int ncol = 0; ncol < numcol; ncol++)  if (iclt_tile[0][ncol] != null){
						color_avg[ncol][i] = iclt_tile[ip][ncol][i];
					}
					port_weights[ip][i] = 1.0; // lt_window[i]; // or use 1.0?
				}
			}
		}
		if (dust_remove && (usedPorts == 4)) {
			dust_remove(port_weights);
		}
		// calculate alpha from channel weights. Start with just a sum of weights?
		for (int i = 0; i < tile_len; i++){
			alpha[i] = 0.0;
			for (int ip = 0; ip < ports; ip++) if ((port_mask & ( 1 << ip)) != 0){
				alpha[i]+=	port_weights[ip][i];
			}
			alpha[i] *= lt_window[i]/usedPorts; // make it configurable?
		}

		for (int ncol = 0; ncol < numcol; ncol++) {
			rgba[ncol] = color_avg[ncol];
		}
		rgba[numcol] = alpha;

		for (int i = 0; i < ports; i++)  rgba[numcol + 1 + i] = port_weights[i];
		if (max_diff != null){
			for (int ip = 0; ip < ports; ip++){
				max_diff[ip] = 0;
				if ((port_mask & ( 1 << ip)) != 0) {
					for (int i = 0; i < tile_len; i++){
						double d2 = 0.0;
						for (int ncol = 0; ncol < numcol; ncol++) if (iclt_tile[0][ncol] != null) {
							double dc = (iclt_tile[ip][ncol][i] - color_avg[ncol][i]);
							d2+=dc*dc*chn_weights[ncol];
						}
						d2 *=lt_window[i];
						if (d2 > max_diff[ip]) max_diff[ip]  = d2;
					}
				}
				max_diff[ip] = Math.sqrt(max_diff[ip]);
			}
		}
		if (ports_rgb != null) {
			for (int ip = 0; ip < ports; ip++){
				for (int ncol = 0; ncol < numcol; ncol++) if (iclt_tile[0][ncol] != null) {
					int indx = ncol*ports+ip;
					ports_rgb[indx]=0;
					for (int i = 0; i < tile_len; i++){
						ports_rgb[indx] += iclt_tile[ip][ncol][i];
					}
					ports_rgb[indx] /= tile_len;
				}
			}
		}
		return rgba;
	}

	public void dust_remove( // redistribute weight between 3 best ports (use only when all 3 are enabled)
			double [][]  port_weights)
	{
		int np = port_weights.length;
		for (int i = 0; i < port_weights[0].length; i++){

			int wi = 0;
			for (int ip = 1; ip < np; ip++) if (port_weights[ip][i] < port_weights[wi][i]) wi = ip;
			double avg = 0;
			for (int ip = 1; ip < np; ip++) if (ip != wi) avg += port_weights[ip][i];
			avg /= (np -1);
			double scale = 1.0 + (avg - port_weights[wi][i])/(avg * (np -1));
			for (int ip = 1; ip < np; ip++) {
				if (ip != wi) port_weights[ip][i] *= scale; // increase weight of non-worst, so if worst == 0.0 sum of 3 (all) ports will be scaled by 4/3, keeping average
			}
			port_weights[wi][i] *= port_weights[wi][i]/avg;
		}
	}


	public double [] tile_debayer_shot_corr(
			boolean   rb,
			double [] tile,
			int tile_size,
			double [] window2, // squared lapping window
			double min_shot,   // 10.0;  // Do not adjust for shot noise if lower than
			double scale_shot,  //3.0;   // scale when dividing by sqrt
			double [] window_back) // re-apply window to the result
	{
		double [] tile_nw = new double [tile.length];
		for (int i = 0; i < tile.length; i++) tile_nw[i] = tile[i]/window2[i]; //unapply squared window
		double [] tile_db = tile_debayer(
				rb,
				tile_nw,
				tile_size);
		if (scale_shot > 0){
			double k = 1.0/Math.sqrt(min_shot);
			for (int i = 0; i < tile.length; i++) tile_db[i] = scale_shot* ((tile_db[i] > min_shot)? Math.sqrt(tile_db[i]) : (k*tile_db[i]));
		}
		if (window_back != null) {
			for (int i = 0; i < tile.length; i++) tile_db[i] = tile_db[i] * window_back[i]; // optionally re-apply window (may be a different one)
		}
		return tile_db;
	}


	public double [] tile_debayer(
			boolean   rb,
			double [] tile,
			int tile_size)
	{
		int [] neib_indices = {-tile_size - 1, -tile_size, -tile_size + 1, -1, 0, 1, tile_size - 1, tile_size, tile_size + 1};
		int tsm1 = tile_size - 1;
		double [] rslt = new double [tile_size*tile_size]; // assuming cleared to 0.0;
		double [] kern = rb ? kern_rb : kern_g;
		double k_corn = rb? (16.0/9.0):(4.0/3.0);
		double k_side = rb? (4.0/3.0):(8.0/7.0);
		// top left
		int indx = 0;
		int side_type = 0;
		for (int dri = 0; dri < corn_side_indices[side_type].length; dri++) {
			int dr = corn_side_indices[side_type][dri]; // middle left
			rslt[indx]+=tile[indx+neib_indices[dr]]*kern[dr];
		}
		rslt[indx] *= k_corn;
		// top middle
		side_type = 1;
		for (int j = 1; j < tsm1; j++){
			indx = j;
			for (int dri = 0; dri < corn_side_indices[side_type].length; dri++) {
				int dr = corn_side_indices[side_type][dri]; // middle left
				rslt[indx]+=tile[indx + neib_indices[dr]] * kern[dr];
			}
			rslt[indx] *= k_side;
		}
		// top right
		indx = tsm1;
		side_type = 2;
		for (int dri = 0; dri < corn_side_indices[side_type].length; dri++) {
			int dr = corn_side_indices[side_type][dri]; // middle left
			rslt[indx]+=tile[indx+neib_indices[dr]]*kern[dr];
		}
		rslt[indx] *= k_corn;
		// middle left
		side_type = 3;
		for (int i = 1; i < tsm1; i++){
			indx = i* tile_size;
			for (int dri = 0; dri < corn_side_indices[side_type].length; dri++) {
				int dr = corn_side_indices[side_type][dri]; // middle left
				rslt[indx]+=tile[indx + neib_indices[dr]] * kern[dr];
			}
			rslt[indx] *= k_side;
		}
		// middle middle
		side_type = 4;
		for (int i = 1; i < tsm1; i++){
			for (int j = 1; j < tsm1; j++){
				indx = i*tile_size+j;
				rslt[indx] = 0.0;
				for (int dr = 0; dr < neib_indices.length; dr++){
					rslt[indx]+=tile[indx+neib_indices[dr]]*kern[dr];
				}
			}
		}

		// middle right
		side_type = 5;
		for (int i = 1; i < tsm1; i++){
			indx = i* tile_size + tsm1;
			for (int dri = 0; dri < corn_side_indices[side_type].length; dri++) {
				int dr = corn_side_indices[side_type][dri]; // middle left
				rslt[indx]+=tile[indx + neib_indices[dr]] * kern[dr];
			}
			rslt[indx] *= k_side;
		}

		// bottom left
		indx = tsm1*tile_size;
		side_type = 6;
		for (int dri = 0; dri < corn_side_indices[side_type].length; dri++) {
			int dr = corn_side_indices[side_type][dri]; // middle left
			rslt[indx]+=tile[indx+neib_indices[dr]]*kern[dr];
		}
		rslt[indx] *= k_corn;
		// bottom middle
		side_type = 7;
//		tsm1*tile_size;
		for (int j = 1; j < tsm1; j++){
			indx++;
			for (int dri = 0; dri < corn_side_indices[side_type].length; dri++) {
				int dr = corn_side_indices[side_type][dri]; // middle left
				rslt[indx]+=tile[indx + neib_indices[dr]] * kern[dr];
			}
			rslt[indx] *= k_side;
		}

		// bottom right
		indx++; // = tile_size*tile_size-1;
		side_type = 8;
		for (int dri = 0; dri < corn_side_indices[side_type].length; dri++) {
			int dr = corn_side_indices[side_type][dri]; // middle left
			rslt[indx]+=tile[indx+neib_indices[dr]]*kern[dr];
		}
		rslt[indx] *= k_corn;
		return rslt;
	}
	//final int    []   neib_indices = {-width-1,-width,-width+1,-1,0,1,width-1,width,width+1};


	// return weights for positive x,y, [(radius+a)*(radius+1)]
	public double [] setMaxXYWeights(
			double sigma,
			int    radius){ // ==3.0, ignore data outside sigma * nSigma
			 //
		double [] weights = new double [(radius + 1)*(radius + 1)];
		int indx = 0;
		for (int i = 0; i <= radius; i ++){
			for (int j = 0; j <= radius; j ++){
				weights[indx++] = Math.exp(-(i*i+j*j)/(2*sigma*sigma));
			}
		}
		return weights;
	}

	// find interpolated location of maximum, return {x,y} or null (if too low or non-existing)

	public int [] getMaxXYInt( // find integer pair or null if below threshold
			double [] data,      // [data_size * data_size]
			int       data_size,
			double    minMax,    // minimal value to consider (at integer location, not interpolated)
			boolean   debug)
	{
		int    imx = 0;
		for (int i = 1; i < data.length; i++){
			if (data[imx] < data[i]){
				imx = i;
			}
		}
		if (data[imx] < minMax){
			if (debug){
					System.out.println("getMaxXYInt() -> null (data["+imx+"] = "+data[imx]+" < "+minMax);
			}
			return null;
		}
		int [] rslt = {imx %  data_size, imx /  data_size};
		if (debug){
			System.out.println("getMaxXYInt() -> "+rslt[0]+"/"+rslt[1]);
		}
		return rslt;
	}

	public double [] getMaxXYCm( // get fractiona center as a "center of mass" inside circle/square from the integer max
			double [] data,      // [data_size * data_size]
			int       data_size,
			int []    icenter, // integer center coordinates (relative to top left)
			double    radius,  // positive - within that distance, negative - within 2*(-radius)+1 square
			boolean   max_corr_double,
			boolean   debug)
	{
		if (icenter == null) {
			double [] rslt = {Double.NaN,Double.NaN};
			return rslt; //gigo
		}
		//calculate as "center of mass"
		int iradius = (int) Math.abs(radius);
		int ir2 = (int) (radius*radius);
		boolean square = radius <0;
		double s0 = 0, sx=0,sy = 0;
		for (int y = - iradius ; y <= iradius; y++){
			int dataY = icenter[1] +y;
			if ((dataY >= 0) && (dataY < data_size)){
				int y2 = y*y;
				for (int x = - iradius ; x <= iradius; x++){
					int dataX = icenter[0] +x;
					double r2 = y2 + x * x;
//					if ((dataX >= 0) && (dataX < data_size) && (square || ((y2 + x * x) <= ir2))){
					if ((dataX >= 0) && (dataX < data_size) && (square || (r2 <= ir2))){
//						double w = max_corr_double? (1.0 - r2/ir2):1.0;
//						double d =  w* data[dataY * data_size + dataX];
						double d =  data[dataY * data_size + dataX];
						s0 += d;
						sx += d * dataX;
						sy += d * dataY;
					}
				}
			}
		}
		double [] rslt = {sx / s0, sy / s0};
		if (debug){
			System.out.println("getMaxXYInt() -> "+rslt[0]+"/"+rslt[1]);
		}
		return rslt;
	}

	public double [] getMaxXSOrtho( // // get fractional center using a quadratic polynomial
			double [] data,            // [data_size * data_size]
			double [] enhortho_scales, // [data_size]
			int       data_size,
//			double    radius,  // positive - within that distance, negative - within 2*(-radius)+1 square
			boolean   debug)
	{
		if (debug) {
			System.out.println("getMaxXSOrtho() old");
			for (int i = 0; i < data_size; i++) {
				System.out.print(String.format("%2d:", i));
				for (int j = 0; j < data_size; j++) {
					System.out.print(String.format(" %8.5f", data[i * data_size + j]));
				}
				System.out.println();
			}
			System.out.println();

		}
		double [] corr_1d = new double [data_size];
		for (int j = 0; j < data_size; j++){
			corr_1d[j] = 0;
			for (int i = 0; i < data_size; i++){
				corr_1d[j] += data[i * data_size + j] * enhortho_scales[i];
			}
		}
		int icenter = 0;
		for (int i = 1; i < data_size; i++){
			if (corr_1d[i] > corr_1d[icenter]) icenter = i;
		}

		if (debug) {
			System.out.println();
			System.out.print ("corr_1d = ");
			for (int j = 0; j < data_size; j++){
				if (debug) System.out.print(String.format(" %8.5f", corr_1d[j]));
			}
			System.out.println();
		}

		double [] coeff = null;
		double xcenter = icenter;
		double [][] pa_data=null;
		// try 3-point parabola
		if ((icenter >0) && (icenter < (data_size - 1))) {
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
		double strength = corr_1d[icenter] / ((data_size+1) / 2);// scale to ~match regular strength
		double [] rslt1 = {xcenter, strength};
		return rslt1;
	}

	// value as weight (boolean), make 0 if <0,
	// select # of points (normally 5 if v_as_w, 3 if not
	// calculate x, f(x) and half width
	// balance strength? Or just assume appropriate window
	// maybe optimized to symmetrical data

	public double [] getMaxXSOrtho2(   // get fractional center using a quadratic polynomial
			double [] data,            // [data_size * data_size]
			double [] vweights,        // [data_size]
			int       data_size,
			int       num_samples,     // number of samples to keep (5?)
			double    vasw_pwr,        // use value to this power as sample weight
			boolean   debug)
	{
		double [] corr_1d = new double [data_size];
		for (int j = 0; j < data_size; j++){
			corr_1d[j] = 0;
			for (int i = 0; i < data_size; i++){
				corr_1d[j] += data[i * data_size + j] * vweights[i];
			}
		}
		int hsize = num_samples/2;
		int i0 = hsize;
		for (int i = hsize; i < (data_size - num_samples + hsize) ; i++){
			if (corr_1d[i] > corr_1d[i0]) i0 = i;
		}
		i0 -= hsize;
		double [] coeff = null;
		double [][] pa_data=new double[num_samples][(vasw_pwr > 0.0)?3:2];
		for (int i = 0; i < num_samples; i++) {
			pa_data[i][0] = i + i0;
			pa_data[i][1] = corr_1d[i + i0];
			if (vasw_pwr > 0.0) {
				pa_data[i][2] =  Math.pow(Math.abs(pa_data[i][1]), vasw_pwr);
			}
		}
		PolynomialApproximation pa = new PolynomialApproximation(debug?5:0); // debugLevel
		coeff = pa.polynomialApproximation1d(pa_data, 2);
		double xcenter, fx, hwidth;
		if ((coeff == null) || (coeff[2] > 0.0)){
			xcenter = i0 + hsize;
			fx = corr_1d[i0+hsize];
			hwidth =  data_size;
		} else {
			xcenter = - coeff[1]/(2* coeff[2]);
			fx = coeff[0] + coeff[1]*xcenter + coeff[2]*xcenter*xcenter;
			hwidth = Math.sqrt(-fx/coeff[2]);
		}
		if (debug){
			System.out.println("getMaxXSOrtho2(): getMaxXYPoly(): xcenter="+xcenter+" fx="+fx+" hwidth="+hwidth+" i0="+i0);
			System.out.println("vweights[i]=");
			for (int i=0; i < vweights.length; i++) System.out.print(vweights[i]+" "); System.out.println();
			System.out.println("corr_1d[i]=");
			for (int i=0; i < corr_1d.length; i++) System.out.print(corr_1d[i]+" "); System.out.println();
			System.out.println("pa_data[i]=");
			for (int i=0; i < pa_data.length; i++) System.out.println(i+": "+pa_data[i][0]+" "+pa_data[i][1]+" "+pa_data[i][2]+" ");
			if (coeff != null) {
				System.out.println("coeff: a = "+coeff[2]+", b="+coeff[1]+", c="+coeff[0]);
			}
			System.out.println("\n----source data ----------");
			for (int i = 0; i < data_size; i++){
				for (int j = 0; j < data_size; j++){
					System.out.print(String.format(" %6.3f", data[i * data_size + j]));
				}
				System.out.println();
			}
			System.out.println("\n---- masked ----------");
			for (int i = 0; i < data_size; i++){
				for (int j = 0; j < data_size; j++){
					System.out.print(String.format(" %6.3f", data[i * data_size + j] * vweights[i]));
				}
				System.out.println();
			}
			for (int j = 0; j < data_size; j++){
				System.out.print(String.format(" %6.3f", corr_1d[j]));
			}
			System.out.println();
		}

		double [] rslt = {xcenter, fx, hwidth};
		return rslt;

	}



	public double [] getMaxXYPoly( // get interpolated maximum coordinates using 2-nd degree polynomial
			PolynomialApproximation pa,
			double [] data,      // [data_size * data_size]
			int       data_size,
			int []    icenter, // integer center coordinates (relative to top left)
			double [] weights,   // [(radius+1) * (radius+1)]
			int       radius,
			double    value_pwr, // raise value to this power (trying to compensate sticking to integer values)
			double    poly_vasw_pwr, // multiply weight by value to this power
			boolean   debug)
	{
		// TODO: make sure it is within 1pxx1px square from the integer maximum? If not - return null and use center of mass instead?
		if (pa == null) pa = new PolynomialApproximation();
		if (icenter == null) return null; //gigo

		double [][] zdata = {{0.0,0.0},{0.0},{0.0}};
//		radius = 1;
		double [][][] mdata = new double[(2 * radius + 1) * (2 * radius + 1)][3][];
		int indx = 0;
		for (int y = - radius ; y <= radius; y++){
			int dataY = icenter[1] +y;
			if ((dataY >= 0) && (dataY < data_size)){
				int ay = (y >= 0)?y:-y;
				for (int x = - radius ; x <= radius; x++){
					int dataX = icenter[0] +x;
					if ((dataX >= 0) && (dataX < data_size)){
						int ax = (x >= 0) ? x: -x;
						mdata[indx][0] = new double [2];
						mdata[indx][0][0] =  dataX;
						mdata[indx][0][1] =  dataY;
						mdata[indx][1] = new double [1];
						mdata[indx][1][0] =  data[dataY * data_size + dataX];
						if (value_pwr != 1.0) {
							if (mdata[indx][1][0] > 0) {
								mdata[indx][1][0] = Math.pow(mdata[indx][1][0], value_pwr);
							} else {
								mdata[indx][1][0] = 0.0;
							}
						}
						mdata[indx][2] = new double [1];
						mdata[indx][2][0] =  weights[ay * (radius + 1) + ax];
						if (poly_vasw_pwr > 0 ) mdata[indx][2][0] *= Math.pow(Math.abs(mdata[indx][1][0]), poly_vasw_pwr);
						indx++;
					}
				}
			}
		}
		for (;indx <  mdata.length; indx++){
			mdata[indx] = zdata;
		}
		if (debug){
			System.out.println("before: getMaxXYPoly(): icenter[0] = "+icenter[0]+" icenter[1] = "+icenter[1]);

			for (int i = 0; i< mdata.length; i++){
				System.out.println(i+": "+mdata[i][0][0]+"/"+mdata[i][0][1]+" z="+mdata[i][1][0]+" w="+mdata[i][2][0]);
			}
		}
//		double [] rslt = pa.quadraticMax2d(
		double [] rslt = pa.quadraticMaxV2dX2Y2XY( // 6 elements - Xc, Yx, f(x,y), A, B, C (from A*x^2 + B*y^2 +C*x*y+...)
				mdata,
				1.0E-30,//25, // 1.0E-15,
				debug? 4:0);
		if (rslt == null) return null;
		// calculate width_x and width_y
		double hwx = Double.NaN, hwy = Double.NaN;
		if ((rslt[2] > 0.0) && (rslt[3] <0.0) && (rslt[4] <0.0)) {
			hwx = Math.sqrt(-rslt[2]/rslt[3]);
			hwy = Math.sqrt(-rslt[2]/rslt[4]);
		}
		double [] xyvwh = {rslt[0], rslt[1], rslt[2], hwx, hwy};
		if (debug){
			System.out.println("after: getMaxXYPoly(): icenter[0] = "+icenter[0]+" icenter[1] = "+icenter[1]);
			for (int i = 0; i< mdata.length; i++){
				System.out.println(i+": "+mdata[i][0][0]+"/"+mdata[i][0][1]+" z="+mdata[i][1][0]+" w="+mdata[i][2][0]);
			}
			System.out.println("quadraticMax2d(mdata) --> "+((rslt==null)?"null":(rslt[0]+"/"+rslt[1])));
		}
		return xyvwh; // rslt;
	}


// perform 2d clt, result is [tileY][tileX][cc_sc_cs_ss][index_in_tile]
	public double [][][][] clt_2d(
			final double [] dpixels,
			final int       width,
			final int       dct_size,
			final int       window_type,
			final int       shiftX, // shift image horizontally (positive - right)
			final int       shiftY, // shift image vertically (positive - down)
			final int       debug_tileX,
			final int       debug_tileY,
			final int       debug_mode,
			final int       threadsMax,  // maximal number of threads to launch
			final int       globalDebugLevel)
	{
		final int height=dpixels.length/width;
		final int tilesX=width/dct_size-1;
		final int tilesY=height/dct_size-1;
		final int nTiles=tilesX*tilesY;
		final double [][][][] dct_data = new double[tilesY][tilesX][4][dct_size*dct_size];
		final Thread[] threads = newThreadArray(threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		for (int tileY = 0; tileY < tilesY; tileY++){
			for (int tileX = 0; tileX < tilesX; tileX++){
				for (int dct_mode = 0; dct_mode < dct_data[tileY][tileX].length; dct_mode++){
					for (int i=0; i<dct_data[tileY][tileX][dct_mode].length;i++) {
						dct_data[tileY][tileX][dct_mode][i]= 0.0; // actually not needed, Java initializes arrays
					}
				}
			}
		}
		double [] dc = new double [dct_size*dct_size];
		for (int i = 0; i<dc.length; i++) dc[i] = 1.0;
		DttRad2 dtt0 = new DttRad2(dct_size);
		dtt0.set_window(window_type);
		if (globalDebugLevel > 0) {
			System.out.println("clt_2d(): width="+width+" height="+height+" dct_size="+dct_size+
					" debug_tileX="+debug_tileX+" debug_tileY="+debug_tileY+" globalDebugLevel="+globalDebugLevel);
		}
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				@Override
				public void run() {
					DttRad2 dtt = new DttRad2(dct_size);
					dtt.set_window(window_type);
					double [] tile_in = new double[4*dct_size * dct_size];
					double [][] tile_folded = new double[4][];
					double [][] tile_out =    new double[4][]; // = new double[dct_size * dct_size];
					int tileY,tileX;
					int n2 = dct_size * 2;
//					showDoubleFloatArrays sdfa_instance = new showDoubleFloatArrays(); // just for debugging?
					for (int nTile = ai.getAndIncrement(); nTile < nTiles; nTile = ai.getAndIncrement()) {
						tileY = nTile/tilesX;
						tileX = nTile - tileY * tilesX;
						if ((shiftX == 0) && (shiftY == 0)){
							for (int i = 0; i < n2;i++){
								System.arraycopy(dpixels, (tileY*width+tileX)*dct_size + i*width, tile_in, i*n2, n2);
							}
						} else {
							int x0 = tileX * dct_size - shiftX;
							if      (x0 < 0)             x0 = 0; // first/last will be incorrect
							else if (x0 >= (width - n2)) x0 = width - n2;
							for (int i = 0; i < n2;i++){
								int y0 = tileY * dct_size + i - shiftY;
								if      (y0 < 0)       y0 = 0;
								else if (y0 >= height) y0 = height -1;
								System.arraycopy(dpixels, y0 * width+ x0, tile_in, i*n2, n2);
							}
						}
						for (int dct_mode = 0; dct_mode <4; dct_mode++) {
							tile_folded[dct_mode] = dtt.fold_tile(tile_in, dct_size, dct_mode); // DCCT, DSCT, DCST, DSST
							if ((debug_mode & 1) != 0) {
								tile_out[dct_mode] = tile_folded[dct_mode];
							} else {
								tile_out[dct_mode] =    dtt.dttt_iv  (tile_folded[dct_mode], dct_mode, dct_size);
							}
							System.arraycopy(tile_out[dct_mode], 0, dct_data[tileY][tileX][dct_mode], 0, tile_out[dct_mode].length);
						}
						if ((globalDebugLevel > 0) && (debug_tileX == tileX) && (debug_tileY == tileY)) {
							ShowDoubleFloatArrays sdfa_instance = new ShowDoubleFloatArrays(); // just for debugging?
							sdfa_instance.showArrays(tile_in,  n2, n2, "tile_in_x"+tileX+"_y"+tileY);
							String [] titles = {"CC","SC","CS","SS"};
							sdfa_instance.showArrays(tile_folded,  dct_size, dct_size, true, "folded_x"+tileX+"_y"+tileY, titles);
							if (globalDebugLevel > 0) {
								sdfa_instance.showArrays(tile_out,     dct_size, dct_size, true, "clt_x"+tileX+"_y"+tileY, titles);
							}
						}
					}
				}
			};
		}
		startAndJoin(threads);
		return dct_data;
	}

	public double [] iclt_2d(
			final double [][][][] dct_data,  // array [tilesY][tilesX][4][dct_size*dct_size]
			final int             dct_size,
			final int             window_type,
			final int             debug_mask, // which transforms to combine
			final int             debug_mode, // skip idct - just unfold
			final int             threadsMax,  // maximal number of threads to launch
			final int             globalDebugLevel)
	{
		final int tilesY=dct_data.length;
		final int tilesX=dct_data[0].length;

//		final int width=  (tilesX+1)*dct_size;
//		final int height= (tilesY+1)*dct_size;
		final int width=  tilesX * dct_size;
		final int height= tilesY * dct_size;
		final double debug_scale = 1.0 /((debug_mask & 1) + ((debug_mask >> 1) & 1) + ((debug_mask >> 2) & 1) + ((debug_mask >> 3) & 1));
		if (globalDebugLevel > 0) {
			System.out.println("iclt_2d():tilesX=        "+tilesX);
			System.out.println("iclt_2d():tilesY=        "+tilesY);
			System.out.println("iclt_2d():width=         "+width);
			System.out.println("iclt_2d():height=        "+height);
			System.out.println("iclt_2d():debug_mask=    "+debug_mask);
			System.out.println("iclt_2d():debug_scale=   "+debug_scale);
		}
		final double [] dpixels = new double[width*height];
		final Thread[] threads = newThreadArray(threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		final AtomicInteger nser = new AtomicInteger(0);
		final int [][][] tiles_list = new int[4][][];
		for (int n=0; n<4; n++){
			int nx = (tilesX + 1 - (n &1)) / 2;
			int ny = (tilesY + 1 - ((n>>1) & 1)) / 2;
			tiles_list[n] = new int [nx*ny][2];
			int indx = 0;
			for (int i = 0;i < ny; i++) for (int j = 0; j < nx; j++){
				tiles_list[n][indx][0]=2*j+(n &1);
				tiles_list[n][indx++][1]=2*i+((n>>1) & 1);
			}
		}
		for (int i=0; i<dpixels.length;i++) dpixels[i]= 0;
		for (int n=0; n<4; n++){
			nser.set(n);
			ai.set(0);
			for (int ithread = 0; ithread < threads.length; ithread++) {
				threads[ithread] = new Thread() {
					@Override
					public void run() {
						DttRad2 dtt = new DttRad2(dct_size);
						dtt.set_window(window_type);
						double [] tile_in =   new double [dct_size * dct_size];
						double [] tile_dct;
						double [] tile_mdct;
						int tileY,tileX;
						int n2 = dct_size * 2;
						int n_half = dct_size / 2;
						int lastY = tilesY-1;
						int lastX = tilesX-1;
						int offset = n_half * (dct_size * tilesX) + n_half;
						for (int nTile = ai.getAndIncrement(); nTile < tiles_list[nser.get()].length; nTile = ai.getAndIncrement()) {
							tileX = tiles_list[nser.get()][nTile][0];
							tileY = tiles_list[nser.get()][nTile][1];
							if (dct_data[tileY][tileX] != null){
								for (int dct_mode = 0; dct_mode < 4; dct_mode++) if (((1 << dct_mode) & debug_mask) != 0) {
									System.arraycopy(dct_data[tileY][tileX][dct_mode], 0, tile_in, 0, tile_in.length);
									if ((debug_mode & 1) != 0) {
										tile_dct = tile_in;
									} else {
										// IDCT-IV should be in reversed order: CC->CC, SC->CS, CS->SC, SS->SS
										int idct_mode = ((dct_mode << 1) & 2) | ((dct_mode >> 1) & 1);
										tile_dct = dtt.dttt_iv  (tile_in, idct_mode, dct_size);
									}
									tile_mdct = dtt.unfold_tile(tile_dct, dct_size, dct_mode); // mode=0 - DCCT
									if ((tileY >0) && (tileX > 0) && (tileY < lastY) && (tileX < lastX)) { // fast, no extra checks
										for (int i = 0; i < n2;i++){
											//									int start_line = ((tileY*dct_size + i) *(tilesX+1) + tileX)*dct_size;
											int start_line = ((tileY*dct_size + i) * tilesX + tileX)*dct_size - offset;
											for (int j = 0; j<n2;j++) {
												dpixels[start_line + j] += debug_scale * tile_mdct[n2 * i + j]; // add (cc+sc+cs+ss)/4
											}
										}
									} else { // be careful with margins
										for (int i = 0; i < n2;i++){
											if (	((tileY > 0) && (tileY < lastY)) ||
													((tileY == 0) && (i >= n_half)) ||
													((tileY == lastY) && (i < (n2 - n_half)))) {
												int start_line = ((tileY*dct_size + i) * tilesX + tileX)*dct_size  - offset;
												for (int j = 0; j<n2;j++) {
													if (	((tileX > 0) && (tileX < lastX)) ||
															((tileX == 0) && (j >= n_half)) ||
															((tileX == lastX) && (j < (n2 - n_half)))) {
														dpixels[start_line + j] += debug_scale * tile_mdct[n2 * i + j]; // add (cc+sc+cs+ss)/4
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
			startAndJoin(threads);
		}
		return dpixels;
	}

	public double [] iclt_2d_debug_gpu(
			final double [][][][] dct_data,  // array [tilesY][tilesX][4][dct_size*dct_size]
			final int             dct_size,
			final int             window_type,
			final int             debug_mask, // which transforms to combine
			final int             debug_mode, // skip idct - just unfold
			final int             threadsMax,  // maximal number of threads to launch
			final int             globalDebugLevel,
			final int             debug_tileX,
			final int             debug_tileY)
	{
		final int tilesY=dct_data.length;
		final int tilesX=dct_data[0].length;
		@SuppressWarnings("unused")
		final double [][] dbg_tile = dct_data[debug_tileY][debug_tileX];

//		final int width=  (tilesX+1)*dct_size;
//		final int height= (tilesY+1)*dct_size;
		final int width=  tilesX * dct_size;
		final int height= tilesY * dct_size;
		final double debug_scale = 1.0 /((debug_mask & 1) + ((debug_mask >> 1) & 1) + ((debug_mask >> 2) & 1) + ((debug_mask >> 3) & 1));
		if (globalDebugLevel > 0) {
			System.out.println("iclt_2d():tilesX=        "+tilesX);
			System.out.println("iclt_2d():tilesY=        "+tilesY);
			System.out.println("iclt_2d():width=         "+width);
			System.out.println("iclt_2d():height=        "+height);
			System.out.println("iclt_2d():debug_mask=    "+debug_mask);
			System.out.println("iclt_2d():debug_scale=   "+debug_scale);
		}
		final double [] dpixels = new double[width*height];
		final Thread[] threads = newThreadArray(threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		final AtomicInteger nser = new AtomicInteger(0);
		final int [][][] tiles_list = new int[4][][];
		for (int n=0; n<4; n++){
			int nx = (tilesX + 1 - (n &1)) / 2;
			int ny = (tilesY + 1 - ((n>>1) & 1)) / 2;
			tiles_list[n] = new int [nx*ny][2];
			int indx = 0;
			for (int i = 0;i < ny; i++) for (int j = 0; j < nx; j++){
				tiles_list[n][indx][0]=2*j+(n &1);
				tiles_list[n][indx++][1]=2*i+((n>>1) & 1);
			}
		}
		for (int i=0; i<dpixels.length;i++) dpixels[i]= 0;
		for (int n=0; n<4; n++){
			nser.set(n);
			ai.set(0);
			for (int ithread = 0; ithread < threads.length; ithread++) {
				threads[ithread] = new Thread() {
					@Override
					public void run() {
						DttRad2 dtt = new DttRad2(dct_size);
						dtt.set_window(window_type);
						double [] tile_in =   new double [dct_size * dct_size];
						double [] tile_dct;
						double [] tile_mdct;
						int tileY,tileX;
						int n2 = dct_size * 2;
						int n_half = dct_size / 2;
						int lastY = tilesY-1;
						int lastX = tilesX-1;
						int offset = n_half * (dct_size * tilesX) + n_half;
						for (int nTile = ai.getAndIncrement(); nTile < tiles_list[nser.get()].length; nTile = ai.getAndIncrement()) {
							if (nTile == 9645) {
								System.out.println("--- nTile="+nTile+", nser.get()="+nser.get());
							}
							tileX = tiles_list[nser.get()][nTile][0];
							tileY = tiles_list[nser.get()][nTile][1];

							boolean dbg_tile =
									(tileX == debug_tileX) && (tileY == debug_tileY) &&
//									(tileX < (debug_tileX+2)) && (tileY < (debug_tileY+2)) &&
									(globalDebugLevel > -10);
							if (dbg_tile) {
								//								unfold_index = new int[4*n*n];
								//								unfold_k = new double[4][4*n*n];
								System.out.println("--- unfold_index:");
								for (int i = 0; i < n2; i++) {
									for (int j = 0; j < n2; j++) {
										System.out.print(String.format("%02x ", dtt.unfold_index[n2 * j + i]));
									}
									System.out.println();
								}
								System.out.println();
								for (int dct_mode = 0; dct_mode < 4; dct_mode++) {
									System.out.println("--- unfold_k["+dct_mode+"]:");
									for (int i = 0; i < n2; i++) {
										for (int j = 0; j < n2; j++) {
											System.out.print(String.format("%8.5f ", dtt.unfold_k[dct_mode][n2 * j + i]));
										}
										System.out.println();
									}
								}
							}


							if (dct_data[tileY][tileX] != null){
								for (int dct_mode = 0; dct_mode < 4; dct_mode++) if (((1 << dct_mode) & debug_mask) != 0) {
									System.arraycopy(dct_data[tileY][tileX][dct_mode], 0, tile_in, 0, tile_in.length);

									if ((debug_mode & 1) != 0) {
										tile_dct = tile_in;
									} else {
										// IDCT-IV should be in reversed order: CC->CC, SC->CS, CS->SC, SS->SS
										int idct_mode = ((dct_mode << 1) & 2) | ((dct_mode >> 1) & 1);

										if (dbg_tile) {
											System.out.println("--- BEFORE idct: dct_mode="+dct_mode+", idct_mode="+idct_mode+" tile_in (trasposed):");
											for (int i = 0; i < dct_size; i++) {
												for (int j = 0; j < dct_size; j++) {
													System.out.print(String.format("%11.5f ", tile_in[dct_size * j + i]));
												}
												System.out.println();
											}
											System.out.println();
										}
										tile_dct = dtt.dttt_iv  (tile_in, idct_mode, dct_size);
										if (dbg_tile) {
											System.out.println("--- AFTER idct: dct_mode="+dct_mode+", idct_mode="+idct_mode+" tile_dct (not transposed):");
											for (int i = 0; i < dct_size; i++) {
												for (int j = 0; j < dct_size; j++) {
													System.out.print(String.format("%11.5f ", tile_dct[dct_size * i + j]));
												}
												System.out.println();
											}
											System.out.println();
										}
									}
									tile_mdct = dtt.unfold_tile(tile_dct, dct_size, dct_mode); // mode=0 - DCCT

									if (dbg_tile) {
										System.out.println("--- MDCT: ");
										for (int i = 0; i < dct_size * 2; i++) {
											for (int j = 0; j < dct_size * 2; j++) {
												System.out.print(String.format("%11.5f ", tile_mdct[dct_size * 2 * i + j]));
											}
											System.out.println();
										}
										System.out.println();
									}


									if ((tileY >0) && (tileX > 0) && (tileY < lastY) && (tileX < lastX)) { // fast, no extra checks
										for (int i = 0; i < n2;i++){
											int start_line = ((tileY*dct_size + i) * tilesX + tileX)*dct_size - offset;
											for (int j = 0; j<n2;j++) {
												dpixels[start_line + j] += debug_scale * tile_mdct[n2 * i + j]; // add (cc+sc+cs+ss)/4
											}
										}
									} else { // be careful with margins
										for (int i = 0; i < n2;i++){
											if (	((tileY > 0) && (tileY < lastY)) ||
													((tileY == 0) && (i >= n_half)) ||
													((tileY == lastY) && (i < (n2 - n_half)))) {
												int start_line = ((tileY*dct_size + i) * tilesX + tileX)*dct_size  - offset;
												for (int j = 0; j<n2;j++) {
													if (	((tileX > 0) && (tileX < lastX)) ||
															((tileX == 0) && (j >= n_half)) ||
															((tileX == lastX) && (j < (n2 - n_half)))) {
														dpixels[start_line + j] += debug_scale * tile_mdct[n2 * i + j]; // add (cc+sc+cs+ss)/4
													}
												}
											}
										}
									}
								}
								if (dbg_tile) {
									System.out.println("--- MDCT, channels combined: ");
									for (int i = 0; i < n2;i++){
										int start_line = ((tileY*dct_size + i) * tilesX + tileX)*dct_size - offset;
										for (int j = 0; j<n2;j++) {
											System.out.print(String.format("%11.5f ", dpixels[start_line + j]));
										}
										System.out.println();
									}
									System.out.println();
								}
							}
						}
					}
				};
			}
			startAndJoin(threads);
		}
		return dpixels;
	}



// in monochrome mode only MONO_CHN == GREEN_CHN is used, R and B are null
	public double [][] combineRBGATiles(
			final double [][][][] texture_tiles,  // array [tilesY][tilesX][4][4*transform_size] or [tilesY][tilesX]{null}
			final int             transform_size,
			final boolean         overlap,    // when false - output each tile as 16x16, true - overlap to make 8x8
			final boolean         sharp_alpha, // combining mode for alpha channel: false - treat as RGB, true - apply center 8x8 only
			final int             threadsMax,  // maximal number of threads to launch
			final int             globalDebugLevel)
	{
		final int tilesY=texture_tiles.length;
		final int tilesX=texture_tiles[0].length;

		final int width=  (overlap?1:2)*tilesX * transform_size;
		final int height=  (overlap?1:2)*tilesY * transform_size;
		if (globalDebugLevel > 0) {
			System.out.println("iclt_2d():tilesX=        "+tilesX);
			System.out.println("iclt_2d():tilesY=        "+tilesY);
			System.out.println("iclt_2d():width=         "+width);
			System.out.println("iclt_2d():height=        "+height);
			System.out.println("iclt_2d():overlap=       "+overlap);
			System.out.println("iclt_2d():sharp_alpha=   "+sharp_alpha);
		}
		boolean has_weights = false;
		boolean set_has_weight = false;
		for (int i = 0; (i < tilesY) && !set_has_weight; i++){
			for (int j = 0; (j < tilesX) && !set_has_weight; j++){
				if (texture_tiles[i][j] != null) {
					set_has_weight = true;
					has_weights = texture_tiles[i][j].length > 4;
				}
			}
		}

//		final double [][] dpixels = new double["RGBA".length()+(has_weights? 4: 0)][width*height]; // assuming java initializes them to 0
		final double [][] dpixels = new double["RGBA".length()+(has_weights? 8: 0)][width*height]; // assuming java initializes them to 0
		final Thread[] threads = newThreadArray(threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		final AtomicInteger nser = new AtomicInteger(0);
		final int [][][] tiles_list = new int[4][][];
		for (int n=0; n<4; n++){
			int nx = (tilesX + 1 - (n &1)) / 2;
			int ny = (tilesY + 1 - ((n>>1) & 1)) / 2;
			tiles_list[n] = new int [nx*ny][2];
			int indx = 0;
			for (int i = 0;i < ny; i++) for (int j = 0; j < nx; j++){
				tiles_list[n][indx][0]=2*j+(n &1);
				tiles_list[n][indx++][1]=2*i+((n>>1) & 1);
			}
		}
		for (int n=0; n<4; n++){
			nser.set(n);
			ai.set(0);
			for (int ithread = 0; ithread < threads.length; ithread++) {
				threads[ithread] = new Thread() {
					@Override
					public void run() {
						int tileY,tileX;
						int n2 = transform_size * 2;
						int n_half = transform_size / 2;
						int lastY = tilesY-1;
						int lastX = tilesX-1;
						int offset = n_half * (transform_size * tilesX) + n_half;
						for (int nTile = ai.getAndIncrement(); nTile < tiles_list[nser.get()].length; nTile = ai.getAndIncrement()) {
							tileX = tiles_list[nser.get()][nTile][0];
							tileY = tiles_list[nser.get()][nTile][1];
							double [][] texture_tile =texture_tiles[tileY][tileX];
							if (texture_tile != null) {
								if (overlap) {
									if ((tileY >0) && (tileX > 0) && (tileY < lastY) && (tileX < lastX)) { // fast, no extra checks
										for (int i = 0; i < n2;i++){
											int start_line = ((tileY*transform_size + i) * tilesX + tileX)*transform_size - offset;
											for (int chn = 0; chn < texture_tile.length; chn++) {
												int schn = chn;
												if (isMonochrome() && (chn<3)) {
													schn = MONO_CHN; // clone green to red and blue output
												}
												if (texture_tile[schn] == null) {
													dpixels[chn] = null;
												} else {
													if ((chn != 3) || !sharp_alpha) {
														for (int j = 0; j<n2;j++) {
															dpixels[chn][start_line + j] += texture_tile[schn][n2 * i + j];
														}
													} else if ((i >= n_half) && (i < (n2-n_half))) {
														for (int j = n_half; j < (n2 - n_half); j++) {
															dpixels[chn][start_line + j] += texture_tile[schn][n2 * i + j];
														}
													}
												}
											}
										}
									} else { // be careful with margins
										for (int i = 0; i < n2;i++){
											if (	((tileY > 0) && (tileY < lastY)) ||
													((tileY == 0) && (i >= n_half)) ||
													((tileY == lastY) && (i < (n2 - n_half)))) {
												int start_line = ((tileY*transform_size + i) * tilesX + tileX)*transform_size  - offset;
												for (int chn = 0; chn < texture_tile.length; chn++) {
													int schn = chn;
													if (isMonochrome() && (chn<3)) {
														schn = MONO_CHN; // clone green to red and blue output
													}
													if (texture_tile[schn] == null) {
														dpixels[chn] = null;
													} else {
														if ((chn != 3) || !sharp_alpha) {
															for (int j = 0; j<n2;j++) {
																if (	((tileX > 0) && (tileX < lastX)) ||
																		((tileX == 0) && (j >= n_half)) ||
																		((tileX == lastX) && (j < (n2 - n_half)))) {
																	dpixels[chn][start_line + j] += texture_tile[schn][n2 * i + j];
																}
															}
														} else if ((i >= n_half) && (i < (n2-n_half))) {
															for (int j = n_half; j < (n2 - n_half); j++) {
																if (	((tileX > 0) && (tileX < lastX)) ||
																		((tileX == 0) && (j >= n_half)) ||
																		((tileX == lastX) && (j < (n2 - n_half)))) {
																	dpixels[chn][start_line + j] += texture_tile[schn][n2 * i + j];
																}
															}
														}
													}
												}
											}
										}

									}
								} else { //if (overlap) - just copy tiles w/o overlapping
									for (int i = 0; i < n2;i++){
										for (int chn = 0; chn < texture_tile.length; chn++) {
											int schn = chn;
											if (isMonochrome() && (chn<3)) {
												schn = MONO_CHN; // clone green to red and blue output
											}
											if (texture_tile[schn] == null) {
												dpixels[chn] = null;
											} else {

												System.arraycopy(
														texture_tile[schn],
														i * n2,
														dpixels[chn],
														(tileY * n2 + i)* width + tileX*n2,
														n2);
											}
										}
									}
								}
							}
						}
					}
				};
			}
			startAndJoin(threads);
		}
		return dpixels;
	}




	public double [][][][] clt_shiftXY(
			final double [][][][] dct_data,  // array [tilesY][tilesX][4][dct_size*dct_size]
			final int             dct_size,
			final double          shiftX,
			final double          shiftY,
			final int             dbg_swap_mode,
			final int             threadsMax,  // maximal number of threads to launch
			final int             globalDebugLevel)
	{
		final int tilesY=dct_data.length;
		final int tilesX=dct_data[0].length;
		final int nTiles = tilesY* tilesX;
		if (globalDebugLevel > 0) {
			System.out.println("clt_shift():tilesX=        "+tilesX);
			System.out.println("clt_shift():tilesY=        "+tilesY);
			System.out.println("clt_shift():shiftX=        "+shiftX);
			System.out.println("clt_shift():shiftY=        "+shiftY);
		}
		final double [] cos_hor =  new double [dct_size*dct_size];
		final double [] sin_hor =  new double [dct_size*dct_size];
		final double [] cos_vert = new double [dct_size*dct_size];
		final double [] sin_vert = new double [dct_size*dct_size];
		for (int i = 0; i < dct_size; i++){
			double ch = Math.cos((i+0.5)*Math.PI*shiftX/dct_size);
			double sh = Math.sin((i+0.5)*Math.PI*shiftX/dct_size);
			double cv = Math.cos((i+0.5)*Math.PI*shiftY/dct_size);
			double sv = Math.sin((i+0.5)*Math.PI*shiftY/dct_size);
			for (int j = 0; j < dct_size; j++){
				int iv = dct_size * j + i; // 2d DTT results are stored transposed!
				int ih = dct_size * i + j;
				cos_hor[ih] = ch;
				sin_hor[ih] = sh;
				cos_vert[iv] = cv;
				sin_vert[iv] = sv;
			}

		}

		if (globalDebugLevel > 1){
			ShowDoubleFloatArrays sdfa_instance = new ShowDoubleFloatArrays(); // just for debugging?
			String [] titles = {"cos_hor","sin_hor","cos_vert","sin_vert"};
			double [][] cs_dbg = {cos_hor, sin_hor, cos_vert, sin_vert};
			sdfa_instance.showArrays(cs_dbg,  dct_size, dct_size, true, "shift_cos_sin", titles);
		}

		final double [][][][] rslt = new double[dct_data.length][dct_data[0].length][dct_data[0][0].length][dct_data[0][0][0].length];
		final Thread[] threads = newThreadArray(threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				@Override
				public void run() {
					int tileY,tileX;
					for (int nTile = ai.getAndIncrement(); nTile < nTiles; nTile = ai.getAndIncrement()) {
						tileY = nTile/tilesX;
						tileX = nTile - tileY * tilesX;
						// Horizontal shift CLT tiled data is stored in transposed way (horizontal - Y, vertical X)
						for (int i = 0; i < cos_hor.length; i++) {
							rslt[tileY][tileX][0][i] = dct_data[tileY][tileX][0][i] * cos_hor[i] - dct_data[tileY][tileX][1][i] * sin_hor[i];
							rslt[tileY][tileX][1][i] = dct_data[tileY][tileX][1][i] * cos_hor[i] + dct_data[tileY][tileX][0][i] * sin_hor[i] ;

							rslt[tileY][tileX][2][i] = dct_data[tileY][tileX][2][i] * cos_hor[i]  - dct_data[tileY][tileX][3][i] * sin_hor[i];
							rslt[tileY][tileX][3][i] = dct_data[tileY][tileX][3][i] * cos_hor[i]  + dct_data[tileY][tileX][2][i] * sin_hor[i] ;
						}
						// Vertical shift (in-place)
						for (int i = 0; i < cos_hor.length; i++) {
							double tmp =               rslt[tileY][tileX][0][i] * cos_vert[i] - rslt[tileY][tileX][2][i] * sin_vert[i];
							rslt[tileY][tileX][2][i] = rslt[tileY][tileX][2][i] * cos_vert[i] + rslt[tileY][tileX][0][i] * sin_vert[i];
							rslt[tileY][tileX][0][i] = tmp;

							tmp =                      rslt[tileY][tileX][1][i] * cos_vert[i] - rslt[tileY][tileX][3][i] * sin_vert[i];
							rslt[tileY][tileX][3][i] = rslt[tileY][tileX][3][i] * cos_vert[i] + rslt[tileY][tileX][1][i] * sin_vert[i];
							rslt[tileY][tileX][1][i] = tmp;
						}
					}
				}
			};
		}
		startAndJoin(threads);
		return rslt;
	}

	public double [][][][] clt_correlate(
			final double [][][][] data1,  // array [tilesY][tilesX][4][dct_size*dct_size]
			final double [][][][] data2,  // array [tilesY][tilesX][4][dct_size*dct_size]
			final int             dct_size,
			final double          fat_zero,    // add to denominator to modify phase correlation (same units as data1, data2)
			final int             debug_tileX,
			final int             debug_tileY,
			final int             threadsMax,  // maximal number of threads to launch
			final int             globalDebugLevel)
	{
		final int tilesY=(data1.length > data2.length)?data2.length:data1.length;
		final int tilesX=(data1[0].length > data2[0].length)?data2[0].length:data1[0].length;
		final int nTiles = tilesY* tilesX;
		if (globalDebugLevel > 0) {
			System.out.println("clt_shift():tilesX= "+tilesX);
			System.out.println("clt_shift():tilesY= "+tilesY);
		}
		/* Direct matrix Z1: X2 ~= Z1 * Shift
		 * {{+cc  -sc  -cs  +ss},
		 *  {+sc  +cc  -ss  -cs},
		 *  {+cs  -ss  +cc  -sc},
		 *  {+ss  +cs  +sc  +cc}}
		 *
		 * T= transp({cc, sc, cs, ss})
		 */
		/*
		final int [][] zi =
			{{ 0, -1, -2,  3},
			 { 1,  0, -3, -2},
			 { 2, -3,  0, -1},
			 { 3,  2,  1,  0}};
		*/
		final int [][] zi =
			{{ 0,  1,  2,  3},
			 {-1,  0, -3,  2},
			 {-2, -3,  0,  1},
			 { 3, -2, -1,  0}};

		final int dct_len = dct_size * dct_size;
		final double [][][][] rslt = new double[tilesY][tilesX][4][dct_len];
		final Thread[] threads = newThreadArray(threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				@Override
				public void run() {
					int tileY,tileX;
					for (int nTile = ai.getAndIncrement(); nTile < nTiles; nTile = ai.getAndIncrement()) {
						tileY = nTile/tilesX;
						tileX = nTile - tileY * tilesX;
						for (int i = 0; i < dct_len; i++) {
							double s1 = 0.0, s2=0.0;
							for (int n = 0; n< 4; n++){
								s1+=data1[tileY][tileX][n][i] * data1[tileY][tileX][n][i];
								s2+=data2[tileY][tileX][n][i] * data2[tileY][tileX][n][i];
							}
							double scale = 1.0 / (Math.sqrt(s1*s2) + fat_zero*fat_zero); // squared to match units
							for (int n = 0; n<4; n++){
								/*
								if (
										(tileY >= rslt.length) ||
										(tileX >= rslt[tileY].length) ||
										(n >= rslt[tileY][tileX].length) ||
										(i >= rslt[tileY][tileX][n].length)) {

									System.out.println("===== tileY="+tileY+" ("+tilesY+") tileX="+tileX+" ("+tilesX+") n="+n+" i="+i);

									System.out.println(
											" rslt.length="+rslt.length+
											" rslt.length[tileY]="+rslt[tileY].length+
											" rslt.length[tileY][tileX]="+rslt[tileY][tileX].length+
											" rslt.length[tileY][tileX][n]="+rslt[tileY][tileX][n].length);
								System.out.println("===== tileY="+tileY+" ("+tilesY+") tileX="+tileX+" ("+tilesX+") n="+n+" i="+i);
								}
								*/
								rslt[tileY][tileX][n][i] = 0;
								for (int k=0; k<4; k++){
									if (zi[n][k] < 0)
										rslt[tileY][tileX][n][i] -=
											data1[tileY][tileX][-zi[n][k]][i] * data2[tileY][tileX][k][i];
									else
										rslt[tileY][tileX][n][i] +=
										data1[tileY][tileX][zi[n][k]][i] * data2[tileY][tileX][k][i];
								}
								rslt[tileY][tileX][n][i] *= scale;
							}
						}
						if ((globalDebugLevel > 0) && (debug_tileX == tileX) && (debug_tileY == tileY)) {
							ShowDoubleFloatArrays sdfa_instance = new ShowDoubleFloatArrays(); // just for debugging?
							String [] titles = {"CC","SC","CS","SS"};
							sdfa_instance.showArrays(data1[tileY][tileX], dct_size, dct_size, true, "data1_x"+tileX+"_y"+tileY, titles);
							sdfa_instance.showArrays(data2[tileY][tileX], dct_size, dct_size, true, "data2_x"+tileX+"_y"+tileY, titles);
							sdfa_instance.showArrays(rslt[tileY][tileX],  dct_size, dct_size, true, "rslt_x"+ tileX+"_y"+tileY, titles);
						}
					}
				}
			};
		}
		startAndJoin(threads);
		return rslt;
	}

	public void clt_lpf(
			final double          sigma,
			final double [][][][] clt_data,
			final int             dct_size,
			final int             threadsMax,     // maximal number of threads to launch
			final int             globalDebugLevel)
	{
		if (clt_data == null) {
			System.out.println("clt_lpf(): clt_data=null");
			return;
		}
		final int tilesY=clt_data.length;
		final int tilesX=clt_data[0].length;
		final int nTiles=tilesX*tilesY;
//		final int dct_size = (int) Math.round(Math.sqrt(clt_data[0][0][0].length));
		final int dct_len = dct_size*dct_size;
		final double [] filter_direct= new double[dct_len];
		if (sigma == 0) {
			filter_direct[0] = 1.0;
			for (int i= 1; i<filter_direct.length;i++) filter_direct[i] =0;
		} else {
			for (int i = 0; i < dct_size; i++){
				for (int j = 0; j < dct_size; j++){
					filter_direct[i*dct_size+j] = Math.exp(-(i*i+j*j)/(2*sigma));
				}
			}
		}
		// normalize
		double sum = 0;
		for (int i = 0; i < dct_size; i++){
			for (int j = 0; j < dct_size; j++){
				double d = 	filter_direct[i*dct_size+j];
				d*=Math.cos(Math.PI*i/(2*dct_size))*Math.cos(Math.PI*j/(2*dct_size));
				if (i > 0) d*= 2.0;
				if (j > 0) d*= 2.0;
				sum +=d;
			}
		}
		for (int i = 0; i<filter_direct.length; i++){
			filter_direct[i] /= sum;
		}

		if (globalDebugLevel > 1) {
			for (int i=0; i<filter_direct.length;i++){
				System.out.println("dct_lpf_psf() "+i+": "+filter_direct[i]);
			}
		}
		DttRad2 dtt = new DttRad2(dct_size);
		final double [] filter= dtt.dttt_iiie(filter_direct);
		final double [] dbg_filter= dtt.dttt_ii(filter);
		for (int i=0; i < filter.length;i++) filter[i] *= 2*dct_size;

		if (globalDebugLevel > 2) {
			System.out.print("__constant__ float lpf_data[64]={");
			for (int i=0; i<filter.length;i++){
				System.out.print(String.format("%5.8ff", filter[i]));
				if (i == 63) {
					System.out.println("};");
				} else {
					System.out.print(", ");
					if ((i % 8) == 7) {
						System.out.print("\n                                 ");
					}
				}
			}
		} else	if (globalDebugLevel > 1) {
			for (int i=0; i<filter.length;i++){
				System.out.println("dct_lpf_psf() "+i+": "+filter[i]);
			}
			ShowDoubleFloatArrays sdfa_instance = new ShowDoubleFloatArrays(); // just for debugging?
			double [][] ff = {filter_direct,filter,dbg_filter};
			sdfa_instance.showArrays(ff,  dct_size,dct_size, true, "filter_lpf");
		}

		final Thread[] threads = newThreadArray(threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);

		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				@Override
				public void run() {
					int tileY,tileX;
					for (int nTile = ai.getAndIncrement(); nTile < nTiles; nTile = ai.getAndIncrement()) {
						tileY = nTile/tilesX;
						tileX = nTile - tileY * tilesX;
						if (clt_data[tileY][tileX] != null) {
							for (int n = 0; n < 4; n++){
								for (int i = 0; i < filter.length; i++){
									clt_data[tileY][tileX][n][i] *= filter[i];
								}
							}
						}
					}
				}
			};
		}
		startAndJoin(threads);
	}

	public void clt_dtt2( // transform dcct2, dsct2, dcst2, dsst2
			final double [][][][] data,
			final boolean         transpose, // when doing inverse transform, the data comes in transposed form, so CS <->SC
			final int             threadsMax,     // maximal number of threads to launch
			final int             globalDebugLevel)
	{
		final int tilesY=data.length;
		final int tilesX=data[0].length;
		final int nTiles=tilesX*tilesY;
		final int dct_size = (int) Math.round(Math.sqrt(data[0][0][0].length));
		final Thread[] threads = newThreadArray(threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);

		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				@Override
				public void run() {
					DttRad2 dtt = new DttRad2(dct_size);
					int tileY,tileX;
					for (int nTile = ai.getAndIncrement(); nTile < nTiles; nTile = ai.getAndIncrement()) {
						tileY = nTile/tilesX;
						tileX = nTile - tileY * tilesX;
						for (int quadrant = 0; quadrant < 4; quadrant++){
							int mode = transpose ? (((quadrant << 1) & 2) | ((quadrant >> 1) & 1)) : quadrant;
							data[tileY][tileX][quadrant] = dtt.dttt_iie(data[tileY][tileX][quadrant], mode, dct_size);
						}
					}
				}
			};
		}
		startAndJoin(threads);
	}

	public double [][][] clt_corr_quad( // combine 4 correlation quadrants after DTT2
			final double [][][][] data,
			final int             threadsMax,     // maximal number of threads to launch
			final int             globalDebugLevel)
	{
		final int tilesY=data.length;
		final int tilesX=data[0].length;
		final int nTiles=tilesX*tilesY;
		final int dct_size = (int) Math.round(Math.sqrt(data[0][0][0].length));
		final int rslt_size=dct_size*2-1;

		final double [][][] rslt = new double[tilesY][tilesX][rslt_size*rslt_size];

		final Thread[] threads = newThreadArray(threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);

		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				@Override
				public void run() {
					int tileY,tileX;
					double scale = 0.25;
					for (int nTile = ai.getAndIncrement(); nTile < nTiles; nTile = ai.getAndIncrement()) {
						tileY = nTile/tilesX;
						tileX = nTile - tileY * tilesX;
						rslt[tileY][tileX][rslt_size*dct_size - dct_size] = scale * data[tileY][tileX][0][0]; // center
						for (int j = 1; j < dct_size; j++) { //  for i == 0
							rslt[tileY][tileX][rslt_size*dct_size - dct_size + j] = scale * (data[tileY][tileX][0][j] + data[tileY][tileX][1][j-1]);
							rslt[tileY][tileX][rslt_size*dct_size - dct_size - j] = scale * (data[tileY][tileX][0][j] - data[tileY][tileX][1][j-1]);
						}
						for (int i = 1; i < dct_size; i++) {
							rslt[tileY][tileX][rslt_size*(dct_size + i) - dct_size] =
									scale * (data[tileY][tileX][0][i*dct_size] + data[tileY][tileX][2][(i-1)*dct_size]);
							rslt[tileY][tileX][rslt_size*(dct_size - i) - dct_size] =
									scale * (data[tileY][tileX][0][i*dct_size] - data[tileY][tileX][2][(i-1)*dct_size]);
							for (int j = 1; j < dct_size; j++) {
								rslt[tileY][tileX][rslt_size*(dct_size + i) - dct_size + j] =
										scale * (data[tileY][tileX][0][i*    dct_size + j] +
												 data[tileY][tileX][1][i*    dct_size + j - 1] +
												 data[tileY][tileX][2][(i-1)*dct_size + j] +
												 data[tileY][tileX][3][(i-1)*dct_size + j - 1]);

								rslt[tileY][tileX][rslt_size*(dct_size + i) - dct_size - j] =
										scale * ( data[tileY][tileX][0][i*    dct_size + j] +
												 -data[tileY][tileX][1][i*    dct_size + j - 1] +
												  data[tileY][tileX][2][(i-1)*dct_size + j] +
												 -data[tileY][tileX][3][(i-1)*dct_size + j - 1]);
								rslt[tileY][tileX][rslt_size*(dct_size - i) - dct_size + j] =
										scale * (data[tileY][tileX][0][i*    dct_size + j] +
												 data[tileY][tileX][1][i*    dct_size + j - 1] +
												 -data[tileY][tileX][2][(i-1)*dct_size + j] +
												 -data[tileY][tileX][3][(i-1)*dct_size + j - 1]);
								rslt[tileY][tileX][rslt_size*(dct_size - i) - dct_size - j] =
										scale * (data[tileY][tileX][0][i*    dct_size + j] +
												 -data[tileY][tileX][1][i*    dct_size + j - 1] +
												 -data[tileY][tileX][2][(i-1)*dct_size + j] +
												 data[tileY][tileX][3][(i-1)*dct_size + j - 1]);
							}
						}
					}
				}
			};
		}
		startAndJoin(threads);
		return rslt;
	}

/*
	private double [] corr_unfold_tile(
			double [][]  qdata, // [4][transform_size*transform_size] data after DCT2 (pixel domain)
			int          transform_size
			)
	{
		int corr_pixsize = transform_size * 2 - 1;
		double corr_pixscale = 0.25;
		double [] rslt = new double [corr_pixsize*corr_pixsize];
		rslt[corr_pixsize*transform_size - transform_size] = corr_pixscale * qdata[0][0]; // center
		for (int j = 1; j < transform_size; j++) { //  for i == 0
			rslt[corr_pixsize*transform_size - transform_size + j] = corr_pixscale * (qdata[0][j] + qdata[1][j-1]);
			rslt[corr_pixsize*transform_size - transform_size - j] = corr_pixscale * (qdata[0][j] - qdata[1][j-1]);
		}
		for (int i = 1; i < transform_size; i++) {
			rslt[corr_pixsize*(transform_size + i) - transform_size] =
					corr_pixscale * (qdata[0][i*transform_size] + qdata[2][(i-1)*transform_size]);
			rslt[corr_pixsize*(transform_size - i) - transform_size] =
					corr_pixscale * (qdata[0][i*transform_size] - qdata[2][(i-1)*transform_size]);
			for (int j = 1; j < transform_size; j++) {
				rslt[corr_pixsize*(transform_size + i) - transform_size + j] =
						corr_pixscale * (qdata[0][i*    transform_size + j] +
								 qdata[1][i*    transform_size + j - 1] +
								 qdata[2][(i-1)*transform_size + j] +
								 qdata[3][(i-1)*transform_size + j - 1]);

				rslt[corr_pixsize*(transform_size + i) - transform_size - j] =
						corr_pixscale * ( qdata[0][i*    transform_size + j] +
								 -qdata[1][i*    transform_size + j - 1] +
								  qdata[2][(i-1)*transform_size + j] +
								 -qdata[3][(i-1)*transform_size + j - 1]);
				rslt[corr_pixsize*(transform_size - i) - transform_size + j] =
						corr_pixscale * (qdata[0][i*    transform_size + j] +
								 qdata[1][i*    transform_size + j - 1] +
								 -qdata[2][(i-1)*transform_size + j] +
								 -qdata[3][(i-1)*transform_size + j - 1]);
				rslt[corr_pixsize*(transform_size - i) - transform_size - j] =
						corr_pixscale * (qdata[0][i*    transform_size + j] +
								 -qdata[1][i*    transform_size + j - 1] +
								 -qdata[2][(i-1)*transform_size + j] +
								 qdata[3][(i-1)*transform_size + j - 1]);
			}
		}


		return rslt;

	}

*/
	// extract correlation result  in linescan order (for visualization)
	public double [] corr_dbg(
			final double [][][] corr_data,
			final int           corr_size,
			final double        border_contrast,
			final int           threadsMax,     // maximal number of threads to launch
			final int           globalDebugLevel)
	{
		final int tilesY=corr_data.length;
		final int tilesX=corr_data[0].length;
		final int nTiles=tilesX*tilesY;
		final int tile_size = corr_size+1;
		final int corr_len = corr_size*corr_size;

		final double [] corr_data_out = new double[tilesY*tilesX*tile_size*tile_size];

		System.out.println("corr_dbg(): tilesY="+tilesY+", tilesX="+tilesX+", corr_size="+corr_size+", corr_len="+corr_len);

		final Thread[] threads = newThreadArray(threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		for (int i=0; i<corr_data_out.length;i++) corr_data_out[i]= 0;

		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				@Override
				public void run() {
					int tileY,tileX;
					for (int nTile = ai.getAndIncrement(); nTile < nTiles; nTile = ai.getAndIncrement()) {
						tileY = nTile/tilesX;
						tileX = nTile - tileY * tilesX;
						if (corr_data[tileY][tileX] != null) {
							for (int i = 0; i < corr_size;i++){
								System.arraycopy(corr_data[tileY][tileX], corr_size* i, corr_data_out, ((tileY*tile_size + i) *tilesX + tileX)*tile_size , corr_size);
								corr_data_out[((tileY*tile_size + i) *tilesX + tileX)*tile_size+corr_size] = border_contrast*((i & 1) - 0.5);
							}
							for (int i = 0; i < tile_size; i++){
								corr_data_out[((tileY*tile_size + corr_size) *tilesX + tileX)*tile_size+i] = border_contrast*((i & 1) - 0.5);
							}
						}
					}
				}
			};
		}
		startAndJoin(threads);
		return corr_data_out;
	}




	// extract correlation result  in linescan order (for visualization)
	public double [][] corr_partial_dbg(
			final double [][][][][] corr_data,
			final int corr_size,
			final int pairs,
			final int colors,
			final double            border_contrast,
			final int               threadsMax,     // maximal number of threads to launch
			final int               globalDebugLevel)
	{
		final int tilesY=corr_data.length;
		final int tilesX=corr_data[0].length;
		final int nTiles=tilesX*tilesY;
		final int tile_size = corr_size+1;
		final int corr_len = corr_size*corr_size;

		System.out.println("corr_partial_dbg(): tilesY="+tilesY+", tilesX="+tilesX+", corr_size="+corr_size+", corr_len="+corr_len+
				" pairs="+pairs +" colors = "+colors+" tile_size="+tile_size);

		final double [][] corr_data_out = new double[pairs*colors][tilesY*tilesX*tile_size*tile_size];
//		final String [] colorNames = {"red","blue","green","composite"};

		final Thread[] threads = newThreadArray(threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		for (int pair = 0; pair< pairs; pair++) {
			for (int nColor = 0; nColor < colors; nColor++) {
				for (int i=0; i<corr_data_out.length;i++) corr_data_out[pair*colors+nColor][i]= 0;
			}
		}

		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				@Override
				public void run() {
					int tileY,tileX;
					for (int nTile = ai.getAndIncrement(); nTile < nTiles; nTile = ai.getAndIncrement()) {
						tileY = nTile/tilesX;
						tileX = nTile - tileY * tilesX;
						if (corr_data[tileY][tileX] != null) {
							for (int pair = 0; pair< pairs; pair++) {
								for (int nColor = 0; nColor < colors; nColor++) {
									int indx = pair*colors+nColor;
									for (int i = 0; i < corr_size;i++){
										System.arraycopy(
												corr_data[tileY][tileX][pair][nColor],
												corr_size* i,
												corr_data_out[indx],
												((tileY*tile_size + i) *tilesX + tileX)*tile_size ,
												corr_size);
										corr_data_out[indx][((tileY*tile_size + i) *tilesX + tileX)*tile_size+corr_size] = border_contrast*((i & 1) - 0.5);
									}
									for (int i = 0; i < tile_size; i++){
										corr_data_out[indx][((tileY*tile_size + corr_size) *tilesX + tileX)*tile_size+i] = border_contrast*((i & 1) - 0.5);
									}
								}
							}
						}
					}
				}
			};
		}
		startAndJoin(threads);
		return corr_data_out;
	}






	public double [][][][][] cltStack(
			final ImageStack                                 imageStack,
			final int                                        subcamera, //
			final CLTParameters   cltParameters, //
			final int                                        shiftX, // shift image horizontally (positive - right)
			final int                                        shiftY, // shift image vertically (positive - down)
			final int                                        threadsMax, // maximal step in pixels on the maxRadius for 1 angular step (i.e. 0.5)
			final int                                        debugLevel,
			final boolean                                    updateStatus) // update status info

	{
	  	  if (imageStack==null) return null;
		  final int imgWidth=imageStack.getWidth();
		  final int nChn=imageStack.getSize();
		  double [][][][][] dct_data = new double [nChn][][][][];
		  float [] fpixels;
		  int i,chn; //tileX,tileY;

		  for (chn=0;chn<nChn;chn++) {
			  fpixels= (float[]) imageStack.getPixels(chn+1);
			  double[] dpixels = new double[fpixels.length];
			  for (i = 0; i <fpixels.length;i++) dpixels[i] = fpixels[i];
			  // convert each to DCT tiles
			  dct_data[chn] = clt_2d(
						dpixels,
						imgWidth,
						cltParameters.transform_size,
						cltParameters.clt_window,
						shiftX,
						shiftY,
						cltParameters.tileX,    //       debug_tileX,
						cltParameters.tileY,    //       debug_tileY,
						cltParameters.dbg_mode, //       debug_mode,
						threadsMax,  // maximal number of threads to launch
						debugLevel);
		  }
		return dct_data;
	}




	// extract DCT transformed parameters in linescan order (for visualization)
	public double [][] clt_dbg(
			final double [][][][] dct_data,
			final int             threadsMax,     // maximal number of threads to launch
			final int             globalDebugLevel)
	{
		final int tilesY=dct_data.length;
		final int tilesX=dct_data[0].length;
		final int nTiles=tilesX*tilesY;
		final int dct_size = (int) Math.round(Math.sqrt(dct_data[0][0][0].length));
		final int dct_len = dct_size*dct_size;
		final double [][] dct_data_out = new double[4][tilesY*tilesX*dct_len];

		System.out.println("clt_dbg(): tilesY="+tilesY+", tilesX="+tilesX+", dct_size="+dct_size+", dct_len="+dct_len+", dct_data_out[0].length="+dct_data_out[0].length);

		final Thread[] threads = newThreadArray(threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		for (int n=0; n<dct_data_out.length;n++) for (int i=0; i<dct_data_out[n].length;i++) dct_data_out[n][i]= 0;

		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				@Override
				public void run() {
					int tileY,tileX;
					for (int nTile = ai.getAndIncrement(); nTile < nTiles; nTile = ai.getAndIncrement()) {
						tileY = nTile/tilesX;
						tileX = nTile - tileY * tilesX;
						for (int n=0; n<dct_data_out.length;n++) {
							for (int i = 0; i < dct_size;i++){
								System.arraycopy(dct_data[tileY][tileX][n], dct_size* i, dct_data_out[n], ((tileY*dct_size + i) *tilesX + tileX)*dct_size , dct_size);
							}
						}
					}
				}
			};
		}
		startAndJoin(threads);
		return dct_data_out;
	}

	void clt_convert_double_kernel( // converts double resolution kernel
			double []   src_kernel, //
			double []   dst_kernel, // should be (2*dtt_size-1) * (2*dtt_size-1) + extra_items size - kernel and dx, dy to the nearest 1/2 pixels + actual full center shift)
			int src_size, // 64
			int dtt_size) // 8
	{

		int [] indices = {0,-src_size,-1,1,src_size,-src_size-1,-src_size+1,src_size-1,src_size+1};
//		double [] weights = {0.25,0.125,0.125,0.125,0.125,0.0625,0.0625,0.0625,0.0625};
		// sum = 4.0, so sum of all kernels is ~ the same
		double [] weights = {1.0, 0.5,  0.5,  0.5,  0.5,  0.25,  0.25,  0.25,  0.25};
		int src_center = src_size / 2; // 32
		// Find center
		double sx=0.0, sy = 0.0, s = 0.0;
		int indx = 0;
		for (int i= -src_center; i < src_center; i++){
			for (int j = -src_center; j < src_center; j++){
				double d = src_kernel[indx++];
				sx+= j*d;
				sy+= i*d;
				s += d;
			}
		}
		int src_x = (int) Math.round(sx / s) + src_center;
		int src_y = (int) Math.round(sy / s) + src_center;
		// make sure selected area (2*dst_size-1) * (2*dst_size-1) fits into src_kernel, move center if not
		if      (src_x < 2 * dtt_size)             src_x = 2 * dtt_size - 1; // 15
		else if (src_x > (src_size - 2* dtt_size)) src_x = src_size - 2* dtt_size;

		if      (src_y < 2 * dtt_size)             src_y = 2 * dtt_size - 1; // 15
		else if (src_y > (src_size - 2* dtt_size)) src_y = src_size - 2* dtt_size;
		indx = 0;
		// downscale, copy
		for (int i = -dtt_size + 1; i < dtt_size; i++){
			int src_i = (src_y + 2 * i) * src_size  + src_x;
			for (int j = -dtt_size + 1; j < dtt_size; j++){
				double d = 0.0;
				for (int k = 0; k < indices.length; k++){
					d += weights[k]*src_kernel[src_i + 2 * j + indices[k]];
				}
				dst_kernel[indx++] = d;
			}
		}
		dst_kernel[indx++] = 0.5*(src_x - src_center);
		dst_kernel[indx++] = 0.5*(src_y - src_center);
		dst_kernel[indx++] = 0.5*(sx / s);             // actual center shift in pixels (to interapolate difference to neighbour regions)
		dst_kernel[indx++] = 0.5*(sy / s);
	}

	void clt_normalize_kernel( //
			double []   kernel, // should be (2*dtt_size-1) * (2*dtt_size-1) + 4 size (last (2*dtt_size-1) are not modified)
			double []   window, // normalizes result kernel * window to have sum of elements == 1.0
			int dtt_size, // 8
			boolean bdebug)
	{
		double s = 0.0;
		int indx = 0;
		for (int i = -dtt_size + 1; i < dtt_size; i++){
			int ai = (i < 0)? -i: i;
			for (int j = -dtt_size + 1; j < dtt_size; j++){
				int aj = (j < 0)? -j: j;
				s += kernel[indx++] * window[ai*dtt_size+aj];
			}
		}
		s = 1.0/s;
		int klen = (2*dtt_size-1) * (2*dtt_size-1);
		if (bdebug)		System.out.println("clt_normalize_kernel(): s="+s);
		for (int i = 0; i < klen; i++) {
//******************** Somewhere scale 16 ? ********************
			kernel[i] *= 16*s;
		}
 	}

	void clt_symmetrize_kernel( //
			double []     kernel,      // should be (2*dtt_size-1) * (2*dtt_size-1) +2 size (last 2 are not modified)
			double [][]   sym_kernels, // set of 4 SS, AS, SA, AA kdernels, each dtt_size * dtt_size (may have 5-th with center shift
			final int     dtt_size) // 8
	{
		int in_size = 2*dtt_size-1;
		int dtt_size_m1 = dtt_size - 1;
		int center = dtt_size_m1 * in_size + dtt_size_m1;

		for (int i = 0; i < dtt_size; i++){
			for (int j = 0; j < dtt_size; j++){
				int indx0 = center - i * in_size - j;
				int indx1 = center - i * in_size + j;
				int indx2 = center + i * in_size - j;
				int indx3 = center + i * in_size + j;
				sym_kernels[0][i*dtt_size+j] =                                 0.25*( kernel[indx0] + kernel[indx1] + kernel[indx2] + kernel[indx3]);
				if (j > 0)              sym_kernels[1][i*dtt_size+j-1] =       0.25*(-kernel[indx0] + kernel[indx1] - kernel[indx2] + kernel[indx3]);
				if (i > 0)              sym_kernels[2][(i-1)*dtt_size+j] =     0.25*(-kernel[indx0] - kernel[indx1] + kernel[indx2] + kernel[indx3]);
				if ((i > 0) && (j > 0)) sym_kernels[3][(i-1)*dtt_size+(j-1)] = 0.25*(-kernel[indx0] + kernel[indx1] - kernel[indx2] + kernel[indx3]);
			}
			sym_kernels[1][i*dtt_size + dtt_size_m1] = 0.0;
			sym_kernels[2][dtt_size_m1*dtt_size + i] = 0.0;
			sym_kernels[3][i*dtt_size + dtt_size_m1] = 0.0;
			sym_kernels[3][dtt_size_m1*dtt_size + i] = 0.0;
		}
 	}

	void clt_dtt3_kernel( //
			double [][]   kernels, // set of 4 SS, AS, SA, AA kdernels, each dtt_size * dtt_size (may have 5-th with center shift
			final int     dtt_size, // 8
			DttRad2       dtt)
	{
		if (dtt == null) dtt = new DttRad2(dtt_size);
		for (int quad = 0; quad < 4; quad ++){
			kernels[quad] = dtt.dttt_iiie(kernels[quad], quad, dtt_size);
		}
 	}
/*
	void clt_fill_coord_corr ( // add 6 more items to extra data:  dxc/dx,dyc/dy, dyc/dx, dyc/dy - pixel shift when applied to different center
			// and x0, y0 (which censor pixel this kernel applies to) ? - not needed
			double [][]   kernels, // set of 4 SS, AS, SA, AA kdernels, each dtt_size * dtt_size (may have 5-th with center shift

			)

*/
	public class CltExtra{
		public double data_x   = 0.0; // kernel data is relative to this displacement X (0.5 pixel increments)
		public double data_y   = 0.0; // kernel data is relative to this displacement Y (0.5 pixel increments)
		public double center_x = 0.0; // actual center X (use to find derivatives)
		public double center_y = 0.0; // actual center X (use to find derivatives)
		public double dxc_dx   = 0.0; // add this to data_x per each pixel X-shift relative to the kernel center location
		public double dxc_dy   = 0.0; // same per each Y-shift pixel
		public double dyc_dx   = 0.0;
		public double dyc_dy   = 0.0;

		public CltExtra(){}
		public CltExtra(double [] data)
		{
			data_x   = data[0]; // kernel data is relative to this displacement X (0.5 pixel increments)
			data_y   = data[1]; // kernel data is relative to this displacement Y (0.5 pixel increments)
			center_x = data[2]; // actual center X (use to find derivatives)
			center_y = data[3]; // actual center X (use to find derivatives)
			dxc_dx   = data[4]; // add this to data_x per each pixel X-shift relative to the kernel centger location
			dxc_dy   = data[5]; // same per each Y-shift pixel
			dyc_dx   = data[6];
			dyc_dy   = data[7];
		}
		public double [] getArray()
		{
			double [] rslt = {
					data_x,
					data_y,
					center_x,
					center_y,
					dxc_dx,
					dxc_dy,
					dyc_dx,
					dyc_dy
			};
			return rslt;
		}
	}

	public void offsetKernelSensor(
			double [][] clt_tile, // clt tile, including [4] - metadata
			double dx,
			double dy) {
		CltExtra ce = new CltExtra (clt_tile[4]);
		ce.center_x += dx;
		ce.center_y += dy;
		ce.data_x +=   dx;
		ce.data_y +=   dy;
		clt_tile[4] = ce.getArray();
	}
	public void clt_fill_coord_corr(
			final int               kern_step, // distance between kernel centers, in pixels.
			final double [][][][][] clt_data,
			final int               threadsMax,     // maximal number of threads to launch
			final int               globalDebugLevel)
	{
		final int nChn=clt_data.length;
		final int tilesY=clt_data[0].length;
		final int tilesX=clt_data[0][0].length;
		final int nTilesInChn=tilesX*tilesY;
		final int nTiles=nTilesInChn*nChn;
		final Thread[] threads = newThreadArray(threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				@Override
				public void run() {
					int tileY,tileX,chn;
					for (int nTile = ai.getAndIncrement(); nTile < nTiles; nTile = ai.getAndIncrement()) {
						chn=nTile/nTilesInChn;
						tileY =(nTile % nTilesInChn)/tilesX;
						tileX = nTile % tilesX;
						double s0=0.0, sx=0.0, sx2= 0.0, sy=0.0, sy2= 0.0, sz=0.0, sxz=0.0,
								syz=0.0, sw=0.0, sxw=0.0, syw=0.0;
						for (int dty = -1; dty < 2; dty++){
							int ty = tileY+dty;
							if ((ty >= 0) && (ty < tilesY)){
								for (int dtx = -1; dtx < 2; dtx++){
									int tx = tileX + dtx;
									if ((tx >= 0) && (tx < tilesX)){
										CltExtra ce = new CltExtra (clt_data[chn][ty][tx][4]);
										s0 +=  1;
										sx +=  dtx;
										sx2 += dtx*dtx;
										sy +=  dty;
										sy2 += dty*dty;
										sz  += ce.center_x;
										sxz += dtx * ce.center_x;
										syz += dty * ce.center_x;
										sw  += ce.center_y;
										sxw += dtx * ce.center_y;
										syw += dty * ce.center_y;
									}
								}
							}
						}
						CltExtra ce = new CltExtra (clt_data[chn][tileY][tileX][4]);
						double denom_x = (sx2*s0-sx*sx)*kern_step;
						double denom_y = (sy2*s0-sy*sy)*kern_step;
						ce.dxc_dx= (sxz*s0 - sz*sx)/denom_x;
						ce.dxc_dy= (syz*s0 - sz*sy)/denom_y;
						ce.dyc_dx= (sxw*s0 - sw*sx)/denom_x;
						ce.dyc_dy= (syw*s0 - sw*sy)/denom_y;
						clt_data[chn][tileY][tileX][4] = ce.getArray();
					}
				}
			};
		}
		startAndJoin(threads);
	}

	public class CltTile{
		public double [][] tile = new double[4][]; // 4 CLT tiles
		public double fract_x; // remaining fractional offset X
		public double fract_y; // remaining fractional offset X
	}



// Extract and correct one image tile using kernel data, required result tile and shifts - x and y
// option - align to Bayer (integer shift by even pixels - no need
// input - RBG stack of sparse data
// return
// kernel [0][0] is centered at  (-kernel_step/2,-kernel_step/2)

	public double [] extract_correct_tile( // return a pair of residual offsets
			double [][]         image_data,
			int                 width,       // image width
			double  [][][][][]  clt_kernels, // [color][tileY][tileX][band][pixel]
			double  [][]        clt_tile,    // should be double [4][];
			int                 kernel_step,
			int                 transform_size,
			DttRad2             dtt,
			int                 chn,     // color channel
			double              centerX, // center of aberration-corrected (common model) tile, X
			double              centerY, //
			int                 debugLevel,
			boolean             dbg_no_deconvolution,
			boolean             dbg_transpose,
			boolean []          saturation_imp, // (near) saturated pixels or null
			int []              overexp_all ) // {number of overexposed,  number of all tiles} or null

	{
//		boolean debug_fpga = debugLevel < -9;
		boolean debug_fpga = (debugLevel < -9); //  || (debugLevel == 2);
		boolean debug_gpu = (debugLevel == 2);
		if (debug_fpga) debugLevel = 1;
		if (debug_gpu) debugLevel = 0; // 1; // skip too many images

		int chn_kernel = chn;
		// use zero kernel channel for monochrome images
		if ((clt_kernels != null) && (clt_kernels.length <= chn_kernel)) {
			chn_kernel = clt_kernels.length - 1;
		}

		int chn_img = chn;
		// use zero kernel channel for monochrome images
		if (image_data.length <= chn_img) {
			chn_img = image_data.length - 1;
		}

		// now for mono both image_data and clt_kernel have outer dimension of 1, but chn == 2 (green)

		boolean use_kernels = (clt_kernels != null) && !dbg_no_deconvolution;
		boolean bdebug0 = debugLevel > 0;
		boolean bdebug =  debugLevel > 1;
		double [] residual_shift = new double[2];
		int height = image_data[chn_img].length/width;
		int transform_size2 = 2* transform_size;
//		if (dtt == null) dtt = new DttRad2(transform_size); should have window set up
		double []   tile_in =  new double [4*transform_size*transform_size];

		double px = centerX - transform_size;
		double py = centerY - transform_size;
		int ktileX = 0;
		int ktileY = 0;
		if (use_kernels) {
			// 1. find closest kernel
			ktileX = (int) Math.round(centerX/kernel_step) + 1;
			ktileY = (int) Math.round(centerY/kernel_step) + 1;
			if      (ktileY < 0)                                ktileY = 0;
			else if (ktileY >= clt_kernels[chn_kernel].length)  ktileY = clt_kernels[chn_kernel].length-1;
			if      (ktileX < 0)                                ktileX = 0;
			else if (ktileX >= clt_kernels[chn_kernel][ktileY].length) ktileX = clt_kernels[chn_kernel][ktileY].length-1;
			// extract center offset data stored with each kernel tile
			CltExtra ce = new CltExtra (clt_kernels[chn_kernel][ktileY][ktileX][4]);
			// 2. calculate correction for center of the kernel offset
			double kdx = centerX - (ktileX -1 +0.5) *  kernel_step; // difference in pixel
			double kdy = centerY - (ktileY -1 +0.5) *  kernel_step;
			// 3. find top-left corner of the
			// check signs, ce.data_x is "-" as if original kernel was shifted to "+" need to take pixel sifted "-"
			// same with extra shift
			px = centerX - transform_size - (ce.data_x + ce.dxc_dx * kdx + ce.dxc_dy * kdy) ; // fractional left corner
			py = centerY - transform_size - (ce.data_y + ce.dyc_dx * kdx + ce.dyc_dy * kdy) ; // fractional top corner
			if (debug_gpu) {
				System.out.println("========= Color channel "+chn+" , kernel channel "+chn_kernel+" input image channel="+chn_img+" =============");
				System.out.println("ce.data_x="+ce.data_x+", ce.data_y="+ce.data_y);
				System.out.println("ce.center_x="+ce.center_x+", ce.center_y="+ce.center_y);
				System.out.println("ce.dxc_dx="+ce.dxc_dx+", ce.dxc_dy="+ce.dxc_dy);
				System.out.println("ce.dyc_dx="+ce.dyc_dx+", ce.dyc_dy="+ce.dyc_dy);
				System.out.println("centerX="+centerX+", centerY="+centerY);
				System.out.println("px="+px+", py="+py);
				System.out.println("ktileX="+ktileX+", ktileY="+ktileY);
				System.out.println("kdx="+kdx+", kdy="+kdy);
			}
		}else {
//			System.out.println("Skipping kernels!!!"); // Happens when using macro_mode, should not happen otherwise
		}
		if (bdebug0){
			System.out.print(px+"\t"+py+"\t");
		}
		// Was wrong rounding, fractional part gets to +0.5
		int ctile_left = (int) -Math.round(-px);
		int ctile_top =  (int) -Math.round(-py);
		residual_shift[0] = -(px - ctile_left);
		residual_shift[1] = -(py - ctile_top);

		if (debug_gpu) {
			System.out.println("ctile_left="+ctile_left+", ctile_top="+ctile_top);
			System.out.println("residual_shift[0]="+residual_shift[0]+", residual_shift[1]="+residual_shift[1]);
		}


		// 4. Verify the tile fits in image and use System.arraycopy(sym_conv, 0, tile_in, 0, n2*n2) to copy data to tile_in
		// if does not fit - extend by duplication? Or just use 0?
		if ((ctile_left >= 0) && (ctile_left < (width - transform_size2)) &&
				(ctile_top >= 0) && (ctile_top < (height - transform_size2))) {
			for (int i = 0; i < transform_size2; i++){
				System.arraycopy(image_data[chn_img], (ctile_top + i) * width + ctile_left, tile_in, transform_size2 * i, transform_size2);
			}
		} else { // copy by 1
			for (int i = 0; i < transform_size2; i++){
				int pi = ctile_top + i;
				if      (pi < 0)       pi &= 1;
				else if (pi >= height) pi = height - 2 + (pi & 1);
				for (int j = 0; j < transform_size2; j++){
					int pj = ctile_left + j;
					if      (pj < 0)      pj &= 1;
					else if (pj >= width) pj = width - 2 + (pj & 1);
					tile_in[transform_size2 * i + j] = image_data[chn_img][pi * width + pj];
				}
			}
		}
		if (debug_gpu) {
			System.out.println("---Image tile for color="+chn_img+"---");
			for (int i = 0; i < transform_size2; i++) {
				for (int j = 0; j < transform_size2; j++) {
					System.out.print(String.format("%10.5f ", tile_in[transform_size2 * i + j]));
				}
				System.out.println();
			}
		}


		if (debug_fpga){ // show extended tile, all colors
			System.out.println("\nFull Bayer fpga tile data");
			int lt = (FPGA_TILE_SIZE - transform_size2)/2;
			double [][] fpga_tile = new double [3][FPGA_TILE_SIZE * FPGA_TILE_SIZE];
			for (int fpga_chn = 0; fpga_chn < 3; fpga_chn++){
				for (int i = 0; i < FPGA_TILE_SIZE; i++){
					System.arraycopy(image_data[fpga_chn], ((ctile_top - lt) + i) * width + (ctile_left - lt), fpga_tile[fpga_chn], FPGA_TILE_SIZE * i, FPGA_TILE_SIZE);
				}
			}
			int id = (1 << (FPGA_PIXEL_BITS - 9)); // 8
			for (int i = 0; i < FPGA_TILE_SIZE*FPGA_TILE_SIZE; i++) {
				double d = 0.0;
				for (int fpga_chn = 0; fpga_chn < 3; fpga_chn++){
					d +=  fpga_tile[fpga_chn][i];
				}
				System.out.print(String.format("%4x ",(int) Math.round(id * d)));
				if (((i+1) %FPGA_TILE_SIZE) == 0) {
					System.out.println();
				}
			}
		}


		if (((chn == GREEN_CHN) || isMonochrome()) && (saturation_imp != null)) {
			//overexp_all
			if ((ctile_left >= 0) && (ctile_left < (width - transform_size2)) &&
					(ctile_top >= 0) && (ctile_top < (height - transform_size2))) {
				for (int i = 0; i < transform_size2; i++){
					int indx = (ctile_top + i) * width + ctile_left;
					for (int j = 0; j < transform_size2; j++) {
						if (saturation_imp[indx++]) {
							overexp_all[0] ++;
						}
					}
				}
				overexp_all[1] += transform_size2 * transform_size2;
			} else { // copy by 1
				for (int i = 0; i < transform_size2; i++){
					int pi = ctile_top + i;
					if      (pi < 0)       pi &= 1;
					else if (pi >= height) pi = height - 2 + (pi & 1);
					for (int j = 0; j < transform_size2; j++){
						int pj = ctile_left + j;
						if      (pj < 0)      pj &= 1;
						else if (pj >= width) pj = width - 2 + (pj & 1);
						if (saturation_imp[pi * width + pj]) {
							overexp_all[0] ++;
						}
					}
				}
				overexp_all[1] += transform_size2 * transform_size2;
			}
		}
		if (debug_fpga) {
			System.out.println("debug_fpga: residual_shift[0]="+residual_shift[0]+", residual_shift[1]="+residual_shift[1]);
			int ishx, ishy;
			ishx = (int) Math.round((1 << (FPGA_SHIFT_BITS)) * residual_shift[0]);
			ishy = (int) Math.round((1 << (FPGA_SHIFT_BITS)) * residual_shift[1]);
			if (ishx >= (1 << (FPGA_SHIFT_BITS-1))) ishx = (1 << (FPGA_SHIFT_BITS-1)) - 1;
			if (ishy >= (1 << (FPGA_SHIFT_BITS-1))) ishy = (1 << (FPGA_SHIFT_BITS-1)) - 1;
			if (ishx < -(1 << (FPGA_SHIFT_BITS-1))) ishx = -(1 << (FPGA_SHIFT_BITS-1));
			if (ishy < -(1 << (FPGA_SHIFT_BITS-1))) ishy = -(1 << (FPGA_SHIFT_BITS-1));
			residual_shift[0] = ishx * (1.0/(1 << (FPGA_SHIFT_BITS)));
			residual_shift[1] = ishy * (1.0/(1 << (FPGA_SHIFT_BITS)));
			System.out.println("rounded: residual_shift[0]="+residual_shift[0]+", residual_shift[1]="+residual_shift[1]);
			double [] fpga_pix_lim = {0.0,0.0};
			for (int i = 0; i < 256; i++){
				if (tile_in[i] > fpga_pix_lim[0]) fpga_pix_lim[0] = tile_in[i];
				if (tile_in[i] < fpga_pix_lim[1]) fpga_pix_lim[1] = tile_in[i];
			}
			System.out.println(String.format("\n// Pixels input range: %f ... %f", fpga_pix_lim[1], fpga_pix_lim[0]));
			System.out.println(String.format("%x // shift_x, %d bits",ishx & ((1 << (FPGA_SHIFT_BITS)) - 1),FPGA_SHIFT_BITS));
			System.out.println(String.format("%x // shift_y, %d bits",ishy & ((1 << (FPGA_SHIFT_BITS)) - 1),FPGA_SHIFT_BITS));
			System.out.println(String.format("%x // bayer",15));
			int id = (1 << (FPGA_PIXEL_BITS - 9)); // 8
			for (int row = 0; row <16; row++){
				for (int col = 0; col <16; col++){
					System.out.print(String.format("%4x ",(int) Math.round(id * tile_in[row*16 + col])));
				}
				System.out.println();
			}
		}

		// Fold and transform
		double [][][] fold_coeff = null;
		if (!dbg_transpose){
			fold_coeff = dtt.get_shifted_fold_2d ( // get_shifted_fold_2d(
					transform_size,
					residual_shift[0],
					residual_shift[1],
					0); // debug level
		}

		if (debug_fpga) {
			System.out.println("debug_fpga: residual_shift[0]="+residual_shift[0]+", residual_shift[1]="+residual_shift[1]);
			System.out.println("Signs table (per mode, per index - bitstring of variants, 0 - positive, 1 - negative");
			for (int dct_mode = 0; dct_mode <4; dct_mode++) {
				for (int i = 0; i < 64; i++){
					int d = 0;
					for (int b = 0; b < 4; b++){
						if (fold_coeff[dct_mode][i][b] < 0){
							d |= (1 << b);
						}
					}
					System.out.print(String.format("%x ",d));
					if ((i % 16) == 15){
						System.out.println();
					}
				}
			}
			System.out.println("Absolute values, shoud be the same for each of 4 modes");
			for (int dct_mode = 0; dct_mode <4; dct_mode++) {
				for (int i = 0; i < 64; i++){
					for (int b = 0; b < 4; b++){
						int d = (int) Math.round(((1 << FPGA_WND_BITS) -1)* Math.abs(fold_coeff[dct_mode][i][b]));
						System.out.print(String.format("%5x ",d & ((1 << (FPGA_WND_BITS)) - 1)));
					}
					if ((i % 4) == 3){
						System.out.println();
					}
				}
				System.out.println();
			}

			double [][][] fold_coeff_direct = dtt.get_shifted_fold_2d_direct ( // get_shifted_fold_2d(
						transform_size,
						residual_shift[0],
						residual_shift[1],
						(1 << FPGA_WND_BITS) -1); // debug level - use as scale

			System.out.println("Direct sin table");
			for (int i = 0; i < 64; i++){
				for (int b = 0; b < 4; b++){
					int d = (int) Math.round(((1 << FPGA_WND_BITS) -1)* Math.abs(fold_coeff_direct[0][i][b])); // dct_mode=0
					System.out.print(String.format("%5x ",d & ((1 << (FPGA_WND_BITS)) - 1)));
				}
				if ((i % 4) == 3){
					System.out.println();
				}
			}
			System.out.println();

			double [][][] fold_coeff_old = dtt.get_shifted_fold_2d ( // get_shifted_fold_2d(
					transform_size,
					residual_shift[0],
					residual_shift[1],
					(1 << FPGA_WND_BITS) -1); // debug level - use as scale

			System.out.println("Direct sin table");
			for (int i = 0; i < 64; i++){
				for (int b = 0; b < 4; b++){
					int d = (int) Math.round(((1 << FPGA_WND_BITS) -1)* Math.abs(fold_coeff_old[0][i][b])); // dct_mode=0
					System.out.print(String.format("%5x ",d & ((1 << (FPGA_WND_BITS)) - 1)));
				}
				if ((i % 4) == 3){
					System.out.println();
				}
			}
			System.out.println();

			System.out.println("Diff: new - old");
			for (int i = 0; i < 64; i++){
				for (int b = 0; b < 4; b++){
					int d0 = (int) Math.round(((1 << FPGA_WND_BITS) -1)* Math.abs(fold_coeff[0][i][b])); // dct_mode=0
					int d = (int) Math.round(((1 << FPGA_WND_BITS) -1)* Math.abs(fold_coeff_direct[0][i][b])); // dct_mode=0
					System.out.print(String.format("%5d ",d - d0));
				}
				if ((i % 4) == 3){
					System.out.println();
				}
			}
			System.out.println();



			System.out.println("\nFold index");
			int [][] fpga_fi = dtt.getFoldIndex();
			for (int i = 0; i < 64; i++){
				for (int k = 0; k <4; k++) {
				        System.out.print(String.format("%02x ", fpga_fi[i][k]));
				}
				System.out.print("  ");
				if (i%8 == 7) System.out.println();
			}
			System.out.println();

			// Show for different Bayer patterns
			int [] bayer_patterns = {0x1, 0x2, 0x4, 0x8, 0x9, 0x6};
			for (int bp:bayer_patterns){
				System.out.println("Pattern (row/col) "+bp+":");
				System.out.println("| "+(((bp & 1) !=0) ? "X ":"  ")+(((bp & 2) !=0) ? "X ":"  ")+"|");
				System.out.println("| "+(((bp & 4) !=0) ? "X ":"  ")+(((bp & 8) !=0) ? "X ":"  ")+"|");
				for (int i = 0; i < 64; i++){
					for (int k = 0; k <4; k++) {
						int row = (fpga_fi[i][k] >> 4);
						int col = (fpga_fi[i][k] & 0xf);
						int indx = (row & 1) + 2 * (col & 1);
						if (((1 << indx) & bp) != 0) {
							System.out.print(String.format("%2x ", fpga_fi[i][k]));
						} else {
							System.out.print(" . ");
						}
					}
					System.out.print("  ");
					if (i%8 == 7) System.out.println();
				}
				System.out.println();
			}

			for (int bp:bayer_patterns){
				System.out.println("Pattern (mode bits) "+bp+":");
				System.out.println("| "+(((bp & 1) !=0) ? "X ":"  ")+(((bp & 2) !=0) ? "X ":"  ")+"|");
				System.out.println("| "+(((bp & 4) !=0) ? "X ":"  ")+(((bp & 8) !=0) ? "X ":"  ")+"|");
				for (int i = 0; i < 64; i++){
					for (int k = 0; k < 4; k++) {
						int row = (fpga_fi[i][k] >> 4);
						int col = (fpga_fi[i][k] & 0xf);
						int indx = (row & 1) + 2 * (col & 1);
						int d = 0;
						for (int dct_mode = 0; dct_mode <4; dct_mode++) {
							if (fold_coeff[dct_mode][i][k] < 0){
								d |= (1 << dct_mode);
							}
						}

						if (((1 << indx) & bp) != 0) {
							System.out.print(String.format("%02x ", d));
						} else {
							System.out.print(" . ");
						}
					}
					System.out.print("  ");
					if (i%8 == 7) System.out.println();
				}
				System.out.println();
			}


			double [][] fpga_w_u = new double [4][256];
			double [][] fpga_w_s = new double [4][256];
			for (int dct_mode = 0; dct_mode <4; dct_mode++) {
				for (int i = 0; i < 64; i++){
					for (int b=0; b < 4; b++){
						fpga_w_u [dct_mode][fpga_fi[i][b]] = Math.abs(fold_coeff[dct_mode][i][b]);
						fpga_w_s [dct_mode][fpga_fi[i][b]] = fold_coeff[dct_mode][i][b];
					}
				}
			}
			ShowDoubleFloatArrays sdfa_instance = new ShowDoubleFloatArrays(); // just for debugging?
			String [] titles = {"CC","SC","CS","SS"};
			sdfa_instance.showArrays(fpga_w_s,  2 * transform_size, 2 * transform_size, true, "fpga_w_s_x"+ctile_left+"_y"+ctile_top, titles);
			sdfa_instance.showArrays(fpga_w_u,  2 * transform_size, 2 * transform_size, true, "fpga_w_u_x"+ctile_left+"_y"+ctile_top, titles);



		} //if (debug_fpga)


		if (!debug_fpga) {
			for (int dct_mode = 0; dct_mode < 4; dct_mode++) {
				if (fold_coeff != null){
					clt_tile[dct_mode] = dtt.fold_tile (tile_in, transform_size, dct_mode, fold_coeff); // DCCT, DSCT, DCST, DSST
				} else {
					clt_tile[dct_mode] = dtt.fold_tile (tile_in, transform_size, dct_mode); // DCCT, DSCT, DCST, DSST
				}
				if (debug_gpu) {
					System.out.println("=== Image tile folded for color="+chn+", dct_mode="+dct_mode+" ===");
					for (int i = 0; i < transform_size; i++) {
						for (int j = 0; j < transform_size; j++) {
							System.out.print(String.format("%10.5f ", clt_tile[dct_mode][transform_size * i + j]));
						}
						System.out.println();
					}
				}

				clt_tile[dct_mode] = dtt.dttt_iv   (clt_tile[dct_mode], dct_mode, transform_size);

				if (debug_gpu) {
					System.out.println("=== Image tile DTT converted for color="+chn+", dct_mode="+dct_mode+" ===");
					for (int i = 0; i < transform_size; i++) {
						for (int j = 0; j < transform_size; j++) {
							System.out.print(String.format("%10.5f ", clt_tile[dct_mode][transform_size * i + j]));
						}
						System.out.println();
					}
				}
			}
		}

		if (debug_fpga){
			for (int dct_mode = 0; dct_mode <4; dct_mode++) {
				if (fold_coeff != null){
					clt_tile[dct_mode] = dtt.fold_tile_debug (tile_in, transform_size, dct_mode, fold_coeff); // DCCT, DSCT, DCST, DSST
				} else {
					clt_tile[dct_mode] = dtt.fold_tile (tile_in, transform_size, dct_mode); // DCCT, DSCT, DCST, DSST
				}
			}



			double [] fpga_dtt_lim = {0.0,0.0};
			for (int dct_mode = 0; dct_mode <4; dct_mode++) {
				for (int i = 0; i < 64; i++){
					if (clt_tile[dct_mode][i] > fpga_dtt_lim[0]) fpga_dtt_lim[0] = clt_tile[dct_mode][i];
					if (clt_tile[dct_mode][i] < fpga_dtt_lim[1]) fpga_dtt_lim[1] = clt_tile[dct_mode][i];
				}
			}
			System.out.println(String.format("// DTT input range: %f ... %f", fpga_dtt_lim[1], fpga_dtt_lim[0]));
			double scale = (1 << (FPGA_DTT_IN - 9)); //  -1;
			scale *= 1.0 *((1 << FPGA_WND_BITS) -1) / (1 << FPGA_WND_BITS);
			scale *= 1.0 *((1 << FPGA_WND_BITS) -1) / (1 << FPGA_WND_BITS);
			for (int dct_mode = 0; dct_mode <4; dct_mode++) {
				for (int i = 0; i < 64; i++){
					int id = (int) Math.round(scale * clt_tile[dct_mode][i]);
					System.out.print(String.format("%7x ", id & ((1 << 25) -1)));
					if ((i % 8) == 7) System.out.println();
				}
				System.out.println();
			}
			System.out.println();

			for (int dct_mode = 0; dct_mode <4; dct_mode++) {
				System.out.println("Color= 2? clt_tile[dct_mode] = dtt.dttt_iv(..., scale="	+scale);
				clt_tile[dct_mode] = dtt.dttt_iv   (clt_tile[dct_mode], dct_mode, transform_size, scale, ((1 << 25) -1)); // debug level
			}
			fpga_dtt_lim[0] = 0.0;
			fpga_dtt_lim[1] = 0.0;
			for (int dct_mode = 0; dct_mode <4; dct_mode++) {
				for (int i = 0; i < 64; i++){
					if (clt_tile[dct_mode][i] > fpga_dtt_lim[0]) fpga_dtt_lim[0] = clt_tile[dct_mode][i];
					if (clt_tile[dct_mode][i] < fpga_dtt_lim[1]) fpga_dtt_lim[1] = clt_tile[dct_mode][i];
				}
			}
			System.out.println(String.format("// DTT output range: %f ... %f", fpga_dtt_lim[1], fpga_dtt_lim[0]));
			// scale = (1 << (FPGA_DTT_IN - 9)); //  -1;
			for (int dct_mode = 0; dct_mode <4; dct_mode++) {
				for (int i = 0; i < 64; i++){
					int id = (int) Math.round(scale * clt_tile[dct_mode][i]);
					System.out.print(String.format("%7x ", id & ((1 << 25) -1)));
					if ((i % 8) == 7) System.out.println();
				}
				System.out.println();
			}
			System.out.println();
			System.out.println("Testing symmetry of checkerboard patterns");
			for (int dct_mode = 0; dct_mode < 2; dct_mode++) {
				for (int i = 0; i < 64; i++){
					int id = (int) Math.round(scale * clt_tile[dct_mode][i]);
					int id1 = (int) Math.round(scale * clt_tile[3-dct_mode][63-i]);
					System.out.print(String.format("%7x ", (id-id1) & ((1 << 25) -1)));
					if ((i % 8) == 7) System.out.println();
				}
				System.out.println();
			}
			System.out.println();
			System.out.println("Testing antisymmetry of checkerboard patterns");
			for (int dct_mode = 0; dct_mode < 2; dct_mode++) {
				for (int i = 0; i < 64; i++){
					int id = (int) Math.round(scale * clt_tile[dct_mode][i]);
					int id1 = (int) Math.round(scale * clt_tile[3-dct_mode][63-i]);
					System.out.print(String.format("%7x ", (id+id1) & ((1 << 25) -1)));
					if ((i % 8) == 7) System.out.println();
				}
				System.out.println();
			}
			System.out.println();
		} // end of if (debug_fpga){



		if (bdebug0) {
			ShowDoubleFloatArrays sdfa_instance = new ShowDoubleFloatArrays(); // just for debugging?
			sdfa_instance.showArrays(tile_in,  transform_size2, transform_size2, "tile_in_x"+ctile_left+"_y"+ctile_top);
			String [] titles = {"CC","SC","CS","SS"};
			sdfa_instance.showArrays(clt_tile,  transform_size, transform_size, true, "clt_x"+ctile_left+"_y"+ctile_top, titles);
		}
		// deconvolve with kernel
		if (use_kernels) {
			double [][] ktile = clt_kernels[chn_kernel][ktileY][ktileX];
			if (debug_gpu) {
				System.out.println("=== kernel tile for color="+chn_kernel+" (color = "+chn+") ===");
				for (int dct_mode = 0; dct_mode < 4; dct_mode++) {
					System.out.println("dct_mode="+dct_mode);
					for (int i = 0; i < transform_size; i++) {
						for (int j = 0; j < transform_size; j++) {
							System.out.print(String.format("%10.5f ", ktile[dct_mode][transform_size * i + j]));
						}
						System.out.println();
					}
				}
			}
			convolve_tile(
					clt_tile,        // double [][]     data,    // array [transform_size*transform_size], will be updated  DTT4 converted
					ktile,           // double [][]     kernel,  // array [4][transform_size*transform_size]  DTT3 converted
					transform_size,
					bdebug);
			if (debug_gpu) {
				System.out.println("=== convolved tile for color="+chn+" ===");
				for (int dct_mode = 0; dct_mode < 4; dct_mode++) {
					System.out.println("dct_mode="+dct_mode);
					for (int i = 0; i < transform_size; i++) {
						for (int j = 0; j < transform_size; j++) {
							System.out.print(String.format("%10.5f ", clt_tile[dct_mode][transform_size * i + j]));
						}
						System.out.println();
					}
				}
			}
		}
		if (bdebug) {
			ShowDoubleFloatArrays sdfa_instance = new ShowDoubleFloatArrays(); // just for debugging?
			String [] titles = {"CC","SC","CS","SS"};
			sdfa_instance.showArrays(clt_tile,  transform_size, transform_size, true, "acorr_x"+ctile_left+"_y"+ctile_top, titles);
		}
		return residual_shift;
	}

//	public
	public void convolve_tile(
			double [][]     data,    // array [transform_size*transform_size], will be updated  DTT4 converted
			double [][]     kernel,  // array [4][transform_size*transform_size]  DTT3 converted
			int             transform_size,
			boolean         bdebug) // externally decoded debug tile
//			boolean         dbg_transpose)

	{
		/* Direct matrix Z1: X2 ~= Z1 * Shift
		 * {{+cc  -sc  -cs  +ss},
		 *  {+sc  +cc  -ss  -cs},
		 *  {+cs  -ss  +cc  -sc},
		 *  {+ss  +cs  +sc  +cc}}
		 *
		 * T= transp({cc, sc, cs, ss})
		 */
		/*
		final int [][] zi =
			{{ 0, -1, -2,  3},
			 { 1,  0, -3, -2},
			 { 2, -3,  0, -1},
			 { 3,  2,  1,  0}};
		final int [][] zi =
			{{ 0,  1,  2,  3},
			 {-1,  0, -3,  2},
			 {-2, -3,  0,  1},
			 { 3, -2, -1,  0}};
		 */
		// opposite sign from correlation
		final int [][] zi =	{ //
				{ 0, -1, -2,  3},
				{ 1,  0, -3, -2},
				{ 2, -3,  0, -1},
				{ 3,  2,  1,  0}};

		final int transform_len = transform_size * transform_size;
		final double [][] rslt = new double[4][transform_len];
		for (int i = 0; i < transform_len; i++) {
			for (int n = 0; n<4; n++){
				rslt[n][i] = 0;
				for (int k=0; k<4; k++){
					if (zi[n][k] < 0)
						rslt[n][i] -= data[-zi[n][k]][i] * kernel[k][i];
					else
						rslt[n][i] += data[ zi[n][k]][i] * kernel[k][i];
				}
			}
		}
		if (bdebug) {
			ShowDoubleFloatArrays sdfa_instance = new ShowDoubleFloatArrays(); // just for debugging?
			String [] titles = {"CC","SC","CS","SS"};
			double [][] dbg_kern = {kernel[0],kernel[1],kernel[2],kernel[3]};
			sdfa_instance.showArrays(data,     transform_size, transform_size, true, "image_data", titles);
			sdfa_instance.showArrays(dbg_kern, transform_size, transform_size, true, "kernel",     titles);
			sdfa_instance.showArrays(rslt,     transform_size, transform_size, true, "aber_corr",  titles);
		}
		for (int n = 0; n<4; n++){
			data[n] = rslt[n];
		}
	}

	public void fract_shift(    // fractional shift in transform domain. Currently uses sin/cos - change to tables with 2? rotations
		double  [][]  clt_tile,
		int           transform_size,
		double        shiftX,
		double        shiftY,
		boolean       bdebug)
	{
        boolean debug_images = false;
		int transform_len = transform_size*transform_size;
		double [] cos_hor =  new double [transform_len];
		double [] sin_hor =  new double [transform_len];
		double [] cos_vert = new double [transform_len];
		double [] sin_vert = new double [transform_len];
		for (int i = 0; i < transform_size; i++){
			double ch = Math.cos((i+0.5)*Math.PI*shiftX/transform_size);
			double sh = Math.sin((i+0.5)*Math.PI*shiftX/transform_size);
			double cv = Math.cos((i+0.5)*Math.PI*shiftY/transform_size);
			double sv = Math.sin((i+0.5)*Math.PI*shiftY/transform_size);
			for (int j = 0; j < transform_size; j++){
				int iv = transform_size * j + i; // 2d DTT results are stored transposed!
				int ih = transform_size * i + j;
				cos_hor[ih] = ch;
				sin_hor[ih] = sh;
				cos_vert[iv] = cv;
				sin_vert[iv] = sv;
			}
		}
		if (bdebug) {
			System.out.println("cos_hor , shift_hor = "+shiftX);
			for (int irow = 0; irow < transform_size; irow++) {
				for (int jcol = 0; jcol < transform_size; jcol++) {
					System.out.print(String.format("%10.5f ", cos_hor[transform_size * irow + jcol]));
				}
				System.out.println();
			}
			System.out.println("\nsin_hor , shift_hor = "+shiftX);
			for (int irow = 0; irow < transform_size; irow++) {
				for (int jcol = 0; jcol < transform_size; jcol++) {
					System.out.print(String.format("%10.5f ", sin_hor[transform_size * irow + jcol]));
				}
				System.out.println();
			}
			System.out.println("cos_vert , shift_vert = "+shiftY);
			for (int irow = 0; irow < transform_size; irow++) {
				for (int jcol = 0; jcol < transform_size; jcol++) {
					System.out.print(String.format("%10.5f ", cos_vert[transform_size * irow + jcol]));
				}
				System.out.println();
			}
			System.out.println("\nsin_vert , shift_vert = "+shiftY);
			for (int irow = 0; irow < transform_size; irow++) {
				for (int jcol = 0; jcol < transform_size; jcol++) {
					System.out.print(String.format("%10.5f ", sin_vert[transform_size * irow + jcol]));
				}
				System.out.println();
			}
			System.out.println();
		}
		if (bdebug && debug_images){
			ShowDoubleFloatArrays sdfa_instance = new ShowDoubleFloatArrays(); // just for debugging?
			String [] titles = {"cos_hor","sin_hor","cos_vert","sin_vert"};
			double [][] cs_dbg = {cos_hor, sin_hor, cos_vert, sin_vert};
			sdfa_instance.showArrays(cs_dbg,  transform_size, transform_size, true, "shift_cos_sin", titles);
		}
		double [][] tmp_tile = new double [4][transform_len];
		// Horizontal shift CLT tiled data is stored in transposed way (horizontal - Y, vertical X)
		for (int i = 0; i < cos_hor.length; i++) {
			tmp_tile[0][i] = clt_tile[0][i] * cos_hor[i] - clt_tile[1][i] * sin_hor[i];
			tmp_tile[1][i] = clt_tile[1][i] * cos_hor[i] + clt_tile[0][i] * sin_hor[i] ;

			tmp_tile[2][i] = clt_tile[2][i] * cos_hor[i]  - clt_tile[3][i] * sin_hor[i];
			tmp_tile[3][i] = clt_tile[3][i] * cos_hor[i]  + clt_tile[2][i] * sin_hor[i] ;
		}
		if (bdebug) {
			System.out.println("---Shifted image tile horizontally, shift_hor = "+shiftX);
			for (int dct_mode=0; dct_mode < 4; dct_mode++ ) {
				System.out.println("dct_mode="+dct_mode);
				for (int irow = 0; irow < transform_size; irow++) {
					for (int jcol = 0; jcol < transform_size; jcol++) {
						System.out.print(String.format("%10.5f ", tmp_tile[dct_mode][transform_size * irow + jcol]));
					}
					System.out.println();
				}
				System.out.println();
			}
		}
		// Vertical shift (back to original array)
		for (int i = 0; i < cos_hor.length; i++) {
			clt_tile[0][i] = tmp_tile[0][i] * cos_vert[i] - tmp_tile[2][i] * sin_vert[i];
			clt_tile[2][i] = tmp_tile[2][i] * cos_vert[i] + tmp_tile[0][i] * sin_vert[i];

			clt_tile[1][i] =                      tmp_tile[1][i] * cos_vert[i] - tmp_tile[3][i] * sin_vert[i];
			clt_tile[3][i] = tmp_tile[3][i] * cos_vert[i] + tmp_tile[1][i] * sin_vert[i];
		}
	}



	public double [][][][] mdctScale(
			final ImageStack                                 imageStack,
			final int                                        subcamera, // not needed
			final EyesisCorrectionParameters.DCTParameters   dctParameters, //
			final int                                        threadsMax, // maximal step in pixels on the maxRadius for 1 angular step (i.e. 0.5)
			final int                                        debugLevel,
			final boolean                                    updateStatus) // update status info

	{
	  	  if (imageStack==null) return null;
		  final int imgWidth=imageStack.getWidth();
		  final int nChn=imageStack.getSize();
		  double [][][][] dct_data = new double [nChn][][][];
		  float [] fpixels;
		  int i,chn; //tileX,tileY;
		  /* find number of the green channel - should be called "green", if none - use last */
		  // Extract float pixels from inage stack, convert each to double


		  for (chn=0;chn<nChn;chn++) {
			  fpixels= (float[]) imageStack.getPixels(chn+1);
			  double[] dpixels = new double[fpixels.length];
			  for (i = 0; i <fpixels.length;i++) dpixels[i] = fpixels[i];
			  // convert each to DCT tiles
			  dct_data[chn] =lapped_dct_scale(
						dpixels,
						imgWidth,
						dctParameters.dct_size,
						(int) Math.round(dctParameters.dbg_src_size),
						dctParameters.dbg_fold_scale,
						dctParameters.dbg_fold_scale,
						0, //     dct_mode,    // 0: dct/dct, 1: dct/dst, 2: dst/dct, 3: dst/dst
						dctParameters.dct_window, // final int       window_type,
						chn,
						dctParameters.tileX,
						dctParameters.tileY,
						dctParameters.dbg_mode,
						threadsMax,  // maximal number of threads to launch
						debugLevel);
		  }
		return dct_data;
	}



	public double [][][] lapped_dct_scale( // scale image to 8/9 size in each direction
			final double [] dpixels,
			final int       width,
			final int       dct_size,
			final int       src_size,    // source step (== dct_size - no scale, == 9 - shrink, ==7 - expand
			final double    scale_hor,
			final double    scale_vert,
			final int       dct_mode,    // 0: dct/dct, 1: dct/dst, 2: dst/dct, 3: dst/dst
			final int       window_type,
			final int       color,
			final int       debug_tileX,
			final int       debug_tileY,
			final int       debug_mode,
			final int       threadsMax,  // maximal number of threads to launch
			final int       globalDebugLevel)
	{
		final int height=dpixels.length/width;
		final int n2 = dct_size * 2;

		final int tilesX = (width - n2) / src_size + 1;
		final int tilesY = (height - n2) / src_size + 1;

		final int nTiles=tilesX*tilesY;
		final double [][][] dct_data = new double[tilesY][tilesX][dct_size*dct_size];
		final Thread[] threads = newThreadArray(threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		for (int tileY = 0; tileY < tilesY; tileY++){
			for (int tileX = 0; tileX < tilesX; tileX++){
				for (int i=0; i<dct_data[tileY][tileX].length;i++) dct_data[tileY][tileX][i]= 0.0; // actually not needed, Java initializes arrays
			}
		}
		double [] dc = new double [dct_size*dct_size];
		for (int i = 0; i<dc.length; i++) dc[i] = 1.0;
//		DttRad2 dtt0 = new DttRad2(dct_size);
//		dtt0.set_window(window_type);
//		final double [] dciii = dtt0.dttt_iii  (dc, dct_size);
//		final double [] dciiie = dtt0.dttt_iiie  (dc, 0, dct_size);


		if (globalDebugLevel > 0) {
			System.out.println("lapped_dctdc(): width="+width+" height="+height);
		}

		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				@Override
				public void run() {
					DttRad2 dtt = new DttRad2(dct_size);
					dtt.set_window(window_type);
					double [] tile_in = new double[4*dct_size * dct_size];
					double [] tile_folded;
					double [] tile_out; // = new double[dct_size * dct_size];
					int tileY,tileX;
					double [][][] fold_k =  dtt.get_fold_2d(
//							int n,
							scale_hor,
							scale_vert
							);
//					double [] tile_out_copy = null;
//					showDoubleFloatArrays sdfa_instance = new showDoubleFloatArrays(); // just for debugging?
					for (int nTile = ai.getAndIncrement(); nTile < nTiles; nTile = ai.getAndIncrement()) {
						tileY = nTile/tilesX;
						tileX = nTile - tileY * tilesX;
						//readDCTKernels() debugLevel = 1 kernels[0].size = 8 kernels[0].img_step = 16 kernels[0].asym_nonzero = 4 nColors = 3 numVert = 123 numHor =  164
						// no aberration correction, just copy data
						for (int i = 0; i < n2;i++){
							System.arraycopy(dpixels, (tileY*width+tileX)*src_size + i*width, tile_in, i*n2, n2);
						}
						tile_folded=dtt.fold_tile(tile_in, dct_size, 0, fold_k); // DCCT
						tile_out=dtt.dttt_iv  (tile_folded, dct_mode, dct_size);

//						if ((tileY == debug_tileY) && (tileX == debug_tileX) && (color == 2)) {
//							tile_out_copy = tile_out.clone();
//						}

						System.arraycopy(tile_out, 0, dct_data[tileY][tileX], 0, tile_out.length);
					}
				}
			};
		}
		startAndJoin(threads);
		return dct_data;
	}

	public void dct_scale(
			final double  scale_hor,  // < 1.0 - enlarge in dct domain (shrink in time/space)
			final double  scale_vert, // < 1.0 - enlarge in dct domain (shrink in time/space)
			final boolean normalize, // preserve weighted dct values
			final double [][][] dct_data,
			final int       debug_tileX,
			final int       debug_tileY,
			final int       threadsMax,     // maximal number of threads to launch
			final int       globalDebugLevel)
	{
		final int tilesY=dct_data.length;
		final int tilesX=dct_data[0].length;
		final int nTiles=tilesX*tilesY;
		final int dct_size = (int) Math.round(Math.sqrt(dct_data[0][0].length));
		final int dct_len = dct_size*dct_size;
		final double [] norm_sym_weights = new double [dct_size*dct_size];
		for (int i = 0; i < dct_size; i++){
			for (int j = 0; j < dct_size; j++){
				double d = 	Math.cos(Math.PI*i/(2*dct_size))*Math.cos(Math.PI*j/(2*dct_size));
				if (i > 0) d*= 2.0;
				if (j > 0) d*= 2.0;
				norm_sym_weights[i*dct_size+j] = d;
			}
		}
		final Thread[] threads = newThreadArray(threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);

		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				@Override
				public void run() {
					int tileY,tileX;
					double [] dct1 = new double [dct_size*dct_size];
					double [] dct;
					double [][] bidata = new double [2][2];
					int dct_m1 = dct_size - 1;
					int dct_m2 = dct_size - 2;
					for (int nTile = ai.getAndIncrement(); nTile < nTiles; nTile = ai.getAndIncrement()) {
						tileY = nTile/tilesX;
						tileX = nTile - tileY * tilesX;
						dct = dct_data[tileY][tileX];
						double sum_orig=0;
						if (normalize) {
							for (int i = 0; i < dct_len; i++){
								sum_orig += dct[i] *norm_sym_weights[i];
							}
						}
						for (int i = 0; i < dct_size; i++){
							double fi = i * scale_vert;
							int i0 = (int) fi;
							fi -= i0;
							for (int j = 0; j < dct_size; j++){
								double fj = j * scale_hor;
								int j0 = (int) fj;
								fj -= j0;
								int indx = i0*dct_size+j0;
								if ((i0 > dct_m1) || (j0 > dct_m1)){
									bidata[0][0] = 0.0;
									bidata[0][1] = 0.0;
									bidata[1][0] = 0.0;
									bidata[1][1] = 0.0;
								} else {
									bidata[0][0] = dct[indx];
									if (i0 > dct_m2) {
										bidata[1][0] = 0.0;
										bidata[1][1] = 0.0;
										if (j0 > dct_m2) {
											bidata[0][1] = 0.0;
										} else {
											bidata[0][1] = dct[indx + 1];
										}
									} else {
										bidata[1][0] = dct[indx+dct_size];
										if (j0 > dct_m2) {
											bidata[0][1] = 0.0;
											bidata[1][1] = 0.0;
										} else {
											bidata[0][1] = dct[indx + 1];
											bidata[1][1] = dct[indx + dct_size + 1];
										}
									}

								}
								// bilinear interpolation
								dct1[i*dct_size+j] =
										bidata[0][0] * (1.0-fi) * (1.0-fj) +
										bidata[0][1] * (1.0-fi) *      fj  +
										bidata[1][0] *      fi  * (1.0-fj) +
										bidata[1][1] *      fi *       fj;
								if ((globalDebugLevel > 0) && (tileY == debug_tileY) && (tileX == debug_tileX)) {
									System.out.println(i+":"+j+" {"+bidata[0][0]+","+bidata[0][1]+","+bidata[1][0]+","+bidata[1][1]+"}, ["+fi+","+fj+"] "+bidata[1][1]);
								}
							}
						}
						if (normalize) {
							double sum=0;
							for (int i = 0; i < dct_len; i++){
								sum += dct1[i] *norm_sym_weights[i];
							}
							if (sum >0.0) {
								double k = sum_orig/sum;
								for (int i = 0; i < dct_len; i++){
									dct1[i] *= k;
								}
							}
						}
//						if ((tileY == debug_tileY) && (tileX == debug_tileX) && (color == 2)) {
						if ((globalDebugLevel > 0) && (tileY == debug_tileY) && (tileX == debug_tileX)) {
							double [][] scaled_tiles = {dct, dct1};
							ShowDoubleFloatArrays sdfa_instance = new ShowDoubleFloatArrays(); // just for debugging?
							String [] titles = {"orig","scaled"};
							sdfa_instance.showArrays(scaled_tiles,  dct_size, dct_size, true, "scaled_tile", titles);
						}

						System.arraycopy(dct1, 0, dct, 0, dct_len); // replace original data
					}
				}
			};
		}
		startAndJoin(threads);
	}





	/* Create a Thread[] array as large as the number of processors available.
	 * From Stephan Preibisch's Multithreading.java class. See:
	 * http://repo.or.cz/w/trakem2.git?a=blob;f=mpi/fruitfly/general/MultiThreading.java;hb=HEAD
	 */
	public static Thread[] newThreadArray(int maxCPUs) {
		int n_cpus = Runtime.getRuntime().availableProcessors();
		if (n_cpus>maxCPUs)n_cpus=maxCPUs;
		return new Thread[n_cpus];
	}
/* Start all given threads and wait on each of them until all are done.
	 * From Stephan Preibisch's Multithreading.java class. See:
	 * http://repo.or.cz/w/trakem2.git?a=blob;f=mpi/fruitfly/general/MultiThreading.java;hb=HEAD
	 */
	public static void startAndJoin(Thread[] threads)
	{
		for (int ithread = 0; ithread < threads.length; ++ithread)
		{
			threads[ithread].setPriority(Thread.NORM_PRIORITY);
			threads[ithread].start();
		}

		try
		{
			for (int ithread = 0; ithread < threads.length; ++ithread)
				threads[ithread].join();
		} catch (InterruptedException ie)
		{
			throw new RuntimeException(ie);
		}
	}

	// temporary switch between old/new implementations

	public double [][][][][][] clt_aberrations_quad_corr(
			final ImageDttParameters  imgdtt_params,   // Now just extra correlation parameters, later will include, most others
			final int                 macro_scale,     // to correlate tile data instead of the pixel data: 1 - pixels, 8 - tiles
			final int [][]            tile_op,         // [tilesY][tilesX] - what to do - 0 - nothing for this tile
			final double [][]         disparity_array, // [tilesY][tilesX] - individual per-tile expected disparity
			final double [][][]       image_data, // first index - number of image in a quad
		    final boolean [][]        saturation_imp, // (near) saturated pixels or null
			 // correlation results - final and partial
			final double [][][][]     clt_corr_combo,  // [type][tilesY][tilesX][(2*transform_size-1)*(2*transform_size-1)] // if null - will not calculate
			                                           // [type][tilesY][tilesX] should be set by caller
													   // types: 0 - selected correlation (product+offset), 1 - sum

			final double [][][][][]   clt_corr_partial,// [tilesY][tilesX][quad]color][(2*transform_size-1)*(2*transform_size-1)] // if null - will not calculate
                                                       // [tilesY][tilesX] should be set by caller
			final double [][]         clt_mismatch,    // [12][tilesY * tilesX] // transpose unapplied. null - do not calculate
			                                           // values in the "main" directions have disparity (*_CM) subtracted, in the perpendicular - as is

			final double [][]         disparity_map,   // [8][tilesY][tilesX], only [6][] is needed on input or null - do not calculate
			                                           // last 2 - contrast, avg/ "geometric average)
			final double [][][][]     texture_tiles,   // [tilesY][tilesX]["RGBA".length()][];  null - will skip images combining

			final int                 width,
			final double              corr_fat_zero,    // add to denominator to modify phase correlation (same units as data1, data2). <0 - pure sum
			final boolean             corr_sym,
			final double              corr_offset,
			final double              corr_red,
			final double              corr_blue,
			final double              corr_sigma,
			final boolean             corr_normalize,  // normalize correlation results by rms
	  		final double              min_corr,        // 0.02; // minimal correlation value to consider valid
			final double              max_corr_sigma,  // 1.2;  // weights of points around global max to find fractional
			final double              max_corr_radius, // 3.9;

//			final int                 enhortho_width,  // 2;    // reduce weight of center correlation pixels from center (0 - none, 1 - center, 2 +/-1 from center)
//			final double              enhortho_scale,  // 0.2;  // multiply center correlation pixels (inside enhortho_width)

			final boolean 			  max_corr_double, //"Double pass when masking center of mass to reduce preference for integer values
			final int                 corr_mode, // Correlation mode: 0 - integer max, 1 - center of mass, 2 - polynomial
			final double              min_shot,        // 10.0;  // Do not adjust for shot noise if lower than
			final double              scale_shot,      // 3.0;   // scale when dividing by sqrt ( <0 - disable correction)
			final double              diff_sigma,      // 5.0;//RMS difference from average to reduce weights (~ 1.0 - 1/255 full scale image)
			final double              diff_threshold,  // 5.0;   // RMS difference from average to discard channel (~ 1.0 - 1/255 full scale image)
			final boolean             diff_gauss,      // true;  // when averaging images, use gaussian around average as weight (false - sharp all/nothing)
			final double              min_agree,       // 3.0;   // minimal number of channels to agree on a point (real number to work with fuzzy averages)
			final boolean             dust_remove,     // Do not reduce average weight when only one image differes much from the average
			final boolean             keep_weights,    // Add port weights to RGBA stack (debug feature)
			final GeometryCorrection  geometryCorrection,
			final GeometryCorrection  geometryCorrection_main, // if not null correct this camera (aux) to the coordinates of the main
			final double [][][][][][] clt_kernels, // [channel_in_quad][color][tileY][tileX][band][pixel] , size should match image (have 1 tile around)
			final int                 kernel_step,
			final int                 transform_size,
			final int                 window_type,
			final double [][]         shiftXY, // [port]{shiftX,shiftY}
			final double              disparity_corr, // disparity at infinity
			final double [][][]       fine_corr, // quadratic coefficients for fine correction (or null)
			final double              corr_magic_scale, // still not understood coefficient that reduces reported disparity value.  Seems to be around 0.85
			final double              shiftX, // shift image horizontally (positive - right) - just for testing
			final double              shiftY, // shift image vertically (positive - down)
			final int                 debug_tileX,
			final int                 debug_tileY,
			final boolean             no_fract_shift,
			final boolean             no_deconvolution,
			final int                 threadsMax,  // maximal number of threads to launch
			final int                 globalDebugLevel)
	{
		if (imgdtt_params.corr_var_cam) { // New correlation mode compatible with 8 subcameras
			return clt_aberrations_quad_corr_new(
					imgdtt_params,   // Now just extra correlation parameters, later will include, most others
					macro_scale,     // to correlate tile data instead of the pixel data: 1 - pixels, 8 - tiles
					tile_op,         // [tilesY][tilesX] - what to do - 0 - nothing for this tile
					disparity_array, // [tilesY][tilesX] - individual per-tile expected disparity
					image_data, // first index - number of image in a quad
				    saturation_imp, // (near) saturated pixels or null
					 // correlation results - final and partial
					clt_corr_combo,  // [type][tilesY][tilesX][(2*transform_size-1)*(2*transform_size-1)] // if null - will not calculate
					                                           // [type][tilesY][tilesX] should be set by caller
															   // types: 0 - selected correlation (product+offset), 1 - sum

					clt_corr_partial,// [tilesY][tilesX][quad]color][(2*transform_size-1)*(2*transform_size-1)] // if null - will not calculate
		                                                       // [tilesY][tilesX] should be set by caller
					clt_mismatch,    // [12][tilesY * tilesX] // transpose unapplied. null - do not calculate
					                                           // values in the "main" directions have disparity (*_CM) subtracted, in the perpendicular - as is

					disparity_map,   // [8][tilesY][tilesX], only [6][] is needed on input or null - do not calculate
					                                           // last 2 - contrast, avg/ "geometric average)
					texture_tiles,   // [tilesY][tilesX]["RGBA".length()][];  null - will skip images combining

					width,
					corr_fat_zero,    // add to denominator to modify phase correlation (same units as data1, data2). <0 - pure sum
					corr_sym,
					corr_offset,
					corr_red,
					corr_blue,
					corr_sigma,
					corr_normalize,  // normalize correlation results by rms
			  		min_corr,        // 0.02; // minimal correlation value to consider valid
					max_corr_sigma,  // 1.2;  // weights of points around global max to find fractional
					max_corr_radius, // 3.9;
					max_corr_double, //"Double pass when masking center of mass to reduce preference for integer values
					corr_mode, // Correlation mode: 0 - integer max, 1 - center of mass, 2 - polynomial
					min_shot,        // 10.0;  // Do not adjust for shot noise if lower than
					scale_shot,      // 3.0;   // scale when dividing by sqrt ( <0 - disable correction)
					diff_sigma,      // 5.0;//RMS difference from average to reduce weights (~ 1.0 - 1/255 full scale image)
					diff_threshold,  // 5.0;   // RMS difference from average to discard channel (~ 1.0 - 1/255 full scale image)
					diff_gauss,      // true;  // when averaging images, use gaussian around average as weight (false - sharp all/nothing)
					min_agree,       // 3.0;   // minimal number of channels to agree on a point (real number to work with fuzzy averages)
					dust_remove,     // Do not reduce average weight when only one image differes much from the average
					keep_weights,    // Add port weights to RGBA stack (debug feature)
					geometryCorrection,
					geometryCorrection_main, // final GeometryCorrection  geometryCorrection_main, // if not null correct this camera (aux) to the coordinates of the main
					clt_kernels, // [channel_in_quad][color][tileY][tileX][band][pixel] , size should match image (have 1 tile around)
					kernel_step,
					transform_size,
					window_type,
					shiftXY, // [port]{shiftX,shiftY}
					disparity_corr, // disparity at infinity
					fine_corr, // quadratic coefficients for fine correction (or null)
					corr_magic_scale, // still not understood coefficient that reduces reported disparity value.  Seems to be around 0.85
					shiftX, // shift image horizontally (positive - right) - just for testing
					shiftY, // shift image vertically (positive - down)
					debug_tileX,
					debug_tileY,
					no_fract_shift,
					no_deconvolution,
					threadsMax,  // maximal number of threads to launch
					globalDebugLevel);
		} else { // old way?
			return clt_aberrations_quad_corr_old(
					imgdtt_params,   // Now just extra correlation parameters, later will include, most others
					macro_scale,     // to correlate tile data instead of the pixel data: 1 - pixels, 8 - tiles
					tile_op,         // [tilesY][tilesX] - what to do - 0 - nothing for this tile
					disparity_array, // [tilesY][tilesX] - individual per-tile expected disparity
					image_data, // first index - number of image in a quad
				    saturation_imp, // (near) saturated pixels or null
					 // correlation results - final and partial
					clt_corr_combo,  // [type][tilesY][tilesX][(2*transform_size-1)*(2*transform_size-1)] // if null - will not calculate
					                                           // [type][tilesY][tilesX] should be set by caller
															   // types: 0 - selected correlation (product+offset), 1 - sum

					clt_corr_partial,// [tilesY][tilesX][quad]color][(2*transform_size-1)*(2*transform_size-1)] // if null - will not calculate
		                                                       // [tilesY][tilesX] should be set by caller
					clt_mismatch,    // [12][tilesY * tilesX] // transpose unapplied. null - do not calculate
					                                           // values in the "main" directions have disparity (*_CM) subtracted, in the perpendicular - as is

					disparity_map,   // [8][tilesY][tilesX], only [6][] is needed on input or null - do not calculate
					                                           // last 2 - contrast, avg/ "geometric average)
					texture_tiles,   // [tilesY][tilesX]["RGBA".length()][];  null - will skip images combining

					width,
					corr_fat_zero,    // add to denominator to modify phase correlation (same units as data1, data2). <0 - pure sum
					corr_sym,
					corr_offset,
					corr_red,
					corr_blue,
					corr_sigma,
					corr_normalize,  // normalize correlation results by rms
			  		min_corr,        // 0.02; // minimal correlation value to consider valid
					max_corr_sigma,  // 1.2;  // weights of points around global max to find fractional
					max_corr_radius, // 3.9;

//					final int                 enhortho_width,  // 2;    // reduce weight of center correlation pixels from center (0 - none, 1 - center, 2 +/-1 from center)
//					final double              enhortho_scale,  // 0.2;  // multiply center correlation pixels (inside enhortho_width)

					max_corr_double, //"Double pass when masking center of mass to reduce preference for integer values
					corr_mode, // Correlation mode: 0 - integer max, 1 - center of mass, 2 - polynomial
					min_shot,        // 10.0;  // Do not adjust for shot noise if lower than
					scale_shot,      // 3.0;   // scale when dividing by sqrt ( <0 - disable correction)
					diff_sigma,      // 5.0;//RMS difference from average to reduce weights (~ 1.0 - 1/255 full scale image)
					diff_threshold,  // 5.0;   // RMS difference from average to discard channel (~ 1.0 - 1/255 full scale image)
					diff_gauss,      // true;  // when averaging images, use gaussian around average as weight (false - sharp all/nothing)
					min_agree,       // 3.0;   // minimal number of channels to agree on a point (real number to work with fuzzy averages)
					dust_remove,     // Do not reduce average weight when only one image differes much from the average
					keep_weights,    // Add port weights to RGBA stack (debug feature)
					geometryCorrection,
					clt_kernels, // [channel_in_quad][color][tileY][tileX][band][pixel] , size should match image (have 1 tile around)
					kernel_step,
					transform_size,
					window_type,
					shiftXY, // [port]{shiftX,shiftY}
					disparity_corr, // disparity at infinity
					fine_corr, // quadratic coefficients for fine correction (or null)
					corr_magic_scale, // still not understood coefficient that reduces reported disparity value.  Seems to be around 0.85
					shiftX, // shift image horizontally (positive - right) - just for testing
					shiftY, // shift image vertically (positive - down)
					debug_tileX,
					debug_tileY,
					no_fract_shift,
					no_deconvolution,
					threadsMax,  // maximal number of threads to launch
					globalDebugLevel);
		}
	}
	public double [][][][][][] clt_aberrations_quad_corr_old(
			final ImageDttParameters  imgdtt_params,   // Now just extra correlation parameters, later will include, most others
			final int                 macro_scale,     // to correlate tile data instead of the pixel data: 1 - pixels, 8 - tiles
			final int [][]            tile_op,         // [tilesY][tilesX] - what to do - 0 - nothing for this tile
			final double [][]         disparity_array, // [tilesY][tilesX] - individual per-tile expected disparity
			final double [][][]       image_data, // first index - number of image in a quad
		    final boolean [][]        saturation_imp, // (near) saturated pixels or null
			 // correlation results - final and partial
			final double [][][][]     clt_corr_combo,  // [type][tilesY][tilesX][(2*transform_size-1)*(2*transform_size-1)] // if null - will not calculate
			                                           // [type][tilesY][tilesX] should be set by caller
													   // types: 0 - selected correlation (product+offset), 1 - sum

			final double [][][][][]   clt_corr_partial,// [tilesY][tilesX][quad]color][(2*transform_size-1)*(2*transform_size-1)] // if null - will not calculate
                                                       // [tilesY][tilesX] should be set by caller
			final double [][]         clt_mismatch,    // [12][tilesY * tilesX] // transpose unapplied. null - do not calculate
			                                           // values in the "main" directions have disparity (*_CM) subtracted, in the perpendicular - as is

			final double [][]         disparity_map,   // [8][tilesY][tilesX], only [6][] is needed on input or null - do not calculate
			                                           // last 2 - contrast, avg/ "geometric average)
			final double [][][][]     texture_tiles,   // [tilesY][tilesX]["RGBA".length()][];  null - will skip images combining

			final int                 width,
			final double              corr_fat_zero,    // add to denominator to modify phase correlation (same units as data1, data2). <0 - pure sum
			final boolean             corr_sym,
			final double              corr_offset,
			final double              corr_red,
			final double              corr_blue,
			final double              corr_sigma,
			final boolean             corr_normalize,  // normalize correlation results by rms
	  		final double              min_corr,        // 0.02; // minimal correlation value to consider valid
			final double              max_corr_sigma,  // 1.2;  // weights of points around global max to find fractional
			final double              max_corr_radius, // 3.9;

//			final int                 enhortho_width,  // 2;    // reduce weight of center correlation pixels from center (0 - none, 1 - center, 2 +/-1 from center)
//			final double              enhortho_scale,  // 0.2;  // multiply center correlation pixels (inside enhortho_width)

			final boolean 			  max_corr_double, //"Double pass when masking center of mass to reduce preference for integer values
			final int                 corr_mode, // Correlation mode: 0 - integer max, 1 - center of mass, 2 - polynomial
			final double              min_shot,        // 10.0;  // Do not adjust for shot noise if lower than
			final double              scale_shot,      // 3.0;   // scale when dividing by sqrt ( <0 - disable correction)
			final double              diff_sigma,      // 5.0;//RMS difference from average to reduce weights (~ 1.0 - 1/255 full scale image)
			final double              diff_threshold,  // 5.0;   // RMS difference from average to discard channel (~ 1.0 - 1/255 full scale image)
			final boolean             diff_gauss,      // true;  // when averaging images, use gaussian around average as weight (false - sharp all/nothing)
			final double              min_agree,       // 3.0;   // minimal number of channels to agree on a point (real number to work with fuzzy averages)
			final boolean             dust_remove,     // Do not reduce average weight when only one image differes much from the average
			final boolean             keep_weights,    // Add port weights to RGBA stack (debug feature)
			final GeometryCorrection  geometryCorrection,
			final double [][][][][][] clt_kernels, // [channel_in_quad][color][tileY][tileX][band][pixel] , size should match image (have 1 tile around)
			final int                 kernel_step,
			final int                 transform_size,
			final int                 window_type,
			final double [][]         shiftXY, // [port]{shiftX,shiftY}
			final double              disparity_corr, // disparity at infinity
			final double [][][]       fine_corr, // quadratic coefficients for fine correction (or null)
			final double              corr_magic_scale, // still not understood coefficient that reduces reported disparity value.  Seems to be around 0.85
			final double              shiftX, // shift image horizontally (positive - right) - just for testing
			final double              shiftY, // shift image vertically (positive - down)
			final int                 debug_tileX,
			final int                 debug_tileY,
			final boolean             no_fract_shift,
			final boolean             no_deconvolution,
			final int                 threadsMax,  // maximal number of threads to launch
			final int                 globalDebugLevel)
	{
//		final boolean debug_ports_coordinates = (debug_tileX == -1234);
		final double poly_corr = imgdtt_params.poly_corr_scale; // maybe add per-tile task bits to select none/near/far
		final boolean macro_mode = macro_scale != 1;      // correlate tile data instead of the pixel data
		final int quad = 4;   // number of subcameras
		final int numcol = 3; // number of colors
		final int nChn = image_data[0].length;
		final int height=image_data[0][0].length/width;
		final int tilesX=width/transform_size;
		final int tilesY=height/transform_size;
		final int nTilesInChn=tilesX*tilesY;
		final double [][][][][][] clt_data = new double[quad][nChn][tilesY][tilesX][][];
		final Thread[] threads = newThreadArray(threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		final double [] col_weights= new double [numcol]; // colors are RBG
		col_weights[2] = 1.0/(1.0 + corr_red + corr_blue);    // green color
		col_weights[0] = corr_red *  col_weights[2];
		col_weights[1] = corr_blue * col_weights[2];
		final int corr_size = transform_size * 2 -1;
		final int [][] transpose_indices = new int [corr_size*(corr_size-1)/2][2];
		int indx = 0;
		if ((globalDebugLevel > -10) && (disparity_corr != 0.0)){
			System.out.println(String.format("Using manual infinity disparity correction of %8.5f pixels",disparity_corr));
		}
		for (int i =0; i < corr_size-1; i++){
			for (int j = i+1; j < corr_size; j++){
				transpose_indices[indx  ][0] = i * corr_size + j;
				transpose_indices[indx++][1] = j * corr_size + i;
			}
		}

		// reducing weight of on-axis correlation values to enhance detection of vertical/horizontal lines
		// multiply correlation results inside the horizontal center strip  2*getEnhOrthoWidth(isAux()) - 1 wide by getEnhOrthoScale(isAux())

		final double [] enh_ortho_scale = new double [corr_size];
		for (int i = 0; i < corr_size; i++){
			if ((i < (transform_size - imgdtt_params.getEnhOrthoWidth(isAux()))) || (i > (transform_size - 2 + imgdtt_params.getEnhOrthoWidth(isAux())))) {
				enh_ortho_scale[i] = 1.0;
			} else {
				enh_ortho_scale[i] = imgdtt_params.getEnhOrthoScale(isAux());
			}
			if (i == (transform_size-1)) enh_ortho_scale[i] = 0.0 ; // hardwired 0 in the center
			enh_ortho_scale[i] *= Math.sin(Math.PI*(i+1.0)/(2*transform_size));
		}
		if (globalDebugLevel > 1){
			System.out.println("getEnhOrthoWidth(isAux())="+ imgdtt_params.getEnhOrthoWidth(isAux())+" getEnhOrthoScale(isAux())="+ imgdtt_params.getEnhOrthoScale(isAux()));
			for (int i = 0; i < corr_size; i++){
				System.out.println(" enh_ortho_scale["+i+"]="+ enh_ortho_scale[i]);

			}
		}

		// Create window  to select center correlation strip using
		// ortho_height - full width of non-zero elements
		// ortho_eff_height - effective height (ratio of the weighted column sum to the center value)
		int wcenter = transform_size - 1;
		final double [] ortho_weights = new double [corr_size]; // [15]
		for (int i = 0; i < corr_size; i++){
			if ((i >= wcenter - imgdtt_params.ortho_height/2) && (i <= wcenter + imgdtt_params.ortho_height/2)) {
				double dx = 1.0*(i-wcenter)/(imgdtt_params.ortho_height/2 + 1);
				ortho_weights[i] = 0.5*(1.0+Math.cos(Math.PI*dx))/imgdtt_params.ortho_eff_height;
			}
		}
		if (globalDebugLevel > 0){
			System.out.println("ortho_height="+ imgdtt_params.ortho_height+" ortho_eff_height="+ imgdtt_params.ortho_eff_height);
			for (int i = 0; i < corr_size; i++){
				System.out.println(" ortho_weights["+i+"]="+ ortho_weights[i]);
			}
		}




		if (globalDebugLevel > 0) {
			System.out.println("clt_aberrations_quad_corr(): width="+width+" height="+height+" transform_size="+transform_size+
					" debug_tileX="+debug_tileX+" debug_tileY="+debug_tileY+" globalDebugLevel="+globalDebugLevel);
		}
		final int [][] zi =
			{{ 0,  1,  2,  3},
			 {-1,  0, -3,  2},
			 {-2, -3,  0,  1},
			 { 3, -2, -1,  0}};
		final int [][] corr_pairs ={ // {first, second, rot} rot: 0 - as is, 1 - swap y,x
				{0,1,0},
				{2,3,0},
				{0,2,1},
				{1,3,1}};

		final double[][] port_offsets = {
				{-0.5, -0.5},
				{ 0.5, -0.5},
				{-0.5,  0.5},
				{ 0.5,  0.5}};
		final int transform_len = transform_size * transform_size;



		final double [] filter_direct= new double[transform_len];
		if (corr_sigma == 0) {
			filter_direct[0] = 1.0;
			for (int i= 1; i<filter_direct.length;i++) filter_direct[i] =0;
		} else {
			for (int i = 0; i < transform_size; i++){
				for (int j = 0; j < transform_size; j++){
					filter_direct[i*transform_size+j] = Math.exp(-(i*i+j*j)/(2*corr_sigma)); // FIXME: should be sigma*sigma !
				}
			}
		}
		// normalize
		double sum = 0;
		for (int i = 0; i < transform_size; i++){
			for (int j = 0; j < transform_size; j++){
				double d = 	filter_direct[i*transform_size+j];
				d*=Math.cos(Math.PI*i/(2*transform_size))*Math.cos(Math.PI*j/(2*transform_size));
				if (i > 0) d*= 2.0;
				if (j > 0) d*= 2.0;
				sum +=d;
			}
		}
		for (int i = 0; i<filter_direct.length; i++){
			filter_direct[i] /= sum;
		}

		DttRad2 dtt = new DttRad2(transform_size);
		final double [] filter= dtt.dttt_iiie(filter_direct);
		for (int i=0; i < filter.length;i++) filter[i] *= 2*transform_size;

		// prepare disparity maps and weights
		final int max_search_radius = (int) Math.abs(max_corr_radius); // use negative max_corr_radius for squares instead of circles?
		final int max_search_radius_poly = 1;
		if (globalDebugLevel > 0){
			System.out.println("max_corr_radius=       "+max_corr_radius);
			System.out.println("max_search_radius=     "+max_search_radius);
			System.out.println("max_search_radius_poly="+max_search_radius_poly);
			System.out.println("corr_fat_zero=         "+corr_fat_zero);
			System.out.println("disparity_array[0][0]= "+disparity_array[0][0]);


		}
		// add optional initialization of debug layers here
		if (disparity_map != null){
			for (int i = 0; i<disparity_map.length;i++){
				if ((i != OVEREXPOSED) || (saturation_imp!= null)){
					disparity_map[i] = new double [tilesY*tilesX];
				}
			}
		}
		if (clt_mismatch != null){
			for (int i = 0; i<clt_mismatch.length;i++){
				clt_mismatch[i] = new double [tilesY*tilesX]; // will use only "center of mass" centers
			}
		}
//		final double [] corr_max_weights =(((max_corr_sigma > 0) && (disparity_map != null))?
//				setMaxXYWeights(max_corr_sigma,max_search_radius): null); // here use square anyway
		final double [] corr_max_weights_poly =(((max_corr_sigma > 0) && (disparity_map != null))?
				setMaxXYWeights(max_corr_sigma,max_search_radius_poly): null); // here use square anyway

		dtt.set_window(window_type);
		final double [] lt_window = dtt.getWin2d();	// [256]
		final double [] lt_window2 = new double [lt_window.length]; // squared
		for (int i = 0; i < lt_window.length; i++) lt_window2[i] = lt_window[i] * lt_window[i];


		if (globalDebugLevel > 1) {
			ShowDoubleFloatArrays sdfa_instance = new ShowDoubleFloatArrays(); // just for debugging?
			sdfa_instance.showArrays(lt_window,  2*transform_size, 2*transform_size, "lt_window");
		}
		if (globalDebugLevel > 0) {
			System.out.println("macro_mode="+macro_mode);
		}



/*
		final double [] corr_wndy = (new Correlation2d(transform_size)).halfFlatTopWindow(
				imgdtt_params.corr_wndy_size, // int     ihwidth,
				imgdtt_params.corr_wndy_hwidth, //  double  hwidth,
				imgdtt_params.corr_wndy_blur, // double  blur,
				true, // boolean normalize,
				1.0); // double  scale);

		final double [] corr_wndx = (new Correlation2d(transform_size)).halfFlatTopWindow(
				imgdtt_params.corr_wndx_size, // int     ihwidth,
				imgdtt_params.corr_wndx_hwidth, //  double  hwidth,
				imgdtt_params.corr_wndx_blur, // double  blur,
				true, // boolean normalize,
				2.0); // double  scale);
		if (globalDebugLevel > -1) {
			System.out.println("\ncorr_wndy:");
			for (int i = 0; i <corr_wndy.length; i++) {
				System.out.println(i+": "+corr_wndy[i]);
			}
			System.out.println("\ncorr_wnxy:");
			for (int i = 0; i <corr_wndx.length; i++) {
				System.out.println(i+": "+corr_wndx[i]);
			}
		}
*/

		final Matrix [] corr_rots = geometryCorrection.getCorrVector().getRotMatrices(); // get array of per-sensor rotation matrices
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				@Override
				public void run() {
					DttRad2 dtt = new DttRad2(transform_size);
					dtt.set_window(window_type);
					int tileY,tileX,tIndex; // , chn;
					//						showDoubleFloatArrays sdfa_instance = new showDoubleFloatArrays(); // just for debugging?
					double centerX; // center of aberration-corrected (common model) tile, X
					double centerY; //
					double [][] fract_shiftsXY = new double[quad][];
					double [][]     tcorr_combo =    null; // [15*15] pixel space
					double [][][]   tcorr_partial =  null; // [quad][numcol+1][15*15]
					double [][][][] tcorr_tpartial = null; // [quad][numcol+1][4][8*8]
					PolynomialApproximation pa =     null;
					Correlation2d corr2d = new Correlation2d(
							imgdtt_params,              // ImageDttParameters  imgdtt_params,
							transform_size,             // int transform_size,
							2.0,                        //  double wndx_scale, // (wndy scale is always 1.0)
							isMonochrome(), // boolean monochrome,
							(globalDebugLevel > -1));   //   boolean debug)

					corr2d.createOrtoNotch(
							imgdtt_params.getEnhOrthoWidth(isAux()), // double enhortho_width,
							imgdtt_params.getEnhOrthoScale(isAux()), //double enhortho_scale,
							false); // true); // boolean debug);
//					public int     enhortho_width =         2;   // reduce weight of center correlation pixels from center (0 - none, 1 - center, 2 +/-1 from center)
//					public double  enhortho_scale =         0.0; // 0.2;  // multiply center correlation pixels (inside enhortho_width)

					if (corr_max_weights_poly !=null)   pa = new PolynomialApproximation(0); // debug level
					for (int nTile = ai.getAndIncrement(); nTile < nTilesInChn; nTile = ai.getAndIncrement()) {

						tileY = nTile /tilesX;
						tileX = nTile % tilesX;
						tIndex = tileY * tilesX + tileX;
						if (tile_op[tileY][tileX] == 0) continue; // nothing to do for this tile
						int                 img_mask = getImgMask(tile_op[tileY][tileX]);         // which images to use
						int                 corr_mask = getPairMask(tile_op[tileY][tileX]);       // which pairs to combine in the combo:  1 - top, 2 bottom, 4 - left, 8 - right
						// mask out pairs that use missing channels
						for (int i = 0; i< corr_pairs.length; i++){
							if ((((1 << corr_pairs[i][0]) & img_mask) == 0) || (((1 << corr_pairs[i][1]) & img_mask) == 0)) {
								corr_mask &= ~ (1 << i);
							}
						}
						boolean debugTile =(tileX == debug_tileX) && (tileY == debug_tileY);

						final int [] overexp_all = (saturation_imp != null) ? ( new int [2]): null;

						// Moved from inside chn loop
						centerX = tileX * transform_size + transform_size/2 - shiftX;
						centerY = tileY * transform_size + transform_size/2 - shiftY;
						// TODO: move port coordinates out of color channel loop
						double [][] centersXY;
						if (macro_mode){
							if ((globalDebugLevel > -1) && (tileX == debug_tileX) && (tileY == debug_tileY)) { // before correction
								System.out.println("\nUsing MACRO mode, centerX="+centerX+", centerY="+centerY);
							}
							centersXY = geometryCorrection.getPortsCoordinatesIdeal(
									macro_scale,
									centerX,
									centerY,
									macro_scale* disparity_array[tileY][tileX] + disparity_corr);

						} else {
							centersXY = geometryCorrection.getPortsCoordinatesAndDerivatives(
									geometryCorrection, //			GeometryCorrection gc_main,
									false,          // boolean use_rig_offsets,
									corr_rots, // Matrix []   rots,
									null,      //  Matrix [][] deriv_rots,
									null,      // double [][] pXYderiv, // if not null, should be double[8][]
									centerX,
									centerY,
									disparity_array[tileY][tileX] + disparity_corr);
/*
							if (dbg_ports_coords != null) {
								double [][] centersXY_shift = geometryCorrection.getPortsCoordinates_old(
										centerX,
										centerY,
										disparity_array[tileY][tileX] + disparity_corr);
								for (int dbg_p = 0; dbg_p < 4; dbg_p++){
									for (int dbg_d = 0; dbg_d < 2; dbg_d++){
										dbg_ports_coords[6 * dbg_p + 3 * dbg_d + 0][nTile] = centersXY[dbg_p][dbg_d] - centersXY_shift[dbg_p][dbg_d];
										dbg_ports_coords[6 * dbg_p + 3 * dbg_d + 1][nTile] = centersXY_shift[dbg_p][dbg_d];
										dbg_ports_coords[6 * dbg_p + 3 * dbg_d + 2][nTile] = centersXY[dbg_p][dbg_d];
									}
								}
//								centersXY = centersXY_rot; // use these
							}
*/
							if ((globalDebugLevel > 0) && (tileX == debug_tileX) && (tileY == debug_tileY)) {
								for (int i = 0; i < quad; i++) {
									System.out.println("clt_aberrations_quad_corr():  tileX="+tileX+", tileY="+tileY+
											" centerX="+centerX+" centerY="+centerY+" disparity="+disparity_array[tileY][tileX]+
											" centersXY["+i+"][0]="+centersXY[i][0]+" centersXY["+i+"][1]="+centersXY[i][1]);
								}
							}

							if ((globalDebugLevel > -1) && (tileX == debug_tileX) && (tileY == debug_tileY)) { // before correction
								System.out.print(disparity_array[tileY][tileX]+"\t"+
										centersXY[0][0]+"\t"+centersXY[0][1]+"\t"+
										centersXY[1][0]+"\t"+centersXY[1][1]+"\t"+
										centersXY[2][0]+"\t"+centersXY[2][1]+"\t"+
										centersXY[3][0]+"\t"+centersXY[3][1]+"\t");
							}

							for (int ip = 0; ip < centersXY.length; ip++){
								centersXY[ip][0] -= shiftXY[ip][0];
								centersXY[ip][1] -= shiftXY[ip][1];
							}
							// TODO: use correction after disparity applied (to work for large disparity values)
							if (fine_corr != null){

								for (int ip = 0; ip < centersXY.length; ip++){
									double [] tXY = geometryCorrection.getRelativeCoords(centersXY[ip]);
									for (int d = 0; d <2; d++) {
										centersXY[ip][d] -= (
												fine_corr[ip][d][0]*tXY[0]*tXY[0]+
												fine_corr[ip][d][1]*tXY[1]*tXY[1]+
												fine_corr[ip][d][2]*tXY[0]*tXY[1]+
												fine_corr[ip][d][3]*tXY[0]+
												fine_corr[ip][d][4]*tXY[1]+
												fine_corr[ip][d][5]);
									}
								}
							}
						} // if (macro_mode) ... else
						if (FPGA_COMPARE_DATA && (globalDebugLevel > 0) && (tileX == debug_tileX) && (tileY == debug_tileY)) {
							final int fpga_cam = 0;
							double [][] manual_offsets={
//									{ 1.3, -2.7},
//									{-1.3,  2.7},
//									{ 0.0,  0.0}};

//							{ 2.3, -2.7},
//							{-0.3,  2.7},
//							{ 0.0,  0.0}};

							{ 2.3, -2.7},
							{-0.3,  2.7},
							{ 1.0,  0.0}};

							double [][] colorCentersXY = {
									{centersXY[fpga_cam][0] + manual_offsets[0][0], centersXY[fpga_cam][1] + manual_offsets[0][1]}, // add manual offsets here
									{centersXY[fpga_cam][0] + manual_offsets[1][0], centersXY[fpga_cam][1] + manual_offsets[1][1]},
									{centersXY[fpga_cam][0] + manual_offsets[2][0], centersXY[fpga_cam][1] + manual_offsets[2][1]}
							};
							generateFPGACompareData(
									image_data[fpga_cam], // final double [][]                  image_data, // for selected subcamera
									colorCentersXY,       // final double [][]         colorCentersXY, // pixel centers per color (2 - green)
									transform_size,       // final int                 transform_size,
									width,                 // final int                 width
									dtt
									);
						}

						for (int chn = 0; chn <numcol; chn++) {
							boolean debug_for_fpga = FPGA_COMPARE_DATA && (globalDebugLevel > 0) && (tileX == debug_tileX) && (tileY == debug_tileY) && (chn == 2);
							if ((globalDebugLevel > -1) && (tileX == debug_tileX) && (tileY == debug_tileY) && (chn == 2)) {
								System.out.println("\nUsing "+(macro_mode?"MACRO":"PIXEL")+" mode, centerX="+centerX+", centerY="+centerY);
								System.out.println(disparity_array[tileY][tileX]+"\t"+
							    centersXY[0][0]+"\t"+centersXY[0][1]+"\t"+
							    centersXY[1][0]+"\t"+centersXY[1][1]+"\t"+
							    centersXY[2][0]+"\t"+centersXY[2][1]+"\t"+
							    centersXY[3][0]+"\t"+centersXY[3][1]+"\t");
							}

							for (int i = 0; i < quad; i++) {
								if (debug_for_fpga && (i==0)){
									double [][] fpga_clt_data = new double [4][];
									double [] fpga_fract_shiftsXY;
									double [] fpga_centersXY = {centersXY[i][0],centersXY[i][1]};
									int fpga_chn = chn; // ==2, green
									// round to nearest 1/128 pix (supported by FPGA code)
									System.out.println(String.format("Center X= %f, center Y = %f", fpga_centersXY[0],fpga_centersXY[1]));
									for (int j=0; j<2;j++){
										fpga_centersXY[j] = Math.round(128*fpga_centersXY[j])/128.0;
									}


									for (int j=0; j<2;j++){
										fpga_centersXY[j] = Math.round(fpga_centersXY[j]);
									}


//									fpga_centersXY[0]+=0.5; // half pixel shift horizontal zero pixel shift vertical
//									fpga_centersXY[1]+=0.5; // half pixel shift vertical, zero pixel shift horizontal

//									fpga_centersXY[0]+=1.0; //

									fpga_chn =           2;


									System.out.println(String.format("Manually changing offset: center X= %f, center Y = %f", fpga_centersXY[0],fpga_centersXY[1]));
									System.out.println(String.format("Manually changing color to %d (was %d)", fpga_chn, chn));



									fpga_fract_shiftsXY = extract_correct_tile( // return a pair of residual offsets
											image_data[i],
											width,       // image width
											null,
											fpga_clt_data, //double  [][]        clt_tile,    // should be double [4][];
											kernel_step,
											transform_size,
											dtt,
											fpga_chn, // chn,
											fpga_centersXY[0], // centersXY[i][0], // centerX, // center of aberration-corrected (common model) tile, X
											fpga_centersXY[1], // centersXY[i][1], // centerY, //
											-10, // globalDebugLevel,
											true, // no_deconvolution,
											false, // ); // transpose);
											null,
											null);
									ShowDoubleFloatArrays sdfa_instance = new ShowDoubleFloatArrays(); // just for debugging?
									String [] titles = {"CC","SC","CS","SS"};
									double [][] dbg_tile = new double [4][];
									for (int im = 0; im < 4; im++) dbg_tile[im]=fpga_clt_data[im];
									sdfa_instance.showArrays(dbg_tile,  transform_size, transform_size, true, "fpre-shifted_x"+tileX+"_y"+tileY+"-z", titles);
									fract_shift(    // fractional shift in transform domain. Currently uses sin/cos - change to tables with 2? rotations
											fpga_clt_data, // double  [][]  clt_tile,
											transform_size,
											fpga_fract_shiftsXY[0],            // double        shiftX,
											fpga_fract_shiftsXY[1],            // double        shiftY,
											true); // debug
									for (int im = 0; im < 4; im++) dbg_tile[im]=fpga_clt_data[im];
									sdfa_instance.showArrays(dbg_tile,  transform_size, transform_size, true, "f-shifted_x"+tileX+"_y"+tileY+"-z", titles);
									System.out.println("Debugging for FPGA data, globalDebugLevel = "+globalDebugLevel+", tileX="+tileX+", tileY="+tileY+", sesnlor="+i+", color="+chn);
									System.out.println("Debugging for FPGA data, fpga_fract_shiftsXY[0] = "+fpga_fract_shiftsXY[0]+", fpga_fract_shiftsXY[1]="+fpga_fract_shiftsXY[1]);
									System.out.println();

									double scale = (1 << (FPGA_DTT_IN - 9)); //  -1;
									// compensate for DTT scale
									scale *= 1.0 *((1 << FPGA_WND_BITS) -1) / (1 << FPGA_WND_BITS);
									scale *= 1.0 *((1 << FPGA_WND_BITS) -1) / (1 << FPGA_WND_BITS);
									// compensate for rotator scale:
									scale *= 1.0 *((1 << FPGA_WND_BITS) -1) / (1 << FPGA_WND_BITS);
									scale *= 1.0 *((1 << FPGA_WND_BITS) -1) / (1 << FPGA_WND_BITS);
									double [] fpga_dtt_lim = {0.0,0.0};
									for (int dct_mode = 0; dct_mode <4; dct_mode++) {
										for (int j = 0; j < 64; j++){
											if (fpga_clt_data[dct_mode][j] > fpga_dtt_lim[0]) fpga_dtt_lim[0] = fpga_clt_data[dct_mode][j];
											if (fpga_clt_data[dct_mode][j] < fpga_dtt_lim[1]) fpga_dtt_lim[1] = fpga_clt_data[dct_mode][j];
										}
									}

									System.out.println(String.format("// DTT rotated, shift_x=%f. shift_y = %f", fpga_fract_shiftsXY[0],fpga_fract_shiftsXY[1]));
									System.out.println(String.format("// DTT rotated  range: %f ... %f", fpga_dtt_lim[1], fpga_dtt_lim[0]));
//									scale = (1 << (FPGA_DTT_IN - 9)); //  -1;
									for (int dct_mode = 0; dct_mode <4; dct_mode++) {
										for (int j = 0; j < 64; j++){
											int id = (int) Math.round(scale * fpga_clt_data[dct_mode][j]);
											System.out.print(String.format("%7x ", id & ((1 << 25) -1)));
											if ((j % 8) == 7) System.out.println();
										}
										System.out.println();
									}
								} // end of debug_for_fpga

								clt_data[i][chn][tileY][tileX] = new double [4][];
								fract_shiftsXY[i] = extract_correct_tile( // return a pair of residual offsets
										image_data[i],
										width,       // image width
										(clt_kernels == null) ? null : clt_kernels[i], // [color][tileY][tileX][band][pixel]
										clt_data[i][chn][tileY][tileX], //double  [][]        clt_tile,    // should be double [4][];
										kernel_step,
										transform_size,
										dtt,
										chn,
										centersXY[i][0], // centerX, // center of aberration-corrected (common model) tile, X
										centersXY[i][1], // centerY, //
//										((globalDebugLevel > -1) && (tileX == debug_tileX) && (tileY == debug_tileY) && (chn == 2)),
										(!FPGA_COMPARE_DATA && (globalDebugLevel > -1) && (tileX == debug_tileX) && (tileY == debug_tileY) && (chn == 2) && (i==0)) ? (globalDebugLevel + 0) : 0, // external tile compare

//										(globalDebugLevel > 0) && (tileX == debug_tileX) && (tileY == debug_tileY) && (chn == 2), // external tile compare
										no_deconvolution,
										false, // ); // transpose);
										((saturation_imp != null) ? saturation_imp[i] : null), //final boolean [][]        saturation_imp, // (near) saturated pixels or null
										((saturation_imp != null) ? overexp_all: null)); // final double [] overexposed)


							}
							if ((globalDebugLevel > -1) && (tileX == debug_tileX) && (tileY == debug_tileY) && (chn == 2)) {
								System.out.println();
							}
							if ((globalDebugLevel > 0) && (debug_tileX == tileX) && (debug_tileY == tileY)  && (chn == 2) && !FPGA_COMPARE_DATA) {
								ShowDoubleFloatArrays sdfa_instance = new ShowDoubleFloatArrays(); // just for debugging?
								String [] titles = {"CC0","SC0","CS0","SS0","CC1","SC1","CS1","SS1","CC2","SC2","CS2","SS2","CC3","SC3","CS3","SS3"};
								double [][] dbg_tile = new double [16][];
								for (int i = 0; i < 16; i++) dbg_tile[i]=clt_data[i>>2][chn][tileY][tileX][i & 3];
								sdfa_instance.showArrays(dbg_tile,  transform_size, transform_size, true, "pre-shifted_x"+tileX+"_y"+tileY, titles);
							}

							if ((globalDebugLevel > 0) && (tileX >= debug_tileX - 2) && (tileX <= debug_tileX + 2) &&
									(tileY >= debug_tileY - 2) && (tileY <= debug_tileY+2)) {
								for (int i = 0; i < quad; i++) {
									System.out.println("clt_aberrations_quad(): color="+chn+", tileX="+tileX+", tileY="+tileY+
											" fract_shiftsXY["+i+"][0]="+fract_shiftsXY[i][0]+" fract_shiftsXY["+i+"][1]="+fract_shiftsXY[i][1]);
								}
							}

//							if (!no_fract_shift && !FPGA_COMPARE_DATA) {
							if (!no_fract_shift) {
								// apply residual shift
								for (int i = 0; i < quad; i++) {
									fract_shift(    // fractional shift in transform domain. Currently uses sin/cos - change to tables with 2? rotations
											clt_data[i][chn][tileY][tileX], // double  [][]  clt_tile,
											transform_size,
											fract_shiftsXY[i][0],            // double        shiftX,
											fract_shiftsXY[i][1],            // double        shiftY,
											//									(globalDebugLevel > 0) && (tileX == debug_tileX) && (tileY == debug_tileY)); // external tile compare
											((globalDebugLevel > 1) && (chn==0) && (tileX >= debug_tileX - 2) && (tileX <= debug_tileX + 2) &&
													(tileY >= debug_tileY - 2) && (tileY <= debug_tileY+2)));
								}
								if ((globalDebugLevel > 0) && (debug_tileX == tileX) && (debug_tileY == tileY)) {
									ShowDoubleFloatArrays sdfa_instance = new ShowDoubleFloatArrays(); // just for debugging?
									String [] titles = {"CC0","SC0","CS0","SS0","CC1","SC1","CS1","SS1","CC2","SC2","CS2","SS2","CC3","SC3","CS3","SS3"};
									double [][] dbg_tile = new double [16][];
									for (int i = 0; i < 16; i++) dbg_tile[i]=clt_data[i>>2][chn][tileY][tileX][i & 3];
									sdfa_instance.showArrays(dbg_tile,  transform_size, transform_size, true, "shifted_x"+tileX+"_y"+tileY+"-z", titles);
								}



							}
						}


						// calculate overexposed fraction
						if (saturation_imp != null) {
							disparity_map[OVEREXPOSED][nTile] = (1.0 * overexp_all[0]) / overexp_all[1];
						}

						// all color channels are done here
						double extra_disparity = 0.0; // if allowed, shift images extra before trying to combine
						if (clt_corr_combo != null){ // not null - calculate correlations

							tcorr_tpartial=new double[corr_pairs.length][numcol+1][4][transform_len];
							tcorr_partial =  new double[quad][numcol+1][];

							for (int pair = 0; pair < corr_pairs.length; pair++){
								for (int chn = 0; chn <numcol; chn++){
									double [][] data1 = clt_data[corr_pairs[pair][0]][chn][tileY][tileX];
									double [][] data2 = clt_data[corr_pairs[pair][1]][chn][tileY][tileX];
									/* for (int i = 0; i < transform_len; i++) {
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
								} */


								double [] a2 = new double[transform_len];
						    	double sa2 = 0.0;
								for (int i = 0; i < transform_len; i++) {
									double s1 = 0.0, s2=0.0;
									for (int n = 0; n< 4; n++){
										s1+=data1[n][i] * data1[n][i];
										s2+=data2[n][i] * data2[n][i];
									}
									a2[i] = Math.sqrt(s1*s2);
									sa2 += a2[i];
								}
								double fz2 = sa2/transform_len * corr_fat_zero * corr_fat_zero; // fat_zero squared to match units
								for (int i = 0; i < transform_len; i++) {
									double scale = 1.0 / (a2[i] + fz2);
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
//									tcorr_partial[pair][chn] = corr_unfold_tile(
									tcorr_partial[pair][chn] = dtt.corr_unfold_tile(
											tcorr_tpartial[pair][chn],
											transform_size);
								}
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
							} // all pairs calculated
							tcorr_combo = new double [TCORR_TITLES.length][corr_size * corr_size];

							int numPairs = 	0, numPairsHor = 0, numPairsVert = 0;
							for (int pair = 0; pair < corr_pairs.length; pair++) if (((corr_mask >> pair) & 1) != 0){
								numPairs++;
								if (corr_pairs[pair][2] == 0) { // horizontal pair)
									numPairsHor++;
								} else {
									numPairsVert++;
								}
							}
							double avScale = 0.0, avScaleHor = 0.0, avScaleVert = 0.0;
							if ((globalDebugLevel > -1) && (tileX == debug_tileX) && (tileY == debug_tileY)){
								System.out.println ("Before combining tiles, numcol="+numcol);
								for (int pair = 0; pair < corr_pairs.length; pair++) if (((corr_mask >> pair) & 1) != 0){
									System.out.println("pair # "+pair);
									for (int i = 0; i < corr_size; i++) {
										System.out.print(String.format("%2d:", i));
										for (int j = 0; j < corr_size; j++) {
											System.out.print(String.format(" %8.5f", tcorr_partial[pair][numcol][i * corr_size + j]));
										}
										System.out.println();
									}
									System.out.println();
								}
							}
//		final int corr_size = transform_size * 2 -1;
							if (numPairs > 0) {
								boolean debugMax = (globalDebugLevel > 1) && (tileX == debug_tileX) && (tileY == debug_tileY);
								avScale = 1.0/numPairs;
								if (numPairsHor > 0)  avScaleHor = 1.0/numPairsHor;
								if (numPairsVert > 0) avScaleVert = 1.0/numPairsVert;
								if (debugMax) {
									System.out.println("avScale = "+avScale+", avScaleHor = "+avScaleHor+", avScaleVert = "+avScaleVert+", corr_offset = "+corr_offset);
								}
								if (corr_offset < 0) { // just add all partial correlations for composite color
									for (int i = 0; i < tcorr_combo[TCORR_COMBO_RSLT].length; i++){
										tcorr_combo[TCORR_COMBO_RSLT][i] = 0.0;
										tcorr_combo[TCORR_COMBO_HOR][i] = 0.0;
										tcorr_combo[TCORR_COMBO_VERT][i] = 0.0;
										for (int pair = 0; pair < corr_pairs.length; pair++) if (((corr_mask >> pair) & 1) != 0){
											tcorr_combo[TCORR_COMBO_RSLT][i] += avScale*tcorr_partial[pair][numcol][i]; // only composite color channel
											if (corr_pairs[pair][2] == 0) { // horizontal pair
												tcorr_combo[TCORR_COMBO_HOR][i] +=  avScaleHor*tcorr_partial[pair][numcol][i]; // only composite color channel
											} else { //vertical pair
												tcorr_combo[TCORR_COMBO_VERT][i] += avScaleVert*tcorr_partial[pair][numcol][i]; // only composite color channel
											}
											if (debugMax) {
												System.out.println("tcorr_combo[TCORR_COMBO_RSLT]["+i+"]="+tcorr_combo[TCORR_COMBO_RSLT][i]+" tcorr_partial["+pair+"]["+numcol+"]["+i+"]="+tcorr_partial[pair][numcol][i]);
											}
										}
//										tcorr_combo[TCORR_COMBO_HOR][i] *=2;   // no have the same scale as tcorr_combo[TCORR_COMBO_RSLT]
//										tcorr_combo[TCORR_COMBO_VERT][i] *=2;
									}
								} else {
									for (int i = 0; i < tcorr_combo[TCORR_COMBO_RSLT].length; i++){
										tcorr_combo[TCORR_COMBO_RSLT][i] = 1.0;
										tcorr_combo[TCORR_COMBO_HOR][i] =  1.0;
										tcorr_combo[TCORR_COMBO_VERT][i] = 1.0;
										for (int pair = 0; pair < corr_pairs.length; pair++) if (((corr_mask >> pair) & 1) != 0){
											tcorr_combo[TCORR_COMBO_RSLT][i] *= (tcorr_partial[pair][numcol][i] + corr_offset); // only composite color channel
											if (corr_pairs[pair][2] == 0) { // horizontal pair
												tcorr_combo[TCORR_COMBO_HOR][i] *= (tcorr_partial[pair][numcol][i] + corr_offset); // only composite color channel
											} else { //vertical pair
												tcorr_combo[TCORR_COMBO_VERT][i] *= (tcorr_partial[pair][numcol][i] + corr_offset); // only composite color channel
											}
											if (debugMax) {
												System.out.println("tcorr_combo[TCORR_COMBO_RSLT]["+i+"]="+tcorr_combo[TCORR_COMBO_RSLT][i]+" tcorr_partial["+pair+"]["+numcol+"]["+i+"]="+tcorr_partial[pair][numcol][i]);
											}
										}
//										tcorr_combo[TCORR_COMBO_HOR][i] *= tcorr_combo[TCORR_COMBO_HOR][i];   // no have the same scale as tcorr_combo[TCORR_COMBO_RSLT]
//										tcorr_combo[TCORR_COMBO_VERT][i] *= tcorr_combo[TCORR_COMBO_VERT][i];
										if (corr_normalize) {
											if (tcorr_combo[TCORR_COMBO_RSLT][i] > 0.0){
												tcorr_combo[TCORR_COMBO_RSLT][i] = Math.pow(tcorr_combo[TCORR_COMBO_RSLT][i],avScale) - corr_offset;
											} else {
												tcorr_combo[TCORR_COMBO_RSLT][i] =  -corr_offset;
											}

											if (tcorr_combo[TCORR_COMBO_HOR][i] > 0.0){
												tcorr_combo[TCORR_COMBO_HOR][i] = Math.pow(tcorr_combo[TCORR_COMBO_HOR][i],avScaleHor) - corr_offset;
											} else {
												tcorr_combo[TCORR_COMBO_HOR][i] =  -corr_offset;
											}

											if (tcorr_combo[TCORR_COMBO_VERT][i] > 0.0){
												tcorr_combo[TCORR_COMBO_VERT][i] = Math.pow(tcorr_combo[TCORR_COMBO_VERT][i],avScaleVert) - corr_offset;
											} else {
												tcorr_combo[TCORR_COMBO_VERT][i] =  -corr_offset;
											}
										}
									}
								}
								if ((globalDebugLevel > -1) && (tileX == debug_tileX) && (tileY == debug_tileY)){
									System.out.println ("After combining tiles, horizontal");
									for (int i = 0; i < corr_size; i++) {
										System.out.print(String.format("%2d:", i));
										for (int j = 0; j < corr_size; j++) {
											System.out.print(String.format(" %8.5f", tcorr_combo[TCORR_COMBO_HOR][i * corr_size + j]));
										}
										System.out.println();
									}
									System.out.println();
								}



								// calculate sum also
								for (int i = 0; i < tcorr_combo[TCORR_COMBO_SUM].length; i++){
									tcorr_combo[TCORR_COMBO_SUM][i] = 0.0;
									for (int pair = 0; pair < corr_pairs.length; pair++) if (((corr_mask >> pair) & 1) != 0){
										tcorr_combo[TCORR_COMBO_SUM][i] += avScale*tcorr_partial[pair][numcol][i]; // only composite color channel
										if (debugMax) {
											System.out.println("tcorr_combo[TCORR_COMBO_SUM]["+i+"]="+tcorr_combo[TCORR_COMBO_SUM][i]+" tcorr_partial["+pair+"]["+numcol+"]["+i+"]="+tcorr_partial[pair][numcol][i]);
										}
									}
								}
/*
								double [] rms = new double [tcorr_combo.length];
								for (int n = 0; n < rms.length; n++) rms[n] = 1.0;
								if (corr_normalize){ // normalize both composite and sum by their RMS
									for (int n = 0; n<tcorr_combo.length; n++){
										rms[n] = 0;
										for (int i = 0; i < tcorr_combo[n].length; i++){
											rms[n] += tcorr_combo[n][i] * tcorr_combo[n][i];
										}
										rms[n] = Math.sqrt(rms[n]/tcorr_combo[n].length);
										if (rms[n] > 0){
											double k = 1.0/rms[n];
											for (int i = 0; i < tcorr_combo[n].length; i++){
												tcorr_combo[n][i] *= k;
											}
										}
									}
								}
*/
								// return results
								for (int n = 0; n < clt_corr_combo.length; n++){ // tcorr_combo now may be longer than clt_corr_combo
									clt_corr_combo[n][tileY][tileX] = tcorr_combo[n];
								}
								if (clt_corr_partial != null){
									clt_corr_partial[tileY][tileX] = tcorr_partial;
								}
								if (disparity_map != null) {
									int [] icorr_max =getMaxXYInt( // find integer pair or null if below threshold
											tcorr_combo[TCORR_COMBO_RSLT],      // [data_size * data_size]
											corr_size,
											min_corr,    // minimal value to consider (at integer location, not interpolated)
											debugMax);
									int max_index = -1;
									if (icorr_max == null){
										disparity_map[DISPARITY_INDEX_INT]          [tIndex] = Double.NaN;
										disparity_map[DISPARITY_INDEX_INT+1]        [tIndex] = Double.NaN;
										disparity_map[DISPARITY_INDEX_CM]           [tIndex] = Double.NaN;
										disparity_map[DISPARITY_INDEX_CM+1]         [tIndex] = Double.NaN;
										disparity_map[DISPARITY_INDEX_HOR]          [tIndex] = Double.NaN;
										disparity_map[DISPARITY_INDEX_HOR_STRENGTH][tIndex] = Double.NaN;
										disparity_map[DISPARITY_INDEX_VERT]         [tIndex] = Double.NaN;
										disparity_map[DISPARITY_INDEX_VERT_STRENGTH][tIndex] = Double.NaN;
										disparity_map[DISPARITY_INDEX_POLY]         [tIndex] = Double.NaN;
										disparity_map[DISPARITY_INDEX_POLY+1]       [tIndex] = Double.NaN;
										if (clt_mismatch != null){
											for (int pair = 0; pair < corr_pairs.length; pair++) if (((corr_mask >> pair) & 1) != 0){
												clt_mismatch[3*pair + 0 ][tIndex] = Double.NaN;
												clt_mismatch[3*pair + 1 ][tIndex] = Double.NaN;
												clt_mismatch[3*pair + 2 ][tIndex] = Double.NaN;
											}
										}

									} else {
										double [] corr_max_XYi = {icorr_max[0],icorr_max[1]};
										disparity_map[DISPARITY_INDEX_INT][tIndex] =  transform_size - 1 -corr_max_XYi[0];
										disparity_map[DISPARITY_INDEX_INT+1][tIndex] = transform_size - 1 -corr_max_XYi[1];
										// for the integer maximum provide contrast and variety
										max_index = icorr_max[1]*corr_size + icorr_max[0];
										disparity_map[DISPARITY_STRENGTH_INDEX][tIndex] = tcorr_combo[TCORR_COMBO_RSLT][max_index]; // correlation combo value at the integer maximum
										// undo scaling caused by optional normalization
//										disparity_map[DISPARITY_VARIATIONS_INDEX][tIndex] = (rms[1]*tcorr_combo[1][max_index])/(rms[0]*tcorr_combo[0][max_index]); // correlation combo value at the integer maximum
										disparity_map[DISPARITY_VARIATIONS_INDEX][tIndex] = (tcorr_combo[TCORR_COMBO_SUM][max_index])/(tcorr_combo[TCORR_COMBO_RSLT][max_index]); // correlation combo value at the integer maximum
										//									Calculate "center of mass" coordinates
										double [] corr_max_XYm = getMaxXYCm( // get fractional center as a "center of mass" inside circle/square from the integer max
												tcorr_combo[TCORR_COMBO_RSLT],      // [data_size * data_size]
												corr_size,
												icorr_max, // integer center coordinates (relative to top left)
												max_corr_radius,  // positive - within that distance, negative - within 2*(-radius)+1 square
												max_corr_double, //"Double pass when masking center of mass to reduce preference for integer values
												debugMax);
										disparity_map[DISPARITY_INDEX_CM][tIndex] = transform_size - 1 -corr_max_XYm[0];
										disparity_map[DISPARITY_INDEX_CM+1][tIndex] = transform_size - 1 -corr_max_XYm[1];
										// returns x and strength, not x,y
										if ((globalDebugLevel > -1) && (tileX == debug_tileX) && (tileY == debug_tileY)){
											System.out.println ("Before combining getMaxXSOrtho() TCORR_COMBO_HOR)");
										}
										double [] corr_max_XS_hor = getMaxXSOrtho( // get fractional center as a "center of mass" inside circle/square from the integer max
												tcorr_combo[TCORR_COMBO_HOR],      // [data_size * data_size]
												enh_ortho_scale, // [data_size]
												corr_size,
//												max_corr_radius,
//												max_corr_double, // reusing, true - just poly for maximum
												(globalDebugLevel > -1) && (tileX == debug_tileX) && (tileY == debug_tileY)); // debugMax);
										disparity_map[DISPARITY_INDEX_HOR][tIndex] = transform_size - 1 - corr_max_XS_hor[0];
										disparity_map[DISPARITY_INDEX_HOR_STRENGTH][tIndex] = corr_max_XS_hor[1];
										double [] corr_max_XS_vert = getMaxXSOrtho( // get fractional center as a "center of mass" inside circle/square from the integer max
												tcorr_combo[TCORR_COMBO_VERT],      // [data_size * data_size]
												enh_ortho_scale, // [data_size]
												corr_size,
//												max_corr_radius,
//												max_corr_double, // reusing, true - just poly for maximum (probably keep it that way)
												(globalDebugLevel > -1) && (tileX == debug_tileX) && (tileY == debug_tileY)); // debugMax);
										disparity_map[DISPARITY_INDEX_VERT][tIndex] = transform_size - 1 - corr_max_XS_vert[0];
										disparity_map[DISPARITY_INDEX_VERT_STRENGTH][tIndex] = corr_max_XS_vert[1];
										//									Calculate polynomial interpolated maximum coordinates
										double [] corr_max_XY = getMaxXYPoly( // get interpolated maximum coordinates using 2-nd degree polynomial
												pa,
												tcorr_combo[TCORR_COMBO_RSLT],      // [data_size * data_size]
												corr_size,
												icorr_max,                          // integer center coordinates (relative to top left)
												corr_max_weights_poly,              // [(radius+1) * (radius+1)]
												max_search_radius_poly,             // max_search_radius, for polynomial - always use 1
												imgdtt_params.poly_pwr,             // double    value_pwr, // raise value to this power (trying to compensate sticking to integer values)
												imgdtt_params.poly_vasw_pwr,        // double   poly_vasw_pwr, // multiply weight by value

												debugMax);
										// Possible to bypass ortho calculation for bad tiles?
										double [][] corr_max_ortho = new double[2][];
										// return argmax, max, half-width
										corr_max_ortho[0] =	getMaxXSOrtho2(   // get fractional center using a quadratic polynomial
												tcorr_combo[TCORR_COMBO_HOR], // double [] data,            // [data_size * data_size]
												ortho_weights,                // double [] enhortho_scales, // [data_size]
												corr_size,                    // int       data_size,
												imgdtt_params.ortho_nsamples, // int       num_samples,     // number of samples to keep (5?)
												imgdtt_params.ortho_vasw_pwr, // double   value_as_weight, // use positive value as sample weight
												false); //(globalDebugLevel > -2) && (tileX == debug_tileX) && (tileY == debug_tileY)); // debugMax);
										corr_max_ortho[1] =	getMaxXSOrtho2(   // get fractional center using a quadratic polynomial
												tcorr_combo[TCORR_COMBO_VERT],// double [] data,            // [data_size * data_size]
												ortho_weights,                // double [] enhortho_scales, // [data_size]
												corr_size,                    // int       data_size,
												imgdtt_params.ortho_nsamples, // int       num_samples,     // number of samples to keep (5?)
												imgdtt_params.ortho_vasw_pwr, // double   value_as_weight, // use positive value as sample weight
												false); //(globalDebugLevel > -2) && (tileX == debug_tileX) && (tileY == debug_tileY)); // debugMax);
//										(globalDebugLevel > 0) && (tileX == debug_tileX) && (tileY == debug_tileY)); // debugMax);
 //										disparity_map[DISPARITY_INDEX_HOR][tIndex] = transform_size - 1 - corr_max_XS_hor[0];
//										disparity_map[DISPARITY_INDEX_HOR_STRENGTH][tIndex] = corr_max_XS_hor[1];




										//double
										if (corr_max_XY != null){
											corr_max_XY[0] = transform_size - 1 -corr_max_XY[0];
											corr_max_XY[1] = transform_size - 1 -corr_max_XY[1];
										} else {
											corr_max_XY = new double[2];
											corr_max_XY[0] = Double.NaN;
											corr_max_XY[1] = Double.NaN;
										}

										if (corr_max_ortho[0] != null){
											corr_max_ortho[0][0] = transform_size - 1 -corr_max_ortho[0][0];
										} else {
											corr_max_ortho[0] = new double[3];
											corr_max_ortho[0][0] = Double.NaN;
											corr_max_ortho[0][1] = Double.NaN;
											corr_max_ortho[0][2] = Double.NaN;
										}

										if (corr_max_ortho[1] != null){
											corr_max_ortho[1][0] = transform_size - 1 -corr_max_ortho[1][0];
										} else {
											corr_max_ortho[1] = new double[3];
											corr_max_ortho[1][0] = Double.NaN;
											corr_max_ortho[1][1] = Double.NaN;
											corr_max_ortho[1][2] = Double.NaN;
										}



										disparity_map[DISPARITY_INDEX_POLY][tIndex] =   corr_max_XY[0];
										disparity_map[DISPARITY_INDEX_POLY+1][tIndex] = corr_max_XY[1];

										// just debug (up to 4 layers)
										// use poly only if half-width y < limit (now it is ~2.0 for good corr)


										if (imgdtt_params.mix_corr_poly) { // regardless of debug
											// apply
											if ((corr_max_XY.length > 4) && !Double.isNaN(corr_max_XY[3]) && !Double.isNaN(corr_max_XY[4]) &&
													(corr_max_XY[3] < imgdtt_params.max_poly_hwidth) &&
													(corr_max_XY[4] < imgdtt_params.max_poly_hwidth) &&
													(disparity_map[DISPARITY_STRENGTH_INDEX][tIndex] > imgdtt_params.min_poly_strength)) { // debug threshold
												// now CM is actually poly!
												double pcorr = (corr_max_XY[3]-corr_max_XY[4]);
												if (pcorr > 0.0) {
													pcorr *= poly_corr; // maybe add per-tile control far/near/none
												} else {
													pcorr = 0.0;
												}

												disparity_map[DISPARITY_INDEX_CM+0][tIndex] = corr_max_XY[0]+pcorr; // only for X
												disparity_map[DISPARITY_INDEX_CM+1][tIndex] = corr_max_XY[1];

												// correct far objects by using hor/vert correlations
												if (imgdtt_params.fo_correct) {
													// check strengths:
													if ((corr_max_XY[2] >= imgdtt_params.fo_min_strength) &&
														(corr_max_ortho[0][1] >= imgdtt_params.fo_min_strength) &&
														(corr_max_ortho[1][1] >= imgdtt_params.fo_min_strength) &&
														// check halw-widths
														(corr_max_XY[3] <=       imgdtt_params.fo_max_hwidth) &&
														(corr_max_ortho[0][2] <= imgdtt_params.fo_max_hwidth) &&
														(corr_max_ortho[1][2] <= imgdtt_params.fo_max_hwidth)) {
														double disp_full = corr_max_XY[0];
														double disp_near =  Math.max(corr_max_ortho[0][0],corr_max_ortho[1][0]);
														double disp_far =   Math.min(corr_max_ortho[0][0],corr_max_ortho[1][0]);
														if ((disp_full < disp_near) && (disp_full > disp_far)) {
															double corr = (disp_near - disp_far) * imgdtt_params.fo_overcorrection + (disp_near -disp_full);
															double lim_corr = (disp_near -disp_full) * imgdtt_params.fo_lim_overcorr;
															corr = Math.min(corr, lim_corr);
															if (corr > 0.0) {
																disparity_map[DISPARITY_INDEX_CM+0][tIndex] += corr;
															}
														}
													}
												}
											}
										}
										if      (corr_mode == 0) extra_disparity = disparity_map[DISPARITY_INDEX_INT][tIndex];
										else if (corr_mode == 1) extra_disparity = disparity_map[DISPARITY_INDEX_CM][tIndex];
										else if (corr_mode == 2) extra_disparity = disparity_map[DISPARITY_INDEX_POLY][tIndex];
										else if (corr_mode == 3) extra_disparity = disparity_map[DISPARITY_INDEX_HOR][tIndex];
										else if (corr_mode == 4) extra_disparity = disparity_map[DISPARITY_INDEX_VERT][tIndex];
										if (Double.isNaN(extra_disparity)) extra_disparity = 0;

										if (clt_mismatch != null){
											for (int pair = 0; pair < corr_pairs.length; pair++) if (((corr_mask >> pair) & 1) != 0){
												icorr_max =getMaxXYInt( // find integer pair or null if below threshold
														tcorr_partial[pair][numcol],      // [data_size * data_size]
														corr_size,
														min_corr,    // minimal value to consider (at integer location, not interpolated)
														debugMax);
												if (icorr_max == null){
													clt_mismatch[3*pair + 0 ][tIndex] = Double.NaN;
													clt_mismatch[3*pair + 1 ][tIndex] = Double.NaN;
													clt_mismatch[3*pair + 2 ][tIndex] = Double.NaN;
												} else {
													double [] corr_max_XYmp = getMaxXYCm( // get fractional center as a "center of mass" inside circle/square from the integer max
															tcorr_partial[pair][numcol],      // [data_size * data_size]
															corr_size,
															icorr_max, // integer center coordinates (relative to top left)
															max_corr_radius,  // positive - within that distance, negative - within 2*(-radius)+1 square
															max_corr_double, //"Double pass when masking center of mass to reduce preference for integer values
															debugMax); // should never return null
													// Only use Y components for pairs 0,1 and X components - for pairs 2,3
													double yp,xp;
													if (corr_pairs[pair][2] > 0){ // transpose - switch x <-> y
														yp = transform_size - 1 -corr_max_XYmp[0] - disparity_map[DISPARITY_INDEX_CM][tIndex];
														xp = transform_size - 1 -corr_max_XYmp[1]; // do not compare to average - it should be 0 anyway

													} else {
														xp = transform_size - 1 -corr_max_XYmp[0] - disparity_map[DISPARITY_INDEX_CM][tIndex];
														yp = transform_size - 1 -corr_max_XYmp[1]; // do not compare to average - it should be 0 anyway
													}
													double strength = tcorr_partial[pair][numcol][max_index]; // using the new location than for combined
													clt_mismatch[3*pair + 0 ][tIndex] = xp;
													clt_mismatch[3*pair + 1 ][tIndex] = yp;
													clt_mismatch[3*pair + 2 ][tIndex] = strength;
												}
											}
										}
									}
								}

							}
						} // end of if (clt_corr_combo != null)

						if (texture_tiles !=null) {

//							if ((extra_disparity != 0) && (((1 << FORCE_DISPARITY_BIT) & tile_op[tileY][tileX]) == 0)){ // 0 - adjust disparity, 1 - use provided
							if ((extra_disparity != 0) && !getForcedDisparity(tile_op[tileY][tileX])){ // 0 - adjust disparity, 1 - use provided
								// shift images by 0.5 * extra disparity in the diagonal direction
								for (int chn = 0; chn <numcol; chn++) { // color
									for (int i = 0; i < quad; i++) {
										fract_shift(    // fractional shift in transform domain. Currently uses sin/cos - change to tables with 2? rotations
												clt_data[i][chn][tileY][tileX], // double  [][]  clt_tile,
												transform_size,
												extra_disparity * port_offsets[i][0] / corr_magic_scale,     // double        shiftX,
												extra_disparity * port_offsets[i][1] / corr_magic_scale,     // double        shiftY,
												//									(globalDebugLevel > 0) && (tileX == debug_tileX) && (tileY == debug_tileY)); // external tile compare
												((globalDebugLevel > 0) && (chn==0) && (tileX >= debug_tileX - 2) && (tileX <= debug_tileX + 2) &&
														(tileY >= debug_tileY - 2) && (tileY <= debug_tileY+2)));
									}
								}
							}
							// lpf tiles (same as images before)
							// iclt tiles
							double [][][] iclt_tile = new double [quad][numcol][];
							double [] clt_tile;
							double scale = 0.25;  // matching iclt_2d
							for (int i = 0; i < quad; i++) {
								for (int chn = 0; chn <numcol; chn++) { // color
									// double [] clt_tile = new double [transform_size*transform_size];
									for (int dct_mode = 0; dct_mode < 4; dct_mode++){
										clt_tile = clt_data[i][chn][tileY][tileX][dct_mode].clone();
										// lpf each of the 4 quadrants before idct
										for (int j = 0; j < filter.length; j++){
											clt_tile[j] *= scale*filter[j];
										}
										// IDCT-IV should be in reversed order: CC->CC, SC->CS, CS->SC, SS->SS
										int idct_mode = ((dct_mode << 1) & 2) | ((dct_mode >> 1) & 1);
										clt_tile = dtt.dttt_iv  (clt_tile, idct_mode, transform_size);
										// iclt_tile[i][chn] = dtt.dttt_iv  (clt_data[i][chn][tileY][tileX][dct_mode], idct_mode, transform_size);
										double [] tile_mdct = dtt.unfold_tile(clt_tile, transform_size, dct_mode); // mode=0 - DCCT 16x16
										// accumulate partial mdct results
										if (dct_mode == 0){
											iclt_tile[i][chn] = tile_mdct;
										} else{
											for (int j = 0; j<tile_mdct.length; j++){
												iclt_tile[i][chn][j] += tile_mdct[j]; // matching iclt_2d
											}
										}
									}
								}
							}
							if ((globalDebugLevel > 0) && debugTile) {
								ShowDoubleFloatArrays sdfa_instance = new ShowDoubleFloatArrays(); // just for debugging?
								String [] titles = {"red0","blue0","green0","red1","blue1","green1","red2","blue2","green2","red3","blue3","green3"};
								double [][] dbg_tile = new double [quad*numcol][];
								for (int i = 0; i < quad; i++) {
									for (int chn = 0; chn <numcol; chn++) { // color
										dbg_tile[i * numcol + chn] = iclt_tile[i][chn];
									}
								}
								sdfa_instance.showArrays(dbg_tile, 2* transform_size, 2* transform_size, true, "iclt_x"+tileX+"_y"+tileY, titles);
							}


							// "de-bayer" tiles for matching, use original data for output
							double [][][] tiles_debayered = new double [quad][numcol][];
							for (int i =0; i<quad; i++){
								for (int chn = 0; chn < numcol; chn++){
									//								tiles_debayered[i][chn] =  tile_debayer(
									//										(chn != 2), // red or blue (false - green)
									//										iclt_tile[i][chn],
									//										2 * transform_size);

									tiles_debayered[i][chn] =  tile_debayer_shot_corr(
											(chn != 2), // red or blue (false - green)
											iclt_tile[i][chn],
											2 * transform_size,
											lt_window2, // squared lapping window
											min_shot,   // 10.0;  // Do not adjust for shot noise if lower than
											scale_shot,  //3.0;   // scale when dividing by sqrt
											lt_window2); // re-apply window to the result
								}
							}
							if ((globalDebugLevel > 0) && debugTile) {
								ShowDoubleFloatArrays sdfa_instance = new ShowDoubleFloatArrays(); // just for debugging?
								String [] titles = {"red0","blue0","green0","red1","blue1","green1","red2","blue2","green2","red3","blue3","green3"};
								double [][] dbg_tile = new double [quad*numcol][];
								for (int i = 0; i < quad; i++) {
									for (int chn = 0; chn <numcol; chn++) { // color
										dbg_tile[i * numcol + chn] = tiles_debayered[i][chn];
									}
								}
								sdfa_instance.showArrays(dbg_tile, 2* transform_size, 2* transform_size, true, "tiles_debayered_x"+tileX+"_y"+tileY, titles);
							}

							double []     max_diff = null;
							if ((disparity_map != null) && (disparity_map.length >= (IMG_DIFF0_INDEX + quad))){
								max_diff = new double[quad];
							}
							texture_tiles[tileY][tileX] =  tile_combine_rgba(
									tiles_debayered, // iclt_tile,      // [port][numcol][256]
									null, // double []     ports_rgb,      // average values of R,G,B for each camera (R0,R1,...,B2,B3)

									max_diff,        // maximal (weighted) deviation of each channel from the average
									lt_window2,      // [256]
									port_offsets,    // [port]{x_off, y_off}
									img_mask,        // which port to use, 0xf - all 4 (will modify as local variable)
									diff_sigma,      // pixel value/pixel change
									diff_threshold,  // pixel value/pixel change
									diff_gauss,      // when averaging images, use gaussian around average as weight (false - sharp all/nothing)
									min_agree,       // minimal number of channels to agree on a point (real number to work with fuzzy averages)
									col_weights,     // color channel weights, sum == 1.0
									dust_remove,     // boolean dust_remove,    // Do not reduce average weight when only one image differes much from the average
									keep_weights,    // keep_weights);   // return channel weights after A in RGBA
									(globalDebugLevel > 0) && debugTile);

							// mix RGB from iclt_tile, mix alpha with - what? correlation strength or 'don't care'? good correlation or all > min?
							for (int i = 0; i < iclt_tile[0][0].length; i++ ) {
								double sw = 0.0;
								for (int ip = 0; ip < quad; ip++) {
									sw += texture_tiles[tileY][tileX][numcol+1+ip][i];
								}
								if (sw != 0 ) sw = 1.0/sw;
								for (int chn = 0; chn <numcol; chn++) { // color
									texture_tiles[tileY][tileX][chn][i] = 0.0; //iclt[tileY][tileX][chn]
									for (int ip = 0; ip < quad; ip++) {
										texture_tiles[tileY][tileX][chn][i] += sw * texture_tiles[tileY][tileX][numcol+1+ip][i] * iclt_tile[ip][chn][i];
									}
								}
							}
							if ((disparity_map != null) && (disparity_map.length >= (IMG_DIFF0_INDEX + quad))){
								for (int i = 0; i < max_diff.length; i++){
									disparity_map[IMG_DIFF0_INDEX + i][tIndex] = max_diff[i];
								}
							}
						}
					}
				}
			};
		}
		startAndJoin(threads);
/*
		if (dbg_ports_coords != null) {
			(new showDoubleFloatArrays()).showArrays(dbg_ports_coords,  tilesX, tilesY, true, "ports_coordinates", dbg_titles);
		}
*/
		return clt_data;
	}

	/**
	 * Calculate disparity and strength for a inter-camera phase correlation of a pair of quad-cameras
	 * @param clt_parameters various configuration parameters
	 * @param fatzero - phase correlation fat zero (higher - ~LPF)
	 * @param corr2d Instance of the 2d correlator class
     * @param clt_data_tile_main aberration-corrected FD CLT data for one tile of the main quad camera  [sub-camera][color][quadrant][index]
     * @param clt_data_tile_aux aberration-corrected FD CLT data for one tile of the auxiliary quad camera  [sub-camera][color][quadrant][index]
	 * @param filter optional low-pass filter
	 * @param col_weights RBG color weights in combined phase correlation
	 * @param ml hwidth Optional to output data for ML: half width of the output tile (0 -> 1x1, 1-> 3x3, 2->5x5). Used only if center_corr != null
	 * @param ml_center_corr Optional to output data for ML:  output array [(2*hwidth+1)*(2*hwidth+1)] or null
	 * @param tcorr_combo if not null then tcorr_combo will contain full 2d correlation for debugging (now 15x15 in line scan order)
	 * @param notch_mode - detect vertical poles, use a notch filter for Y, normal - for X
	 * @param tileX debug tile X
	 * @param tileY debug tile X
	 * @param debugLevel debug level
	 * @return
	 */

	public double [] tileInterCamCorrs(
			final boolean                                   no_int_x0, // do not offset window to integer - used when averaging low textures to avoid "jumps" for very wide
			final CLTParameters  clt_parameters,
			final double                                    fatzero,         // May use correlation fat zero from 2 different parameters - fat_zero and rig.ml_fatzero
			final Correlation2d                             corr2d,
    		final double [][][][]                           clt_data_tile_main,
    		final double [][][][]                           clt_data_tile_aux,
			final double []       							filter,
			final double []       							col_weights,
    		final int                                       ml_hwidth,
    		final double []                                 ml_center_corr,
			final double [][]                               tcorr_combo,
			final boolean                                   notch_mode,
			final int             							tileX, // only used in debug output
			final int             							tileY,
			final int             							debugLevel) {
		double [] inter_cam_corr = corr2d.correlateInterCamerasFD(
				clt_data_tile_main,       // double [][][][]     clt_data_tile_main,
				clt_data_tile_aux,        // double [][][][]     clt_data_tile_aux,
	    		filter,                   // double []           lpf,
	    		scale_strengths,
	    		col_weights,              // double []           col_weights,
	    		fatzero);                 // double              fat_zero)

		return tileInterCamCorrs(
				no_int_x0, // do not offset window to integer - used when averaging low textures to avoid "jumps" for very wide
				clt_parameters, // final EyesisCorrectionParameters.CLTParameters  clt_parameters,
				inter_cam_corr, // final double []                                 inter_cam_corr,
				corr2d,         // final Correlation2d                             corr2d,
				ml_hwidth,      // final int                                       ml_hwidth,
				ml_center_corr, // final double []                                 ml_center_corr,
				tcorr_combo,    // final double [][]                               tcorr_combo,
				notch_mode,     // final boolean                                   notch_mode,
				tileX,          // final int             							tileX, // only used in debug output
				tileY,          // final int             							tileY,
				debugLevel);    // final int             							debugLevel);

	}

	public double [] tileInterCamCorrs(
			final boolean                                   no_int_x0, // do not offset window to integer - used when averaging low textures to avoid "jumps" for very wide
			// maximums. That reduces the residual disparity, but works continuously when it is known the maximum should be near zero
			final CLTParameters  clt_parameters,
			final double []                                 inter_cam_corr,
			final Correlation2d                             corr2d,
    		final int                                       ml_hwidth,
    		final double []                                 ml_center_corr,
			final double [][]                               tcorr_combo,
			final boolean                                   notch_mode,
			final int             							tileX, // only used in debug output
			final int             							tileY,
			final int             							debugLevel) {
//		int strip_hight = notch_mode? clt_parameters.img_dtt.corr_strip_notch : clt_parameters.img_dtt.corr_strip_hight;
		double [] result = {Double.NaN,  0.0, Double.NaN, Double.NaN};

// using exacltly as for main/aux

		int all_pairs = 1;
		double [][] corrs = { inter_cam_corr}; // single-channel
	    double [][] strips = corr2d.scaleRotateInterpoateCorrelations(
	    		corrs,                          // double [][] correlations,
	    		all_pairs,                      // int         pairs_mask,
	    		clt_parameters.img_dtt.corr_strip_hight, //);    // int         hwidth);
	    		(debugLevel > 0) ? all_pairs:0); // debugMax);

	    // Combine strips for selected pairs. Now using only for all available pairs.
	    // Other combinations are used only if requested (clt_corr_partial != null)

	    double [] stripe_combo = corr2d.combineInterpolatedCorrelations(
	    		strips,                                 // double [][] strips,
	    		all_pairs,                              // int         pairs_mask,
	    		clt_parameters.img_dtt.corr_offset,     // double      offset);
	    		clt_parameters.img_dtt.twice_diagonal); //    		boolean     twice_diagonal)


	    // calculate CM maximums for all mixed channels
	    // First get integer correlation center, relative to the center
		int [] ixy_combo =  corr2d.getMaxXYInt( // find integer pair or null if below threshold
				stripe_combo,              // double [] data,
				true,                     // boolean   axis_only,
// reduce minimal weight when averaging
				(no_int_x0?0.5:1.0) * clt_parameters.img_dtt.min_corr,   //  double    minMax,    // minimal value to consider (at integer location, not interpolated)
				debugLevel > 0); // boolean   debug);

		double [] stripe_inter = stripe_combo;
		int [] ixy =  ixy_combo;
/*
		double [] stripe_inter = corr2d. scaleRotateInterpoateSingleCorrelation(
				inter_cam_corr,                           // double []   corr,
				strip_hight,                              // int         hwidth,
				Correlation2d.PAIR_HORIZONTAL,            // int         dir, // 0 - hor, 1 - vert, 2 - parallel to row = col (main) diagonal (0->3), 3 -2->1
				1,                                        // int         ss,
				(debugLevel > 0));                        // boolean     debug
*/
		if (ml_center_corr != null) {
			corr2d.corrCenterValues(
					ml_hwidth,
					inter_cam_corr,
					ml_center_corr);
		}
		if (tcorr_combo != null) {

			tcorr_combo[0] = inter_cam_corr;
			if (tcorr_combo.length > 1) {
				tcorr_combo[1] = corr2d.debugStrip(stripe_inter);
			}
			if (tcorr_combo.length > 2) {
				tcorr_combo[2] = corr2d.debugStrip2(stripe_inter);
			}
		}


	    // First get integer correlation center, relative to the center
// TODO: multiply/acummulate by Y window?
/*
		int [] ixy =  corr2d.getMaxXYInt(            // find integer pair or null if below threshold
				stripe_inter,                        // double [] data,
				true,                                // boolean   axis_only, for strip it is always true
				clt_parameters.img_dtt.min_corr,     //  double    minMax,    // minimal value to consider (at integer location, not interpolated)
				debugLevel > 0); // boolean   debug);
*/
		double [] corr_stat = null;
		// if integer argmax was strong enough, calculate CM argmax
		// will not fill out DISPARITY_INDEX_INT+1, DISPARITY_INDEX_CM+1, DISPARITY_INDEX_POLY+1
		// use clt_mismatch for that
		double strength = 0.0;
		double disparity = 0.0;
		if (ixy != null) {
			strength = stripe_inter[ixy[0]+clt_parameters.transform_size-1]; // strength at integer max on axis
			disparity =      -ixy[0];
			result[INDEX_STRENGTH] =  strength;
			result[INDEX_DISP] =      disparity;
			if (Double.isNaN(strength)) {
				System.out.println("BUG: 1. strength should not be NaN");
			}
//getMaxXCmNotch			if
			if (notch_mode) {
				corr_stat = corr2d.getMaxXCmNotch( // get fractional center as a "center of mass" inside circle/square from the integer max
						stripe_inter,              // double [] data,      // [data_size * data_size]
						0, // ixy[0],                    // int       ixcenter,  // integer center x
						// corr_wndy,              // double [] window_y,  // (half) window function in y-direction(perpendicular to disparity: for row0  ==1
						// corr_wndx,              // double [] window_x,  // half of a window function in x (disparity) direction
						(debugLevel > 0));         // boolean   debug);
			} else {
				corr_stat = corr2d.getMaxXCm(      // get fractional center as a "center of mass" inside circle/square from the integer max
						stripe_inter,              // double [] data,      // [data_size * data_size]
						(no_int_x0 ? 0:ixy[0]), // 0, // ixy[0],                    // int       ixcenter,  // integer center x
						// corr_wndy,              // double [] window_y,  // (half) window function in y-direction(perpendicular to disparity: for row0  ==1
						// corr_wndx,              // double [] window_x,  // half of a window function in x (disparity) direction
						(debugLevel > 0));         // boolean   debug);
			}
		}

		if (corr_stat != null) {
			if (debugLevel > 0){
				System.out.println(String.format("Tile X/Y = %d/%d corr_stat= {%f, %f, %f} ",tileX, tileY,corr_stat[0],corr_stat[1],corr_stat[2]));
				System.out.println("notch_mode ="+notch_mode);
			}


			disparity = -corr_stat[0]; // yes, uses this value
			result[INDEX_DISP] = disparity;
			double eff_radius = corr_stat[2] * clt_parameters.img_dtt.cm_max_normalization;
			strength = corr_stat[1]/(eff_radius * eff_radius); // total mass by square effective radius
			result[INDEX_STRENGTH] = strength;

			// see if strength is enough to proceed with LMA/poly (otherwise keep disp/strength
			if ((strength > clt_parameters.img_dtt.min_poly_strength) &&  !notch_mode) {  // for now notch mode is only for CM
				// create LMA instance, calculate LMA composite argmax
		    	// Create 2 groups: ortho & diag
		    	Correlations2dLMA lma = corr2d.corrLMA(
		    			clt_parameters.img_dtt,           // ImageDttParameters  clt_parameters.img_dtt,
		    			inter_cam_corr,                   // double []         corr,
		        		false,                            // boolean             run_poly_instead, // true - run LMA, false - run 2d polynomial approximation
		    			corr_stat[0],                     // double    xcenter,   // preliminary center x in pixels for largest baseline
		    			clt_parameters.img_dtt.ortho_vasw_pwr, // double    vasw_pwr,  // value as weight to this power,
		    			debugLevel,         // int                 debug_level,
		        		tileX,                        // int                 tileX, // just for debug output
		        		tileY );                      // int                 tileY
		    	double [] lma_disparity_strength = null;
		    	if (lma != null) {
			    	lma_disparity_strength = lma.getDisparityStrength();
		    		if (debugLevel > 0){
		    				System.out.println(String.format("Tile X/Y = %d/%d LMA disparity = %7.4f, strength = %7.4f",
		    						tileX, tileY,
		    						lma_disparity_strength[0],lma_disparity_strength[1]));
		    		}
					// if enabled overwrite - replace  DISPARITY_INDEX_CM and DISPARITY_STRENGTH_INDEX
					if (clt_parameters.rig.use_poly) {
						disparity = lma_disparity_strength[0];
						strength =  lma_disparity_strength[1];
						result[INDEX_STRENGTH] = strength;
						result[INDEX_DISP] =     disparity;
						if (Double.isNaN(strength)) {
							System.out.println("BUG: 2. strength should not be NaN");
						}
					}

					double [] corrXY = null;
					if (clt_parameters.rig.use_xy_poly) {
						corrXY= corr2d.single2dPoly( // returns x-xcenter, y, strength (sign same as disparity)
							clt_parameters.img_dtt, //  ImageDttParameters  imgdtt_params,
							inter_cam_corr, // double []           corr,
							// See if -disparity should be here as in mismatch
							0.0, // double              xcenter,   // -disparity to compare. use 0?
							clt_parameters.img_dtt.ortho_vasw_pwr, // double              vasw_pwr,  // value as weight to this power,
							debugLevel, // int                 debug_level,
							tileX, // int                 tileX,
							tileY); // int                 tileY
					} else {
						corrXY= corr2d.single2dCM( // returns x-xcenter, y, strength (sign same as disparity)
								clt_parameters.img_dtt, //  ImageDttParameters  imgdtt_params,
								inter_cam_corr, // double []           corr,
								// See if -disparity should be here as in mismatch
								0.0, // double              xcenter,   // -disparity to compare. use 0? //-ixy[0], //
								clt_parameters.max_corr_radius, // double              vasw_pwr,  // value as weight to this power,
								debugLevel, // int                 debug_level,
								tileX, // int                 tileX,
								tileY); // int                 tileY
					}
					if (corrXY != null) {
						result[INDEX_DX] =  corrXY[0];
						result[INDEX_DY] =  corrXY[1];
					}
		    	}
			}
		} // end of if (corr_stat != null)
		return result;
	}



	/**
	 * Calculate correlation/strength, start with center of mass (CM) for all available pairs, proceed with LMA
	 * if strength is sufficient. Calculate 4 directional correlation/strengths if requested and strong enough
	 * @param clt_parameters various configuration parameters
	 * @param fatzero phaase correlation fat zero (higher ~LPF)
	 * @param get4dirs request 4 directional correlations (horizontal, vertical main diagonal, other diagonal)
	 * @param corr2d Instance of the 2d correlator class
	 * @param clt_data aberration-corrected FD CLT data [camera][color][quadrant][index]
	 * @param filter optional low-pass filter
	 * @param col_weights RBG color weights in combined phase correlation
	 * @param ml hwidth Optional to output data for ML: half width of the output tile (0 -> 1x1, 1-> 3x3, 2->5x5). Used only if center_corr != null
	 * @param ml_center_corr Optional to output data for ML:  output array [(2*hwidth+1)*(2*hwidth+1)] or null
	 * @param tileX debug tile X
	 * @param tileY debug tile X
	 * @param debugLevel debug level
	 * @return {disparity, disp_hor, disp_vert, disp_diagm, disp_diago, strength, str_hor, str_vert, str_diagm, str_diago}
	 * indexed by DISP_*_INDEX and STR_*_INDEX constants
	 */
	public double [] tileCorrs(
			final CLTParameters       clt_parameters,
			final double          fatzero,         // May use correlation fat zero from 2 different parameters - fat_zero and rig.ml_fatzero
			final boolean         get4dirs, // calculate disparity/strength for each of the 4 directions
			final Correlation2d   corr2d,
			final double [][][][] clt_data,
			final double []       filter,
			final double []       col_weights,
    		final int             ml_hwidth,
    		final double [][]     ml_center_corr,
			final int             tileX, // only used in debug output
			final int             tileY,
			final int             debugLevel
			) {
		double [] result = new double [SNGL_DISPARITY_TITLES.length];
		for (int i:SNGL_DISPARITY_NAN) {
			result[i] = Double.NaN;
		}

		// calculate all selected pairs correlations
		int all_pairs = clt_parameters.img_dtt.dbg_pair_mask; //TODO: use tile tasks

	    double [][]  corrs = corr2d.correlateCompositeFD(
	    		clt_data,       // double [][][][] clt_data, // aberration-corrected FD CLT data for one tile [camera][color][quadrant][index]
	    		all_pairs,      // int                 pairs_mask,
	    		filter,         // double []           lpf,
	    		scale_strengths,
	    		col_weights,    // double []           col_weights,
	    		fatzero);       // double              fat_zero)

	    // calculate interpolated "strips" to match different scales and orientations (ortho/diagonal) on the
	    // fine (0.5 pix) grid. ortho for scale == 1 provide even/even samples (1/4 of all), diagonal ones -
	    // checkerboard pattern
	    if (ml_center_corr != null) {
	    	corr2d.corrCenterValues(
	    			ml_hwidth,                           // int         hwidth,
	    			clt_parameters.img_dtt.corr_offset,  //double      offset,
	        		corrs,                               // double [][] full_corr,
	        		ml_center_corr);                     // double [][] center_corr)
	    }

	    double [][] strips = corr2d.scaleRotateInterpoateCorrelations(
	    		corrs,                          // double [][] correlations,
	    		all_pairs,                      // int         pairs_mask,
	    		clt_parameters.img_dtt.corr_strip_hight, //);    // int         hwidth);
	    		(debugLevel > 0) ? all_pairs:0); // debugMax);

	    // Combine strips for selected pairs. Now using only for all available pairs.
	    // Other combinations are used only if requested (clt_corr_partial != null)

	    double [] strip_combo = corr2d.combineInterpolatedCorrelations(
	    		strips,                                 // double [][] strips,
	    		all_pairs,                              // int         pairs_mask,
	    		clt_parameters.img_dtt.corr_offset,     // double      offset);
	    		clt_parameters.img_dtt.twice_diagonal); //    		boolean     twice_diagonal)


	    // calculate CM maximums for all mixed channels
	    // First get integer correlation center, relative to the center
		int [] ixy =  corr2d.getMaxXYInt( // find integer pair or null if below threshold
				strip_combo,              // double [] data,
				true,                     // boolean   axis_only,
				clt_parameters.img_dtt.min_corr,   //  double    minMax,    // minimal value to consider (at integer location, not interpolated)
				debugLevel > 0); // boolean   debug);
		double [] corr_stat = null;

		// if integer argmax was strong enough, calculate CM argmax
		// will not fill out DISPARITY_INDEX_INT+1, DISPARITY_INDEX_CM+1, DISPARITY_INDEX_POLY+1
		// use clt_mismatch for that
		double strength = 0.0;
		double disparity = 0.0;
		if (ixy != null) {
			strength = strip_combo[ixy[0]+clt_parameters.transform_size-1]; // strength at integer max on axis
			disparity =      -ixy[0];
			result[STR_FULL_INDEX] = strength;
			result[DISP_FULL_INDEX] = disparity;
			if (Double.isNaN(strength)) {
				System.out.println("BUG: 1. strength should not be NaN");
			}
			// Ignores negative values!
			corr_stat = corr2d.getMaxXCm(   // get fractional center as a "center of mass" inside circle/square from the integer max
					strip_combo,            // double [] data,      // [data_size * data_size]
					ixy[0],                 // int       ixcenter,  // integer center x
					// corr_wndy,           // double [] window_y,  // (half) window function in y-direction(perpendicular to disparity: for row0  ==1
					// corr_wndx,           // double [] window_x,  // half of a window function in x (disparity) direction
					(debugLevel > 0));      // boolean   debug);
		}
// removed HOR/VERT
		// proceed only if CM correlation result is non-null // for compatibility with old code we need it to run regardless of the strength of the normal correlation
		if (corr_stat != null) {
			disparity = -corr_stat[0];
			result[DISP_FULL_INDEX] = disparity;
			// see if strength is enough to proceed with LMA/poly (otherwise keep disp/strength
			if (strength > clt_parameters.img_dtt.min_poly_strength) {
				// create LMA instance, calculate LMA composite argmax
		    	// Create 2 groups: ortho & diag
		    	Correlations2dLMA lma = corr2d.corrLMA(
		    			clt_parameters.img_dtt,                // ImageDttParameters  clt_parameters.img_dtt,
		    			corrs,                        // double [][]         corrs,
		    			clt_parameters.img_dtt.dbg_pair_mask,  // int                 pair_mask, // which pairs to process
		        		false,                        // boolean             run_poly_instead, // true - run LMA, false - run 2d polynomial approximation
		    			corr_stat[0],                 // double    xcenter,   // preliminary center x in pixels for largest baseline
		    			clt_parameters.img_dtt.ortho_vasw_pwr, // double    vasw_pwr,  // value as weight to this power,
		    			debugLevel,         // int                 debug_level,
		        		tileX,                        // int                 tileX, // just for debug output
		        		tileY );                      // int                 tileY
		    	double [] lma_disparity_strength = null;
		    	if (lma != null) {
//		    		double []   mod_disparity_diff = null;
		    		double [][] dir_corr_strength =  null;
			    	lma_disparity_strength = lma.getDisparityStrength();
		    		if (debugLevel > 0){
		    				System.out.println(String.format("Tile X/Y = %d/%d LMA disparity = %7.4f, strength = %7.4f",
		    						tileX, tileY,
		    						lma_disparity_strength[0],lma_disparity_strength[1]));
		    		}
					// if enabled overwrite - replace  DISPARITY_INDEX_CM and DISPARITY_STRENGTH_INDEX
					if (clt_parameters.rig.use_poly) {
						disparity = lma_disparity_strength[0];
						strength =  lma_disparity_strength[1];
						result[STR_FULL_INDEX] = strength;
						result[DISP_FULL_INDEX] = disparity;
						if (Double.isNaN(strength)) {
							System.out.println("BUG: 2. strength should not be NaN");
						}
					}

                    // Correction for far foreground objects
//					if ((clt_parameters.img_dtt.fo_correct && (strength > 0 * clt_parameters.img_dtt.fo_min_strength)) || get4dirs) {
					// no fo_correct for the rig!
					if (get4dirs) {
			    		// try all dirs:
						dir_corr_strength = corr2d.corr4dirsLMA(
			    				clt_parameters.img_dtt,                // ImageDttParameters  clt_parameters.img_dtt,
			    				corrs,                        // double [][]         corrs,
			    				clt_parameters.img_dtt.dbg_pair_mask,  // int                 pair_mask, // which pairs to process
			    				-disparity,                   // double    xcenter,   // preliminary center x in pixels for largest baseline
			    				clt_parameters.img_dtt.ortho_vasw_pwr, // double    vasw_pwr,  // value as weight to this power,
			    				debugLevel,         // int                 debug_level,
			    				tileX,                        // int                 tileX, // just for debug output
			    				tileY );                      // int                 tileY
			    		if ((debugLevel > 0) && (dir_corr_strength != null)) {
			    			double [] nan2 = {Double.NaN, Double.NaN, Double.NaN, Double.NaN};
			    			for (int ii = 0; ii < dir_corr_strength.length; ii++) {
			    				if (dir_corr_strength[ii] == null) dir_corr_strength[ii] = nan2;
			    			}
			    			System.out.println(String.format("corr4dirsLMA -> ↔: %7.4f (%7.4f) ↕:%7.4f (%7.4f) ⤡:%7.4f (%7.4f) ⤢:%7.4f (%7.4f)",
			    					dir_corr_strength[0][0],dir_corr_strength[0][1],dir_corr_strength[1][0],dir_corr_strength[1][1],
			    					dir_corr_strength[2][0],dir_corr_strength[2][1],dir_corr_strength[3][0],dir_corr_strength[3][1]));
			    		}
			    		if (dir_corr_strength[0] != null) {
			    			result[DISP_HOR_INDEX] =   dir_corr_strength[0][0];
			    			result[STR_HOR_INDEX] =    dir_corr_strength[0][1];
			    		} else {
			    			result[DISP_HOR_INDEX] =   Double.NaN;
			    			result[STR_HOR_INDEX] =    0.0;
			    		}
			    		if (dir_corr_strength[1] != null) {
			    			result[DISP_VERT_INDEX] =  dir_corr_strength[1][0];
			    			result[STR_VERT_INDEX] =   dir_corr_strength[1][1];
			    		} else {
			    			result[DISP_VERT_INDEX] =  Double.NaN;
			    			result[STR_VERT_INDEX] =   0.0;
			    		}
			    		if (dir_corr_strength[2] != null) {
			    			result[DISP_DIAGM_INDEX] = dir_corr_strength[2][0];
			    			result[STR_DIAGM_INDEX] =  dir_corr_strength[2][1];
			    		} else {
			    			result[DISP_DIAGM_INDEX] = Double.NaN;
			    			result[STR_DIAGM_INDEX] =  0.0;
			    		}
			    		if (dir_corr_strength[3] != null) {
			    			result[DISP_DIAGO_INDEX] = dir_corr_strength[3][0];
			    			result[STR_DIAGO_INDEX] =  dir_corr_strength[3][1];
			    		} else {
			    			result[DISP_DIAGO_INDEX] = Double.NaN;
			    			result[STR_DIAGO_INDEX] =  0.0;

			    		}
/*
			    		mod_disparity_diff =     corr2d.foregroundCorrect(
			    				clt_parameters.img_dtt.fo_far,            // boolean   bg,
			    				clt_parameters.img_dtt.fo_ortho,          // boolean   ortho,
			    				dir_corr_strength,               // double [] dir_disp,
			    	    		disparity,                       // double    full_disp,
			    	    		clt_parameters.img_dtt.fo_min_strength,   // double      min_strength,
			    	    		clt_parameters.img_dtt.fo_min_eff,        // double      min_eff,
			    	    		clt_parameters.img_dtt.fo_min_eff_ratio,  // double      min_eff_ratio,
			    	    		clt_parameters.img_dtt.fo_max_hwidth,    // double      max_hwidth, //  =          3.0;  // maximal half-width disparity  direction to try to correct far objects
			    	    		clt_parameters.img_dtt.fo_min_diff,       // double    fo_min_diff,
			    	    		clt_parameters.img_dtt.fo_overcorrection, // double    fo_overcorrection,
			    	    		clt_parameters.img_dtt.fo_lim_overcorr,   // double    fo_lim_overcorr,
			    	    		(debugLevel > 0) );    // boolean debug);

						if (mod_disparity_diff[0] != disparity){ // if it changed
							if (clt_parameters.img_dtt.fo_correct && (strength > clt_parameters.img_dtt.fo_min_strength)) { // always
								disparity = mod_disparity_diff[0];
								result[DISP_FULL_INDEX] = disparity;
							}
						}
*/
					}
		    	}
			}
		} // end of if (corr_stat != null)
		return result;
	}


	public void generateTextureTiles(
			final CLTParameters       clt_parameters,
			final double           extra_disparity,
			final int              quad,      // number of subcameras
			final int              numcol, // number of colors
			int                    img_mask,
			final int [][]         tile_op,         // [tilesY][tilesX] - what to do - 0 - nothing for this tile
			final double [][][][]  clt_data,
			final double [][][][]  texture_tiles,   // [tilesY][tilesX]["RGBA".length()][];  null - will skip images combining
			final double []        filter,
			final double []        lt_window2,
			final double[][]       port_offsets,
			final double []        col_weights,
			final DttRad2          dtt,
			final int              tileX, // only used in debug output
			final int              tileY,
			final int              debugLevel
			) {

		if ((extra_disparity != 0) && !getForcedDisparity(tile_op[tileY][tileX])){ // 0 - adjust disparity, 1 - use provided
			// shift images by 0.5 * extra disparity in the diagonal direction
			for (int chn = 0; chn <numcol; chn++) { // color
				for (int i = 0; i < quad; i++) {
					fract_shift(    // fractional shift in transform domain. Currently uses sin/cos - change to tables with 2? rotations
							clt_data[i][chn], // [tileY][tileX], // double  [][]  clt_tile,
							clt_parameters.transform_size,
							extra_disparity * port_offsets[i][0] / clt_parameters.corr_magic_scale,     // double        shiftX,
							extra_disparity * port_offsets[i][1] / clt_parameters.corr_magic_scale,     // double        shiftY,
							debugLevel > 0);
				}
			}
		}
		// lpf tiles (same as images before)
		// iclt tiles
		double [][][] iclt_tile = new double [quad][numcol][];
		double [] clt_tile;
		double scale = 0.25;  // matching iclt_2d
		for (int i = 0; i < quad; i++) {
			for (int chn = 0; chn <numcol; chn++) { // color
				// double [] clt_tile = new double [transform_size*transform_size];
				for (int dct_mode = 0; dct_mode < 4; dct_mode++){
					clt_tile = clt_data[i][chn][dct_mode].clone();
					// lpf each of the 4 quadrants before idct
					for (int j = 0; j < filter.length; j++){
						clt_tile[j] *= scale*filter[j];
					}
					// IDCT-IV should be in reversed order: CC->CC, SC->CS, CS->SC, SS->SS
					int idct_mode = ((dct_mode << 1) & 2) | ((dct_mode >> 1) & 1);
					clt_tile = dtt.dttt_iv  (clt_tile, idct_mode, clt_parameters.transform_size);
					double [] tile_mdct = dtt.unfold_tile(clt_tile, clt_parameters.transform_size, dct_mode); // mode=0 - DCCT 16x16
					// accumulate partial mdct results
					if (dct_mode == 0){
						iclt_tile[i][chn] = tile_mdct;
					} else{
						for (int j = 0; j<tile_mdct.length; j++){
							iclt_tile[i][chn][j] += tile_mdct[j]; // matching iclt_2d
						}
					}
				}
			}
		}
		if (debugLevel > 0) {
			ShowDoubleFloatArrays sdfa_instance = new ShowDoubleFloatArrays(); // just for debugging?
			String [] titles = {"red0","blue0","green0","red1","blue1","green1","red2","blue2","green2","red3","blue3","green3"};
			double [][] dbg_tile = new double [quad*numcol][];
			for (int i = 0; i < quad; i++) {
				for (int chn = 0; chn <numcol; chn++) { // color
					dbg_tile[i * numcol + chn] = iclt_tile[i][chn];
				}
			}
			sdfa_instance.showArrays(dbg_tile, 2* clt_parameters.transform_size, 2* clt_parameters.transform_size, true, "iclt_x"+tileX+"_y"+tileY, titles);
		}


		// "de-bayer" tiles for matching, use original data for output
		double [][][] tiles_debayered = new double [quad][numcol][];
		for (int i =0; i<quad; i++){
			for (int chn = 0; chn < numcol; chn++){

				tiles_debayered[i][chn] =  tile_debayer_shot_corr(
						(chn != 2),                 // red or blue (false - green)
						iclt_tile[i][chn],
						2 * clt_parameters.transform_size,
						lt_window2,                 // squared lapping window
						clt_parameters.min_shot,    // 10.0;  // Do not adjust for shot noise if lower than
						clt_parameters.scale_shot,  //3.0;   // scale when dividing by sqrt
						lt_window2);                // re-apply window to the result
			}
		}
		if (debugLevel > 0) {
			ShowDoubleFloatArrays sdfa_instance = new ShowDoubleFloatArrays(); // just for debugging?
			String [] titles = {"red0","blue0","green0","red1","blue1","green1","red2","blue2","green2","red3","blue3","green3"};
			double [][] dbg_tile = new double [quad*numcol][];
			for (int i = 0; i < quad; i++) {
				for (int chn = 0; chn <numcol; chn++) { // color
					dbg_tile[i * numcol + chn] = tiles_debayered[i][chn];
				}
			}
			sdfa_instance.showArrays(dbg_tile, 2* clt_parameters.transform_size, 2* clt_parameters.transform_size, true, "tiles_debayered_x"+tileX+"_y"+tileY, titles);
		}

//		double []     max_diff = null;
//		if ((disparity_map != null) && (disparity_map.length >= (IMG_DIFF0_INDEX + quad))){
//			max_diff = new double[quad];
//		}
		texture_tiles[tileY][tileX] =  tile_combine_rgba(
				tiles_debayered,                // iclt_tile,      // [port][numcol][256]
				null,                           // double []     ports_rgb,      // average values of R,G,B for each camera (R0,R1,...,B2,B3)
				null,                           // max_diff,        // maximal (weighted) deviation of each channel from the average
				lt_window2,                     // [256]
				port_offsets,                   // [port]{x_off, y_off}
				img_mask,                       // which port to use, 0xf - all 4 (will modify as local variable)
				clt_parameters.diff_sigma,      // pixel value/pixel change
				clt_parameters.diff_threshold,  // pixel value/pixel change
				clt_parameters.diff_gauss,      // when averaging images, use gaussian around average as weight (false - sharp all/nothing)
				clt_parameters.min_agree,       // minimal number of channels to agree on a point (real number to work with fuzzy averages)
				col_weights,                    // color channel weights, sum == 1.0
				clt_parameters.dust_remove,     // boolean dust_remove,    // Do not reduce average weight when only one image differes much from the average
				clt_parameters.keep_weights,    // keep_weights);   // return channel weights after A in RGBA
				(debugLevel > 0));

		// mix RGB from iclt_tile, mix alpha with - what? correlation strength or 'don't care'? good correlation or all > min?
		for (int i = 0; i < iclt_tile[0][0].length; i++ ) {
			double sw = 0.0;
			for (int ip = 0; ip < quad; ip++) {
				sw += texture_tiles[tileY][tileX][numcol+1+ip][i];
			}
			if (sw != 0 ) sw = 1.0/sw;
			for (int chn = 0; chn <numcol; chn++) { // color
				texture_tiles[tileY][tileX][chn][i] = 0.0; //iclt[tileY][tileX][chn]
				for (int ip = 0; ip < quad; ip++) {
					texture_tiles[tileY][tileX][chn][i] += sw * texture_tiles[tileY][tileX][numcol+1+ip][i] * iclt_tile[ip][chn][i];
				}
			}
		}
	}



//	public double [][][][][][] clt_bi_quad(

	public double [][][][][][][]  clt_bi_quad_dbg(
			final CLTParameters       clt_parameters,
			final double              fatzero,         // May use correlation fat zero from 2 different parameters - fat_zero and rig.ml_fatzero
			final boolean             notch_mode,      // use notch filter for inter-camera correlation to detect poles
			final int                 lt_rad,          // low texture mode - inter-correlation is averaged between the neighbors before argmax-ing, using (2*notch_mode+1)^2 square
			final boolean             no_int_x0,       // do not offset window to integer maximum - used when averaging low textures to avoid "jumps" for very wide
			final int [][]            tile_op,         // [tilesY][tilesX] - what to do - 0 - nothing for this tile
			final double [][]         disparity_array, // [tilesY][tilesX] - individual per-tile expected disparity
			final double [][][]       image_data_main, // first index - number of image in a quad
			final double [][][]       image_data_aux,  // first index - number of image in a quad
		    final boolean [][]        saturation_main, // (near) saturated pixels or null
		    final boolean [][]        saturation_aux,  // (near) saturated pixels or null
			 // correlation results - combo will be for the correlation between two quad cameras
			final double [][][][]     clt_corr_combo,  // [type][tilesY][tilesX][(2*transform_size-1)*(2*transform_size-1)] // if null - will not calculate
			                                           // [type][tilesY][tilesX] should be set by caller
													   // types: 0 - selected correlation (product+offset), 1 - sum
			final double [][]         disparity_bimap, // [23][tilesY][tilesX], only [6][] is needed on input or null - do not calculate
			                                           // last 2 - contrast, avg/ "geometric average)
			final double [][]         ml_data,         // data for ML - 18 layers - 4 center areas (3x3, 5x5,..) per camera-per direction, 1 - composite, and 1 with just 1 data (target disparity)
			final double [][][][]     texture_tiles_main, // [tilesY][tilesX]["RGBA".length()][];  null - will skip images combining
			final double [][][][]     texture_tiles_aux,  // [tilesY][tilesX]["RGBA".length()][];  null - will skip images combining
			final int                 width,           // may be not multiple of 8, same for the height

			final GeometryCorrection  geometryCorrection_main,
			final GeometryCorrection  geometryCorrection_aux,
			final double [][][][][][] clt_kernels_main, // [channel_in_quad][color][tileY][tileX][band][pixel] , size should match image (have 1 tile around)
			final double [][][][][][] clt_kernels_aux,  // [channel_in_quad][color][tileY][tileX][band][pixel] , size should match image (have 1 tile around)
			final double              corr_magic_scale, // still not understood coefficient that reduces reported disparity value.  Seems to be around 0.85
			final boolean             keep_clt_data,
//			final int [][]            woi_tops,
			final double [][][]       ers_delay,        // if not null - fill with tile center acquisition delay
			final int                 threadsMax,  // maximal number of threads to launch
			final int                 debugLevel,
			final double [][][]       port_xy_main_dbg, // for each tile/port save x,y pixel coordinates (gpu code development)
			final double [][][]       port_xy_aux_dbg) // for each tile/port save x,y pixel coordinates (gpu code development)
	{
		final int                 globalDebugLevel = clt_parameters.rig.rig_mode_debug?debugLevel:-2;
		final int                 debug_tileX = clt_parameters.tileX;
		final int                 debug_tileY = clt_parameters.tileY;
		final int quad_main = image_data_main.length;   // number of subcameras
		final int quad_aux =  image_data_aux.length;   // number of subcameras
		final int numcol = 3; // number of colors
		final int nChn = image_data_main[0].length;
		final int height=image_data_main[0][0].length/width;
		final int tilesX=width/clt_parameters.transform_size;
		final int tilesY=height/clt_parameters.transform_size;
		final int nTilesInChn=tilesX*tilesY;
		// clt_data does not need to be for the whole image (no, it is used for textures)
		final double [][][][][][][] clt_bidata = (keep_clt_data)? (new double[2][][][][][][]):null;
		if (clt_bidata != null) {
			clt_bidata[0] = new double[quad_main][nChn][tilesY][tilesX][][];
			clt_bidata[1] = new double[quad_aux][nChn][tilesY][tilesX][][];
		}

		final double [][] lt_corr = (lt_rad > 0)? (new double [nTilesInChn][]):null; // will keep inter-camera combo correlation, later combined in a separate multi-thread run

		final Thread[] threads = newThreadArray(threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		final double [] col_weights= new double [numcol]; // colors are RBG
		col_weights[2] = 1.0/(1.0 + clt_parameters.corr_red + clt_parameters.corr_blue);    // green color
		col_weights[0] = clt_parameters.corr_red *  col_weights[2];
		col_weights[1] = clt_parameters.corr_blue * col_weights[2];
		final int corr_size = clt_parameters.transform_size * 2 -1;
		final int [][] transpose_indices = new int [corr_size*(corr_size-1)/2][2];
		int indx = 0;
		for (int i =0; i < corr_size-1; i++){
			for (int j = i+1; j < corr_size; j++){
				transpose_indices[indx  ][0] = i * corr_size + j;
				transpose_indices[indx++][1] = j * corr_size + i;
			}
		}

		// get ml_data half width
		final int ml_hwidth = (ml_data != null)?(((int) Math.round(Math.sqrt(ml_data[0].length/nTilesInChn)) - 1) / 2):0;

		// Create window  to select center correlation strip using
		// ortho_height - full width of non-zero elements
		// ortho_eff_height - effective height (ration of the weighted column sum to the center value)
		int wcenter = clt_parameters.transform_size - 1;
		final double [] ortho_weights = new double [corr_size]; // [15]
		for (int i = 0; i < corr_size; i++){
			if ((i >= wcenter - clt_parameters.img_dtt.ortho_height/2) && (i <= wcenter + clt_parameters.img_dtt.ortho_height/2)) {
				double dx = 1.0*(i-wcenter)/(clt_parameters.img_dtt.ortho_height/2 + 1);
				ortho_weights[i] = 0.5*(1.0+Math.cos(Math.PI*dx))/clt_parameters.img_dtt.ortho_eff_height;
			}
		}
		if (globalDebugLevel > 0){
			System.out.println("ortho_height="+ clt_parameters.img_dtt.ortho_height+" ortho_eff_height="+ clt_parameters.img_dtt.ortho_eff_height);
			for (int i = 0; i < corr_size; i++){
				System.out.println(" ortho_weights["+i+"]="+ ortho_weights[i]);
			}
		}

		if (globalDebugLevel > -2) {
			System.out.println("clt_aberrations_quad_corr(): width="+width+" height="+height+" transform_size="+clt_parameters.transform_size+
					" debug_tileX="+debug_tileX+" debug_tileY="+debug_tileY+" globalDebugLevel="+globalDebugLevel);
		}
///		final int [][] corr_pairs ={ // {first, second, rot} rot: 0 - as is, 1 - swap y,x
///				{0,1,0},
///				{2,3,0},
///				{0,2,1},
///				{1,3,1}};

		final double[][] port_offsets = {
				{-0.5, -0.5},
				{ 0.5, -0.5},
				{-0.5,  0.5},
				{ 0.5,  0.5}};
		final int transform_len = clt_parameters.transform_size * clt_parameters.transform_size;



		final double [] filter_direct= new double[transform_len];
		if (clt_parameters.getCorrSigma(isMonochrome()) == 0) {
			filter_direct[0] = 1.0;
			for (int i= 1; i<filter_direct.length;i++) filter_direct[i] =0;
		} else {
			for (int i = 0; i < clt_parameters.transform_size; i++){
				for (int j = 0; j < clt_parameters.transform_size; j++){
					filter_direct[i * clt_parameters.transform_size+j] = Math.exp(-(i*i+j*j)/(2*clt_parameters.getCorrSigma(isMonochrome()))); // FIXME: should be sigma*sigma !
				}
			}
		}
		// normalize
		double sum = 0;
		for (int i = 0; i < clt_parameters.transform_size; i++){
			for (int j = 0; j < clt_parameters.transform_size; j++){
				double d = 	filter_direct[i*clt_parameters.transform_size+j];
				d*=Math.cos(Math.PI*i/(2*clt_parameters.transform_size))*Math.cos(Math.PI*j/(2*clt_parameters.transform_size));
				if (i > 0) d*= 2.0;
				if (j > 0) d*= 2.0;
				sum +=d;
			}
		}
		for (int i = 0; i<filter_direct.length; i++){
			filter_direct[i] /= sum;
		}

		DttRad2 dtt = new DttRad2(clt_parameters.transform_size);
		final double [] filter= dtt.dttt_iiie(filter_direct);
		for (int i=0; i < filter.length;i++) filter[i] *= 2*clt_parameters.transform_size;

		// prepare disparity maps and weights
		final int max_search_radius = (int) Math.abs(clt_parameters.max_corr_radius); // use negative max_corr_radius for squares instead of circles?
		final int max_search_radius_poly = 1;
		if (globalDebugLevel > 0){
			System.out.println("max_corr_radius=       "+clt_parameters.max_corr_radius);
			System.out.println("max_search_radius=     "+max_search_radius);
			System.out.println("max_search_radius_poly="+max_search_radius_poly);
			System.out.println("corr_fat_zero=         "+fatzero);
			System.out.println("disparity_array[0][0]= "+disparity_array[0][0]);


		}
		// add optional initialization of debug layers here
		if (disparity_bimap != null){
			for (int i = 0; i < disparity_bimap.length;i++){
				disparity_bimap[i] = new double [tilesY*tilesX];
			}
		}

		if (ers_delay != null) {
			ers_delay[0] = new double [quad_main][];
			for (int i = 0; i < quad_main; i++) ers_delay[0][i] = new double [tilesX*tilesY];
			ers_delay[1] = new double [quad_aux][];
			for (int i = 0; i < quad_aux; i++)  ers_delay[1][i] = new double [tilesX*tilesY];
		}

		dtt.set_window(clt_parameters.clt_window);
		final double [] lt_window = dtt.getWin2d();	// [256]
		final double [] lt_window2 = new double [lt_window.length]; // squared
		for (int i = 0; i < lt_window.length; i++) lt_window2[i] = lt_window[i] * lt_window[i];


		if (globalDebugLevel > 1) {
			ShowDoubleFloatArrays sdfa_instance = new ShowDoubleFloatArrays(); // just for debugging?
			sdfa_instance.showArrays(lt_window,  2*clt_parameters.transform_size, 2*clt_parameters.transform_size, "lt_window");
		}

		final Matrix [] corr_rots_main = geometryCorrection_main.getCorrVector().getRotMatrices(); // get array of per-sensor rotation matrices
		final Matrix rigMatrix = geometryCorrection_aux.getRotMatrix(true);
		final Matrix [] corr_rots_aux =  geometryCorrection_aux.getCorrVector().getRotMatrices(rigMatrix); // get array of per-sensor rotation matrices

		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				@Override
				public void run() {
					DttRad2 dtt = new DttRad2(clt_parameters.transform_size);
					dtt.set_window(clt_parameters.clt_window);
					int tileY,tileX,tIndex; // , chn;
					//						showDoubleFloatArrays sdfa_instance = new showDoubleFloatArrays(); // just for debugging?
					double centerX; // center of aberration-corrected (common model) tile, X
					double centerY; //
					double [][] fract_shiftsXY_main = new double[quad_main][];
					double [][] fract_shiftsXY_aux =  new double[quad_aux][];
					double [][]     tcorr_combo =     null; // [15*15] pixel space
					double [][][][] clt_data_main =   new double[quad_main][nChn][][];
					double [][][][] clt_data_aux =    new double[quad_aux][nChn][][];
					double [][] ml_data_main =  (ml_data != null)? new double [ML_TOP_AUX_INDEX][(2*ml_hwidth +1)*(2*ml_hwidth +1)]:null;
					double [][] ml_data_aux =   (ml_data != null)? new double [ML_TOP_AUX_INDEX][(2*ml_hwidth +1)*(2*ml_hwidth +1)]:null;
					double []   ml_data_inter = (ml_data != null)? new double [(2*ml_hwidth +1)*(2*ml_hwidth +1)]:null;
//					double []   ml_data_other = (ml_data != null)? new double [(2*ml_hwidth +1)*(2*ml_hwidth +1)]:null;
					double []   ml_data_dbg1 =  (ml_data != null)? new double [(2*ml_hwidth +1)*(2*ml_hwidth +1)]:null;


					Correlation2d corr2d = new Correlation2d(
							clt_parameters.img_dtt,              // ImageDttParameters  imgdtt_params,
							clt_parameters.transform_size,             // int transform_size,
							2.0,                        //  double wndx_scale, // (wndy scale is always 1.0)
							isMonochrome(), // boolean monochrome,
							(globalDebugLevel > -1));   //   boolean debug)

					for (int nTile = ai.getAndIncrement(); nTile < nTilesInChn; nTile = ai.getAndIncrement()) {

						tileY = nTile /tilesX;
						tileX = nTile % tilesX;
						tIndex = tileY * tilesX + tileX;
						if (tile_op[tileY][tileX] == 0) {
							if (disparity_bimap != null){
								disparity_bimap[BI_TARGET_INDEX][tIndex] = Double.NaN;
							}
							 continue; // nothing to do for this tile
						}
						int                 img_mask =       getImgMask(tile_op[tileY][tileX]);         // which images to use
///						int                 corr_mask =      getPairMask(tile_op[tileY][tileX]);       // which pairs to combine in the combo:  1 - top, 2 bottom, 4 - left, 8 - right
						// mask out pairs that use missing channels

// Is it currently used with diagonals?
						// TODO: use masks from tile task
///						for (int i = 0; i< corr_pairs.length; i++){
///							if ((((1 << corr_pairs[i][0]) & img_mask) == 0) || (((1 << corr_pairs[i][1]) & img_mask) == 0)) {
///								corr_mask &= ~ (1 << i);
///							}
///						}
//						boolean debugTile =(tileX == debug_tileX) && (tileY == debug_tileY);

						final int [] overexp_main = (saturation_main != null) ? ( new int [2]): null;
						final int [] overexp_aux =  (saturation_aux != null) ? ( new int [2]): null;

						// Moved from inside chn loop
						centerX = tileX * clt_parameters.transform_size + clt_parameters.transform_size/2; //  - shiftX;
						centerY = tileY * clt_parameters.transform_size + clt_parameters.transform_size/2; //  - shiftY;
						// TODO: move port coordinates out of color channel loop
						double [][] centersXY_main;
						double [][] centersXY_aux;
						double disparity_main = disparity_array[tileY][tileX];
						double disparity_aux =  disparity_main * geometryCorrection_aux.getDisparityRadius()/geometryCorrection_main.getDisparityRadius();
						if (disparity_bimap != null){
							disparity_bimap[BI_TARGET_INDEX][tIndex] = disparity_main;
						}
						centersXY_main = geometryCorrection_main.getPortsCoordinatesAndDerivatives(
								geometryCorrection_main, //			GeometryCorrection gc_main,
								false,          // boolean use_rig_offsets,
								corr_rots_main, // Matrix []   rots,
								null,           //  Matrix [][] deriv_rots,
								null,           // double [][] pXYderiv, // if not null, should be double[8][]
								centerX,
								centerY,
								disparity_main); //  + disparity_corr);

						centersXY_aux =  geometryCorrection_aux.getPortsCoordinatesAndDerivatives(
								geometryCorrection_main, //			GeometryCorrection gc_main,
								true,            // boolean use_rig_offsets,
								corr_rots_aux,   // Matrix []   rots,
								null,            //  Matrix [][] deriv_rots,
								null,            // double [][] pXYderiv, // if not null, should be double[8][]
								centerX,
								centerY,
								disparity_aux); //  + disparity_corr);
						// acquisition time of the tiles centers in scanline times
						if (ers_delay != null) {
							for (int i = 0; i < quad_main; i++) ers_delay[0][i][nTile] = centersXY_main[i][1]-geometryCorrection_main.woi_tops[i];
							for (int i = 0; i < quad_aux; i++)  ers_delay[1][i][nTile] = centersXY_aux[i][1]- geometryCorrection_aux.woi_tops[i];
						}

						if ((globalDebugLevel > 0) && (tileX == debug_tileX) && (tileY == debug_tileY)) {
							for (int i = 0; i < quad_main; i++) {
								System.out.println("clt_aberrations_quad_corr():  tileX="+tileX+", tileY="+tileY+
										" centerX="+centerX+" centerY="+centerY+" disparity="+disparity_array[tileY][tileX]+
										" centersXY_main["+i+"][0]="+centersXY_main[i][0]+" centersXY_main["+i+"][1]="+centersXY_main[i][1]);
							}
							for (int i = 0; i < quad_aux; i++) {
								System.out.println("clt_aberrations_quad_corr():  tileX="+tileX+", tileY="+tileY+
										" centerX="+centerX+" centerY="+centerY+" disparity="+disparity_array[tileY][tileX]+
										" centersXY_aux["+i+"][0]="+centersXY_aux[i][0]+" centersXY_aux["+i+"][1]="+centersXY_aux[i][1]);
							}
						}
						if (port_xy_main_dbg != null) {
							port_xy_main_dbg[nTile] = centersXY_main;
						}
						if (port_xy_aux_dbg != null) {
							port_xy_aux_dbg [nTile] = centersXY_main;
						}
						if ((globalDebugLevel > -1) && (tileX == debug_tileX) && (tileY == debug_tileY)) { // before correction
							System.out.print(disparity_array[tileY][tileX]+"\t"+
									centersXY_main[0][0]+"\t"+centersXY_main[0][1]+"\t"+
									centersXY_main[1][0]+"\t"+centersXY_main[1][1]+"\t"+
									centersXY_main[2][0]+"\t"+centersXY_main[2][1]+"\t"+
									centersXY_main[3][0]+"\t"+centersXY_main[3][1]+"\t");
							System.out.print(disparity_array[tileY][tileX]+"\t"+
									centersXY_aux[0][0]+"\t"+centersXY_aux[0][1]+"\t"+
									centersXY_aux[1][0]+"\t"+centersXY_aux[1][1]+"\t"+
									centersXY_aux[2][0]+"\t"+centersXY_aux[2][1]+"\t"+
									centersXY_aux[3][0]+"\t"+centersXY_aux[3][1]+"\t");
						}

						for (int chn = 0; chn <numcol; chn++) {
							if ((globalDebugLevel > -1) && (tileX == debug_tileX) && (tileY == debug_tileY) && (chn == 2)) {
								System.out.println("\nMain camera, centerX="+centerX+", centerY="+centerY);
								System.out.println(disparity_array[tileY][tileX]+"\t"+
							    centersXY_main[0][0]+"\t"+centersXY_main[0][1]+"\t"+
							    centersXY_main[1][0]+"\t"+centersXY_main[1][1]+"\t"+
							    centersXY_main[2][0]+"\t"+centersXY_main[2][1]+"\t"+
							    centersXY_main[3][0]+"\t"+centersXY_main[3][1]+"\t");
								System.out.println("\nAux camera, centerX="+centerX+", centerY="+centerY);
								System.out.println(disparity_array[tileY][tileX]+"\t"+
							    centersXY_aux[0][0]+"\t"+centersXY_aux[0][1]+"\t"+
							    centersXY_aux[1][0]+"\t"+centersXY_aux[1][1]+"\t"+
							    centersXY_aux[2][0]+"\t"+centersXY_aux[2][1]+"\t"+
							    centersXY_aux[3][0]+"\t"+centersXY_aux[3][1]+"\t");
							}

							for (int i = 0; i < quad_main; i++) {
								clt_data_main[i][chn] = new double [4][];
								fract_shiftsXY_main[i] = extract_correct_tile( // return a pair of residual offsets
										image_data_main[i],
										width,       // image width
										(clt_kernels_main == null) ? null : clt_kernels_main[i], // [color][tileY][tileX][band][pixel]
										clt_data_main[i][chn], //double  [][]        clt_tile,    // should be double [4][];
										clt_parameters.kernel_step,
										clt_parameters.transform_size,
										dtt,
										chn,
										centersXY_main[i][0], // centerX, // center of aberration-corrected (common model) tile, X
										centersXY_main[i][1], // centerY, //
//										0, // (globalDebugLevel > -2) && (tileX == debug_tileX) && (tileY == debug_tileY)? 2:0, // external tile compare
										(globalDebugLevel > -2) && (tileX == debug_tileX) && (tileY == debug_tileY)? 2:0, // external tile compare
//										(globalDebugLevel > -2) && (tileX == debug_tileX) && (tileY == debug_tileY) && (i == 0)? 2:0, // external tile compare
										false,// no_deconvolution,
										false, // ); // transpose);
										((saturation_main != null) ? saturation_main[i] : null), //final boolean [][]        saturation_imp, // (near) saturated pixels or null
										((saturation_main != null) ? overexp_main: null)); // final double [] overexposed)


							}
							for (int i = 0; i < quad_aux; i++) {
								clt_data_aux[i][chn] = new double [4][];
								fract_shiftsXY_aux[i] = extract_correct_tile( // return a pair of residual offsets
										image_data_aux[i],
										width,       // image width
										(clt_kernels_aux == null) ? null : clt_kernels_aux[i], // [color][tileY][tileX][band][pixel]
										clt_data_aux[i][chn], //double  [][]        clt_tile,    // should be double [4][];
										clt_parameters.kernel_step,
										clt_parameters.transform_size,
										dtt,
										chn,
										centersXY_aux[i][0], // centerX, // center of aberration-corrected (common model) tile, X
										centersXY_aux[i][1], // centerY, //
										 0, // external tile compare
										false,// no_deconvolution,
										false, // ); // transpose);
										((saturation_aux != null) ? saturation_aux[i] : null), //final boolean [][]        saturation_imp, // (near) saturated pixels or null
										((saturation_aux != null) ? overexp_aux: null)); // final double [] overexposed)


							}
							if ((globalDebugLevel > -1) && (tileX == debug_tileX) && (tileY == debug_tileY) && (chn == 2)) {
								System.out.println();
							}
							if ((globalDebugLevel > 0) && (debug_tileX == tileX) && (debug_tileY == tileY)  && (chn == 2)) {
								ShowDoubleFloatArrays sdfa_instance = new ShowDoubleFloatArrays(); // just for debugging?
								String [] titles = {"CC0","SC0","CS0","SS0","CC1","SC1","CS1","SS1","CC2","SC2","CS2","SS2","CC3","SC3","CS3","SS3"};
								double [][] dbg_tile = new double [16][];
								for (int i = 0; i < 16; i++) dbg_tile[i]=clt_data_main[i>>2][chn][i & 3];
								sdfa_instance.showArrays(dbg_tile,  clt_parameters.transform_size, clt_parameters.transform_size, true, "MAIN_pre-shifted_x"+tileX+"_y"+tileY, titles);

								double [][] dbg_tile_aux = new double [16][];
								for (int i = 0; i < 16; i++) dbg_tile[i]=clt_data_aux[i>>2][chn][i & 3];
								sdfa_instance.showArrays(dbg_tile_aux,  clt_parameters.transform_size, clt_parameters.transform_size, true, "AUX_pre-shifted_x"+tileX+"_y"+tileY, titles);
							}

							if ((globalDebugLevel > 0) && (tileX >= debug_tileX - 2) && (tileX <= debug_tileX + 2) &&
									(tileY >= debug_tileY - 2) && (tileY <= debug_tileY+2)) {
								for (int i = 0; i < quad_main; i++) {
									System.out.println("clt_bi_quad(): color="+chn+", tileX="+tileX+", tileY="+tileY+
											" fract_shiftsXY_main["+i+"][0]="+fract_shiftsXY_main[i][0]+" fract_shiftsXY_main["+i+"][1]="+fract_shiftsXY_main[i][1]);
								}
								for (int i = 0; i < quad_aux; i++) {
									System.out.println("clt_bi_quad(): color="+chn+", tileX="+tileX+", tileY="+tileY+
											" fract_shiftsXY_aux["+i+"][0]="+fract_shiftsXY_aux[i][0]+" fract_shiftsXY_aux["+i+"][1]="+fract_shiftsXY_aux[i][1]);
								}
							}

							// apply residual shift
							boolean debug_gpu = (globalDebugLevel > -2) && (tileX == debug_tileX) && (tileY == debug_tileY);
							for (int i = 0; i < quad_main; i++) {
								fract_shift(    // fractional shift in transform domain. Currently uses sin/cos - change to tables with 2? rotations
										clt_data_main[i][chn], // double  [][]  clt_tile,
										clt_parameters.transform_size,
										fract_shiftsXY_main[i][0],            // double        shiftX,
										fract_shiftsXY_main[i][1],            // double        shiftY,
//										((globalDebugLevel > 1) && (chn==0) && (tileX >= debug_tileX - 2) && (tileX <= debug_tileX + 2) &&
//												(tileY >= debug_tileY - 2) && (tileY <= debug_tileY+2))
										false //debug_gpu
										);
								// (globalDebugLevel > -2) && (tileX == debug_tileX) && (tileY == debug_tileY)? 2:0, // external tile compare
								if (debug_gpu) {
									System.out.println("---Shifted image tile for quad="+i+" color="+chn+", shift_hor = "+fract_shiftsXY_main[i][0]+", shift_vert = "+fract_shiftsXY_main[i][1]+"---");
									for (int dct_mode=0; dct_mode < 4; dct_mode++ ) {
										System.out.println("dct_mode="+dct_mode);
										for (int irow = 0; irow < clt_parameters.transform_size; irow++) {
											for (int jcol = 0; jcol < clt_parameters.transform_size; jcol++) {
												System.out.print(String.format("%10.5f ", clt_data_main[i][chn][dct_mode][clt_parameters.transform_size * irow + jcol]));
											}
											System.out.println();
										}
										System.out.println();
									}
								}
							}

							for (int i = 0; i < quad_aux; i++) {
								fract_shift(    // fractional shift in transform domain. Currently uses sin/cos - change to tables with 2? rotations
										clt_data_aux[i][chn], // double  [][]  clt_tile,
										clt_parameters.transform_size,
										fract_shiftsXY_aux[i][0],            // double        shiftX,
										fract_shiftsXY_aux[i][1],            // double        shiftY,
										((globalDebugLevel > 1) && (chn==0) && (tileX >= debug_tileX - 2) && (tileX <= debug_tileX + 2) &&
												(tileY >= debug_tileY - 2) && (tileY <= debug_tileY+2)));
							}



							if ((globalDebugLevel > 0) && (debug_tileX == tileX) && (debug_tileY == tileY)) {
								ShowDoubleFloatArrays sdfa_instance = new ShowDoubleFloatArrays(); // just for debugging?
								String [] titles = {"CC0","SC0","CS0","SS0","CC1","SC1","CS1","SS1","CC2","SC2","CS2","SS2","CC3","SC3","CS3","SS3"};
								double [][] dbg_tile = new double [16][];
								for (int i = 0; i < 16; i++) dbg_tile[i]=clt_data_main[i>>2][chn][i & 3];
								sdfa_instance.showArrays(dbg_tile,  clt_parameters.transform_size, clt_parameters.transform_size, true, "MAIN_shifted_x"+tileX+"_y"+tileY+"-z", titles);
								double [][] dbg_tile_aux = new double [16][];
								for (int i = 0; i < 16; i++) dbg_tile[i]=clt_data_aux[i>>2][chn][i & 3];
								sdfa_instance.showArrays(dbg_tile_aux,  clt_parameters.transform_size, clt_parameters.transform_size, true, "AUX_shifted_x"+tileX+"_y"+tileY+"-z", titles);
							}

						}

						int tile_lma_debug_level =  ((tileX == debug_tileX) && (tileY == debug_tileY))? clt_parameters.img_dtt.lma_debug_level : -1;

						// all color channels are done here
						double extra_disparity_main = 0.0; // used for textures:  if allowed, shift images extra before trying to combine
						double extra_disparity_aux = 0.0; // used for textures:  if allowed, shift images extra before trying to combine

						// fill clt_corr_combo if it exists
						if (disparity_bimap != null){ // not null - calculate correlations

							// calculate overexposed fraction - remove ?
//							if (saturation_imp != null){
//								disparity_map[OVEREXPOSED][nTile] = (1.0 * overexp_all[0]) / overexp_all[1];
//							}


							double [] tile_corrs_main = tileCorrs(
									clt_parameters,        // final EyesisCorrectionParameters.CLTParameters       clt_parameters,
									fatzero,               // final double          fatzero,         // May use correlation fat zero from 2 different parameters - fat_zero and rig.ml_fatzero
									true,                  // final boolean         get4dirs, // calculate disparity/strength for each of the 4 directions
									corr2d,                // final Correlation2d  corr2d,
									clt_data_main,         // final double [][][][] clt_data,
									filter,                // final double []       filter,
									col_weights,           // final double []        col_weights,
									ml_hwidth,             // final int             ml_hwidth,
						    		ml_data_main,          // final double [][]     ml_center_corr,
									tileX,                 // final int              tileX, // only used in debug output
									tileY,                 // final int              tileY,
									tile_lma_debug_level); // final int              debugLevel)

							double [] tile_corrs_aux  = tileCorrs(
									clt_parameters,        // final EyesisCorrectionParameters.CLTParameters       clt_parameters,
									fatzero,               // final double          fatzero,         // May use correlation fat zero from 2 different parameters - fat_zero and rig.ml_fatzero
									true,                  // final boolean         get4dirs, // calculate disparity/strength for each of the 4 directions
									corr2d,                // final Correlation2d  corr2d,
									clt_data_aux,          // final double [][][][] clt_data,
									filter,                // final double []       filter,
									col_weights,           // final double []        col_weights,
									ml_hwidth,             // final int             ml_hwidth,
						    		ml_data_aux,          // final double [][]     ml_center_corr,
									tileX,                 // final int              tileX, // only used in debug output
									tileY,                 // final int              tileY,
									tile_lma_debug_level); // final int              debugLevel)

							extra_disparity_main = tile_corrs_main[DISP_FULL_INDEX];
							extra_disparity_aux =  tile_corrs_aux[DISP_FULL_INDEX];
							if (Double.isNaN(extra_disparity_main)) extra_disparity_main = 0;
							if (Double.isNaN(extra_disparity_aux))  extra_disparity_aux = 0;

							for (int i = 0; i < tile_corrs_main.length; i++) {
								int dest = SNGL_TO_BI[0][i];
								if (disparity_bimap[dest] != null) disparity_bimap[dest][nTile] = tile_corrs_main[i];
							}
							for (int i = 0; i < tile_corrs_aux.length; i++) {
								int dest = SNGL_TO_BI[1][i];
								if (disparity_bimap[dest] != null) disparity_bimap[dest][nTile] = tile_corrs_aux[i];
							}
							if (Double.isNaN(disparity_bimap[BI_STR_FULL_INDEX][tIndex])) {
								System.out.println("BUG: 3. disparity_map[BI_STR_FULL_INDEX][tIndex] should not be NaN");
							}
							if (Double.isNaN(disparity_bimap[BI_ASTR_FULL_INDEX][tIndex])) {
								System.out.println("BUG: 3a. disparity_map[BI_ASTR_FULL_INDEX][tIndex] should not be NaN");
							}

							if (clt_corr_combo != null) { // [type][tilesY][tilesX][(2*transform_size-1)*(2*transform_size-1)] // if null - will not calculate
								tcorr_combo = new double [TCORR_TITLES.length][corr_size * corr_size];
							}
							double [] inter_cam_corr = corr2d.correlateInterCamerasFD(
									clt_data_main,       // double [][][][]     clt_data_tile_main,
									clt_data_aux,        // double [][][][]     clt_data_tile_aux,
						    		filter,                   // double []           lpf,
						    		scale_strengths,
						    		col_weights,              // double []           col_weights,
						    		fatzero);                 // double              fat_zero)

							double [] inter_corrs_dxy = tileInterCamCorrs(
									no_int_x0,             // final boolean                                   no_int_x0, // do not offset window to integer - used when averaging low textures to avoid "jumps" for very wide
									clt_parameters,        // final EyesisCorrectionParameters.CLTParameters  clt_parameters,
									inter_cam_corr,        // final double []                                 inter_cam_corr,
									corr2d,                // final Correlation2d                             corr2d,
									ml_hwidth,             // final int                                       ml_hwidth,
						    		ml_data_inter,         // final double []                                 ml_center_corr,
									tcorr_combo,           // double [][]                                     tcorr_combo,
									notch_mode,            // final boolean                                   notch_mode,
									tileX,                 // final int              tileX, // only used in debug output
									tileY,                 // final int              tileY,
									tile_lma_debug_level); // final int              debugLevel)

							if (lt_corr != null) {
								lt_corr[tIndex] = inter_cam_corr;
							}


							if (inter_corrs_dxy != null) {
								disparity_bimap[BI_DISP_CROSS_INDEX][nTile] =    inter_corrs_dxy[INDEX_DISP];
								disparity_bimap[BI_STR_CROSS_INDEX][nTile] =     inter_corrs_dxy[INDEX_STRENGTH];
								disparity_bimap[BI_DISP_CROSS_DX_INDEX][nTile] = inter_corrs_dxy[INDEX_DX];
								disparity_bimap[BI_DISP_CROSS_DY_INDEX][nTile] = inter_corrs_dxy[INDEX_DY];
								// TODO: Use strength for the same residual disparity
								disparity_bimap[BI_STR_ALL_INDEX][nTile] =
										Math.pow(inter_corrs_dxy[INDEX_STRENGTH]*
										disparity_bimap[BI_STR_FULL_INDEX][nTile]*
										disparity_bimap[BI_ASTR_FULL_INDEX][nTile], 1.0/3);
							}
							// finalize ML stuff
							if (ml_data != null) {
								// save data for the main camera
								for (int nlayer = 0; nlayer <  ML_TOP_AUX_INDEX; nlayer++) {
									// save main camera data
									corr2d.saveMlTile(
								    		tileX,                     // int         tileX,
								    		tileY,                     // int         tileY,
								    		ml_hwidth,                 // int         ml_hwidth,
								    		ml_data,                   // double [][] ml_data,
								    		nlayer + 0,                // int         ml_layer,
								    		ml_data_main[nlayer],      // double []   ml_tile,
								    		tilesX);                   // int         tilesX);
									// save aux_camera data
									corr2d.saveMlTile(
								    		tileX,                     // int         tileX,
								    		tileY,                     // int         tileY,
								    		ml_hwidth,                 // int         ml_hwidth,
								    		ml_data,                   // double [][] ml_data,
								    		nlayer + ML_TOP_AUX_INDEX, // int         ml_layer,
								    		ml_data_aux[nlayer],       // double []   ml_tile,
								    		tilesX);                   // int         tilesX);
								}
								// save inter-camera correlation
								corr2d.saveMlTile(
							    		tileX,                         // int         tileX,
							    		tileY,                         // int         tileY,
							    		ml_hwidth,                     // int         ml_hwidth,
							    		ml_data,                       // double [][] ml_data,
							    		ML_INTER_INDEX,                // int         ml_layer,
							    		ml_data_inter,                 // double []   ml_tile,
							    		tilesX);                       // int         tilesX);
								// save other data (just 1 value)
/*
								corr2d.saveMlTile(
							    		tileX,                         // int         tileX,
							    		tileY,                         // int         tileY,
							    		ml_hwidth,                     // int         ml_hwidth,
							    		ml_data,                       // double [][] ml_data,
							    		ML_OTHER_INDEX,                // int         ml_layer,
							    		ml_data_other,                 // double []   ml_tile,
							    		tilesX);                       // int         tilesX);
*/
								corr2d.saveMlTilePixel(
							    		tileX,                         // int         tileX,
							    		tileY,                         // int         tileY,
							    		ml_hwidth,                     // int         ml_hwidth,
							    		ml_data,                       // double [][] ml_data,
							    		ML_OTHER_INDEX,                // int         ml_layer,
							    		ML_OTHER_TARGET ,              // int         ml_index,
							    		disparity_main,                // target disparitydouble      ml_value,
							    		tilesX);                       // int         tilesX);

								if (ml_data_dbg1 != null) {
									tileInterCamCorrs(
											false,             // final boolean                                   no_int_x0, // do not offset window to integer - used when averaging low textures to avoid "jumps" for very wide
											clt_parameters,    // final EyesisCorrectionParameters.CLTParameters  clt_parameters,
											fatzero,           // final double fatzero,         // May use correlation fat zero from 2 different parameters - fat_zero and rig.ml_fatzero

											corr2d,            // final Correlation2d                             corr2d,
											clt_data_main,     // double [][][][]                                 clt_data_tile_main,
											clt_data_main,     // double [][][][]                                 clt_data_tile_aux,
											filter,            // final double []       							filter,
											col_weights,       // final double []       							col_weights,
											ml_hwidth,         // final int                                        ml_hwidth,
											ml_data_dbg1,      // final double []                                 ml_center_corr,
											null,              // double [][]                                     tcorr_combo,
											false,             // final boolean                                   notch_mode,
											tileX,                 // final int              tileX, // only used in debug output
											tileY,                 // final int              tileY,
											tile_lma_debug_level); // final int              debugLevel)
									corr2d.saveMlTile(
								    		tileX,                         // int         tileX,
								    		tileY,                         // int         tileY,
								    		ml_hwidth,                     // int         ml_hwidth,
								    		ml_data,                       // double [][] ml_data,
								    		ML_DBG1_INDEX,                 // int         ml_layer,
								    		ml_data_dbg1,                  // double []   ml_tile,
								    		tilesX);                       // int         tilesX);
								}

							}


						} // if (disparity_map != null){ // not null - calculate correlations

						if (tcorr_combo != null) { // [type][tilesY][tilesX][(2*transform_size-1)*(2*transform_size-1)] // if null - will not calculate
							for (int i = 0; i < tcorr_combo.length; i++) {
								clt_corr_combo[i][tileY][tileX] = tcorr_combo[i];
							}
						}

						if (texture_tiles_main !=null) {
							generateTextureTiles (
									clt_parameters,        // final EyesisCorrectionParameters.CLTParameters       clt_parameters,
									extra_disparity_main,  // final double           extra_disparity,
									quad_main,                  // final int              quad,      // number of subcameras
									numcol,                // final int              numcol, // number of colors
									img_mask, // _main,         // int                    img_mask,
									tile_op, // _main,          // final int [][]         tile_op,         // [tilesY][tilesX] - what to do - 0 - nothing for this tile
									clt_data_main,         // final double [][][][]  clt_data,
									texture_tiles_main,    // final double [][][][]  texture_tiles,   // [tilesY][tilesX]["RGBA".length()][];  null - will skip images combining
									filter,                // final double []        filter,
									lt_window2,            // final double []        lt_window2,
									port_offsets,          // final double[][]       port_offsets,
									col_weights,           // final double []        col_weights,
									dtt,                   // final DttRad2          dtt,
									tileX,                 // final int              tileX, // only used in debug output
									tileY,                 // final int              tileY,
									tile_lma_debug_level); // final int              debugLevel);
						}
						if (texture_tiles_aux !=null) {
							generateTextureTiles (
									clt_parameters,        // final EyesisCorrectionParameters.CLTParameters       clt_parameters,
									extra_disparity_aux,   // final double           extra_disparity,
									quad_aux,              // final int              quad,      // number of subcameras
									numcol,                // final int              numcol, // number of colors
									img_mask, // _aux,         // int                    img_mask,
									tile_op, // _aux,          // final int [][]         tile_op,         // [tilesY][tilesX] - what to do - 0 - nothing for this tile
									clt_data_aux,         // final double [][][][]  clt_data,
									texture_tiles_aux,    // final double [][][][]  texture_tiles,   // [tilesY][tilesX]["RGBA".length()][];  null - will skip images combining
									filter,                // final double []        filter,
									lt_window2,            // final double []        lt_window2,
									port_offsets,          // final double[][]       port_offsets,
									col_weights,           // final double []        col_weights,
									dtt,                   // final DttRad2          dtt,
									tileX,                 // final int              tileX, // only used in debug output
									tileY,                 // final int              tileY,
									tile_lma_debug_level); // final int              debugLevel);
						}
						// Save channel tiles to result
						if (clt_bidata != null) {
							for (int i = 0; i < quad_main; i++) for (int j = 0; j < numcol; j++){
								clt_bidata[0][i][j][tileY][tileX] = clt_data_main[i][j];
							}
							for (int i = 0; i < quad_aux; i++) for (int j = 0; j < numcol; j++){
								clt_bidata[1][i][j][tileY][tileX] = clt_data_aux[i][j];
							}
						}
					}
				}
			};
		}
		startAndJoin(threads);

// If it was low-texture mode, 	use lt_corr to average bi-quad inter-correlation between neighbor tiles and then calculate disparity/strength
		if (lt_corr != null) {
			//notch_mode
			// prepare weights for neighbors
			final int lt_rad_x = notch_mode? 0: lt_rad;
			final double [][] neib_weights = new double[lt_rad+1][lt_rad_x+1];
			for (int i = 0; i <= lt_rad; i++) {
				for (int j = 0; j <= lt_rad_x; j++) {
					neib_weights[i][j] = Math.cos(Math.PI * i /(2 * lt_rad + 1)) * Math.cos(Math.PI * j /(2 * lt_rad + 1)); // no need to normalize - it will need to skip empty tiles anyway
				}
			}
			//		final int corr_size = clt_parameters.transform_size * 2 -1;

			final TileNeibs tnImage  = new TileNeibs(tilesX, tilesY);
			ai.set(0);
			for (int ithread = 0; ithread < threads.length; ithread++) {
				threads[ithread] = new Thread() {
					@Override
					public void run() {
						double []    ml_data_inter = (ml_data != null)? new double [(2*ml_hwidth +1)*(2*ml_hwidth +1)]:null;
						double [][]  tcorr_combo =     null; // [15*15] pixel space

						Correlation2d corr2d = new Correlation2d(
								clt_parameters.img_dtt,              // ImageDttParameters  imgdtt_params,
								clt_parameters.transform_size,             // int transform_size,
								2.0,                        //  double wndx_scale, // (wndy scale is always 1.0)
								isMonochrome(), // boolean monochrome,
								(globalDebugLevel > -1));   //   boolean debug)
						for (int nTile = ai.getAndIncrement(); nTile < nTilesInChn; nTile = ai.getAndIncrement()) if (lt_corr[nTile] != null){ // center must be non-null (from tile_op)
							int tileX = nTile % tilesX;
							int tileY = nTile / tilesX;
							double [] tile_corrs = new double [corr_size * corr_size];
							double sw = 0.0;
							for (int dy = -lt_rad; dy <= lt_rad; dy++) {
								int ady = (dy > 0)? dy:(-dy);
								for (int dx = -lt_rad_x; dx <= lt_rad_x; dx++) {
									int nTile1 =  tnImage.getNeibIndex(nTile, dx, dy);
									if (nTile1 >= 0) {
										double [] ot_corr = lt_corr[nTile1];
										if (ot_corr != null) {
											int adx = (dx > 0)? dx:(-dx);
											double nw = neib_weights[ady][adx];
											for (int i = 0; i < tile_corrs.length; i++) {
												tile_corrs[i] += nw * ot_corr[i];
											}
											sw+=nw;
										}
									}
								}
							}
							if (sw > 0.0) { // with the current window should always be so, as te center tile is non-null
								double s = 1.0/sw;
								for (int i = 0; i < tile_corrs.length; i++) {
									tile_corrs[i] *= s;
								}
								if (clt_corr_combo != null) { // [type][tilesY][tilesX][(2*transform_size-1)*(2*transform_size-1)] // if null - will not calculate
									tcorr_combo = new double [TCORR_TITLES.length][corr_size * corr_size];
								}
								int tile_lma_debug_level =  ((tileX == debug_tileX) && (tileY == debug_tileY))? clt_parameters.img_dtt.lma_debug_level : -1;

	  							double [] inter_corrs_dxy = tileInterCamCorrs(
										true,                  // final boolean                                   no_int_x0, // do not offset window to integer - used when averaging low textures to avoid "jumps" for very wide
										clt_parameters,        // final EyesisCorrectionParameters.CLTParameters  clt_parameters,
										tile_corrs,            // final double []                                 inter_cam_corr,
										corr2d,                // final Correlation2d                             corr2d,
										ml_hwidth,             // final int                                       ml_hwidth,
							    		ml_data_inter,         // final double []                                 ml_center_corr,
										tcorr_combo,           // double [][]                                     tcorr_combo,
										notch_mode,            // final boolean                                   notch_mode,
										tileX,                 // final int              tileX, // only used in debug output
										tileY,                 // final int              tileY,
										tile_lma_debug_level); // final int              debugLevel)
								if (inter_corrs_dxy != null) {
									disparity_bimap[BI_DISP_CROSS_INDEX][nTile] =    inter_corrs_dxy[INDEX_DISP];
									disparity_bimap[BI_STR_CROSS_INDEX][nTile] =     inter_corrs_dxy[INDEX_STRENGTH];
									disparity_bimap[BI_DISP_CROSS_DX_INDEX][nTile] = inter_corrs_dxy[INDEX_DX];
									disparity_bimap[BI_DISP_CROSS_DY_INDEX][nTile] = inter_corrs_dxy[INDEX_DY];
									// TODO: Use strength for the same residual disparity
									disparity_bimap[BI_STR_ALL_INDEX][nTile] =
											Math.pow(inter_corrs_dxy[INDEX_STRENGTH]*
											disparity_bimap[BI_STR_FULL_INDEX][nTile]*
											disparity_bimap[BI_ASTR_FULL_INDEX][nTile], 1.0/3);
								}
								if (tcorr_combo != null) { // [type][tilesY][tilesX][(2*transform_size-1)*(2*transform_size-1)] // if null - will not calculate
									for (int i = 0; i < tcorr_combo.length; i++) {
										clt_corr_combo[i][tileY][tileX] = tcorr_combo[i];
									}
								}
								if (ml_data != null) {
									// save inter-camera correlation
									corr2d.saveMlTile(
											tileX,                         // int         tileX,
											tileY,                         // int         tileY,
											ml_hwidth,                     // int         ml_hwidth,
											ml_data,                       // double [][] ml_data,
											ML_INTER_INDEX,                // int         ml_layer,
											ml_data_inter,                 // double []   ml_tile,
											tilesX);                       // int         tilesX);
								}
								// TODO: save ml_data_inter, tcorr_combo
							}
						}
					}
				};
			}
			startAndJoin(threads);
		}
		return  clt_bidata;
	}

	public double [][][][][][][]  clt_bi_quad(
			final CLTParameters       clt_parameters,
			final double              fatzero,         // May use correlation fat zero from 2 different parameters - fat_zero and rig.ml_fatzero
			final boolean             notch_mode,      // use notch filter for inter-camera correlation to detect poles
			final int                 lt_rad,          // low texture mode - inter-correlation is averaged between the neighbors before argmax-ing, using (2*notch_mode+1)^2 square
			final boolean             no_int_x0,       // do not offset window to integer maximum - used when averaging low textures to avoid "jumps" for very wide
			final int [][]            tile_op,         // [tilesY][tilesX] - what to do - 0 - nothing for this tile
			final double [][]         disparity_array, // [tilesY][tilesX] - individual per-tile expected disparity
			final double [][][]       image_data_main, // first index - number of image in a quad
			final double [][][]       image_data_aux,  // first index - number of image in a quad
		    final boolean [][]        saturation_main, // (near) saturated pixels or null
		    final boolean [][]        saturation_aux,  // (near) saturated pixels or null
			 // correlation results - combo will be for the correlation between two quad cameras
			final double [][][][]     clt_corr_combo,  // [type][tilesY][tilesX][(2*transform_size-1)*(2*transform_size-1)] // if null - will not calculate
			                                           // [type][tilesY][tilesX] should be set by caller
													   // types: 0 - selected correlation (product+offset), 1 - sum
			final double [][]         disparity_bimap, // [23][tilesY][tilesX], only [6][] is needed on input or null - do not calculate
			                                           // last 2 - contrast, avg/ "geometric average)
			final double [][]         ml_data,         // data for ML - 18 layers - 4 center areas (3x3, 5x5,..) per camera-per direction, 1 - composite, and 1 with just 1 data (target disparity)
			final double [][][][]     texture_tiles_main, // [tilesY][tilesX]["RGBA".length()][];  null - will skip images combining
			final double [][][][]     texture_tiles_aux,  // [tilesY][tilesX]["RGBA".length()][];  null - will skip images combining
			final int                 width,           // may be not multiple of 8, same for the height

			final GeometryCorrection  geometryCorrection_main,
			final GeometryCorrection  geometryCorrection_aux,
			final double [][][][][][] clt_kernels_main, // [channel_in_quad][color][tileY][tileX][band][pixel] , size should match image (have 1 tile around)
			final double [][][][][][] clt_kernels_aux,  // [channel_in_quad][color][tileY][tileX][band][pixel] , size should match image (have 1 tile around)
			final double              corr_magic_scale, // still not understood coefficient that reduces reported disparity value.  Seems to be around 0.85
			final boolean             keep_clt_data,
//			final int [][]            woi_tops,
			final double [][][]       ers_delay,        // if not null - fill with tile center acquisition delay
			final int                 threadsMax,  // maximal number of threads to launch
			final int                 debugLevel)
	{
		final int                 globalDebugLevel = clt_parameters.rig.rig_mode_debug?debugLevel:-2;
		final int                 debug_tileX = clt_parameters.tileX;
		final int                 debug_tileY = clt_parameters.tileY;
		final int quad_main = image_data_main.length;   // number of subcameras
		final int quad_aux =  image_data_aux.length;   // number of subcameras
		final int numcol = 3; // number of colors
		final int nChn = image_data_main[0].length;
		final int height=image_data_main[0][0].length/width;
		final int tilesX=width/clt_parameters.transform_size;
		final int tilesY=height/clt_parameters.transform_size;
		final int nTilesInChn=tilesX*tilesY;
		// clt_data does not need to be for the whole image (no, it is used for textures)
		final double [][][][][][][] clt_bidata = (keep_clt_data)? (new double[2][][][][][][]):null;
		if (clt_bidata != null) {
			clt_bidata[0] = new double[quad_main][nChn][tilesY][tilesX][][];
			clt_bidata[1] = new double[quad_aux][nChn][tilesY][tilesX][][];
		}

		final double [][] lt_corr = (lt_rad > 0)? (new double [nTilesInChn][]):null; // will keep inter-camera combo correlation, later combined in a separate multi-thread run

		final Thread[] threads = newThreadArray(threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		final double [] col_weights= new double [numcol]; // colors are RBG
		col_weights[2] = 1.0/(1.0 + clt_parameters.corr_red + clt_parameters.corr_blue);    // green color
		col_weights[0] = clt_parameters.corr_red *  col_weights[2];
		col_weights[1] = clt_parameters.corr_blue * col_weights[2];
		final int corr_size = clt_parameters.transform_size * 2 -1;
		final int [][] transpose_indices = new int [corr_size*(corr_size-1)/2][2];
		int indx = 0;
		for (int i =0; i < corr_size-1; i++){
			for (int j = i+1; j < corr_size; j++){
				transpose_indices[indx  ][0] = i * corr_size + j;
				transpose_indices[indx++][1] = j * corr_size + i;
			}
		}

		// get ml_data half width
		final int ml_hwidth = (ml_data != null)?(((int) Math.round(Math.sqrt(ml_data[0].length/nTilesInChn)) - 1) / 2):0;

		// Create window  to select center correlation strip using
		// ortho_height - full width of non-zero elements
		// ortho_eff_height - effective height (ration of the weighted column sum to the center value)
		int wcenter = clt_parameters.transform_size - 1;
		final double [] ortho_weights = new double [corr_size]; // [15]
		for (int i = 0; i < corr_size; i++){
			if ((i >= wcenter - clt_parameters.img_dtt.ortho_height/2) && (i <= wcenter + clt_parameters.img_dtt.ortho_height/2)) {
				double dx = 1.0*(i-wcenter)/(clt_parameters.img_dtt.ortho_height/2 + 1);
				ortho_weights[i] = 0.5*(1.0+Math.cos(Math.PI*dx))/clt_parameters.img_dtt.ortho_eff_height;
			}
		}
		if (globalDebugLevel > 0){
			System.out.println("ortho_height="+ clt_parameters.img_dtt.ortho_height+" ortho_eff_height="+ clt_parameters.img_dtt.ortho_eff_height);
			for (int i = 0; i < corr_size; i++){
				System.out.println(" ortho_weights["+i+"]="+ ortho_weights[i]);
			}
		}

		if (globalDebugLevel > -2) {
			System.out.println("clt_aberrations_quad_corr(): width="+width+" height="+height+" transform_size="+clt_parameters.transform_size+
					" debug_tileX="+debug_tileX+" debug_tileY="+debug_tileY+" globalDebugLevel="+globalDebugLevel);
		}
///		final int [][] corr_pairs ={ // {first, second, rot} rot: 0 - as is, 1 - swap y,x
///				{0,1,0},
///				{2,3,0},
///				{0,2,1},
///				{1,3,1}};

		final double[][] port_offsets = {
				{-0.5, -0.5},
				{ 0.5, -0.5},
				{-0.5,  0.5},
				{ 0.5,  0.5}};
		final int transform_len = clt_parameters.transform_size * clt_parameters.transform_size;



		final double [] filter_direct= new double[transform_len];
		if (clt_parameters.getCorrSigma(isMonochrome()) == 0) {
			filter_direct[0] = 1.0;
			for (int i= 1; i<filter_direct.length;i++) filter_direct[i] =0;
		} else {
			for (int i = 0; i < clt_parameters.transform_size; i++){
				for (int j = 0; j < clt_parameters.transform_size; j++){
					filter_direct[i * clt_parameters.transform_size+j] = Math.exp(-(i*i+j*j)/(2*clt_parameters.getCorrSigma(isMonochrome()))); // FIXME: should be sigma*sigma !
				}
			}
		}
		// normalize
		double sum = 0;
		for (int i = 0; i < clt_parameters.transform_size; i++){
			for (int j = 0; j < clt_parameters.transform_size; j++){
				double d = 	filter_direct[i*clt_parameters.transform_size+j];
				d*=Math.cos(Math.PI*i/(2*clt_parameters.transform_size))*Math.cos(Math.PI*j/(2*clt_parameters.transform_size));
				if (i > 0) d*= 2.0;
				if (j > 0) d*= 2.0;
				sum +=d;
			}
		}
		for (int i = 0; i<filter_direct.length; i++){
			filter_direct[i] /= sum;
		}

		DttRad2 dtt = new DttRad2(clt_parameters.transform_size);
		final double [] filter= dtt.dttt_iiie(filter_direct);
		for (int i=0; i < filter.length;i++) filter[i] *= 2*clt_parameters.transform_size;

		// prepare disparity maps and weights
		final int max_search_radius = (int) Math.abs(clt_parameters.max_corr_radius); // use negative max_corr_radius for squares instead of circles?
		final int max_search_radius_poly = 1;
		if (globalDebugLevel > 0){
			System.out.println("max_corr_radius=       "+clt_parameters.max_corr_radius);
			System.out.println("max_search_radius=     "+max_search_radius);
			System.out.println("max_search_radius_poly="+max_search_radius_poly);
			System.out.println("corr_fat_zero=         "+fatzero);
			System.out.println("disparity_array[0][0]= "+disparity_array[0][0]);


		}
		// add optional initialization of debug layers here
		if (disparity_bimap != null){
			for (int i = 0; i < disparity_bimap.length;i++){
				disparity_bimap[i] = new double [tilesY*tilesX];
			}
		}

		if (ers_delay != null) {
			ers_delay[0] = new double [quad_main][];
			for (int i = 0; i < quad_main; i++) ers_delay[0][i] = new double [tilesX*tilesY];
			ers_delay[1] = new double [quad_aux][];
			for (int i = 0; i < quad_aux; i++)  ers_delay[1][i] = new double [tilesX*tilesY];
		}

		dtt.set_window(clt_parameters.clt_window);
		final double [] lt_window = dtt.getWin2d();	// [256]
		final double [] lt_window2 = new double [lt_window.length]; // squared
		for (int i = 0; i < lt_window.length; i++) lt_window2[i] = lt_window[i] * lt_window[i];


		if (globalDebugLevel > 1) {
			ShowDoubleFloatArrays sdfa_instance = new ShowDoubleFloatArrays(); // just for debugging?
			sdfa_instance.showArrays(lt_window,  2*clt_parameters.transform_size, 2*clt_parameters.transform_size, "lt_window");
		}

		final Matrix [] corr_rots_main = geometryCorrection_main.getCorrVector().getRotMatrices(); // get array of per-sensor rotation matrices
		final Matrix rigMatrix = geometryCorrection_aux.getRotMatrix(true);
		final Matrix [] corr_rots_aux =  geometryCorrection_aux.getCorrVector().getRotMatrices(rigMatrix); // get array of per-sensor rotation matrices

		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				@Override
				public void run() {
					DttRad2 dtt = new DttRad2(clt_parameters.transform_size);
					dtt.set_window(clt_parameters.clt_window);
					int tileY,tileX,tIndex; // , chn;
					//						showDoubleFloatArrays sdfa_instance = new showDoubleFloatArrays(); // just for debugging?
					double centerX; // center of aberration-corrected (common model) tile, X
					double centerY; //
					double [][] fract_shiftsXY_main = new double[quad_main][];
					double [][] fract_shiftsXY_aux =  new double[quad_aux][];
					double [][]     tcorr_combo =     null; // [15*15] pixel space
					double [][][][] clt_data_main =   new double[quad_main][nChn][][];
					double [][][][] clt_data_aux =    new double[quad_aux][nChn][][];
					double [][] ml_data_main =  (ml_data != null)? new double [ML_TOP_AUX_INDEX][(2*ml_hwidth +1)*(2*ml_hwidth +1)]:null;
					double [][] ml_data_aux =   (ml_data != null)? new double [ML_TOP_AUX_INDEX][(2*ml_hwidth +1)*(2*ml_hwidth +1)]:null;
					double []   ml_data_inter = (ml_data != null)? new double [(2*ml_hwidth +1)*(2*ml_hwidth +1)]:null;
//					double []   ml_data_other = (ml_data != null)? new double [(2*ml_hwidth +1)*(2*ml_hwidth +1)]:null;
					double []   ml_data_dbg1 =  (ml_data != null)? new double [(2*ml_hwidth +1)*(2*ml_hwidth +1)]:null;


					Correlation2d corr2d = new Correlation2d(
							clt_parameters.img_dtt,              // ImageDttParameters  imgdtt_params,
							clt_parameters.transform_size,             // int transform_size,
							2.0,                        //  double wndx_scale, // (wndy scale is always 1.0)
							isMonochrome(), // boolean monochrome,
							(globalDebugLevel > -1));   //   boolean debug)

					for (int nTile = ai.getAndIncrement(); nTile < nTilesInChn; nTile = ai.getAndIncrement()) {

						tileY = nTile /tilesX;
						tileX = nTile % tilesX;
						tIndex = tileY * tilesX + tileX;
						if (tile_op[tileY][tileX] == 0) {
							if (disparity_bimap != null){
								disparity_bimap[BI_TARGET_INDEX][tIndex] = Double.NaN;
							}
							 continue; // nothing to do for this tile
						}
						int                 img_mask =       getImgMask(tile_op[tileY][tileX]);         // which images to use
///						int                 corr_mask =      getPairMask(tile_op[tileY][tileX]);       // which pairs to combine in the combo:  1 - top, 2 bottom, 4 - left, 8 - right
						// mask out pairs that use missing channels

// Is it currently used with diagonals?
						// TODO: use masks from tile task
///						for (int i = 0; i< corr_pairs.length; i++){
///							if ((((1 << corr_pairs[i][0]) & img_mask) == 0) || (((1 << corr_pairs[i][1]) & img_mask) == 0)) {
///								corr_mask &= ~ (1 << i);
///							}
///						}
//						boolean debugTile =(tileX == debug_tileX) && (tileY == debug_tileY);

						final int [] overexp_main = (saturation_main != null) ? ( new int [2]): null;
						final int [] overexp_aux =  (saturation_aux != null) ? ( new int [2]): null;

						// Moved from inside chn loop
						centerX = tileX * clt_parameters.transform_size + clt_parameters.transform_size/2; //  - shiftX;
						centerY = tileY * clt_parameters.transform_size + clt_parameters.transform_size/2; //  - shiftY;
						// TODO: move port coordinates out of color channel loop
						double [][] centersXY_main;
						double [][] centersXY_aux;
						double disparity_main = disparity_array[tileY][tileX];
						double disparity_aux =  disparity_main * geometryCorrection_aux.getDisparityRadius()/geometryCorrection_main.getDisparityRadius();
						if (disparity_bimap != null){
							disparity_bimap[BI_TARGET_INDEX][tIndex] = disparity_main;
						}
						centersXY_main = geometryCorrection_main.getPortsCoordinatesAndDerivatives(
								geometryCorrection_main, //			GeometryCorrection gc_main,
								false,          // boolean use_rig_offsets,
								corr_rots_main, // Matrix []   rots,
								null,           //  Matrix [][] deriv_rots,
								null,           // double [][] pXYderiv, // if not null, should be double[8][]
								centerX,
								centerY,
								disparity_main); //  + disparity_corr);

						centersXY_aux =  geometryCorrection_aux.getPortsCoordinatesAndDerivatives(
								geometryCorrection_main, //			GeometryCorrection gc_main,
								true,            // boolean use_rig_offsets,
								corr_rots_aux,   // Matrix []   rots,
								null,            //  Matrix [][] deriv_rots,
								null,            // double [][] pXYderiv, // if not null, should be double[8][]
								centerX,
								centerY,
								disparity_aux); //  + disparity_corr);
						// acquisition time of the tiles centers in scanline times
						if (ers_delay != null) {
							for (int i = 0; i < quad_main; i++) ers_delay[0][i][nTile] = centersXY_main[i][1]-geometryCorrection_main.woi_tops[i];
							for (int i = 0; i < quad_aux; i++)  ers_delay[1][i][nTile] = centersXY_aux[i][1]- geometryCorrection_aux.woi_tops[i];
						}

						if ((globalDebugLevel > 0) && (tileX == debug_tileX) && (tileY == debug_tileY)) {
							for (int i = 0; i < quad_main; i++) {
								System.out.println("clt_aberrations_quad_corr():  tileX="+tileX+", tileY="+tileY+
										" centerX="+centerX+" centerY="+centerY+" disparity="+disparity_array[tileY][tileX]+
										" centersXY_main["+i+"][0]="+centersXY_main[i][0]+" centersXY_main["+i+"][1]="+centersXY_main[i][1]);
							}
							for (int i = 0; i < quad_aux; i++) {
								System.out.println("clt_aberrations_quad_corr():  tileX="+tileX+", tileY="+tileY+
										" centerX="+centerX+" centerY="+centerY+" disparity="+disparity_array[tileY][tileX]+
										" centersXY_aux["+i+"][0]="+centersXY_aux[i][0]+" centersXY_aux["+i+"][1]="+centersXY_aux[i][1]);
							}
						}

						if ((globalDebugLevel > -1) && (tileX == debug_tileX) && (tileY == debug_tileY)) { // before correction
							System.out.print(disparity_array[tileY][tileX]+"\t"+
									centersXY_main[0][0]+"\t"+centersXY_main[0][1]+"\t"+
									centersXY_main[1][0]+"\t"+centersXY_main[1][1]+"\t"+
									centersXY_main[2][0]+"\t"+centersXY_main[2][1]+"\t"+
									centersXY_main[3][0]+"\t"+centersXY_main[3][1]+"\t");
							System.out.print(disparity_array[tileY][tileX]+"\t"+
									centersXY_aux[0][0]+"\t"+centersXY_aux[0][1]+"\t"+
									centersXY_aux[1][0]+"\t"+centersXY_aux[1][1]+"\t"+
									centersXY_aux[2][0]+"\t"+centersXY_aux[2][1]+"\t"+
									centersXY_aux[3][0]+"\t"+centersXY_aux[3][1]+"\t");
						}

						for (int chn = 0; chn <numcol; chn++) {
							if ((globalDebugLevel > -1) && (tileX == debug_tileX) && (tileY == debug_tileY) && (chn == 2)) {
								System.out.println("\nMain camera, centerX="+centerX+", centerY="+centerY);
								System.out.println(disparity_array[tileY][tileX]+"\t"+
							    centersXY_main[0][0]+"\t"+centersXY_main[0][1]+"\t"+
							    centersXY_main[1][0]+"\t"+centersXY_main[1][1]+"\t"+
							    centersXY_main[2][0]+"\t"+centersXY_main[2][1]+"\t"+
							    centersXY_main[3][0]+"\t"+centersXY_main[3][1]+"\t");
								System.out.println("\nAux camera, centerX="+centerX+", centerY="+centerY);
								System.out.println(disparity_array[tileY][tileX]+"\t"+
							    centersXY_aux[0][0]+"\t"+centersXY_aux[0][1]+"\t"+
							    centersXY_aux[1][0]+"\t"+centersXY_aux[1][1]+"\t"+
							    centersXY_aux[2][0]+"\t"+centersXY_aux[2][1]+"\t"+
							    centersXY_aux[3][0]+"\t"+centersXY_aux[3][1]+"\t");
							}

							for (int i = 0; i < quad_main; i++) {
								clt_data_main[i][chn] = new double [4][];
								fract_shiftsXY_main[i] = extract_correct_tile( // return a pair of residual offsets
										image_data_main[i],
										width,       // image width
										(clt_kernels_main == null) ? null : clt_kernels_main[i], // [color][tileY][tileX][band][pixel]
										clt_data_main[i][chn], //double  [][]        clt_tile,    // should be double [4][];
										clt_parameters.kernel_step,
										clt_parameters.transform_size,
										dtt,
										chn,
										centersXY_main[i][0], // centerX, // center of aberration-corrected (common model) tile, X
										centersXY_main[i][1], // centerY, //
										(globalDebugLevel > -2) && (tileX == debug_tileX) && (tileY == debug_tileY)? 2:0, // external tile compare
										false,// no_deconvolution,
										false, // ); // transpose);
										((saturation_main != null) ? saturation_main[i] : null), //final boolean [][]        saturation_imp, // (near) saturated pixels or null
										((saturation_main != null) ? overexp_main: null)); // final double [] overexposed)


							}
							for (int i = 0; i < quad_aux; i++) {
								clt_data_aux[i][chn] = new double [4][];
								fract_shiftsXY_aux[i] = extract_correct_tile( // return a pair of residual offsets
										image_data_aux[i],
										width,       // image width
										(clt_kernels_aux == null) ? null : clt_kernels_aux[i], // [color][tileY][tileX][band][pixel]
										clt_data_aux[i][chn], //double  [][]        clt_tile,    // should be double [4][];
										clt_parameters.kernel_step,
										clt_parameters.transform_size,
										dtt,
										chn,
										centersXY_aux[i][0], // centerX, // center of aberration-corrected (common model) tile, X
										centersXY_aux[i][1], // centerY, //
										 0, // external tile compare
										false,// no_deconvolution,
										false, // ); // transpose);
										((saturation_aux != null) ? saturation_aux[i] : null), //final boolean [][]        saturation_imp, // (near) saturated pixels or null
										((saturation_aux != null) ? overexp_aux: null)); // final double [] overexposed)


							}
							if ((globalDebugLevel > -1) && (tileX == debug_tileX) && (tileY == debug_tileY) && (chn == 2)) {
								System.out.println();
							}
							if ((globalDebugLevel > 0) && (debug_tileX == tileX) && (debug_tileY == tileY)  && (chn == 2)) {
								ShowDoubleFloatArrays sdfa_instance = new ShowDoubleFloatArrays(); // just for debugging?
								String [] titles = {"CC0","SC0","CS0","SS0","CC1","SC1","CS1","SS1","CC2","SC2","CS2","SS2","CC3","SC3","CS3","SS3"};
								double [][] dbg_tile = new double [16][];
								for (int i = 0; i < 16; i++) dbg_tile[i]=clt_data_main[i>>2][chn][i & 3];
								sdfa_instance.showArrays(dbg_tile,  clt_parameters.transform_size, clt_parameters.transform_size, true, "MAIN_pre-shifted_x"+tileX+"_y"+tileY, titles);

								double [][] dbg_tile_aux = new double [16][];
								for (int i = 0; i < 16; i++) dbg_tile[i]=clt_data_aux[i>>2][chn][i & 3];
								sdfa_instance.showArrays(dbg_tile_aux,  clt_parameters.transform_size, clt_parameters.transform_size, true, "AUX_pre-shifted_x"+tileX+"_y"+tileY, titles);
							}

							if ((globalDebugLevel > 0) && (tileX >= debug_tileX - 2) && (tileX <= debug_tileX + 2) &&
									(tileY >= debug_tileY - 2) && (tileY <= debug_tileY+2)) {
								for (int i = 0; i < quad_main; i++) {
									System.out.println("clt_bi_quad(): color="+chn+", tileX="+tileX+", tileY="+tileY+
											" fract_shiftsXY_main["+i+"][0]="+fract_shiftsXY_main[i][0]+" fract_shiftsXY_main["+i+"][1]="+fract_shiftsXY_main[i][1]);
								}
								for (int i = 0; i < quad_aux; i++) {
									System.out.println("clt_bi_quad(): color="+chn+", tileX="+tileX+", tileY="+tileY+
											" fract_shiftsXY_aux["+i+"][0]="+fract_shiftsXY_aux[i][0]+" fract_shiftsXY_aux["+i+"][1]="+fract_shiftsXY_aux[i][1]);
								}
							}

							// apply residual shift
							for (int i = 0; i < quad_main; i++) {
								fract_shift(    // fractional shift in transform domain. Currently uses sin/cos - change to tables with 2? rotations
										clt_data_main[i][chn], // double  [][]  clt_tile,
										clt_parameters.transform_size,
										fract_shiftsXY_main[i][0],            // double        shiftX,
										fract_shiftsXY_main[i][1],            // double        shiftY,
										((globalDebugLevel > 1) && (chn==0) && (tileX >= debug_tileX - 2) && (tileX <= debug_tileX + 2) &&
												(tileY >= debug_tileY - 2) && (tileY <= debug_tileY+2)));
							}
							for (int i = 0; i < quad_aux; i++) {
								fract_shift(    // fractional shift in transform domain. Currently uses sin/cos - change to tables with 2? rotations
										clt_data_aux[i][chn], // double  [][]  clt_tile,
										clt_parameters.transform_size,
										fract_shiftsXY_aux[i][0],            // double        shiftX,
										fract_shiftsXY_aux[i][1],            // double        shiftY,
										((globalDebugLevel > 1) && (chn==0) && (tileX >= debug_tileX - 2) && (tileX <= debug_tileX + 2) &&
												(tileY >= debug_tileY - 2) && (tileY <= debug_tileY+2)));
							}



							if ((globalDebugLevel > 0) && (debug_tileX == tileX) && (debug_tileY == tileY)) {
								ShowDoubleFloatArrays sdfa_instance = new ShowDoubleFloatArrays(); // just for debugging?
								String [] titles = {"CC0","SC0","CS0","SS0","CC1","SC1","CS1","SS1","CC2","SC2","CS2","SS2","CC3","SC3","CS3","SS3"};
								double [][] dbg_tile = new double [16][];
								for (int i = 0; i < 16; i++) dbg_tile[i]=clt_data_main[i>>2][chn][i & 3];
								sdfa_instance.showArrays(dbg_tile,  clt_parameters.transform_size, clt_parameters.transform_size, true, "MAIN_shifted_x"+tileX+"_y"+tileY+"-z", titles);
								double [][] dbg_tile_aux = new double [16][];
								for (int i = 0; i < 16; i++) dbg_tile[i]=clt_data_aux[i>>2][chn][i & 3];
								sdfa_instance.showArrays(dbg_tile_aux,  clt_parameters.transform_size, clt_parameters.transform_size, true, "AUX_shifted_x"+tileX+"_y"+tileY+"-z", titles);
							}

						}

						int tile_lma_debug_level =  ((tileX == debug_tileX) && (tileY == debug_tileY))? clt_parameters.img_dtt.lma_debug_level : -1;

						// all color channels are done here
						double extra_disparity_main = 0.0; // used for textures:  if allowed, shift images extra before trying to combine
						double extra_disparity_aux = 0.0; // used for textures:  if allowed, shift images extra before trying to combine

						// fill clt_corr_combo if it exists
						if (disparity_bimap != null){ // not null - calculate correlations

							// calculate overexposed fraction - remove ?
//							if (saturation_imp != null){
//								disparity_map[OVEREXPOSED][nTile] = (1.0 * overexp_all[0]) / overexp_all[1];
//							}


							double [] tile_corrs_main = tileCorrs(
									clt_parameters,        // final EyesisCorrectionParameters.CLTParameters       clt_parameters,
									fatzero,               // final double          fatzero,         // May use correlation fat zero from 2 different parameters - fat_zero and rig.ml_fatzero
									true,                  // final boolean         get4dirs, // calculate disparity/strength for each of the 4 directions
									corr2d,                // final Correlation2d  corr2d,
									clt_data_main,         // final double [][][][] clt_data,
									filter,                // final double []       filter,
									col_weights,           // final double []        col_weights,
									ml_hwidth,             // final int             ml_hwidth,
						    		ml_data_main,          // final double [][]     ml_center_corr,
									tileX,                 // final int              tileX, // only used in debug output
									tileY,                 // final int              tileY,
									tile_lma_debug_level); // final int              debugLevel)

							double [] tile_corrs_aux  = tileCorrs(
									clt_parameters,        // final EyesisCorrectionParameters.CLTParameters       clt_parameters,
									fatzero,               // final double          fatzero,         // May use correlation fat zero from 2 different parameters - fat_zero and rig.ml_fatzero
									true,                  // final boolean         get4dirs, // calculate disparity/strength for each of the 4 directions
									corr2d,                // final Correlation2d  corr2d,
									clt_data_aux,          // final double [][][][] clt_data,
									filter,                // final double []       filter,
									col_weights,           // final double []        col_weights,
									ml_hwidth,             // final int             ml_hwidth,
						    		ml_data_aux,          // final double [][]     ml_center_corr,
									tileX,                 // final int              tileX, // only used in debug output
									tileY,                 // final int              tileY,
									tile_lma_debug_level); // final int              debugLevel)

							extra_disparity_main = tile_corrs_main[DISP_FULL_INDEX];
							extra_disparity_aux =  tile_corrs_aux[DISP_FULL_INDEX];
							if (Double.isNaN(extra_disparity_main)) extra_disparity_main = 0;
							if (Double.isNaN(extra_disparity_aux))  extra_disparity_aux = 0;

							for (int i = 0; i < tile_corrs_main.length; i++) {
								int dest = SNGL_TO_BI[0][i];
								if (disparity_bimap[dest] != null) disparity_bimap[dest][nTile] = tile_corrs_main[i];
							}
							for (int i = 0; i < tile_corrs_aux.length; i++) {
								int dest = SNGL_TO_BI[1][i];
								if (disparity_bimap[dest] != null) disparity_bimap[dest][nTile] = tile_corrs_aux[i];
							}
							if (Double.isNaN(disparity_bimap[BI_STR_FULL_INDEX][tIndex])) {
								System.out.println("BUG: 3. disparity_map[BI_STR_FULL_INDEX][tIndex] should not be NaN");
							}
							if (Double.isNaN(disparity_bimap[BI_ASTR_FULL_INDEX][tIndex])) {
								System.out.println("BUG: 3a. disparity_map[BI_ASTR_FULL_INDEX][tIndex] should not be NaN");
							}

							if (clt_corr_combo != null) { // [type][tilesY][tilesX][(2*transform_size-1)*(2*transform_size-1)] // if null - will not calculate
								tcorr_combo = new double [TCORR_TITLES.length][corr_size * corr_size];
							}
							double [] inter_cam_corr = corr2d.correlateInterCamerasFD(
									clt_data_main,       // double [][][][]     clt_data_tile_main,
									clt_data_aux,        // double [][][][]     clt_data_tile_aux,
						    		filter,                   // double []           lpf,
						    		scale_strengths,
						    		col_weights,              // double []           col_weights,
						    		fatzero);                 // double              fat_zero)

							double [] inter_corrs_dxy = tileInterCamCorrs(
									no_int_x0,             // final boolean                                   no_int_x0, // do not offset window to integer - used when averaging low textures to avoid "jumps" for very wide
									clt_parameters,        // final EyesisCorrectionParameters.CLTParameters  clt_parameters,
									inter_cam_corr,        // final double []                                 inter_cam_corr,
									corr2d,                // final Correlation2d                             corr2d,
									ml_hwidth,             // final int                                       ml_hwidth,
						    		ml_data_inter,         // final double []                                 ml_center_corr,
									tcorr_combo,           // double [][]                                     tcorr_combo,
									notch_mode,            // final boolean                                   notch_mode,
									tileX,                 // final int              tileX, // only used in debug output
									tileY,                 // final int              tileY,
									tile_lma_debug_level); // final int              debugLevel)

							if (lt_corr != null) {
								lt_corr[tIndex] = inter_cam_corr;
							}


							if (inter_corrs_dxy != null) {
								disparity_bimap[BI_DISP_CROSS_INDEX][nTile] =    inter_corrs_dxy[INDEX_DISP];
								disparity_bimap[BI_STR_CROSS_INDEX][nTile] =     inter_corrs_dxy[INDEX_STRENGTH];
								disparity_bimap[BI_DISP_CROSS_DX_INDEX][nTile] = inter_corrs_dxy[INDEX_DX];
								disparity_bimap[BI_DISP_CROSS_DY_INDEX][nTile] = inter_corrs_dxy[INDEX_DY];
								// TODO: Use strength for the same residual disparity
								disparity_bimap[BI_STR_ALL_INDEX][nTile] =
										Math.pow(inter_corrs_dxy[INDEX_STRENGTH]*
										disparity_bimap[BI_STR_FULL_INDEX][nTile]*
										disparity_bimap[BI_ASTR_FULL_INDEX][nTile], 1.0/3);
							}
							// finalize ML stuff
							if (ml_data != null) {
								// save data for the main camera
								for (int nlayer = 0; nlayer <  ML_TOP_AUX_INDEX; nlayer++) {
									// save main camera data
									corr2d.saveMlTile(
								    		tileX,                     // int         tileX,
								    		tileY,                     // int         tileY,
								    		ml_hwidth,                 // int         ml_hwidth,
								    		ml_data,                   // double [][] ml_data,
								    		nlayer + 0,                // int         ml_layer,
								    		ml_data_main[nlayer],      // double []   ml_tile,
								    		tilesX);                   // int         tilesX);
									// save aux_camera data
									corr2d.saveMlTile(
								    		tileX,                     // int         tileX,
								    		tileY,                     // int         tileY,
								    		ml_hwidth,                 // int         ml_hwidth,
								    		ml_data,                   // double [][] ml_data,
								    		nlayer + ML_TOP_AUX_INDEX, // int         ml_layer,
								    		ml_data_aux[nlayer],       // double []   ml_tile,
								    		tilesX);                   // int         tilesX);
								}
								// save inter-camera correlation
								corr2d.saveMlTile(
							    		tileX,                         // int         tileX,
							    		tileY,                         // int         tileY,
							    		ml_hwidth,                     // int         ml_hwidth,
							    		ml_data,                       // double [][] ml_data,
							    		ML_INTER_INDEX,                // int         ml_layer,
							    		ml_data_inter,                 // double []   ml_tile,
							    		tilesX);                       // int         tilesX);
								// save other data (just 1 value)
/*
								corr2d.saveMlTile(
							    		tileX,                         // int         tileX,
							    		tileY,                         // int         tileY,
							    		ml_hwidth,                     // int         ml_hwidth,
							    		ml_data,                       // double [][] ml_data,
							    		ML_OTHER_INDEX,                // int         ml_layer,
							    		ml_data_other,                 // double []   ml_tile,
							    		tilesX);                       // int         tilesX);
*/
								corr2d.saveMlTilePixel(
							    		tileX,                         // int         tileX,
							    		tileY,                         // int         tileY,
							    		ml_hwidth,                     // int         ml_hwidth,
							    		ml_data,                       // double [][] ml_data,
							    		ML_OTHER_INDEX,                // int         ml_layer,
							    		ML_OTHER_TARGET ,              // int         ml_index,
							    		disparity_main,                // target disparitydouble      ml_value,
							    		tilesX);                       // int         tilesX);

								if (ml_data_dbg1 != null) {
									tileInterCamCorrs(
											false,             // final boolean                                   no_int_x0, // do not offset window to integer - used when averaging low textures to avoid "jumps" for very wide
											clt_parameters,    // final EyesisCorrectionParameters.CLTParameters  clt_parameters,
											fatzero,           // final double fatzero,         // May use correlation fat zero from 2 different parameters - fat_zero and rig.ml_fatzero

											corr2d,            // final Correlation2d                             corr2d,
											clt_data_main,     // double [][][][]                                 clt_data_tile_main,
											clt_data_main,     // double [][][][]                                 clt_data_tile_aux,
											filter,            // final double []       							filter,
											col_weights,       // final double []       							col_weights,
											ml_hwidth,         // final int                                        ml_hwidth,
											ml_data_dbg1,      // final double []                                 ml_center_corr,
											null,              // double [][]                                     tcorr_combo,
											false,             // final boolean                                   notch_mode,
											tileX,                 // final int              tileX, // only used in debug output
											tileY,                 // final int              tileY,
											tile_lma_debug_level); // final int              debugLevel)
									corr2d.saveMlTile(
								    		tileX,                         // int         tileX,
								    		tileY,                         // int         tileY,
								    		ml_hwidth,                     // int         ml_hwidth,
								    		ml_data,                       // double [][] ml_data,
								    		ML_DBG1_INDEX,                 // int         ml_layer,
								    		ml_data_dbg1,                  // double []   ml_tile,
								    		tilesX);                       // int         tilesX);
								}

							}


						} // if (disparity_map != null){ // not null - calculate correlations

						if (tcorr_combo != null) { // [type][tilesY][tilesX][(2*transform_size-1)*(2*transform_size-1)] // if null - will not calculate
							for (int i = 0; i < tcorr_combo.length; i++) {
								clt_corr_combo[i][tileY][tileX] = tcorr_combo[i];
							}
						}

						if (texture_tiles_main !=null) {
							generateTextureTiles (
									clt_parameters,        // final EyesisCorrectionParameters.CLTParameters       clt_parameters,
									extra_disparity_main,  // final double           extra_disparity,
									quad_main,                  // final int              quad,      // number of subcameras
									numcol,                // final int              numcol, // number of colors
									img_mask, // _main,         // int                    img_mask,
									tile_op, // _main,          // final int [][]         tile_op,         // [tilesY][tilesX] - what to do - 0 - nothing for this tile
									clt_data_main,         // final double [][][][]  clt_data,
									texture_tiles_main,    // final double [][][][]  texture_tiles,   // [tilesY][tilesX]["RGBA".length()][];  null - will skip images combining
									filter,                // final double []        filter,
									lt_window2,            // final double []        lt_window2,
									port_offsets,          // final double[][]       port_offsets,
									col_weights,           // final double []        col_weights,
									dtt,                   // final DttRad2          dtt,
									tileX,                 // final int              tileX, // only used in debug output
									tileY,                 // final int              tileY,
									tile_lma_debug_level); // final int              debugLevel);
						}
						if (texture_tiles_aux !=null) {
							generateTextureTiles (
									clt_parameters,        // final EyesisCorrectionParameters.CLTParameters       clt_parameters,
									extra_disparity_aux,   // final double           extra_disparity,
									quad_aux,              // final int              quad,      // number of subcameras
									numcol,                // final int              numcol, // number of colors
									img_mask, // _aux,         // int                    img_mask,
									tile_op, // _aux,          // final int [][]         tile_op,         // [tilesY][tilesX] - what to do - 0 - nothing for this tile
									clt_data_aux,         // final double [][][][]  clt_data,
									texture_tiles_aux,    // final double [][][][]  texture_tiles,   // [tilesY][tilesX]["RGBA".length()][];  null - will skip images combining
									filter,                // final double []        filter,
									lt_window2,            // final double []        lt_window2,
									port_offsets,          // final double[][]       port_offsets,
									col_weights,           // final double []        col_weights,
									dtt,                   // final DttRad2          dtt,
									tileX,                 // final int              tileX, // only used in debug output
									tileY,                 // final int              tileY,
									tile_lma_debug_level); // final int              debugLevel);
						}
						// Save channel tiles to result
						if (clt_bidata != null) {
							for (int i = 0; i < quad_main; i++) for (int j = 0; j < numcol; j++){
								clt_bidata[0][i][j][tileY][tileX] = clt_data_main[i][j];
							}
							for (int i = 0; i < quad_aux; i++) for (int j = 0; j < numcol; j++){
								clt_bidata[1][i][j][tileY][tileX] = clt_data_aux[i][j];
							}
						}
					}
				}
			};
		}
		startAndJoin(threads);

// If it was low-texture mode, 	use lt_corr to average bi-quad inter-correlation between neighbor tiles and then calculate disparity/strength
		if (lt_corr != null) {
			//notch_mode
			// prepare weights for neighbors
			final int lt_rad_x = notch_mode? 0: lt_rad;
			final double [][] neib_weights = new double[lt_rad+1][lt_rad_x+1];
			for (int i = 0; i <= lt_rad; i++) {
				for (int j = 0; j <= lt_rad_x; j++) {
					neib_weights[i][j] = Math.cos(Math.PI * i /(2 * lt_rad + 1)) * Math.cos(Math.PI * j /(2 * lt_rad + 1)); // no need to normalize - it will need to skip empty tiles anyway
				}
			}
			//		final int corr_size = clt_parameters.transform_size * 2 -1;

			final TileNeibs tnImage  = new TileNeibs(tilesX, tilesY);
			ai.set(0);
			for (int ithread = 0; ithread < threads.length; ithread++) {
				threads[ithread] = new Thread() {
					@Override
					public void run() {
						double []    ml_data_inter = (ml_data != null)? new double [(2*ml_hwidth +1)*(2*ml_hwidth +1)]:null;
						double [][]  tcorr_combo =     null; // [15*15] pixel space

						Correlation2d corr2d = new Correlation2d(
								clt_parameters.img_dtt,              // ImageDttParameters  imgdtt_params,
								clt_parameters.transform_size,             // int transform_size,
								2.0,                        //  double wndx_scale, // (wndy scale is always 1.0)
								isMonochrome(), // boolean monochrome,
								(globalDebugLevel > -1));   //   boolean debug)
						for (int nTile = ai.getAndIncrement(); nTile < nTilesInChn; nTile = ai.getAndIncrement()) if (lt_corr[nTile] != null){ // center must be non-null (from tile_op)
							int tileX = nTile % tilesX;
							int tileY = nTile / tilesX;
							double [] tile_corrs = new double [corr_size * corr_size];
							double sw = 0.0;
							for (int dy = -lt_rad; dy <= lt_rad; dy++) {
								int ady = (dy > 0)? dy:(-dy);
								for (int dx = -lt_rad_x; dx <= lt_rad_x; dx++) {
									int nTile1 =  tnImage.getNeibIndex(nTile, dx, dy);
									if (nTile1 >= 0) {
										double [] ot_corr = lt_corr[nTile1];
										if (ot_corr != null) {
											int adx = (dx > 0)? dx:(-dx);
											double nw = neib_weights[ady][adx];
											for (int i = 0; i < tile_corrs.length; i++) {
												tile_corrs[i] += nw * ot_corr[i];
											}
											sw+=nw;
										}
									}
								}
							}
							if (sw > 0.0) { // with the current window should always be so, as te center tile is non-null
								double s = 1.0/sw;
								for (int i = 0; i < tile_corrs.length; i++) {
									tile_corrs[i] *= s;
								}
								if (clt_corr_combo != null) { // [type][tilesY][tilesX][(2*transform_size-1)*(2*transform_size-1)] // if null - will not calculate
									tcorr_combo = new double [TCORR_TITLES.length][corr_size * corr_size];
								}
								int tile_lma_debug_level =  ((tileX == debug_tileX) && (tileY == debug_tileY))? clt_parameters.img_dtt.lma_debug_level : -1;

	  							double [] inter_corrs_dxy = tileInterCamCorrs(
										true,                  // final boolean                                   no_int_x0, // do not offset window to integer - used when averaging low textures to avoid "jumps" for very wide
										clt_parameters,        // final EyesisCorrectionParameters.CLTParameters  clt_parameters,
										tile_corrs,            // final double []                                 inter_cam_corr,
										corr2d,                // final Correlation2d                             corr2d,
										ml_hwidth,             // final int                                       ml_hwidth,
							    		ml_data_inter,         // final double []                                 ml_center_corr,
										tcorr_combo,           // double [][]                                     tcorr_combo,
										notch_mode,            // final boolean                                   notch_mode,
										tileX,                 // final int              tileX, // only used in debug output
										tileY,                 // final int              tileY,
										tile_lma_debug_level); // final int              debugLevel)
								if (inter_corrs_dxy != null) {
									disparity_bimap[BI_DISP_CROSS_INDEX][nTile] =    inter_corrs_dxy[INDEX_DISP];
									disparity_bimap[BI_STR_CROSS_INDEX][nTile] =     inter_corrs_dxy[INDEX_STRENGTH];
									disparity_bimap[BI_DISP_CROSS_DX_INDEX][nTile] = inter_corrs_dxy[INDEX_DX];
									disparity_bimap[BI_DISP_CROSS_DY_INDEX][nTile] = inter_corrs_dxy[INDEX_DY];
									// TODO: Use strength for the same residual disparity
									disparity_bimap[BI_STR_ALL_INDEX][nTile] =
											Math.pow(inter_corrs_dxy[INDEX_STRENGTH]*
											disparity_bimap[BI_STR_FULL_INDEX][nTile]*
											disparity_bimap[BI_ASTR_FULL_INDEX][nTile], 1.0/3);
								}
								if (tcorr_combo != null) { // [type][tilesY][tilesX][(2*transform_size-1)*(2*transform_size-1)] // if null - will not calculate
									for (int i = 0; i < tcorr_combo.length; i++) {
										clt_corr_combo[i][tileY][tileX] = tcorr_combo[i];
									}
								}
								if (ml_data != null) {
									// save inter-camera correlation
									corr2d.saveMlTile(
											tileX,                         // int         tileX,
											tileY,                         // int         tileY,
											ml_hwidth,                     // int         ml_hwidth,
											ml_data,                       // double [][] ml_data,
											ML_INTER_INDEX,                // int         ml_layer,
											ml_data_inter,                 // double []   ml_tile,
											tilesX);                       // int         tilesX);
								}
								// TODO: save ml_data_inter, tcorr_combo
							}
						}
					}
				};
			}
			startAndJoin(threads);
		}
		return  clt_bidata;
	}

	public void  clt_bi_macro( // not yet operational
			final CLTParameters       clt_parameters,
			final double              fatzero,         // May use correlation fat zero from 2 different parameters - fat_zero and rig.ml_fatzero
			final int                 macro_scale,
			final int [][]            tile_op,         // [tilesY][tilesX] - what to do - 0 - nothing for this tile
			final double [][]         disparity_array, // [tilesY][tilesX] - individual per-tile expected disparity
			final double [][][]       image_rig_data,  // [2][1][pixels] (single color channel)
			final double [][]         disparity_bimap, // [23][tilesY][tilesX], only [6][] is needed on input or null - do not calculate
			final int                 width,
			final GeometryCorrection  geometryCorrection_main,
			final GeometryCorrection  geometryCorrection_aux, // it has rig offset)
			final double              corr_magic_scale, // still not understood coefficient that reduces reported disparity value.  Seems to be around 0.85
			final int                 threadsMax,  // maximal number of threads to launch
			final int                 debugLevel)
	{
		final int                 globalDebugLevel = clt_parameters.rig.rig_mode_debug?debugLevel:-2;
		final int                 debug_tileX = clt_parameters.tileX;
		final int                 debug_tileY = clt_parameters.tileY;
		final int height=         image_rig_data[0][0].length/width;
		final int tilesX=width/clt_parameters.transform_size;
		final int tilesY=height/clt_parameters.transform_size;
		final int nTilesInChn=tilesX*tilesY;

		// clt_data does not need to be for the whole image (no, it is used for textures)

		final Thread[] threads = newThreadArray(threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		final int corr_size = clt_parameters.transform_size * 2 -1;
		final int [][] transpose_indices = new int [corr_size*(corr_size-1)/2][2];
		int indx = 0;
		for (int i =0; i < corr_size-1; i++){
			for (int j = i+1; j < corr_size; j++){
				transpose_indices[indx  ][0] = i * corr_size + j;
				transpose_indices[indx++][1] = j * corr_size + i;
			}
		}


		// Create window  to select center correlation strip using
		// ortho_height - full width of non-zero elements
		// ortho_eff_height - effective height (ration of the weighted column sum to the center value)
		int wcenter = clt_parameters.transform_size - 1;
		final double [] ortho_weights = new double [corr_size]; // [15]
		for (int i = 0; i < corr_size; i++){
			if ((i >= wcenter - clt_parameters.img_dtt.ortho_height/2) && (i <= wcenter + clt_parameters.img_dtt.ortho_height/2)) {
				double dx = 1.0*(i-wcenter)/(clt_parameters.img_dtt.ortho_height/2 + 1);
				ortho_weights[i] = 0.5*(1.0+Math.cos(Math.PI*dx))/clt_parameters.img_dtt.ortho_eff_height;
			}
		}
		if (globalDebugLevel > 0){
			System.out.println("ortho_height="+ clt_parameters.img_dtt.ortho_height+" ortho_eff_height="+ clt_parameters.img_dtt.ortho_eff_height);
			for (int i = 0; i < corr_size; i++){
				System.out.println(" ortho_weights["+i+"]="+ ortho_weights[i]);
			}
		}

		if (globalDebugLevel > 0) {
			System.out.println("clt_aberrations_quad_corr(): width="+width+" height="+height+" transform_size="+clt_parameters.transform_size+
					" debug_tileX="+debug_tileX+" debug_tileY="+debug_tileY+" globalDebugLevel="+globalDebugLevel);
		}
		final int transform_len = clt_parameters.transform_size * clt_parameters.transform_size;



		final double [] filter_direct= new double[transform_len];
		if (clt_parameters.getCorrSigma(isMonochrome()) == 0) {
			filter_direct[0] = 1.0;
			for (int i= 1; i<filter_direct.length;i++) filter_direct[i] =0;
		} else {
			for (int i = 0; i < clt_parameters.transform_size; i++){
				for (int j = 0; j < clt_parameters.transform_size; j++){
					filter_direct[i * clt_parameters.transform_size+j] = Math.exp(-(i*i+j*j)/(2*clt_parameters.getCorrSigma(isMonochrome()))); // FIXME: should be sigma*sigma !
				}
			}
		}
		// normalize
		double sum = 0;
		for (int i = 0; i < clt_parameters.transform_size; i++){
			for (int j = 0; j < clt_parameters.transform_size; j++){
				double d = 	filter_direct[i*clt_parameters.transform_size+j];
				d*=Math.cos(Math.PI*i/(2*clt_parameters.transform_size))*Math.cos(Math.PI*j/(2*clt_parameters.transform_size));
				if (i > 0) d*= 2.0;
				if (j > 0) d*= 2.0;
				sum +=d;
			}
		}
		for (int i = 0; i<filter_direct.length; i++){
			filter_direct[i] /= sum;
		}

		DttRad2 dtt = new DttRad2(clt_parameters.transform_size);
		final double [] filter= dtt.dttt_iiie(filter_direct);
		for (int i=0; i < filter.length;i++) filter[i] *= 2*clt_parameters.transform_size;

		// prepare disparity maps and weights
		final int max_search_radius = (int) Math.abs(clt_parameters.max_corr_radius); // use negative max_corr_radius for squares instead of circles?
		final int max_search_radius_poly = 1;
		if (globalDebugLevel > 0){
			System.out.println("max_corr_radius=       "+clt_parameters.max_corr_radius);
			System.out.println("max_search_radius=     "+max_search_radius);
			System.out.println("max_search_radius_poly="+max_search_radius_poly);
			System.out.println("corr_fat_zero=         "+clt_parameters.getFatZero(isMonochrome()));
			System.out.println("disparity_array[0][0]= "+disparity_array[0][0]);


		}
		// add optional initialization of debug layers here
		// FIXME: init only relevant, leave others null?
		if (disparity_bimap != null){
			for (int i = 0; i < disparity_bimap.length;i++){
				disparity_bimap[i] = new double [tilesY*tilesX];
			}
		}

//		final double [] corr_max_weights_poly =(((clt_parameters.max_corr_sigma > 0) && (disparity_bimap != null))?
//				setMaxXYWeights(clt_parameters.max_corr_sigma,max_search_radius_poly): null); // here use square anyway

		dtt.set_window(clt_parameters.clt_window);
		final double [] lt_window = dtt.getWin2d();	// [256]
		final double [] lt_window2 = new double [lt_window.length]; // squared
		for (int i = 0; i < lt_window.length; i++) lt_window2[i] = lt_window[i] * lt_window[i];


		if (globalDebugLevel > 1) {
			ShowDoubleFloatArrays sdfa_instance = new ShowDoubleFloatArrays(); // just for debugging?
			sdfa_instance.showArrays(lt_window,  2*clt_parameters.transform_size, 2*clt_parameters.transform_size, "lt_window");
		}

		//FIXME: no rotation matrices, no rig disparity - it is all included in source data in macro mode. Distortion?

//		final Matrix rigMatrix = geometryCorrection_aux.getRotMatrix(true);
		final int nChn = image_rig_data[0].length; // 1; ==1
		final double [] col_weights= {1.0};
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				@Override
				public void run() {
					DttRad2 dtt = new DttRad2(clt_parameters.transform_size);
					dtt.set_window(clt_parameters.clt_window);
					int tileY,tileX,tIndex; // , chn;
					//						showDoubleFloatArrays sdfa_instance = new showDoubleFloatArrays(); // just for debugging?
					double centerX; // center of aberration-corrected (common model) tile, X
					double centerY; //
					double [] fract_shiftsXY_main = null; //new double[quad_main][];
					double [] fract_shiftsXY_aux =  null; // new double[quad_aux][];
					double [][][][] clt_data_main =   new double[1][nChn][][]; // single channel
					double [][][][] clt_data_aux =    new double[1][nChn][][];


					Correlation2d corr2d = new Correlation2d(
							clt_parameters.img_dtt,              // ImageDttParameters  imgdtt_params,
							clt_parameters.transform_size,             // int transform_size,
							2.0,                        //  double wndx_scale, // (wndy scale is always 1.0)
							isMonochrome(), // boolean monochrome,
							(globalDebugLevel > -1));   //   boolean debug)

					for (int nTile = ai.getAndIncrement(); nTile < nTilesInChn; nTile = ai.getAndIncrement()) {

						tileY = nTile /tilesX;
						tileX = nTile % tilesX;
						tIndex = tileY * tilesX + tileX;
						if (tile_op[tileY][tileX] == 0) {
							if (disparity_bimap != null){
								disparity_bimap[BI_TARGET_INDEX][tIndex] = Double.NaN;
							}
							 continue; // nothing to do for this tile
						}


						// Moved from inside chn loop
						centerX = tileX * clt_parameters.transform_size + clt_parameters.transform_size/2; //  - shiftX;
						centerY = tileY * clt_parameters.transform_size + clt_parameters.transform_size/2; //  - shiftY;
						// TODO: move port coordinates out of color channel loop
						double disparity = disparity_array[tileY][tileX];
						if (disparity_bimap != null){
							disparity_bimap[BI_TARGET_INDEX][tIndex] = disparity;
						}

						// here it needs disparity_aux, as the aux camera knows nothing about the main

						double [] centersXY_main = {centerX, centerY};
						double [] centersXY_aux = geometryCorrection_aux.getRigAuxCoordinatesIdeal(
								macro_scale,              // int  macro_scale, // 1 for pixels, 8 - for tiles when correlating tiles instead of the pixels
								geometryCorrection_main,  // GeometryCorrection gc_main,
								null, // rigMatrix,       // Matrix             aux_rot,
								centerX,                  // double             px,
								centerY,                  // double             py,
								disparity*macro_scale);   // double             disparity)


						if ((globalDebugLevel > -1) && (tileX == debug_tileX) && (tileY == debug_tileY)) { // before correction
							System.out.println(disparity_array[tileY][tileX]+"\t"+
									centersXY_main[0]+"\t"+centersXY_main[1]+"\t"+
									centersXY_aux[0]+" \t"+centersXY_aux[1]);
						}

						for (int chn = 0; chn < nChn; chn++) { // now single color channel in each (

							clt_data_main[0][chn] = new double [4][];
							fract_shiftsXY_main = extract_correct_tile( // return a pair of residual offsets
									image_rig_data[0],
									width,       // image width
									null,      // [color][tileY][tileX][band][pixel]
									clt_data_main[0][chn], //double  [][]        clt_tile,    // should be double [4][];
									clt_parameters.kernel_step,
									clt_parameters.transform_size,
									dtt,
									chn,
									centersXY_main[0], // centerX, // center of aberration-corrected (common model) tile, X
									centersXY_main[1], // centerY, //
									0, // external tile compare
									true,// no_deconvolution,
									false, // ); // transpose);
									null, //final boolean [][]        saturation_imp, // (near) saturated pixels or null
									null); // final double [] overexposed)


							clt_data_aux[0][chn] = new double [4][];
							fract_shiftsXY_aux = extract_correct_tile( // return a pair of residual offsets
									image_rig_data[1],
									width,       // image width
									null, // [color][tileY][tileX][band][pixel]
									clt_data_aux[0][chn], //double  [][]        clt_tile,    // should be double [4][];
									clt_parameters.kernel_step,
									clt_parameters.transform_size,
									dtt,
									chn,
									centersXY_aux[0], // centerX, // center of aberration-corrected (common model) tile, X
									centersXY_aux[1], // centerY, //
									0, // external tile compare
									true,// no_deconvolution,
									false, // ); // transpose);
									null, //final boolean [][]        saturation_imp, // (near) saturated pixels or null
									null); // final double [] overexposed)


							if ((globalDebugLevel > -1) && (tileX == debug_tileX) && (tileY == debug_tileY) && (chn == 2)) {
								System.out.println();
							}

							if ((globalDebugLevel > 0) && (tileX >= debug_tileX - 2) && (tileX <= debug_tileX + 2) &&
									(tileY >= debug_tileY - 2) && (tileY <= debug_tileY+2)) {
									System.out.println("clt_bi_quad(): color="+chn+", tileX="+tileX+", tileY="+tileY+
											" fract_shiftsXY_main[0]="+fract_shiftsXY_main[0]+" fract_shiftsXY_main[1]="+fract_shiftsXY_main[1]);
									System.out.println("clt_bi_quad(): color="+chn+", tileX="+tileX+", tileY="+tileY+
											" fract_shiftsXY_aux[0]="+fract_shiftsXY_aux[0]+" fract_shiftsXY_aux[1]="+fract_shiftsXY_aux[1]);
							}

							// apply residual shift
							fract_shift(    // fractional shift in transform domain. Currently uses sin/cos - change to tables with 2? rotations
									clt_data_main[0][chn], // double  [][]  clt_tile,
									clt_parameters.transform_size,
									fract_shiftsXY_main[0],            // double        shiftX,
									fract_shiftsXY_main[1],            // double        shiftY,
									((globalDebugLevel > 1) && (chn==0) && (tileX >= debug_tileX - 2) && (tileX <= debug_tileX + 2) &&
											(tileY >= debug_tileY - 2) && (tileY <= debug_tileY+2)));
							fract_shift(    // fractional shift in transform domain. Currently uses sin/cos - change to tables with 2? rotations
									clt_data_aux[0][chn], // double  [][]  clt_tile,
									clt_parameters.transform_size,
									fract_shiftsXY_aux[0],            // double        shiftX,
									fract_shiftsXY_aux[1],            // double        shiftY,
									((globalDebugLevel > 1) && (chn==0) && (tileX >= debug_tileX - 2) && (tileX <= debug_tileX + 2) &&
											(tileY >= debug_tileY - 2) && (tileY <= debug_tileY+2)));



							if ((globalDebugLevel > 100) && (debug_tileX == tileX) && (debug_tileY == tileY)) {
								ShowDoubleFloatArrays sdfa_instance = new ShowDoubleFloatArrays(); // just for debugging?
								String [] titles = {"CC0","SC0","CS0","SS0"};
								double [][] dbg_tile = new double [4][];
								for (int i = 0; i < 4; i++) dbg_tile[i]=clt_data_main[0][chn][i & 3];
								sdfa_instance.showArrays(dbg_tile,  clt_parameters.transform_size, clt_parameters.transform_size, true, "MAIN_shifted_x"+tileX+"_y"+tileY+"-z", titles);
								double [][] dbg_tile_aux = new double [16][];
								for (int i = 0; i < 4; i++) dbg_tile[i]=clt_data_aux[0][chn][i & 3];
								sdfa_instance.showArrays(dbg_tile_aux,  clt_parameters.transform_size, clt_parameters.transform_size, true, "AUX_shifted_x"+tileX+"_y"+tileY+"-z", titles);
							}
						}

						int tile_lma_debug_level =  ((tileX == debug_tileX) && (tileY == debug_tileY))? clt_parameters.img_dtt.lma_debug_level : -1;

						if (disparity_bimap != null){ // not null - calculate correlations
//     * @param clt_data_tile_main aberration-corrected FD CLT data for one tile of the main quad camera  [sub-camera][color][quadrant][index]

							double [] inter_corrs_dxy = tileInterCamCorrs(
									false,             // final boolean                                   no_int_x0, // do not offset window to integer - used when averaging low textures to avoid "jumps" for very wide
									clt_parameters,    // final EyesisCorrectionParameters.CLTParameters  clt_parameters,
									fatzero,           // final double                                    fatzero,         // May use correlation fat zero from 2 different parameters - fat_zero and rig.ml_fatzero
									corr2d,            // final Correlation2d                             corr2d,
									clt_data_main,     // double [][][][]                                 clt_data_tile_main,
									clt_data_aux,      // double [][][][]                                 clt_data_tile_aux,
									filter,            // final double []       							filter,
									col_weights,       // final double []       							col_weights,
									0, // ml_hwidth,         // final int                                        ml_hwidth,
						    		null, // ml_data_inter,     // final double []                                 ml_center_corr,
									null, // tcorr_combo,       // double [][]                                     tcorr_combo,
									false,             // final boolean                                   notch_mode,
									tileX,                 // final int              tileX, // only used in debug output
									tileY,                 // final int              tileY,
									tile_lma_debug_level); // final int              debugLevel)
							if (inter_corrs_dxy != null) {
								disparity_bimap[BI_DISP_CROSS_INDEX][nTile] =    inter_corrs_dxy[INDEX_DISP];
								disparity_bimap[BI_STR_CROSS_INDEX][nTile] =     inter_corrs_dxy[INDEX_STRENGTH];
								disparity_bimap[BI_DISP_CROSS_DX_INDEX][nTile] = inter_corrs_dxy[INDEX_DX];
								disparity_bimap[BI_DISP_CROSS_DY_INDEX][nTile] = inter_corrs_dxy[INDEX_DY];
								// TODO: Use strength for the same residual disparity
								disparity_bimap[BI_STR_ALL_INDEX][nTile] =
										Math.pow(inter_corrs_dxy[INDEX_STRENGTH]*
										disparity_bimap[BI_STR_FULL_INDEX][nTile]*
										disparity_bimap[BI_ASTR_FULL_INDEX][nTile], 1.0/3);
							}

						} // if (disparity_map != null){ // not null - calculate correlations

					}
				}
			};
		}
		startAndJoin(threads);
	}
}
