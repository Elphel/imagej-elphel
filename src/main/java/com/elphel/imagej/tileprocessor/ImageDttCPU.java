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
// ← → ↑ ↓ ⇖ ⇗ ⇘ ⇙ ↔ ↕ 

import java.awt.Rectangle;
import java.util.Arrays;
import java.util.concurrent.atomic.AtomicInteger;

import com.elphel.imagej.cameras.CLTParameters;
import com.elphel.imagej.cameras.EyesisCorrectionParameters;
import com.elphel.imagej.common.PolynomialApproximation;
import com.elphel.imagej.common.ShowDoubleFloatArrays;
import com.elphel.imagej.correction.EyesisDCT;
import com.elphel.imagej.gpu.GPUTileProcessor;
import com.elphel.imagej.gpu.TpTask;

import Jama.Matrix;
import ij.ImageStack;

public class ImageDttCPU {
	  static boolean FPGA_COMPARE_DATA= false; // true; // false; //
	  static int     FPGA_SHIFT_BITS =  7; // number of bits for fractional pixel shift
	  static int     FPGA_PIXEL_BITS = 15; // bits to represent pixel data (positive)
	  static int     FPGA_WND_BITS =   17; // bits to represent mclt window (positive for 18-bit signed mpy input)
	  static int     FPGA_DTT_IN =     22; // bits to represent maximal value after folding (input to DTT)
	  static int     FPGA_TILE_SIZE =  22; // size of square side for the composite colors tile (16..22)
	public static  int [][] ZI =
			{{ 0,  1,  2,  3},
			 {-1,  0, -3,  2},
			 {-2, -3,  0,  1},
			 { 3, -2, -1,  0}};
	  
	public static int [][] CORR_PAIRS ={ // {first, second, rot} rot: 0 - as is, 1 - swap y,x
			{0,1,0},
			{2,3,0},
			{0,2,1},
			{1,3,1}};

	public static double[][] PORT_OFFSETS4 = {
			{-0.5, -0.5},
			{ 0.5, -0.5},
			{-0.5,  0.5},
			{ 0.5,  0.5}};

	  // kernels ar designed to have sum = 1.0 and completely reject Bayer modulation for each color
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

//	  static int  QUAD =                           4; // number of cameras in camera
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
//	  static int  IMG_DIFF0_INDEX =               12; // index of noise- normalized image difference for port 0 in disparity map
//	  static int  OVEREXPOSED =                   16; // index of overexposed fraction of all pixels
	  static int  OVEREXPOSED =                   12; // index of overexposed fraction of all pixels
	  static int  IMG_DIFF0_INDEX =               13; // index of noise- normalized image difference for port 0 in disparity map
//	  static int  IMG_TONE_RGB =                  17; // 12 entries of r0,r1,r2,r3,g0,g1,g2,g3,b0,b1,b2,b3
	  
	  static int  IMG_DIFFS =                     16;
	  static int  IMG_TONES_RGB =                 17;
	  
	  static boolean isCorrBit (int indx) {
		  return indx <= DISPARITY_STRENGTH_INDEX; // maybe use DISPARITY_VARIATIONS_INDEX?
	  }
	  static boolean isSliceBit(int indx) {
		  return indx < IMG_DIFFS;
	  }
	  static int  getImgToneRGB(int numSensors) {
		  return IMG_DIFF0_INDEX + numSensors;
	  }
	  static boolean needImgDiffs(int mode) {
		  return ((mode >> IMG_DIFFS) & 1) != 0;
	  }
	  static boolean needTonesRGB(int mode) {
		  return ((mode >> IMG_TONES_RGB) & 1) != 0;
	  }
	  public boolean isDiffIndex(int indx) {
		  return (indx >= IMG_DIFF0_INDEX) && (indx < (IMG_DIFF0_INDEX + numSensors));
	  }
	  public boolean isToneRGBIndex(int indx) {
		  return (indx >= (IMG_DIFF0_INDEX + numSensors)) && (indx < (IMG_DIFF0_INDEX + 4 * numSensors));
	  }
	  public int getImgToneRGB() {
		  return getImgToneRGB(numSensors);
	  }

// remove when not needed
	  static int BITS_ALL_DISPARITIES = (
			  (3 << DISPARITY_INDEX_INT) |           // 0 - disparity from correlation integer pixels, 1 - ortho
			  (3 << DISPARITY_INDEX_CM) |            // 2 - disparity from correlation "center mass", 3 - ortho (only used for fine correction)
			  (1 << DISPARITY_INDEX_HOR) |           // 4 - disparity from correlation of the horizontal pairs with center suppressed
			  (1 << DISPARITY_INDEX_HOR_STRENGTH) |  // 5; // strength for hor mode (emphasis on vertical lines)
			  (1 << DISPARITY_INDEX_VERT) |          // 6; // disparity from correlation of the vertical pairs with center suppressed
			  (1 << DISPARITY_INDEX_VERT_STRENGTH) | // 7; // strength in vert mode (horizontal lines detection)
			  (3 << DISPARITY_INDEX_POLY) |          // 8; // index of disparity value in disparity_map == 2 (0,2 or 4)
			  (1 << DISPARITY_STRENGTH_INDEX) |      // 10; // index of strength data in disparity map ==6
			  (1 << DISPARITY_VARIATIONS_INDEX));    // 11; // index of strength data in disparity map ==6
	  
	  // TODO - use 1 bit per type 
	  static int BITS_OVEREXPOSED = (1     << OVEREXPOSED);
      /*
	  static int BITS_ALL_DIFFS =   (0xf   << IMG_DIFF0_INDEX);
	  static int BITS_TONE_RGB =    (0xfff << IMG_TONE_RGB);
	  */
	  
	  // TODO: modify decoding of bits - turns on/off all channels simultaneously
	  static int BITS_ALL_DIFFS =   (0x1 << IMG_DIFFS);
	  static int BITS_TONE_RGB =    (0x1 << IMG_TONES_RGB);
	  static int BITS_FROM_GPU = BITS_ALL_DIFFS | BITS_TONE_RGB;
	  
	  static String [] DISPARITY_TITLES4 = {
			  "int_disp","int_y_disp","cm_disp","cm_y_disp","hor_disp","hor_strength","vert_disp","vert_strength",
			  "poly_disp", "poly_y_disp", "strength_disp", "vary_disp","overexp","diff0","diff1","diff2","diff3",
			  "r0","r1","r2","r3",
			  "g0","g1","g2","g3",
			  "b0","b1","b2","b3",
			  };
	  static String [] getDisparityTitles(int numSensors) {
		  String [] disparity_titles0 = {
				  "int_disp","int_y_disp","cm_disp","cm_y_disp","hor_disp","hor_strength","vert_disp","vert_strength",
				  "poly_disp", "poly_y_disp", "strength_disp", "vary_disp","overexp"};
		  String [] disparity_titles = new String [disparity_titles0.length + 4* numSensors];
		  int indx = 0;
		  for (String s: disparity_titles0) {
			  disparity_titles[indx++] = s;
		  }
		  for (int i = 0; i < numSensors; i++) disparity_titles[indx++] = "diff"+i;
		  for (int i = 0; i < numSensors; i++) disparity_titles[indx++] = "r"+i;
		  for (int i = 0; i < numSensors; i++) disparity_titles[indx++] = "b"+i;
		  for (int i = 0; i < numSensors; i++) disparity_titles[indx++] = "g"+i;
		  return disparity_titles;
	  }
	  String [] getDisparityTitles() {
		  return getDisparityTitles(numSensors);
	  }
	  
	  static public String[] CORR_TITLES = {
			  "top","bottom","left","right","diag-m","diag-o",
			  "quad","cross","hor","vert",
			  "s-hor","s-vert","s-quad","s-cross","s-quad-cross","s-combo"}; 
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

	  static String [] GEOM_TITLES_DBG ={
			  "px0","py0","px1","py1","px2","py2","px3","py3",
			  "dd0-0","dd0-1","dd0-2","dd0-3",
			  "dd1-0","dd1-1","dd1-2","dd1-3",
			  "dd2-0","dd2-1","dd2-2","dd2-3",
			  "dd3-0","dd3-1","dd3-2","dd3-3"};
	  
	  public static int  ML_OTHER_TARGET =            0;  // Offset to target disparity data in  ML_OTHER_INDEX layer tile
	  public static int  ML_OTHER_GTRUTH =            2;  // Offset to ground truth disparity data in ML_OTHER_INDEX layer tile
	  public static int  ML_OTHER_GTRUTH_STRENGTH =   4;  // Offset to ground truth confidence data in  ML_OTHER_INDEX layer tile
	  public static int  ML_OTHER_GTRUTH_RMS =        6;  // Offset to ground truth RMS in  ML_OTHER_INDEX layer tile
	  public static int  ML_OTHER_GTRUTH_RMS_SPLIT =  8;  // Offset to ground truth combined FG/BG RMS in  ML_OTHER_INDEX layer tile
	  public static int  ML_OTHER_GTRUTH_FG_DISP =   10;  // Offset to ground truth FG disparity in  ML_OTHER_INDEX layer tile
	  public static int  ML_OTHER_GTRUTH_FG_STR =    12;  // Offset to ground truth FG strength in  ML_OTHER_INDEX layer tile
	  public static int  ML_OTHER_GTRUTH_BG_DISP =   14;  // Offset to ground truth BG disparity in  ML_OTHER_INDEX layer tile
	  public static int  ML_OTHER_GTRUTH_BG_STR =    16;  // Offset to ground truth BG strength in  ML_OTHER_INDEX layer tile
	  public static int  ML_OTHER_AUX_DISP =         18;  // Offset to AUX heuristic disparity in  ML_OTHER_INDEX layer tile
	  public static int  ML_OTHER_AUX_STR =          20;  // Offset to AUX heuristic strength in  ML_OTHER_INDEX layer tile


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

	  static int CORR_SEL_BIT_ALL =   0; 
	  static int CORR_SEL_BIT_DIA =   1; 
	  static int CORR_SEL_BIT_SQ =    2;
	  static int CORR_SEL_BIT_NEIB =  3;
	  static int CORR_SEL_BIT_HOR =   4;
	  static int CORR_SEL_BIT_VERT =  5;
	  
	  public static int corrSelEncode(
			  boolean sel_all,
			  boolean sel_dia,
			  boolean sel_sq,
			  boolean sel_neib,
			  boolean sel_hor,
			  boolean sel_vert) {
		  return  (sel_all ?  (1 << CORR_SEL_BIT_ALL):  0) |
				  (sel_dia ?  (1 << CORR_SEL_BIT_DIA):  0) |
				  (sel_sq ?   (1 << CORR_SEL_BIT_SQ):   0) |
				  (sel_neib ? (1 << CORR_SEL_BIT_NEIB): 0) |
				  (sel_hor ?  (1 << CORR_SEL_BIT_HOR):  0) |
				  (sel_vert ? (1 << CORR_SEL_BIT_VERT): 0);
	  }
	  
	  public static int corrSelEncode(ImageDttParameters  img_dtt, int num_sensors) {
		  return corrSelEncode(
				  img_dtt.getMcorrAll  (num_sensors),  // boolean sel_all,
				  img_dtt.getMcorrDia  (num_sensors),  // boolean sel_dia,
				  img_dtt.getMcorrSq   (num_sensors),  // boolean sel_sq,
				  img_dtt.getMcorrNeib (num_sensors),  // boolean sel_neib,
				  img_dtt.getMcorrHor  (num_sensors),  // boolean sel_hor,
				  img_dtt.getMcorrVert (num_sensors)); // boolean sel_vert);
	  }

/*	  
	  int mcorr_sel = ImageDtt. corrSelEncode( // maybe update?
			  clt_parameters.img_dtt.getMcorrAll  (getNumSensors()),  // boolean sel_all,
			  clt_parameters.img_dtt.getMcorrDia  (getNumSensors()),  // boolean sel_dia,
			  clt_parameters.img_dtt.getMcorrSq   (getNumSensors()),  // boolean sel_sq,
			  clt_parameters.img_dtt.getMcorrNeib (getNumSensors()),  // boolean sel_neib,
			  clt_parameters.img_dtt.getMcorrHor  (getNumSensors()),  // boolean sel_hor,
			  clt_parameters.img_dtt.getMcorrVert (getNumSensors())); // boolean sel_vert);
*/
	  
	  public static boolean isCorrAll (int sel) { return ((sel >> CORR_SEL_BIT_ALL)  & 1) != 0;}
	  public static boolean isCorrDia (int sel) { return ((sel >> CORR_SEL_BIT_DIA)  & 1) != 0;}
	  public static boolean isCorrSq  (int sel) { return ((sel >> CORR_SEL_BIT_SQ)   & 1) != 0;}
	  public static boolean isCorrNeib(int sel) { return ((sel >> CORR_SEL_BIT_NEIB) & 1) != 0;}
	  public static boolean isCorrHor (int sel) { return ((sel >> CORR_SEL_BIT_HOR)  & 1) != 0;}
	  public static boolean isCorrVert(int sel) { return ((sel >> CORR_SEL_BIT_VERT) & 1) != 0;}
	  
	  private final boolean monochrome;
	  private final boolean lwir;       // means that no sqrt correction
	  private final boolean aux;
	  private final double scale_strengths; // scale all correlation strengths (to compensate for LPF sigma changes)
	  public final int transform_size;
	  public final int numSensors;
	  public Correlation2d correlation2d = null;
	  final ImageDttParameters imgdtt_params;
	  final double [][] port_offsets;
	  // Save correlation parameters to instantiate  Correlation2d when first needed (called getCorrelation2d(),
	  // do it before threads to prevent races
	  

	  // save lpf arrays (rgb and corr) to use in C-program for debugging
	 public double []   dbg_filter_corr;
	 public double [][] dbg_filter_rbg;

     public static int getImgMask  (int data){ return (data & 0xf);}      // which images to use
     public static int getPairMask (int data){ return ((data >> 4) & 0xf);} // which pairs to combine in the combo:  1 - top, 2 bottom, 4 - left, 8 - right
     public static int setImgMask  (int data, int mask) {return (data & ~0xf) | (mask & 0xf);}
     public static int setPairMask (int data, int mask) {return (data & ~0xf0) | ((mask & 0xf) << 4);}
     public static boolean getForcedDisparity (int data){return (data & 0x100) != 0;}
     public static int     setForcedDisparity (int data, boolean force) {return (data & ~0x100) | (force?0x100:0);}
     public static boolean getOrthoLines (int data){return (data & 0x200) != 0;} // not used in lwir
     public static int     setOrthoLines (int data, boolean force) {return (data & ~0x200) | (force?0x200:0);} // not used in lwir
     
     public int getNumSensors() {
    	 return numSensors;
     }

     // different scale for quad and multi: Multi - diameter == 1.0, quad - side = 1.0
     public static double [][]getPortOffsets(int num_sensors, boolean quad_sequence, boolean topis0) {
    	 if ((num_sensors==4) && quad_sequence) {
    		 return PORT_OFFSETS4;
    	 } else {
    		 double [][] port_offsets = new double [num_sensors][2];
    		 for (int i = 0; i < num_sensors; i++) {
    			 double alpha = 2 * Math.PI * (i + (topis0 ? 0 : 0.5))/num_sensors;
    			 port_offsets[i][0] =  0.5 * Math.sin((alpha));
    			 port_offsets[i][1] = -0.5 * Math.cos((alpha));
    		 }
        	 return port_offsets;
    	 }
     }
     
     
     public Correlation2d getCorrelation2d() {
    	 if (correlation2d == null) {
    		 // instantiate
    			this.correlation2d = new Correlation2d (
    					this.numSensors,      // int                numSensors,
    					this.imgdtt_params,   // ImageDttParameters imgdtt_params,
    					this.transform_size,  // int                transform_size,
    		    		1.0, // 2.0,                  // double             wndx_scale, // (wndy scale is always 1.0)
    		    		this.monochrome,      // boolean            monochrome,
    		    		false);               // boolean            debug) 
    	 }
    	 return correlation2d;
     }
     
	public ImageDttCPU(
			int numSensors,
			int transform_size,
			ImageDttParameters imgdtt_params,
			boolean aux,
			boolean mono,
			boolean lwir,
			double scale_strengths){
		this.numSensors = numSensors;
		this.transform_size = transform_size;
		this.aux = aux;
		this.monochrome = mono || lwir;
		this.lwir =       lwir;
		this.scale_strengths = scale_strengths;
		this.imgdtt_params = imgdtt_params;
		this.port_offsets = getPortOffsets(
				numSensors, // int num_sensors,
				imgdtt_params.mcorr_quad_sequence, // boolean quad_sequence, 
				imgdtt_params.getTopIs0(numSensors));// boolean topis0
	}

	public double getScaleStrengths() {
		return scale_strengths;
	}

	public boolean isLwir() {
		return lwir;
	}

	public boolean isMonochrome() {
		return monochrome;
	}
	// maybe change in the future
	public boolean isAux() {
		return aux;
	}


	public double [][][][] mdctStack( // not used in lwir
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

	public double [][][] lapped_dct( // not used in lwir
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
	public double [] lapped_dct_dbg( // not used in lwir
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

	public void dct_lpf( // not used in lwir
			final double sigma,
			final double [][][] dct_data,
			final int       threadsMax,     // maximal number of threads to launch
			final int       globalDebugLevel)
	{
		final int tilesY=dct_data.length;
		final int tilesX=dct_data[0].length;
		final int nTiles=tilesX*tilesY;
		final double [] filter =  doubleGetCltLpfFd(sigma);
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

	public double [][][][] dct_color_convert( // not used in lwir
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






	public double [] lapped_idct( // not used in lwir
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
	public double [][][][][] clt_aberrations( // not used in lwir
			final double [][]       image_data,
			final int               width,
			final double [][][][][] clt_kernels, // [color][tileY][tileX][band][pixel] , size should match image (have 1 tile around)
			final int               kernel_step,
//			final int               transform_size,
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
//								transform_size,
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
//									transform_size,
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

	@Deprecated
	public double [][][][][][] clt_aberrations_quad( // not used in lwir
			final double              disparity,
			final double [][][]       image_data, // first index - number of image in a quad
			final int                 width,
			final GeometryCorrection  geometryCorrection,
			final double [][][][][][] clt_kernels, // [channel_in_quad][color][tileY][tileX][band][pixel] , size should match image (have 1 tile around)
			final int                 kernel_step,
//			final int                 transform_size,
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
								null,      // double [][] disp_dist used to correct 3D correlations
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
//										transform_size,
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
//											transform_size,
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

	public void printSignsFPGA ( // not used in lwir
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

	public void generateFPGACompareData( // not used in lwir
			final double [][]         image_data, // for selected subcamera
			final double [][]         colorCentersXY, // pixel centers per color (2 - green)
//			final int                 transform_size,
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
//					transform_size,
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

// removing macro and FPGA modes
	
	public double[][] cltMeasureLazyEye ( // returns d,s lazy eye parameters 
			final ImageDttParameters  imgdtt_params,   // Now just extra correlation parameters, later will include, most others
			final int [][]            tile_op,         // [tilesY][tilesX] - what to do - 0 - nothing for this tile
			final double [][]         disparity_array, // [tilesY][tilesX] - individual per-tile expected disparity
			final double [][][]       image_data, // first index - number of image in a quad
		    final boolean [][]        saturation_imp, // (near) saturated pixels or null
			 // correlation results - final and partial

			final int                 width,
			final double              corr_fat_zero,    // add to denominator to modify phase correlation (same units as data1, data2). <0 - pure sum
			final double              corr_red,
			final double              corr_blue,
			final double              corr_sigma,
	  		final double              min_corr,        // 0.02; // minimal correlation value to consider valid
			final GeometryCorrection  geometryCorrection,
			final GeometryCorrection  geometryCorrection_main, // if not null correct this camera (aux) to the coordinates of the main
			final double [][][][][][] clt_kernels, // [channel_in_quad][color][tileY][tileX][band][pixel] , size should match image (have 1 tile around)
			final int                 kernel_step,
			final int                 window_type,
			final double [][]         shiftXY, // [port]{shiftX,shiftY}
			final double              disparity_corr, // disparity at infinity
			final double              shiftX, // shift image horizontally (positive - right) - just for testing
			final double              shiftY, // shift image vertically (positive - down)
			final int                 tileStep, // process tileStep x tileStep cluster of tiles when adjusting lazy eye parameters
			
			final int                 mcorr_sel, // +1 - all, +2 - dia, +4 - sq, +8 - neibs, +16 - hor + 32 - vert 
			final int                 mcorr_comb_width,  // combined correlation tile width
			final int                 mcorr_comb_height, // combined correlation tile full height
			final int                 mcorr_comb_offset, // combined correlation tile height offset: 0 - centered (-height/2 to height/2), height/2 - only positive (0 to height)
			final double              mcorr_comb_disp,   // Combined tile per-pixel disparity for baseline == side of a square
			
			final int                 debug_tileX,
			final int                 debug_tileY,
			final int                 threadsMax,  // maximal number of threads to launch
			final int                 globalDebugLevel)
	{
		final boolean debug_distort= (globalDebugLevel >0); // .false; // true;
		final double [][] debug_offsets = null;
		/*
		final double [][] debug_offsets = new double[getNumSensors()][2];  
		for (int i = 0; i < imgdtt_params.lma_dbg_offset.length; i++) for (int j = 0; j < debug_offsets[i].length; j++) {
			debug_offsets[i][j] = imgdtt_params.lma_dbg_offset[i][j]*imgdtt_params.lma_dbg_scale;
		}
	   */
		final int numcol = 3; // number of colors // keep the same, just do not use [0] and [1], [2] - green

		final int height=image_data[0][0].length/width;
		final int tilesX=width/transform_size;
		final int tilesY=height/transform_size;

		final int clustersX= (tilesX + tileStep - 1) / tileStep;
		final int clustersY= (tilesY + tileStep - 1) / tileStep;
		final double [][] lazy_eye_data = new double [clustersY*clustersX][];
		final double [][] dbg_num_good_tiles = new double [2][clustersY*clustersX];

		final int nClustersInChn=clustersX * clustersY;
		final int clustSize = tileStep*tileStep;

		final int debug_clustX = debug_tileX / tileStep;
		final int debug_clustY = debug_tileY / tileStep;
///tileStep
		// this.correlation2d should be non-null
		boolean [] corr_calculate = null;
		{
			if (isCorrAll  (mcorr_sel)) corr_calculate = correlation2d.selectAll();
			if (isCorrDia  (mcorr_sel)) corr_calculate = correlation2d.selectDiameters  (corr_calculate);
			if (isCorrSq   (mcorr_sel)) corr_calculate = correlation2d.selectSquares    (corr_calculate);
			if (isCorrNeib (mcorr_sel)) corr_calculate = correlation2d.selectNeibs      (corr_calculate);
			if (isCorrHor  (mcorr_sel)) corr_calculate = correlation2d.selectHorizontal (corr_calculate);
			if (isCorrVert (mcorr_sel)) corr_calculate = correlation2d.selectVertical   (corr_calculate);
			correlation2d.setCorrPairs(corr_calculate); // will limit correlation pairs calculation
			correlation2d.generateResample( // should be called before
					mcorr_comb_width,  // combined correlation tile width
					mcorr_comb_height, // combined correlation tile full height
					mcorr_comb_offset, // combined correlation tile height offset: 0 - centered (-height/2 to height/2), height/2 - only positive (0 to height)
					mcorr_comb_disp);
			
		}
		final double [][][][][][] clt_data = new double[numSensors][numcol][tilesY][tilesX][][];
		final Thread[] threads = newThreadArray(threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		final double [] col_weights= new double [numcol]; // colors are RBG
		final double [][] dbg_distort = debug_distort? (new double [4*numSensors][tilesX*tilesY]) : null;

		final double [][] corr_wnd = Corr2dLMA.getCorrWnd(
				transform_size,
				imgdtt_params.lma_wnd);
		final double [] corr_wnd_inv_limited = (imgdtt_params.lma_min_wnd <= 1.0)?  new double [corr_wnd.length * corr_wnd[0].length]: null;
		if (corr_wnd_inv_limited != null) {
			double inv_pwr = imgdtt_params.lma_wnd_pwr - (imgdtt_params.lma_wnd - 1.0); // compensate for lma_wnd
			for (int i = imgdtt_params.lma_hard_marg; i < (corr_wnd.length - imgdtt_params.lma_hard_marg); i++) {
				for (int j = imgdtt_params.lma_hard_marg; j < (corr_wnd.length - imgdtt_params.lma_hard_marg); j++) {
					corr_wnd_inv_limited[i * (corr_wnd.length) + j] = 1.0/Math.max(Math.pow(corr_wnd[i][j],
							inv_pwr),
							imgdtt_params.lma_min_wnd);
				}
			}
		}
		// keep for now for mono, find out  what do they mean for macro mode
		if (isMonochrome()) {
			col_weights[2] = 1.0;// green color/mono
			col_weights[0] = 0;
			col_weights[1] = 0;
		} else {
			col_weights[2] = 1.0/(1.0 + corr_red + corr_blue);    // green color
			col_weights[0] = corr_red *  col_weights[2];
			col_weights[1] = corr_blue * col_weights[2];
		}

		final int corr_size = transform_size * 2 -1;
		if ((globalDebugLevel > -10) && (disparity_corr != 0.0)){
			System.out.println(String.format("Using manual infinity disparity correction of %8.5f pixels",disparity_corr));
		}

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

		final double [] filter =  doubleGetCltLpfFd(corr_sigma);

		// prepare disparity maps and weights
		if (globalDebugLevel > 0){
			System.out.println("corr_fat_zero=         "+corr_fat_zero);
			System.out.println("disparity_array[0][0]= "+disparity_array[0][0]);


		}
		DttRad2 dtt = new DttRad2(transform_size);
		dtt.set_window(window_type);
		final double [] lt_window = dtt.getWin2d();	// [256]
		final double [] lt_window2 = new double [lt_window.length]; // squared
		for (int i = 0; i < lt_window.length; i++) lt_window2[i] = lt_window[i] * lt_window[i];

		if (globalDebugLevel > 1) {
			ShowDoubleFloatArrays sdfa_instance = new ShowDoubleFloatArrays(); // just for debugging?
			sdfa_instance.showArrays(lt_window,  2*transform_size, 2*transform_size, "lt_window");
		}
		Matrix [] corr_rots_aux = null;
		Matrix [][] deriv_rots_aux = null;
		if (geometryCorrection_main != null) {
			corr_rots_aux = geometryCorrection.getCorrVector().getRotMatrices(geometryCorrection.getRotMatrix(true));
			deriv_rots_aux = geometryCorrection.getCorrVector().getRotDeriveMatrices();
		}

		final boolean use_main = corr_rots_aux != null;
		final Matrix [] corr_rots = use_main ? corr_rots_aux : geometryCorrection.getCorrVector().getRotMatrices(); // get array of per-sensor rotation matrices
		final Matrix [][] deriv_rots = use_main ? deriv_rots_aux : geometryCorrection.getCorrVector().getRotDeriveMatrices();

		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				@Override
				public void run() {
					DttRad2 dtt = new DttRad2(transform_size);
					dtt.set_window(window_type);
					int tileY,tileX, clustX, clustY, cTile; // , tIndex; // , chn;
					//						showDoubleFloatArrays sdfa_instance = new showDoubleFloatArrays(); // just for debugging?
					double centerX; // center of aberration-corrected (common model) tile, X
					double centerY; //
					double [][] fract_shiftsXY = new double[numSensors][];
					/*
					Correlation2d corr2d = new Correlation2d(
							numSensors,
							imgdtt_params,              // ImageDttParameters  imgdtt_params,
							transform_size,             // int transform_size,
							2.0,                        //  double wndx_scale, // (wndy scale is always 1.0)
							isMonochrome(), // boolean monochrome,
							(globalDebugLevel > -1));   //   boolean debug)
					corr2d.createOrtoNotch(
							imgdtt_params.getEnhOrthoWidth(isAux()), // double getEnhOrthoWidth(isAux()),
							imgdtt_params.getEnhOrthoScale(isAux()), //double getEnhOrthoScale(isAux()),
							(imgdtt_params.lma_debug_level > 1)); // boolean debug);
					 */							

					double [][] rXY;
					if (use_main) {
						rXY = geometryCorrection.getRXY(true); // boolean use_rig_offsets,
					} else  {
						rXY = geometryCorrection.getRXY(false); // boolean use_rig_offsets,
					}

					//					for (int nTile = ai.getAndIncrement(); nTile < nTilesInChn; nTile = ai.getAndIncrement()) {
					for (int nCluster = ai.getAndIncrement(); nCluster < nClustersInChn; nCluster = ai.getAndIncrement()) {
						clustY = nCluster / clustersX;
						clustX = nCluster % clustersX;
						double [][][] centersXY = new double [clustSize][][];
						double [][][] disp_dist = new double[clustSize][numSensors][]; // used to correct 3D correlations
						double [][][] corrs =     new double [clustSize][][];
						double [][]   disp_str =  new double [clustSize][];
						double [][]   pxpy =      new double [clustSize][2];
//						double [] tile_weights =  new double [clustSize];

						boolean debugCluster =  (clustX == debug_clustX) && (clustY == debug_clustY);
						if (debugCluster) {
							System.out.println("debugCluster");
						}
						int clust_lma_debug_level =  debugCluster? imgdtt_params.lma_debug_level : -5;
// filter only tiles with similar disparity to enable lazy eye for the ERS.
						int num_good_tiles = 0;
						while (true) {
							num_good_tiles = 0; // FIXME: Was missing - uncomment?
							int mnTx = -1, mnTy = -1, mxTx = -1, mxTy = -1;
							double mn = Double.NaN;
							double mx = Double.NaN;
							double avg= 0.0;
							for (int cTileY = 0; cTileY < tileStep; cTileY++) {
								tileY = clustY * tileStep + cTileY ;
								if (tileY < tilesY) {
									for (int cTileX = 0; cTileX < tileStep; cTileX++) {
										tileX = clustX * tileStep + cTileX ;
										if ((tileX < tilesX) && (tile_op[tileY][tileX] != 0)) {
											double d = disparity_array [tileY][tileX];
											avg += d;
											if (!(d <= mx)) {
												mx = d;
												mxTx = tileX;
												mxTy = tileY;
											}
											if (!(d >= mn)) {
												mn = d;
												mnTx = tileX;
												mnTy = tileY;
											}
											num_good_tiles++;
										}
									}
								}
							}
							if (num_good_tiles ==0) {
								break;
							}
							if ((mx-mn) <= imgdtt_params.lma_disp_range ) {
								break;
							}
							avg /= num_good_tiles;
							if ((mx-avg) > (avg-mn)) {
								tile_op[mxTy][mxTx] = 0;
							} else {
								tile_op[mnTy][mnTx] = 0;
							}
						}
						if (num_good_tiles == 0) {
							continue;
						}
						
						dbg_num_good_tiles[0][nCluster] = num_good_tiles;
						int dbg_num_good_lma= 0;
						// num_good_tiles in a cluster > 0
						for (int cTileY = 0; cTileY < tileStep; cTileY++) {
							tileY = clustY * tileStep + cTileY ;
							if (tileY < tilesY) {
								for (int cTileX = 0; cTileX < tileStep; cTileX++) {
									tileX = clustX * tileStep + cTileX ;
									if (tileX < tilesX) { 
										cTile = cTileY * tileStep + cTileX;
//										tIndex = tileY * tilesX + tileX;
										int nTile = tileY * tilesX + tileX; // how is it different from tIndex?
										if (tile_op[tileY][tileX] == 0) {
											disp_str[cTile] = null;
											continue; // nothing to do for this tile
										}
										boolean debugTile =(tileX == debug_tileX) && (tileY == debug_tileY);
										final int [] overexp_all = (saturation_imp != null) ? ( new int [2]): null;
										// Moved from inside chn loop
										centerX = tileX * transform_size + transform_size/2 - shiftX;
										centerY = tileY * transform_size + transform_size/2 - shiftY;
										// TODO: move port coordinates out of color channel loop
										if (    (disparity_array == null) ||
												(disparity_array[tileY] == null) ||
												(Double.isNaN(disparity_array[tileY][tileX]))) {
											System.out.println("Bug with disparity_array !!! - 1, tileX="+tileX+", tileY="+tileY);
											continue; // nothing to do for this tile
										}
										if ((globalDebugLevel > -2) && debugTile) {
											System.out.println ("Calculating offsets");
										}

										if (use_main) { // this is AUX camera that uses main coordinates  // not used in lwir
											centersXY[cTile] =  geometryCorrection.getPortsCoordinatesAndDerivatives(
													geometryCorrection_main, //			GeometryCorrection gc_main,
													true,            // boolean use_rig_offsets,
													corr_rots,       // Matrix []   rots,
													deriv_rots,      //  Matrix [][] deriv_rots,
													null,            // double [][] pXYderiv, // if not null, should be double[8][]
													disp_dist[cTile],       // used to correct 3D correlations
													centerX,
													centerY,
													disparity_array[tileY][tileX] + disparity_corr); // _aux); //  + disparity_corr);
										} else {  // used in lwir
											centersXY[cTile] = geometryCorrection.getPortsCoordinatesAndDerivatives(
													geometryCorrection, //			GeometryCorrection gc_main,
													false,           // boolean use_rig_offsets,
													corr_rots,       // Matrix []   rots,
													deriv_rots,      //  Matrix [][] deriv_rots,
													null,            // double [][] pXYderiv, // if not null, should be double[8][]
													disp_dist[cTile],       // used to correct 3D correlations
													centerX,
													centerY,
													disparity_array[tileY][tileX] + disparity_corr);
										}
										pxpy[cTile][0] = centerX;
										pxpy[cTile][1] = centerY;
										if ((globalDebugLevel > 0) || 	debug_distort ||
												(debugTile && (globalDebugLevel > -2))) {
											for (int i = 0; i < numSensors; i++) {
												System.out.println("clt_aberrations_quad_corr():  tileX="+tileX+", tileY="+tileY+
														" centerX="+centerX+" centerY="+centerY+" disparity="+disparity_array[tileY][tileX]+
														" centersXY["+cTile+"]["+i+"][0]="+centersXY[cTile][i][0]+" centersXY["+cTile+"]["+i+"][1]="+centersXY[cTile][i][1]);
											}
										}
										if (debug_offsets != null) {
											double [][] debug_offsets_xy = new double [debug_offsets.length][2];
											for (int i = 0; i < debug_offsets.length; i++) {
												debug_offsets_xy[i][0] = disp_dist[cTile][i][0] * debug_offsets[i][0] + disp_dist[cTile][i][1] * debug_offsets[i][1];
												debug_offsets_xy[i][1] = disp_dist[cTile][i][2] * debug_offsets[i][0] + disp_dist[cTile][i][3] * debug_offsets[i][1];
											}
											for (int i = 0; i < debug_offsets.length; i++) {
												centersXY[cTile][i][0] += debug_offsets_xy[i][0];
												centersXY[cTile][i][1] += debug_offsets_xy[i][1];
											}
											if ((debug_distort && debugCluster) || debugTile) {
												for (int i = 0; i < numSensors; i++) {
													System.out.println(String.format("%d: {%8.3f, %8.3f}",i,debug_offsets_xy[i][0],debug_offsets_xy[i][1]));
												}
												for (int i = 0; i < numSensors; i++) {
													System.out.println("Corrected clt_aberrations_quad_corr():  tileX="+tileX+", tileY="+tileY+
															" centerX="+centerX+" centerY="+centerY+" disparity="+disparity_array[tileY][tileX]+
															" centersXY["+cTile+"]["+i+"][0]="+centersXY[cTile][i][0]+" centersXY["+cTile+"]["+i+"][1]="+centersXY[cTile][i][1]);
												}
											}
										}

										if ((globalDebugLevel > -1) && (tileX == debug_tileX) && (tileY == debug_tileY)) { // before correction
											System.out.print(disparity_array[tileY][tileX]+"\t"+
													centersXY[cTile][0][0]+"\t"+centersXY[cTile][0][1]+"\t"+
													centersXY[cTile][1][0]+"\t"+centersXY[cTile][1][1]+"\t"+
													centersXY[cTile][2][0]+"\t"+centersXY[cTile][2][1]+"\t"+
													centersXY[cTile][3][0]+"\t"+centersXY[cTile][3][1]+"\t");
										}

										for (int ip = 0; ip < centersXY[cTile].length; ip++){
											centersXY[cTile][ip][0] -= shiftXY[ip][0];
											centersXY[cTile][ip][1] -= shiftXY[ip][1];
										}
										// save disparity distortions for visualization:
										if (dbg_distort != null) {
											for (int cam = 0; cam <numSensors; cam++) {
												dbg_distort[cam * 4 + 0 ][nTile] = disp_dist[cTile][cam][0];
												dbg_distort[cam * 4 + 1 ][nTile] = disp_dist[cTile][cam][1];
												dbg_distort[cam * 4 + 2 ][nTile] = disp_dist[cTile][cam][2];
												dbg_distort[cam * 4 + 3 ][nTile] = disp_dist[cTile][cam][3];
											}
										}

										for (int ncol = 0; ncol <numcol; ncol++) {
											if (!isMonochrome() || (ncol == MONO_CHN)) { // in monochrome mode skip all non-mono (green) channels  // used in lwir (5 of 6 branches)
												if ((globalDebugLevel > -1) && (tileX == debug_tileX) && (tileY == debug_tileY) && (ncol == 2)) {
													System.out.println("\nUsing "+"PIXEL"+" mode, centerX="+centerX+", centerY="+centerY);
													System.out.println(disparity_array[tileY][tileX]+"\t"+
															centersXY[cTile][0][0]+"\t"+centersXY[cTile][0][1]+"\t"+
															centersXY[cTile][1][0]+"\t"+centersXY[cTile][1][1]+"\t"+
															centersXY[cTile][2][0]+"\t"+centersXY[cTile][2][1]+"\t"+
															centersXY[cTile][3][0]+"\t"+centersXY[cTile][3][1]+"\t");
												}

												for (int i = 0; i < numSensors; i++) {
													clt_data[i][ncol][tileY][tileX] = new double [4][];
													// Extract image tiles and kernels, correct aberrations, return (ut do not apply) fractional shifts
													fract_shiftsXY[i] = extract_correct_tile( // return a pair of residual offsets
															image_data[i],
															width,       // image width
															((clt_kernels == null) ? null : clt_kernels[i]), // [color][tileY][tileX][band][pixel]
															clt_data[i][ncol][tileY][tileX], //double  [][]        clt_tile,    // should be double [4][];
															kernel_step,
															dtt,
															ncol,
															centersXY[cTile][i][0], // centerX, // center of aberration-corrected (common model) tile, X
															centersXY[cTile][i][1], // centerY, //
															((!FPGA_COMPARE_DATA && (globalDebugLevel > -1) && (tileX == debug_tileX) && (tileY == debug_tileY) && (ncol == 2) && (i==0)) ? (globalDebugLevel + 0) : 0), // external tile compare
															false, // no_deconvolution,
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
													for (int i = 0; i < numSensors; i++) {
														System.out.println("clt_aberrations_quad(): color="+ncol+", tileX="+tileX+", tileY="+tileY+
																" fract_shiftsXY["+i+"][0]="+fract_shiftsXY[i][0]+" fract_shiftsXY["+i+"][1]="+fract_shiftsXY[i][1]);
													}
												}

												// apply residual shift
												for (int i = 0; i < numSensors; i++) {
													fract_shift(    // fractional shift in transform domain. Currently uses sin/cos - change to tables with 2? rotations
															clt_data[i][ncol][tileY][tileX], // double  [][]  clt_tile,
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

											} else { // if (!isMonochrome() || (chn == MONO_CHN) || macro_mode) { // in monochrome mode skip all non-mono (green) channels
												for (int i = 0; i < numSensors; i++) {  // used in lwir
													clt_data[i][ncol] = null; // erase unused clt_data
												}
											}
										}// end of for (int chn = 0; chn <numcol; chn++)
										
										int tile_lma_debug_level =  ((tileX == debug_tileX) && (tileY == debug_tileY))? imgdtt_params.lma_debug_level : -1;
										// all color channels are done here
										// calculate all selected pairs correlations
										corrs[cTile] = correlation2d.correlateCompositeFD(
									    		clt_data,                     // double [][][][][][] clt_data,
									    		tileX,                        // int                 tileX,
									    		tileY,                        // int                 tileY,
									    		correlation2d.getCorrPairs(), // boolean[]           pairs_mask,
									    		filter,                       // double []           lpf,
									    		getScaleStrengths(),          // double              scale_value, // scale correlation value
									    		col_weights,                  // double []           col_weights,
									    		corr_fat_zero);               // double              fat_zero)

										disp_str[cTile] = new double[2];
										if (tile_lma_debug_level > 0) {
											System.out.println("Will run getMaxXSOrtho( ) for tileX="+tileX+", tileY="+tileY);
										}
										// get CM disparity/strength
										if (true) {
							                double [] corr_combo_tile = correlation2d.accumulateInit();
							                double sumw = correlation2d.accummulatePairs(
							                		corr_combo_tile,            // double []   accum_tile,
							                        corrs[cTile],              // double [][] corr_tiles,
							                        correlation2d.selectAll(), // boolean []  selection,
							                        1.0);             // double      weight);
							                correlation2d.normalizeAccumulatedPairs(
							                		corr_combo_tile,
							                        sumw);
							                //tile_weights[cTile] = sumw;
							                // double [] disp_str_combo = new double[2];
							                int [] ixy =   correlation2d.getMaxXYInt( // find integer pair or null if below threshold // USED in lwir
							                		corr_combo_tile, // double [] data,      // [data_size * data_size]
							                		null, // disp_str_combo,
							                        correlation2d.getCombWidth(), //       data_width,
							                        correlation2d.getCombHeight()/2 - correlation2d.getCombOffset(), // int       center_row, ??????????????
							                        true, // boolean   axis_only,
							                        -1.0, // imgdtt_params.min_corr,  // ???? double    minMax,    // minimal value to consider (at integer location, not interpolated)
							                        false); // debugCluster); // tile_lma_debug_level > 0); // boolean   debug);
											double [] corr_stat = correlation2d.getMaxXCm(         // get fractional center as a "center of mass" inside circle/square from the integer max
													corr_combo_tile,                      // double [] data,      // [data_size * data_size]
													correlation2d.getCombWidth(),        // int       data_width,      //  = 2 * transform_size - 1;
													correlation2d.getCombHeight()/2 - correlation2d.getCombOffset(),// int       center_row,
													ixy[0],                              // int       ixcenter,  // integer center x
													false); // debugCluster); // (tile_lma_debug_level > 0)); // boolean   debug);
											if (corr_stat != null) { // almost always
												disp_str[cTile] = new double [] {-corr_stat[0], corr_stat[1]};
											}
										}
										
										

										// debug new LMA correlations
										int tdl = debugCluster ? tile_lma_debug_level : -3;
										// find disp_str for each tile in a cluster
										if (globalDebugLevel > 1000) { // true) { // debugCluster1) {
											if (debugCluster && (globalDebugLevel > -1)) { // -2)) {
												System.out.println("Will run new LMA for tileX="+tileX+", tileY="+tileY);
											}
											double [] poly_disp = {Double.NaN, 0.0};
											/*
											Corr2dLMA lma2 = corr2d.corrLMA2Single(
													imgdtt_params,                // ImageDttParameters  imgdtt_params,
										    		false,                        // boolean             adjust_ly, // adjust Lazy Eye
													corr_wnd,                     // double [][]         corr_wnd, // correlation window to save on re-calculation of the window
													corr_wnd_inv_limited,         // corr_wnd_limited, // correlation window, limited not to be smaller than threshold - used for finding max/convex areas (or null)
													corrs[cTile],                        // double [][]         corrs,
													disp_dist[cTile],
													rXY,                          // double [][]         rXY, // non-distorted X,Y offset per nominal pixel of disparity
													corr2d.longToArray(imgdtt_params.dbg_pair_mask), // imgdtt_params.dbg_pair_mask,  // int                 pair_mask, // which pairs to process
													null,                         // disp_str[cTile],  //corr_stat[0],                 // double    xcenter,   // preliminary center x in pixels for largest baseline
													poly_disp,                    // double[]            poly_ds,    // null or pair of disparity/strength
													imgdtt_params.ortho_vasw_pwr, // double    vasw_pwr,  // value as weight to this power,
													tdl, // tile_lma_debug_level, //+2,         // int                 debug_level,
													tileX,                        // int                 tileX, // just for debug output
													tileY);                       // int                 tileY
											*/
											if (debugCluster) {
												System.out.println("Before single: clustX="+clustX+", clustY="+clustY+", tileX="+tileX+
														", tileY="+tileY+", cTile="+cTile+", cTileX="+cTileX+", cTileY="+cTileY);
											}
											
											Corr2dLMA lma2 = correlation2d.corrLMA2Single(
													imgdtt_params,                // ImageDttParameters  imgdtt_params,
													imgdtt_params.lmas_LY_single_LY, // false, // false,                        // boolean             adjust_ly, // adjust Lazy Eye
													corr_wnd,                     // double [][]         corr_wnd, // correlation window to save on re-calculation of the window
													corr_wnd_inv_limited,         // corr_wnd_limited, // correlation window, limited not to be smaller than threshold - used for finding max/convex areas (or null)
													corrs[cTile], // corrs,          // double [][]         corrs,
													disp_dist[cTile],
													rXY,                          // double [][]         rXY, // non-distorted X,Y offset per nominal pixel of disparity
													// all that are not null in corr_tiles
													correlation2d.selectAll(),   // longToArray(imgdtt_params.dbg_pair_mask),  // int                 pair_mask, // which pairs to process
													null, // disp_str,  //corr_stat[0],                 // double    xcenter,   // preliminary center x in pixels for largest baseline
													poly_disp,                    // double[]            poly_ds,    // null or pair of disparity/strength
													imgdtt_params.ortho_vasw_pwr, // double    vasw_pwr,  // value as weight to this power,
													tdl, // tile_lma_debug_level, // +2,         // int                 debug_level,
													tileX,                        // int                 tileX, // just for debug output
													tileY );                      // int                 tileY
											
											disp_str[cTile] = null;
											if (lma2 != null) {
												dbg_num_good_lma ++;
												disp_str[cTile] = lma2.lmaDisparityStrength(
									    				imgdtt_params.lmas_min_amp,      //  minimal ratio of minimal pair correlation amplitude to maximal pair correlation amplitude
														imgdtt_params.lmas_max_rel_rms,  // maximal relative (to average max/min amplitude LMA RMS) // May be up to 0.3)
														imgdtt_params.lmas_min_strength, // minimal composite strength (sqrt(average amp squared over absolute RMS)
														imgdtt_params.lmas_min_ac,       // minimal of A and C coefficients maximum (measures sharpest point/line)
														imgdtt_params.lmas_min_min_ac,   // minimal of A and C coefficients minimum (measures sharpest point)
														imgdtt_params.lmas_max_area,     // double  lma_max_area,     // maximal half-area (if > 0.0)
														imgdtt_params.lma_str_scale,     // convert lma-generated strength to match previous ones - scale
														imgdtt_params.lma_str_offset     // convert lma-generated strength to match previous ones - add to result
														)[0];
													if ((disp_str[cTile]!=null) && Double.isNaN(disp_str[cTile][1])) {
														System.out.println();
													}
												if (debugCluster || (tile_lma_debug_level > 0)) {
													double [][] ds_dbg = {disp_str[cTile]};
													System.out.println("multi: clustX="+clustX+", clustY="+clustY+", tileX="+tileX+", tileY="+tileY+", dbg_num_good_lma="+dbg_num_good_lma);
													lma2.printStats(ds_dbg,1);
												}
											} else {
												if (debugCluster || (tile_lma_debug_level > 0)) {
													System.out.println("LMA failed, dbg_num_good_lma="+dbg_num_good_lma);
												}
											}
										} // if (globalDebugLevel > 1000) - skip single-tile LMA, similar to GPU version
									}
								}
							}
						} //for (int cTileY = 0; cTileY < tileStep; cTileY++) {
						dbg_num_good_tiles[1][nCluster] = dbg_num_good_lma;

						// above - scanned each tile in a cluster
						if (true) { //debugCluster1) {
							if (debugCluster || (globalDebugLevel > 0)) {
								System.out.println("Will run new LMA for clustX="+clustX+", clustY="+clustY);
							}
							Corr2dLMA lma2 = null;
							if (imgdtt_params.lma_multi_cons) {
								double [][] corrs_cons = new double [correlation2d.getNumPairs()][];
								double [][] disp_dist_cons = new double[numSensors][]; // used to correct 3D correlations
								double []   disp_str_cons = {0.0,0.0};
								int num_tiles = 0;
								double sum_w = 0.0;
								double sum_wd = 0.0;
								for (int nTile = 0; nTile < disp_str.length; nTile++) if (corrs[nTile] != null) { //(disp_str[nTile] != null) {
									// Init if it is a first non-null tile
									if (num_tiles == 0) {
										for (int np = 0; np < corrs_cons.length; np++) if (corrs[nTile][np] != null){
											corrs_cons[np] = new double [corrs[nTile][np].length];
										}
										for (int nsens = 0; nsens < numSensors; nsens++) {
											for (int i = 0; i < disp_dist[nTile][nsens].length; i++) {
												disp_dist_cons[nsens] = new double [disp_dist[nTile][nsens].length];
											}
										}
									}
									
									double w = disp_str[nTile][1];
									sum_w += w;
									sum_wd += w * disp_str[nTile][0];
									for (int np = 0; np < corrs_cons.length; np++) if (corrs[nTile][np] != null) {
										for (int i = 0; i < corrs_cons[np].length; i++) {
											corrs_cons[np][i] += w * corrs[nTile][np][i];
										}
									}
									for (int nsens = 0; nsens < numSensors; nsens++) {
										for (int i = 0; i < disp_dist_cons[nsens].length; i++) {
											disp_dist_cons[nsens][i] += w * disp_dist[nTile][nsens][i];
										}
									}
									num_tiles ++;
								}
								if ((num_tiles > 0) && (sum_w > 0.0)) {
//									if (sum_w > 0.0) {
										disp_str_cons[0] = sum_wd/sum_w;
										disp_str_cons[1] = sum_w/num_tiles; 
//									}
									double s = 1.0/sum_w; // num_tiles;
									for (int np = 0; np < corrs_cons.length; np++) {
										if (corrs_cons[np] != null) {
											for (int i = 0; i < corrs_cons[np].length; i++) {
												corrs_cons[np][i] *= s;
											}
										}
									}
									for (int nsens = 0; nsens < numSensors; nsens++) {
										for (int i = 0; i < disp_dist_cons[nsens].length; i++) {
											disp_dist_cons[nsens][i] *=s;
										}
									}
									
									// Estimate disparity - consolidate all correlation pairs and find maximum
					                double [] corr_combo_all = correlation2d.accumulateInit();
					                double sumw = correlation2d.accummulatePairs(
					                        corr_combo_all,            // double []   accum_tile,
					                        corrs_cons,                // double [][] corr_tiles,
					                        correlation2d.selectAll(), // boolean []  selection,
					                        1.0);             // double      weight);
					                correlation2d.normalizeAccumulatedPairs(
					                        corr_combo_all,
					                        sumw);
									
					                // find argmax on y==0
					                double [] disp_str_combo = new double[2];
					                int [] ixy =   correlation2d.getMaxXYInt( // find integer pair or null if below threshold // USED in lwir
					                		corr_combo_all, // double [] data,      // [data_size * data_size]
					                		disp_str_combo,
					                        correlation2d.getCombWidth(), //       data_width,
					                        correlation2d.getCombHeight()/2 - correlation2d.getCombOffset(), // int       center_row, ??????????????
					                        true, // boolean   axis_only,
					                        imgdtt_params.min_corr,  // ???? double    minMax,    // minimal value to consider (at integer location, not interpolated)
					                        false); // debugCluster); // tile_lma_debug_level > 0); // boolean   debug);
									
									if (debugCluster) { // && (globalDebugLevel > -1)) { // -2)) {
										System.out.println("Will run new LMA for clustX="+clustX+", clustY="+clustY);
										if (disp_dist_cons != null) {
											System.out.println("disp_dist_cons[0]="+disp_dist_cons[0]+", disp_dist_cons[0]="+disp_dist_cons[0]+" (this tile - weighted average will be discarded)");
										}
//										(new ShowDoubleFloatArrays()).showArrays(corrs_cons,  15, 15, true, "corrs_cons_CX"+clustX+"-CY"+clustY,correlation2d.getCorrTitles());
										// 					mcorr_comb_width,  // combined correlation tile width
										//mcorr_comb_height, // combined correlation tile full height
//										(new ShowDoubleFloatArrays()).showArrays(corr_combo_all,  mcorr_comb_width, mcorr_comb_height, "corr_combo_all_CX"+clustX+"-CY"+clustY);
									}
									
									if (ixy != null) { //TODO - for CM use magic!
										double [] corr_stat = correlation2d.getMaxXCm(         // get fractional center as a "center of mass" inside circle/square from the integer max
												corr_combo_all,                      // double [] data,      // [data_size * data_size]
												correlation2d.getCombWidth(),        // int       data_width,      //  = 2 * transform_size - 1;
												correlation2d.getCombHeight()/2 - correlation2d.getCombOffset(),// int       center_row,
												ixy[0],                              // int       ixcenter,  // integer center x
												// corr_wndy,                        // double [] window_y,  // (half) window function in y-direction(perpendicular to disparity: for row0  ==1
												// corr_wndx,                        // double [] window_x,  // half of a window function in x (disparity) direction
												false); // debugCluster); // (tile_lma_debug_level > 0)); // boolean   debug);
										if (corr_stat != null) {
											disp_str_combo[0] = -corr_stat[0]; // -ixy[0]; for CM use magic???
											disp_str_combo[1] =  corr_stat[1];
										}
										if (debugCluster) { // && (globalDebugLevel > -1)) { // -2)) {
											System.out.println("Will run new LMA for clustX="+clustX+", clustY="+clustY);
											if (disp_str_combo != null) {
												System.out.println("disp_str_combo[0]="+disp_str_combo[0]+", disp_str_combo[1]="+disp_str_combo[1]);
											}
//											(new ShowDoubleFloatArrays()).showArrays(corrs_cons,  15, 15, true, "corrs_cons_CX"+clustX+"-CY"+clustY,correlation2d.getCorrTitles());
										}
										
										lma2 = correlation2d.corrLMA2Single(
												imgdtt_params,                  // ImageDttParameters  imgdtt_params,
												true,// false, // false,                        // boolean             adjust_ly, // adjust Lazy Eye
												corr_wnd,                       // double [][]         corr_wnd, // correlation window to save on re-calculation of the window
												corr_wnd_inv_limited,           // corr_wnd_limited, // correlation window, limited not to be smaller than threshold - used for finding max/convex areas (or null)
												corrs_cons,                     // corrs,          // double [][]         corrs,
												disp_dist_cons,
												rXY,                          // double [][]         rXY, // non-distorted X,Y offset per nominal pixel of disparity
												// all that are not null in corr_tiles
												correlation2d.selectAll(),   // longToArray(imgdtt_params.dbg_pair_mask),  // int                 pair_mask, // which pairs to process
												// TODO: Verify sign is correct (disparity, not X0)?
												disp_str_combo, // null, // disp_str,  //corr_stat[0],                 // double    xcenter,   // preliminary center x in pixels for largest baseline
												disp_str_combo, // disp_str_cons,                // double[]            poly_ds,    // null or pair of disparity/strength
												imgdtt_params.ortho_vasw_pwr, // double    vasw_pwr,  // value as weight to this power,
												clust_lma_debug_level + 0, // tdl, // tile_lma_debug_level, // +2,         // int                 debug_level,
												clustX,                        // int                 tileX, // just for debug output
												clustY );                      // int                 tileY
									}
								}
							} else {
								lma2 = correlation2d.corrLMA2Multi( // slow
										imgdtt_params,                // ImageDttParameters  imgdtt_params,
										tileStep,                     // int                 clust_width,
										corr_wnd,                     // double [][]         corr_wnd, // correlation window to save on re-calculation of the window
										corr_wnd_inv_limited,         // corr_wnd_limited, // correlation window, limited not to be smaller than threshold - used for finding max/convex areas (or null)
										corrs, // corrs,          // double [][][]         corrs,
										disp_dist,
										rXY,                          // double [][]         rXY, // non-distorted X,Y offset per nominal pixel of disparity
										// all that are not null in corr_tiles
										correlation2d.selectAll(),   // longToArray(imgdtt_params.dbg_pair_mask),  // int                 pair_mask, // which pairs to process
										disp_str,  //corr_stat[0],                 // double    xcenter,   // preliminary center x in pixels for largest baseline
										imgdtt_params.ortho_vasw_pwr, // double    vasw_pwr,  // value as weight to this power,
										clust_lma_debug_level + 0, // // +2,         // int                 debug_level,
										clustX,                        // int                 tileX, // just for debug output
										clustY );                      // int                 tileY
							}
							if (lma2 != null) {
								double [][] ddnd = lma2.getDdNd();
								double [] stats  = lma2.getStats(num_good_tiles);
								double [][] lma_ds = lma2.lmaDisparityStrength(
					    				imgdtt_params.lmas_min_amp,      //  minimal ratio of minimal pair correlation amplitude to maximal pair correlation amplitude
					    				imgdtt_params.lma_max_rel_rms,  // maximal relative (to average max/min amplitude LMA RMS) // May be up to 0.3)
					    				imgdtt_params.lma_min_strength, // minimal composite strength (sqrt(average amp squared over absolute RMS)
					    				imgdtt_params.lma_min_ac,       // minimal of A and C coefficients maximum (measures sharpest point/line)
										imgdtt_params.lma_min_min_ac,   // minimal of A and C coefficients minimum (measures sharpest point)
					    				imgdtt_params.lma_max_area,      //double  lma_max_area,     // maximal half-area (if > 0.0)
					    				1.0, // imgdtt_params.lma_str_scale,    // convert lma-generated strength to match previous ones - scale
					    				0.0); // imgdtt_params.lma_str_offset);  // convert lma-generated strength to match previous ones - add to result
//								double [][] extra_stats = lma2.getTileStats();
								
								if (debugCluster) {
									System.out.println("ClustX="+clustX+", clustY="+clustY);
							        lma2.printStats(lma_ds, 1);
							        if (ddnd != null) {
							            double [][] dxy= new double [ddnd.length][2];
							            for (int i = 0; i < dxy.length; i++) {
							                dxy[i][0] = ddnd[i][0] * rXY[i][0] - ddnd[i][1] * rXY[i][1];
							                dxy[i][1] = ddnd[i][0] * rXY[i][1] + ddnd[i][1] * rXY[i][0];
							            }
							            System.out.print("       Port:  ");
							            for (int i = 0; i < dxy.length; i++) System.out.print(String.format("   %2d   ", i)); System.out.println();
							            System.out.print("Radial_in =  [");
							            for (int i = 0; i < dxy.length; i++) System.out.print(String.format(" %6.3f,", ddnd[i][0])); System.out.println("]");
							            System.out.print("Tangent_CW = [");
							            for (int i = 0; i < dxy.length; i++) System.out.print(String.format(" %6.3f,", ddnd[i][1])); System.out.println("]");
							            System.out.print("X =          [");
							            for (int i = 0; i < dxy.length; i++) System.out.print(String.format(" %6.3f,", dxy[i][0])); System.out.println("]");
							            System.out.print("Y =          [");
							            for (int i = 0; i < dxy.length; i++) System.out.print(String.format(" %6.3f,", dxy[i][1])); System.out.println("]");
							            System.out.println();
							        }
								} 
								
								
								//		final double [][] lazy_eye_data = new double [clustersY*clustersX][];
								// calculate average disparity per cluster using a sum of the disparity_array and the result of the LMA
								lazy_eye_data[nCluster] = new double [ExtrinsicAdjustment.get_INDX_LENGTH(numSensors)];
								double sum_w = 0;
								if (imgdtt_params.lma_multi_cons) {
									if ((lma_ds[0] != null) && (lma_ds[0][1]> 0.0)) {
										for (int cTileY = 0; cTileY < tileStep; cTileY++) {
											tileY = clustY * tileStep + cTileY ;
											if (tileY < tilesY) {
												for (int cTileX = 0; cTileX < tileStep; cTileX++) {
													tileX = clustX * tileStep + cTileX ;
													if (tileX < tilesX) {
														cTile = cTileY * tileStep + cTileX;
//														if ((lma_ds[cTile] != null) && (lma_ds[cTile][1]> 0.0)) {
														if (disp_str[cTile] != null) {
//															double w = lma_ds[cTile][1];
															double w = disp_str[cTile][1];
															lazy_eye_data[nCluster][ExtrinsicAdjustment.INDX_DISP] += (lma_ds[0][0] + disparity_array[tileY][tileX] + disparity_corr) * w;
															lazy_eye_data[nCluster][ExtrinsicAdjustment.INDX_TARGET] += (disparity_array[tileY][tileX] + disparity_corr) * w;
															lazy_eye_data[nCluster][ExtrinsicAdjustment.INDX_DIFF] += lma_ds[0][0] * w;
															lazy_eye_data[nCluster][ExtrinsicAdjustment.INDX_PX + 0] += pxpy[cTile][0] * w;
															lazy_eye_data[nCluster][ExtrinsicAdjustment.INDX_PX + 1] += pxpy[cTile][1] * w;
															for (int cam = 0; cam < numSensors; cam++) {
																lazy_eye_data[nCluster][ExtrinsicAdjustment.get_INDX_DYDDISP0(numSensors) + cam] += disp_dist[cTile][cam][2] * w;
																lazy_eye_data[nCluster][ExtrinsicAdjustment.get_INDX_PYDIST(numSensors) + cam] += centersXY[cTile][cam][1] * w;
															}
															sum_w += w;
														}
													}
												}
											}
										}
									}
								} else {
									for (int cTileY = 0; cTileY < tileStep; cTileY++) {
										tileY = clustY * tileStep + cTileY ;
										if (tileY < tilesY) {
											for (int cTileX = 0; cTileX < tileStep; cTileX++) {
												tileX = clustX * tileStep + cTileX ;
												if (tileX < tilesX) {
													cTile = cTileY * tileStep + cTileX;
													if ((lma_ds[cTile] != null) && (lma_ds[cTile][1]> 0.0)) {
														double w = lma_ds[cTile][1];
														lazy_eye_data[nCluster][ExtrinsicAdjustment.INDX_DISP] += (lma_ds[cTile][0] + disparity_array[tileY][tileX] + disparity_corr) * w;
														lazy_eye_data[nCluster][ExtrinsicAdjustment.INDX_TARGET] += (disparity_array[tileY][tileX] + disparity_corr) * w;
														lazy_eye_data[nCluster][ExtrinsicAdjustment.INDX_DIFF] += lma_ds[cTile][0] * w;
														lazy_eye_data[nCluster][ExtrinsicAdjustment.INDX_PX + 0] += pxpy[cTile][0] * w;
														lazy_eye_data[nCluster][ExtrinsicAdjustment.INDX_PX + 1] += pxpy[cTile][1] * w;
														for (int cam = 0; cam < numSensors; cam++) {
															lazy_eye_data[nCluster][ExtrinsicAdjustment.get_INDX_DYDDISP0(numSensors) + cam] += disp_dist[cTile][cam][2] * w;
															lazy_eye_data[nCluster][ExtrinsicAdjustment.get_INDX_PYDIST(numSensors) + cam] += centersXY[cTile][cam][1] * w;
														}
														sum_w += w;
													}
												}
											}
										}
									}
								}
								if (sum_w > 0.0) {
									lazy_eye_data[nCluster][ExtrinsicAdjustment.INDX_STRENGTH] = stats[0];
									lazy_eye_data[nCluster][ExtrinsicAdjustment.INDX_DISP]   /= sum_w;
									lazy_eye_data[nCluster][ExtrinsicAdjustment.INDX_TARGET] /= sum_w;
									lazy_eye_data[nCluster][ExtrinsicAdjustment.INDX_DIFF]   /= sum_w;
									lazy_eye_data[nCluster][ExtrinsicAdjustment.INDX_PX + 0] /= sum_w;
									lazy_eye_data[nCluster][ExtrinsicAdjustment.INDX_PX + 1] /= sum_w;
									for (int cam = 0; cam < numSensors; cam++) {
										lazy_eye_data[nCluster][ExtrinsicAdjustment.get_INDX_DYDDISP0(numSensors) + cam] /= sum_w;
										lazy_eye_data[nCluster][ExtrinsicAdjustment.get_INDX_PYDIST(numSensors) +   cam] /= sum_w;
									}
									for (int cam = 0; cam < ddnd.length; cam++) {
										if (ddnd[cam] != null) { //convert to x,y from dd/nd
											lazy_eye_data[nCluster][2 * cam + ExtrinsicAdjustment.INDX_X0 + 0] = ddnd[cam][0] * rXY[cam][0] - ddnd[cam][1] * rXY[cam][1];
											lazy_eye_data[nCluster][2 * cam + ExtrinsicAdjustment.INDX_X0 + 1] = ddnd[cam][0] * rXY[cam][1] + ddnd[cam][1] * rXY[cam][0];
											lazy_eye_data[nCluster][ExtrinsicAdjustment.get_INDX_DD0(numSensors) + cam] = ddnd[cam][0];
											lazy_eye_data[nCluster][ExtrinsicAdjustment.get_INDX_ND0(numSensors) + cam] = ddnd[cam][1];
										} else {
											lazy_eye_data[nCluster][2 * cam + ExtrinsicAdjustment.INDX_X0 + 0] = Double.NaN;
											lazy_eye_data[nCluster][2 * cam + ExtrinsicAdjustment.INDX_X0 + 1] = Double.NaN;
											lazy_eye_data[nCluster][ExtrinsicAdjustment.get_INDX_DD0(numSensors) + cam] = Double.NaN;
											lazy_eye_data[nCluster][ExtrinsicAdjustment.get_INDX_ND0(numSensors) + cam] = Double.NaN;
										}
									}
								} else {
									lazy_eye_data[nCluster] = null;
								}
							}
						}
					}
				}
			};
		}
		startAndJoin(threads);
		if (dbg_num_good_tiles != null) {
//			(new ShowDoubleFloatArrays()).showArrays(dbg_num_good_tiles,  clustersX, clustersY, true, "num_good_tiles"); // , dbg_titles);
		}
		if ((dbg_distort != null) &&(globalDebugLevel >=0)) {
			(new ShowDoubleFloatArrays()).showArrays(dbg_distort,  tilesX, tilesY, true, "disparity_distortions"); // , dbg_titles);
		}
		return lazy_eye_data; // clt_data;
	}
	
	
	public double [][][][][][] clt_aberrations_quad_corr( // 08/28/21 version
			final ImageDttParameters  imgdtt_params,   // Now just extra correlation parameters, later will include, most others
			final int                 macro_scale,     // to correlate tile data instead of the pixel data: 1 - pixels, 8 - tiles
			final int [][]            tile_op,         // [tilesY][tilesX] - what to do - 0 - nothing for this tile
			final double [][]         disparity_array, // [tilesY][tilesX] - individual per-tile expected disparity
			final double [][][]       image_data,      // first index - number of image in a quad
		    final boolean [][]        saturation_imp,  // (near) saturated pixels or null
			 // correlation results
		    final double [][][][]     clt_corr_out,   // sparse (by the first index) [type][tilesY][tilesX][(2*transform_size-1)*(2*transform_size-1)] or null
		    final double [][][][]     clt_combo_out,  // sparse (by the first index) [type][tilesY][tilesX][(combo_tile_size] or null
			final double [][][][]     clt_combo_dbg,  // generate sparse  partial rotated/scaled pairs
			// When clt_mismatch is non-zero, no far objects extraction will be attempted
			// clt_mismatch is used in older code, not supported in GPU - there is cltMeasureLazyEye for that purpose
			// this.correlation2d should be not null if disparity_map != null
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
			final int                 window_type,
			final double [][]         shiftXY, // [port]{shiftX,shiftY}
			final double              disparity_corr, // disparity at infinity
			final double [][][]       fine_corr, // quadratic coefficients for fine correction (or null)
			final double              corr_magic_scale, // still not understood coefficient that reduces reported disparity value.  Seems to be around 0.85
			final double              shiftX, // shift image horizontally (positive - right) - just for testing
			final double              shiftY, // shift image vertically (positive - down)
			
			final int                 mcorr_sel,    // Which pairs to correlate // +1 - all, +2 - dia, +4 - sq, +8 - neibs, +16 - hor + 32 - vert
			final int                 mcorr_comb_width,  // combined correlation tile width
			final int                 mcorr_comb_height, // combined correlation tile full height
			final int                 mcorr_comb_offset, // combined correlation tile height offset: 0 - centered (-height/2 to height/2), height/2 - only positive (0 to height)
			final double              mcorr_comb_disp,   // Combined tile per-pixel disparity for baseline == side of a square
			
			final int                 debug_tileX,
			final int                 debug_tileY,
			final boolean             no_fract_shift,
			final boolean             no_deconvolution,
			final int                 threadsMax,  // maximal number of threads to launch
			final int                 globalDebugLevel)
	{
		// Initiate Correlation2d instance
//		getCorrelation2d(); // should actually be done by caller to determine size of clt_corr_out (first index)
		
		final int height=image_data[0][0].length/width;
		final int tilesX=width/transform_size;
		final int tilesY=height/transform_size;
		boolean [][] pcombo_sels = null;
		if (correlation2d != null){
			// Initialize correlation pairs selection to be used by all threads
			boolean [] corr_calculate = null;
			if (isCorrAll  (mcorr_sel)) corr_calculate = correlation2d.selectAll();
			if (isCorrDia  (mcorr_sel)) corr_calculate = correlation2d.selectDiameters  (corr_calculate);
			if (isCorrSq   (mcorr_sel)) corr_calculate = correlation2d.selectSquares    (corr_calculate);
			if (isCorrNeib (mcorr_sel)) corr_calculate = correlation2d.selectNeibs      (corr_calculate);
			if (isCorrHor  (mcorr_sel)) corr_calculate = correlation2d.selectHorizontal (corr_calculate);
			if (isCorrVert (mcorr_sel)) corr_calculate = correlation2d.selectVertical   (corr_calculate);
			correlation2d.setCorrPairs(corr_calculate); // will limit correlation pairs calculation
			correlation2d.generateResample( // should be called before
					mcorr_comb_width,  // combined correlation tile width
					mcorr_comb_height, // combined correlation tile full height
					mcorr_comb_offset, // combined correlation tile height offset: 0 - centered (-height/2 to height/2), height/2 - only positive (0 to height)
					mcorr_comb_disp);
			pcombo_sels = new boolean [Correlation2d.MCORR_COMB.values().length][];
			if (imgdtt_params.mcorr_cons_all) {
				int indx = Correlation2d.MCORR_COMB.ALL.ordinal();
				pcombo_sels[indx] = correlation2d.selectAll(); 
			}
			if (imgdtt_params.mcorr_cons_dia) {
				int indx = Correlation2d.MCORR_COMB.DIA.ordinal();
				pcombo_sels[indx] = correlation2d.selectDiameters(null); 
			}
			if (imgdtt_params.mcorr_cons_sq) {
				int indx = Correlation2d.MCORR_COMB.SQ.ordinal();
				pcombo_sels[indx] = correlation2d.selectSquares(null); 
			}
			if (imgdtt_params.mcorr_cons_neib) {
				int indx = Correlation2d.MCORR_COMB.NEIB.ordinal();
				pcombo_sels[indx] = correlation2d.selectNeibs(null); 
			}
			if (imgdtt_params.mcorr_cons_hor) {
				int indx = Correlation2d.MCORR_COMB.HOR.ordinal();
				pcombo_sels[indx] = correlation2d.selectHorizontal(null); 
			}
			if (imgdtt_params.mcorr_cons_vert) {
				int indx = Correlation2d.MCORR_COMB.VERT.ordinal();
				pcombo_sels[indx] = correlation2d.selectVertical(null); 
			}
			if (clt_corr_out != null) {
				boolean [] calc_corr_pairs = correlation2d.getCorrPairs();
				for (int i = 0; i < calc_corr_pairs.length; i++) if (calc_corr_pairs[i]){
					clt_corr_out[i] = new double[tilesY][tilesX][]; 
				}
			}
			if (clt_combo_out != null) {
				for (int i = 0; i < pcombo_sels.length; i++) if (pcombo_sels[i] != null){
					clt_combo_out[i] = new double[tilesY][tilesX][]; 
				}
			}
			if (clt_combo_dbg != null) {
				boolean [] calc_corr_pairs = correlation2d.getCorrPairs();
				for (int i = 0; i < calc_corr_pairs.length; i++) if (calc_corr_pairs[i]){
					clt_combo_dbg[i] = new double[tilesY][tilesX][]; 
				}
			}
			
		}
		final boolean [][] combo_sels = pcombo_sels; 
		final boolean debug_distort= globalDebugLevel > 0; ///false; // true;
		final double [][] debug_offsets = null;

		final boolean macro_mode = macro_scale != 1;      // correlate tile data instead of the pixel data
		final int numcol = 3; // number of colors // keep the same, just do not use [0] and [1], [2] - green
		final int nTilesInChn=tilesX*tilesY;
		final double [][][][][][] clt_data = new double[numSensors][numcol][tilesY][tilesX][][];
		final Thread[] threads = newThreadArray(threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		final double [] col_weights= new double [numcol]; // colors are RBG
		final double [][] dbg_distort = debug_distort? (new double [4*numSensors][tilesX*tilesY]) : null;
		final double [][] corr_wnd = Corr2dLMA.getCorrWnd(
				transform_size,
				imgdtt_params.lma_wnd);
		final double [] corr_wnd_inv_limited = (imgdtt_params.lma_min_wnd <= 1.0)?  new double [corr_wnd.length * corr_wnd[0].length]: null;
		if (corr_wnd_inv_limited != null) {
			double inv_pwr = imgdtt_params.lma_wnd_pwr - (imgdtt_params.lma_wnd - 1.0); // compensate for lma_wnd
			for (int i = imgdtt_params.lma_hard_marg; i < (corr_wnd.length - imgdtt_params.lma_hard_marg); i++) {
				for (int j = imgdtt_params.lma_hard_marg; j < (corr_wnd.length - imgdtt_params.lma_hard_marg); j++) {
					corr_wnd_inv_limited[i * (corr_wnd.length) + j] = 1.0/Math.max(Math.pow(corr_wnd[i][j],inv_pwr), imgdtt_params.lma_min_wnd);
				}
			}
		}

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

		// TODO: move transpose_indices to Correlation2d class
		final int first_color = isMonochrome()? MONO_CHN : 0; // color that is non-zero
		

		if (globalDebugLevel > 0) {
			System.out.println("clt_aberrations_quad_corr(): width="+width+" height="+height+" transform_size="+transform_size+
					" debug_tileX="+debug_tileX+" debug_tileY="+debug_tileY+" globalDebugLevel="+globalDebugLevel);
		}
		final double [] filter =     doubleGetCltLpfFd(corr_sigma);
		

		if (disparity_map != null){
			for (int i = 0; i<disparity_map.length;i++){
				if (i == OVEREXPOSED) {
					if (saturation_imp!= null) {
						disparity_map[i] = new double [tilesY*tilesX];
					}
				} else if (isSliceBit(i)) {
					disparity_map[i] = new double [tilesY*tilesX];
				} else if (isDiffIndex(i)) {
					disparity_map[i] = new double [tilesY*tilesX];
				} else if (isToneRGBIndex(i)) {
					if (texture_tiles != null) { // for now - enable 12 tone layers only together with texture tiles
						disparity_map[i] = new double [tilesY*tilesX];
					}
				}
			}
		}
		DttRad2 dtt = new DttRad2(transform_size);
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
					double [][] fract_shiftsXY = new double[numSensors][];
					double [] ports_rgb = null;
					double [][] rXY;

					if (use_main) {
						rXY = geometryCorrection.getRXY(true); // boolean use_rig_offsets,
					} else  {
						rXY = geometryCorrection.getRXY(false); // boolean use_rig_offsets,
					}
					for (int nTile = ai.getAndIncrement(); nTile < nTilesInChn; nTile = ai.getAndIncrement()) {
						tileY = nTile /tilesX;
						tileX = nTile % tilesX;
						tIndex = tileY * tilesX + tileX;
						if (tile_op[tileY][tileX] == 0) continue; // nothing to do for this tile
						int                 img_mask =  getImgMask(tile_op[tileY][tileX]);         // which images to use
						if (numSensors > 4) {
							if (img_mask == 0xf) {
								for (int i = 0; i < numSensors; i++){
									img_mask |= 1 << i;
								}
							}
						}
						boolean debugTile =(tileX == debug_tileX) && (tileY == debug_tileY) && (globalDebugLevel > -1);
						boolean debugTile0 =(tileX == debug_tileX) && (tileY == debug_tileY) && (globalDebugLevel > -3);
						final int [] overexp_all = (saturation_imp != null) ? ( new int [2]): null;

						// Moved from inside chn loop
						centerX = tileX * transform_size + transform_size/2 - shiftX;
						centerY = tileY * transform_size + transform_size/2 - shiftY;
						// TODO: move port coordinates out of color channel loop
						double [][] centersXY;
						double [][] disp_dist = new double[numSensors][]; // used to correct 3D correlations

						if ((disparity_array == null) || (disparity_array[tileY] == null) || (Double.isNaN(disparity_array[tileY][tileX]))) {
							System.out.println("Bug with disparity_array !!! - 2, tileX="+tileX+", tileY="+tileY);
							continue; // nothing to do for this tile
						}

						if (macro_mode){
							if ((globalDebugLevel > -1) && debugTile) { // before correction
								System.out.println("\nUsing MACRO mode, centerX="+centerX+", centerY="+centerY);
							}
							centersXY = geometryCorrection.getPortsCoordinatesIdeal(
									macro_scale,
									centerX,
									centerY,
									macro_scale* disparity_array[tileY][tileX] + disparity_corr);

						} else {
							if (use_main) { // this is AUX camera that uses main coordinates  // not used in lwir
								centersXY =  geometryCorrection.getPortsCoordinatesAndDerivatives(
										geometryCorrection_main, //			GeometryCorrection gc_main,
										true,            // boolean use_rig_offsets,
										corr_rots,       // Matrix []   rots,
										null,            //  Matrix [][] deriv_rots,
										null,            // double [][] pXYderiv, // if not null, should be double[8][]
										disp_dist,       // used to correct 3D correlations
										centerX,
										centerY,
										disparity_array[tileY][tileX] + disparity_corr); // _aux); //  + disparity_corr);


							} else {  // used in lwir
								centersXY = geometryCorrection.getPortsCoordinatesAndDerivatives(
										geometryCorrection, //			GeometryCorrection gc_main,
										false,           // boolean use_rig_offsets,
										corr_rots,       // Matrix []   rots,
										null,            //  Matrix [][] deriv_rots,
										null,            // double [][] pXYderiv, // if not null, should be double[8][]
										disp_dist,       // used to correct 3D correlations
										centerX,
										centerY,
										disparity_array[tileY][tileX] + disparity_corr);
							}

							if (((globalDebugLevel > 0) || debug_distort) && debugTile) {
								for (int i = 0; i < numSensors; i++) {
									System.out.println("clt_aberrations_quad_corr():  tileX="+tileX+", tileY="+tileY+
											" centerX="+centerX+" centerY="+centerY+" disparity="+disparity_array[tileY][tileX]+
											" centersXY["+i+"][0]="+centersXY[i][0]+" centersXY["+i+"][1]="+centersXY[i][1]);
								}
							}
							if (debug_distort && debugTile && (debug_offsets != null)) {
								double [][] debug_offsets_xy = new double [debug_offsets.length][2];
								for (int i = 0; i < debug_offsets.length; i++) {
									debug_offsets_xy[i][0] = disp_dist[i][0] * debug_offsets[i][0] + disp_dist[i][1] * debug_offsets[i][1];
									debug_offsets_xy[i][1] = disp_dist[i][2] * debug_offsets[i][0] + disp_dist[i][3] * debug_offsets[i][1];
								}
								for (int i = 0; i < numSensors; i++) {
									System.out.println(String.format("%d: {%8.3f, %8.3f}",i,debug_offsets_xy[i][0],debug_offsets_xy[i][1]));
								}

								for (int i = 0; i < debug_offsets.length; i++) {
									centersXY[i][0] += debug_offsets_xy[i][0];
									centersXY[i][1] += debug_offsets_xy[i][1];
								}
								for (int i = 0; i < numSensors; i++) {
									System.out.println("Corrected clt_aberrations_quad_corr():  tileX="+tileX+", tileY="+tileY+
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

							for (int ip = 0; ip < centersXY.length; ip++){ // 4 OOB
								centersXY[ip][0] -= shiftXY[ip][0];
								centersXY[ip][1] -= shiftXY[ip][1];
							}
							// save disparity distortions for visualization:
							if (dbg_distort != null) {
								for (int cam = 0; cam <numSensors; cam++) {
									dbg_distort[cam * 4 + 0 ][nTile] = disp_dist[cam][0];
									dbg_distort[cam * 4 + 1 ][nTile] = disp_dist[cam][1];
									dbg_distort[cam * 4 + 2 ][nTile] = disp_dist[cam][2];
									dbg_distort[cam * 4 + 3 ][nTile] = disp_dist[cam][3];
								}
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
// See if macro_mode uses color channels for non-color?
						for (int ncol = 0; ncol <numcol; ncol++) {
							if (!isMonochrome() || (ncol == MONO_CHN) || macro_mode) { // in monochrome mode skip all non-mono (green) channels  // used in lwir (5 of 6 branches)
								if ((globalDebugLevel > -1) && (tileX == debug_tileX) && (tileY == debug_tileY) && (ncol == 2)) {
									System.out.println("\nUsing "+(macro_mode?"MACRO":"PIXEL")+" mode, centerX="+centerX+", centerY="+centerY);
									System.out.println(disparity_array[tileY][tileX]+"\t"+
											centersXY[0][0]+"\t"+centersXY[0][1]+"\t"+
											centersXY[1][0]+"\t"+centersXY[1][1]+"\t"+
											centersXY[2][0]+"\t"+centersXY[2][1]+"\t"+
											centersXY[3][0]+"\t"+centersXY[3][1]+"\t");
								}

								for (int i = 0; i < numSensors; i++) {
									clt_data[i][ncol][tileY][tileX] = new double [4][];
									// Extract image tiles and kernels, correct aberrations, return (ut do not apply) fractional shifts
									fract_shiftsXY[i] = extract_correct_tile( // return a pair of residual offsets
											image_data[i],
											width,       // image width
											((clt_kernels == null) ? null : clt_kernels[i]), // [color][tileY][tileX][band][pixel]
											clt_data[i][ncol][tileY][tileX], //double  [][]        clt_tile,    // should be double [4][];
											kernel_step,
											dtt,
											ncol,
											centersXY[i][0], // centerX, // center of aberration-corrected (common model) tile, X
											centersXY[i][1], // centerY, //
											((!FPGA_COMPARE_DATA && (globalDebugLevel > -1) && (tileX == debug_tileX) && (tileY == debug_tileY) && (ncol == 2) && (i==0)) ? (globalDebugLevel + 0) : 0), // external tile compare
											no_deconvolution,
											false, // ); // transpose);
											((saturation_imp != null) ? saturation_imp[i] : null), //final boolean [][]        saturation_imp, // (near) saturated pixels or null
											((saturation_imp != null) ? overexp_all: null)); // final double [] overexposed)
								} // for (int i = 0; i < quad; i++)
								if ((globalDebugLevel > -1) && (tileX == debug_tileX) && (tileY == debug_tileY) && (ncol == 2)) {
									System.out.println();
								}

								if (!no_fract_shift) {  // USED in lwir
									// apply residual shift
									for (int i = 0; i < numSensors; i++) {
										fract_shift(    // fractional shift in transform domain. Currently uses sin/cos - change to tables with 2? rotations
												clt_data[i][ncol][tileY][tileX], // double  [][]  clt_tile,
												fract_shiftsXY[i][0],            // double        shiftX,
												fract_shiftsXY[i][1],            // double        shiftY,
												//									(globalDebugLevel > 0) && (tileX == debug_tileX) && (tileY == debug_tileY)); // external tile compare
												((globalDebugLevel > 1) &&
														((ncol==0) || isMonochrome()) &&
														(tileX >= debug_tileX - 2) && (tileX <= debug_tileX + 2) &&
														(tileY >= debug_tileY - 2) && (tileY <= debug_tileY+2)));
									}
								}
							} else { // if (!isMonochrome() || (chn == MONO_CHN) || macro_mode) { // in monochrome mode skip all non-mono (green) channels
								for (int i = 0; i < numSensors; i++) {  // used in lwir
									clt_data[i][ncol] = null; // erase unused clt_data
								}
							}
						}// end of for (int chn = 0; chn <numcol; chn++)
						int tile_lma_debug_level =  ((tileX == debug_tileX) && (tileY == debug_tileY))? (imgdtt_params.lma_debug_level-1) : -2;
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
							// calculate overexposed fraction
							if (saturation_imp != null){
								disparity_map[OVEREXPOSED][nTile] = (1.0 * overexp_all[0]) / overexp_all[1];
							}
							// calculate all selected pairs correlations
							double [][]  corr_tiles = correlation2d.correlateCompositeFD(
						    		clt_data,                     // double [][][][][][] clt_data,
						    		tileX,                        // int                 tileX,
						    		tileY,                        // int                 tileY,
						    		correlation2d.getCorrPairs(), // boolean[]           pairs_mask,
						    		filter,                       // double []           lpf,
						    		getScaleStrengths(),          // double              scale_value, // scale correlation value
						    		col_weights,                  // double []           col_weights,
						    		corr_fat_zero);               // double              fat_zero)
							
							double [][] corr_combo = new double [combo_sels.length][];
							for (int nsel = 0; nsel < combo_sels.length; nsel++) if (combo_sels[nsel] != null){
								corr_combo[nsel] = correlation2d.accumulateInit();
								double sumw = correlation2d.accummulatePairs(
										corr_combo[nsel], // double []   accum_tile,
										corr_tiles,       // double [][] corr_tiles,
										combo_sels[nsel], // boolean []  selection,
							    		1.0);             // double      weight);
								correlation2d.normalizeAccumulatedPairs(
										corr_combo[nsel],
										sumw);
							}
							// save for output
							if (clt_corr_out != null) {
								for (int num_pair = 0; num_pair < corr_tiles.length; num_pair++) if (corr_tiles[num_pair] != null){
									clt_corr_out[num_pair][tileY][tileX] = corr_tiles[num_pair];
								}
							}
							if (clt_combo_out != null) {
								for (int num_comb = 0; num_comb < corr_combo.length; num_comb++) if (corr_combo[num_comb] != null){
									clt_combo_out[num_comb][tileY][tileX] = corr_combo[num_comb];
								}
							}
						    if (clt_combo_dbg != null) { // debug feature, will re-calculate scaled/rotated pairs
								for (int num_pair = 0; num_pair < corr_tiles.length; num_pair++) if (corr_tiles[num_pair] != null){
									clt_combo_dbg[num_pair][tileY][tileX] = correlation2d.accumulateInit();
									correlation2d.accummulatePair(
											clt_combo_dbg[num_pair][tileY][tileX], // double [] accum_tile,
											corr_tiles[num_pair], // double [] corr_tile,
											num_pair, // int       num_pair,
								    		1.0); // double    weight)
								}
						    }
							
						    // calculate CM maximums for all mixed channels
						    // First get integer correlation center, relative to the center
							double [] corr_combo_max = corr_combo[Correlation2d.MCORR_COMB.ALL.ordinal()];
							if (debugTile0) {
								System.out.println("tileX = "+tileX+", tileY = "+tileY);
							}
							double [] disp_str_int = new double[2];
							int [] ixy =   correlation2d.getMaxXYInt( // find integer pair or null if below threshold // USED in lwir
									corr_combo_max, // double [] data,      // [data_size * data_size]
									disp_str_int,
									correlation2d.getCombWidth(), //       data_width,
									correlation2d.getCombHeight()/2 - correlation2d.getCombOffset(), // int       center_row, ??????????????
									true, // boolean   axis_only,
									imgdtt_params.min_corr,  // double    minMax,    // minimal value to consider (at integer location, not interpolated)
									tile_lma_debug_level > 0); // boolean   debug);
							double [] corr_stat = null;
							

							// if integer argmax was strong enough, calculate CM argmax
							// will not fill out DISPARITY_INDEX_INT+1, DISPARITY_INDEX_CM+1, DISPARITY_INDEX_POLY+1
							// use clt_mismatch for that
							double [] disp_str = new double[2];
							if (ixy != null) { //TODO - for CM use magic!
								disp_str = disp_str_int;
								disparity_map[DISPARITY_INDEX_INT][tIndex] =      disp_str[0]; // -ixy[0];
								disparity_map[DISPARITY_INDEX_CM+1][tIndex] =     disp_str[1]; // reusing old fine-corr
								disparity_map[DISPARITY_STRENGTH_INDEX][tIndex] = disp_str[1];
								
								if (Double.isNaN(disparity_map[DISPARITY_STRENGTH_INDEX][tIndex])) {
									System.out.println("BUG: 1. disparity_map[DISPARITY_STRENGTH_INDEX]["+tIndex+"] should not be NaN");
								}
								corr_stat = correlation2d.getMaxXCm(         // get fractional center as a "center of mass" inside circle/square from the integer max
										corr_combo_max,                      // double [] data,      // [data_size * data_size]
										correlation2d.getCombWidth(),        // int       data_width,      //  = 2 * transform_size - 1;
										correlation2d.getCombHeight()/2 - correlation2d.getCombOffset(),// int       center_row,
										ixy[0],                              // int       ixcenter,  // integer center x
										// corr_wndy,                        // double [] window_y,  // (half) window function in y-direction(perpendicular to disparity: for row0  ==1
										// corr_wndx,                        // double [] window_x,  // half of a window function in x (disparity) direction
										(tile_lma_debug_level > 0)); // boolean   debug);
								if (corr_stat != null) {
									disparity_map[DISPARITY_INDEX_CM][tIndex] =      -corr_stat[0]; // -ixy[0];
									disparity_map[DISPARITY_INDEX_CM+1][tIndex] =     corr_stat[1]; // reusing old fine-corr
									disparity_map[DISPARITY_STRENGTH_INDEX][tIndex] = corr_stat[1];
								}
							}
							
							double [] corr_combo_hor =  corr_combo[Correlation2d.MCORR_COMB.HOR.ordinal()];
							double [] corr_combo_vert = corr_combo[Correlation2d.MCORR_COMB.VERT.ordinal()];
							if (imgdtt_params.pcorr_use_hv) { // use combined horizontal and vertical pairs (vertical not transposed) 
								// for compatibility with old code executed unconditionally. TODO: Move to if (corr_stat != null) ... condition below
								if ((corr_combo_hor != null) && (ixy != null)) {
									double [] hor_pair1 = correlation2d.getMaxXCmNotch(
											corr_combo_hor,                      // double [] data,      // [data_size * data_size]
											correlation2d.getCombWidth(),        // int       data_width,      //  = 2 * transform_size - 1;
											correlation2d.getCombHeight()/2 - correlation2d.getCombOffset(),// int       center_row,
											ixy[0],                              // int       ixcenter,  // integer center x
											(tile_lma_debug_level > 0)); // boolean   debug);
									if (hor_pair1 != null) {
										disparity_map[DISPARITY_INDEX_HOR][tIndex] =          -hor_pair1[0];
										disparity_map[DISPARITY_INDEX_HOR_STRENGTH][tIndex] =  hor_pair1[1];
									}

								}
								if ((corr_combo_vert != null) && (ixy != null)) {
									double [] vert_pair1 = correlation2d.getMaxXCmNotch(
											corr_combo_vert,                     // double [] data,      // [data_size * data_size]
											correlation2d.getCombWidth(),        // int       data_width,      //  = 2 * transform_size - 1;
											correlation2d.getCombHeight()/2 - correlation2d.getCombOffset(),// int       center_row,
											ixy[0],                              // int       ixcenter,  // integer center x
											(tile_lma_debug_level > 0)); // boolean   debug);
									if (vert_pair1 != null) {
										disparity_map[DISPARITY_INDEX_VERT][tIndex] =         -vert_pair1[0];
										disparity_map[DISPARITY_INDEX_VERT_STRENGTH][tIndex] = vert_pair1[1];
									}
								}
							}
							// proceed only if CM correlation result is non-null // for compatibility with old code we need it to run regardless of the strength of the normal correlation
							if (corr_stat != null) {
// skipping DISPARITY_VARIATIONS_INDEX - it was not used
								disp_str[0] = -corr_stat[0];
								disp_str[1] =  corr_stat[1];
//								disparity_map[DISPARITY_INDEX_CM][tIndex] = disp_str[0]; // disparity is negative X. Already done above
								if (tile_lma_debug_level > 0) {
									System.out.println("Will run getMaxXSOrtho( ) for tileX="+tileX+", tileY="+tileY);
								}

								// debug new LMA correlations
								if (debugTile0) { // should be debugTile
									System.out.println("Will run new LMA for tileX="+tileX+", tileY="+tileY);
								}
								double [] poly_disp = {Double.NaN, 0.0};
								Corr2dLMA lma2 = correlation2d.corrLMA2Single(
										imgdtt_params,                // ImageDttParameters  imgdtt_params,
										imgdtt_params.lmas_LY_single, // false,                        // boolean             adjust_ly, // adjust Lazy Eye
										corr_wnd,                     // double [][]         corr_wnd, // correlation window to save on re-calculation of the window
										corr_wnd_inv_limited,         // corr_wnd_limited, // correlation window, limited not to be smaller than threshold - used for finding max/convex areas (or null)
										corr_tiles, // corrs,          // double [][]         corrs,
										disp_dist,
										rXY,                          // double [][]         rXY, // non-distorted X,Y offset per nominal pixel of disparity
										// all that are not null in corr_tiles
										correlation2d.selectAll(),   // longToArray(imgdtt_params.dbg_pair_mask),  // int                 pair_mask, // which pairs to process
										disp_str,  //corr_stat[0],                 // double    xcenter,   // preliminary center x in pixels for largest baseline
										poly_disp,                    // double[]            poly_ds,    // null or pair of disparity/strength
										imgdtt_params.ortho_vasw_pwr, // double    vasw_pwr,  // value as weight to this power,
										tile_lma_debug_level, // +2,         // int                 debug_level,
										tileX,                        // int                 tileX, // just for debug output
										tileY );                      // int                 tileY
								if (debugTile0) { // should be debugTile
									System.out.println("Ran LMA for tileX="+tileX+", tileY="+tileY);
								}
								double [][] ds = null;
								if (lma2 != null) {
									ds = lma2.lmaDisparityStrength(
						    				imgdtt_params.lmas_min_amp,      //  minimal ratio of minimal pair correlation amplitude to maximal pair correlation amplitude
											imgdtt_params.lmas_max_rel_rms,  // maximal relative (to average max/min amplitude LMA RMS) // May be up to 0.3)
											imgdtt_params.lmas_min_strength, // minimal composite strength (sqrt(average amp squared over absolute RMS)
											imgdtt_params.lmas_min_ac,       // minimal of A and C coefficients maximum (measures sharpest point/line)
											imgdtt_params.lmas_min_min_ac,   // minimal of A and C coefficients minimum (measures sharpest point)
											imgdtt_params.lmas_max_area,      //double  lma_max_area,     // maximal half-area (if > 0.0)
											imgdtt_params.lma_str_scale,    // convert lma-generated strength to match previous ones - scale
											imgdtt_params.lma_str_offset    // convert lma-generated strength to match previous ones - add to result
											);
									if (ds != null) { // always true
										disparity_map[DISPARITY_INDEX_POLY][tIndex] =   ds[0][0];
										disparity_map[DISPARITY_INDEX_POLY+1][tIndex] = ds[0][1];
										disparity_map[DISPARITY_STRENGTH_INDEX][tIndex] = ds[0][1]; // overwrite with LMA strength

										if (debugTile0) {
											lma2.printStats(ds,1);
											if (imgdtt_params.lmas_LY_single) {
												double [][] ddnd = lma2.getDdNd();
												if (ddnd != null) {
													double [][] dxy= new double [ddnd.length][2];
													for (int i = 0; i < dxy.length; i++) {
														dxy[i][0] = ddnd[i][0] * rXY[i][0] - ddnd[i][1] * rXY[i][1];
														dxy[i][1] = ddnd[i][0] * rXY[i][1] + ddnd[i][1] * rXY[i][0];
													}
													System.out.print("       Port:  ");
													for (int i = 0; i < dxy.length; i++) System.out.print(String.format("   %2d   ", i)); System.out.println();
													System.out.print("Radial_in =  [");
													for (int i = 0; i < dxy.length; i++) System.out.print(String.format(" %6.3f,", ddnd[i][0])); System.out.println("]");
													System.out.print("Tangent_CW = [");
													for (int i = 0; i < dxy.length; i++) System.out.print(String.format(" %6.3f,", ddnd[i][1])); System.out.println("]");
													System.out.print("X =          [");
													for (int i = 0; i < dxy.length; i++) System.out.print(String.format(" %6.3f,", dxy[i][0])); System.out.println("]");
													System.out.print("Y =          [");
													for (int i = 0; i < dxy.length; i++) System.out.print(String.format(" %6.3f,", dxy[i][1])); System.out.println("]");
													System.out.println();
												} else {
													System.out.println("No dd/nd and x/y offsets data is available ");
												}
											} else {
												System.out.println("LY offsets are not measured");
											}
										}
									}
								}
//								} //if (debugTile0) { // should be debugTile
							} // end of if (corr_stat != null)
							
						} // if (disparity_map != null){ // not null - calculate correlations

						if (texture_tiles !=null) {
							if ((extra_disparity != 0) && !getForcedDisparity(tile_op[tileY][tileX])){ // 0 - adjust disparity, 1 - use provided
								// shift images by 0.5 * extra disparity in the diagonal direction  // not used in lwir
								for (int ncol = 0; ncol <numcol; ncol++) { // color
									for (int i = 0; i < numSensors; i++) {
										if (clt_data[i][ncol] != null) {
											fract_shift(    // fractional shift in transform domain. Currently uses sin/cos - change to tables with 2? rotations
													clt_data[i][ncol][tileY][tileX], // double  [][]  clt_tile,
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
							double [][][] iclt_tile = new double [numSensors][numcol][]; // in mono some may remain null
							double [] clt_tile;
							double scale = 0.25;  // matching iclt_2d
							for (int i = 0; i < numSensors; i++) {  // USED in lwir
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

										if ((globalDebugLevel > 10) && debugTile) {
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
							if ((globalDebugLevel > 0) && debugTile0) {
								ShowDoubleFloatArrays sdfa_instance = new ShowDoubleFloatArrays(); // just for debugging?
								String [] col_names= {"red","blue","green"};
								String [] titles = new String[numSensors * numcol];
//								String [] titles = {"red0","blue0","green0","red1","blue1","green1","red2","blue2","green2","red3","blue3","green3"};
								double [][] dbg_tile = new double [numSensors*numcol][];
								for (int i = 0; i < numSensors; i++) {
									for (int ncol = 0; ncol <numcol; ncol++) if (iclt_tile[i][ncol] != null) { // color
										String col_name = (ncol < col_names.length)?col_names[ncol]:("c<"+ncol+">");
										titles[i * numcol + ncol] = col_name + i; 
										dbg_tile[i * numcol + ncol] = iclt_tile[i][ncol];
									}
								}
								sdfa_instance.showArrays(dbg_tile, 2* transform_size, 2* transform_size, true, "iclt_x"+tileX+"_y"+tileY, titles);
							}


							// "de-bayer" tiles for matching, use original data for output
							double [][][] tiles_debayered = new double [numSensors][numcol][];
							for (int i =0; i<numSensors; i++){
								for (int ncol = 0; ncol < numcol; ncol++) if (iclt_tile[i][ncol] != null) {
									if (isMonochrome()) {  // used in lwir
										tiles_debayered[i][ncol] =  iclt_tile[i][ncol];
									} else {  // used in lwir
										tiles_debayered[i][ncol] =  tile_debayer_shot_corr(
												(ncol != 2), // red or blue (false - green)
												iclt_tile[i][ncol],
												2 * transform_size,
												lt_window2, // squared lapping window
												min_shot,   // 10.0;  // Do not adjust for shot noise if lower than
												scale_shot,  //3.0;   // scale when dividing by sqrt
												lt_window2, // re-apply window to the result
												false); //boolean  debug_gpu);

									}
								}
							}
							if ((globalDebugLevel > 0) && debugTile0) {
								ShowDoubleFloatArrays sdfa_instance = new ShowDoubleFloatArrays(); // just for debugging?
								String [] col_names= {"red","blue","green"};
								String [] titles = new String[numSensors * numcol];
								
								double [][] dbg_tile = new double [numSensors*numcol][];
								for (int i = 0; i < numSensors; i++) {
									for (int chn = 0; chn <numcol; chn++) { // color
										String col_name = (chn< col_names.length)?col_names[chn]:("c<"+chn+">");
										titles[i * numcol + chn] = col_name + i; 
										dbg_tile[i * numcol + chn] = tiles_debayered[i][chn];
									}
								}
								sdfa_instance.showArrays(dbg_tile, 2* transform_size, 2* transform_size, true, "tiles_debayered_x"+tileX+"_y"+tileY, titles);
							}
							// ... used in lwir
							double []     max_diff = null;
							if ((disparity_map != null) && (disparity_map.length >= (IMG_DIFF0_INDEX + numSensors))){
								max_diff = new double[numSensors];
							}
							int ports_rgb_len = numSensors*numcol;  // 12
							if ((disparity_map != null) && (disparity_map.length >= (getImgToneRGB() + ports_rgb_len))) {
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
									(globalDebugLevel > 0) && debugTile,
									false);          // boolean       debug_gpu)      // generate output fro matching with GPU processing

							if ((globalDebugLevel > 0) && debugTile0) {
								ShowDoubleFloatArrays sdfa_instance = new ShowDoubleFloatArrays(); // just for debugging?
								sdfa_instance.showArrays(texture_tiles[tileY][tileX], 2* transform_size, 2* transform_size, true, "tile_combine_rgba_x"+tileX+"_y"+tileY);
							}
							// mix RGB from iclt_tile, mix alpha with - what? correlation strength or 'don't care'? good correlation or all > min?
							for (int i = 0; i < iclt_tile[0][first_color].length; i++ ) {
								double sw = 0.0;
								for (int ip = 0; ip < numSensors; ip++) {
									sw += texture_tiles[tileY][tileX][numcol+1+ip][i];
								}
								if (sw != 0 ) sw = 1.0/sw;
								for (int ncol = 0; ncol < numcol; ncol++) if (iclt_tile[0][ncol] !=null){ // color
									texture_tiles[tileY][tileX][ncol][i] = 0.0; //iclt[tileY][tileX][chn]
									for (int ip = 0; ip < numSensors; ip++) {
										texture_tiles[tileY][tileX][ncol][i] += sw * texture_tiles[tileY][tileX][numcol+1+ip][i] * iclt_tile[ip][ncol][i];
									}
								}
							}
							if ((globalDebugLevel > 0) && debugTile0) { // same as previous for mono
								ShowDoubleFloatArrays sdfa_instance = new ShowDoubleFloatArrays(); // just for debugging?
								sdfa_instance.showArrays(texture_tiles[tileY][tileX], 2* transform_size, 2* transform_size, true, "tile_mixed_rgba_x"+tileX+"_y"+tileY);
							}
							if ((disparity_map != null) && (disparity_map.length >= (IMG_DIFF0_INDEX + numSensors))){
								for (int i = 0; i < max_diff.length; i++){
									disparity_map[IMG_DIFF0_INDEX + i][tIndex] = max_diff[i];
								}
							}
							if (ports_rgb != null) {
								for (int i = 0; i < ports_rgb.length; i++){
									disparity_map[getImgToneRGB() + i][tIndex] = ports_rgb[i];
								}
							}
						}
					}
				}
			};
		}
		startAndJoin(threads);

		if ((dbg_distort != null) &&(globalDebugLevel >=0)) {
			(new ShowDoubleFloatArrays()).showArrays(dbg_distort,  tilesX, tilesY, true, "disparity_distortions"); // , dbg_titles);
		}

		return clt_data;
	}

	// correlate multiple tiles around each selected one
	// generates non-zero disparities (cm and lma) pair with zero strength if that particular tile fails to provide argmax,
	// but the residual disparity still should be applied during refine sequence
	public double [][][][][][] clt_aberrations_quad_corr_multi(
			final ImageDttParameters  imgdtt_params,   // Now just extra correlation parameters, later will include, most others
			final int                 macro_scale,     // to correlate tile data instead of the pixel data: 1 - pixels, 8 - tiles
			final int [][]            tile_op,         // [tilesY][tilesX] - what to do - 0 - nothing for this tile
			final double [][]         disparity_array, // [tilesY][tilesX] - individual per-tile expected disparity
			final double [][][]       image_data, // first index - number of image in a quad
		    final boolean [][]        saturation_imp, // (near) saturated pixels or null
			 // correlation results
		    // clt_corr_out is needed to store intermediate correlation results
		    final double [][][][]     clt_corr_out0,   // sparse (by the first index) [type][tilesY][tilesX][(2*transform_size-1)*(2*transform_size-1)] or null
		    final double [][][][]     clt_combo_out,  // sparse (by the first index) [type][tilesY][tilesX][(combo_tile_size] or null
			final double [][][][]     clt_combo_dbg,  // generate sparse  partial rotated/scaled pairs
			// When clt_mismatch is non-zero, no far objects extraction will be attempted
			// clt_mismatch is used in older code, not supported in GPU - there is cltMeasureLazyEye for that purpose
			// this.correlation2d should be not null if disparity_map != null
			final double [][]         disparity_map,   // [8][tilesY][tilesX], only [6][] is needed on input or null - do not calculate
			                                           // last 2 - contrast, avg/ "geometric average)
			
			final int                 width,
			final double              corr_fat_zero,    // add to denominator to modify phase correlation (same units as data1, data2). <0 - pure sum
			final boolean             corr_sym,
			final double              corr_offset,
			final double              corr_red,
			final double              corr_blue,
			final double              corr_sigma,

			final GeometryCorrection  geometryCorrection,
			final GeometryCorrection  geometryCorrection_main, // if not null correct this camera (aux) to the coordinates of the main
			final double [][][][][][] clt_kernels, // [channel_in_quad][color][tileY][tileX][band][pixel] , size should match image (have 1 tile around)
			final int                 kernel_step,
			final int                 window_type,
			final double [][]         shiftXY, // [port]{shiftX,shiftY}
			final double              disparity_corr, // disparity at infinity
			final double [][][]       fine_corr, // quadratic coefficients for fine correction (or null)

			final double              shiftX, // shift image horizontally (positive - right) - just for testing
			final double              shiftY, // shift image vertically (positive - down)
			
			final int                 clustRadius,  // 1 - single tile, 2 - 3x3, 3 - 5x5, ...
			final double              arange, // absolute disparity range to consolidate
			final double              rrange, // relative disparity range to consolidate
			final int                 mcorr_sel,    // Which pairs to correlate // +1 - all, +2 - dia, +4 - sq, +8 - neibs, +16 - hor + 32 - vert
			final int                 mcorr_comb_width,  // combined correlation tile width
			final int                 mcorr_comb_height, // combined correlation tile full height
			final int                 mcorr_comb_offset, // combined correlation tile height offset: 0 - centered (-height/2 to height/2), height/2 - only positive (0 to height)
			final double              mcorr_comb_disp,   // Combined tile per-pixel disparity for baseline == side of a square
			
			final int                 debug_tileX,
			final int                 debug_tileY,
			final boolean             no_fract_shift,
			final boolean             no_deconvolution,
			final int                 threadsMax,  // maximal number of threads to launch
			final int                 globalDebugLevel)
	{
		final int clustDiameter = 2 * clustRadius - 1;
		final int center_indx = (clustRadius - 1) * (clustDiameter + 1);
		final double [] wnd_neib = new double [clustDiameter * clustDiameter];
		for (int iy = 0; iy < clustDiameter; iy++) {
			double wy = Math.sin(Math.PI *(iy+1) / (clustDiameter + 1));
			for (int ix = 0; ix < clustDiameter; ix++) {
				wnd_neib[iy * clustDiameter + ix] = wy * Math.sin(Math.PI *(ix+1) / (clustDiameter + 1));
			}
		}
		final int height=image_data[0][0].length/width;
		final int tilesX=width/transform_size;
		final int tilesY=height/transform_size;
		boolean [][] pcombo_sels = null;
		if (correlation2d == null) {
			  throw new IllegalArgumentException ("clt_aberrations_quad_corr_multi(): correlation2d == null!");
		}
		if (disparity_map == null) {
			  throw new IllegalArgumentException ("clt_aberrations_quad_corr_multi(): disparity_map == null!");
		}
		// Initialize correlation pairs selection to be used by all threads
		boolean [] corr_calculate = null;
		if (isCorrAll  (mcorr_sel)) corr_calculate = correlation2d.selectAll();
		if (isCorrDia  (mcorr_sel)) corr_calculate = correlation2d.selectDiameters  (corr_calculate);
		if (isCorrSq   (mcorr_sel)) corr_calculate = correlation2d.selectSquares    (corr_calculate);
		if (isCorrNeib (mcorr_sel)) corr_calculate = correlation2d.selectNeibs      (corr_calculate);
		if (isCorrHor  (mcorr_sel)) corr_calculate = correlation2d.selectHorizontal (corr_calculate);
		if (isCorrVert (mcorr_sel)) corr_calculate = correlation2d.selectVertical   (corr_calculate);
		correlation2d.setCorrPairs(corr_calculate); // will limit correlation pairs calculation
		correlation2d.generateResample( // should be called before
				mcorr_comb_width,  // combined correlation tile width
				mcorr_comb_height, // combined correlation tile full height
				mcorr_comb_offset, // combined correlation tile height offset: 0 - centered (-height/2 to height/2), height/2 - only positive (0 to height)
				mcorr_comb_disp);
		pcombo_sels = new boolean [Correlation2d.MCORR_COMB.values().length][];
		if (imgdtt_params.mcorr_cons_all) {
			int indx = Correlation2d.MCORR_COMB.ALL.ordinal();
			pcombo_sels[indx] = correlation2d.selectAll(); 
		}
		if (imgdtt_params.mcorr_cons_dia) {
			int indx = Correlation2d.MCORR_COMB.DIA.ordinal();
			pcombo_sels[indx] = correlation2d.selectDiameters(null); 
		}
		if (imgdtt_params.mcorr_cons_sq) {
			int indx = Correlation2d.MCORR_COMB.SQ.ordinal();
			pcombo_sels[indx] = correlation2d.selectSquares(null); 
		}
		if (imgdtt_params.mcorr_cons_neib) {
			int indx = Correlation2d.MCORR_COMB.NEIB.ordinal();
			pcombo_sels[indx] = correlation2d.selectNeibs(null); 
		}
		if (imgdtt_params.mcorr_cons_hor) {
			int indx = Correlation2d.MCORR_COMB.HOR.ordinal();
			pcombo_sels[indx] = correlation2d.selectHorizontal(null); 
		}
		if (imgdtt_params.mcorr_cons_vert) {
			int indx = Correlation2d.MCORR_COMB.VERT.ordinal();
			pcombo_sels[indx] = correlation2d.selectVertical(null); 
		}
		// if it was null, it will not be visible form outside, only used internally
		final double [][][][] clt_corr_out = (clt_corr_out0 != null) ? clt_corr_out0 : new double [correlation2d.getCorrPairs().length][][][]; 
		{
			boolean [] calc_corr_pairs = correlation2d.getCorrPairs();
			for (int i = 0; i < calc_corr_pairs.length; i++) if (calc_corr_pairs[i]){
				clt_corr_out[i] = new double[tilesY][tilesX][]; 
			}
		}

		if (clt_combo_out != null) {
			for (int i = 0; i < pcombo_sels.length; i++) if (pcombo_sels[i] != null){
				clt_combo_out[i] = new double[tilesY][tilesX][]; 
			}
		}
		if (clt_combo_dbg != null) {
			boolean [] calc_corr_pairs = correlation2d.getCorrPairs();
			for (int i = 0; i < calc_corr_pairs.length; i++) if (calc_corr_pairs[i]){
				clt_combo_dbg[i] = new double[tilesY][tilesX][]; 
			}
		}
		final boolean [][] combo_sels = pcombo_sels; 
		final boolean debug_distort= globalDebugLevel > 0; ///false; // true;
		final double [][] debug_offsets = null;
		final boolean macro_mode = macro_scale != 1;      // correlate tile data instead of the pixel data
		final int numcol = 3; // number of colors // keep the same, just do not use [0] and [1], [2] - green
		final int nTilesInChn=tilesX*tilesY;
		final double [][][][][][] clt_data = new double[numSensors][numcol][tilesY][tilesX][][];
		final Thread[] threads = newThreadArray(threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		final double [] col_weights= new double [numcol]; // colors are RBG
		final double [][] dbg_distort = debug_distort? (new double [4*numSensors][tilesX*tilesY]) : null;
		final double [][] corr_wnd = Corr2dLMA.getCorrWnd(
				transform_size,
				imgdtt_params.lma_wnd);
		final double [] corr_wnd_inv_limited = (imgdtt_params.lma_min_wnd <= 1.0)?  new double [corr_wnd.length * corr_wnd[0].length]: null;
		if (corr_wnd_inv_limited != null) {
			double inv_pwr = imgdtt_params.lma_wnd_pwr - (imgdtt_params.lma_wnd - 1.0); // compensate for lma_wnd
			for (int i = imgdtt_params.lma_hard_marg; i < (corr_wnd.length - imgdtt_params.lma_hard_marg); i++) {
				for (int j = imgdtt_params.lma_hard_marg; j < (corr_wnd.length - imgdtt_params.lma_hard_marg); j++) {
					corr_wnd_inv_limited[i * (corr_wnd.length) + j] = 1.0/Math.max(Math.pow(corr_wnd[i][j],inv_pwr), imgdtt_params.lma_min_wnd);
				}
			}
		}

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

		// TODO: move transpose_indices to Correlation2d class
		
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

		final double [] filter =     doubleGetCltLpfFd(corr_sigma);
		
		if (disparity_map != null){
			for (int i = 0; i<disparity_map.length;i++){
				if (i == OVEREXPOSED) {
					if (saturation_imp!= null) {
						disparity_map[i] = new double [tilesY*tilesX];
					}
				} else if (isSliceBit(i)) {
					disparity_map[i] = new double [tilesY*tilesX];
				} else if (isDiffIndex(i)) {
					disparity_map[i] = new double [tilesY*tilesX];
				}
			}
		}
		DttRad2 dtt = new DttRad2(transform_size);
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
					int tileY,tileX; // ,tIndex; // , chn;
					double centerX; // center of aberration-corrected (common model) tile, X
					double centerY; //
					double [][] fract_shiftsXY = new double[numSensors][];

					for (int nTile = ai.getAndIncrement(); nTile < nTilesInChn; nTile = ai.getAndIncrement()) {
						tileY = nTile /tilesX;
						tileX = nTile % tilesX;
						// now only procress selected tiles
						if (tile_op[tileY][tileX] == 0) { // see if it is still needed
							continue;
						}

						boolean debugTile =(tileX == debug_tileX) && (tileY == debug_tileY) && (globalDebugLevel > -1);
						boolean debugTile0 =(tileX == debug_tileX) && (tileY == debug_tileY) && (globalDebugLevel > -3);
						final int [] overexp_all = (saturation_imp != null) ? ( new int [2]): null;

						// Moved from inside chn loop
						centerX = tileX * transform_size + transform_size/2 - shiftX;
						centerY = tileY * transform_size + transform_size/2 - shiftY;
						// TODO: move port coordinates out of color channel loop
						double [][] centersXY;
						double [][] disp_dist = new double[numSensors][]; // used to correct 3D correlations

						if ((disparity_array == null) || (disparity_array[tileY] == null) || (Double.isNaN(disparity_array[tileY][tileX]))) {
							System.out.println("Bug with disparity_array !!! - 3");
							continue; // nothing to do for this tile
						}

						if (macro_mode){
							if ((globalDebugLevel > -1) && debugTile) { // before correction
								System.out.println("\nUsing MACRO mode, centerX="+centerX+", centerY="+centerY);
							}
							centersXY = geometryCorrection.getPortsCoordinatesIdeal(
									macro_scale,
									centerX,
									centerY,
									macro_scale* disparity_array[tileY][tileX] + disparity_corr);

						} else {
							if (use_main) { // this is AUX camera that uses main coordinates  // not used in lwir
								centersXY =  geometryCorrection.getPortsCoordinatesAndDerivatives(
										geometryCorrection_main, //			GeometryCorrection gc_main,
										true,            // boolean use_rig_offsets,
										corr_rots,       // Matrix []   rots,
										null,            //  Matrix [][] deriv_rots,
										null,            // double [][] pXYderiv, // if not null, should be double[8][]
										disp_dist,       // used to correct 3D correlations
										centerX,
										centerY,
										disparity_array[tileY][tileX] + disparity_corr); // _aux); //  + disparity_corr);


							} else {  // used in lwir
								centersXY = geometryCorrection.getPortsCoordinatesAndDerivatives(
										geometryCorrection, //			GeometryCorrection gc_main,
										false,           // boolean use_rig_offsets,
										corr_rots,       // Matrix []   rots,
										null,            //  Matrix [][] deriv_rots,
										null,            // double [][] pXYderiv, // if not null, should be double[8][]
										disp_dist,       // used to correct 3D correlations
										centerX,
										centerY,
										disparity_array[tileY][tileX] + disparity_corr);
							}

							if (((globalDebugLevel > 0) || debug_distort) && debugTile) {
								for (int i = 0; i < numSensors; i++) {
									System.out.println("clt_aberrations_quad_corr():  tileX="+tileX+", tileY="+tileY+
											" centerX="+centerX+" centerY="+centerY+" disparity="+disparity_array[tileY][tileX]+
											" centersXY["+i+"][0]="+centersXY[i][0]+" centersXY["+i+"][1]="+centersXY[i][1]);
								}
							}
							if (debug_distort && debugTile && (debug_offsets != null)) {
								double [][] debug_offsets_xy = new double [debug_offsets.length][2];
								for (int i = 0; i < debug_offsets.length; i++) {
									debug_offsets_xy[i][0] = disp_dist[i][0] * debug_offsets[i][0] + disp_dist[i][1] * debug_offsets[i][1];
									debug_offsets_xy[i][1] = disp_dist[i][2] * debug_offsets[i][0] + disp_dist[i][3] * debug_offsets[i][1];
								}
								for (int i = 0; i < numSensors; i++) {
									System.out.println(String.format("%d: {%8.3f, %8.3f}",i,debug_offsets_xy[i][0],debug_offsets_xy[i][1]));
								}

								for (int i = 0; i < debug_offsets.length; i++) {
									centersXY[i][0] += debug_offsets_xy[i][0];
									centersXY[i][1] += debug_offsets_xy[i][1];
								}
								for (int i = 0; i < numSensors; i++) {
									System.out.println("Corrected clt_aberrations_quad_corr():  tileX="+tileX+", tileY="+tileY+
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

							for (int ip = 0; ip < centersXY.length; ip++){ // 4 OOB
								centersXY[ip][0] -= shiftXY[ip][0];
								centersXY[ip][1] -= shiftXY[ip][1];
							}
							// save disparity distortions for visualization:
							if (dbg_distort != null) {
								for (int cam = 0; cam <numSensors; cam++) {
									dbg_distort[cam * 4 + 0 ][nTile] = disp_dist[cam][0];
									dbg_distort[cam * 4 + 1 ][nTile] = disp_dist[cam][1];
									dbg_distort[cam * 4 + 2 ][nTile] = disp_dist[cam][2];
									dbg_distort[cam * 4 + 3 ][nTile] = disp_dist[cam][3];
								}
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
						if (FPGA_COMPARE_DATA && (globalDebugLevel > 0) && (tileX == debug_tileX) && (tileY == debug_tileY)) { // not used in lwir
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
//									transform_size,       // final int                 transform_size,
									width,                 // final int                 width
									dtt
									);
						}
// See if macro_mode uses color channels for non-color?
						for (int ncol = 0; ncol <numcol; ncol++) {
							if (!isMonochrome() || (ncol == MONO_CHN) || macro_mode) { // in monochrome mode skip all non-mono (green) channels  // used in lwir (5 of 6 branches)
								boolean debug_for_fpga = FPGA_COMPARE_DATA && (globalDebugLevel > 0) && (tileX == debug_tileX) && (tileY == debug_tileY) && (ncol == 2);
								if ((globalDebugLevel > -1) && (tileX == debug_tileX) && (tileY == debug_tileY) && (ncol == 2)) {
									System.out.println("\nUsing "+(macro_mode?"MACRO":"PIXEL")+" mode, centerX="+centerX+", centerY="+centerY);
									System.out.println(disparity_array[tileY][tileX]+"\t"+
											centersXY[0][0]+"\t"+centersXY[0][1]+"\t"+
											centersXY[1][0]+"\t"+centersXY[1][1]+"\t"+
											centersXY[2][0]+"\t"+centersXY[2][1]+"\t"+
											centersXY[3][0]+"\t"+centersXY[3][1]+"\t");
								}

								for (int i = 0; i < numSensors; i++) {
									if (debug_for_fpga && (i==0)){ // not used in lwir
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

										fpga_chn =           2;
										System.out.println(String.format("Manually changing offset: center X= %f, center Y = %f", fpga_centersXY[0],fpga_centersXY[1]));
										System.out.println(String.format("Manually changing color to %d (was %d)", fpga_chn, ncol));
										fpga_fract_shiftsXY = extract_correct_tile( // return a pair of residual offsets
												image_data[i],
												width,       // image width
												null,
												fpga_clt_data, //double  [][]        clt_tile,    // should be double [4][];
												kernel_step,
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
									// Extract image tiles and kernels, correct aberrations, return (ut do not apply) fractional shifts
									fract_shiftsXY[i] = extract_correct_tile( // return a pair of residual offsets
											image_data[i],
											width,       // image width
											((clt_kernels == null) ? null : clt_kernels[i]), // [color][tileY][tileX][band][pixel]
											clt_data[i][ncol][tileY][tileX], //double  [][]        clt_tile,    // should be double [4][];
											kernel_step,
											dtt,
											ncol,
											centersXY[i][0], // centerX, // center of aberration-corrected (common model) tile, X
											centersXY[i][1], // centerY, //
											((!FPGA_COMPARE_DATA && (globalDebugLevel > -1) && (tileX == debug_tileX) && (tileY == debug_tileY) && (ncol == 2) && (i==0)) ? (globalDebugLevel + 0) : 0), // external tile compare
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
									for (int i = 0; i < numSensors; i++) {
										System.out.println("clt_aberrations_quad(): color="+ncol+", tileX="+tileX+", tileY="+tileY+
												" fract_shiftsXY["+i+"][0]="+fract_shiftsXY[i][0]+" fract_shiftsXY["+i+"][1]="+fract_shiftsXY[i][1]);
									}
								}

								if (!no_fract_shift) {  // USED in lwir
									// apply residual shift
									for (int i = 0; i < numSensors; i++) {
										fract_shift(    // fractional shift in transform domain. Currently uses sin/cos - change to tables with 2? rotations
												clt_data[i][ncol][tileY][tileX], // double  [][]  clt_tile,
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
								for (int i = 0; i < numSensors; i++) {  // used in lwir
									clt_data[i][ncol] = null; // erase unused clt_data
								}
							}
						}// end of for (int chn = 0; chn <numcol; chn++)
						// all color channels are done here
						// fill clt_corr_combo if it exists
						if (debugTile0) {
							System.out.println("tileX = "+tileX+", tileY = "+tileY);
						}
						if (disparity_map != null){ // not null - calculate correlations 
							for (int i = 0; i < disparity_map.length; i++) {
								if (disparity_map[i] != null) disparity_map[i][nTile] = (
										(i == DISPARITY_STRENGTH_INDEX) ||
										(i == DISPARITY_INDEX_HOR_STRENGTH) ||
										(i == DISPARITY_INDEX_VERT_STRENGTH)) ? 0.0 : Double.NaN; // once and for all
							}
							//clt_mismatch should only be used with disparity_map != null;
							// calculate overexposed fraction
							if (saturation_imp != null){
								disparity_map[OVEREXPOSED][nTile] = (1.0 * overexp_all[0]) / overexp_all[1];
							}
							// calculate all selected pairs correlations
							double [][]  corr_tiles = correlation2d.correlateCompositeFD(
									clt_data,                     // double [][][][][][] clt_data,
									tileX,                        // int                 tileX,
									tileY,                        // int                 tileY,
									correlation2d.getCorrPairs(), // boolean[]           pairs_mask,
									filter,                       // double []           lpf,
									getScaleStrengths(),          // double              scale_value, // scale correlation value
									col_weights,                  // double []           col_weights,
									corr_fat_zero);               // double              fat_zero)
							for (int num_pair = 0; num_pair < corr_tiles.length; num_pair++) if (corr_tiles[num_pair] != null){
								clt_corr_out[num_pair][tileY][tileX] = corr_tiles[num_pair];
							}
						}
					}
				}
			};
		}
		startAndJoin(threads);

		final boolean [][] used_neibs = new boolean [tilesX*tilesY][];
		final double [] disp_lma =      new double [tilesX*tilesY];
		final double [] str_lma =       new double [tilesX*tilesY];
		final double [] disp_cm =       new double [tilesX*tilesY];
		final double [] str_cm =        new double [tilesX*tilesY];

		ai.set(0);
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				@Override
				public void run() {
					DttRad2 dtt = new DttRad2(transform_size);
					dtt.set_window(window_type);
					int tileY,tileX; // , chn;
					double centerX; // center of aberration-corrected (common model) tile, X
					double centerY; //
					double [][] rXY;
					TileNeibs tn = new TileNeibs(tilesX,tilesY);
					if (use_main) {
						rXY = geometryCorrection.getRXY(true); // boolean use_rig_offsets,
					} else  {
						rXY = geometryCorrection.getRXY(false); // boolean use_rig_offsets,
					}
					for (int nTile = ai.getAndIncrement(); nTile < nTilesInChn; nTile = ai.getAndIncrement()) {
						tileY = nTile /tilesX;
						tileX = nTile % tilesX;
						if (tile_op[tileY][tileX] == 0) continue; // nothing to do for this tile
						boolean debugTile =(tileX == debug_tileX) && (tileY == debug_tileY) && (globalDebugLevel > -1);
						boolean debugTile0 =(tileX == debug_tileX) && (tileY == debug_tileY) && (globalDebugLevel > -3);
						centerX = tileX * transform_size + transform_size/2 - shiftX;
						centerY = tileY * transform_size + transform_size/2 - shiftY;
						// TODO: move port coordinates out of color channel loop
//						double [][] centersXY;
						double [][] disp_dist = new double[numSensors][]; // used to correct 3D correlations
						int tile_lma_debug_level =  ((tileX == debug_tileX) && (tileY == debug_tileY))? (imgdtt_params.lma_debug_level-1) : -2;
						
						if ((disparity_array == null) || (disparity_array[tileY] == null) || (Double.isNaN(disparity_array[tileY][tileX]))) {
							System.out.println("Bug with disparity_array !!! - 4");
							continue; // nothing to do for this tile
						}

						if (macro_mode){
							if ((globalDebugLevel > -1) && debugTile) { // before correction
								System.out.println("\nUsing MACRO mode, centerX="+centerX+", centerY="+centerY);
							}
//							centersXY = 
									geometryCorrection.getPortsCoordinatesIdeal(
									macro_scale,
									centerX,
									centerY,
									macro_scale* disparity_array[tileY][tileX] + disparity_corr);

						} else {
							if (use_main) { // this is AUX camera that uses main coordinates  // not used in lwir
//								centersXY =  
										geometryCorrection.getPortsCoordinatesAndDerivatives(
										geometryCorrection_main, //			GeometryCorrection gc_main,
										true,            // boolean use_rig_offsets,
										corr_rots,       // Matrix []   rots,
										null,            //  Matrix [][] deriv_rots,
										null,            // double [][] pXYderiv, // if not null, should be double[8][]
										disp_dist,       // used to correct 3D correlations
										centerX,
										centerY,
										disparity_array[tileY][tileX] + disparity_corr); // _aux); //  + disparity_corr);


							} else {  // used in lwir
//								centersXY = 
										geometryCorrection.getPortsCoordinatesAndDerivatives(
										geometryCorrection, //			GeometryCorrection gc_main,
										false,           // boolean use_rig_offsets,
										corr_rots,       // Matrix []   rots,
										null,            //  Matrix [][] deriv_rots,
										null,            // double [][] pXYderiv, // if not null, should be double[8][]
										disp_dist,       // used to correct 3D correlations
										centerX,
										centerY,
										disparity_array[tileY][tileX] + disparity_corr);
							}
                        }
						// determine which of the neighbors to use
						used_neibs[nTile] = new boolean [clustDiameter * clustDiameter];
						double disparity0 = disparity_array[tileY][tileX];
						double disp_min = disparity0 - arange - rrange * Math.abs(disparity0);
						double disp_max = disparity0 + arange + rrange * Math.abs(disparity0);
						double sum_w = 0.0;
						for (int dty = -clustRadius+1; dty < clustRadius; dty++) {
							for (int dtx = -clustRadius+1; dtx < clustRadius; dtx++) {
								int nTile1 = tn.getNeibIndex(nTile, dtx, dty);
								if (nTile1 >= 0){
									int tileX1 = nTile1 % tilesX;
									int tileY1 = nTile1 / tilesX;
									if ((tile_op[tileY1][tileX1] != 0) && (disparity_array[tileY1][tileX1] >= disp_min)  && (disparity_array[tileY1][tileX1] <= disp_max)){
										int windx = (dty * clustDiameter + dtx + center_indx);
										used_neibs[nTile][windx] = true;
										sum_w += wnd_neib[windx];
									}
								}
							}
						}
						double [][]  corr_tiles = new double [clt_corr_out.length][];
						if (sum_w > 0) {
							double s = 1.0/sum_w;
							// combine correlation tiles
							for (int num_pair = 0; num_pair < corr_tiles.length; num_pair ++) if (clt_corr_out[num_pair] != null) {
								corr_tiles[num_pair] = new double[(2 * transform_size - 1) * (2 * transform_size - 1)];
								for (int dty = -clustRadius+1; dty < clustRadius; dty++) {
									for (int dtx = -clustRadius+1; dtx < clustRadius; dtx++) {
										int nTile1 = tn.getNeibIndex(nTile, dtx, dty);
										if (nTile1 >= 0){
											int tileX1 = nTile1 % tilesX;
											int tileY1 = nTile1 / tilesX;
											if (clt_corr_out[num_pair][tileY1][tileX1] != null) {
												int windx = (dty * clustDiameter + dtx + center_indx);
												double w = s * wnd_neib[windx];
												for (int i = 0; i < corr_tiles[num_pair].length; i++) {
													corr_tiles[num_pair][i] += w * clt_corr_out[num_pair][tileY1][tileX1][i];
												}
											}

										}

									}
								}
							}
						}

						double [][] corr_combo = new double [combo_sels.length][];
						for (int nsel = 0; nsel < combo_sels.length; nsel++) if (combo_sels[nsel] != null){
							corr_combo[nsel] = correlation2d.accumulateInit();
							double sumw = correlation2d.accummulatePairs(
									corr_combo[nsel], // double []   accum_tile,
									corr_tiles,       // double [][] corr_tiles,
									combo_sels[nsel], // boolean []  selection,
									1.0);             // double      weight);
							correlation2d.normalizeAccumulatedPairs(
									corr_combo[nsel],
									sumw);
						}
						// save for output
						if (clt_corr_out != null) {
							for (int num_pair = 0; num_pair < corr_tiles.length; num_pair++) if (corr_tiles[num_pair] != null){
								clt_corr_out[num_pair][tileY][tileX] = corr_tiles[num_pair];
							}
						}
						if (clt_combo_out != null) {
							for (int num_comb = 0; num_comb < corr_combo.length; num_comb++) if (corr_combo[num_comb] != null){
								clt_combo_out[num_comb][tileY][tileX] = corr_combo[num_comb];
							}
						}
						if (clt_combo_dbg != null) { // debug feature, will re-calculate scaled/rotated pairs
							for (int num_pair = 0; num_pair < corr_tiles.length; num_pair++) if (corr_tiles[num_pair] != null){
								clt_combo_dbg[num_pair][tileY][tileX] = correlation2d.accumulateInit();
								correlation2d.accummulatePair(
										clt_combo_dbg[num_pair][tileY][tileX], // double [] accum_tile,
										corr_tiles[num_pair], // double [] corr_tile,
										num_pair, // int       num_pair,
										1.0); // double    weight)
							}
						}

						// calculate CM maximums for all mixed channels
						// First get integer correlation center, relative to the center
						double [] corr_combo_max = corr_combo[Correlation2d.MCORR_COMB.ALL.ordinal()];
						if (debugTile0) {
							System.out.println("tileX = "+tileX+", tileY = "+tileY);
						}
						double [] disp_str_int = new double[2];
						int [] ixy =   correlation2d.getMaxXYInt( // find integer pair or null if below threshold // USED in lwir
								corr_combo_max, // double [] data,      // [data_size * data_size]
								disp_str_int,
								correlation2d.getCombWidth(), //       data_width,
								correlation2d.getCombHeight()/2 - correlation2d.getCombOffset(), // int       center_row, ??????????????
								true, // boolean   axis_only,
								imgdtt_params.min_corr,  // double    minMax,    // minimal value to consider (at integer location, not interpolated)
								tile_lma_debug_level > 0); // boolean   debug);
						double [] corr_stat = null;


						// if integer argmax was strong enough, calculate CM argmax
						// will not fill out DISPARITY_INDEX_INT+1, DISPARITY_INDEX_CM+1, DISPARITY_INDEX_POLY+1
						// use clt_mismatch for that
						double [] disp_str = new double[2];
						if (ixy != null) { //TODO - for CM use magic!
							disp_str = disp_str_int;
							corr_stat = correlation2d.getMaxXCm(         // get fractional center as a "center of mass" inside circle/square from the integer max
									corr_combo_max,                      // double [] data,      // [data_size * data_size]
									correlation2d.getCombWidth(),        // int       data_width,      //  = 2 * transform_size - 1;
									correlation2d.getCombHeight()/2 - correlation2d.getCombOffset(),// int       center_row,
									ixy[0],                              // int       ixcenter,  // integer center x
									// corr_wndy,                        // double [] window_y,  // (half) window function in y-direction(perpendicular to disparity: for row0  ==1
									// corr_wndx,                        // double [] window_x,  // half of a window function in x (disparity) direction
									(tile_lma_debug_level > 0)); // boolean   debug);
							if (corr_stat != null) {
								disp_cm[nTile] =    -corr_stat[0]; // -ixy[0];
								str_cm [nTile] =     corr_stat[1]; // reusing old fine-corr
							}
						}

						// proceed only if CM correlation result is non-null // for compatibility with old code we need it to run regardless of the strength of the normal correlation
						if (corr_stat != null) {
							// skipping DISPARITY_VARIATIONS_INDEX - it was not used
							disp_str[0] = -corr_stat[0];
							disp_str[1] =  corr_stat[1];
							//								disparity_map[DISPARITY_INDEX_CM][tIndex] = disp_str[0]; // disparity is negative X. Already done above
							if (tile_lma_debug_level > 0) {
								System.out.println("Will run getMaxXSOrtho( ) for tileX="+tileX+", tileY="+tileY);
							}

							// debug new LMA correlations
							if (debugTile0) { // should be debugTile
								System.out.println("Will run new LMA for tileX="+tileX+", tileY="+tileY);
							}
							double [] poly_disp = {Double.NaN, 0.0};
							Corr2dLMA lma2 = correlation2d.corrLMA2Single(
									imgdtt_params,                // ImageDttParameters  imgdtt_params,
									imgdtt_params.lmas_LY_single, // false,                        // boolean             adjust_ly, // adjust Lazy Eye
									corr_wnd,                     // double [][]         corr_wnd, // correlation window to save on re-calculation of the window
									corr_wnd_inv_limited,         // corr_wnd_limited, // correlation window, limited not to be smaller than threshold - used for finding max/convex areas (or null)
									corr_tiles, // corrs,          // double [][]         corrs,
									disp_dist,
									rXY,                          // double [][]         rXY, // non-distorted X,Y offset per nominal pixel of disparity
									// all that are not null in corr_tiles
									correlation2d.selectAll(),   // longToArray(imgdtt_params.dbg_pair_mask),  // int                 pair_mask, // which pairs to process
									disp_str,  //corr_stat[0],                 // double    xcenter,   // preliminary center x in pixels for largest baseline
									poly_disp,                    // double[]            poly_ds,    // null or pair of disparity/strength
									imgdtt_params.ortho_vasw_pwr, // double    vasw_pwr,  // value as weight to this power,
									tile_lma_debug_level, // +2,         // int                 debug_level,
									tileX,                        // int                 tileX, // just for debug output
									tileY );                      // int                 tileY
							if (debugTile0) { // should be debugTile
								System.out.println("Ran LMA for tileX="+tileX+", tileY="+tileY);
							}
							double [][] ds = null;
							if (lma2 != null) {
								ds = lma2.lmaDisparityStrength(
										imgdtt_params.lmas_min_amp,      //  minimal ratio of minimal pair correlation amplitude to maximal pair correlation amplitude
										imgdtt_params.lmas_max_rel_rms,  // maximal relative (to average max/min amplitude LMA RMS) // May be up to 0.3)
										imgdtt_params.lmas_min_strength, // minimal composite strength (sqrt(average amp squared over absolute RMS)
										imgdtt_params.lmas_min_ac,       // minimal of A and C coefficients maximum (measures sharpest point/line)
										imgdtt_params.lmas_min_min_ac,   // minimal of A and C coefficients minimum (measures sharpest point)
										imgdtt_params.lmas_max_area,      //double  lma_max_area,     // maximal half-area (if > 0.0)
										imgdtt_params.lma_str_scale,    // convert lma-generated strength to match previous ones - scale
										imgdtt_params.lma_str_offset    // convert lma-generated strength to match previous ones - add to result
										);
								if (ds != null) { // always true
									disp_lma[nTile] =   ds[0][0];
									str_lma[nTile] =    ds[0][1];
									if (debugTile0) {
										lma2.printStats(ds,1);
										if (imgdtt_params.lmas_LY_single) {
											double [][] ddnd = lma2.getDdNd();
											if (ddnd != null) {
												double [][] dxy= new double [ddnd.length][2];
												for (int i = 0; i < dxy.length; i++) {
													dxy[i][0] = ddnd[i][0] * rXY[i][0] - ddnd[i][1] * rXY[i][1];
													dxy[i][1] = ddnd[i][0] * rXY[i][1] + ddnd[i][1] * rXY[i][0];
												}
												System.out.print("       Port:  ");
												for (int i = 0; i < dxy.length; i++) System.out.print(String.format("   %2d   ", i)); System.out.println();
												System.out.print("Radial_in =  [");
												for (int i = 0; i < dxy.length; i++) System.out.print(String.format(" %6.3f,", ddnd[i][0])); System.out.println("]");
												System.out.print("Tangent_CW = [");
												for (int i = 0; i < dxy.length; i++) System.out.print(String.format(" %6.3f,", ddnd[i][1])); System.out.println("]");
												System.out.print("X =          [");
												for (int i = 0; i < dxy.length; i++) System.out.print(String.format(" %6.3f,", dxy[i][0])); System.out.println("]");
												System.out.print("Y =          [");
												for (int i = 0; i < dxy.length; i++) System.out.print(String.format(" %6.3f,", dxy[i][1])); System.out.println("]");
												System.out.println();
											} else {
												System.out.println("No dd/nd and x/y offsets data is available ");
											}
										} else {
											System.out.println("LY offsets are not measured");
										}
									}
								}
							}
							//								} //if (debugTile0) { // should be debugTile
						} // end of if (corr_stat != null)
						// REMOVED if (texture_tiles !=null)
					}
				}
			};
		}
		startAndJoin(threads);
		
		// third pass apply corrections to neighbors, even if they did not produce argmax-es 
		ai.set(0);
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				@Override
				public void run() {
					DttRad2 dtt = new DttRad2(transform_size);
					dtt.set_window(window_type);
					int tileY,tileX; // , chn;
					TileNeibs tn = new TileNeibs(tilesX,tilesY);
					for (int nTile = ai.getAndIncrement(); nTile < nTilesInChn; nTile = ai.getAndIncrement()) {
						tileY = nTile /tilesX;
						tileX = nTile % tilesX;
						if (tile_op[tileY][tileX] == 0) continue; // nothing to do for this tile
						double sum_w_lma = 0.0;
						double sum_w_cm = 0.0;
						for (int dty = -clustRadius+1; dty < clustRadius; dty++) {
							for (int dtx = -clustRadius+1; dtx < clustRadius; dtx++) {
								int windx = (center_indx - dty * clustDiameter - dtx);
								if (used_neibs[nTile][windx]) {
									double w = wnd_neib[windx];
									int nTile1 = tn.getNeibIndex(nTile, -dtx, -dty); // no need to check for >=0
									sum_w_lma += (str_lma[nTile1] > 0) ? w : 0;
									sum_w_cm +=  (str_cm[nTile1] > 0) ? w : 0;
								}
							}
						}
						double s_lma = (sum_w_lma > 0.0) ? (1.0/sum_w_lma) : 0.0;
						double s_cm =  (sum_w_cm > 0.0) ?  (1.0/sum_w_cm) : 0.0;
						double d_lma = 0.0;
						double d_cm = 0.0;
						if ((sum_w_lma > 0.0) || (sum_w_cm > 0.0)) {
							for (int dty = -clustRadius+1; dty < clustRadius; dty++) {
								for (int dtx = -clustRadius+1; dtx < clustRadius; dtx++) {
									int windx = (center_indx - dty * clustDiameter - dtx);
									if (used_neibs[nTile][windx]) {
										double w = wnd_neib[windx];
										int nTile1 = tn.getNeibIndex(nTile, -dtx, -dty); // no need to check for >=0
										d_lma += (str_lma[nTile1] > 0) ? (s_lma * w * disp_lma[nTile1]) : 0;
										d_cm +=  (str_cm[nTile1] > 0) ?  (s_cm * w * disp_cm[nTile1]) : 0;
									}
								}
							}
						}
						disparity_map[DISPARITY_INDEX_CM+1][nTile] =     str_cm[nTile];
						disparity_map[DISPARITY_INDEX_POLY+1][nTile] =   str_lma[nTile];
						disparity_map[DISPARITY_STRENGTH_INDEX][nTile] = str_cm[nTile];
						if (str_lma[nTile] > 0) {
							disparity_map[DISPARITY_STRENGTH_INDEX][nTile] = str_lma[nTile];
						}
						if (sum_w_cm > 0.0) {
							disparity_map[DISPARITY_INDEX_POLY][nTile] =   d_cm;
						} else {
							disparity_map[DISPARITY_INDEX_CM][nTile] =   Double.NaN;
						}
						if (sum_w_lma > 0.0) {
							disparity_map[DISPARITY_INDEX_POLY][nTile] =   d_lma;
						} else {
							disparity_map[DISPARITY_INDEX_POLY][nTile] =   Double.NaN;
						}
					}
				}
			};
		}
		startAndJoin(threads);
		
		if ((dbg_distort != null) &&(globalDebugLevel >=0)) {
			(new ShowDoubleFloatArrays()).showArrays(dbg_distort,  tilesX, tilesY, true, "disparity_distortions"); // , dbg_titles);
		}
		return clt_data;
	}

	
//	public double [][][][][][] clt_aberrations_quad_corr_tilted(
	public void clt_aberrations_quad_corr_tilted(
			final ImageDttParameters  imgdtt_params,   // Now just extra correlation parameters, later will include, most others
			final int                 macro_scale,     // to correlate tile data instead of the pixel data: 1 - pixels, 8 - tiles
			final int [][]            tile_op,         // [tilesY][tilesX] - what to do - 0 - nothing for this tile
			final double [][]         disparity_array, // [tilesY][tilesX] - individual per-tile expected disparity
			final double [][][]       image_data, // first index - number of image in a quad
		    final boolean [][]        saturation_imp, // (near) saturated pixels or null
			 // correlation results
		    // clt_corr_out is needed to store intermediate correlation results
		    final double [][][][]     clt_corr_out0,   // sparse (by the first index) [type][tilesY][tilesX][(2*transform_size-1)*(2*transform_size-1)] or null
		    final double [][][][]     clt_combo_out,  // sparse (by the first index) [type][tilesY][tilesX][(combo_tile_size] or null
			final double [][][][]     clt_combo_dbg,  // generate sparse  partial rotated/scaled pairs
			// When clt_mismatch is non-zero, no far objects extraction will be attempted
			// clt_mismatch is used in older code, not supported in GPU - there is cltMeasureLazyEye for that purpose
			// this.correlation2d should be not null if disparity_map != null
			final double [][]         disparity_map,   // [8][tilesY][tilesX], only [6][] is needed on input or null - do not calculate
			                                           // last 2 - contrast, avg/ "geometric average)
			
			final int                 width,
			final double              corr_fat_zero,    // add to denominator to modify phase correlation (same units as data1, data2). <0 - pure sum
			final boolean             corr_sym,
			final double              corr_offset,
			final double              corr_red,
			final double              corr_blue,
			final double              corr_sigma,

			final GeometryCorrection  geometryCorrection,
			final GeometryCorrection  geometryCorrection_main, // if not null correct this camera (aux) to the coordinates of the main
			final double [][][][][][] clt_kernels, // [channel_in_quad][color][tileY][tileX][band][pixel] , size should match image (have 1 tile around)
			final int                 kernel_step,
			final int                 window_type,
			final double [][]         shiftXY, // [port]{shiftX,shiftY}
			final double              disparity_corr, // disparity at infinity
			final double [][][]       fine_corr, // quadratic coefficients for fine correction (or null)

			final double              shiftX, // shift image horizontally (positive - right) - just for testing
			final double              shiftY, // shift image vertically (positive - down)
			
			final int                 clustRadius,  // 1 - single tile, 2 - 3x3, 3 - 5x5, ...
			final double              arange,  // absolute disparity range to consolidate
			final double              rrange,  // relative disparity range to consolidate
			final double              no_tilt, // no tilt if center disparity is lower
			final double              damp_tilt,    // 0.1?
			final int                 mcorr_sel,    // Which pairs to correlate // +1 - all, +2 - dia, +4 - sq, +8 - neibs, +16 - hor + 32 - vert
			final int                 mcorr_comb_width,  // combined correlation tile width
			final int                 mcorr_comb_height, // combined correlation tile full height
			final int                 mcorr_comb_offset, // combined correlation tile height offset: 0 - centered (-height/2 to height/2), height/2 - only positive (0 to height)
			final double              mcorr_comb_disp,   // Combined tile per-pixel disparity for baseline == side of a square
			
			final int                 debug_tileX,
			final int                 debug_tileY,
			final boolean             no_fract_shift,
			final boolean             no_deconvolution,
			final int                 threadsMax,  // maximal number of threads to launch
			final int                 globalDebugLevel)
	{
		final double [] damping = {damp_tilt, damp_tilt, 0.0}; // 0.0 will be applied to average value, tilt_cost - to both tilts
		final int clustDiameter = 2 * clustRadius - 1;
		final int center_indx = (clustRadius - 1) * (clustDiameter + 1);
		final double [] wnd_neib = new double [clustDiameter * clustDiameter];
		for (int iy = 0; iy < clustDiameter; iy++) {
			double wy = Math.sin(Math.PI *(iy+1) / (clustDiameter + 1));
			for (int ix = 0; ix < clustDiameter; ix++) {
				wnd_neib[iy * clustDiameter + ix] = wy * Math.sin(Math.PI *(ix+1) / (clustDiameter + 1));
			}
		}
		final int height=image_data[0][0].length/width;
		final int tilesX=width/transform_size;
		final int tilesY=height/transform_size;
		boolean [][] pcombo_sels = null;
		if (correlation2d == null) {
			  throw new IllegalArgumentException ("clt_aberrations_quad_corr_multi(): correlation2d == null!");
		}
		if (disparity_map == null) {
			  throw new IllegalArgumentException ("clt_aberrations_quad_corr_multi(): disparity_map == null!");
		}
		// Initialize correlation pairs selection to be used by all threads
		boolean [] corr_calculate = null;
		if (isCorrAll  (mcorr_sel)) corr_calculate = correlation2d.selectAll();
		if (isCorrDia  (mcorr_sel)) corr_calculate = correlation2d.selectDiameters  (corr_calculate);
		if (isCorrSq   (mcorr_sel)) corr_calculate = correlation2d.selectSquares    (corr_calculate);
		if (isCorrNeib (mcorr_sel)) corr_calculate = correlation2d.selectNeibs      (corr_calculate);
		if (isCorrHor  (mcorr_sel)) corr_calculate = correlation2d.selectHorizontal (corr_calculate);
		if (isCorrVert (mcorr_sel)) corr_calculate = correlation2d.selectVertical   (corr_calculate);
		correlation2d.setCorrPairs(corr_calculate); // will limit correlation pairs calculation
		correlation2d.generateResample( // should be called before
				mcorr_comb_width,  // combined correlation tile width
				mcorr_comb_height, // combined correlation tile full height
				mcorr_comb_offset, // combined correlation tile height offset: 0 - centered (-height/2 to height/2), height/2 - only positive (0 to height)
				mcorr_comb_disp);
		pcombo_sels = new boolean [Correlation2d.MCORR_COMB.values().length][];
		if (imgdtt_params.mcorr_cons_all) {
			int indx = Correlation2d.MCORR_COMB.ALL.ordinal();
			pcombo_sels[indx] = correlation2d.selectAll(); 
		}
		if (imgdtt_params.mcorr_cons_dia) {
			int indx = Correlation2d.MCORR_COMB.DIA.ordinal();
			pcombo_sels[indx] = correlation2d.selectDiameters(null); 
		}
		if (imgdtt_params.mcorr_cons_sq) {
			int indx = Correlation2d.MCORR_COMB.SQ.ordinal();
			pcombo_sels[indx] = correlation2d.selectSquares(null); 
		}
		if (imgdtt_params.mcorr_cons_neib) {
			int indx = Correlation2d.MCORR_COMB.NEIB.ordinal();
			pcombo_sels[indx] = correlation2d.selectNeibs(null); 
		}
		if (imgdtt_params.mcorr_cons_hor) {
			int indx = Correlation2d.MCORR_COMB.HOR.ordinal();
			pcombo_sels[indx] = correlation2d.selectHorizontal(null); 
		}
		if (imgdtt_params.mcorr_cons_vert) {
			int indx = Correlation2d.MCORR_COMB.VERT.ordinal();
			pcombo_sels[indx] = correlation2d.selectVertical(null); 
		}
		// if it was null, it will not be visible form outside, only used internally
		final double [][][][] clt_corr_out = (clt_corr_out0 != null) ? clt_corr_out0 : new double [correlation2d.getCorrPairs().length][][][]; 
		{
			boolean [] calc_corr_pairs = correlation2d.getCorrPairs();
			for (int i = 0; i < calc_corr_pairs.length; i++) if (calc_corr_pairs[i]){
				clt_corr_out[i] = new double[tilesY][tilesX][]; 
			}
		}

		if (clt_combo_out != null) {
			for (int i = 0; i < pcombo_sels.length; i++) if (pcombo_sels[i] != null){
				clt_combo_out[i] = new double[tilesY][tilesX][]; 
			}
		}
		if (clt_combo_dbg != null) {
			boolean [] calc_corr_pairs = correlation2d.getCorrPairs();
			for (int i = 0; i < calc_corr_pairs.length; i++) if (calc_corr_pairs[i]){
				clt_combo_dbg[i] = new double[tilesY][tilesX][]; 
			}
		}
//		final boolean [][] combo_sels = pcombo_sels; 
		final boolean debug_distort= globalDebugLevel > 0; ///false; // true;
		final double [][] debug_offsets = null;
		final boolean macro_mode = macro_scale != 1;      // correlate tile data instead of the pixel data
		final int numcol = 3; // number of colors // keep the same, just do not use [0] and [1], [2] - green
		final int nTilesInChn=tilesX*tilesY;
//		final double [][][][][][] clt_data = new double[numSensors][numcol][tilesY][tilesX][][];
		final Thread[] threads = newThreadArray(threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		final double [] col_weights= new double [numcol]; // colors are RBG
		final double [][] dbg_distort = debug_distort? (new double [4*numSensors][tilesX*tilesY]) : null;
		final double [][] corr_wnd = Corr2dLMA.getCorrWnd(
				transform_size,
				imgdtt_params.lma_wnd);
		final double [] corr_wnd_inv_limited = (imgdtt_params.lma_min_wnd <= 1.0)?  new double [corr_wnd.length * corr_wnd[0].length]: null;
		if (corr_wnd_inv_limited != null) {
			double inv_pwr = imgdtt_params.lma_wnd_pwr - (imgdtt_params.lma_wnd - 1.0); // compensate for lma_wnd
			for (int i = imgdtt_params.lma_hard_marg; i < (corr_wnd.length - imgdtt_params.lma_hard_marg); i++) {
				for (int j = imgdtt_params.lma_hard_marg; j < (corr_wnd.length - imgdtt_params.lma_hard_marg); j++) {
					corr_wnd_inv_limited[i * (corr_wnd.length) + j] = 1.0/Math.max(Math.pow(corr_wnd[i][j],inv_pwr), imgdtt_params.lma_min_wnd);
				}
			}
		}

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

		// TODO: move transpose_indices to Correlation2d class
		
		final int corr_size = transform_size * 2 -1;

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

		final double [] filter =     doubleGetCltLpfFd(corr_sigma);
		
		if (disparity_map != null){
			for (int i = 0; i<disparity_map.length;i++){
				if (i == OVEREXPOSED) {
					if (saturation_imp!= null) {
						disparity_map[i] = new double [tilesY*tilesX];
					}
				} else if (isSliceBit(i)) {
					disparity_map[i] = new double [tilesY*tilesX];
				} else if (isDiffIndex(i)) {
					disparity_map[i] = new double [tilesY*tilesX];
				}
			}
		}
		DttRad2 dtt = new DttRad2(transform_size);
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
					int tileYC,tileXC; // ,tIndex; // , chn;
					double centerX; // center of aberration-corrected (common model) tile, X
					double centerY; //
					double [][] fract_shiftsXY = new double[numSensors][];
					PolynomialApproximation pa = new PolynomialApproximation();
					double [][] rXY;
					if (use_main) {
						rXY = geometryCorrection.getRXY(true); // boolean use_rig_offsets,
					} else  {
						rXY = geometryCorrection.getRXY(false); // boolean use_rig_offsets,
					}
					double [][][][] clt_data_tile = new double[numSensors][numcol][][];
					for (int nTileC = ai.getAndIncrement(); nTileC < nTilesInChn; nTileC = ai.getAndIncrement()) {
						tileYC = nTileC /tilesX;
						tileXC = nTileC % tilesX;
						// now only process selected tiles
						if (tile_op[tileYC][tileXC] == 0) { // see if it is still needed
							continue;
						}

						boolean debugTile =(tileXC == debug_tileX) && (tileYC == debug_tileY) && (globalDebugLevel > -1);
						boolean debugTile0 =(tileXC == debug_tileX) && (tileYC == debug_tileY) && (globalDebugLevel > -3);
						final int [] overexp_all = (saturation_imp != null) ? ( new int [2]): null;
						if (debugTile0) {
							System.out.println("tileXC = "+tileXC+", tileYC = "+tileYC);
						}
						for (int i = 0; i < disparity_map.length; i++) {
							if (disparity_map[i] != null) disparity_map[i][nTileC] = (
									(i == DISPARITY_STRENGTH_INDEX) ||
									(i == DISPARITY_INDEX_HOR_STRENGTH) ||
									(i == DISPARITY_INDEX_VERT_STRENGTH)) ? 0.0 : Double.NaN; // once and for all
						}

						//calculate tilts
						TileNeibs tn = new TileNeibs(tilesX,tilesY);
						boolean [] used_neibs = new boolean [clustDiameter * clustDiameter];
						double disparity_center = disparity_array[tileYC][tileXC];
						double disp_min = disparity_center - arange - rrange * Math.abs(disparity_center);
						double disp_max = disparity_center + arange + rrange * Math.abs(disparity_center);
						double tiltY = 0, tiltX=0;
						double [][][] mdata = new double [wnd_neib.length][][]; // now empty lines will be skipped [3][];
						double sum_w = 0.0;
						for (int dty = -clustRadius+1; dty < clustRadius; dty++) {
							for (int dtx = -clustRadius+1; dtx < clustRadius; dtx++) {
								int nTile1 = tn.getNeibIndex(nTileC, dtx, dty);
								if (nTile1 >= 0){
									int tileX1 = nTile1 % tilesX;
									int tileY1 = nTile1 / tilesX;
									if ((tile_op[tileY1][tileX1] != 0) && (disparity_array[tileY1][tileX1] >= disp_min)  && (disparity_array[tileY1][tileX1] <= disp_max)){
										int mindx = (dty * clustDiameter + dtx + center_indx);
										used_neibs[mindx] = true;
										double w = wnd_neib[mindx];
										mdata[mindx] = new double[3][];
										mdata[mindx][0] = new double [2];
										mdata[mindx][0][0] =  dtx;
										mdata[mindx][0][1] =  dty;
										mdata[mindx][1] = new double [1];
										mdata[mindx][1][0] =  disparity_array[tileY1][tileX1]; 
										mdata[mindx][2] = new double [1];
										mdata[mindx][2][0] =  w;
										sum_w += w;
									}
								}
							}
						}
						double scale_weighths = 1.0/sum_w;
						if (disparity_center > no_tilt) {
							double[][] approx2d = pa.quadraticApproximation(
									mdata,
									true,          // boolean forceLinear,  // use linear approximation
									damping,       // double [] damping,
									-1);           // debug level
							if (approx2d != null){
								tiltX = approx2d[0][0];
								tiltY = approx2d[0][1]; // approx2d[0][2] - const C (A*x+B*y+C), C is not used here 
							}
//							if (num_tiles == 0) {
//								continue; // should never happen anyway
//							}
						}
						
						double [][] corrs_cons = new double [correlation2d.getNumPairs()][];
						double [][] disp_dist_cons = new double[numSensors][]; // used to correct 3D correlations
						int num_tiles = 0;

						for (int dty = -clustRadius+1; dty < clustRadius; dty++) {
							for (int dtx = -clustRadius+1; dtx < clustRadius; dtx++) {
								int mindx = (dty * clustDiameter + dtx + center_indx);
								if (used_neibs[mindx]) {
									int nTile = tn.getNeibIndex(nTileC, dtx, dty);
									if (nTile >= 0){
										int tileX = nTile % tilesX;
										int tileY = nTile / tilesX;
										
										double disparity_target = disparity_center + dty * tiltY + dtx * tiltX;

										// Moved from inside chn loop
										centerX = tileX * transform_size + transform_size/2 - shiftX;
										centerY = tileY * transform_size + transform_size/2 - shiftY;
										// TODO: move port coordinates out of color channel loop
										double [][] centersXY;
										double [][] disp_dist = new double[numSensors][]; // used to correct 3D correlations

										if (macro_mode){
											if ((globalDebugLevel > -1) && debugTile) { // before correction
												System.out.println("\nUsing MACRO mode, centerX="+centerX+", centerY="+centerY);
											}
											centersXY = geometryCorrection.getPortsCoordinatesIdeal(
													macro_scale,
													centerX,
													centerY,
													macro_scale* disparity_target + disparity_corr);

										} else {
											if (use_main) { // this is AUX camera that uses main coordinates  // not used in lwir
												centersXY =  geometryCorrection.getPortsCoordinatesAndDerivatives(
														geometryCorrection_main, //			GeometryCorrection gc_main,
														true,            // boolean use_rig_offsets,
														corr_rots,       // Matrix []   rots,
														null,            //  Matrix [][] deriv_rots,
														null,            // double [][] pXYderiv, // if not null, should be double[8][]
														disp_dist,       // used to correct 3D correlations
														centerX,
														centerY,
														disparity_target + disparity_corr); // _aux); //  + disparity_corr);


											} else {  // used in lwir
												centersXY = geometryCorrection.getPortsCoordinatesAndDerivatives(
														geometryCorrection, //			GeometryCorrection gc_main,
														false,           // boolean use_rig_offsets,
														corr_rots,       // Matrix []   rots,
														null,            //  Matrix [][] deriv_rots,
														null,            // double [][] pXYderiv, // if not null, should be double[8][]
														disp_dist,       // used to correct 3D correlations
														centerX,
														centerY,
														disparity_target + disparity_corr);
											}

											if (((globalDebugLevel > 0) || debug_distort) && debugTile) {
												for (int i = 0; i < numSensors; i++) {
													System.out.println("clt_aberrations_quad_corr():  tileX="+tileX+", tileY="+tileY+
															" centerX="+centerX+" centerY="+centerY+" disparity="+disparity_target+
															" centersXY["+i+"][0]="+centersXY[i][0]+" centersXY["+i+"][1]="+centersXY[i][1]);
												}
											}
											if (debug_distort && debugTile && (debug_offsets != null)) {
												double [][] debug_offsets_xy = new double [debug_offsets.length][2];
												for (int i = 0; i < debug_offsets.length; i++) {
													debug_offsets_xy[i][0] = disp_dist[i][0] * debug_offsets[i][0] + disp_dist[i][1] * debug_offsets[i][1];
													debug_offsets_xy[i][1] = disp_dist[i][2] * debug_offsets[i][0] + disp_dist[i][3] * debug_offsets[i][1];
												}
												for (int i = 0; i < numSensors; i++) {
													System.out.println(String.format("%d: {%8.3f, %8.3f}",i,debug_offsets_xy[i][0],debug_offsets_xy[i][1]));
												}

												for (int i = 0; i < debug_offsets.length; i++) {
													centersXY[i][0] += debug_offsets_xy[i][0];
													centersXY[i][1] += debug_offsets_xy[i][1];
												}
												for (int i = 0; i < numSensors; i++) {
													System.out.println("Corrected clt_aberrations_quad_corr():  tileX="+tileX+", tileY="+tileY+
															" centerX="+centerX+" centerY="+centerY+" disparity="+disparity_target+
															" centersXY["+i+"][0]="+centersXY[i][0]+" centersXY["+i+"][1]="+centersXY[i][1]);
												}
											}

											if ((globalDebugLevel > -1) && (tileX == debug_tileX) && (tileY == debug_tileY)) { // before correction
												System.out.print(disparity_target+"\t"+
														centersXY[0][0]+"\t"+centersXY[0][1]+"\t"+
														centersXY[1][0]+"\t"+centersXY[1][1]+"\t"+
														centersXY[2][0]+"\t"+centersXY[2][1]+"\t"+
														centersXY[3][0]+"\t"+centersXY[3][1]+"\t");
											}

											for (int ip = 0; ip < centersXY.length; ip++){ // 4 OOB
												centersXY[ip][0] -= shiftXY[ip][0];
												centersXY[ip][1] -= shiftXY[ip][1];
											}
											// save disparity distortions for visualization:
											if (dbg_distort != null) {
												for (int cam = 0; cam <numSensors; cam++) {
													dbg_distort[cam * 4 + 0 ][nTile] = disp_dist[cam][0];
													dbg_distort[cam * 4 + 1 ][nTile] = disp_dist[cam][1];
													dbg_distort[cam * 4 + 2 ][nTile] = disp_dist[cam][2];
													dbg_distort[cam * 4 + 3 ][nTile] = disp_dist[cam][3];
												}
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
										if (FPGA_COMPARE_DATA && (globalDebugLevel > 0) && (tileX == debug_tileX) && (tileY == debug_tileY)) { // not used in lwir
											final int fpga_cam = 0;
											double [][] manual_offsets={

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
													width,                 // final int                 width
													dtt
													);
										}
										// See if macro_mode uses color channels for non-color?
										for (int ncol = 0; ncol <numcol; ncol++) {
											if (!isMonochrome() || (ncol == MONO_CHN) || macro_mode) { // in monochrome mode skip all non-mono (green) channels  // used in lwir (5 of 6 branches)
												boolean debug_for_fpga = FPGA_COMPARE_DATA && (globalDebugLevel > 0) && (tileX == debug_tileX) && (tileY == debug_tileY) && (ncol == 2);
												if ((globalDebugLevel > 0) && (tileX == debug_tileX) && (tileY == debug_tileY) && (ncol == 2)) {
													System.out.println("\nUsing "+(macro_mode?"MACRO":"PIXEL")+" mode, centerX="+centerX+", centerY="+centerY);
													System.out.println(disparity_target+"\t"+
															centersXY[0][0]+"\t"+centersXY[0][1]+"\t"+
															centersXY[1][0]+"\t"+centersXY[1][1]+"\t"+
															centersXY[2][0]+"\t"+centersXY[2][1]+"\t"+
															centersXY[3][0]+"\t"+centersXY[3][1]+"\t");
												}

												for (int i = 0; i < numSensors; i++) {
													if (debug_for_fpga && (i==0)){ // not used in lwir
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

														fpga_chn =           2;
														System.out.println(String.format("Manually changing offset: center X= %f, center Y = %f", fpga_centersXY[0],fpga_centersXY[1]));
														System.out.println(String.format("Manually changing color to %d (was %d)", fpga_chn, ncol));
														fpga_fract_shiftsXY = extract_correct_tile( // return a pair of residual offsets
																image_data[i],
																width,       // image width
																null,
																fpga_clt_data, //double  [][]        clt_tile,    // should be double [4][];
																kernel_step,
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
													clt_data_tile[i][ncol] = new double [4][];
													// Extract image tiles and kernels, correct aberrations, return (ut do not apply) fractional shifts
													fract_shiftsXY[i] = extract_correct_tile( // return a pair of residual offsets
															image_data[i],
															width,       // image width
															((clt_kernels == null) ? null : clt_kernels[i]), // [color][tileY][tileX][band][pixel]
															clt_data_tile[i][ncol], //  clt_data[i][ncol][tileY][tileX], //double  [][]        clt_tile,    // should be double [4][];
															kernel_step,
															dtt,
															ncol,
															centersXY[i][0], // centerX, // center of aberration-corrected (common model) tile, X
															centersXY[i][1], // centerY, //
															((!FPGA_COMPARE_DATA && (globalDebugLevel > -1) && (tileX == debug_tileX) && (tileY == debug_tileY) && (ncol == 2) && (i==0)) ? (globalDebugLevel + 0) : 0), // external tile compare
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
//													for (int i = 0; i < 16; i++) dbg_tile[i]=clt_data[i>>2][ncol][tileY][tileX][i & 3];
													for (int i = 0; i < 16; i++) dbg_tile[i]=clt_data_tile[i>>2][ncol][i & 3];
													sdfa_instance.showArrays(dbg_tile,  transform_size, transform_size, true, "pre-shifted_x"+tileX+"_y"+tileY, titles);
												}

												if ((globalDebugLevel > 0) && (tileX >= debug_tileX - 2) && (tileX <= debug_tileX + 2) &&
														(tileY >= debug_tileY - 2) && (tileY <= debug_tileY+2)) {
													for (int i = 0; i < numSensors; i++) {
														System.out.println("clt_aberrations_quad(): color="+ncol+", tileX="+tileX+", tileY="+tileY+
																" fract_shiftsXY["+i+"][0]="+fract_shiftsXY[i][0]+" fract_shiftsXY["+i+"][1]="+fract_shiftsXY[i][1]);
													}
												}

												if (!no_fract_shift) {  // USED in lwir
													// apply residual shift
													for (int i = 0; i < numSensors; i++) {
														fract_shift(    // fractional shift in transform domain. Currently uses sin/cos - change to tables with 2? rotations
																clt_data_tile[i][ncol], // clt_data[i][ncol][tileY][tileX], // double  [][]  clt_tile,
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
														for (int i = 0; i < 16; i++) {
															dbg_tile[i]=clt_data_tile[i>>2][ncol][i & 3];
														}
														sdfa_instance.showArrays(dbg_tile,  transform_size, transform_size, true, "shifted_x"+tileX+"_y"+tileY+"-z", titles);
													}
												}
											} else { // if (!isMonochrome() || (chn == MONO_CHN) || macro_mode) { // in monochrome mode skip all non-mono (green) channels
												for (int i = 0; i < numSensors; i++) {  // used in lwir
													clt_data_tile[i][ncol] = null; // clt_data[i][ncol] = null; // erase unused clt_data
												}
											}
										}// end of for (int chn = 0; chn <numcol; chn++)





										// all color channels are done here
										// fill clt_corr_combo if it exists
										//clt_mismatch should only be used with disparity_map != null;
										// calculate overexposed fraction
										
										
										if (saturation_imp != null){ // TODO: Remove??
											disparity_map[OVEREXPOSED][nTile] = (1.0 * overexp_all[0]) / overexp_all[1];
										}
										// calculate all selected pairs correlations
										double [][]  corr_tiles = correlation2d.correlateCompositeFD(
												clt_data_tile,                // double [][][][][][] clt_data,
//												tileX,                        // int                 tileX,
//												tileY,                        // int                 tileY,
												correlation2d.getCorrPairs(), // boolean[]           pairs_mask,
												filter,                       // double []           lpf,
												getScaleStrengths(),          // double              scale_value, // scale correlation value
												col_weights,                  // double []           col_weights,
												corr_fat_zero);               // double              fat_zero)
										for (int num_pair = 0; num_pair < corr_tiles.length; num_pair++) if (corr_tiles[num_pair] != null){
											clt_corr_out[num_pair][tileY][tileX] = corr_tiles[num_pair];
										}
										// Init if it is the first non-null tile
										if (num_tiles == 0) {
											for (int np = 0; np < corrs_cons.length; np++) if (corr_tiles[np] != null){
												corrs_cons[np] = new double [corr_tiles[np].length];
											}
											for (int nsens = 0; nsens < numSensors; nsens++) {
												for (int i = 0; i < disp_dist[nsens].length; i++) {
													disp_dist_cons[nsens] = new double [disp_dist[nsens].length];
												}
											}
										}
										
										double w = scale_weighths * wnd_neib[mindx];
									    for (int np = 0; np < corrs_cons.length; np++) if (corr_tiles[np] != null) {
									        for (int i = 0; i < corrs_cons[np].length; i++) {
									            corrs_cons[np][i] += w * corr_tiles[np][i];
									        }
									    }
									    for (int nsens = 0; nsens < numSensors; nsens++) {
									        for (int i = 0; i < disp_dist_cons[nsens].length; i++) {
									            disp_dist_cons[nsens][i] += w * disp_dist[nsens][i];
									        }
									    }
										num_tiles++;
									}
								}
							}
						} // for (int dty = -clustRadius+1; dty < clustRadius; dty++)  - done with preparing cluster data
						int tile_lma_debug_level =  ((tileXC == debug_tileX) && (tileYC == debug_tileY))? (imgdtt_params.lma_debug_level-1) : -2;
						centerX = tileXC * transform_size + transform_size/2 - shiftX;
						centerY = tileYC * transform_size + transform_size/2 - shiftY;
						
					    double [] corr_combo_all = correlation2d.accumulateInit();
					    double sumw = correlation2d.accummulatePairs(
					            corr_combo_all,            // double []   accum_tile,
					            corrs_cons,                // double [][] corr_tiles,
					            correlation2d.selectAll(), // boolean []  selection,
					            1.0);             // double      weight);
					    correlation2d.normalizeAccumulatedPairs(
					            corr_combo_all,
					            sumw);
					    // save for output
						if (clt_corr_out != null) {
							for (int num_pair = 0; num_pair < corrs_cons.length; num_pair++) if (corrs_cons[num_pair] != null){
								clt_corr_out[num_pair][tileYC][tileXC] = corrs_cons[num_pair];
							}
						}
						if (clt_combo_out != null) { // only all
								clt_combo_out[0][tileYC][tileXC] = corr_combo_all;
						}

						// calculate CM maximums for all mixed channels
						// First get integer correlation center, relative to the center
						if (debugTile0) {
							System.out.println("tileXC = "+tileXC+", tileYC = "+tileYC);
						}
					    double [] disp_str_combo = new double[2];
					    int [] ixy =   correlation2d.getMaxXYInt( // find integer pair or null if below threshold // USED in lwir
					            corr_combo_all, // double [] data,      // [data_size * data_size]
					            disp_str_combo,
					            correlation2d.getCombWidth(), //       data_width,
					            correlation2d.getCombHeight()/2 - correlation2d.getCombOffset(), // int       center_row, ??????????????
					            true, // boolean   axis_only,
					            imgdtt_params.min_corr,  // ???? double    minMax,    // minimal value to consider (at integer location, not interpolated)
					            false); // debugCluster); // tile_lma_debug_level > 0); // boolean   debug);

						// if integer argmax was strong enough, calculate CM argmax
						// will not fill out DISPARITY_INDEX_INT+1, DISPARITY_INDEX_CM+1, DISPARITY_INDEX_POLY+1
						// use clt_mismatch for that
//						double [] disp_str = new double[2];
					    if (ixy != null) {
							disparity_map[DISPARITY_INDEX_INT][nTileC] =      disp_str_combo[0]; // -ixy[0];
							disparity_map[DISPARITY_INDEX_CM+1][nTileC] =     disp_str_combo[1]; // reusing old fine-corr
							disparity_map[DISPARITY_STRENGTH_INDEX][nTileC] = disp_str_combo[1];
					    	
					        double [] corr_stat = correlation2d.getMaxXCm(         // get fractional center as a "center of mass" inside circle/square from the integer max
					                corr_combo_all,                      // double [] data,      // [data_size * data_size]
					                correlation2d.getCombWidth(),        // int       data_width,      //  = 2 * transform_size - 1;
					                correlation2d.getCombHeight()/2 - correlation2d.getCombOffset(),// int       center_row,
					                ixy[0],                              // int       ixcenter,  // integer center x
					                // corr_wndy,                        // double [] window_y,  // (half) window function in y-direction(perpendicular to disparity: for row0  ==1
					                // corr_wndx,                        // double [] window_x,  // half of a window function in x (disparity) direction
					                false); // debugCluster); // (tile_lma_debug_level > 0)); // boolean   debug);
					        if (corr_stat != null) {
					            disp_str_combo[0] = -corr_stat[0]; // -ixy[0]; for CM use magic???
					            disp_str_combo[1] =  corr_stat[1];
								disparity_map[DISPARITY_INDEX_CM][nTileC] =      -corr_stat[0]; // -ixy[0];
								disparity_map[DISPARITY_INDEX_CM+1][nTileC] =     corr_stat[1]; // reusing old fine-corr
								disparity_map[DISPARITY_STRENGTH_INDEX][nTileC] = corr_stat[1];
					        }
					        if (debugTile0) { // && (globalDebugLevel > -1)) { // -2)) {
					            System.out.println("Will run new LMA for tileXC="+tileXC+", tileYC="+tileXC);
					            if (disp_str_combo != null) {
					                System.out.println("disp_str_combo[0]="+disp_str_combo[0]+", disp_str_combo[1]="+disp_str_combo[1]);
					            }
//								(new ShowDoubleFloatArrays()).showArrays(corrs_cons,  15, 15, true, "corrs_cons_CX"+clustX+"-CY"+clustY,correlation2d.getCorrTitles());
					        }
					        
					        Corr2dLMA lma2 = correlation2d.corrLMA2Single(
					                imgdtt_params,                  // ImageDttParameters  imgdtt_params,
					                false,                          // imgdtt_params.lmas_LY_single, //                   // boolean             adjust_ly, // adjust Lazy Eye
					                corr_wnd,                       // double [][]         corr_wnd, // correlation window to save on re-calculation of the window
					                corr_wnd_inv_limited,           // corr_wnd_limited, // correlation window, limited not to be smaller than threshold - used for finding max/convex areas (or null)
					                corrs_cons,                     // corrs,          // double [][]         corrs,
					                disp_dist_cons,
					                rXY,                          // double [][]         rXY, // non-distorted X,Y offset per nominal pixel of disparity
					                // all that are not null in corr_tiles
					                correlation2d.selectAll(),   // longToArray(imgdtt_params.dbg_pair_mask),  // int                 pair_mask, // which pairs to process
					                // TODO: Verify sign is correct (disparity, not X0)?
					                disp_str_combo, // null, // disp_str,  //corr_stat[0],                 // double    xcenter,   // preliminary center x in pixels for largest baseline
					                disp_str_combo, // disp_str_cons,                // double[]            poly_ds,    // null or pair of disparity/strength
					                imgdtt_params.ortho_vasw_pwr, // double    vasw_pwr,  // value as weight to this power,
					                tile_lma_debug_level, // +2,         // int                 debug_level,
					                tileXC,                        // int                 tileX, // just for debug output
					                tileYC );                      // int                 tileY
							if (debugTile0) { // should be debugTile
								System.out.println("Ran LMA for tileXC="+tileXC+", tileYC="+tileYC);
							}
							double [][] ds = null;
							if (lma2 != null) {
								ds = lma2.lmaDisparityStrength(
										imgdtt_params.lmas_min_amp,      //  minimal ratio of minimal pair correlation amplitude to maximal pair correlation amplitude
										imgdtt_params.lmas_max_rel_rms,  // maximal relative (to average max/min amplitude LMA RMS) // May be up to 0.3)
										imgdtt_params.lmas_min_strength, // minimal composite strength (sqrt(average amp squared over absolute RMS)
										imgdtt_params.lmas_min_ac,       // minimal of A and C coefficients maximum (measures sharpest point/line)
										imgdtt_params.lmas_min_min_ac,   // minimal of A and C coefficients minimum (measures sharpest point)
										imgdtt_params.lmas_max_area,      //double  lma_max_area,     // maximal half-area (if > 0.0)
										imgdtt_params.lma_str_scale,    // convert lma-generated strength to match previous ones - scale
										imgdtt_params.lma_str_offset    // convert lma-generated strength to match previous ones - add to result
										);
								if (ds != null) { // always true
									disparity_map[DISPARITY_INDEX_POLY][nTileC] =   ds[0][0];
									disparity_map[DISPARITY_INDEX_POLY+1][nTileC] = ds[0][1];
									disparity_map[DISPARITY_STRENGTH_INDEX][nTileC] = ds[0][1]; // overwrite with LMA strength
									if (debugTile0) {
										lma2.printStats(ds,1);
										if (imgdtt_params.lmas_LY_single) {
											double [][] ddnd = lma2.getDdNd();
											if (ddnd != null) {
												double [][] dxy= new double [ddnd.length][2];
												for (int i = 0; i < dxy.length; i++) {
													dxy[i][0] = ddnd[i][0] * rXY[i][0] - ddnd[i][1] * rXY[i][1];
													dxy[i][1] = ddnd[i][0] * rXY[i][1] + ddnd[i][1] * rXY[i][0];
												}
												System.out.print("       Port:  ");
												for (int i = 0; i < dxy.length; i++) System.out.print(String.format("   %2d   ", i)); System.out.println();
												System.out.print("Radial_in =  [");
												for (int i = 0; i < dxy.length; i++) System.out.print(String.format(" %6.3f,", ddnd[i][0])); System.out.println("]");
												System.out.print("Tangent_CW = [");
												for (int i = 0; i < dxy.length; i++) System.out.print(String.format(" %6.3f,", ddnd[i][1])); System.out.println("]");
												System.out.print("X =          [");
												for (int i = 0; i < dxy.length; i++) System.out.print(String.format(" %6.3f,", dxy[i][0])); System.out.println("]");
												System.out.print("Y =          [");
												for (int i = 0; i < dxy.length; i++) System.out.print(String.format(" %6.3f,", dxy[i][1])); System.out.println("]");
												System.out.println();
											} else {
												System.out.println("No dd/nd and x/y offsets data is available ");
											}
										} else {
											System.out.println("LY offsets are not measured");
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
		if ((dbg_distort != null) &&(globalDebugLevel >=0)) {
			(new ShowDoubleFloatArrays()).showArrays(dbg_distort,  tilesX, tilesY, true, "disparity_distortions"); // , dbg_titles);
		}
//		return clt_data;
	}
	
	
	public boolean dmExists(double [][] dm, int indx) {  // not used in lwir or else
		return (dm != null) && (dm.length > indx) && (dm[indx]!= null);
	}

	public double [][] tile_combine_rgba(  // used in lwir
			double [][][] iclt_tile,     // [port][numcol][256] // in mono some are null
			double []     ports_rgb,     // average values of R,G,B for each camera (R0,R1,...,B2,B3) // null
			double []     max_diff,      // maximal (weighted) deviation of each channel from the average /null
			double []     lt_window,     // [256]
			double [][]   port_offsets,  // [port]{x_off, y_off} - just to scale pixel value differences
			int           port_mask,      // which port to use, 0xf - all 4 (will modify as local variable)
			double        diff_sigma,     // pixel value/pixel change
			double        diff_threshold, // pixel value/pixel change
			// next not used
			boolean       diff_gauss,     // when averaging images, use gaussian around average as weight (false - sharp all/nothing)
			double        min_agree,      // minimal number of channels to agree on a point (real number to work with fuzzy averages)
			double []     chn_weights,     // color channel weights, sum == 1.0
			boolean       dust_remove,    // Do not reduce average weight when only one image differes much from the average
			boolean       keep_weights,   // return channel weights after A in RGBA
			boolean       debug,
			boolean       debug_gpu)      // generate output fro matching with GPU processing
	{
		if (debug_gpu) {
			System.out.println("tile_combine_rgba(): generatinhg GPU debug data");
			System.out.println(String.format(
					"diff_sigma =     %f\n" + // 1.5
					"diff_threshold = %f\n" + // 10.0
					"diff_gauss =     %b\n" + // false never used?
					"min_agree =      %f\n" + // 3.0
					"dust_remove =    %b\n" + // true
					"keep_weights =   %b\n",  // true
					diff_sigma, diff_threshold, diff_gauss, min_agree, dust_remove, keep_weights));
			for (int i = 0; i < chn_weights.length; i++) {
				System.out.println(String.format("Channel%d weight = %f", i, chn_weights[i]));
			}
			// Channel0 weight = 0.294118
			// Channel1 weight = 0.117647
			// Channel2 weight = 0.588235

		}

		int ports =   iclt_tile.length;    // number of cameras
		int numcol =  iclt_tile[0].length; // number of color channels
		int tile_len = 0;
		for (int ncol = 0; ncol < numcol; ncol++) {
			if (iclt_tile[0][ncol] != null) {
				tile_len = iclt_tile[0][ncol].length;
				break;
			}
		}

//		int usedPorts = ((port_mask >> 0) & 1) + ((port_mask >> 1) & 1) + ((port_mask >> 2) & 1) + ((port_mask >> 3) & 1);
		int usedPorts = 0;
		for (int i = port_mask; i != 0; i >>= 1) {
			if (( i & 1) !=0) {
				usedPorts++;
			}
		}

		double [][] port_weights = new double[ports][tile_len];
		double [][] color_avg =    new double[numcol][tile_len];
//		double [][] rgba = new double[numcol + 1 + (keep_weights?(ports + numcol + 1):0)][];
// need to pass keep_weights to the caller
		double [][] rgba = new double[numcol + 1 + ports + (keep_weights?(numcol + 1):0)][];
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
			if (debug_gpu) {
				for (int ncol = 0; ncol < numcol; ncol ++) {
					System.out.println("--- rms_col["+ncol+"] ---");
					for (int ii = 0; ii < 2 * transform_size; ii++) {
						for (int jj = 0; jj < 2 * transform_size; jj++) {
							System.out.print(String.format("%10.4f ", rgba[rms_start+ncol][2* transform_size * ii + jj]));
						}
						System.out.println();
					}
				}
				System.out.println("--- rms_col combo ---");
				for (int ii = 0; ii < 2 * transform_size; ii++) {
					for (int jj = 0; jj < 2 * transform_size; jj++) {
						System.out.print(String.format("%10.4f ", rgba[rms_start+numcol][2* transform_size * ii + jj]));
					}
					System.out.println();
				}
			}
		}


		double []  alpha = new double[tile_len];
		double threshold2 = diff_sigma * diff_threshold;
		threshold2 *= threshold2; // squared to compare with diff^2
		if (usedPorts > 1) { // lwir16:
			double [] pair_dist2r =  new double [ports*(ports-1)/2]; // reversed squared distance between images - to be used with gaussian. Can be calculated once !
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
					d = Math.exp(-pair_dist2r[ip]*d)+GPUTileProcessor.FAT_ZERO_WEIGHT; // 0.5 for exact match, lower for mismatch. Add this weight to both ports involved
					// Add weight to both channels in a pair
					port_weights[pair_ports[ip][0]][i] +=d;
					port_weights[pair_ports[ip][1]][i] +=d;
				}
// TODO: Merge loops when debug is done
			}
			if (debug_gpu) {
				for (int ccam = 0; ccam < port_weights.length; ccam ++) {
					System.out.println("--- port_weights["+ccam+"] ---");
					for (int ii = 0; ii < 2 * transform_size; ii++) {
						for (int jj = 0; jj < 2 * transform_size; jj++) {
							System.out.print(String.format("%10.4f ", port_weights[ccam][2* transform_size * ii + jj]));
						}
						System.out.println();
					}
				}
			}
			int [][] dbg_bestPorts = new int [2][tile_len];
			for (int i = 0; i < tile_len; i++){
				// find 2 best ports (resolving 2 pairs of close values)
				int bestPort1=0;
				for (int ip = bestPort1+1; ip < ports; ip++) if (port_weights[ip][i] > port_weights[bestPort1][i]) bestPort1 = ip;
				int bestPort2 = (bestPort1 == 0)?1:0;
				for (int ip = bestPort2+1; ip < ports; ip++) if ((ip != bestPort1) && (port_weights[ip][i] > port_weights[bestPort2][i])) bestPort2 = ip;
				dbg_bestPorts[0][i] = bestPort1;
				dbg_bestPorts[1][i] = bestPort2;
				// find weighted average between these 2 ports
				double w1 = port_weights[bestPort1][i]/(port_weights[bestPort1][i]+port_weights[bestPort2][i]);
				double w2 = 1.0 - w1;
				for (int ncol = 0; ncol < numcol; ncol++) if (iclt_tile[0][ncol] != null) {
					color_avg[ncol][i] = w1 * iclt_tile[bestPort1][ncol][i] + w2 * iclt_tile[bestPort2][ncol][i];
				}
			}
			if (debug_gpu) {
				System.out.println("--- bestPorts ---");
				for (int ii = 0; ii < 2 * transform_size; ii++) {
					for (int jj = 0; jj < 2 * transform_size; jj++) {
						System.out.print(String.format("%1d[%1d] ", dbg_bestPorts[0][2* transform_size * ii + jj],dbg_bestPorts[1][2* transform_size * ii + jj]));
					}
					System.out.println();
				}
				for (int ncol = 0; ncol < numcol; ncol++) if (iclt_tile[0][ncol] != null) {
					System.out.println("--- color_avg["+ncol+" ---");
					for (int ii = 0; ii < 2 * transform_size; ii++) {
						for (int jj = 0; jj < 2 * transform_size; jj++) {
							System.out.print(String.format("%10.4f ", color_avg[ncol][2* transform_size * ii + jj]));
						}
						System.out.println();
					}
				}
			}

			for (int i = 0; i < tile_len; i++){
				// recalculate all weights using difference from this average of the best pair
				double [] d2 = new double [ports]; //weighted squared differences
				for (int ip = 0; ip < ports; ip++) if ((port_mask & ( 1 << ip)) != 0){
					d2[ip] = 0;
					for (int ncol = 0; ncol < numcol; ncol++)  if (iclt_tile[0][ncol] != null){
						double dc = iclt_tile[ip][ncol][i] - color_avg[ncol][i];
						dc /= lt_window[i]; // to compensate fading near the edges
						d2[ip]+= chn_weights[ncol]*dc*dc;
					}
					port_weights[ip][i] = Math.exp(-ksigma * d2[ip]) + GPUTileProcessor.FAT_ZERO_WEIGHT;
				}
			}
			// dust remove here - otherwise it seems to do nothing
			if (dust_remove && (usedPorts == 4)) {
				dust_remove(port_weights);
			}
			for (int i = 0; i < tile_len; i++){
				// and now make a new average with those weights
				double k = 0.0;
				for (int ip = 0; ip < ports; ip++) k+=port_weights[ip][i];
				k = 1.0/k;
				for (int ncol = 0; ncol < numcol; ncol++) if (iclt_tile[0][ncol] != null) {
					color_avg[ncol][i] = 0;
					for (int ip = 0; ip < ports; ip++) {
						color_avg[ncol][i] += k* port_weights[ip][i] * iclt_tile[ip][ncol][i]; // sum of k* port_weights[ip][i] == 1.0
					}
				}

			} // or (int i = 0; i < tile_len; i++){

			if (debug_gpu) {
				for (int ccam = 0; ccam < port_weights.length; ccam ++) {
					System.out.println("--- UPDATED port_weights["+ccam+"] ---");
					for (int ii = 0; ii < 2 * transform_size; ii++) {
						for (int jj = 0; jj < 2 * transform_size; jj++) {
							System.out.print(String.format("%10.4f ", port_weights[ccam][2* transform_size * ii + jj]));
						}
						System.out.println();
					}
				}
				for (int ncol = 0; ncol < numcol; ncol++) if (iclt_tile[0][ncol] != null) {
					System.out.println("--- UPDATED color_avg["+ncol+" ---");
					for (int ii = 0; ii < 2 * transform_size; ii++) {
						for (int jj = 0; jj < 2 * transform_size; jj++) {
							System.out.print(String.format("%10.4f ", color_avg[ncol][2* transform_size * ii + jj]));
						}
						System.out.println();
					}
				}
			}


		} else if (usedPorts > 0){ // just copy from a single channel  // not used in lwir
			for (int ip = 0; ip < ports; ip++) if ((port_mask & ( 1 << ip)) != 0){
				for (int i = 0; i < tile_len; i++){
					for (int ncol = 0; ncol < numcol; ncol++)  if (iclt_tile[0][ncol] != null){
						color_avg[ncol][i] = iclt_tile[ip][ncol][i];
					}
					port_weights[ip][i] = 1.0; // lt_window[i]; // or use 1.0?
				}
			}
		}
		// Original location of dust_remove = was it too late?
///		if (dust_remove && (usedPorts == 4)) {
///			dust_remove(port_weights);
///		}
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
//		if (keep_weights){
			for (int i = 0; i < ports; i++) {
				rgba[numcol + 1 + i] = port_weights[i];
			}
//		}
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
		if (debug_gpu) {
			System.out.println("--- color_avg  ---");
			for (int ncol = 0; ncol < numcol; ncol++) if (iclt_tile[0][ncol] != null) {
				System.out.println("\n --- color_avg ["+ncol+"] ---");
				for (int ii = 0; ii < 2 * transform_size; ii++) {
					for (int jj = 0; jj < 2 * transform_size; jj++) {
						System.out.print(String.format("%10.4f ", color_avg[ncol][2* transform_size * ii + jj]));
					}
					System.out.println();
				}
				for (int ccam = 0; ccam < port_weights.length; ccam ++) {
					System.out.println("--- max_diff ["+ccam+"] ---");
					System.out.println("\n --- imclt ["+ccam+"]["+ncol+"] ---");
					for (int ii = 0; ii < 2 * transform_size; ii++) {
						for (int jj = 0; jj < 2 * transform_size; jj++) {
							System.out.print(String.format("%10.4f ", iclt_tile[ccam][ncol][2* transform_size * ii + jj]));
						}
						System.out.println();
					}
				}
			}
			System.out.println("--- max_diff: "+max_diff[0]+","+max_diff[1]+","+max_diff[2]+","+max_diff[3]);
			System.out.println("---        R: "+ports_rgb[ 0]+", "+ports_rgb[ 1]+", "+ports_rgb[ 2]+", "+ports_rgb[ 3]);
			System.out.println("---        B: "+ports_rgb[ 4]+", "+ports_rgb[ 5]+", "+ports_rgb[ 6]+", "+ports_rgb[ 7]);
			System.out.println("---        G: "+ports_rgb[ 8]+", "+ports_rgb[ 9]+", "+ports_rgb[10]+", "+ports_rgb[11]);
		}
		return rgba;
	}

	public void dust_remove( // redistribute weight between 3 best ports (use only when all 3 are enabled)  // USED in lwir
			double [][]  port_weights)
	{
		int np = port_weights.length;
		for (int i = 0; i < port_weights[0].length; i++){

			int wi = 0;
			for (int ip = 1; ip < np; ip++) if (port_weights[ip][i] < port_weights[wi][i]) wi = ip;
			double avg = 0;
			for (int ip = 0; ip < np; ip++) if (ip != wi) avg += port_weights[ip][i];
			avg /= (np -1);
			double worst = port_weights[wi][i];
//			double scale = 1.0 + (avg - port_weights[wi][i])/(avg * (np -1));
			double scale = 1.0 + worst * (avg - worst)/(avg * avg * (np -1));
			for (int ip = 0; ip < np; ip++) {
				if (ip != wi) port_weights[ip][i] *= scale; // increase weight of non-worst, so if worst == 0.0 sum of 3 (all) ports will be scaled by 4/3, keeping average
			}
			port_weights[wi][i] *= worst / avg;
		}
	}


	public double [] tile_debayer_shot_corr(  // USED in lwir
			boolean   rb,
			double [] tile,
			int tile_size,
			double [] window2, // squared lapping window
			double min_shot,   // 10.0;  // Do not adjust for shot noise if lower than
			double scale_shot,  //3.0;   // scale when dividing by sqrt
			double [] window_back, // re-apply window to the result
			boolean  debug_gpu)

	{
		double [] tile_nw = new double [tile.length];
		for (int i = 0; i < tile.length; i++) tile_nw[i] = tile[i]/window2[i]; //unapply squared window
		if (debug_gpu) {
			System.out.println("=== tile_debayer_shot_corr after window2 ===");
			for (int i = 0; i < 2 * transform_size; i++) {
				for (int j = 0; j < 2 * transform_size; j++) {
					System.out.print(String.format("%10.4f ", tile_nw[2* transform_size * i + j]));
				}
				System.out.println();
			}
		}
		double [] tile_db = tile_debayer(
				rb,
				tile_nw,
				tile_size,
				debug_gpu);
		if (debug_gpu) {
			System.out.println("=== tile_debayer_shot_corr rb="+rb+", min_shot="+min_shot+", scale_shot="+scale_shot+" ===");
			for (int i = 0; i < 2 * transform_size; i++) {
				for (int j = 0; j < 2 * transform_size; j++) {
					System.out.print(String.format("%10.4f ", tile_db[2* transform_size * i + j]));
				}
				System.out.println();
			}
		}

		if (!isLwir() && (scale_shot > 0)){
			double k = 1.0/Math.sqrt(min_shot);
			for (int i = 0; i < tile.length; i++) tile_db[i] = scale_shot* ((tile_db[i] > min_shot)? Math.sqrt(tile_db[i]) : (k*tile_db[i]));
			if (debug_gpu) {
				System.out.println("=== tile_debayer_shot_corr after scale_sho ===");
				for (int i = 0; i < 2 * transform_size; i++) {
					for (int j = 0; j < 2 * transform_size; j++) {
						System.out.print(String.format("%10.4f ", tile_db[2* transform_size * i + j]));
					}
					System.out.println();
				}
			}
		}

		if (window_back != null) {
			for (int i = 0; i < tile.length; i++) tile_db[i] = tile_db[i] * window_back[i]; // optionally re-apply window (may be a different one)
			if (debug_gpu) {
				System.out.println("=== tile_debayer_shot_corr after window_back ===");
				for (int i = 0; i < 2 * transform_size; i++) {
					for (int j = 0; j < 2 * transform_size; j++) {
						System.out.print(String.format("%10.4f ", tile_db[2* transform_size * i + j]));
					}
					System.out.println();
				}
			}
		}
		return tile_db;
	}


	public double [] tile_debayer( // USED in lwir
			boolean   rb,
			double [] tile,
			int tile_size,
			boolean debug)
	{
		int [] neib_indices = {-tile_size - 1, -tile_size, -tile_size + 1, -1, 0, 1, tile_size - 1, tile_size, tile_size + 1};
		int tsm1 = tile_size - 1;
		double [] rslt = new double [tile_size*tile_size]; // assuming cleared to 0.0;
		double [] kern = rb ? kern_rb : kern_g;
		double k_corn = rb? (16.0/9.0):(4.0/3.0); // reciprocal of sum kern_* inside respective corner
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
	public double [] setMaxXYWeights(  // not used in lwir
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

	public int [] getMaxXYInt( // find integer pair or null if below threshold  // not used in lwir
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

	public double [] getMaxXYCm( // get fractiona center as a "center of mass" inside circle/square from the integer max  // not used in lwir
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

	public double [] getMaxXSOrtho( // // get fractional center using a quadratic polynomial  // not used in lwir
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

	public double [] getMaxXSOrtho2(   // get fractional center using a quadratic polynomial  // not used in lwir
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



	public double [] getMaxXYPoly( // get interpolated maximum coordinates using 2-nd degree polynomial  // not used in lwir
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
	public double [][][][] clt_2d(  // not used in lwir
			final double [] dpixels,
			final int       width,
//			final int       transform_size,
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
		final int tilesX=width/transform_size-1;
		final int tilesY=height/transform_size-1;
		final int nTiles=tilesX*tilesY;
		final double [][][][] dct_data = new double[tilesY][tilesX][4][transform_size*transform_size];
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
		double [] dc = new double [transform_size*transform_size];
		for (int i = 0; i<dc.length; i++) dc[i] = 1.0;
		DttRad2 dtt0 = new DttRad2(transform_size);
		dtt0.set_window(window_type);
		if (globalDebugLevel > 0) {
			System.out.println("clt_2d(): width="+width+" height="+height+" transform_size="+transform_size+
					" debug_tileX="+debug_tileX+" debug_tileY="+debug_tileY+" globalDebugLevel="+globalDebugLevel);
		}
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				@Override
				public void run() {
					DttRad2 dtt = new DttRad2(transform_size);
					dtt.set_window(window_type);
					double [] tile_in = new double[4*transform_size * transform_size];
					double [][] tile_folded = new double[4][];
					double [][] tile_out =    new double[4][]; // = new double[dct_size * dct_size];
					int tileY,tileX;
					int n2 = transform_size * 2;
//					showDoubleFloatArrays sdfa_instance = new showDoubleFloatArrays(); // just for debugging?
					for (int nTile = ai.getAndIncrement(); nTile < nTiles; nTile = ai.getAndIncrement()) {
						tileY = nTile/tilesX;
						tileX = nTile - tileY * tilesX;
						if ((shiftX == 0) && (shiftY == 0)){
							for (int i = 0; i < n2;i++){
								System.arraycopy(dpixels, (tileY*width+tileX)*transform_size + i*width, tile_in, i*n2, n2);
							}
						} else {
							int x0 = tileX * transform_size - shiftX;
							if      (x0 < 0)             x0 = 0; // first/last will be incorrect
							else if (x0 >= (width - n2)) x0 = width - n2;
							for (int i = 0; i < n2;i++){
								int y0 = tileY * transform_size + i - shiftY;
								if      (y0 < 0)       y0 = 0;
								else if (y0 >= height) y0 = height -1;
								System.arraycopy(dpixels, y0 * width+ x0, tile_in, i*n2, n2);
							}
						}
						for (int dct_mode = 0; dct_mode <4; dct_mode++) {
							tile_folded[dct_mode] = dtt.fold_tile(tile_in, transform_size, dct_mode); // DCCT, DSCT, DCST, DSST
							if ((debug_mode & 1) != 0) {
								tile_out[dct_mode] = tile_folded[dct_mode];
							} else {
								tile_out[dct_mode] =    dtt.dttt_iv  (tile_folded[dct_mode], dct_mode, transform_size);
							}
							System.arraycopy(tile_out[dct_mode], 0, dct_data[tileY][tileX][dct_mode], 0, tile_out[dct_mode].length);
						}
						if ((globalDebugLevel > 0) && (debug_tileX == tileX) && (debug_tileY == tileY)) {
							ShowDoubleFloatArrays sdfa_instance = new ShowDoubleFloatArrays(); // just for debugging?
							sdfa_instance.showArrays(tile_in,  n2, n2, "tile_in_x"+tileX+"_y"+tileY);
							String [] titles = {"CC","SC","CS","SS"};
							sdfa_instance.showArrays(tile_folded,  transform_size, transform_size, true, "folded_x"+tileX+"_y"+tileY, titles);
							if (globalDebugLevel > 0) {
								sdfa_instance.showArrays(tile_out,     transform_size, transform_size, true, "clt_x"+tileX+"_y"+tileY, titles);
							}
						}
					}
				}
			};
		}
		startAndJoin(threads);
		return dct_data;
	}

	public double [] iclt_2d(  // USED in lwir
			final double [][][][] dct_data,  // array [tilesY][tilesX][4][dct_size*dct_size]
//			final int             transform_size,
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
		final int width=  tilesX * transform_size;
		final int height= tilesY * transform_size;
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
						DttRad2 dtt = new DttRad2(transform_size);
						dtt.set_window(window_type);
						double [] tile_in =   new double [transform_size * transform_size];
						double [] tile_dct;
						double [] tile_mdct;
						int tileY,tileX;
						int n2 = transform_size * 2;
						int n_half = transform_size / 2;
						int lastY = tilesY-1;
						int lastX = tilesX-1;
						int offset = n_half * (transform_size * tilesX) + n_half;
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
										tile_dct = dtt.dttt_iv  (tile_in, idct_mode, transform_size);
									}
									tile_mdct = dtt.unfold_tile(tile_dct, transform_size, dct_mode); // mode=0 - DCCT
									if ((tileY >0) && (tileX > 0) && (tileY < lastY) && (tileX < lastX)) { // fast, no extra checks
										for (int i = 0; i < n2;i++){
											//									int start_line = ((tileY*dct_size + i) *(tilesX+1) + tileX)*dct_size;
											int start_line = ((tileY*transform_size + i) * tilesX + tileX)*transform_size - offset;
											for (int j = 0; j<n2;j++) {
												dpixels[start_line + j] += debug_scale * tile_mdct[n2 * i + j]; // add (cc+sc+cs+ss)/4
											}
										}
									} else { // be careful with margins
										for (int i = 0; i < n2;i++){
											if (	((tileY > 0) && (tileY < lastY)) ||
													((tileY == 0) && (i >= n_half)) ||
													((tileY == lastY) && (i < (n2 - n_half)))) {
												int start_line = ((tileY*transform_size + i) * tilesX + tileX)*transform_size  - offset;
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

	public double [] iclt_2d_debug_gpu( // not used in lwir
			final double [][][][] dct_data,  // array [tilesY][tilesX][4][dct_size*dct_size]
//			final int             transform_size,
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
		final double [][] dbg_tile = ((debug_tileY >=0) && (debug_tileX >=0))? dct_data[debug_tileY][debug_tileX]:null;

//		final int width=  (tilesX+1)*dct_size;
//		final int height= (tilesY+1)*dct_size;
		final int width=  tilesX * transform_size;
		final int height= tilesY * transform_size;
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
						DttRad2 dtt = new DttRad2(transform_size);
						dtt.set_window(window_type);
						double [] tile_in =   new double [transform_size * transform_size];
						double [] tile_dct;
						double [] tile_mdct;
						int tileY,tileX;
						int n2 = transform_size * 2;
						int n_half = transform_size / 2;
						int lastY = tilesY-1;
						int lastX = tilesX-1;
						int offset = n_half * (transform_size * tilesX) + n_half;
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
											for (int i = 0; i < transform_size; i++) {
												for (int j = 0; j < transform_size; j++) {
													System.out.print(String.format("%11.5f ", tile_in[transform_size * j + i]));
												}
												System.out.println();
											}
											System.out.println();
										}
										tile_dct = dtt.dttt_iv  (tile_in, idct_mode, transform_size);
										if (dbg_tile) {
											System.out.println("--- AFTER idct: dct_mode="+dct_mode+", idct_mode="+idct_mode+" tile_dct (not transposed):");
											for (int i = 0; i < transform_size; i++) {
												for (int j = 0; j < transform_size; j++) {
													System.out.print(String.format("%11.5f ", tile_dct[transform_size * i + j]));
												}
												System.out.println();
											}
											System.out.println();
										}
									}
									tile_mdct = dtt.unfold_tile(tile_dct, transform_size, dct_mode); // mode=0 - DCCT

									if (dbg_tile) {
										System.out.println("--- MDCT: ");
										for (int i = 0; i < transform_size * 2; i++) {
											for (int j = 0; j < transform_size * 2; j++) {
												System.out.print(String.format("%11.5f ", tile_mdct[transform_size * 2 * i + j]));
											}
											System.out.println();
										}
										System.out.println();
									}


									if ((tileY >0) && (tileX > 0) && (tileY < lastY) && (tileX < lastX)) { // fast, no extra checks
										for (int i = 0; i < n2;i++){
											int start_line = ((tileY*transform_size + i) * tilesX + tileX)*transform_size - offset;
											for (int j = 0; j<n2;j++) {
												dpixels[start_line + j] += debug_scale * tile_mdct[n2 * i + j]; // add (cc+sc+cs+ss)/4
											}
										}
									} else { // be careful with margins
										for (int i = 0; i < n2;i++){
											if (	((tileY > 0) && (tileY < lastY)) ||
													((tileY == 0) && (i >= n_half)) ||
													((tileY == lastY) && (i < (n2 - n_half)))) {
												int start_line = ((tileY*transform_size + i) * tilesX + tileX)*transform_size  - offset;
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
										int start_line = ((tileY*transform_size + i) * tilesX + tileX)*transform_size - offset;
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
	public double [][] combineRBGATiles(  // USED in lwir
			final double [][][][] texture_tiles,  // array [tilesY][tilesX][4][4*transform_size] or [tilesY][tilesX]{null}
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
		boolean set_has_weight = false; // not used
		for (int i = 0; (i < tilesY) && !set_has_weight; i++){
			for (int j = 0; (j < tilesX) && !set_has_weight; j++){
				if (texture_tiles[i][j] != null) {
					set_has_weight = true;
					has_weights = texture_tiles[i][j].length > 4;
				}
			}
		}

		final double [][] dpixels = new double["RGBA".length()+(has_weights? (numSensors+4): 0)][width*height]; // assuming java initializes them to 0
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
						int offset = n_half * (transform_size * tilesX) + n_half; // 4 pixels left and down (right/up when subtracted below)
						for (int nTile = ai.getAndIncrement(); nTile < tiles_list[nser.get()].length; nTile = ai.getAndIncrement()) {
							tileX = tiles_list[nser.get()][nTile][0];
							tileY = tiles_list[nser.get()][nTile][1];
							double [][] texture_tile =texture_tiles[tileY][tileX];
							if (texture_tile != null) {
								if (overlap) {
									if ((tileY >0) && (tileX > 0) && (tileY < lastY) && (tileX < lastX)) { // fast, no extra checks - ignore first/last rows and columns
										for (int i = 0; i < n2; i++){
											int start_line = ((tileY*transform_size + i) * tilesX + tileX)*transform_size - offset;
											for (int chn = 0; chn < texture_tile.length; chn++) {
												int schn = chn;
												if (isMonochrome() && (chn<3)) {
													schn = MONO_CHN; // clone green to red and blue output
												}
												if (texture_tile[schn] == null) {
													dpixels[chn] = null;
												} else {
													// should it be better to multiply each color by alpha before accumulating? No, it is already windowed!
													if ((chn != 3) || !sharp_alpha) {
														for (int j = 0; j<n2;j++) {
															dpixels[chn][start_line + j] += texture_tile[schn][n2 * i + j];
														}
													} else if ((i >= n_half) && (i < (n2-n_half))) {
														for (int j = n_half; j < (n2 - n_half); j++) { // not used in lwir
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
															for (int j = n_half; j < (n2 - n_half); j++) { // not used in lwir
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
									for (int i = 0; i < n2;i++){ // not used in lwir
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




	public double [][][][] clt_shiftXY( // not used in lwir
			final double [][][][] dct_data,  // array [tilesY][tilesX][4][dct_size*dct_size]
//			final int             transform_size,
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
		final double [] cos_hor =  new double [transform_size*transform_size];
		final double [] sin_hor =  new double [transform_size*transform_size];
		final double [] cos_vert = new double [transform_size*transform_size];
		final double [] sin_vert = new double [transform_size*transform_size];
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

		if (globalDebugLevel > 1){
			ShowDoubleFloatArrays sdfa_instance = new ShowDoubleFloatArrays(); // just for debugging?
			String [] titles = {"cos_hor","sin_hor","cos_vert","sin_vert"};
			double [][] cs_dbg = {cos_hor, sin_hor, cos_vert, sin_vert};
			sdfa_instance.showArrays(cs_dbg,  transform_size, transform_size, true, "shift_cos_sin", titles);
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

	public double [][][][] clt_correlate( // not used in lwir
			final double [][][][] data1,  // array [tilesY][tilesX][4][dct_size*dct_size]
			final double [][][][] data2,  // array [tilesY][tilesX][4][dct_size*dct_size]
//			final int             transform_size,
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
		final int dct_len = transform_size * transform_size;
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
								rslt[tileY][tileX][n][i] = 0;
								for (int k=0; k<4; k++){
									if (ZI[n][k] < 0)
										rslt[tileY][tileX][n][i] -=
											data1[tileY][tileX][-ZI[n][k]][i] * data2[tileY][tileX][k][i];
									else
										rslt[tileY][tileX][n][i] +=
										data1[tileY][tileX][ZI[n][k]][i] * data2[tileY][tileX][k][i];
								}
								rslt[tileY][tileX][n][i] *= scale;
							}
						}
						if ((globalDebugLevel > 0) && (debug_tileX == tileX) && (debug_tileY == tileY)) {
							ShowDoubleFloatArrays sdfa_instance = new ShowDoubleFloatArrays(); // just for debugging?
							String [] titles = {"CC","SC","CS","SS"};
							sdfa_instance.showArrays(data1[tileY][tileX], transform_size, transform_size, true, "data1_x"+tileX+"_y"+tileY, titles);
							sdfa_instance.showArrays(data2[tileY][tileX], transform_size, transform_size, true, "data2_x"+tileX+"_y"+tileY, titles);
							sdfa_instance.showArrays(rslt[tileY][tileX],  transform_size, transform_size, true, "rslt_x"+ tileX+"_y"+tileY, titles);
						}
					}
				}
			};
		}
		startAndJoin(threads);
		return rslt;
	}

	/**
	 * Get frequency-domain representation of the LPF
	 * @param sigma blurring in pixels
	 * @return double array of the filter, 64 long for 8-pixel DTT
	 */
	public double [] doubleGetCltLpfFd(
			double   sigma) {
		DttRad2 dtt = new DttRad2(transform_size);
		double [] clt_fd = dtt.dttt_iiie(getLpf(sigma));
		int l = clt_fd.length;
		double []   lpf_flat = new double [l];
		for (int j = 0; j < l; j++) {
			lpf_flat[j] = (clt_fd[j]*2*transform_size);
		}
		return lpf_flat;
	}
	
	public static double [] doubleGetCltLpfFd(
			double   sigma,
			int transform_size) {
		DttRad2 dtt = new DttRad2(transform_size);
		double [] clt_fd = dtt.dttt_iiie(getLpf(sigma, transform_size));
		int l = clt_fd.length;
		double []   lpf_flat = new double [l];
		for (int j = 0; j < l; j++) {
			lpf_flat[j] = (clt_fd[j]*2*transform_size);
		}
		return lpf_flat;
	}
	
	
	/**
	 * Get frequency-domain representation of the LPF (version for the GPU, in floats)
	 * @param sigma2 squared Gaussian sigma in pixels
	 * @return float array of the filter, 64 long for 8-pixel DTT
	 */
	public float [] floatGetCltLpfFd(
			double   sigma2) {
		DttRad2 dtt = new DttRad2(transform_size);
		double [] clt_fd = dtt.dttt_iiie(getLpf(sigma2));
		int l = clt_fd.length;
		float []   lpf_flat = new float [l];
		for (int j = 0; j < l; j++) {
			lpf_flat[j] = (float) (clt_fd[j]*2*transform_size);
		}
		return lpf_flat;
	}

	public static float [] floatGetCltLpfFd(
			double   sigma2,
			int      transform_size) {
		DttRad2 dtt = new DttRad2(transform_size);
		double [] clt_fd = dtt.dttt_iiie(getLpf(sigma2, transform_size));
		int l = clt_fd.length;
		float []   lpf_flat = new float [l];
		for (int j = 0; j < l; j++) {
			lpf_flat[j] = (float) (clt_fd[j]*2*transform_size);
		}
		return lpf_flat;
	}
	
	
	/**
	 * Get frequency-domain representation of the LPF (version for the GPU, in floats)
	 * @param sigma Gaussian sigma in pixels
	 * @return float array of the filter, 64 long for 8-pixel DTT
	 */
	public float [] floatGetCltHpfFd(
			double   sigma) {
		DttRad2 dtt = new DttRad2(transform_size);
		double [] clt_fd = (sigma == 0.0)? (new double[transform_size*transform_size]) : dtt.dttt_iiie(getLpf(sigma * sigma));
		int l = clt_fd.length;
		float []   hpf_flat = new float [l];
		for (int j = 0; j < l; j++) {
			hpf_flat[j] = (float) (1.0 - clt_fd[j]*2*transform_size);
		}
		return hpf_flat;
	}
	
	
	/**
	 * Get pixel-domain representation of the LPF
	 * @param sigma2  squared Gaussian sigma in pixels
	 * @return double array of the filter, 64 long for 8-pixel DTT
	 */

	public double [] getLpf(
			double   sigma2) // sigma squared
	{
		int transform_len = transform_size * transform_size;
		final double [] filter_direct= new double[transform_len];
		if (sigma2 == 0) {
			filter_direct[0] = 1.0;
			for (int i= 1; i<filter_direct.length;i++) {
                filter_direct[i] =0;
			}
		} else {
			for (int i = 0; i < transform_size; i++){
				for (int j = 0; j < transform_size; j++){
					filter_direct[i*transform_size+j] = Math.exp(-(i*i+j*j)/(2*sigma2));
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
		return filter_direct;
	}

	public static double [] getLpf(
			double   sigma2,
			int      transform_size) // sigma squared
	{
		int transform_len = transform_size * transform_size;
		final double [] filter_direct= new double[transform_len];
		if (sigma2 == 0) {
			filter_direct[0] = 1.0;
			for (int i= 1; i<filter_direct.length;i++) {
                filter_direct[i] =0;
			}
		} else {
			for (int i = 0; i < transform_size; i++){
				for (int j = 0; j < transform_size; j++){
					filter_direct[i*transform_size+j] = Math.exp(-(i*i+j*j)/(2*sigma2));
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
		return filter_direct;
	}
	
	
	
	
	/**
	 * Get frequency-domain representation of the LoG (version for the GPU, in floats)
	 * @param sigma Gaussian sigma in pixels
	 * @return float array of the filter, 64 long for 8-pixel DTT
	 */
	public float [] floatGetCltLoGFd(
			double   sigma) {
		DttRad2 dtt = new DttRad2(transform_size);
		double [] clt_fd = dtt.dttt_iiie(getLoG(sigma));
		int l = clt_fd.length;
		float []   log_flat = new float [l];
		for (int j = 0; j < l; j++) {
			log_flat[j] = (float) (clt_fd[j]*2*transform_size);
		}
		return log_flat;
	}
	
	/**
	 * Get pixel-domain representation of the LoG
	 * @param sigma Gaussian sigma in pixels
	 * @return double array of the filter, 64 long for 8-pixel DTT
	 */
	public double [] getLoG(
			double   sigma)
	{
		int transform_len = transform_size * transform_size;
		final double sigma2 = sigma*sigma;
		final double sigma4 = sigma2*sigma2;
		final double [] filter_direct= new double[transform_len];
		if (sigma == 0) {
			filter_direct[0] = 1.0;
			for (int i= 1; i<filter_direct.length;i++) {
                filter_direct[i] =0;
			}
		} else {
			for (int i = 0; i < transform_size; i++){
				for (int j = 0; j < transform_size; j++){
//https://homepages.inf.ed.ac.uk/rbf/HIPR2/log.htm					
					filter_direct[i*transform_size+j] =
							-1.0/(Math.PI * sigma4)*(1.0 - (i*i+j*j)/(2*sigma2))*
							Math.exp(-(i*i+j*j)/(2*sigma2));
				}
			}
		}
		(new ShowDoubleFloatArrays()).showArrays(
				filter_direct,
				8,
				8,
				"log_direct-"+sigma);
		// normalize
		double sum2 = 0;
		for (int i = 0; i < transform_size; i++){
			for (int j = 0; j < transform_size; j++){
				double d = 	filter_direct[i*transform_size+j];
				d*=d;
				d*=Math.cos(Math.PI*i/(2*transform_size))*Math.cos(Math.PI*j/(2*transform_size));
				if (i > 0) d*= 2.0;
				if (j > 0) d*= 2.0;
				sum2 +=d;
			}
		}
		double sum = Math.sqrt(sum2);
		for (int i = 0; i<filter_direct.length; i++){
			filter_direct[i] /= sum;
		}
		System.out.println("getLoG("+sigma+") sum="+sum);
		/*
		sum2 = 0;
		for (int i = 0; i < transform_size; i++){
			for (int j = 0; j < transform_size; j++){
				double d = 	filter_direct[i*transform_size+j];
				d*=d;
				d*=Math.cos(Math.PI*i/(2*transform_size))*Math.cos(Math.PI*j/(2*transform_size));
				if (i > 0) d*= 2.0;
				if (j > 0) d*= 2.0;
				sum2 +=d;
			}
		}
		*/
		(new ShowDoubleFloatArrays()).showArrays(
				filter_direct,
				8,
				8,
				"log_direct_norm-"+sigma);
		
		
		
		return filter_direct;
	}


	public void clt_lpf(  // USED in lwir
			final double          sigma,
			final double [][][][] clt_data,
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

		final double [] filter =  doubleGetCltLpfFd(sigma);


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

			// code to simulate old commented out
			double [] filter_direct = getLpf(sigma);
			double [] dbg_filter0 = new double [filter.length];
			for (int i=0; i < filter.length;i++) {
				dbg_filter0[i]= filter[i] / (2*transform_size); //simulating old debug data
			}
			double [] dbg_filter= (new DttRad2(transform_size)).dttt_ii(dbg_filter0);
			// end of code to simulate old commented out

			ShowDoubleFloatArrays sdfa_instance = new ShowDoubleFloatArrays(); // just for debugging?
			double [][] ff = {filter_direct,filter,dbg_filter};
			sdfa_instance.showArrays(ff,  transform_size,transform_size, true, "filter_lpf");
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


	public double [][][][] clt_lpf_copy(  // USED in lwir
			final double          sigma,
			final double [][][][] clt_data,
			final int             threadsMax,     // maximal number of threads to launch
			final int             globalDebugLevel)
	{
		if (clt_data == null) {
			System.out.println("clt_lpf(): clt_data=null");
			return null;
		}
		final int tilesY=clt_data.length;
		final int tilesX=clt_data[0].length;
		final int nTiles=tilesX*tilesY;
		final double [][][][] clt_out = new double [tilesY][tilesX][][];

		final double [] filter =  doubleGetCltLpfFd(sigma);


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

			// code to simulate old commented out
			double [] filter_direct = getLpf(sigma);
			double [] dbg_filter0 = new double [filter.length];
			for (int i=0; i < filter.length;i++) {
				dbg_filter0[i]= filter[i] / (2*transform_size); //simulating old debug data
			}
			double [] dbg_filter= (new DttRad2(transform_size)).dttt_ii(dbg_filter0);
			// end of code to simulate old commented out

			ShowDoubleFloatArrays sdfa_instance = new ShowDoubleFloatArrays(); // just for debugging?
			double [][] ff = {filter_direct,filter,dbg_filter};
			sdfa_instance.showArrays(ff,  transform_size,transform_size, true, "filter_lpf");
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
							clt_out[tileY][tileX] = new double [4][filter.length];
							for (int n = 0; n < 4; n++){
								for (int i = 0; i < filter.length; i++){
									clt_out[tileY][tileX][n][i]= clt_data[tileY][tileX][n][i] * filter[i];
								}
							}
						}
					}
				}
			};
		}
		startAndJoin(threads);
		return clt_out;
	}




	public void clt_dtt2( // transform dcct2, dsct2, dcst2, dsst2 // not used in lwir
			final double [][][][] data,
			final boolean         transpose, // when doing inverse transform, the data comes in transposed form, so CS <->SC
			final int             threadsMax,     // maximal number of threads to launch
			final int             globalDebugLevel)
	{
		final int tilesY=data.length;
		final int tilesX=data[0].length;
		final int nTiles=tilesX*tilesY;
///		final int transform_size = (int) Math.round(Math.sqrt(data[0][0][0].length));
		final Thread[] threads = newThreadArray(threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);

		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				@Override
				public void run() {
					DttRad2 dtt = new DttRad2(transform_size);
					int tileY,tileX;
					for (int nTile = ai.getAndIncrement(); nTile < nTiles; nTile = ai.getAndIncrement()) {
						tileY = nTile/tilesX;
						tileX = nTile - tileY * tilesX;
						for (int quadrant = 0; quadrant < 4; quadrant++){
							int mode = transpose ? (((quadrant << 1) & 2) | ((quadrant >> 1) & 1)) : quadrant;
							data[tileY][tileX][quadrant] = dtt.dttt_iie(data[tileY][tileX][quadrant], mode, transform_size);
						}
					}
				}
			};
		}
		startAndJoin(threads);
	}

	public double [][][] clt_corr_quad( // combine 4 correlation quadrants after DTT2  // not used in lwir
			final double [][][][] data,
			final int             threadsMax,     // maximal number of threads to launch
			final int             globalDebugLevel)
	{
		final int tilesY=data.length;
		final int tilesX=data[0].length;
		final int nTiles=tilesX*tilesY;
///		final int transform_size = (int) Math.round(Math.sqrt(data[0][0][0].length));
		final int rslt_size=transform_size*2-1;

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
						rslt[tileY][tileX][rslt_size*transform_size - transform_size] = scale * data[tileY][tileX][0][0]; // center
						for (int j = 1; j < transform_size; j++) { //  for i == 0
							rslt[tileY][tileX][rslt_size*transform_size - transform_size + j] = scale * (data[tileY][tileX][0][j] + data[tileY][tileX][1][j-1]);
							rslt[tileY][tileX][rslt_size*transform_size - transform_size - j] = scale * (data[tileY][tileX][0][j] - data[tileY][tileX][1][j-1]);
						}
						for (int i = 1; i < transform_size; i++) {
							rslt[tileY][tileX][rslt_size*(transform_size + i) - transform_size] =
									scale * (data[tileY][tileX][0][i*transform_size] + data[tileY][tileX][2][(i-1)*transform_size]);
							rslt[tileY][tileX][rslt_size*(transform_size - i) - transform_size] =
									scale * (data[tileY][tileX][0][i*transform_size] - data[tileY][tileX][2][(i-1)*transform_size]);
							for (int j = 1; j < transform_size; j++) {
								rslt[tileY][tileX][rslt_size*(transform_size + i) - transform_size + j] =
										scale * (data[tileY][tileX][0][i*    transform_size + j] +
												 data[tileY][tileX][1][i*    transform_size + j - 1] +
												 data[tileY][tileX][2][(i-1)*transform_size + j] +
												 data[tileY][tileX][3][(i-1)*transform_size + j - 1]);

								rslt[tileY][tileX][rslt_size*(transform_size + i) - transform_size - j] =
										scale * ( data[tileY][tileX][0][i*    transform_size + j] +
												 -data[tileY][tileX][1][i*    transform_size + j - 1] +
												  data[tileY][tileX][2][(i-1)*transform_size + j] +
												 -data[tileY][tileX][3][(i-1)*transform_size + j - 1]);
								rslt[tileY][tileX][rslt_size*(transform_size - i) - transform_size + j] =
										scale * (data[tileY][tileX][0][i*    transform_size + j] +
												 data[tileY][tileX][1][i*    transform_size + j - 1] +
												 -data[tileY][tileX][2][(i-1)*transform_size + j] +
												 -data[tileY][tileX][3][(i-1)*transform_size + j - 1]);
								rslt[tileY][tileX][rslt_size*(transform_size - i) - transform_size - j] =
										scale * (data[tileY][tileX][0][i*    transform_size + j] +
												 -data[tileY][tileX][1][i*    transform_size + j - 1] +
												 -data[tileY][tileX][2][(i-1)*transform_size + j] +
												 data[tileY][tileX][3][(i-1)*transform_size + j - 1]);
							}
						}
					}
				}
			};
		}
		startAndJoin(threads);
		return rslt;
	}

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

	public static float [][] corr_td_dbg(
			final float [][][][] fcorr_td,
			// if 0 - fcorr_combo_td = new float[4][tilesY][tilesX][];
			// if > 0 - fcorr_td =       new float[tilesY][tilesX][num_slices][];
			final int            num_slices,
			final int            transform_size,
			final int []         wh, // should be initialized as int[2];
			final int            threadsMax)     // maximal number of threads to launch
	{
		final int tilesY = (num_slices == 0) ? fcorr_td[0].length : fcorr_td.length;
		final int tilesX = (num_slices == 0) ? fcorr_td[0][0].length : fcorr_td[0].length;
		final int nTiles = tilesX*tilesY;
		final int width =  tilesX * 2 * transform_size;
		final int height = tilesY * 2 * transform_size;
		if (wh != null) {
			wh[0] = width;
			wh[1] = height;
		}
		final int fnum_slices = (num_slices == 0) ? fcorr_td.length : num_slices;
		final int transform_len = transform_size*transform_size; // 64
		float [][] dbg_img = new float [fnum_slices][width * height];
		for (int i = 0; i < dbg_img.length; i++) {
			Arrays.fill(dbg_img[i], Float.NaN);
		}
		final Thread[] threads = newThreadArray(threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				@Override
				public void run() {
					for (int nTile = ai.getAndIncrement(); nTile < nTiles; nTile = ai.getAndIncrement()) {
						int tileY = nTile/tilesX;
						int tileX = nTile - tileY * tilesX;
						if ((num_slices == 0) || (fcorr_td[tileY][tileX] != null)) {
							for (int slice = 0; slice < fnum_slices; slice ++) {
								float [] ftile = (num_slices > 0) ? fcorr_td[tileY][tileX][slice] : fcorr_td[slice][tileY][tileX];
								if (ftile != null) {
									for (int qy = 0; qy < 2; qy++) {
										for (int qx = 0; qx < 2; qx++) {
											for (int ty = 0; ty < transform_size; ty++) {
												int py = (tileY * 2 + qy) * transform_size + ty;
												int px = (tileX * 2 + qx) * transform_size;
												System.arraycopy(
														ftile,
														(2 * qy + qx) * transform_len + transform_size * ty,
														dbg_img[slice],
														py * width + px,
														transform_size);
												/*
												for (int tilesY = 0; tx < transform_size; tx++) {
													int px = (tileX * 2 + qx) * transform_size + tx;
													dbg_img[slice][py * width + px] = ftile[(2 * qy + qx) * transform_len + transform_size * ty + tx];
												}
												*/
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
		return dbg_img;
	}
	
	
	
//	final float  [][][][]     fcorr_td =       new float[tilesY][tilesX][][];
//	final float  [][][][]     fcorr_combo_td = new float[4][tilesY][tilesX][];

	public static void corr_td_normalize(
			final float [][][][] fcorr_td, // will be updated
			// if 0 - fcorr_combo_td = new float[4][tilesY][tilesX][];
			// if > 0 - fcorr_td =       new float[tilesY][tilesX][num_slices][];
			final int            num_slices,
			final int            transform_size,
			final double         fat_zero_abs,
			final double         output_amplitude,
			final int            threadsMax)     // maximal number of threads to launch
	{
		final double fat_zero_abs2 = fat_zero_abs * fat_zero_abs; 
		final int tilesY = (num_slices == 0) ? fcorr_td[0].length : fcorr_td.length;
		final int tilesX = (num_slices == 0) ? fcorr_td[0][0].length : fcorr_td[0].length;
		final int nTiles = tilesX*tilesY;
		final int fnum_slices = (num_slices == 0) ? fcorr_td.length : num_slices;
		final int transform_len = transform_size*transform_size; // 64
		final Thread[] threads = newThreadArray(threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				@Override
				public void run() {
					for (int nTile = ai.getAndIncrement(); nTile < nTiles; nTile = ai.getAndIncrement()) {
						int tileY = nTile/tilesX;
						int tileX = nTile - tileY * tilesX;
						if ((num_slices == 0) || (fcorr_td[tileY][tileX] != null)) {
							for (int slice = 0; slice < fnum_slices; slice ++) {
								float [] ftile = (num_slices > 0) ? fcorr_td[tileY][tileX][slice] : fcorr_td[slice][tileY][tileX];
								if (ftile != null) {
									for (int i = 0; i < transform_len; i++) {
										double s2 = fat_zero_abs2;
										for (int q = 0; q < 4; q++) {
											double d = ftile[q * transform_len + i]; 
											s2 += d*d;
										}
										double k = output_amplitude/Math.sqrt(s2);
										for (int q = 0; q < 4; q++) {
											ftile[q * transform_len + i] *= k;
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

	public static void corr_td_lpf(
			final float [][][][] fcorr_td, // will be updated
			// if 0 - fcorr_combo_td = new float[4][tilesY][tilesX][];
			// if > 0 - fcorr_td =       new float[tilesY][tilesX][num_slices][];
			final int            num_slices,
			final int            transform_size,
			final double         corr_sigma,
			final int            threadsMax)     // maximal number of threads to launch
	{
		final float [] filter =     ImageDtt.floatGetCltLpfFd(corr_sigma, transform_size);
		final int tilesY = (num_slices == 0) ? fcorr_td[0].length : fcorr_td.length;
		final int tilesX = (num_slices == 0) ? fcorr_td[0][0].length : fcorr_td[0].length;
		final int nTiles = tilesX*tilesY;
		final int fnum_slices = (num_slices == 0) ? fcorr_td.length : num_slices;
		final int transform_len = transform_size*transform_size; // 64
		final Thread[] threads = newThreadArray(threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				@Override
				public void run() {
					for (int nTile = ai.getAndIncrement(); nTile < nTiles; nTile = ai.getAndIncrement()) {
						int tileY = nTile/tilesX;
						int tileX = nTile - tileY * tilesX;
						if ((num_slices == 0) || (fcorr_td[tileY][tileX] != null)) {
							for (int slice = 0; slice < fnum_slices; slice ++) {
								float [] ftile = (num_slices > 0) ? fcorr_td[tileY][tileX][slice] : fcorr_td[slice][tileY][tileX];
								if (ftile != null) {
									for (int q = 0; q < 4; q++) {
										for (int i = 0; i < transform_len; i++) {
											ftile[q * transform_len + i] *= filter[i];
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
	
	
	
	/**
	 * Convert normalized individual pairs correlations or composite (float, from GPU) from
	 * transform domain to pixel domain (debugging irregularities in correlations)
	 * @param fcorr_td Individual or composite correlations (different format)
	 * @param num_slices  If 0 - fcorr_td is float[4][tilesY][tilesX][256], if > 0:float[tilesY][tilesX][num_slices][256]
	 * @param transform_size 8
	 * @param threadsMax max number of threads
	 * @return tiles of [4][64] pixel-domain correlations, first indices same as in fcorr_td
	 */
	public static double [][][][][] corr_td_dct2(
			final float [][][][] fcorr_td, // normalized TD representation
			// if 0 - fcorr_combo_td = new float[4][tilesY][tilesX][]; // last index - [4 * 64]
			// if > 0 - fcorr_td =       new float[tilesY][tilesX][num_slices][];
			final int            num_slices,
			final int            transform_size,
			final int            threadsMax)     // maximal number of threads to launch
	{
		final int tilesY = (num_slices == 0) ? fcorr_td[0].length : fcorr_td.length;
		final int tilesX = (num_slices == 0) ? fcorr_td[0][0].length : fcorr_td[0].length;
		final int nTiles = tilesX*tilesY;
		final int fnum_slices = (num_slices == 0) ? fcorr_td.length : num_slices;
		final int transform_len = transform_size*transform_size; // 64
		final double [][][][][] pd_data = (num_slices == 0) ? (new double[4][tilesY][tilesX][][]):(new double [tilesY][tilesX][][][]);
		
		final Thread[] threads = newThreadArray(threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				@Override
				public void run() {
					DttRad2 dtt = new DttRad2(transform_size);
					for (int nTile = ai.getAndIncrement(); nTile < nTiles; nTile = ai.getAndIncrement()) {
						int tileY = nTile/tilesX;
						int tileX = nTile - tileY * tilesX;
						if ((num_slices == 0) || (fcorr_td[tileY][tileX] != null)) {
							if (num_slices > 0) {
								pd_data[tileY][tileX] = new double[num_slices][][];
							}
							for (int slice = 0; slice < fnum_slices; slice ++) {
								float [] ftile = (num_slices > 0) ? fcorr_td[tileY][tileX][slice] : fcorr_td[slice][tileY][tileX];
								if (ftile != null) {
									// convert to double
									double [][] dtile = new double [4][transform_len];
									for (int q = 0; q < dtile.length; q++) {
										for (int i = 0; i < transform_len; i++) {
											dtile[q][i] = ftile[q*transform_len + i];
										}
									}
									for (int quadrant = 0; quadrant < 4; quadrant++){
										int mode = ((quadrant << 1) & 2) | ((quadrant >> 1) & 1); // transpose
										dtile[quadrant] = dtt.dttt_iie(dtile[quadrant], mode, transform_size);
									}
									if (num_slices == 0) {
										pd_data[slice][tileY][tileX] = dtile;
									} else {
										pd_data[tileY][tileX][slice] = dtile;
									}
								}
							}
						}
					}
				}
			};
		}
		startAndJoin(threads);
		return pd_data;
	}

	
	public static float [][][][] corr_transpose(
			final float [][][][] fcorr_td, // normalized TD representation
			// if 0 - fcorr_combo_td = new float[4][tilesY][tilesX][]; // last index - [4 * 64]
			// if > 0 - fcorr_td =       new float[tilesY][tilesX][num_slices][];
			final int            num_slices,
			final int            transform_size,
			final int            threadsMax)     // maximal number of threads to launch
	{
		final int tilesY = (num_slices == 0) ? fcorr_td[0].length : fcorr_td.length;
		final int tilesX = (num_slices == 0) ? fcorr_td[0][0].length : fcorr_td[0].length;
		final int nTiles = tilesX*tilesY;
		final int fnum_slices = (num_slices == 0) ? fcorr_td.length : num_slices;
		final int transform_len = transform_size*transform_size; // 64
		final float [][][][] fcorr_transposed= (num_slices == 0) ? (new float[4][tilesY][tilesX][]):(new float [tilesY][tilesX][][]);
		
		final Thread[] threads = newThreadArray(threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				@Override
				public void run() {
					for (int nTile = ai.getAndIncrement(); nTile < nTiles; nTile = ai.getAndIncrement()) {
						int tileY = nTile/tilesX;
						int tileX = nTile - tileY * tilesX;
						if ((num_slices == 0) || (fcorr_td[tileY][tileX] != null)) {
							if (num_slices > 0) {
								fcorr_transposed[tileY][tileX] = new float[num_slices][];
							}
							for (int slice = 0; slice < fnum_slices; slice ++) {
								float [] ftile = (num_slices > 0) ? fcorr_td[tileY][tileX][slice] : fcorr_td[slice][tileY][tileX];
								if (ftile != null) {
									float [] ftile_transposed = new float [4*transform_len];
									for (int q = 0; q < 4; q++) {
										for (int i = 0; i < transform_size; i++) {
											for (int j = 0; j < transform_size; j++) {
												ftile_transposed[q*transform_len + i * transform_size + j] =	
											ftile[q*transform_len + j * transform_size + i];
											}
										}
									}
									if (num_slices == 0) {
										fcorr_transposed[slice][tileY][tileX] = ftile_transposed;
									} else {
										fcorr_transposed[tileY][tileX][slice] = ftile_transposed;
									}
								}
							}
						}
					}
				}
			};
		}
		startAndJoin(threads);
		return fcorr_transposed;
	}
	
	
	
	/**
	 * Convert pixel-domain CLT correlation (after DTT-II) to 15x15 2D correlation tiles.
	 * @param pd_data
	 * @param num_slices
	 * @param transform_size
	 * @param threadsMax
	 * @return
	 */
	public static double [][][][] corr_pd_unfold(
			final double [][][][][] pd_data, // pixel-domain correlation representation
			// if 0 - pd_data = new double[4][tilesY][tilesX][4][64]; // last index - [4 * 64]
			// if > 0 - pd_data =       new double[tilesY][tilesX][num_slices][4][64];
			final int            num_slices,
			final int            transform_size,
			final int            threadsMax)     // maximal number of threads to launch
	{
		final int tilesY = (num_slices == 0) ? pd_data[0].length : pd_data.length;
		final int tilesX = (num_slices == 0) ? pd_data[0][0].length : pd_data[0].length;
		final int nTiles = tilesX*tilesY;
		final int fnum_slices = (num_slices == 0) ? pd_data.length : num_slices;
//		final int transform_len = transform_size*transform_size; // 64
		final double [][][][] corr = (num_slices == 0) ? (new double[4][tilesY][tilesX][]):(new double [tilesY][tilesX][][]);
		final Thread[] threads = newThreadArray(threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				@Override
				public void run() {
					DttRad2 dtt = new DttRad2(transform_size);
					for (int nTile = ai.getAndIncrement(); nTile < nTiles; nTile = ai.getAndIncrement()) {
						int tileY = nTile/tilesX;
						int tileX = nTile - tileY * tilesX;
						if (num_slices == 0) {
							for (int slice = 0; slice < fnum_slices; slice ++) {
								double [][] dtile = pd_data[slice][tileY][tileX];
								if (dtile != null) {
									corr[slice][tileY][tileX] = dtt.corr_unfold_tile(dtile,	transform_size);
								}
							}
						} else if (pd_data[tileY][tileX] != null) {
							corr[tileY][tileX] = new double [num_slices][];
							for (int slice = 0; slice < fnum_slices; slice ++) {
								double [][] dtile = pd_data[tileY][tileX][slice];
								if (dtile != null) {
									corr[tileY][tileX][slice] = dtt.corr_unfold_tile(dtile,	transform_size);
								}
							}
						}
					}
				}
			};
		}
		startAndJoin(threads);
		return corr;
	}
	
	//TODO: remove dependence on quad, remove @Deprecated annotation 
	@Deprecated
	public static void clt_process_pd_correlations( // process correlations already prepared in dcorr and/or dcorr_combo
			final QuadCLT             quadCLT,
			final ImageDttParameters  imgdtt_params,   // Now just extra correlation parameters, later will include, most others
			// both arrays should have same non-null tiles
			final double  [][][][]    dcorr,           // [tilesY][tilesX][pair][225] pixel domain representation of 6 corr pairs SHOULD not == null
			final double  [][][][]    dcorr_combo,     // [4][tilesY][tilesX][225] PD of combo corrs: qud, cross, hor,vert MAY BE null
			// each of the top elements may be null to skip particular combo type
			final double [][][][]     corr_tiles,      // [tilesY][tilesX][pair][] ([(2*gpu_corr_rad+1)*(2*gpu_corr_rad+1)]) or null
			final double [][][][][]   clt_corr_partial,// [tilesY][tilesX][quad]color][(2*transform_size-1)*(2*transform_size-1)] // if null - will not calculate
			// 											  [tilesY][tilesX] should be set by caller
			final double  [][][]       dcorr_tiles,     // [tile][index][(2*transform_size-1)*(2*transform_size-1)] // if null - will not calculate
			// When clt_mismatch is non-zero, no far objects extraction will be attempted
			final double [][]         clt_mismatch,    // [12][tilesY * tilesX] // ***** transpose unapplied ***** ?. null - do not calculate
			// values in the "main" directions have disparity (*_CM) subtracted, in the perpendicular - as is
			final double [][]         disparity_map,   // [8][tilesY][tilesX], only [6][] is needed on input or null - do not calculate
			// last 2 - contrast, avg/ "geometric average)
			final int                 disparity_modes, // bit mask of disparity_map slices to calculate/return
			final int                 corr_pairs_mask, // which pairs to use in LMA
			final double              max_corr_radius, // 3.9;
			final int                 debug_tileX,
			final int                 debug_tileY,
			final int                 threadsMax,      // maximal number of threads to launch
			final int                 globalDebugLevel)
	{
		final int transform_size = quadCLT.getTileProcessor().getTileSize();
		final boolean is_aux = quadCLT.isAux();
		final boolean is_monochrome = quadCLT.isMonochrome();
		final int [] dbg_tiles = (globalDebugLevel > -1)? new int []{37721, 39022, 40975, 38042, 40963, 42253, 38694, 39343, 36443} : null;
		final boolean [][] saturation_imp = quadCLT.saturation_imp;               // boolean [][] saturation_imp, // (near) saturated pixels or null
		final boolean debug_distort= globalDebugLevel > 0; ///false; // true;

//		final double [][] debug_offsets = new double[imgdtt_params.lma_dbg_offset.length][2];
		final double [][] debug_offsets = new double[quadCLT.getNumSensors()][2];  
		for (int i = 0; i < imgdtt_params.lma_dbg_offset.length; i++) for (int j = 0; j < debug_offsets[i].length; j++) {
			debug_offsets[i][j] = imgdtt_params.lma_dbg_offset[i][j]*imgdtt_params.lma_dbg_scale;
		}
		final int quad = 4;   // number of subcameras
		final int tilesX = quadCLT.getTileProcessor().getTilesX(); // gpuQuad.getTilesX(); // width/transform_size;
		final int tilesY = quadCLT.getTileProcessor().getTilesY(); // gpuQuad.getTilesY(); // final int tilesY=height/transform_size;
		final int tiles = tilesY * tilesX;
		final Thread[] threads = newThreadArray(threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		final double [][] dbg_distort = debug_distort? (new double [4*quad][tilesX*tilesY]) : null;
		final int corr_size = transform_size * 2 - 1;
		// reducing weight of on-axis correlation values to enhance detection of vertical/horizontal lines
		// multiply correlation results inside the horizontal center strip  2*enhortho_width - 1 wide by enhortho_scale
		final double [] enh_ortho_scale = new double [corr_size];
		for (int i = 0; i < corr_size; i++){
			if ((i < (transform_size - imgdtt_params.getEnhOrthoWidth(is_aux))) || (i > (transform_size - 2 + imgdtt_params.getEnhOrthoWidth(is_aux)))) {
				enh_ortho_scale[i] = 1.0;
			} else {
				enh_ortho_scale[i] = imgdtt_params.getEnhOrthoScale(is_aux);
			}
			if (i == (transform_size-1)) enh_ortho_scale[i] = 0.0 ; // hardwired 0 in the center
			enh_ortho_scale[i] *= Math.sin(Math.PI*(i+1.0)/(2*transform_size));
		}
		if (globalDebugLevel > 1){
			System.out.println("getEnhOrthoWidth(isAux())="+ imgdtt_params.getEnhOrthoWidth(is_aux)+" getEnhOrthoScale(isAux())="+ imgdtt_params.getEnhOrthoScale(is_aux));
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

		if (disparity_map != null){
			for (int i = 0; i<disparity_map.length;i++) if ((disparity_modes & (1 << i)) != 0){
				if ((i == OVEREXPOSED) && (saturation_imp == null)) {
					continue;
				}
				disparity_map[i] = new double [tilesY*tilesX];
			}
		}

		if (clt_mismatch != null){
			for (int i = 0; i<clt_mismatch.length;i++){
				clt_mismatch[i] = new double [tilesY*tilesX]; // will use only "center of mass" centers
			}
		}
		final int num_combo = (dcorr_combo == null)? 0 : dcorr_combo.length;
		final int corr_length = (2 * transform_size - 1) * (2 * transform_size - 1);
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				@Override
				public void run() {

					Correlation2d corr2d = new Correlation2d(
							quadCLT.getNumSensors(),
							imgdtt_params,              // ImageDttParameters  imgdtt_params,
							transform_size,             // int transform_size,
							2.0,                        //  double wndx_scale, // (wndy scale is always 1.0)
							is_monochrome,              // boolean monochrome,
							(globalDebugLevel > -1));   //   boolean debug)
					corr2d.createOrtoNotch(
							imgdtt_params.getEnhOrthoWidth(is_aux), // double getEnhOrthoWidth(isAux()),
							imgdtt_params.getEnhOrthoScale(is_aux), //double getEnhOrthoScale(isAux()),
							(imgdtt_params.lma_debug_level > 1)); // boolean debug);
					for (int tIndex = ai.getAndIncrement(); tIndex < tiles; tIndex = ai.getAndIncrement()) {
						int tileY = tIndex / tilesX;
						int tileX = tIndex % tilesX;
						if (dcorr[tileY][tileX] != null) {
							// added quad and cross combos
							double [][]  corrs = new double [GPUTileProcessor.NUM_PAIRS + num_combo][corr_length]; // 225-long (15x15)
							// Prepare the same (currently 10-layer) corrs as double [][], as in CPU version
							int pair_mask = 0;
							if (dcorr != null) {
								for (int pair = 0; pair < dcorr[tileY][tileX].length; pair++) {
									corrs[pair] = 	dcorr[tileY][tileX][pair];
									if (corrs[pair] != null) {
										pair_mask |= (1 << pair);
									}
								}
							}
							// add 4 combo layers : quad, cross, hor, vert

							if (num_combo > 0) {
								for (int pair_combo = 0; pair_combo < dcorr_combo.length; pair_combo++) {
									corrs[GPUTileProcessor.NUM_PAIRS + pair_combo] = dcorr_combo[pair_combo][tileY][tileX];
								}
							}
							if (corr_tiles != null) {
								corr_tiles[tileY][tileX] = corrs; 
							}
							if (dcorr_tiles != null) {
								dcorr_tiles[tIndex] = corrs; // does not require corr_common_GPU()
							}
							if ((disparity_map != null) || (clt_corr_partial != null) || (clt_mismatch != null)) {
								//								int used_pairs = pair_mask; // imgdtt_params.dbg_pair_mask; //TODO: use tile tasks
								int used_pairs = pair_mask & corr_pairs_mask; // imgdtt_params.dbg_pair_mask; // imgdtt_params.dbg_pair_mask; //TODO: use tile tasks
								int tile_lma_debug_level =  ((tileX == debug_tileX) && (tileY == debug_tileY))? (imgdtt_params.lma_debug_level-1) : -2;
								boolean debugTile =(tileX == debug_tileX) && (tileY == debug_tileY) && (globalDebugLevel > -1);
								if (dbg_tiles != null) for (int ii : dbg_tiles){
									if (ii == tIndex) {
										debugTile=true;
										System.out.println("Debugging tile "+tIndex+" ( "+tileX+" / "+tileY+" )");
										break;
									}
								}
								corr_common( // _GPU(
										imgdtt_params,        // final ImageDttParameters  imgdtt_params,
										clt_corr_partial,     // final double [][][][][]   clt_corr_partial,
										used_pairs,           // final int           used_pairs,
										disparity_map,        // final double [][]   disparity_map,
										clt_mismatch,         // final double [][]   clt_mismatch,
										saturation_imp,       // final boolean [][]  saturation_imp,
										false, 				  // final boolean       fneed_macro,
										corr2d,               // final Correlation2d corr2d,
										corrs,                // final double [][]   corrs,
										tilesX,               // final int           tilesX, // new from corr_common_GPU
										tileX,                // final int           tileX,
										tileY,                // final int           tileY,
										max_corr_radius,      // final double        max_corr_radius, // 3.9;
										is_monochrome,        // final boolean       is_monochrome,
										transform_size,       // final int           transform_size,
										tile_lma_debug_level, // int                 tile_lma_debug_level,
										debugTile,            // boolean             debugTile,
										globalDebugLevel);    // final int           globalDebugLevel)							
							}
							if (Double.isNaN(disparity_map[DISPARITY_STRENGTH_INDEX][tIndex])) {
								System.out.println("BUG: 3. disparity_map[DISPARITY_STRENGTH_INDEX][tIndex] should not be NaN");
							}
						}
					} // end of tile
				}
			};
		}
		startAndJoin(threads);
		if ((dbg_distort != null) &&(globalDebugLevel >= 0)) {
			(new ShowDoubleFloatArrays()).showArrays(dbg_distort,  tilesX, tilesY, true, "disparity_distortions"); // , dbg_titles);
		}
	}

	
	//TODO: remove dependence on quad, remove @Deprecated annotation 
	@Deprecated
	public static void corr_common( // remove corr_common_GPU in ImageDtt ?
			
			final ImageDttParameters  imgdtt_params,
			final double [][][][][]   clt_corr_partial,
			final int           used_pairs,
			final double [][]   disparity_map,
			final double [][]   clt_mismatch,
			final boolean [][]  saturation_imp,
			final boolean       fneed_macro,
			final Correlation2d corr2d,
			final double [][]   corrs,
			final int           tilesX, // new from corr_common_GPU
			final int           tileX,
			final int           tileY,
			final double        max_corr_radius, // 3.9; only used with clt_mismatch
			final boolean       is_monochrome,
			final int           transform_size,
			int                 tile_lma_debug_level,
			boolean             debugTile,
			final int           globalDebugLevel)
	{
		final int quad = 4;   // number of subcameras
		final int numcol = is_monochrome?1:3;
		// does not include combo
//		int tile_lma_debug_level =  ((tileX == debug_tileX) && (tileY == debug_tileY))? (imgdtt_params.lma_debug_level-1) : -2;
//		boolean debugTile =(tileX == debug_tileX) && (tileY == debug_tileY) && (globalDebugLevel > -1);
		// non-GPU initialization of the data structures
//		final int tilesX=gpuQuad.getTilesX(); // width/transform_size;
//		final int tilesY=gpuQuad.getTilesY(); // final int tilesY=height/transform_size;
		int tIndex = tileY * tilesX + tileX;
//		final int [] overexp_all = (saturation_imp != null) ? ( new int [2]): null;
		if (disparity_map != null) {
			for (int i = 0; i < disparity_map.length; i++) {
				if (disparity_map[i] != null) {
					if ((((1 << i) & BITS_FROM_GPU) == 0) || !fneed_macro) { // do not touch data already set up
						disparity_map[i][tIndex] = (
								(i == DISPARITY_STRENGTH_INDEX) ||
								(i == DISPARITY_INDEX_HOR_STRENGTH) ||
								(i == DISPARITY_INDEX_VERT_STRENGTH)) ? 0.0 : Double.NaN; // once and for all
					}
				}
			}
			// calculate overexposed fraction
			if (saturation_imp != null){
				// not yet implemented in GPU
				disparity_map[OVEREXPOSED][tIndex] = 0.0; // (1.0 * overexp_all[0]) / overexp_all[1];
			}
			//clt_mismatch should only be used with disparity_map != null;
			if (clt_mismatch != null) {
				for (int np = 0; np < clt_mismatch.length/3; np++) {
					clt_mismatch[3 * np + 0 ][tIndex] = Double.NaN;
					clt_mismatch[3 * np + 1 ][tIndex] = Double.NaN;
					clt_mismatch[3 * np + 2 ][tIndex] = 0;
				}
			}
		}


		// calculate all selected pairs correlations
		//int all_pairs = imgdtt_params.dbg_pair_mask; //TODO: use tile tasks
		// Code that was after correlations calculation

		double [][] strips = corr2d.scaleRotateInterpoateCorrelations(
				corrs,                          // double [][] correlations,
				used_pairs,                      // int         pairs_mask,
				imgdtt_params.corr_strip_hight, //);    // int         hwidth);
				(tile_lma_debug_level > 0) ? used_pairs:0); // debugMax);

		// Combine strips for selected pairs. Now using only for all available pairs.
		// Other combinations are used only if requested (clt_corr_partial != null)

		double [] strip_combo = corr2d.combineInterpolatedCorrelations(
				strips,                        // double [][] strips,
				used_pairs,                     // int         pairs_mask,
				imgdtt_params.corr_offset,     // double      offset);
				imgdtt_params.twice_diagonal); //    		boolean     twice_diagonal)

		double [][] strips_intra = null;
		double [] strip_combo_intra = null;
		if (corrs.length > 6) {
			strips_intra = new double[2][];
			strips_intra[0] = corr2d.scaleRotateInterpoateSingleCorrelation(
					corrs[6], // double []   corr,  // quad 
					imgdtt_params.corr_strip_hight, // int         hwidth,
					Correlation2d.PAIR_HORIZONTAL,  // int         dir, // 0 - hor, 1 - vert, 2 - parallel to row = col (main) diagonal (0->3), 3 -2->1
					1,                              // int         ss, // 1
					false // boolean     debug
					); 
			strips_intra[1] = corr2d.scaleRotateInterpoateSingleCorrelation(
					corrs[7], // double []   corr,  // quad 
					imgdtt_params.corr_strip_hight, // int         hwidth,
					Correlation2d.PAIR_DIAGONAL_MAIN,  // int         dir, // 0 - hor, 1 - vert, 2 - parallel to row = col (main) diagonal (0->3), 3 -2->1
					1,                              // int         ss, // 1
					false // boolean     debug
					); 
			strip_combo_intra = corr2d.combineInterpolatedCorrelations(
					strips_intra,                  // double [][] strips,
					3,                     // int         pairs_mask,
					imgdtt_params.corr_offset,     // double      offset);
					imgdtt_params.twice_diagonal); //    		boolean     twice_diagonal)
		}
		// Debug feature - only calculated if requested
		if ((clt_corr_partial != null) && (imgdtt_params.corr_mode_debug || imgdtt_params.gpu_mode_debug)) {
			@SuppressWarnings("unused")
			double [] strip_ortho = corr2d.combineInterpolatedCorrelations(
					strips,                         // double [][] strips,
					0x0f,                           // int         pairs_mask,
					imgdtt_params.corr_offset,      // double      offset);
					imgdtt_params.twice_diagonal);  //    		boolean     twice_diagonal)
			@SuppressWarnings("unused")
			double [] strip_diag = corr2d.combineInterpolatedCorrelations(
					strips,                         // double [][] strips,
					0x30,                           // int         pairs_mask,
					imgdtt_params.corr_offset,      // double      offset);
					imgdtt_params.twice_diagonal);  //    		boolean     twice_diagonal)
			@SuppressWarnings("unused")
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
			clt_corr_partial[tileY][tileX][1][2] = corrs[6];                        // 7
			clt_corr_partial[tileY][tileX][1][3] = corrs[7];                        // 8
			//												    	clt_corr_partial[tileY][tileX][1][2] = corrs_ortho;                     // 7
			//												    	clt_corr_partial[tileY][tileX][1][3] = corrs_cross;                     // 8
			//												    	clt_corr_partial[tileY][tileX][1][2] = corr2d.debugStrip(strip_hor);    // 7
			//												    	clt_corr_partial[tileY][tileX][1][3] = corr2d.debugStrip(strip_vert);   // 8
			//strip_combo_intra						    	
			clt_corr_partial[tileY][tileX][2][0] = corrs[8];                        // 9
			clt_corr_partial[tileY][tileX][2][1] = corrs[9];                        // 10
			//												    	clt_corr_partial[tileY][tileX][2][0] = corr2d.debugStrip(strips[4]);    // 9
			//												    	clt_corr_partial[tileY][tileX][2][1] = corr2d.debugStrip(strips[5]);    // 10
			clt_corr_partial[tileY][tileX][2][2] = corr2d.debugStrip2(strip_hor);   // 11
			clt_corr_partial[tileY][tileX][2][3] = corr2d.debugStrip2(strip_vert);  // 12
			if (strips_intra != null) {
				clt_corr_partial[tileY][tileX][3][0] = corr2d.debugStrip2(strips_intra[0]); // 13
				clt_corr_partial[tileY][tileX][3][1] = corr2d.debugStrip2(strips_intra[1]);  // 14
			}
			//												    	clt_corr_partial[tileY][tileX][3][0] = corr2d.debugStrip2(strip_ortho); // 13
			//												    	clt_corr_partial[tileY][tileX][3][1] = corr2d.debugStrip2(strip_diag);  // 14
			if (strip_combo_intra != null) {
				clt_corr_partial[tileY][tileX][3][2] = corr2d.debugStrip2(strip_combo_intra);    // 15
			}
			//												    	clt_corr_partial[tileY][tileX][3][2] = corr2d.debugStrip(strip_all);    // 15
			clt_corr_partial[tileY][tileX][3][3] = corr2d.debugStrip2(strip_combo); // 16
		}
		if (imgdtt_params.pcorr_use && (strip_combo_intra != null)) {
			strip_combo = strip_combo_intra;
		}
		// calculate CM maximums for all mixed channels
		// First get integer correlation center, relative to the center
		for (int i = 0; i < strip_combo.length; i++) if (Double.isNaN(strip_combo[i])){
			strip_combo[i] = 0.0; // ????
		}
		int [] ixy =  corr2d.getMaxXYInt( // find integer pair or null if below threshold
				strip_combo,              // double [] data,
				true,                     // boolean   axis_only,
				imgdtt_params.min_corr,   //  double    minMax,    // minimal value to consider (at integer location, not interpolated)
				tile_lma_debug_level > 0); // boolean   debug);

		double [] corr_stat = null;

		// if integer argmax was strong enough, calculate CM argmax
		// will not fill out DISPARITY_INDEX_INT+1, DISPARITY_INDEX_CM+1, DISPARITY_INDEX_POLY+1
		// use clt_mismatch for that
		//							double strength = 0.0;
		//							double disparity = 0.0;
		double [] disp_str = new double[2];
		if (ixy != null) {
			disp_str[1] = strip_combo[ixy[0]+transform_size-1]; // strength at integer max on axis
			if (disparity_map != null) {
				disparity_map[DISPARITY_INDEX_INT][tIndex] =      -ixy[0];
				disparity_map[DISPARITY_STRENGTH_INDEX][tIndex] = disp_str[1];
				if (Double.isNaN(disparity_map[DISPARITY_STRENGTH_INDEX][tIndex])) {
					System.out.println("BUG: 1. disparity_map[DISPARITY_STRENGTH_INDEX]["+tIndex+"] should not be NaN");
				}
			}
			corr_stat = corr2d.getMaxXCm(   // get fractional center as a "center of mass" inside circle/square from the integer max
					strip_combo,                      // double [] data,      // [data_size * data_size]
					ixy[0],                           // int       ixcenter,  // integer center x
					// corr_wndy,                        // double [] window_y,  // (half) window function in y-direction(perpendicular to disparity: for row0  ==1
					// corr_wndx,                        // double [] window_x,  // half of a window function in x (disparity) direction
					(tile_lma_debug_level > 0)); // boolean   debug);
		}
		if (disparity_map != null) {
			if (imgdtt_params.pcorr_use_hv) {
				// for compatibility with old code executed unconditionally. TODO: Move to if (corr_stat != null) ... condition below
				double [] hor_pair1 = corr2d.getMaxXSOrtho(
						corrs,                              // double [][] correlations,
						0x100, //  corrs[8] Correlation2d.getMaskHorizontal(1), // int         pairs_mask,
						imgdtt_params.corr_offset,          // double      corr_offset,
						true,                               // boolean     symmetric,   // for comparing with old implementation average with symmetrical before multiplication
						false,                              // boolean     is_vert,      // transpose X/Y
						debugTile); // tile_lma_debug_level > 0);          // boolean   debug);
				if (hor_pair1 != null){
					disparity_map[DISPARITY_INDEX_HOR][tIndex] =          -hor_pair1[0];
					disparity_map[DISPARITY_INDEX_HOR_STRENGTH][tIndex] =  hor_pair1[1];
				}

				double [] vert_pair1 = corr2d.getMaxXSOrtho(
						corrs,                              // double [][] correlations,
						0x200, // corrs[9] Correlation2d.getMaskVertical(1), // int         pairs_mask,
						imgdtt_params.corr_offset,        // double      corr_offset,
						true,                             // boolean     symmetric,   // for comparing with old implementation average with symmetrical before multiplication
						true, // not anymore transposed false, // already transposed  // true,                             // boolean     is_vert,      // transpose X/Y
						debugTile); // tile_lma_debug_level > 0); // boolean   debug);
				if (vert_pair1 != null) {
					disparity_map[DISPARITY_INDEX_VERT][tIndex] =         -vert_pair1[0];
					disparity_map[DISPARITY_INDEX_VERT_STRENGTH][tIndex] = vert_pair1[1];
				}
			} else  {								
				// for compatibility with old code executed unconditionally. TODO: Move to if (corr_stat != null) ... condition below
				double [] hor_pair1 = corr2d.getMaxXSOrtho(
						corrs,                              // double [][] correlations,
						Correlation2d.getMaskHorizontal(1), // int         pairs_mask,
						imgdtt_params.corr_offset,          // double      corr_offset,
						true,                               // boolean     symmetric,   // for comparing with old implementation average with symmetrical before multiplication
						false,                              // boolean     is_vert,      // transpose X/Y
						debugTile); // tile_lma_debug_level > 0);          // boolean   debug);
				if ((hor_pair1 != null)) {
					disparity_map[DISPARITY_INDEX_HOR][tIndex] =          -hor_pair1[0];
					disparity_map[DISPARITY_INDEX_HOR_STRENGTH][tIndex] =  hor_pair1[1];
				}

				double [] vert_pair1 = corr2d.getMaxXSOrtho(
						corrs,                              // double [][] correlations,
						Correlation2d.getMaskVertical(1), // int         pairs_mask,
						imgdtt_params.corr_offset,        // double      corr_offset,
						true,                             // boolean     symmetric,   // for comparing with old implementation average with symmetrical before multiplication
						true,                             // boolean     is_vert,      // transpose X/Y
						debugTile); // tile_lma_debug_level > 0); // boolean   debug);
				if (vert_pair1 != null) {
					disparity_map[DISPARITY_INDEX_VERT][tIndex] =         -vert_pair1[0];
					disparity_map[DISPARITY_INDEX_VERT_STRENGTH][tIndex] = vert_pair1[1];
				}
			}
		}
		// proceed only if CM correlation result is non-null // for compatibility with old code we need it to run regardless of the strength of the normal correlation
		if (corr_stat != null) {
			// skipping DISPARITY_VARIATIONS_INDEX - it was not used
			disp_str[0] = -corr_stat[0];
			if (disparity_map != null) {
				disparity_map[DISPARITY_INDEX_CM][tIndex] = disp_str[0]; // disparity is negative X
			}
			if (tile_lma_debug_level > 0) {
				System.out.println("Will run getMaxXSOrtho( ) for tileX="+tileX+", tileY="+tileY);
			}

			// debug new LMA correlations
			if (debugTile) {
				System.out.println("Will run new LMA for tileX="+tileX+", tileY="+tileY+"\n\n\n*** Disabled for the GPU ***\n\n\n");
				// temporarily disabling for the GPU (can be restored as disp_dist is available)
				/*
				double [] poly_disp = {Double.NaN, 0.0};
				Corr2dLMA lma2 = corr2d.corrLMA2(
						imgdtt_params,                // ImageDttParameters  imgdtt_params,
						corr_wnd,                     // double [][]         corr_wnd, // correlation window to save on re-calculation of the window
						corr_wnd_inv_limited,         // corr_wnd_limited, // correlation window, limited not to be smaller than threshold - used for finding max/convex areas (or null)
						corrs,                        // double [][]         corrs,
						disp_dist,
						rXY,                          // double [][]         rXY, // non-distorted X,Y offset per nominal pixel of disparity
						imgdtt_params.dbg_pair_mask,  // int                 pair_mask, // which pairs to process
						disp_str,  //corr_stat[0],                 // double    xcenter,   // preliminary center x in pixels for largest baseline
						poly_disp,                    // double[]            poly_ds,    // null or pair of disparity/strength
						imgdtt_params.ortho_vasw_pwr, // double    vasw_pwr,  // value as weight to this power,
						tile_lma_debug_level+2,         // int                 debug_level,
						tileX,                        // int                 tileX, // just for debug output
						tileY );                      // int                 tileY
				double [][] ds = null;
				if (lma2 != null) {
					ds = lma2.lmaDisparityStrength(
							imgdtt_params.lmas_max_rel_rms,  // maximal relative (to average max/min amplitude LMA RMS) // May be up to 0.3)
							imgdtt_params.lmas_min_strength, // minimal composite strength (sqrt(average amp squared over absolute RMS)
							imgdtt_params.lmas_min_ac,       // minimal of A and C coefficients maximum (measures sharpest point/line)
							imgdtt_params.lmas_max_area,      //double  lma_max_area,     // maximal half-area (if > 0.0)
							imgdtt_params.lma_str_scale,    // convert lma-generated strength to match previous ones - scale
							imgdtt_params.lma_str_offset    // convert lma-generated strength to match previous ones - add to result
							);
					lma2.printStats(ds,1);
				}
				 */
			}

			//								disparity_map[DISPARITY_INDEX_CM + 1][tIndex] = // y not available here
			// calculate/fill out hor and vert
			// convert to multi-baseline combining results from several integer scales

			// see if strength is enough to proceed with LMA/poly (otherwise keep disp/strength
			if (disp_str[1] > imgdtt_params.min_poly_strength) {
				// create LMA instance, calculate LMA composite argmax
				// Create 2 groups: ortho & diag
				Correlations2dLMA lma;
				int num_pairs = (quad * (quad-1)) / 2;
				if (imgdtt_params.pcorr_use) { // new group phase correlation
					double [][] fake_corrs = {corrs[6],null,null,null,corrs[7],null};
					lma = corr2d.corrLMA(
							imgdtt_params,                // ImageDttParameters  imgdtt_params,
							fake_corrs,                   // double [][]         corrs,
							corr2d.longToArray(0x11), // imgdtt_params.dbg_pair_mask,  // int                 pair_mask, // which pairs to process
							false,                        // boolean             run_poly_instead, // true - run LMA, false - run 2d polynomial approximation
							corr_stat[0],                 // double    xcenter,   // preliminary center x in pixels for largest baseline
							imgdtt_params.ortho_vasw_pwr, // double    vasw_pwr,  // value as weight to this power,
							tile_lma_debug_level,         // int                 debug_level,
							tileX,                        // int                 tileX, // just for debug output
							tileY );                      // int                 tileY
				} else {
					lma = corr2d.corrLMA(
							imgdtt_params,                // ImageDttParameters  imgdtt_params,
							corrs,                        // double [][]         corrs,
							corr2d.longToArray(used_pairs), // used_pairs, // imgdtt_params.dbg_pair_mask,  // int                 pair_mask, // which pairs to process
							false,                        // boolean             run_poly_instead, // true - run LMA, false - run 2d polynomial approximation
							corr_stat[0],                 // double    xcenter,   // preliminary center x in pixels for largest baseline
							imgdtt_params.ortho_vasw_pwr, // double    vasw_pwr,  // value as weight to this power,
							tile_lma_debug_level,         // int                 debug_level,
							tileX,                        // int                 tileX, // just for debug output
							tileY );                      // int                 tileY
				}
				double [] lma_disparity_strength = null;
				double max_disp_diff_lma = 3.0;
				if (lma != null) {
					lma_disparity_strength = lma.getDisparityStrength();
					if (Math.abs(lma_disparity_strength[0] - disp_str[0] ) > max_disp_diff_lma) {
						if (globalDebugLevel > -1) {
							System.out.println("Crazy LMA for tileX="+tileX+", tileY="+tileY+": disparity="+disp_str[0]+",lma_disparity_strength[0]="+lma_disparity_strength[0]);
						}
						lma = null;
					}
				}
				if (lma != null) {
					double []   mod_disparity_diff = null;
					double [][] dir_corr_strength =  null;
					lma_disparity_strength = lma.getDisparityStrength();
					if (tile_lma_debug_level > 0){
						System.out.println(String.format("Tile X/Y = %d/%d LMA disparity = %7.4f, strength = %7.4f",
								tileX, tileY,
								lma_disparity_strength[0],lma_disparity_strength[1]));
					}
					if (disparity_map != null) {
						disparity_map[DISPARITY_INDEX_POLY]         [tIndex] = lma_disparity_strength[0];
					}

					// if enabled overwrite - replace  DISPARITY_INDEX_CM and DISPARITY_STRENGTH_INDEX
					if (imgdtt_params.mix_corr_poly) { //true
						disp_str[0] =  lma_disparity_strength[0];
						disp_str[1] =  lma_disparity_strength[1];
						if (disparity_map != null) {
							disparity_map[DISPARITY_INDEX_CM]       [tIndex] = disp_str[0];
							disparity_map[DISPARITY_STRENGTH_INDEX] [tIndex] = disp_str[1];
							if (Double.isNaN(disparity_map[DISPARITY_STRENGTH_INDEX][tIndex])) {
								System.out.println("BUG: 2. disparity_map[DISPARITY_STRENGTH_INDEX][tIndex] should not be NaN");
							}
						}
					}
					// store debug data
					// if strong enough and enabled - try to improve far objects
					// temporarily removed strength requirement for debugging, restore later to make faster
					///										if ((imgdtt_params.fo_correct && (strength > imgdtt_params.fo_min_strength)) || (clt_mismatch != null)) {

					if ((imgdtt_params.fo_correct && (disp_str[1] > 0 * imgdtt_params.fo_min_strength)) || (clt_mismatch != null)) {
						// try all dirs:
						dir_corr_strength = corr2d.corr4dirsLMA(
								imgdtt_params,                // ImageDttParameters  imgdtt_params,
								corrs,                        // double [][]         corrs,
								used_pairs, // imgdtt_params.dbg_pair_mask,  // int                 pair_mask, // which pairs to process
								-disp_str[0],                 // double    xcenter,   // preliminary center x in pixels for largest baseline
								imgdtt_params.ortho_vasw_pwr, // double    vasw_pwr,  // value as weight to this power,
								tile_lma_debug_level,         // int                 debug_level,
								tileX,                        // int                 tileX, // just for debug output
								tileY );                      // int                 tileY
						if ((tile_lma_debug_level > 0) && (dir_corr_strength != null)) {
							double [] nan2 = {Double.NaN, Double.NaN, Double.NaN, Double.NaN};
							for (int ii = 0; ii < dir_corr_strength.length; ii++) {
								if (dir_corr_strength[ii] == null) dir_corr_strength[ii] = nan2;
							}
							System.out.println(String.format("corr4dirsLMA -> ↔: %7.4f (%7.4f) ↕:%7.4f (%7.4f) ⇖:%7.4f (%7.4f) ⇗:%7.4f (%7.4f)",
									dir_corr_strength[0][0],dir_corr_strength[0][1],dir_corr_strength[1][0],dir_corr_strength[1][1],
									dir_corr_strength[2][0],dir_corr_strength[2][1],dir_corr_strength[3][0],dir_corr_strength[3][1]));
						}

						mod_disparity_diff =     corr2d.foregroundCorrect(
								imgdtt_params.fo_far,            // boolean   bg,
								imgdtt_params.fo_ortho,          // boolean   ortho,
								dir_corr_strength,               // double [] dir_disp,
								disp_str[0],                     // double    full_disp,
								imgdtt_params.fo_min_strength,   // double      min_strength,
								imgdtt_params.fo_min_eff,        // double      min_eff,
								imgdtt_params.fo_min_eff_ratio,  // double      min_eff_ratio,
								imgdtt_params.fo_max_hwidth,    // double      max_hwidth, //  =          3.0;  // maximal half-width disparity  direction to try to correct far objects
								imgdtt_params.fo_min_diff,       // double    fo_min_diff,
								imgdtt_params.fo_overcorrection, // double    fo_overcorrection,
								imgdtt_params.fo_lim_overcorr,   // double    fo_lim_overcorr,
								(tile_lma_debug_level > 0) );    // boolean debug);

						// Do not use modified far object distance when mismatch is measured
						if ((mod_disparity_diff[0] != disp_str[0]) && (clt_mismatch == null)){ // if it changed
							if (imgdtt_params.fo_correct && (disp_str[1] > imgdtt_params.fo_min_strength)) { // always
								disp_str[0] = mod_disparity_diff[0];
								if (disparity_map != null) {
									disparity_map[DISPARITY_INDEX_CM]       [tIndex] = disp_str[0];
								}
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
							if (imgdtt_params.ly_poly) {  // not used in lwir
								mismatch_result = corr2d.mismatchPairs( // returns x-xcenter, y, strength (sign same as disparity)
										imgdtt_params,                // ImageDttParameters  imgdtt_params,
										corrs,                        // double [][]         corrs,
										used_pairs,                    // int                 pair_mask, // which pairs to process
										-disp_str[0],                   // double    xcenter,   // preliminary center x in pixels for largest baseline
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
										used_pairs,                    // int                 pair_mask, // which pairs to process
										-disp_str[0],                 // double    xcenter,   // preliminary center x in pixels for largest baseline
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
				} else { // if (lma != null)
					if (imgdtt_params.corr_poly_only) { // discard tile if LMA failed
						disp_str[0] =  Double.NaN;
						disp_str[1] =  0.0;
						if (disparity_map != null) {
							disparity_map[DISPARITY_INDEX_CM]       [tIndex] = disp_str[0];
							disparity_map[DISPARITY_STRENGTH_INDEX] [tIndex] = disp_str[1];
						}
					}
				}
			}
		} // end of if (corr_stat != null)
		// double extra_disparity = 0.0; // used for textures:  if allowed, shift images extra before trying to combine
		/*
//Disabled for GPU
			if      (corr_mode == 0) extra_disparity = disparity_map[DISPARITY_INDEX_INT][tIndex];
			else if (corr_mode == 1) extra_disparity = disparity_map[DISPARITY_INDEX_CM][tIndex];
			else if (corr_mode == 2) extra_disparity = disparity_map[DISPARITY_INDEX_POLY][tIndex];
			else if (corr_mode == 3) extra_disparity = disparity_map[DISPARITY_INDEX_HOR][tIndex];  // not used in lwir
			else if (corr_mode == 4) extra_disparity = disparity_map[DISPARITY_INDEX_VERT][tIndex];  // not used in lwir
			if (Double.isNaN(extra_disparity)) extra_disparity = 0;  // used in lwir
		 */
	} 
	
	
	
	

	// extract correlation result  in linescan order (for visualization)
	@Deprecated
	public static double [][] corr_partial_dbg( // not used in lwir
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

		final double [][] corr_data_out = new double[pairs*colors][tilesY*tilesX*tile_size*tile_size];

		final Thread[] threads = newThreadArray(threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		for (int pair = 0; pair< pairs; pair++) {
			for (int nColor = 0; nColor < colors; nColor++) {
				Arrays.fill(corr_data_out[pair*colors+nColor], Double.NaN);
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
//										corr_data_out[indx][((tileY*tile_size + i) *tilesX + tileX)*tile_size+corr_size] = border_contrast*((i & 1) - 0.5);
									}
									for (int i = 0; i < tile_size; i++){
//										corr_data_out[indx][((tileY*tile_size + corr_size) *tilesX + tileX)*tile_size+i] = border_contrast*((i & 1) - 0.5);
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
	
	public static double [][] corr_partial_dbg(
			final double [][][][] corr_data,
			final int             corr_width,
			final int             corr_height,
			final int             threadsMax,     // maximal number of threads to launch
			final int             globalDebugLevel)
	{
		int w=0, h=0;
		for (int i = 0; i < corr_data.length; i++) {
			if (corr_data[i]!=null) {
				h = corr_data[i].length;
				w = corr_data[i][0].length;
				break;
			}
		}
		final int tilesY = h;
		final int tilesX = w;
		final int tile_width =  corr_width + 1;
		final int tile_height = corr_height + 1;

		final double [][] corr_data_out = new double[corr_data.length][]; // tilesY*tilesX*tile_height*tile_width];

		final Thread[] threads = newThreadArray(threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		int num_pairs = 0;
		for (int pair = 0; pair< corr_data.length; pair++) {
			if (corr_data[pair] != null) {
				corr_data_out[pair] = new double [tilesY*tilesX*tile_height*tile_width];
				Arrays.fill(corr_data_out[pair], Double.NaN);
				num_pairs++;
			}
		}
		final int [] pair_indx = new int [num_pairs];
		num_pairs = 0;
		for (int pair = 0; pair< corr_data.length; pair++) {
			if (corr_data[pair] != null) {
				pair_indx[num_pairs++] = pair;
			}
		}
		final int nTiles = tilesX * tilesY;
		final int nAllTiles = nTiles * pair_indx.length;

		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				@Override
				public void run() {
					int tileY,tileX, pair;
					for (int nTile = ai.getAndIncrement(); nTile < nAllTiles; nTile = ai.getAndIncrement()) {
						pair = pair_indx[nTile / nTiles];
						int tile = nTile % nTiles;
						tileY = tile / tilesX;
						tileX = tile % tilesX;
						if (corr_data[pair][tileY][tileX] != null) {
							for (int i = 0; i < corr_height; i++){
								System.arraycopy(
										corr_data[pair][tileY][tileX],
										corr_width * i,
										corr_data_out[pair],
										((tileY*tile_height + i) * tilesX + tileX) * tile_width ,
										corr_width);
							}
						}
					}
				}
			};
		}
		startAndJoin(threads);
		return corr_data_out;
	}
		
	/**
	 * Extract selected window from fcorr [tile][pair][]
	 * @param copy try copy, false - reference
	 * @param fcorr correlation tiles: first index tile number (assuming linescan order), second index - pair index, third - pixel in a tile
	 * @param woi Rectangle x,y,width,height (in tiles)
	 * @param tilesX fcorr width in tiles
	 * @param threadsMax
	 * @return
	 */
	@Deprecated
	public static float [][][] extract_corr_woi(
			final boolean      copy, // copy tiles stack, not reference
			final float [][][] fcorr,
			final Rectangle    woi,
			final int          tilesX,
			final int          threadsMax)     // maximal number of threads to launch

	{
		final int tilesY = fcorr.length/tilesX;
		if ((woi.width + woi.x) >= tilesX) {
			int ww = woi.width; 
			woi.width = tilesX - woi.x;
			if (woi.width <= 0) {
				if (ww > tilesX) ww = tilesX;
				woi.width = ww;
				woi.x = tilesX - woi.width; 
			}
		}
		if ((woi.height + woi.y) >= tilesY) {
			int wndh = woi.height; 
			woi.height = tilesY - woi.y;
			if (woi.height <= 0) {
				if (wndh > tilesY) wndh = tilesY;
				woi.height = wndh;
				woi.y = tilesY - woi.height; 
			}
		}
		final int nTiles=woi.width * woi.height;
		final float [][][] fcorr_out = new float [fcorr.length][][];
		final Thread[] threads = newThreadArray(threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				@Override
				public void run() {
					for (int nTile = ai.getAndIncrement(); nTile < nTiles; nTile = ai.getAndIncrement()) {
						int tileY = nTile / woi.width + woi.y;
						int tileX = nTile % woi.width + woi.x;
						int tile = tileY * tilesX + tileX;
						
						if (copy && (fcorr[tile] != null)) {
							fcorr_out[tile] = fcorr[tile].clone();	// seems a bug, should be fcorr_out[nTile] = fcorr[tile].clone();
						} else {
							fcorr_out[tile] = fcorr[tile];	       // seems a bug, should be fcorr_out[nTile] = fcorr[tile];
						}
					}
				}
			};
		}
		startAndJoin(threads);
		return fcorr_out;
	}

	// 2021
	public static float [][][] extract_corr_woi(
			final boolean      copy, // copy tiles stack, not reference
			final float [][][] fcorr,
			final Rectangle    woi,
			final TpTask []    tp_tasks,        //
			final int          tilesX,
			final int          tilesY,
			final int          threadsMax)     // maximal number of threads to launch

	{
//		int []    tXY = getCorrSize(tp_tasks);
//		final int tilesX = tXY[0];
//		final int tilesY = tXY[1];
		if ((woi.width + woi.x) >= tilesX) {
			int ww = woi.width; 
			woi.width = tilesX - woi.x;
			if (woi.width <= 0) {
				if (ww > tilesX) ww = tilesX;
				woi.width = ww;
				woi.x = tilesX - woi.width; 
			}
		}
		if ((woi.height + woi.y) >= tilesY) {
			int wndh = woi.height; 
			woi.height = tilesY - woi.y;
			if (woi.height <= 0) {
				if (wndh > tilesY) wndh = tilesY;
				woi.height = wndh;
				woi.y = tilesY - woi.height; 
			}
		}
		final float [][][] fcorr_out = new float [fcorr.length][][];
		final Thread[] threads = newThreadArray(threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				@Override
				public void run() {
					for (int nTile = ai.getAndIncrement(); nTile < tp_tasks.length; nTile = ai.getAndIncrement()) {
						int stileY = tp_tasks[nTile].getTileY(); //  nTile/tilesX;
						int stileX = tp_tasks[nTile].getTileX(); //nTile - tileY * tilesX;
						int tileY = stileY - woi.y; //  nTile/tilesX;
						int tileX = stileX - woi.x; //nTile - tileY * tilesX;
						if ((tileY >=0) && (tileX >=0) && (tileY < woi.height) && (tileX < woi.width)) {

//							int tileY = nTile / woi.width + woi.y;
//							int tileX = nTile % woi.width + woi.x;
							int tile = tileY * tilesX + tileX;

							if (copy && (fcorr[tile] != null)) {
								fcorr_out[nTile] = fcorr[tile].clone();	
							} else {
								fcorr_out[nTile] = fcorr[tile];
							}
						}
					}
				}
			};
		}
		startAndJoin(threads);
		return fcorr_out;
	}
	
	
	
	@Deprecated
	public static float [][] corr_partial_dbg( // not used in lwir
			final float  [][][]     fcorr_data,       // [tile][index][(2*transform_size-1)*(2*transform_size-1)] // if null - will not calculate
			final int               tilesX,
			final int               corr_size,
			final int               layers,
			final double            border_contrast,
			final int               threadsMax,     // maximal number of threads to launch
			final int               globalDebugLevel)
	{
		final int tilesY=fcorr_data.length/tilesX;
		final int nTiles=tilesX*tilesY;
		final int tile_size = corr_size+1;
//		final int corr_len = corr_size*corr_size;
		final float [][] fcorr_data_out = new float[layers][tilesY*tilesX*tile_size*tile_size];
		final Thread[] threads = newThreadArray(threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		for (int layer = 0; layer < layers; layer++) {
			Arrays.fill(fcorr_data_out[layer], Float.NaN);
		}
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				@Override
				public void run() {
					int tileY,tileX;
					for (int nTile = ai.getAndIncrement(); nTile < nTiles; nTile = ai.getAndIncrement()) {
						tileY = nTile/tilesX;
						tileX = nTile - tileY * tilesX;
						if (fcorr_data[nTile] != null) {
							for (int layer = 0; layer < layers; layer++) {
								for (int i = 0; i < corr_size;i++){
									System.arraycopy(
											fcorr_data[nTile][layer],
											corr_size* i,
											fcorr_data_out[layer],
											((tileY*tile_size + i) *tilesX + tileX)*tile_size ,
											corr_size);
								}
							}
						}
					}
				}
			};
		}
		startAndJoin(threads);
		return fcorr_data_out;
	}

	// 2021 version, uses tp_tasks as an index to fcorr_data
	public static float [][] corr_partial_dbg( // not used in lwir
			final float  [][][]     fcorr_data,       // [tile][pair][(2*transform_size-1)*(2*transform_size-1)] // if null - will not calculate
			final TpTask []         tp_tasks,        //
			final int               tilesX,
			final int               tilesY,
			final int               corr_size,
			final int               layers0,
			final double            border_contrast,
			final int               threadsMax,     // maximal number of threads to launch
			final int               globalDebugLevel)
	{
		if ((fcorr_data == null) || (fcorr_data.length == 0)) { // assuming each defined tile has the same number of layers and non-null layers
			System.out.println("corr_partial_dbg(): empty fcorr_data");
			return null;
		}
//		int []    tXY = getCorrSize(tp_tasks);
//		final int tilesX = tXY[0];
//		final int tilesY = tXY[1];
		final int tile_size = corr_size+1;
		final int layers = (layers0 < fcorr_data[0].length) ? layers0 : fcorr_data[0].length;
	
		final float [][] fcorr_data_out = new float[layers][]; //  [tilesY*tilesX*tile_size*tile_size];
		final Thread[] threads = newThreadArray(threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		
		for (int layer = 0; layer < layers; layer++) if (fcorr_data[0][layer] != null){
			fcorr_data_out[layer] = new float[tilesY*tilesX*tile_size*tile_size];
			Arrays.fill(fcorr_data_out[layer], Float.NaN);
		}
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				@Override
				public void run() {
					int tileY,tileX;
					for (int nTile = ai.getAndIncrement(); nTile < tp_tasks.length; nTile = ai.getAndIncrement()) {
						tileY = tp_tasks[nTile].getTileY(); //  nTile/tilesX;
						tileX = tp_tasks[nTile].getTileX(); //nTile - tileY * tilesX;
						if (fcorr_data[nTile] != null) {
							for (int layer = 0; layer < layers; layer++) if (fcorr_data[nTile][layer] != null){ // sparse
								for (int i = 0; i < corr_size;i++){
									System.arraycopy(
											fcorr_data[nTile][layer],
											corr_size* i,
											fcorr_data_out[layer],
											((tileY*tile_size + i) *tilesX + tileX)*tile_size ,
											corr_size);
								}
							}
						}
					}
				}
			};
		}
		startAndJoin(threads);
		return fcorr_data_out;
	}

	@Deprecated
	public static float [] corr_partial_wnd( // not used in lwir
			final float  [][][]     fcorr_data,       // [tile][index][(2*transform_size-1)*(2*transform_size-1)] // if null - will not calculate
			final int               tilesX,
			final int               corr_size,
			final Rectangle         woi,
			final int               gap,
			final int []            wh,
			final int               threadsMax)     // maximal number of threads to launch
	{
		final int tile_size = corr_size+1;
		final int tilesY = fcorr_data.length/tilesX;
		final int [][] layout = {{0,0,0},{1,1,0},{2,0,1},{3,1,1},{4,0,2},{5,1,2},{6,0,3},{7,1,3},{8,0,4},{9,1,4}}; // {source_index, row, col};
		if ((woi.width + woi.x) >= tilesX) {
			int ww = woi.width; 
			woi.width = tilesX - woi.x;
			if (woi.width <= 0) {
				if (ww > tilesX) ww = tilesX;
				woi.width = ww;
				woi.x = tilesX - woi.width; 
			}
		}
		if ((woi.height + woi.y) >= tilesY) {
			int wndh = woi.height; 
			woi.height = tilesY - woi.y;
			if (woi.height <= 0) {
				if (wndh > tilesY) wndh = tilesY;
				woi.height = wndh;
				woi.y = tilesY - woi.height; 
			}
		}
		final int nTiles=woi.width * woi.height;
		final int clust_width =  5 * tile_size + gap;
		final int clust_height = 2 * tile_size + gap;
		final int width =  woi.width* clust_width - gap;
		final int height = woi.height*clust_height - gap;
		if (wh != null) {
			wh[0] = width;
			wh[1] = height;
		}
		final float [] corr_data_out = new float[width * height];
		Arrays.fill(corr_data_out, Float.NaN);
		final Thread[] threads = newThreadArray(threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				@Override
				public void run() {
					for (int nTile = ai.getAndIncrement(); nTile < nTiles; nTile = ai.getAndIncrement()) {
						int tileY = nTile / woi.width; // relative to woi
						int tileX = nTile % woi.width;
						int stileY = tileY + woi.y;    // absolute in the corr_data
						int stileX = tileX + woi.x;
						int stile = stileY * tilesX + stileX;
						if (fcorr_data[stile] != null) {
							for (int n = 0; n < layout.length; n++) {
								int src_layer = layout[n][0];
								int v_tile = layout[n][1];
								int h_tile = layout[n][2];
								float [] fcorr_tile = fcorr_data[stile][src_layer]; // tiles were organized as 4x4
								int out_x = tileX * clust_width +  h_tile * tile_size;
								int out_y = tileY * clust_height + v_tile * tile_size;
								for (int i = 0; i < corr_size;i++){
									System.arraycopy(
											fcorr_tile,
											corr_size* i,
											corr_data_out,
											(out_y + i) * width + out_x,
											corr_size);
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
	
	@Deprecated
	public static int [] getCorrSize( // not used in lwir
			final TpTask []         tp_tasks) {
		int ty = -1,tx = -1;
		for (int i = 0; i < tp_tasks.length; i++) {
			if (tp_tasks[i].getTileX() > tx) tx = tp_tasks[i].getTileX(); 
			if (tp_tasks[i].getTileY() > ty) ty = tp_tasks[i].getTileY(); 
		}
		return new int [] {tx,ty};
	}

	public static int getCorrPairs( // not used in lwir
			final float  [][][]     fcorr_data) {       // [tile][pair][(2*transform_size-1)*(2*transform_size-1)] // if null - will not calculate
		for (int i = 0; i < fcorr_data.length; i++) {
			if (fcorr_data[i] != null) {
				return fcorr_data[i].length;
			}
		}
		return 0;
	}
	
	
	
	
	public static int [][] getCorrLayout(   // {source_index, row, col};
			Object corr_pairs,
			double max_aspect){ // 3.0
		if (Double.isNaN(max_aspect)) {
			max_aspect = 3.0;
		}
		boolean [] sel_pairs = null;
		int num_pairs = -1;
		if (corr_pairs instanceof double[][][]) {
			double [][][] parray = (double [][][]) corr_pairs;
			for (int nTile = 0; nTile < parray.length; nTile++) {
				if (parray[nTile] != null) {
					num_pairs = parray[nTile].length;
					sel_pairs = new boolean [num_pairs];
					for (int npair = 0; npair < num_pairs; npair++) {
						sel_pairs[npair] = parray[nTile][npair] != null;
					}
					break;
				}
			}
		} else if (corr_pairs instanceof float [][][]) {
			float [][][] parray = (float [][][]) corr_pairs;
			for (int nTile = 0; nTile < parray.length; nTile++) {
				if (parray[nTile] != null) {
					num_pairs = parray[nTile].length;
					sel_pairs = new boolean [num_pairs];
					for (int npair = 0; npair < num_pairs; npair++) {
						sel_pairs[npair] = parray[nTile][npair] != null;
					}
					break;
				}
			}
		} else {
			  throw new IllegalArgumentException ("getCorrLayout(): only double [][][] and float [][][] are supported");
		}
		if (sel_pairs == null) {
			System.out.println("getCorrLayout(): empty");
			return null;
		}
		int num_used = 0;
		for (int i = 0; i < sel_pairs.length; i++) {
			if (sel_pairs[i]) num_used++;
		}
		if (num_used == 0) {
			System.out.println("getCorrLayout(): Number of define correlation pairs = 0");
			return null;
		}
		int max_height = (int) Math.sqrt(num_used);
		int waste = num_used, best_width = -1;
		
		for (int height = max_height; (height > 0) && (Math.ceil(1.0*num_used/height) <  height *max_aspect ); height--) {
			int width = num_used/ height;
			if (width * height < num_used) width ++;
			int new_waste = width * height - num_used;
			if (new_waste < waste) {
				waste = new_waste;
				best_width = width;
			}
		}
		int [][] rslt = new int [num_used][3];
		int indx = 0, row = 0, col = 0;
		for (int i = 0; i < num_pairs; i++) {
			if (sel_pairs[i]) {
				rslt[indx  ][0] = i;
				rslt[indx  ][1] = row;
				rslt[indx++][2] = col++;
				if (col >= best_width) {
					row++;
					col = 0;
				}
			}
		}
		return rslt;     
	}
	
	// 2021 version, uses tp_tasks as an index to fcorr_data
	public static float [] corr_partial_wnd( // not used in lwir
			final float  [][][]     fcorr_data,       // [tile][index][(2*transform_size-1)*(2*transform_size-1)] // if null - will not calculate
			final TpTask []         tp_tasks,        //
			final int               tilesX,
			final int               tilesY,
			final int               corr_size,
			final Rectangle         woi,
			final int               gap,
			final int []            wh,
			final int               threadsMax)     // maximal number of threads to launch
	{
		if ((fcorr_data == null) || (fcorr_data.length == 0)) { // assuming each defined tile has the same number of layers and non-null layers
			System.out.println("corr_partial_dbg(): empty fcorr_data");
			return null;
		}
//		int []    tXY = getCorrSize(tp_tasks);
//		final int tilesX = tXY[0];
//		final int tilesY = tXY[1];
		
		final int tile_size = corr_size+1;
//		final int tilesY = fcorr_data.length/tilesX;
//		final int [][] layout = {{0,0,0},{1,1,0},{2,0,1},{3,1,1},{4,0,2},{5,1,2},{6,0,3},{7,1,3},{8,0,4},{9,1,4}}; // {source_index, row, col};
		final int [][] layout = getCorrLayout(  // {source_index, row, col};
				fcorr_data, // Object corr_pairs,
				3.0); // double max_aspect);
		
		if ((woi.width + woi.x) >= tilesX) {
			int ww = woi.width; 
			woi.width = tilesX - woi.x;
			if (woi.width <= 0) {
				if (ww > tilesX) ww = tilesX;
				woi.width = ww;
				woi.x = tilesX - woi.width; 
			}
		}
		if ((woi.height + woi.y) >= tilesY) {
			int wndh = woi.height; 
			woi.height = tilesY - woi.y;
			if (woi.height <= 0) {
				if (wndh > tilesY) wndh = tilesY;
				woi.height = wndh;
				woi.y = tilesY - woi.height; 
			}
		}
//		final int nTiles=woi.width * woi.height;
		final int clust_width =  5 * tile_size + gap;
		final int clust_height = 2 * tile_size + gap;
		final int width =  woi.width* clust_width - gap;
		final int height = woi.height*clust_height - gap;
		if (wh != null) {
			wh[0] = width;
			wh[1] = height;
		}
		final float [] corr_data_out = new float[width * height];
		Arrays.fill(corr_data_out, Float.NaN);
		final Thread[] threads = newThreadArray(threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				@Override
				public void run() {
					for (int nTile = ai.getAndIncrement(); nTile < tp_tasks.length; nTile = ai.getAndIncrement()) {
						int stileY = tp_tasks[nTile].getTileY(); //  nTile/tilesX;
						int stileX = tp_tasks[nTile].getTileX(); //nTile - tileY * tilesX;
						int tileY = stileY - woi.y; //  nTile/tilesX;
						int tileX = stileX - woi.x; //nTile - tileY * tilesX;
						if ((tileY >=0) && (tileX >=0) && (tileY < woi.height) && (tileX < woi.width)) {
//							int tileY = nTile / woi.width; // relative to woi
//							int tileX = nTile % woi.width;
//							int stileY = tileY + woi.y;    // absolute in the corr_data
//							int stileX = tileX + woi.x;
//							int stile = stileY * tilesX + stileX;
							if (fcorr_data[nTile] != null) {
								for (int n = 0; n < layout.length; n++) {
									int src_layer = layout[n][0];
									int v_tile = layout[n][1];
									int h_tile = layout[n][2];
									float [] fcorr_tile = fcorr_data[nTile][src_layer]; // tiles were organized as 4x4 ??
									int out_x = tileX * clust_width +  h_tile * tile_size;
									int out_y = tileY * clust_height + v_tile * tile_size;
									for (int i = 0; i < corr_size;i++){
										int out_start = (out_y + i) * width + out_x;

										System.arraycopy(
												fcorr_tile,
												corr_size* i,
												corr_data_out,
												(out_y + i) * width + out_x,
												corr_size);
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

	@Deprecated
	public static double [] corr_partial_wnd( // not used in lwir
			final double  [][][]    dcorr_data,       // [tile][index][(2*transform_size-1)*(2*transform_size-1)] // if null - will not calculate
			final int               tilesX,
			final int               corr_size,
			final Rectangle         woi,
			final int               gap,
			final int []            wh,
			final int               threadsMax)     // maximal number of threads to launch
	{
		final int tile_size = corr_size+1;
		final int tilesY = dcorr_data.length/tilesX;
		final int [][] layout = {{0,0,0},{1,1,0},{2,0,1},{3,1,1},{4,0,2},{5,1,2},{6,0,3},{7,1,3},{8,0,4},{9,1,4}}; // {source_index, row, col};
		if ((woi.width + woi.x) >= tilesX) {
			int ww = woi.width; 
			woi.width = tilesX - woi.x;
			if (woi.width <= 0) {
				if (ww > tilesX) ww = tilesX;
				woi.width = ww;
				woi.x = tilesX - woi.width; 
			}
		}
		if ((woi.height + woi.y) >= tilesY) {
			int wndh = woi.height; 
			woi.height = tilesY - woi.y;
			if (woi.height <= 0) {
				if (wndh > tilesY) wndh = tilesY;
				woi.height = wndh;
				woi.y = tilesY - woi.height; 
			}
		}
		final int nTiles=woi.width * woi.height;
		final int clust_width =  5 * tile_size + gap;
		final int clust_height = 2 * tile_size + gap;
		final int width =  woi.width* clust_width - gap;
		final int height = woi.height*clust_height - gap;
		if (wh != null) {
			wh[0] = width;
			wh[1] = height;
		}
		final double [] corr_data_out = new double[width * height];
		Arrays.fill(corr_data_out, Float.NaN);
		final Thread[] threads = newThreadArray(threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				@Override
				public void run() {
					for (int nTile = ai.getAndIncrement(); nTile < nTiles; nTile = ai.getAndIncrement()) {
						int tileY = nTile / woi.width; // relative to woi
						int tileX = nTile % woi.width;
						int stileY = tileY + woi.y;    // absolute in the corr_data
						int stileX = tileX + woi.x;
						int stile = stileY * tilesX + stileX;
						if (dcorr_data[stile] != null) {
							for (int n = 0; n < layout.length; n++) {
								int src_layer = layout[n][0];
								if (dcorr_data[stile].length > src_layer) {
									int v_tile = layout[n][1];
									int h_tile = layout[n][2];
									double [] dcorr_tile = dcorr_data[stile][src_layer]; // tiles were organized as 4x4
									int out_x = tileX * clust_width +  h_tile * tile_size;
									int out_y = tileY * clust_height + v_tile * tile_size;
									for (int i = 0; i < corr_size;i++){
										System.arraycopy(
												dcorr_tile,
												corr_size* i,
												corr_data_out,
												(out_y + i) * width + out_x,
												corr_size);
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
	
	// 2021 
	public static double [] corr_partial_wnd( // not used in lwir
			final double  [][][]    dcorr_data,      // [tile][index][(2*transform_size-1)*(2*transform_size-1)] // if null - will not calculate
			final TpTask []         tp_tasks,        //
			final int               tilesX,
			final int               tilesY,
			final int               corr_size,
			final Rectangle         woi,
			final int               gap,
			final int []            wh,
			final int               threadsMax)     // maximal number of threads to launch
	{
//		int []    tXY = getCorrSize(tp_tasks);
//		final int tilesX = tXY[0];
//		final int tilesY = tXY[1];
		final int tile_size = corr_size+1;

//		final int [][] layout = {{0,0,0},{1,1,0},{2,0,1},{3,1,1},{4,0,2},{5,1,2},{6,0,3},{7,1,3},{8,0,4},{9,1,4}}; // {source_index, row, col};
		final int [][] layout = getCorrLayout(  // {source_index, row, col};
				dcorr_data, // Object corr_pairs,
				3.0); // double max_aspect);
		
		if ((woi.width + woi.x) >= tilesX) {
			int ww = woi.width; 
			woi.width = tilesX - woi.x;
			if (woi.width <= 0) {
				if (ww > tilesX) ww = tilesX;
				woi.width = ww;
				woi.x = tilesX - woi.width; 
			}
		}
		if ((woi.height + woi.y) >= tilesY) {
			int wndh = woi.height; 
			woi.height = tilesY - woi.y;
			if (woi.height <= 0) {
				if (wndh > tilesY) wndh = tilesY;
				woi.height = wndh;
				woi.y = tilesY - woi.height; 
			}
		}
//		final int nTiles=woi.width * woi.height;
		final int clust_width =  5 * tile_size + gap;
		final int clust_height = 2 * tile_size + gap;
		final int width =  woi.width* clust_width - gap;
		final int height = woi.height*clust_height - gap;
		if (wh != null) {
			wh[0] = width;
			wh[1] = height;
		}
		final double [] corr_data_out = new double[width * height];
		Arrays.fill(corr_data_out, Float.NaN);
		final Thread[] threads = newThreadArray(threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				@Override
				public void run() {
					for (int nTile = ai.getAndIncrement(); nTile < tp_tasks.length; nTile = ai.getAndIncrement()) {
						int stileY = tp_tasks[nTile].getTileY(); //  nTile/tilesX;
						int stileX = tp_tasks[nTile].getTileX(); //nTile - tileY * tilesX;
						int tileY = stileY - woi.y; //  nTile/tilesX;
						int tileX = stileX - woi.x; //nTile - tileY * tilesX;
						if ((tileY >=0) && (tileX >=0) && (tileY < woi.height) && (tileX < woi.width)) {
//							int tileY = nTile / woi.width; // relative to woi
//							int tileX = nTile % woi.width;
//							int stileY = tileY + woi.y;    // absolute in the corr_data
//							int stileX = tileX + woi.x;
//							int stile = stileY * tilesX + stileX;
							if (dcorr_data[nTile] != null) {
								for (int n = 0; n < layout.length; n++) {
									int src_layer = layout[n][0];
									if (dcorr_data[nTile].length > src_layer) {
										int v_tile = layout[n][1];
										int h_tile = layout[n][2];
										double [] dcorr_tile = dcorr_data[nTile][src_layer]; // tiles were organized as 4x4
										int out_x = tileX * clust_width +  h_tile * tile_size;
										int out_y = tileY * clust_height + v_tile * tile_size;
										for (int i = 0; i < corr_size;i++){
											System.arraycopy(
													dcorr_tile,
													corr_size* i,
													corr_data_out,
													(out_y + i) * width + out_x,
													corr_size);
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
		return corr_data_out;
	}
	
	@Deprecated
	public static float [] corr_td_wnd(
			final float [][][][] fcorr_td,       // float[tilesY][tilesX][num_slices][];
			final float [][][][] fcorr_combo_td, // float[4][tilesY][tilesX][];
			final Rectangle      woi,
			final int            gap,
			final int []         wh,
			final int            transform_size,
			final int            threadsMax)     // maximal number of threads to launch
	{
		final int tile_size = 2 * transform_size;
		final int [][] layout1 = {{0,0,0},{1,1,0},{2,0,1},{3,1,1},{4,0,2},{5,1,2}}; // {source_index, row, col};
		final int [][] layout2 = {{0,0,3},{1,1,3},{2,0,4},{3,1,4}}; // {source_index, row, col};
		if ((woi.width + woi.x) >= fcorr_td[0].length) {
			int ww = woi.width; 
			woi.width = fcorr_td[0].length - woi.x;
			if (woi.width <= 0) {
				if (ww > fcorr_td[0].length) ww = fcorr_td[0].length;
				woi.width = ww;
				woi.x = fcorr_td[0].length - woi.width; 
			}
		}
		if ((woi.height + woi.y) >= fcorr_td.length) {
			int wndh = woi.height; 
			woi.height = fcorr_td.length - woi.y;
			if (woi.height <= 0) {
				if (wndh > fcorr_td.length) wndh = fcorr_td.length;
				woi.height = wndh;
				woi.y = fcorr_td.length - woi.height; 
			}
		}
		
		final int nTiles=woi.width * woi.height;
		
		final int clust_width =  5 * tile_size + gap;
		final int clust_height = 2 * tile_size + gap;
		
		final int width =  woi.width* clust_width - gap;
		final int height = woi.height*clust_height - gap;
		if (wh != null) {
			wh[0] = width;
			wh[1] = height;
		}
		final int transform_len = transform_size * transform_size;
		final float [] fcorr_data_out = new float[width * height];
		Arrays.fill(fcorr_data_out, Float.NaN);
		
		final Thread[] threads = newThreadArray(threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				@Override
				public void run() {
					for (int nTile = ai.getAndIncrement(); nTile < nTiles; nTile = ai.getAndIncrement()) {
						int tileY = nTile / woi.width; // relative to woi
						int tileX = nTile % woi.width;
						int stileY = tileY + woi.y;    // absolute in the corr_data
						int stileX = tileX + woi.x;
						// first 6 tiles from fcorr_td =       new float[tilesY][tilesX][num_slices][];
						if (fcorr_td[stileY][stileX] != null) {
							for (int n = 0; n < layout1.length; n++) {
								int src_layer = layout1[n][0];
								int v_tile = layout1[n][1];
								int h_tile = layout1[n][2];
								float [] fcorr_tile = fcorr_td[stileY][stileX][src_layer];
								for (int qy = 0; qy < 2; qy++) {
									for (int qx = 0; qx < 2; qx++) {
										for (int i = 0; i < transform_size;i++){
											System.arraycopy(
													fcorr_tile,
													transform_len * (2 *qy + qx) + transform_size * i,
													fcorr_data_out,
													(tileY * clust_height + v_tile * tile_size + qy * transform_size + i) * width +
													(tileX * clust_width +  h_tile * tile_size + qx * transform_size),
													transform_size);
										}
									}
								}
							}
						}
						// last 4 tiles from fcorr_combo_td = new float[4][tilesY][tilesX][];
						if (fcorr_combo_td != null) {
							for (int n = 0; n < layout2.length; n++) {
								if (fcorr_combo_td[n][stileY][stileX] != null) {
									int src_layer = layout2[n][0];
									int v_tile = layout2[n][1];
									int h_tile = layout2[n][2];
									float [] fcorr_tile = fcorr_combo_td[src_layer][stileY][stileX];
									for (int qy = 0; qy < 2; qy++) {
										for (int qx = 0; qx < 2; qx++) {
											for (int i = 0; i < transform_size;i++){
												System.arraycopy(
														fcorr_tile,
														transform_len * (2 *qy + qx) + transform_size * i,
														fcorr_data_out,
														(tileY * clust_height + v_tile * tile_size + qy * transform_size + i) * width +
														(tileX * clust_width +  h_tile * tile_size + qx * transform_size),
														transform_size);
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
		return fcorr_data_out;
	}

	@Deprecated
	public static double [] corr_td_wnd( // not used
			final double [][][][] dcorr_td,       // float[tilesY][tilesX][num_slices][];
			final double [][][][] dcorr_combo_td, // float[4][tilesY][tilesX][];
			final Rectangle      woi,
			final int            gap,
			final int []         wh,
			final int            transform_size,
			final int            threadsMax)     // maximal number of threads to launch
	{
		final int tile_size = 2 * transform_size;
		final int [][] layout1 = {{0,0,0},{1,1,0},{2,0,1},{3,1,1},{4,0,2},{5,1,2}}; // {source_index, row, col};
		final int [][] layout2 = {{0,0,3},{1,1,3},{2,0,4},{3,1,4}}; // {source_index, row, col};
		if ((woi.width + woi.x) >= dcorr_td[0].length) {
			int ww = woi.width; 
			woi.width = dcorr_td[0].length - woi.x;
			if (woi.width <= 0) {
				if (ww > dcorr_td[0].length) ww = dcorr_td[0].length;
				woi.width = ww;
				woi.x = dcorr_td[0].length - woi.width; 
			}
		}
		if ((woi.height + woi.y) >= dcorr_td.length) {
			int wndh = woi.height; 
			woi.height = dcorr_td.length - woi.y;
			if (woi.height <= 0) {
				if (wndh > dcorr_td.length) wndh = dcorr_td.length;
				woi.height = wndh;
				woi.y = dcorr_td.length - woi.height; 
			}
		}
		
		final int nTiles=woi.width * woi.height;
		
		final int clust_width =  5 * tile_size + gap;
		final int clust_height = 2 * tile_size + gap;
		
		final int width =  woi.width* clust_width - gap;
		final int height = woi.height*clust_height - gap;
		if (wh != null) {
			wh[0] = width;
			wh[1] = height;
		}
		final int transform_len = transform_size * transform_size;
		final double [] dcorr_data_out = new double[width * height];
		Arrays.fill(dcorr_data_out, Float.NaN);
		
		final Thread[] threads = newThreadArray(threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				@Override
				public void run() {
					for (int nTile = ai.getAndIncrement(); nTile < nTiles; nTile = ai.getAndIncrement()) {
						int tileY = nTile / woi.width; // relative to woi
						int tileX = nTile % woi.width;
						int stileY = tileY + woi.y;    // absolute in the corr_data
						int stileX = tileX + woi.x;
						// first 6 tiles from fcorr_td =       new float[tilesY][tilesX][num_slices][];
						if (dcorr_td[stileY][stileX] != null) {
							for (int n = 0; n < layout1.length; n++) {
								int src_layer = layout1[n][0];
								int v_tile = layout1[n][1];
								int h_tile = layout1[n][2];
								double [] dcorr_tile = dcorr_td[stileY][stileX][src_layer];
								for (int qy = 0; qy < 2; qy++) {
									for (int qx = 0; qx < 2; qx++) {
										for (int i = 0; i < transform_size;i++){
											System.arraycopy(
													dcorr_tile,
													transform_len * (2 *qy + qx) + transform_size * i,
													dcorr_data_out,
													(tileY * clust_height + v_tile * tile_size + qy * transform_size + i) * width +
													(tileX * clust_width +  h_tile * tile_size + qx * transform_size),
													transform_size);
										}
									}
								}
							}
						}
						// last 4 tiles from fcorr_combo_td = new float[4][tilesY][tilesX][];
						if (dcorr_combo_td != null) {
							for (int n = 0; n < layout2.length; n++) {
								int src_layer = layout2[n][0];
								int v_tile = layout2[n][1];
								int h_tile = layout2[n][2];
								double [] dcorr_tile = dcorr_combo_td[src_layer][stileY][stileX];
								for (int qy = 0; qy < 2; qy++) {
									for (int qx = 0; qx < 2; qx++) {
										for (int i = 0; i < transform_size;i++){
											System.arraycopy(
													dcorr_tile,
													transform_len * (2 *qy + qx) + transform_size * i,
													dcorr_data_out,
													(tileY * clust_height + v_tile * tile_size + qy * transform_size + i) * width +
													(tileX * clust_width +  h_tile * tile_size + qx * transform_size),
													transform_size);
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
		return dcorr_data_out;
	}

	
	@Deprecated
	public static double [] corr_td_wnd( // only used in debug version
			final double [][][][][] dcorr_td,       // float[tilesY][tilesX][num_slices][];
			final double [][][][][] dcorr_combo_td, // float[4][tilesY][tilesX][];
			final Rectangle         woi,
			final int               gap,
			final int []            wh,
			final int               transform_size,
			final int               threadsMax)     // maximal number of threads to launch
	{
		final int tile_size = 2 * transform_size;
		final int [][] layout1 = {{0,0,0},{1,1,0},{2,0,1},{3,1,1},{4,0,2},{5,1,2}}; // {source_index, row, col};
		final int [][] layout2 = {{0,0,3},{1,1,3},{2,0,4},{3,1,4}}; // {source_index, row, col};
		if ((woi.width + woi.x) >= dcorr_td[0].length) {
			int ww = woi.width; 
			woi.width = dcorr_td[0].length - woi.x;
			if (woi.width <= 0) {
				if (ww > dcorr_td[0].length) ww = dcorr_td[0].length;
				woi.width = ww;
				woi.x = dcorr_td[0].length - woi.width; 
			}
		}
		if ((woi.height + woi.y) >= dcorr_td.length) {
			int wndh = woi.height; 
			woi.height = dcorr_td.length - woi.y;
			if (woi.height <= 0) {
				if (wndh > dcorr_td.length) wndh = dcorr_td.length;
				woi.height = wndh;
				woi.y = dcorr_td.length - woi.height; 
			}
		}
		
		final int nTiles=woi.width * woi.height;
		
		final int clust_width =  5 * tile_size + gap;
		final int clust_height = 2 * tile_size + gap;
		
		final int width =  woi.width* clust_width - gap;
		final int height = woi.height*clust_height - gap;
		if (wh != null) {
			wh[0] = width;
			wh[1] = height;
		}
		final int transform_len = transform_size * transform_size;
		final double [] dcorr_data_out = new double[width * height];
		Arrays.fill(dcorr_data_out, Float.NaN);
		
		final Thread[] threads = newThreadArray(threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				@Override
				public void run() {
					for (int nTile = ai.getAndIncrement(); nTile < nTiles; nTile = ai.getAndIncrement()) {
						int tileY = nTile / woi.width; // relative to woi
						int tileX = nTile % woi.width;
						int stileY = tileY + woi.y;    // absolute in the corr_data
						int stileX = tileX + woi.x;
						// first 6 tiles from fcorr_td =       new float[tilesY][tilesX][num_slices][];
						if (dcorr_td[stileY][stileX] != null) {
							for (int n = 0; n < layout1.length; n++) {
								int src_layer = layout1[n][0];
								int v_tile = layout1[n][1];
								int h_tile = layout1[n][2];
								double [][] dcorr_tile = dcorr_td[stileY][stileX][src_layer];
								for (int qy = 0; qy < 2; qy++) {
									for (int qx = 0; qx < 2; qx++) {
										for (int i = 0; i < transform_size;i++){
											System.arraycopy(
													dcorr_tile[2*qy+qx],
													transform_size * i,
													dcorr_data_out,
													(tileY * clust_height + v_tile * tile_size + qy * transform_size + i) * width +
													(tileX * clust_width +  h_tile * tile_size + qx * transform_size),
													transform_size);
										}
									}
								}
							}
						}
						// last 4 tiles from fcorr_combo_td = new float[4][tilesY][tilesX][];
						if (dcorr_combo_td != null) {
							for (int n = 0; n < layout2.length; n++) {
								int src_layer = layout2[n][0];
								int v_tile = layout2[n][1];
								int h_tile = layout2[n][2];
								double [][] dcorr_tile = dcorr_combo_td[src_layer][stileY][stileX];
								for (int qy = 0; qy < 2; qy++) {
									for (int qx = 0; qx < 2; qx++) {
										for (int i = 0; i < transform_size;i++){
											System.arraycopy(
													dcorr_tile[2*qy+qx],
													transform_size * i,
													dcorr_data_out,
													(tileY * clust_height + v_tile * tile_size + qy * transform_size + i) * width +
													(tileX * clust_width +  h_tile * tile_size + qx * transform_size),
													transform_size);
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
		return dcorr_data_out;
	}
	
	
		
	public static float [][][] corr_get_extra(
			final float [][][][] fcorrs,
			final int            tilesX,
			final int            ncombo,
			final int            slices,
			final int            threadsMax)     // maximal number of threads to launch
	{
		final int tiles = fcorrs[ncombo].length;
		final int ncorrs = fcorrs.length - 1;
		float [][][] fcorr_extra = new float [tiles][][];
		
		final Thread[] threads = newThreadArray(threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				@Override
				public void run() {
					for (int nTile = ai.getAndIncrement(); nTile < tiles; nTile = ai.getAndIncrement()) if (fcorrs[ncombo][nTile] != null){
						fcorr_extra[nTile] = new float [slices][ncorrs];
						for (int indx0 = 0; indx0 < ncorrs; indx0++) {
							int indx = indx0 + ((indx0 < ncombo) ? 0 : 1);
							if (fcorrs[indx][nTile] != null) {
								for (int slice = 0; slice < slices; slice++) {
									if ((fcorrs[ncombo][nTile][slice] != null) && (fcorrs[indx][nTile][slice] != null)) {
										float [] t0 = fcorrs[ncombo][nTile][slice];
										float [] t1 = fcorrs[indx][nTile][slice];
										float s00 = 0.0f, s11 = 0.0f, s01 = 0.0f;
										for (int i = 0; i < t0.length; i++) {
											s00 += t0[i] * t0[i];
											s11 += t1[i] * t1[i];
											s01 += t0[i] * t1[i];
										}
										fcorr_extra[nTile][slice][indx] = (float) (s01/Math.sqrt(s00*s11));
									}
								}
							}
						}
					}
				}
			};
		}
		startAndJoin(threads);
		return fcorr_extra;
	}
	
	// prepare tile-to-tile correlation data as an extra layer, after layers by corr_td_wnd()
	@Deprecated
	public static float [] corr_show_extra(
			final float [][][] fcorr_extra,  // float[tile][slices][extra];
			final int          tilesX,
			final Rectangle    woi,
			final int          gap,
			final int          step, // 3
			final int          size, //2
			final int []       wh,
			final int          transform_size,
			final int          threadsMax)     // maximal number of threads to launch
	{
		final int tilesY = fcorr_extra.length / tilesX;
		final int tile_size = 2 * transform_size;
		final int per_row = (tile_size + (step-size)) /step;
		final int [][] layout = {{0,0,0},{1,1,0},{2,0,1},{3,1,1},{4,0,2},{5,1,2},{6,0,3},{7,1,3},{8,0,4},{9,1,4}}; // {source_index, row, col};
		if ((woi.width + woi.x) >= tilesX) {
			int ww = woi.width; 
			woi.width = tilesX - woi.x;
			if (woi.width <= 0) {
				if (ww > tilesX) ww = tilesX;
				woi.width = ww;
				woi.x = tilesX - woi.width; 
			}
		}
		if ((woi.height + woi.y) >= tilesY) {
			int wndh = woi.height; 
			woi.height = tilesY - woi.y;
			if (woi.height <= 0) {
				if (wndh > tilesY) wndh = fcorr_extra.length;
				woi.height = wndh;
				woi.y = tilesY - woi.height; 
			}
		}
		final int nTiles=woi.width * woi.height;
		final int clust_width =  5 * tile_size + gap;
		final int clust_height = 2 * tile_size + gap;
		final int width =  woi.width* clust_width - gap;
		final int height = woi.height*clust_height - gap;
		if (wh != null) {
			wh[0] = width;
			wh[1] = height;
		}
		final float [] fcorr_data_out = new float[width * height];
		Arrays.fill(fcorr_data_out, Float.NaN);
		final Thread[] threads = newThreadArray(threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				@Override
				public void run() {
					for (int nTile = ai.getAndIncrement(); nTile < nTiles; nTile = ai.getAndIncrement()) {
						int tileY = nTile / woi.width; // relative to woi
						int tileX = nTile % woi.width;
						int stileY = tileY + woi.y;    // absolute in the corr_data
						int stileX = tileX + woi.x;
						int stile = stileY * tilesX + stileX; 
						// first 6 tiles from fcorr_td =       new float[tilesY][tilesX][num_slices][];
						if (fcorr_extra[stile] != null) {
							for (int n = 0; n < layout.length; n++) {
								int src_layer = layout[n][0];
								int v_tile = layout[n][1];
								int h_tile = layout[n][2];
								float [] extra_data = fcorr_extra[stile][src_layer];
								for (int i = 0; i < extra_data.length; i++) {
									int extra_row = i / per_row;
									int extra_col = i % per_row;
									int extra_0 =
											(tileY * clust_height + v_tile * tile_size + extra_row * step) * width +
											(tileX * clust_width +  h_tile * tile_size + extra_col * step);
									for (int sy = 0; sy < size; sy++) {
										for (int sx = 0; sx < size; sx++) {
											fcorr_data_out[extra_0 + sy* width + sx] = extra_data[i];
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
		return fcorr_data_out;
	}

	// calculate inter-tile correlation from the data already converted to the debug images
	// (and so only for the selected woi)




	public double [][][][][] cltStack( // not used in lwir
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
//						cltParameters.transform_size,
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
	public double [][] clt_dbg( // not used in lwir
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

	void clt_convert_double_kernel( // converts double resolution kernel // not used in lwir
			double []   src_kernel, //
			double []   dst_kernel, // should be (2*dtt_size-1) * (2*dtt_size-1) + extra_items size - kernel and dx, dy to the nearest 1/2 pixels + actual full center shift)
			int src_size) //, // 64
//			int transform_size) // 8
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
		if      (src_x < 2 * transform_size)             src_x = 2 * transform_size - 1; // 15
		else if (src_x > (src_size - 2* transform_size)) src_x = src_size - 2* transform_size;

		if      (src_y < 2 * transform_size)             src_y = 2 * transform_size - 1; // 15
		else if (src_y > (src_size - 2* transform_size)) src_y = src_size - 2* transform_size;
		indx = 0;
		// downscale, copy
		for (int i = -transform_size + 1; i < transform_size; i++){
			int src_i = (src_y + 2 * i) * src_size  + src_x;
			for (int j = -transform_size + 1; j < transform_size; j++){
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

	void clt_normalize_kernel(  // not used in lwir
			double []   kernel, // should be (2*dtt_size-1) * (2*dtt_size-1) + 4 size (last (2*dtt_size-1) are not modified)
			double []   window, // normalizes result kernel * window to have sum of elements == 1.0
//			int transform_size, // 8
			boolean bdebug)
	{
		double s = 0.0;
		int indx = 0;
		for (int i = -transform_size + 1; i < transform_size; i++){
			int ai = (i < 0)? -i: i;
			for (int j = -transform_size + 1; j < transform_size; j++){
				int aj = (j < 0)? -j: j;
				s += kernel[indx++] * window[ai*transform_size+aj];
			}
		}
		s = 1.0/s;
		int klen = (2*transform_size-1) * (2*transform_size-1);
		if (bdebug)		System.out.println("clt_normalize_kernel(): s="+s);
		for (int i = 0; i < klen; i++) {
//******************** Somewhere scale 16 ? ********************
			kernel[i] *= 16*s;
		}
 	}

	void clt_symmetrize_kernel(  // not used in lwir
			double []     kernel,      // should be (2*dtt_size-1) * (2*dtt_size-1) +2 size (last 2 are not modified)
			double [][]   sym_kernels) //, // set of 4 SS, AS, SA, AA kdernels, each dtt_size * dtt_size (may have 5-th with center shift
///			final int     transform_size) // 8
	{
		int in_size = 2*transform_size-1;
		int dtt_size_m1 = transform_size - 1;
		int center = dtt_size_m1 * in_size + dtt_size_m1;

		for (int i = 0; i < transform_size; i++){
			for (int j = 0; j < transform_size; j++){
				int indx0 = center - i * in_size - j;
				int indx1 = center - i * in_size + j;
				int indx2 = center + i * in_size - j;
				int indx3 = center + i * in_size + j;
				sym_kernels[0][i*transform_size+j] =                                 0.25*( kernel[indx0] + kernel[indx1] + kernel[indx2] + kernel[indx3]);
				if (j > 0)              sym_kernels[1][i*transform_size+j-1] =       0.25*(-kernel[indx0] + kernel[indx1] - kernel[indx2] + kernel[indx3]);
				if (i > 0)              sym_kernels[2][(i-1)*transform_size+j] =     0.25*(-kernel[indx0] - kernel[indx1] + kernel[indx2] + kernel[indx3]);
				if ((i > 0) && (j > 0)) sym_kernels[3][(i-1)*transform_size+(j-1)] = 0.25*(-kernel[indx0] + kernel[indx1] - kernel[indx2] + kernel[indx3]);
			}
			sym_kernels[1][i*transform_size + dtt_size_m1] = 0.0;
			sym_kernels[2][dtt_size_m1*transform_size + i] = 0.0;
			sym_kernels[3][i*transform_size + dtt_size_m1] = 0.0;
			sym_kernels[3][dtt_size_m1*transform_size + i] = 0.0;
		}
 	}

	void clt_dtt3_kernel( // not used in lwir
			double [][]   kernels, // set of 4 SS, AS, SA, AA kdernels, each dtt_size * dtt_size (may have 5-th with center shift
			final int     dtt_size, // 8
			DttRad2       dtt)
	{
		if (dtt == null) dtt = new DttRad2(dtt_size);
		for (int quadrant = 0; quadrant < 4; quadrant ++){
			kernels[quadrant] = dtt.dttt_iiie(kernels[quadrant], quadrant, dtt_size);
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

		public CltExtra(){} // not used in lwir
		public CltExtra(double [] data)  // USED in lwir
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
		public double [] getArray() // not used in lwir
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

	public void offsetKernelSensor( // not used in lwir
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
	public void clt_fill_coord_corr( // not used in lwir
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

	public class CltTile{ // not used in lwir
		public double [][] tile = new double[4][]; // 4 CLT tiles
		public double fract_x; // remaining fractional offset X
		public double fract_y; // remaining fractional offset X
	}



// Extract and correct one image tile using kernel data, required result tile and shifts - x and y
// option - align to Bayer (integer shift by even pixels - no need
// input - RBG stack of sparse data
// return
// kernel [0][0] is centered at  (-kernel_step/2,-kernel_step/2)

	public double [] extract_correct_tile( // return a pair of residual offsets // USED in lwir (except debug, fpga, gpu)
			double [][]         image_data,
			int                 width,       // image width
			double  [][][][][]  clt_kernels, // [color][tileY][tileX][band][pixel]
			double  [][]        clt_tile,    // should be double [4][];
			int                 kernel_step,
//			int                 transform_size,
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
///					transform_size,
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
	public void convolve_tile( // USED in lwir
			double [][]     data,    // array [transform_size*transform_size], will be updated  DTT4 converted
			double [][]     kernel,  // array [4][transform_size*transform_size]  DTT3 converted
///			int             transform_size,
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

		final int transform_len = transform_size * transform_size;
		final double [][] rslt = new double[4][transform_len];
		for (int i = 0; i < transform_len; i++) {
			for (int n = 0; n<4; n++){
				rslt[n][i] = 0;
				for (int k=0; k<4; k++){
//					if (ZI[n][k] < 0)
//						rslt[n][i] -= data[-ZI[n][k]][i] * kernel[k][i];
//					else
//						rslt[n][i] += data[ ZI[n][k]][i] * kernel[k][i];
					if (ZI[k][n] < 0)
					    rslt[n][i] -= data[-ZI[k][n]][i] * kernel[k][i]; // data[-ZI[n][k]][i] * kernel[k][i];
                    else
                        rslt[n][i] += data[ ZI[k][n]][i] * kernel[k][i]; // data[ ZI[n][k]][i] * kernel[k][i];
					
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

	public void fract_shift(    // fractional shift in transform domain. Currently uses sin/cos - change to tables with 2? rotations // USED in lwir
		double  [][]  clt_tile,
//		int           transform_size,
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



	public double [][][][] mdctScale( // not used in lwir
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



	public double [][][] lapped_dct_scale( // scale image to 8/9 size in each direction // not used in lwir
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

	public void dct_scale( // not used in lwir
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
	public static Thread[] newThreadArray(int maxCPUs) { // USED in lwir
		int n_cpus = Runtime.getRuntime().availableProcessors();
		if (n_cpus>maxCPUs)n_cpus=maxCPUs;
		return new Thread[n_cpus];
	}
/* Start all given threads and wait on each of them until all are done.
	 * From Stephan Preibisch's Multithreading.java class. See:
	 * http://repo.or.cz/w/trakem2.git?a=blob;f=mpi/fruitfly/general/MultiThreading.java;hb=HEAD
	 */
	public static void startAndJoin(Thread[] threads) // USED in lwir
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

	public double [] tileInterCamCorrs( // not used in lwir
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
	    		getScaleStrengths(),
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

	public double [] tileInterCamCorrs( // not used in lwir
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
			strength = stripe_inter[ixy[0]+transform_size-1]; // strength at integer max on axis
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
	public double [] tileCorrs( // USED in lwir
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
	    		getScaleStrengths(),
	    		col_weights,    // double []           col_weights,
	    		fatzero);       // double              fat_zero)

	    if (ml_center_corr != null) { // USED in lwir
	    	corr2d.corrCenterValues(
	    			ml_hwidth,                           // int         hwidth,
	    			clt_parameters.img_dtt.corr_offset,  //double      offset,
	        		corrs,                               // double [][] full_corr,
	        		ml_center_corr);                     // double [][] center_corr)
	    }

	    // calculate interpolated "strips" to match different scales and orientations (ortho/diagonal) on the
	    // fine (0.5 pix) grid. ortho for scale == 1 provide even/even samples (1/4 of all), diagonal ones -
	    // checkerboard pattern
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
		if (ixy != null) { // USED in lwir
			strength = strip_combo[ixy[0]+transform_size-1]; // strength at integer max on axis
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
		if (corr_stat != null) { // USED in lwir
			disparity = -corr_stat[0];
			result[DISP_FULL_INDEX] = disparity;
			// see if strength is enough to proceed with LMA/poly (otherwise keep disp/strength
			if (strength > clt_parameters.img_dtt.min_poly_strength) {
				// create LMA instance, calculate LMA composite argmax
		    	// Create 2 groups: ortho & diag
		    	Correlations2dLMA lma = corr2d.corrLMA(
		    			clt_parameters.img_dtt,                // ImageDttParameters  clt_parameters.img_dtt,
		    			corrs,                        // double [][]         corrs,
		    			corr2d.longToArray(clt_parameters.img_dtt.dbg_pair_mask), // clt_parameters.img_dtt.dbg_pair_mask,  // int                 pair_mask, // which pairs to process
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
					if (clt_parameters.rig.use_poly) { // not used in lwir
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
					if (get4dirs) { // USED in lwir
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
			    			System.out.println(String.format("corr4dirsLMA -> ↔: %7.4f (%7.4f) ↕:%7.4f (%7.4f) ⇖:%7.4f (%7.4f) ⇗:%7.4f (%7.4f)",
			    					dir_corr_strength[0][0],dir_corr_strength[0][1],dir_corr_strength[1][0],dir_corr_strength[1][1],
			    					dir_corr_strength[2][0],dir_corr_strength[2][1],dir_corr_strength[3][0],dir_corr_strength[3][1]));
			    		}
			    		if (dir_corr_strength[0] != null) { // USED in lwir
			    			result[DISP_HOR_INDEX] =   dir_corr_strength[0][0];
			    			result[STR_HOR_INDEX] =    dir_corr_strength[0][1];
			    		} else { // USED in lwir
			    			result[DISP_HOR_INDEX] =   Double.NaN;
			    			result[STR_HOR_INDEX] =    0.0;
			    		}
			    		if (dir_corr_strength[1] != null) { // USED in lwir
			    			result[DISP_VERT_INDEX] =  dir_corr_strength[1][0];
			    			result[STR_VERT_INDEX] =   dir_corr_strength[1][1];
			    		} else { // USED in lwir
			    			result[DISP_VERT_INDEX] =  Double.NaN;
			    			result[STR_VERT_INDEX] =   0.0;
			    		}
			    		if (dir_corr_strength[2] != null) { // USED in lwir
			    			result[DISP_DIAGM_INDEX] = dir_corr_strength[2][0];
			    			result[STR_DIAGM_INDEX] =  dir_corr_strength[2][1];
			    		} else { // USED in lwir
			    			result[DISP_DIAGM_INDEX] = Double.NaN;
			    			result[STR_DIAGM_INDEX] =  0.0;
			    		}
			    		if (dir_corr_strength[3] != null) { // USED in lwir
			    			result[DISP_DIAGO_INDEX] = dir_corr_strength[3][0];
			    			result[STR_DIAGO_INDEX] =  dir_corr_strength[3][1];
			    		} else { // USED in lwir
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


	// old version with common lpf filter
	public void generateTextureTiles(// not used in lwir
			final CLTParameters    clt_parameters,
			final double           extra_disparity,
			final int              nSens,      // number of subcameras
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
			final int              debugLevel,
			final boolean          debug_gpu
			) {
		// old version with common lpf filter
		final double [][] lpf_rgb = new double [numcol][];
		for (int i = 0; i < numcol; i++) {
			lpf_rgb[i] = filter;
		}
		generateTextureTiles(// not used in lwir
				null, // final double []        ports_rgb,      // average values of R,G,B for each camera (R0,R1,...,B2,B3)
				null, // final double []        max_diff,       // maximal (weighted) deviation of each channel from the average
				clt_parameters,
				extra_disparity,
				nSens,      // number of subcameras
				numcol, // number of colors
				img_mask,
				tile_op,         // [tilesY][tilesX] - what to do - 0 - nothing for this tile
				clt_data,
				texture_tiles,   // [tilesY][tilesX]["RGBA".length()][];  null - will skip images combining
				lpf_rgb,// final double []        filter,
				lt_window2,
				port_offsets,
				col_weights,
				dtt,
				tileX, // only used in debug output
				tileY,
				debugLevel,
				debug_gpu);
	}

	public void generateTextureTiles(// not used in lwir
			final double []        ports_rgb,      // average values of R,G,B for each camera (R0,R1,...,B2,B3)
			final double []        max_diff,       // maximal (weighted) deviation of each channel from the average
			final CLTParameters    clt_parameters,
			final double           extra_disparity,
			final int              nSens,      // number of subcameras
			final int              numcol, // number of colors
			int                    img_mask,
			final int [][]         tile_op,         // [tilesY][tilesX] - what to do - 0 - nothing for this tile
			final double [][][][]  clt_data,
			final double [][][][]  texture_tiles,   // [tilesY][tilesX]["RGBA".length()][];  null - will skip images combining
			final double [][]      lpf_rgb,// final double []        filter,
			final double []        lt_window2,
			final double[][]       port_offsets,
			final double []        col_weights,
			final DttRad2          dtt,
			final int              tileX, // only used in debug output
			final int              tileY,
			final int              debugLevel,
			final boolean          debug_gpu
			) {
		if (debug_gpu) {
			System.out.println("Generating GPU debug data for tileX="+tileX+", tileY="+tileY);
		}
		if ((extra_disparity != 0) && !getForcedDisparity(tile_op[tileY][tileX])){ // 0 - adjust disparity, 1 - use provided (now do not adjust!)
			// shift images by 0.5 * extra disparity in the diagonal direction
			for (int chn = 0; chn <numcol; chn++) { // color
				for (int i = 0; i < nSens; i++) {
					fract_shift(    // fractional shift in transform domain. Currently uses sin/cos - change to tables with 2? rotations
							clt_data[i][chn], // [tileY][tileX], // double  [][]  clt_tile,
//							transform_size,
							extra_disparity * port_offsets[i][0] / clt_parameters.corr_magic_scale,     // double        shiftX,
							extra_disparity * port_offsets[i][1] / clt_parameters.corr_magic_scale,     // double        shiftY,
							debugLevel > 0);
				}
			}
		}
		// lpf tiles (same as images before)
		// iclt tiles

		if (debug_gpu) {
			// just one camera, all colors
			for (int cam = 0; cam < 1; cam++) { // nSens; i++) {
				for (int chn = 0; chn < numcol; chn++) { // color
					for (int dct_mode = 0; dct_mode < 4; dct_mode++){
					System.out.println("=== CLT camera="+cam+", color="+chn+" mode="+dct_mode+" ===");
						for (int i = 0; i < transform_size; i++) {
							for (int j = 0; j < transform_size; j++) {
								System.out.print(String.format("%10.5f ", clt_data[cam][chn][dct_mode][transform_size * i + j]));
							}
							System.out.println();
						}
					}
				}
			}
		}



		double [][][] iclt_tile = new double [nSens][numcol][];
		double [] clt_tile;
		double scale = 0.25;  // matching iclt_2d ???
		for (int i = 0; i < nSens; i++) {
			for (int chn = 0; chn <numcol; chn++) { // color
				// double [] clt_tile = new double [transform_size*transform_size];
				for (int dct_mode = 0; dct_mode < 4; dct_mode++){
					clt_tile = clt_data[i][chn][dct_mode].clone();
					// lpf each of the 4 quadrants before idct
					if (lpf_rgb == null) {
						for (int j = 0; j < clt_tile.length; j++){
							clt_tile[j] *= scale;
						}
					} else {
						for (int j = 0; j < clt_tile.length; j++){
							clt_tile[j] *= scale*lpf_rgb[chn][j]; // filter[j];
						}
					}

					if (debug_gpu) {
						if ((i ==0) && (dct_mode == 0)) {
							System.out.println("=== Color LPF , color="+chn+" ===");
							for (int ii = 0; ii < transform_size; ii++) {
								for (int j = 0; j < transform_size; j++) {
									System.out.print(String.format("%10.5f ", lpf_rgb[chn][transform_size * ii + j]));
								}
								System.out.println();
							}
						}

						// just one camera, all colors
						System.out.println("=== CLT-scaledLPF camera="+i+", color="+chn+" mode="+dct_mode+" ===");
						for (int ii = 0; ii < transform_size; ii++) {
							for (int j = 0; j < transform_size; j++) {
								System.out.print(String.format("%10.5f ", clt_tile[transform_size * ii + j]));
							}
							System.out.println();
						}
					}



					// IDCT-IV should be in reversed order: CC->CC, SC->CS, CS->SC, SS->SS
					int idct_mode = ((dct_mode << 1) & 2) | ((dct_mode >> 1) & 1);
					clt_tile = dtt.dttt_iv  (clt_tile, idct_mode, transform_size);
					if (debug_gpu) {
						System.out.println("=== IDTT camera="+i+", color="+chn+" mode="+dct_mode+" idct_mode="+idct_mode+" ===");
						for (int ii = 0; ii < transform_size; ii++) {
							for (int j = 0; j < transform_size; j++) {
								System.out.print(String.format("%10.5f ", clt_tile[transform_size * ii + j]));
							}
							System.out.println();
						}
					}

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
		if (debug_gpu) {
			System.out.println("=============== IMCLT done  ===============");
			// just one camera, all colors
			for (int cam = 0; cam < nSens; cam++) { // nSens; i++) {
				for (int chn = 0; chn < numcol; chn++) { // color
					System.out.println("=== IMCLT camera="+cam+", color="+chn+" ===");
						for (int i = 0; i < 2 * transform_size; i++) {
							for (int j = 0; j < 2 * transform_size; j++) {
								System.out.print(String.format("%10.4f ", iclt_tile[cam][chn][2* transform_size * i + j]));
							}
							System.out.println();
						}
				}
			}
		}
		if ((debugLevel > 0) || debug_gpu) {
			ShowDoubleFloatArrays sdfa_instance = new ShowDoubleFloatArrays(); // just for debugging?
			String [] titles = {"red0","blue0","green0","red1","blue1","green1","red2","blue2","green2","red3","blue3","green3"};
			double [][] dbg_tile = new double [nSens*numcol][];
			for (int i = 0; i < nSens; i++) {
				for (int chn = 0; chn <numcol; chn++) { // color
					dbg_tile[i * numcol + chn] = iclt_tile[i][chn];
				}
			}
			sdfa_instance.showArrays(dbg_tile, 2* transform_size, 2* transform_size, true, "iclt_x"+tileX+"_y"+tileY, titles);
		}


		// "de-bayer" tiles for matching, use original data for output
		double [][][] tiles_debayered = new double [nSens][numcol][];
		for (int i =0; i<nSens; i++){
			for (int chn = 0; chn < numcol; chn++){
				tiles_debayered[i][chn] =  tile_debayer_shot_corr(
						(chn != 2),                 // red or blue (false - green)
						iclt_tile[i][chn],
						2 * transform_size,
						lt_window2,                 // squared lapping window
						clt_parameters.min_shot,    // 10.0;  // Do not adjust for shot noise if lower than
						clt_parameters.scale_shot,  //3.0;   // scale when dividing by sqrt
						lt_window2,                 // re-apply window to the result
						debug_gpu); // false); // debug_gpu && (i == 0));     //boolean  debug_gpu);
			}
		}
		if ((debugLevel > 0) || debug_gpu) {
			ShowDoubleFloatArrays sdfa_instance = new ShowDoubleFloatArrays(); // just for debugging?
			String [] titles = {"red0","blue0","green0","red1","blue1","green1","red2","blue2","green2","red3","blue3","green3"};
			double [][] dbg_tile = new double [nSens*numcol][];
			for (int i = 0; i < nSens; i++) {
				for (int chn = 0; chn <numcol; chn++) { // color
					dbg_tile[i * numcol + chn] = tiles_debayered[i][chn];
				}
			}
			sdfa_instance.showArrays(dbg_tile, 2* transform_size, 2* transform_size, true, "tiles_debayered_x"+tileX+"_y"+tileY, titles);
		}
		if (debug_gpu) {
			System.out.println("\n=== BEFORE tile_combine_rgba ===");
			for (int ccam = 0; ccam < nSens; ccam++) {
				for (int nncol = 0; nncol < numcol; nncol++) {
					System.out.println(String.format("--- before tile_combine_rgba cam = %d, color=%d ---", ccam, nncol));
					for (int i = 0; i < 2 * transform_size; i++) {
						for (int j = 0; j < 2 * transform_size; j++) {
							System.out.print(String.format("%10.4f ", tiles_debayered[ccam][nncol][2* transform_size * i + j]));
						}
						System.out.println();
					}
				}
			}
		}

		texture_tiles[tileY][tileX] =  tile_combine_rgba(
				tiles_debayered,                // iclt_tile,      // [port][numcol][256]
				ports_rgb, // null,                           // double []     ports_rgb,      // average values of R,G,B for each camera (R0,R1,...,B2,B3)
				max_diff,  // null,                           // max_diff,        // maximal (weighted) deviation of each channel from the average
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
				(debugLevel > 0),
				debug_gpu);          // boolean       debug_gpu)      // generate output for matching with GPU processing

		if (debug_gpu) {
			ShowDoubleFloatArrays sdfa_instance = new ShowDoubleFloatArrays(); // just for debugging?
//			String [] titles = {"red0","blue0","green0","red1","blue1","green1","red2","blue2","green2","red3","blue3","green3"};
			sdfa_instance.showArrays(texture_tiles[tileY][tileX], 2* transform_size, 2* transform_size, true, "tile_combine_rgba_x"+tileX+"_y"+tileY); //, titles);


			// just one camera, all colors
			for (int chn = 0; chn < texture_tiles[tileY][tileX].length; chn++) { // color
				System.out.println("=== AFTER tile_combine_rgba, chn="+chn+" ===");
				for (int i = 0; i < 2 * transform_size; i++) {
					for (int j = 0; j < 2 * transform_size; j++) {
						System.out.print(String.format("%10.4f ", texture_tiles[tileY][tileX][chn][2* transform_size * i + j]));
					}
					System.out.println();
				}
			}
		}

		// mix RGB from iclt_tile, mix alpha with - what? correlation strength or 'don't care'? good correlation or all > min?
		if (lpf_rgb != null) {
			for (int i = 0; i < iclt_tile[0][0].length; i++ ) {
				double sw = 0.0;
				for (int ip = 0; ip < nSens; ip++) {
					sw += texture_tiles[tileY][tileX][numcol+1+ip][i];
				}
				if (sw != 0 ) sw = 1.0/sw;
				for (int chn = 0; chn <numcol; chn++) { // color
					texture_tiles[tileY][tileX][chn][i] = 0.0; //iclt[tileY][tileX][chn]
					for (int ip = 0; ip < nSens; ip++) {
						texture_tiles[tileY][tileX][chn][i] += sw * texture_tiles[tileY][tileX][numcol+1+ip][i] * iclt_tile[ip][chn][i];
					}
				}
			}
		}
		if (debug_gpu) {
			ShowDoubleFloatArrays sdfa_instance = new ShowDoubleFloatArrays(); // just for debugging?
//			String [] titles = {"red0","blue0","green0","red1","blue1","green1","red2","blue2","green2","red3","blue3","green3"};
			sdfa_instance.showArrays(texture_tiles[tileY][tileX], 2* transform_size, 2* transform_size, true, "final_tile_combine_rgba_x"+tileX+"_y"+tileY); //, titles);
			for (int chn = 0; chn < texture_tiles[tileY][tileX].length; chn++) { // color
				System.out.println("=== AFTER combining, chn="+chn+" ===");
				for (int i = 0; i < 2 * transform_size; i++) {
					for (int j = 0; j < 2 * transform_size; j++) {
						System.out.print(String.format("%10.4f ", texture_tiles[tileY][tileX][chn][2* transform_size * i + j]));
					}
					System.out.println();
				}
			}
		}
		// fix: removing extra slices
		if (!clt_parameters.keep_weights && (texture_tiles[tileY][tileX]!=null)) {
			if (numcol == 3 ) {
				texture_tiles[tileY][tileX] = new double[][] {texture_tiles[tileY][tileX][0],texture_tiles[tileY][tileX][1],texture_tiles[tileY][tileX][2],texture_tiles[tileY][tileX][3]};
			} else {
				texture_tiles[tileY][tileX] = new double[][] {texture_tiles[tileY][tileX][0],texture_tiles[tileY][tileX][1]};
			}
		}
	}

	public double [][] get2DCorrs(
			final CLTParameters       clt_parameters,
			final double [][][][][][] clt_data, // [channel_in_quad][color][tileY][tileX][band][pixel];
			final int    []           pairs_list,
			final int                 threadsMax,  // maximal number of threads to launch
			final int                 debugLevel
			){
		final int debug_threshold = 0;
		final boolean debug_gpu = (debugLevel > -2);// (debugLevel == 2);
		final int num_corrs = pairs_list.length;
//		final int num_cams = clt_data.length;
		int   nc = 0;
		int   mc = -1;
		for (int i = 0; i < clt_data[0].length; i++){
			if (clt_data[0][i] != null) {
				nc++;
				if (mc<0) mc = i;
			}
		}
		final int colors = nc;
		final int mono_color = mc;

//		final int tilesY = clt_data[0][mono_color].length;
		final int tilesX = clt_data[0][mono_color][0].length;
		final boolean is_mono = (colors < 3); //  isMonochrome(); //
		final double [][] corrs2d = new double [num_corrs][];
		final Thread[] threads = newThreadArray(threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
///		final int transform_size = transform_size;
		final int transform_len =  transform_size*transform_size;
		final double afat_zero2 = clt_parameters.getGpuFatZero(is_mono); // ? clt_parameters.gpu_fatz_m : clt_parameters.gpu_fatz;
		final double sigma_corr = clt_parameters.getGpuCorrSigma(is_mono); // ? clt_parameters.gpu_sigma_corr_m : clt_parameters.gpu_sigma_corr;
		final double sigma_rb_corr = clt_parameters.getGpuCorrRBSigma(is_mono);
		final int corr_radius = clt_parameters.gpu_corr_rad;
//		final double[] lpf =  getLpf(sigma_corr);

		final double[] lpf_fd =  doubleGetCltLpfFd(
				sigma_corr);
		final double[] lpf_rb_fd =  doubleGetCltLpfFd(
				sigma_rb_corr);

		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				@Override
				public void run() {
					double [][] tcorr = null;
					Correlation2d correlation2d = new  Correlation2d (
							numSensors,
							clt_parameters.img_dtt,
							transform_size,
							1.0, // double wndx_scale, // (wndy scale is always 1.0)
							is_mono,
							debugLevel > debug_threshold); //  boolean debug) {
					for (int n_corr = ai.getAndIncrement(); n_corr < num_corrs; n_corr = ai.getAndIncrement()) {
						int nTile = pairs_list[n_corr] >> GPUTileProcessor.CORR_NTILE_SHIFT;
					int n_pair = pairs_list[n_corr] & GPUTileProcessor.CORR_PAIRS_MASK; // ((1 <<  GPUTileProcessor.CORR_NTILE_SHIFT) -1);
						int tileY = nTile /tilesX;
						int tileX = nTile % tilesX;

						boolean debug_tile = (tileY == clt_parameters.tileY) && (tileX == clt_parameters.tileX) && (n_pair == 0);
						if (debug_tile) {
							System.out.println("Debugging tile "+nTile+", tileX="+tileX+", tileY="+tileY);
						}
						int cam1 = Correlation2d.PAIRS[n_pair][0];
						int cam2 = Correlation2d.PAIRS[n_pair][1];
						if (is_mono) {
							tcorr = correlation2d.correlateSingleColorFD(
									clt_data[cam1][mono_color][tileY][tileX], //  double [][] clt_data1,
									clt_data[cam2][mono_color][tileY][tileX], // double [][] clt_data2,
						      		null);
						} else {
							tcorr = new double [4][transform_len];
							double [] weights = {
									clt_parameters.gpu_weight_r,
									clt_parameters.gpu_weight_b,
									1.0 - clt_parameters.gpu_weight_r - clt_parameters.gpu_weight_b};
							for (int c = 0; c < colors; c++) {
								if (debug_gpu && debug_tile) {
									System.out.println("\n=== CLT CAMERA1="+cam1+" color="+c+" ===");
									for (int dct_mode = 0; dct_mode < 4; dct_mode++) {
										System.out.println("------dct_mode="+dct_mode);
										for (int i = 0; i < transform_size; i++) {
											for (int j = 0; j < transform_size; j++) {
												System.out.print(String.format("%10.5f ", clt_data[cam1][c][tileY][tileX][dct_mode][transform_size * i + j]));
											}
											System.out.println();
										}
									}
									System.out.println("\n=== CLT CAMERA2="+cam2+" color="+c+" ===");
									for (int dct_mode = 0; dct_mode < 4; dct_mode++) {
										System.out.println("------dct_mode="+dct_mode);
										for (int i = 0; i < transform_size; i++) {
											for (int j = 0; j < transform_size; j++) {
												System.out.print(String.format("%10.5f ", clt_data[cam2][c][tileY][tileX][dct_mode][transform_size * i + j]));
											}
											System.out.println();
										}
									}

								}
								double [][] tcorri = correlation2d.correlateSingleColorFD(
										clt_data[cam1][c][tileY][tileX], //  double [][] clt_data1,
										clt_data[cam2][c][tileY][tileX], // double [][] clt_data2,
							      		null);
								if (tcorri != null) {
									for (int n = 0; n < 4; n++) {
										for (int i = 0; i < transform_len; i++) {
											tcorr[n][i]+= tcorri[n][i]* weights[c];
										}
									}
								}
								if (debug_gpu && debug_tile) {
									System.out.println("=== ACCUMMULATED CORRELATION ===");
									for (int dct_mode = 0; dct_mode < 4; dct_mode++) {
										System.out.println("------dct_mode="+dct_mode);
										for (int i = 0; i < transform_size; i++) {
											for (int j = 0; j < transform_size; j++) {
												System.out.print(String.format("%10.3f ", tcorr[dct_mode][transform_size * i + j]));
											}
											System.out.println();
										}
									}
								}
								//c
								if (c == 1) { // after red and blue, color only
									if (lpf_rb_fd != null) {
										for (int n = 0; n<4; n++) {
											for (int i = 0; i < transform_len; i++) {
												tcorr[n][i] *= lpf_rb_fd[i];
											}
										}
										if (debug_gpu && debug_tile) {
											System.out.println("=== LPF-RB for CORRELATION ===");
											System.out.print("__constant__ float lpf_rb_corr[64]={");
											for (int i=0; i<lpf_rb_fd.length;i++){
												System.out.print(String.format("%10.8ff", lpf_rb_fd[i]));
												if (i == 63) {
													System.out.println("};");
												} else {
													System.out.print(", ");
													if ((i % 8) == 7) {
														System.out.print("\n                                 ");
													}
												}
											}

											System.out.println("=== R+B LPF-ed CORRELATION ===");
											for (int dct_mode = 0; dct_mode < 4; dct_mode++) {
												System.out.println("------dct_mode="+dct_mode);
												for (int i = 0; i < transform_size; i++) {
													for (int j = 0; j < transform_size; j++) {
														System.out.print(String.format("%10.3f ", tcorr[dct_mode][transform_size * i + j]));
													}
													System.out.println();
												}
											}
										}
									}
								}
							}
						}
						corrs2d[n_corr] = correlation2d.normalizeConvertCorr(
					    		  tcorr, // null or initialized to [4][transform_len]
					    		  lpf_fd,
					    		  afat_zero2, // double      afat_zero2, // absolute fat zero, same units as components squared values
					    		  corr_radius, // int corr_radius);
					    		  debug_gpu && debug_tile); // boolean debug_gpu); //
					}
				}
			};
		}
		startAndJoin(threads);
		return corrs2d;
	}

	// Used quad
	@Deprecated
	public double [][][][][][][]  clt_bi_quad_dbg(// not used in lwir
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
// added in debug version
			final double [][]         disparity_map,   // [8][tilesY][tilesX], only [6][] is needed on input or null - do not calculate

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
		final int quad = getNumSensors(); // 4;   // number of subcameras
		final int numcol = 3; // number of colors
		final int nChn = image_data_main[0].length;
		final int height=image_data_main[0][0].length/width;
		final int tilesX=width/transform_size;
		final int tilesY=height/transform_size;
		final int nTilesInChn=tilesX*tilesY;
		final double [][] geom_dbg = new double [24][nTilesInChn];
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
		final int corr_size = transform_size * 2 -1;
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
		int wcenter = transform_size - 1;
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
			System.out.println("clt_aberrations_quad_corr(): width="+width+" height="+height+" transform_size="+transform_size+
					" debug_tileX="+debug_tileX+" debug_tileY="+debug_tileY+" globalDebugLevel="+globalDebugLevel);
		}

		final double [] filter =  doubleGetCltLpfFd(clt_parameters.getCorrSigma(isMonochrome()));
		dbg_filter_corr = filter;


		// add optional initialization of debug layers here
		/*
		if (disparity_map != null){
			for (int i = 0; i<disparity_map.length;i++){
				if (i < OVEREXPOSED) {
					disparity_map[i] = new double [tilesY*tilesX];
				} else if (i == OVEREXPOSED) {
//					if (saturation_imp!= null) {
//					disparity_map[i] = new double [tilesY*tilesX];
//					}
				} else if (i >= IMG_TONE_RGB) {
//					if (texture_tiles != null) { // for now - enable 12 tone layers only together with texture tiles
						disparity_map[i] = new double [tilesY*tilesX];
//					}
				}
			}
		}
	 */
		if (disparity_map != null){
			for (int i = 0; i<disparity_map.length;i++){
				if (i == OVEREXPOSED) {
//					if (saturation_imp!= null) {
//						disparity_map[i] = new double [tilesY*tilesX];
//					}
				} else if (isSliceBit(i)) {
					disparity_map[i] = new double [tilesY*tilesX];
				} else if (isDiffIndex(i)) {
					disparity_map[i] = new double [tilesY*tilesX];
				} else if (isToneRGBIndex(i)) {
//					if (texture_tiles != null) { // for now - enable 12 tone layers only together with texture tiles
						disparity_map[i] = new double [tilesY*tilesX];
//					}
				}
			}
		}
		
		
		



		final double [][] lpf_rgb = new double[isMonochrome()?1:3][];
		if (isMonochrome()) {
			lpf_rgb[0] = doubleGetCltLpfFd(clt_parameters.gpu_sigma_m);
		} else {
			String [] comments = {"red", "blue", "green", "mono"};
			lpf_rgb[0] = doubleGetCltLpfFd(clt_parameters.gpu_sigma_r);
			lpf_rgb[1] = doubleGetCltLpfFd(clt_parameters.gpu_sigma_b);
			lpf_rgb[2] = doubleGetCltLpfFd(clt_parameters.gpu_sigma_g);
			if (globalDebugLevel > -2) {
				double [][] lpf_rbgm = {
						lpf_rgb[0],
						lpf_rgb[1],
						lpf_rgb[2],
						doubleGetCltLpfFd(clt_parameters.gpu_sigma_m)
				};


				System.out.println("__constant__ float lpf_data["+lpf_rbgm.length+"][64]={");
				for (int ncol = 0; ncol < lpf_rbgm.length; ncol++) {
					if (ncol==0) {
						System.out.print("		{");
					} else {
						System.out.print("		},{");
					}
					System.out.println(" // "+ comments[ncol]);
					for (int i=0; i<lpf_rbgm[ncol].length;i++){
						if ((i % 8) == 0) {
							System.out.print("				");
						}
						System.out.print(String.format("%5.8ff", lpf_rbgm[ncol][i]));
						if (i == (lpf_rbgm[ncol].length - 1)) {
							System.out.println();
						} else {
							System.out.print(", ");
							if ((i % 8) == 7) {
								System.out.print("\n");
							}
						}
					}
				}
				System.out.println("		}};");
			}
		}

		// prepare disparity maps and weights
		final int max_search_radius = (int) Math.abs(clt_parameters.max_corr_radius); // use negative max_corr_radius for squares instead of circles?
		final int max_search_radius_poly = 1;
		if (globalDebugLevel > -3) { // 0){
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

		DttRad2 dtt = new DttRad2(transform_size);
		dtt.set_window(clt_parameters.clt_window);
		final double [] lt_window = dtt.getWin2d();	// [256]
		final double [] lt_window2 = new double [lt_window.length]; // squared
		for (int i = 0; i < lt_window.length; i++) lt_window2[i] = lt_window[i] * lt_window[i];


		if (globalDebugLevel > 1) {
			ShowDoubleFloatArrays sdfa_instance = new ShowDoubleFloatArrays(); // just for debugging?
			sdfa_instance.showArrays(lt_window,  2*transform_size, 2*transform_size, "lt_window");
		}

		final Matrix [] corr_rots_main = geometryCorrection_main.getCorrVector().getRotMatrices(); // get array of per-sensor rotation matrices
		final Matrix rigMatrix = geometryCorrection_aux.getRotMatrix(true);
		final Matrix [] corr_rots_aux =  geometryCorrection_aux.getCorrVector().getRotMatrices(rigMatrix); // get array of per-sensor rotation matrices

		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				@Override
				public void run() {
					DttRad2 dtt = new DttRad2(transform_size);
					dtt.set_window(clt_parameters.clt_window);
					int tileY,tileX,tIndex; // , chn;
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
					double []   ml_data_dbg1 =  (ml_data != null)? new double [(2*ml_hwidth +1)*(2*ml_hwidth +1)]:null;


					Correlation2d corr2d = new Correlation2d(
							numSensors,
							clt_parameters.img_dtt,              // ImageDttParameters  imgdtt_params,
							transform_size,             // int transform_size,
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

						final int [] overexp_main = (saturation_main != null) ? ( new int [2]): null;
						final int [] overexp_aux =  (saturation_aux != null) ? ( new int [2]): null;

						// Moved from inside chn loop
						centerX = tileX * transform_size + transform_size/2; //  - shiftX;
						centerY = tileY * transform_size + transform_size/2; //  - shiftY;
						// TODO: move port coordinates out of color channel loop
						double [][] centersXY_main;
						double [][] centersXY_aux;
						double disparity_main = disparity_array[tileY][tileX];
						double disparity_aux =  disparity_main * geometryCorrection_aux.getDisparityRadius()/geometryCorrection_main.getDisparityRadius();
						if (disparity_bimap != null){
							disparity_bimap[BI_TARGET_INDEX][tIndex] = disparity_main;
						}
						double [][] disp_dist_main = new double[quad_main][]; // used to correct 3D correlations
						double [][] disp_dist_aux =  new double[quad_aux][]; // used to correct 3D correlations

						centersXY_main = geometryCorrection_main.getPortsCoordinatesAndDerivatives(
								geometryCorrection_main, //			GeometryCorrection gc_main,
								false,          // boolean use_rig_offsets,
								corr_rots_main, // Matrix []   rots,
								null,           //  Matrix [][] deriv_rots,
								null,           // double [][] pXYderiv, // if not null, should be double[8][]
								disp_dist_main, // used to correct 3D correlations
								centerX,
								centerY,
								disparity_main); //  + disparity_corr);
						if (geom_dbg != null) {
							for (int i = 0; i < quad_main; i++) {
								geom_dbg[2 * i + 0][nTile] = centersXY_main[i][0]; // x
								geom_dbg[2 * i + 1][nTile] = centersXY_main[i][1]; // y
								for (int j = 0; j < 4; j++) {
									geom_dbg[2 * quad_main + 4 * i + j][nTile] = disp_dist_main[i][j];
								}
							}
						}
						centersXY_aux =  geometryCorrection_aux.getPortsCoordinatesAndDerivatives(
								geometryCorrection_main, //			GeometryCorrection gc_main,
								true,            // boolean use_rig_offsets,
								corr_rots_aux,   // Matrix []   rots,
								null,            //  Matrix [][] deriv_rots,
								null,            // double [][] pXYderiv, // if not null, should be double[8][]
								disp_dist_aux,   // used to correct 3D correlations
								centerX,
								centerY,
								disparity_aux); //  + disparity_corr);

						if ((tileX == debug_tileX ) && (tileY == debug_tileY )) {
							// will just print debug data
							geometryCorrection_main.getPortsCoordinatesAndDerivativesDbg(
									geometryCorrection_main, //			GeometryCorrection gc_main,
									false,          // boolean use_rig_offsets,
									corr_rots_main, // Matrix []   rots,
									null,           //  Matrix [][] deriv_rots,
									null,           // double [][] pXYderiv, // if not null, should be double[8][]
									disp_dist_main,       // used to correct 3D correlations
									centerX,
									centerY,
									disparity_main); //  + disparity_corr);
						}



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

						boolean debug_gpu = (globalDebugLevel > -2) && (tileX == debug_tileX) && (tileY == debug_tileY);

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
//										transform_size,
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
//										transform_size,
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
								sdfa_instance.showArrays(dbg_tile,  transform_size, transform_size, true, "MAIN_pre-shifted_x"+tileX+"_y"+tileY, titles);

								double [][] dbg_tile_aux = new double [16][];
								for (int i = 0; i < 16; i++) dbg_tile[i]=clt_data_aux[i>>2][chn][i & 3];
								sdfa_instance.showArrays(dbg_tile_aux,  transform_size, transform_size, true, "AUX_pre-shifted_x"+tileX+"_y"+tileY, titles);
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
//							boolean debug_gpu = (globalDebugLevel > -2) && (tileX == debug_tileX) && (tileY == debug_tileY);
							for (int i = 0; i < quad_main; i++) {
								fract_shift(    // fractional shift in transform domain. Currently uses sin/cos - change to tables with 2? rotations
										clt_data_main[i][chn], // double  [][]  clt_tile,
										fract_shiftsXY_main[i][0],            // double        shiftX,
										fract_shiftsXY_main[i][1],            // double        shiftY,
										false //debug_gpu
										);
								if (debug_gpu) {
									System.out.println("---Shifted image tile for quad="+i+" color="+chn+", shift_hor = "+fract_shiftsXY_main[i][0]+", shift_vert = "+fract_shiftsXY_main[i][1]+"---");
									for (int dct_mode=0; dct_mode < 4; dct_mode++ ) {
										System.out.println("dct_mode="+dct_mode);
										for (int irow = 0; irow < transform_size; irow++) {
											for (int jcol = 0; jcol < transform_size; jcol++) {
												System.out.print(String.format("%10.5f ", clt_data_main[i][chn][dct_mode][transform_size * irow + jcol]));
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
								sdfa_instance.showArrays(dbg_tile,  transform_size, transform_size, true, "MAIN_shifted_x"+tileX+"_y"+tileY+"-z", titles);
								double [][] dbg_tile_aux = new double [16][];
								for (int i = 0; i < 16; i++) dbg_tile[i]=clt_data_aux[i>>2][chn][i & 3];
								sdfa_instance.showArrays(dbg_tile_aux,  transform_size, transform_size, true, "AUX_shifted_x"+tileX+"_y"+tileY+"-z", titles);
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
									col_weights,           // final double []       col_weights,
									ml_hwidth,             // final int             ml_hwidth,
						    		ml_data_main,          // final double [][]     ml_center_corr,
									tileX,                 // final int             tileX, // only used in debug output
									tileY,                 // final int             tileY,
									tile_lma_debug_level); // final int             debugLevel)

							double [] tile_corrs_aux  = tileCorrs(
									clt_parameters,        // final EyesisCorrectionParameters.CLTParameters       clt_parameters,
									fatzero,               // final double          fatzero,         // May use correlation fat zero from 2 different parameters - fat_zero and rig.ml_fatzero
									true,                  // final boolean         get4dirs, // calculate disparity/strength for each of the 4 directions
									corr2d,                // final Correlation2d   corr2d,
									clt_data_aux,          // final double [][][][] clt_data,
									filter,                // final double []       filter,
									col_weights,           // final double []       col_weights,
									ml_hwidth,             // final int             ml_hwidth,
						    		ml_data_aux,           // final double [][]     ml_center_corr,
									tileX,                 // final int             tileX, // only used in debug output
									tileY,                 // final int             tileY,
									tile_lma_debug_level); // final int             debugLevel)

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
						    		getScaleStrengths(),
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

						double []     max_diff = null;
						if ((disparity_map != null) && (disparity_map.length >= (IMG_DIFF0_INDEX + quad))){
							max_diff = new double[quad];
						}
						double [] ports_rgb = null;

						int ports_rgb_len = quad*numcol;  // 12
						if ((disparity_map != null) && (disparity_map.length >= (getImgToneRGB() + ports_rgb_len))) {
							ports_rgb = new double[ports_rgb_len];
						}



						if (texture_tiles_main !=null) {
							generateTextureTiles (
									ports_rgb,      //final double []        ports_rgb,      // average values of R,G,B for each camera (R0,R1,...,B2,B3)
									max_diff,       // final double []        max_diff,       // maximal (weighted) deviation of each channel from the average
									clt_parameters,        // final EyesisCorrectionParameters.CLTParameters       clt_parameters,
									extra_disparity_main,  // final double           extra_disparity,
									quad_main,                  // final int              quad,      // number of subcameras
									numcol,                // final int              numcol, // number of colors
									img_mask, // _main,         // int                    img_mask,
									tile_op, // _main,          // final int [][]         tile_op,         // [tilesY][tilesX] - what to do - 0 - nothing for this tile
									clt_data_main,         // final double [][][][]  clt_data,
									texture_tiles_main,    // final double [][][][]  texture_tiles,   // [tilesY][tilesX]["RGBA".length()][];  null - will skip images combining
									lpf_rgb, //null, // filter, //null               // final double []        filter,
									lt_window2,            // final double []        lt_window2,
									port_offsets,          // final double[][]       port_offsets,
									col_weights,           // final double []        col_weights,
									dtt,                   // final DttRad2          dtt,
									tileX,                 // final int              tileX, // only used in debug output
									tileY,                 // final int              tileY,
									tile_lma_debug_level,  // final int              debugLevel);
									debug_gpu); //boolean debug_gpu
						}
						if (texture_tiles_aux !=null) {
							generateTextureTiles (
									null,                  // final double []        ports_rgb,      // average values of R,G,B for each camera (R0,R1,...,B2,B3)
									null,                  // final double []        max_diff,       // maximal (weighted) deviation of each channel from the average
									clt_parameters,        // final EyesisCorrectionParameters.CLTParameters       clt_parameters,
									extra_disparity_aux,   // final double           extra_disparity,
									quad_aux,              // final int              quad,      // number of subcameras
									numcol,                // final int              numcol, // number of colors
									img_mask, // _aux,         // int                    img_mask,
									tile_op, // _aux,          // final int [][]         tile_op,         // [tilesY][tilesX] - what to do - 0 - nothing for this tile
									clt_data_aux,         // final double [][][][]  clt_data,
									texture_tiles_aux,    // final double [][][][]  texture_tiles,   // [tilesY][tilesX]["RGBA".length()][];  null - will skip images combining
									lpf_rgb, //  null, // filter,                // final double []        filter,
									lt_window2,            // final double []        lt_window2,
									port_offsets,          // final double[][]       port_offsets,
									col_weights,           // final double []        col_weights,
									dtt,                   // final DttRad2          dtt,
									tileX,                 // final int              tileX, // only used in debug output
									tileY,                 // final int              tileY,
									tile_lma_debug_level,  // final int              debugLevel);
									false); // debug_gpu); //boolean debug_gpu
						}


						if (debug_gpu) {
							double [][] texture_tile = texture_tiles_main[clt_parameters.tileY][clt_parameters.tileX];
							int tile = +clt_parameters.tileY * tilesX + +clt_parameters.tileX;
			    			System.out.println("=== inside clt_bi_quad_dbg(): tileX= "+clt_parameters.tileX+" tileY= "+clt_parameters.tileY+" tile="+tile+" ===");

			    			for (int slice =0; slice < texture_tile.length; slice++) {
			    				System.out.println("\n=== Slice (inside) ="+slice+" ===");
			    				for (int i = 0; i < 2 * GPUTileProcessor.DTT_SIZE; i++) {
			    					for (int j = 0; j < 2 * GPUTileProcessor.DTT_SIZE; j++) {
			    						System.out.print(String.format("%10.4f ",
			    								texture_tile[slice][2 * GPUTileProcessor.DTT_SIZE * i + j]));
			    					}
			    					System.out.println();
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
								disparity_map[getImgToneRGB() + i][tIndex] = ports_rgb[i];
							}
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
		if (geom_dbg != null) {
			(new ShowDoubleFloatArrays()).showArrays(
					geom_dbg,
					tilesX,
					tilesY,
					true,
					"geom_dbg",
					GEOM_TITLES_DBG);
		}
		
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
			//		final int corr_size = transform_size * 2 -1;

			final TileNeibs tnImage  = new TileNeibs(tilesX, tilesY);
			ai.set(0);
			for (int ithread = 0; ithread < threads.length; ithread++) {
				threads[ithread] = new Thread() {
					@Override
					public void run() {
						double []    ml_data_inter = (ml_data != null)? new double [(2*ml_hwidth +1)*(2*ml_hwidth +1)]:null;
						double [][]  tcorr_combo =     null; // [15*15] pixel space

						Correlation2d corr2d = new Correlation2d(
								numSensors,
								clt_parameters.img_dtt,              // ImageDttParameters  imgdtt_params,
								transform_size,             // int transform_size,
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

	public double [][][][][][][]  clt_bi_quad(// USED in lwir
			final CLTParameters       clt_parameters,
			final double              fatzero,         // May use correlation fat zero from 2 different parameters - fat_zero and rig.ml_fatzero
			final boolean             notch_mode,      // use notch filter for inter-camera correlation to detect poles
			// lt_rad = 0 for lwir
			final int                 lt_rad,          // low texture mode - inter-correlation is averaged between the neighbors before argmax-ing, using (2*notch_mode+1)^2 square
			final boolean             no_int_x0,       // do not offset window to integer maximum - used when averaging low textures to avoid "jumps" for very wide
			final int [][]            tile_op,         // [tilesY][tilesX] - what to do - 0 - nothing for this tile
			final double [][]         disparity_array, // [tilesY][tilesX] - individual per-tile expected disparity. If main camera is not calculated - use for AUX directly
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
		final boolean calc_main = (clt_kernels_main != null) && (image_data_main != null) && (geometryCorrection_main != null);
		final boolean calc_aux =  (clt_kernels_aux != null)  && (image_data_aux != null)  && (geometryCorrection_aux  != null);
		final boolean calc_both = calc_main && calc_aux;
		final int                 globalDebugLevel = clt_parameters.rig.rig_mode_debug?debugLevel:-2;
		final int                 debug_tileX = clt_parameters.tileX;
		final int                 debug_tileY = clt_parameters.tileY;
		final int quad_main = calc_main ? image_data_main.length: 0;   // number of subcameras
		final int quad_aux =  calc_aux ? image_data_aux.length : 0;   // number of subcameras
		final int numcol = 3; // number of colors
//		final int nColors = calc_main ? image_data_main[0].length : image_data_aux[0].length;
		final int height = (calc_main ? image_data_main[0][0].length: image_data_aux[0][0].length)/width;
		final int tilesX=width/transform_size;
		final int tilesY=height/transform_size;
		final int nTilesInChn=tilesX*tilesY;
		// clt_data does not need to be for the whole image (no, it is used for textures)
		final double [][][][][][][] clt_bidata = (keep_clt_data)? (new double[2][][][][][][]):null;
		if (clt_bidata != null) {// not used in lwir
			clt_bidata[0] = new double[quad_main][numcol][tilesY][tilesX][][];
			clt_bidata[1] = new double[quad_aux][numcol][tilesY][tilesX][][];
		}

		final double [][] lt_corr = (lt_rad > 0)? (new double [nTilesInChn][]):null; // will keep inter-camera combo correlation, later combined in a separate multi-thread run

		final Thread[] threads = newThreadArray(threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		final double [] col_weights= new double [numcol]; // colors are RBG

		if (isMonochrome()) {
			col_weights[2] = 1.0;// green color/mono
			col_weights[0] = 0;
			col_weights[1] = 0;
		} else {// not used in lwir
			col_weights[2] = 1.0/(1.0 + clt_parameters.corr_red + clt_parameters.corr_blue);    // green color
			col_weights[0] = clt_parameters.corr_red *  col_weights[2];
			col_weights[1] = clt_parameters.corr_blue * col_weights[2];
		}


		final int corr_size = transform_size * 2 -1;
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
		int wcenter = transform_size - 1;
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
			System.out.println("clt_aberrations_quad_corr(): width="+width+" height="+height+" transform_size="+transform_size+
					" debug_tileX="+debug_tileX+" debug_tileY="+debug_tileY+" globalDebugLevel="+globalDebugLevel);
		}

		final double[][] port_offsets = {
				{-0.5, -0.5},
				{ 0.5, -0.5},
				{-0.5,  0.5},
				{ 0.5,  0.5}};
		final double [] filter =  doubleGetCltLpfFd(clt_parameters.getCorrSigma(isMonochrome()));
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
			for (int i = 0; i < disparity_bimap.length;i++){// USED in lwir
				disparity_bimap[i] = new double [tilesY*tilesX];
			}
		}

		if (ers_delay != null) {
			ers_delay[0] = new double [quad_main][]; // not used in lwir
			for (int i = 0; i < quad_main; i++) ers_delay[0][i] = new double [tilesX*tilesY];
			ers_delay[1] = new double [quad_aux][];
			for (int i = 0; i < quad_aux; i++)  ers_delay[1][i] = new double [tilesX*tilesY];
		}

		DttRad2 dtt = new DttRad2(transform_size);
		dtt.set_window(clt_parameters.clt_window);
		final double [] lt_window = dtt.getWin2d();	// [256]
		final double [] lt_window2 = new double [lt_window.length]; // squared
		for (int i = 0; i < lt_window.length; i++) lt_window2[i] = lt_window[i] * lt_window[i];


		if (globalDebugLevel > 1) {
			ShowDoubleFloatArrays sdfa_instance = new ShowDoubleFloatArrays(); // just for debugging?
			sdfa_instance.showArrays(lt_window,  2*transform_size, 2*transform_size, "lt_window");
		}

		final Matrix [] corr_rots_main = calc_main?geometryCorrection_main.getCorrVector().getRotMatrices():null; // get array of per-sensor rotation matrices
		final Matrix rigMatrix = calc_aux? geometryCorrection_aux.getRotMatrix(true):null;
		final Matrix [] corr_rots_aux =  calc_aux? geometryCorrection_aux.getCorrVector().getRotMatrices(rigMatrix):null; // get array of per-sensor rotation matrices (will be OK if (rigMatrix==null)

		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				@Override
				public void run() {
					DttRad2 dtt = new DttRad2(transform_size);
					dtt.set_window(clt_parameters.clt_window);
					int tileY,tileX,tIndex; // , chn;
					//						showDoubleFloatArrays sdfa_instance = new showDoubleFloatArrays(); // just for debugging?
					double centerX; // center of aberration-corrected (common model) tile, X
					double centerY; //
					double [][] fract_shiftsXY_main = calc_main ? new double[quad_main][]:null;
					double [][] fract_shiftsXY_aux =  calc_aux ?  new double[quad_aux][]:null;
					double [][]     tcorr_combo =     null; // [15*15] pixel space
					double [][][][] clt_data_main =   calc_main ? new double[quad_main][numcol][][] : null;
					double [][][][] clt_data_aux =    calc_aux?  new double[quad_aux][numcol][][]:null;
					double [][] ml_data_main =  (calc_main && (ml_data != null)) ? new double [ML_VERT_INDEX+1][(2*ml_hwidth +1)*(2*ml_hwidth +1)]:null;
					double [][] ml_data_aux =   (calc_aux &&  (ml_data != null)) ? new double [ML_VERT_INDEX+1][(2*ml_hwidth +1)*(2*ml_hwidth +1)]:null;
					double []   ml_data_inter = (calc_both && (ml_data != null)) ? new double [(2*ml_hwidth +1)*(2*ml_hwidth +1)]:null;
					double []   ml_data_dbg1 =  (ml_data != null)? new double [(2*ml_hwidth +1)*(2*ml_hwidth +1)]:null;


					Correlation2d corr2d = new Correlation2d(
							numSensors,
							clt_parameters.img_dtt,              // ImageDttParameters  imgdtt_params,
							transform_size,             // int transform_size,
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
						int img_mask =       getImgMask(tile_op[tileY][tileX]);         // which images to use

						final int [] overexp_main = (calc_main && (saturation_main != null)) ? ( new int [2]): null;
						final int [] overexp_aux =  (calc_aux && (saturation_aux != null)) ?   ( new int [2]): null;

						// Moved from inside chn loop
						centerX = tileX * transform_size + transform_size/2; //  - shiftX;
						centerY = tileY * transform_size + transform_size/2; //  - shiftY;
						// TODO: move port coordinates out of color channel loop
						double [][] centersXY_main = null;
						double [][] centersXY_aux = null;
						double disparity_main = disparity_array[tileY][tileX]; // will be used for AUX directly if !calc_main
						double disparity_aux =  calc_both ? (disparity_main * geometryCorrection_aux.getDisparityRadius()/geometryCorrection_main.getDisparityRadius()) : 0.0;
						if (disparity_bimap != null){
							disparity_bimap[BI_TARGET_INDEX][tIndex] = disparity_main;
						}
						double [][] disp_dist_main = new double[quad_main][]; // used to correct 3D correlations
						double [][] disp_dist_aux =  new double[quad_aux][]; // used to correct 3D correlations
						if (calc_main) {// not used in lwir
							centersXY_main = geometryCorrection_main.getPortsCoordinatesAndDerivatives(
									geometryCorrection_main, //			GeometryCorrection gc_main,
									false,          // boolean use_rig_offsets,
									corr_rots_main, // Matrix []   rots,
									null,           //  Matrix [][] deriv_rots,
									null,           // double [][] pXYderiv, // if not null, should be double[8][]
									disp_dist_main, // used to correct 3D correlations
									centerX,
									centerY,
									disparity_main); //  + disparity_corr);
						}
						if (calc_aux) {
							if (calc_main) { // use rig and main coordinates // not used in lwir
								centersXY_aux =  geometryCorrection_aux.getPortsCoordinatesAndDerivatives(
										geometryCorrection_main, //			GeometryCorrection gc_main,
										true,            // boolean use_rig_offsets,
										corr_rots_aux,   // Matrix []   rots,
										null,            //  Matrix [][] deriv_rots,
										null,            // double [][] pXYderiv, // if not null, should be double[8][]
										disp_dist_aux, // used to correct 3D correlations
										centerX,
										centerY,
										disparity_aux); //  + disparity_corr);
							} else {// USED in lwir
								centersXY_aux =  geometryCorrection_aux.getPortsCoordinatesAndDerivatives(
										geometryCorrection_aux, //			GeometryCorrection gc_main,
										false,            // boolean use_rig_offsets,
										corr_rots_aux,   // Matrix []   rots,
										null,            //  Matrix [][] deriv_rots,
										null,            // double [][] pXYderiv, // if not null, should be double[8][]
										disp_dist_aux, // used to correct 3D correlations
										centerX,
										centerY,
										disparity_main); //  + disparity_corr);

							}
						}
						// acquisition time of the tiles centers in scanline times
						if (ers_delay != null) {// not used in lwir
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
							if (calc_both && !isMonochrome() && (globalDebugLevel > -1) && (tileX == debug_tileX) && (tileY == debug_tileY) && (chn == 2)) {
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

							for (int i = 0; i < quad_main; i++) { // quad_main == 0 if not calculated // not used in lwir
								clt_data_main[i][chn] = new double [4][];
								fract_shiftsXY_main[i] = extract_correct_tile( // return a pair of residual offsets
										image_data_main[i],
										width,       // image width
										(clt_kernels_main == null) ? null : clt_kernels_main[i], // [color][tileY][tileX][band][pixel]
										clt_data_main[i][chn], //double  [][]        clt_tile,    // should be double [4][];
										clt_parameters.kernel_step,
//										transform_size,
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
							for (int i = 0; i < quad_aux; i++) { // quad_aux == 0 if not calculated // not used in lwir
								clt_data_aux[i][chn] = new double [4][];
								fract_shiftsXY_aux[i] = extract_correct_tile( // return a pair of residual offsets
										image_data_aux[i],
										width,       // image width
										(clt_kernels_aux == null) ? null : clt_kernels_aux[i], // [color][tileY][tileX][band][pixel]
										clt_data_aux[i][chn], //double  [][]        clt_tile,    // should be double [4][];
										clt_parameters.kernel_step,
//										transform_size,
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
								if (calc_main) {
									double [][] dbg_tile = new double [16][];
									for (int i = 0; i < 16; i++) dbg_tile[i]=clt_data_main[i>>2][chn][i & 3];
									sdfa_instance.showArrays(dbg_tile,  transform_size, transform_size, true, "MAIN_pre-shifted_x"+tileX+"_y"+tileY, titles);
								}
								if (calc_aux) {
									double [][] dbg_tile_aux = new double [16][];
									for (int i = 0; i < 16; i++) dbg_tile_aux[i]=clt_data_aux[i>>2][chn][i & 3];
									sdfa_instance.showArrays(dbg_tile_aux,  transform_size, transform_size, true, "AUX_pre-shifted_x"+tileX+"_y"+tileY, titles);
								}
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
							for (int i = 0; i < quad_main; i++) {// not used in lwir
								fract_shift(    // fractional shift in transform domain. Currently uses sin/cos - change to tables with 2? rotations
										clt_data_main[i][chn], // double  [][]  clt_tile,
//										transform_size,
										fract_shiftsXY_main[i][0],            // double        shiftX,
										fract_shiftsXY_main[i][1],            // double        shiftY,
										((globalDebugLevel > 1) && (chn==0) && (tileX >= debug_tileX - 2) && (tileX <= debug_tileX + 2) &&
												(tileY >= debug_tileY - 2) && (tileY <= debug_tileY+2)));
							}
							for (int i = 0; i < quad_aux; i++) {// USED in lwir
								fract_shift(    // fractional shift in transform domain. Currently uses sin/cos - change to tables with 2? rotations
										clt_data_aux[i][chn], // double  [][]  clt_tile,
//										transform_size,
										fract_shiftsXY_aux[i][0],            // double        shiftX,
										fract_shiftsXY_aux[i][1],            // double        shiftY,
										((globalDebugLevel > 1) && (chn==0) && (tileX >= debug_tileX - 2) && (tileX <= debug_tileX + 2) &&
												(tileY >= debug_tileY - 2) && (tileY <= debug_tileY+2)));
							}

							if ((globalDebugLevel > 0) && (debug_tileX == tileX) && (debug_tileY == tileY)) {
								ShowDoubleFloatArrays sdfa_instance = new ShowDoubleFloatArrays(); // just for debugging?
								String [] titles = {"CC0","SC0","CS0","SS0","CC1","SC1","CS1","SS1","CC2","SC2","CS2","SS2","CC3","SC3","CS3","SS3"};
								if (calc_main) {
									double [][] dbg_tile = new double [16][];
									for (int i = 0; i < 16; i++) dbg_tile[i]=clt_data_main[i>>2][chn][i & 3];
									sdfa_instance.showArrays(dbg_tile,  transform_size, transform_size, true, "MAIN_shifted_x"+tileX+"_y"+tileY+"-z", titles);
								}
								if (calc_aux) {
									double [][] dbg_tile_aux = new double [16][];
									for (int i = 0; i < 16; i++) dbg_tile_aux[i]=clt_data_aux[i>>2][chn][i & 3];
									sdfa_instance.showArrays(dbg_tile_aux,  transform_size, transform_size, true, "AUX_shifted_x"+tileX+"_y"+tileY+"-z", titles);
								}
							}
						}

						int tile_lma_debug_level =  ((tileX == debug_tileX) && (tileY == debug_tileY))? clt_parameters.img_dtt.lma_debug_level : -1;

						// all color channels are done here
						double extra_disparity_main = 0.0; // used for textures:  if allowed, shift images extra before trying to combine
						double extra_disparity_aux = 0.0; // used for textures:  if allowed, shift images extra before trying to combine

						// fill clt_corr_combo if it exists
						if (disparity_bimap != null){ // not null - calculate correlations
							double [] tile_corrs_main = new double[0]; //null;
							double [] tile_corrs_aux  = new double[0]; //null;
							if (calc_main) {// not used in lwir
								tile_corrs_main = tileCorrs(
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
								extra_disparity_main = tile_corrs_main[DISP_FULL_INDEX];
								if (Double.isNaN(extra_disparity_main)) extra_disparity_main = 0;
							}
							if (calc_aux) {// USED in lwir
								tile_corrs_aux  = tileCorrs(
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
								extra_disparity_aux =  tile_corrs_aux[DISP_FULL_INDEX];
								if (Double.isNaN(extra_disparity_aux))  extra_disparity_aux = 0;
							}

							for (int i = 0; i < tile_corrs_main.length; i++) {// not used in lwir
								int dest = SNGL_TO_BI[0][i];
								if (disparity_bimap[dest] != null) disparity_bimap[dest][nTile] = tile_corrs_main[i];
							}
							for (int i = 0; i < tile_corrs_aux.length; i++) {
								int dest = SNGL_TO_BI[1][i];
								if (disparity_bimap[dest] != null) disparity_bimap[dest][nTile] = tile_corrs_aux[i];
							}
							if (calc_main && Double.isNaN(disparity_bimap[BI_STR_FULL_INDEX][tIndex])) {
								System.out.println("BUG: 3. disparity_map[BI_STR_FULL_INDEX][tIndex] should not be NaN");
							}
							if (calc_aux && Double.isNaN(disparity_bimap[BI_ASTR_FULL_INDEX][tIndex])) {
								System.out.println("BUG: 3a. disparity_map[BI_ASTR_FULL_INDEX][tIndex] should not be NaN");
							}

							if (clt_corr_combo != null) { // [type][tilesY][tilesX][(2*transform_size-1)*(2*transform_size-1)] // if null - will not calculate
								tcorr_combo = new double [TCORR_TITLES.length][corr_size * corr_size];
							}
							double [] inter_cam_corr = null;
							double [] inter_corrs_dxy = null;
							if (calc_both) {// not used in lwir
								inter_cam_corr = corr2d.correlateInterCamerasFD(
										clt_data_main,       // double [][][][]     clt_data_tile_main,
										clt_data_aux,        // double [][][][]     clt_data_tile_aux,
										filter,                   // double []           lpf,
										getScaleStrengths(),
										col_weights,              // double []           col_weights,
										fatzero);                 // double              fat_zero)

								inter_corrs_dxy = tileInterCamCorrs(
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
							}
							// finalize ML stuff
							if (ml_data != null) {
								// save data for the main camera
								for (int nlayer = 0; nlayer <  ML_TOP_AUX_INDEX; nlayer++) {
									// save main camera data
									if (calc_main) {// not used in lwir
										corr2d.saveMlTile(
												tileX,                     // int         tileX,
												tileY,                     // int         tileY,
												ml_hwidth,                 // int         ml_hwidth,
												ml_data,                   // double [][] ml_data,
												nlayer + 0,                // int         ml_layer,
												ml_data_main[nlayer],      // double []   ml_tile,
												tilesX);                   // int         tilesX);
									}
									// save aux_camera data
									if (calc_aux) {// USED in lwir
										corr2d.saveMlTile(
												tileX,                     // int         tileX,
												tileY,                     // int         tileY,
												ml_hwidth,                 // int         ml_hwidth,
												ml_data,                   // double [][] ml_data,
												nlayer + ML_TOP_AUX_INDEX, // int         ml_layer,
												ml_data_aux[nlayer],       // double []   ml_tile,
												tilesX);                   // int         tilesX);
									}
								}
								// save inter-camera correlation
								if (calc_both ) {// not used in lwir
									corr2d.saveMlTile(
											tileX,                         // int         tileX,
											tileY,                         // int         tileY,
											ml_hwidth,                     // int         ml_hwidth,
											ml_data,                       // double [][] ml_data,
											ML_INTER_INDEX,                // int         ml_layer,
											ml_data_inter,                 // double []   ml_tile,
											tilesX);                       // int         tilesX);
								}
								// save other data (just 1 value)

								corr2d.saveMlTilePixel(// USED in lwir
							    		tileX,                         // int         tileX,
							    		tileY,                         // int         tileY,
							    		ml_hwidth,                     // int         ml_hwidth,
							    		ml_data,                       // double [][] ml_data,
							    		ML_OTHER_INDEX,                // int         ml_layer,
							    		ML_OTHER_TARGET ,              // int         ml_index,
							    		disparity_main,                // target disparitydouble      ml_value,
							    		tilesX);                       // int         tilesX);

								if (calc_both && (ml_data_dbg1 != null)) { // not used in lwir
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
							for (int i = 0; i < tcorr_combo.length; i++) {// not used in lwir
								clt_corr_combo[i][tileY][tileX] = tcorr_combo[i];
							}
						}

						if (calc_main && (texture_tiles_main !=null)) {// not used in lwir
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
									tile_lma_debug_level,  // final int              debugLevel);
									false); //boolean debug_gpu
						}
						if (calc_aux && (texture_tiles_aux !=null)) {// not used in lwir
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
									tile_lma_debug_level,  // final int              debugLevel);
									false); //boolean debug_gpu
						}
						// Save channel tiles to result
						if (clt_bidata != null) {// not used in lwir
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
		if (calc_both && (lt_corr != null)) {// not used in lwir
			//notch_mode
			// prepare weights for neighbors
			final int lt_rad_x = notch_mode? 0: lt_rad;
			final double [][] neib_weights = new double[lt_rad+1][lt_rad_x+1];
			for (int i = 0; i <= lt_rad; i++) {
				for (int j = 0; j <= lt_rad_x; j++) {
					neib_weights[i][j] = Math.cos(Math.PI * i /(2 * lt_rad + 1)) * Math.cos(Math.PI * j /(2 * lt_rad + 1)); // no need to normalize - it will need to skip empty tiles anyway
				}
			}
			//		final int corr_size = transform_size * 2 -1;

			final TileNeibs tnImage  = new TileNeibs(tilesX, tilesY);
			ai.set(0);
			for (int ithread = 0; ithread < threads.length; ithread++) {
				threads[ithread] = new Thread() {
					@Override
					public void run() {
						double []    ml_data_inter = (ml_data != null)? new double [(2*ml_hwidth +1)*(2*ml_hwidth +1)]:null;
						double [][]  tcorr_combo =     null; // [15*15] pixel space

						Correlation2d corr2d = new Correlation2d(
								numSensors,
								clt_parameters.img_dtt,              // ImageDttParameters  imgdtt_params,
								transform_size,             // int transform_size,
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

	public void  clt_bi_macro( // not yet operational // not used in lwir
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
		final int tilesX=width/transform_size;
		final int tilesY=height/transform_size;
		final int nTilesInChn=tilesX*tilesY;

		// clt_data does not need to be for the whole image (no, it is used for textures)

		final Thread[] threads = newThreadArray(threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		final int corr_size = transform_size * 2 -1;
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
		int wcenter = transform_size - 1;
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
			System.out.println("clt_aberrations_quad_corr(): width="+width+" height="+height+" transform_size="+transform_size+
					" debug_tileX="+debug_tileX+" debug_tileY="+debug_tileY+" globalDebugLevel="+globalDebugLevel);
		}
		final double [] filter =  doubleGetCltLpfFd(clt_parameters.getCorrSigma(isMonochrome()));

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

		DttRad2 dtt = new DttRad2(transform_size);
		dtt.set_window(clt_parameters.clt_window);
		final double [] lt_window = dtt.getWin2d();	// [256]
		final double [] lt_window2 = new double [lt_window.length]; // squared
		for (int i = 0; i < lt_window.length; i++) lt_window2[i] = lt_window[i] * lt_window[i];


		if (globalDebugLevel > 1) {
			ShowDoubleFloatArrays sdfa_instance = new ShowDoubleFloatArrays(); // just for debugging?
			sdfa_instance.showArrays(lt_window,  2*transform_size, 2*transform_size, "lt_window");
		}

		//FIXME: no rotation matrices, no rig disparity - it is all included in source data in macro mode. Distortion?

//		final Matrix rigMatrix = geometryCorrection_aux.getRotMatrix(true);
		final int nChn = image_rig_data[0].length; // 1; ==1
		final double [] col_weights= {1.0};
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				@Override
				public void run() {
					DttRad2 dtt = new DttRad2(transform_size);
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
							numSensors,
							clt_parameters.img_dtt,              // ImageDttParameters  imgdtt_params,
							transform_size,             // int transform_size,
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
						centerX = tileX * transform_size + transform_size/2; //  - shiftX;
						centerY = tileY * transform_size + transform_size/2; //  - shiftY;
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
//									transform_size,
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
//									transform_size,
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
//									transform_size,
									fract_shiftsXY_main[0],            // double        shiftX,
									fract_shiftsXY_main[1],            // double        shiftY,
									((globalDebugLevel > 1) && (chn==0) && (tileX >= debug_tileX - 2) && (tileX <= debug_tileX + 2) &&
											(tileY >= debug_tileY - 2) && (tileY <= debug_tileY+2)));
							fract_shift(    // fractional shift in transform domain. Currently uses sin/cos - change to tables with 2? rotations
									clt_data_aux[0][chn], // double  [][]  clt_tile,
//									transform_size,
									fract_shiftsXY_aux[0],            // double        shiftX,
									fract_shiftsXY_aux[1],            // double        shiftY,
									((globalDebugLevel > 1) && (chn==0) && (tileX >= debug_tileX - 2) && (tileX <= debug_tileX + 2) &&
											(tileY >= debug_tileY - 2) && (tileY <= debug_tileY+2)));



							if ((globalDebugLevel > 100) && (debug_tileX == tileX) && (debug_tileY == tileY)) {
								ShowDoubleFloatArrays sdfa_instance = new ShowDoubleFloatArrays(); // just for debugging?
								String [] titles = {"CC0","SC0","CS0","SS0"};
								double [][] dbg_tile = new double [4][];
								for (int i = 0; i < 4; i++) dbg_tile[i]=clt_data_main[0][chn][i & 3];
								sdfa_instance.showArrays(dbg_tile,  transform_size, transform_size, true, "MAIN_shifted_x"+tileX+"_y"+tileY+"-z", titles);
								double [][] dbg_tile_aux = new double [16][];
								for (int i = 0; i < 4; i++) dbg_tile[i]=clt_data_aux[0][chn][i & 3];
								sdfa_instance.showArrays(dbg_tile_aux,  transform_size, transform_size, true, "AUX_shifted_x"+tileX+"_y"+tileY+"-z", titles);
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
	
	@Deprecated
	public double [][][][][][] clt_aberrations_quad_corr( // USED in LWIR
			final ImageDttParameters  imgdtt_params,   // Now just extra correlation parameters, later will include, most others
			final int                 macro_scale,     // to correlate tile data instead of the pixel data: 1 - pixels, 8 - tiles
			final int [][]            tile_op,         // [tilesY][tilesX] - what to do - 0 - nothing for this tile
			final double [][]         disparity_array, // [tilesY][tilesX] - individual per-tile expected disparity
			final double [][][]       image_data, // first index - number of image in a quad
		    final boolean [][]        saturation_imp, // (near) saturated pixels or null
			 // correlation results - final and partial
			final double [][][][]     clt_corr_combo_in,  // [type][tilesY][tilesX][(2*transform_size-1)*(2*transform_size-1)] // if null - will not calculate
			                                           // [type][tilesY][tilesX] should be set by caller
													   // types: 0 - selected correlation (product+offset), 1 - sum
			final double [][][][][]   clt_corr_partial,// [tilesY][tilesX][quad]color][(2*transform_size-1)*(2*transform_size-1)] // if null - will not calculate
                                                       // [tilesY][tilesX] should be set by caller
			// When clt_mismatch is non-zero, no far objects extraction will be attempted
			// clt_mismatch is used in older code, not supported in GPU - there is cltMeasureLazyEye for that purpose
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
	  		final double              min_corr,        // 0.02; // minimal correlation value to consider valid - not used
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
//			final int                 transform_size,
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
		// ****** FIXME tries to use color == 3, should be disabled!
		if (clt_corr_combo_in != null) {
			System.out.println("clt_corr_combo != null");
			System.out.println("clt_corr_combo != null");
			System.out.println("clt_corr_combo != null");
		}
		final double [][][][]     clt_corr_combo = (globalDebugLevel >-10000)?null:clt_corr_combo_in;
		
		final boolean debug_distort= globalDebugLevel > 0; ///false; // true;

		final double [][] debug_offsets = null; // new double[imgdtt_params.lma_dbg_offset.length][2];
		/*
		final double [][] debug_offsets = new double[getNumSensors()][2];  
		for (int i = 0; i < imgdtt_params.lma_dbg_offset.length; i++) for (int j = 0; j < debug_offsets[i].length; j++) {
			debug_offsets[i][j] = imgdtt_params.lma_dbg_offset[i][j]*imgdtt_params.lma_dbg_scale;
		}
		*/
		final double [] dbg_corr_shift =((imgdtt_params.pcorr_dbg_offsx != 0.0) || (imgdtt_params.pcorr_dbg_offsy != 0.0))?
				(new double [] {imgdtt_params.pcorr_dbg_offsx,imgdtt_params.pcorr_dbg_offsy}):null;

		final boolean macro_mode = macro_scale != 1;      // correlate tile data instead of the pixel data
		final int numcol = 3; // number of colors // keep the same, just do not use [0] and [1], [2] - green
		final int height=image_data[0][0].length/width;
		final int tilesX=width/transform_size;
		final int tilesY=height/transform_size;
		final int nTilesInChn=tilesX*tilesY;
		final double [][][][][][] clt_data = new double[numSensors][numcol][tilesY][tilesX][][];
		final Thread[] threads = newThreadArray(threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		final double [] col_weights= new double [numcol]; // colors are RBG
		final double [][] dbg_distort = debug_distort? (new double [4*numSensors][tilesX*tilesY]) : null;
		final double [][] corr_wnd = Corr2dLMA.getCorrWnd(
				transform_size,
				imgdtt_params.lma_wnd);
		final double [] corr_wnd_inv_limited = (imgdtt_params.lma_min_wnd <= 1.0)?  new double [corr_wnd.length * corr_wnd[0].length]: null;
		if (corr_wnd_inv_limited != null) {
			double inv_pwr = imgdtt_params.lma_wnd_pwr - (imgdtt_params.lma_wnd - 1.0); // compensate for lma_wnd
			for (int i = imgdtt_params.lma_hard_marg; i < (corr_wnd.length - imgdtt_params.lma_hard_marg); i++) {
				for (int j = imgdtt_params.lma_hard_marg; j < (corr_wnd.length - imgdtt_params.lma_hard_marg); j++) {
					corr_wnd_inv_limited[i * (corr_wnd.length) + j] = 1.0/Math.max(Math.pow(corr_wnd[i][j],inv_pwr), imgdtt_params.lma_min_wnd);
				}
			}
		}

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
		final int transform_len = transform_size * transform_size;

		final double [] filter =     doubleGetCltLpfFd(corr_sigma);
		
		final double [] filter_rb =      doubleGetCltLpfFd(imgdtt_params.pcorr_sigma_rb);
		final double [] filter_common =  doubleGetCltLpfFd(imgdtt_params.getCorrSigma(isMonochrome()));
		final double pcorr_fat_zero =    imgdtt_params.getFatZero(isMonochrome());

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
		/*
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
		*/
		if (disparity_map != null){
			for (int i = 0; i<disparity_map.length;i++){
				if (i == OVEREXPOSED) {
					if (saturation_imp!= null) {
						disparity_map[i] = new double [tilesY*tilesX];
					}
				} else if (isSliceBit(i)) {
					disparity_map[i] = new double [tilesY*tilesX];
				} else if (isDiffIndex(i)) {
					disparity_map[i] = new double [tilesY*tilesX];
				} else if (isToneRGBIndex(i)) {
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

		DttRad2 dtt = new DttRad2(transform_size);
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
					double [][] fract_shiftsXY = new double[numSensors][];
					double [][]     tcorr_combo =    null; // [15*15] pixel space
					double [][][]   tcorr_partial =  null; // [quad][numcol+1][15*15]
					double [][][][] tcorr_tpartial = null; // [quad][numcol+1][4][8*8]
					double [] ports_rgb = null;
					double [][] rXY;

					if (use_main) {
						rXY = geometryCorrection.getRXY(true); // boolean use_rig_offsets,
					} else  {
						rXY = geometryCorrection.getRXY(false); // boolean use_rig_offsets,
					}
					Correlation2d corr2d = new Correlation2d(
							numSensors,
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
						for (int i = 0; i< CORR_PAIRS.length; i++){
							if ((((1 << CORR_PAIRS[i][0]) & img_mask) == 0) || (((1 << CORR_PAIRS[i][1]) & img_mask) == 0)) {
								corr_mask &= ~ (1 << i);
							}
						}
						boolean debugTile =(tileX == debug_tileX) && (tileY == debug_tileY) && (globalDebugLevel > -1);
						boolean debugTile0 =(tileX == debug_tileX) && (tileY == debug_tileY) && (globalDebugLevel > -3);
						final int [] overexp_all = (saturation_imp != null) ? ( new int [2]): null;

						// Moved from inside chn loop
						centerX = tileX * transform_size + transform_size/2 - shiftX;
						centerY = tileY * transform_size + transform_size/2 - shiftY;
						// TODO: move port coordinates out of color channel loop
						double [][] centersXY;
						double [][] disp_dist = new double[numSensors][]; // used to correct 3D correlations

						if ((disparity_array == null) || (disparity_array[tileY] == null) || (Double.isNaN(disparity_array[tileY][tileX]))) {
							System.out.println("Bug with disparity_array !!! - 5");
							continue; // nothing to do for this tile
						}

						if (macro_mode){
							if ((globalDebugLevel > -1) && debugTile) { // before correction
								System.out.println("\nUsing MACRO mode, centerX="+centerX+", centerY="+centerY);
							}
							centersXY = geometryCorrection.getPortsCoordinatesIdeal(
									macro_scale,
									centerX,
									centerY,
									macro_scale* disparity_array[tileY][tileX] + disparity_corr);

						} else {
							if (use_main) { // this is AUX camera that uses main coordinates  // not used in lwir
								centersXY =  geometryCorrection.getPortsCoordinatesAndDerivatives(
										geometryCorrection_main, //			GeometryCorrection gc_main,
										true,            // boolean use_rig_offsets,
										corr_rots,       // Matrix []   rots,
										null,            //  Matrix [][] deriv_rots,
										null,            // double [][] pXYderiv, // if not null, should be double[8][]
										disp_dist,       // used to correct 3D correlations
										centerX,
										centerY,
										disparity_array[tileY][tileX] + disparity_corr); // _aux); //  + disparity_corr);


							} else {  // used in lwir
								centersXY = geometryCorrection.getPortsCoordinatesAndDerivatives(
										geometryCorrection, //			GeometryCorrection gc_main,
										false,           // boolean use_rig_offsets,
										corr_rots,       // Matrix []   rots,
										null,            //  Matrix [][] deriv_rots,
										null,            // double [][] pXYderiv, // if not null, should be double[8][]
										disp_dist,       // used to correct 3D correlations
										centerX,
										centerY,
										disparity_array[tileY][tileX] + disparity_corr);
							}

							if (((globalDebugLevel > 0) || debug_distort) && debugTile) {
								for (int i = 0; i < numSensors; i++) {
									System.out.println("clt_aberrations_quad_corr():  tileX="+tileX+", tileY="+tileY+
											" centerX="+centerX+" centerY="+centerY+" disparity="+disparity_array[tileY][tileX]+
											" centersXY["+i+"][0]="+centersXY[i][0]+" centersXY["+i+"][1]="+centersXY[i][1]);
								}
							}
							if (debug_distort && debugTile && (debug_offsets != null)) {
								double [][] debug_offsets_xy = new double [debug_offsets.length][2];
								for (int i = 0; i < debug_offsets.length; i++) {
									debug_offsets_xy[i][0] = disp_dist[i][0] * debug_offsets[i][0] + disp_dist[i][1] * debug_offsets[i][1];
									debug_offsets_xy[i][1] = disp_dist[i][2] * debug_offsets[i][0] + disp_dist[i][3] * debug_offsets[i][1];
								}
								for (int i = 0; i < numSensors; i++) {
									System.out.println(String.format("%d: {%8.3f, %8.3f}",i,debug_offsets_xy[i][0],debug_offsets_xy[i][1]));
								}

								for (int i = 0; i < debug_offsets.length; i++) {
									centersXY[i][0] += debug_offsets_xy[i][0];
									centersXY[i][1] += debug_offsets_xy[i][1];
								}
								for (int i = 0; i < numSensors; i++) {
									System.out.println("Corrected clt_aberrations_quad_corr():  tileX="+tileX+", tileY="+tileY+
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

							for (int ip = 0; ip < centersXY.length; ip++){ // 4 OOB
								centersXY[ip][0] -= shiftXY[ip][0];
								centersXY[ip][1] -= shiftXY[ip][1];
							}
							// save disparity distortions for visualization:
							if (dbg_distort != null) {
								for (int cam = 0; cam <numSensors; cam++) {
									dbg_distort[cam * 4 + 0 ][nTile] = disp_dist[cam][0];
									dbg_distort[cam * 4 + 1 ][nTile] = disp_dist[cam][1];
									dbg_distort[cam * 4 + 2 ][nTile] = disp_dist[cam][2];
									dbg_distort[cam * 4 + 3 ][nTile] = disp_dist[cam][3];
								}
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
						if (FPGA_COMPARE_DATA && (globalDebugLevel > 0) && (tileX == debug_tileX) && (tileY == debug_tileY)) { // not used in lwir
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
//									transform_size,       // final int                 transform_size,
									width,                 // final int                 width
									dtt
									);
						}
// See if macro_mode uses color channels for non-color?
						for (int ncol = 0; ncol <numcol; ncol++) {
							if (!isMonochrome() || (ncol == MONO_CHN) || macro_mode) { // in monochrome mode skip all non-mono (green) channels  // used in lwir (5 of 6 branches)
								boolean debug_for_fpga = FPGA_COMPARE_DATA && (globalDebugLevel > 0) && (tileX == debug_tileX) && (tileY == debug_tileY) && (ncol == 2);
								if ((globalDebugLevel > -1) && (tileX == debug_tileX) && (tileY == debug_tileY) && (ncol == 2)) {
									System.out.println("\nUsing "+(macro_mode?"MACRO":"PIXEL")+" mode, centerX="+centerX+", centerY="+centerY);
									System.out.println(disparity_array[tileY][tileX]+"\t"+
											centersXY[0][0]+"\t"+centersXY[0][1]+"\t"+
											centersXY[1][0]+"\t"+centersXY[1][1]+"\t"+
											centersXY[2][0]+"\t"+centersXY[2][1]+"\t"+
											centersXY[3][0]+"\t"+centersXY[3][1]+"\t");
								}

								for (int i = 0; i < numSensors; i++) {
									if (debug_for_fpga && (i==0)){ // not used in lwir
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
//												transform_size,
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
//												transform_size,
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
									// Extract image tiles and kernels, correct aberrations, return (ut do not apply) fractional shifts
									fract_shiftsXY[i] = extract_correct_tile( // return a pair of residual offsets
											image_data[i],
											width,       // image width
											((clt_kernels == null) ? null : clt_kernels[i]), // [color][tileY][tileX][band][pixel]
											clt_data[i][ncol][tileY][tileX], //double  [][]        clt_tile,    // should be double [4][];
											kernel_step,
//											transform_size,
											dtt,
											ncol,
											centersXY[i][0], // centerX, // center of aberration-corrected (common model) tile, X
											centersXY[i][1], // centerY, //
											((!FPGA_COMPARE_DATA && (globalDebugLevel > -1) && (tileX == debug_tileX) && (tileY == debug_tileY) && (ncol == 2) && (i==0)) ? (globalDebugLevel + 0) : 0), // external tile compare
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
									for (int i = 0; i < numSensors; i++) {
										System.out.println("clt_aberrations_quad(): color="+ncol+", tileX="+tileX+", tileY="+tileY+
												" fract_shiftsXY["+i+"][0]="+fract_shiftsXY[i][0]+" fract_shiftsXY["+i+"][1]="+fract_shiftsXY[i][1]);
									}
								}

								if (!no_fract_shift) {  // USED in lwir
									// apply residual shift
									for (int i = 0; i < numSensors; i++) {
										fract_shift(    // fractional shift in transform domain. Currently uses sin/cos - change to tables with 2? rotations
												clt_data[i][ncol][tileY][tileX], // double  [][]  clt_tile,
//												transform_size,
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
								for (int i = 0; i < numSensors; i++) {  // used in lwir
									clt_data[i][ncol] = null; // erase unused clt_data
								}
							}
						}// end of for (int chn = 0; chn <numcol; chn++)
						// used in lwir
//						int tile_lma_debug_level =  ((tileX == debug_tileX) && (tileY == debug_tileY))? imgdtt_params.lma_debug_level : -1;
						int tile_lma_debug_level =  ((tileX == debug_tileX) && (tileY == debug_tileY))? (imgdtt_params.lma_debug_level-1) : -2;
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

							double [][]  corrs6 = corr2d.correlateCompositeFD( // now works with nulls for some clt_data colors
						    		clt_data,        // double [][][][][][] clt_data,
						    		tileX,           // int                 tileX,
						    		tileY,           // int                 tileY,
						    		all_pairs,       // int                 pairs_mask,
						    		filter,          // double []           lpf,
						    		getScaleStrengths(), // double              scale_value, // scale correlation value
						    		col_weights,     // double []           col_weights,
						    		corr_fat_zero);  // double              fat_zero)
							double [][]  corrs = {corrs6[0],corrs6[1],corrs6[2],corrs6[3],corrs6[4],corrs6[5],null,null,null,null};

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
//						    if ((clt_corr_partial != null) && imgdtt_params.corr_mode_debug) {
//						    if ((clt_corr_partial != null) && (imgdtt_params.corr_mode_debug || imgdtt_params.gpu_mode_debug)) {
						    if (debugTile0) {
						    	System.out.println("Debugging tileY = "+tileY+"  tileX = " + tileX);
//						    	System.out.println("Debugging tileY = "+tileY+"  tileX = " + tileX);
//						    	System.out.println("Debugging tileY = "+tileY+"  tileX = " + tileX);
						    }

						    double [][][] corr_pairs_td = corr2d.correlateCompositeTD(
						    		clt_data,            // double [][][][][][] clt_data,
						    		tileX,               // int                 tileX,
						    		tileY,               // int                 tileY,
						    		0x3f,                // int                 pairs_mask,
						    		filter_rb,           // double []           lpf_rb, // extra lpf for red and blue (unused for mono) or null
						    		getScaleStrengths(), // double              scale_value, // scale correlation value
						    		col_weights);        // double []           col_weights);

						    double [][] corrs_ortho_td = corr2d.combineCorrelationsTD(
						    		corr_pairs_td,       // double [][][] corr_combo_td, // vertical will be transposed, other diagonal flipped vertically (should be separate)
						    		0x0f);               // int           pairs_mask); // vertical
						    if (dbg_corr_shift!= null) {
						    	fract_shift(    // fractional shift in transform domain. Currently uses sin/cos - change to tables with 2? rotations // USED in lwir
						    			corrs_ortho_td, // double  [][]  clt_tile,
						    			dbg_corr_shift[0], // double        shiftX,
						    			dbg_corr_shift[1], // double        shiftY,
						    			false); // boolean       bdebug);
						    }

						    corr2d.normalize_TD(
						    		corrs_ortho_td,      // double [][] td,
						    		filter_common,       // double []   lpf, // or null
						    		pcorr_fat_zero);     // double      fat_zero);

						    //						    	double [] corrs_ortho = corr2d.convertCorrToPD(
						    corrs[6] = corr2d.convertCorrToPD(
						    		corrs_ortho_td);     // double [][] td);

						    double [][] corrs_cross_td = corr2d.combineCorrelationsTD(
						    		corr_pairs_td,       // double [][][] corr_combo_td, // vertical will be transposed, other diagonal flipped vertically (should be separate)
						    		0x30);               // int           pairs_mask); // vertical
						    if (dbg_corr_shift!= null) {
						    	fract_shift(    // fractional shift in transform domain. Currently uses sin/cos - change to tables with 2? rotations // USED in lwir
						    			corrs_cross_td, // double  [][]  clt_tile,
						    			dbg_corr_shift[0], // double        shiftX,
						    			dbg_corr_shift[1], // double        shiftY,
						    			false); // boolean       bdebug);
						    }

						    corr2d.normalize_TD(
						    		corrs_cross_td,      // double [][] td,
						    		filter_common,       // double []   lpf, // or null
						    		pcorr_fat_zero);     // double      fat_zero);

						    //						    	double [] corrs_cross = corr2d.convertCorrToPD(
						    corrs[7] = corr2d.convertCorrToPD(
						    		corrs_cross_td);     // double [][] td);
						    
						    double [][] corrs_hor_td = corr2d.combineCorrelationsTD(
						    		corr_pairs_td,       // double [][][] corr_combo_td, // vertical will be transposed, other diagonal flipped vertically (should be separate)
						    		0x03);               // int           pairs_mask); // horizontal
						    corr2d.normalize_TD(
						    		corrs_hor_td,      // double [][] td,
						    		filter_common,       // double []   lpf, // or null
						    		pcorr_fat_zero);     // double      fat_zero);
						    corrs[8] = corr2d.convertCorrToPD(
						    		corrs_hor_td);     // double [][] td);
						    

						    double [][] corrs_vert_td = corr2d.combineCorrelationsTD(
						    		corr_pairs_td,       // double [][][] corr_combo_td, // vertical will be transposed, other diagonal flipped vertically (should be separate)
						    		0x0c,                // int           pairs_mask); // vertical
						    		true);               // boolean       no_transpose_vertical
						    corr2d.normalize_TD(
						    		corrs_vert_td,      // double [][] td,
						    		filter_common,       // double []   lpf, // or null
						    		pcorr_fat_zero);     // double      fat_zero);
						    corrs[9] = corr2d.convertCorrToPD( // ** Now - not transposed anymore!
						    		corrs_vert_td);     // double [][] td);
						    
						    
						    
						    double [][] strips_intra = new double[2][];
						    strips_intra[0] = corr2d.scaleRotateInterpoateSingleCorrelation(
						    		corrs[6], // double []   corr,  // quad 
						    		imgdtt_params.corr_strip_hight, // int         hwidth,
						    		Correlation2d.PAIR_HORIZONTAL,  // int         dir, // 0 - hor, 1 - vert, 2 - parallel to row = col (main) diagonal (0->3), 3 -2->1
						    		1,                              // int         ss, // 1
						    		false // boolean     debug
						    		); 
						    strips_intra[1] = corr2d.scaleRotateInterpoateSingleCorrelation(
						    		corrs[7], // double []   corr,  // quad 
						    		imgdtt_params.corr_strip_hight, // int         hwidth,
						    		Correlation2d.PAIR_DIAGONAL_MAIN,  // int         dir, // 0 - hor, 1 - vert, 2 - parallel to row = col (main) diagonal (0->3), 3 -2->1
						    		1,                              // int         ss, // 1
						    		false // boolean     debug
						    		); 
						    double [] strip_combo_intra = corr2d.combineInterpolatedCorrelations(
						    		strips_intra,                  // double [][] strips,
						    		3,                     // int         pairs_mask,
						    		imgdtt_params.corr_offset,     // double      offset);
						    		imgdtt_params.twice_diagonal); //    		boolean     twice_diagonal)
						    
						    
						    if ((clt_corr_partial != null) && (imgdtt_params.corr_mode_debug || imgdtt_params.gpu_mode_debug)) {

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
							    clt_corr_partial[tileY][tileX] = new double[numSensors][numcol+1][];
						    	clt_corr_partial[tileY][tileX][0][0] = corrs[0];                        // 1
						    	clt_corr_partial[tileY][tileX][0][1] = corrs[1];                        // 2
						    	clt_corr_partial[tileY][tileX][0][2] = corrs[2];                        // 3
						    	clt_corr_partial[tileY][tileX][0][3] = corrs[3];                        // 4
						    	clt_corr_partial[tileY][tileX][1][0] = corrs[4];                        // 5
						    	clt_corr_partial[tileY][tileX][1][1] = corrs[5];                        // 6
						    	clt_corr_partial[tileY][tileX][1][2] = corrs[6];                        // 7
						    	clt_corr_partial[tileY][tileX][1][3] = corrs[7];                        // 8
//						    	clt_corr_partial[tileY][tileX][1][2] = corrs_ortho;                     // 7
//						    	clt_corr_partial[tileY][tileX][1][3] = corrs_cross;                     // 8
//						    	clt_corr_partial[tileY][tileX][1][2] = corr2d.debugStrip(strip_hor);    // 7
//						    	clt_corr_partial[tileY][tileX][1][3] = corr2d.debugStrip(strip_vert);   // 8
//strip_combo_intra						    	
						    	clt_corr_partial[tileY][tileX][2][0] = corrs[8];                        // 9
						    	clt_corr_partial[tileY][tileX][2][1] = corrs[9];                        // 10
//						    	clt_corr_partial[tileY][tileX][2][0] = corr2d.debugStrip(strips[4]);    // 9
//						    	clt_corr_partial[tileY][tileX][2][1] = corr2d.debugStrip(strips[5]);    // 10
						    	clt_corr_partial[tileY][tileX][2][2] = corr2d.debugStrip2(strip_hor);   // 11
						    	clt_corr_partial[tileY][tileX][2][3] = corr2d.debugStrip2(strip_vert);  // 12
						    	clt_corr_partial[tileY][tileX][3][0] = corr2d.debugStrip2(strips_intra[0]); // 13
						    	clt_corr_partial[tileY][tileX][3][1] = corr2d.debugStrip2(strips_intra[1]);  // 14
//						    	clt_corr_partial[tileY][tileX][3][0] = corr2d.debugStrip2(strip_ortho); // 13
//						    	clt_corr_partial[tileY][tileX][3][1] = corr2d.debugStrip2(strip_diag);  // 14
						    	clt_corr_partial[tileY][tileX][3][2] = corr2d.debugStrip2(strip_combo_intra);    // 15
//						    	clt_corr_partial[tileY][tileX][3][2] = corr2d.debugStrip(strip_all);    // 15
						    	clt_corr_partial[tileY][tileX][3][3] = corr2d.debugStrip2(strip_combo); // 16
						    }
						    if (imgdtt_params.pcorr_use) {
						    	strip_combo = strip_combo_intra;
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
//							double strength = 0.0;
//							double disparity = 0.0;
							double [] disp_str = new double[2];
							if (ixy != null) { //TODO - for CM use magic!
								disp_str[1] = strip_combo[ixy[0]+transform_size-1]; // strength at integer max on axis
								disparity_map[DISPARITY_INDEX_INT][tIndex] =      -ixy[0];
//								disparity_map[DISPARITY_INDEX_INT + 1][tIndex] =
								disparity_map[DISPARITY_STRENGTH_INDEX][tIndex] = disp_str[1];
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
							
							if (imgdtt_params.pcorr_use_hv) { // use combined horizontal and vertical pairs (vertical not transposed) 
								// for compatibility with old code executed unconditionally. TODO: Move to if (corr_stat != null) ... condition below
								double [] hor_pair1 = corr2d.getMaxXSOrtho(
										corrs,                              // double [][] correlations,
										0x100, //  corrs[8] Correlation2d.getMaskHorizontal(1), // int         pairs_mask,
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
										0x200, // corrs[9] Correlation2d.getMaskVertical(1), // int         pairs_mask,
										imgdtt_params.corr_offset,        // double      corr_offset,
										true,                             // boolean     symmetric,   // for comparing with old implementation average with symmetrical before multiplication
										// change to true, un-rotate source
										true, // false, // already transposed  // true,                             // boolean     is_vert,      // transpose X/Y
										tile_lma_debug_level > 0); // boolean   debug);
								if (vert_pair1 != null) {
									disparity_map[DISPARITY_INDEX_VERT][tIndex] =         -vert_pair1[0];
									disparity_map[DISPARITY_INDEX_VERT_STRENGTH][tIndex] = vert_pair1[1];
								}
							} else {
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
							}


							// proceed only if CM correlation result is non-null // for compatibility with old code we need it to run regardless of the strength of the normal correlation
							if (corr_stat != null) {
// skipping DISPARITY_VARIATIONS_INDEX - it was not used
								disp_str[0] = -corr_stat[0];
								disparity_map[DISPARITY_INDEX_CM][tIndex] = disp_str[0]; // disparity is negative X
								if (tile_lma_debug_level > 0) {
									System.out.println("Will run getMaxXSOrtho( ) for tileX="+tileX+", tileY="+tileY);
								}

								// debug new LMA correlations
								if (debugTile) {
									System.out.println("Will run new LMA for tileX="+tileX+", tileY="+tileY);

									double [] poly_disp = {Double.NaN, 0.0};
							    	Corr2dLMA lma2 = corr2d.corrLMA2Single(
							    			imgdtt_params,                // ImageDttParameters  imgdtt_params,
								    		false,                        // boolean             adjust_ly, // adjust Lazy Eye
							    			corr_wnd,                     // double [][]         corr_wnd, // correlation window to save on re-calculation of the window
							    			corr_wnd_inv_limited,         // corr_wnd_limited, // correlation window, limited not to be smaller than threshold - used for finding max/convex areas (or null)
							    			corrs,                        // double [][]         corrs,
							    			disp_dist,
											rXY,                          // double [][]         rXY, // non-distorted X,Y offset per nominal pixel of disparity
											corr2d.longToArray(imgdtt_params.dbg_pair_mask), // imgdtt_params.dbg_pair_mask  // int                 pair_mask, // which pairs to process
							    			disp_str,  //corr_stat[0],                 // double    xcenter,   // preliminary center x in pixels for largest baseline
											poly_disp,                    // double[]            poly_ds,    // null or pair of disparity/strength
							    			imgdtt_params.ortho_vasw_pwr, // double    vasw_pwr,  // value as weight to this power,
							    			tile_lma_debug_level+2,         // int                 debug_level,
							        		tileX,                        // int                 tileX, // just for debug output
							        		tileY );                      // int                 tileY
							    	double [][] ds = null;
							    	if (lma2 != null) {
							    		ds = lma2.lmaDisparityStrength(
							    				imgdtt_params.lmas_min_amp,      //  minimal ratio of minimal pair correlation amplitude to maximal pair correlation amplitude
							    				imgdtt_params.lmas_max_rel_rms,  // maximal relative (to average max/min amplitude LMA RMS) // May be up to 0.3)
							    				imgdtt_params.lmas_min_strength, // minimal composite strength (sqrt(average amp squared over absolute RMS)
							    				imgdtt_params.lmas_min_ac,       // minimal of A and C coefficients maximum (measures sharpest point/line)
												imgdtt_params.lmas_min_min_ac,   // minimal of A and C coefficients minimum (measures sharpest point)
							    				imgdtt_params.lmas_max_area,      //double  lma_max_area,     // maximal half-area (if > 0.0)
							    				imgdtt_params.lma_str_scale,    // convert lma-generated strength to match previous ones - scale
							    				imgdtt_params.lma_str_offset    // convert lma-generated strength to match previous ones - add to result
							    				);
							    		lma2.printStats(ds,1);
							    	}
								}

//								disparity_map[DISPARITY_INDEX_CM + 1][tIndex] = // y not available here
								// calculate/fill out hor and vert
								// convert to multi-baseline combining results from several integer scales

								// see if strength is enough to proceed with LMA/poly (otherwise keep disp/strength
								if (disp_str[1] > imgdtt_params.min_poly_strength) {
									// create LMA instance, calculate LMA composite argmax
							    	// Create 2 groups: ortho & diag
									Correlations2dLMA lma;
									if (imgdtt_params.pcorr_use && !imgdtt_params.pcorr_use_hv) { // incorrect as the ellipsoid is rotated!
										double [][] fake_corrs = {corrs[6],null,null,null,corrs[7],null}; // 
										lma = corr2d.corrLMA(
								    			imgdtt_params,                // ImageDttParameters  imgdtt_params,
								    			fake_corrs,                   // double [][]         corrs,
								    			corr2d.longToArray(0x11),     // 0x11, // imgdtt_params.dbg_pair_mask,  // int                 pair_mask, // which pairs to process
								        		false,                        // boolean             run_poly_instead, // true - run LMA, false - run 2d polynomial approximation
								    			corr_stat[0],                 // double    xcenter,   // preliminary center x in pixels for largest baseline
								    			imgdtt_params.ortho_vasw_pwr, // double    vasw_pwr,  // value as weight to this power,
								    			tile_lma_debug_level,         // int                 debug_level,
								        		tileX,                        // int                 tileX, // just for debug output
								        		tileY );                      // int                 tileY
									} else 	if (imgdtt_params.pcorr_use_hv) { // correct way to combine
										// use hor pairs as one (top), use vert pairs as one (left)
										double [][] fake_corrs = {corrs[8],null,corrs[9],null,corrs[4],corrs[4]}; // hor, null, vert, diagm, diago
										lma = corr2d.corrLMA(
								    			imgdtt_params,                // ImageDttParameters  imgdtt_params,
								    			fake_corrs,                   // double [][]         corrs,
								    			corr2d.longToArray(0x11),     // 0x11, // imgdtt_params.dbg_pair_mask,  // int                 pair_mask, // which pairs to process
								        		false,                        // boolean             run_poly_instead, // true - run LMA, false - run 2d polynomial approximation
								    			corr_stat[0],                 // double    xcenter,   // preliminary center x in pixels for largest baseline
								    			imgdtt_params.ortho_vasw_pwr, // double    vasw_pwr,  // value as weight to this power,
								    			tile_lma_debug_level,         // int                 debug_level,
								        		tileX,                        // int                 tileX, // just for debug output
								        		tileY );                      // int                 tileY
 
										
									} else { //  old way, 6 pairs top, bottom, left, right, diagm, diago
										lma = corr2d.corrLMA(
								    			imgdtt_params,                // ImageDttParameters  imgdtt_params,
								    			corrs,                        // double [][]         corrs,
								    			corr2d.longToArray(imgdtt_params.dbg_pair_mask),     // imgdtt_params.dbg_pair_mask,  // int                 pair_mask, // which pairs to process
								        		false,                        // boolean             run_poly_instead, // true - run LMA, false - run 2d polynomial approximation
								    			corr_stat[0],                 // double    xcenter,   // preliminary center x in pixels for largest baseline
								    			imgdtt_params.ortho_vasw_pwr, // double    vasw_pwr,  // value as weight to this power,
								    			tile_lma_debug_level,         // int                 debug_level,
								        		tileX,                        // int                 tileX, // just for debug output
								        		tileY );                      // int                 tileY
									}
									
							    	
							    	double [] lma_disparity_strength = null;
							    	double max_disp_diff_lma = 3.0;
							    	if (lma != null) {
							    		lma_disparity_strength = lma.getDisparityStrength();
							    		if (Math.abs(lma_disparity_strength[0] - disp_str[0] ) > max_disp_diff_lma) {
							    			if (globalDebugLevel > -1) {
							    				System.out.println("Crazy LMA for tileX="+tileX+", tileY="+tileY+": disparity="+disp_str[0]+",lma_disparity_strength[0]="+lma_disparity_strength[0]);
							    			}
							    			lma = null;
							    		}
							    	}
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
										if (imgdtt_params.mix_corr_poly) { //true
											disp_str[0] =  lma_disparity_strength[0];
											disp_str[1] =  lma_disparity_strength[1];
											disparity_map[DISPARITY_INDEX_CM]       [tIndex] = disp_str[0];
											disparity_map[DISPARITY_STRENGTH_INDEX] [tIndex] = disp_str[1];
											if (Double.isNaN(disparity_map[DISPARITY_STRENGTH_INDEX][tIndex])) {
												System.out.println("BUG: 2. disparity_map[DISPARITY_STRENGTH_INDEX][tIndex] should not be NaN");
											}
										}
										// store debug data
										// if strong enough and enabled - try to improve far objects
										// temporarily removed strength requirement for debugging, restore later to make faster
///										if ((imgdtt_params.fo_correct && (strength > imgdtt_params.fo_min_strength)) || (clt_mismatch != null)) {

										if ((imgdtt_params.fo_correct && (disp_str[1] > 0 * imgdtt_params.fo_min_strength)) || (clt_mismatch != null)) {
								    		// try all dirs:
											dir_corr_strength = corr2d.corr4dirsLMA(
								    				imgdtt_params,                // ImageDttParameters  imgdtt_params,
								    				corrs,                        // double [][]         corrs,
								    				imgdtt_params.dbg_pair_mask,  // int                 pair_mask, // which pairs to process
								    				-disp_str[0],                 // double    xcenter,   // preliminary center x in pixels for largest baseline
								    				imgdtt_params.ortho_vasw_pwr, // double    vasw_pwr,  // value as weight to this power,
								    				tile_lma_debug_level,         // int                 debug_level,
								    				tileX,                        // int                 tileX, // just for debug output
								    				tileY );                      // int                 tileY
								    		if ((tile_lma_debug_level > 0) && (dir_corr_strength != null)) {
								    			double [] nan2 = {Double.NaN, Double.NaN, Double.NaN, Double.NaN};
								    			for (int ii = 0; ii < dir_corr_strength.length; ii++) {
								    				if (dir_corr_strength[ii] == null) dir_corr_strength[ii] = nan2;
								    			}
								    			System.out.println(String.format("corr4dirsLMA -> ↔: %7.4f (%7.4f) ↕:%7.4f (%7.4f) ⇖:%7.4f (%7.4f) ⇗:%7.4f (%7.4f)",
								    					dir_corr_strength[0][0],dir_corr_strength[0][1],dir_corr_strength[1][0],dir_corr_strength[1][1],
								    					dir_corr_strength[2][0],dir_corr_strength[2][1],dir_corr_strength[3][0],dir_corr_strength[3][1]));
								    		}

								    		mod_disparity_diff =     corr2d.foregroundCorrect(
								    				imgdtt_params.fo_far,            // boolean   bg,
								    				imgdtt_params.fo_ortho,          // boolean   ortho,
								    				dir_corr_strength,               // double [] dir_disp,
								    				disp_str[0],                     // double    full_disp,
								    	    		imgdtt_params.fo_min_strength,   // double      min_strength,
								    	    		imgdtt_params.fo_min_eff,        // double      min_eff,
								    	    		imgdtt_params.fo_min_eff_ratio,  // double      min_eff_ratio,
								    	    		imgdtt_params.fo_max_hwidth,    // double      max_hwidth, //  =          3.0;  // maximal half-width disparity  direction to try to correct far objects
								    	    		imgdtt_params.fo_min_diff,       // double    fo_min_diff,
								    	    		imgdtt_params.fo_overcorrection, // double    fo_overcorrection,
								    	    		imgdtt_params.fo_lim_overcorr,   // double    fo_lim_overcorr,
								    	    		(tile_lma_debug_level > 0) );    // boolean debug);

// Do not use modified far object distance when mismatch is measured
											if ((mod_disparity_diff[0] != disp_str[0]) && (clt_mismatch == null)){ // if it changed
												if (imgdtt_params.fo_correct && (disp_str[1] > imgdtt_params.fo_min_strength)) { // always
													disp_str[0] = mod_disparity_diff[0];
													disparity_map[DISPARITY_INDEX_CM]       [tIndex] = disp_str[0];
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
												if (imgdtt_params.ly_poly) {  // not used in lwir
													mismatch_result = corr2d.mismatchPairs( // returns x-xcenter, y, strength (sign same as disparity)
															imgdtt_params,                // ImageDttParameters  imgdtt_params,
															corrs,                        // double [][]         corrs,
															all_pairs,                    // int                 pair_mask, // which pairs to process
															-disp_str[0],                   // double    xcenter,   // preliminary center x in pixels for largest baseline
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
															-disp_str[0],                 // double    xcenter,   // preliminary center x in pixels for largest baseline
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
							else if (corr_mode == 3) extra_disparity = disparity_map[DISPARITY_INDEX_HOR][tIndex];  // not used in lwir
							else if (corr_mode == 4) extra_disparity = disparity_map[DISPARITY_INDEX_VERT][tIndex];  // not used in lwir
							if (Double.isNaN(extra_disparity)) extra_disparity = 0;  // used in lwir

							if (Double.isNaN(disparity_map[DISPARITY_STRENGTH_INDEX][tIndex])) {
								System.out.println("BUG: 3. disparity_map[DISPARITY_STRENGTH_INDEX][tIndex] should not be NaN");
							}
						} // if (disparity_map != null){ // not null - calculate correlations
						// only debug is left
						// old (per-color correlation)
						// ****** FIXME tries to use color == 3, should be disabled!
						if ((clt_corr_combo != null)  && !imgdtt_params.corr_mode_debug){ // not null - calculate correlations  // not used in lwir
							tcorr_tpartial=  new double[CORR_PAIRS.length][numcol+1][4][transform_len];
							tcorr_partial =  new double[numSensors][numcol+1][];

							for (int pair = 0; pair < CORR_PAIRS.length; pair++){
								for (int ncol = 0; ncol <numcol; ncol++) if (clt_data[ncol] != null){
									double [][] data1 = clt_data[CORR_PAIRS[pair][0]][ncol][tileY][tileX];
									double [][] data2 = clt_data[CORR_PAIRS[pair][1]][ncol][tileY][tileX];
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
													if (ZI[n][k] < 0)
														tcorr_tpartial[pair][ncol][n][i] -=
														data1[-ZI[n][k]][i] * data2[k][i];
													else
														tcorr_tpartial[pair][ncol][n][i] +=
														data1[ZI[n][k]][i] * data2[k][i];
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
										for (int ncol= 0; ncol < tcorr_tpartial[pair].length; ncol++) { // tcorr_tpartial[pair].length = 4 > colors
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
								if (CORR_PAIRS[pair][2] != 0) {
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
							for (int pair = 0; pair < CORR_PAIRS.length; pair++) if (((corr_mask >> pair) & 1) != 0){
								numPairs++;
								if (CORR_PAIRS[pair][2] == 0) { // horizontal pair)
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
										for (int pair = 0; pair < CORR_PAIRS.length; pair++) if (((corr_mask >> pair) & 1) != 0){
											tcorr_combo[TCORR_COMBO_RSLT][i] += avScale*tcorr_partial[pair][numcol][i]; // only composite color channel
											if (CORR_PAIRS[pair][2] == 0) { // horizontal pair
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
										for (int pair = 0; pair < CORR_PAIRS.length; pair++) if (((corr_mask >> pair) & 1) != 0){
											tcorr_combo[TCORR_COMBO_RSLT][i] *= (tcorr_partial[pair][numcol][i] + corr_offset); // only composite color channel
											if (CORR_PAIRS[pair][2] == 0) { // horizontal pair
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
									for (int pair = 0; pair < CORR_PAIRS.length; pair++) if (((corr_mask >> pair) & 1) != 0){
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
								// shift images by 0.5 * extra disparity in the diagonal direction  // not used in lwir
								for (int ncol = 0; ncol <numcol; ncol++) { // color
									for (int i = 0; i < numSensors; i++) {
										if (clt_data[i][ncol] != null) {
											fract_shift(    // fractional shift in transform domain. Currently uses sin/cos - change to tables with 2? rotations
													clt_data[i][ncol][tileY][tileX], // double  [][]  clt_tile,
//													transform_size,
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
							double [][][] iclt_tile = new double [numSensors][numcol][]; // in mono some may remain null
							double [] clt_tile;
							double scale = 0.25;  // matching iclt_2d
							for (int i = 0; i < numSensors; i++) {  // USED in lwir
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
								double [][] dbg_tile = new double [numSensors*numcol][];
								for (int i = 0; i < numSensors; i++) {
									for (int ncol = 0; ncol <numcol; ncol++) if (iclt_tile[i][ncol] != null) { // color
										dbg_tile[i * numcol + ncol] = iclt_tile[i][ncol];
									}
								}
								sdfa_instance.showArrays(dbg_tile, 2* transform_size, 2* transform_size, true, "iclt_x"+tileX+"_y"+tileY, titles);
							}


							// "de-bayer" tiles for matching, use original data for output
							double [][][] tiles_debayered = new double [numSensors][numcol][];
							for (int i =0; i<numSensors; i++){
								for (int ncol = 0; ncol < numcol; ncol++) if (iclt_tile[i][ncol] != null) {
									if (isMonochrome()) {  // used in lwir
										tiles_debayered[i][ncol] =  iclt_tile[i][ncol];
									} else {  // used in lwir
										tiles_debayered[i][ncol] =  tile_debayer_shot_corr(
												(ncol != 2), // red or blue (false - green)
												iclt_tile[i][ncol],
												2 * transform_size,
												lt_window2, // squared lapping window
												min_shot,   // 10.0;  // Do not adjust for shot noise if lower than
												scale_shot,  //3.0;   // scale when dividing by sqrt
												lt_window2, // re-apply window to the result
												false); //boolean  debug_gpu);

									}
								}
							}
							if ((globalDebugLevel > 0) && debugTile) {
								ShowDoubleFloatArrays sdfa_instance = new ShowDoubleFloatArrays(); // just for debugging?
								String [] titles = {"red0","blue0","green0","red1","blue1","green1","red2","blue2","green2","red3","blue3","green3"};
								double [][] dbg_tile = new double [numSensors*numcol][];
								for (int i = 0; i < numSensors; i++) {
									for (int chn = 0; chn <numcol; chn++) { // color
										dbg_tile[i * numcol + chn] = tiles_debayered[i][chn];
									}
								}
								sdfa_instance.showArrays(dbg_tile, 2* transform_size, 2* transform_size, true, "tiles_debayered_x"+tileX+"_y"+tileY, titles);
							}
							// ... used in lwir
							double []     max_diff = null;
							if ((disparity_map != null) && (disparity_map.length >= (IMG_DIFF0_INDEX + numSensors))){
								max_diff = new double[numSensors];
							}
							int ports_rgb_len = numSensors*numcol;  // 12
							if ((disparity_map != null) && (disparity_map.length >= (getImgToneRGB() + ports_rgb_len))) {
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
									(globalDebugLevel > 0) && debugTile,
									false);          // boolean       debug_gpu)      // generate output fro matching with GPU processing


							// mix RGB from iclt_tile, mix alpha with - what? correlation strength or 'don't care'? good correlation or all > min?
							for (int i = 0; i < iclt_tile[0][first_color].length; i++ ) {
								double sw = 0.0;
								for (int ip = 0; ip < numSensors; ip++) {
									sw += texture_tiles[tileY][tileX][numcol+1+ip][i];
								}
								if (sw != 0 ) sw = 1.0/sw;
								for (int ncol = 0; ncol < numcol; ncol++) if (iclt_tile[0][ncol] !=null){ // color
									texture_tiles[tileY][tileX][ncol][i] = 0.0; //iclt[tileY][tileX][chn]
									for (int ip = 0; ip < numSensors; ip++) {
										texture_tiles[tileY][tileX][ncol][i] += sw * texture_tiles[tileY][tileX][numcol+1+ip][i] * iclt_tile[ip][ncol][i];
									}
								}
							}
							if ((disparity_map != null) && (disparity_map.length >= (IMG_DIFF0_INDEX + numSensors))){
								for (int i = 0; i < max_diff.length; i++){
									disparity_map[IMG_DIFF0_INDEX + i][tIndex] = max_diff[i];
								}
							}
							if (ports_rgb != null) {
								for (int i = 0; i < ports_rgb.length; i++){
									disparity_map[getImgToneRGB() + i][tIndex] = ports_rgb[i];
								}
							}
						}
					}
				}
			};
		}
		startAndJoin(threads);

		if ((dbg_distort != null) &&(globalDebugLevel >=0)) {
			(new ShowDoubleFloatArrays()).showArrays(dbg_distort,  tilesX, tilesY, true, "disparity_distortions"); // , dbg_titles);
		}

		return clt_data;
	}

	// reimplementing from GPU version (will also need upgrade for multi-sensor > 4)
	
	public double [][][][][][] quadCorrTD(
			final double [][][]       image_data,      // first index - number of image in a quad
			final int                 width,
			final TpTask []           tp_tasks,
			final ImageDttParameters  imgdtt_params,    // Now just extra correlation parameters, later will include, most others
			final double [][][][][]   dcorr_td,        // [pair][tilesY][tilesX][4][64] sparse transform domain representation of corr pairs
			// no combo here - rotate, combine in pixel domain after interframe
			final double [][][][][][] clt_kernels,     // [channel_in_quad][color][tileY][tileX][band][pixel] , size should match image (have 1 tile around)
			final int                 kernel_step,
			final int                 window_type,
			final double              corr_red,
			final double              corr_blue,

			final int                 mcorr_sel,    // Which pairs to correlate // +1 - all, +2 - dia, +4 - sq, +8 - neibs, +16 - hor + 32 - vert
			final int                 debug_tileX,
			final int                 debug_tileY,
			final int                 threadsMax,       // maximal number of threads to launch
			final int                 globalDebugLevel)
	{
		final int height=image_data[0][0].length/width;
		final int tilesX=width/transform_size;
		final int tilesY=height/transform_size;
//		if (correlation2d == null){
//			  throw new IllegalArgumentException ("quadCorrTD(): correlation2d == null!");
//		}
		// Initialize correlation pairs selection to be used by all threads
		if (correlation2d != null){
			boolean [] corr_calculate = null;
			if (isCorrAll  (mcorr_sel)) corr_calculate = correlation2d.selectAll();
			if (isCorrDia  (mcorr_sel)) corr_calculate = correlation2d.selectDiameters  (corr_calculate);
			if (isCorrSq   (mcorr_sel)) corr_calculate = correlation2d.selectSquares    (corr_calculate);
			if (isCorrNeib (mcorr_sel)) corr_calculate = correlation2d.selectNeibs      (corr_calculate);
			if (isCorrHor  (mcorr_sel)) corr_calculate = correlation2d.selectHorizontal (corr_calculate);
			if (isCorrVert (mcorr_sel)) corr_calculate = correlation2d.selectVertical   (corr_calculate);
			correlation2d.setCorrPairs(corr_calculate); // will limit correlation pairs calculation

			boolean [] calc_corr_pairs = correlation2d.getCorrPairs();
			for (int i = 0; i < calc_corr_pairs.length; i++) if (calc_corr_pairs[i]){
				dcorr_td[i] = new double[tilesY][tilesX][][]; 
			}
		}
//		final int numcol = isMonochrome()?1:3;
		final int numcol = 3; // number of colors // keep the same, just do not use [0] and [1], [2] - green

		final double [] col_weights= new double [numcol]; // colors are RBG
		if (isMonochrome()) {
			col_weights[2] = 1.0;// green color/mono
			col_weights[0] = 0;
			col_weights[1] = 0;
		} else {
			col_weights[2] = 1.0/(1.0 + corr_red + corr_blue);    // green color
			col_weights[0] = corr_red *  col_weights[2];
			col_weights[1] = corr_blue * col_weights[2];
		}
		
		final double [] filter_rb =  isMonochrome() ? null: doubleGetCltLpfFd(imgdtt_params.pcorr_sigma_rb);
		
		DttRad2 dtt = new DttRad2(transform_size);
		dtt.set_window(window_type);
		final double [] lt_window = dtt.getWin2d();	// [256]

		final double [] lt_window2 = new double [lt_window.length]; // squared
		for (int i = 0; i < lt_window.length; i++) lt_window2[i] = lt_window[i] * lt_window[i];

		if (globalDebugLevel > 1) {
			ShowDoubleFloatArrays sdfa_instance = new ShowDoubleFloatArrays(); // just for debugging?
			sdfa_instance.showArrays(lt_window,  2*transform_size, 2*transform_size, "lt_window");
		}
		
		final double [][][][][][] clt_data = new double[numSensors][numcol][tilesY][tilesX][][];
		final Thread[] threads = newThreadArray(threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				@Override
				public void run() {
					DttRad2 dtt = new DttRad2(transform_size);
					dtt.set_window(window_type);
					int tileY,tileX; // ,tIndex; // , chn;
					double [][] fract_shiftsXY = new double[numSensors][];
					for (int iTile = ai.getAndIncrement(); iTile < tp_tasks.length; iTile = ai.getAndIncrement()) if (tp_tasks[iTile].getTask() != 0) {
						tileY = tp_tasks[iTile].getTileY(); //  /tilesX;
						tileX = tp_tasks[iTile].getTileX(); //nTile % tilesX;
						/*
						int                 img_mask =  0xf; // getImgMask(tile_op[tileY][tileX]);         // which images to use
						if (numSensors > 4) {
							if (img_mask == 0xf) {
								for (int i = 0; i < numSensors; i++){
									img_mask |= 1 << i;
								}
							}
						}
						*/
//						boolean debugTile =(tileX == debug_tileX) && (tileY == debug_tileY) && (globalDebugLevel > -1);
//						boolean debugTile0 =(tileX == debug_tileX) && (tileY == debug_tileY) && (globalDebugLevel > -3);
						// TODO: move port coordinates out of color channel loop
						double [][] centersXY = tp_tasks[iTile].getDoubleXY();// isAux());

						// save disparity distortions for visualization:
						// TODO: use correction after disparity applied (to work for large disparity values)
						// See if macro_mode uses color channels for non-color?
						for (int ncol = 0; ncol <numcol; ncol++) {
							if (!isMonochrome() || (ncol == MONO_CHN)) { // in monochrome mode skip all non-mono (green) channels  // used in lwir (5 of 6 branches)

								for (int i = 0; i < numSensors; i++) {
									clt_data[i][ncol][tileY][tileX] = new double [4][];
									// Extract image tiles and kernels, correct aberrations, return (ut do not apply) fractional shifts
									fract_shiftsXY[i] = extract_correct_tile( // return a pair of residual offsets
											image_data[i],
											width,       // image width
											((clt_kernels == null) ? null : clt_kernels[i]), // [color][tileY][tileX][band][pixel]
											clt_data[i][ncol][tileY][tileX], //double  [][]        clt_tile,    // should be double [4][];
											kernel_step,
											dtt,
											ncol,
											centersXY[i][0], // centerX, // center of aberration-corrected (common model) tile, X
											centersXY[i][1], // centerY, //
											0, // ((!FPGA_COMPARE_DATA && (globalDebugLevel > -1) && (tileX == debug_tileX) && (tileY == debug_tileY) && (ncol == 2) && (i==0)) ? (globalDebugLevel + 0) : 0), // external tile compare
											false, // no_deconvolution,
											false, // ); // transpose);
											null, //final boolean [][]        saturation_imp, // (near) saturated pixels or null
											null); // final double [] overexposed)
								} // for (int i = 0; i < quad; i++)
								if ((globalDebugLevel > -1) && (tileX == debug_tileX) && (tileY == debug_tileY) && (ncol == 2)) {
									System.out.println();
								}

								// apply residual shift
								for (int i = 0; i < numSensors; i++) {
									fract_shift(    // fractional shift in transform domain. Currently uses sin/cos - change to tables with 2? rotations
											clt_data[i][ncol][tileY][tileX], // double  [][]  clt_tile,
											fract_shiftsXY[i][0],            // double        shiftX,
											fract_shiftsXY[i][1],            // double        shiftY,
											false);
								}
							} else { // if (!isMonochrome() || (chn == MONO_CHN) || macro_mode) { // in monochrome mode skip all non-mono (green) channels
								for (int i = 0; i < numSensors; i++) {  // used in lwir
									clt_data[i][ncol] = null; // erase unused clt_data
								}
							}
						}// end of for (int chn = 0; chn <numcol; chn++)
						// all color channels are done here
						if (correlation2d != null) { // will only calculate clt_data
							// calculate all selected pairs correlations
							// change filter for lpf_rb (null for mono)
							double [][][] corr_tiles_td =  correlation2d.correlateCompositeTD(
									clt_data,                     // double [][][][][][] clt_data,
									tileX,                        // int                 tileX,
									tileY,                        // int                 tileY,
									correlation2d.getCorrPairs(), // boolean []          pairs_mask,
									filter_rb,                    // double []           lpf_rb, // extra lpf for red and blue (unused for mono) or null
									getScaleStrengths(),          // double              scale_value, // scale correlation value
									col_weights);                 // double []           col_weights)
							for (int pair = 0; pair < corr_tiles_td.length; pair++) if (corr_tiles_td[pair] != null) {
								dcorr_td[pair][tileY][tileX] = corr_tiles_td[pair]; 
							}
						}
					}
				}
			};
		}
		startAndJoin(threads);
		return clt_data;
	}
	
	
	/*
	public double [][] quadCorrTD_tilted( // process tilted multiframe, returns per-tile weights (for fat zero application)
			final double [][][]       image_data,      // first index - number of image in a quad
			final int                 width,
			final TpTask []           tp_tasks,        // should exclude strong tiles (e.g. having disparity_lma)
			final ImageDttParameters  imgdtt_params,   // Now just extra correlation parameters, later will include, most others
			final double [][][][][]   dcorr_td,        // [pair][tilesY][tilesX][4][64] sparse transform domain representation of corr pairs
			final double [][][][][][] clt_kernels,     // [channel_in_quad][color][tileY][tileX][band][pixel] , size should match image (have 1 tile around)
			final int                 kernel_step,
			final int                 window_type,
			final double              corr_red,
			final double              corr_blue,
			
			final int                 clustRadius,  // 1 - single tile, 2 - 3x3, 3 - 5x5, ...
			final double              arange,  // absolute disparity range to consolidate
			final double              rrange,  // relative disparity range to consolidate
			final double              no_tilt, // no tilt if center disparity is lower
			final double              damp_tilt,    // 0.1?

			final int                 mcorr_sel,    // Which pairs to correlate // +1 - all, +2 - dia, +4 - sq, +8 - neibs, +16 - hor + 32 - vert
			final int                 debug_tileX,
			final int                 debug_tileY,
			final int                 threadsMax,       // maximal number of threads to launch
			final int                 globalDebugLevel)
	{
		final int height=image_data[0][0].length/width;
		final int tilesX=width/transform_size;
		final int tilesY=height/transform_size;
		final double [][] tile_weights = new double [tilesY][tilesX];
		final TpTask [][] tp_tasks_full = new TpTask[tilesY][tilesX]; 
		final double [] damping = {damp_tilt, damp_tilt, 0.0}; // 0.0 will be applied to average value, tilt_cost - to both tilts
		final int clustDiameter = 2 * clustRadius - 1;
		final int center_indx = (clustRadius - 1) * (clustDiameter + 1);
		final double [] wnd_neib = new double [clustDiameter * clustDiameter];
		for (int iy = 0; iy < clustDiameter; iy++) {
			double wy = Math.sin(Math.PI *(iy+1) / (clustDiameter + 1));
			for (int ix = 0; ix < clustDiameter; ix++) {
				wnd_neib[iy * clustDiameter + ix] = wy * Math.sin(Math.PI *(ix+1) / (clustDiameter + 1));
			}
		}
		
//		if (correlation2d == null){
//			  throw new IllegalArgumentException ("quadCorrTD(): correlation2d == null!");
//		}
		// Initialize correlation pairs selection to be used by all threads
		if (correlation2d != null){
			boolean [] corr_calculate = null;
			if (isCorrAll  (mcorr_sel)) corr_calculate = correlation2d.selectAll();
			if (isCorrDia  (mcorr_sel)) corr_calculate = correlation2d.selectDiameters  (corr_calculate);
			if (isCorrSq   (mcorr_sel)) corr_calculate = correlation2d.selectSquares    (corr_calculate);
			if (isCorrNeib (mcorr_sel)) corr_calculate = correlation2d.selectNeibs      (corr_calculate);
			if (isCorrHor  (mcorr_sel)) corr_calculate = correlation2d.selectHorizontal (corr_calculate);
			if (isCorrVert (mcorr_sel)) corr_calculate = correlation2d.selectVertical   (corr_calculate);
			correlation2d.setCorrPairs(corr_calculate); // will limit correlation pairs calculation

			boolean [] calc_corr_pairs = correlation2d.getCorrPairs();
			for (int i = 0; i < calc_corr_pairs.length; i++) if (calc_corr_pairs[i]){
				dcorr_td[i] = new double[tilesY][tilesX][][]; 
			}
		}
//		final int numcol = isMonochrome()?1:3;
		final int numcol = 3; // number of colors // keep the same, just do not use [0] and [1], [2] - green

		final double [] col_weights= new double [numcol]; // colors are RBG
		if (isMonochrome()) {
			col_weights[2] = 1.0;// green color/mono
			col_weights[0] = 0;
			col_weights[1] = 0;
		} else {
			col_weights[2] = 1.0/(1.0 + corr_red + corr_blue);    // green color
			col_weights[0] = corr_red *  col_weights[2];
			col_weights[1] = corr_blue * col_weights[2];
		}
		
		final double [] filter_rb =  isMonochrome() ? null: doubleGetCltLpfFd(imgdtt_params.pcorr_sigma_rb);
		
		DttRad2 dtt = new DttRad2(transform_size);
		dtt.set_window(window_type);
		final double [] lt_window = dtt.getWin2d();	// [256]

		final double [] lt_window2 = new double [lt_window.length]; // squared
		for (int i = 0; i < lt_window.length; i++) lt_window2[i] = lt_window[i] * lt_window[i];

		if (globalDebugLevel > 1) {
			ShowDoubleFloatArrays sdfa_instance = new ShowDoubleFloatArrays(); // just for debugging?
			sdfa_instance.showArrays(lt_window,  2*transform_size, 2*transform_size, "lt_window");
		}
		
		final double [][][][][][] clt_data = new double[numSensors][numcol][tilesY][tilesX][][];
		final Thread[] threads = newThreadArray(threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				@Override
				public void run() {
					DttRad2 dtt = new DttRad2(transform_size);
					dtt.set_window(window_type);
					int tileY,tileX; // ,tIndex; // , chn;
					for (int iTile = ai.getAndIncrement(); iTile < tp_tasks.length; iTile = ai.getAndIncrement()) if (tp_tasks[iTile].getTask() != 0){
						tileY = tp_tasks[iTile].getTileY(); //  /tilesX;
						tileX = tp_tasks[iTile].getTileX(); //nTile % tilesX;
						tp_tasks_full[tileY][tileX] = tp_tasks[iTile];
					}
				}
			};
		}
		startAndJoin(threads);
		
		ai.set(0);
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				@Override
				public void run() {
					DttRad2 dtt = new DttRad2(transform_size);
					dtt.set_window(window_type);
					int tileYC,tileXC; // ,tIndex; // , chn;
					double [][] fract_shiftsXY = new double[numSensors][];
					PolynomialApproximation pa = new PolynomialApproximation();
					TileNeibs tn = new TileNeibs(tilesX,tilesY);
					double [][][][] clt_data_tile = new double[numSensors][numcol][][];
					for (int iTile = ai.getAndIncrement(); iTile < tp_tasks.length; iTile = ai.getAndIncrement()) if (tp_tasks[iTile].getTask() != 0){
						tileYC = tp_tasks[iTile].getTileY(); //  /tilesX;
						tileXC = tp_tasks[iTile].getTileX(); //nTile % tilesX;
						boolean debugTile =(tileXC == debug_tileX) && (tileYC == debug_tileY) && (globalDebugLevel > -1);
						boolean debugTile0 =(tileXC == debug_tileX) && (tileYC == debug_tileY) && (globalDebugLevel > -3);
						double [][] centersXY = tp_tasks[iTile].getDoubleXY();// isAux());
						
						if (debugTile0) {
							System.out.println("tileXC = "+tileXC+", tileYC = "+tileYC);
						}
						//calculate tilts
						TileNeibs tn = new TileNeibs(tilesX,tilesY);
						boolean [] used_neibs = new boolean [clustDiameter * clustDiameter];
						double disparity_center = disparity_array[tileYC][tileXC];
						double disp_min = disparity_center - arange - rrange * Math.abs(disparity_center);
						double disp_max = disparity_center + arange + rrange * Math.abs(disparity_center);
						double tiltY = 0, tiltX=0;
						double [][][] mdata = new double [wnd_neib.length][][]; // now empty lines will be skipped [3][];
						double sum_w = 0.0;
						for (int dty = -clustRadius+1; dty < clustRadius; dty++) {
							for (int dtx = -clustRadius+1; dtx < clustRadius; dtx++) {
								int nTile1 = tn.getNeibIndex(nTileC, dtx, dty);
								if (nTile1 >= 0){
									int tileX1 = nTile1 % tilesX;
									int tileY1 = nTile1 / tilesX;
									if ((tile_op[tileY1][tileX1] != 0) && (disparity_array[tileY1][tileX1] >= disp_min)  && (disparity_array[tileY1][tileX1] <= disp_max)){
										int mindx = (dty * clustDiameter + dtx + center_indx);
										used_neibs[mindx] = true;
										double w = wnd_neib[mindx];
										mdata[mindx] = new double[3][];
										mdata[mindx][0] = new double [2];
										mdata[mindx][0][0] =  dtx;
										mdata[mindx][0][1] =  dty;
										mdata[mindx][1] = new double [1];
										mdata[mindx][1][0] =  disparity_array[tileY1][tileX1]; 
										mdata[mindx][2] = new double [1];
										mdata[mindx][2][0] =  w;
										sum_w += w;
									}
								}
							}
						}
						double scale_weighths = 1.0/sum_w;
						if (disparity_center > no_tilt) {
							double[][] approx2d = pa.quadraticApproximation(
									mdata,
									true,          // boolean forceLinear,  // use linear approximation
									damping,       // double [] damping,
									-1);           // debug level
							if (approx2d != null){
								tiltX = approx2d[0][0];
								tiltY = approx2d[0][1]; // approx2d[0][2] - const C (A*x+B*y+C), C is not used here 
							}
//							if (num_tiles == 0) {
//								continue; // should never happen anyway
//							}
						}
						
						
						
//						boolean debugTile =(tileX == debug_tileX) && (tileY == debug_tileY) && (globalDebugLevel > -1);
//						boolean debugTile0 =(tileX == debug_tileX) && (tileY == debug_tileY) && (globalDebugLevel > -3);
						// TODO: move port coordinates out of color channel loop
						double [][] centersXY = tp_tasks[iTile].getDoubleXY();// isAux());

						// save disparity distortions for visualization:
						// TODO: use correction after disparity applied (to work for large disparity values)
						// See if macro_mode uses color channels for non-color?
						for (int ncol = 0; ncol <numcol; ncol++) {
							if (!isMonochrome() || (ncol == MONO_CHN)) { // in monochrome mode skip all non-mono (green) channels  // used in lwir (5 of 6 branches)

								for (int i = 0; i < numSensors; i++) {
									clt_data[i][ncol][tileY][tileX] = new double [4][];
									// Extract image tiles and kernels, correct aberrations, return (ut do not apply) fractional shifts
									fract_shiftsXY[i] = extract_correct_tile( // return a pair of residual offsets
											image_data[i],
											width,       // image width
											((clt_kernels == null) ? null : clt_kernels[i]), // [color][tileY][tileX][band][pixel]
											clt_data[i][ncol][tileY][tileX], //double  [][]        clt_tile,    // should be double [4][];
											kernel_step,
											dtt,
											ncol,
											centersXY[i][0], // centerX, // center of aberration-corrected (common model) tile, X
											centersXY[i][1], // centerY, //
											0, // ((!FPGA_COMPARE_DATA && (globalDebugLevel > -1) && (tileX == debug_tileX) && (tileY == debug_tileY) && (ncol == 2) && (i==0)) ? (globalDebugLevel + 0) : 0), // external tile compare
											false, // no_deconvolution,
											false, // ); // transpose);
											null, //final boolean [][]        saturation_imp, // (near) saturated pixels or null
											null); // final double [] overexposed)
								} // for (int i = 0; i < quad; i++)
								if ((globalDebugLevel > -1) && (tileX == debug_tileX) && (tileY == debug_tileY) && (ncol == 2)) {
									System.out.println();
								}

								// apply residual shift
								for (int i = 0; i < numSensors; i++) {
									fract_shift(    // fractional shift in transform domain. Currently uses sin/cos - change to tables with 2? rotations
											clt_data[i][ncol][tileY][tileX], // double  [][]  clt_tile,
											fract_shiftsXY[i][0],            // double        shiftX,
											fract_shiftsXY[i][1],            // double        shiftY,
											//									(globalDebugLevel > 0) && (tileX == debug_tileX) && (tileY == debug_tileY)); // external tile compare
											((globalDebugLevel > 1) &&
													((ncol==0) || isMonochrome()) &&
													(tileX >= debug_tileX - 2) && (tileX <= debug_tileX + 2) &&
													(tileY >= debug_tileY - 2) && (tileY <= debug_tileY+2)));
								}
							} else { // if (!isMonochrome() || (chn == MONO_CHN) || macro_mode) { // in monochrome mode skip all non-mono (green) channels
								for (int i = 0; i < numSensors; i++) {  // used in lwir
									clt_data[i][ncol] = null; // erase unused clt_data
								}
							}
						}// end of for (int chn = 0; chn <numcol; chn++)
						// all color channels are done here
						if (correlation2d != null) { // will only calculate clt_data
							// calculate all selected pairs correlations
							// change filter for lpf_rb (null for mono)
							double [][][] corr_tiles_td =  correlation2d.correlateCompositeTD(
									clt_data,                     // double [][][][][][] clt_data,
									tileX,                        // int                 tileX,
									tileY,                        // int                 tileY,
									correlation2d.getCorrPairs(), // boolean []          pairs_mask,
									filter_rb,                    // double []           lpf_rb, // extra lpf for red and blue (unused for mono) or null
									getScaleStrengths(),          // double              scale_value, // scale correlation value
									col_weights);                 // double []           col_weights)
							for (int pair = 0; pair < corr_tiles_td.length; pair++) if (corr_tiles_td[pair] != null) {
								dcorr_td[pair][tileY][tileX] = corr_tiles_td[pair]; 
							}
						}
						
						
						
					}
				}
			};
		}
		startAndJoin(threads);
		
		
		
		return tile_weights;
	}
	
*/	
	
	
	
	
	public void clt_process_tl_correlations( // convert to pixel domain and process correlations already prepared in fcorr_td and/or fcorr_combo_td
			final ImageDttParameters  imgdtt_params,   // Now just extra correlation parameters, later will include, most others
			final TpTask []           tp_tasks,        // data from the reference frame - will be applied to LMW for the integrated correlations
			                                           // only listed tiles will be processed
			final double [][]         rXY,             // from geometryCorrection
			// no fcorr_combo_td here both arrays should have same non-null tiles
			final double [][][][][]   dcorr_td,      // [pair][tilesY][tilesX][4][64] sparse transform domain representation of corr pairs
			// next both can be nulls
		    final double [][][][]     clt_corr_out,   // sparse (by the first index) [type][tilesY][tilesX][(2*transform_size-1)*(2*transform_size-1)] or null
  		    final double [][][][]     clt_combo_out,  // sparse (by the first index) [>=1][tilesY][tilesX][(combo_tile_size] or null
			// to be converted to float
			final double  [][][]      dcorr_tiles,     // [tile][pair][(2*transform_size-1)*(2*transform_size-1)] // if null - will not calculate
			// When clt_mismatch is non-zero, no far objects extraction will be attempted
			//optional, may be null
			final double [][]         disparity_map,   // [8][tilesY][tilesX], only [6][] is needed on input or null - do not calculate
			                                           // last 2 - contrast, avg/ "geometric average)
  		    final double              afat_zero2,      // gpu_fat_zero ==30? clt_parameters.getGpuFatZero(is_mono); absolute fat zero, same units as components squared values
  		    final double              corr_sigma,      //
  		    // define combining of all 2D correlation pairs for CM (LMA does not use them)
			final int                 mcorr_comb_width,  // combined correlation tile width
			final int                 mcorr_comb_height, // combined correlation tile full height
			final int                 mcorr_comb_offset, // combined correlation tile height offset: 0 - centered (-height/2 to height/2), height/2 - only positive (0 to height)
			final double              mcorr_comb_disp,   // Combined tile per-pixel disparity for baseline == side of a square
  		    
			final int                 window_type,     // GPU: will not be used
			final int                 debug_tileX,
			final int                 debug_tileY,
			final int                 threadsMax,      // maximal number of threads to launch
			final int                 globalDebugLevel)
	{
		final boolean combine_corrs = (mcorr_comb_width > 0); // use 0 to disable combining (use LMA only)
		int tx=-1, ty=-1;
		for (int np = 0; np < dcorr_td.length; np++) if (dcorr_td[np] != null){
			for (int ity = 0; ity < dcorr_td[np].length; ity++) if (dcorr_td[np][ity] != null){
				ty = dcorr_td[np].length;;
				tx = dcorr_td[np][ity].length;
				ity =  dcorr_td[np].length;
				np = dcorr_td.length;
				break;
			}
		}
		if (tx < 0) {
			System.out.println("clt_process_tl_correlations(): dcorr_td is empty");
			return;
		}
		final int tilesX=tx;
		final int tilesY=ty;
		final double [] corr_lpf =     doubleGetCltLpfFd(corr_sigma);

		if (clt_corr_out != null) {
			boolean [] calc_corr_pairs = correlation2d.getCorrPairs();
			for (int i = 0; i < calc_corr_pairs.length; i++) if (calc_corr_pairs[i]){
				clt_corr_out[i] = new double[tilesY][tilesX][]; 
			}
		}
		
		if (combine_corrs) {
			correlation2d.generateResample( // should be called before
					mcorr_comb_width,  // combined correlation tile width
					mcorr_comb_height, // combined correlation tile full height
					mcorr_comb_offset, // combined correlation tile height offset: 0 - centered (-height/2 to height/2), height/2 - only positive (0 to height)
					mcorr_comb_disp);
			if (clt_combo_out != null) {
				clt_combo_out[0] = new double[tilesY][tilesX][]; // only supports [0] - all pairs available 
			}
		}
		final double [][] corr_wnd = Corr2dLMA.getCorrWnd(
				transform_size,
				imgdtt_params.lma_wnd);
		final double [] corr_wnd_inv_limited = (imgdtt_params.lma_min_wnd <= 1.0)?  new double [corr_wnd.length * corr_wnd[0].length]: null;
		if (corr_wnd_inv_limited != null) {
			double inv_pwr = imgdtt_params.lma_wnd_pwr - (imgdtt_params.lma_wnd - 1.0); // compensate for lma_wnd
			for (int i = imgdtt_params.lma_hard_marg; i < (corr_wnd.length - imgdtt_params.lma_hard_marg); i++) {
				for (int j = imgdtt_params.lma_hard_marg; j < (corr_wnd.length - imgdtt_params.lma_hard_marg); j++) {
					corr_wnd_inv_limited[i * (corr_wnd.length) + j] = 1.0/Math.max(Math.pow(corr_wnd[i][j],inv_pwr), imgdtt_params.lma_min_wnd);
				}
			}
		}
		
		if (disparity_map != null){
			int [] used_slices = {DISPARITY_INDEX_CM, DISPARITY_INDEX_CM+1,DISPARITY_INDEX_POLY,DISPARITY_INDEX_POLY+1,DISPARITY_STRENGTH_INDEX};
			int [] nan_slices =  {DISPARITY_INDEX_CM, DISPARITY_INDEX_POLY};
			for (int indx:used_slices) {
				disparity_map[indx] = new double[tilesY*tilesX];
			}
			for (int indx:nan_slices) {
				Arrays.fill(disparity_map[indx], Double.NaN);
			}
		}
		
		final Thread[] threads = newThreadArray(threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				@Override
				public void run() {
					DttRad2 dtt = new DttRad2(transform_size);
					dtt.set_window(window_type);
					int tileY,tileX,nTile; // , chn;
					for (int iTile = ai.getAndIncrement(); iTile < tp_tasks.length; iTile = ai.getAndIncrement()) {
						tileY = tp_tasks[iTile].getTileY(); //  /tilesX;
						tileX = tp_tasks[iTile].getTileX(); //nTile % tilesX;
						nTile = tileY * tilesX + tileX;
						if (tp_tasks[iTile].getTask() == 0) continue; // nothing to do for this tile
						int                 img_mask =  0xf; // getImgMask(tile_op[tileY][tileX]);         // which images to use
						if (numSensors > 4) {
							if (img_mask == 0xf) {
								for (int i = 0; i < numSensors; i++){
									img_mask |= 1 << i;
								}
							}
						}
						boolean debugTile0 =(tileX == debug_tileX) && (tileY == debug_tileY) && (globalDebugLevel > -3);
						// TODO: move port coordinates out of color channel loop
						double [][] centersXY = tp_tasks[iTile].getDoubleXY(); // isAux());
						double [][] disp_dist = tp_tasks[iTile].getDoubleDispDist();
						double [][] corrs =      new double [dcorr_td.length][];
						for (int i = 0; i < dcorr_td.length; i++) if (dcorr_td[i] != null) {
							corrs[i] = correlation2d. normalizeConvertCorr(
									dcorr_td[i][tileY][tileX], // double [][] tcorr, // null or initialized to [4][transform_len]
									corr_lpf, //   double []   lpf,
									afat_zero2, //   double      afat_zero2, // absolute fat zero, same units as components squared values
						    		0, //   int corr_radius,
						    		false); //   boolean debug_gpu)							
						}
						
						if (dcorr_tiles != null) {
							dcorr_tiles[iTile] =  corrs;
//							for (int i = 0; i < corrs.length; i++) { //  if (corrs[i] != null) {
//								dcorr_tiles[iTile][i] =  corrs[i];
//							}
						}
						if (clt_corr_out != null) {
							for (int num_pair = 0; num_pair < corrs.length; num_pair++) if (corrs[num_pair] != null){
								clt_corr_out[num_pair][tileY][tileX] = corrs[num_pair];
							}
						}
						// get CM disparity/strength
						double [] disp_str = {0.0, 0.0}; // diaprity = 0 will be initial approximation for LMA if no averaging
						if (combine_corrs) {
							double [] corr_combo_tile = correlation2d.accumulateInit(); // combine all available pairs
							double sumw = correlation2d.accummulatePairs(
									corr_combo_tile,           // double []   accum_tile,
									corrs,                     // double [][] corr_tiles,
									correlation2d.selectAll(), // boolean []  selection,
									1.0);             // double      weight);
							correlation2d.normalizeAccumulatedPairs(
									corr_combo_tile,
									sumw);
							// copy to output for monitoring if non-null;
							if (clt_combo_out != null) {
								clt_combo_out[0][tileY][tileX] = corr_combo_tile;
							}

							int [] ixy =   correlation2d.getMaxXYInt( // find integer pair or null if below threshold // USED in lwir
									corr_combo_tile, // double [] data,      // [data_size * data_size]
									null, // disp_str_combo,
									correlation2d.getCombWidth(), //       data_width,
									correlation2d.getCombHeight()/2 - correlation2d.getCombOffset(), // int       center_row, ??????????????
									true, // boolean   axis_only,
									-1.0, // imgdtt_params.min_corr,  // ???? double    minMax,    // minimal value to consider (at integer location, not interpolated)
									false); // debugCluster); // tile_lma_debug_level > 0); // boolean   debug);
							double [] corr_stat = correlation2d.getMaxXCm(         // get fractional center as a "center of mass" inside circle/square from the integer max
									corr_combo_tile,                      // double [] data,      // [data_size * data_size]
									correlation2d.getCombWidth(),        // int       data_width,      //  = 2 * transform_size - 1;
									correlation2d.getCombHeight()/2 - correlation2d.getCombOffset(),// int       center_row,
									ixy[0],                              // int       ixcenter,  // integer center x
									false); // debugCluster); // (tile_lma_debug_level > 0)); // boolean   debug);
							if (corr_stat != null) { // almost always
								disp_str = new double [] {-corr_stat[0], corr_stat[1]};
								if (disparity_map!=null) {
									disparity_map[DISPARITY_INDEX_CM      ][nTile] = disp_str[0];
									disparity_map[DISPARITY_INDEX_CM + 1  ][nTile] = disp_str[1];
									disparity_map[DISPARITY_STRENGTH_INDEX][nTile] = disp_str[1];
								}
							}
						}

						// debug new LMA correlations
						if (debugTile0) { // should be debugTile
							System.out.println("Will run new LMA for tileX="+tileX+", tileY="+tileY);
						}
						double [] poly_disp = {Double.NaN, 0.0};
						Corr2dLMA lma2 = correlation2d.corrLMA2Single(
								imgdtt_params,                // ImageDttParameters  imgdtt_params,
								imgdtt_params.lmas_LY_single, // false,                        // boolean             adjust_ly, // adjust Lazy Eye
								corr_wnd,                     // double [][]         corr_wnd, // correlation window to save on re-calculation of the window
								corr_wnd_inv_limited,         // corr_wnd_limited, // correlation window, limited not to be smaller than threshold - used for finding max/convex areas (or null)
								corrs,                        // corrs,          // double [][]         corrs,
								disp_dist,
								rXY,                          // double [][]         rXY, // non-distorted X,Y offset per nominal pixel of disparity
								// all that are not null in corr_tiles
								correlation2d.selectAll(),    // longToArray(imgdtt_params.dbg_pair_mask),  // int                 pair_mask, // which pairs to process
								disp_str,  //corr_stat[0],                 // double    xcenter,   // preliminary center x in pixels for largest baseline
								poly_disp,                    // double[]            poly_ds,    // null or pair of disparity/strength
								imgdtt_params.ortho_vasw_pwr, // double    vasw_pwr,  // value as weight to this power,
								-2, //0,                            // tile_lma_debug_level, // +2,         // int                 debug_level,
								tileX,                        // int                 tileX, // just for debug output
								tileY );                      // int                 tileY
						if (debugTile0) { // should be debugTile
							System.out.println("Ran LMA for tileX="+tileX+", tileY="+tileY);
						}
						double [][] ds = null;
						if (lma2 != null) {
							ds = lma2.lmaDisparityStrength(
									imgdtt_params.lmas_min_amp,      //  minimal ratio of minimal pair correlation amplitude to maximal pair correlation amplitude
									imgdtt_params.lmas_max_rel_rms,  // maximal relative (to average max/min amplitude LMA RMS) // May be up to 0.3)
									imgdtt_params.lmas_min_strength, // minimal composite strength (sqrt(average amp squared over absolute RMS)
									imgdtt_params.lmas_min_ac,       // minimal of A and C coefficients maximum (measures sharpest point/line)
									imgdtt_params.lmas_min_min_ac,   // minimal of A and C coefficients minimum (measures sharpest point)
									imgdtt_params.lmas_max_area,      //double  lma_max_area,     // maximal half-area (if > 0.0)
									imgdtt_params.lma_str_scale,    // convert lma-generated strength to match previous ones - scale
									imgdtt_params.lma_str_offset    // convert lma-generated strength to match previous ones - add to result
									);
							if (ds != null) { // always true
								if (disparity_map!=null) {
									disparity_map[DISPARITY_INDEX_POLY    ][nTile] = ds[0][0];
									disparity_map[DISPARITY_INDEX_POLY + 1][nTile] = ds[0][1];
									disparity_map[DISPARITY_STRENGTH_INDEX][nTile] = ds[0][1]; // overwrite with LMA strength
								}
								if (debugTile0) {
									lma2.printStats(ds,1);
									if (imgdtt_params.lmas_LY_single) {
										double [][] ddnd = lma2.getDdNd(); // will not be used here
										if (ddnd != null) {
											double [][] dxy= new double [ddnd.length][2];
											for (int i = 0; i < dxy.length; i++) {
												dxy[i][0] = ddnd[i][0] * rXY[i][0] - ddnd[i][1] * rXY[i][1];
												dxy[i][1] = ddnd[i][0] * rXY[i][1] + ddnd[i][1] * rXY[i][0];
											}
											System.out.print("       Port:  ");
											for (int i = 0; i < dxy.length; i++) System.out.print(String.format("   %2d   ", i)); System.out.println();
											System.out.print("Radial_in =  [");
											for (int i = 0; i < dxy.length; i++) System.out.print(String.format(" %6.3f,", ddnd[i][0])); System.out.println("]");
											System.out.print("Tangent_CW = [");
											for (int i = 0; i < dxy.length; i++) System.out.print(String.format(" %6.3f,", ddnd[i][1])); System.out.println("]");
											System.out.print("X =          [");
											for (int i = 0; i < dxy.length; i++) System.out.print(String.format(" %6.3f,", dxy[i][0])); System.out.println("]");
											System.out.print("Y =          [");
											for (int i = 0; i < dxy.length; i++) System.out.print(String.format(" %6.3f,", dxy[i][1])); System.out.println("]");
											System.out.println();
										} else {
											System.out.println("No dd/nd and x/y offsets data is available ");
										}
									} else {
										System.out.println("LY offsets are not measured");
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
	
	// convertions between double and float for compartibility CPU/GPU
	
	// float  [][][][][]   fcorrs_td =       new float[num_scenes][tilesY][tilesX][pair][256];
	// double [][][][][]   dcorr_td,        // [pair][tilesY][tilesX][4][64] sparse transform domain representation of corr pairs
	// change? dcorr_td to  [tilesY][tilesX][pair][4][64] sparse transform domain representation of corr pairs
    public float [][][][] convertCorrTd(
			double [][][][][]   dcorr_td,
			float [][][][] fcorr_td) // may be null
	{
		int tilesX=-1, tilesY=-1; 
		for (int np = 0; np < dcorr_td.length; np++) if (dcorr_td[np] != null){
			for (int ity = 0; ity < dcorr_td[np].length; ity++) if (dcorr_td[np][ity] != null){
				tilesY = dcorr_td[np].length;
				tilesX = dcorr_td[np][ity].length;
				np = dcorr_td.length;
				break;
			}
		}
		if (tilesX < 0) {
			System.out.println("convertCorrTd(): dcorr_td is empty");
			return null;
		}
		if (fcorr_td == null) {
			fcorr_td = new float[tilesY][tilesX][][];
		}
		int num_pairs = dcorr_td.length;
//		float [][][][] fcorr_td = new float [tilesY][tilesX][][];
		//for (int npair = 0;)
		for (int tileY = 0; tileY < tilesY; tileY++) {
			for (int tileX = 0; tileX < tilesX; tileX++) {
//				fcorr_td[tileY][tileX] = new float [dcorr_td[tileY][tileX].length][];
				fcorr_td[tileY][tileX] = new float [num_pairs][];
				for (int npair = 0; npair < num_pairs; npair++) {
					if ((dcorr_td[npair] != null) && (dcorr_td[npair][tileY] != null)) {
						double [][] dtile = dcorr_td[npair][tileY][tileX];
						if (dtile != null) {
							if (fcorr_td[tileY][tileX] == null) {
								fcorr_td[tileY][tileX] = new float[num_pairs][]; // if it was not created before
							}
							fcorr_td[tileY][tileX][npair] = new float [dtile.length * dtile[0].length]; // 256
							int indx = 0;
							for (int quadrant = 0; quadrant < dtile.length; quadrant++) {
								for (int i = 0; i < dtile[quadrant].length; i++) {
									fcorr_td[tileY][tileX][npair][indx++] = (float) dtile[quadrant][i];
								}
							}
						}
					}
				}
			}
		}
		return fcorr_td;
	}
	
    public double [][][][][] convertCorrTd(
			float [][][][] fcorr_td) {
    	int tilesY = fcorr_td.length;
    	int tilesX = -1;
    	int num_pairs = -1;
    	for (int tileY = 0; tileY < fcorr_td.length; tileY++) if (fcorr_td[tileY]!=null){
    		tilesX = fcorr_td[tileY].length;
    		for (int tileX = 0; tileX < fcorr_td[tileY].length; tileX++) if (fcorr_td[tileY][tileX]!=null){
    			num_pairs=fcorr_td[tileY][tileX].length;
    			tileY = fcorr_td.length; // exit from the outer loop
    			break;
    		}
    	}
    	if (num_pairs < 0) {
			System.out.println("convertCorrTd(): fcorr_td is empty");
			return null;
    	}
    	int tile_len = transform_size*transform_size;
//    	double [][][][][] dcorr_td = new double [num_pairs][tilesY][tilesX][][];
    	double [][][][][] dcorr_td = new double [num_pairs][][][][];
    	for (int tileY = 0; tileY < tilesY; tileY++) if (fcorr_td [tileY] != null) {
    		for (int tileX = 0; tileX < tilesX; tileX++)  if (fcorr_td [tileY][tileX] != null){
    			for (int npair = 0; npair < num_pairs; npair++) {
    				if (fcorr_td [tileY][tileX][npair] != null) {
    					if (dcorr_td[npair] == null) {
    						dcorr_td[npair] = new double [tilesY][tilesX][][];
    					}
    					dcorr_td[npair][tileY][tileX] = new double [4][tile_len];
						int indx = 0;
						for (int quadrant = 0; quadrant < 4; quadrant++) {
							for (int i = 0; i < tile_len; i++) {
								dcorr_td[npair][tileY][tileX][quadrant][i] = fcorr_td[tileY][tileX][npair][indx++]; 
							}
						}
    				}
    			}			
    		}
    	}
    	return dcorr_td;
    }
    // conver from double[][][] sequence of correlation tiles from CPU implementation to the same float format to match the GPU one
    public float [][][] convertFcltCorr(
    		double [][][] dcorr_tiles,// [tile][sparse, correlation pair][(2*transform_size-1)*(2*transform_size-1)] // if null - will not calculate
    		float  [][][] fclt_corr) //  new float [tilesX * tilesY][][] or null
    {
    	if (dcorr_tiles == null) {
    		return null;
    	}
    	if (fclt_corr == null) {
    		fclt_corr = new float [dcorr_tiles.length][][];
    	}
    	for (int nTile = 0; nTile < dcorr_tiles.length; nTile++) if (dcorr_tiles[nTile] != null) { // should not be null if the length is not extended
    		fclt_corr[nTile] = new float [dcorr_tiles[nTile].length][];
    		for (int npair = 0; npair < dcorr_tiles[nTile].length; npair++) if (dcorr_tiles[nTile][npair] != null) {
        		fclt_corr[nTile][npair] = new float [dcorr_tiles[nTile][npair].length];
    			for (int i = 0; i < dcorr_tiles[nTile][npair].length; i++) {
    				fclt_corr[nTile][npair][i] = (float) dcorr_tiles[nTile][npair][i];
    			}
    		}
    	}
    	return fclt_corr;
    }
	
	
}
