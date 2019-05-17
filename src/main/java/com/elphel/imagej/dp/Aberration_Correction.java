package com.elphel.imagej.dp;
/**
** -----------------------------------------------------------------------------**
** Aberration_Correction.java
**
** De-mosaic and correct aberrations using array of inverted PSF kernels
** 
**
** Copyright (C) 2010-2011 Elphel, Inc.
**
** -----------------------------------------------------------------------------**
**  
**  Aberration_Correction.java is free software: you can redistribute it and/or modify
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

import ij.*;
import ij.process.*;
import ij.gui.*;
import ij.io.FileInfo;
import ij.io.FileSaver;
import ij.io.OpenDialog;
import ij.io.Opener;
import ij.plugin.frame.*;

import java.awt.*;
import java.awt.event.*;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
//import java.io.FileFilter;
//import java.io.FilenameFilter;
import java.util.Properties;
import java.util.concurrent.atomic.AtomicInteger;

import javax.swing.JFileChooser;
import javax.swing.filechooser.FileFilter;

import com.elphel.imagej.calibration.Aberration_Calibration;
import com.elphel.imagej.calibration.DeBayerScissors;
import com.elphel.imagej.common.DoubleFHT;
import com.elphel.imagej.common.DoubleGaussianBlur;
import com.elphel.imagej.common.ShowDoubleFloatArrays;
import com.elphel.imagej.jp4.JP46_Reader_camera;

//import javax.swing.SwingUtilities;
//import javax.swing.UIManager;

public class Aberration_Correction extends PlugInFrame implements ActionListener {
   /**
	 * 
	 */
	private static final long serialVersionUID = -1507307664341265263L;
private Panel panel1,panel2,panel3,panel4,panel5,panel5a, panel6;
   JP46_Reader_camera JP4_INSTANCE=null;

//   private deBayerScissors debayer_instance;
   private ShowDoubleFloatArrays SDFA_INSTANCE; // just for debugging?
   private DoubleFHT             FHT_INSTANCE;
   static Frame instance;
   static Properties PROPERTIES=new Properties();
   public static int DEBUG_LEVEL = 1;
   public static int MASTER_DEBUG_LEVEL = 1;
   public static boolean UPDATE_STATUS= true; // update ImageJ status info
   public static SplitParameters SPLIT_PARAMETERS = new SplitParameters(
		   2,  // oversample;
		   32, // addLeft
		   32, // addTop
		   32, // addRight
		   32  // addBottom
   );		   
   public static DebayerParameters DEBAYER_PARAMETERS = new DebayerParameters(
		   64,    // size //128;
		   0.5,   // polarStep -  size of largest polar cell to Cartesian one;
		   0.35,  // debayerThreshold - Measuring maximal spectral component amplitude in mid-frequencies. If below this - use default de-bayer mask
		   1.5,   //debayerRelativeWidthGreen - result green mask mpy by scaled default (diamond))
		   1.5,   //debayerRelativeWidthRedblue - result red/blue mask mpy by scaled default (square)
		   1.3,   // debayerRelativeWidthRedblueMain - green mask when applied to red/blue, main (center)
		   2.0,   // debayerRelativeWidthRedblueClones - green mask when applied to red/blue, clones
		   0.3,   // public double debayerGamma; - gamma applied to spectrum amplitude
		   0.5,   // public double debayerBonus; //0.5 //0.8 ; //0.4;  // increase far pixels value  (relative to distance to nearest alias)
		   0.6,   // mainToAlias - relative main/alias amplitudes to enable pixels (i.e. 0.5 means that if alias is >0.5*main, the pixel will be masked out)
		   1.6,   // debayerMaskBlur - sigma for gaussian blur of the green and red/blue masks
		   true,  // debayerUseScissors - use "scissors", if false - just apply "diamond" ands "square" with DEBAYER_PARAMETERS.debayerRelativeWidthGreen and DEBAYER_PARAMETERS.debayerRelativeWidthRedblue
//		   false, // showEnergy - plot debayer high frequency energy (use to select between "scissors" and default uniform

		   false, // debug - display internal data for the tile around selected pixel
		   2592,  // xDebug - near the center  
		   1936,  // yDebug - near the center
		   true  // debayerStacks - show debayer debug images as stack (false - individual)
   );

   public static NonlinParameters NONLIN_PARAMETERS=new NonlinParameters(
		   false, //true, // useRejectBlocksFilter - apply blocking artifacts rejection filter to the mask 
		   true, // combineBothModes - combine (sqrt(A*B) mask with rejected blocks and the one without 
		   256,  // maskFFTSize- FFT size of the sliding filter for the mask
		    32,  // blockPeriod - size of the JPEG macroblock (16) scaled for oversampling
		   1.0,  // rejectFreqSigma - width of the rejection filter zeros (in frequency counts)
		   2.5,  //5.0  // lowPassSigma -  sigma for the nonlinear filtering (higher the sigma, farther from the edges is the "sharpened" image used)
		   1.0,  //0.005, // filtMin -   minimal low-pass filtered squared difference between the corrected and original pixels to trigger sharpness enhancement
		   5.0, //0.025, // filtMax -   squared low-pass filtered difference between the corrected and original pixels, so above that level 100% corrected image is used
 		   1.0, // thresholdCorrection_11 - multiply filtMin and filtMax for channel 1-1 
 		   1.0, // thresholdCorrection_12,
 		   1.0, // thresholdCorrection_13,
 		   1.0, // thresholdCorrection_21,
 		   1.0, // thresholdCorrection_22,
 		   1.0, // thresholdCorrection_23,
 		   1.0, // thresholdCorrection_31,
 		   1.0, // thresholdCorrection_32,
 		   1.0, // thresholdCorrection_33,
		   0.25, // threshold - when blurred intensity is below this value, use it as a denominator
		   true, // useDiffNoiseGains;
			0.0, // noiseGainWeights_0, // r weights used to calculate noise gains from individual color. Currently mask is calculated for green only
			0.0, // noiseGainWeights_1, // b
			1.0, // noiseGainWeights_2, // g
			1.0,  // blurSigma,     // blur sigma for mask calculation (blur convolution kernels for noise gain calculation 
			3.0,  // noiseGainPower;
			// ring filter
		  	true,  // useRingFilter;    // filter out spots on denoise mask
		    0.3,   // minMaxValue;       //  minimal value (relative to filtMax) of the local maximum to be processed
		    1.2,   // overRingThreshold; // ratio of local max. and maximal value in the surrounding ring to trigger filter 
		    1.1,   //  overRingLimit;     // limit values in the center circle to scaled maximum in a ring
		    6.0,   //  ringIR;            // ring inner radius (center circle radius)
		    9.0    //  ringOR;            // ring outer radius

   );
   public static ColorProcParameters COLOR_PROC_PARAMETERS=new ColorProcParameters (
		1.245,   // balanceRed   - manual color balance, gain 1.0 matches 0.0.255.0 range of the input Bayer data
		1.34,   // balanceBlue;
		1.45,   // gain;
		1.0,   // weightScaleR; - (now not used) additional correction for the weights of colors for different sub-pixels in a Bayer cell
		1.0,   // weightScaleB;
//		2.0,   // sigma;
		0.5,   // gamma;        -  Gaussian sigma to low-pass color components when calculating "smooth" color
		0.003, // minLin;
		0.299, // kr;
		0.114, // kb;
		2.5,   // saturationRed;                              
		2.5,   // saturationBlue;
		true,  // useFirstY;
        3.0,   // maskSigma,        3);
	    0.1,  // maskMin,        3);
	    0.35,  //maskMax,        3);
		true,  // combineWithSharpnessMask, // combine chroma mask with sharpness mask to reduce color leak through borders
	    2.0,   //chromaBrightSigma,        3);
	    10.0    //chromaDarkSigma,        3);
		);
   
//ColorCalibrationParameters
   public static ColorCalibParameters COLOR_CALIB_PARAMETERS= new ColorCalibParameters(
 			1.05, // gain_11,
  			1.05, // gain_12,
  			1.04, // gain_13,
  			1.05, // gain_21,
  			1.00, // gain_22,
  			0.97, // gain_23,
  			2.05, // gain_31,
  			1.10, // gain_32,
  			1.15, // gain_33,
  			0.97, // balanceRed_11,
  			1.02, // balanceRed_12,
  			1.04, // balanceRed_13,
  			1.09, // balanceRed_21,
  			0.98, // balanceRed_22,
  			0.98, // balanceRed_23,
  			1.10, // balanceRed_31,
  			1.005,// balanceRed_32,
  			0.96, // balanceRed_33,
  			0.95, // balanceBlue_11,
  			1.0,  // balanceBlue_12,
  			0.97, // balanceBlue_13,
  			1.0,  // balanceBlue_21,
  			0.96, // balanceBlue_22,
  			1.0,  // balanceBlue_23,
  			1.25, // balanceBlue_31,
  			1.0,  // balanceBlue_32,
  			0.96  // balanceBlue_33){
		   );
   public static RGBParameters RGB_PARAMETERS = new  RGBParameters (
		0.075, // r_min;
		0.075, //  g_min;
		0.075, //  b_min;
		1.0, //  r_max;
		1.0, //  g_max;
		1.0,  //  b_max;
		0.0, //  alpha_min;
		1.0 //  alpha_max;
   );
   
   public static FilesParameters FILE_PARAMETERS= new FilesParameters (
		   "RPSF_11.tiff",    "RPSF_12.tiff",    "RPSF_13.tiff",
		   "RPSF_21.tiff",    "RPSF_22.tiff",    "RPSF_23.tiff",
		   "RPSF_31.tiff",    "RPSF_32.tiff",    "RPSF_33.tiff",
		   "GAUSSIAN_11.tiff","GAUSSIAN_12.tiff","GAUSSIAN_13.tiff",
		   "GAUSSIAN_21.tiff","GAUSSIAN_22.tiff","GAUSSIAN_23.tiff",
		   "GAUSSIAN_31.tiff","GAUSSIAN_32.tiff","GAUSSIAN_33.tiff",
		   "", //./kernels",
		   "", //./results"
		   true //useXML
   );

   public static ProcessParameters PROCESS_PARAMETERS= new ProcessParameters (
		   true, //eyesisMode
		   true,true,true, // first camera channel
		   true,true,true, // second camera channel
		   true,true,true, // third camera channel
		    true,  // selectFile
			false, // thisFileOnly;
			1,    //  subChannelToProcess=subChannelToProcess (1..3)
			true, // split;
			true, // debayer;
			false,// showDebayerEnergy;
			true, // saveDebayerEnergy;
			true, // deconvolve;
			true, // combine;
			false,// showDenoiseMask;
			true, // saveDenoiseMask;
			false,// showChromaDenoiseMask;
			true, // saveChromaDenoiseMask;
			false,// showNoiseGains;
			false,// saveNoiseGains;
			true, // colorProc;
			true, // toRGB;
			true, // rotate;
			true, // crop - crop image to the sensor size
			true, // jpeg - convert to 8-bit RGB and save jpeg (if save is true)
			true, // save - save result
			false, //true, // save16 - save 16-bit tiff also if the end result is 8 bit
			false, //true, // save32 - save 32-bit tiff also if the end result is 8 or 16 bit
			true, // show;
			95,   // JPEG_quality
			0.5,   // JPEG_scale
			true
   );   
   public static int     CONVOLVE_FFT_SIZE= 128; // FFT size for sliding convolution with kernel
   public static int     THREADS_MAX=      100; // testing multi-threading, limit maximal number of threads

   public double GAUSS_WIDTH=0.4; //0 - use Hamming window

/* replace */
   public static int     PSF_SUBPIXEL_SHOULD_BE_4=4;         // sub-pixel decimation 

//   public float [][] OUT_PIXELS=null; // (global, used with threads to accumulate output result
   public double []  DEBAYER_ENERGY=null; // (global, used with threads to accumulate output result)
   public int        DEBAYER_ENERGY_WIDTH; // width of the DEBAYER_ENERGY image 

   public double []  DENOISE_MASK=null; // (global, used to return denoise mask to save/show
   public int        DENOISE_MASK_WIDTH; // width of the DENOISE_MASK image 
 
   public double []  DENOISE_MASK_CHROMA=null; // (global, used to return denoise mask to save/show
   public int        DENOISE_MASK_CHROMA_WIDTH; // width of the DENOISE_MASK_CHROMA image 

   
   public double []  MASK_LOHIRES=null; // (global, used with threads to accumulate output result)
//   public float []  MASK_LOHIRES=null; // (global, used with threads to accumulate output result)


//   public static String [] stackColorNames={"red","green","blue"};
   public static String [] stackColorNames= {"Red","Green","Blue"};
  
   public static ImageStack convolutionKernelStack=null; // select to use for image convolution
   public static ImageStack convolutionKernelStack2=null; // select to use for image convolution
   public static ImagePlus imp_gaussian=null;
   public static String DEFAULT_DIRECTORY=null;
   public static boolean ADVANCED_MODE=true; // show all buttons

   public Aberration_Correction() {
	   super("Aberration_Correction");
	   if (IJ.versionLessThan("1.43q")) return;
	   if (instance!=null) {
		   instance.toFront();
		   return;
	   }
	   instance = this;
	   addKeyListener(IJ.getInstance());
	   //    setLayout(new FlowLayout());
	   //    setLayout(new GridLayout(4, 1, 50, 5));
	   setLayout(new GridLayout(ADVANCED_MODE?7:3, 1));

	   panel6 = new Panel();
	   panel6.setLayout(new GridLayout(1, 0, 5, 5));
	   addButton("Process",panel6);
	   addButton("Save",panel6);
	   addButton("Restore",panel6);
	   add(panel6);
	   
	   panel5 = new Panel();
	   panel5.setLayout(new GridLayout(1, 0, 5, 5));
	   addButton("Kernels directory",panel5);
	   addButton("Results directory",panel5);
	   addButton("Configure paths",panel5);
	   addButton("Configure spilt",    panel5);
	   addButton("Configure demosaic", panel5);
	   add(panel5);
	   panel5a = new Panel();
	   panel5a.setLayout(new GridLayout(1, 0, 5, 5));
	   addButton("Configure convolution", panel5a);
	   addButton("Configure denoise", panel5a);
	   addButton("Configure color", panel5a);
	   addButton("Channel gains", panel5a);
	   addButton("Configure RGB", panel5a);
	   add(panel5a);


	   // Debug/development options 

	   if (ADVANCED_MODE) {
		   panel1 = new Panel();
		   panel1.setLayout(new GridLayout(1, 0, 5, 5)); // rows, columns, vgap, hgap
		   addButton("Split Image",panel1);
		   addButton("Debayer Image",panel1);
		   add(panel1);

		   panel2 = new Panel();
		   panel2.setLayout(new GridLayout(1, 0, 5, 5));
		   addButton("Select kernel stack",panel2);
		   addButton("Convolve with stack",panel2);
		   add(panel2);

		   panel3 = new Panel();
		   panel3.setLayout(new GridLayout(1, 0, 5, 5));
		   addButton("Select lo-res",panel3);
		   addButton("Combine pair",panel3);
		   addButton("Colors",panel3);
		   addButton("RGB",panel3);
		   add(panel3);


		   panel4 = new Panel();
		   panel4.setLayout(new GridLayout(1, 0, 5, 5));
		   addButton("Test",panel4);
		   addButton("Test Debayer",panel4);
		   add(panel4);
	   }
	   pack();

	   GUI.center(this);
	   setVisible(true);
	   FHT_INSTANCE=       new DoubleFHT();
	   SDFA_INSTANCE=      new ShowDoubleFloatArrays();
   }
   void addButton(String label, Panel panel) {
    Button b = new Button(label);
    b.addActionListener(this);
    b.addKeyListener(IJ.getInstance());
    panel.add(b);
  }
  public void processWindowEvent(WindowEvent e) {
    super.processWindowEvent(e);
    if (e.getID()==WindowEvent.WINDOW_CLOSING) {
      instance = null;	
    }
  }

  public void actionPerformed(ActionEvent e) {
    int i, j; //,l,iq;
    String label = e.getActionCommand();

    if (label==null) return;
/* ======================================================================== */
    if (label.equals("Configure spilt")) {
      showSplitBayerToStackDialog(SPLIT_PARAMETERS);
      return;
    }  
/* ======================================================================== */
    if (label.equals("Configure demosaic")) {
        showDeBayerDialog(DEBAYER_PARAMETERS, PROCESS_PARAMETERS);
        return;
    }  
/* ======================================================================== */
    if (label.equals("Configure convolution")) {
    	showStackConvolutionDialog();
        return;
    }  
/* ======================================================================== */
    if (label.equals("Configure denoise")) {
    	showCombinePairDialog(NONLIN_PARAMETERS, PROCESS_PARAMETERS);
        return;
    }  
/* ======================================================================== */
    if (label.equals("Configure color")) {
    	showColorProcessDialog(COLOR_PROC_PARAMETERS);
        return;
    } 
/* ======================================================================== */
    if (label.equals("Channel gains")) {
    	showColorCalibDialog(COLOR_CALIB_PARAMETERS);
        return;
    } 
/* ======================================================================== */
    if (label.equals("Configure RGB")) {
    	showRGBProcessDialog(RGB_PARAMETERS);
        return;
    }  
    
//    
/* ======================================================================== */
    if (label.equals("Split Image")) {
      if (!showSplitBayerToStackDialog(SPLIT_PARAMETERS)) return;
      DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
      ImagePlus imp_src = WindowManager.getCurrentImage();
      if (imp_src==null){
        IJ.showMessage("Error","Bayer Image required");
        return;
      }
      ImageStack sourceStack= bayerToStack(imp_src, // source Bayer image, linearized, 32-bit (float))
    		                               SPLIT_PARAMETERS);
      ImagePlus imp_srcStack = new ImagePlus(imp_src.getTitle()+"-EXP", sourceStack);
      imp_srcStack.getProcessor().resetMinAndMax();
      imp_srcStack.show();
      return;
/* ======================================================================== */
    } else if (label.equals("Debayer Image")) {
      if (! showDeBayerDialog(DEBAYER_PARAMETERS, PROCESS_PARAMETERS)) return;
      DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
      ImagePlus imp_src = WindowManager.getCurrentImage();
      if ((imp_src==null) || (imp_src.getStackSize()<3)){
        IJ.showMessage("Error","Bayer image stack required");
        return;
      }
    
      ImageStack imageStack= aliasScissorsStack(imp_src.getStack(),  // stack with 3 colors/slices with the image
    		                                    DEBAYER_PARAMETERS,
    		                                    PROCESS_PARAMETERS.showDebayerEnergy,
                                                THREADS_MAX, // number of image pixels/ sensor pixels (each direction) == 2
                                                UPDATE_STATUS);// update status info
	  if (PROCESS_PARAMETERS.showDebayerEnergy) {
		  SDFA_INSTANCE.showArrays (DEBAYER_ENERGY,DEBAYER_ENERGY_WIDTH,DEBAYER_ENERGY.length/DEBAYER_ENERGY_WIDTH, "Debayer-Energy");
	  }

      ImagePlus imp_debyrStack = new ImagePlus(imp_src.getTitle()+"-DeBayer", imageStack);
      imp_debyrStack.getProcessor().resetMinAndMax();
      imp_debyrStack.show();
      return;

/* ======================================================================== */
    } else if (label.equals("Select kernel stack")) {
//      int loop_debug_level=1;
      DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
      ImagePlus imp_kernels = WindowManager.getCurrentImage();
      if (imp_kernels==null){
        IJ.showMessage("Error","There is no kernel stack to process");
        return;
      }
      if (imp_kernels.getStackSize()<3) {
        IJ.showMessage("Error","Need a 3-layer stack with kernels");
        return;
      }
      convolutionKernelStack=imp_kernels.getStack();
      System.out.println("Select kernel stack "+imp_kernels.getTitle());
      return;

/* ======================================================================== */
    } else if (label.equals("Convolve with stack")) {
      if (convolutionKernelStack==null) {
        IJ.showMessage("Error","Convolution kernel stack needed, please select one with 'Select kernel stack' command");
        return;
      }
//      int loop_debug_level=1;
      if (!showStackConvolutionDialog()) return;
      DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
      ImagePlus imp_src = WindowManager.getCurrentImage();
      if (imp_src==null){
        IJ.showMessage("Error","There is no kernel stack to process");
        return;
      }
      if (imp_src.getStackSize()<3) {
        IJ.showMessage("Error","Need a 3-layer stack with kernels");
        return;
      }
      ImageStack convolvedStack= convolveStackWithKernelStack (imp_src.getStack(),  // stack with 3 colors/slices with the image
                                                           convolutionKernelStack, // stack with 3 colors/slices convolution kernels
                                                                CONVOLVE_FFT_SIZE, // 128 - fft size, kernel size should be size/2 
                                                                      THREADS_MAX,
                                                                    UPDATE_STATUS); // update status info

      ImagePlus imp_convolvedStack = new ImagePlus(imp_src.getTitle()+"-convolved", convolvedStack);
      imp_convolvedStack.getProcessor().resetMinAndMax();
      imp_convolvedStack.show();
      return;
/* ======================================================================== */
    } else if (label.equals("Select lo-res")) {
      imp_gaussian=WindowManager.getCurrentImage();
      if ((imp_gaussian==null) || (imp_gaussian.getStackSize()<3)) {
        IJ.showMessage("Error","Please select image stack of 3 (float) slices (r,g,b), resulted from convolving input bayer image with Gaussian kernels");
        return;
      }
      if (DEBUG_LEVEL>1) System.out.println ( "Read image "+imp_gaussian.getTitle()+" as the image convolved with Gaussian kernels");

      return;
/* ======================================================================== */
    } else if (label.equals("Combine pair")) {
      if ((imp_gaussian==null) || (imp_gaussian.getStackSize()<3)) {
        IJ.showMessage("Error","Wrong/empty Gaussian convolved image, please use 'Read Gaussian' command first");
        return;
      }
      ImagePlus imp_convolved=WindowManager.getCurrentImage();
      if ((imp_convolved==null) || (imp_convolved.getStackSize()<3)) {
        IJ.showMessage("Error","Please select image stack of 3 (float) slices (r,g,b), resulted from convolving input bayer image with the reversed PSF kernels");
        return;
      }
      if (!showCombinePairDialog(NONLIN_PARAMETERS, PROCESS_PARAMETERS)) return;
      DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
      if (DEBUG_LEVEL>1) System.out.println ( "Read image "+imp_convolved.getTitle()+" as the image convolved with the reversed PSF kernels");
      NONLIN_PARAMETERS.showMask=PROCESS_PARAMETERS.showDenoiseMask;
      ImageStack stack_combo=  combineLoHiStacks(imp_convolved.getStack(), // ImageStack with the image, convolved with the reversed PSF (sharp but with high noise)
                                                  imp_gaussian.getStack(),  // ImageStack with the image, convolved with the Gaussian (just lateral compensated)  (blurred, but low noise)
                                                  -1, // do not scale filtMin and filtMax for individual channel/
                                                  -1,
                                                  NONLIN_PARAMETERS, // show mask generated and used
                                                  null, //    noiseMask, // 2-d array of kernelsNoiseGain (divide mask by it)
                                                  32,   //    noiseStep, // linear pixels per noiseMask pixels (32)
                                                  
                    							  THREADS_MAX,
                    							  UPDATE_STATUS); // update status info
                                                  
	  if (PROCESS_PARAMETERS.showDenoiseMask) {
		  SDFA_INSTANCE.showArrays (DENOISE_MASK,DENOISE_MASK_WIDTH,DENOISE_MASK.length/DENOISE_MASK_WIDTH, "mask");
	  }

      ImagePlus imp_stack_combo = new ImagePlus(imp_convolved.getTitle()+"combo-rgb", stack_combo);
      imp_stack_combo.getProcessor().resetMinAndMax();
      imp_stack_combo.show();
      return;
/* ======================================================================== */
    } else if (label.equals("Colors")) {

      ImagePlus imp_convolved=WindowManager.getCurrentImage();
      if ((imp_convolved==null) || (imp_convolved.getStackSize()<3)) {
        IJ.showMessage("Error","Please select image stack of 5 (float) slices (r,g,b) and 2 weight slices");
        return;
      }
      if (!showColorProcessDialog(COLOR_PROC_PARAMETERS)) return;


      DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
      if (DEBUG_LEVEL>1) System.out.println ( "Read image "+imp_convolved.getTitle()+" as the image convolved with the reversed PSF kernels");
      ImageStack stack_convolved=imp_convolved.getStack();

/* Set stack sequence to r-g-b */
// public static String [] stackColorNames={"red","green","blue"};

/* find number of the green channel - should be called "green", if none - use last */
     if (!fixSliceSequence (stack_convolved)) {
    	 return;     
     }
     processColorsWeights(stack_convolved,
                               255.0/PSF_SUBPIXEL_SHOULD_BE_4/PSF_SUBPIXEL_SHOULD_BE_4,
                               COLOR_PROC_PARAMETERS,
                               null, // no channel calibration
                               0,
                               0
                           );

//COLOR_PROC_PARAMETERS.useFirstY
     YPrPbToRGB(stack_convolved,
                       COLOR_PROC_PARAMETERS.kr,        // 0.299;
                       COLOR_PROC_PARAMETERS.kb,        // 0.114;
                COLOR_PROC_PARAMETERS.useFirstY?9:8,        //  int sliceY,
                              6, // int slicePr,
                              7// int slicePb
                           );

      imp_convolved.getProcessor().resetMinAndMax();
      imp_convolved.updateAndDraw();

      SDFA_INSTANCE.showImageStackThree(stack_convolved, imp_convolved.getTitle()+"restored_rgb");
      return;
/* ======================================================================== */
    } else if (label.equals("Test")) {
        IJ.showMessage("Error","Nothing here, just a place holder");
     return;
/* ======================================================================== */
    } else if (label.equals("Test Debayer")) {

/**TODO - use stacks? Parameter? */


      ImagePlus imp_debayer=WindowManager.getCurrentImage();
      if (imp_debayer==null) {
        IJ.showMessage("Error","Please select monochrome 32-bit image");
        return;
      }
      if (!showDeBayerDialog(DEBAYER_PARAMETERS, PROCESS_PARAMETERS)) return;
      DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
      if ((imp_debayer.getWidth()<DEBAYER_PARAMETERS.size) || (imp_debayer.getHeight()<DEBAYER_PARAMETERS.size)) {
        IJ.showMessage("Error","Please select monochrome image at least "+DEBAYER_PARAMETERS.size+"x"+DEBAYER_PARAMETERS.size+" pixels");
        return;
      }
      Roi roi_debayer= imp_debayer.getRoi();
      if (roi_debayer==null){
        imp_debayer.setRoi(0, 0, imp_debayer.getWidth(), imp_debayer.getHeight());
        roi_debayer= imp_debayer.getRoi();
      }
      Rectangle rroi=roi_debayer.getBounds();
/* center to the selection center */
      rroi.x= (rroi.x+rroi.width/2)- (DEBAYER_PARAMETERS.size/2);
      rroi.y= (rroi.y+rroi.height/2)-(DEBAYER_PARAMETERS.size/2);
      rroi.width=DEBAYER_PARAMETERS.size;
      rroi.height=DEBAYER_PARAMETERS.size;
      if (rroi.x<0) rroi.x=0;
      else if (rroi.x>(imp_debayer.getWidth()-DEBAYER_PARAMETERS.size)) rroi.x= imp_debayer.getWidth()-DEBAYER_PARAMETERS.size;
      if (rroi.y<0) rroi.y=0;
      else if (rroi.y>(imp_debayer.getHeight()-DEBAYER_PARAMETERS.size)) rroi.y= imp_debayer.getHeight()-DEBAYER_PARAMETERS.size;
      float [] fpixels_debayer=(float[]) imp_debayer.getProcessor().getPixels();
      double [][] pixels_debayer=new double [3][];
      pixels_debayer[0]=new double [DEBAYER_PARAMETERS.size*DEBAYER_PARAMETERS.size];
      int debayer_width=imp_debayer.getWidth();
      int debayer_base=rroi.y*debayer_width+rroi.x;;
      for (i=0;i<DEBAYER_PARAMETERS.size;i++) for (j=0;j<DEBAYER_PARAMETERS.size;j++) pixels_debayer[0][i*DEBAYER_PARAMETERS.size+j]=fpixels_debayer[debayer_base+i*debayer_width+j];
      pixels_debayer[1]=pixels_debayer[0].clone();
      pixels_debayer[2]=pixels_debayer[0].clone();
      SDFA_INSTANCE.showArrays(pixels_debayer[0],  "original");

      pixels_debayer= normalizeAndWindow (pixels_debayer, initWindowFunction((int)Math.sqrt(pixels_debayer[1].length)),false);
      pixels_debayer=extendFFTInputTo (pixels_debayer, DEBAYER_PARAMETERS.size);
      if (DEBUG_LEVEL>2)   SDFA_INSTANCE.showArrays(pixels_debayer[0],  "wnd-orig");
      for (i=0;i<DEBAYER_PARAMETERS.size;i++) for (j=0;j<DEBAYER_PARAMETERS.size;j++) {
        pixels_debayer[1][i*DEBAYER_PARAMETERS.size+j]*= ((((i % PSF_SUBPIXEL_SHOULD_BE_4)==0) && ((j % PSF_SUBPIXEL_SHOULD_BE_4)==0 )) || (((i % PSF_SUBPIXEL_SHOULD_BE_4)==(PSF_SUBPIXEL_SHOULD_BE_4/2)) && ((j % PSF_SUBPIXEL_SHOULD_BE_4)==(PSF_SUBPIXEL_SHOULD_BE_4/2) )))?1.0:0.0;
        pixels_debayer[0][i*DEBAYER_PARAMETERS.size+j]*=  (((i % PSF_SUBPIXEL_SHOULD_BE_4)==0) && ((j % PSF_SUBPIXEL_SHOULD_BE_4)==(PSF_SUBPIXEL_SHOULD_BE_4/2) ))?1.0:0.0;
        pixels_debayer[2][i*DEBAYER_PARAMETERS.size+j]*=  (((i % PSF_SUBPIXEL_SHOULD_BE_4)==(PSF_SUBPIXEL_SHOULD_BE_4/2)) && ((j % PSF_SUBPIXEL_SHOULD_BE_4)==0 ))?1.0:0.0;
      }
      if (DEBUG_LEVEL>1)   SDFA_INSTANCE.showArrays(pixels_debayer,  DEBAYER_PARAMETERS.debayerStacks, "mosaic");
/* swap quadrants and perform FHT */
      double [][] amps=null;
      if (DEBUG_LEVEL>1)  amps=new double[3][];
      for (i=0;i<3;i++) {
        FHT_INSTANCE.swapQuadrants(pixels_debayer[i]);
        FHT_INSTANCE.transform(pixels_debayer[i]);
        if (amps!=null) amps[i]=FHT_INSTANCE.calculateAmplitude(pixels_debayer[i]);

      }
      DeBayerScissors debayer_instance=new DeBayerScissors(DEBAYER_PARAMETERS.size, // size of the square array, center is at size/2, size/2, only top half+line will be used
                                                         DEBAYER_PARAMETERS.polarStep, // maximal step in pixels on the maxRadius for 1 angular step (i.e. 0.5)
                                                        DEBAYER_PARAMETERS.debayerRelativeWidthGreen, // result green mask mpy by scaled default (diamond)
                                                      DEBAYER_PARAMETERS.debayerRelativeWidthRedblue, // result red/blue mask mpy by scaled default (square)
                                                 DEBAYER_PARAMETERS.debayerRelativeWidthRedblueMain, // green mask when applied to red/blue, main (center)
                                               DEBAYER_PARAMETERS.debayerRelativeWidthRedblueClones); // green mask when applied to red/blue, clones 

      double [][] both_masks= debayer_instance.aliasScissors(pixels_debayer[1], // fht array for green, will be masked in-place
                                                             DEBAYER_PARAMETERS.debayerThreshold, // no high frequencies - use default uniform filter
                                                                 DEBAYER_PARAMETERS.debayerGamma, // power function applied to the amplitudes before generating spectral masks
                                                                 DEBAYER_PARAMETERS.debayerBonus, // for green mask - rays start radius from center, relative to distance to the first alias
                                                           DEBAYER_PARAMETERS.mainToAlias, // relative main/alias amplitudes to enable lixels (i.e. 0.5 means that if alias is >0.5*main, the pixel will be masked out)
                                                             DEBAYER_PARAMETERS.debayerMaskBlur, // for both masks  sigma for gaussian blur of the binary masks (<0 -do not use "scissors")
                                                          DEBAYER_PARAMETERS.debayerUseScissors, // use "scissors", if false - just apply "diamond" ands "square" with DEBAYER_PARAMETERS.debayerRelativeWidthGreen and DEBAYER_PARAMETERS.debayerRelativeWidthRedblue
                                                                   DEBUG_LEVEL); // internal debug level

      double [] green_mask=   both_masks[0];
      double [] red_blue_mask=both_masks[1];
        
      pixels_debayer[1]=FHT_INSTANCE.multiply(pixels_debayer[1],green_mask,false);
      pixels_debayer[0]=FHT_INSTANCE.multiply(pixels_debayer[0],red_blue_mask,false);
      pixels_debayer[2]=FHT_INSTANCE.multiply(pixels_debayer[2],red_blue_mask,false);

      FHT_INSTANCE.inverseTransform(pixels_debayer[1]);
      FHT_INSTANCE.swapQuadrants(pixels_debayer[1]);

      FHT_INSTANCE.inverseTransform(pixels_debayer[0]);
      FHT_INSTANCE.swapQuadrants(pixels_debayer[0]);

      FHT_INSTANCE.inverseTransform(pixels_debayer[2]);
      FHT_INSTANCE.swapQuadrants(pixels_debayer[2]);

      SDFA_INSTANCE.showArrays(pixels_debayer,  DEBAYER_PARAMETERS.debayerStacks, "filtered");

/* Swap quadrants in masks and display them */
      if (DEBUG_LEVEL>1){
        FHT_INSTANCE.swapQuadrants(green_mask);
        FHT_INSTANCE.swapQuadrants(red_blue_mask);
        SDFA_INSTANCE.showArrays(green_mask,    "G-mask");
        SDFA_INSTANCE.showArrays(red_blue_mask, "RB-mask");
        if (amps!=null) {
/**normalize amplitudes, apply gamma */
          double dmax=0.0;
          for (i=0;i<amps.length;i++) {
            for (j=0;j<amps[i].length;j++) if (amps[i][j]>dmax) dmax=amps[i][j];
            dmax=1.0/dmax;
            for (j=0;j<amps[i].length;j++) amps[i][j]= Math.pow(amps[i][j]*dmax,DEBAYER_PARAMETERS.debayerGamma);
          }
          SDFA_INSTANCE.showArrays(amps,  DEBAYER_PARAMETERS.debayerStacks,"PWR");
          for(i=0;i<amps[1].length;i++){
             amps[1][i]*=green_mask[i];
             amps[0][i]*=red_blue_mask[i];
             amps[2][i]*=red_blue_mask[i];
          }
          SDFA_INSTANCE.showArrays(amps,  DEBAYER_PARAMETERS.debayerStacks, "PWR-MSK");
       }
      }
      return;
/* ======================================================================== */
    } else if (label.equals("Configure paths")) {
    	if (!showFilesDialog(FILE_PARAMETERS,PROCESS_PARAMETERS)) return;
     return;
/* ======================================================================== */

    } else if (label.equals("Kernels directory")) {
    	String fileName= selectKernelsDirectory(FILE_PARAMETERS.kernelDirectory);
        if (fileName!=null) FILE_PARAMETERS.kernelDirectory=fileName;
     return;
/* ======================================================================== */
    } else if (label.equals("Results directory")) {
    	String fileName= selectResultsDirectory(FILE_PARAMETERS.resultsDirectory);
        if (fileName!=null) FILE_PARAMETERS.resultsDirectory=fileName;
     return;
     
//* ======================================================================== */
    } else if (label.equals("Source file(s)")) {
    	String [] fileNames= selectSourceFiles(FILE_PARAMETERS.sourceFiles);
        if (fileNames!=null) FILE_PARAMETERS.sourceFiles=fileNames;
     return;
     
/* ======================================================================== */
     
    } else if (label.equals("Save")) {
    	saveProperties(null,FILE_PARAMETERS.resultsDirectory,FILE_PARAMETERS.useXML, PROPERTIES);
    	return;
/* ======================================================================== */
        
    } else if (label.equals("Restore")) {
    	loadProperties(null,FILE_PARAMETERS.resultsDirectory,FILE_PARAMETERS.useXML, PROPERTIES);
    	return;
/* ======================================================================== */
     
    } else if (label.equals("Process")) {
    	if (JP4_INSTANCE==null) JP4_INSTANCE= new JP46_Reader_camera(false);
    	if (!showProcessDialog(PROCESS_PARAMETERS)) return;
        DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
    	String [] paths= FILE_PARAMETERS.sourceFiles;
    	if (PROCESS_PARAMETERS.selectFile || (paths==null) || ((paths.length==1) && (paths[0].length()==0))) {
    		paths=selectSourceFiles(paths);
    		if (paths==null) {
    	          IJ.showMessage("Error","No files selected");
    	          return;
    		}
        	FILE_PARAMETERS.sourceFiles=paths;
    	}
        File [] files =new File [FILE_PARAMETERS.sourceFiles.length];
        String [] names=   new String[files.length];
        String [] dirNames=new String[files.length];
        for (int numFile=0;numFile<files.length;numFile++){
        	files[numFile]= new File (FILE_PARAMETERS.sourceFiles[numFile]);
            names[numFile]=   files[numFile].getName();
            dirNames[numFile]=files[numFile].getParent()+Prefs.getFileSeparator();
            if (DEBUG_LEVEL>1) System.out.println(">> FILE_PARAMETERS.sourceFiles["+numFile+ "]="+FILE_PARAMETERS.sourceFiles[numFile]+
            		               " dirNames["+numFile+ "]="+dirNames[numFile]+
            		               " names["+numFile+ "]="+names[numFile]);
        }
/*
        if ((FILE_PARAMETERS.kernelDirectory==null) || (FILE_PARAMETERS.kernelDirectory.length()==0)){
        	String fileName= selectKernelsDirectory(FILE_PARAMETERS.kernelDirectory);
            if (fileName!=null) FILE_PARAMETERS.kernelDirectory=fileName;
        }
 */       	
    	processChannelFiles(
    			  false, //	saveResult, // free memory, return empty results if false
    			  JP4_INSTANCE,
    			  dirNames,
    			  names,
    			  PROCESS_PARAMETERS,
    			  SPLIT_PARAMETERS,
    			  DEBAYER_PARAMETERS,
    			  FILE_PARAMETERS,
    			  NONLIN_PARAMETERS,
    			  COLOR_PROC_PARAMETERS,
    			  COLOR_CALIB_PARAMETERS,
    			  RGB_PARAMETERS,
    			  CONVOLVE_FFT_SIZE, // FFT size for sliding convolution with kernel
    			  THREADS_MAX,// maximal number of threads
    			  UPDATE_STATUS); // update ImageJ status
//    	saveProperties(null,FILE_PARAMETERS.resultsDirectory,FILE_PARAMETERS.useXML, properties);
if (PROCESS_PARAMETERS.saveSettings) saveProperties(FILE_PARAMETERS.resultsDirectory+Prefs.getFileSeparator()+"settings",null,FILE_PARAMETERS.useXML, PROPERTIES);
    	return;
/* ======================================================================== */
    } else if (label.equals("RGB" )) {
        ImagePlus imp_colorStack=WindowManager.getCurrentImage();
        if ((imp_colorStack==null) || (imp_colorStack.getStackSize()<3)) {
          IJ.showMessage("Error","Please select image stack of 3 (float) slices (r,g,b)");
          return;
        }
        
        if (!showRGBProcessDialog(RGB_PARAMETERS)) return;

        ImageStack stack=imp_colorStack.getStack();
        fixSliceSequence (stack);
        ImageStack stack_crop=cropStack32(stack,SPLIT_PARAMETERS);
        ImageStack stack_rot=rotateStack32CW(stack_crop);
        ImageStack stack16=convertRGB32toRGB16Stack(
        	stack_rot,
      		RGB_PARAMETERS); 
        ImagePlus imp_stack16 = new ImagePlus(imp_colorStack.getTitle()+"-16b", stack16);
        CompositeImage compositeImage=convertToComposite(imp_stack16);
        compositeImage.show();
	    ImagePlus imp_RGB24=convertRGB48toRGB24(
	    		  stack16,
   	    		  imp_colorStack.getTitle()+"-RGB24",
	    		  0, 65536, // r range 0->0, 65536->256
	    		  0, 65536, // g range
	    		  0, 65536, // b range
	    		  0, 65536);// alpha range
	    imp_RGB24.setTitle(imp_colorStack.getTitle()+"rgb24");
	    imp_RGB24.show();
        return;
    }
/* ======================================================================== */
    DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
  }

/* ======================================================================== */
  private String [] selectSourceFiles(String [] defaultPaths) {
	  String []patterns={".jp4",".jp46",".tiff",".tif"};
	  return selectFiles(false, // save  
			  "Source file(s) selection", // title
			  "Select source file", // button
			  new MultipleExtensionsFileFilter(patterns,"JP4/TIFF files"), // filter
			  defaultPaths);
  }
  public String selectSourceFile(String defaultPath) {
	  String []patterns={".jp4",".jp46",".tiff",".tif"};
	  return selectFile(false, // save  
			  "Source file(s) selection", // title
			  "Select source file", // button
			  new MultipleExtensionsFileFilter(patterns,"JP4/TIFF files"), // filter
			  defaultPath);
  }
  private String selectResultsDirectory(String defaultPath) {
	  return  selectDirectory(
			  true, // save
			  "Results directory selection", // title
			  "Select results directory", // button
			  null, //FileFIlter
			  defaultPath);
  }
  private String selectKernelsDirectory(String defaultPath) {
	  return  selectDirectory(
			  false, // save
			  "Convolution kernels directory selection",
			  "Select kernels directory",
			  null, //FileFIlter
			  defaultPath);

  }
  public String selectDirectory(boolean save, String title, String button, FileFilter filter, String defaultPath) {
	  return selectDirectoryOrFile(save,true, title, button, filter, defaultPath);
  }
  public String [] selectDirectories(boolean save, String title, String button, FileFilter filter, String [] defaultPaths) {
	  return selectDirectoriesOrFiles(save,true, title, button, filter, defaultPaths);
  }
  public String selectFile(boolean save,  String title, String button, FileFilter filter, String defaultPath) {
	  return selectDirectoryOrFile(save,false, title, button, filter, defaultPath );
  }
  public String [] selectFiles(boolean save,  String title, String button, FileFilter filter, String [] defaultPaths) {
	  return selectDirectoriesOrFiles(save,false, title, button, filter, defaultPaths );
  }
  public String [] selectDirectoriesOrFiles(boolean save,
		  boolean directory,
		  String title,
		  String button,
		  FileFilter filter,
		  String [] defaultPaths) {
	  File dir=null;
	  String defaultPath=null;
	  File [] files=null;
	  int fileNum;
	  if ((defaultPaths!=null) && (defaultPaths.length>0)) {
		  File [] tfiles=new File [defaultPaths.length];
		  int nf=defaultPaths.length;
		  for (fileNum=0;fileNum<defaultPaths.length; fileNum++) {
			  tfiles[fileNum]=new File(defaultPaths[fileNum]);
			  if ((!tfiles[fileNum].exists()) ||(!tfiles[fileNum].isFile())) {
				  tfiles[fileNum]=null;
				  nf--;
			  }
		  }
		  files=new File[nf];
		  nf=0;
		  for (fileNum=0;fileNum<defaultPaths.length; fileNum++) if (tfiles[fileNum]!=null){
			  files[nf++]=tfiles[fileNum];
		  }
	  }
	  if ((defaultPaths!=null) && (defaultPaths.length>0) &&  (!defaultPaths[0].equals(""))) {
		  defaultPath=defaultPaths[0];
		  dir = new File(defaultPath);
	  }
	  if ((dir==null) || (!dir.exists())) {
		  if (DEFAULT_DIRECTORY!=null) {
			  defaultPath = DEFAULT_DIRECTORY; 
			  dir = new File(defaultPath);
		  }
	  }
	  if ((dir==null) || (!dir.exists())) {
		  defaultPath = OpenDialog.getDefaultDirectory();
		  if (defaultPath!=null) dir = new File(defaultPath);
	  }
	  if ((dir!=null) && (!dir.exists())) dir=null;
	  if ((dir!=null) && (!dir.isDirectory())){
		  dir=dir.getParentFile();
	  }
//getSelectedFiles

	  JFileChooser fc= new JFileChooser();
	  fc.setFileSelectionMode(directory?JFileChooser.DIRECTORIES_ONLY:JFileChooser.FILES_ONLY);
	  fc.setMultiSelectionEnabled(true);
	  if ((title!=null)  && (title.length()>0)) fc.setDialogTitle(title); 
	  if ((button!=null) && (button.length()>0)) fc.setApproveButtonText(button); 
	  if (filter!=null) fc.setFileFilter(filter) ; 
	  if (dir!=null) 	fc.setCurrentDirectory(dir);
	  fc.setSelectedFiles(files);
	  int returnVal = save?(fc.showSaveDialog(IJ.getInstance())):(fc.showOpenDialog(IJ.getInstance()));
	  if (returnVal!=JFileChooser.APPROVE_OPTION)	return null;
	  DEFAULT_DIRECTORY=fc.getCurrentDirectory().getPath();
	  files=fc.getSelectedFiles();
	  if (files.length<1) return null;
	  String [] filenames=new String[files.length];
//	  for (int nFile=0;nFile<files.length;nFile++) filenames[nFile]= files[nFile].getName();
	  for (int nFile=0;nFile<files.length;nFile++) filenames[nFile]= files[nFile].getPath();
	  return filenames;
  }
  public String selectDirectoryOrFile(boolean save,
		  boolean directory,
		  String title,
		  String button,
		  FileFilter filter,
		  String defaultPath) {
	  File dir=null;
	  if ((defaultPath!=null) &&  (!defaultPath.equals(""))) {
		  dir = new File(defaultPath);
	  }
	  if ((dir==null) || (!dir.exists())) {
		  if (DEFAULT_DIRECTORY!=null) {
			  defaultPath = DEFAULT_DIRECTORY; 
			  dir = new File(defaultPath);
		  }
	  }
	  if ((dir==null) || (!dir.exists())) {
		  defaultPath = OpenDialog.getDefaultDirectory();
		  if (defaultPath!=null) dir = new File(defaultPath);
	  }
	  if ((dir!=null) && (!dir.exists())) dir=null;
	  if ((dir!=null) && (!dir.isDirectory())){
		  dir=dir.getParentFile();
	  }


	  JFileChooser fc= new JFileChooser();
	  fc.setFileSelectionMode(directory?JFileChooser.DIRECTORIES_ONLY:JFileChooser.FILES_ONLY);
	  fc.setMultiSelectionEnabled(false);
	  if ((title!=null)  && (title.length()>0)) fc.setDialogTitle(title); 
	  if ((button!=null) && (button.length()>0)) fc.setApproveButtonText(button); 
	  if (filter!=null) fc.setFileFilter(filter) ; 
	  if (dir!=null) 	fc.setCurrentDirectory(dir);
	  int returnVal = save?(fc.showSaveDialog(IJ.getInstance())):(fc.showOpenDialog(IJ.getInstance()));
	  if (returnVal!=JFileChooser.APPROVE_OPTION)	return null;
	  DEFAULT_DIRECTORY=fc.getCurrentDirectory().getPath();
	  return fc.getSelectedFile().getPath();
  }
  class MultipleExtensionsFileFilter extends FileFilter {
	  protected String [] patterns;
	  protected String    description="JP4 files";
	  
	  public MultipleExtensionsFileFilter (String [] patterns,String description) {
		  this.description=description;
		  this.patterns=patterns.clone();
	  }
	  public MultipleExtensionsFileFilter (String [] patterns) {
		  this.patterns=patterns.clone();
	  }
	  public boolean accept (File file) {
		  int i;  
		  String name=file.getName();
		  if (file.isDirectory()) return true;
		  for (i=0;i<patterns.length;i++) {
			  if (name.toLowerCase().endsWith(patterns[i].toLowerCase())) return true;
		  }
		  return false; 
	  }
	  public String getDescription() {
		  return description;
	  }
  }

  
 /* ======================================================================== */
  
  
  private boolean fixSliceSequence (ImageStack stack){
	  int i,j;
	  int [] rgbNumbers= {0,0,0};
	  for (j=0;j<3;j++) {
		  for (i=1;i<=3;i++) if (stack.getSliceLabel(i).toLowerCase().equals(stackColorNames[j].toLowerCase())){
// fix case (capitalized)
//			  System.out.println ( "stack.getSliceLabel("+i+")="+stack.getSliceLabel(i));
//			  System.out.println ( "stackColorNames["+j+"]="+stackColorNames[j]);
			  stack.setSliceLabel(stackColorNames[j],i);
			  rgbNumbers[j]=i;
//			  System.out.println ( "rgbNumbers["+j+"]="+rgbNumbers[j]);
			  break;
		  }
	  }
	  if (DEBUG_LEVEL>2) {
		  System.out.println ( "Input file color slice numbers:");
		  System.out.println ( "  Red -   slice "+((rgbNumbers[0]>0)?rgbNumbers[0]:"missing"));
		  System.out.println ( "  Green - slice "+((rgbNumbers[1]>0)?rgbNumbers[1]:"missing"));
		  System.out.println ( "  Blue -  slice "+((rgbNumbers[2]>0)?rgbNumbers[2]:"missing"));
	  }

	  for (i=0;i<3;i++) if (rgbNumbers[i]<=0) {
		  System.out.println ( stackColorNames[i]+ "  slice is missing in the input file. Please check slice names");
		  return false;
	  }
	  while ((rgbNumbers[0]!=1) || (rgbNumbers[1]!=2) ||(rgbNumbers[2]!=3)) {
		  if      (rgbNumbers[0]==1) swapStackSlices(stack,2,3);
		  else if (rgbNumbers[2]==3) swapStackSlices(stack,1,2);
		  else                       swapStackSlices(stack,1,3);
		  for (j=0;j<3;j++) {
			  for (i=1;i<=3;i++) if (stack.getSliceLabel(i).equals(stackColorNames[j])){
				  rgbNumbers[j]=i;
				  break;
			  }
		  }
	  }
	  return true;
  }
/* ======================================================================== */
  public void swapStackSlices(ImageStack stack,
                              int slice1,
                              int slice2) {
    String label=stack.getSliceLabel(slice1);
    stack.setSliceLabel(stack.getSliceLabel(slice2), slice1);
    stack.setSliceLabel(label,                       slice2);
    Object pixels=stack.getPixels(slice1);
    stack.setPixels   (stack.getPixels(slice2),      slice1);
    stack.setPixels   (pixels,                       slice2);
  }
 

/* ======================================================================== */
   public ImageStack cropStack32(
		  ImageStack stack,
		  SplitParameters splitParameters) {
	  int size=stack.getSize();
	  int iWidth=stack.getWidth();
	  int height=stack.getHeight()-splitParameters.addTop-splitParameters.addBottom;
	  int width= stack.getWidth()-splitParameters.addLeft-splitParameters.addRight;
	  int length=width*height;
	  ImageStack stack_crop = new ImageStack(width,height);
	  int i,x,y,index,base;
	  float [] ipixels;
	  float [] opixels;
	  for (i=0;i<size;i++) {
		 ipixels= (float[])stack.getPixels(i+1);
		 opixels=new float [length];
		 index=0;
		 for (y=0;y< height;y++) {
			  base=iWidth*(y+splitParameters.addTop)+splitParameters.addLeft;
			  for (x=0;x<width;x++) opixels[index++]=ipixels[base++];
		  }
		 stack_crop.addSlice(stack.getSliceLabel(i+1), opixels);
	  }
	  return stack_crop;
  }

/* ======================================================================== */
  public ImageStack rotateStack32CW(
		  ImageStack stack) {
	  int size=stack.getSize();
	  int height=stack.getHeight();
	  int width= stack.getWidth();
	  int length=width*height;
	  ImageStack stack_rot = new ImageStack(height, width);
	  int i,x,y,index;
	  float [] ipixels;
	  float [] opixels;
	  for (i=0;i<size;i++) {
		 ipixels= (float[])stack.getPixels(i+1);
		 opixels=new float [length];
		 index=0;
		 for (x=0;x<width;x++)for(y=height-1;y>=0;y--) opixels[index++]=ipixels[y*width+x];
		 stack_rot.addSlice(stack.getSliceLabel(i+1), opixels);
	  }
	  return stack_rot;
	  
  }
  /* ======================================================================== */
  public ImagePlus cropImage32(
		  ImagePlus imp,
		  SplitParameters splitParameters) {
	  int iWidth=imp.getWidth();
	  int height=imp.getHeight()-splitParameters.addTop-splitParameters.addBottom;
	  int width= imp.getWidth()-splitParameters.addLeft-splitParameters.addRight;
	  int length=width*height;
	  int x,y,index,base;
	  float [] ipixels=(float[])imp.getProcessor().getPixels();
	  float [] opixels=new float [length];
	  index=0;
	  for (y=0;y< height;y++) {
		base=iWidth*(y+splitParameters.addTop)+splitParameters.addLeft;
		for (x=0;x<width;x++) opixels[index++]=ipixels[base++];
	  }
	  ImageProcessor ip_crop=new FloatProcessor(width,height,opixels,null);
	  ImagePlus imp_crop = new ImagePlus(imp.getTitle(),ip_crop); // same title?
	  return imp_crop;
 }

/* ======================================================================== */
 public ImagePlus rotateImage32CW(
		 ImagePlus imp) {
	  int width=imp.getWidth();
	  int height=imp.getHeight();
	  int length=width*height;
	  int x,y,index;
	  float [] ipixels=(float[])imp.getProcessor().getPixels();
	  float [] opixels=new float [length];
	  index=0;
	  for (x=0;x<width;x++)for(y=height-1;y>=0;y--) opixels[index++]=ipixels[y*width+x];
	  ImageProcessor ip_rot=new FloatProcessor(height,width,opixels,null);
	  ImagePlus imp_rot = new ImagePlus(imp.getTitle(),ip_rot); // same title?
	  return imp_rot;
 }


/* ======================================================================== */
  public CompositeImage convertToComposite( ImagePlus imp) {
//	  if (imp.isComposite()) return imp;
	  if (imp.isComposite()) return null;
	  if (imp.getNChannels()>1) {
		  return null; // number of channels should be just 1
	  }
	  int c = imp.getStackSize();
	  imp.setDimensions(c, 1, 1);
	  CompositeImage ci = new CompositeImage(imp, CompositeImage.COMPOSITE);
//	  ci.show();
//	  imp.hide();
	  return ci;
  }

/* ======================================================================== */
  public ImageStack convertRGB32toRGB16Stack(
		  ImageStack stack32,
		  RGBParameters rgbParameters) {
	  ImageStack stack16 = new ImageStack(stack32.getWidth(), stack32.getHeight());
	  int length=stack32.getWidth()*stack32.getHeight();
	  int i,j;
	  float [] fpixels;
	  short [] spixels;
	  double [] mins= {rgbParameters.r_min,rgbParameters.g_min,rgbParameters.b_min};
	  double [] maxs= {rgbParameters.r_max,rgbParameters.g_max,rgbParameters.b_max};
	  if (stack32.getSize()<3) return null;
	  double value;
	  double scale;
	  for (i=0;i<3;i++) {
		 fpixels= (float[])stack32.getPixels(i+1);
		 scale=65535.0/(maxs[i]-mins[i]);
		 spixels=new short [length];
		 for (j=0;j<length;j++) {
			 value=(fpixels[j]-mins[i])*scale;
			 if      (value<0.0) value=0.0;
			 else if (value>65535.0) value=65535.0;
			 spixels[j]=(short)(value+0.5);
		 }
		 stack16.addSlice(stack32.getSliceLabel(i+1), spixels);
	  }
	  return stack16;
	  
  }
  
  public ImagePlus convertRGB48toRGB24(
		  ImageStack stack16,
		  String title,
		  int r_min,
		  int r_max,
		  int g_min,
		  int g_max,
		  int b_min,
		  int b_max,
		  int alpha_min,
		  int alpha_max){
	  int [] mins= {r_min,g_min,b_min,alpha_min};
	  int [] maxs= {r_max,g_max,b_max,alpha_max};
      int i;
	  int length=stack16.getWidth()*stack16.getHeight();
	  int numSlices  = stack16.getSize();
	  if (numSlices > 4) numSlices = 4;
	  short [][] spixels=new short[numSlices][];
	  int [] pixels=new int[length];
	  int c,d;
	  double [] scale=new double[numSlices];
	  for (c = 0; c < numSlices; c++) {
		  scale[c]=256.0/(maxs[c]-mins[c]);
	  }
	  for (i = 0; i < numSlices; i++) spixels[i]= (short[])stack16.getPixels(i+1);
	  for (i = 0; i < length; i++) {
		  pixels[i]=0;
		  for (c=0;c < numSlices;c++) {
			  d=(int)(((spixels[c][i]& 0xffff)-mins[c])*scale[c]);
			  if (d > 255) d=255;
			  else if (d < 0) d=0;
			  pixels[i]= d | (pixels[i] << 8);
		  }
	  }
	  ColorProcessor cp=new ColorProcessor(stack16.getWidth(),stack16.getHeight());
	  cp.setPixels(pixels);
	  ImagePlus imp=new ImagePlus(title,cp);
	  return imp;
  }
/* ======================================================================== */
  public ImagePlus Image32toGreyRGB24(
		  ImagePlus imp){
	  int width=imp.getWidth();
	  int height=imp.getHeight();
	  int length=width*height;
	  int i;
	  float [] ipixels=(float[])imp.getProcessor().getPixels();
	  int [] pixels=new int[length];
      float min=ipixels[0];
      float max=ipixels[0];
      for (i=0;i<length;i++) {
    	  if (min>ipixels[i])min=ipixels[i];
    	  if (max<ipixels[i])max=ipixels[i];
      }
      double d= 256.0/(max-min);
      int c;
      for (i=0;i<length;i++) {
    	  c=(int)((ipixels[i]-min)*d);
    	  if (c>255) c=255;
    	  pixels[i]=c | (c <<8) | (c<<16);
      }
	  ColorProcessor cp=new ColorProcessor(width,height,pixels);
      ImagePlus imp_rgb=new ImagePlus(imp.getTitle(),cp);
	  return imp_rgb;
  }
/* ======================================================================== */

	public void saveProperties(
			String path,      // full path or null
			String directory, // use as default directory if path==null 
			boolean useXML,
			Properties properties ){
   	    String [] XMLPatterns= {".xml"};
   	    String [] confPatterns={".conf"};
   	    String [] patterns=useXML?XMLPatterns:confPatterns;
     if (path==null) {
	    path= selectFile(true, // save  
			  "Save configuration selection", // title
			  "Select configuration file", // button
			  new MultipleExtensionsFileFilter(patterns, (useXML?"XML ":"")+"Configuration files"), // filter
			  directory); // may be ""
     } else path+=patterns[0];
     if (path==null) return;
     setAllProperties(properties);
    
     OutputStream os;
	try {
		os = new FileOutputStream(path);
	} catch (FileNotFoundException e1) {
   	 IJ.showMessage("Error","Failed to open configuration file: "+path);
	 return;
	}
    if (useXML) {
         try {
     		properties.storeToXML(os,
     		 "last updated " + new java.util.Date(), "UTF8");
     		
     	 } catch (IOException e) {
         	 IJ.showMessage("Error","Failed to write XML configuration file: "+path);
         	 return;
     	 }
     } else {
         try {
      		properties.store(os,
      		 "last updated " + new java.util.Date());
      	 } catch (IOException e) {
          	 IJ.showMessage("Error","Failed to write configuration file: "+path);
          	 return;
      	 }
     }
     try {
		os.close();
	 } catch (IOException e) {
		// TODO Auto-generated catch block
		e.printStackTrace();
	 }
	 if (DEBUG_LEVEL>0) System.out.println("Configuration parameters are saved to "+path);

	}
/* ======================================================================== */
	public void loadProperties(
			String path,
			String directory,
			boolean useXML,
			Properties properties){
	    String [] XMLPatterns= {".xml"};
	    String [] confPatterns={".conf"};
	    String [] patterns=useXML?XMLPatterns:confPatterns;
	     if (path==null) {
	 	    path= selectFile(false, // save  
	 			  "Configuration file selection", // title
	 			  "Read configuration file", // button
	 			  new MultipleExtensionsFileFilter(patterns,(useXML?"XML ":"")+"Configuration files"), // filter
	 			  directory); // may be ""
	      }  else path+=patterns[0];
	     if (path==null) return;
	     InputStream is;
		try {
			is = new FileInputStream(path);
		} catch (FileNotFoundException e) {
        	 IJ.showMessage("Error","Failed to open configuration file: "+path);
         	 return;
		}

	     if (useXML) {
	         try {
	     		properties.loadFromXML(is);
	     		
	     	 } catch (IOException e) {
	         	 IJ.showMessage("Error","Failed to read XML configuration file: "+path);
	         	 return;
	     	 }
	     } else {
	         try {
	      		properties.load(is);
	      	 } catch (IOException e) {
	          	 IJ.showMessage("Error","Failed to read configuration file: "+path);
	          	 return;
	      	 }
	     }
	     try {
	 		is.close();
	 	 } catch (IOException e) {
	 		// TODO Auto-generated catch block
	 		e.printStackTrace();
	 	 }      
	     getAllProperties(properties);
		 if (DEBUG_LEVEL>0) System.out.println("Configuration parameters are restored from "+path);
	}
/* ======================================================================== */
    public void setAllProperties(Properties properties){
    	properties.setProperty("MASTER_DEBUG_LEVEL",MASTER_DEBUG_LEVEL+"");
    	properties.setProperty("UPDATE_STATUS",     UPDATE_STATUS+     "");
    	SPLIT_PARAMETERS.setProperties("SPLIT_PARAMETERS.", properties);
    	DEBAYER_PARAMETERS.setProperties("DEBAYER_PARAMETERS.", properties);
    	NONLIN_PARAMETERS.setProperties("NONLIN_PARAMETERS.", properties);
    	COLOR_PROC_PARAMETERS.setProperties("COLOR_PROC_PARAMETERS.", properties);
    	COLOR_CALIB_PARAMETERS.setProperties("COLOR_CALIB_PARAMETERS.", properties);
    	RGB_PARAMETERS.setProperties("RGB_PARAMETERS.", properties);
    	FILE_PARAMETERS.setProperties("FILE_PARAMETERS.", properties);
    	PROCESS_PARAMETERS.setProperties("PROCESS_PARAMETERS.", properties);
    	properties.setProperty("CONVOLVE_FFT_SIZE",CONVOLVE_FFT_SIZE+""); //128,  FFT size for sliding convolution with kernel
    	properties.setProperty("THREADS_MAX",THREADS_MAX+""); // 100, testing multi-threading, limit maximal number of threads
    	properties.setProperty("GAUSS_WIDTH",GAUSS_WIDTH+""); // 0.4 (0 - use Hamming window)
    	properties.setProperty("PSF_SUBPIXEL_SHOULD_BE_4",PSF_SUBPIXEL_SHOULD_BE_4+""); // 4, sub-pixel decimation 
    	
    }
/* ======================================================================== */
    public void getAllProperties(Properties properties){
       MASTER_DEBUG_LEVEL = Integer.parseInt(properties.getProperty("MASTER_DEBUG_LEVEL"));
       UPDATE_STATUS= Boolean.parseBoolean(properties.getProperty("UPDATE_STATUS"));
   	   SPLIT_PARAMETERS.getProperties("SPLIT_PARAMETERS.", properties);
       DEBAYER_PARAMETERS.getProperties("DEBAYER_PARAMETERS.", properties);
   	   NONLIN_PARAMETERS.getProperties("NONLIN_PARAMETERS.", properties);
       COLOR_PROC_PARAMETERS.getProperties("COLOR_PROC_PARAMETERS.", properties);
       COLOR_CALIB_PARAMETERS.getProperties("COLOR_CALIB_PARAMETERS.", properties);
	   RGB_PARAMETERS.getProperties("RGB_PARAMETERS.", properties);
   	   FILE_PARAMETERS.getProperties("FILE_PARAMETERS.", properties);
   	   PROCESS_PARAMETERS.getProperties("PROCESS_PARAMETERS.", properties);
	   CONVOLVE_FFT_SIZE=Integer.parseInt(properties.getProperty("CONVOLVE_FFT_SIZE"));
	   THREADS_MAX=Integer.parseInt(properties.getProperty("THREADS_MAX"));
	   GAUSS_WIDTH=Double.parseDouble(properties.getProperty("GAUSS_WIDTH"));
	   PSF_SUBPIXEL_SHOULD_BE_4=Integer.parseInt(properties.getProperty("PSF_SUBPIXEL_SHOULD_BE_4")); 
    }

/* ======================================================================== */
  public ImagePlus [][][] processChannelFiles(
		  boolean                     saveResult, // free memory, return empty results if false
		  JP46_Reader_camera        JP4_instance,
		  String []                 idirectories,
		  String []                   ifilenames,
		  ProcessParameters     processParameters,
		  SplitParameters         splitParameters,
		  DebayerParameters     debayerParameters,
		  FilesParameters         filesParameters,
		  NonlinParameters       nonlinParameters,
		  ColorProcParameters colorProcParameters,
		  ColorCalibParameters colorCalibParameters,
		  RGBParameters             rgbParameters,
		  int                     convolveFFTSize, // 128; // FFT size for sliding convolution with kernel
		  int                          threadsMax,     //   100; // testing multi-threading, limit maximal number of threads
		  boolean                    updateStatus){
	  Runtime runtime = Runtime.getRuntime();
	  long 	  startTime=System.nanoTime();
	  ImagePlus [][][] result=null;
	  ImagePlus [][]   ImageNoiseGains=null;
	  int nFile,nChn,nSubChn;
	  int numFiles=ifilenames.length;
	  if (numFiles==0) return null;
	  double [][] noiseMask=null;
	  for (nFile=0;nFile<ifilenames.length;nFile++) {
		  if (DEBUG_LEVEL>1) System.out.println (nFile+": directory="+idirectories[nFile]+" filename="+ ifilenames[nFile] );
	  }
	  boolean eyesisMode=processParameters.eyesisMode;
	  // turn eyesisMode off, if all the selected files do not have the proper signature	(*_1-*, *_2-*, *_3-*)
	  if (eyesisMode) {
		  boolean allEyesis=true;
		  boolean noneEyesis=true;
		  boolean thisEyesis;
		  for (nFile=0;nFile<numFiles;nFile++) {
			  thisEyesis=false;
			  for (nChn=0; nChn<processParameters.numEyesisChannels; nChn++) if  (ifilenames[nFile].indexOf("_"+(nChn+1)+"-")>=0) {
				  thisEyesis=true;
				  break;
			  }
			  if (thisEyesis) noneEyesis=false;
			  else            allEyesis=false;
		  }
		  if (!noneEyesis && !allEyesis) {
			  IJ.showMessage("Error","Some of the selected files are Eyesis (with _1-, _2-, _3-), some - not\n Please select only one type or turn off the eysisMode");
			  return null;
		  }
		  if (noneEyesis) {
			  if (DEBUG_LEVEL>0) System.out.println(">> eyesisMode was on, but non files have channel number, turning eyesisMode off.");
			  eyesisMode=false;
		  }
	  }
	  int numChannels=   eyesisMode?processParameters.numEyesisChannels:1;
	  int numSubChannels=eyesisMode?processParameters.numEyesisSubChannels:1;

	  // Add same shot files (if enabled) - duplicates are possible
	  String [][]filenames = new String [numFiles][numChannels];
	  boolean [][][] channels = new boolean [ifilenames.length][numChannels][numSubChannels];
	  String [] directories=idirectories.clone();
	  int channel;
	  int indexOfDot;
	  String extension;
	  String filename;
	  String [] prefix = new String[numFiles];
	  String [] suffix = new String[numFiles];
	  for (nFile=0;nFile<numFiles;nFile++) {
		  indexOfDot=ifilenames[nFile].toLowerCase().lastIndexOf(".");
		  if (indexOfDot<0) {
			  IJ.showMessage("Error","Valid source filename should have extension");
			  return null;
		  }
		  extension=ifilenames[nFile].substring(indexOfDot);
		  if (eyesisMode){
			  filename= ifilenames[nFile].substring(0,indexOfDot);
			  channel=-1;
			  for (nChn=1;nChn<=numChannels;nChn++) if (filename.indexOf("_"+nChn+"-")>=0) {
				  channel=nChn;
				  break;
			  }
			  if (channel<0) { // should not get here
				  return null;
			  }
			  prefix[nFile]=filename.substring(0,filename.indexOf("_"+channel+"-"));
			  suffix[nFile]=filename.substring(filename.indexOf("_"+channel+"-")+3);
			  if (DEBUG_LEVEL>1) {
				  System.out.println ("channel="+channel);
				  System.out.println ("prefix="+prefix[nFile]);
				  System.out.println ("suffix="+ suffix[nFile]);
			  }
			  for (nChn=0;nChn<numChannels;nChn++)  {
				  if ((!processParameters.thisFileOnly) || (nChn==(channel-1))) { 
					  filenames[nFile][nChn]=prefix[nFile]+"_"+(nChn+1)+"-"+suffix[nFile]+extension;
					  for (nSubChn=0;nSubChn<numSubChannels;nSubChn++) {
						  if ((!processParameters.thisFileOnly) || (nSubChn==(processParameters.subChannelToProcess-1))) 
							  channels[nFile][nChn][nSubChn]=true;  
						  else channels[nFile][nChn][nSubChn]=false;
					  }
				  } else{
					  filenames[nFile][nChn]=null;
					  channels[nFile][nChn]=null;
				  }
			  }		  
		  } else {
			  filenames[nFile][0]=ifilenames[nFile];
			  channels[nFile][0][0]=true;
			  prefix[nFile]= ifilenames[nFile].substring(0,indexOfDot);
			  suffix[nFile]="";
		  }
	  }
	  int nFile1;
	  for (nFile=1;nFile<numFiles;nFile++) {
		  for (nChn=0;(nChn<numChannels) && (filenames[nFile]!=null);nChn++) if (filenames[nFile][nChn]!=null) {
			  for (nFile1=0;nFile1<nFile;nFile1++)
				  if ((filenames[nFile1]!=null) &&
						  (filenames[nFile1][nChn]!=null) && 
						  (directories[nFile1].equals(directories[nFile]))&&
						  (filenames[nFile1][nChn].equals(filenames[nFile][nChn]))) {
					  filenames[nFile]=null;
					  //						  nf--;
					  break;
				  }
		  }

	  }
	  File file;
	  result=new ImagePlus [numFiles][numChannels][numSubChannels] ;
	  for (nFile=1;nFile<numFiles;nFile++)
		  for (nChn=0;nChn<numChannels;nChn++)
			  for (nSubChn=0;nSubChn<numSubChannels;nSubChn++) result[nFile][nChn][nSubChn]=null;

	  ImagePlus imp_composite=null;
	  ImageStack stack,stack_d,stack_g;
	  ImageStack [][] kernelsNoise=new ImageStack[numChannels][numSubChannels];
	  String title, titleFull;
	  String kernelPath;
	  ImagePlus imp_kernels= null;
	  ImagePlus imp_kernels2=null;
	  CompositeImage compositeImage=null;

	  for (nChn=0;nChn<numChannels;nChn++)
		  for (nSubChn=0;nSubChn<numSubChannels;nSubChn++) kernelsNoise[nChn][nSubChn]=null;
  
// == Calculate noise gains of convolution kernels ==
	  //processParameters.frames[i][j]
	  ImageNoiseGains=new ImagePlus [numChannels][numSubChannels] ;

	  if (processParameters.deconvolve && (nonlinParameters.noiseGainPower!=0)) for (nChn=0;nChn<numChannels;nChn++) for (nSubChn=0;nSubChn<numSubChannels;nSubChn++)
	      if (processParameters.frames[nChn][nSubChn]){
		  //Ask for the kernel directory if it is undefined
		  if ((filesParameters.kernelDirectory==null) || (filesParameters.kernelDirectory.length()==0)) {
			  kernelPath= selectKernelsDirectory(filesParameters.kernelDirectory);
			  if (kernelPath!=null) filesParameters.kernelDirectory=kernelPath;
		  }
		  // Read deconvolution kernels
		  kernelPath=filesParameters.kernelDirectory+Prefs.getFileSeparator()+filesParameters.rPSFNames[nChn][nSubChn];
		  file=new File(kernelPath);
		  if (!file.exists()) {
			  System.out.println("Kernel stack file "+kernelPath+" does not exist");
			  continue;
		  }
		  imp_kernels=new ImagePlus(kernelPath);
		  if (imp_kernels.getStackSize()<3) {
			  System.out.println("Need a 3-layer stack with kernels");
			  continue;
		  }
		  convolutionKernelStack=imp_kernels.getStack();
		  if (DEBUG_LEVEL>1) System.out.println("Using kernel stack "+kernelPath+" for noise gain calculation");
		  if (nonlinParameters.useDiffNoiseGains) {
			  // Read Gaussian kernels
			  kernelPath=filesParameters.kernelDirectory+Prefs.getFileSeparator()+filesParameters.gaussianNames[nChn][nSubChn];
			  file=new File(kernelPath);
			  if (!file.exists()) {
				  System.out.println("Gaussian stack file "+kernelPath+" does not exist");
				  continue;
			  }
			  imp_kernels2=new ImagePlus(kernelPath);
			  if (imp_kernels.getStackSize()<3) {
				  System.out.println("Need a 3-layer stack with gaussian kernels");
				  continue;
			  }
			  convolutionKernelStack2=imp_kernels2.getStack();
			  if (DEBUG_LEVEL>1) System.out.println("Using second (lowpass) stack "+kernelPath+" for noise gain calculation");
		  }  else convolutionKernelStack2=null; 
		  kernelsNoise[nChn][nSubChn]=calculateKernelsNoiseGains (
				  convolutionKernelStack,  // first stack with 3 colors/slices convolution kernels
				  convolutionKernelStack2, // second stack with 3 colors/slices convolution kernels (or null)
				  convolveFFTSize, // 128 - fft size, kernel size should be size/2
				  nonlinParameters.blurSigma,
				  threadsMax,
				  updateStatus); // update status info
		  convolutionKernelStack=null;
		  convolutionKernelStack2=null;
		  runtime.gc();
/*  ImageNoiseGains is used , maybe change below to use just the array	kernelsNoise[nChn][nSubChn] */	  
     	  title="noiseGains_"+(nonlinParameters.useDiffNoiseGains?"diff_":"")+(nChn+1)+""+(nSubChn+1);
		  ImageNoiseGains[nChn][nSubChn]= new ImagePlus(title, kernelsNoise[nChn][nSubChn]);
		  if (processParameters.saveNoiseGains || processParameters.showNoiseGains) {
//			  title="noiseGains_"+(nonlinParameters.useDiffNoiseGains?"diff_":"")+(nChn+1)+""+(nSubChn+1);
//			  ImageNoiseGains[nChn][nSubChn]= new ImagePlus(title, kernelsNoise[nChn][nSubChn]);
			  saveAndShow(ImageNoiseGains[nChn][nSubChn],
					  filesParameters,
					  processParameters.saveNoiseGains,
					  processParameters.showNoiseGains
			  );
//			  if (!processParameters.showNoiseGains) ImageNoiseGains[nChn][nSubChn]=null; // it is needed!
		  }		  
	  }  
	  convolutionKernelStack=null;
	  convolutionKernelStack2=null;
	  runtime.gc();
	  
	  
//==================================================	  
	  //	  int i,nFile,nChn,nSubChn;
	  for (nFile=0;nFile<numFiles;nFile++) if (filenames[nFile]!=null) 
		  for (nChn=0;nChn<numChannels;nChn++) if (filenames[nFile][nChn]!=null) {
			  if (DEBUG_LEVEL>1) System.out.println(">> using: "+directories[nFile]+filenames[nFile][nChn]+ ", eyesisMode="+eyesisMode);
			  file=new File(directories[nFile]+filenames[nFile][nChn]);
			  if (!file.exists()) continue;
			  String thisExtension=filenames[nFile][nChn].substring(filenames[nFile][nChn].toLowerCase().lastIndexOf("."));
			  if (thisExtension.equals(".tiff") || thisExtension.equals(".tif")){
				  String fullPath=directories[nFile]+Prefs.getFileSeparator()+filenames[nFile][nChn];
		        	Opener opener=new Opener();
		        	imp_composite=opener.openImage("", fullPath);
		        	if (imp_composite==null) {
		        		String msg="Failed to read map file "+fullPath;
		        		IJ.showMessage("Error",msg);
		        		throw new IllegalArgumentException (msg);
		        	}
		        	(new JP46_Reader_camera(false)).decodeProperiesFromInfo(imp_composite);
		        	imp_composite.setProperty("isTiff", true);
			  } else {
				  imp_composite=JP4_instance.open(
						  directories[nFile], // path,
						  filenames[nFile][nChn],
						  "",  //arg - not used in JP46 reader
						  true, // un-apply camera color gains
						  null, // new window
						  false); // do not show
				  if (imp_composite==null) continue;
			  }
			  for (nSubChn=0;nSubChn<numSubChannels;nSubChn++) if (channels[nFile][nChn][nSubChn])  if (processParameters.frames[nChn][nSubChn]) {
				  DENOISE_MASK=null; // color processing might try to use it, prevent using wrong one
				  // Extract sub-frame
				  if (eyesisMode) {
					  result[nFile][nChn][nSubChn]= JP4_instance.demuxImage(imp_composite, nSubChn);
					  if (result[nFile][nChn][nSubChn]==null) continue;
					  title=prefix[nFile]+"-"+suffix[nFile]+"__"+(nChn+1)+(nSubChn+1);
				  } else {
					  result[nFile][nChn][nSubChn]= imp_composite; // single channel
					  title=prefix[nFile]+suffix[nFile];
				  }
				  if (DEBUG_LEVEL>1) System.out.println("processing: "+title+ ", eyesisMode="+eyesisMode);
				  if (result[nFile][nChn][nSubChn]==null) continue;
				  result[nFile][nChn][nSubChn].setTitle(title+"RAW");
				  if (!processParameters.split){
					  saveAndShow(result[nFile][nChn][nSubChn], processParameters, filesParameters);
					  if(!saveResult) {
						  result[nFile][nChn][nSubChn]=null; // erase results to save memory
						  runtime.gc();
						  if (DEBUG_LEVEL>1) System.out.println("--- Free memory="+runtime.freeMemory()+" (of "+runtime.totalMemory()+")");
					  }
					  continue;
				  }
				  // Split into Bayer components, oversample, increase canvas    		  
				  stack= bayerToStack(
						  result[nFile][nChn][nSubChn], // source Bayer image, linearized, 32-bit (float))
						  splitParameters);
				  titleFull=title+"-SPLIT";
				  if (!processParameters.debayer) {
					  result[nFile][nChn][nSubChn]= new ImagePlus(titleFull, stack);    			  
					  saveAndShow(result[nFile][nChn][nSubChn], processParameters, filesParameters);
					  if(!saveResult) {
						  result[nFile][nChn][nSubChn]=null; // erase results to save memory
						  runtime.gc();
						  if (DEBUG_LEVEL>1) System.out.println("--- Free memory="+runtime.freeMemory()+" (of "+runtime.totalMemory()+")");
					  }
					  continue;
				  }
				  // Demosaic image
				  stack= aliasScissorsStack(stack,  // stack with 3 colors/slices with the image
						  debayerParameters,
						  (processParameters.saveDebayerEnergy || processParameters.showDebayerEnergy),
						  threadsMax, // number of image pixels/ sensor pixels (each direction) == 2
						  updateStatus);// update status info
				  if (processParameters.saveDebayerEnergy || processParameters.showDebayerEnergy) {
					  if (DEBAYER_ENERGY!=null) {
						  ImagePlus debayerMask=SDFA_INSTANCE.makeArrays (DEBAYER_ENERGY,
								  DEBAYER_ENERGY_WIDTH,
								  DEBAYER_ENERGY.length/DEBAYER_ENERGY_WIDTH,
								  title+"-DEBAYER-ENERGY");
						  saveAndShow(debayerMask,
								  filesParameters,
								  processParameters.saveDebayerEnergy,
								  processParameters.showDebayerEnergy
						  );
					  }
				  }
				  titleFull=title+"-DEMOSAIC";
				  result[nFile][nChn][nSubChn]= new ImagePlus(titleFull, stack);
				  if (processParameters.deconvolve) {
					  //Ask for the kernel directory if it is undefined
					  if ((filesParameters.kernelDirectory==null) || (filesParameters.kernelDirectory.length()==0)) {
						  kernelPath= selectKernelsDirectory(filesParameters.kernelDirectory);
						  if (kernelPath!=null) filesParameters.kernelDirectory=kernelPath;
					  }
					  // Read deconvolution kernels
					  kernelPath=filesParameters.kernelDirectory+Prefs.getFileSeparator()+filesParameters.rPSFNames[nChn][nSubChn];
					  file=new File(kernelPath);
					  if (!file.exists()) {
						  System.out.println("Kernel stack file "+kernelPath+" does not exist");
						  if(!saveResult) {
							  if(!saveResult) {
								  result[nFile][nChn][nSubChn]=null; // erase results to save memory
								  runtime.gc();
								  if (DEBUG_LEVEL>1) System.out.println("--- Free memory="+runtime.freeMemory()+" (of "+runtime.totalMemory()+")");
							  }
							  continue;
						  }
					  }
					  imp_kernels=new ImagePlus(kernelPath);
					  if (imp_kernels.getStackSize()<3) {
						  System.out.println("Need a 3-layer stack with kernels");
						  if(!saveResult) {
							  result[nFile][nChn][nSubChn]=null; // erase results to save memory
							  runtime.gc();
							  if (DEBUG_LEVEL>1) System.out.println("--- Free memory="+runtime.freeMemory()+" (of "+runtime.totalMemory()+")");
						  }
						  continue;
					  }
					  convolutionKernelStack=imp_kernels.getStack();
					  if (DEBUG_LEVEL>1) System.out.println("Using kernel stack "+kernelPath+" for convolution with "+result[nFile][nChn][nSubChn].getTitle());
					  stack_d= convolveStackWithKernelStack(
							  stack,  // stack with 3 colors/slices with the image
							  convolutionKernelStack, // stack with 3 colors/slices convolution kernels
							  convolveFFTSize, // 128 - fft size, kernel size should be size/2 
							  threadsMax,
							  updateStatus); // update status info
					  titleFull=title+"-DECONV";
					  if (processParameters.combine) {
						  // Read Gaussian kernels
						  kernelPath=filesParameters.kernelDirectory+Prefs.getFileSeparator()+filesParameters.gaussianNames[nChn][nSubChn];
						  file=new File(kernelPath);
						  if (!file.exists()) {
							  System.out.println("Gaussian stack file "+kernelPath+" does not exist");
							  if(!saveResult) {
								  result[nFile][nChn][nSubChn]=null; // erase results to save memory
								  runtime.gc();
								  if (DEBUG_LEVEL>1) System.out.println("--- Free memory="+runtime.freeMemory()+" (of "+runtime.totalMemory()+")");
							  }
							  continue;
						  }
						  imp_kernels=new ImagePlus(kernelPath);
						  if (imp_kernels.getStackSize()<3) {
							  System.out.println("Need a 3-layer stack with gaussian kernels");
							  if(!saveResult)  result[nFile][nChn][nSubChn]=null; // erase results to save memory
							  continue;
						  }
						  convolutionKernelStack=imp_kernels.getStack();
						  if (DEBUG_LEVEL>1) System.out.println("Using Gaussian stack "+kernelPath+" for convolution with "+result[nFile][nChn][nSubChn].getTitle());
						  stack_g= convolveStackWithKernelStack(
								  stack,  // stack with 3 colors/slices with the image
								  convolutionKernelStack, // stack with 3 colors/slices convolution kernels
								  convolveFFTSize, // 128 - fft size, kernel size should be size/2 
								  threadsMax,
								  updateStatus); // update status info
						  // Combine Gaussian and Deconvolved
						  noiseMask=extractNoiseMask(
								     ImageNoiseGains[nChn][nSubChn],// contains 3-slice stack (r,b,g)
								     nonlinParameters.noiseGainWeights[0], // coefficient for slice 0 (r)
								     nonlinParameters.noiseGainWeights[1], // coefficient for slice 1 (b)
								     nonlinParameters.noiseGainWeights[2], // coefficient for slice 2 (g)
								     1,     // decimate result (not yet supported)
								     nonlinParameters.noiseGainPower
								  );
// show noise mask here?						  
						  nonlinParameters.showMask=processParameters.showDenoiseMask;
//				          if (DEBUG_LEVEL>1) System.out.println ( " noiseMask.length="+((noiseMask==null)?"null":(noiseMask.length+" noiseMask[0].length="+noiseMask[0].length)));
						  
						  stack=  combineLoHiStacks(stack_d, // ImageStack with the image, convolved with the reversed PSF (sharp but with high noise)
								  stack_g,  // ImageStack with the image, convolved with the Gaussian (just lateral compensated)  (blurred, but low noise)
								  nChn,
								  nSubChn,
								  nonlinParameters, // show mask generated and used
                                  noiseMask, // 2-d array of kernelsNoiseGain (divide mask by it)
                                  32,        // linear pixels per noiseMask pixels (32)
								  threadsMax,
								  updateStatus); // update status info
						  if (processParameters.saveDenoiseMask || processParameters.showDenoiseMask) {
							  ImagePlus denoiseMask=SDFA_INSTANCE.makeArrays (DENOISE_MASK,
									  DENOISE_MASK_WIDTH,
									  DENOISE_MASK.length/DENOISE_MASK_WIDTH,
									  title+"-MASK");
							  if (processParameters.jpeg) {
// crop Mask to original image size
								  if (processParameters.crop){
									  denoiseMask=cropImage32(denoiseMask,splitParameters);
								  }
// rotate the result			  
								  if (processParameters.rotate){
									  denoiseMask=rotateImage32CW(denoiseMask);
								  }
// scale the result
								  if (processParameters.JPEG_scale!=1.0){
									  ImageProcessor ip=denoiseMask.getProcessor();
									  ip.setInterpolationMethod(ImageProcessor.BICUBIC);
									  ip=ip.resize((int)(ip.getWidth()*processParameters.JPEG_scale),(int) (ip.getHeight()*processParameters.JPEG_scale));
									  denoiseMask= new ImagePlus(denoiseMask.getTitle(),ip);
									  denoiseMask.updateAndDraw();
								  }
								  if (processParameters.showDenoiseMask) denoiseMask.show(); 
// public ImagePlus Image32toGreyRGB24(ImagePlus imp);
								  if (processParameters.saveDenoiseMask) {
									  ImagePlus denoiseMaskRGB24=Image32toGreyRGB24(denoiseMask);
								    saveAndShow(denoiseMaskRGB24,
									  filesParameters,
									  processParameters.saveDenoiseMask,
									  false, //processParameters.showDenoiseMask,
									  processParameters.JPEG_quality);
								    denoiseMaskRGB24=null;
								  }
							  } else {
							    saveAndShow(denoiseMask,
									  filesParameters,
									  processParameters.saveDenoiseMask,
									  processParameters.showDenoiseMask
							    );
							  }
						  }
						  stack_g=null;
						  stack_d=null;
						  titleFull=title+"-COMBO";
					  } else { // end of if (processParameters.combine)
						  stack=stack_d;
					  } // end of else if (processParameters.combine)
				  }  else if (processParameters.combine) { // just create convolution with Gaussians
					  // Read Gaussian kernels
					  kernelPath=filesParameters.kernelDirectory+Prefs.getFileSeparator()+filesParameters.gaussianNames[nChn][nSubChn];
					  file=new File(kernelPath);
					  if (!file.exists()) {
						  System.out.println("Gaussian stack file "+kernelPath+" does not exist");
						  if(!saveResult) {
							  result[nFile][nChn][nSubChn]=null; // erase results to save memory
							  runtime.gc();
							  if (DEBUG_LEVEL>1) System.out.println("--- Free memory="+runtime.freeMemory()+" (of "+runtime.totalMemory()+")");
						  }
						  continue;
					  }
					  imp_kernels=new ImagePlus(kernelPath);
					  if (imp_kernels.getStackSize()<3) {
						  System.out.println("Need a 3-layer stack with gaussian kernels");
						  if(!saveResult) {
							  result[nFile][nChn][nSubChn]=null; // erase results to save memory
							  runtime.gc();
							  if (DEBUG_LEVEL>1) System.out.println("--- Free memory="+runtime.freeMemory()+" (of "+runtime.totalMemory()+")");
						  }
						  continue;
					  }
					  convolutionKernelStack=imp_kernels.getStack();
					  if (DEBUG_LEVEL>1) System.out.println("Using Gaussian stack "+kernelPath+" for convolution with "+result[nFile][nChn][nSubChn].getTitle());
					  stack_g= convolveStackWithKernelStack(
							  stack,  // stack with 3 colors/slices with the image
							  convolutionKernelStack, // stack with 3 colors/slices convolution kernels
							  convolveFFTSize, // 128 - fft size, kernel size should be size/2 
							  threadsMax,
							  updateStatus); // update status info
					  titleFull=title+"-LOWRES";
				  }// end of if (processParameters.deconvolve)
				  //stack now has the result, titleFull - correct title for the image 
				  if (!processParameters.colorProc){
					  result[nFile][nChn][nSubChn]= new ImagePlus(titleFull, stack);    			  
					  saveAndShow(result[nFile][nChn][nSubChn], processParameters, filesParameters);
					  if(!saveResult) {
						  result[nFile][nChn][nSubChn]=null; // erase results to save memory
						  runtime.gc();
						  if (DEBUG_LEVEL>1) System.out.println("--- Free memory="+runtime.freeMemory()+" (of "+runtime.totalMemory()+")");
					  }
					  continue;
				  }
      SDFA_INSTANCE.showImageStack(stack, "processColorsWeights-RBG");
				  //Processing colors - changing stack sequence to r-g-b (was r-b-g)
				  if (!fixSliceSequence (stack)){
					  if(!saveResult) {
						  result[nFile][nChn][nSubChn]=null; // erase results to save memory
						  runtime.gc();
						  if (DEBUG_LEVEL>1) System.out.println("--- Free memory="+runtime.freeMemory()+" (of "+runtime.totalMemory()+")");
					  }
					  continue;
				  }
				  
// Debugging				  
//      SDFA_INSTANCE.showImageStackThree(stack, "processColorsWeights");
      SDFA_INSTANCE.showImageStack(stack, "processColorsWeights-RGB");
				  
				  processColorsWeights(stack,
						  255.0/PSF_SUBPIXEL_SHOULD_BE_4/PSF_SUBPIXEL_SHOULD_BE_4,
						  colorProcParameters,
						  colorCalibParameters,
						  nChn,
						  nSubChn
				  );
				  if (DEBUG_LEVEL>1) System.out.println("Processed colors to YPbPr, total number of slices="+stack.getSize());
// Show/save color denoise mask				  
				  if ((processParameters.saveChromaDenoiseMask || processParameters.showChromaDenoiseMask) && (DENOISE_MASK_CHROMA!=null)) {
					  ImagePlus chromaDenoiseMask=SDFA_INSTANCE.makeArrays (DENOISE_MASK_CHROMA,
							  DENOISE_MASK_CHROMA_WIDTH,
							  DENOISE_MASK_CHROMA.length/DENOISE_MASK_CHROMA_WIDTH,
							  title+"-MASK_CHROMA");
					  if (processParameters.jpeg) {
//crop Mask to original image size
						  if (processParameters.crop){
							  chromaDenoiseMask=cropImage32(chromaDenoiseMask,splitParameters);
						  }
//rotate the result			  
						  if (processParameters.rotate){
							  chromaDenoiseMask=rotateImage32CW(chromaDenoiseMask);
						  }
//scale the result
						  if (processParameters.JPEG_scale!=1.0){
							  ImageProcessor ip=chromaDenoiseMask.getProcessor();
							  ip.setInterpolationMethod(ImageProcessor.BICUBIC);
							  ip=ip.resize((int)(ip.getWidth()*processParameters.JPEG_scale),(int) (ip.getHeight()*processParameters.JPEG_scale));
							  chromaDenoiseMask= new ImagePlus(chromaDenoiseMask.getTitle(),ip);
							  chromaDenoiseMask.updateAndDraw();
						  }
						  if (processParameters.showChromaDenoiseMask) chromaDenoiseMask.show(); 
//public ImagePlus Image32toGreyRGB24(ImagePlus imp);
						  if (processParameters.saveChromaDenoiseMask) {
							  ImagePlus chromaDenoiseMaskRGB24=Image32toGreyRGB24(chromaDenoiseMask);
						    saveAndShow(chromaDenoiseMaskRGB24,
							  filesParameters,
							  processParameters.saveChromaDenoiseMask,
							  false, //processParameters.showChromaDenoiseMask,
							  processParameters.JPEG_quality);
						    chromaDenoiseMaskRGB24=null;
						  }
					  } else {
					    saveAndShow(chromaDenoiseMask,
							  filesParameters,
							  processParameters.saveChromaDenoiseMask,
							  processParameters.showChromaDenoiseMask
					    );
					  }
				  }
				  
				  if (processParameters.toRGB) {
					  YPrPbToRGB(stack,
							  colorProcParameters.kr,        // 0.299;
							  colorProcParameters.kb,        // 0.114;
							  colorProcParameters.useFirstY?9:8,        //  int sliceY,
									  6, // int slicePr,
									  7// int slicePb
					  );
					  title=titleFull; // including "-DECONV" or "-COMBO"
					  titleFull=title+"-RGB-float";
					  //Trim stack to just first 3 slices
					  while (stack.getSize()>3) stack.deleteLastSlice();
					  if (DEBUG_LEVEL>1) System.out.println("Trimming color stack");
				  } else {
					  title=titleFull; // including "-DECONV" or "-COMBO"
					  titleFull=title+"-YPrPb"; // including "-DECONV" or "-COMBO"
					  if (DEBUG_LEVEL>1) System.out.println("Using full stack, including YPbPr");
				  }
				  result[nFile][nChn][nSubChn]= new ImagePlus(titleFull, stack);    			  
				  // Crop image to match original one (scaled to oversampling)
				  if (processParameters.crop){
					  stack=cropStack32(stack,splitParameters);
				  }
				  // rotate the result			  
				  if (processParameters.rotate){
					  stack=rotateStack32CW(stack);
				  }
				  if (!processParameters.toRGB && !processParameters.jpeg){
					  saveAndShow(result[nFile][nChn][nSubChn], processParameters, filesParameters);
					  if(!saveResult) {
						  result[nFile][nChn][nSubChn]=null; // erase results to save memory
						  runtime.gc();
						  if (DEBUG_LEVEL>1) System.out.println("--- Free memory="+runtime.freeMemory()+" (of "+runtime.totalMemory()+")");
					  }
					  continue;
				  } else { // that's not the end result, save if required
					  saveAndShow(result[nFile][nChn][nSubChn], filesParameters, processParameters.save32, false,processParameters.JPEG_quality); // save, no show
				  }
				  // convert to RGB48 (16 bits per color component)			  
				  stack=convertRGB32toRGB16Stack(
						  stack,
						  rgbParameters); 
				  titleFull=title+"-RGB48";
				  result[nFile][nChn][nSubChn]= new ImagePlus(titleFull, stack);    			  
				  compositeImage=convertToComposite(result[nFile][nChn][nSubChn]);
				  if (!processParameters.jpeg){ // RGB48 was the end result
					  saveAndShow(compositeImage, processParameters, filesParameters);
					  if(!saveResult) {
						  result[nFile][nChn][nSubChn]=null; // erase results to save memory
						  runtime.gc();
						  if (DEBUG_LEVEL>1) System.out.println("--- Free memory="+runtime.freeMemory()+" (of "+runtime.totalMemory()+")");
					  }
					  continue;
				  } else { // that's not the end result, save if required
					  saveAndShow(compositeImage, filesParameters, processParameters.save16, false); // save, no show
				  }
				  ImagePlus imp_RGB24=convertRGB48toRGB24(
						  stack,
						  title+"-RGB24",
						  0, 65536, // r range 0->0, 65536->256
						  0, 65536, // g range
			    		  0, 65536, // b range
			    		  0, 65536);// alpha range
				  if (processParameters.JPEG_scale!=1.0){
					  ImageProcessor ip=imp_RGB24.getProcessor();
					  ip.setInterpolationMethod(ImageProcessor.BICUBIC);
					  ip=ip.resize((int)(ip.getWidth()*processParameters.JPEG_scale),(int) (ip.getHeight()*processParameters.JPEG_scale));
					  imp_RGB24= new ImagePlus(imp_RGB24.getTitle(),ip);
					  imp_RGB24.updateAndDraw();

				  }
				  saveAndShow(imp_RGB24, processParameters, filesParameters);
				  if(!saveResult) {
					  result[nFile][nChn][nSubChn]=null; // erase results to save memory
					  runtime.gc();
					  if (DEBUG_LEVEL>1) System.out.println("---- Free memory="+runtime.freeMemory()+" (of "+runtime.totalMemory()+")");
				  }
			  }  // next nSubChn
			  if(!saveResult) {
				  result[nFile][nChn]=null;
				  runtime.gc();
				  if (DEBUG_LEVEL>0) System.out.println("At "+IJ.d2s(0.000000001*(System.nanoTime()-startTime),3)+
						  " seconds  Free memory="+IJ.d2s(runtime.freeMemory()/(1024.0*1024.0*1024.0),3)+" GB (of "+
						  IJ.d2s(runtime.totalMemory()/(1024.0*1024.0*1024.0),3)+" GB), used "+
						  IJ.d2s((runtime.totalMemory()-runtime.freeMemory())/(1024.0*1024.0*1024.0),3)+" GB");
			  }

		  } // next nChn, next nFile
	  System.out.println("Processing done in "+IJ.d2s(0.000000001*(System.nanoTime()-startTime),3)+ " seconds");
	  return result;
  }
  /* ======================================================================== */
  private void saveAndShow(
		  ImagePlus             imp,
		  ProcessParameters     processParameters,
		  FilesParameters       filesParameters){
	  saveAndShow( imp,    processParameters, filesParameters , true, true);
  }

  private void saveAndShow(
		  ImagePlus             imp,
		  ProcessParameters     processParameters,
		  FilesParameters       filesParameters,
		  boolean               enableSave,
		  boolean               enableShow){
	  saveAndShow(
			  imp,
			  filesParameters,
			  processParameters.save && enableSave,
			  processParameters.show && enableShow,
			  processParameters.JPEG_quality);
  }

  private void saveAndShow(
		  ImagePlus             imp,
		  FilesParameters       filesParameters,
		  boolean               save,
		  boolean               show){
	  saveAndShow(imp, filesParameters,  save,  show, -1);
  } 
  
  private void saveAndShow(
		  ImagePlus             imp,
		  FilesParameters       filesParameters,
		  boolean               save,
		  boolean               show,
		  int                   jpegQuality){
	  String path;
	  if (save) {
		//Ask for the destination directory if it is undefined
		  if ((filesParameters.resultsDirectory==null) || (filesParameters.resultsDirectory.length()==0)) {
			  path= selectResultsDirectory(filesParameters.resultsDirectory);
			  if (path!=null) filesParameters.resultsDirectory=path;
		  }
		  File destDir= new File (filesParameters.resultsDirectory);
		  if (!destDir.exists()){
			  if (!destDir.mkdirs()) {
				  IJ.showMessage("Error","Failed to create results directory "+filesParameters.resultsDirectory);
				  save=false;
			  }
		  }
	  }
	  if (save) {
		  path=filesParameters.resultsDirectory+Prefs.getFileSeparator()+imp.getTitle();
		  if (((imp.getStackSize()==1)) && ((imp.getFileInfo().fileType== FileInfo.RGB) || (jpegQuality>0))) {
			  if (DEBUG_LEVEL>0) System.out.println("Saving result to "+path+".jpeg");
			  FileSaver fs=new FileSaver(imp);
			  if (jpegQuality>0) FileSaver.setJpegQuality(jpegQuality);
			  fs.saveAsJpeg(path+".jpeg");
		  }
		  else {
			  if (DEBUG_LEVEL>0) System.out.println("Saving result to "+path+".tiff");
			  FileSaver fs=new FileSaver(imp);
			  if (imp.getStackSize()>1)  fs.saveAsTiffStack(path+".tiff");
			  else fs.saveAsTiff(path+".tiff");
		  }
	  }
	  if (show) {
		  imp.getProcessor().resetMinAndMax(); // probably not needed
		  imp.show();
	  }
  }
 
  private void saveAndShow(
		  CompositeImage        compositeImage,
		  ProcessParameters     processParameters,
		  FilesParameters       filesParameters){
	  saveAndShow(compositeImage,    processParameters, filesParameters , true, true);
  }

  private void saveAndShow(
		  CompositeImage        compositeImage,
		  ProcessParameters     processParameters,
		  FilesParameters       filesParameters,
		  boolean               enableSave,
		  boolean               enableShow){

	  saveAndShow(
			  compositeImage,
			  filesParameters,
			  processParameters.save && enableSave,
			  processParameters.show && enableShow);
  }

  private void saveAndShow(
		  CompositeImage        compositeImage,
		  FilesParameters       filesParameters,
		  boolean               save,
		  boolean               show){
	  String path;
	  if (save) {
			//Ask for the destination directory if it is undefined
			  if ((filesParameters.resultsDirectory==null) || (filesParameters.resultsDirectory.length()==0)) {
				  path= selectResultsDirectory(filesParameters.resultsDirectory);
				  if (path!=null) filesParameters.resultsDirectory=path;
			  }
			  File destDir= new File (filesParameters.resultsDirectory);
			  if (!destDir.exists()){
				  if (!destDir.mkdirs()) {
					  IJ.showMessage("Error","Failed to create results directory "+filesParameters.resultsDirectory);
					  save=false;
				  }
			  }
		  }

	  if (save) {
		  path=filesParameters.resultsDirectory+Prefs.getFileSeparator()+compositeImage.getTitle();
		  if (DEBUG_LEVEL>0) System.out.println("Saving result to "+path+".tiff");
		  FileSaver fs=new FileSaver(compositeImage);
		  if (compositeImage.getStackSize()>1)  fs.saveAsTiffStack(path+".tiff");
		  else fs.saveAsTiff(path+".tiff");
//		  IJ.saveAs(compositeImage,"tif",path);
	  }

	  if (show) {
		  compositeImage.show();
	  }
  }
  
  
  
  /* ======================================================================== */
	/* Create a Thread[] array as large as the number of processors available.
	 * From Stephan Preibisch's Multithreading.java class. See:
	 * http://repo.or.cz/w/trakem2.git?a=blob;f=mpi/fruitfly/general/MultiThreading.java;hb=HEAD
	 */
	private Thread[] newThreadArray(int maxCPUs) {
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

/* ======================================================================== */
/* ======================================================================== */
/* ======================================================================== */
/* ======================================================================== */
/* ======================================================================== */
/* ======================================================================== */
/* ======================================================================== */
/* ======================================================================== */
/* ======================================================================== */
  public double[] initWindowFunction(int size) {
    double [] windowFunction =new double [size*size];
    double [] windowFunction_line=new double [size];
    double a,k;
    int i,j;
    if (GAUSS_WIDTH<=0) {
      for (i=0; i<size; i++) windowFunction_line[i]= (0.54-0.46*Math.cos((i*2.0*Math.PI)/size));
    } else {
      k=2.0/(size*GAUSS_WIDTH);
      for (i=0; i<size; i++) {
         a=(i-size/2)*k;
         windowFunction_line[i]= Math.exp( - a*a);
      }
    }
    for (i=0; i<size; i++) for (j=0; j<size; j++){
       windowFunction[size*i+j]=windowFunction_line[i]*windowFunction_line[j];
    }
    return windowFunction;
  }

/* ======================================================================== */
  /* convolve image stack with the kernel stack using FHT. kernels ahould be (size/2)*(size/2) - currently 64x64, then image will be split into same 
      (size/2)*(size/2) overlapping by step=size/4 segments. Both are zero-padded to size x size, so after convolution the result will not roll over, and
      processed 128x128 result arrays are accumulated in the output stack.
      The input image should be properly extended by size/4 in each direction (and so the kernel arrays should match it) - that would minimize border effects.*/
 
  /* ======================================================================== */
  public ImageStack convolveStackWithKernelStack (
		  final ImageStack  imageStack,  // stack with 3 colors/slices with the image
		  final ImageStack kernelStack, // stack with 3 colors/slices convolution kernels
		  final int               size, // 128 - fft size, kernel size should be size/2 
		  final int          threadsMax,  // maximal number of threads to launch                         
		  final boolean    updateStatus) // update status info
  {
	  if ((imageStack==null) || (kernelStack==null)) return null;
	  final int imgWidth=imageStack.getWidth();
	  final int imgHeight=imageStack.getHeight();
	  final int length=imgWidth*imgHeight;
	  final int step=size/4;
	  final int kernelSize=size/2;
	  final int tilesX=imgWidth/step-1; // horizontal number of overlapping tiles in the source image (should be expanded from the registerd one by "step" in each direction)
	  final int tilesY=imgHeight/step-1; // vertical number of overlapping tiles in the source image (should be expanded from the registerd one by "step" in each direction)
	  final int kernelWidth=kernelStack.getWidth();
	  final int kernelNumHor=kernelWidth/(size/2);

	  final int nChn=imageStack.getSize();
	  final float [][] outPixels=new float[nChn][length]; // GLOBAL same as input
	  //	   float [][] outPixels=new float[nChn][length]; // same as input
	  int i,j;
	  for (i=0;i<nChn;i++) for (j=0;j<length;j++) outPixels[i][j]=0.0f;
	  final double [] slidingWindow=getSlidingMask(kernelSize); // 64x64
	  final Thread[] threads = newThreadArray(threadsMax);
	  final AtomicInteger ai = new AtomicInteger(0);
	  final int numberOfKernels=     tilesY*tilesX*nChn;
	  final int numberOfKernelsInChn=tilesY*tilesX;
	  final long startTime = System.nanoTime();
	  for (int ithread = 0; ithread < threads.length; ithread++) {
		  threads[ithread] = new Thread() {
			  public void run() {
				  float [] pixels=null;       // will be initialized at first use
				  float [] kernelPixels=null; // will be initialized at first use
				  double [] kernel=       new double[kernelSize*kernelSize];
				  double [] inTile=       new double[kernelSize*kernelSize];
				  double [] outTile=      new double[size * size];
				  double [] doubleKernel= new double[size * size];
				  int chn,tileY,tileX;
				  int chn0=-1;
//				  double debug_sum;
//				  int i;
				  DoubleFHT fht_instance =new DoubleFHT(); // provide DoubleFHT instance to save on initializations (or null)
				  for (int nTile = ai.getAndIncrement(); nTile < numberOfKernels; nTile = ai.getAndIncrement()) {
					  chn=nTile/numberOfKernelsInChn;
					  tileY =(nTile % numberOfKernelsInChn)/tilesX;
					  tileX = nTile % tilesX;
					  if (tileX==0) {
						  if (updateStatus) IJ.showStatus("Convolving image with kernels, channel "+(chn+1)+" of "+nChn+", row "+(tileY+1)+" of "+tilesY);
						  if (MASTER_DEBUG_LEVEL>2) System.out.println("Processing kernels, channel "+(chn+1)+" of "+nChn+", row "+(tileY+1)+" of "+tilesY+" : "+IJ.d2s(0.000000001*(System.nanoTime()-startTime),3));
					  }
					  
					  if (chn!=chn0) {
						  pixels=      (float[]) imageStack.getPixels(chn+1);
						  kernelPixels=(float[]) kernelStack.getPixels(chn+1);
						  chn0=chn;
					  }
					  /* Read source image tile */
					  extractSquareTile( pixels, // source pixel array,
							  inTile, // will be filled, should have correct size before call
							  slidingWindow, // window (same size as the kernel)
							  imgWidth, // width of pixels array
							  tileX*step, // left corner X
							  tileY*step); // top corner Y
					  /* zero pad twice the original size*/
					  outTile=extendFFTInputTo (inTile, size);
					  /* FHT transform of the source image data*/
					  fht_instance.swapQuadrants(outTile);
					  fht_instance.transform(    outTile);
					  /* read convolution kernel */
					  extractOneKernel(kernelPixels, //  array of combined square kernels, each 
							  kernel, // will be filled, should have correct size before call
							  kernelNumHor, // number of kernels in a row
							  //tileX*kernelSize, // horizontal number of kernel to extract
							  //tileY*kernelSize); // vertical number of kernel to extract
							  tileX, // horizontal number of kernel to extract
							  tileY); // vertical number of kernel to extract
					  /* zero pad twice the original size*/
					  doubleKernel=extendFFTInputTo (kernel, size);
//					  debug_sum=0;
//					  for (i=0;i<doubleKernel.length;i++) debug_sum+=doubleKernel[i];
//					  if (MASTER_DEBUG_LEVEL>1) System.out.println("kernel sum="+debug_sum);
					  
					  //if ((tileY==tilesY/2) && (tileX==tilesX/2))  SDFA_INSTANCE.showArrays(doubleKernel,size,size, "doubleKernel-"+chn);
					  /* FHT transform of the kernel */
					  fht_instance.swapQuadrants(doubleKernel);
					  fht_instance.transform(    doubleKernel);
					  /* multiply in frequency domain */
					  outTile=     fht_instance.multiply(outTile, doubleKernel, false);
					  /* FHT inverse transform of the product - back to space domain */
					  fht_instance.inverseTransform(outTile);
					  fht_instance.swapQuadrants(outTile);
					  /* accumulate result */
					  //if ((tileY==tilesY/2) && (tileX==tilesX/2))  SDFA_INSTANCE.showArrays(outTile,size,size, "out-"+chn);
					  /*This is synchronized method. It is possible to make threads to write to non-overlapping regions of the outPixels, but as the accumulation
					   * takes just small fraction of severtal FHTs, it should be OK - reasonable number of threads will spread and not "stay in line"
					   */
					  accumulateSquareTile(outPixels[chn], //  float pixels array to accumulate tile
							  outTile, // data to accumulate to the pixels array
							  imgWidth, // width of pixels array
							  (tileX-1)*step, // left corner X
							  (tileY-1)*step); // top corner Y
				  }
			  }
		  };
	  }		      
	  startAndJoin(threads);
	  if (DEBUG_LEVEL > 1) System.out.println("Threads done at "+IJ.d2s(0.000000001*(System.nanoTime()-startTime),3));

	  /* prepare result stack to return */
	  ImageStack outStack=new ImageStack(imgWidth,imgHeight);
	  for (i=0;i<nChn;i++) {
		  outStack.addSlice(imageStack.getSliceLabel(i+1), outPixels[i]);
	  }
	  return outStack;
  }
 /* Extract noise mask (proportional to noise gain of the kernels), the denoise mask should be divided by this
  *  
  */
  public double [][] extractNoiseMask(
		     ImagePlus imp,// contains 3-slice stack (r,b,g)
		     double k0,    // coefficient for slice 0 (r)
		     double k1,    // coefficient for slice 1 (b)
		     double k2,    // coefficient for slice 2 (g)
		     int decim,     // decimate result (not yet supported)
		     double gainPower
		  ){
	  if (imp==null) return null;
	  if (gainPower==0.0) return null; // do not use noise gain correction
	  ImageStack stack=imp.getStack();
	  int width=stack.getWidth();
	  int height=stack.getHeight();
	  double [][]mask=new double[height*decim][width*decim];
	  float [][] pixels= new float [3][];
	  pixels[0]=(float[]) stack.getPixels(1);
	  pixels[1]=(float[]) stack.getPixels(2);
	  pixels[2]=(float[]) stack.getPixels(3);
	  int i,j,x,y,indx;
	  for (y=0;y<height;y++) for (x=0;x<width;x++) {
		  indx=x+y*width;
		  for (i=0;i<decim;i++) for (j=0;j<decim;j++) {
		     mask[y*decim+i][x*decim+j]=Math.pow(k0*pixels[0][indx]+k1*pixels[1][indx]+k2*pixels[2][indx],gainPower);
		  }
	  }
	  return mask;
  }
/* ======================================================================== */ 
/* Calculate deconvolution kernel (or difference of the two) noise gain
 *  to be used when calculating mask that selects between deconvolved with
 *  different kernels
 */
  public ImageStack calculateKernelsNoiseGains (
		  final ImageStack kernelStack1, // first stack with 3 colors/slices convolution kernels
		  final ImageStack kernelStack2, // second stack with 3 colors/slices convolution kernels (or null)
		  final int               size, // 128 - fft size, kernel size should be size/2
		  final double       blurSigma,
		  final int          threadsMax,  // maximal number of threads to launch                         
		  final boolean    updateStatus) // update status info
  {
	  if (kernelStack1==null) return null;
	  final boolean useDiff= (kernelStack2 != null);
	  final int kernelSize=size/2;
	  final int kernelWidth=kernelStack1.getWidth();
	  final int kernelNumHor=kernelWidth/(size/2);
	  final int kernelNumVert=kernelStack1.getHeight()/(size/2);
	  final int length=kernelNumHor*kernelNumVert;
	  final int nChn=kernelStack1.getSize();
	  final float [][] outPixles=new float[nChn][length]; // GLOBAL same as input
	  int i,j;
	  for (i=0;i<nChn;i++) for (j=0;j<length;j++) outPixles[i][j]=0.0f;
	  final Thread[] threads = newThreadArray(threadsMax);
	  final AtomicInteger ai = new AtomicInteger(0);
	  final int numberOfKernels=     kernelNumHor*kernelNumVert*nChn;
	  final int numberOfKernelsInChn=kernelNumHor*kernelNumVert;
	  final long startTime = System.nanoTime();
	  for (int ithread = 0; ithread < threads.length; ithread++) {
		  threads[ithread] = new Thread() {
			  public void run() {
				  DoubleGaussianBlur gb=null;
				  if (blurSigma>0)	 gb=new DoubleGaussianBlur();
				  float [] kernelPixels1= null; // will be initialized at first use
				  float [] kernelPixels2= null; // will be initialized at first use
				  double [] kernel1=      new double[kernelSize*kernelSize];
				  double [] kernel2=      new double[kernelSize*kernelSize];
				  int chn,tileY,tileX;
				  int chn0=-1;
				  int i;
				  double sum;
				  for (int nTile = ai.getAndIncrement(); nTile < numberOfKernels; nTile = ai.getAndIncrement()) {
					  chn=nTile/numberOfKernelsInChn;
					  tileY =(nTile % numberOfKernelsInChn)/kernelNumHor;
					  tileX = nTile % kernelNumHor;
					  if (tileX==0) {
						  if (updateStatus) IJ.showStatus("Processing kernels, channel "+(chn+1)+" of "+nChn+", row "+(tileY+1)+" of "+kernelNumVert);
						  if (MASTER_DEBUG_LEVEL>2) System.out.println("Processing kernels, channel "+(chn+1)+" of "+nChn+", row "+(tileY+1)+" of "+kernelNumVert+" : "+IJ.d2s(0.000000001*(System.nanoTime()-startTime),3));
					  }
					  
					  if (chn!=chn0) {
						  kernelPixels1=(float[]) kernelStack1.getPixels(chn+1);
						  if (useDiff) kernelPixels2=(float[]) kernelStack2.getPixels(chn+1);
						  chn0=chn;
					  }
					  /* read convolution kernel */
					  extractOneKernel(kernelPixels1, //  array of combined square kernels, each 
							  kernel1, // will be filled, should have correct size before call
							  kernelNumHor, // number of kernels in a row
							  tileX, // horizontal number of kernel to extract
							  tileY); // vertical number of kernel to extract
					  /* optionally read the second convolution kernel */
					  if (useDiff) {extractOneKernel(kernelPixels2, //  array of combined square kernels, each 
							  kernel2, // will be filled, should have correct size before call
							  kernelNumHor, // number of kernels in a row
							  tileX, // horizontal number of kernel to extract
							  tileY); // vertical number of kernel to extract
					     for (i=0; i<kernel1.length;i++) kernel1[i]-=kernel2[i];
					  }
					  if (blurSigma>0) gb.blurDouble(kernel1, kernelSize, kernelSize, blurSigma, blurSigma, 0.01);
					  /* Calculate sum of squared kernel1  elements */
					  sum=0.0;
					  for (i=0; i<kernel1.length;i++) sum+=kernel1[i]*kernel1[i];
					  outPixles[chn][tileY*kernelNumHor+tileX]= (float) (Math.sqrt(sum));
//					  System.out.println("Processing kernels, channel "+(chn+1)+" of "+nChn+", row "+(tileY+1)+" of "+kernelNumVert+" sum="+sum);
				  }
			  }
		  };
	  }		      
	  startAndJoin(threads);
	  if (DEBUG_LEVEL > 1) System.out.println("Threads done at "+IJ.d2s(0.000000001*(System.nanoTime()-startTime),3));
	  /* prepare result stack to return */
	  ImageStack outStack=new ImageStack(kernelNumHor,kernelNumVert);
	  for (i=0;i<nChn;i++) {
		  outStack.addSlice(kernelStack1.getSliceLabel(i+1), outPixles[i]);
	  }
	  return outStack;
  }
 
  
  /* ======================================================================== */ 
  // uses global OUT_PIXELS to accumulate results
  public ImageStack aliasScissorsStack (
		  final ImageStack                        imageStack,  // stack with 3 colors/slices with the image
		  final DebayerParameters          debayerParameters, // 64 - fft size
		  final boolean                generateDebayerEnergy,
		  final int                               threadsMax, // maximal step in pixels on the maxRadius for 1 angular step (i.e. 0.5)
		  final boolean                         updateStatus) // update status info

  {
	  final int wasDebugLevel=DEBUG_LEVEL;
	  // this debug of a specific tile will work in a single-threading only, because it uses global DEBUG_LEVEL
	  final int xTileDebug, yTileDebug;
      if (!debayerParameters.debug || (debayerParameters.xDebug<0) || (debayerParameters.xDebug<0)) {
    	  xTileDebug=-1;
    	  yTileDebug=-1;
      } else {
    	  xTileDebug=(debayerParameters.xDebug>=(debayerParameters.size/4))?((debayerParameters.xDebug-debayerParameters.size/4)/(debayerParameters.size/2)):0;
    	  yTileDebug=(debayerParameters.yDebug>=(debayerParameters.size/4))?((debayerParameters.yDebug-debayerParameters.size/4)/(debayerParameters.size/2)):0;
      }


	  if (imageStack==null) return null;
	  final int imgWidth=imageStack.getWidth();
	  final int imgHeight=imageStack.getHeight();
	  final int length=imgWidth*imgHeight;
	  final int step=debayerParameters.size/2;
	  final int tilesX=imgWidth/step-1; // horizontal number of overlapping tiles in the source image (should be expanded from the registerd one by "step" in each direction)
	  final int tilesY=imgHeight/step-1; // vertical number of overlapping tiles in the source image (should be expanded from the registerd one by "step" in each direction)
	  final int nChn=imageStack.getSize();
	  int i,chn; //tileX,tileY;
	  /* find number of the green channel - should be called "green", if none - use last */
	  i=nChn-1;
	  for (chn=0;chn<nChn;chn++) if (imageStack.getSliceLabel(chn+1).equals("green")){
		  i=chn;
		  break;
	  }
	  final int greenChn=i;
	  final float [][] outPixles=new float[nChn][length]; // same as input
	  DEBAYER_ENERGY=null;
	  if (generateDebayerEnergy) {
		  DEBAYER_ENERGY=new double[tilesY*tilesX];
	  }

	  for (chn=0;chn<nChn;chn++) for (i=0;i<length;i++) outPixles[chn][i]=0.0f;
	  final double [] slidingWindow= getSlidingMask(debayerParameters.size); // 64x64

//	  outPixles=new float[nChn][length]; // GLOBAL same as input
	  final Thread[] threads = newThreadArray(threadsMax);
	  final AtomicInteger ai = new AtomicInteger(0);
	  final int numberOfKernels=tilesY*tilesX;
	  final long startTime = System.nanoTime();
	  for (int ithread = 0; ithread < threads.length; ithread++) {
		  threads[ithread] = new Thread() {
			  public void run() {
				  double [][] tile=        new double[nChn][debayerParameters.size * debayerParameters.size ];
				  double [][] both_masks;
				  float [][] pixels=       new float[nChn][];
				  int chn,tileY,tileX,i;
				  for (chn=0;chn<nChn;chn++) pixels[chn]= (float[]) imageStack.getPixels(chn+1);
				  DoubleFHT       fht_instance =   new DoubleFHT(); // provide DoubleFHT instance to save on initializations (or null)
				  ShowDoubleFloatArrays SDFA_instance=null; // just for debugging?

				  DeBayerScissors debayer_instance=new DeBayerScissors( debayerParameters.size, // size of the square array, centar is at size/2, size/2, only top half+line will be used
						  debayerParameters.polarStep, // maximal step in pixels on the maxRadius for 1 angular step (i.e. 0.5)
						  debayerParameters.debayerRelativeWidthGreen, // result green mask mpy by scaled default (diamond)
						  debayerParameters.debayerRelativeWidthRedblue, // result red/blue mask mpy by scaled default (square)
						  debayerParameters.debayerRelativeWidthRedblueMain, // green mask when applied to red/blue, main (center)
						  debayerParameters.debayerRelativeWidthRedblueClones);// green mask when applied to red/blue, clones 
				  
				  for (int nTile = ai.getAndIncrement(); nTile < numberOfKernels; nTile = ai.getAndIncrement()) {
					  tileY = nTile /tilesX;
					  tileX = nTile % tilesX;
					  if (tileX==0) {
						  if (updateStatus) IJ.showStatus("(3)Reducing sampling aliases, row "+(tileY+1)+" of "+tilesY);
						  if (MASTER_DEBUG_LEVEL>2) System.out.println("(3)Reducing sampling aliases, row "+(tileY+1)+" of "+tilesY+" : "+IJ.d2s(0.000000001*(System.nanoTime()-startTime),3));
					  }
				  
					  if ((tileY==yTileDebug) && (tileX==xTileDebug)) DEBUG_LEVEL=4;
					  else DEBUG_LEVEL=wasDebugLevel;
					  for (chn=0;chn<nChn;chn++){
						  extractSquareTile( pixels[chn], // source pixel array,
								  tile[chn], // will be filled, should have correct size before call
								  slidingWindow, // window (same size as the kernel)
								  imgWidth, // width of pixels array
								  tileX*step, // left corner X
								  tileY*step); // top corner Y
					  }

					  /* Scale green channel x0.5 as there are twice more pixels there as in red or blue. Or move it somewhere else and multiply to original range ? */
					  for (i=0;i<tile[greenChn].length;i++) tile[greenChn][i]*=0.5;
					  if ((tileY==yTileDebug) && (tileX==xTileDebug)) {
						  if (SDFA_instance==null) SDFA_instance=      new ShowDoubleFloatArrays();
						  SDFA_instance.showArrays (tile.clone(),debayerParameters.size,debayerParameters.size, "x"+(tileX*step)+"_y"+(tileY*step));
					  }
					  for (chn=0;chn<nChn;chn++){
						  fht_instance.swapQuadrants(tile[chn]);
						  fht_instance.transform(tile[chn]);
					  }
					  if ((tileY==yTileDebug) && (tileX==xTileDebug) && (SDFA_instance!=null)) SDFA_instance.showArrays (tile.clone(),debayerParameters.size,debayerParameters.size, "tile-fht");
					  both_masks= debayer_instance.aliasScissors(tile[greenChn], // fht array for green, will be masked in-place
							  debayerParameters.debayerThreshold, // no high frequencies - use default uniform filter
							  debayerParameters.debayerGamma, // power function applied to the amplitudes before generating spectral masks
							  debayerParameters.debayerBonus, // scale far pixels as (1.0+bonus*r/rmax)
							  debayerParameters.mainToAlias,// relative main/alias amplitudes to enable lixels (i.e. 0.5 means that if alias is >0.5*main, the pixel will be masked out)
							  debayerParameters.debayerMaskBlur, // for both masks  sigma for gaussian blur of the binary masks (<0 -do not use "scissors")
							  debayerParameters.debayerUseScissors, // use "scissors", if false - just apply "diamond" ands "square" with DEBAYER_PARAMETERS.debayerRelativeWidthGreen and DEBAYER_PARAMETERS.debayerRelativeWidthRedblue
							  ((tileY==yTileDebug) && (tileX==xTileDebug))?4:1);
					  //                                               1); // internal debug level ((DEBUG_LEVEL>2) && (yTile==yTile0) && (xTile==xTile0))?3:1;
					  if ((tileY==yTileDebug) && (tileX==xTileDebug) && (SDFA_instance!=null)) {
						  SDFA_instance.showArrays (tile.clone(),debayerParameters.size,debayerParameters.size, "A00");
						  SDFA_instance.showArrays (both_masks.clone(),debayerParameters.size,debayerParameters.size, "masks");
					  }
					  if (DEBAYER_ENERGY!=null) {
						  DEBAYER_ENERGY[tileY*tilesX+tileX]=debayer_instance.getMidEnergy();
					  }
					  for (chn=0;chn<nChn;chn++) {
						  tile[chn]=fht_instance.multiply(tile[chn],both_masks[(chn==greenChn)?0:1],false);
						  fht_instance.inverseTransform(tile[chn]);
						  fht_instance.swapQuadrants(tile[chn]);
						  /* accumulate result */
						  /*This is synchronized method. It is possible to make threads to write to non-overlapping regions of the outPixles, but as the accumulation
						   * takes just small fraction of severtal FHTs, it should be OK - reasonable number of threads will spread and not "stay in line"
						   */

						  accumulateSquareTile(outPixles[chn], //  float pixels array to accumulate tile
								  tile[chn], // data to accumulate to the pixels array
								  imgWidth, // width of pixels array
								  tileX*step, // left corner X
								  tileY*step); // top corner Y
					  }
					  if ((tileY==yTileDebug) && (tileX==xTileDebug) && (SDFA_instance!=null)) SDFA_instance.showArrays (tile.clone(),debayerParameters.size,debayerParameters.size, "B00");
					  
				  }
			  }
		  };
	  }		      
	  startAndJoin(threads);
	  DEBUG_LEVEL=wasDebugLevel;
	  /* prepare result stack to return */
	  ImageStack outStack=new ImageStack(imgWidth,imgHeight);
	  for (chn=0;chn<nChn;chn++) {
		  outStack.addSlice(imageStack.getSliceLabel(chn+1), outPixles[chn]);
	  }
	  DEBAYER_ENERGY_WIDTH=	 (DEBAYER_ENERGY!=null)?tilesX:0; // for the image to be displayed externally 
//	  if (debayerParameters.showEnergy) {
//		  SDFA_INSTANCE.showArrays (DEBAYER_ENERGY,tilesX,tilesY, "Debayer-Energy");
//	  }

	  return outStack;
  }


  /* ======================================================================== */
  /* Filters mask that selects between hi-res/high-noise deconvolved image and lo-res/lo-noise image convolved with Gaussian
   * by rejecting frequencies that correspond to multiples of JPEG blocks (here with the current settings it is 32 pixels - twice 16x16 macroblock
   */
  // uses global MASK_LOHIRES to accumulate results
  public double [] createFilterForBlockArtifacts(
		  final int             size, // size of square FHT
		  final int     rejectPeriod, // period (pixels) of the block artifacts to reject (32)
		  final double   rejectSigma, // sigma of the rejection spots ( 0.0 - just zero a single point)
		  final double   lopassSigma) // sigma of the low pass filter (frequency units, 0.0 - do not filter)
  {
	  double [] maskFHT=new double[size*size];
	  int freqPeriod=size/rejectPeriod;
	  int i,j;
	  for (i=0;i<size*size;i++) maskFHT[i]=1.0; // Initialize mask
	  double []rejGauss;
	  double k;
	  int pointX,pointY,x,y;
	  if (rejectSigma>=0.0) {
		  int rejSize= (int) (rejectSigma*4)+1;
		  if (rejSize>(freqPeriod/2)) rejSize=(freqPeriod/2);
		  rejGauss=new double[rejSize];
		  rejGauss[0]=1.0;
		  k=1.0/(2*rejectSigma*rejectSigma);
		  for (i=1;i<rejSize;i++) {
			  rejGauss[i]=Math.exp(-(k*i*i));
		  }
		  for (pointY=0;pointY<size;pointY+=freqPeriod) for (pointX=0;pointX<size;pointX+=freqPeriod) if ((pointY>0) || (pointX>0)) {
			  for (y=-rejSize+1;y<rejSize;y++) for (x=-rejSize+1;x<rejSize;x++) {
				  i=(pointY+y+size)%size;
				  j=(pointX+x+size)%size;
				  maskFHT[i*size+j]=1.0-rejGauss[Math.abs(y)]*rejGauss[Math.abs(x)];

			  }
		  }
	  }

	  if (lopassSigma>0) {
		  double [] maskTmp=maskFHT;
		  maskFHT=new double[size*size];
		  for (i=0;i<size*size;i++) maskFHT[i]=0.0; // Initialize mask

		  int lopassSize= (int) (lopassSigma*4)+1;
		  if (lopassSize>(size/2)) lopassSize=(size/2);
		  rejGauss=new double[lopassSize];
		  rejGauss[0]=1.0;
		  k=1.0/(2*lopassSigma*lopassSigma);
		  for (i=1;i<lopassSize;i++) {
			  rejGauss[i]=Math.exp(-(k*i*i));
		  }
		  for (y=-lopassSize+1;y<lopassSize;y++) for (x=-lopassSize+1;x<lopassSize;x++) {
			  i=(y+size)%size;
			  j=(x+size)%size;
			  maskFHT[i*size+j]=maskTmp[i*size+j]*rejGauss[Math.abs(y)]*rejGauss[Math.abs(x)];
		  }
	  }
	  return maskFHT;
  }
  
  private int extendDimension(
		  int dimension,
		  int step) {
	  return step*((int) Math.ceil(dimension/step) +2);
  }
  public double [] extendDoubleArrayForSlidingWindow(
		  double []     ipixels, // input pixel array
		  int           imgWidth,  // width of the image
		  int             step){ // size of sliding step (half of the sliding window size)

	  int imgHeight= ipixels.length/imgWidth;  // width of the image
	  int width=extendDimension(imgWidth,step);
	  int height=extendDimension(imgHeight,step);
	  double [] pixels=new double[width*height];
	  int i,j,k;
	  int index=0;
	  for (i=0;i<height;i++) {
		  if ((i<step) || (i>=(imgHeight+step))) {
			  k=(i+1)*width;
			  for (j=i*width;j<k;j++) pixels[j]=0.0;
		  } else {
			  k=i*width+step;
			  for (j=i*width;j<k;j++) pixels[j]=0.0;
			  k=(i+1)*width;
			  for (j=i*width+step+imgWidth;j<k;j++) pixels[j]=0.0;
			  k=i*width+step+imgWidth;
			  for (j=i*width+step;j<k;j++) pixels[j]=ipixels[index++];

		  }
	  }
	  return pixels;
  }

  public double [] reducedDoubleArrayAfterSlidingWindow(
		  double []     ipixels, // input pixel array
		  int           imgWidth,  // width of the image
		  int           imgHeight,
		  int             step){ // size of sliding step (half of the sliding window size)

	  int width=extendDimension(imgWidth,step);
	  double [] pixels=new double[imgWidth*imgHeight];
	  int i,j,base;
	  int index=0;
	  for (i=0;i< imgHeight;i++) {
		  base=width*(i+step)+step;
		  for (j=0;j<imgWidth;j++) pixels[index++]=ipixels[base++];
	  }
	  return pixels;
  }

  
  
  
  public double [] filterMaskFromBlockArtifacts(
//  public float [] filterMaskFromBlockArtifacts(
		  final double []     pixels, // input pixel array
//		  final float []     pixels, // input pixel array
		  final int         imgWidth, // width of the image
		  final int        imgHeight, // width of the image
		  final int             size, // size of sliding FHT
		  final double []     filter, // filter to multiply FHT (created once for the whole filter mask)
		  final int       threadsMax, // maximal step in pixels on the maxRadius for 1 angular step (i.e. 0.5)
		  final boolean updateStatus) // update status info
  {

	  if (MASTER_DEBUG_LEVEL>1) System.out.println("filterMaskFromBlockArtifacts, imgWidth="+imgWidth);
	  if (MASTER_DEBUG_LEVEL>1) System.out.println("filterMaskFromBlockArtifacts, imgHeight="+imgHeight);

	  if (pixels==null) return null;
	  final int length=imgWidth*imgHeight;
	  final int step=size/2;
	  final int tilesX=imgWidth/step-1; // horizontal number of overlapping tiles in the source image (should be expanded from the registered one by "step" in each direction)
	  final int tilesY=imgHeight/step-1; // vertical number of overlapping tiles in the source image (should be expanded from the registered one by "step" in each direction)
	  if (MASTER_DEBUG_LEVEL>1) System.out.println("filterMaskFromBlockArtifacts, tilesX="+tilesX);
	  if (MASTER_DEBUG_LEVEL>1) System.out.println("filterMaskFromBlockArtifacts, tilesY="+tilesY);
	  
	  int i; //tileX,tileY;

//	  for (i=0;i<length;i++) MASK_LOHIRES[i]=0.0;
//	  MASK_LOHIRES=new float[length];
	  MASK_LOHIRES=new double[length];
	  for (i=0;i<length;i++) MASK_LOHIRES[i]=0.0;
	  final double [] slidingWindow= getSlidingMask(size); // 256x256?
	  final Thread[] threads = newThreadArray(threadsMax);
	  final AtomicInteger ai = new AtomicInteger(0);
	  final int numberOfKernels=tilesY*tilesX;
	  final long startTime = System.nanoTime();
	  for (int ithread = 0; ithread < threads.length; ithread++) {
		  threads[ithread] = new Thread() {
			  public void run() {
				  double [] tile=        new double[size * size ];
				  int tileY,tileX;
				  DoubleFHT       fht_instance =   new DoubleFHT(); // provide DoubleFHT instance to save on initializations (or null)
//				  showDoubleFloatArrays SDFA_instance=null; // just for debugging?
				  for (int nTile = ai.getAndIncrement(); nTile < numberOfKernels; nTile = ai.getAndIncrement()) {
					  tileY = nTile /tilesX;
					  tileX = nTile % tilesX;
					  if (tileX==0) {
						  if (updateStatus) IJ.showStatus("Filtering noise rejection mask, row "+(tileY+1)+" of "+tilesY);
						  if (MASTER_DEBUG_LEVEL>1) System.out.println("Filtering noise rejection mask, row "+(tileY+1)+" of "+tilesY+" : "+IJ.d2s(0.000000001*(System.nanoTime()-startTime),3));
					  }
					  extractSquareTile( pixels, // source pixel array,
								  tile, // will be filled, should have correct size before call
								  slidingWindow, // window (same size as the kernel)
								  imgWidth, // width of pixels array
								  tileX*step, // left corner X
								  tileY*step); // top corner Y
					  fht_instance.swapQuadrants(tile);
					  fht_instance.transform(tile);
					  tile=fht_instance.multiply(tile,filter,false);
					  fht_instance.inverseTransform(tile);
					  fht_instance.swapQuadrants(tile);
					  /* accumulate result */
					  /*This is synchronized method. It is possible to make threads to write to non-overlapping regions of the OUT_PIXELS, but as the accumulation
					   * takes just small fraction of several FHTs, it should be OK - reasonable number of threads will spread and not "stay in line"
					   */
					  accumulateSquareTile(MASK_LOHIRES, //  float pixels array to accumulate tile
								  tile, // data to accumulate to the pixels array
								  imgWidth, // width of pixels array
								  tileX*step, // left corner X
								  tileY*step); // top corner Y
					  
				  }
			  }
		  };
	  }		      
	  startAndJoin(threads);
	  return MASK_LOHIRES;
  }
 
  
  
  
/* ======================================================================== */
    public double [] ringFilter(double [] dmask,       // mask to be filtered
    		                    int       width,       // mask width
    		                    double    minMaxValue, // min value for the local maximum to be processed (absolute, not relative)
    		                    double    overRingThreshold, // Ratio of the local maximum to maximal value in a ring to trigger filter
    		                    double    overRingLimit, // limit for the pixels in the center ring relative to the maximum in a ring
    		                    double    ringIR,        // ring inner radius
    		                    double    ringOR) {      // ring outer radius
    	if (dmask==null) return null;
    	int margins=(int)Math.ceil(ringIR);
    	int height=dmask.length/width;
    	double [] result= dmask.clone();
    	int i,j,nr,nc,y,x,n;
    	double ringIR2=ringIR*ringIR;
    	double ringOR2=ringOR*ringOR;
    	double r2;
    	nc=0;
    	nr=0;
    	for (i=-margins+1;i<margins;i++) for (j=-margins+1;j<margins;j++) {
    		r2=i*i+j*j;
    		if (r2<ringIR2) nc++;
    		else if (r2<=ringOR2) nr++;
    	}
    	if ((nc==0) || (nr==0)) return result; // do not filter; 
    	int [] indxc=new int[nc];
    	int [] indxr=new int[nr];
    	nc=0;
    	nr=0;
    	for (i=-margins+1;i<margins;i++) for (j=-margins+1;j<margins;j++) {
    		r2=i*i+j*j;
    		if (r2<ringIR2) {
    			indxc[nc++]=j+width*i;
    		}
    		else if (r2<=ringOR2) {
    			indxr[nr++]=j+width*i;
    		}
    	}
    	int [] neighb={-width,-width+1,1,width+1,width,width-1,-1,-width-1};
    	int index;
    	boolean isMax;
    	double d,maxInRing;
    	for (y=margins; y<height-margins;y++) for (x=margins; x<width-margins;x++) {
    		index=y*width+x;
    		d=dmask[index];
    		if (d<minMaxValue) continue; // too small value - don't bother to filter
    		isMax=true;
    		for (n=0;n<neighb.length;n++) if (dmask[index+neighb[n]]>d) {
    			isMax=false;
    			break;
    		}
    		if (!isMax) continue; // only process local maximums
    		maxInRing=dmask[index+indxr[0]];
    		for (n=1;n<nr;n++) if (dmask[index+indxr[n]]>maxInRing) maxInRing=dmask[index+indxr[n]];
    		if (d< (maxInRing*overRingThreshold)) continue; // under threshold, nop
// limit values in the circle
    		maxInRing*=overRingLimit;
    		for (n=0; n<nc;n++) if (dmask[index+indxc[n]]>maxInRing) result[index+indxc[n]]=maxInRing;
    	}
    	return result;
    }
/* ======================================================================== */
/* Combine two  3-slice image stacks generated from the same source image - one high-res/high noise, other low-res/low noise 
 * @param nonlinParameters TODO*/
  public ImageStack combineLoHiStacks(ImageStack        stack_convolved, // ImageStack with the image, convolved with the reversed PSF (sharp but with high noise)
                                      ImageStack         stack_gaussian, // ImageStack with the image, convolved with the Gaussian (just lateral compensated)  (blurred, but low noise)
                                      int                          nChn, // number of channel to apply to the min/max. If <0 - do not apply
                                      int                       nSubChn, // number of sub channel to apply to the min/max. If <0 - do not apply
                                      NonlinParameters nonlinParameters, // show mask generated and used
                                      final double [][]       noiseMask, // 2-d array of kernelsNoiseGain (divide mask by it)
                                      final int               noiseStep, // linear pixels per noiseMask pixels (32)
                                      final int              threadsMax, // maximal step in pixels on the maxRadius for 1 angular step (i.e. 0.5)
                                      final boolean        updateStatus){ // update status info

      int i,j,k;
      int imgWidth= stack_convolved.getWidth();
      int imgHeight=stack_convolved.getHeight();
      double [] diffGreens= new double [imgWidth*imgHeight];
      double [] diffGreens1;
      double filtMin=nonlinParameters.filtMin;
      double filtMax=nonlinParameters.filtMax;
      if ((nChn>=0) && (nSubChn>=0)){
    	  filtMin*=nonlinParameters.thresholdCorrection[nChn][nSubChn];
    	  filtMax*=nonlinParameters.thresholdCorrection[nChn][nSubChn];
      }
/* find number of the green channel - should be called "green", if none - use last */
      int greenChn=2;
      for (i=0;i<3;i++) if (stack_convolved.getSliceLabel(i+1).equals("green")){
        greenChn=i;
        break;
      }
      double d;
      double max=0.0f;
//      double average=0.0f;
	  DoubleGaussianBlur gb=new DoubleGaussianBlur();

      float [] hipassPixels=(float[]) stack_convolved.getPixels(greenChn+1);
      float [] lopassPixels=(float[]) stack_gaussian.getPixels(greenChn+1);
/*
      for (i=0;i<lopassPixels.length;i++) {
 
        d=hipassPixels[i]-lopassPixels[i];
        diffGreens[i]=d*d;
        if (max<lopassPixels[i]) max=lopassPixels[i];
      }
*/      
      for (i=0;i<lopassPixels.length;i++) {
    	  
//          d=hipassPixels[i]-lopassPixels[i];
//          diffGreens[i]=d*d;
    	  diffGreens[i]=hipassPixels[i]-lopassPixels[i];
      }
      if ((DEBUG_LEVEL>3) && nonlinParameters.showMask) SDFA_INSTANCE.showArrays(diffGreens.clone(), imgWidth, imgHeight,"diffGreens-nofilter");
      if (nonlinParameters.blurSigma>0)	{
	   	  if (DEBUG_LEVEL>1) System.out.println ( "Applying gaussian blur to difference hi/lo pass, blurSigma="+nonlinParameters.blurSigma);
		  gb.blurDouble(diffGreens, imgWidth, imgHeight, nonlinParameters.blurSigma, nonlinParameters.blurSigma, 0.01);
	  }
      if ((DEBUG_LEVEL>3) && nonlinParameters.showMask) SDFA_INSTANCE.showArrays(diffGreens.clone(), imgWidth, imgHeight,"diffGreens-blurred");
      for (i=0;i<lopassPixels.length;i++) {
    	  diffGreens[i]=diffGreens[i]*diffGreens[i];
      }
      if ((DEBUG_LEVEL>2) && nonlinParameters.showMask) SDFA_INSTANCE.showArrays(lopassPixels.clone(), imgWidth, imgHeight,"lopassPixels");
      
      for (i=0;i<lopassPixels.length;i++) {
    	  if (max<lopassPixels[i]) max=lopassPixels[i];
      }
   	  if (DEBUG_LEVEL>1) System.out.println ( "max(lopassPixels)="+max);
//      max*=((float) NONLIN_PARAMETERS.threshold);
// Make threshold absolute - when (blured) intensity is below thershold, the divisor is not decreasing
      max=((float) NONLIN_PARAMETERS.threshold);
      if ((DEBUG_LEVEL>3) && nonlinParameters.showMask) SDFA_INSTANCE.showArrays(diffGreens.clone(), imgWidth, imgHeight,"diffGreens-squared");
      for (i=0;i<lopassPixels.length;i++) {
        diffGreens[i]/=(float) Math.max(max,lopassPixels[i]);
      }
//      if ((DEBUG_LEVEL>3) && nonlinParameters.showMask) SDFA_INSTANCE.showArrays(diffGreens.clone(), imgWidth, imgHeight,"diffG-norm-limited");
      if ((DEBUG_LEVEL>1) && nonlinParameters.showMask) SDFA_INSTANCE.showArrays(diffGreens.clone(), imgWidth, imgHeight,"diffG-norm-limited");
      if (nonlinParameters.useRejectBlocksFilter) { // use frequency domain filtering
    	  
    	  double lowpassSigmaFreq=1.0*nonlinParameters.maskFFTSize/(2*Math.PI*nonlinParameters.lowPassSigma); // low pass sigma in frequency domain
    	  double [] filterFHT = createFilterForBlockArtifacts(
    			  nonlinParameters.maskFFTSize, // size of square FHT
    			  nonlinParameters.blockPeriod, // period (pixels) of the block artifacts to reject (32)
    			  nonlinParameters.rejectFreqSigma, // sigma of the rejection spots ( 0.0 - just zero a single point)
    			  lowpassSigmaFreq); // sigma of the low pass filter (frequency units, 0.0 - do not filter)
    	  if ((DEBUG_LEVEL>3) && nonlinParameters.showMask) SDFA_INSTANCE.showArrays(filterFHT,"filterFHT");
// Extend at least by half sliding window in each direction to reduce border effect 	  
    	  diffGreens1=extendDoubleArrayForSlidingWindow(
    			  diffGreens, // input pixel array
    			  imgWidth,  // width of the image
    			  nonlinParameters.maskFFTSize/2); // size of sliding step (half of the sliding window size)
    	  int extendedWidth=  extendDimension(imgWidth, (nonlinParameters.maskFFTSize/2));
          int extendedHeight= extendDimension(imgHeight,(nonlinParameters.maskFFTSize/2));
         
          if ((DEBUG_LEVEL>3) && nonlinParameters.showMask) SDFA_INSTANCE.showArrays(diffGreens1.clone(),extendedWidth,  extendedHeight,"diffGreens-extended");
    	  
// run block rejection filter
    	  diffGreens1=filterMaskFromBlockArtifacts(
    			  diffGreens1, // input pixel array
    			  extendedWidth, // width of the image
    			  extendedHeight, // width of the image
    			  nonlinParameters.maskFFTSize, // size of sliding FHT
    			  filterFHT, // filter to multiply FHT (created once for the whole filter mask)
    			  threadsMax, // maximal step in pixels on the maxRadius for 1 angular step (i.e. 0.5)
    			  updateStatus); // update status info
          if ((DEBUG_LEVEL>3) && nonlinParameters.showMask) SDFA_INSTANCE.showArrays(diffGreens1.clone(),extendedWidth,
        		  extendedHeight,"diffGreens-filtered-extended");/**/
// cut extra margins, crop to original size
    	  diffGreens1=reducedDoubleArrayAfterSlidingWindow(
    			  diffGreens1, // input pixel array
    			  imgWidth,  // width of the image
    			  imgHeight,
    			  nonlinParameters.maskFFTSize/2); // size of sliding step (half of the sliding window size)
          if ((DEBUG_LEVEL>2) && nonlinParameters.showMask) SDFA_INSTANCE.showArrays(diffGreens1.clone(),imgWidth,
        		  imgHeight,"diffGreens-filtered");
    	  if (nonlinParameters.combineBothModes) {
//      		DoubleGaussianBlur gb=new DoubleGaussianBlur();
    		gb.blurDouble(diffGreens, imgWidth, imgHeight, nonlinParameters.lowPassSigma, nonlinParameters.lowPassSigma, 0.01);
    		for (i=0;i<diffGreens.length;i++){
    			d=diffGreens[i]*diffGreens1[i];
    			diffGreens[i]=(d>0)?Math.sqrt(diffGreens[i]*diffGreens1[i]):0.0;
    		}
    	  } else {
    		  diffGreens=diffGreens1; 
    	  }
      } else { // just apply low-pass filter to the mask
//    		DoubleGaussianBlur gb=new DoubleGaussianBlur();
    		gb.blurDouble(diffGreens, imgWidth, imgHeight, nonlinParameters.lowPassSigma, nonlinParameters.lowPassSigma, 0.01);
      }
      if ((DEBUG_LEVEL>1) && nonlinParameters.showMask) SDFA_INSTANCE.showArrays(diffGreens.clone(), imgWidth, imgHeight,"diffGreens-filtered");
      
//      final double [][]       noiseMask, // 2-d array of kernelsNoiseGain (divide mask by it)
//      final int               noiseStep, // linear pixels per noiseMask pixels (32)
/* divide mask by noiseMask, if defined */
      if (noiseMask!=null) {
    	  if (DEBUG_LEVEL>1) System.out.println ( "diffGreens.length="+diffGreens.length+" imgWidth="+imgWidth+" noiseMask.length="+noiseMask.length+" noiseMask[0].length="+noiseMask[0].length);
    	  
          for (i=0;i<diffGreens.length;i++) {
        	  j=(i/imgWidth)/noiseStep;
        	  k=(i%imgWidth)/noiseStep;
        	  if (j>=noiseMask.length)    j=noiseMask.length-1;
        	  if (k>=noiseMask[j].length) k=noiseMask[j].length-1;
        	  diffGreens[i]/=noiseMask[j][k];
          }    	  
      }
      if ((DEBUG_LEVEL>1) && nonlinParameters.showMask) SDFA_INSTANCE.showArrays(diffGreens.clone(), imgWidth, imgHeight,"diffGreens-noise");
      if (nonlinParameters.useRingFilter) {
    	  diffGreens=ringFilter(diffGreens,            // mask to be filtered
    			  imgWidth,                            // mask width
    			  nonlinParameters.minMaxValue*nonlinParameters.filtMax, // min value for the local maximum to be processed (absolute, not relative)
    			  nonlinParameters.overRingThreshold,  // Ratio of the local maximum to maximal value in a ring to trigger filter
    			  nonlinParameters.overRingLimit,      // limit for the pixels in the center ring relative to the maximum in a ring
    			  nonlinParameters.ringIR,             // ring inner radius
    			  nonlinParameters.ringOR);            // ring outer radius
          if ((DEBUG_LEVEL>1) && nonlinParameters.showMask) SDFA_INSTANCE.showArrays(diffGreens.clone(), imgWidth, imgHeight,"diffGreens-ring");
      }
   	  if (DEBUG_LEVEL>1) System.out.println ( "filtMax="+filtMax+" filtMin="+filtMin);
      d= (float) ( 1.0/(filtMax-filtMin));
      if (filtMax>filtMin) {
        for (i=0;i<diffGreens.length;i++) {
          if (diffGreens[i]<filtMin) diffGreens[i]=0.0f;
          else if (diffGreens[i]>filtMax) diffGreens[i]=1.0f;
          else diffGreens[i]=d*(diffGreens[i]- (float) filtMin);
        }
      }
//      if (nonlinParameters.showMask) {
//    	  SDFA_INSTANCE.showArrays(diffGreens, imgWidth, imgHeight,"mask");
//      }
     DENOISE_MASK=diffGreens;
     DENOISE_MASK_WIDTH=imgWidth; 

/* Combine 2 stacks and a mask */
      return combineStacksWithMask (stack_gaussian,
                                   stack_convolved, 
                                        diffGreens);
   }

/* ======================================================================== */

  /* Combine 2 stacks and a mask */
  public ImageStack combineStacksWithMask (ImageStack stack_bg,
		  ImageStack stack_fg, 
		  //                                                 float [] mask ) {
		  double [] mask ) {

	  ImageStack stack=new ImageStack(stack_bg.getWidth(),stack_bg.getHeight());
	  int slice,i;
	  float [] fpixels;
	  float [] fpixels_bg;
	  float [] fpixels_fg;
	  for (slice=1; slice <=stack_bg.getSize(); slice++) {
		  fpixels_bg= (float[])stack_bg.getPixels(slice);
		  fpixels_fg= (float[])stack_fg.getPixels(slice);
		  fpixels=new float [fpixels_bg.length];
		  for (i=0;i<fpixels_bg.length;i++) fpixels[i]= (float) (mask[i]*fpixels_fg[i]+(1.0f-mask[i])*fpixels_bg[i]);
		  stack.addSlice(stack_fg.getSliceLabel(slice), fpixels);
	  }
	  return stack;
  }


/* ======================================================================== */
/* Convert source Bayer pattern (GR/BG) image to higher resolution, add margins by duplicating pattern around */
  public ImageStack  bayerToStack(ImagePlus imp, // source bayer image, linearized, 32-bit (float))
                                 SplitParameters splitParameters){

    if (imp==null) return null;
//    String [] chnNames={"red","blue","green"};
    String [] chnNames={"Red","Blue","Green"}; //Different sequence than RGB!!
    int nChn=chnNames.length;
    ImageProcessor ip=imp.getProcessor();
    int inWidth=imp.getWidth();
    int inHeight=imp.getHeight();
    int outHeight=inHeight*splitParameters.oversample+splitParameters.addTop+splitParameters.addBottom;
    int outWidth=inWidth*splitParameters.oversample+splitParameters.addLeft+splitParameters.addRight;
    int outLength=outWidth*outHeight;

    float [][] outPixels=new float[nChn][outLength];
    float [] pixels = (float[]) ip.getPixels();
    int chn,y,x,i,index;
    int bayerPeriod=2*splitParameters.oversample;
    int ovrWidth= inWidth*splitParameters.oversample;
    int ovrHeight=inHeight*splitParameters.oversample;
    for (chn=0;chn<nChn;chn++) for (i=0;i<outPixels[chn].length;i++) outPixels[chn][i]=0.0f;
/* Can be optimized - now it calculate input address for all those 0-es */
    for (index=0; index<outLength; index++) {
      y=(index / outWidth)-splitParameters.addTop;
      x=(index % outWidth)-splitParameters.addLeft;
      if (y<0) y= (bayerPeriod-((-y) % bayerPeriod))%bayerPeriod;
      else if (y>=ovrHeight) y= ovrHeight-bayerPeriod +((y-ovrHeight) % bayerPeriod);
      if (x<0) x= (bayerPeriod-((-x) % bayerPeriod))%bayerPeriod;
      else  if (x>=ovrWidth) x= ovrWidth-bayerPeriod +((x-ovrWidth) % bayerPeriod);
      if (((y% splitParameters.oversample)==0) && ((x% splitParameters.oversample)==0)) {
        x/=splitParameters.oversample;
        y/=splitParameters.oversample;
        chn=((x&1)==(y&1))?2:(((x&1)!=0)?0:1);
        outPixels[chn][index]=pixels[y*inWidth+x];
      }
    }
/* prepare result stack to return */
    ImageStack outStack=new ImageStack(outWidth,outHeight);
    for (chn=0;chn<nChn;chn++) {
      outStack.addSlice(chnNames[chn], outPixels[chn]);
    }
    return outStack;
}
/* ======================================================================== */
  void extractOneKernel(float [] pixels, //  array of combined square kernels, each 
                        double [] kernel, // will be filled, should have correct size before call
                              int numHor, // number of kernels in a row
                               int xTile, // horizontal number of kernel to extract
                               int yTile) { // vertical number of kernel to extract
    int length=kernel.length;
    int size=(int) Math.sqrt(length);
    int i,j;
    int pixelsWidth=numHor*size;
    int pixelsHeight=pixels.length/pixelsWidth;
    int numVert=pixelsHeight/size;
/* limit tile numbers - effectively add margins around the known kernels */
    if (xTile<0) xTile=0;
    else if (xTile>=numHor) xTile=numHor-1;
    if (yTile<0) yTile=0;
    else if (yTile>=numVert) yTile=numVert-1;
    int base=(yTile*pixelsWidth+xTile)*size;
    for (i=0;i<size;i++) for (j=0;j<size;j++) kernel [i*size+j]=pixels[base+i*pixelsWidth+j];
  }
/* ======================================================================== */
 /**extract and multiply by window function (same size as kernel itself) */
  void extractSquareTile(float [] pixels, // source pixel array,
                          double [] tile, // will be filled, should have correct size before call
                        double [] window, // window (same size as the kernel)
                               int width, // width of pixels array
                                  int x0, // left corner X
                                  int y0) { // top corner Y
    int length=tile.length;
    int size=(int) Math.sqrt(length);
    int i,j,x,y;
    int height=pixels.length/width;
    int index=0;
    for (i=0;i<size;i++) {
      y=y0+i;
      if ((y>=0) && (y<height)) {
        index=i*size;
        for (j=0;j<size;j++) {
         x=x0+j;
         if ((x>=0) && (x<width)) tile [index]=pixels[y*width+x]*window[index];
         index++;
        }
      }
    }
  }
/* ======================================================================== */
  void extractSquareTile(double [] pixels, // source pixel array,
		  double [] tile, // will be filled, should have correct size before call
		  double [] window, // window (same size as the kernel)
		  int width, // width of pixels array
		  int x0, // left corner X
		  int y0) { // top corner Y
	  int length=tile.length;
	  int size=(int) Math.sqrt(length);
	  int i,j,x,y;
	  int height=pixels.length/width;
	  int index=0;
	  for (i=0;i<size;i++) {
		  y=y0+i;
		  if ((y>=0) && (y<height)) {
			  index=i*size;
			  for (j=0;j<size;j++) {
				  x=x0+j;
				  if ((x>=0) && (x<width)) tile [index]=pixels[y*width+x]*window[index];
				  index++;
			  }
		  }
	  }
  }

  
/* ======================================================================== */
/* accumulate square tile to the pixel array (tile may extend beyond the array, will be cropped) */
  synchronized void  accumulateSquareTile(
		  float [] pixels, //  float pixels array to accumulate tile
		  double []  tile, // data to accumulate to the pixels array
		  int       width, // width of pixels array
		  int          x0, // left corner X
		  int          y0) { // top corner Y
	  int length=tile.length;
	  int size=(int) Math.sqrt(length);
	  int i,j,x,y;
	  int height=pixels.length/width;
	  int index=0;
	  for (i=0;i<size;i++) {
		  y=y0+i;
		  if ((y>=0) && (y<height)) {
			  index=i*size;
			  for (j=0;j<size;j++) {
				  x=x0+j;
				  if ((x>=0) && (x<width)) pixels[y*width+x]+=tile [index];
				  index++;
			  }
		  }
	  }
  }
  synchronized void  accumulateSquareTile(
		  double [] pixels, //  float pixels array to accumulate tile
		  double []  tile, // data to accumulate to the pixels array
		  int       width, // width of pixels array
		  int          x0, // left corner X
		  int          y0) { // top corner Y
	  int length=tile.length;
	  int size=(int) Math.sqrt(length);
	  int i,j,x,y;
	  int height=pixels.length/width;
	  int index=0;
	  for (i=0;i<size;i++) {
		  y=y0+i;
		  if ((y>=0) && (y<height)) {
			  index=i*size;
			  for (j=0;j<size;j++) {
				  x=x0+j;
				  if ((x>=0) && (x<width)) pixels[y*width+x]+=tile [index];
				  index++;
			  }
		  }
	  }
  }

/* ======================================================================== */
  double [] getSlidingMask(int size) {
    double [] mask = new double [size*size];
    double [] maskLine=new double [size];
    double k=2*Math.PI/size;
    int i,j,index;
    for (i=0;i<size;i++) maskLine[i]= 0.5*(1.0-Math.cos(i*k));
    index=0;
    for (i=0;i<size;i++) for (j=0;j<size;j++) mask[index++]=maskLine[i]*maskLine[j];
    return mask;
  }
/* ======================================================================== */

 /* Adds zero pixels around the image, "extending canvas" */

  public double [][] extendFFTInputTo (double[][] input_pixels,
                                              int newSize) {
    double [][] pixels=new double[input_pixels.length][];
    int i;
    for (i=0;i<pixels.length;i++) pixels[i]= extendFFTInputTo (input_pixels[i], newSize);
    return pixels;
  }
  public double [][] extendFFTInput (double[][] input_pixels,
                                              int subDivFreq) {
    double [][] pixels=new double[input_pixels.length][];
    int i;
    for (i=0;i<pixels.length;i++) pixels[i]= extendFFTInput (input_pixels[i], subDivFreq);
    return pixels;
  }
  public double [] extendFFTInputTo (double[] input_pixels,
                                               int newSize) {
    int subDivFreq=newSize/((int)Math.sqrt (input_pixels.length));
    return extendFFTInput (input_pixels,
                             subDivFreq);

  }
  public double [] extendFFTInput (double[] input_pixels,
                                          int subDivFreq) {
    if (input_pixels==null) return null;
    int width=(int) Math.sqrt(input_pixels.length);
    return extendFFTInput (input_pixels,
                                  width,   // width of the image
                             subDivFreq);
  }

  public double [] extendFFTInput (double[] input_pixels,
                                               int width,   // width of the image
                                          int subDivFreq) {
    if (input_pixels==null) return null;
    double [] pixels=new double[input_pixels.length*subDivFreq*subDivFreq];
    int j,base,x,y;
    int height=input_pixels.length/width;
    for (j=0;j<pixels.length;j++) pixels[j]=0.0;
    j=0;
    for (y=0;y<height;y++) {
      base=width*(subDivFreq-1)*(width*subDivFreq +1)/2+y*width*subDivFreq;
      for (x=0;x<width;x++) pixels[base+x]=input_pixels[j++];
    }
    return pixels;
  }

/* ======================================================================== */
  public double[][] normalizeAndWindow (double [][] pixels, double [] windowFunction) {
    return normalizeAndWindow (pixels, windowFunction, true);
  }

  public double[] normalizeAndWindow (double [] pixels, double [] windowFunction) {
    return normalizeAndWindow (pixels, windowFunction, true);
  }


  public double[][] normalizeAndWindow (double [][] pixels, double [] windowFunction, boolean removeDC) {
    int i;
    for (i=0;i<pixels.length;i++)  if (pixels[i]!=null) pixels[i]=normalizeAndWindow (pixels[i],  windowFunction, removeDC);
    return pixels;
  }


  public double[] normalizeAndWindow (double [] pixels, double [] windowFunction, boolean removeDC) {
    int j;
    double s=0.0;
    if (pixels==null) return null;
    if (removeDC) {
      for (j=0;j<pixels.length;j++) s+=pixels[j];
      s/=pixels.length;
    }
    for (j=0;j<pixels.length;j++) pixels[j]=(pixels[j]-s)*windowFunction[j];
    return pixels;
  }


/* ======================================================================== */
  public void  YPrPbToRGB(ImageStack stack,
                                 double Kr,        // 0.299;
                                 double Kb,        // 0.114;
                                int sliceY,
                               int slicePr,
                               int slicePb
                           ) {
      float [] fpixels_r= (float[]) stack.getPixels(1);
      float [] fpixels_g= (float[]) stack.getPixels(2);
      float [] fpixels_b= (float[]) stack.getPixels(3);
      float [] fpixels_Y= (float[]) stack.getPixels(sliceY);
      float [] fpixels_Pr=(float[]) stack.getPixels(slicePr);
      float [] fpixels_Pb=(float[]) stack.getPixels(slicePb);
      int length=fpixels_r.length;
      double Kg=1.0-Kr-Kb;
      int i;
/**
R= Y+ Pr*2.0*(1-Kr)
B= Y+ Pb*2.0*(1-Kb)
G= Y  +Pr*(- 2*Kr*(1-Kr))/Kg + Pb*(-2*Kb*(1-Kb))/Kg

*/
      double KPrR=  2.0*(1-Kr);
      double KPbB=  2.0*(1-Kb);
      double KPrG= -2.0*Kr*(1-Kr)/Kg;
      double KPbG= -2.0*Kb*(1-Kb)/Kg;
      double Y,Pr,Pb;
      for (i=0;i<length;i++) {
        Pb=fpixels_Pb[i];
        Pr=fpixels_Pr[i];
        Y =fpixels_Y [i];
        fpixels_r[i]=(float) (Y+ Pr*KPrR);
        fpixels_b[i]=(float) (Y+ Pb*KPbB);
        fpixels_g[i]=(float) (Y+ Pr*KPrG + Pb*KPbG);
      }
  }


/* ======================================================================== */
/* ======================================================================== */
/* ============== Dialogs ========================================================== */
  /*
			true, // crop - crop image to the sensor size
			true, // jpeg - convert to 8-bit RGB and save jpeg (if save is true)
			true, // save - save result
			true, // save16 - save 16-bit tiff also if the end result is 8 bit
			true, // save32 - save 32-bit tiff also if the end result is 8 or 16 bit

*/
  
  public boolean showProcessDialog(ProcessParameters processParameters) {
	  GenericDialog gd = new GenericDialog("Process parameters");
	  //	    gd.addCheckbox    ("Use first approximation for Y",           colorProcParameters.useFirstY);
	  int i,j;
      gd.addCheckbox ("Eyesis camera mode (3 composite frames)",          processParameters.eyesisMode);
	  for (i=0; i<processParameters.frames.length;i++) for (j=0;j<processParameters.frames[i].length;j++)
		  gd.addCheckbox ("Enable processing channel "+(i+1)+", subframe "+(j+1), processParameters.frames[i][j]);
      gd.addCheckbox ("Open file selection dialog (false - use selected)",processParameters.selectFile);
	  gd.addCheckbox ("Process only selected file (false - all 3)",       processParameters.thisFileOnly);
	  gd.addNumericField("Channel to process (if only selected file, 1..3)", processParameters.subChannelToProcess,0);
	  gd.addCheckbox ("Splt into Bayer stack (if false will exit)",       processParameters.split);
	  gd.addCheckbox ("De-mosaic (if false will exit)",                   processParameters.debayer);
	  gd.addCheckbox ("Show de-mosaic middle-frequency 'energy",          processParameters.showDebayerEnergy);
	  gd.addCheckbox ("Save de-mosaic middle-frequency 'energy",          processParameters.saveDebayerEnergy);
	  gd.addCheckbox ("Sharpen (convolve with calibration kernels)",      processParameters.deconvolve);
	  gd.addCheckbox ("Denoise (convolve with Gaussian in smooth areas)", processParameters.combine);
	  gd.addCheckbox ("Show denoise mask (white - use hi-res, black - low-res)", processParameters.showDenoiseMask);
	  gd.addCheckbox ("Save denoise mask (white - use hi-res, black - low-res)", processParameters.saveDenoiseMask);
	  gd.addCheckbox ("Show kernel noise gains",                          processParameters.showNoiseGains);
	  gd.addCheckbox ("Save kernel noise gains",                          processParameters.saveNoiseGains);
	  gd.addCheckbox ("Convert colors",                                   processParameters.colorProc);
	  gd.addCheckbox ("Show chroma denoise mask (white - use hi-res, black - low-res)", processParameters.showChromaDenoiseMask);
	  gd.addCheckbox ("Save chroma denoise mask (white - use hi-res, black - low-res)", processParameters.saveChromaDenoiseMask);
	  gd.addCheckbox ("Rotate result image",                              processParameters.rotate);
	  gd.addCheckbox ("Crop result image to the original size",           processParameters.crop);
	  gd.addCheckbox ("Convert to RGB48",                                 processParameters.toRGB);
	  gd.addCheckbox ("Convert to 8 bit RGB (and save JPEG if save is enabled)", processParameters.jpeg);
	  gd.addCheckbox ("Save the result to file system",                   processParameters.save);
	  gd.addCheckbox ("Save 16-bit tiff if the result is 8 bit",          processParameters.save16);
	  gd.addCheckbox ("Save 32-bit tiff if the result is 8 or 16 bit",    processParameters.save32);
	  gd.addCheckbox ("Show the result image",                            processParameters.show);
	  gd.addNumericField("JPEG quality (%)",                              processParameters.JPEG_quality,0);
	  gd.addNumericField("JPEG scale   (%)",                         100* processParameters.JPEG_scale,0);
      gd.addCheckbox ("Save current settings with results",               processParameters.saveSettings);
      gd.addCheckbox ("Update ImageJ status",                             UPDATE_STATUS);
      WindowTools.addScrollBars(gd);
	  gd.addNumericField("Debug Level:",                          MASTER_DEBUG_LEVEL,      0);
	  gd.showDialog();
	  if (gd.wasCanceled()) return false;
	  processParameters.eyesisMode=        gd.getNextBoolean();
	  for (i=0; i<processParameters.frames.length;i++) for (j=0;j<processParameters.frames[i].length;j++) 
		  processParameters.frames[i][j]=  gd.getNextBoolean();
	  processParameters.selectFile=        gd.getNextBoolean();
	  processParameters.thisFileOnly=      gd.getNextBoolean();
	  processParameters.subChannelToProcess=(int) gd.getNextNumber();
	  processParameters.split=             gd.getNextBoolean();
	  processParameters.debayer=           gd.getNextBoolean();
      processParameters.showDebayerEnergy= gd.getNextBoolean();
      processParameters.saveDebayerEnergy= gd.getNextBoolean();
	  processParameters.deconvolve=        gd.getNextBoolean();
	  processParameters.combine=           gd.getNextBoolean();
      processParameters.showDenoiseMask=   gd.getNextBoolean();
      processParameters.saveDenoiseMask=   gd.getNextBoolean();
      processParameters.showNoiseGains=    gd.getNextBoolean();
      processParameters.saveNoiseGains=    gd.getNextBoolean();
	  processParameters.colorProc=         gd.getNextBoolean();
      processParameters.showChromaDenoiseMask=   gd.getNextBoolean();
      processParameters.saveChromaDenoiseMask=   gd.getNextBoolean();
	  processParameters.rotate=            gd.getNextBoolean();
	  processParameters.crop=              gd.getNextBoolean();
	  processParameters.toRGB=             gd.getNextBoolean();
	  processParameters.jpeg=              gd.getNextBoolean();
	  processParameters.save=              gd.getNextBoolean();
	  processParameters.save16=            gd.getNextBoolean();
	  processParameters.save32=            gd.getNextBoolean();
	  processParameters.show=              gd.getNextBoolean();
	  processParameters.JPEG_quality=(int) gd.getNextNumber();
	  processParameters.JPEG_scale=   0.01*gd.getNextNumber();
	  processParameters.saveSettings=      gd.getNextBoolean();
      UPDATE_STATUS=                       gd.getNextBoolean();
	  MASTER_DEBUG_LEVEL=            (int) gd.getNextNumber();
	  return true;
  }  
/* ======================================================================== */
  
   public boolean showFilesDialog(FilesParameters filesParameters, ProcessParameters processParameters) {
	    GenericDialog gd = new GenericDialog("Kernel paths");
	    int i,j;
        gd.addCheckbox ("Eyesis camera mode (3 composite frames)",          processParameters.eyesisMode);
	    if (processParameters.eyesisMode) {
	    	for (i=0; i<filesParameters.rPSFNames.length;i++)
	    		for (j=0;j<filesParameters.rPSFNames[i].length;j++) if (filesParameters.rPSFNames[i][j]!=null) {
	    			gd.addStringField ("Deconvolution kernel name, channel "+(i+1)+", subframe "+(j+1), filesParameters.rPSFNames[i][j], 40);

	    		}
	    	for (i=0; i<filesParameters.rPSFNames.length;i++)
	    		for (j=0;j<filesParameters.gaussianNames[i].length;j++) if (filesParameters.gaussianNames[i][j]!=null) {
	    			gd.addStringField ("Gaussian kernel name, channel "+(i+1)+", subframe "+(j+1), filesParameters.gaussianNames[i][j], 40);

	    		}
	    } else {
			gd.addStringField ("Deconvolution kernel name ", filesParameters.rPSFNames[0][0], 40);
			gd.addStringField ("Gaussian kernel name",       filesParameters.gaussianNames[0][0], 40);
	    }
//sourceFile
	    if (filesParameters.sourceFiles!=null) {
	    	if (filesParameters.sourceFiles.length==1) {
	  		  gd.addStringField ("Source file(s) path (select one - will include other channels):    ", filesParameters.sourceFiles[0], 80);	    		
	        } else if (filesParameters.sourceFiles.length>1) {
	          gd.addMessage("Multiple source files selected");
	        }

	    } else {
	        gd.addMessage("Source file list is empty");
	    }
		gd.addStringField ("Kernel  directory path:    ", filesParameters.kernelDirectory, 80);
		gd.addStringField ("Results directory path:    ", filesParameters.resultsDirectory, 80);
        gd.addCheckbox    ("Use XML format to save/restore settings",  filesParameters.useXML);
        WindowTools.addScrollBars(gd);
	    gd.showDialog();
	    if (gd.wasCanceled()) return false;
		boolean newEyesisMode=           gd.getNextBoolean();

		
	    if (processParameters.eyesisMode) {
	    	for (i=0; i<filesParameters.rPSFNames.length;i++)
	    		for (j=0;j<filesParameters.rPSFNames[i].length;j++) if (filesParameters.rPSFNames[i][j]!=null) {
	    			filesParameters.rPSFNames[i][j]=gd.getNextString();

	    		}
	    	for (i=0; i<filesParameters.rPSFNames.length;i++)
	    		for (j=0;j<filesParameters.gaussianNames[i].length;j++) if (filesParameters.gaussianNames[i][j]!=null) {
	    			filesParameters.gaussianNames[i][j]=gd.getNextString();

	    		}
	    } else {
			filesParameters.rPSFNames[0][0]=            gd.getNextString();
			filesParameters.gaussianNames[0][0]=        gd.getNextString();
	    }
	    if ((filesParameters.sourceFiles!=null) && (filesParameters.sourceFiles.length==1)) {
		    filesParameters.sourceFiles[0]=             gd.getNextString();
	    }
		filesParameters.kernelDirectory=                gd.getNextString();
		filesParameters.resultsDirectory=               gd.getNextString();
		filesParameters.useXML=                         gd.getNextBoolean();
	    if (newEyesisMode!=processParameters.eyesisMode) {
	    	processParameters.eyesisMode=newEyesisMode;
	    	return showFilesDialog(filesParameters,processParameters);
	    }
	  return true;
  }

  /* ======================================================================== */
  
  public boolean showRGBProcessDialog(RGBParameters rgbParameters) {
    GenericDialog gd = new GenericDialog("RGB conversion parameters");

    gd.addNumericField("Red   color black level",             rgbParameters.r_min,     3);
    gd.addNumericField("Green color black level",             rgbParameters.g_min,     3);
    gd.addNumericField("Blue  color black level",             rgbParameters.b_min,     3);
    gd.addNumericField("Red   color white level",             rgbParameters.r_max,     3);
    gd.addNumericField("Green color white level",             rgbParameters.g_max,     3);
    gd.addNumericField("Blue  color white level",             rgbParameters.b_max,     3);
    gd.addNumericField("Alpha channel min",                   rgbParameters.alpha_min, 3);
    gd.addNumericField("Alpha channel max",                   rgbParameters.alpha_max, 3);
    
    gd.showDialog();
    if (gd.wasCanceled()) return false;
    rgbParameters.r_min=                gd.getNextNumber();
    rgbParameters.g_min=                gd.getNextNumber();
    rgbParameters.b_min=                gd.getNextNumber();
    rgbParameters.r_max=                gd.getNextNumber();
    rgbParameters.g_max=                gd.getNextNumber();
    rgbParameters.b_max=                gd.getNextNumber();
    rgbParameters.alpha_min=            gd.getNextNumber();
    rgbParameters.alpha_max=            gd.getNextNumber();
    return true;
  }

/* ======================================================================== */
  public boolean showColorProcessDialog(ColorProcParameters colorProcParameters) {
	    GenericDialog gd = new GenericDialog("Color processing parameters");
	    gd.addNumericField("YCbCr gamma",                             colorProcParameters.gamma,     3); //0.53
	    gd.addNumericField("Color saturation, red.green",             colorProcParameters.saturationRed,3); //2.0
	    gd.addNumericField("Color saturation, blue/green",            colorProcParameters.saturationBlue,3); //2.0
	    gd.addNumericField("Color balance, red-to-green",             colorProcParameters.balanceRed,     3); //1.8
	    gd.addNumericField("Color balance, blue-to-green",            colorProcParameters.balanceBlue,    3); //1.8
	    gd.addNumericField("Gain green",                              colorProcParameters.gain,      3); //1.8
	    gd.addNumericField("Weight scale RED  (which color to use)",  colorProcParameters.weightScaleR,  3); //1.8
	    gd.addNumericField("Weight scale BLUE (which color to use)",  colorProcParameters.weightScaleB,  3); //1.8
	    gd.addNumericField("Minimal linear value to apply gammaY",    colorProcParameters.minLin,    3); //0.53
	    gd.addNumericField("YCbCr Kb",                                colorProcParameters.kb,        3); //0.114
	    gd.addNumericField("YCbCr Kr",                                colorProcParameters.kr,        3); //0.299
	    gd.addCheckbox    ("Use first approximation for Y",           colorProcParameters.useFirstY);
	    gd.addMessage("Color denoise parameters");
	    gd.addNumericField("Low-pass sigma for luma to calculate mask",colorProcParameters.maskSigma,        3);
	    gd.addNumericField("Filtered luma minimal level for mask transition",colorProcParameters.maskMin,        3);
	    gd.addNumericField("Filtered luma maximal level for mask transition",colorProcParameters.maskMax,        3);
	    gd.addCheckbox    ("Combine chroma mask with sharpness one (reduce color leak)",  colorProcParameters.combineWithSharpnessMask);
	    gd.addNumericField("Low-pass sigma for chroma in bright areas",colorProcParameters.chromaBrightSigma,        3);
	    gd.addNumericField("Low-pass sigma for chroma in dark areas",colorProcParameters.chromaDarkSigma,        3);
	    gd.addNumericField("Debug Level:",                          MASTER_DEBUG_LEVEL,      0);
	    WindowTools.addScrollBars(gd);
	    gd.showDialog();
	    if (gd.wasCanceled()) return false;

	    colorProcParameters.gamma=                    gd.getNextNumber();
	    colorProcParameters.saturationRed=            gd.getNextNumber();
	    colorProcParameters.saturationBlue=           gd.getNextNumber();
	    colorProcParameters.balanceRed=               gd.getNextNumber();
	    colorProcParameters.balanceBlue=              gd.getNextNumber();
	    colorProcParameters.gain=                     gd.getNextNumber();
	    colorProcParameters.weightScaleR=             gd.getNextNumber();
	    colorProcParameters.weightScaleB=             gd.getNextNumber();
	    colorProcParameters.minLin=                   gd.getNextNumber();
	    colorProcParameters.kb=                       gd.getNextNumber(); //---
	    colorProcParameters.kr=                       gd.getNextNumber();
	    colorProcParameters.useFirstY=                gd.getNextBoolean();
	    colorProcParameters.maskSigma=                gd.getNextNumber();
	    colorProcParameters.maskMin=                  gd.getNextNumber();
	    colorProcParameters.maskMax=                  gd.getNextNumber();
	    colorProcParameters.combineWithSharpnessMask= gd.getNextBoolean();
	    colorProcParameters.chromaBrightSigma=        gd.getNextNumber();
	    colorProcParameters.chromaDarkSigma=          gd.getNextNumber();
	    
	    MASTER_DEBUG_LEVEL=      (int) gd.getNextNumber();
	    return true;
	  }
  /* ======================================================================== */
  public boolean showColorCalibDialog(ColorCalibParameters colorCalibParameters) {
		int i,j;
	    GenericDialog gd = new GenericDialog("Individual channels colors/gains");
	    for (i=0; i<colorCalibParameters.gain.length;i++) for (j=0;j<colorCalibParameters.gain[i].length;j++)
	    	   gd.addNumericField("Gain (brightness) for channel "+(i+1)+" subframe "+(j+1), colorCalibParameters.gain[i][j],  3);
	    for (i=0; i<colorCalibParameters.balanceRed.length;i++) for (j=0;j<colorCalibParameters.balanceRed[i].length;j++)
	    	   gd.addNumericField("Balance Red/Green for channel "+(i+1)+" subframe "+(j+1), colorCalibParameters.balanceRed[i][j],  3);
	    for (i=0; i<colorCalibParameters.balanceBlue.length;i++) for (j=0;j<colorCalibParameters.balanceBlue[i].length;j++)
	    	   gd.addNumericField("Balance Blue/Green for channel "+(i+1)+" subframe "+(j+1), colorCalibParameters.balanceBlue[i][j],  3);
	    WindowTools.addScrollBars(gd);
	    gd.showDialog();
	    if (gd.wasCanceled()) return false;

	    for (i=0; i<colorCalibParameters.gain.length;i++) for (j=0;j<colorCalibParameters.gain[i].length;j++)
	    	   colorCalibParameters.gain[i][j]=gd.getNextNumber();
	    for (i=0; i<colorCalibParameters.balanceRed.length;i++) for (j=0;j<colorCalibParameters.balanceRed[i].length;j++)
	    	   colorCalibParameters.balanceRed[i][j]=gd.getNextNumber();
	    for (i=0; i<colorCalibParameters.balanceBlue.length;i++) for (j=0;j<colorCalibParameters.balanceBlue[i].length;j++)
	    	   colorCalibParameters.balanceBlue[i][j]=gd.getNextNumber();
	    return true;
	  }
///    	showColorCalibDialog(COLOR_CALIB_PARAMETERS);

/* ======================================================================== */
  
  public boolean showCombinePairDialog(NonlinParameters nonlinParameters, ProcessParameters processParameters) {
	int i,j;
    GenericDialog gd = new GenericDialog("Combining high-re and low-res images");
    gd.addNumericField("Mask low-pass filter mask sigma (spatial domain)" ,  nonlinParameters.lowPassSigma,     3); // 5.0-  sigma for the nonlinear filtering
    gd.addCheckbox    ("Apply block rejection filter to the mask",           nonlinParameters.useRejectBlocksFilter); 
    gd.addCheckbox    ("Combine masks withe blocks rejected and not",        nonlinParameters.combineBothModes); 
    gd.addNumericField("Sliding FFT size for block filtering",               nonlinParameters.maskFFTSize,      0);
    gd.addNumericField("JPEG block period adjusted to oversampling (32)",    nonlinParameters.blockPeriod,      0);
    gd.addNumericField("Mask block rejection filter sigma (freq. domain)" ,  nonlinParameters.rejectFreqSigma,     3); // 1.0-
    gd.addNumericField("Nonlinear filter mask min. level",      nonlinParameters.filtMin,              3); // 0.01  minimal low-pass filtered squared difference between the corrected and original pixels to trigger sharpness enhancement
    gd.addNumericField("Nonlinear filter mask max. level",      nonlinParameters.filtMax,              3); // 0.15 squared low-pass filtered difference between the corrected and original pixels, so abopve that level 100% corrected image is used
    for (i=0; i<nonlinParameters.thresholdCorrection.length;i++) for (j=0;j<nonlinParameters.thresholdCorrection[i].length;j++)
    	   gd.addNumericField("Filter Min/Max scale for channel "+(i+1)+" subframe "+(j+1), nonlinParameters.thresholdCorrection[i][j],  3);
    gd.addNumericField("Nonlinear filter threshold",            nonlinParameters.threshold,        3); // 0.01 when blurred intensity is below this value, use it as a denominator
    gd.addCheckbox    ("Show generated/used mask",              processParameters.showDenoiseMask);   
    gd.addCheckbox    ("Save generated/used mask",              processParameters.saveDenoiseMask);
	gd.addCheckbox    ("Use differntial noise gains",           nonlinParameters.useDiffNoiseGains);
    gd.addNumericField("Noise gain weight 0 (red)",             nonlinParameters.noiseGainWeights[0], 3); // Weight of red component noise gain (0.0)
    gd.addNumericField("Noise gain weight 1 (blue)",             nonlinParameters.noiseGainWeights[1], 3); // Weight of blue component noise gain (0.0)
    gd.addNumericField("Noise gain weight 2 (green)",             nonlinParameters.noiseGainWeights[2], 3); // Weight of green component noise gain (1.0)
    gd.addNumericField("Blur kernels for noise calculation",    nonlinParameters.blurSigma, 3); // // blur sigma for mask calculation (blur convolution kernels for noise gain calculation
	gd.addNumericField("Noise Gain power",                      nonlinParameters.noiseGainPower,3);
	
	gd.addCheckbox    ("Filter out spots on denoise mask",      nonlinParameters.useRingFilter);    // filter out spots on denoise mask
	gd.addNumericField("Minimal relative value of local max",   nonlinParameters.minMaxValue,3);       // minimal value (relative to filtMax) of the local maximum to be processed
	gd.addNumericField("Relative max over ring threshold",      nonlinParameters.overRingThreshold,3); // ratio of local max. and maximal value in the surrounding ring to trigger filter 
	gd.addNumericField("Limit center circle over ring max",     nonlinParameters.overRingLimit,3);     // limit values in the center circle to scaled maximum in a ring
	gd.addNumericField("Ring inner radius",                     nonlinParameters.ringIR,3);            // ring inner radius (center circle radius)
	gd.addNumericField("Ring outer radius",                     nonlinParameters.ringOR,3);            // ring outer radius
    gd.addNumericField("Debug Level:",                          MASTER_DEBUG_LEVEL,      0);
	WindowTools.addScrollBars(gd);
    gd.showDialog();
    if (gd.wasCanceled()) return false;
    nonlinParameters.lowPassSigma=                gd.getNextNumber();
    nonlinParameters.useRejectBlocksFilter=       gd.getNextBoolean();
    nonlinParameters.combineBothModes=            gd.getNextBoolean();
    nonlinParameters.maskFFTSize=           (int) gd.getNextNumber();
    nonlinParameters.blockPeriod=           (int) gd.getNextNumber();
    nonlinParameters.rejectFreqSigma=             gd.getNextNumber();
    nonlinParameters.filtMin=                     gd.getNextNumber();
    nonlinParameters.filtMax=                     gd.getNextNumber();
    for (i=0; i<nonlinParameters.thresholdCorrection.length;i++) for (j=0;j<nonlinParameters.thresholdCorrection[i].length;j++)
 	   nonlinParameters.thresholdCorrection[i][j]=gd.getNextNumber();
    nonlinParameters.threshold=                   gd.getNextNumber();
    processParameters.showDenoiseMask=            gd.getNextBoolean();   
    processParameters.saveDenoiseMask=            gd.getNextBoolean();
    nonlinParameters.useDiffNoiseGains=           gd.getNextBoolean();
    nonlinParameters.noiseGainWeights[0]=         gd.getNextNumber();
    nonlinParameters.noiseGainWeights[1]=         gd.getNextNumber();
    nonlinParameters.noiseGainWeights[2]=         gd.getNextNumber();
    nonlinParameters.blurSigma=                   gd.getNextNumber();
    nonlinParameters.noiseGainPower=              gd.getNextNumber();
    nonlinParameters.useRingFilter=               gd.getNextBoolean();
	nonlinParameters.minMaxValue=                 gd.getNextNumber();
	nonlinParameters.overRingThreshold=           gd.getNextNumber(); 
	nonlinParameters.overRingLimit=               gd.getNextNumber();
	nonlinParameters.ringIR=                      gd.getNextNumber();
	nonlinParameters.ringOR=                      gd.getNextNumber();
    MASTER_DEBUG_LEVEL=                     (int) gd.getNextNumber();
    return true;
  }

/* ======================================================================== */
  public boolean showStackConvolutionDialog() {
    int i;
    GenericDialog gd = new GenericDialog("Stack convolution parameters");
    gd.addNumericField("Convolution FFT size (twice the kernel size)",                             CONVOLVE_FFT_SIZE,     0); // 128
    gd.addNumericField("Maximal number of concurrent threads",                                     THREADS_MAX, 0); //   100
    gd.addCheckbox    ("Update ImageJ status",                                                         UPDATE_STATUS);

    gd.addNumericField("Debug Level:",                          MASTER_DEBUG_LEVEL,      0);
    gd.showDialog();
    if (gd.wasCanceled()) return false;
    CONVOLVE_FFT_SIZE=1;
    for (i=(int) gd.getNextNumber(); i >1; i>>=1) CONVOLVE_FFT_SIZE <<=1; /* make it to be power of 2 */
    UPDATE_STATUS=               gd.getNextBoolean();
    THREADS_MAX=           (int) gd.getNextNumber();    
    MASTER_DEBUG_LEVEL=    (int) gd.getNextNumber();
    return true;
 }  
/* ======================================================================== */
  public boolean showDeBayerDialog(DebayerParameters debayerParameters, ProcessParameters processParameters) {
    int i;
    GenericDialog gd = new GenericDialog("De-bayer parameters");
    gd.addNumericField("Debayer threshold (lower use default filtering)",     debayerParameters.debayerThreshold,3) ; //=0.2; Measuring maximal spectral component amplitude in mid-frequencies. If below this - use default de-bayer mask
    gd.addNumericField("Debayer lo-pass relative width for green",            debayerParameters.debayerRelativeWidthGreen,   3); //1.5 result green mask mpy by scaled default (diamond)
    gd.addNumericField("Debayer lo-pass relative width for red/blue",         debayerParameters.debayerRelativeWidthRedblue, 3); //1.5 result red/blue mask mpy by scaled default (square)
    gd.addNumericField("Debayer lo-pass relative width for red/blue (main)",  debayerParameters.debayerRelativeWidthRedblueMain, 3); //1.3 green mask when applied to red/blue, main (center)
    gd.addNumericField("Debayer lo-pass relative width for red/blue (clones)",debayerParameters.debayerRelativeWidthRedblueClones, 3); //2.0 green mask when applied to red/blue, clones 
    gd.addNumericField("Debayer mask - power for the ampliutude",             debayerParameters.debayerGamma,   3); //0.3  
    gd.addNumericField("relative increase value with radius",                 debayerParameters.debayerBonus,   3); //0.5 - scale far pixels as (1.0+bonus*r/rmax)
    gd.addNumericField("Relative alais strength to mask out point",           debayerParameters.mainToAlias,3); //0.5; // relative main/alias amplitudes to enable lixels (i.e. 0.5 means that if alias is >0.5*main, the pixel will be masked out)
    gd.addNumericField("Gaussian blur sigma for the alias-rejecting masks",   debayerParameters.debayerMaskBlur,   3); //2.0

    gd.addCheckbox    ("Use 'scissors' filter",                               debayerParameters.debayerUseScissors);    // true; // use "scissors", if false - just apply "diamond" ands "square" with debayerParameters.debayerRelativeWidthGreen and debayerParameters.debayerRelativeWidthRedblue
    gd.addCheckbox    ("Show mid frequency components energy plots",         processParameters.showDebayerEnergy);    // false - plot debayer high frequency energy (use to select between "scissors" and default uniform
    gd.addCheckbox    ("Save mid frequency components energy plots",         processParameters.saveDebayerEnergy);    // false - plot debayer high frequency energy (use to select between "scissors" and default uniform
    
    
//	   

    gd.addNumericField("Polar grid largest cell size to cartesian one ratio", debayerParameters.polarStep, 3); //   0.5;// size of largest polar cell to cartesian one
    gd.addNumericField("Debayer FFT Size (64)",                               debayerParameters.size,   0); // 64
    
    gd.addCheckbox    ("Debug: show data for selected tile",                  debayerParameters.debug); // false
    gd.addNumericField("Debug: X-coordinate of the point of interest  (oversampled)",   debayerParameters.xDebug, 0); 
    gd.addNumericField("Debug: Y-coordinate of the point of interest  (oversampled)",   debayerParameters.yDebug, 0);
    
    gd.addNumericField("Maximal number of concurrent threads",                THREADS_MAX, 0); //   100
    gd.addCheckbox    ("Update ImageJ status",                                UPDATE_STATUS);
    gd.addCheckbox    ("Show debayer debug images as stacks (false - individual)", debayerParameters.debayerStacks); // true

    gd.addNumericField("Debug Level:",                                         MASTER_DEBUG_LEVEL, 0);
    WindowTools.addScrollBars(gd);
    gd.showDialog();
    if (gd.wasCanceled()) return false;
    debayerParameters.debayerThreshold=            gd.getNextNumber();
    debayerParameters.debayerRelativeWidthGreen=          gd.getNextNumber();
    debayerParameters.debayerRelativeWidthRedblue=        gd.getNextNumber();
    debayerParameters.debayerRelativeWidthRedblueMain=   gd.getNextNumber();
    debayerParameters.debayerRelativeWidthRedblueClones= gd.getNextNumber();

    debayerParameters.debayerGamma=                gd.getNextNumber();
    debayerParameters.debayerBonus=                gd.getNextNumber();
    debayerParameters.mainToAlias=                 gd.getNextNumber();
    debayerParameters.debayerMaskBlur=             gd.getNextNumber();

    debayerParameters.debayerUseScissors=          gd.getNextBoolean();
//    debayerParameters.showEnergy=                  gd.getNextBoolean();    
    processParameters.showDebayerEnergy=           gd.getNextBoolean();
    processParameters.saveDebayerEnergy=           gd.getNextBoolean();
    debayerParameters.polarStep=                   gd.getNextNumber();

    debayerParameters.size=1;
    for (i=(int) gd.getNextNumber(); i >1; i>>=1) debayerParameters.size <<=1; /* make it to be power of 2 */
    
    debayerParameters.debug=                gd.getNextBoolean();
    debayerParameters.xDebug=        (int) gd.getNextNumber(); 
    debayerParameters.yDebug=        (int) gd.getNextNumber();

    
    THREADS_MAX=            (int) gd.getNextNumber();    
    UPDATE_STATUS=                gd.getNextBoolean();
    debayerParameters.debayerStacks=               gd.getNextBoolean();
    MASTER_DEBUG_LEVEL=     (int) gd.getNextNumber();
    return true;
  }
/* ======================================================================== */
  public boolean showSplitBayerToStackDialog(SplitParameters splitParameters) {
    GenericDialog gd = new GenericDialog("Interpolate kernels parameters");
    gd.addNumericField("Interpolation step between original kernels",                                  splitParameters.oversample, 0); //2
    gd.addNumericField("Interpolation add on the top (in output, subdivided steps)",                   splitParameters.addTop,     0); //32
    gd.addNumericField("Interpolation add on the left (in output, subdivided steps)",                  splitParameters.addLeft,    0); //32
    gd.addNumericField("Interpolation add on the right (in output, subdivided steps)",                 splitParameters.addRight,   0); //32
    gd.addNumericField("Interpolation add on the bottom (in output, subdivided steps)",                splitParameters.addBottom,  0); //32

    gd.addNumericField("Debug Level:",                          MASTER_DEBUG_LEVEL,      0);
    gd.showDialog();
    if (gd.wasCanceled()) return false;
    splitParameters.oversample=   (int) gd.getNextNumber();
    splitParameters.addTop=       (int) gd.getNextNumber();
    splitParameters.addLeft=      (int) gd.getNextNumber();
    splitParameters.addRight=     (int) gd.getNextNumber();
    splitParameters.addBottom=    (int) gd.getNextNumber();
    MASTER_DEBUG_LEVEL= (int) gd.getNextNumber();
    return true;
 }  

/* ======================================================================== */
/* ======================================================================== */
/* ======================================================================== */
/* ======================================================================== */
/* unfinished/not used */
/* ======================================================================== 
 * @param colorProcParameters TODO*/
  public void  processColorsWeights(ImageStack stack,
                                        double scale,     // initila maximal pixel value (16))
                                  ColorProcParameters  colorProcParameters,
                                  ColorCalibParameters colorCalibParameters, // if null - assume all 1-s
        						  int nChn,
        						  int nSubChn
                           ) {
	  double thisGain=       colorProcParameters.gain;
	  double thisBalanceRed= colorProcParameters.balanceRed;
	  double thisBalanceBlue=colorProcParameters.balanceBlue;
	  if (colorCalibParameters!=null) {
		  thisGain*=       colorCalibParameters.gain[nChn][nSubChn];
		  thisBalanceRed*= colorCalibParameters.balanceRed[nChn][nSubChn];
		  thisBalanceBlue*=colorCalibParameters.balanceBlue[nChn][nSubChn];
	  }
      float [] fpixels_r= (float[]) stack.getPixels(1);
      float [] fpixels_g= (float[]) stack.getPixels(2);
      float [] fpixels_b= (float[]) stack.getPixels(3);
      boolean useWeights=(stack.getSize()>=5);
      if (!useWeights) {
        stack.addSlice("dummy1",  fpixels_r);
        stack.addSlice("dummy2",  fpixels_g);
      }
      float [] fpixels_wr=(float[]) stack.getPixels(4);
      float [] fpixels_wb=(float[]) stack.getPixels(5);
      int length=fpixels_r.length;
      int width= stack.getWidth();
      int height=stack.getHeight();
/* Scale colors, gamma-convert */
      int i;
      double gain_red= thisBalanceRed* thisGain/scale;
      double gain_blue=thisBalanceBlue*thisGain/scale;
      double gain_green=thisGain/scale;
      double gamma_a=Math.pow(colorProcParameters.minLin,colorProcParameters.gamma)*(1.0-colorProcParameters.gamma);
      gamma_a=gamma_a/(1.0-gamma_a);
      double gamma_linK=(1.0+gamma_a)*colorProcParameters.gamma*Math.pow(colorProcParameters.minLin,colorProcParameters.gamma)/colorProcParameters.minLin;



      for (i=0;i<length;i++) {
        fpixels_r[i]=(float) linGamma(colorProcParameters.gamma, gamma_a, gamma_linK, colorProcParameters.minLin, fpixels_r[i]*gain_red);
        fpixels_g[i]=(float) linGamma(colorProcParameters.gamma, gamma_a, gamma_linK, colorProcParameters.minLin, fpixels_g[i]*gain_green);
        fpixels_b[i]=(float) linGamma(colorProcParameters.gamma, gamma_a, gamma_linK, colorProcParameters.minLin, fpixels_b[i]*gain_blue);
      }
/* Convert to YPbPr */
      double Y,Pb,Pr;
      double Kg=1.0-colorProcParameters.kr-colorProcParameters.kb;
      double Sb=0.5/(1.0-colorProcParameters.kb)*colorProcParameters.saturationBlue;
      double Sr=0.5/(1.0-colorProcParameters.kr)*colorProcParameters.saturationRed;
      double Yr,Yg,Yb,Wr,Wg,Wb,S;
/* coefficients to find Y from Pb, Pr and a color (R,G or B)
 Yr = R- Pr*KPrR
 Yb = B- Pb*KPbB
 Yg = G+ Pr*KPrG  + Pb*KPbG
 */
      double KPrR= -(2.0*(1-colorProcParameters.kr))/colorProcParameters.saturationRed;
      double KPbB= -(2.0*(1-colorProcParameters.kb))/colorProcParameters.saturationBlue;
      double KPrG=  2.0*colorProcParameters.kr*(1-colorProcParameters.kr)/Kg/colorProcParameters.saturationRed;
      double KPbG=  2.0*colorProcParameters.kb*(1-colorProcParameters.kb)/Kg/colorProcParameters.saturationBlue;
      if (DEBUG_LEVEL>1) {
          System.out.println ( " processColorsWeights() gain_red="+gain_red+" gain_green="+gain_green+" gain_blue="+gain_blue);
          System.out.println ( " processColorsWeights() gamma="+colorProcParameters.gamma+      " minLin="+colorProcParameters.minLin+" gamma_a="+gamma_a+" gamma_linK="+gamma_linK);
          System.out.println ( " processColorsWeights() Kr="+colorProcParameters.kr+" Kg="+Kg+" Kb="+colorProcParameters.kb+" Sr="+Sr+" Sb="+Sb);
          System.out.println ( " processColorsWeights() KPrR="+KPrR+" KPbB="+KPbB+" KPrG="+KPrG+" KPbG="+KPbG);
      }

      float [] fpixels_pb= new float [length];
      float [] fpixels_pr= new float [length];
      float [] fpixels_y0= new float [length];
      float [] fpixels_y= fpixels_y0;

      float [] fpixels_yR=null;
      float [] fpixels_yG=null;
      float [] fpixels_yB=null;

      if (DEBUG_LEVEL>2) {
        fpixels_yR= new float [length];
        fpixels_yG= new float [length];
        fpixels_yB= new float [length];
      }
      for (i=0;i<length;i++) {
        Y=colorProcParameters.kr*fpixels_r[i]+Kg*fpixels_g[i]+colorProcParameters.kb*fpixels_b[i];
        fpixels_pb[i] = (float) (Sb*(fpixels_b[i]-Y));
        fpixels_pr[i] = (float) (Sr*(fpixels_r[i]-Y));
        fpixels_y0[i]=(float) Y;
      }
/* calculate Y from weighted colors, weights derived from how good each color component predicts signal in each subpixel of Bayer pattern */
      if (useWeights) {
        fpixels_y=  new float [length];
        for (i=0;i<length;i++) {
          Pb=fpixels_pb[i];
          Pr=fpixels_pr[i];
          Yr = fpixels_r[i]- Pr*KPrR;
          Yb = fpixels_b[i]- Pb*KPbB;
          Yg = fpixels_g[i]+ Pr*KPrG  + Pb*KPbG;
          Wr=fpixels_wr[i];
          Wb=fpixels_wb[i];
          Wg=1.0-Wr-Wb;
          S=1.0/(Wr*(colorProcParameters.weightScaleR-1.0)+Wb*(colorProcParameters.weightScaleB-1.0)+1.0);
          Wr*=S*colorProcParameters.weightScaleR;
          Wb*=S*colorProcParameters.weightScaleB;
          Wg*=S;
          Y=Yr*Wr+Yg*Wg+Yb*Wb;
          fpixels_y[i]=(float) Y;
          if (DEBUG_LEVEL>2) {
            fpixels_yR[i]= (float) Yr;
            fpixels_yG[i]= (float) Yg;
            fpixels_yB[i]= (float) Yb;
          }
        }
      }
/* Low-pass filter Pb and Pr */
	  DoubleGaussianBlur gb=new DoubleGaussianBlur();
      double [] dpixels_pr=new double[fpixels_pr.length];
      double [] dpixels_pb=new double[fpixels_pb.length];
      for (i=0;i<dpixels_pr.length;i++) {
    	  dpixels_pr[i]=fpixels_pr[i];
    	  dpixels_pb[i]=fpixels_pb[i];
      }
      if (colorProcParameters.maskMax>0.0) {
          double [] dmask=new double[fpixels_y0.length];
          for (i=0;i<dpixels_pr.length;i++)  dmask[i]=fpixels_y0[i];
          double [] dpixels_pr_dark=dpixels_pr.clone();
          double [] dpixels_pb_dark=dpixels_pb.clone();
    	  gb.blurDouble(dmask, width, height, colorProcParameters.maskSigma, colorProcParameters.maskSigma, 0.01);
    	  gb.blurDouble(dpixels_pr, width, height, colorProcParameters.chromaBrightSigma, colorProcParameters.chromaBrightSigma, 0.01);
    	  gb.blurDouble(dpixels_pb, width, height, colorProcParameters.chromaBrightSigma, colorProcParameters.chromaBrightSigma, 0.01);
    	  gb.blurDouble(dpixels_pr_dark, width, height, colorProcParameters.chromaDarkSigma, colorProcParameters.chromaDarkSigma, 0.01);
    	  gb.blurDouble(dpixels_pb_dark, width, height, colorProcParameters.chromaDarkSigma, colorProcParameters.chromaDarkSigma, 0.01);
          if (DEBUG_LEVEL>2) {
        	  SDFA_INSTANCE.showArrays(dmask, width, height,"dmask");          
        	  SDFA_INSTANCE.showArrays(dpixels_pr, width, height,"dpixels_pr");          
        	  SDFA_INSTANCE.showArrays(dpixels_pb, width, height,"dpixels_pb");          
        	  SDFA_INSTANCE.showArrays(dpixels_pr_dark, width, height,"dpixels_pr_dark");          
        	  SDFA_INSTANCE.showArrays(dpixels_pb_dark, width, height,"dpixels_pb_dark");          
          }
          double mp;
          double k =1.0/(colorProcParameters.maskMax-colorProcParameters.maskMin);
          for (i=0;i<dmask.length;i++) {
        	  mp=dmask[i];
        	  if (mp < colorProcParameters.maskMin) {
        		  dmask[i]=0.0;
        	  } else if (mp< colorProcParameters.maskMax) {
        		  dmask[i]= k*(mp-colorProcParameters.maskMin);
        	  } else dmask[i]=1.0;
          }
//TODO: null DENOISE_MASK if it is not calculated
          if (colorProcParameters.combineWithSharpnessMask) {
        	  if (DENOISE_MASK==null) {
                  System.out.println ( "Can not combine masks as DENOISE_MASK is null (i.e. no denoise was performed)");
        	  } else if (DENOISE_MASK.length!=dmask.length) {
                  System.out.println ( "Can not combine masks as DENOISE_MASK length is different from that of dmask");
        	  } else {
                  for (i=0;i<dmask.length;i++) {
                	  dmask[i]+=DENOISE_MASK[i];
                	  if (dmask[i]>1.0) dmask[i]=1.0;
                  }
        	  }
        	  
          }
          for (i=0;i<dmask.length;i++) {
        	  mp=dmask[i];
      		  dpixels_pb[i]= (1.0-mp)*dpixels_pb_dark[i]+ mp* dpixels_pb[i];
      		  dpixels_pr[i]= (1.0-mp)*dpixels_pr_dark[i]+ mp* dpixels_pr[i];
          }
          DENOISE_MASK_CHROMA=dmask; // (global, used to return denoise mask to save/show
          DENOISE_MASK_CHROMA_WIDTH=width; // width of the DENOISE_MASK_CHROMA image
      } else {
    	  gb.blurDouble(dpixels_pr, width, height, colorProcParameters.chromaBrightSigma, colorProcParameters.chromaBrightSigma, 0.01);
    	  gb.blurDouble(dpixels_pb, width, height, colorProcParameters.chromaBrightSigma, colorProcParameters.chromaBrightSigma, 0.01);
   	      DENOISE_MASK_CHROMA=null; // (global, used to return denoise mask to save/show
      }
      for (i=0;i<dpixels_pr.length;i++) {
    	  fpixels_pr[i]=(float) dpixels_pr[i];
    	  fpixels_pb[i]=(float) dpixels_pb[i];
      }
      stack.addSlice("Pr", fpixels_pr);
      stack.addSlice("Pb", fpixels_pb);
      stack.addSlice("Y",  fpixels_y);
      stack.addSlice("Y0", fpixels_y0); // not filtered by low-pass, preliminary (for comaprison only)
      if (DEBUG_LEVEL>2) {
        stack.addSlice("Yr",fpixels_yR);
        stack.addSlice("Yg",fpixels_yG);
        stack.addSlice("Yb",fpixels_yB);
      }

  }
/* ======================================================================== */
  public double linGamma(double gamma, double a, double k, double x0, double x) {
    if (x<0) return 0.0;

    if (x<=x0) return k*x;
    return (1.0+a)*Math.pow(x,gamma)-a;
//  return x;
  }
/* ======================================================================== */
/* ======================================================================== */

  
/* === Parameter classes === */
  public static class ProcessParameters {
	    public int numEyesisChannels=3;
	    public int numEyesisSubChannels=3;
	    public boolean eyesisMode;
		public boolean [][] frames=new boolean[3][3];
		public boolean selectFile;
		public boolean thisFileOnly;
		public int     subChannelToProcess;
		public boolean split;
		public boolean debayer;
		public boolean showDebayerEnergy;
		public boolean saveDebayerEnergy;
		public boolean deconvolve;
		public boolean combine;
		public boolean showDenoiseMask;
		public boolean saveDenoiseMask;
		public boolean showChromaDenoiseMask;
		public boolean saveChromaDenoiseMask;
		public boolean showNoiseGains;
		public boolean saveNoiseGains;
		public boolean colorProc;
		public boolean toRGB;
		public boolean rotate;
		public boolean crop;   // crop to the sennor size 
		public boolean jpeg;   // convert to RGB and save jpeg (if save is true)
		public boolean save;
		public boolean save16; // save 16-bit tiff also if the end result is 8 bit 
		public boolean save32; // save 32-bit tiff also if the end result is 8 or 16 bit
		public boolean show;
		public int     JPEG_quality;
		public double  JPEG_scale;
		public boolean saveSettings;

		public ProcessParameters(
			boolean eyesisMode,	
			boolean frames_11,
			boolean frames_12,
			boolean frames_13,
			boolean frames_21,
			boolean frames_22,
			boolean frames_23,
			boolean frames_31,
			boolean frames_32,
			boolean frames_33,
			boolean selectFile, // ask for file(s) to process
			boolean thisFileOnly,
			int     subChannelToProcess,
			boolean split,
			boolean debayer,
			boolean showDebayerEnergy,
			boolean saveDebayerEnergy,
			boolean deconvolve,
			boolean combine,
			boolean showDenoiseMask,
			boolean saveDenoiseMask,
			boolean showChromaDenoiseMask,
			boolean saveChromaDenoiseMask,
			boolean showNoiseGains,
			boolean saveNoiseGains,
			boolean colorProc,
			boolean toRGB,
			boolean rotate,
			boolean crop,   // crop to the sennor size 
			boolean jpeg,   // convert to RGB and save jpeg (if save is true)
			boolean save,
			boolean save16, // save 16-bit tiff also if the end result is 8 bit 
			boolean save32, // save 32-bit tiff also if the end result is 8 or 16 bit
			boolean show,
			int     JPEG_quality,
			double  JPEG_scale,
			boolean saveSettings
		) {
			this.eyesisMode=eyesisMode;
			this.frames[0][0]=frames_11;
			this.frames[0][1]=frames_12;
			this.frames[0][2]=frames_13;
			this.frames[1][0]=frames_21;
			this.frames[1][1]=frames_22;
			this.frames[1][2]=frames_23;
			this.frames[2][0]=frames_31;
			this.frames[2][1]=frames_32;
			this.frames[2][2]=frames_33;
			this.selectFile=selectFile;
			this.thisFileOnly=thisFileOnly;
			this.subChannelToProcess=subChannelToProcess;
			this.split=split;
			this.debayer=debayer;
			this.showDebayerEnergy=showDebayerEnergy;
			this.saveDebayerEnergy=saveDebayerEnergy;
			this.deconvolve=deconvolve;
			this.combine=combine;
			this.showDenoiseMask=showDenoiseMask;
			this.saveDenoiseMask=saveDenoiseMask;
			this.showNoiseGains=showNoiseGains;
			this.saveNoiseGains=saveNoiseGains;
			this.showChromaDenoiseMask=showChromaDenoiseMask;
			this.saveChromaDenoiseMask=saveChromaDenoiseMask;
			this.colorProc=colorProc;
			this.toRGB=toRGB;
			this.rotate=rotate;
			this.crop=crop;
			this.jpeg=jpeg;
			this.save=save;
			this.save16=save16;
			this.save32=save32;
			this.show=show;
			this.JPEG_quality=JPEG_quality;
			this.JPEG_scale=  JPEG_scale;
			this.saveSettings=saveSettings;
		}
		public void setProperties(String prefix,Properties properties){
			int i,j;
			properties.setProperty(prefix+"numEyesisChannels",this.numEyesisChannels+"");
			properties.setProperty(prefix+"numEyesisSubChannels",this.numEyesisSubChannels+"");
			properties.setProperty(prefix+"eyesisMode",this.eyesisMode+"");
		    for (i=0;i<this.frames.length;i++) for (j=0;j<this.frames[0].length;j++)
				properties.setProperty(prefix+"frames_"+i+"_"+j,this.frames[i][j]+"");
			properties.setProperty(prefix+"selectFile",this.selectFile+"");
			properties.setProperty(prefix+"thisFileOnly",this.thisFileOnly+"");
			properties.setProperty(prefix+"subChannelToProcess",this.subChannelToProcess+"");
			properties.setProperty(prefix+"split",this.split+"");
			properties.setProperty(prefix+"debayer",this.debayer+"");
			properties.setProperty(prefix+"showDebayerEnergy",this.showDebayerEnergy+"");
			properties.setProperty(prefix+"saveDebayerEnergy",this.saveDebayerEnergy+"");
			properties.setProperty(prefix+"deconvolve",this.deconvolve+"");
			properties.setProperty(prefix+"combine",this.combine+"");
			properties.setProperty(prefix+"showDenoiseMask",this.showDenoiseMask+"");
			properties.setProperty(prefix+"saveDenoiseMask",this.saveDenoiseMask+"");
			properties.setProperty(prefix+"showChromaDenoiseMask",this.showChromaDenoiseMask+"");
			properties.setProperty(prefix+"saveChromaDenoiseMask",this.saveChromaDenoiseMask+"");
			properties.setProperty(prefix+"showNoiseGains",this.showNoiseGains+"");
			properties.setProperty(prefix+"saveNoiseGains",this.saveNoiseGains+"");
			properties.setProperty(prefix+"colorProc",this.colorProc+"");
			properties.setProperty(prefix+"toRGB",this.toRGB+"");
			properties.setProperty(prefix+"rotate",this.rotate+"");
			properties.setProperty(prefix+"crop",this.crop+"");
			properties.setProperty(prefix+"jpeg",this.jpeg+"");
			properties.setProperty(prefix+"save",this.save+"");
			properties.setProperty(prefix+"save16",this.save16+"");
			properties.setProperty(prefix+"save32",this.save32+"");
			properties.setProperty(prefix+"show",this.show+"");
			properties.setProperty(prefix+"JPEG_quality",this.JPEG_quality+"");
			properties.setProperty(prefix+"JPEG_scale",this.JPEG_scale+"");
			properties.setProperty(prefix+"saveSettings",this.saveSettings+"");
		}
		public void getProperties(String prefix,Properties properties){
			int i,j;
			if (properties.getProperty(prefix+"numEyesisChannels")!=null) this.numEyesisChannels=Integer.parseInt(properties.getProperty(prefix+"numEyesisChannels"));
			if (properties.getProperty(prefix+"numEyesisSubChannels")!=null) this.numEyesisSubChannels=Integer.parseInt(properties.getProperty(prefix+"numEyesisSubChannels"));
			if (properties.getProperty(prefix+"eyesisMode")!=null) this.eyesisMode=Boolean.parseBoolean(properties.getProperty(prefix+"eyesisMode"));
		    for (i=0;i<this.frames.length;i++) for (j=0;j<this.frames[0].length;j++)
		    	if (properties.getProperty(prefix+"frames_"+i+"_"+j)!=null) this.frames[i][j]=Boolean.parseBoolean(properties.getProperty(prefix+"frames_"+i+"_"+j));
		    if (properties.getProperty(prefix+"selectFile")!=null) this.selectFile=Boolean.parseBoolean(properties.getProperty(prefix+"selectFile"));
		    if (properties.getProperty(prefix+"thisFileOnly")!=null) this.thisFileOnly=Boolean.parseBoolean(properties.getProperty(prefix+"thisFileOnly"));
		    if (properties.getProperty(prefix+"subChannelToProcess")!=null) this.subChannelToProcess=Integer.parseInt(properties.getProperty(prefix+"subChannelToProcess"));
		    if (properties.getProperty(prefix+"split")!=null) this.split=Boolean.parseBoolean(properties.getProperty(prefix+"split"));
		    if (properties.getProperty(prefix+"debayer")!=null) this.debayer=Boolean.parseBoolean(properties.getProperty(prefix+"debayer"));
		    if (properties.getProperty(prefix+"showDebayerEnergy")!=null) this.showDebayerEnergy=Boolean.parseBoolean(properties.getProperty(prefix+"showDebayerEnergy"));
		    if (properties.getProperty(prefix+"saveDebayerEnergy")!=null) this.saveDebayerEnergy=Boolean.parseBoolean(properties.getProperty(prefix+"saveDebayerEnergy"));
		    if (properties.getProperty(prefix+"deconvolve")!=null) this.deconvolve=Boolean.parseBoolean(properties.getProperty(prefix+"deconvolve"));
		    if (properties.getProperty(prefix+"combine")!=null) this.combine=Boolean.parseBoolean(properties.getProperty(prefix+"combine"));
		    if (properties.getProperty(prefix+"showDenoiseMask")!=null) this.showDenoiseMask=Boolean.parseBoolean(properties.getProperty(prefix+"showDenoiseMask"));
		    if (properties.getProperty(prefix+"saveDenoiseMask")!=null) this.saveDenoiseMask=Boolean.parseBoolean(properties.getProperty(prefix+"saveDenoiseMask"));
		    if (properties.getProperty(prefix+"showChromaDenoiseMask")!=null) this.showChromaDenoiseMask=Boolean.parseBoolean(properties.getProperty(prefix+"showChromaDenoiseMask"));
		    if (properties.getProperty(prefix+"saveChromaDenoiseMask")!=null) this.saveChromaDenoiseMask=Boolean.parseBoolean(properties.getProperty(prefix+"saveChromaDenoiseMask"));
		    if (properties.getProperty(prefix+"showNoiseGains")!=null) this.showNoiseGains=Boolean.parseBoolean(properties.getProperty(prefix+"showNoiseGains"));
		    if (properties.getProperty(prefix+"saveNoiseGains")!=null) this.saveNoiseGains=Boolean.parseBoolean(properties.getProperty(prefix+"saveNoiseGains"));
		    if (properties.getProperty(prefix+"colorProc")!=null) this.colorProc=Boolean.parseBoolean(properties.getProperty(prefix+"colorProc"));
		    if (properties.getProperty(prefix+"toRGB")!=null) this.toRGB=Boolean.parseBoolean(properties.getProperty(prefix+"toRGB"));
		    if (properties.getProperty(prefix+"rotate")!=null) this.rotate=Boolean.parseBoolean(properties.getProperty(prefix+"rotate"));
		    if (properties.getProperty(prefix+"crop")!=null) this.crop=Boolean.parseBoolean(properties.getProperty(prefix+"crop"));   // crop to the sensor size 
		    if (properties.getProperty(prefix+"jpeg")!=null) this.jpeg=Boolean.parseBoolean(properties.getProperty(prefix+"jpeg"));   // convert to RGB and save jpeg (if save is true)
		    if (properties.getProperty(prefix+"save")!=null) this.save=Boolean.parseBoolean(properties.getProperty(prefix+"save"));
		    if (properties.getProperty(prefix+"save16")!=null) this.save16=Boolean.parseBoolean(properties.getProperty(prefix+"save16")); // save 16-bit tiff also if the end result is 8 bit 
		    if (properties.getProperty(prefix+"save32")!=null) this.save32=Boolean.parseBoolean(properties.getProperty(prefix+"save32")); // save 32-bit tiff also if the end result is 8 or 16 bit
		    if (properties.getProperty(prefix+"show")!=null) this.show=Boolean.parseBoolean(properties.getProperty(prefix+"show"));
		    if (properties.getProperty(prefix+"JPEG_quality")!=null) this.JPEG_quality=Integer.parseInt(properties.getProperty(prefix+"JPEG_quality"));
		    if (properties.getProperty(prefix+"JPEG_scale")!=null) this.JPEG_scale=Double.parseDouble(properties.getProperty(prefix+"JPEG_scale"));
		    if (properties.getProperty(prefix+"saveSettings")!=null) this.saveSettings=Boolean.parseBoolean(properties.getProperty(prefix+"saveSettings"));
		}
	}
  
/* ======================================================================== */
  public static class FilesParameters {
	  public String [][] rPSFNames=new String [3][3];
	  public String [][] gaussianNames=new String [3][3];
	  public String kernelDirectory;
	  public String resultsDirectory;
	  public String [] sourceFiles;
	  public boolean useXML;
	  public FilesParameters(
			  String rPSFNames_11,
			  String rPSFNames_12,
			  String rPSFNames_13,
			  String rPSFNames_21,
			  String rPSFNames_22,
			  String rPSFNames_23,
			  String rPSFNames_31,
			  String rPSFNames_32,
			  String rPSFNames_33,
			  String gaussianNames_11,
			  String gaussianNames_12,
			  String gaussianNames_13,
			  String gaussianNames_21,
			  String gaussianNames_22,
			  String gaussianNames_23,
			  String gaussianNames_31,
			  String gaussianNames_32,
			  String gaussianNames_33,
			  String kernelDirectory,
  		      String resultsDirectory,
  		      boolean useXML
  		      ){
		  this.rPSFNames[0][0]=rPSFNames_11;
		  this.rPSFNames[0][1]=rPSFNames_12;
		  this.rPSFNames[0][2]=rPSFNames_13;
		  this.rPSFNames[1][0]=rPSFNames_21;
		  this.rPSFNames[1][1]=rPSFNames_22;
		  this.rPSFNames[1][2]=rPSFNames_23;
		  this.rPSFNames[2][0]=rPSFNames_31;
		  this.rPSFNames[2][1]=rPSFNames_32;
		  this.rPSFNames[2][2]=rPSFNames_33;
		  this.gaussianNames[0][0]=gaussianNames_11;
		  this.gaussianNames[0][1]=gaussianNames_12;
		  this.gaussianNames[0][2]=gaussianNames_13;
		  this.gaussianNames[1][0]=gaussianNames_21;
		  this.gaussianNames[1][1]=gaussianNames_22;
		  this.gaussianNames[1][2]=gaussianNames_23;
		  this.gaussianNames[2][0]=gaussianNames_31;
		  this.gaussianNames[2][1]=gaussianNames_32;
		  this.gaussianNames[2][2]=gaussianNames_33;
		  this.kernelDirectory=    kernelDirectory;
		  this.resultsDirectory=   resultsDirectory;
		  this.useXML=useXML;
	  }
		public void setProperties(String prefix,Properties properties){
//			properties.setProperty(prefix+"",this.+"");
			int i,j;
			for (i=0;i<this.rPSFNames.length;i++) for (j=0;j<this.rPSFNames[i].length;j++)
				properties.setProperty(prefix+"rPSFNames_"+i+"_"+j,this.rPSFNames[i][j]);				
			for (i=0;i<this.gaussianNames.length;i++) for (j=0;j<this.gaussianNames[i].length;j++)
				properties.setProperty(prefix+"gaussianNames_"+i+"_"+j,this.gaussianNames[i][j]);
			properties.setProperty(prefix+"kernelDirectory",this.kernelDirectory);
			properties.setProperty(prefix+"resultsDirectory",this.resultsDirectory);	
			properties.setProperty(prefix+"useXML",this.useXML+"");
			j=(this.sourceFiles==null)?0:this.sourceFiles.length;
			properties.setProperty(prefix+"sourceFiles_length",j+"");
			for (i=0;i<j;i++)
				properties.setProperty(prefix+"sourceFiles_"+i,this.sourceFiles[i]);
		}
		public void getProperties(String prefix,Properties properties){
			int i,j;
			for (i=0;i<this.rPSFNames.length;i++) for (j=0;j<this.rPSFNames[i].length;j++)
				this.rPSFNames[i][j]=properties.getProperty(prefix+"rPSFNames_"+i+"_"+j);
			for (i=0;i<this.gaussianNames.length;i++) for (j=0;j<this.gaussianNames[i].length;j++)
				this.gaussianNames[i][j]=properties.getProperty(prefix+"gaussianNames_"+i+"_"+j);
			this.kernelDirectory=properties.getProperty(prefix+"kernelDirectory");
			this.resultsDirectory=properties.getProperty(prefix+"resultsDirectory");
			this.useXML=Boolean.parseBoolean(properties.getProperty(prefix+"useXML"));
			j=Integer.parseInt(properties.getProperty(prefix+"sourceFiles_length"));
			this.sourceFiles=new String[j];
			for (i=0;i<j;i++)
				this.sourceFiles[i]=properties.getProperty(prefix+"sourceFiles_"+i);
		}
  }

/* ======================================================================== */
  public static class RGBParameters {
		public double r_min = 0.075;
		public double g_min = 0.075;
		public double b_min = 0.075;
		public double r_max = 1.0;
		public double g_max = 1.0;
		public double b_max = 1.0;
		public double alpha_min = 0.0;
		public double alpha_max = 1.0;
/*
		public RGBParameters(double r_min, double g_min, double b_min, double r_max, double g_max, double b_max) {
			this.r_min = r_min;
			this.g_min = g_min;
			this.b_min = b_min;
			this.r_max = r_max;
			this.g_max = g_max;
			this.b_max = b_max;
		}
*/		
  		public RGBParameters(double r_min, double g_min, double b_min, double r_max, double g_max, double b_max, double alpha_min, double alpha_max) {
  			this.r_min = r_min;
  			this.g_min = g_min;
  			this.b_min = b_min;
  			this.r_max = r_max;
  			this.g_max = g_max;
  			this.b_max = b_max;
  			this.alpha_min = alpha_min;
  			this.alpha_max = alpha_max;
  		}
		
		public void setProperties(String prefix,Properties properties){
			properties.setProperty(prefix+"r_min",this.r_min+"");
			properties.setProperty(prefix+"g_min",this.g_min+"");
			properties.setProperty(prefix+"b_min",this.b_min+"");
			properties.setProperty(prefix+"r_max",this.r_max+"");
			properties.setProperty(prefix+"g_max",this.g_max+"");
			properties.setProperty(prefix+"b_max",this.b_max+"");
  			properties.setProperty(prefix+"alpha_min",this.alpha_min+"");
  			properties.setProperty(prefix+"alpha_max",this.alpha_max+"");
			
		}
		public void getProperties(String prefix,Properties properties){
			this.r_min=Double.parseDouble(properties.getProperty(prefix+"r_min"));
			this.g_min=Double.parseDouble(properties.getProperty(prefix+"g_min"));
			this.b_min=Double.parseDouble(properties.getProperty(prefix+"b_min"));
			this.r_max=Double.parseDouble(properties.getProperty(prefix+"r_max"));
			this.g_max=Double.parseDouble(properties.getProperty(prefix+"g_max"));
			this.b_max=Double.parseDouble(properties.getProperty(prefix+"b_max"));
  			if (properties.getProperty(prefix+"alpha_min")!=null)  this.alpha_min=Double.parseDouble(properties.getProperty(prefix+"alpha_min"));
  			if (properties.getProperty(prefix+"alpha_max")!=null)  this.alpha_max=Double.parseDouble(properties.getProperty(prefix+"alpha_max"));
		}
		
	}
/* ======================================================================== */
  
  public static class ColorProcParameters {
		public double balanceRed;
		public double balanceBlue;
		public double gain;
		public double weightScaleR;
		public double weightScaleB;
//		public double sigma;
		public double gamma;
		public double minLin;
		public double kr;
		public double kb;
		public double saturationRed;
		public double saturationBlue;
		public boolean useFirstY;
		public double maskSigma; // LPF for luma to calculate chroma mask
		public double maskMin; // minimal level for the mask (absolute luma values)
		public double maskMax; // maximal level for the mask (absolute luma values)
		public boolean combineWithSharpnessMask; // combine chroma mask with sharpness mask to reduce color leak through borders
		public double chromaBrightSigma; // LPF for chroma in the bright areas (determined by the mask)
		public double chromaDarkSigma;   // LPF for chroma in the dark areas (determined by the mask)
		
		public ColorProcParameters(
				double balanceRed,
				double balanceBlue,
				double gain,
				double weightScaleR,
				double weightScaleB,
//				double sigma,
				double gamma,
				double minLin,
				double kr,
				double kb,
				double saturationRed,
				double saturationBlue,
				boolean useFirstY,
				double maskSigma, // LPF for luma to calculate chroma mask
				double maskMin, // minimal level for the mask (absolute luma values)
				double maskMax, // maximal level for the mask (absolute luma values)
				boolean combineWithSharpnessMask, // combine chroma mask with sharpness mask to reduce color leak through borders
				double chromaBrightSigma, // LPF for chroma in the bright areas (determined by the mask)
				double chromaDarkSigma   // LPF for chroma in the dark areas (determined by the mask)
             ) {
			this.balanceRed = balanceRed;
			this.balanceBlue = balanceBlue;
			this.gain = gain;
			this.weightScaleR = weightScaleR;
			this.weightScaleB = weightScaleB;
//			this.sigma = sigma;
			this.gamma = gamma;
			this.minLin = minLin;
			this.kr = kr;
			this.kb = kb;
			this.saturationRed = saturationRed;
			this.saturationBlue = saturationBlue;
			this.useFirstY=useFirstY;
			this.maskSigma=maskSigma;
			this.maskMin=maskMin;
			this.maskMax=maskMax;
			this.combineWithSharpnessMask=combineWithSharpnessMask;
			this.chromaBrightSigma=chromaBrightSigma;
			this.chromaDarkSigma=chromaDarkSigma;
		}
		public void setProperties(String prefix,Properties properties){
			properties.setProperty(prefix+"balanceRed",this.balanceRed+"");
			properties.setProperty(prefix+"balanceBlue",this.balanceBlue+"");
			properties.setProperty(prefix+"gain",this.gain+"");
			properties.setProperty(prefix+"weightScaleR",this.weightScaleR+"");
			properties.setProperty(prefix+"weightScaleB",this.weightScaleB+"");
//			properties.setProperty(prefix+"sigma",this.sigma+"");
			properties.setProperty(prefix+"gamma",this.gamma+"");
			properties.setProperty(prefix+"minLin",this.minLin+"");
			properties.setProperty(prefix+"kr",this.kr+"");
			properties.setProperty(prefix+"kb",this.kb+"");
			properties.setProperty(prefix+"saturationRed",this.saturationRed+"");
			properties.setProperty(prefix+"saturationBlue",this.saturationBlue+"");
			properties.setProperty(prefix+"useFirstY",this.useFirstY+"");
			properties.setProperty(prefix+"maskSigma",this.maskSigma+"");
			properties.setProperty(prefix+"maskMin",this.maskMin+"");
			properties.setProperty(prefix+"maskMax",this.maskMax+"");
			properties.setProperty(prefix+"combineWithSharpnessMask",this.combineWithSharpnessMask+"");
			properties.setProperty(prefix+"chromaBrightSigma",this.chromaBrightSigma+"");
			properties.setProperty(prefix+"chromaDarkSigma",this.chromaDarkSigma+"");
		}
		public void getProperties(String prefix,Properties properties){
			if (properties.getProperty(prefix+"balanceRed")!=null) this.balanceRed=Double.parseDouble(properties.getProperty(prefix+"balanceRed"));
			if (properties.getProperty(prefix+"balanceBlue")!=null) this.balanceBlue=Double.parseDouble(properties.getProperty(prefix+"balanceBlue"));
			if (properties.getProperty(prefix+"gain")!=null) this.gain=Double.parseDouble(properties.getProperty(prefix+"gain"));
			if (properties.getProperty(prefix+"weightScaleR")!=null) this.weightScaleR=Double.parseDouble(properties.getProperty(prefix+"weightScaleR"));
			if (properties.getProperty(prefix+"weightScaleB")!=null) this.weightScaleB=Double.parseDouble(properties.getProperty(prefix+"weightScaleB"));
//			this.sigma=Double.parseDouble(properties.getProperty(prefix+"sigma"));
			if (properties.getProperty(prefix+"gamma")!=null) this.gamma=Double.parseDouble(properties.getProperty(prefix+"gamma"));
			if (properties.getProperty(prefix+"minLin")!=null) this.minLin=Double.parseDouble(properties.getProperty(prefix+"minLin"));
			if (properties.getProperty(prefix+"kr")!=null) this.kr=Double.parseDouble(properties.getProperty(prefix+"kr"));
			if (properties.getProperty(prefix+"kb")!=null) this.kb=Double.parseDouble(properties.getProperty(prefix+"kb"));
			if (properties.getProperty(prefix+"saturationRed")!=null) this.saturationRed=Double.parseDouble(properties.getProperty(prefix+"saturationRed"));
			if (properties.getProperty(prefix+"saturationBlue")!=null) this.saturationBlue=Double.parseDouble(properties.getProperty(prefix+"saturationBlue"));
			if (properties.getProperty(prefix+"useFirstY")!=null) this.useFirstY=Boolean.parseBoolean(properties.getProperty(prefix+"useFirstY"));
			if (properties.getProperty(prefix+"maskSigma")!=null) this.maskSigma=Double.parseDouble(properties.getProperty(prefix+"maskSigma"));
			if (properties.getProperty(prefix+"maskMin")!=null) this.maskMin=Double.parseDouble(properties.getProperty(prefix+"maskMin"));
			if (properties.getProperty(prefix+"maskMax")!=null) this.maskMax=Double.parseDouble(properties.getProperty(prefix+"maskMax"));
			if (properties.getProperty(prefix+"combineWithSharpnessMask")!=null) this.combineWithSharpnessMask=Boolean.parseBoolean(properties.getProperty(prefix+"combineWithSharpnessMask"));
			if (properties.getProperty(prefix+"chromaBrightSigma")!=null) this.chromaBrightSigma=Double.parseDouble(properties.getProperty(prefix+"chromaBrightSigma"));
			if (properties.getProperty(prefix+"chromaDarkSigma")!=null) this.chromaDarkSigma=Double.parseDouble(properties.getProperty(prefix+"chromaDarkSigma"));
		}
	}
  /* ======================================================================== */
// individual per-channel color balance and gain
  public static class ColorCalibParameters {
		public double[][] gain=new double[3][3];
		public double[][] balanceRed=new double[3][3];
		public double[][] balanceBlue=new double[3][3];
	  	public ColorCalibParameters(
	  			double gain_11,
	  			double gain_12,
	  			double gain_13,
	  			double gain_21,
	  			double gain_22,
	  			double gain_23,
	  			double gain_31,
	  			double gain_32,
	  			double gain_33,
	  			double balanceRed_11,
	  			double balanceRed_12,
	  			double balanceRed_13,
	  			double balanceRed_21,
	  			double balanceRed_22,
	  			double balanceRed_23,
	  			double balanceRed_31,
	  			double balanceRed_32,
	  			double balanceRed_33,
	  			double balanceBlue_11,
	  			double balanceBlue_12,
	  			double balanceBlue_13,
	  			double balanceBlue_21,
	  			double balanceBlue_22,
	  			double balanceBlue_23,
	  			double balanceBlue_31,
	  			double balanceBlue_32,
	  			double balanceBlue_33){
	  		this.gain[0][0]=gain_11;
	  		this.gain[0][1]=gain_12;
	  		this.gain[0][2]=gain_13;
	  		this.gain[1][0]=gain_21;
	  		this.gain[1][1]=gain_22;
	  		this.gain[1][2]=gain_23;
	  		this.gain[2][0]=gain_31;
	  		this.gain[2][1]=gain_32;
	  		this.gain[2][2]=gain_33;
	  		this.balanceRed[0][0]=balanceRed_11;
	  		this.balanceRed[0][1]=balanceRed_12;
	  		this.balanceRed[0][2]=balanceRed_13;
	  		this.balanceRed[1][0]=balanceRed_21;
	  		this.balanceRed[1][1]=balanceRed_22;
	  		this.balanceRed[1][2]=balanceRed_23;
	  		this.balanceRed[2][0]=balanceRed_31;
	  		this.balanceRed[2][1]=balanceRed_32;
	  		this.balanceRed[2][2]=balanceRed_33;
	  		this.balanceBlue[0][0]=balanceBlue_11;
	  		this.balanceBlue[0][1]=balanceBlue_12;
	  		this.balanceBlue[0][2]=balanceBlue_13;
	  		this.balanceBlue[1][0]=balanceBlue_21;
	  		this.balanceBlue[1][1]=balanceBlue_22;
	  		this.balanceBlue[1][2]=balanceBlue_23;
	  		this.balanceBlue[2][0]=balanceBlue_31;
	  		this.balanceBlue[2][1]=balanceBlue_32;
	  		this.balanceBlue[2][2]=balanceBlue_33;
	    }
		public void setProperties(String prefix,Properties properties){
			int i,j;
			for (i=0;i<this.gain.length;i++) for (j=0;j<this.gain[i].length;j++)
			  properties.setProperty(prefix+"gain_"+i+"_"+j,this.gain[i][j]+"");
			for (i=0;i<this.balanceRed.length;i++) for (j=0;j<this.balanceRed[i].length;j++)
				  properties.setProperty(prefix+"balanceRed_"+i+"_"+j,this.balanceRed[i][j]+"");
			for (i=0;i<this.balanceBlue.length;i++) for (j=0;j<this.balanceBlue[i].length;j++)
				  properties.setProperty(prefix+"balanceBlue_"+i+"_"+j,this.balanceBlue[i][j]+"");
		}
		public void getProperties(String prefix,Properties properties){
//			this.oversample=Integer.parseInt(properties.getProperty(prefix+"oversample"));
			int i,j;
			String s;
			for (i=0;i<this.gain.length;i++) for (j=0;j<this.gain[i].length;j++) {
				s=properties.getProperty(prefix+"gain_"+i+"_"+j);
				if (s!=null) this.gain[i][j]=Double.parseDouble(s);
			}
			for (i=0;i<this.balanceRed.length;i++) for (j=0;j<this.balanceRed[i].length;j++) {
				s=properties.getProperty(prefix+"balanceRed_"+i+"_"+j);
				if (s!=null) this.balanceRed[i][j]=Double.parseDouble(s);
			}
			for (i=0;i<this.balanceBlue.length;i++) for (j=0;j<this.balanceBlue[i].length;j++) {
				s=properties.getProperty(prefix+"balanceBlue_"+i+"_"+j);
				if (s!=null) this.balanceBlue[i][j]=Double.parseDouble(s);
			}
		}

  }
  /* ======================================================================== */
  public static class NonlinParameters {
  	public boolean useRejectBlocksFilter;
  	public boolean combineBothModes; 
 	public int    maskFFTSize; // 256
 	public int    blockPeriod; // 32
  	public double rejectFreqSigma; // 1.0, frequency domain
  	public double lowPassSigma;    // 5.0, spatial domain
  	public double filtMin;
  	public double filtMax;
	public double [][] thresholdCorrection=new double[3][3]; // apply to filtMin and filtMax
  	public double threshold;
	public boolean useDiffNoiseGains;
	public double [] noiseGainWeights=new double[3]; 
	double blurSigma;     // blur sigma for mask calculation (blur convolution kernels for noise gain calculation 
	public double  noiseGainPower;
  	public boolean showMask;
// ring filter
  	public boolean useRingFilter;    // filter out spots on denoise mask
    public double minMaxValue;       // minimal value (relative to filtMax) of the local maximum to be processed
    public double overRingThreshold; // ratio of local max. and maximal value in the surrounding ring to trigger filter 
    public double overRingLimit;     // limit values in the center circle to scaled maximum in a ring
    public double ringIR;            // ring inner radius (center circle radius)
    public double ringOR;            // ring outer radius
    
  	public NonlinParameters(
  			boolean useRejectBlocksFilter,
  			boolean combineBothModes,
  			int maskFFTSize,
  			int blockPeriod,
  			double rejectFreqSigma,
  			double lowPassSigma,
  			double filtMin,
  			double filtMax,
  			double thresholdCorrection_11,
  			double thresholdCorrection_12,
  			double thresholdCorrection_13,
  			double thresholdCorrection_21,
  			double thresholdCorrection_22,
  			double thresholdCorrection_23,
  			double thresholdCorrection_31,
  			double thresholdCorrection_32,
  			double thresholdCorrection_33,
  			double threshold,
			boolean useDiffNoiseGains,
			double noiseGainWeights_0, // r
			double noiseGainWeights_1, // b
			double noiseGainWeights_2, // g
			double blurSigma,     // blur sigma for mask calculation (blur convolution kernels for noise gain calculation 
			double noiseGainPower,
			boolean useRingFilter,    // filter out spots on denoise mask
  		    double minMaxValue,       // minimal value (relative to filtMax) of the local maximum to be processed
  		    double overRingThreshold, // ratio of local max. and maximal value in the surrounding ring to trigger filter 
  		    double overRingLimit,     // limit values in the center circle to scaled maximum in a ring
  		    double ringIR,            // ring inner radius (center circle radius)
  		    double ringOR             // ring outer radius
			
) {
  		this.useRejectBlocksFilter=useRejectBlocksFilter;
  		this.combineBothModes=combineBothModes;
  		this.maskFFTSize = maskFFTSize;
  		this.blockPeriod = blockPeriod;
  		this.rejectFreqSigma = rejectFreqSigma;
  		this.lowPassSigma = lowPassSigma;
  		this.filtMin = filtMin;
  		this.filtMax = filtMax;
  		this.thresholdCorrection[0][0]=thresholdCorrection_11;
  		this.thresholdCorrection[0][1]=thresholdCorrection_12;
  		this.thresholdCorrection[0][2]=thresholdCorrection_13;
  		this.thresholdCorrection[1][0]=thresholdCorrection_21;
  		this.thresholdCorrection[1][1]=thresholdCorrection_22;
  		this.thresholdCorrection[1][2]=thresholdCorrection_23;
  		this.thresholdCorrection[2][0]=thresholdCorrection_31;
  		this.thresholdCorrection[2][1]=thresholdCorrection_32;
  		this.thresholdCorrection[2][2]=thresholdCorrection_33;
  		this.threshold = threshold;
		this.useDiffNoiseGains=useDiffNoiseGains;
		this.noiseGainWeights[0]=noiseGainWeights_0;
		this.noiseGainWeights[1]=noiseGainWeights_1;
		this.noiseGainWeights[2]=noiseGainWeights_2;
		this.blurSigma=blurSigma;
		this.noiseGainPower=noiseGainPower;
		this.useRingFilter=useRingFilter;
		this.minMaxValue=minMaxValue;
		this.overRingThreshold=overRingThreshold; 
		this.overRingLimit=overRingLimit;
		this.ringIR=ringIR;
		this.ringOR=ringOR;
		
  	}
	public void setProperties(String prefix,Properties properties){
		int i,j;
//		properties.setProperty(prefix+"oversample",this.oversample+"");
		properties.setProperty(prefix+"useRejectBlocksFilter",this.useRejectBlocksFilter+"");
		properties.setProperty(prefix+"combineBothModes",this.combineBothModes+"");
		properties.setProperty(prefix+"maskFFTSize",this.maskFFTSize+"");
		properties.setProperty(prefix+"blockPeriod",this.blockPeriod+"");
		properties.setProperty(prefix+"rejectFreqSigma",this.rejectFreqSigma+"");
		properties.setProperty(prefix+"lowPassSigma",this.lowPassSigma+"");
		properties.setProperty(prefix+"filtMin",this.filtMin+"");
		properties.setProperty(prefix+"filtMax",this.filtMax+"");
		for (i=0;i<this.thresholdCorrection.length;i++) for (j=0;j<this.thresholdCorrection[i].length;j++)
		  properties.setProperty(prefix+"thresholdCorrection_"+i+"_"+j,this.thresholdCorrection[i][j]+"");
		properties.setProperty(prefix+"threshold",this.threshold+"");
		properties.setProperty(prefix+"useDiffNoiseGains",this.useDiffNoiseGains+"");
		for (i=0;i<this.noiseGainWeights.length;i++)
		   properties.setProperty(prefix+"noiseGainWeights_"+i,this.noiseGainWeights[i]+"");
		properties.setProperty(prefix+"blurSigma",this.blurSigma+"");
		properties.setProperty(prefix+"noiseGainPower",this.noiseGainPower+"");
		properties.setProperty(prefix+"useRingFilter",this.useRingFilter+"");
		properties.setProperty(prefix+"minMaxValue",this.minMaxValue+"");
		properties.setProperty(prefix+"overRingThreshold",this.overRingThreshold+""); 
		properties.setProperty(prefix+"overRingLimit",this.overRingLimit+"");
		properties.setProperty(prefix+"ringIR",this.ringIR+"");
		properties.setProperty(prefix+"ringOR",this.ringOR+"");
	}
	public void getProperties(String prefix,Properties properties){
//		this.oversample=Integer.parseInt(properties.getProperty(prefix+"oversample"));
		int i,j;
		String s;
		this.useRejectBlocksFilter=Boolean.parseBoolean(properties.getProperty(prefix+"useRejectBlocksFilter"));
		this.combineBothModes=Boolean.parseBoolean(properties.getProperty(prefix+"combineBothModes")); 
		this.maskFFTSize=Integer.parseInt(properties.getProperty(prefix+"maskFFTSize"));
		this.blockPeriod=Integer.parseInt(properties.getProperty(prefix+"blockPeriod"));
		this.rejectFreqSigma=Double.parseDouble(properties.getProperty(prefix+"rejectFreqSigma"));
		this.lowPassSigma=Double.parseDouble(properties.getProperty(prefix+"lowPassSigma"));
		this.filtMin=Double.parseDouble(properties.getProperty(prefix+"filtMin"));
		this.filtMax=Double.parseDouble(properties.getProperty(prefix+"filtMax"));
		for (i=0;i<this.thresholdCorrection.length;i++) for (j=0;j<this.thresholdCorrection[i].length;j++)
		    this.thresholdCorrection[i][j]=Double.parseDouble(properties.getProperty(prefix+"thresholdCorrection_"+i+"_"+j));
		this.threshold=Double.parseDouble(properties.getProperty(prefix+"threshold"));
		this.useDiffNoiseGains=Boolean.parseBoolean(properties.getProperty(prefix+"useDiffNoiseGains"));
		for (i=0;i<this.noiseGainWeights.length;i++)
			this.noiseGainWeights[i]=Double.parseDouble(properties.getProperty(prefix+"noiseGainWeights_"+i)); 
		this.blurSigma=Double.parseDouble(properties.getProperty(prefix+"blurSigma"));
		this.noiseGainPower=Double.parseDouble(properties.getProperty(prefix+"noiseGainPower"));
		s=properties.getProperty(prefix+"useRingFilter");
		if ((s==null) || (s=="")) return; // earlier revision
		this.useRingFilter=Boolean.parseBoolean(properties.getProperty(prefix+"useRingFilter"));
		this.minMaxValue=Double.parseDouble(properties.getProperty(prefix+"minMaxValue"));
		this.overRingThreshold=Double.parseDouble(properties.getProperty(prefix+"overRingThreshold")); 
		this.overRingLimit=Double.parseDouble(properties.getProperty(prefix+"overRingLimit"));
		this.ringIR=Double.parseDouble(properties.getProperty(prefix+"ringIR"));
		this.ringOR=Double.parseDouble(properties.getProperty(prefix+"ringOR"));
	}
  	
  }
/* ======================================================================== */
  public static class SplitParameters {
		public int oversample;
		public int addLeft;
		public int addTop;
		public int addRight;
		public int addBottom;

		public SplitParameters(int oversample, int addLeft, int addTop,
				int addRight, int addBottom) {
			this.oversample = oversample;
			this.addLeft = addLeft;
			this.addTop = addTop;
			this.addRight = addRight;
			this.addBottom = addBottom;
		}
		public void setProperties(String prefix,Properties properties){
			properties.setProperty(prefix+"oversample",this.oversample+"");
			properties.setProperty(prefix+"addLeft",   this.addLeft+"");
			properties.setProperty(prefix+"addTop",    this.addTop+"");
			properties.setProperty(prefix+"addRight",  this.addRight+"");
			properties.setProperty(prefix+"addBottom", this.addBottom+"");
		}
		public void getProperties(String prefix,Properties properties){
			this.oversample=Integer.parseInt(properties.getProperty(prefix+"oversample"));
			this.addLeft=Integer.parseInt(properties.getProperty(prefix+"addLeft"));
			this.addTop=Integer.parseInt(properties.getProperty(prefix+"addTop"));
			this.addRight=Integer.parseInt(properties.getProperty(prefix+"addRight"));
			this.addBottom=Integer.parseInt(properties.getProperty(prefix+"addBottom"));
		}
	}
/* ======================================================================== */
  public static class DebayerParameters {
		public int size;
		public double polarStep;
		public double debayerThreshold;
		public double debayerRelativeWidthGreen;
		public double debayerRelativeWidthRedblue;
		public double debayerRelativeWidthRedblueMain;
		public double debayerRelativeWidthRedblueClones;
		public double debayerGamma;
		public double debayerBonus;
		public double mainToAlias;
		public double debayerMaskBlur;
		public boolean debayerUseScissors;
		public boolean debug;
		public int xDebug;
		public int yDebug;
		public boolean debayerStacks;
		public DebayerParameters(int size, double polarStep,
				double debayerThreshold, double debayerRelativeWidthGreen,
				double debayerRelativeWidthRedblue,
				double debayerRelativeWidthRedblueMain,
				double debayerRelativeWidthRedblueClones, double debayerGamma,
				double debayerBonus, double mainToAlias, double debayerMaskBlur,
				boolean debayerUseScissors, 
				boolean debug, int xDebug, int yDebug,
				boolean debayerStacks) {
			this.size = size;
			this.polarStep = polarStep;
			this.debayerThreshold = debayerThreshold;
			this.debayerRelativeWidthGreen = debayerRelativeWidthGreen;
			this.debayerRelativeWidthRedblue = debayerRelativeWidthRedblue;
			this.debayerRelativeWidthRedblueMain = debayerRelativeWidthRedblueMain;
			this.debayerRelativeWidthRedblueClones = debayerRelativeWidthRedblueClones;
			this.debayerGamma = debayerGamma;
			this.debayerBonus = debayerBonus;
			this.mainToAlias = mainToAlias;
			this.debayerMaskBlur = debayerMaskBlur;
			this.debayerUseScissors = debayerUseScissors;
			this.debug = debug;
			this.xDebug = xDebug;
			this.yDebug = yDebug;
			this.debayerStacks = debayerStacks;
		}
		public void setProperties(String prefix,Properties properties){
//			properties.setProperty(prefix+"oversample",this.oversample+"");
			properties.setProperty(prefix+"size",this.size+"");
			properties.setProperty(prefix+"polarStep",this.polarStep+"");
			properties.setProperty(prefix+"debayerThreshold",this.debayerThreshold+"");
			properties.setProperty(prefix+"debayerRelativeWidthGreen",this.debayerRelativeWidthGreen+"");
			properties.setProperty(prefix+"debayerRelativeWidthRedblue",this.debayerRelativeWidthRedblue+"");
			properties.setProperty(prefix+"debayerRelativeWidthRedblueMain",this.debayerRelativeWidthRedblueMain+"");
			properties.setProperty(prefix+"debayerRelativeWidthRedblueClones",this.debayerRelativeWidthRedblueClones+"");
			properties.setProperty(prefix+"debayerGamma",this.debayerGamma+"");
			properties.setProperty(prefix+"debayerBonus",this.debayerBonus+"");
			properties.setProperty(prefix+"mainToAlias",this.mainToAlias+"");
			properties.setProperty(prefix+"debayerMaskBlur",this.debayerMaskBlur+"");
			properties.setProperty(prefix+"debayerUseScissors",this.debayerUseScissors+"");
			properties.setProperty(prefix+"debug",this.debug+"");
			properties.setProperty(prefix+"xDebug",this.xDebug+"");
			properties.setProperty(prefix+"yDebug",this.yDebug+"");
			properties.setProperty(prefix+"debayerStacks",this.debayerStacks+"");
		}
		public void getProperties(String prefix,Properties properties){
			this.size=                             Integer.parseInt(properties.getProperty(prefix+"size"));
			this.polarStep=                        Double.parseDouble(properties.getProperty(prefix+"polarStep"));
			this.debayerThreshold=                 Double.parseDouble(properties.getProperty(prefix+"debayerThreshold"));
			this.debayerRelativeWidthGreen=        Double.parseDouble(properties.getProperty(prefix+"debayerRelativeWidthGreen"));
			this.debayerRelativeWidthRedblue=      Double.parseDouble(properties.getProperty(prefix+"debayerRelativeWidthRedblue"));
			this.debayerRelativeWidthRedblueMain=  Double.parseDouble(properties.getProperty(prefix+"debayerRelativeWidthRedblueMain"));
			this.debayerRelativeWidthRedblueClones=Double.parseDouble(properties.getProperty(prefix+"debayerRelativeWidthRedblueClones"));
			this.debayerGamma=                     Double.parseDouble(properties.getProperty(prefix+"debayerGamma"));
			this.debayerBonus=                     Double.parseDouble(properties.getProperty(prefix+"debayerBonus"));
			this.mainToAlias=                      Double.parseDouble(properties.getProperty(prefix+"mainToAlias"));
			this.debayerMaskBlur=                  Double.parseDouble(properties.getProperty(prefix+"debayerMaskBlur"));
			this.debayerUseScissors=               Boolean.parseBoolean(properties.getProperty(prefix+"debayerUseScissors"));
			this.debug=                            Boolean.parseBoolean(properties.getProperty(prefix+"debug"));
			this.xDebug=                           Integer.parseInt(properties.getProperty(prefix+"xDebug"));
			this.yDebug=                           Integer.parseInt(properties.getProperty(prefix+"yDebug"));
			this.debayerStacks=                    Boolean.parseBoolean(properties.getProperty(prefix+"debayerStacks"));
		}
	}

/* ======================================================================== */
 	/**
	 * Main method for debugging.
	 *
	 * For debugging, it is convenient to have a method that starts ImageJ, loads an
	 * image and calls the plugin, e.g. after setting breakpoints.
	 * Grabbed from https://github.com/imagej/minimal-ij1-plugin
	 * @param args unused
	 */
	public static void main(String[] args) {
		// set the plugins.dir property to make the plugin appear in the Plugins menu
		Class<?> clazz = Aberration_Calibration.class;
		String url = clazz.getResource("/" + clazz.getName().replace('.', '/') + ".class").toString();
		String pluginsDir = url.substring(5, url.length() - clazz.getName().length() - 6);
		System.setProperty("plugins.dir", pluginsDir);
		// start ImageJ
		new ImageJ();
		// run the plugin
		IJ.runPlugIn(clazz.getName(), "");
	}

}

